//! Interval-based annotation features.

extern crate bio;

use std::cmp::{max, min};
use std::collections::HashMap;
use std::iter;
use std::ops::Deref;

use bio::utils::{self, Interval, IntervalError, Strand};
use sliding_windows::{IterExt, Storage};

use self::error::FeatureError;


pub mod error {

    use super::*;

    quick_error! {
        #[derive(Debug)]
        pub enum FeatureError {
            IntervalError(err: utils::IntervalError) {
                description(
                    match err {
                        &IntervalError::InvalidRange =>
                            "interval start coordinate larger than its end coordinate",
                        ref otherwise => otherwise.description(),
                    })
                from()
            }
            StrandCharError(err: utils::StrandError) {
                description(err.description())
                from()
            }
            ConflictingStrandError {
                description("conflicting strand inputs specified")
            }
            UnspecifiedStrandError {
                description("strand not specified")
            }
            SubFeatureIntervalError {
                description("subfeature interval is not completely enveloped by parent")
            }
            IncompleteTranscriptError {
                description("transcript annotation is incomplete")
            }
        }
    }
}

pub trait Annotation {
    fn seq_name(&self) -> &str;
    fn interval(&self) -> &Interval<u64>;
    fn attributes(&self) -> &HashMap<String, String>;
    fn attribute(&self, key: &str) -> Option<&str>;
    fn strand(&self) -> &Strand;
    fn id(&self) -> Option<&str>;

    #[inline]
    fn span(&self) -> u64 {
        self.interval().end - self.interval().start
    }
}

macro_rules! impl_annotation {
    ($struct_ty:ty) => (

        impl Annotation for $struct_ty {

            fn seq_name(&self) -> &str {
                self.seq_name.as_str()
            }

            fn interval(&self) -> &Interval<u64> {
                &self.interval
            }

            fn attributes(&self) -> &HashMap<String, String> {
                &self.attributes
            }

            fn attribute(&self, key: &str) -> Option<&str> {
                self.attributes.get(key).map(|n| n.as_str())
            }

            fn strand(&self) -> &Strand {
                &self.strand
            }

            fn id(&self) -> Option<&str> {
                self.id.as_ref().map(|ref n| n.as_str())
            }
        }

    );
}

#[derive(Debug, Clone, PartialEq)]
pub enum TxFeature {
    Exon,
    UTR,
    UTR5,
    UTR3,
    CDS,
    StartCodon,
    StopCodon,
    Any(String),
}

fn resolve_strand_input(strand: Option<Strand>, strand_char: Option<char>) -> Result<Strand, FeatureError> {
    match (strand, strand_char) {
        (None, None) => Err(FeatureError::UnspecifiedStrandError),
        (Some(sv), None) => Ok(sv),
        (None, Some(ref scv)) => Strand::from_char(scv).map_err(FeatureError::from),
        (Some(sv), Some(ref scv)) => {
            let sv_from_char = Strand::from_char(scv).map_err(FeatureError::from)?;
            if sv == sv_from_char {
                Ok(sv)
            } else {
                Err(FeatureError::ConflictingStrandError)
            }
        }
    }
}

fn coords_to_interval(start: u64, end: u64) -> Result<Interval<u64>, FeatureError> {
    Interval::new(start..end).map_err(FeatureError::from)
}

fn resolve_transcript_features(
    transcript_seqname: &String,
    transcript_interval: &Interval<u64>,
    transcript_strand: &Strand,
    features: Option<Vec<TranscriptFeature>>,
    exon_coords: Option<&Vec<(u64, u64)>>,
    coding_coord: Option<(u64, u64)>
) -> Result<Vec<TranscriptFeature>, FeatureError>
{
    // Deliberately not handling all possible input types to avoid
    // overcomplicating code. The inputs are expected to come from
    // either GTF or refFlat after all.

    match (features, exon_coords, coding_coord) {
        // nothing defined -> the transcript doesn't have any known subfeatures
        (None, None, None) => Ok(Vec::new()),

        // only CDS defined -> must be an error
        (None, None, Some(_)) => Err(FeatureError::IncompleteTranscriptError),

        // features defined ~ takes precedence over coords (GTF input, since we need
        // to construct the tx features first to store its annotations)
        // TODO: Maybe do some checks to ensure the given features are correct?
        (Some(fxs), _, _) => Ok(fxs.into_iter().collect()),

        // exon defined & coords possibly defined (refFlat input)
        (None, Some(raw_exon_coords), raw_coding_coord) =>
            infer_features(transcript_seqname, transcript_interval,
                           transcript_strand, raw_exon_coords, raw_coding_coord),
    }
}

// Helper function for backtracking and adding possibly skipped features
fn backtrack_and_push<F>(
    features: &mut Vec<TranscriptFeature>,
    tx_feature_kind: TxFeature,
    window: &[Option<(u64, u64)>],
    codon_rem: &mut u64,
    feature_maker: &F)
where F: Fn(TxFeature, u64, u64) -> TranscriptFeature
{
    let mut backtrack_count = 1;
    let window_size = window.len();
    while *codon_rem > 0 && backtrack_count < (window_size-1) {
        // Get the nth previous item, if it exists
        if let Some(prev) = window[window_size-(backtrack_count+1)] {
            let fx = feature_maker(tx_feature_kind.clone(),
                                   max(prev.0, prev.1 - *codon_rem), prev.1);
            *codon_rem = *codon_rem - fx.span();
            backtrack_count += 1;
            features.push(fx);
        }
    };
    // Ensure the features vec is sorted if any backtrack was done
    if backtrack_count > 1 {
        features.sort_by_key(|fx| fx.interval().start)
    }
}

fn infer_features(
    transcript_seqname: &String,
    transcript_interval: &Interval<u64>,
    transcript_strand: &Strand,
    exon_coords: &Vec<(u64, u64)>,
    coding_coord: Option<(u64, u64)>
) -> Result<Vec<TranscriptFeature>, FeatureError>
{

    if exon_coords.len() == 0 {
        return Err(FeatureError::IncompleteTranscriptError);
    }

    let mut m_exon_coords = Vec::with_capacity(exon_coords.len());
    for &(a, b) in exon_coords.iter() {
        if a >= b {
            return Err(FeatureError::SubFeatureIntervalError)
        }
        m_exon_coords.push((a, b));
    }

    let exon_r = (m_exon_coords.first().unwrap().0, m_exon_coords.last().unwrap().1);

    if exon_r.0 != transcript_interval.start || exon_r.1 != transcript_interval.end {
        return Err(FeatureError::SubFeatureIntervalError);
    }

    match coding_coord {

        // Improper coding region is an error
        Some(coding_r) if coding_r.0 > coding_r.1 =>
            Err(FeatureError::SubFeatureIntervalError),

        // Presence of proper coding region means we can resolve UTRs
        Some(coding_r) if coding_r.0 < coding_r.1 => {
            // coding coord must be fully enveloped by exon max-min
            if coding_r.0 < exon_r.0 || coding_r.1 > exon_r.1 {
                return Err(FeatureError::SubFeatureIntervalError);
            }
            // TODO: no exon start == CDS end and no exon end == CDS start
            // Rough heuristic: num of features (incl exons) ~ 2 * num of exons + 4
            let mut features: Vec<TranscriptFeature> =
                Vec::with_capacity(m_exon_coords.len() * 2 + 4);

            let (utr1, utr2) = match transcript_strand {
                &Strand::Forward => (TxFeature::UTR5, TxFeature::UTR3),
                &Strand::Reverse => (TxFeature::UTR3, TxFeature::UTR5),
                &Strand::Unknown => (TxFeature::UTR, TxFeature::UTR),
            };

            let make_feature = |kind, start, end| {
                TranscriptFeature {
                    seq_name: transcript_seqname.clone(),
                    kind: kind,
                    interval: Interval::new(start..end).unwrap(),
                    strand: *transcript_strand,
                    attributes: HashMap::new(),
                    id: None,
                }
            };

            // TODO: Look at using proper trees instead of these manual steps.
            // 4 is 1 (current block) + possible # of backtracks (3, which is the codon size)
            let window_size = 4;
            // prepend `window_size-1` None so the block of interest is always at end of the window
            let iter = iter::repeat(None).take(window_size-1)
                .chain(m_exon_coords.iter().map(|&v| Some(v)));
            let mut window_storage: Storage<Option<(u64, u64)>> =
                Storage::optimized(&iter, window_size);

            // variables for tracking how much we have consumed the 5' or 3' codon
            let (mut codon1_rem, mut codon2_rem) = (3, 3);
            for window in iter.sliding_windows(&mut window_storage) {

                let (start, end) = window[window_size-1].unwrap();

                // Whole UTR exon blocks
                if end < coding_r.0 {
                    features.push(make_feature(TxFeature::Exon, start, end));
                    features.push(make_feature(utr1.clone(), start, end));

                // Whole UTR exon block with potential stop codon
                } else if end == coding_r.0 {
                    features.push(make_feature(TxFeature::Exon, start, end));
                    features.push(make_feature(utr1.clone(), start, end));
                    if let &Strand::Reverse = transcript_strand {
                        let fx = make_feature(TxFeature::StopCodon,
                                                max(start, coding_r.0 - codon1_rem),
                                                coding_r.0);
                        codon1_rem -= fx.span();
                        features.push(fx);
                        backtrack_and_push(&mut features, TxFeature::StopCodon,
                                           window.deref(), &mut codon1_rem, &make_feature);
                    }

                // UTR-CDS exon block
                } else if start < coding_r.0 {
                    features.push(make_feature(TxFeature::Exon, start, end));
                    features.push(make_feature(utr1.clone(), start, coding_r.0));

                    match transcript_strand {
                        &Strand::Forward => {
                            let fx = make_feature(TxFeature::StartCodon,
                                                  coding_r.0,
                                                  min(coding_r.0 + codon1_rem, end));
                            codon1_rem -= fx.span();
                            features.push(fx);
                        },
                        &Strand::Reverse => {
                            // most complex case: stop codon split over 3 previous 1-bp exons
                            // biologically unexpected but technically possible
                            let fx = make_feature(TxFeature::StopCodon,
                                                  max(start, coding_r.0 - codon1_rem),
                                                  coding_r.0);
                            codon1_rem -= fx.span();
                            features.push(fx);
                            backtrack_and_push(&mut features, TxFeature::StopCodon,
                                               window.deref(), &mut codon1_rem, &make_feature);
                        },
                        &Strand::Unknown => {},
                    }
                    features.push(make_feature(TxFeature::CDS, coding_r.0, end));

                // Whole CDS exon blocks
                } else if end < coding_r.1 {
                    features.push(make_feature(TxFeature::Exon, start, end));

                    if codon1_rem > 0 {
                        match transcript_strand {
                            &Strand::Forward => {
                                // Ensure the start codon coordinate is not in an intron
                                let start_ok = window.iter()
                                    .flat_map(|item| item.iter())
                                    .any(|&(a, b)| (a <= coding_r.0) && (coding_r.0 <= b));
                                if !start_ok {
                                    return Err(FeatureError::SubFeatureIntervalError);
                                }

                                let fx = make_feature(TxFeature::StartCodon,
                                                      start,
                                                      min(start + codon1_rem, end));
                                codon1_rem -= fx.span();
                                features.push(fx);
                            },
                            &Strand::Reverse if start == coding_r.0 => {
                                backtrack_and_push(&mut features, TxFeature::StopCodon,
                                                   window.deref(), &mut codon1_rem, &make_feature);
                            },
                            _ => {},
                        }
                    }
                    features.push(make_feature(TxFeature::CDS, start, end));

                // Whole CDS exon blocks at the end
                } else if end == coding_r.1 {
                    features.push(make_feature(TxFeature::Exon, start, end));
                    features.push(make_feature(TxFeature::CDS, start, end));

                    if let &Strand::Reverse = transcript_strand {
                        let fx = make_feature(TxFeature::StartCodon,
                                              max(start, coding_r.1 - codon2_rem),
                                              coding_r.1);
                        codon2_rem -= fx.span();
                        features.push(fx);
                        backtrack_and_push(&mut features, TxFeature::StartCodon,
                                           window.deref(), &mut codon2_rem, &make_feature);
                    }

                // CDS-UTR exon block
                } else if start < coding_r.1 {
                    features.push(make_feature(TxFeature::Exon, start, end));
                    features.push(make_feature(TxFeature::CDS, start, coding_r.1));

                    match transcript_strand {
                        &Strand::Forward => {
                            let fx = make_feature(TxFeature::StopCodon,
                                                  coding_r.1,
                                                  min(coding_r.1 + codon2_rem, end));
                            codon2_rem -= fx.span();
                            features.push(fx);
                        },
                        &Strand::Reverse => {
                            let fx = make_feature(TxFeature::StartCodon,
                                                  max(start, coding_r.1 - codon2_rem),
                                                  coding_r.1);
                            codon2_rem -= fx.span();
                            features.push(fx);
                            backtrack_and_push(&mut features, TxFeature::StartCodon,
                                               window.deref(), &mut codon2_rem, &make_feature);
                        },
                        &Strand::Unknown => {},
                    }
                    features.push(make_feature(utr2.clone(), coding_r.1, end));

                // Whole UTR exon blocks
                } else {
                    features.push(make_feature(TxFeature::Exon, start, end));
                    if codon2_rem > 0 {
                        match transcript_strand {
                            &Strand::Forward => {
                                let fx = make_feature(TxFeature::StopCodon,
                                                    start,
                                                    min(start + codon2_rem, end));
                                codon2_rem -= fx.span();
                                features.push(fx);
                            },
                            &Strand::Reverse => backtrack_and_push(
                                &mut features, TxFeature::StartCodon,
                                window.deref(), &mut codon2_rem, &make_feature),
                            _ => {},
                        }
                    }
                    features.push(make_feature(utr2.clone(), start, end));
                }
            }

            Ok(features)
        },

        // No CDS intervals mean we just sort the coordinates and create the exons
        _ => {
            let mut features = Vec::with_capacity(m_exon_coords.len());
            for &(start, end) in m_exon_coords.iter() {
                if start > end {
                    return Err(FeatureError::SubFeatureIntervalError)
                }
                features.push(
                    TranscriptFeature {
                        seq_name: transcript_seqname.clone(),
                        kind: TxFeature::Exon,
                        interval: Interval::new(start..end).unwrap(),
                        strand: *transcript_strand,
                        id: None,
                        attributes: HashMap::new(),
                    });
            }
            Ok(features)
        }
    }
}

#[cfg(test)]
mod test_transcript {
    use super::*;
    use self::TxFeature::*;

    fn get_coords_by_feature(fxs: &Vec<TranscriptFeature>, kind: TxFeature) -> Vec<(u64, u64)> {
        fxs.iter()
            .filter(|fx| *fx.kind() == kind)
            .map(|fx| (fx.interval().start, fx.interval().end))
            .collect()
    }

    fn get_features(fxs: &Vec<TranscriptFeature>) -> Vec<TxFeature> {
        fxs.iter().map(|fx| fx.kind().clone()).collect()
    }

    #[test]
    fn builder_coords_fwd() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Forward)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], None)
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs), vec![Exon, Exon, Exon]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
    }

    #[test]
    fn builder_coords_fwd_coding() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Forward)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 900)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs),
                   vec![Exon, UTR5, StartCodon, CDS, Exon, CDS, Exon, CDS, StopCodon, UTR3]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(100, 200)]);
        assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(200, 203)]);
        assert_eq!(get_coords_by_feature(fxs, CDS), vec![(200, 300), (400, 500), (700, 900)]);
        assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(900, 903)]);
        assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(900, 1000)]);
    }

    #[test]
    fn builder_coords_fwd_coding_from_exon5end() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Forward)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((400, 900)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs),
                   vec![Exon, UTR5, Exon, StartCodon, CDS, Exon, CDS, StopCodon, UTR3]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(100, 300)]);
        assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(400, 403)]);
        assert_eq!(get_coords_by_feature(fxs, CDS), vec![(400, 500), (700, 900)]);
        assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(900, 903)]);
        assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(900, 1000)]);
    }

    #[test]
    fn builder_coords_fwd_coding_from_exon3end() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Forward)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((300, 900)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs),
                   vec![Exon, UTR5, Exon, StartCodon, CDS, Exon, CDS, StopCodon, UTR3]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(100, 300)]);
        assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(400, 403)]);
        assert_eq!(get_coords_by_feature(fxs, CDS), vec![(400, 500), (700, 900)]);
        assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(900, 903)]);
        assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(900, 1000)]);
    }

    #[test]
    fn builder_coords_fwd_coding_from_near_exon3end() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Forward)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((297, 900)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs),
                   vec![Exon, UTR5, StartCodon, CDS, Exon, CDS, Exon, CDS, StopCodon, UTR3]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(100, 297)]);
        assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(297, 300)]);
        assert_eq!(get_coords_by_feature(fxs, CDS), vec![(297, 300), (400, 500), (700, 900)]);
        assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(900, 903)]);
        assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(900, 1000)]);
    }

     #[test]
     fn builder_coords_fwd_coding_from_split() {
         let tm = TranscriptBuilder::new("chrT", 100, 1000)
             .strand(Strand::Forward)
             .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((298, 900)))
             .build();
         assert!(tm.is_ok(), "{:?}", tm);
         let t = tm.unwrap();
         let fxs = t.features();
         assert_eq!(get_features(fxs),
                    vec![Exon, UTR5, StartCodon, CDS, Exon, StartCodon, CDS, Exon, CDS, StopCodon,
                         UTR3]);
         assert_eq!(get_coords_by_feature(fxs, Exon),
                    vec![(100, 300), (400, 500), (700, 1000)]);
         assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(100, 298)]);
         assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(298, 300), (400, 401)]);
         assert_eq!(get_coords_by_feature(fxs, CDS), vec![(298, 300), (400, 500), (700, 900)]);
         assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(900, 903)]);
         assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(900, 1000)]);
     }

     #[test]
     fn builder_coords_fwd_coding_to_exon3end() {
         let tm = TranscriptBuilder::new("chrT", 100, 1000)
             .strand(Strand::Forward)
             .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 500)))
             .build();
         assert!(tm.is_ok(), "{:?}", tm);
         let t = tm.unwrap();
         let fxs = t.features();
         assert_eq!(get_features(fxs),
                    vec![Exon, UTR5, StartCodon, CDS, Exon, CDS, Exon, StopCodon, UTR3]);
         assert_eq!(get_coords_by_feature(fxs, Exon),
                    vec![(100, 300), (400, 500), (700, 1000)]);
         assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(100, 200)]);
         assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(200, 203)]);
         assert_eq!(get_coords_by_feature(fxs, CDS), vec![(200, 300), (400, 500)]);
         assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(700, 703)]);
         assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(700, 1000)]);
     }

     #[test]
     fn builder_coords_fwd_coding_to_exon5end() {
         let tm = TranscriptBuilder::new("chrT", 100, 1000)
             .strand(Strand::Forward)
             .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 700)))
             .build();
         assert!(tm.is_ok(), "{:?}", tm);
         let t = tm.unwrap();
         let fxs = t.features();
         assert_eq!(get_features(fxs),
                    vec![Exon, UTR5, StartCodon, CDS, Exon, CDS, Exon, StopCodon, UTR3]);
         assert_eq!(get_coords_by_feature(fxs, Exon),
                    vec![(100, 300), (400, 500), (700, 1000)]);
         assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(100, 200)]);
         assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(200, 203)]);
         assert_eq!(get_coords_by_feature(fxs, CDS), vec![(200, 300), (400, 500)]);
         assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(700, 703)]);
         assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(700, 1000)]);
     }

     #[test]
     fn builder_coords_fwd_coding_to_split() {
         let tm = TranscriptBuilder::new("chrT", 100, 1000)
             .strand(Strand::Forward)
             .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 499)))
             .build();
         assert!(tm.is_ok(), "{:?}", tm);
         let t = tm.unwrap();
         let fxs = t.features();
         assert_eq!(get_features(fxs),
                    vec![Exon, UTR5, StartCodon, CDS, Exon, CDS, StopCodon, UTR3,
                         Exon, StopCodon, UTR3]);
         assert_eq!(get_coords_by_feature(fxs, Exon),
                    vec![(100, 300), (400, 500), (700, 1000)]);
         assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(100, 200)]);
         assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(200, 203)]);
         assert_eq!(get_coords_by_feature(fxs, CDS), vec![(200, 300), (400, 499)]);
         assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(499, 500), (700, 702)]);
         assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(499, 500), (700, 1000)]);
     }

    #[test]
    fn builder_coords_rev() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Reverse)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], None)
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs), vec![Exon, Exon, Exon]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
    }

    #[test]
    fn builder_coords_rev_coding() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Reverse)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 900)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs),
                   vec![Exon, UTR3, StopCodon, CDS, Exon, CDS, Exon, CDS, StartCodon, UTR5]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(100, 200)]);
        assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(197, 200)]);
        assert_eq!(get_coords_by_feature(fxs, CDS), vec![(200, 300), (400, 500), (700, 900)]);
        assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(897, 900)]);
        assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(900, 1000)]);
    }

    #[test]
    fn builder_coords_rev_coding_from_exon3end() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Reverse)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((300, 900)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs),
                   vec![Exon, UTR3, StopCodon, Exon, CDS, Exon, CDS, StartCodon, UTR5]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(900, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(897, 900)]);
        assert_eq!(get_coords_by_feature(fxs, CDS), vec![(400, 500), (700, 900)]);
        assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(297, 300)]);
        assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(100, 300)]);
    }

    #[test]
    fn builder_coords_rev_coding_from_exon5end() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Reverse)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((400, 900)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs),
                   vec![Exon, UTR3, StopCodon, Exon, CDS, Exon, CDS, StartCodon, UTR5]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(900, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(897, 900)]);
        assert_eq!(get_coords_by_feature(fxs, CDS), vec![(400, 500), (700, 900)]);
        assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(297, 300)]);
        assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(100, 300)]);
    }

    #[test]
    fn builder_coords_rev_coding_to_exon5end() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Reverse)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 700)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs),
                   vec![Exon, UTR3, StopCodon, CDS, Exon, CDS, StartCodon, Exon, UTR5]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(497, 500)]);
        assert_eq!(get_coords_by_feature(fxs, CDS), vec![(200, 300), (400, 500)]);
        assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(197, 200)]);
        assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(100, 200)]);
    }

    #[test]
    fn builder_coords_rev_coding_to_near_exon5end() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Reverse)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 703)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs),
                   vec![Exon, UTR3, StopCodon, CDS, Exon, CDS, Exon, CDS, StartCodon, UTR5]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(703, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(700, 703)]);
        assert_eq!(get_coords_by_feature(fxs, CDS), vec![(200, 300), (400, 500), (700, 703)]);
        assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(197, 200)]);
        assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(100, 200)]);
    }

    #[test]
    fn builder_coords_rev_coding_to_exon3end() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Reverse)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 500)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs),
                   vec![Exon, UTR3, StopCodon, CDS, Exon, CDS, StartCodon, Exon, UTR5]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, UTR5), vec![(700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, StartCodon), vec![(497, 500)]);
        assert_eq!(get_coords_by_feature(fxs, CDS), vec![(200, 300), (400, 500)]);
        assert_eq!(get_coords_by_feature(fxs, StopCodon), vec![(197, 200)]);
        assert_eq!(get_coords_by_feature(fxs, UTR3), vec![(100, 200)]);
    }

    #[test]
    fn builder_coords_coding_unk() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Unknown)
            .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 900)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        let fxs = t.features();
        assert_eq!(get_features(fxs),
                   vec![Exon, UTR, CDS, Exon, CDS, Exon, CDS, UTR]);
        assert_eq!(get_coords_by_feature(fxs, Exon),
                   vec![(100, 300), (400, 500), (700, 1000)]);
        assert_eq!(get_coords_by_feature(fxs, CDS), vec![(200, 300), (400, 500), (700, 900)]);
        assert_eq!(get_coords_by_feature(fxs, UTR), vec![(100, 200), (900, 1000)]);
    }
}

#[derive(Debug)]
pub struct TranscriptFeature {
    kind: TxFeature,
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    id: Option<String>,
    attributes: HashMap<String, String>
}

impl_annotation!(TranscriptFeature);

impl TranscriptFeature {

    pub fn kind(&self) -> &TxFeature {
        &self.kind
    }
}

pub struct TxFeatureBuilder {
    seq_name: String,
    start: u64,
    end: u64,
    attributes: HashMap<String, String>,
    kind: TxFeature,
    strand: Option<Strand>,
    strand_char: Option<char>,
    id: Option<String>,
}

impl TxFeatureBuilder {

    pub fn new<T>(seq_name: T, start: u64, end: u64, kind: TxFeature) -> TxFeatureBuilder
        where T: Into<String>
    {
        TxFeatureBuilder {
            seq_name: seq_name.into(),
            start: start, end: end,
            kind: kind,
            attributes: HashMap::new(),
            strand: None,
            strand_char: None,
            id: None,
        }
    }

    pub fn strand(mut self, strand: Strand) -> TxFeatureBuilder {
        self.strand = Some(strand);
        self
    }

    pub fn strand_char(mut self, strand_char: char) -> TxFeatureBuilder {
        self.strand_char = Some(strand_char);
        self
    }

    pub fn attribute<K, V>(mut self, key: K, value: V) -> TxFeatureBuilder
        where K: Into<String>, V: Into<String>
    {
        self.attributes.insert(key.into(), value.into());
        self
    }

    pub fn id<T>(mut self, id: T) -> TxFeatureBuilder
        where T: Into<String>
    {
        self.id = Some(id.into());
        self
    }

    pub fn build(self) -> Result<TranscriptFeature, FeatureError> {
        let interval = coords_to_interval(self.start, self.end)?;
        let strand = resolve_strand_input(self.strand, self.strand_char)?;
        let feature = TranscriptFeature {
            seq_name: self.seq_name, kind: self.kind, interval: interval,
            strand: strand, attributes: self.attributes, id: self.id,
        };
        Ok(feature)
    }
}

#[derive(Debug)]
pub struct Transcript {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    id: Option<String>,
    attributes: HashMap<String, String>,
    features: Vec<TranscriptFeature>,
}

impl_annotation!(Transcript);

impl Transcript {

    pub fn features(&self) -> &Vec<TranscriptFeature> {
        &self.features
    }
}

pub struct TranscriptBuilder {
    seq_name: String,
    start: u64, end: u64,
    strand: Option<Strand>,
    strand_char: Option<char>,
    // Input can be a vector of pre-made features ...
    features: Option<Vec<TranscriptFeature>>,
    // Or exon coordinates, possibly coupled with cds coord
    // NOTE: Can we instead of using Vec<_> here keep it as an unconsumed iterator?
    exon_coords: Option<Vec<(u64, u64)>>,
    coding_coord: Option<(u64, u64)>,
    id: Option<String>,
    attributes: HashMap<String, String>,
}

impl TranscriptBuilder {

    pub fn new<T>(seq_name: T, start: u64, end: u64) -> TranscriptBuilder
        where T: Into<String>
    {
        TranscriptBuilder {
            seq_name: seq_name.into(),
            start: start, end: end,
            strand: None,
            strand_char: None,
            features: None,
            exon_coords: None,
            coding_coord: None,
            id: None,
            attributes: HashMap::new(),
        }
    }

    pub fn strand(mut self, strand: Strand) -> TranscriptBuilder {
        self.strand = Some(strand);
        self
    }

    pub fn strand_char(mut self, strand_char: char) -> TranscriptBuilder {
        self.strand_char = Some(strand_char);
        self
    }

    pub fn attribute<K, V>(mut self, key: K, value: V) -> TranscriptBuilder
        where K: Into<String>, V: Into<String>
    {
        self.attributes.insert(key.into(), value.into());
        self
    }

    pub fn id<T>(mut self, id: T) -> TranscriptBuilder
        where T: Into<String>
    {
        self.id = Some(id.into());
        self
    }

    pub fn feature_coords<E>(
        mut self,
        exon_coords: E,
        coding_coord: Option<(u64, u64)>
    )-> TranscriptBuilder
        where E: IntoIterator<Item=(u64, u64)>
    {
        self.exon_coords = Some(exon_coords.into_iter().collect());
        self.coding_coord = coding_coord;
        self
    }

    pub fn build(self) -> Result<Transcript, FeatureError> {
        let interval = coords_to_interval(self.start, self.end)?;
        let strand = resolve_strand_input(self.strand, self.strand_char)?;
        let features = resolve_transcript_features(
            &self.seq_name, &interval, &strand,
            self.features, self.exon_coords.as_ref(), self.coding_coord)?;

        let transcript = Transcript {
            seq_name: self.seq_name, interval: interval, strand: strand,
            features: features, attributes: self.attributes, id: self.id,
        };
        Ok(transcript)
    }
}

#[derive(Debug)]
pub struct Gene {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    id: Option<String>,
    attributes: HashMap<String, String>,
    transcripts: HashMap<String, Transcript>,
}

impl_annotation!(Gene);


#[cfg(test)]
mod test_transcript_feature {
    use super::*;
    use self::TxFeature::*;

    #[test]
    fn builder() {
        let tfm1 = TxFeatureBuilder::new("chrT", 10, 20, Exon)
            .strand(Strand::Forward)
            .attribute("name", "ex1")
            .id("ex1.1")
            .build();
        assert!(tfm1.is_ok());
        let tf = tfm1.unwrap();
        assert_eq!(tf.seq_name(), "chrT");
        assert_eq!(tf.kind(), &TxFeature::Exon);
        assert_eq!(tf.strand(), &Strand::Forward);
        assert_eq!(tf.id(), Some("ex1.1"));
        assert_eq!(tf.attribute("name"), Some("ex1"));
        assert_eq!(tf.attributes.len(), 1);

        let tfm2 = TxFeatureBuilder::new("chrO", 10, 10, Exon)
            .strand_char('-')
            .strand(Strand::Reverse)
            .build();
        assert!(tfm2.is_ok());
    }

    #[test]
    fn builder_interval_invalid() {
        let tfm = TxFeatureBuilder::new("chrE", 20, 10, Exon)
            .build();
        assert!(tfm.is_err());
        assert!(matches!(tfm.unwrap_err(),
                         FeatureError::IntervalError(utils::IntervalError::InvalidRange)));
    }

    #[test]
    fn builder_strand_unspecified() {
        let tfm = TxFeatureBuilder::new("chrT", 20, 30, Exon)
            .build();
        assert!(tfm.is_err());
        assert!(matches!(tfm.unwrap_err(), FeatureError::UnspecifiedStrandError));
    }

    #[test]
    fn builder_strand_char_unexpected() {
        let tfm = TxFeatureBuilder::new("chrE", 10, 20, Exon)
            .strand_char('w')
            .build();
        assert!(tfm.is_err());
        assert!(matches!(tfm.unwrap_err(),
                         FeatureError::StrandCharError(utils::StrandError::InvalidChar(_))));
    }

    #[test]
    fn builder_strand_char_conflicting() {
        let tfm = TxFeatureBuilder::new("chrE", 10, 20, Exon)
            .strand_char('-')
            .strand(Strand::Reverse)
            .build();
        assert!(tfm.is_ok());
        let tf = tfm.unwrap();
        assert_eq!(tf.strand(), &Strand::Reverse);
    }
}
