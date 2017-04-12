extern crate bio;
#[macro_use]
#[cfg(test)]
extern crate matches;
#[macro_use]
extern crate quick_error;
extern crate sliding_windows;

use std::cmp::{max, min};
use std::collections::HashMap;
use std::iter;
use std::ops::Deref;

use bio::utils::{self, Interval, IntervalError, Strand};
use sliding_windows::{IterExt, Storage};

use self::error::FeatureError;
use self::TranscriptFeatureKind::*;

#[cfg(test)]
mod test;


macro_rules! impl_feature {
    ($struct_ty:ty) => (

        impl $struct_ty {

            fn seq_name(&self) -> &str {
                self.seq_name.as_str()
            }

            fn interval(&self) -> &Interval<u64> {
                &self.interval
            }

            fn strand(&self) -> &Strand {
                &self.strand
            }

            #[inline]
            fn span(&self) -> u64 {
                self.interval().end - self.interval().start
            }
        }

    );
}

#[derive(Debug, Clone, PartialEq)]
pub enum TranscriptFeatureKind {
    Exon,
    UTR,
    UTR5,
    UTR3,
    CDS,
    StartCodon,
    StopCodon,
    Intron,
    Any(String),
}

#[derive(Debug)]
pub struct TranscriptFeature {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    pub id: Option<String>,
    pub attributes: HashMap<String, String>,
    kind: TranscriptFeatureKind,
    frame: Option<u64>,
}

impl_feature!(TranscriptFeature);

impl TranscriptFeature {

    pub fn kind(&self) -> &TranscriptFeatureKind {
        &self.kind
    }

    pub fn frame(&self) -> Option<u64> {
        self.frame
    }
}

pub struct TFBuilder {
    seq_name: String,
    start: u64,
    end: u64,
    strand: Option<Strand>,
    strand_char: Option<char>,
    pub id: Option<String>,
    pub attributes: HashMap<String, String>,
    kind: TranscriptFeatureKind,
    frame: Option<u64>,
}

impl TFBuilder {

    pub fn new<T>(seq_name: T, start: u64, end: u64, kind: TranscriptFeatureKind) -> TFBuilder
        where T: Into<String>
    {
        TFBuilder {
            seq_name: seq_name.into(),
            start: start, end: end,
            kind: kind,
            attributes: HashMap::new(),
            strand: None,
            strand_char: None,
            id: None,
            frame: None,
        }
    }

    pub fn strand(mut self, strand: Strand) -> TFBuilder {
        self.strand = Some(strand);
        self
    }

    pub fn strand_char(mut self, strand_char: char) -> TFBuilder {
        self.strand_char = Some(strand_char);
        self
    }

    pub fn attribute<K, V>(mut self, key: K, value: V) -> TFBuilder
        where K: Into<String>, V: Into<String>
    {
        self.attributes.insert(key.into(), value.into());
        self
    }

    pub fn id<T>(mut self, id: T) -> TFBuilder
        where T: Into<String>
    {
        self.id = Some(id.into());
        self
    }

    pub fn frame(mut self, frame: u64) -> TFBuilder {
        self.frame = Some(frame);
        self
    }

    pub fn build(self) -> Result<TranscriptFeature, FeatureError> {
        let interval = coords_to_interval(self.start, self.end)?;
        let strand = resolve_strand_input(self.strand, self.strand_char)?;
        let feature = TranscriptFeature {
            seq_name: self.seq_name, kind: self.kind, interval: interval,
            strand: strand, attributes: self.attributes, id: self.id,
            frame: self.frame,
        };
        Ok(feature)
    }
}

#[derive(Debug)]
pub struct Transcript {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    pub id: Option<String>,
    pub attributes: HashMap<String, String>,
    features: Vec<TranscriptFeature>,
}

impl_feature!(Transcript);

impl Transcript {

    pub fn features(&self) -> &Vec<TranscriptFeature> {
        &self.features
    }
}

pub struct TBuilder {
    seq_name: String,
    start: u64,
    end: u64,
    strand: Option<Strand>,
    strand_char: Option<char>,
    id: Option<String>,
    attributes: HashMap<String, String>,
    // Input can be a vector of pre-made features ...
    features: Option<Vec<TranscriptFeature>>,
    // Or exon coordinates, possibly coupled with cds coord
    // NOTE: Can we instead of using Vec<_> here keep it as an unconsumed iterator?
    exon_coords: Option<Vec<(u64, u64)>>,
    coding_coord: Option<(u64, u64)>,
}

impl TBuilder {

    pub fn new<T>(seq_name: T, start: u64, end: u64) -> TBuilder
        where T: Into<String>
    {
        TBuilder {
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

    pub fn strand(mut self, strand: Strand) -> TBuilder {
        self.strand = Some(strand);
        self
    }

    pub fn strand_char(mut self, strand_char: char) -> TBuilder {
        self.strand_char = Some(strand_char);
        self
    }

    pub fn attribute<K, V>(mut self, key: K, value: V) -> TBuilder
        where K: Into<String>, V: Into<String>
    {
        self.attributes.insert(key.into(), value.into());
        self
    }

    pub fn id<T>(mut self, id: T) -> TBuilder
        where T: Into<String>
    {
        self.id = Some(id.into());
        self
    }

    pub fn feature_coords<E>(
        mut self,
        exon_coords: E,
        coding_coord: Option<(u64, u64)>
    )-> TBuilder
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

impl_feature!(Gene);

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

fn coords_to_interval(start: u64, end: u64) -> Result<Interval<u64>, FeatureError> {
    Interval::new(start..end).map_err(FeatureError::from)
}

fn resolve_strand_input(
    strand: Option<Strand>,
    strand_char: Option<char>)
-> Result<Strand, FeatureError>
{
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

        Some(coding_r) => {
            // Improper coding region is an error
            if coding_r.0 > coding_r.1 {
                return Err(FeatureError::SubFeatureIntervalError);
            }
            // Coding coord must be fully enveloped by exon max-min
            if coding_r.0 < exon_r.0 || coding_r.1 > exon_r.1 {
                return Err(FeatureError::SubFeatureIntervalError);
            }
            // There must be room for stop codons (which is not inclusive in coding_coord)
            let stop_codon_ok = match transcript_strand {
                &Strand::Forward => coding_r.1 + 3 <= exon_r.1,
                &Strand::Reverse => coding_r.0 - 3 >= exon_r.0,
                &Strand::Unknown => true,
            };
            if !stop_codon_ok {
                return Err(FeatureError::SubFeatureIntervalError);
            }
            infer_coding_features(&m_exon_coords, coding_r, &transcript_seqname, transcript_strand)
        }

        // No CDS intervals mean we just sort the coordinates and create the exons
        None => {
            let mut features = Vec::with_capacity(m_exon_coords.len());
            for &(start, end) in m_exon_coords.iter() {
                features.push(
                    TranscriptFeature {
                        seq_name: transcript_seqname.clone(),
                        kind: Exon,
                        interval: Interval::new(start..end).unwrap(),
                        strand: *transcript_strand,
                        id: None,
                        attributes: HashMap::new(),
                        frame: None,
                    });
            }
            Ok(features)
        }
    }
}

// Window size for iteration to make CDS features.
// 4 is 1 (current block) + possible # of backtracks (3, which is the codon size)
const WINDOW_SIZE: usize = 4;

// requirements:
//  - exon coords sorted and nonempty
//  - coding coord within exon span
//  - coding_coord.0 < coding_coord.1
fn infer_coding_features(
    exon_coords: &Vec<(u64, u64)>,
    coding_r: (u64, u64),
    transcript_seqname: &String,
    transcript_strand: &Strand)
-> Result<Vec<TranscriptFeature>, FeatureError> {

    let mut features: Vec<TranscriptFeature> =
        Vec::with_capacity(exon_coords.len() * 2 + 4);

    let (utr1, utr2) = match transcript_strand {
        &Strand::Forward => (UTR5, UTR3),
        &Strand::Reverse => (UTR3, UTR5),
        &Strand::Unknown => (UTR, UTR),
    };

    let feat = |kind, start, end| {
        TranscriptFeature {
            seq_name: transcript_seqname.clone(),
            kind: kind,
            interval: Interval::new(start..end).unwrap(),
            strand: *transcript_strand,
            attributes: HashMap::new(),
            id: None,
            frame: None,
        }
    };

    // prepend `window_size-1` of `None`s so the block of interest is always at end of the window
    let iter = iter::repeat(None).take(WINDOW_SIZE-1)
        .chain(exon_coords.iter().map(|&v| Some(v)));
    let mut window_storage: Storage<Option<(u64, u64)>> = Storage::optimized(&iter, WINDOW_SIZE);

    // how much we have consumed the 5' or 3' codon
    let (mut codon1_rem, mut codon2_rem) = (3, 3);

    for window in iter.sliding_windows(&mut window_storage) {

        let (start, end) = window[WINDOW_SIZE-1].unwrap();

        if start < coding_r.0 {

            if end < coding_r.0 {
                features.push(feat(Exon, start, end));
                features.push(feat(utr1.clone(), start, end));

            } else if end == coding_r.0 {
                features.push(feat(Exon, start, end));
                features.push(feat(utr1.clone(), start, end));
                if let &Strand::Reverse = transcript_strand {
                    let fx = feat(StopCodon, max(start, coding_r.0 - codon1_rem), coding_r.0);
                    codon1_rem -= fx.span();
                    features.push(fx);
                    codon1_rem = backtrack_and_push(
                        &mut features, StopCodon, window.deref(), codon1_rem, &feat);
                }

            } else if end > coding_r.0 && end < coding_r.1 {
                features.push(feat(Exon, start, end));
                features.push(feat(utr1.clone(), start, coding_r.0));
                match transcript_strand {
                    &Strand::Forward => {
                        let fx = feat(StartCodon, coding_r.0, min(end, coding_r.0 + 3));
                        codon1_rem -= fx.span();
                        features.push(fx);
                    },
                    &Strand::Reverse => {
                        let fx = feat(StopCodon, max(start, coding_r.0 - codon1_rem), coding_r.0);
                        codon1_rem -= fx.span();
                        features.push(fx);
                        codon1_rem = backtrack_and_push(
                            &mut features, StopCodon, window.deref(), codon1_rem, &feat);
                    },
                    &Strand::Unknown => {},
                }
                features.push(feat(CDS, coding_r.0, end));

            } else if end == coding_r.1 {
                if coding_r.1 - coding_r.0 < 3 {
                    // a coding region must have at least 3 bases for the start codon
                    return Err(FeatureError::SubFeatureIntervalError);
                }
                features.push(feat(Exon, start, end));
                features.push(feat(utr1.clone(), start, coding_r.0));
                match transcript_strand {
                    &Strand::Forward => {
                        let fx = feat(StartCodon, coding_r.0, coding_r.0 + 3);
                        codon1_rem -= fx.span();
                        features.push(fx);
                    },
                    &Strand::Reverse => {
                        let fx = feat(StopCodon, start, coding_r.0 - codon1_rem);
                        codon1_rem -= fx.span();
                        features.push(fx);
                        codon1_rem = backtrack_and_push(
                            &mut features, StopCodon, window.deref(), codon1_rem, &feat);
                    },
                    &Strand::Unknown => {},
                }
                features.push(feat(CDS, coding_r.0, coding_r.1));

            } else if end > coding_r.1 {
                if coding_r.1 - coding_r.0 < 3 {
                    // coding region must have at least 3 bases
                    return Err(FeatureError::SubFeatureIntervalError);
                }
                features.push(feat(Exon, start, end));
                features.push(feat(utr1.clone(), start, coding_r.0));
                match transcript_strand {
                    &Strand::Forward => {
                        let fx = feat(StartCodon, coding_r.0, coding_r.0 + 3);
                        codon1_rem -= 3;
                        features.push(fx);
                        features.push(feat(CDS, coding_r.0, coding_r.1));
                        let fx = feat(StopCodon, coding_r.1, min(end, coding_r.1 + codon2_rem));
                        codon2_rem -= fx.span();
                        features.push(fx);
                    },
                    &Strand::Reverse => {
                        let fx = feat(StopCodon, max(start, coding_r.0 - codon2_rem), coding_r.0);
                        codon1_rem -= fx.span();
                        features.push(fx);
                        codon1_rem = backtrack_and_push(
                            &mut features, StopCodon, window.deref(), codon1_rem, &feat);
                        features.push(feat(CDS, coding_r.0, coding_r.1));
                        let fx = feat(StartCodon, coding_r.1 - codon2_rem, coding_r.1);
                        codon2_rem -= fx.span();
                        features.push(fx);
                    },
                    &Strand::Unknown => {
                        features.push(feat(CDS, coding_r.0, coding_r.1));
                    },
                }
                features.push(feat(utr2.clone(), coding_r.1, end));

            } else {
                assert!(false, "unexpected: exon=[{},{}) cds=[{},{})",
                        start, end, coding_r.0, coding_r.1)
            }

        } else if start == coding_r.0 {

            if end < coding_r.1 {
                features.push(feat(Exon, start, end));
                match transcript_strand {
                    &Strand::Forward => {
                        let fx = feat(StartCodon, start, min(start + 3, end));
                        codon1_rem -= fx.span();
                        features.push(fx);
                    },
                    &Strand::Reverse => {
                        codon1_rem = backtrack_and_push(
                            &mut features, StopCodon, window.deref(), codon1_rem, &feat);
                    },
                    &Strand::Unknown => {},
                }
                features.push(feat(CDS, start, end));

            } else if end == coding_r.1 {
                if coding_r.1 - coding_r.0 < 3 {
                    // coding region must have at least 3 bases
                    return Err(FeatureError::SubFeatureIntervalError);
                }
                features.push(feat(Exon, start, end));
                match transcript_strand {
                    &Strand::Forward => {
                        features.push(feat(StartCodon, start, start + 3));
                        codon1_rem = 0;
                        features.push(feat(CDS, start, end));
                    },
                    &Strand::Reverse => {
                        codon1_rem = backtrack_and_push(
                            &mut features, StopCodon, window.deref(), codon1_rem, &feat);
                        features.push(feat(CDS, start, end));
                        features.push(feat(StartCodon, end - 3, end));
                        codon2_rem = 0;
                    },
                    &Strand::Unknown => {
                        features.push(feat(CDS, start, end));
                    },
                }

            } else if end > coding_r.1 {
                if coding_r.1 - coding_r.0 < 3 {
                    // coding region must have at least 3 bases
                    return Err(FeatureError::SubFeatureIntervalError);
                }
                features.push(feat(Exon, start, end));
                match transcript_strand {
                    &Strand::Forward => {
                        features.push(feat(StartCodon, start, start + 3));
                        codon1_rem -= 3;
                        features.push(feat(CDS, start, coding_r.1));
                        let fx = feat(StopCodon, coding_r.1, min(end, coding_r.1 + codon2_rem));
                        codon2_rem -= fx.span();
                        features.push(fx);
                    },
                    &Strand::Reverse => {
                        codon1_rem = backtrack_and_push(
                            &mut features, StopCodon, window.deref(), codon1_rem, &feat);
                        features.push(feat(CDS, start, end));
                        features.push(feat(StartCodon, end - 3, end));
                        codon2_rem -= 3;
                    },
                    &Strand::Unknown => {
                        features.push(feat(CDS, start, coding_r.1));
                    },
                }
                features.push(feat(utr2.clone(), coding_r.1, end));

            } else {
                assert!(false, "unexpected: exon=[{},{}) cds=[{},{})",
                        start, end, coding_r.0, coding_r.1)
            }

        } else if start > coding_r.0 && start < coding_r.1 {

            if end < coding_r.1 {
                features.push(feat(Exon, start, end));
                if let &Strand::Forward = transcript_strand {
                    if codon1_rem > 0 {
                        let fx = feat(StartCodon, start, min(end, start + codon1_rem));
                        codon1_rem -= fx.span();
                        features.push(fx);
                    }
                }
                features.push(feat(CDS, start, end));

            } else if end == coding_r.1 {
                features.push(feat(Exon, start, end));
                match transcript_strand {
                    &Strand::Forward => {
                        if codon1_rem > 0 {
                            let fx = feat(StartCodon, start, min(end, start + codon1_rem));
                            codon1_rem -= fx.span();
                            features.push(fx);
                        }
                        features.push(feat(CDS, start, end));
                    },
                    &Strand::Reverse => {
                        features.push(feat(CDS, start, end));
                        let fx = feat(StartCodon, max(start, coding_r.1 - codon2_rem), coding_r.1);
                        codon2_rem -= fx.span();
                        features.push(fx);
                        codon2_rem = backtrack_and_push(
                            &mut features, StartCodon, window.deref(), codon2_rem, &feat);
                    },
                    &Strand::Unknown => {
                        features.push(feat(CDS, start, end));
                    },
                }


            } else if end > coding_r.1 {
                features.push(feat(Exon, start, end));
                match transcript_strand {
                    &Strand::Forward => {
                        if codon1_rem > 0 {
                            let fx = feat(StartCodon, start, min(coding_r.1, start + codon1_rem));
                            codon1_rem -= fx.span();
                            if codon1_rem > 0 {
                                // codon1 bases must be exhausted at this point
                                return Err(FeatureError::SubFeatureIntervalError);
                            }
                            features.push(fx);
                        }
                        features.push(feat(CDS, start, coding_r.1));
                        let fx = feat(StopCodon, coding_r.1, min(coding_r.1 + codon2_rem, end));
                        codon2_rem -= fx.span();
                        features.push(fx);
                        features.push(feat(utr2.clone(), coding_r.1, end));
                    },
                    &Strand::Reverse => {
                        features.push(feat(CDS, start, coding_r.1));
                        let fx = feat(StartCodon, max(start, coding_r.1 - codon2_rem), coding_r.1);
                        codon2_rem -= fx.span();
                        features.push(fx);
                        codon2_rem = backtrack_and_push(
                            &mut features, StartCodon, window.deref(), codon2_rem, &feat);
                        features.push(feat(utr2.clone(), coding_r.1, end));
                    },
                    &Strand::Unknown => {
                        features.push(feat(CDS, start, coding_r.1));
                        features.push(feat(utr2.clone(), coding_r.1, end));
                    },
                }

            } else {
                assert!(false, "unexpected: exon=[{},{}) cds=[{},{})",
                        start, end, coding_r.0, coding_r.1)
            }

        } else if start >= coding_r.1 {
            features.push(feat(Exon, start, end));
            match transcript_strand {
                &Strand::Forward => {
                    if codon2_rem > 0 {
                        let fx = feat(StopCodon, start, min(start + codon2_rem, end));
                        codon2_rem -= fx.span();
                        features.push(fx);
                    }
                },
                &Strand::Reverse => {
                    codon2_rem = backtrack_and_push(
                        &mut features, StartCodon, window.deref(), codon2_rem, &feat);
                },
                &Strand::Unknown => {},
            }
            features.push(feat(utr2.clone(), start, end));

        } else {
            assert!(false, "unexpected: exon=[{},{}) cds=[{},{})",
                    start, end, coding_r.0, coding_r.1)
        }
    }

    match transcript_strand {
        &Strand::Forward => set_coding_frames(features.iter_mut()),
        &Strand::Reverse => set_coding_frames(features.iter_mut().rev()),
        _ => {}
    }

    Ok(features)
}

// Helper function for backtracking and adding possibly skipped features
fn backtrack_and_push<F>(
    features: &mut Vec<TranscriptFeature>,
    tfk: TranscriptFeatureKind,
    window: &[Option<(u64, u64)>],
    mut codon_rem: u64,
    feature_maker: &F) -> u64
where F: Fn(TranscriptFeatureKind, u64, u64) -> TranscriptFeature
{
    let mut backtrack_count = 1;
    while codon_rem > 0 && backtrack_count < (WINDOW_SIZE-1) {
        // Get the nth previous item, if it exists
        if let Some(prev) = window[WINDOW_SIZE-(backtrack_count+1)] {
            let fx = feature_maker(tfk.clone(), max(prev.0, prev.1 - codon_rem), prev.1);
            codon_rem -= fx.span();
            features.push(fx);
        }
        backtrack_count += 1;
    };
    // Ensure the features vec is sorted if any backtrack was done
    if backtrack_count > 1 {
        features.sort_by_key(|fx| fx.interval().start)
    };
    codon_rem
}

fn set_coding_frames<'a, T>(features_miter: T)
where T: Iterator<Item=&'a mut TranscriptFeature>
{
    let (mut sc_frame, mut cds_frame) = (0, 0);
    for coding_fx in features_miter {
        match coding_fx.kind() {
            &StartCodon => {
                coding_fx.frame = Some(sc_frame);
                sc_frame = calc_next_frame(coding_fx.span(), sc_frame);
            },
            &CDS => {
                coding_fx.frame = Some(cds_frame);
                cds_frame = calc_next_frame(coding_fx.span(), cds_frame);
            },
            _ => {},
        }
    }
}

#[inline(always)]
// Adapted from: http://mblab.wustl.edu/GTF22.html
fn calc_next_frame(cur_span: u64, cur_frame: u64) -> u64 {
    if cur_span >= cur_frame {
        (3 - ((cur_span - cur_frame) % 3)) % 3
    } else {
        (3 - (cur_frame - cur_span) % 3)
    }
}
