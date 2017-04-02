//! Interval-based annotation features.

extern crate bio;

use std::cmp::{max, min};
use std::error::Error;
use std::collections::HashMap;
use std::fmt::{self, Display};
use std::iter::FromIterator;

use bio::data_structures::interval_tree::IntervalTree;
use bio::utils::{Interval, IntervalError, Strand};

use self::error::FeatureError;


pub mod error {

    use super::*;

    #[derive(Debug, PartialEq)]
    pub enum FeatureError {
        IntervalError,
        StrandCharError,
        ConflictingStrandError,
        UnspecifiedStrandError,
        SubFeatureIntervalError,
        IncompleteTranscriptError,
    }

    impl Error for FeatureError {

        fn description(&self) -> &str {
            match *self {
                FeatureError::IntervalError =>
                    "interval start coordinate larger than its end coordinate",
                FeatureError::StrandCharError =>
                    "strand character is invalid",
                FeatureError::ConflictingStrandError =>
                    "conflicting strand inputs specified",
                FeatureError::UnspecifiedStrandError =>
                    "strand not specified",
                FeatureError::SubFeatureIntervalError =>
                    "subfeature interval is not completely enveloped by parent",
                FeatureError::IncompleteTranscriptError =>
                    "transcript annotation is incomplete",
            }
        }

        fn cause(&self) -> Option<&Error> {
            None
        }
    }

    impl Display for FeatureError {

        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "FeatureError: {}", self.description())
        }
    }

    impl From<IntervalError> for FeatureError {
        fn from(err: IntervalError) -> FeatureError {
            FeatureError::IntervalError
        }
    }
}

trait Annotation {
    fn seq_name(&self) -> &str;
    fn interval(&self) -> &Interval<u64>;
    fn attributes(&self) -> &HashMap<String, String>;
    fn attribute(&self, key: &str) -> Option<&str>;
    fn strand(&self) -> &Strand;
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
        }

    );
}

fn char_to_strand(strand_char: char) -> Result<Strand, FeatureError> {
    match strand_char {
        '+' | 'f' | 'F' => Ok(Strand::Forward),
        '-' | 'r' | 'R' => Ok(Strand::Reverse),
        '.' | '?' => Ok(Strand::Unknown),
        _ => Err(FeatureError::StrandCharError),
    }
}

fn resolve_strand_input(strand: Option<Strand>, strand_char: Option<char>) -> Result<Strand, FeatureError> {
    match (strand, strand_char) {
        (None, None) => Err(FeatureError::UnspecifiedStrandError),
        (Some(sv), None) => Ok(sv),
        (None, Some(scv)) => char_to_strand(scv),
        (Some(sv), Some(scv)) => {
            let sv_from_char = char_to_strand(scv)?;
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
    cds_coord: Option<(u64, u64)>
) -> Result<Vec<TranscriptFeature>, FeatureError>
{
    // Deliberately not handling all possible input types to avoid
    // overcomplicating code. The inputs are expected to come from
    // either GTF or refFlat after all.

    match (features, exon_coords, cds_coord) {
        // nothing defined -> the transcript doesn't have any known subfeatures
        (None, None, None) => Ok(Vec::new()),

        // only CDS defined -> must be an error
        (None, None, Some(_)) => Err(FeatureError::IncompleteTranscriptError),

        // features defined ~ takes precedence over coords (GTF input, since we need
        // to construct the tx features first to store its annotations)
        // TODO: Maybe do some checks to ensure the given features are correct?
        (Some(fxs), _, _) => Ok(fxs.into_iter().collect()),

        // exon defined & coords possibly defined (refFlat input)
        (None, Some(raw_exon_coords), cdsc) =>
            infer_features(transcript_seqname, transcript_interval,
                           transcript_strand, raw_exon_coords, cdsc),
    }
}

#[cfg(test)]
mod test_transcript {
    use super::*;

    #[test]
    fn builder_from_coords() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Forward)
            .exon_and_cds_coords(vec![(100, 300), (400, 500), (700, 1000)], None)
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        assert_eq!(t.features().len(), 3, "{:?}", t.features());
    }

    #[test]
    fn builder_from_coords_with_cds() {
        let tm = TranscriptBuilder::new("chrT", 100, 1000)
            .strand(Strand::Forward)
            .exon_and_cds_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 900)))
            .build();
        assert!(tm.is_ok(), "{:?}", tm);
        let t = tm.unwrap();
        assert_eq!(t.features().len(), 10, "{:?}", t.features());
    }
}

fn infer_features(
    transcript_seqname: &String,
    transcript_interval: &Interval<u64>,
    transcript_strand: &Strand,
    exon_coords: &Vec<(u64, u64)>,
    cds_coord: Option<(u64, u64)>
) -> Result<Vec<TranscriptFeature>, FeatureError>
{

    // A transcript must have an exon at this point
    let first_exon_coord = exon_coords.first().ok_or(FeatureError::IncompleteTranscriptError)?;
    // Get min exon start and max exon end for some sanity checks
    let exon_iv: (u64, u64) = exon_coords.iter().skip(1)
        .fold(*first_exon_coord, |acc, &(a, b)| (min(acc.0, a), max(acc.1, b)));

    if exon_iv.0 != transcript_interval.start || exon_iv.1 != transcript_interval.end {
        return Err(FeatureError::SubFeatureIntervalError);
    }

    match exon_coords.len() {
        0 => Err(FeatureError::IncompleteTranscriptError),
        _ => match cds_coord {

            // Improper CDS interval is an error
            Some(coord_cds) if coord_cds.0 > coord_cds.1 =>
                Err(FeatureError::SubFeatureIntervalError),

            // Presence of proper CDS interval means we can resolve UTRs and start/stop codons
            Some(coord_cds) if coord_cds.0 < coord_cds.1 => {
                // CDS coords must be fully enveloped by exon max-min
                if coord_cds.0 < exon_iv.0 || coord_cds.1 > exon_iv.1 {
                    return Err(FeatureError::SubFeatureIntervalError);
                }
                let nodes = exon_coords.iter().map(|&(a, b)| (a..b, ()));
                let tree = IntervalTree::from_iter(nodes);
                let mut features: Vec<TranscriptFeature> = Vec::new();

                // Features at the 5' of cds interval
                for fx5 in tree.find(transcript_interval.start..coord_cds.0) {
                    let end = min(fx5.interval().end, coord_cds.0);
                    let fx_interval = Interval::new(fx5.interval().start..end).unwrap();
                    let utr = match *transcript_strand {
                        Strand::Forward => TxFeature::UTR5,
                        Strand::Reverse => TxFeature::UTR3,
                        Strand::Unknown => TxFeature::UTR,
                    };

                    features.push(TranscriptFeature {
                        seq_name: transcript_seqname.clone(),
                        kind: TxFeature::Exon,
                        interval: fx_interval.clone(),
                        strand: *transcript_strand,
                        attributes: HashMap::new(),
                    });
                    features.push(TranscriptFeature {
                        seq_name: transcript_seqname.clone(),
                        kind: utr,
                        interval: fx_interval.clone(),
                        strand: *transcript_strand,
                        attributes: HashMap::new(),
                    });
                }
                // CDS features
                // TODO: Codons and frames
                for fx in tree.find(coord_cds.0..coord_cds.1) {
                    let start = max(coord_cds.0, fx.interval().start);
                    let end = min(coord_cds.1, fx.interval().end);
                    let fx_interval = Interval::new(start..end).unwrap();

                    features.push(TranscriptFeature {
                        seq_name: transcript_seqname.clone(),
                        kind: TxFeature::Exon,
                        interval: fx_interval.clone(),
                        strand: *transcript_strand,
                        attributes: HashMap::new(),
                    });
                    features.push(TranscriptFeature {
                        seq_name: transcript_seqname.clone(),
                        kind: TxFeature::CDS,
                        interval: fx_interval.clone(),
                        strand: *transcript_strand,
                        attributes: HashMap::new(),
                    });
                }
                // Features at the 3' of cds interval
                for fx3 in tree.find(coord_cds.1..transcript_interval.end) {
                    let start = max(fx3.interval().start, coord_cds.1);
                    let fx_interval = Interval::new(start..fx3.interval().end).unwrap();
                    let utr = match *transcript_strand {
                        Strand::Forward => TxFeature::UTR3,
                        Strand::Reverse => TxFeature::UTR5,
                        Strand::Unknown => TxFeature::UTR,
                    };

                    features.push(TranscriptFeature {
                        seq_name: transcript_seqname.clone(),
                        kind: TxFeature::Exon,
                        interval: fx_interval.clone(),
                        strand: *transcript_strand,
                        attributes: HashMap::new(),
                    });
                    features.push(TranscriptFeature {
                        seq_name: transcript_seqname.clone(),
                        kind: utr,
                        interval: fx_interval.clone(),
                        strand: *transcript_strand,
                        attributes: HashMap::new(),
                    });
                }

                Ok(features)
            },

            // No CDS intervals mean we just sort the coordinates and create the exons
            _ => {
                let mut exon_coords_s = exon_coords.clone();
                exon_coords_s.sort();
                if exon_coords_s.first().unwrap().0 < transcript_interval.start ||
                    exon_coords_s.last().unwrap().1 > transcript_interval.end {
                    return Err(FeatureError::SubFeatureIntervalError)
                }

                exon_coords_s.into_iter()
                    .map(|(start, end)| {
                        TxFeatureBuilder::new(transcript_seqname.clone(), start, end)
                            .kind(TxFeature::Exon)
                            .strand(*transcript_strand)
                            .build()
                    })
                    .collect::<Result<Vec<TranscriptFeature>, FeatureError>>()
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TxFeature {
    Exon,
    UTR,
    UTR5,
    UTR3,
    CDS,
    StartCodon,
    StopCodon,
    Any,
}

#[derive(Debug)]
pub struct TranscriptFeature {
    kind: TxFeature,
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    attributes: HashMap<String, String>
}

impl_annotation!(TranscriptFeature);

impl TranscriptFeature {

    fn kind(&self) -> &TxFeature {
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
}

impl TxFeatureBuilder {

    fn new<T>(seq_name: T, start: u64, end: u64) -> TxFeatureBuilder
        where T: Into<String>
    {
        TxFeatureBuilder {
            seq_name: seq_name.into(),
            start: start, end: end,
            kind: TxFeature::Any,
            attributes: HashMap::new(),
            strand: None,
            strand_char: None,
        }
    }

    fn kind(mut self, kind: TxFeature) -> TxFeatureBuilder {
        self.kind = kind;
        self
    }

    fn strand(mut self, strand: Strand) -> TxFeatureBuilder {
        self.strand = Some(strand);
        self
    }

    fn strand_char(mut self, strand_char: char) -> TxFeatureBuilder {
        self.strand_char = Some(strand_char);
        self
    }

    fn attribute<K, V>(mut self, key: K, value: V) -> TxFeatureBuilder
        where K: Into<String>, V: Into<String>
    {
        self.attributes.insert(key.into(), value.into());
        self
    }

    fn build(self) -> Result<TranscriptFeature, FeatureError> {
        let interval = coords_to_interval(self.start, self.end)?;
        let strand = resolve_strand_input(self.strand, self.strand_char)?;
        let feature = TranscriptFeature {
            seq_name: self.seq_name, kind: self.kind, interval: interval,
            strand: strand, attributes: self.attributes,
        };
        Ok(feature)
    }
}

#[derive(Debug)]
pub struct Transcript {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    attributes: HashMap<String, String>,
    features: Vec<TranscriptFeature>,
}

impl_annotation!(Transcript);

impl Transcript {

    fn features(&self) -> &Vec<TranscriptFeature> {
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
    cds_coord: Option<(u64, u64)>,
    attributes: HashMap<String, String>,
}

impl TranscriptBuilder {

    fn new<T>(seq_name: T, start: u64, end: u64) -> TranscriptBuilder
        where T: Into<String>
    {
        TranscriptBuilder {
            seq_name: seq_name.into(),
            start: start, end: end,
            strand: None,
            strand_char: None,
            features: None,
            exon_coords: None,
            cds_coord: None,
            attributes: HashMap::new(),
        }
    }

    fn strand(mut self, strand: Strand) -> TranscriptBuilder {
        self.strand = Some(strand);
        self
    }

    fn strand_char(mut self, strand_char: char) -> TranscriptBuilder {
        self.strand_char = Some(strand_char);
        self
    }

    fn attribute<K, V>(mut self, key: K, value: V) -> TranscriptBuilder
        where K: Into<String>, V: Into<String>
    {
        self.attributes.insert(key.into(), value.into());
        self
    }

    fn exon_and_cds_coords<E>(mut self, exon_coords: E, cds_coord: Option<(u64, u64)>)-> TranscriptBuilder
        where E: IntoIterator<Item=(u64, u64)>
    {
        self.exon_coords = Some(exon_coords.into_iter().collect());
        self.cds_coord = cds_coord;
        self
    }

    fn build(mut self) -> Result<Transcript, FeatureError> {
        let interval = coords_to_interval(self.start, self.end)?;
        let strand = resolve_strand_input(self.strand, self.strand_char)?;
        let features = resolve_transcript_features(
            &self.seq_name, &interval, &strand,
            self.features, self.exon_coords.as_ref(), self.cds_coord)?;

        let transcript = Transcript {
            seq_name: self.seq_name, interval: interval, strand: strand,
            features: features, attributes: self.attributes,
        };
        Ok(transcript)
    }
}

#[derive(Debug)]
pub struct Gene {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    attributes: HashMap<String, String>,
    transcripts: HashMap<String, Transcript>,
}

impl_annotation!(Gene);


#[cfg(test)]
mod test_transcript_feature {
    use super::*;

    #[test]
    fn builder() {
        let tfm1 = TxFeatureBuilder::new("chrT", 10, 20)
            .kind(TxFeature::Exon)
            .strand(Strand::Forward)
            .attribute("name", "ex1")
            .build();
        assert!(tfm1.is_ok());
        let tf = tfm1.unwrap();
        assert_eq!(tf.seq_name(), "chrT");
        assert_eq!(tf.kind(), &TxFeature::Exon);
        assert_eq!(tf.strand(), &Strand::Forward);
        assert_eq!(tf.attribute("name"), Some("ex1"));
        assert_eq!(tf.attributes.len(), 1);

        let tfm2 = TxFeatureBuilder::new("chrO", 10, 10)
            .strand_char('-')
            .strand(Strand::Reverse)
            .build();
        assert!(tfm2.is_ok());
    }

    #[test]
    fn builder_interval_invalid() {
        let tfm = TxFeatureBuilder::new("chrE", 20, 10).build();
        assert!(tfm.is_err());
        assert_eq!(tfm.unwrap_err(), FeatureError::IntervalError);
    }

    #[test]
    fn builder_strand_unspecified() {
        let tfm = TxFeatureBuilder::new("chrT", 20, 30)
            .build();
        assert!(tfm.is_err());
        assert_eq!(tfm.unwrap_err(), FeatureError::UnspecifiedStrandError);
    }

    #[test]
    fn builder_strand_char_unexpected() {
        let tfm = TxFeatureBuilder::new("chrE", 10, 20)
            .strand_char('w')
            .build();
        assert!(tfm.is_err());
        assert_eq!(tfm.unwrap_err(), FeatureError::StrandCharError);
    }

    #[test]
    fn builder_strand_char_conflicting() {
        let tfm = TxFeatureBuilder::new("chrE", 10, 20)
            .strand_char('-')
            .strand(Strand::Reverse)
            .build();
        assert!(tfm.is_ok());
        let tf = tfm.unwrap();
        assert_eq!(tf.strand(), &Strand::Reverse);
    }
}
