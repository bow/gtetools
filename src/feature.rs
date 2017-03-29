//! Interval-based annotation features.

use std::error::Error;
use std::collections::HashMap;
use std::fmt::{self, Display};

use bio::utils::{Interval, IntervalError};

use self::error::FeatureError;

// TODO: PR for rust-bio so that it implements Copy
//       for now, we use our own
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown
}


pub mod error {

    use super::*;

    #[derive(Debug, PartialEq)]
    pub enum FeatureError {
        IntervalError,
        StrandCharError,
        ConflictingStrandError,
        UnspecifiedStrandError,
    }

    impl Error for FeatureError {

        fn description(&self) -> &str {
            match *self {
                FeatureError::IntervalError => "interval start coordinate larger than its end coordinate",
                FeatureError::StrandCharError => "strand character is invalid",
                FeatureError::ConflictingStrandError => "conflicting strand inputs specified",
                FeatureError::UnspecifiedStrandError => "strand not specified",
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
    Interval::new(start..end)
        .map_err(|err| match err {
            IntervalError::InvalidRange => FeatureError::IntervalError
        })
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
