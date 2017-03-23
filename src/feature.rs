//! Features and intervals for gene annotations.

use std::collections::HashMap;

use bio::io::Strand;

use self::error::FeatureError;


pub mod error {

    use std::error::Error;
    use std::fmt;

    #[derive(Debug)]
    pub enum FeatureError {
        InvalidCoords,
    }

    impl Error for FeatureError {
        fn description(&self) -> &str {
            match *self {
                FeatureError::InvalidCoords => "interval start coordinate larger than its end coordinate"
            }
        }

        fn cause(&self) -> Option<&Error> {
            None
        }
    }

    impl fmt::Display for FeatureError {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "FeatureError: {}", self.description())
        }
    }
}

/// Interface for features with start and end coordinates.
///
/// The intervals represented here are zero-based, half open intervals.
/// Coordinates start from zero and start coordinates must be less than
/// or equal to the end coordinates.
pub trait Interval: Sized {

    /// Start coordinate of the interval.
    fn start(&self) -> u64;

    /// End coordinate of the interval.
    fn end(&self) -> u64;

    /// Name of the interval.
    fn name(&self) -> Option<&str>;

    /// Name setter that returns the implementor itself.
    ///
    /// This function is expected to mutate the implementing type.
    fn with_name<T: Into<String>>(self, name: T) -> Self;

    /// Coordinate setter that returns the implementor itself.
    ///
    /// This function is expected to mutate the implementing type.
    fn with_coords(self, start: u64, end: u64) -> Self;

    /// Performs various validation of the interval.
    fn validate(self) -> Result<Self, FeatureError> {
        self.validate_coords()
    }

    /// The number of bases covered by the interval.
    fn span(&self) -> u64 {
        self.end() - self.start()
    }

    /// Whether two intervals have an intersection or not.
    fn overlaps(&self, other: &Self) -> bool {
        (other.start() <= self.start() && self.start() <= other.end()) ||
        (other.start() <= self.end() && self.end() <= other.end())
    }

    /// Whether one interval completely contains the other.
    fn envelops(&self, other: &Self) -> bool {
        self.start() <= other.start() && self.end() >= other.end()
    }

    /// Whether two intervals cover a contiguous region without any overlaps.
    fn adjacent(&self, other: &Self) -> bool {
        self.end() == other.start() || self.start() == other.end()
    }

    /// Performs validation of the interval coordinates.
    ///
    /// If the validation fails, an error message is returned within
    /// the `Result` type.
    fn validate_coords(self) -> Result<Self, FeatureError> {
        if self.start() > self.end() {
            return Err(FeatureError::InvalidCoords)
        }

        Ok(self)
    }
}

/// Macro for default function implementations of interval types.
macro_rules! impl_interval {
    ($struct_ty:ty) => (

        impl Interval for $struct_ty {

            /// Name of the interval.
            fn name(&self) -> Option<&str> {
                self.name.as_ref().map(|n| n.as_str())
            }

            /// Start coordinate of the interval.
            fn start(&self) -> u64 {
                self.start
            }

            /// End coordinate of the interval.
            fn end(&self) -> u64 {
                self.end
            }

            fn with_name<T>(mut self, name: T) -> $struct_ty
                where T: Into<String>
            {
                self.name = Some(name.into());
                self
            }

            fn with_coords(mut self, start: u64, end: u64) -> $struct_ty {
                self.start = start;
                self.end = end;
                self
            }
        }

    );
}

/// Default implementation of the `Interval` trait.
///
/// This struct also provides static methods for creating exons, transcripts, and genes.
#[derive(Debug, Default)]
pub struct Feature {
    start: u64,
    end: u64,
    name: Option<String>,
}

impl Feature {

    /// Creates a gene interval with default values.
    ///
    /// A gene interval is a container for transcript intervals.
    ///
    /// # Examples
    ///
    /// ```
    /// let gene = Feature::gene();
    ///
    /// assert_eq!(gene.transcript().len(), 0);
    /// assert_eq!(gene.start(), 0);
    /// assert_eq!(gene.end(), 0);
    /// assert_eq!(gene.name(), None);
    /// ```
    pub fn gene() -> Gene {
        Gene::default()
    }

    /// Creates a transcript interval with default values.
    ///
    /// A transcript interval is a container for exon intervals.
    ///
    /// # Examples
    ///
    /// ```
    /// use bio::io::Strand;
    ///
    /// let transcript = Feature::transcript();
    ///
    /// assert_eq!(transcript.exons().len(), 0);
    /// assert_eq!(transcript.strand(), &Strand::Unknown)
    /// assert_eq!(transcript.start(), 0);
    /// assert_eq!(transcript.end(), 0);
    /// assert_eq!(transcript.name(), None);
    /// ```
    pub fn transcript() -> Transcript {
        Transcript::default()
    }

    /// Creates an exon interval with default values.
    ///
    /// # Examples
    ///
    /// ```
    /// let exon = Feature::exon();
    ///
    /// assert_eq!(exon.start(), 0);
    /// assert_eq!(exon.end(), 0);
    /// assert_eq!(exon.name(), None);
    /// ```
    pub fn exon() -> Exon {
        Exon::default()
    }

    pub fn with_coords(mut self, start: u64, end: u64) -> Feature {
        self.start = start;
        self.end = end;
        self
    }
}

impl_interval!(Feature);

#[cfg(test)]
mod test_feature {
    use super::*;

    #[test]
    fn default() {
        let fx = Feature::default();
        assert_eq!(fx.start(), 0);
        assert_eq!(fx.end(), 0);
        assert!(fx.validate().is_ok());
    }

    #[test]
    fn with_name() {
        let fx1 = Feature::default()
            .with_name("fx1");
        assert_eq!(fx1.start(), 0);
        assert_eq!(fx1.end(), 0);
        assert_eq!(fx1.name(), Some("fx1"));

        let fx2 = Feature::default()
            .with_name("fx2".to_owned());
        assert_eq!(fx2.start(), 0);
        assert_eq!(fx2.end(), 0);
        assert_eq!(fx2.name(), Some("fx2"));
    }

    #[test]
    fn with_coords() {
        let fx = Feature::default()
            .with_coords(1, 3);
        assert_eq!(fx.start(), 1);
        assert_eq!(fx.end(), 3);
        assert_eq!(fx.name(), None);
    }

    #[test]
    fn with_multiples() {
        let fx = Feature::default()
            .with_name("fx")
            .with_coords(20, 30);
        assert_eq!(fx.start(), 20);
        assert_eq!(fx.end(), 30);
        assert_eq!(fx.name(), Some("fx"));
    }
}

/// Gene annotation.
#[derive(Debug, Default)]
pub struct Gene {
    name: Option<String>,
    start: u64,
    end: u64,
    transcripts: HashMap<String, Transcript>,
}

impl Gene {

    pub fn transcripts(&self) -> &HashMap<String, Transcript> {
        &self.transcripts
    }
}

impl_interval!(Gene);

#[cfg(test)]
mod test_gene {
    use super::*;

    #[test]
    fn default() {
        let gene = Feature::gene();
        assert_eq!(gene.start(), 0);
        assert_eq!(gene.end(), 0);
        assert_eq!(gene.name(), None);
        assert_eq!(gene.transcripts().len(), 0);
    }
}

/// Transcript annotation.
#[derive(Debug)]
pub struct Transcript {
    name: Option<String>,
    strand: Strand,
    start: u64,
    end: u64,
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    exons: Vec<Exon>,
}

impl Transcript {

    pub fn strand(&self) -> &Strand {
        &self.strand
    }

    pub fn with_strand(mut self, strand: Strand) -> Transcript {
        self.strand = strand;
        self
    }

    pub fn cds_start(&self) -> Option<u64> {
        self.cds_start
    }

    pub fn with_cds_start(mut self, cds_start: u64) -> Transcript {
        self.cds_start = Some(cds_start);
        self
    }

    pub fn cds_end(&self) -> Option<u64> {
        self.cds_end
    }

    pub fn with_cds_end(mut self, cds_end: u64) -> Transcript {
        self.cds_end = Some(cds_end);
        self
    }

    pub fn exons(&self) -> &Vec<Exon> {
        &self.exons
    }
}

impl Default for Transcript {

    fn default() -> Transcript {
        Transcript {
            name: None,
            strand: Strand::Unknown,
            start: 0,
            end: 0,
            cds_start: None,
            cds_end: None,
            exons: Vec::new(),
        }
    }
}

impl_interval!(Transcript);

#[cfg(test)]
mod test_transcript {
    use super::*;

    #[test]
    fn default() {
        let trx = Feature::transcript();
        assert_eq!(trx.start(), 0);
        assert_eq!(trx.end(), 0);
        assert_eq!(trx.name(), None);
        assert_eq!(trx.strand(), &Strand::Unknown);
        assert_eq!(trx.cds_start(), None);
        assert_eq!(trx.cds_end(), None);
        assert_eq!(trx.exons().len(), 0);
    }

    #[test]
    fn with_strand() {
        let trx = Feature::transcript()
            .with_strand(Strand::Forward);
        assert_eq!(trx.strand(), &Strand::Forward);
    }

    #[test]
    fn with_cds_start() {
        let trx = Feature::transcript()
            .with_cds_start(20);
        assert_eq!(trx.cds_start(), Some(20));
    }

    #[test]
    fn with_cds_end() {
        let trx = Feature::transcript()
            .with_cds_end(40);
        assert_eq!(trx.cds_end(), Some(40))
    }
}

/// Exon annotation.
#[derive(Debug, Default)]
pub struct Exon {
    name: Option<String>,
    start: u64,
    end: u64,
}

impl_interval!(Exon);

#[cfg(test)]
mod test_exon {
    use super::*;

    #[test]
    fn default() {
        let exon = Feature::exon();
        assert_eq!(exon.start(), 0);
        assert_eq!(exon.end(), 0);
        assert_eq!(exon.name(), None);
    }
}
