//! Features and intervals for gene annotations.

use std::collections::HashMap;

use bio::io::Strand;


/// Builder for exons, transcripts, or genes.
///
/// Use this struct to initialize annotation features with the default values.
pub struct Feature;

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

    /// Performs various validation of the interval.
    fn validate(self) -> Result<Self, &'static str>;

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
    fn validate_coords(self) -> Result<Self, &'static str> {
        if self.start() > self.end() {
            return Err("interval start coordinate larger than its end coordinate")
        }

        Ok(self)
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

    pub fn name(&self) -> Option<&str> {
        self.name.as_ref().map(|n| n.as_str())
    }

    pub fn transcripts(&self) -> &HashMap<String, Transcript> {
        &self.transcripts
    }
}

impl Interval for Gene {

    fn start(&self) -> u64 {
        self.start
    }

    fn end(&self) -> u64 {
        self.end
    }

    fn validate(self) -> Result<Self, &'static str> {
        self.validate_coords()
    }
}

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

impl Transcript {

    pub fn name(&self) -> Option<&str> {
        self.name.as_ref().map(|n| n.as_str())
    }

    pub fn strand(&self) -> &Strand {
        &self.strand
    }

    pub fn cds_start(&self) -> Option<u64> {
        self.cds_start
    }

    pub fn cds_end(&self) -> Option<u64> {
        self.cds_end
    }

    pub fn exons(&self) -> &Vec<Exon> {
        &self.exons
    }
}

impl Interval for Transcript {

    fn start(&self) -> u64 {
        self.start
    }

    fn end(&self) -> u64 {
        self.end
    }

    fn validate(self) -> Result<Self, &'static str> {
        self.validate_coords()
    }
}

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
}

/// Exon annotation.
#[derive(Debug, Default)]
pub struct Exon {
    name: Option<String>,
    start: u64,
    end: u64,
}

impl Exon {

    pub fn name(&self) -> Option<&str> {
        self.name.as_ref().map(|n| n.as_str())
    }

    pub fn with_name<T>(mut self, n: T) -> Exon
        where T: Into<String>
    {
        self.name = Some(n.into());
        self
    }
}

impl Interval for Exon {

    fn start(&self) -> u64 {
        self.start
    }

    fn end(&self) -> u64 {
        self.end
    }

    fn validate(self) -> Result<Self, &'static str> {
        self.validate_coords()
    }
}

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

    #[test]
    fn with_name() {
        let exon1 = Feature::exon().with_name("ex1");
        assert_eq!(exon1.start(), 0);
        assert_eq!(exon1.end(), 0);
        assert_eq!(exon1.name(), Some("ex1"));

        let exon2 = Feature::exon().with_name("ex2".to_owned());
        assert_eq!(exon2.start(), 0);
        assert_eq!(exon2.end(), 0);
        assert_eq!(exon2.name(), Some("ex2"));
    }
}
