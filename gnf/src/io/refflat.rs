//! Reader and writer for the refFlat format.
//!
//! The refFlat format is a transcript-oriented format in which each
//! transcript is denoted in a single line. It is used most prominently
//! by the [picard suite tools](https://broadinstitute.github.io/picard/)
//!
//! A minimum specification of the columns can be found here:
//! https://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat


use std::io;
use std::fs;
use std::convert::AsRef;
use std::path::Path;
use std::str::FromStr;

use csv;

use error::FeatureError;
use {Strand, Transcript, TBuilder};
use super::error::ParseError;


pub type RefFlatRow = (String, String, String, char, u64, u64, u64, u64, usize, String, String);

pub struct RefFlatRecord {
    pub gene_id: String,
    pub transcript_name: String,
    pub seq_name: String,
    pub strand_char: char,
    pub trx_start: u64,
    pub trx_end: u64,
    pub coding_start: u64,
    pub coding_end: u64,
    pub num_exons: usize,
    pub exon_starts: String,
    pub exon_ends: String,
}

impl RefFlatRecord {
    pub fn new(row: RefFlatRow) -> RefFlatRecord {
        RefFlatRecord {
            gene_id: row.0,
            transcript_name: row.1,
            seq_name: row.2,
            strand_char: row.3,
            trx_start: row.4,
            trx_end: row.5,
            coding_start: row.6,
            coding_end: row.7,
            num_exons: row.8,
            exon_starts: row.9,
            exon_ends: row.10,
        }
    }

    pub fn to_transcript(self) -> Result<Transcript, ParseError> {

        let exon_coords = self.zip_raw_exon_coords();
        if exon_coords.len() != self.num_exons {
            return Err(ParseError::FormatSpecific(
                "number of exon and exon coordinates mismatch"));
        }

        let strand = Strand::from_char(&self.strand_char)
            .map_err(FeatureError::from)
            .map_err(ParseError::from)?;
        let mut tb = TBuilder::new(self.seq_name, self.trx_start, self.trx_end)
            .id(self.transcript_name)
            .strand(strand)
            .exon_coords(exon_coords)
            .coding_incl_stop(true)
            .attribute("gene_id", self.gene_id);

        if self.coding_start != self.coding_end {
            tb = tb.coding_coord(self.coding_start, self.coding_end)
        }

        tb.build().map_err(ParseError::from)
    }

    #[inline]
    fn zip_raw_exon_coords(&self) -> Vec<(u64, u64)> {
        let exon_starts = self.exon_starts
            .split(',').filter_map(|item| u64::from_str(item).ok());
        let exon_ends = self.exon_ends
            .split(',').filter_map(|item| u64::from_str(item).ok());
        exon_starts.zip(exon_ends).collect()
    }
}

pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
}

impl<R: io::Read> Reader<R> {

    pub fn from_reader(in_reader: R) -> Reader<R> {
        Reader {
            inner: csv::Reader::from_reader(in_reader)
                .delimiter(b'\t')
                .has_headers(false)
        }
    }

    pub fn transcripts<'a>(&'a mut self) -> RefFlatTranscripts<'a, R> {
        RefFlatTranscripts { inner: self.inner.decode() }
    }
}

impl Reader<fs::File> {
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::from_reader)
    }
}

pub struct RefFlatTranscripts<'a, R: 'a> where R: io::Read {
    inner: csv::DecodedRecords<'a, R, RefFlatRow>,
}

impl<'a, R> Iterator for RefFlatTranscripts<'a, R> where R: io::Read {

    type Item = Result<Transcript, ParseError>;

    fn next(&mut self) -> Option<Result<Transcript, ParseError>> {
        self.inner.next().map(|row| {
            row.or_else(|err| Err(ParseError::from(err)))
                .map(RefFlatRecord::new)
                .and_then(|record| record.to_transcript())
        })
    }
}
