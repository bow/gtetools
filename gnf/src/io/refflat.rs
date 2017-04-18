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
        self.inner.next().map(|row_res| {
            let row = row_res.map_err(ParseError::from)?;
            let (gene_id,
                 transcript_name,
                 seq_name,
                 strand_char,
                 trx_start,
                 trx_end,
                 coding_start,
                 coding_end,
                 num_exons,
                 raw_exon_starts,
                 raw_exon_ends): RefFlatRow = row;

            let exon_coords = zip_raw_exon_coords(raw_exon_starts, raw_exon_ends);
            if exon_coords.len() != num_exons {
                return Err(ParseError::FormatSpecific(
                    "number of exon and exon coordinates mismatch"));
            }

            let strand = Strand::from_char(&strand_char)
                .map_err(FeatureError::from)
                .map_err(ParseError::from)?;
            let mut tb = TBuilder::new(seq_name, trx_start, trx_end)
                .id(transcript_name)
                .strand(strand)
                .exon_coords(exon_coords)
                .coding_incl_stop(true)
                .attribute("gene_id", gene_id);

            if coding_start != coding_end {
                tb = tb.coding_coord(coding_start, coding_end)
            }

            tb.build().map_err(ParseError::from)
        })
    }
}

#[inline]
fn zip_raw_exon_coords(starts: String, ends: String) -> Vec<(u64, u64)> {
    let exon_starts = starts
        .split(',').filter_map(|item| u64::from_str(item).ok());
    let exon_ends = ends
        .split(',').filter_map(|item| u64::from_str(item).ok());
    exon_starts.zip(exon_ends).collect()
}
