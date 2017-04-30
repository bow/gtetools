//! Reader and writer for the refFlat format.
//!
//! The refFlat format is a transcript-oriented format in which each
//! transcript is denoted in a single line. It is used most prominently
//! by the [picard suite tools](https://broadinstitute.github.io/picard/)
//!
//! A minimum specification of the columns can be found here:
//! https://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat


use std::cmp::{max, min};
use std::collections::HashMap;
use std::convert::AsRef;
use std::io;
use std::fs;
use std::path::Path;
use std::str::FromStr;

use csv;
use itertools::{GroupBy, Itertools};

use feature::FeatureError;
use {Coord, Gene, GBuilder, Strand, Transcript, TBuilder, Error};


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

    pub fn to_transcript(self) -> Result<Transcript, Error> {

        let exon_coords = self.zip_raw_exon_coords();
        if exon_coords.len() != self.num_exons {
            return Err(Error::RefFlat(
                "number of exon and exon coordinates mismatch"));
        }

        let strand = Strand::from_char(&self.strand_char)
            .map_err(FeatureError::from)
            .map_err(Error::from)?;

        let coding_coord =
            if self.coding_start != self.coding_end {
                Some((self.coding_start, self.coding_end))
            } else {
                None
            };

        let mut attribs = HashMap::new();
        attribs.insert("gene_id".to_owned(), self.gene_id);
        TBuilder::new(self.seq_name, self.trx_start, self.trx_end)
            .id(self.transcript_name)
            .strand(strand)
            .coords(exon_coords, coding_coord)
            .coding_incl_stop(true)
            .attributes(attribs)
            .build()
            .map_err(Error::from)
    }

    #[inline]
    fn zip_raw_exon_coords(&self) -> Vec<Coord<u64>> {
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

    pub fn records<'a>(&'a mut self) -> RefFlatRecords<'a, R> {
        RefFlatRecords {
            inner: self.inner.decode()
        }
    }

    pub fn transcripts<'a>(&'a mut self) -> RefFlatTranscripts<'a, R> {
        RefFlatTranscripts {
            inner: self.records()
        }
    }

    pub fn genes<'a>(&'a mut self) -> RefFlatGenes<'a, R> {
        RefFlatGenes {
            inner: self.records().group_by(groupf),
        }
    }
}

impl Reader<fs::File> {
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::from_reader)
    }
}

pub struct RefFlatRecords<'a, R: 'a> where R: io::Read {
    inner: csv::DecodedRecords<'a, R, RefFlatRow>,
}

impl<'a, R> Iterator for RefFlatRecords<'a, R> where R: io::Read {

    type Item = Result<RefFlatRecord, Error>;

    fn next(&mut self) -> Option<Result<RefFlatRecord, Error>> {
        self.inner.next()
            .map(|row| {
                row.or_else(|err| Err(Error::from(err)))
                    .map(RefFlatRecord::new)
            })
    }
}

pub struct RefFlatTranscripts<'a, R: 'a> where R: io::Read {
    inner: RefFlatRecords<'a, R>,
}

impl<'a, R> Iterator for RefFlatTranscripts<'a, R> where R: io::Read {

    type Item = Result<Transcript, Error>;

    fn next(&mut self) -> Option<Result<Transcript, Error>> {
        self.inner.next()
            .map(|record| record.and_then(|rec| rec.to_transcript()))
    }
}

pub struct RefFlatGenes<'a, R: 'a> where R: io::Read, {
    inner: GroupBy<GroupKey, RefFlatRecords<'a, R>, GroupFunc>,
}

impl<'a, R> Iterator for RefFlatGenes<'a, R> where R: io::Read {

    type Item = Result<Gene, Error>;

    fn next(&mut self) -> Option<Result<Gene, Error>> {
        self.inner.into_iter()
            .map(|(okey, records)| {
                match okey {
                    Some((gid, seq_name, strand_char)) => {
                        let mut transcripts = HashMap::new();
                        let (mut gene_start, mut gene_end) = (u64::max_value(), u64::min_value());
                        for record in records {
                            let transcript = record.and_then(|rec| rec.to_transcript())?;
                            gene_start = min(gene_start, transcript.interval().start);
                            gene_end = max(gene_end, transcript.interval().end);
                            let tid = transcript.id.clone()
                                .ok_or(Error::RefFlat("transcript does not have ID"))?;
                            transcripts.insert(tid, transcript);
                        }
                        GBuilder::new(seq_name, gene_start, gene_end)
                            .id(gid)
                            .strand_char(strand_char)
                            .transcripts(transcripts)
                            .build()
                    },
                    None => Err(records.filter_map(|x| x.err()).next().unwrap())
                }
            }).next()
    }
}

type GroupKey = Option<(String, String, char)>;

type GroupFunc = fn(&Result<RefFlatRecord, Error>) -> GroupKey;

fn groupf(result: &Result<RefFlatRecord, Error>) -> GroupKey {
    result.as_ref().ok()
        .map(|ref res| (res.gene_id.clone(), res.seq_name.clone(), res.strand_char.clone()))
}
