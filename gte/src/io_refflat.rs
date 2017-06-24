//! Reader and writer for the refFlat format.
//!
//! The refFlat format is a transcript-oriented format in which each
//! transcript is denoted in a single line. It is used most prominently
//! by the [picard suite tools](https://broadinstitute.github.io/picard/)
//!
//! A minimum specification of the columns can be found here:
//! https://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat


use std::cmp::{max, min};
use std::convert::AsRef;
use std::error::Error;
use std::io;
use std::num::ParseIntError;
use std::fs;
use std::path::Path;
use std::str::FromStr;

use csv;
use itertools::{GroupBy, Group, Itertools};
use linked_hash_map::LinkedHashMap;

use {Coord, Gene, GBuilder, Strand, Transcript, TBuilder, DEF_ID, INIT_COORD};
use utils::{OptionDeref, update_contig};


quick_error! {
    #[derive(Debug)]
    pub enum RefFlatError {
        ExonCountMismatch(tid: Option<String>) {
            description("number of exons and number of exon coordinates are not equal")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        ExonCoordsMismatch(tid: Option<String>) {
            description("number of exon start and stop coordinates are not equal")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        DuplicateTranscriptId(gid: Option<String>) {
            description("gene has multiple transcripts with the same identifier")
            display(self_) -> ("{}, gene ID: {}",
                               self_.description(), gid.as_deref().unwrap_or(DEF_ID))
        }
        MissingGeneId {
            description("gene identifier attribute not found")
        }
        MissingTranscriptId {
            description("transcript identifier attribute not found")
        }
        ParseInt(err: ParseIntError, tid: Option<String>) {
            description(err.description())
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
            cause(err)
        }
        Csv(err: csv::Error) {
            description(err.description())
            from()
            cause(err)
        }
    }
}

pub type RefFlatRow = (String, String, String, char, u64, u64, u64, u64, usize, String, String);

#[derive(Debug)]
pub struct RefFlatRecord {
    gene_id: String,
    transcript_id: String,
    seq_name: String,
    strand: char,
    transcript_start: u64,
    transcript_end: u64,
    coding_start: u64,
    coding_end: u64,
    exon_starts: Vec<u64>,
    exon_ends: Vec<u64>,
}

impl RefFlatRecord {

    pub fn gene_id(&self) -> &str {
        self.gene_id.as_str()
    }

    pub fn set_gene_id<T>(&mut self, gene_id: T)
        where T: Into<String>
    {
        self.gene_id = gene_id.into();
    }

    pub fn transcript_id(&self) -> &str {
        self.transcript_id.as_str()
    }

    pub fn set_transcript_id<T>(&mut self, transcript_id: T)
        where T: Into<String> + AsRef<str>
    {
        self.transcript_id = transcript_id.into();
    }

    pub fn seq_name(&self) -> &str {
        self.seq_name.as_str()
    }

    pub fn set_seq_name<T>(&mut self, seq_name: T)
        where T: Into<String>
    {
        self.seq_name = seq_name.into();
    }

    pub fn strand(&self) -> char {
        self.strand
    }

    pub fn set_strand(&mut self, strand_char: char) {
        self.strand = strand_char;
    }

    pub fn transcript_start(&self) -> u64 {
        self.transcript_start
    }

    pub fn set_transcript_start(&mut self, coord: u64) {
        self.transcript_start = coord;
    }

    pub fn transcript_end(&self) -> u64 {
        self.transcript_end
    }

    pub fn set_transcript_end(&mut self, coord: u64) {
        self.transcript_end = coord;
    }

    pub fn coding_start(&self) -> u64 {
        self.coding_start
    }

    pub fn set_coding_start(&mut self, coord: u64) {
        self.coding_start = coord;
    }

    pub fn coding_end(&self) -> u64 {
        self.coding_end
    }

    pub fn set_coding_end(&mut self, coord: u64) {
        self.coding_end = coord;
    }

    pub fn num_exons(&self) -> usize {
        self.exon_starts.len() // must be the same as exon_ends
    }

    pub fn exon_starts(&self) -> &[u64] {
        self.exon_starts.as_slice()
    }

    pub fn exon_ends(&self) -> &[u64] {
        self.exon_ends.as_slice()
    }

    pub fn set_exon_coords(&mut self, starts: Vec<u64>, ends: Vec<u64>) -> ::Result<()> {
        if starts.len() != ends.len() {
            let tid = self.transcript_id.clone();
            let err = ::Error::from(RefFlatError::ExonCountMismatch(Some(tid)));
            return Err(err);
        }
        self.exon_starts = starts;
        self.exon_ends = ends;
        Ok(())
    }

    pub fn try_from_row(row: RefFlatRow) -> ::Result<Self> {

        let exon_starts = Self::parse_coords(row.9.as_str(), row.1.as_str())
            .map_err(::Error::from)?;
        let exon_ends = Self::parse_coords(row.10.as_str(), row.1.as_str())
            .map_err(::Error::from)?;
        if exon_starts.len() != row.8 {
            let err = RefFlatError::ExonCountMismatch(Some(row.1.clone()));
            return Err(::Error::RefFlat(err));
        }
        if exon_starts.len() != exon_ends.len() {
            let err = RefFlatError::ExonCoordsMismatch(Some(row.1.clone()));
            return Err(::Error::RefFlat(err));
        }

        Ok(RefFlatRecord {
            gene_id: row.0,
            transcript_id: row.1,
            seq_name: row.2,
            strand: row.3,
            transcript_start: row.4,
            transcript_end: row.5,
            coding_start: row.6,
            coding_end: row.7,
            exon_starts: exon_starts,
            exon_ends: exon_ends,
        })
    }

    pub fn into_transcript(self) -> ::Result<Transcript> {

        if self.transcript_id.is_empty() {
            return Err(::Error::from(::RefFlatError::MissingTranscriptId));
        }
        if self.gene_id.is_empty() {
            return Err(::Error::from(::RefFlatError::MissingGeneId));
        }
        let coding_interval =
            if self.coding_start == self.coding_end {
                None
            } else {
                Some((self.coding_start, self.coding_end))
            };

        let exon_coords = self.exon_starts.into_iter().zip(self.exon_ends.into_iter())
            .collect::<Vec<Coord<u64>>>();

        TBuilder::new(self.seq_name, self.transcript_start, self.transcript_end)
            .id(self.transcript_id)
            .gene_id(self.gene_id)
            .strand_char(self.strand)
            .coords(exon_coords, coding_interval)
            .coding_incl_stop(true)
            .build()
            .map_err(::Error::from)
    }

    #[inline]
    fn parse_coords(raw_coords: &str, tid: &str) -> Result<Vec<u64>, RefFlatError> {
        let rcoords = raw_coords
            .trim_matches(',')
            .split(',')
            .map(|item| u64::from_str(item)
                 .map_err(|e| RefFlatError::ParseInt(e, Some(tid.to_owned()))));

        let mut res = vec![];
        for rcoord in rcoords {
            res.push(rcoord?);
        }
        Ok(res)
    }
}

pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
    contig_prefix: Option<String>,
    contig_lstrip: Option<String>,
}

impl<R: io::Read> Reader<R> {

    pub fn from_reader(in_reader: R) -> Reader<R> {
        Reader {
            inner: csv::Reader::from_reader(in_reader)
                .delimiter(b'\t')
                .has_headers(false),
            contig_prefix: None,
            contig_lstrip: None,
        }
    }

    pub fn contig_prefix<T>(&mut self, prefix: T) -> &mut Self
        where T: Into<String>
    {
        self.contig_prefix = Some(prefix.into());
        self
    }

    pub fn contig_lstrip<T>(&mut self, lstrip: T) -> &mut Self
        where T: Into<String>
    {
        self.contig_lstrip = Some(lstrip.into());
        self
    }

    pub fn records_stream(&mut self) -> RefFlatRecordsStream<R> {
        RefFlatRecordsStream {
            inner: self.inner.decode(),
            contig_prefix: self.contig_prefix.as_deref(),
            contig_lstrip: self.contig_lstrip.as_deref(),
        }
    }

    pub fn transcripts_stream(&mut self) -> RefFlatTranscriptsStream<R> {
        RefFlatTranscriptsStream {
            inner: self.records_stream()
        }
    }

    pub fn genes_stream(&mut self) -> RefFlatGenesStream<R> {
        RefFlatGenesStream {
            inner: self.records_stream()
                .group_by(RefFlatGenesStream::<R>::group_func),
        }
    }
}

impl Reader<fs::File> {
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::from_reader)
    }
}

pub struct RefFlatRecordsStream<'a, R: 'a> where R: io::Read {
    inner: csv::DecodedRecords<'a, R, RefFlatRow>,
    contig_prefix: Option<&'a str>,
    contig_lstrip: Option<&'a str>,
}

impl<'a, R> Iterator for RefFlatRecordsStream<'a, R> where R: io::Read {

    type Item = ::Result<RefFlatRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let lstrip = self.contig_lstrip.map(|v| (v, v.len()));
        let prefix = self.contig_prefix;
        self.inner.next()
            .map(|row| {
                row
                    .or_else(|err| Err(::Error::from(RefFlatError::from(err))))
                    .map(|mut row| {
                        update_contig(&mut row.2, prefix, lstrip);
                        row
                    })
                    .and_then(RefFlatRecord::try_from_row)
            })
    }
}

pub struct RefFlatTranscriptsStream<'a, R: 'a> where R: io::Read {
    inner: RefFlatRecordsStream<'a, R>,
}

impl<'a, R> Iterator for RefFlatTranscriptsStream<'a, R> where R: io::Read {

    type Item = ::Result<Transcript>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
            .map(|record| record.and_then(|rec| rec.into_transcript()))
    }
}

type GroupKey = Option<(String, String, char)>;

type GroupFunc = fn(&::Result<RefFlatRecord>) -> GroupKey;

type GroupedRecords<'a, 'b, R> = Group<'b, GroupKey, RefFlatRecordsStream<'a, R>, GroupFunc>;

pub struct RefFlatGenesStream<'a, R: 'a> where R: io::Read, {
    inner: GroupBy<GroupKey, RefFlatRecordsStream<'a, R>, GroupFunc>,
}

impl<'a, R> RefFlatGenesStream<'a, R> where R: io::Read {

    fn group_func(result: &::Result<RefFlatRecord>) -> GroupKey {
        result.as_ref().ok()
            .map(|ref res| (res.gene_id.clone(), res.seq_name.clone(), res.strand.clone()))
    }

    fn group_to_gene<'b>(group: (GroupKey, GroupedRecords<'a, 'b, R>)) -> ::Result<Gene> {
        let (group_key, records) = group;
        match group_key {

            None => Err(records.filter_map(|x| x.err()).next().unwrap()),

            Some((gid, seq_name, strand_char)) => {
                let mut transcripts = LinkedHashMap::new();
                let (mut gene_start, mut gene_end) = INIT_COORD;
                for record in records {
                    let transcript = record.and_then(|rec| rec.into_transcript())?;
                    gene_start = min(gene_start, transcript.start());
                    gene_end = max(gene_end, transcript.end());
                    let tid = transcript.id().map(|id| id.to_owned())
                        .ok_or(::Error::from(RefFlatError::MissingTranscriptId))?;
                    let existing_trx = transcripts.insert(tid, transcript);
                    if existing_trx.is_some() {
                        let err = RefFlatError::DuplicateTranscriptId(Some(gid));
                        return Err(::Error::from(err));
                    }
                }
                GBuilder::new(seq_name, gene_start, gene_end)
                    .id(gid)
                    .strand_char(strand_char)
                    .transcripts(transcripts)
                    .transcript_coding_incl_stop(true)
                    .build()
            },
        }
    }
}

impl<'a, R> Iterator for RefFlatGenesStream<'a, R> where R: io::Read {

    type Item = ::Result<Gene>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.into_iter().map(Self::group_to_gene).next()
    }
}

pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
}

impl<W: io::Write> Writer<W> {

    pub fn from_writer(in_writer: W) -> Writer<W> {
        Writer {
            inner: csv::Writer::from_writer(in_writer)
                .delimiter(b'\t')
                .quote_style(csv::QuoteStyle::Never)
        }
    }

    pub fn write(&mut self, row: &RefFlatRow) -> ::Result<()> {
        self.inner
            .encode((&row.0, &row.1, &row.2, row.3, row.4, row.5, row.6, row.7, row.8,
                     &row.9, &row.10))
            .map_err(|e| ::Error::from(RefFlatError::from(e)))
    }

    pub fn write_record(&mut self, record: &RefFlatRecord) -> ::Result<()> {
        let mut exon_starts = record.exon_starts.iter().join(",");
        exon_starts.push(',');
        let mut exon_ends = record.exon_ends.iter().join(",");
        exon_ends.push(',');
        self.inner
            .encode((&record.gene_id, &record.transcript_id, &record.seq_name,
                     record.strand, record.transcript_start, record.transcript_end,
                     record.coding_start, record.coding_end, record.num_exons(),
                     exon_starts, exon_ends))
            .map_err(|e| ::Error::from(RefFlatError::from(e)))
    }

    pub fn write_transcript(&mut self, transcript: &Transcript) -> ::Result<()> {
        let transcript_name = transcript.id()
            .ok_or(::Error::RefFlat(RefFlatError::MissingTranscriptId))?;
        let strand_char = match transcript.strand() {
            &Strand::Forward => '+',
            &Strand::Reverse => '-',
            &Strand::Unknown => '.',
        };

        let (coding_start, coding_end) = transcript.coding_coord(true)
            .unwrap_or((transcript.end(), transcript.end()));
        let (exon_starts, exon_ends) = transcript.coords_field();

        self.inner
            .encode((transcript.gene_id(), transcript_name, transcript.seq_name(), strand_char,
                     transcript.start(), transcript.end(),
                     coding_start, coding_end, transcript.exons().len(),
                     exon_starts, exon_ends))
            .map_err(|e| ::Error::from(RefFlatError::from(e)))
    }

    pub fn write_gene(&mut self, gene: &Gene) -> ::Result<()> {
        for transcript in gene.transcripts().values() {
            self.write_transcript(&transcript)?;
        }
        Ok(())
    }
}

impl Transcript {

    #[inline(always)]
    fn coords_field(&self) -> (String, String) {
        let mut starts = self.exons().iter().map(|exon| exon.start()).join(",");
        starts.push(',');
        let mut ends = self.exons().iter().map(|exon| exon.end()).join(",");
        ends.push(',');
        (starts, ends)
    }
}

impl Writer<fs::File> {

    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let f = fs::File::create(path)?;
        Ok(Writer::from_writer(f))
    }
}

impl Writer<Vec<u8>> {

    pub fn from_memory() -> Writer<Vec<u8>> {
        Writer::from_writer(Vec::with_capacity(1024 * 64))
    }

    pub fn as_string(&mut self) -> &str {
        self.inner.as_string()
    }
}
