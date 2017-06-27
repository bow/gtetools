/*! Reader and writer for the refFlat format.

The refFlat format is a transcript-oriented format in which each transcript is denoted in a single
line. It is used most prominently by the
[picard suite tools](https://broadinstitute.github.io/picard/)

A minimum specification of the columns can be found on
[this page](https://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat).
*/
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
use utils::{OptionDeref, update_seq_name};


quick_error! {
    /// Errors that occur when reading or writing refFlat files.
    #[derive(Debug)]
    pub enum RefFlatError {
        /// Occurs when the value of the number of exons column, the number of exon start
        /// coordinates, and/or the number of exon end coordinates are not the same.
        ExonCountMismatch(tid: Option<String>) {
            description("number of exons and number of exon coordinates are not equal")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        /// Indicates a duplicate transcript identifier with the same gene identifier.
        DuplicateTranscriptId(gid: Option<String>) {
            description("gene has multiple transcripts with the same identifier")
            display(self_) -> ("{}, gene ID: {}",
                               self_.description(), gid.as_deref().unwrap_or(DEF_ID))
        }
        /// Occurs when the gene identifier column is empty.
        MissingGeneId {
            description("gene identifier column has no value")
        }
        /// Occurs when the transcript identifier column is empty.
        MissingTranscriptId {
            description("transcript identifier column has no value")
        }
        /// Occurs when any of the exon start or end coordinates is not a valid u64 value.
        InvalidExonCoord(err: ParseIntError, tid: Option<String>) {
            description(err.description())
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
            cause(err)
        }
        /// Errors propagated from the underlying `csv` crate.
        Csv(err: csv::Error) {
            description(err.description())
            from()
            cause(err)
        }
    }
}

/// Raw refFlat row type.
///
/// This type represents the simplest value types that compose a refFlat row. The provided reader
/// does not provide this type, but it is rather the type that is used for creating the refFlat
/// record type.
///
/// Each tuple element represents a refFlat column:
///
/// 1.  gene identifier
/// 2.  transcript identifier
/// 3.  sequence name
/// 4.  strand
/// 5.  transcript 5' coordinate
/// 6.  transcript 3' coordinate
/// 7.  coding region 5' coordinate
/// 8.  coding region 3' coordinate
/// 9.  number of exons
/// 10. exon 5' coordinates (as a comma-separated string)
/// 11. exon 3' coordinates (as a comma-separated string)
///
/// All coordinates are zero-based, half-open.
pub type RefFlatRow = (String, String, String, char, u64, u64, u64, u64, usize, String, String);

/// RefFlat record type.
///
/// This type represents the essential information present in a refFlat record. The main
/// differences between this record type and the row type are:
///
/// * The exon coordinates are represented here as `Vec<u64>`, as opposed to just `String` in the
///   row type.
/// * The number of exon start and end coordinates are guaranteed to be equal in this type.
#[derive(Debug, Clone, PartialEq)]
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

    /// Returns the gene identifier.
    pub fn gene_id(&self) -> &str {
        self.gene_id.as_str()
    }

    /// Sets the gene identifier.
    pub fn set_gene_id<T>(&mut self, gene_id: T)
        where T: Into<String>
    {
        self.gene_id = gene_id.into();
    }

    /// Returns the transcript identifier.
    pub fn transcript_id(&self) -> &str {
        self.transcript_id.as_str()
    }

    /// Sets the transcript identifier.
    pub fn set_transcript_id<T>(&mut self, transcript_id: T)
        where T: Into<String> + AsRef<str>
    {
        self.transcript_id = transcript_id.into();
    }

    /// Returns the sequence name.
    pub fn seq_name(&self) -> &str {
        self.seq_name.as_str()
    }

    /// Sets the sequence name.
    pub fn set_seq_name<T>(&mut self, seq_name: T)
        where T: Into<String>
    {
        self.seq_name = seq_name.into();
    }

    /// Returns the strand.
    pub fn strand(&self) -> char {
        self.strand
    }

    /// Sets the strand.
    pub fn set_strand(&mut self, strand_char: char) {
        self.strand = strand_char;
    }

    /// Returns the genome-wise 5'-most transcript coordinate of the record.
    pub fn transcript_start(&self) -> u64 {
        self.transcript_start
    }

    /// Sets the genome-wise 5'-most transcript coordinate of the record.
    pub fn set_transcript_start(&mut self, coord: u64) {
        self.transcript_start = coord;
    }

    /// Returns the genome-wise 3'-most transcript coordinate of the record.
    pub fn transcript_end(&self) -> u64 {
        self.transcript_end
    }

    /// Sets the genome-wise 3'-most transcript coordinate of the record.
    pub fn set_transcript_end(&mut self, coord: u64) {
        self.transcript_end = coord;
    }

    /// Returns the genome-wise 5'-most coding region coordinate of the record.
    ///
    /// This includes the stop codon coordinate on minus strand records.
    pub fn coding_start(&self) -> u64 {
        self.coding_start
    }

    /// Sets the genome-wise 5'-most coding region coordinate of the record.
    ///
    /// This must include the stop codon coordinate on minus strand records.
    pub fn set_coding_start(&mut self, coord: u64) {
        self.coding_start = coord;
    }

    /// Returns the genome-wise 3'-most coding region coordinate of the record.
    ///
    /// This includes the stop codon coordinate on plus strand records.
    pub fn coding_end(&self) -> u64 {
        self.coding_end
    }

    /// Sets the genome-wise 3'-most coding region coordinate of the record.
    ///
    /// This must include the stop codon coordinate on plus strand records.
    pub fn set_coding_end(&mut self, coord: u64) {
        self.coding_end = coord;
    }

    /// Returns the number of exons contained within the record.
    pub fn num_exons(&self) -> usize {
        self.exon_starts.len() // must be the same as exon_ends
    }

    /// Returns a slice of the genome-wise 5'-most coordinates of exons in the record.
    pub fn exon_starts(&self) -> &[u64] {
        self.exon_starts.as_slice()
    }

    /// Returns a slice of the genome-wise 3'-most coordinates of the exons in the record.
    pub fn exon_ends(&self) -> &[u64] {
        self.exon_ends.as_slice()
    }

    /// Sets the exon coordinates of the record.
    ///
    /// An error type will be returned if the number of coordinates differ.
    pub fn set_exon_coords(&mut self, coord_starts: Vec<u64>, coord_ends: Vec<u64>) -> ::Result<()> {
        if coord_starts.len() != coord_ends.len() {
            let tid = self.transcript_id.clone();
            let err = ::Error::from(RefFlatError::ExonCountMismatch(Some(tid)));
            return Err(err);
        }
        self.exon_starts = coord_starts;
        self.exon_ends = coord_ends;
        Ok(())
    }

    /// Creates a record from a row.
    ///
    /// This method will return an error if:
    /// * any of the exon coordinates are not valid u64 values, or
    /// * the number of exon coordinates and the number of exons column value are not equal
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
            let err = RefFlatError::ExonCountMismatch(Some(row.1.clone()));
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

    /// Transforms the record into a transcript.
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

    /// Parses the given raw coordinate string into a vector of u64s.
    ///
    /// The transcript identifier argument is required for when an error type is returned.
    #[inline]
    fn parse_coords(raw_coords: &str, tid: &str) -> Result<Vec<u64>, RefFlatError> {
        let rcoords = raw_coords
            .trim_matches(',')
            .split(',')
            .map(|item| u64::from_str(item)
                 .map_err(|e| RefFlatError::InvalidExonCoord(e, Some(tid.to_owned()))));

        let mut res = vec![];
        for rcoord in rcoords {
            res.push(rcoord?);
        }
        Ok(res)
    }
}

/// RefFlat reader.
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
    seq_name_prefix: Option<String>,
    seq_name_lstrip: Option<String>,
}

impl<R: io::Read> Reader<R> {

    /// Creates a refFlat reader from another reader.
    pub fn from_reader(in_reader: R) -> Reader<R> {
        Reader {
            inner: csv::Reader::from_reader(in_reader)
                .delimiter(b'\t')
                .has_headers(false),
            seq_name_prefix: None,
            seq_name_lstrip: None,
        }
    }

    /// Sets the reader to add the given prefix to all sequence names.
    pub fn seq_name_prefix<T>(&mut self, prefix: T) -> &mut Self
        where T: Into<String>
    {
        self.seq_name_prefix = Some(prefix.into());
        self
    }

    /// Sets the reader to trim the given string from all sequence names if present at the
    /// beginning.
    pub fn seq_name_lstrip<T>(&mut self, lstrip: T) -> &mut Self
        where T: Into<String>
    {
        self.seq_name_lstrip = Some(lstrip.into());
        self
    }

    /// Creates an iterator of refFlat records.
    pub fn records_stream(&mut self) -> RefFlatRecordsStream<R> {
        RefFlatRecordsStream {
            inner: self.inner.decode(),
            seq_name_prefix: self.seq_name_prefix.as_deref(),
            seq_name_lstrip: self.seq_name_lstrip.as_deref(),
        }
    }

    /// Creates an iterator of transcripts.
    pub fn transcripts_stream(&mut self) -> RefFlatTranscriptsStream<R> {
        RefFlatTranscriptsStream {
            inner: self.records_stream()
        }
    }

    /// Creates an iterator of genes.
    ///
    /// This iterator groups consecutive records based on their gene identifiers into genes.
    pub fn genes_stream(&mut self) -> RefFlatGenesStream<R> {
        RefFlatGenesStream {
            inner: self.records_stream()
                .group_by(RefFlatGenesStream::<R>::group_func),
        }
    }
}

impl Reader<fs::File> {

    /// Creates a refFlat reader that reads from the given path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::from_reader)
    }
}

/// Iterator over refFlat records.
pub struct RefFlatRecordsStream<'a, R: 'a> where R: io::Read {
    inner: csv::DecodedRecords<'a, R, RefFlatRow>,
    seq_name_prefix: Option<&'a str>,
    seq_name_lstrip: Option<&'a str>,
}

impl<'a, R> Iterator for RefFlatRecordsStream<'a, R> where R: io::Read {

    type Item = ::Result<RefFlatRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let lstrip = self.seq_name_lstrip.map(|v| (v, v.len()));
        let prefix = self.seq_name_prefix;
        self.inner.next()
            .map(|row| {
                row
                    .or_else(|err| Err(::Error::from(RefFlatError::from(err))))
                    .map(|mut row| {
                        update_seq_name(&mut row.2, prefix, lstrip);
                        row
                    })
                    .and_then(RefFlatRecord::try_from_row)
            })
    }
}

/// Iterator over transcripts created from refFlat records.
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

/// The type used for grouping records into genes.
///
/// The tuple elements represent gene identifier, sequence name, and strand.
type GroupKey = Option<(String, String, char)>;

/// The type of the function used for creating record-grouping keys for genes.
type GroupFunc = fn(&::Result<RefFlatRecord>) -> GroupKey;

/// The type of the grouped records for creating genes.
type GroupedRecords<'a, 'b, R> = Group<'b, GroupKey, RefFlatRecordsStream<'a, R>, GroupFunc>;

/// Iterator over genes created from refFlat records.
pub struct RefFlatGenesStream<'a, R: 'a> where R: io::Read, {
    inner: GroupBy<GroupKey, RefFlatRecordsStream<'a, R>, GroupFunc>,
}

impl<'a, R> RefFlatGenesStream<'a, R> where R: io::Read {

    /// Creates the group key from the given refFlat record result.
    fn group_func(result: &::Result<RefFlatRecord>) -> GroupKey {
        result.as_ref().ok()
            .map(|ref res| (res.gene_id.clone(), res.seq_name.clone(), res.strand.clone()))
    }

    /// Creates genes from the given grouped records.
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

/// RefFlat writer.
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
}

impl<W: io::Write> Writer<W> {

    /// Creates a refFlat writer from another writer.
    pub fn from_writer(in_writer: W) -> Writer<W> {
        Writer {
            inner: csv::Writer::from_writer(in_writer)
                .delimiter(b'\t')
                .quote_style(csv::QuoteStyle::Never)
        }
    }

    /// Writes the given row.
    pub fn write(&mut self, row: &RefFlatRow) -> ::Result<()> {
        self.inner
            .encode((&row.0, &row.1, &row.2, row.3, row.4, row.5, row.6, row.7, row.8,
                     &row.9, &row.10))
            .map_err(|e| ::Error::from(RefFlatError::from(e)))
    }

    /// Writes the given record.
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

    /// Writes the given transcript as a single row.
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

    /// Writes the given gene as multiple rows.
    pub fn write_gene(&mut self, gene: &Gene) -> ::Result<()> {
        for transcript in gene.transcripts().values() {
            self.write_transcript(&transcript)?;
        }
        Ok(())
    }
}

impl Transcript {

    /// Returns the string values of the exon coordinate columns.
    #[inline(always)]
    fn coords_field(&self) -> (String, String) {
        let mut coord_starts = self.exons().iter().map(|exon| exon.start()).join(",");
        coord_starts.push(',');
        let mut coord_ends = self.exons().iter().map(|exon| exon.end()).join(",");
        coord_ends.push(',');
        (coord_starts, coord_ends)
    }
}

impl Writer<fs::File> {

    /// Creates a refFlat writer that writes to the given path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let f = fs::File::create(path)?;
        Ok(Writer::from_writer(f))
    }
}

impl Writer<Vec<u8>> {

    /// Creates a refFlat writer that writes to an in-memory buffer.
    ///
    /// The initial capacity of the buffer is 64 KiB.
    pub fn from_memory() -> Writer<Vec<u8>> {
        Writer::from_writer(Vec::with_capacity(1024 * 64))
    }

    /// Returns the values of the in-memory buffer as a string.
    pub fn as_string(&mut self) -> &str {
        self.inner.as_string()
    }
}
