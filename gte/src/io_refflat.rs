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
use std::io;
use std::fs;
use std::path::Path;
use std::str::FromStr;

use csv;
use itertools::{GroupBy, Group, Itertools};
use linked_hash_map::LinkedHashMap;

use {Coord, FeatureError, Gene, GBuilder, Strand, Transcript, TBuilder, Error, consts};


pub type RefFlatRow = (String, String, String, char, u64, u64, u64, u64, usize, String, String);

pub struct RefFlatRecord {
    pub gene_id: String,
    pub transcript_name: String,
    pub seq_name: String,
    pub strand: char,
    pub trx_start: u64,
    pub trx_end: u64,
    pub coding_start: u64,
    pub coding_end: u64,
    pub num_exons: usize,
    pub exon_starts: String,
    pub exon_ends: String,
}

impl RefFlatRecord {

    pub fn into_transcript(self) -> ::Result<Transcript> {

        let exon_coords = self.zip_raw_exon_coords()?;
        if exon_coords.len() != self.num_exons {
            return Err(Error::RefFlat(
                "number of exon and exon coordinates mismatch"));
        }

        let strand = Strand::from_char(&self.strand)
            .map_err(FeatureError::from)
            .map_err(Error::from)?;

        let coding_coord =
            if self.coding_start != self.coding_end {
                Some((self.coding_start, self.coding_end))
            } else {
                None
            };

        TBuilder::new(self.seq_name, self.trx_start, self.trx_end)
            .id(self.transcript_name)
            .gene_id(self.gene_id)
            .strand(strand)
            .coords(exon_coords, coding_coord)
            .coding_incl_stop(true)
            .build()
            .map_err(Error::from)
    }

    #[inline]
    fn zip_raw_exon_coords(&self) -> Result<Vec<Coord<u64>>, Error> {
        let exon_starts = self.exon_starts
            .trim_matches(',').split(',').map(|item| u64::from_str(item).map_err(Error::from));
        let exon_ends = self.exon_ends
            .trim_matches(',').split(',').map(|item| u64::from_str(item).map_err(Error::from));

        let mut res = vec![];
        for (rstart, rend) in exon_starts.zip(exon_ends) {
            let start = rstart?;
            let end = rend?;
            res.push((start, end));
        }
        Ok(res)
    }
}

impl From<RefFlatRow> for RefFlatRecord {

    fn from(row: RefFlatRow) -> Self {
        RefFlatRecord {
            gene_id: row.0,
            transcript_name: row.1,
            seq_name: row.2,
            strand: row.3,
            trx_start: row.4,
            trx_end: row.5,
            coding_start: row.6,
            coding_end: row.7,
            num_exons: row.8,
            exon_starts: row.9,
            exon_ends: row.10,
        }
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

    pub fn records_stream<'a>(&'a mut self) -> RefFlatRecordsStream<'a, R> {
        RefFlatRecordsStream {
            inner: self.inner.decode()
        }
    }

    pub fn transcripts_stream<'a>(&'a mut self) -> RefFlatTranscriptsStream<'a, R> {
        RefFlatTranscriptsStream {
            inner: self.records_stream()
        }
    }

    pub fn genes_stream<'a>(&'a mut self) -> RefFlatGenesStream<'a, R> {
        RefFlatGenesStream {
            inner: self.records_stream().group_by(RefFlatGenesStream::<'a, R>::group_func),
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
}

impl<'a, R> Iterator for RefFlatRecordsStream<'a, R> where R: io::Read {

    type Item = ::Result<RefFlatRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
            .map(|row| {
                row.or_else(|err| Err(Error::from(err)))
                    .map(RefFlatRecord::from)
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
                let (mut gene_start, mut gene_end) = consts::INIT_COORD;
                for record in records {
                    let transcript = record.and_then(|rec| rec.into_transcript())?;
                    gene_start = min(gene_start, transcript.start());
                    gene_end = max(gene_end, transcript.end());
                    let tid = transcript.id().map(|id| id.to_owned())
                        .ok_or(Error::RefFlat("transcript does not have ID"))?;
                    transcripts.insert(tid, transcript);
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
            .map_err(Error::from)
    }

    pub fn write_record(&mut self, record: &RefFlatRecord) -> ::Result<()> {
        self.inner
            .encode((&record.gene_id, &record.transcript_name, &record.seq_name,
                     record.strand, record.trx_start, record.trx_end,
                     record.coding_start, record.coding_end, record.num_exons,
                     &record.exon_starts, &record.exon_ends))
            .map_err(Error::from)
    }

    pub fn write_transcript(&mut self, transcript: &Transcript) -> ::Result<()> {
        let transcript_name = transcript.id()
            .map(|tn| tn.as_ref()).unwrap_or("");
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
            .map_err(Error::from)
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

    pub fn from_file<P: AsRef<Path>>(path: P) -> ::Result<Self> {
        let f = fs::File::create(path).map_err(Error::from)?;
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
