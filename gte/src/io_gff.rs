use std::cmp::{max, min};
use std::convert::AsRef;
use std::io;
use std::iter::Filter;
use std::fs;
use std::path::Path;
use std::vec;

use bio::io::gff::{self, GffType};
use itertools::{GroupBy, Group, Itertools};
use linked_hash_map::LinkedHashMap;
use multimap::MultiMap;
use regex::Regex;

use {Coord, Error, Exon, ExonFeatureKind as EFK, Gene, GBuilder, Strand, TBuilder, Transcript,
     RawTrxCoord};
use consts::*;


quick_error! {
    #[derive(Debug)]
    pub enum GffError {
        MissingGeneId {
            description("gene identifier attribute not found")
        }
        MissingTranscriptId {
            description("transcript identifier attribute not found")
        }
        MultipleTranscriptIds {
            description("more than one 'transcript_id' found")
        }
        StopCodonInCds {
            description("'stop_codon' feature intersects cds")
        }
        MissingTranscript {
            description("no 'transcript' feature present")
        }
        MultipleTranscripts {
            description("multiple 'transcript' features present")
        }
        OrphanStop {
            description("stop codon exists without start codon")
        }
        OrphanStart {
            description("start codon exists without stop codon")
        }
        OrphanCds {
            description("cds exists without start and/or stop codon")
        }
        UnsupportedGffType {
            description("unsupported gff type")
        }
    }
}

pub struct Reader<R: io::Read> {
    inner: gff::Reader<R>,
    pub(crate) gff_type: GffType,
}

impl<R: io::Read> Reader<R> {

    pub fn from_reader(in_reader: R, gff_type: GffType) -> Reader<R> {
        Reader {
            inner: gff::Reader::new(in_reader, gff_type),
            gff_type: gff_type.clone(),
        }
    }

    pub fn genes_stream(&mut self) -> GffGenesStream<R> {
        GffGenesStream {
            inner: self.records_stream()
                .filter(GffGenesStream::<R>::gene_filter_func as GxFilterFunc)
                .group_by(GffGenesStream::<R>::gene_group_func),
        }
    }

    pub fn transcripts_stream(&mut self) -> GffTranscriptsStream<R> {
        GffTranscriptsStream {
            inner: self.records_stream()
                .filter(GffTranscriptsStream::<R>::transcript_filter_func as TrxFilterFunc)
                .group_by(GffTranscriptsStream::<R>::transcript_group_func),
        }
    }

    pub fn transcripts<'a>(
        &mut self,
        gene_id_attr: &'a str,
        transcript_id_attr: &'a str,
        contig_prefix: Option<&'a str>,
        contig_lstrip: Option<&'a str>,
        strict_mode: bool,
    ) -> ::Result<GffTranscripts> {

        let gid_regex = make_gff_id_regex(gene_id_attr, self.gff_type)?;
        let tid_regex = make_gff_id_regex(transcript_id_attr, self.gff_type)?;

        let lstrip = contig_lstrip.map(|v| (v, v.len()));

        let mut parts = Vec::new();
        for result in self.raw_rows_stream() {
            let mut row = result.map_err(Error::from)?;
            if let Some(ref pre) = contig_prefix {
                row.0 = format!("{}{}", pre, row.0);
            }
            if let Some((ref lstr, lstr_len)) = lstrip {
                if row.0.starts_with(lstr) {
                    let _ = row.0.drain(..lstr_len);
                }
            }
            match row.2.as_str() {
                TRANSCRIPT_STR | EXON_STR | CDS_STR | START_CODON_STR | STOP_CODON_STR => {
                    let rf = TranscriptPart::try_from_row(row, &gid_regex, &tid_regex)
                        .map_err(Error::from)?;
                    parts.push(rf);
                },
                _ => {},
            }
        }
        parts.sort_by_key(|ref elem| elem.sort_key());

        Ok(GffTranscripts {
            groups: parts.into_iter().group_by(TranscriptPart::group_key),
            strict_mode: strict_mode,
        })
    }

    pub(crate) fn records_stream(&mut self) -> GffRecords<R> {
        GffRecords {
            inner: self.inner.records()
        }
    }

    pub(crate) fn raw_rows_stream(&mut self) -> GffRawRows<R> {
        GffRawRows {
            inner: self.inner.raw_rows()
        }
    }
}

impl Reader<fs::File> {
    pub fn from_file<P: AsRef<Path>>(path: P, gff_type: GffType) -> io::Result<Self> {
        fs::File::open(path).map(|file| Reader::from_reader(file, gff_type))
    }
}

pub(crate) struct GffRecords<'a, R: 'a> where R: io::Read {
    inner: gff::Records<'a, R>,
}

impl<'a, R> Iterator for GffRecords<'a, R> where R: io::Read {

    type Item = ::Result<gff::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
            .map(|row| row.map_err(Error::from))
    }
}

pub(crate) struct GffRawRows<'a, R: 'a> where R: io::Read {
    inner: gff::RawRows<'a, R>,
}

impl<'a, R> Iterator for GffRawRows<'a, R> where R: io::Read {

    type Item = ::Result<gff::RawRow>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
            .map(|row| row.map_err(Error::from))
    }
}

pub struct GffGenesStream<'a, R: 'a> where R: io::Read, {
    inner: GroupBy<GxGroupKey, GxRecords<'a, R>, GxGroupFunc>,
}

type GxFilterFunc = fn(&::Result<gff::Record>) -> bool;

type GxRecords<'a, R> = Filter<GffRecords<'a, R>, GxFilterFunc>;

type GxGroupKey = Option<(String, String, Strand)>;

type GxGroupFunc = fn(&::Result<gff::Record>) -> GxGroupKey;

type GxGroupedRecs<'a, 'b, R> = Group<'b, GxGroupKey, GxRecords<'a, R>, GxGroupFunc>;

impl<'a, R> GffGenesStream<'a, R> where R: io::Read {

    fn gene_filter_func(result: &::Result<gff::Record>) -> bool {
        match result {
            &Ok(ref rec) => rec.attributes().contains_key(GENE_ID_STR),
            &Err(_) => true,  // Err(_) is supposed to be handled elsewhere
        }
    }

    fn gene_group_func(result: &::Result<gff::Record>) -> GxGroupKey {
        result.as_ref().ok()
            .map(|res| {
                let gene_id = res.attributes().get(GENE_ID_STR)
                    .expect("'gene_id' attribute not found")
                    .clone();
                let seq_name = res.seqname().to_owned();
                let strand = res.strand();

                (seq_name, gene_id, strand)
            })
    }

    fn group_to_gene<'b>(group: (GxGroupKey, GxGroupedRecs<'a, 'b, R>)) -> ::Result<Gene> {
        let (group_key, records) = group;
        match group_key {

            None => Err(records.filter_map(|x| x.err()).next().unwrap()),

            // TODO: Create features of the parsed records instead of just capturing coordinates.
            Some((seq_name, gid, strand)) => {

                let mut gene_coord = INIT_COORD;
                let mut trx_coords: LinkedHashMap<String, RawTrxCoord> = LinkedHashMap::new();

                for record in records {

                    let mut rec = record?;
                    rec.attributes_mut().remove(GENE_ID_STR);

                    let rtrx_entry = rec.attributes_mut()
                        .remove(TRANSCRIPT_ID_STR)
                        .ok_or(GffError::MissingTranscriptId)
                        .and_then(|mut ids| {
                            if ids.len() > 1 {
                                Err(GffError::MultipleTranscriptIds)
                            } else {
                                Ok(ids.pop().unwrap())
                            }
                        })
                        .map(|id| trx_coords.entry(id).or_insert((INIT_COORD, vec![], None)));

                    match rec.feature_type() {
                        Some(GENE_STR) => {
                            gene_coord = adjust_coord(gene_coord, &rec);
                        },
                        Some(TRANSCRIPT_STR) => {
                            let trx_entry = rtrx_entry?;
                            trx_entry.0 = adjust_coord(trx_entry.0, &rec);
                        },
                        Some(EXON_STR) => {
                            let trx_entry = rtrx_entry?;
                            (trx_entry.1).push((rec.start(), rec.end()));
                        },
                        Some(CDS_STR) => {
                            let trx_entry = rtrx_entry?;
                            trx_entry.2 = (trx_entry.2).or(Some(INIT_COORD))
                                .map(|coord| adjust_coord(coord, &rec));
                        },
                        _ => {},
                    }
                }

                GBuilder::new(seq_name, gene_coord.0, gene_coord.1)
                    .id(gid)
                    .strand(strand)
                    .transcript_coords(trx_coords)
                    .transcript_coding_incl_stop(false)
                    .build()
            },
        }
    }
}

impl<'a, R> Iterator for GffGenesStream<'a, R> where R: io::Read {

    type Item = ::Result<Gene>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.into_iter().map(Self::group_to_gene).next()
            .map(|res| res.map_err(Error::from))
    }
}

pub struct GffTranscriptsStream<'a, R: 'a> where R: io::Read {
    inner: GroupBy<TrxGroupKey, TrxRecords<'a, R>, TrxGroupFunc>,
}

type TrxFilterFunc = fn(&::Result<gff::Record>) -> bool;

type TrxRecords<'a, R> = Filter<GffRecords<'a, R>, TrxFilterFunc>;

type TrxGroupKey = Option<(String, String, String, Strand)>;

type TrxGroupFunc = fn(&::Result<gff::Record>) -> TrxGroupKey;

type TrxGroupedRecs<'a, 'b, R> = Group<'b, TrxGroupKey, TrxRecords<'a, R>, TrxGroupFunc>;

impl<'a, R> GffTranscriptsStream<'a, R> where R: io::Read {

    fn transcript_filter_func(result: &::Result<gff::Record>) -> bool {
        match result {
            &Ok(ref rec) =>
                rec.attributes().contains_key(GENE_ID_STR) &&
                rec.attributes().contains_key(TRANSCRIPT_ID_STR),
            &Err(_) => true,
        }
    }

    fn transcript_group_func(result: &::Result<gff::Record>) -> TrxGroupKey {
        result.as_ref().ok()
            .map(|res| {
                let gene_id = res.attributes().get(GENE_ID_STR)
                    .expect("'gene_id' attribute not found")
                    .clone();
                let transcript_id = res.attributes().get(TRANSCRIPT_ID_STR)
                    .expect("'transcript_id' attribute not found")
                    .clone();
                let seq_name = res.seqname().to_owned();
                let strand = res.strand();

                (seq_name, gene_id, transcript_id, strand)
            })
    }

    fn group_to_transcript<'b>(group: (TrxGroupKey, TrxGroupedRecs<'a, 'b, R>)
    ) -> ::Result<Transcript>
    {
        let (group_key, records) = group;
        match group_key {

            None => Err(records.filter_map(|x| x.err()).next().unwrap()),

            // TODO: Create features of the parsed records instead of just capturing coordinates.
            Some((seq_name, gid, tid, strand)) => {

                let mut trx_coord = INIT_COORD;
                let mut exn_coords = Vec::<Coord<u64>>::new();
                let mut coding_coord: Option<Coord<u64>> = None;

                for record in records {

                    let mut rec = record?;
                    rec.attributes_mut().remove(GENE_ID_STR);
                    rec.attributes_mut().remove(TRANSCRIPT_ID_STR);

                    match rec.feature_type() {
                        Some(TRANSCRIPT_STR) => {
                            trx_coord = adjust_coord(trx_coord, &rec);
                        },
                        Some(EXON_STR) => {
                            exn_coords.push((rec.start(), rec.end()));
                        },
                        Some(CDS_STR) => {
                            coding_coord = (coding_coord).or(Some(INIT_COORD))
                                .map(|coord| adjust_coord(coord, &rec));
                        },
                        _ => {},
                    }
                }

                TBuilder::new(seq_name, trx_coord.0, trx_coord.1)
                    .id(tid)
                    .gene_id(gid.clone())
                    .strand(strand)
                    .coords(exn_coords, coding_coord)
                    .coding_incl_stop(false)
                    .build()
            }
        }
    }
}

impl<'a, R> Iterator for GffTranscriptsStream<'a, R> where R: io::Read {

    type Item = ::Result<Transcript>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.into_iter().map(Self::group_to_transcript).next()
            .map(|res| res.map_err(Error::from))
    }
}


#[derive(Debug, PartialEq)]
struct TranscriptPart {
    feature: String,
    chrom: String,
    coord: Coord<u64>,
    strand: Strand,
    transcript_id: String,
    gene_id: String,
}

// sort key: gene identifier, transcript identifier, contig name, start, end, strand num
type TPSortKey = (String, String, String, u64, u64, u8);

// group key: gene identifier, transcript identifier, contig name, strand
type TPGroupKey = (String, String, String, Strand);

impl TranscriptPart {

    fn try_from_row(
        row: gff::RawRow,
        gx_regex: &Regex,
        trx_regex: &Regex,
    ) -> Result<Self, GffError> {

        let gx_id = gx_regex.captures(&row.8)
            .and_then(|cap| cap.name("value"))
            .map(|v| v.as_str().to_owned())
            .ok_or(GffError::MissingGeneId)?;

        let trx_id = trx_regex.captures(&row.8)
            .and_then(|cap| cap.name("value"))
            .map(|v| v.as_str().to_owned())
            .ok_or(GffError::MissingTranscriptId)?;

        Ok(TranscriptPart {
            feature: row.2,
            chrom: row.0,
            coord: (row.3 - 1, row.4),
            strand: Strand::from_char(&row.6).unwrap(),
            transcript_id: trx_id,
            gene_id: gx_id,
        })
    }

    fn sort_key(&self) -> TPSortKey {
        (self.gene_id.clone(), self.transcript_id.clone(),
         self.chrom.clone(), self.coord.0, self.coord.1, self.strand_ord())
    }

    fn group_key(&self) -> TPGroupKey {
        (self.gene_id.clone(), self.transcript_id.clone(), self.chrom.clone(), self.strand)
    }

    fn strand_ord(&self) -> u8 {
        match &self.strand {
            &Strand::Unknown => 0,
            &Strand::Forward => 1,
            &Strand::Reverse => 2,
        }
    }
}

pub struct GffTranscripts {
    groups: GroupBy<TPGroupKey, vec::IntoIter<TranscriptPart>, TPGroupFunc>,
    strict_mode: bool,
}

type TPGroupFunc = fn(&TranscriptPart) -> TPGroupKey;

type TPGroup<'a> = Group<'a, TPGroupKey, vec::IntoIter<TranscriptPart>, TPGroupFunc>;

impl Iterator for GffTranscripts {

    type Item = ::Result<Transcript>;

    fn next(&mut self) -> Option<Self::Item> {
        let group_to_transcript = |(key, rfs): (TPGroupKey, TPGroup)| {
            let (gid, tid, chrom, strand) = key;
            let mut exn_coords = Vec::new();
            let mut trx_coord = None;
            let mut cds_min = None;
            let mut cds_max = None;
            let mut startc_min = None;
            let mut stopc_max = None;

            for rf in rfs {
                match rf.feature.as_str() {
                    TRANSCRIPT_STR => {
                        if let None = trx_coord {
                            trx_coord = Some(rf.coord);
                        } else {
                            let err = Error::Gff(GffError::MultipleTranscripts);
                            return Err(err);
                        }
                    },
                    EXON_STR => {
                        exn_coords.push(rf.coord);
                    },
                    CDS_STR => {
                        cds_min = (cds_min).or(Some(INIT_START)).map(|c| min(c, rf.coord.0));
                        cds_max = (cds_max).or(Some(INIT_END)).map(|c| max(c, rf.coord.1));
                    },
                    START_CODON_STR => {
                        startc_min = (startc_min).or(Some(INIT_START)).map(|c| min(c, rf.coord.0))
                    },
                    STOP_CODON_STR => {
                        stopc_max = (stopc_max).or(Some(INIT_END)).map(|c| max(c, rf.coord.1))
                    },
                    _ => {},
                }
            }

            let (trx_start, trx_end) = trx_coord.ok_or(GffError::MissingTranscript)
                .map_err(Error::from)?;

            let coding_coord = resolve_coding(startc_min, stopc_max, cds_min, cds_max, &strand,
                                              self.strict_mode)
                .map_err(Error::from)?;

            TBuilder::new(chrom, trx_start, trx_end)
                .id(tid)
                .gene_id(gid)
                .strand(strand)
                .coords(exn_coords, coding_coord)
                .coding_incl_stop(true)
                .build()
        };

        self.groups.into_iter().map(group_to_transcript).next()
    }
}

fn resolve_coding(
    startc_min: Option<u64>,
    stopc_max: Option<u64>,
    cds_min: Option<u64>,
    cds_max: Option<u64>,
    strand: &Strand,
    strict_mode: bool,
) -> Result<Option<Coord<u64>>, GffError>
{
    let coding_coord = match (startc_min, stopc_max) {
        // common case: stop and start codon defined
        (Some(c0), Some(c1)) => Some((min(c0, c1), max(c0, c1))),
        // expected case: no stop and start codon defined
        (None, None) => None,
        // error case: only stop or start codon defined
        (a, b) => {
            if strict_mode {
                let start = a.ok_or(GffError::OrphanStop)?;
                let end = b.ok_or(GffError::OrphanStart)?;
                Some((min(start, end), max(start, end)))
            } else {
                let start = a.or(cds_min)
                    .ok_or(GffError::OrphanCds)?;
                let end = b.or(cds_max)
                    .ok_or(GffError::OrphanCds)?;
                Some((min(start, end), max(start, end)))
            }
        }
    };

    if let Some((start, end)) = coding_coord {
        match strand {
            &Strand::Forward if end == cds_max.unwrap() => {
                return Err(GffError::StopCodonInCds);
            },
            &Strand::Reverse if start > cds_min.unwrap() => {
                return Err(GffError::StopCodonInCds);
            },
            _ => {},
        }
    }

    Ok(coding_coord)
}

fn make_gff_id_regex(attr_name: &str, gff_type: GffType) -> ::Result<Regex> {
    let fmts = match gff_type {
        GffType::GFF2 | GffType::GTF2 => Ok((" ", ";", r#"""#)),
        GffType::GFF3 => Ok(("=", ",", "")),
        _ => Err(Error::from(GffError::UnsupportedGffType)),
    };
    fmts.and_then(|(delim, term, nest)| {
        let pat = format!(
            r#"{attr_name}{delim}{nest}(?P<value>[^{delim}{term}\t]+){nest}{term}?"#,
            attr_name=attr_name, delim=delim, term=term, nest=nest);
        Regex::new(&pat).map_err(Error::from)
    })

}

impl Gene {

    #[inline(always)]
    fn num_records(&self) -> usize {
        1 + self.transcripts().values()
            .map(|ref trx| trx.num_records())
            .fold(0, |acc, x| acc + x)
    }

    // TODO: also handle gene-level features
    pub fn into_gff_records(mut self) -> Result<Vec<gff::Record>, Error> {

        let mut attribs = self.take_attributes();

        self.id()
            .ok_or(GffError::MissingGeneId)
            .map(|gid| attribs.insert(GENE_ID_STR.to_owned(), gid.to_owned()))?;

        let (source, score) = extract_source_score(&mut attribs);

        let mut recs = Vec::with_capacity(self.num_records());

        let gx_record = gff::RecordBuilder::new(self.seq_name(), self.start(), self.end())
            .source(source)
            .feature_type(GENE_STR)
            .score(score)
            .strand(strand_to_char(&self.strand()))
            .frame(UNK_CHAR)
            .attributes(attribs)
            .build()?;
        recs.push(gx_record);

        for (_, transcript) in self.take_transcripts() {
            recs.append(&mut transcript.into_gff_records()?);
        }

        Ok(recs)
    }
}

impl Transcript {

    #[inline(always)]
    fn num_records(&self) -> usize {
        1 + self.exons().iter()
            .map(|ref exn| 1 + exn.features().len())
            .fold(0, |acc, x| acc + x)
    }

    // TODO: also handle transcript-level features
    pub fn into_gff_records(mut self) -> Result<Vec<gff::Record>, Error> {

        let mut attribs = self.take_attributes();

        self.gene_id()
            .ok_or(GffError::MissingGeneId)
            .map(|gid| attribs.insert(GENE_ID_STR.to_owned(), gid.to_owned()))?;

        self.id()
            .ok_or(GffError::MissingTranscriptId)
            .map(|tid| attribs.insert(TRANSCRIPT_ID_STR.to_owned(), tid.to_owned()))?;

        let (source, score) = extract_source_score(&mut attribs);

        let mut recs = Vec::with_capacity(self.num_records());

        let trx_record = gff::RecordBuilder::new(self.seq_name(), self.start(), self.end())
            .source(source)
            .feature_type(TRANSCRIPT_STR)
            .score(score)
            .strand(strand_to_char(&self.strand()))
            .frame(UNK_CHAR)
            .attributes(attribs)
            .build()?;
        recs.push(trx_record);

        for exon in self.take_exons() {
            recs.append(&mut exon.into_gff_records()?);
        }

        Ok(recs)
    }
}

impl EFK {

    #[inline(always)]
    fn get_feature_frame(&self) -> (String, char) {
        let (feature, frame) = match self {
            &EFK::UTR => (UTR_STR, UNK_CHAR),
            &EFK::UTR5 => (UTR5_STR, UNK_CHAR),
            &EFK::UTR3 => (UTR3_STR, UNK_CHAR),
            &EFK::CDS { frame: ref f } => (CDS_STR, frame_to_char(f)),
            &EFK::StopCodon { frame: ref f } => (STOP_CODON_STR, frame_to_char(f)),
            &EFK::StartCodon { frame: ref f } => (START_CODON_STR, frame_to_char(f)),
            &EFK::Any(ref s) => (s.as_str(), UNK_CHAR),
        };
        (feature.to_owned(), frame)
    }
}

impl Exon {

    pub fn into_gff_records(mut self) -> Result<Vec<gff::Record>, Error> {

        let mut attribs = self.take_attributes();

        self.gene_id()
            .ok_or(GffError::MissingGeneId)
            .map(|gid| attribs.insert(GENE_ID_STR.to_owned(), gid.to_owned()))?;

        self.transcript_id()
            .ok_or(GffError::MissingTranscriptId)
            .map(|tid| attribs.insert(TRANSCRIPT_ID_STR.to_owned(), tid.to_owned()))?;

        let (source, score) = extract_source_score(&mut attribs);

        let mut recs = Vec::with_capacity(1 + self.features().len());

        for (idx, fx) in self.features().iter().enumerate() {
            let (feature, frame) = fx.kind().get_feature_frame();
            let fx_record = gff::RecordBuilder::new(self.seq_name(), fx.start(), fx.end())
                .source(source.as_str())
                .feature_type(feature.as_str())
                .score(score.as_str())
                .strand(strand_to_char(&self.strand()))
                .frame(frame)
                .attributes(attribs.clone())
                .build()?;

            recs[1 + idx] = fx_record;
        }

        let exn_record = gff::RecordBuilder::new(self.seq_name(), self.start(), self.end())
            .source(source.as_str())
            .feature_type(EXON_STR)
            .score(score.as_str())
            .strand(strand_to_char(&self.strand()))
            .frame(UNK_CHAR)
            .attributes(attribs)
            .build()?;
        recs[0] = exn_record;

        Ok(recs)
    }
}

#[inline(always)]
fn adjust_coord(cur_coord: Coord<u64>, record: &gff::Record) -> Coord<u64> {
    (min(cur_coord.0, record.start()),
     max(cur_coord.1, record.end()))
}

#[inline(always)]
fn extract_source_score(attributes: &mut MultiMap<String, String>) -> (String, String) {
    let source = attributes.remove("source")
        .and_then(|mut vec| vec.pop())
        .unwrap_or(UNK_STR.to_owned());
    let score = attributes.remove("score")
        .and_then(|mut vec| vec.pop())
        .unwrap_or(UNK_STR.to_owned());
    (source, score)
}

#[inline(always)]
fn strand_to_char(strand: &Strand) -> char {
    match strand {
        &Strand::Forward => '+',
        &Strand::Reverse => '-',
        &Strand::Unknown => UNK_CHAR,
    }
}

#[inline(always)]
fn frame_to_char(frame: &Option<u8>) -> char {
    match frame {
        &Some(0) => '0',
        &Some(1) => '1',
        &Some(2) => '2',
        _ => UNK_CHAR,
    }
}
