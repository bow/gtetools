use std::cmp::{max, min};
use std::convert::AsRef;
use std::error::Error;
use std::io;
use std::fs;
use std::path::Path;
use std::vec;

use bio::io::gff::{self, GffType};
use itertools::{GroupBy, Group, Itertools};
use multimap::MultiMap;
use regex::Regex;

use {Coord, Exon, ExonFeatureKind as EFK, Gene, Strand, TBuilder, Transcript,
     RawTrxCoords};
use consts::*;
use utils::{OptionDeref, update_contig};


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
        StopCodonInCds(tid: Option<String>) {
            description("'stop_codon' feature intersects cds")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        MissingTranscript(tid: Option<String>) {
            description("no 'transcript' feature present")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        MultipleTranscripts {
            description("multiple 'transcript' features present")
        }
        OrphanStop(tid: Option<String>) {
            description("stop codon exists without start codon")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        OrphanStart(tid: Option<String>) {
            description("start codon exists without stop codon")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        OrphanCodon(tid: Option<String>) {
            description("start and stop codon exists without cds")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        OrphanCds(tid: Option<String>) {
            description("cds exists without start and/or stop codon")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        UnsupportedGffType {
            description("unsupported gff type")
        }
        Bio(err: gff::GffError) {
            description(err.description())
            display("{}", err)
            from()
            cause(err)
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

    pub fn transcripts<'a>(
        &mut self,
        gene_id_attr: Option<&'a str>,
        transcript_id_attr: Option<&'a str>,
        contig_prefix: Option<&'a str>,
        contig_lstrip: Option<&'a str>,
        loose_codons: bool,
    ) -> ::Result<GffTranscripts> {

        let gid_regex = make_gff_id_regex(
            gene_id_attr.unwrap_or(GENE_ID_STR), self.gff_type)?;
        let tid_regex = make_gff_id_regex(
            transcript_id_attr.unwrap_or(TRANSCRIPT_ID_STR), self.gff_type)?;

        let lstrip = contig_lstrip.map(|v| (v, v.len()));

        let mut parts = Vec::new();
        for result in self.raw_rows_stream() {
            let mut row = result.map_err(::Error::from)?;
            update_contig(&mut row.0, contig_prefix, lstrip);
            match row.2.as_str() {
                TRANSCRIPT_STR | EXON_STR | CDS_STR | START_CODON_STR | STOP_CODON_STR => {
                    let rf = TrxPart::try_from_row(row, &gid_regex, &tid_regex)
                        .map_err(::Error::from)?;
                    parts.push(rf);
                },
                _ => {},
            }
        }
        parts.sort_by_key(|ref elem| elem.sort_key());

        Ok(GffTranscripts {
            groups: parts.into_iter().group_by(TrxPart::transcript_group_key),
            loose_codons: loose_codons,
        })
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

pub(crate) struct GffRawRows<'a, R: 'a> where R: io::Read {
    inner: gff::RawRows<'a, R>,
}

impl<'a, R> Iterator for GffRawRows<'a, R> where R: io::Read {

    type Item = ::Result<gff::RawRow>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
            .map(|row| row.map_err(|e| ::Error::from(GffError::from(e))))
    }
}

#[derive(Debug, PartialEq)]
struct TrxPart {
    feature: String,
    chrom: String,
    coord: Coord<u64>,
    strand: Strand,
    transcript_id: String,
    gene_id: String,
}

// sort key: gene identifier, transcript identifier, contig name, start, end, strand num
type TrxSortKey = (String, String, String, u64, u64, u8);

impl TrxPart {

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

        Ok(TrxPart {
            feature: row.2,
            chrom: row.0,
            coord: (row.3 - 1, row.4),
            strand: Strand::from_char(&row.6).unwrap(),
            transcript_id: trx_id,
            gene_id: gx_id,
        })
    }

    fn sort_key(&self) -> TrxSortKey {
        (self.gene_id.clone(), self.transcript_id.clone(),
         self.chrom.clone(), self.coord.0, self.coord.1, self.strand_ord())
    }

    fn transcript_group_key(&self) -> TrxGroupKey {
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

#[derive(Debug, Default)]
struct TrxCoords {
    trx_coord: Option<Coord<u64>>,
    exon_coords: Vec<Coord<u64>>,
    cds_coord: Option<Coord<u64>>,
    codon_5: Option<u64>,
    codon_3: Option<u64>,
}

impl TrxCoords {

    fn set_trx_coord(&mut self, coord: Coord<u64>) -> Result<(), GffError> {
        if let None = self.trx_coord {
            self.trx_coord = Some(coord);
        } else {
            return Err(GffError::MultipleTranscripts);
        }
        Ok(())
    }

    fn add_exon_coord(&mut self, coord: Coord<u64>) {
        self.exon_coords.push(coord);
    }

    fn include_cds_coord(&mut self, coord: Coord<u64>) {
        self.cds_coord = (self.cds_coord).or(Some(INIT_COORD))
            .map(|(a, b)| (min(a, coord.0), max(b, coord.1)));
    }

    fn include_codon_5(&mut self, coord_5: u64) {
        self.codon_5 = (self.codon_5).or(Some(INIT_START))
            .map(|c| min(c, coord_5));
    }

    fn include_codon_3(&mut self, coord_3: u64) {
        self.codon_3 = (self.codon_3).or(Some(INIT_END))
            .map(|c| max(c, coord_3));
    }

    fn resolve<'a>(
        self,
        strand: Strand,
        loose_codons: bool,
        tid: Option<&'a str>
    ) -> Result<RawTrxCoords, GffError> {

        let trx_coord = self.trx_coord
            .ok_or(GffError::MissingTranscript(tid.map(|v| v.to_owned())))?;

        let coding_coord = match (self.codon_5, self.codon_3) {
            // common case: stop and start codon defined
            (Some(c5), Some(c3)) => Some((c5, c3)),
            // expected case: no stop and start codon defined
            (None, None) => None,
            // error case: only stop or start codon defined
            (a, b) => {
                if !loose_codons {
                    let start = a
                        .ok_or(GffError::OrphanStop(tid.map(|v| v.to_owned())))?;
                    let end = b
                        .ok_or(GffError::OrphanStart(tid.map(|v| v.to_owned())))?;
                    Some((start, end))
                } else {
                    let start = a.or(self.cds_coord.map(|c| c.0))
                        .ok_or(GffError::OrphanCds(tid.map(|v| v.to_owned())))?;
                    let end = b.or(self.cds_coord.map(|c| c.1))
                        .ok_or(GffError::OrphanCds(tid.map(|v| v.to_owned())))?;
                    Some((start, end))
                }
            }
        };

        if !loose_codons {
            if let Some((start, end)) = coding_coord {
                let cdsc = self.cds_coord
                    .ok_or(GffError::OrphanCodon(tid.map(|v| v.to_owned())))?;
                match strand {
                    Strand::Forward if end == cdsc.1 => {
                        return Err(GffError::StopCodonInCds(tid.map(|v| v.to_owned())));
                    },
                    Strand::Reverse if start > cdsc.0 => {
                        return Err(GffError::StopCodonInCds(tid.map(|v| v.to_owned())));
                    },
                    _ => {},
                }
            }
        }

        Ok((trx_coord, self.exon_coords, coding_coord))
    }
}

pub struct GffTranscripts {
    groups: GroupBy<TrxGroupKey, vec::IntoIter<TrxPart>, TrxGroupFunc>,
    loose_codons: bool,
}

type TrxGroupKey = (String, String, String, Strand);

type TrxGroupFunc = fn(&TrxPart) -> TrxGroupKey;

type TrxGroup<'a> = Group<'a, TrxGroupKey, vec::IntoIter<TrxPart>, TrxGroupFunc>;

impl Iterator for GffTranscripts {

    type Item = ::Result<Transcript>;

    fn next(&mut self) -> Option<Self::Item> {
        let group_to_transcript = |(key, tps): (TrxGroupKey, TrxGroup)| {
            let (gid, tid, chrom, strand) = key;
            let mut tc = TrxCoords::default();

            for tp in tps {
                match (tp.feature.as_str(), strand) {
                    (TRANSCRIPT_STR, _) => {
                        tc.set_trx_coord(tp.coord)
                            .map_err(::Error::from)?;
                    },
                    (EXON_STR, _) => {
                        tc.add_exon_coord(tp.coord);
                    },
                    (CDS_STR, _) => {
                        tc.include_cds_coord(tp.coord);
                    },
                    (START_CODON_STR, Strand::Forward) | (STOP_CODON_STR, Strand::Reverse) => {
                        tc.include_codon_5(tp.coord.0);
                    },
                    (STOP_CODON_STR, Strand::Forward) | (START_CODON_STR, Strand::Reverse) => {
                        tc.include_codon_3(tp.coord.1);
                    },
                    _ => {},
                }
            }

            let ((trx_start, trx_end), exn_coords, coding_coord) =
                tc.resolve(strand, self.loose_codons, Some(tid.as_str()))
                    .map_err(::Error::from)?;

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

fn make_gff_id_regex(attr_name: &str, gff_type: GffType) -> ::Result<Regex> {
    let fmts = match gff_type {
        GffType::GFF2 | GffType::GTF2 => Ok((" ", ";", r#"""#)),
        GffType::GFF3 => Ok(("=", ",", "")),
        _ => Err(::Error::from(GffError::UnsupportedGffType)),
    };
    fmts.and_then(|(delim, term, nest)| {
        let pat = format!(
            r#"{attr_name}{delim}{nest}(?P<value>[^{delim}{term}\t]+){nest}{term}?"#,
            attr_name=attr_name, delim=delim, term=term, nest=nest);
        Regex::new(&pat).map_err(::Error::from)
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
    pub fn into_gff_records(mut self) -> ::Result<Vec<gff::Record>> {

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
            .build()
            .map_err(|e| ::Error::from(GffError::from(e)))?;
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
    pub fn into_gff_records(mut self) -> ::Result<Vec<gff::Record>> {

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
            .build()
            .map_err(|e| ::Error::from(GffError::from(e)))?;
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

    pub fn into_gff_records(mut self) -> ::Result<Vec<gff::Record>> {

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
                .build()
                .map_err(|e| ::Error::from(GffError::from(e)))?;

            recs[1 + idx] = fx_record;
        }

        let exn_record = gff::RecordBuilder::new(self.seq_name(), self.start(), self.end())
            .source(source.as_str())
            .feature_type(EXON_STR)
            .score(score.as_str())
            .strand(strand_to_char(&self.strand()))
            .frame(UNK_CHAR)
            .attributes(attribs)
            .build()
            .map_err(|e| ::Error::from(GffError::from(e)))?;
        recs[0] = exn_record;

        Ok(recs)
    }
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
