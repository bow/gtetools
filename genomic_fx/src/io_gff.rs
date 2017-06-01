use std::collections::HashMap;
use std::cmp::{max, min};
use std::convert::AsRef;
use std::io;
use std::iter::Filter;
use std::fs;
use std::path::Path;

use bio::io::gff::{self, GffType};
use itertools::{GroupBy, Group, Itertools};
use linked_hash_map::LinkedHashMap;

use {Coord, Exon, ExonFeatureKind as EFK, Gene, GBuilder, Strand, Transcript, TBuilder, Error,
     RawTrxCoord};

// Various commonly-used feature column values
const GENE_STR: &'static str = "gene";
const TRANSCRIPT_STR: &'static str = "transcript";
const EXON_STR: &'static str = "exon";
const UTR_STR: &'static str = "UTR";
const UTR5_STR: &'static str = "UTR5";
const UTR3_STR: &'static str = "UTR3";
const CDS_STR: &'static str = "CDS";
const START_CODON_STR: &'static str = "start_codon";
const STOP_CODON_STR: &'static str = "stop_codon";

// Value for unknown columns.
const UNK_STR: &'static str = ".";

// Commonly-used attribute keys.
const GENE_ID_STR: &'static str = "gene_id";
const TRANSCRIPT_ID_STR: &'static str = "transcript_id";

// Initial coordinate for features.
const INIT_COORD: (u64, u64) = (::std::u64::MAX, ::std::u64::MIN);


pub struct Reader<R: io::Read> {
    inner: gff::Reader<R>,
}

impl<R: io::Read> Reader<R> {

    pub fn from_reader(in_reader: R, gff_type: GffType) -> Reader<R> {
        Reader {
            inner: gff::Reader::new(in_reader, gff_type),
        }
    }

    pub fn genes(&mut self) -> GffGenes<R> {
        GffGenes {
            inner: self.records()
                .filter(GffGenes::<R>::gene_filter_func as GxFilterFunc)
                .group_by(GffGenes::<R>::gene_group_func),
        }
    }

    pub fn transcripts(&mut self) -> GffTranscripts<R> {
        GffTranscripts {
            inner: self.records()
                .filter(GffTranscripts::<R>::transcript_filter_func as TrxFilterFunc)
                .group_by(GffTranscripts::<R>::transcript_group_func),
        }
    }

    fn records(&mut self) -> GffRecords<R> {
        GffRecords {
            inner: self.inner.records()
        }
    }
}

impl Reader<fs::File> {
    pub fn from_file<P: AsRef<Path>>(path: P, gff_type: GffType) -> io::Result<Self> {
        fs::File::open(path).map(|file| Reader::from_reader(file, gff_type))
    }
}

struct GffRecords<'a, R: 'a> where R: io::Read {
    inner: gff::Records<'a, R>,
}

impl<'a, R> Iterator for GffRecords<'a, R> where R: io::Read {

    type Item = Result<gff::Record, Error>;

    fn next(&mut self) -> Option<Result<gff::Record, Error>> {
        self.inner.next()
            .map(|row| row.map_err(Error::from))
    }
}

pub struct GffGenes<'a, R: 'a> where R: io::Read, {
    inner: GroupBy<GxGroupKey, GxRecords<'a, R>, GxGroupFunc>,
}

type GxFilterFunc = fn(&Result<gff::Record, Error>) -> bool;

type GxRecords<'a, R> = Filter<GffRecords<'a, R>, GxFilterFunc>;

type GxGroupKey = Option<(String, String, Strand)>;

type GxGroupFunc = fn(&Result<gff::Record, Error>) -> GxGroupKey;

type GxGroupedRecs<'a, 'b, R> = Group<'b, GxGroupKey, GxRecords<'a, R>, GxGroupFunc>;

impl<'a, R> GffGenes<'a, R> where R: io::Read {

    fn gene_filter_func(result: &Result<gff::Record, Error>) -> bool {
        match result {
            &Ok(ref rec) => rec.attributes().contains_key(GENE_ID_STR),
            &Err(_) => true,  // Err(_) is supposed to be handled elsewhere
        }
    }

    fn gene_group_func(result: &Result<gff::Record, Error>) -> GxGroupKey {
        result.as_ref().ok()
            .map(|res| {
                let gene_id = res.attributes().get(GENE_ID_STR)
                    .expect("'gene_id' attribute not found")
                    .clone();
                let seq_name = res.seqname().to_owned();
                let strand = res.strand().unwrap_or(Strand::Unknown);

                (seq_name, gene_id, strand)
            })
    }

    fn group_to_gene<'b>(group: (GxGroupKey, GxGroupedRecs<'a, 'b, R>)) -> Result<Gene, Error> {
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
                        .ok_or(Error::Gff("required 'transcript_id' attribute not found"))
                        .map(|id| trx_coords.entry(id).or_insert((INIT_COORD, vec![], None)));

                    match rec.feature_type() {
                        GENE_STR => {
                            gene_coord = adjust_coord(gene_coord, &rec);
                        },
                        TRANSCRIPT_STR => {
                            let trx_entry = rtrx_entry?;
                            trx_entry.0 = adjust_coord(trx_entry.0, &rec);
                        },
                        EXON_STR => {
                            let trx_entry = rtrx_entry?;
                            (trx_entry.1).push((*rec.start() - 1, *rec.end()));
                        },
                        CDS_STR => {
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

impl<'a, R> Iterator for GffGenes<'a, R> where R: io::Read {

    type Item = Result<Gene, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.into_iter().map(Self::group_to_gene).next()
    }
}

pub struct GffTranscripts<'a, R: 'a> where R: io::Read {
    inner: GroupBy<TrxGroupKey, TrxRecords<'a, R>, TrxGroupFunc>,
}

type TrxFilterFunc = fn(&Result<gff::Record, Error>) -> bool;

type TrxRecords<'a, R> = Filter<GffRecords<'a, R>, TrxFilterFunc>;

type TrxGroupKey = Option<(String, String, String, Strand)>;

type TrxGroupFunc = fn(&Result<gff::Record, Error>) -> TrxGroupKey;

type TrxGroupedRecs<'a, 'b, R> = Group<'b, TrxGroupKey, TrxRecords<'a, R>, TrxGroupFunc>;

impl<'a, R> GffTranscripts<'a, R> where R: io::Read {

    fn transcript_filter_func(result: &Result<gff::Record, Error>) -> bool {
        match result {
            &Ok(ref rec) =>
                rec.attributes().contains_key(GENE_ID_STR) &&
                rec.attributes().contains_key(TRANSCRIPT_ID_STR),
            &Err(_) => true,
        }
    }

    fn transcript_group_func(result: &Result<gff::Record, Error>) -> TrxGroupKey {
        result.as_ref().ok()
            .map(|res| {
                let gene_id = res.attributes().get(GENE_ID_STR)
                    .expect("'gene_id' attribute not found")
                    .clone();
                let transcript_id = res.attributes().get(TRANSCRIPT_ID_STR)
                    .expect("'transcript_id' attribute not found")
                    .clone();
                let seq_name = res.seqname().to_owned();
                let strand = res.strand().unwrap_or(Strand::Unknown);

                (seq_name, gene_id, transcript_id, strand)
            })
    }

    fn group_to_transcript<'b>(group: (TrxGroupKey, TrxGroupedRecs<'a, 'b, R>)
    ) -> Result<Transcript, Error>
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
                    rec.attributes_mut().remove(TRANSCRIPT_ID_STR);

                    match rec.feature_type() {
                        TRANSCRIPT_STR => {
                            trx_coord = adjust_coord(trx_coord, &rec);
                        },
                        EXON_STR => {
                            exn_coords.push((*rec.start() - 1, *rec.end()));
                        },
                        CDS_STR => {
                            (coding_coord).or(Some(INIT_COORD))
                                .map(|coord| adjust_coord(coord, &rec));
                        },
                        _ => {},
                    }
                }

                TBuilder::new(seq_name, trx_coord.0, trx_coord.1)
                    .id(tid)
                    .strand(strand)
                    .coords(exn_coords, coding_coord)
                    .coding_incl_stop(false)
                    .build()
            }
        }
    }
}

impl<'a, R> Iterator for GffTranscripts<'a, R> where R: io::Read {

    type Item = Result<Transcript, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.into_iter().map(Self::group_to_transcript).next()
    }
}

struct GffRow(String, String, String, u64, u64, String, String, String, HashMap<String, String>);


impl From<GffRow> for gff::Record {
    fn from(row: GffRow) -> Self {
        let mut rec = gff::Record::default();
        *rec.seqname_mut() = row.0;
        *rec.source_mut() = row.1;
        *rec.feature_type_mut() = row.2;
        *rec.start_mut() = row.3 + 1;
        *rec.end_mut() = row.4;
        *rec.score_mut() = row.5;
        *rec.strand_mut() = row.6;
        *rec.frame_mut() = row.7;
        *rec.attributes_mut() = row.8;
        rec
    }
}

impl Gene {

    #[inline(always)]
    fn num_records(&self) -> usize {
        1 + self.transcripts().values()
            .map(|ref trx| trx.num_records())
            .fold(0, |acc, x| acc + x)
    }

    // TODO: also handle gene-level features
    fn to_gff_records(mut self) -> Result<Vec<gff::Record>, Error> {

        let (source, score) = extract_source_score(&mut self.attributes);
        let strand = strand_to_string(self.strand());
        let gene_id = match self.id.take() {
            Some(gid) => {
                self.attributes.insert(GENE_ID_STR.to_owned(), gid.clone());
                gid
            },
            None =>  {
                self.attributes.get(GENE_ID_STR)
                    .ok_or(Error::Gff("required 'gene_id' attribute not found"))
                    .map(|gid| gid.clone())?
            },
        };

        let mut recs = Vec::with_capacity(self.num_records());

        let row = GffRow(
            self.seq_name().to_owned(), source, GENE_STR.to_owned(),
            self.start(), self.end(), score, strand.to_owned(),
            UNK_STR.to_owned(), self.attributes.clone());
        recs.push(gff::Record::from(row));

        for (_, mut transcript) in self.into_transcripts() {
            let mut trx_recs = transcript.to_gff_records(gene_id.as_str())?;
            recs.append(&mut trx_recs);
        }

        Ok(recs)
    }
}

impl Transcript {

    #[inline(always)]
    fn num_records(&self) -> usize {
        1 + self.exons().iter()
            .map(|ref exn| 1 + exn.features.len())
            .fold(0, |acc, x| acc + x)
    }

    // TODO: also handle transcript-level features
    fn to_gff_records(mut self, gene_id: &str) -> Result<Vec<gff::Record>, Error> {

        self.attributes.insert(GENE_ID_STR.to_owned(), gene_id.to_owned());
        let (source, score) = extract_source_score(&mut self.attributes);
        let strand = strand_to_string(self.strand());
        let transcript_id = match self.id.take() {
            Some(tid) => {
                self.attributes.insert(TRANSCRIPT_ID_STR.to_owned(), tid.clone());
                tid
            },
            None =>  {
                self.attributes.get(TRANSCRIPT_ID_STR)
                    .ok_or(Error::Gff("required 'transcript_id' attribute not found"))
                    .map(|tid| tid.clone())?
            },
        };

        let mut recs = Vec::with_capacity(self.num_records());

        let trx_row =
            GffRow(self.seq_name().to_owned(), source, TRANSCRIPT_STR.to_owned(),
                   self.start(), self.end(), score,
                   strand.to_owned(), UNK_STR.to_owned(), self.attributes.clone());
        recs.push(gff::Record::from(trx_row));

        for exon in self.into_exons() {
            let mut exon_recs = exon.to_gff_records(gene_id, transcript_id.as_str())?;
            recs.append(&mut exon_recs);
        }

        Ok(recs)
    }
}

impl EFK {

    #[inline(always)]
    fn get_feature_frame(&self) -> (String, String) {
        let (feature, frame) = match self {
            &EFK::UTR => (UTR_STR, UNK_STR),
            &EFK::UTR5 => (UTR5_STR, UNK_STR),
            &EFK::UTR3 => (UTR3_STR, UNK_STR),
            &EFK::CDS { frame: ref f } => (CDS_STR, frame_to_str(f)),
            &EFK::StopCodon { frame: ref f } => (STOP_CODON_STR, frame_to_str(f)),
            &EFK::StartCodon { frame: ref f } => (START_CODON_STR, frame_to_str(f)),
            &EFK::Any(ref s) => (s.as_str(), UNK_STR),
        };
        (feature.to_owned(), frame.to_owned())
    }
}

impl Exon {

    fn to_gff_records(mut self, gene_id: &str, transcript_id: &str
    ) -> Result<Vec<gff::Record>, Error> {

        self.attributes.insert(GENE_ID_STR.to_owned(), gene_id.to_owned());
        self.attributes.insert(TRANSCRIPT_ID_STR.to_owned(), transcript_id.to_owned());
        let (source, score) = extract_source_score(&mut self.attributes);
        let strand = strand_to_string(self.strand());

        let mut recs = Vec::with_capacity(1 + self.features.len());

        let exon_row =
            GffRow(self.seq_name().to_owned(), source.clone(), EXON_STR.to_owned(),
                   self.start(), self.end(), score, strand.to_owned(),
                   UNK_STR.to_owned(), self.attributes.clone());
        recs.push(gff::Record::from(exon_row));

        for fx in self.features.iter() {
            let (feature, frame) = fx.kind.get_feature_frame();
            let fx_row =
                GffRow(self.seq_name().to_owned(), source.clone(), feature,
                       self.start(), self.end(), UNK_STR.to_owned(),
                       strand.to_owned(), UNK_STR.to_owned(), self.attributes.clone());
            recs.push(gff::Record::from(fx_row));
        }

        Ok(recs)
    }
}

pub struct Writer<W: io::Write> {
    inner: gff::Writer<W>,
}

impl<W: io::Write> Writer<W> {

    pub fn from_writer(in_writer: W, gff_type: GffType) -> Writer<W> {
        Writer {
            inner: gff::Writer::new(in_writer, gff_type),
        }
    }

    pub fn write_gene(&mut self, mut gene: Gene) -> Result<(), Error> {
        for record in gene.to_gff_records()?.iter() {
            self.inner.write(&record)?;
        }
        Ok(())
    }
}

#[inline(always)]
fn adjust_coord(cur_coord: Coord<u64>, record: &gff::Record) -> Coord<u64> {
    (min(cur_coord.0, *record.start() - 1),
     max(cur_coord.1, *record.end()))
}

#[inline(always)]
fn extract_source_score(attributes: &mut HashMap<String, String>) -> (String, String) {
    let source = attributes.remove("source").unwrap_or(UNK_STR.to_owned());
    let score = attributes.remove("score").unwrap_or(UNK_STR.to_owned());
    (source, score)
}

#[inline(always)]
fn strand_to_string(strand: &Strand) -> String {
    let value = match strand {
        &Strand::Forward => "+",
        &Strand::Reverse => "-",
        &Strand::Unknown => UNK_STR,
    };
    value.to_owned()
}

#[inline(always)]
fn frame_to_str(frame: &Option<u8>) -> &str {
    match frame {
        &Some(0) => "0",
        &Some(1) => "1",
        &Some(2) => "2",
        _ => UNK_STR,
    }
}
