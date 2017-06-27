/*! Reader for GFF format variants.

The GFF format is a feature-oriented format that is commonly used to store gene annotation data.

An unofficial specification of the formats can be found [here](http://mblab.wustl.edu/GTF22.html)
or [here](http://www.ensembl.org/info/website/upload/gff.html).

The reader provided by this module is based on a modified version of the GFF reader provided by
the [rust-bio](https://github.com/rust-bio/rust-bio) library.
*/
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
use regex::{Error as RegexError, Regex};

use {Coord, Exon, ExonFeatureKind as EFK, Gene, Strand, TBuilder, Transcript,
     RawTrxCoords, INIT_START, INIT_END, INIT_COORD, DEF_ID};
use utils::{OptionDeref, update_seq_name};


/// Name for gene features.
const GENE_STR: &'static str = "gene";

/// Name for transcript features.
const TRANSCRIPT_STR: &'static str = "transcript";

/// Name for exon features.
const EXON_STR: &'static str = "exon";

/// Name for generic UTR features.
const UTR_STR: &'static str = "UTR";

/// Name for 5'UTR features.
const UTR5_STR: &'static str = "UTR5";

/// Name for 3'UTR features.
const UTR3_STR: &'static str = "UTR3";

/// Name for CDS features.
const CDS_STR: &'static str = "CDS";

/// Name for start codon features.
const START_CODON_STR: &'static str = "start_codon";

/// Name for stop codon features.
const STOP_CODON_STR: &'static str = "stop_codon";

/// Name for attribute key of gene identifiers.
const GENE_ID_STR: &'static str = "gene_id";

/// Name for attribute key of transcript identifiers.
const TRANSCRIPT_ID_STR: &'static str = "transcript_id";

/// Value for columns that are undefined, as a string.
const UNK_STR: &'static str = ".";

/// Value for columns that are undefined, as a char.
const UNK_CHAR: char = '.';

quick_error! {
    /// Errors that occur when reading GFF file variants.
    #[derive(Debug)]
    pub enum GffError {
        /// Occurs when a record does not have any gene identifier attribute.
        MissingGeneId {
            description("gene identifier attribute not found")
        }
        /// Occurs when a record does not have an expected transcript identifier attribute.
        MissingTranscriptId {
            description("transcript identifier attribute not found")
        }
        /// Occurs when a record contains multiple transcript identifier attributes.
        MultipleTranscriptIds {
            description("more than one 'transcript_id' found")
        }
        /// Occurs when a stop codon feature intersects a CDS feature.
        StopCodonInCds(tid: Option<String>) {
            description("'stop_codon' feature intersects cds")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        /// Occurs when an expected transcript feature record is not found.
        MissingTranscript(tid: Option<String>) {
            description("no 'transcript' feature present")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        /// Occurs when more than one expected transcript features are found.
        MultipleTranscripts {
            description("multiple 'transcript' features present")
        }
        /// Occurs when a stop codon exists in a transcript but no start codons are found.
        OrphanStop(tid: Option<String>) {
            description("stop codon exists without start codon")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        /// Occurs when a start codon exists in a transcript but no stop codons are found.
        OrphanStart(tid: Option<String>) {
            description("start codon exists without stop codon")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        /// Occurs when start and/or stop codons exists in a transcript without any CDS.
        OrphanCodon(tid: Option<String>) {
            description("start and stop codon exists without cds")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        /// Occurs when one or more CDSes exist in a transcript without any start and/or stop
        /// codons.
        OrphanCds(tid: Option<String>) {
            description("cds exists without start and/or stop codon")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        /// Occurs when an unsupported GFF variant is used.
        UnsupportedGffType {
            description("unsupported gff type")
        }
        /// Generic wrapper type for errors from the regex crate.
        Regex(err: RegexError) {
            description(err.description())
            from()
            cause(err)
        }
        /// Generic wrapper for GFF errors from the rust-bio crate.
        Bio(err: gff::GffError) {
            description(err.description())
            display("{}", err)
            from()
            cause(err)
        }
    }
}

/// GFF reader.
pub struct Reader<R: io::Read> {
    inner: gff::Reader<R>,
    gene_id_attr: String,
    transcript_id_attr: String,
    seq_name_prefix: Option<String>,
    seq_name_lstrip: Option<String>,
    loose_codons: bool,
    pub(crate) gff_type: GffType,
}

impl<R: io::Read> Reader<R> {

    /// Creates a GFF reader of the given variant from another reader.
    pub fn from_reader(in_reader: R, gff_type: GffType) -> Reader<R> {
        Reader {
            inner: gff::Reader::new(in_reader, gff_type),
            gene_id_attr: GENE_ID_STR.to_owned(),
            transcript_id_attr: TRANSCRIPT_ID_STR.to_owned(),
            seq_name_prefix: None,
            seq_name_lstrip: None,
            loose_codons: false,
            gff_type: gff_type.clone(),
        }
    }

    /// Sets the reader to use the given attribute key for getting gene identifiers.
    pub fn gene_id_attr<T>(&mut self, gene_id_attr: T) -> &mut Self
        where T: Into<String>
    {
        self.gene_id_attr = gene_id_attr.into();
        self
    }

    /// Sets the reader to use the given attribute key for getting transcript identifiers.
    pub fn transcript_id_attr<T>(&mut self, transcript_id_attr: T) -> &mut Self
        where T: Into<String>
    {
        self.transcript_id_attr = transcript_id_attr.into();
        self
    }

    /// Sets the reader to add the given prefix to all sequence names.
    pub fn seq_name_prefix<T>(&mut self, prefix: Option<T>) -> &mut Self
        where T: Into<String>
    {
        self.seq_name_prefix = prefix.map(|v| v.into());
        self
    }

    /// Sets the reader to trim the given string from all sequence names if present at the
    /// beginning.
    pub fn seq_name_lstrip<T>(&mut self, lstrip: Option<T>) -> &mut Self
        where T: Into<String>
    {
        self.seq_name_lstrip = lstrip.map(|v| v.into());
        self
    }

    /// Sets the reader to use CDS coordinates when start and/or stop codons for transcripts
    /// can not be found.
    pub fn loose_codons(&mut self, loose_codons: bool) -> &mut Self {
        self.loose_codons = loose_codons;
        self
    }

    /// Creates an iterator of transcripts.
    ///
    /// This iterator reads all GFF records into memory first, before sorting and grouping them
    /// into transcripts. This is because features of a transcript may be interspersed with
    /// features from another transcript.
    pub fn transcripts(&mut self) -> ::Result<GffTranscripts> {

        let gid_regex = make_gff_id_regex(self.gene_id_attr.as_str(), self.gff_type)?;
        let tid_regex = make_gff_id_regex(self.transcript_id_attr.as_str(), self.gff_type)?;
        let prefix = self.seq_name_prefix.clone();
        let lstrip = self.seq_name_lstrip.clone();

        let mut parts = Vec::new();
        for result in self.raw_rows_stream() {
            let mut row = result.map_err(::Error::from)?;
            update_seq_name(&mut row.0, prefix.as_deref(),
                            lstrip.as_deref().map(|v| (v, v.len())));
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
            loose_codons: self.loose_codons,
        })
    }

    /// Creates an iterator of GFF rows.
    pub(crate) fn raw_rows_stream(&mut self) -> GffRawRows<R> {
        GffRawRows {
            inner: self.inner.raw_rows()
        }
    }
}

impl Reader<fs::File> {

    /// Creates a refFlat reader that reads from the given path.
    pub fn from_file<P: AsRef<Path>>(path: P, gff_type: GffType) -> io::Result<Self> {
        fs::File::open(path).map(|file| Reader::from_reader(file, gff_type))
    }
}

/// Iterator over GFF rows.
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

/// Helper struct for creating transcripts.
///
/// This struct is meant to be used when complete parsing of the GFF attribute column is not
/// required. Instead, only the minimum required values (gene and transcript identifiers) are
/// parsed.
#[derive(Debug, PartialEq)]
struct TrxPart {
    feature: String,
    chrom: String,
    coord: Coord<u64>,
    strand: Strand,
    transcript_id: String,
    gene_id: String,
}

/// The type used for sorting GFF records.
///
/// The tuple elements represent gene identifier, transcript identifier, sequence, start
/// coordinate, end coordinate, end strand.
type TrxSortKey = (String, String, String, u64, u64, u8);

impl TrxPart {

    /// Creates a `TrxPart` from the given GFF row and the gene and transcript identifier regexes.
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

    /// Returns a tuple of sorting key.
    fn sort_key(&self) -> TrxSortKey {
        (self.gene_id.clone(), self.transcript_id.clone(),
         self.chrom.clone(), self.coord.0, self.coord.1, self.strand_ord())
    }

    /// Returns a tuple of grouping key.
    fn transcript_group_key(&self) -> TrxGroupKey {
        (self.gene_id.clone(), self.transcript_id.clone(), self.chrom.clone(), self.strand)
    }

    /// Returns the u8 value for distinguishing strands.
    fn strand_ord(&self) -> u8 {
        match &self.strand {
            &Strand::Unknown => 0,
            &Strand::Forward => 1,
            &Strand::Reverse => 2,
        }
    }
}

/// Helper container of raw coordinate values.
///
/// This struct is meant to be updated as features of a transcript are parsed.
#[derive(Debug, Default)]
struct TrxCoords {
    trx_coord: Option<Coord<u64>>,
    exon_coords: Vec<Coord<u64>>,
    cds_coord: Option<Coord<u64>>,
    codon_5: Option<u64>,
    codon_3: Option<u64>,
}

impl TrxCoords {

    /// Sets the transcript 5'-most and 3'-most coordinates.
    ///
    /// If this is set more than once, an error will be returned.
    fn set_trx_coord(&mut self, coord: Coord<u64>) -> Result<(), GffError> {
        if let None = self.trx_coord {
            self.trx_coord = Some(coord);
        } else {
            return Err(GffError::MultipleTranscripts);
        }
        Ok(())
    }

    /// Adds an exon coordinate of the transcript.
    fn add_exon_coord(&mut self, coord: Coord<u64>) {
        self.exon_coords.push(coord);
    }

    /// Adds a CDS coordinate.
    ///
    /// This will update the 5' and 3'-most CDS coordinates.
    fn include_cds_coord(&mut self, coord: Coord<u64>) {
        self.cds_coord = (self.cds_coord).or(Some(INIT_COORD))
            .map(|(a, b)| (min(a, coord.0), max(b, coord.1)));
    }

    /// Adds a 5'-most codon coordinate.
    fn include_codon_5(&mut self, coord_5: u64) {
        self.codon_5 = (self.codon_5).or(Some(INIT_START))
            .map(|c| min(c, coord_5));
    }

    /// Ads a 3'-most codon coordinate.
    fn include_codon_3(&mut self, coord_3: u64) {
        self.codon_3 = (self.codon_3).or(Some(INIT_END))
            .map(|c| max(c, coord_3));
    }

    /// Returns coordinates required to create a transcript.
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

/// Iterator over transcripts created from GFF records.
pub struct GffTranscripts {
    groups: GroupBy<TrxGroupKey, vec::IntoIter<TrxPart>, TrxGroupFunc>,
    loose_codons: bool,
}

/// The type used for grouping records into transcripts.
///
/// The tuple elements represent gene identifier, transcript identifier, sequence name, and strand.
type TrxGroupKey = (String, String, String, Strand);

/// The type of the function used for creating record-grouping keys for transcripts.
type TrxGroupFunc = fn(&TrxPart) -> TrxGroupKey;

/// The type of the grouped records for creating transcripts.
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

/// Helper function to create regex for parsing GFF identifiers.
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
        Regex::new(&pat)
            .map_err(|e| ::Error::from(GffError::from(e)))
    })

}

impl Gene {

    /// Returns the number of GFF records the gene has.
    #[inline(always)]
    fn num_records(&self) -> usize {
        1 + self.transcripts().values()
            .map(|ref trx| trx.num_records())
            .fold(0, |acc, x| acc + x)
    }

    // TODO: also handle gene-level features
    /// Transforms the gene into GFF records.
    pub fn into_gff_records(mut self) -> ::Result<Vec<gff::Record>> {

        let mut attribs = self.set_attributes(MultiMap::new());

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

    /// Returns the number of GFF records the transcript has.
    #[inline(always)]
    fn num_records(&self) -> usize {
        1 + self.exons().iter()
            .map(|ref exn| 1 + exn.features().len())
            .fold(0, |acc, x| acc + x)
    }

    // TODO: also handle transcript-level features
    /// Transforms the transcript into GFF records.
    pub fn into_gff_records(mut self) -> ::Result<Vec<gff::Record>> {

        let mut attribs = self.set_attributes(MultiMap::new());

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

    /// Returns the feature name and the frame of the exon feature kind.
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

    /// Transforms the exon into GFF records.
    pub fn into_gff_records(mut self) -> ::Result<Vec<gff::Record>> {

        let mut attribs = self.set_attributes(MultiMap::new());

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

/// Helper function to extract source and score attributes.
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

/// Helper function to create a char given a strand reference.
#[inline(always)]
fn strand_to_char(strand: &Strand) -> char {
    match strand {
        &Strand::Forward => '+',
        &Strand::Reverse => '-',
        &Strand::Unknown => UNK_CHAR,
    }
}

/// Helper function to create a char given an optional frame.
#[inline(always)]
fn frame_to_char(frame: &Option<u8>) -> char {
    match frame {
        &Some(0) => '0',
        &Some(1) => '1',
        &Some(2) => '2',
        _ => UNK_CHAR,
    }
}
