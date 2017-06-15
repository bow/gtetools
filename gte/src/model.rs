use std::cmp::{max, min};
use std::mem;
use std::error::Error;

use bio::utils::{self as bio_utils, Interval, IntervalError};
use bio::utils::Strand;
use linked_hash_map::LinkedHashMap;
use multimap::MultiMap;

use {Coord, RawTrxCoord, consts};
use consts::DEF_ID;
use utils::{OptionDeref, coord_to_interval};

use self::ExonFeatureKind::*;


macro_rules! impl_common {
    ($struct_ty:ty) => (

        impl $struct_ty {

            pub fn seq_name(&self) -> &str {
                self.seq_name.as_str()
            }

            pub fn set_seq_name<T>(&mut self, name: T)
                where T: Into<String>
            {
                self.seq_name = name.into()
            }

            pub fn id(&self) -> Option<&str> {
                self.id.as_deref()
            }

            pub fn strand(&self) -> &Strand {
                &self.strand
            }

            pub fn set_strand(&mut self, strand: Strand) {
                self.strand = strand
            }

            pub fn attributes(&self) -> &MultiMap<String, String> {
                &self.attributes
            }

            pub fn attributes_mut(&mut self) -> &mut MultiMap<String, String> {
                &mut self.attributes
            }

            pub fn set_attributes(&mut self, attributes: MultiMap<String, String>) {
                self.attributes = attributes
            }

            pub fn take_attributes(&mut self) -> MultiMap<String, String> {
                mem::replace(&mut self.attributes, MultiMap::new())
            }

            pub fn interval(&self) -> &Interval<u64> {
                &self.interval
            }

            pub fn start(&self) -> u64 {
                self.interval.start
            }

            pub fn end(&self) -> u64 {
                self.interval.end
            }

            #[inline]
            pub fn span(&self) -> u64 {
                self.end() - self.start()
            }
        }

    );
}

#[derive(Debug, Clone, PartialEq)]
pub struct Feature<K: FeatureKind> {
    interval: Interval<u64>,
    kind: K,
}

impl<K: FeatureKind> Feature<K> {

    pub fn new(interval: Interval<u64>, kind: K) -> Self {
        Feature {
            interval: interval,
            kind: kind,
        }
    }

    pub fn kind(&self) -> &K {
        &self.kind
    }

    pub fn interval(&self) -> &Interval<u64> {
        &self.interval
    }

    pub fn start(&self) -> u64 {
        self.interval.start
    }

    pub fn end(&self) -> u64 {
        self.interval.end
    }

    #[inline]
    pub fn span(&self) -> u64 {
        self.end() - self.start()
    }
}

pub trait FeatureKind {}

#[derive(Debug, Clone, PartialEq)]
pub enum ExonFeatureKind {
    UTR,
    UTR5,
    UTR3,
    CDS { frame: Option<u8> },
    StartCodon { frame: Option<u8> },
    StopCodon { frame: Option<u8> },
    Any(String),
}

impl FeatureKind for ExonFeatureKind {}

pub type ExonFeature = Feature<ExonFeatureKind>;

#[derive(Debug, Clone, PartialEq)]
pub enum TranscriptFeatureKind {
    Intron,
    Any(String),
}

impl FeatureKind for TranscriptFeatureKind {}

pub type TranscriptFeature = Feature<TranscriptFeatureKind>;

#[derive(Debug, Clone, PartialEq)]
pub struct GeneFeatureKind(String);

impl FeatureKind for GeneFeatureKind {}

pub type GeneFeature = Feature<GeneFeatureKind>;

#[derive(Debug, Clone)]
pub struct Exon {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    id: Option<String>,
    gene_id: Option<String>,
    transcript_id: Option<String>,
    attributes: MultiMap<String, String>,
    features: Vec<ExonFeature>,
}

impl_common!(Exon);

impl Exon {

    pub fn set_id<T>(&mut self, id: Option<T>)
        where T: Into<String>
    {
        self.id = id.map(|v| v.into())
    }

    pub fn transcript_id(&self) -> Option<&str> {
        self.transcript_id.as_deref()
    }

    pub fn set_transcript_id<T>(&mut self, transcript_id: Option<T>)
        where T: Into<String>
    {
        self.transcript_id = transcript_id.map(|v| v.into())
    }

    pub fn gene_id(&self) -> Option<&str> {
        self.gene_id.as_deref()
    }

    pub fn set_gene_id<T>(&mut self, gene_id: Option<T>)
        where T: Into<String>
    {
        self.gene_id = gene_id.map(|v| v.into())
    }

    pub fn features(&self) -> &[ExonFeature] {
        self.features.as_slice()
    }

    pub fn features_mut(&mut self) -> &mut [ExonFeature] {
        self.features.as_mut_slice()
    }

    pub fn set_features(&mut self, features: Vec<ExonFeature>) -> Result<(), FeatureError> {
        let (new_start, new_end) = features.iter()
            .fold(consts::INIT_COORD,
                  |acc, x| (min(acc.0, x.start()),
                            max(acc.1, x.end())));
        self.interval = coord_to_interval(new_start, new_end)?;
        self.features = features;
        Ok(())
    }

    pub fn take_features(self) -> Vec<ExonFeature> {
        self.features
    }
}

pub struct EBuilder {
    seq_name: String,
    start: u64,
    end: u64,
    strand: Option<Strand>,
    strand_char: Option<char>,
    id: Option<String>,
    transcript_id: Option<String>,
    gene_id: Option<String>,
    attributes: MultiMap<String, String>,
    features: Vec<ExonFeature>,
}

impl EBuilder {

    pub fn new<T>(seq_name: T, start: u64, end: u64) -> Self
        where T: Into<String>
    {
        EBuilder {
            seq_name: seq_name.into(),
            start: start,
            end: end,
            strand: None,
            strand_char: None,
            id: None,
            transcript_id: None,
            gene_id: None,
            attributes: MultiMap::new(),
            features: Vec::new(),
        }
    }

    pub fn strand(mut self, strand: Strand) -> Self {
        self.strand = Some(strand);
        self
    }

    pub fn strand_char(mut self, strand_char: char) -> Self {
        self.strand_char = Some(strand_char);
        self
    }

    pub fn id<T>(mut self, id: T) -> Self
        where T: Into<String>
    {
        self.id = Some(id.into());
        self
    }

    pub fn transcript_id<T>(mut self, transcript_id: T) -> Self
        where T: Into<String>
    {
        self.transcript_id = Some(transcript_id.into());
        self
    }

    pub fn gene_id<T>(mut self, gene_id: T) -> Self
        where T: Into<String>
    {
        self.gene_id = Some(gene_id.into());
        self
    }

    pub fn attribute<K, V>(mut self, key: K, value: V) -> Self
        where K: Into<String>, V: Into<String>
    {
        self.attributes.insert(key.into(), value.into());
        self
    }

    pub fn attributes(mut self, attributes: MultiMap<String, String>) -> Self {
        self.attributes = attributes;
        self
    }

    pub fn feature(mut self, feature: ExonFeature) -> Self {
        self.features.push(feature);
        self
    }

    pub fn features(mut self, features: Vec<ExonFeature>) -> Self {
        self.features = features;
        self
    }

    pub fn build(self) -> ::Result<Exon> {
        let interval = coord_to_interval(self.start, self.end)
            .map_err(::Error::Feature)?;
        let strand = resolve_strand_input(self.strand, self.strand_char)
            .map_err(::Error::Feature)?;
        let feature = Exon {
            seq_name: self.seq_name,
            interval: interval,
            strand: strand,
            id: self.id,
            transcript_id: self.transcript_id,
            gene_id: self.gene_id,
            attributes: self.attributes,
            features: self.features,
        };
        Ok(feature)
    }
}

#[derive(Debug, Clone)]
pub struct Transcript {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    id: Option<String>,
    gene_id: Option<String>,
    attributes: MultiMap<String, String>,
    exons: Vec<Exon>,
}

impl_common!(Transcript);

impl Transcript {

    pub fn set_id<T>(&mut self, id: Option<T>)
        where T: Into<String> + Clone
    {
        for exon in self.exons.iter_mut() {
            exon.set_transcript_id(id.clone())
        }
        self.id = id.map(|id| id.into())
    }

    pub fn gene_id(&self) -> Option<&str> {
        self.gene_id.as_deref()
    }

    pub fn set_gene_id<T>(&mut self, gene_id: Option<T>)
        where T: Into<String> + Clone
    {
        for exon in self.exons.iter_mut() {
            exon.set_gene_id(gene_id.clone())
        }
        self.gene_id = gene_id.map(|v| v.into())
    }

    pub fn exons(&self) -> &[Exon] {
        self.exons.as_slice()
    }

    pub fn take_exons(self) -> Vec<Exon> {
        self.exons
    }

    pub fn coding_coord(&self, incl_stop: bool) -> Option<Coord<u64>> {
        let start = self.coding_start_coord(incl_stop);
        let end = self.coding_end_coord(incl_stop);
        match (start, end) {
            (Some(s), Some(e)) => Some((s, e)),
            _ => None,
        }
    }

    fn coding_start_coord(&self, incl_stop: bool) -> Option<u64> {
        match &self.strand {
            &Strand::Forward => {
                for exon in self.exons.iter() {
                    for fx in exon.features.iter() {
                        if let StartCodon { .. } = fx.kind {
                            return Some(fx.interval.start);
                        }
                    }
                }
                None
            },
            &Strand::Reverse => {
                let mut codon_rem = if incl_stop { 0 } else { 3 };
                for exon in self.exons.iter() {
                    for fx in exon.features.iter() {
                        if let StopCodon { .. } = fx.kind {
                            if incl_stop {
                                return Some(fx.interval.start);
                            }
                            codon_rem -= fx.span();
                            if codon_rem == 0 {
                                return Some(fx.interval.end);
                            }
                        }
                    }
                }
                None
            },
            &Strand::Unknown if incl_stop => {
                for exon in self.exons.iter() {
                    for fx in exon.features.iter() {
                        if let CDS { .. } = fx.kind {
                            return Some(fx.interval.start)
                        }
                    }
                }
                None
            },
            _ => None
        }
    }

    fn coding_end_coord(&self, incl_stop: bool) -> Option<u64> {
        match &self.strand {
            &Strand::Forward => {
                let mut codon_rem = if incl_stop { 0 } else { 3 };
                for exon in self.exons.iter().rev() {
                    for fx in exon.features.iter().rev() {
                        if let StopCodon { .. } = fx.kind {
                            if incl_stop {
                                return Some(fx.interval.end);
                            }
                            codon_rem -= fx.span();
                            if codon_rem == 0 {
                                return Some(fx.interval.start);
                            }
                        }
                    }
                }
                None
            },
            &Strand::Reverse => {
                for exon in self.exons.iter().rev() {
                    for fx in exon.features.iter().rev() {
                        if let StartCodon { .. } = fx.kind {
                            return Some(fx.interval.end);
                        }
                    }
                }
                None
            },
            &Strand::Unknown if incl_stop => {
                for exon in self.exons.iter().rev() {
                    for fx in exon.features.iter().rev() {
                        if let CDS { .. } = fx.kind {
                            return Some(fx.interval.end);
                        }
                    }
                }
                None
            },
            _ => None,
        }
    }

}

pub struct TBuilder {
    seq_name: String,
    start: u64,
    end: u64,
    strand: Option<Strand>,
    strand_char: Option<char>,
    id: Option<String>,
    gene_id: Option<String>,
    attributes: MultiMap<String, String>,
    // Input can be a vector of pre-made features ...
    exons: Option<Vec<Exon>>,
    // Or exon coordinates, possibly coupled with cds coord
    // NOTE: Can we instead of using Vec<_> here keep it as an unconsumed iterator?
    exon_coords: Option<Vec<Coord<u64>>>,
    coding_coord: Option<Coord<u64>>,
    coding_incl_stop: bool,
}

impl TBuilder {

    pub fn new<T>(seq_name: T, start: u64, end: u64) -> Self
        where T: Into<String>
    {
        TBuilder {
            seq_name: seq_name.into(),
            start: start,
            end: end,
            strand: None,
            strand_char: None,
            id: None,
            gene_id: None,
            attributes: MultiMap::new(),
            exons: None,
            exon_coords: None,
            coding_coord: None,
            coding_incl_stop: false,
        }
    }

    pub fn strand(mut self, strand: Strand) -> Self {
        self.strand = Some(strand);
        self
    }

    pub fn strand_char(mut self, strand_char: char) -> Self {
        self.strand_char = Some(strand_char);
        self
    }

    pub fn id<T>(mut self, id: T) -> Self
        where T: Into<String>
    {
        self.id = Some(id.into());
        self
    }

    pub fn gene_id<T>(mut self, gene_id: T) -> Self
        where T: Into<String>
    {
        self.gene_id = Some(gene_id.into());
        self
    }

    pub fn attribute<K, V>(mut self, key: K, value: V) -> Self
        where K: Into<String>, V: Into<String>
    {
        self.attributes.insert(key.into(), value.into());
        self
    }

    pub fn attributes(mut self, attributes: MultiMap<String, String>) -> Self {
        self.attributes = attributes;
        self
    }

    pub fn exons(mut self, exons: Vec<Exon>) -> Self {
        self.exons =
            if exons.is_empty() {
                None
            } else {
                Some(exons)
            };
        self
    }

    pub fn coords<E>(mut self, exon_coords: E, coding_coord: Option<Coord<u64>>) -> Self
        where E: IntoIterator<Item=Coord<u64>>
    {
        self.exon_coords = Some(exon_coords.into_iter().collect());
        self.coding_coord = coding_coord;
        self
    }

    pub fn coding_incl_stop(mut self, incl_stop: bool) -> Self {
        self.coding_incl_stop = incl_stop;
        self
    }

    pub fn build(self) -> ::Result<Transcript> {
        let interval = coord_to_interval(self.start, self.end)
            .map_err(::Error::Feature)?;
        let strand = resolve_strand_input(self.strand, self.strand_char)
            .map_err(::Error::Feature)?;
        let exons = resolve_exons_input(
            &self.seq_name, &interval, &strand, self.id.as_deref(),
            self.gene_id.as_deref(), None, // TODO: allow for exon IDs here
            self.exons, self.exon_coords.as_ref(), self.coding_coord,
            self.coding_incl_stop).map_err(::Error::Feature)?;

        let transcript = Transcript {
            seq_name: self.seq_name,
            interval: interval,
            strand: strand,
            id: self.id,
            gene_id: self.gene_id,
            attributes: self.attributes,
            exons: exons,
        };
        Ok(transcript)
    }
}

#[derive(Debug, Clone)]
pub struct Gene {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    id: Option<String>,
    attributes: MultiMap<String, String>,
    transcripts: LinkedHashMap<String, Transcript>,
}

impl_common!(Gene);

impl Gene {

    pub fn set_id<T>(&mut self, id: Option<T>)
        where T: Into<String> + Clone
    {
        for (_, transcript) in self.transcripts.iter_mut() {
            transcript.set_id(id.clone())
        }
        self.id = id.map(|v| v.into())
    }

    pub fn transcripts(&self) -> &LinkedHashMap<String, Transcript> {
        &self.transcripts
    }

    pub fn take_transcripts(self) -> LinkedHashMap<String, Transcript> {
        self.transcripts
    }
}

pub struct GBuilder {
    seq_name: String,
    start: u64,
    end: u64,
    strand: Option<Strand>,
    strand_char: Option<char>,
    id: Option<String>,
    attributes: MultiMap<String, String>,
    transcripts: Option<LinkedHashMap<String, Transcript>>,
    transcript_coords: Option<LinkedHashMap<String, RawTrxCoord>>,
    transcript_coding_incl_stop: bool,
}

impl GBuilder {

    pub fn new<T>(seq_name: T, start: u64, end: u64) -> Self
        where T: Into<String>
    {
        GBuilder {
            seq_name: seq_name.into(),
            start: start,
            end: end,
            strand: None,
            strand_char: None,
            id: None,
            attributes: MultiMap::new(),
            transcripts: None,
            transcript_coords: None,
            transcript_coding_incl_stop: false,
        }
    }

    pub fn strand(mut self, strand: Strand) -> Self {
        self.strand = Some(strand);
        self
    }

    pub fn strand_char(mut self, strand_char: char) -> Self {
        self.strand_char = Some(strand_char);
        self
    }

    pub fn id<T>(mut self, id: T) -> Self
        where T: Into<String>
    {
        self.id = Some(id.into());
        self
    }

    pub fn attribute<K, V>(mut self, key: K, value: V) -> Self
        where K: Into<String>, V: Into<String>
    {
        self.attributes.insert(key.into(), value.into());
        self
    }

    pub fn attributes(mut self, attributes: MultiMap<String, String>) -> Self {
        self.attributes = attributes;
        self
    }

    pub fn transcripts(mut self, transcripts: LinkedHashMap<String, Transcript>) -> Self {
        self.transcripts = Some(transcripts);
        self
    }

    pub fn transcript_coords(mut self, coords: LinkedHashMap<String, RawTrxCoord>)-> Self {
        self.transcript_coords = Some(coords);
        self
    }

    pub fn transcript_coding_incl_stop(mut self, incl_stop: bool) -> Self {
        self.transcript_coding_incl_stop = incl_stop;
        self
    }

    pub fn build(self) -> ::Result<Gene> {
        let interval = coord_to_interval(self.start, self.end)
            .map_err(::Error::Feature)?;
        let strand = resolve_strand_input(self.strand, self.strand_char)
            .map_err(::Error::Feature)?;
        let transcripts = resolve_transcripts_input(
            &self.seq_name, &interval, &strand, self.id.as_deref(),
            self.transcripts, self.transcript_coords, self.transcript_coding_incl_stop)?;

        let gene = Gene {
            seq_name: self.seq_name,
            interval: interval,
            strand: strand,
            id: self.id,
            attributes: self.attributes,
            transcripts: transcripts,
        };
        Ok(gene)
    }
}

quick_error! {
    #[derive(Debug)]
    pub enum FeatureError {
        InvalidInterval(err: bio_utils::IntervalError) {
            description(
                match err {
                    &IntervalError::InvalidRange =>
                        "interval start coordinate larger than its end coordinate",
                    ref otherwise => otherwise.description(),
                })
            from()
        }
        InvalidStrandChar(err: bio_utils::StrandError) {
            description(err.description())
            from()
        }
        ConflictingStrand {
            description("conflicting strand inputs specified")
        }
        UnspecifiedStrand {
            description("strand not specified")
        }
        InvalidExonInterval(tid: Option<String>) {
            description("exon has larger start than end coordinate")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        InvalidCodingInterval(tid: Option<String>) {
            description("coding region has larger start than end coordinate")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        UnspecifiedExons(tid: Option<String>) {
            description("transcript is defined without exons")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        UnmatchedExons(tid: Option<String>) {
            description("first and/or last exon coordinates do not match transcript \
                         start and/or end coordinates")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        CodingTooLarge(tid: Option<String>) {
            description("coding region leaves no room for stop codon in transcript")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        CodingTooSmall(tid: Option<String>) {
            description("coding region leaves no room for start codon")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        CodingNotFullyEnveloped(tid: Option<String>) {
            description("coding region not fully enveloped by exons")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        CodingInIntron(tid: Option<String>) {
            description("coding start and/or end lies in introns")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
        TranscriptNotFullyEnveloped(tid: Option<String>) {
            description("transcript coordinate not fully enveloped in gene coordinate")
            display(self_) -> ("{}, transcript ID: {}",
                               self_.description(), tid.as_deref().unwrap_or(DEF_ID))
        }
    }
}

fn resolve_strand_input(
    strand: Option<Strand>,
    strand_char: Option<char>)
-> Result<Strand, FeatureError>
{
    match (strand, strand_char) {
        (None, None) => Err(FeatureError::UnspecifiedStrand),
        (Some(sv), None) => Ok(sv),
        (None, Some(ref scv)) => Strand::from_char(scv).map_err(FeatureError::from),
        (Some(sv), Some(ref scv)) => {
            let sv_from_char = Strand::from_char(scv).map_err(FeatureError::from)?;
            if sv == sv_from_char {
                Ok(sv)
            } else {
                Err(FeatureError::ConflictingStrand)
            }
        }
    }
}

fn resolve_exons_input(
    transcript_seqname: &String,
    transcript_interval: &Interval<u64>,
    transcript_strand: &Strand,
    transcript_id: Option<&str>,
    gene_id: Option<&str>,
    exon_id: Option<&str>,
    exons: Option<Vec<Exon>>,
    exon_coords: Option<&Vec<Coord<u64>>>,
    coding_coord: Option<Coord<u64>>,
    coding_incl_stop: bool
) -> Result<Vec<Exon>, FeatureError>
{
    match (exons, exon_coords, coding_coord) {
        // nothing defined -> the transcript doesn't have any known exons
        (None, None, None) => Ok(Vec::new()),

        // only CDS defined -> must be an error
        (None, None, Some(_)) => Err(
            FeatureError::UnspecifiedExons(transcript_id.map(|tid| tid.to_owned()))),

        // features defined ~ takes precedence over coords (GTF input, since we need
        // to construct the tx features first to store its annotations)
        // TODO: Maybe do some checks to ensure the given features are correct?
        (Some(exns), _, _) => Ok(exns.into_iter().collect()),

        // exon defined & coords possibly defined (refFlat input)
        (None, Some(raw_exon_coords), raw_coding_coord) =>
            infer_exons(transcript_seqname, transcript_interval, transcript_strand, transcript_id,
                        gene_id, exon_id, raw_exon_coords, raw_coding_coord, coding_incl_stop),
    }
}

fn resolve_transcripts_input(
    gene_seqname: &String,
    gene_interval: &Interval<u64>,
    gene_strand: &Strand,
    gene_id: Option<&str>,
    transcripts: Option<LinkedHashMap<String, Transcript>>,
    transcript_coords: Option<LinkedHashMap<String, RawTrxCoord>>,
    transcript_coding_incl_stop: bool
) -> ::Result<LinkedHashMap<String, Transcript>>
{
    match (transcripts, transcript_coords) {
        // nothing defined -> the gene doesn't have any known transcripts
        (None, None) => Ok(LinkedHashMap::new()),

        // transcript defined, return it
        (Some(trxs), _) => Ok(trxs),

        // transcripts coords defined, create transcripts
        (None, Some(trxs_coords)) => {
            let mut trxs = LinkedHashMap::new();
            for (trx_id, (trx_coord, exon_coords, coding_coord)) in trxs_coords.into_iter() {

                if trx_coord.0 < gene_interval.start || trx_coord.1 > gene_interval.end {
                    let tid = Some(trx_id);
                    return Err(::Error::Feature(FeatureError::TranscriptNotFullyEnveloped(tid)));
                }

                let btrx = TBuilder::new(gene_seqname.clone(), trx_coord.0, trx_coord.1)
                    .strand(*gene_strand)
                    .id(trx_id.clone())
                    .coords(exon_coords, coding_coord)
                    .coding_incl_stop(transcript_coding_incl_stop);
                let trx = match gene_id {
                    Some(ref gid) => btrx
                        .gene_id(gid.to_owned())
                        .build()?,
                    None => btrx.build()?
                };
                trxs.insert(trx_id, trx);
            }
            Ok(trxs)
        },
    }
}

fn infer_exons(
    transcript_seqname: &String,
    transcript_interval: &Interval<u64>,
    transcript_strand: &Strand,
    transcript_id: Option<&str>,
    gene_id: Option<&str>,
    exon_id: Option<&str>,
    exon_coords: &Vec<Coord<u64>>,
    coding_coord: Option<Coord<u64>>,
    coding_incl_stop: bool
) -> Result<Vec<Exon>, FeatureError>
{

    let tid = transcript_id.map(|id| id.to_owned());

    if exon_coords.len() == 0 {
        return Err(FeatureError::UnspecifiedExons(tid));
    }

    let mut m_exon_coords = Vec::with_capacity(exon_coords.len());
    for &(a, b) in exon_coords.iter() {
        if a >= b {
            return Err(FeatureError::InvalidExonInterval(tid))
        }
        m_exon_coords.push((a, b));
    }
    m_exon_coords.sort();

    let adj_coding_coord =
        if coding_incl_stop {
            coding_coord.and_then(|(a, b)| {
                adjust_coding_coord(a, b, &transcript_strand, &m_exon_coords)
            })
        } else {
            coding_coord
        };

    let exon_r = (m_exon_coords.first().unwrap().0, m_exon_coords.last().unwrap().1);

    if exon_r.0 != transcript_interval.start || exon_r.1 != transcript_interval.end {
        return Err(FeatureError::UnmatchedExons(tid));
    }

    match adj_coding_coord {

        Some(coding_r) => {
            // Improper coding region is an error
            if coding_r.0 >= coding_r.1 {
                return Err(FeatureError::InvalidCodingInterval(tid));
            }
            // Coding coord must be fully enveloped by exon max-min
            if coding_r.0 < exon_r.0 || coding_r.1 > exon_r.1 {
                return Err(FeatureError::CodingNotFullyEnveloped(tid));
            }
            // Coding start and end must be in exons
            let cine = m_exon_coords.iter()
                .fold((false, false), |acc, c| {
                    (acc.0 || (c.0 <= coding_r.0 && coding_r.0 <= c.1),
                     acc.1 || (c.0 <= coding_r.1 && coding_r.1 <= c.1))
                });
            if !cine.0 || !cine.1 {
                return Err(FeatureError::CodingInIntron(tid));
            }
            // There must be room for stop codons (which is not inclusive in coding_coord)
            let stop_codon_ok = match transcript_strand {
                &Strand::Forward => coding_r.1 + 3 <= exon_r.1,
                &Strand::Reverse => coding_r.0 - 3 >= exon_r.0,
                &Strand::Unknown =>
                    coding_r.0 - 3 >= exon_r.0 && coding_r.1 + 3 <= exon_r.1,
            };
            if !stop_codon_ok {
                return Err(FeatureError::CodingTooLarge(tid));
            }
            infer_exon_features(&m_exon_coords, coding_r, &transcript_seqname, transcript_strand,
                                transcript_id, gene_id, exon_id)
        }

        // No CDS intervals mean we just sort the coordinates and create the exons
        None => {
            let mut features = Vec::with_capacity(m_exon_coords.len());
            for &(start, end) in m_exon_coords.iter() {
                features.push(
                    Exon {
                        seq_name: transcript_seqname.clone(),
                        interval: Interval::new(start..end).unwrap(),
                        strand: *transcript_strand,
                        id: exon_id.map(|id| id.to_owned()),
                        transcript_id: tid.clone(),
                        gene_id: gene_id.map(|id| id.to_owned()),
                        attributes: MultiMap::new(),
                        features: Vec::new(),
                    });
            }
            Ok(features)
        }
    }
}

fn adjust_coding_coord(mut start: u64, mut end: u64,
                       strand: &Strand, exon_coords: &Vec<Coord<u64>>
) -> Option<Coord<u64>>
{
    let mut codon_rem = 3;
    match strand {
        &Strand::Forward => {
            for &(exon_start, exon_end) in exon_coords.iter().rev() {
                if exon_start <= end && end <= exon_end {
                    let adj_end = max(end - codon_rem, exon_start);
                    codon_rem -= end - adj_end;
                    end = adj_end;
                    if codon_rem == 0 {
                        break;
                    }
                }
            }
        },
        &Strand::Reverse => {
            for &(exon_start, exon_end) in exon_coords.iter() {
                if exon_start <= start  && start <= exon_end {
                    let adj_start = min(start + codon_rem, exon_end);
                    codon_rem -= adj_start - start;
                    start = adj_start;
                    if codon_rem == 0 {
                        break;
                    }
                }
            }
        },
        &Strand::Unknown => {},
    }
    Some((start, end))
}

// requirements:
//  - exon coords sorted and nonempty
//  - coding coord within exon span
//  - coding_coord.0 < coding_coord.1
fn infer_exon_features(
    exon_coords: &Vec<Coord<u64>>,
    coding_r: Coord<u64>,
    transcript_seqname: &String,
    transcript_strand: &Strand,
    transcript_id: Option<&str>,
    gene_id: Option<&str>,
    exon_id: Option<&str>,
) -> Result<Vec<Exon>, FeatureError> {

    let mut exons: Vec<Exon> = Vec::with_capacity(exon_coords.len() * 2 + 4);
    let (utr1, utr2) = match transcript_strand {
        &Strand::Forward => (UTR5, UTR3),
        &Strand::Reverse => (UTR3, UTR5),
        &Strand::Unknown => (UTR, UTR),
    };
    let exn = |start, end, features| {
        Exon {
            seq_name: transcript_seqname.clone(),
            interval: Interval::new(start..end).unwrap(),
            strand: *transcript_strand,
            id: exon_id.map(|v| v.to_owned()),
            transcript_id: transcript_id.map(|v| v.to_owned()),
            gene_id: gene_id.map(|v| v.to_owned()),
            attributes: MultiMap::new(),
            features: features,
        }
    };
    let feat = |start, end, kind| {
        ExonFeature {
            interval: Interval::new(start..end).unwrap(),
            kind: kind,
        }
    };

    let tid = transcript_id.map(|id| id.to_owned());

    // how much we have consumed the 5' or 3' codon
    let (mut codon1_rem, mut codon2_rem) = (3, 3);
    for &(start, end) in exon_coords.iter() {

        if start < coding_r.0 {

            let mut exon = exn(start, end, vec![]);
            let utr_end =
                if let &Strand::Reverse = transcript_strand {
                    min(end, coding_r.0 - codon1_rem)
                } else {
                    min(end, coding_r.0)
                };
            if start < utr_end {
                exon.features.push(feat(start, utr_end, utr1.clone()));
            }

            if end < coding_r.0 {
                exons.push(exon);

            } else if end == coding_r.0 {
                if let &Strand::Reverse = transcript_strand {
                    let fx = feat(max(start, coding_r.0 - codon1_rem), coding_r.0,
                                  StopCodon { frame: None });
                    codon1_rem -= fx.span();
                    exon.features.push(fx);
                    codon1_rem = backtrack_and_push(&mut exons, StopCodon { frame: None },
                                                    codon1_rem, &feat);
                }
                exons.push(exon);

            } else if end > coding_r.0 && end < coding_r.1 {
                match transcript_strand {
                    &Strand::Forward => {
                        let fx = feat(coding_r.0, min(end, coding_r.0 + 3),
                                      StartCodon { frame: None });
                        codon1_rem -= fx.span();
                        exon.features.push(fx);
                    },
                    &Strand::Reverse => {
                        let fx = feat(max(start, coding_r.0 - codon1_rem), coding_r.0,
                                      StopCodon { frame: None });
                        codon1_rem -= fx.span();
                        exon.features.push(fx);
                        codon1_rem = backtrack_and_push(&mut exons, StopCodon { frame: None },
                                                        codon1_rem, &feat);
                    },
                    &Strand::Unknown => {},
                }
                exon.features.push(feat(coding_r.0, end, CDS { frame: None }));
                exons.push(exon);

            } else if end == coding_r.1 {
                if coding_r.1 - coding_r.0 < 3 {
                    // a coding region must have at least 3 bases for the start codon
                    return Err(FeatureError::CodingTooSmall(tid));
                }
                match transcript_strand {
                    &Strand::Forward => {
                        let fx = feat(coding_r.0, coding_r.0 + 3, StartCodon { frame: None });
                        codon1_rem -= fx.span();
                        exon.features.push(fx);
                        exon.features.push(feat(coding_r.0, coding_r.1, CDS { frame: None }));
                    },
                    &Strand::Reverse => {
                        let fx = feat(max(coding_r.0 - codon1_rem, start), coding_r.0,
                                      StopCodon { frame: None });
                        codon1_rem -= fx.span();
                        exon.features.push(fx);
                        codon1_rem = backtrack_and_push(&mut exons, StopCodon { frame: None },
                                                        codon1_rem, &feat);
                        exon.features.push(feat(coding_r.0, coding_r.1, CDS { frame: None }));
                        let fx = feat(max(start, coding_r.1 - codon2_rem),
                                      coding_r.1, StartCodon { frame: None });
                        codon2_rem -= fx.span();
                        exon.features.push(fx);
                        codon2_rem = backtrack_and_push(&mut exons, StartCodon { frame: None },
                                                        codon2_rem, &feat);
                    },
                    &Strand::Unknown => {
                        exon.features.push(feat(coding_r.0, coding_r.1, CDS { frame: None }));
                    },
                }
                exons.push(exon);

            } else if end > coding_r.1 {
                if coding_r.1 - coding_r.0 < 3 {
                    // coding region must have at least 3 bases
                    return Err(FeatureError::CodingTooSmall(tid));
                }
                match transcript_strand {
                    &Strand::Forward => {
                        exon.features.push(feat(coding_r.0, coding_r.0 + 3,
                                                StartCodon { frame: None }));
                        codon1_rem -= 3;
                        exon.features.push(feat(coding_r.0, coding_r.1, CDS { frame: None }));
                        let fx = feat(coding_r.1, min(end, coding_r.1 + codon2_rem),
                                      StopCodon { frame: None });
                        codon2_rem -= fx.span();
                        let codon2_end = fx.end();
                        exon.features.push(fx);
                        if codon2_end < end {
                            exon.features.push(feat(codon2_end, end, utr2.clone()));
                        }
                    },
                    &Strand::Reverse => {
                        let fx = feat(max(start, coding_r.0 - codon2_rem), coding_r.0,
                                      StopCodon { frame: None });
                        codon1_rem -= fx.span();
                        exon.features.push(fx);
                        codon1_rem = backtrack_and_push(&mut exons, StopCodon { frame: None },
                                                        codon1_rem, &feat);
                        exon.features.push(feat(coding_r.0, coding_r.1, CDS { frame: None }));
                        let fx = feat(coding_r.1 - codon2_rem, coding_r.1,
                                      StartCodon { frame: None });
                        codon2_rem -= fx.span();
                        exon.features.push(fx);
                        exon.features.push(feat(coding_r.1, end, utr2.clone()));
                    },
                    &Strand::Unknown => {
                        exon.features.push(feat(coding_r.0, coding_r.1, CDS { frame: None }));
                        exon.features.push(feat(coding_r.1, end, utr2.clone()));
                    },
                }
                exons.push(exon);

            } else {
                assert!(false, "unexpected: exon=[{},{}) cds=[{},{})",
                        start, end, coding_r.0, coding_r.1)
            }

        } else if start == coding_r.0 {

            if end < coding_r.1 {
                let mut exon = exn(start, end, vec![]);
                match transcript_strand {
                    &Strand::Forward => {
                        let fx = feat(start, min(start + codon1_rem, end),
                                      StartCodon { frame: None });
                        codon1_rem -= fx.span();
                        exon.features.push(fx);
                    },
                    &Strand::Reverse => {
                        codon1_rem = backtrack_and_push(&mut exons, StopCodon { frame: None },
                                                        codon1_rem, &feat);
                    },
                    &Strand::Unknown => {},
                }
                exon.features.push(feat(start, end, CDS { frame: None }));
                exons.push(exon);

            } else if end == coding_r.1 {
                if coding_r.1 - coding_r.0 < 3 {
                    // coding region must have at least 3 bases
                    return Err(FeatureError::CodingTooSmall(tid));
                }
                let mut exon = exn(start, end, vec![]);
                match transcript_strand {
                    &Strand::Forward => {
                        exon.features.push(feat(start, start + 3, StartCodon { frame: None }));
                        codon1_rem = 0;
                        exon.features.push(feat(start, end, CDS { frame: None }));
                    },
                    &Strand::Reverse => {
                        codon1_rem = backtrack_and_push(&mut exons, StopCodon { frame: None },
                                                        codon1_rem, &feat);
                        exon.features.push(feat(start, end, CDS { frame: None }));
                        exon.features.push(feat(end - 3, end, StartCodon { frame: None }));
                        codon2_rem = 0;
                    },
                    &Strand::Unknown => {
                        exon.features.push(feat(start, end, CDS { frame: None }));
                    },
                }
                exons.push(exon);

            } else if end > coding_r.1 {
                if coding_r.1 - coding_r.0 < 3 {
                    // coding region must have at least 3 bases
                    return Err(FeatureError::CodingTooSmall(tid));
                }
                let mut exon = exn(start, end, vec![]);
                match transcript_strand {
                    &Strand::Forward => {
                        exon.features.push(feat(start, start + 3, StartCodon { frame: None }));
                        codon1_rem -= 3;
                        exon.features.push(feat(start, coding_r.1, CDS { frame: None }));
                        let fx = feat(coding_r.1, min(end, coding_r.1 + codon2_rem),
                                      StopCodon { frame: None });
                        codon2_rem -= fx.span();
                        let codon2_end = fx.end();
                        exon.features.push(fx);
                        if codon2_end < end {
                            exon.features.push(feat(codon2_end, end, utr2.clone()));
                        }
                    },
                    &Strand::Reverse => {
                        codon1_rem = backtrack_and_push(&mut exons, StopCodon { frame: None },
                                                        codon1_rem, &feat);
                        exon.features.push(feat(start, coding_r.1, CDS { frame: None }));
                        exon.features.push(feat(coding_r.1 - 3, coding_r.1,
                                                StartCodon { frame: None }));
                        codon2_rem -= 3;
                        exon.features.push(feat(coding_r.1, end, utr2.clone()));
                    },
                    &Strand::Unknown => {
                        exon.features.push(feat(start, coding_r.1, CDS { frame: None }));
                        exon.features.push(feat(coding_r.1, end, utr2.clone()));
                    },
                }
                exons.push(exon);

            } else {
                assert!(false, "unexpected: exon=[{},{}) cds=[{},{})",
                        start, end, coding_r.0, coding_r.1)
            }

        } else if start > coding_r.0 && start < coding_r.1 {

            if end < coding_r.1 {
                let mut exon = exn(start, end, vec![]);
                if let &Strand::Forward = transcript_strand {
                    if codon1_rem > 0 {
                        let fx = feat(start, min(end, start + codon1_rem),
                                      StartCodon { frame: None });
                        codon1_rem -= fx.span();
                        exon.features.push(fx);
                    }
                }
                exon.features.push(feat(start, end, CDS { frame: None }));
                exons.push(exon);

            } else if end == coding_r.1 {
                let mut exon = exn(start, end, vec![]);
                match transcript_strand {
                    &Strand::Forward => {
                        if codon1_rem > 0 {
                            let fx = feat(start, min(end, start + codon1_rem),
                                          StartCodon { frame: None });
                            codon1_rem -= fx.span();
                            exon.features.push(fx);
                        }
                        exon.features.push(feat(start, end, CDS { frame: None }));
                    },
                    &Strand::Reverse => {
                        exon.features.push(feat(start, end, CDS { frame: None }));
                        let fx = feat(max(start, coding_r.1 - codon2_rem), coding_r.1,
                                      StartCodon { frame: None });
                        codon2_rem -= fx.span();
                        exon.features.push(fx);
                        codon2_rem = backtrack_and_push(
                            &mut exons, StartCodon { frame: None }, codon2_rem, &feat);
                    },
                    &Strand::Unknown => {
                        exon.features.push(feat(start, end, CDS { frame: None }));
                    },
                }
                exons.push(exon);

            } else if end > coding_r.1 {
                let mut exon = exn(start, end, vec![]);
                match transcript_strand {
                    &Strand::Forward => {
                        if codon1_rem > 0 {
                            let fx = feat(start, min(coding_r.1, start + codon1_rem),
                                          StartCodon { frame: None });
                            codon1_rem -= fx.span();
                            exon.features.push(fx);
                        }
                        exon.features.push(feat(start, coding_r.1, CDS { frame: None }));
                        let fx = feat(coding_r.1, min(coding_r.1 + codon2_rem, end),
                                      StopCodon { frame: None });
                        codon2_rem -= fx.span();
                        let codon2_end = fx.end();
                        exon.features.push(fx);
                        if codon2_end < end {
                            exon.features.push(feat(codon2_end, end, utr2.clone()));
                        }
                    },
                    &Strand::Reverse => {
                        exon.features.push(feat(start, coding_r.1, CDS { frame: None }));
                        let fx = feat(max(start, coding_r.1 - codon2_rem), coding_r.1,
                                      StartCodon { frame: None });
                        codon2_rem -= fx.span();
                        exon.features.push(fx);
                        codon2_rem = backtrack_and_push(
                            &mut exons, StartCodon { frame: None }, codon2_rem, &feat);
                        exon.features.push(feat(coding_r.1, end, utr2.clone()));
                    },
                    &Strand::Unknown => {
                        exon.features.push(feat(start, coding_r.1, CDS { frame: None }));
                        exon.features.push(feat(coding_r.1, end, utr2.clone()));
                    },
                }
                exons.push(exon);

            } else {
                assert!(false, "unexpected: exon=[{},{}) cds=[{},{})",
                        start, end, coding_r.0, coding_r.1)
            }

        } else if start >= coding_r.1 {
            let mut exon = exn(start, end, vec![]);
            match transcript_strand {
                &Strand::Forward => {
                    if codon2_rem > 0 {
                        let fx = feat(start, min(start + codon2_rem, end),
                                      StopCodon { frame: None });
                        codon2_rem -= fx.span();
                        let codon2_end = fx.end();
                        exon.features.push(fx);
                        if codon2_end < end {
                            exon.features.push(feat(codon2_end, end, utr2.clone()));
                        }
                    } else {
                        exon.features.push(feat(start, end, utr2.clone()));
                    }
                },
                &Strand::Reverse => {
                    codon2_rem = backtrack_and_push(
                        &mut exons, StartCodon { frame: None }, codon2_rem, &feat);
                    exon.features.push(feat(start, end, utr2.clone()));
                },
                &Strand::Unknown => {
                    exon.features.push(feat(start, end, utr2.clone()));
                },
            }
            exons.push(exon);

        } else {
            assert!(false, "unexpected: exon=[{},{}) cds=[{},{})",
                    start, end, coding_r.0, coding_r.1)
        }
    }

    match transcript_strand {
        &Strand::Forward => set_coding_frames(exons.iter_mut()),
        &Strand::Reverse => set_coding_frames(exons.iter_mut().rev()),
        _ => {}
    }

    Ok(exons)
}

// Helper function for backtracking and adding possibly skipped features
fn backtrack_and_push<F>(
    exons: &mut Vec<Exon>,
    efk: ExonFeatureKind,
    mut codon_rem: u64,
    feature_maker: &F) -> u64
where F: Fn(u64, u64, ExonFeatureKind) -> ExonFeature
{
    for mut exon in exons.iter_mut().rev() {
        if codon_rem == 0 {
            break;
        };
        let fx = feature_maker(max(exon.start(), exon.end() - codon_rem),
                               exon.end(), efk.clone());
        codon_rem -= fx.span();
        let ofxp_start = exon.features.last()
            .and_then(|fxp| {
                match (fxp.kind(), fx.kind()) {
                    (&UTR3, &StopCodon { .. }) => Some(fxp.start()),
                    _ => None,
                }
            });
        if let Some(fxp_start) = ofxp_start {
            let adj_fxp_end = fx.start();
            if fxp_start == adj_fxp_end {
                exon.features.pop();
            } else {
                let n_fxs = exon.features.len();
                let new_fxp_interval = Interval::new(fxp_start..adj_fxp_end).unwrap();
                exon.features[n_fxs-1].interval = new_fxp_interval;
            }
        }
        exon.features.push(fx);
    }
    codon_rem
}

fn set_coding_frames<'a, T>(exons_miter: T)
where T: Iterator<Item=&'a mut Exon>
{
    let (mut startc_frame, mut cds_frame, mut stopc_frame) = (0, 0, 0);
    for mut exon in exons_miter {
        for coding_fx in exon.features.iter_mut() {
            match coding_fx.kind {
                StartCodon { .. } => {
                    coding_fx.kind = StartCodon { frame: Some(startc_frame) };
                    startc_frame = calc_next_frame(coding_fx.span(), startc_frame);
                },
                CDS { .. } => {
                    coding_fx.kind = CDS { frame: Some(cds_frame) };
                    cds_frame = calc_next_frame(coding_fx.span(), cds_frame);
                },
                StopCodon { .. } => {
                    coding_fx.kind = StopCodon { frame: Some(stopc_frame) };
                    stopc_frame = calc_next_frame(coding_fx.span(), stopc_frame);
                },
                _ => {}
            };
        }
    }
}

#[inline(always)]
// Adapted from: http://mblab.wustl.edu/GTF22.html
fn calc_next_frame(cur_span: u64, cur_frame: u8) -> u8 {
    let cast_cur_frame = cur_frame as u64;
    let result =
        if cur_span >= cast_cur_frame {
            (3 - ((cur_span - cast_cur_frame) % 3)) % 3
        } else {
            (3 - (cast_cur_frame - cur_span) % 3)
        };
    result as u8
}
