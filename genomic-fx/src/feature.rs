use std::cmp::{max, min};
use std::collections::HashMap;

use bio::utils::{self, Interval, IntervalError};
use bio::utils::Strand;

use Error;
use self::ExonFeatureKind::*;


macro_rules! impl_common {
    ($struct_ty:ty) => (

        impl $struct_ty {

            pub fn seq_name(&self) -> &str {
                self.seq_name.as_str()
            }

            pub fn interval(&self) -> &Interval<u64> {
                &self.interval
            }

            pub fn strand(&self) -> &Strand {
                &self.strand
            }

            #[inline]
            pub fn span(&self) -> u64 {
                self.interval().end - self.interval().start
            }
        }

    );
}

#[derive(Debug, Clone, PartialEq)]
pub struct Feature<K: FeatureKind> {
    pub interval: Interval<u64>,
    pub kind: K,
}

impl<K: FeatureKind> Feature<K> {

    #[inline]
    fn span(&self) -> u64 {
        self.interval.end - self.interval.start
    }
}

pub trait FeatureKind {}

#[derive(Debug, Clone, PartialEq)]
pub enum ExonFeatureKind {
    UTR,
    UTR5,
    UTR3,
    CDS { frame: Option<u64> },
    StartCodon { frame: Option<u64> },
    StopCodon { frame: Option<u64> },
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

#[derive(Debug)]
pub struct Exon {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    pub id: Option<String>,
    pub attributes: HashMap<String, String>,
    pub features: Vec<ExonFeature>,
}

impl_common!(Exon);

pub struct EBuilder {
    seq_name: String,
    start: u64,
    end: u64,
    strand: Option<Strand>,
    strand_char: Option<char>,
    pub id: Option<String>,
    pub attributes: HashMap<String, String>,
    pub features: Vec<ExonFeature>,
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
            attributes: HashMap::new(),
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

    pub fn attribute<K, V>(mut self, key: K, value: V) -> Self
        where K: Into<String>, V: Into<String>
    {
        self.attributes.insert(key.into(), value.into());
        self
    }

    pub fn attributes(mut self, attributes: HashMap<String, String>) -> Self {
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

    pub fn build(self) -> Result<Exon, Error> {
        let interval = coords_to_interval(self.start, self.end)
            .map_err(Error::Feature)?;
        let strand = resolve_strand_input(self.strand, self.strand_char)
            .map_err(Error::Feature)?;
        let feature = Exon {
            seq_name: self.seq_name,
            interval: interval,
            strand: strand,
            id: self.id,
            attributes: self.attributes,
            features: self.features,
        };
        Ok(feature)
    }
}

#[derive(Debug)]
pub struct Transcript {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    pub id: Option<String>,
    pub attributes: HashMap<String, String>,
    exons: Vec<Exon>,
}

impl_common!(Transcript);

impl Transcript {

    pub fn exons(&self) -> &[Exon] {
        self.exons.as_slice()
    }
}

pub struct TBuilder {
    seq_name: String,
    start: u64,
    end: u64,
    strand: Option<Strand>,
    strand_char: Option<char>,
    id: Option<String>,
    attributes: HashMap<String, String>,
    // Input can be a vector of pre-made features ...
    exons: Option<Vec<Exon>>,
    // Or exon coordinates, possibly coupled with cds coord
    // NOTE: Can we instead of using Vec<_> here keep it as an unconsumed iterator?
    exon_coords: Option<Vec<(u64, u64)>>,
    coding_coord: Option<(u64, u64)>,
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
            attributes: HashMap::new(),
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

    pub fn attribute<K, V>(mut self, key: K, value: V) -> Self
        where K: Into<String>, V: Into<String>
    {
        self.attributes.insert(key.into(), value.into());
        self
    }

    pub fn attributes(mut self, attributes: HashMap<String, String>) -> Self {
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

    pub fn exon_coords<E>(mut self, exon_coords: E)-> Self
        where E: IntoIterator<Item=(u64, u64)>
    {
        self.exon_coords = Some(exon_coords.into_iter().collect());
        self
    }

    pub fn coding_coord(mut self, start: u64, end: u64) -> Self {
        self.coding_coord = Some((start, end));
        self
    }

    pub fn coding_incl_stop(mut self, incl_stop: bool) -> Self {
        self.coding_incl_stop = incl_stop;
        self
    }

    pub fn build(self) -> Result<Transcript, Error> {
        let interval = coords_to_interval(self.start, self.end)
            .map_err(Error::Feature)?;
        let strand = resolve_strand_input(self.strand, self.strand_char)
            .map_err(Error::Feature)?;
        let exons = resolve_exons_input(
            &self.seq_name, &interval, &strand,
            self.exons, self.exon_coords.as_ref(), self.coding_coord,
            self.coding_incl_stop).map_err(Error::Feature)?;

        let transcript = Transcript {
            seq_name: self.seq_name,
            interval: interval,
            strand: strand,
            id: self.id,
            attributes: self.attributes,
            exons: exons,
        };
        Ok(transcript)
    }
}

#[derive(Debug)]
pub struct Gene {
    seq_name: String,
    interval: Interval<u64>,
    strand: Strand,
    id: Option<String>,
    attributes: HashMap<String, String>,
    transcripts: HashMap<String, Transcript>,
}

impl_common!(Gene);

quick_error! {
    #[derive(Debug)]
    pub enum FeatureError {
        InvalidInterval(err: utils::IntervalError) {
            description(
                match err {
                    &IntervalError::InvalidRange =>
                        "interval start coordinate larger than its end coordinate",
                    ref otherwise => otherwise.description(),
                })
            from()
        }
        InvalidStrandChar(err: utils::StrandError) {
            description(err.description())
            from()
        }
        ConflictingStrand {
            description("conflicting strand inputs specified")
        }
        UnspecifiedStrand {
            description("strand not specified")
        }
        InvalidExonInterval {
            description("exon has larger start than end coordinate")
        }
        InvalidCodingInterval {
            description("coding region has larger start than end coordinate")
        }
        UnspecifiedExons {
            description("transcript is defined without exons")
        }
        UnmatchedExons {
            description("first and/or last exon coordinates do not match transcript \
                         start and/or end coordinates")
        }
        CodingTooLarge {
            description("coding region leaves no room for stop codon in transcript")
        }
        CodingTooSmall {
            description("coding region leaves no room for start codon")
        }
        CodingNotFullyEnveloped {
            description("coding region not fully enveloped by exons")
        }
        CodingInIntron {
            description("coding start and/or end lies in introns")
        }
    }
}

fn coords_to_interval(start: u64, end: u64) -> Result<Interval<u64>, FeatureError> {
    Interval::new(start..end).map_err(FeatureError::from)
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
    exons: Option<Vec<Exon>>,
    exon_coords: Option<&Vec<(u64, u64)>>,
    coding_coord: Option<(u64, u64)>,
    coding_incl_stop: bool
) -> Result<Vec<Exon>, FeatureError>
{
    // Deliberately not handling all possible input types to avoid
    // overcomplicating code. The inputs are expected to come from
    // either GTF or refFlat after all.

    match (exons, exon_coords, coding_coord) {
        // nothing defined -> the transcript doesn't have any known exons
        (None, None, None) => Ok(Vec::new()),

        // only CDS defined -> must be an error
        (None, None, Some(_)) => Err(FeatureError::UnspecifiedExons),

        // features defined ~ takes precedence over coords (GTF input, since we need
        // to construct the tx features first to store its annotations)
        // TODO: Maybe do some checks to ensure the given features are correct?
        (Some(exns), _, _) => Ok(exns.into_iter().collect()),

        // exon defined & coords possibly defined (refFlat input)
        (None, Some(raw_exon_coords), raw_coding_coord) =>
            infer_exons(transcript_seqname, transcript_interval,
                        transcript_strand, raw_exon_coords, raw_coding_coord,
                        coding_incl_stop),
    }
}

fn infer_exons(
    transcript_seqname: &String,
    transcript_interval: &Interval<u64>,
    transcript_strand: &Strand,
    exon_coords: &Vec<(u64, u64)>,
    coding_coord: Option<(u64, u64)>,
    coding_incl_stop: bool
) -> Result<Vec<Exon>, FeatureError>
{

    if exon_coords.len() == 0 {
        return Err(FeatureError::UnspecifiedExons);
    }

    let mut m_exon_coords = Vec::with_capacity(exon_coords.len());
    for &(a, b) in exon_coords.iter() {
        if a >= b {
            return Err(FeatureError::InvalidExonInterval)
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
        return Err(FeatureError::UnmatchedExons);
    }

    match adj_coding_coord {

        Some(coding_r) => {
            // Improper coding region is an error
            if coding_r.0 >= coding_r.1 {
                return Err(FeatureError::InvalidCodingInterval);
            }
            // Coding coord must be fully enveloped by exon max-min
            if coding_r.0 < exon_r.0 || coding_r.1 > exon_r.1 {
                return Err(FeatureError::CodingNotFullyEnveloped);
            }
            // Coding start and end must be in exons
            let cine = m_exon_coords.iter()
                .fold((false, false), |acc, c| {
                    (acc.0 || (c.0 <= coding_r.0 && coding_r.0 <= c.1),
                     acc.1 || (c.0 <= coding_r.1 && coding_r.1 <= c.1))
                });
            if !cine.0 || !cine.1 {
                return Err(FeatureError::CodingInIntron);
            }
            // There must be room for stop codons (which is not inclusive in coding_coord)
            let stop_codon_ok = match transcript_strand {
                &Strand::Forward => coding_r.1 + 3 <= exon_r.1,
                &Strand::Reverse => coding_r.0 - 3 >= exon_r.0,
                &Strand::Unknown =>
                    coding_r.0 - 3 >= exon_r.0 && coding_r.1 + 3 <= exon_r.1,
            };
            if !stop_codon_ok {
                return Err(FeatureError::CodingTooLarge);
            }
            infer_exon_features(&m_exon_coords, coding_r, &transcript_seqname, transcript_strand)
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
                        id: None,
                        attributes: HashMap::new(),
                        features: Vec::new(),
                    });
            }
            Ok(features)
        }
    }
}

fn adjust_coding_coord(mut start: u64, mut end: u64,
                       strand: &Strand, exon_coords: &Vec<(u64, u64)>
) -> Option<(u64, u64)>
{
    let mut codon_rem = 3;
    match strand {
        &Strand::Forward => {
            for &(exon_start, exon_end) in exon_coords.iter().rev() {
                if end <= exon_end {
                    let adj_end = max(end - codon_rem, exon_start);
                    codon_rem -= end - adj_end;
                    end = adj_end;
                }
            }
        },
        &Strand::Reverse => {
            for &(exon_start, exon_end) in exon_coords.iter() {
                if exon_start <= start {
                    let adj_start = min(start + codon_rem, exon_end);
                    codon_rem -= adj_start - start;
                    start = adj_start;
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
    exon_coords: &Vec<(u64, u64)>,
    coding_r: (u64, u64),
    transcript_seqname: &String,
    transcript_strand: &Strand)
-> Result<Vec<Exon>, FeatureError> {

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
            id: None,
            attributes: HashMap::new(),
            features: features,
        }
    };
    let feat = |start, end, kind| {
        ExonFeature {
            interval: Interval::new(start..end).unwrap(),
            kind: kind,
        }
    };

    // how much we have consumed the 5' or 3' codon
    let (mut codon1_rem, mut codon2_rem) = (3, 3);
    for &(start, end) in exon_coords.iter() {

        if start < coding_r.0 {

            if end < coding_r.0 {
                let exon = exn(start, end, vec![feat(start, end, utr1.clone())]);
                exons.push(exon);

            } else if end == coding_r.0 {
                let mut exon = exn(start, end, vec![feat(start, end, utr1.clone())]);

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
                let mut exon = exn(start, end, vec![feat(start, coding_r.0, utr1.clone())]);
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
                    return Err(FeatureError::CodingTooSmall);
                }
                let mut exon = exn(start, end, vec![feat(start, coding_r.0, utr1.clone())]);
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
                    return Err(FeatureError::CodingTooSmall);
                }
                let mut exon = exn(start, end, vec![feat(start, coding_r.0, utr1.clone())]);
                match transcript_strand {
                    &Strand::Forward => {
                        exon.features.push(feat(coding_r.0, coding_r.0 + 3,
                                                StartCodon { frame: None }));
                        codon1_rem -= 3;
                        exon.features.push(feat(coding_r.0, coding_r.1, CDS { frame: None }));
                        let fx = feat(coding_r.1, min(end, coding_r.1 + codon2_rem),
                                      StopCodon { frame: None });
                        codon2_rem -= fx.span();
                        exon.features.push(fx);
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
                    },
                    &Strand::Unknown => {
                        exon.features.push(feat(coding_r.0, coding_r.1, CDS { frame: None }));
                    },
                }
                exon.features.push(feat(coding_r.1, end, utr2.clone()));
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
                    return Err(FeatureError::CodingTooSmall);
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
                    return Err(FeatureError::CodingTooSmall);
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
                        exon.features.push(fx);
                    },
                    &Strand::Reverse => {
                        codon1_rem = backtrack_and_push(&mut exons, StopCodon { frame: None },
                                                        codon1_rem, &feat);
                        exon.features.push(feat(start, coding_r.1, CDS { frame: None }));
                        exon.features.push(feat(coding_r.1 - 3, coding_r.1,
                                                StartCodon { frame: None }));
                        codon2_rem -= 3;
                    },
                    &Strand::Unknown => {
                        exon.features.push(feat(start, coding_r.1, CDS { frame: None }));
                    },
                }
                exon.features.push(feat(coding_r.1, end, utr2.clone()));
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
                        exon.features.push(fx);
                        exon.features.push(feat(coding_r.1, end, utr2.clone()));
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
                        exon.features.push(fx);
                    }
                },
                &Strand::Reverse => {
                    codon2_rem = backtrack_and_push(
                        &mut exons, StartCodon { frame: None }, codon2_rem, &feat);
                },
                &Strand::Unknown => {},
            }
            exon.features.push(feat(start, end, utr2.clone()));
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
        let fx = feature_maker(max(exon.interval().start, exon.interval().end - codon_rem),
                               exon.interval().end, efk.clone());
        codon_rem -= fx.span();
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
fn calc_next_frame(cur_span: u64, cur_frame: u64) -> u64 {
    if cur_span >= cur_frame {
        (3 - ((cur_span - cur_frame) % 3)) % 3
    } else {
        (3 - (cur_frame - cur_span) % 3)
    }
}
