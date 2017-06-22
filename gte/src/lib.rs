#![deny(
        trivial_casts, trivial_numeric_casts,
        unsafe_code,
        unstable_features,
        unused_extern_crates, unused_import_braces, unused_qualifications)]
#![warn(unused_results)]

extern crate bio;
extern crate csv;
extern crate itertools;
extern crate linked_hash_map;
extern crate multimap;
#[macro_use]
extern crate quick_error;
extern crate regex;

use std::num::{ParseFloatError, ParseIntError};

pub use bio::utils::Strand;
pub use bio::io::gff::GffType;
use csv::Error as CsvError;
use regex::Error as RegexError;

mod model;
pub use model::{Feature, ModelError, FeatureKind,
                EBuilder, Exon, ExonFeature, ExonFeatureKind,
                TBuilder, Transcript, TranscriptFeature, TranscriptFeatureKind,
                GBuilder, Gene, GeneFeature, GeneFeatureKind};

mod io_refflat;
pub use io_refflat::{Reader as RefFlatReader, Writer as RefFlatWriter,
                     RefFlatError, RefFlatRow, RefFlatRecord,
                     RefFlatRecordsStream, RefFlatTranscriptsStream, RefFlatGenesStream};

mod io_gff;
pub use io_gff::{Reader as GffReader, GffError, GffTranscripts};


quick_error! {
    #[derive(Debug)]
    pub enum Error {
        Model(err: ModelError) {
            description(err.description())
            display("{}", err)
            from()
            cause(err)
        }
        RefFlat(err: RefFlatError) {
            description(err.description())
            display("{}", err)
            from()
            cause(err)
        }
        Gff(err: GffError) {
            description(err.description())
            display("{}", err)
            from()
            cause(err)
        }
        Regex(err: RegexError) {
            description(err.description())
            from()
            cause(err)
        }
        Csv(err: CsvError) {
            description(err.description())
            from()
            cause(err)
        }
        ParseInt(err: ParseIntError) {
            description(err.description())
            from()
            cause(err)
        }
        ParseFloat(err: ParseFloatError) {
            description(err.description())
            from()
            cause(err)
        }
    }
}

pub type Result<T> = ::std::result::Result<T, Error>;

// Helper type for raw coordinates
type Coord<T> = (T, T);

// Helper type alias for raw transcript coordinate inputs
type RawTrxCoords = (Coord<u64>, Vec<Coord<u64>>, Option<Coord<u64>>);

// Crate-wide constants
mod consts {
    // Initial coordinate for features.
    pub(crate) const INIT_START: u64 = ::std::u64::MAX;
    pub(crate) const INIT_END: u64 = ::std::u64::MIN;
    pub(crate) const INIT_COORD: (u64, u64) = (INIT_START, INIT_END);

    // Various commonly-used feature column values
    pub(crate) const GENE_STR: &'static str = "gene";
    pub(crate) const TRANSCRIPT_STR: &'static str = "transcript";
    pub(crate) const EXON_STR: &'static str = "exon";
    pub(crate) const UTR_STR: &'static str = "UTR";
    pub(crate) const UTR5_STR: &'static str = "UTR5";
    pub(crate) const UTR3_STR: &'static str = "UTR3";
    pub(crate) const CDS_STR: &'static str = "CDS";
    pub(crate) const START_CODON_STR: &'static str = "start_codon";
    pub(crate) const STOP_CODON_STR: &'static str = "stop_codon";

    // Value for unknown columns.
    pub(crate) const UNK_STR: &'static str = ".";
    pub(crate) const UNK_CHAR: char = '.';

    // Commonly-used attribute keys.
    pub(crate) const GENE_ID_STR: &'static str = "gene_id";
    pub(crate) const TRANSCRIPT_ID_STR: &'static str = "transcript_id";

    // Value for optionally known strings.
    pub(crate) const DEF_ID: &'static str = "<unknown>";
}

// Generic utilities
mod utils {
    use std::ops::Deref;

    // taken from: https://stackoverflow.com/q/31233938/243058
    pub(crate) trait OptionDeref<T: Deref> {
        fn as_deref(&self) -> Option<&T::Target>;
    }

    impl<T: Deref> OptionDeref<T> for Option<T> {
        fn as_deref(&self) -> Option<&T::Target> {
            self.as_ref().map(Deref::deref)
        }
    }

    #[inline]
    pub(crate) fn update_contig<'a>(
        value: &'a mut String,
        prefix: Option<&'a str>,
        lstrip: Option<(&'a str, usize)>,
    ) {
        if let Some(ref pre) = prefix {
            value.insert_str(0, pre);
        }
        if let Some((ref lstr, lstr_len)) = lstrip {
            if value.starts_with(lstr) {
                let _ = value.drain(..lstr_len);
            }
        }
    }
}
