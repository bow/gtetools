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

pub use bio::utils::Strand;
pub use bio::io::gff::GffType;

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
    }
}

pub type Result<T> = ::std::result::Result<T, Error>;

// Helper type for raw coordinates
type Coord<T> = (T, T);

// Helper type alias for raw transcript coordinate inputs
type RawTrxCoords = (Coord<u64>, Vec<Coord<u64>>, Option<Coord<u64>>);

/// Fallback identifier for when a string value is required.
const DEF_ID: &'static str = "<unknown>";

/// Initial start coordinate value.
///
/// This is meant to be used with `std::cmp::min` as coordinates from an input is parsed.
const INIT_START: u64 = ::std::u64::MAX;

/// Initial end coordinate value.
///
/// This is meant to be used with `std::cmp::max` as coordinates from an input is parsed.
const INIT_END: u64 = ::std::u64::MIN;

/// Initial start and end coordinates values.
const INIT_COORD: (u64, u64) = (INIT_START, INIT_END);

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
    pub(crate) fn update_seq_name<'a>(
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
