/*!
The `gte` crate provides functionalities for working with genes, transcripts, and exons.

It defines simple structs for genes, transcripts, and exons, along with builders of these
structures that accept a flexible range of arguments. You can create these structs on your own
or from formats such as GFF and refFlat which are commonly used for storing gene annotations.

*/
#![deny(missing_docs,
        trivial_casts, trivial_numeric_casts,
        unsafe_code,
        unstable_features,
        unused_extern_crates, unused_import_braces, unused_qualifications)]
#![warn(unused_results)]
#![recursion_limit="128"]

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
    /// The error type returned by the `gte` crate.
    #[derive(Debug)]
    pub enum Error {
        /// Errors that occur when building the exon, transcript, or gene model.
        Model(err: ModelError) {
            description(err.description())
            display("{}", err)
            from()
            cause(err)
        }
        /// Errors that occur when reading or writing refFlat files.
        RefFlat(err: RefFlatError) {
            description(err.description())
            display("{}", err)
            from()
            cause(err)
        }
        /// Errors that occur when reading or writing GFF file variants.
        Gff(err: GffError) {
            description(err.description())
            display("{}", err)
            from()
            cause(err)
        }
    }
}

/// Result type whose error variant is bound to [`gte::Error`].
///
/// [`gte::Error`]: enum.Error.html
pub type Result<T> = ::std::result::Result<T, Error>;

/// Helper type alias for raw coordinates.
type Coord<T> = (T, T);

/// Helper type alias for raw transcript coordinate inputs.
///
/// This is an alias for a three-element tuple which contains:
/// * The coordinate of the full transcript (as Coord<u64>).
/// * The coordinates of all the transcript's exons (as (Vec<Coord<u64>>)).
/// * The optional coordinate of the coding region, which may or may not include
///   the stop codon.
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

/// Utility functions.
mod utils {
    use std::ops::Deref;

    // taken from: https://stackoverflow.com/q/31233938/243058
    /// Helper trait for dereferencing wrapped option values.
    ///
    /// This is mainly used by `gte` to read `Option<String>` values as `Option<&str>`.
    pub(crate) trait OptionDeref<T: Deref> {
        fn as_deref(&self) -> Option<&T::Target>;
    }

    impl<T: Deref> OptionDeref<T> for Option<T> {
        fn as_deref(&self) -> Option<&T::Target> {
            self.as_ref().map(Deref::deref)
        }
    }

    /// Helper function for in-place prefixing and/or left-stripping of owned strings.
    ///
    /// This is used by `gte` to prepend and/or remove strings to sequence names.
    ///
    /// Arguments:
    /// * `value`: Mutable reference to `String` to update.
    /// * `prefix`: Optional value of `&str` to prefix. If set to `None`, no prefix will be added.
    /// * `lstrip`: Optional value of `&str` to left-strip. If set to `None`, no left-stripping
    ///             will be done.
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
