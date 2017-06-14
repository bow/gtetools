extern crate bio;
extern crate csv;
extern crate itertools;
extern crate linked_hash_map;
extern crate multimap;
#[macro_use]
extern crate quick_error;
extern crate regex;

use std::io::Error as StdIoError;
use std::num::{ParseFloatError, ParseIntError};

pub use bio::utils::Strand;
pub use bio::io::gff::GffType;
use bio::io::gff;
use csv::Error as CsvError;
use regex::Error as RegexError;

mod model;
pub use model::{Feature, FeatureError, FeatureKind,
                EBuilder, Exon, ExonFeature, ExonFeatureKind,
                TBuilder, Transcript, TranscriptFeature, TranscriptFeatureKind,
                GBuilder, Gene, GeneFeature, GeneFeatureKind};

mod io_refflat;
pub use io_refflat::{Reader as RefFlatReader, Writer as RefFlatWriter,
                     RefFlatRow, RefFlatRecord,
                     RefFlatRecordsStream, RefFlatTranscriptsStream, RefFlatGenesStream};

mod io_gff;
pub use io_gff::{Reader as GffReader,
                 GffError, GffGenesStream, GffTranscriptsStream, GffTranscripts};


quick_error! {
    #[derive(Debug)]
    pub enum Error {
        Feature(err: FeatureError) {
            description(err.description())
            from()
            cause(err)
        }
        RefFlat(string: &'static str) {
            description(string)
        }
        Gff(err: GffError) {
            description(err.description())
            from()
            cause(err)
        }
        BioGff(err: gff::GffError) {
            description(err.description())
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
        Io(err: StdIoError) {
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
type RawTrxCoord = (Coord<u64>, Vec<Coord<u64>>, Option<Coord<u64>>);

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
}
