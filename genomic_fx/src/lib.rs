extern crate bio;
extern crate csv;
extern crate itertools;
extern crate linked_hash_map;
#[macro_use]
extern crate quick_error;

use std::io::Error as StdIoError;
use std::num::ParseIntError;

pub use bio::utils::Strand;
pub use bio::io::gff::GffType;
use csv::Error as CsvError;

mod feature;
pub use feature::{Feature, FeatureError, FeatureKind,
                  EBuilder, Exon, ExonFeature, ExonFeatureKind,
                  TBuilder, Transcript, TranscriptFeature, TranscriptFeatureKind,
                  GBuilder, Gene, GeneFeature, GeneFeatureKind};

mod io_refflat;
pub use io_refflat::{Reader as RefFlatReader, Writer as RefFlatWriter,
                     RefFlatRow, RefFlatRecord,
                     RefFlatRecords, RefFlatTranscripts, RefFlatGenes};

mod io_gff;
pub use io_gff::{Reader as GffReader,
                 GffGenes, GffTranscripts};

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
        Gff(string: &'static str) {
            description(string)
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
    }
}

// Helper type for raw coordinates
type Coord<T> = (T, T);

// Helper type alias for raw transcript coordinate inputs
type RawTrxCoord = (Coord<u64>, Vec<Coord<u64>>, Option<Coord<u64>>);
