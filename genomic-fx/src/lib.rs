extern crate bio;
extern crate csv;
extern crate itertools;
#[macro_use]
extern crate quick_error;

use std::io::Error as StdIoError;

pub use bio::utils::Strand;
use csv::Error as CsvError;

mod feature;
pub use feature::{Feature, FeatureError, FeatureKind,
                  EBuilder, Exon, ExonFeature, ExonFeatureKind,
                  TBuilder, Transcript, TranscriptFeature, TranscriptFeatureKind,
                  GBuilder, Gene, GeneFeature, GeneFeatureKind};

mod io_refflat;
pub use io_refflat::{Reader as RefFlatReader, RefFlatRow, RefFlatRecord,
                     RefFlatRecords, RefFlatTranscripts, RefFlatGenes};

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
    }
}

// Helper type for raw coordinates
type Coord<T> = (T, T);
