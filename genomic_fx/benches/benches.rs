#![feature(test)]

extern crate bio;
extern crate genomic_fx;
extern crate linked_hash_map;
extern crate test;

use test::Bencher;

use bio::utils::{Interval, Strand};
use linked_hash_map::LinkedHashMap;

use genomic_fx::{ExonFeature, ExonFeatureKind, EBuilder, TBuilder, GBuilder, RefFlatReader};
use ExonFeatureKind::*;


#[cfg(test)]
mod tests {
    use super::*;

    #[bench]
    fn ebuilder_no_features(b: &mut Bencher) {
        b.iter(|| {
            EBuilder::new("chrT", 100, 500)
                .strand(Strand::Forward)
                .build()
        });
    }

    #[bench]
    fn ebuiler_5_features(b: &mut Bencher) {
        b.iter(|| {
            EBuilder::new("chrT", 100, 300)
                .strand(Strand::Forward)
                .features(vec![
                    ExonFeature::new(Interval::new(100..150).unwrap(), UTR5),
                    ExonFeature::new(Interval::new(150..153).unwrap(),
                                     StartCodon { frame: Some(0) }),
                    ExonFeature::new(Interval::new(150..300).unwrap(),
                                     CDS { frame: Some(0) }),
                    ExonFeature::new(Interval::new(300..303).unwrap(),
                                     StopCodon { frame: Some(0) }),
                    ExonFeature::new(Interval::new(303..500).unwrap(),
                                     UTR3),
                ])
                .build()
        });
    }

    #[bench]
    fn tbuilder_1_exon_no_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .coords(vec![(100, 10000)], None)
                .build()
        });
    }

    #[bench]
    fn tbuilder_1_exon_1_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .coords(vec![(100, 10000)], Some((200, 9500)))
                .build()
        });
    }

    #[bench]
    fn tbuilder_12_exons_no_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .coords(
                    vec![
                        (100, 300), (400, 500), (700, 1000),
                        (1100, 1300), (1400, 1500), (1700, 2000),
                        (2100, 2300), (2400, 2500), (2700, 3000),
                        (3000, 6000), (7000, 8000), (9000, 10000),
                    ], None)
                .build()
        });
    }

    #[bench]
    fn tbuilder_12_exons_with_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .coords(
                    vec![
                        (100, 300), (400, 500), (700, 1000),
                        (1100, 1300), (1400, 1500), (1700, 2000),
                        (2100, 2300), (2400, 2500), (2700, 3000),
                        (3000, 6000), (7000, 8000), (9000, 10000),
                    ], Some((200, 9500)))
                .build()
        });
    }

    #[bench]
    fn gbuilder_no_transcripts(b: &mut Bencher) {
        b.iter(|| {
            GBuilder::new("chrT", 100, 20000)
                .strand(Strand::Forward)
                .id("gx01")
                .build()
        });
    }

    #[bench]
    fn gbuilder_3_transcripts_with_cds(b: &mut Bencher) {
        b.iter(|| {
            let mut coords = LinkedHashMap::new();
            coords.insert(
                "trx01".to_owned(),
                ((100, 10000),
                vec![
                    (100, 300), (400, 500), (700, 1000),
                    (1100, 1300), (1400, 1500), (1700, 2000),
                    (2100, 2300), (2400, 2500), (2700, 3000),
                    (3000, 6000), (7000, 8000), (9000, 10000)],
                Some((200, 9500))));
            coords.insert(
                "trx02".to_owned(),
                ((100, 10000),
                vec![
                    (100, 300), (400, 500), (700, 1000),
                    (1100, 1300), (1400, 1500), (1700, 2000),
                    (2100, 2300), (2400, 2500), (2700, 3000),
                    (3000, 6000), (7000, 8000), (9000, 10000)],
                Some((200, 9500))));
            coords.insert(
                "trx03".to_owned(),
                ((100, 10000),
                vec![
                    (100, 300), (400, 500), (700, 1000),
                    (1100, 1300), (1400, 1500), (1700, 2000),
                    (2100, 2300), (2400, 2500), (2700, 3000),
                    (3000, 6000), (7000, 8000), (9000, 10000)],
                Some((200, 9500))));
            GBuilder::new("chrT", 100, 20000)
                .strand(Strand::Forward)
                .id("gx01")
                .transcript_coords(coords)
                .build()
        });
    }

    const REFFLAT_MULT_ROWS: &'static [u8] = b"TNFRSF14\tNM_001297605\tchr1\t+\t2556364\t2565622\t2556664\t2562868\t7\t2556364,2557725,2558342,2559822,2560623,2562864,2563147,\t2556733,2557834,2558468,2559978,2560714,2562896,2565622,
TNFRSF14\tNM_003820\tchr1\t+\t2556364\t2565622\t2556664\t2563273\t8\t2556364,2557725,2558342,2559822,2560623,2561672,2562864,2563147,\t2556733,2557834,2558468,2559978,2560714,2561815,2562896,2565622,
SMIM12\tNM_001164824\tchr1\t-\t34850361\t34859045\t34855698\t34855977\t3\t34850361,34856555,34858839,\t34855982,34856739,34859045,
SMIM12\tNM_001164825\tchr1\t-\t34850361\t34859737\t34855698\t34855977\t2\t34850361,34859454,\t34855982,34859737,
SMIM12\tNM_138428\tchr1\t-\t34850361\t34859816\t34855698\t34855977\t2\t34850361,34859676,\t34855982,34859816,";

    #[bench]
    fn refflat_reader_mult_rows_records(b: &mut Bencher) {
        b.iter(|| {
            let mut count = 0;
            let mut reader = RefFlatReader::from_reader(REFFLAT_MULT_ROWS);
            for rec in reader.records_stream() {
                if let Ok(_) = rec {
                    count += 1;
                }
            }
            assert_eq!(count, 5);
        });
    }

    #[bench]
    fn refflat_reader_mult_rows_transcripts(b: &mut Bencher) {
        b.iter(|| {
            let mut count = 0;
            let mut reader = RefFlatReader::from_reader(REFFLAT_MULT_ROWS);
            for trx in reader.transcripts_stream() {
                if let Ok(_) = trx {
                    count += 1;
                }
            }
            assert_eq!(count, 5);
        });
    }

    #[bench]
    fn refflat_reader_mult_rows_genes(b: &mut Bencher) {
        b.iter(|| {
            let mut count = 0;
            let mut reader = RefFlatReader::from_reader(REFFLAT_MULT_ROWS);
            for gx in reader.genes_stream() {
                if let Ok(_) = gx {
                    count += 1;
                }
            }
            assert_eq!(count, 2);
        });
    }
}
