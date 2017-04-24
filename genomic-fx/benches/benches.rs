#![feature(test)]

extern crate test;
extern crate genomic_fx;
extern crate bio;

use std::collections::HashMap;

use bio::utils::{Interval, Strand};
use test::Bencher;

use genomic_fx::{ExonFeature, ExonFeatureKind, EBuilder, TBuilder, GBuilder};
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
                    ExonFeature {
                        interval: Interval::new(100..150).unwrap(),
                        kind: UTR5,
                    },
                    ExonFeature {
                        interval: Interval::new(150..153).unwrap(),
                        kind: StartCodon { frame: Some(0) },
                    },
                    ExonFeature {
                        interval: Interval::new(150..300).unwrap(),
                        kind: CDS { frame: Some(0) },
                    },
                    ExonFeature {
                        interval: Interval::new(300..303).unwrap(),
                        kind: StopCodon { frame: Some(0) },
                    },
                    ExonFeature {
                        interval: Interval::new(300..500).unwrap(),
                        kind: UTR3,
                    },
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
                .exon_coords(vec![(100, 10000)], None)
                .build()
        });
    }

    #[bench]
    fn tbuilder_1_exon_1_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .exon_coords(vec![(100, 10000)], Some(200, 9500))
                .build()
        });
    }

    #[bench]
    fn tbuilder_12_exons_no_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .exon_coords(
                    vec![
                        (100, 300), (400, 500), (700, 1000),
                        (1100, 1300), (1400, 1500), (1700, 2000),
                        (2100, 2300), (2400, 2500), (2700, 3000),
                        (3000, 6000), (7000, 8000), (9000, 10000),
                    ])
                .build()
        });
    }

    #[bench]
    fn tbuilder_12_exons_with_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .exon_coords(
                    vec![
                        (100, 300), (400, 500), (700, 1000),
                        (1100, 1300), (1400, 1500), (1700, 2000),
                        (2100, 2300), (2400, 2500), (2700, 3000),
                        (3000, 6000), (7000, 8000), (9000, 10000),
                    ], Some(200, 9500))
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
            let mut coords = HashMap::new();
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
}
