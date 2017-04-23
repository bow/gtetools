#![feature(test)]

extern crate test;
extern crate genomic_fx;
extern crate bio;

use bio::utils::{Interval, Strand};
use test::Bencher;

use genomic_fx::{EBuilder, TBuilder, ExonFeature, ExonFeatureKind};
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
    fn ebuiler_mult_features(b: &mut Bencher) {
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
    fn tbuilder_single_exon_no_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .exon_coords(vec![(100, 10000)])
                .build()
        });
    }

    #[bench]
    fn tbuilder_single_exon_single_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .exon_coords(vec![(100, 10000)])
                .coding_coord(200, 9500)
                .build()
        });
    }

    #[bench]
    fn tbuilder_mult_exons_no_cds(b: &mut Bencher) {
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
    fn tbuilder_mult_exons_mult_cds(b: &mut Bencher) {
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
                .coding_coord(200, 9500)
                .build()
        });
    }
}
