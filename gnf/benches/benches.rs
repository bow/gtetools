#![feature(test)]

extern crate test;
extern crate gnf;
extern crate bio;

use bio::utils::Strand;
use gnf::TBuilder;
use test::Bencher;


#[cfg(test)]
mod tests {
    use super::*;

    #[bench]
    fn tbuilder_single_exon_no_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .feature_coords(
                    vec![(100, 10000)], None)
                .build()
        });
    }

    #[bench]
    fn tbuilder_single_exon_single_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .feature_coords(
                    vec![(100, 10000)], Some((200, 9500)))
                .build()
        });
    }

    #[bench]
    fn tbuilder_mult_exons_no_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .feature_coords(
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
    fn tbuilder_mult_exons_mult_cds(b: &mut Bencher) {
        b.iter(|| {
            TBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .id("trx01")
                .feature_coords(
                    vec![
                        (100, 300), (400, 500), (700, 1000),
                        (1100, 1300), (1400, 1500), (1700, 2000),
                        (2100, 2300), (2400, 2500), (2700, 3000),
                        (3000, 6000), (7000, 8000), (9000, 10000),
                    ],
                    Some((200, 9500)))
                .build()
        });
    }
}
