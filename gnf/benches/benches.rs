#![feature(test)]

extern crate test;
extern crate gnf;
extern crate bio;

use bio::utils::Strand;
use gnf::feature::TranscriptBuilder;
use test::Bencher;


#[cfg(test)]
mod tests {
    use super::*;

    #[bench]
    fn transcript_builder_s_exon_0_cds(b: &mut Bencher) {
        b.iter(|| {
            TranscriptBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .exon_and_cds_coords(
                    vec![(100, 10000)], None)
                .build()
        });
    }

    #[bench]
    fn transcript_builder_s_exon_s_cds(b: &mut Bencher) {
        b.iter(|| {
            TranscriptBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .exon_and_cds_coords(
                    vec![(100, 10000)], Some((200, 9500)))
                .build()
        });
    }

    #[bench]
    fn transcript_builder_m_exons_0_cds(b: &mut Bencher) {
        b.iter(|| {
            TranscriptBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .exon_and_cds_coords(
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
    fn transcript_builder_m_exons_m_cds(b: &mut Bencher) {
        b.iter(|| {
            TranscriptBuilder::new("chrT", 100, 10000)
                .strand(Strand::Forward)
                .exon_and_cds_coords(
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
