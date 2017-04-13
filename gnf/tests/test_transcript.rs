extern crate bio;
extern crate gnf;

use gnf::{ExonFeatureKind, Strand, TBuilder, Transcript};
use ExonFeatureKind::*;

fn exon_coords(transcript: &Transcript) -> Vec<(u64, u64)> {
    transcript.exons().iter()
        .map(|exn| (exn.interval().start, exn.interval().end))
        .collect()
}

fn exonf_coords(transcript: &Transcript) -> Vec<Vec<(u64, u64, ExonFeatureKind)>> {
    transcript.exons().iter()
        .map(|exn| {
            exn.features.iter()
                .map(|fx| (fx.interval.start, fx.interval.end, fx.kind.clone()))
                .collect()
        })
        .collect()
}

#[test]
fn tbuilder_coords_fwd() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert!(exonf_coords(&trx).iter().all(|fx| fx.len() == 0));
}

#[test]
fn tbuilder_coords_fwd_coding() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(200, 800)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 200, UTR5), (200, 203, StartCodon { frame: Some(0) }),
                            (200, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(2) })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: Some(1) }),
                            (800, 803, StopCodon { frame: Some(0) }), (800, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_from_exon5end() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(400, 900)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 403, StartCodon { frame: Some(0) }),
                            (400, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 900, CDS { frame: Some(2) }),
                            (900, 903, StopCodon { frame: Some(0) }), (900, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_from_exon3end() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(300, 900)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 403, StartCodon { frame: Some(0) }),
                            (400, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 900, CDS { frame: Some(2) }),
                            (900, 903, StopCodon { frame: Some(0) }), (900, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_from_near_exon3end() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(297, 800)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 297, UTR5), (297, 300, StartCodon { frame: Some(0) }),
                            (297, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: Some(2) }),
                            (800, 803, StopCodon { frame: Some(0) }), (800, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_from_split() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(298, 901)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 298, UTR5), (298, 300, StartCodon { frame: Some(0) }),
                            (298, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 401, StartCodon { frame: Some(1) }),
                            (400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 901, CDS { frame: Some(0) }),
                            (901, 904, StopCodon { frame: Some(0) }), (901, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_to_exon3end() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(190, 500)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 190, UTR5), (190, 193, StartCodon { frame: Some(0) }),
                            (190, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0)}), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_to_exon5end() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(190, 700)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 190, UTR5), (190, 193, StartCodon { frame: Some(0) }),
                            (190, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0) }), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_to_split() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(199, 499)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 199, UTR5), (199, 202, StartCodon { frame: Some(0) }),
                            (199, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 499, CDS { frame: Some(1) }),
                            (499, 500, StopCodon { frame: Some(0) }), (499, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 702, StopCodon { frame: Some(2) }), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_rev() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert!(exonf_coords(&trx).iter().all(|fx| fx.len() == 0));
}

#[test]
fn tbuilder_coords_rev_coding() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(200, 800)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 200, UTR3), (197, 200, StopCodon { frame: Some(0) }),
                            (200, 300, CDS { frame: Some(1) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(2) })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: Some(0) }),
                            (797, 800, StartCodon { frame: Some(0) }),
                            (800, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_from_exon3end() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(300, 900)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 900, CDS { frame: Some(0) }),
                            (897, 900, StartCodon { frame: Some(0) }), (900, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_from_exon5end() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(400, 900)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 900, CDS { frame: Some(0) }),
                            (897, 900, StartCodon { frame: Some(0) }), (900, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_from_split() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(401, 901)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (298, 300, StopCodon { frame: Some(2) })]);
    assert_eq!(fxs[1], vec![(400, 401, UTR3), (400, 401, StopCodon { frame: Some(0) }),
                            (401, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 901, CDS { frame: Some(0) }),
                            (898, 901, StartCodon { frame: Some(0) }), (901, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_to_exon5end() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(190, 700)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 190, UTR3), (187, 190, StopCodon { frame: Some(0) }),
                            (190, 300, CDS { frame: Some(2) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(0) }),
                            (497, 500, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_to_near_exon5end() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(200, 703)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 200, UTR3), (197, 200, StopCodon { frame: Some(0) }),
                            (200, 300, CDS { frame: Some(2) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 703, CDS { frame: Some(0) }),
                            (700, 703, StartCodon { frame: Some(0) }), (703, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_to_exon3end() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(190, 500)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 190, UTR3), (187, 190, StopCodon { frame: Some(0) }),
                            (190, 300, CDS { frame: Some(2) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(0) }),
                            (497, 500, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_to_split() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(199, 702)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 199, UTR3), (196, 199, StopCodon { frame: Some(0) }),
                            (199, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) }),
                            (499, 500, StartCodon { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 702, CDS { frame: Some(0) }),
                            (700, 702, StartCodon { frame: Some(0) }), (702, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_coding_unk() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Unknown)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .coding_coord(200, 800)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    let fxs = exonf_coords(&trx);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 200, UTR), (200, 300, CDS { frame: None })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: None })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: None }), (800, 1000, UTR)]);
}
