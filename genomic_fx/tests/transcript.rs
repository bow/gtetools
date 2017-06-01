extern crate bio;
extern crate genomic_fx;

use std::collections::HashMap;

use genomic_fx::{ExonFeatureKind, Strand, TBuilder, Transcript};
use ExonFeatureKind::*;
use Strand::*;

fn exon_coords(transcript: &Transcript) -> Vec<(u64, u64)> {
    transcript.exons().iter()
        .map(|exn| (exn.start(), exn.end()))
        .collect()
}

fn exon_fxs_coords(transcript: &Transcript) -> Vec<Vec<(u64, u64, ExonFeatureKind)>> {
    transcript.exons().iter()
        .map(|exn| {
            exn.features.iter()
                .map(|fx| (fx.start(), fx.end(), fx.kind().clone()))
                .collect()
        })
        .collect()
}

fn trx_fxs<T>(start: u64, end: u64, strand: Strand, exon_coords: T,
               coding_coord: Option<(u64, u64)>
) -> (Transcript, Vec<Vec<(u64, u64, ExonFeatureKind)>>)
where T: IntoIterator<Item=(u64, u64)>
{
    let trx = TBuilder::new("chrT", start, end)
            .strand(strand)
            .coords(exon_coords, coding_coord)
            .build().unwrap();
    let fxs = exon_fxs_coords(&trx);
    (trx, fxs)
}

#[test]
fn tbuilder_basic() {
    let mut attribs = HashMap::new();
    attribs.insert("key1".to_owned(), "value1".to_owned());
    attribs.insert("key2".to_owned(), "value2".to_owned());
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Reverse)
        .coords(vec![(100, 300), (400, 500), (700, 1000)], None)
        .attributes(attribs)
        .id("transcript-1")
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(trx.span(), 900);
    assert_eq!(trx.seq_name(), "chrT");
    assert_eq!(trx.strand(), &Strand::Reverse);
    assert_eq!(trx.id, Some("transcript-1".to_owned()));
    assert_eq!(trx.attributes().get("key1"), Some(&"value1".to_owned()));
    assert_eq!(trx.attributes().get("key2"), Some(&"value2".to_owned()));
    assert_eq!(trx.attributes().len(), 2);
    assert_eq!(trx.exons().len(), 3);
}

#[test]
fn tbuilder_alt1() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Reverse)
        .coords(vec![(100, 300), (400, 500), (700, 1000)], None)
        .attribute("tag", "basic")
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(trx.span(), 900);
    assert_eq!(trx.seq_name(), "chrT");
    assert_eq!(trx.strand(), &Strand::Reverse);
    assert_eq!(trx.id, None);
    assert_eq!(trx.attributes().get("tag"), Some(&"basic".to_owned()));
    assert_eq!(trx.attributes().len(), 1);
    assert_eq!(trx.exons().len(), 3);
}

// Forward strand cases

#[test]
fn tbuilder_coords_fwd() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)], None);
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert!(fxs.iter().all(|fx| fx.len() == 0));
}

#[test]
fn tbuilder_coords_fwd_coding_1_1() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((150, 210)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 150, UTR5),
                            (150, 153, StartCodon { frame: Some(0) }),
                            (150, 210, CDS { frame: Some(0) }),
                            (210, 213, StopCodon { frame: Some(0) }),
                            (210, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_1_from_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((100, 160)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 103, StartCodon { frame: Some(0) }),
                            (100, 160, CDS { frame: Some(0) }),
                            (160, 163, StopCodon { frame: Some(0) }),
                            (160, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_2_to_near_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((150, 297)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 150, UTR5), (150, 153, StartCodon { frame: Some(0) }),
                            (150, 297, CDS { frame: Some(0) }),
                            (297, 300, StopCodon { frame: Some(0) }), (297, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_2_from_trx5end_to_near_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 301), (400, 500), (700, 1000)],
                             Some((100, 298)));
    assert_eq!(exon_coords(&trx), vec![(100, 301), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 103, StartCodon { frame: Some(0) }),
                            (100, 298, CDS { frame: Some(0) }),
                            (298, 301, StopCodon { frame: Some(0) }), (298, 301, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_2_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((150, 300)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 150, UTR5), (150, 153, StartCodon { frame: Some(0) }),
                            (150, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 403, StopCodon { frame: Some(0) }), (400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_2_from_trx5end_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 280), (400, 500), (700, 1000)],
                             Some((100, 280)));
    assert_eq!(exon_coords(&trx), vec![(100, 280), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 103, StartCodon { frame: Some(0) }),
                            (100, 280, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 403, StopCodon { frame: Some(0) }), (400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_2_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((149, 299)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 149, UTR5), (149, 152, StartCodon { frame: Some(0) }),
                            (149, 299, CDS { frame: Some(0) }),
                            (299, 300, StopCodon { frame: Some(0) }), (299, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 402, StopCodon { frame: Some(2) }), (400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_2_from_trx5end_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 280), (400, 500), (700, 1000)],
                             Some((100, 279)));
    assert_eq!(exon_coords(&trx), vec![(100, 280), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 103, StartCodon { frame: Some(0) }),
                            (100, 279, CDS { frame: Some(0) }),
                            (279, 280, StopCodon { frame: Some(0) }), (279, 280, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 402, StopCodon { frame: Some(2) }), (400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_2_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((150, 400)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 150, UTR5), (150, 153, StartCodon { frame: Some(0) }),
                            (150, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 403, StopCodon { frame: Some(0) }), (400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_2_from_trx5end_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 280), (400, 500), (700, 1000)],
                             Some((100, 400)));
    assert_eq!(exon_coords(&trx), vec![(100, 280), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 103, StartCodon { frame: Some(0) }),
                            (100, 280, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 403, StopCodon { frame: Some(0) }), (400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_3() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((150, 460)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 150, UTR5), (150, 153, StartCodon { frame: Some(0) }),
                            (150, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 460, CDS { frame: Some(0) }),
                            (460, 463, StopCodon { frame: Some(0) }), (460, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_3_from_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((100, 470)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 103, StartCodon { frame: Some(0) }),
                            (100, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 470, CDS { frame: Some(1) }),
                            (470, 473, StopCodon { frame: Some(0) }), (470, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_4_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((190, 500)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 190, UTR5), (190, 193, StartCodon { frame: Some(0) }),
                            (190, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0)}), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_4_from_trx5end_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((100, 500)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 103, StartCodon { frame: Some(0) }),
                            (100, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0)}), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_4_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((199, 499)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 199, UTR5), (199, 202, StartCodon { frame: Some(0) }),
                            (199, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 499, CDS { frame: Some(1) }),
                            (499, 500, StopCodon { frame: Some(0) }), (499, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 702, StopCodon { frame: Some(2) }), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_4_from_trx5end_to_split() {
    let (trx, fxs) = trx_fxs(101, 1000, Forward, vec![(101, 300), (400, 500), (700, 1000)],
                             Some((101, 499)));
    assert_eq!(exon_coords(&trx), vec![(101, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(101, 104, StartCodon { frame: Some(0) }),
                            (101, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 499, CDS { frame: Some(2) }),
                            (499, 500, StopCodon { frame: Some(0) }), (499, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 702, StopCodon { frame: Some(2) }), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_4_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((190, 700)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 190, UTR5), (190, 193, StartCodon { frame: Some(0) }),
                            (190, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0) }), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_4_from_trx5end_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((100, 700)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 103, StartCodon { frame: Some(0) }),
                            (100, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0)}), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_5() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((200, 800)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 200, UTR5), (200, 203, StartCodon { frame: Some(0) }),
                            (200, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(2) })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: Some(1) }),
                            (800, 803, StopCodon { frame: Some(0) }), (800, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_1_5_from_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((100, 760)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 103, StartCodon { frame: Some(0) }),
                            (100, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 760, CDS { frame: Some(0) }),
                            (760, 763, StopCodon { frame: Some(0) }), (760, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_4_from_exon3end_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((300, 520)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 403, StartCodon { frame: Some(0) }),
                            (400, 520, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0) }),
                            (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_4_from_exon3end_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 298), (400, 520), (700, 1000)],
                             Some((298, 518)));
    assert_eq!(exon_coords(&trx), vec![(100, 298), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 298, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 403, StartCodon { frame: Some(0) }),
                            (400, 518, CDS { frame: Some(0) }),
                            (518, 520, StopCodon { frame: Some(0) }), (518, 520, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 701, StopCodon { frame: Some(1) }), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_4_from_exon3end_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((300, 700)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 403, StartCodon { frame: Some(0) }),
                            (400, 520, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0) }),
                            (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_4_from_split_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 518), (700, 1000)],
                             Some((298, 518)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 518), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 298, UTR5), (298, 300, StartCodon { frame: Some(0) }),
                            (298, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 401, StartCodon { frame: Some(1) }),
                            (400, 518, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0) }), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_4_from_split_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((298, 518)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 298, UTR5), (298, 300, StartCodon { frame: Some(0) }),
                            (298, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 401, StartCodon { frame: Some(1) }),
                            (400, 518, CDS { frame: Some(1) }),
                            (518, 520, StopCodon { frame: Some(0) }), (518, 520, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 701, StopCodon { frame: Some(1) }), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_4_from_split_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 518), (700, 1000)],
                             Some((298, 700)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 518), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 298, UTR5), (298, 300, StartCodon { frame: Some(0) }),
                            (298, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 401, StartCodon { frame: Some(1) }),
                            (400, 518, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0) }), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_4_from_exon5end_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((400, 520)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 403, StartCodon { frame: Some(0) }),
                            (400, 520, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0) }),
                            (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_4_from_exon5end_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((400, 499)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 403, StartCodon { frame: Some(0) }),
                            (400, 499, CDS { frame: Some(0) }),
                            (499, 500, StopCodon { frame: Some(0) }), (499, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 702, StopCodon { frame: Some(2) }), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_4_from_exon5end_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((400, 700)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 403, StartCodon { frame: Some(0) }),
                            (400, 520, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0) }),
                            (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_5_from_near_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((297, 800)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 297, UTR5), (297, 300, StartCodon { frame: Some(0) }),
                            (297, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: Some(2) }),
                            (800, 803, StopCodon { frame: Some(0) }), (800, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_5_from_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((300, 900)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 403, StartCodon { frame: Some(0) }),
                            (400, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 900, CDS { frame: Some(2) }),
                            (900, 903, StopCodon { frame: Some(0) }), (900, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_5_from_split_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((298, 901)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 298, UTR5), (298, 300, StartCodon { frame: Some(0) }),
                            (298, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 401, StartCodon { frame: Some(1) }),
                            (400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 901, CDS { frame: Some(0) }),
                            (901, 904, StopCodon { frame: Some(0) }), (901, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_2_5_from_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((400, 900)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 403, StartCodon { frame: Some(0) }),
                            (400, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 900, CDS { frame: Some(2) }),
                            (900, 903, StopCodon { frame: Some(0) }), (900, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_3_3() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((450, 480)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 450, UTR5), (450, 453, StartCodon { frame: Some(0) }),
                            (450, 480, CDS { frame: Some(0) }),
                            (480, 483, StopCodon { frame: Some(0) }), (480, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_3_4_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((440, 500)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 440, UTR5), (440, 443, StartCodon { frame: Some(0) }),
                            (440, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0) }),
                            (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_3_4_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((447, 499)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 447, UTR5), (447, 450, StartCodon { frame: Some(0) }),
                            (447, 499, CDS { frame: Some(0) }),
                            (499, 500, StopCodon { frame: Some(0) }), (499, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 702, StopCodon { frame: Some(2) }), (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_3_4_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((440, 700)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 440, UTR5), (440, 443, StartCodon { frame: Some(0) }),
                            (440, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 703, StopCodon { frame: Some(0) }),
                            (700, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_3_5() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((450, 800)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 450, UTR5), (450, 453, StartCodon { frame: Some(0) }),
                            (450, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: Some(1) }),
                            (800, 803, StopCodon { frame: Some(0) }), (800, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_4_5_from_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((500, 950)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 703, StartCodon { frame: Some(0) }),
                            (700, 950, CDS { frame: Some(0) }),
                            (950, 953, StopCodon { frame: Some(0) }), (950, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_4_5_from_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((499, 939)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 499, UTR5), (499, 500, StartCodon { frame: Some(0) }),
                            (499, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 702, StartCodon { frame: Some(2) }),
                            (700, 939, CDS { frame: Some(2) }),
                            (939, 942, StopCodon { frame: Some(0) }), (939, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_4_5_from_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((700, 950)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 703, StartCodon { frame: Some(0) }),
                            (700, 950, CDS { frame: Some(0) }),
                            (950, 953, StopCodon { frame: Some(0) }), (950, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_5_5() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((800, 950)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 800, UTR5), (800, 803, StartCodon { frame: Some(0) }),
                            (800, 950, CDS { frame: Some(0) }),
                            (950, 953, StopCodon { frame: Some(0) }), (950, 1000, UTR3)]);
}

#[test]
fn tbuilder_coords_fwd_coding_incl_stop() {
    let btrxe = TBuilder::new("chrT", 100, 1000)
        .strand(Forward)
        .coords(vec![(100, 400), (700, 1000)], Some((100, 1000)))
        .build();
    assert!(btrxe.is_err());

    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Forward)
        .coords(vec![(100, 400), (700, 1000)], Some((100, 1000)))
        .coding_incl_stop(true)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    let fxs = exon_fxs_coords(&trx);
    assert_eq!(fxs.len(), 2);
    assert_eq!(fxs[0], vec![(100, 103, StartCodon { frame: Some(0) }),
                            (100, 400, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(700, 997, CDS { frame: Some(0) }),
                            (997, 1000, StopCodon { frame: Some(0) }), (997, 1000, UTR3)]);
}

// Reverse strand cases

#[test]
fn tbuilder_coords_rev() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)], None);
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert!(fxs.iter().all(|fx| fx.len() == 0));
}

#[test]
fn tbuilder_coords_rev_coding_1_1() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((150, 210)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 150, UTR3), (147, 150, StopCodon { frame: Some(0) }),
                            (150, 210, CDS { frame: Some(0) }),
                            (207, 210, StartCodon { frame: Some(0) }), (210, 300, UTR5)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_1_2_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((150, 300)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 150, UTR3), (147, 150, StopCodon { frame: Some(0) }),
                            (150, 300, CDS { frame: Some(0) }),
                            (297, 300, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_1_2_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((150, 400)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 150, UTR3), (147, 150, StopCodon { frame: Some(0) }),
                            (150, 300, CDS { frame: Some(0) }),
                            (297, 300, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_1_2_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((148, 402)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 148, UTR3), (145, 148, StopCodon { frame: Some(0) }),
                            (148, 300, CDS { frame: Some(1) }),
                            (299, 300, StartCodon { frame: Some(1) })]);
    assert_eq!(fxs[1], vec![(400, 402, CDS { frame: Some(0) }),
                            (400, 402, StartCodon { frame: Some(0) }), (402, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_1_3() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((200, 450)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 200, UTR3), (197, 200, StopCodon { frame: Some(0) }),
                            (200, 300, CDS { frame: Some(1) })]);
    assert_eq!(fxs[1], vec![(400, 450, CDS { frame: Some(0) }),
                            (447, 450, StartCodon { frame: Some(0) }), (450, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_1_4_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((190, 500)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 190, UTR3), (187, 190, StopCodon { frame: Some(0) }),
                            (190, 300, CDS { frame: Some(2) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(0) }),
                            (497, 500, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_1_4_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((199, 702)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 199, UTR3), (196, 199, StopCodon { frame: Some(0) }),
                            (199, 300, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) }),
                            (499, 500, StartCodon { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 702, CDS { frame: Some(0) }),
                            (700, 702, StartCodon { frame: Some(0) }), (702, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_1_4_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((190, 700)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 190, UTR3), (187, 190, StopCodon { frame: Some(0) }),
                            (190, 300, CDS { frame: Some(2) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(0) }),
                            (497, 500, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_1_5() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((200, 800)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 200, UTR3), (197, 200, StopCodon { frame: Some(0) }),
                            (200, 300, CDS { frame: Some(1) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(2) })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: Some(0) }),
                            (797, 800, StartCodon { frame: Some(0) }),
                            (800, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_1_5_to_near_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((200, 703)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 200, UTR3), (197, 200, StopCodon { frame: Some(0) }),
                            (200, 300, CDS { frame: Some(2) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 703, CDS { frame: Some(0) }),
                            (700, 703, StartCodon { frame: Some(0) }), (703, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_1_5_to_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((250, 1000)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 250, UTR3), (247, 250, StopCodon { frame: Some(0) }),
                            (250, 300, CDS { frame: Some(2) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, CDS { frame: Some(0) }),
                            (997, 1000, StartCodon { frame: Some(0) })]);
}

#[test]
fn tbuilder_coords_rev_coding_2_3_from_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((300, 460)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 460, CDS { frame: Some(0) }),
                            (457, 460, StartCodon { frame: Some(0) }), (460, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_3_from_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((400, 460)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 460, CDS { frame: Some(0) }),
                            (457, 460, StartCodon { frame: Some(0) }), (460, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_3_from_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((402, 462)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (299, 300, StopCodon { frame: Some(1) })]);
    assert_eq!(fxs[1], vec![(400, 402, UTR3), (400, 402, StopCodon { frame: Some(0) }),
                            (402, 462, CDS { frame: Some(0) }),
                            (459, 462, StartCodon { frame: Some(0) }), (462, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_4_from_exon5end_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((300, 520)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 520, CDS { frame: Some(0) }),
                            (517, 520, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_4_from_exon5end_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((300, 700)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 520, CDS { frame: Some(0) }),
                            (517, 520, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_4_from_exon5end_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((300, 702)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) }),
                            (499, 500, StartCodon { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 702, CDS { frame: Some(0) }),
                            (700, 702, StartCodon { frame: Some(0) }), (702, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_4_from_exon3end_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((400, 520)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 520, CDS { frame: Some(0) }),
                            (517, 520, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_4_from_exon3end_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((400, 700)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 520, CDS { frame: Some(0) }),
                            (517, 520, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_4_from_exon3end_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 519), (700, 1000)],
                             Some((400, 701)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 519), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 519, CDS { frame: Some(2) }),
                            (517, 519, StartCodon { frame: Some(2) })]);
    assert_eq!(fxs[2], vec![(700, 701, CDS { frame: Some(0) }),
                            (700, 701, StartCodon { frame: Some(0) }), (701, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_4_from_split_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 522), (700, 1000)],
                             Some((402, 522)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 522), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (299, 300, StopCodon { frame: Some(1) })]);
    assert_eq!(fxs[1], vec![(400, 402, UTR3), (400, 402, StopCodon { frame: Some(0) }),
                            (402, 522, CDS { frame: Some(0) }),
                            (519, 522, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_4_from_split_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 522), (700, 1000)],
                             Some((402, 700)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 522), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (299, 300, StopCodon { frame: Some(1) })]);
    assert_eq!(fxs[1], vec![(400, 402, UTR3), (400, 402, StopCodon { frame: Some(0) }),
                            (402, 522, CDS { frame: Some(0) }),
                            (519, 522, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_4_from_split_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((402, 702)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (299, 300, StopCodon { frame: Some(1) })]);
    assert_eq!(fxs[1], vec![(400, 402, UTR3), (400, 402, StopCodon { frame: Some(0) }),
                            (402, 520, CDS { frame: Some(1) }),
                            (519, 520, StartCodon { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 702, CDS { frame: Some(0) }),
                            (700, 702, StartCodon { frame: Some(0) }), (702, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_5_from_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((300, 900)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 900, CDS { frame: Some(0) }),
                            (897, 900, StartCodon { frame: Some(0) }), (900, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_5_from_exon5end_to_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((300, 1000)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 520, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, CDS { frame: Some(0) }),
                            (997, 1000, StartCodon { frame: Some(0) })]);
}

#[test]
fn tbuilder_coords_rev_coding_2_5_from_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((401, 901)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (298, 300, StopCodon { frame: Some(2) })]);
    assert_eq!(fxs[1], vec![(400, 401, UTR3), (400, 401, StopCodon { frame: Some(0) }),
                            (401, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 901, CDS { frame: Some(0) }),
                            (898, 901, StartCodon { frame: Some(0) }), (901, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_5_from_split_to_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 552), (700, 1000)],
                             Some((402, 1000)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 552), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (299, 300, StopCodon { frame: Some(1) })]);
    assert_eq!(fxs[1], vec![(400, 402, UTR3), (400, 402, StopCodon { frame: Some(0) }),
                            (402, 552, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, CDS { frame: Some(0) }),
                            (997, 1000, StartCodon { frame: Some(0) })]);
}

#[test]
fn tbuilder_coords_rev_coding_2_5_from_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((400, 900)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: Some(1) })]);
    assert_eq!(fxs[2], vec![(700, 900, CDS { frame: Some(0) }),
                            (897, 900, StartCodon { frame: Some(0) }), (900, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_2_5_from_exon3end_to_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 520), (700, 1000)],
                             Some((400, 1000)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 520), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3), (297, 300, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(400, 520, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, CDS { frame: Some(0) }),
                            (997, 1000, StartCodon { frame: Some(0) })]);
}

#[test]
fn tbuilder_coords_rev_coding_3_3() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((420, 480)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 420, UTR3), (417, 420, StopCodon { frame: Some(0) }),
                            (420, 480, CDS { frame: Some(0) }),
                            (477, 480, StartCodon { frame: Some(0) }), (480, 500, UTR5)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_3_4_to_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((440, 500)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 440, UTR3), (437, 440, StopCodon { frame: Some(0) }),
                            (440, 500, CDS { frame: Some(0) }),
                            (497, 500, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_3_4_to_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((420, 701)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 420, UTR3), (417, 420, StopCodon { frame: Some(0) }),
                            (420, 500, CDS { frame: Some(2) }),
                            (498, 500, StartCodon { frame: Some(2) })]);
    assert_eq!(fxs[2], vec![(700, 701, CDS { frame: Some(0) }),
                            (700, 701, StartCodon { frame: Some(0) }), (701, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_3_4_to_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((440, 700)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 440, UTR3), (437, 440, StopCodon { frame: Some(0) }),
                            (440, 500, CDS { frame: Some(0) }),
                            (497, 500, StartCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_3_5() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((450, 800)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 450, UTR3), (447, 450, StopCodon { frame: Some(0) }),
                            (450, 500, CDS { frame: Some(2) })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: Some(0) }),
                            (797, 800, StartCodon { frame: Some(0) }), (800, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_3_5_to_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((460, 1000)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 460, UTR3), (457, 460, StopCodon { frame: Some(0) }),
                            (460, 500, CDS { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, CDS { frame: Some(0) }),
                            (997, 1000, StartCodon { frame: Some(0) })]);
}

#[test]
fn tbuilder_coords_rev_coding_4_5_from_exon5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((500, 850)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3), (497, 500, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 850, CDS { frame: Some(0) }),
                            (847, 850, StartCodon { frame: Some(0) }), (850, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_4_5_from_exon5end_to_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((500, 1000)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3), (497, 500, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, CDS { frame: Some(0) }),
                            (997, 1000, StartCodon { frame: Some(0) })]);
}

#[test]
fn tbuilder_coords_rev_coding_4_5_from_exon3end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((700, 850)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3), (497, 500, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 850, CDS { frame: Some(0) }),
                            (847, 850, StartCodon { frame: Some(0) }), (850, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_4_5_from_exon3end_to_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((700, 1000)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3), (497, 500, StopCodon { frame: Some(0) })]);
    assert_eq!(fxs[2], vec![(700, 1000, CDS { frame: Some(0) }),
                            (997, 1000, StartCodon { frame: Some(0) })]);
}

#[test]
fn tbuilder_coords_rev_coding_4_5_from_split() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((701, 851)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3), (498, 500, StopCodon { frame: Some(2) })]);
    assert_eq!(fxs[2], vec![(700, 701, UTR3), (700, 701, StopCodon { frame: Some(0) }),
                            (701, 851, CDS { frame: Some(0) }),
                            (848, 851, StartCodon { frame: Some(0) }), (851, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_4_5_from_split_to_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1001, Reverse, vec![(100, 300), (400, 500), (700, 1001)],
                             Some((701, 1001)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1001)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3), (498, 500, StopCodon { frame: Some(2) })]);
    assert_eq!(fxs[2], vec![(700, 701, UTR3), (700, 701, StopCodon { frame: Some(0) }),
                            (701, 1001, CDS { frame: Some(0) }),
                            (998, 1001, StartCodon { frame: Some(0) })]);
}

#[test]
fn tbuilder_coords_rev_coding_5_5() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((720, 870)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 720, UTR3), (717, 720, StopCodon { frame: Some(0) }),
                            (720, 870, CDS { frame: Some(0) }),
                            (867, 870, StartCodon { frame: Some(0) }), (870, 1000, UTR5)]);
}

#[test]
fn tbuilder_coords_rev_coding_5_5_to_trx5end() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((850, 1000)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR3)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR3)]);
    assert_eq!(fxs[2], vec![(700, 850, UTR3), (847, 850, StopCodon { frame: Some(0) }),
                            (850, 1000, CDS { frame: Some(0) }),
                            (997, 1000, StartCodon { frame: Some(0) })]);
}

#[test]
fn tbuilder_coords_rev_coding_incl_stop() {
    let btrxe = TBuilder::new("chrT", 100, 1000)
        .strand(Reverse)
        .coords(vec![(100, 400), (700, 1000)], Some((100, 1000)))
        .build();
    assert!(btrxe.is_err());

    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Reverse)
        .coords(vec![(100, 400), (700, 1000)], Some((100, 1000)))
        .coding_incl_stop(true)
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    let fxs = exon_fxs_coords(&trx);
    assert_eq!(fxs.len(), 2);
    assert_eq!(fxs[0], vec![(100, 103, UTR3), (100, 103, StopCodon { frame: Some(0) }),
                            (103, 400, CDS { frame: Some(0) })]);
    assert_eq!(fxs[1], vec![(700, 1000, CDS { frame: Some(0) }),
                            (997, 1000, StartCodon { frame: Some(0) })]);
}

// Unknown strand cases

#[test]
fn tbuilder_coords_unk() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)], None);
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert!(fxs.iter().all(|fx| fx.len() == 0));
}

#[test]
fn tbuilder_coords_unk_coding_1_1() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((150, 210)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 150, UTR), (150, 210, CDS { frame: None }), (210, 300, UTR)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_1_2() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((150, 300)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 150, UTR), (150, 300, CDS { frame: None })]);
    assert_eq!(fxs[1], vec![(400, 500, UTR)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_1_3() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((180, 430)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 180, UTR), (180, 300, CDS { frame: None })]);
    assert_eq!(fxs[1], vec![(400, 430, CDS { frame: None} ), (430, 500, UTR)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_1_4() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((250, 500)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 250, UTR), (250, 300, CDS { frame: None })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: None })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_1_5() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((200, 800)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 200, UTR), (200, 300, CDS { frame: None })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: None })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: None }), (800, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_2_3() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((300, 460)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR)]);
    assert_eq!(fxs[1], vec![(400, 460, CDS { frame: None }), (460, 500, UTR)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_2_4() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((300, 500)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR)]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: None })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_2_5() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((300, 800)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR)]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: None })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: None }), (800, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_3_3() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((420, 480)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR)]);
    assert_eq!(fxs[1], vec![(400, 420, UTR), (420, 480, CDS { frame: None }), (480, 500, UTR)]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_3_4() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((440, 500)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR)]);
    assert_eq!(fxs[1], vec![(400, 440, UTR), (440, 500, CDS { frame: None })]);
    assert_eq!(fxs[2], vec![(700, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_3_5() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((450, 800)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR)]);
    assert_eq!(fxs[1], vec![(400, 450, UTR), (450, 500, CDS { frame: None })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: None }), (800, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_4_5() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((500, 850)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR)]);
    assert_eq!(fxs[2], vec![(700, 850, CDS { frame: None }), (850, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_5_5() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((750, 900)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 300, UTR)]);
    assert_eq!(fxs[1], vec![(400, 500, UTR)]);
    assert_eq!(fxs[2], vec![(700, 750, UTR), (750, 900, CDS { frame: None }), (900, 1000, UTR)]);
}

#[test]
fn tbuilder_coords_unk_coding_incl_stop() {
    let btrxe = TBuilder::new("chrT", 100, 1000)
        .strand(Unknown)
        .coords(vec![(100, 400), (700, 1000)], Some((100, 1000)))
        .build();
    assert!(btrxe.is_err());

    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Unknown)
        .coords(vec![(100, 400), (700, 1000)], Some((100, 1000)))
        .coding_incl_stop(true)
        .build();
    assert!(btrx.is_err());
}
