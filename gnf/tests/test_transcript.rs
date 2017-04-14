extern crate bio;
extern crate gnf;

use gnf::{ExonFeatureKind, Strand, TBuilder, Transcript};
use ExonFeatureKind::*;
use Strand::*;

fn exon_coords(transcript: &Transcript) -> Vec<(u64, u64)> {
    transcript.exons().iter()
        .map(|exn| (exn.interval().start, exn.interval().end))
        .collect()
}

fn exon_fxs_coords(transcript: &Transcript) -> Vec<Vec<(u64, u64, ExonFeatureKind)>> {
    transcript.exons().iter()
        .map(|exn| {
            exn.features.iter()
                .map(|fx| (fx.interval.start, fx.interval.end, fx.kind.clone()))
                .collect()
        })
        .collect()
}

fn trx_fxs<T>(start: u64, end: u64, strand: Strand, exon_coords: T,
               coding_coord: Option<(u64, u64)>
) -> (Transcript, Vec<Vec<(u64, u64, ExonFeatureKind)>>)
where T: IntoIterator<Item=(u64, u64)>
{
    let mut btrx = TBuilder::new("chrT", start, end)
            .strand(strand)
            .exon_coords(exon_coords);
    if let Some((a, b)) = coding_coord {
        btrx = btrx.coding_coord(a, b);
    }
    let rtrx = btrx.build();
    let trx = rtrx.unwrap();
    let fxs = exon_fxs_coords(&trx);
    (trx, fxs)
}

#[test]
fn tbuilder_basic() {
    let btrx = TBuilder::new("chrT", 100, 1000)
        .strand(Reverse)
        .exon_coords(vec![(100, 300), (400, 500), (700, 1000)])
        .id("transcript-1")
        .attribute("tag", "basic")
        .build();
    assert!(btrx.is_ok(), "{:?}", btrx);
    let trx = btrx.unwrap();
    assert_eq!(trx.span(), 900);
    assert_eq!(trx.seq_name(), "chrT");
    assert_eq!(trx.strand(), &Strand::Reverse);
    assert_eq!(trx.id, Some("transcript-1".to_owned()));
    assert_eq!(trx.attributes.get("tag"), Some(&"basic".to_owned()));
    assert_eq!(trx.attributes.len(), 1);
    assert_eq!(trx.exons().len(), 3);
}

#[test]
fn tbuilder_coords_fwd() {
    let (trx, fxs) = trx_fxs(100, 1000, Forward, vec![(100, 300), (400, 500), (700, 1000)], None);
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert!(fxs.iter().all(|fx| fx.len() == 0));
}

#[test]
fn tbuilder_coords_fwd_coding() {
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
fn tbuilder_coords_fwd_coding_from_exon5end() {
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
fn tbuilder_coords_fwd_coding_from_exon3end() {
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
fn tbuilder_coords_fwd_coding_from_near_exon3end() {
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
fn tbuilder_coords_fwd_coding_from_split() {
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
fn tbuilder_coords_fwd_coding_to_exon3end() {
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
fn tbuilder_coords_fwd_coding_to_exon5end() {
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
fn tbuilder_coords_fwd_coding_to_split() {
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
fn tbuilder_coords_rev() {
    let (trx, fxs) = trx_fxs(100, 1000, Reverse, vec![(100, 300), (400, 500), (700, 1000)], None);
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert!(fxs.iter().all(|fx| fx.len() == 0));
}

#[test]
fn tbuilder_coords_rev_coding() {
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
fn tbuilder_coords_rev_coding_from_exon3end() {
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
fn tbuilder_coords_rev_coding_from_exon5end() {
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
fn tbuilder_coords_rev_coding_from_split() {
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
fn tbuilder_coords_rev_coding_to_exon5end() {
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
fn tbuilder_coords_rev_coding_to_near_exon5end() {
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
fn tbuilder_coords_rev_coding_to_exon3end() {
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
fn tbuilder_coords_rev_coding_to_split() {
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
fn tbuilder_coords_coding_unk() {
    let (trx, fxs) = trx_fxs(100, 1000, Unknown, vec![(100, 300), (400, 500), (700, 1000)],
                             Some((200, 800)));
    assert_eq!(exon_coords(&trx), vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(fxs.len(), 3);
    assert_eq!(fxs[0], vec![(100, 200, UTR), (200, 300, CDS { frame: None })]);
    assert_eq!(fxs[1], vec![(400, 500, CDS { frame: None })]);
    assert_eq!(fxs[2], vec![(700, 800, CDS { frame: None }), (800, 1000, UTR)]);
}
