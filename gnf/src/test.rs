use super::*;
use TranscriptFeatureKind as TFK;

#[test]
fn tfbuilder() {
    let tfm1 = TFBuilder::new("chrT", 10, 20, Exon)
        .strand(Strand::Forward)
        .attribute("name", "ex1")
        .id("ex1.1")
        .build();
    assert!(tfm1.is_ok());
    let tf = tfm1.unwrap();
    assert_eq!(tf.seq_name(), "chrT");
    assert_eq!(tf.kind(), &Exon);
    assert_eq!(tf.strand(), &Strand::Forward);
    assert_eq!(tf.id, Some("ex1.1".to_owned()));
    assert_eq!(tf.attributes.get("name"), Some(&"ex1".to_owned()));
    assert_eq!(tf.attributes.len(), 1);

    let tfm2 = TFBuilder::new("chrO", 10, 10, Exon)
        .strand_char('-')
        .strand(Strand::Reverse)
        .build();
    assert!(tfm2.is_ok());
}

#[test]
fn tfbuilder_interval_invalid() {
    let tfm = TFBuilder::new("chrE", 20, 10, Exon)
        .build();
    assert!(tfm.is_err());
    assert!(matches!(tfm.unwrap_err(),
                     FeatureError::IntervalError(utils::IntervalError::InvalidRange)));
}

#[test]
fn tfbuilder_strand_unspecified() {
    let tfm = TFBuilder::new("chrT", 20, 30, Exon)
        .build();
    assert!(tfm.is_err());
    assert!(matches!(tfm.unwrap_err(), FeatureError::UnspecifiedStrandError));
}

#[test]
fn tfbuilder_strand_char_unexpected() {
    let tfm = TFBuilder::new("chrE", 10, 20, Exon)
        .strand_char('w')
        .build();
    assert!(tfm.is_err());
    assert!(matches!(tfm.unwrap_err(),
                     FeatureError::StrandCharError(utils::StrandError::InvalidChar(_))));
}

#[test]
fn tfbuilder_strand_char_conflicting() {
    let tfm = TFBuilder::new("chrE", 10, 20, Exon)
        .strand_char('-')
        .strand(Strand::Reverse)
        .build();
    assert!(tfm.is_ok());
    let tf = tfm.unwrap();
    assert_eq!(tf.strand(), &Strand::Reverse);
}

fn get_trfk_coords(fxs: &Vec<TranscriptFeature>, kind: TFK) -> Vec<(u64, u64)> {
    fxs.iter()
        .filter(|fx| *fx.kind() == kind)
        .map(|fx| (fx.interval().start, fx.interval().end))
        .collect()
}

fn get_trfks(fxs: &Vec<TranscriptFeature>) -> Vec<TFK> {
    fxs.iter().map(|fx| fx.kind().clone()).collect()
}

fn get_framed_trfk(fxs: &Vec<TranscriptFeature>) -> Vec<(TFK, u64)> {
    fxs.iter().filter_map(|fx| {
        match (fx.frame(), fx.kind()) {
            (Some(frame), fxk @ &StartCodon) => Some((fxk.clone(), frame)),
            (Some(frame), fxk @ &CDS) => Some((fxk.clone(), frame)),
            _ => None,
        }
    }).collect()
}

#[test]
fn tbuilder_coords_fwd() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], None)
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs), vec![Exon, Exon, Exon]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
}

#[test]
fn tbuilder_coords_fwd_coding() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 800)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR5, StartCodon, CDS, Exon, CDS, Exon, CDS, StopCodon, UTR3]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(100, 200)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(200, 203)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(200, 300), (400, 500), (700, 800)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(800, 803)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(800, 1000)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(StartCodon, 0), (CDS, 0), (CDS, 2), (CDS, 1)]);
}

#[test]
fn tbuilder_coords_fwd_coding_from_exon5end() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((400, 900)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR5, Exon, StartCodon, CDS, Exon, CDS, StopCodon, UTR3]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(100, 300)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(400, 403)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(400, 500), (700, 900)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(900, 903)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(900, 1000)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(StartCodon, 0), (CDS, 0), (CDS, 2)]);
}

#[test]
fn tbuilder_coords_fwd_coding_from_exon3end() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((300, 900)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR5, Exon, StartCodon, CDS, Exon, CDS, StopCodon, UTR3]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(100, 300)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(400, 403)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(400, 500), (700, 900)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(900, 903)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(900, 1000)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(StartCodon, 0), (CDS, 0), (CDS, 2)]);
}

#[test]
fn tbuilder_coords_fwd_coding_from_near_exon3end() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((297, 800)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR5, StartCodon, CDS, Exon, CDS, Exon, CDS, StopCodon, UTR3]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(100, 297)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(297, 300)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(297, 300), (400, 500), (700, 800)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(800, 803)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(800, 1000)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(StartCodon, 0), (CDS, 0), (CDS, 0), (CDS, 2)]);
}

#[test]
fn tbuilder_coords_fwd_coding_from_split() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((298, 901)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR5, StartCodon, CDS, Exon, StartCodon, CDS, Exon, CDS, StopCodon,
                    UTR3]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(100, 298)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(298, 300), (400, 401)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(298, 300), (400, 500), (700, 901)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(901, 904)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(901, 1000)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(StartCodon, 0), (CDS, 0),
                    (StartCodon, 1), (CDS, 1), (CDS, 0)]);
}

#[test]
fn tbuilder_coords_fwd_coding_to_exon3end() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((190, 500)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR5, StartCodon, CDS, Exon, CDS, Exon, StopCodon, UTR3]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(100, 190)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(190, 193)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(190, 300), (400, 500)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(700, 703)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(700, 1000)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(StartCodon, 0), (CDS, 0), (CDS, 1)]);
}

#[test]
fn tbuilder_coords_fwd_coding_to_exon5end() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((190, 700)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR5, StartCodon, CDS, Exon, CDS, Exon, StopCodon, UTR3]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(100, 190)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(190, 193)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(190, 300), (400, 500)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(700, 703)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(700, 1000)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(StartCodon, 0), (CDS, 0), (CDS, 1)]);
}

#[test]
fn tbuilder_coords_fwd_coding_to_split() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Forward)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((199, 499)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR5, StartCodon, CDS, Exon, CDS, StopCodon, UTR3,
                    Exon, StopCodon, UTR3]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(100, 199)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(199, 202)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(199, 300), (400, 499)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(499, 500), (700, 702)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(499, 500), (700, 1000)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(StartCodon, 0), (CDS, 0), (CDS, 1)]);
}

#[test]
fn tbuilder_coords_rev() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], None)
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs), vec![Exon, Exon, Exon]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
}

#[test]
fn tbuilder_coords_rev_coding() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 800)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR3, StopCodon, CDS, Exon, CDS, Exon, CDS, StartCodon, UTR5]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(100, 200)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(197, 200)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(200, 300), (400, 500), (700, 800)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(797, 800)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(800, 1000)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(CDS, 1), (CDS, 2), (CDS, 0), (StartCodon, 0)]);
}

#[test]
fn tbuilder_coords_rev_coding_from_exon3end() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((300, 900)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR3, StopCodon, Exon, CDS, Exon, CDS, StartCodon, UTR5]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(900, 1000)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(897, 900)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(400, 500), (700, 900)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(297, 300)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(100, 300)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(CDS, 1), (CDS, 0), (StartCodon, 0)]);
}

#[test]
fn tbuilder_coords_rev_coding_from_exon5end() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((400, 900)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR3, StopCodon, Exon, CDS, Exon, CDS, StartCodon, UTR5]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(900, 1000)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(897, 900)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(400, 500), (700, 900)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(297, 300)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(100, 300)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(CDS, 1), (CDS, 0), (StartCodon, 0)]);
}

#[test]
fn tbuilder_coords_rev_coding_from_split() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((401, 901)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR3, StopCodon, Exon, UTR3, StopCodon, CDS, Exon, CDS,
                    StartCodon, UTR5]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(901, 1000)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(898, 901)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(401, 500), (700, 901)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(298, 300), (400, 401)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(100, 300), (400, 401)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(CDS, 0), (CDS, 0), (StartCodon, 0)]);
}

#[test]
fn tbuilder_coords_rev_coding_to_exon5end() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((190, 700)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR3, StopCodon, CDS, Exon, CDS, StartCodon, Exon, UTR5]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(497, 500)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(190, 300), (400, 500)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(187, 190)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(100, 190)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(CDS, 2), (CDS, 0), (StartCodon, 0)]);
}

#[test]
fn tbuilder_coords_rev_coding_to_near_exon5end() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 703)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR3, StopCodon, CDS, Exon, CDS, Exon, CDS, StartCodon, UTR5]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(703, 1000)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(700, 703)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(200, 300), (400, 500), (700, 703)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(197, 200)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(100, 200)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(CDS, 2), (CDS, 0), (CDS, 0), (StartCodon, 0)]);
}

#[test]
fn tbuilder_coords_rev_coding_to_exon3end() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((190, 500)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR3, StopCodon, CDS, Exon, CDS, StartCodon, Exon, UTR5]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(497, 500)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(190, 300), (400, 500)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(187, 190)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(100, 190)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(CDS, 2), (CDS, 0), (StartCodon, 0)]);
}

#[test]
fn tbuilder_coords_rev_coding_to_split() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Reverse)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((199, 702)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR3, StopCodon, CDS, Exon, CDS, StartCodon,
                    Exon, CDS, StartCodon, UTR5]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, UTR5), vec![(702, 1000)]);
    assert_eq!(get_trfk_coords(fxs, StartCodon), vec![(499, 500), (700, 702)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(199, 300), (400, 500), (700, 702)]);
    assert_eq!(get_trfk_coords(fxs, StopCodon), vec![(196, 199)]);
    assert_eq!(get_trfk_coords(fxs, UTR3), vec![(100, 199)]);
    assert_eq!(get_framed_trfk(fxs),
               vec![(CDS, 0), (CDS, 1), (StartCodon, 1), (CDS, 0), (StartCodon, 0)]);
}

#[test]
fn tbuilder_coords_coding_unk() {
    let tm = TBuilder::new("chrT", 100, 1000)
        .strand(Strand::Unknown)
        .feature_coords(vec![(100, 300), (400, 500), (700, 1000)], Some((200, 800)))
        .build();
    assert!(tm.is_ok(), "{:?}", tm);
    let t = tm.unwrap();
    let fxs = t.features();
    assert_eq!(get_trfks(fxs),
               vec![Exon, UTR, CDS, Exon, CDS, Exon, CDS, UTR]);
    assert_eq!(get_trfk_coords(fxs, Exon),
               vec![(100, 300), (400, 500), (700, 1000)]);
    assert_eq!(get_trfk_coords(fxs, CDS), vec![(200, 300), (400, 500), (700, 800)]);
    assert_eq!(get_trfk_coords(fxs, UTR), vec![(100, 200), (800, 1000)]);
    assert_eq!(get_framed_trfk(fxs).len(), 0);
}
