extern crate bio;
#[macro_use]
extern crate matches;
extern crate gnf;

use bio::utils::{self, Interval, Strand};

use gnf::{EBuilder, ExonFeature, ExonFeatureKind};
use gnf::error::FeatureError;
use ExonFeatureKind::*;

fn make_feat(start: u64, end: u64, kind: ExonFeatureKind) -> ExonFeature {
    ExonFeature {
        interval: Interval::new(start..end).unwrap(),
        kind: kind,
    }
}

#[test]
fn ebuilder_basic() {
    let exonb = EBuilder::new("chrT", 10, 20)
        .strand(Strand::Forward)
        .id("ex1.1")
        .attribute("name", "ex1")
        .feature(make_feat(10, 15, UTR5))
        .build();
    assert!(exonb.is_ok());
    let exon = exonb.unwrap();
    assert_eq!(exon.seq_name(), "chrT");
    assert_eq!(exon.strand(), &Strand::Forward);
    assert_eq!(exon.id, Some("ex1.1".to_owned()));
    assert_eq!(exon.attributes.get("name"), Some(&"ex1".to_owned()));
    assert_eq!(exon.attributes.len(), 1);
    assert_eq!(exon.features, vec![make_feat(10, 15, UTR5)]);
}

#[test]
fn ebuilder_multiple_strand_input() {
    let exonb = EBuilder::new("chrO", 10, 10)
        .strand_char('-')
        .strand(Strand::Reverse)
        .build();
    assert!(exonb.is_ok());
    let exon = exonb.unwrap();
    assert_eq!(exon.seq_name(), "chrO");
    assert_eq!(exon.strand(), &Strand::Reverse);
    assert_eq!(exon.id, None);
    assert_eq!(exon.attributes.len(), 0);
    assert_eq!(exon.features, vec![]);
}

#[test]
fn ebuilder_interval_invalid() {
    let exonb = EBuilder::new("chrE", 20, 10).build();
    assert!(exonb.is_err());
    assert!(matches!(exonb.unwrap_err(),
                     FeatureError::IntervalError(utils::IntervalError::InvalidRange)));
}

#[test]
fn ebuilder_strand_unspecified() {
    let exonb = EBuilder::new("chrT", 20, 30).build();
    assert!(exonb.is_err());
    assert!(matches!(exonb.unwrap_err(), FeatureError::UnspecifiedStrandError));
}

#[test]
fn ebuilder_strand_char_unexpected() {
    let exonb = EBuilder::new("chrE", 10, 20)
        .strand_char('w')
        .build();
    assert!(exonb.is_err());
    assert!(matches!(exonb.unwrap_err(),
                     FeatureError::StrandCharError(utils::StrandError::InvalidChar(_))));
}

#[test]
fn ebuilder_strand_char_conflicting() {
    let exonb = EBuilder::new("chrE", 10, 20)
        .strand_char('-')
        .strand(Strand::Reverse)
        .build();
    assert!(exonb.is_ok());
    let exon = exonb.unwrap();
    assert_eq!(exon.strand(), &Strand::Reverse);
}
