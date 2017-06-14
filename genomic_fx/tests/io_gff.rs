extern crate bio;
extern crate genomic_fx;

use std::io;

use genomic_fx::{GffType, GffReader,
                 GffGenesStream, GffTranscriptsStream,
                 Gene, Transcript, ExonFeatureKind as EFK, Strand};
use Strand::*;


static SINGLE_GENE_GTF: &'static str = include_str!("data/single_gene.gtf");


fn next_gx<R>(rg: &mut GffGenesStream<R>) -> Gene where R: io::Read {
    rg.next().expect("a gene result").expect("a gene")
}

#[test]
fn gtf_reader_single_gene() {
    let mut reader = GffReader::from_reader(SINGLE_GENE_GTF.as_bytes(), GffType::GTF2);
    let mut genes = reader.genes_stream();

    let gx = next_gx(&mut genes);
    assert_eq!(gx.seq_name(), "chr2");
    assert_eq!(gx.id(), Some("ENSG00000128645.13"));
    assert_eq!(gx.start(), 176188578);
    assert_eq!(gx.end(), 176190907);
    assert_eq!(gx.strand(), &Forward);
    assert_eq!(gx.transcripts().len(), 2);

    let trx1 = gx.transcripts().get("ENST00000331462.5").expect("a transcript");
    assert_eq!(trx1.seq_name(), "chr2");
    assert_eq!(trx1.start(), 176188578);
    assert_eq!(trx1.end(), 176190907);
    assert_eq!(trx1.strand(), &Forward);
    assert_eq!(trx1.exons().len(), 2);
    let trx1_exon1 = [
        ((176188578, 176188801), EFK::UTR5),
        ((176188801, 176188804), EFK::StartCodon { frame: Some(0) }),
        ((176188801, 176189453), EFK::CDS { frame: Some(0) })];
    for (eidx, feat) in trx1.exons()[0].features().iter().enumerate() {
        assert_eq!(feat.start(), (trx1_exon1[eidx].0).0);
        assert_eq!(feat.end(), (trx1_exon1[eidx].0).1);
        assert_eq!(feat.kind(), &trx1_exon1[eidx].1)
    }
    let trx1_exon2 = [
        ((176189807, 176190139), EFK::CDS { frame: Some(2) }),
        ((176190139, 176190142), EFK::StopCodon { frame: Some(0) }),
        ((176190142, 176190907), EFK::UTR3)];
    for (eidx, feat) in trx1.exons()[1].features().iter().enumerate() {
        assert_eq!(feat.start(), (trx1_exon2[eidx].0).0);
        assert_eq!(feat.end(), (trx1_exon2[eidx].0).1);
        assert_eq!(feat.kind(), &trx1_exon2[eidx].1)
    }

    let trx2 = gx.transcripts().get("ENST00000610524.1").expect("a transcript");
    assert_eq!(trx2.seq_name(), "chr2");
    assert_eq!(trx2.start(), 176188842);
    assert_eq!(trx2.end(), 176188901);
    assert_eq!(trx2.strand(), &Forward);
    assert_eq!(trx2.exons().len(), 1);
    assert_eq!(trx2.exons()[0].features().len(), 0);

    assert!(genes.next().is_none());
}

fn next_trx<R>(rt: &mut GffTranscriptsStream<R>) -> Transcript where R: io::Read {
    rt.next().expect("a transcript result").expect("a transcript")
}

#[test]
fn gtf_reader_multiple_transcripts_stream() {
    let mut reader = GffReader::from_reader(SINGLE_GENE_GTF.as_bytes(), GffType::GTF2);
    let mut transcripts = reader.transcripts_stream();

    let trx1 = next_trx(&mut transcripts);
    assert_eq!(trx1.seq_name(), "chr2");
    assert_eq!(trx1.start(), 176188578);
    assert_eq!(trx1.end(), 176190907);
    assert_eq!(trx1.strand(), &Forward);
    assert_eq!(trx1.exons().len(), 2);
    let trx1_exon1 = [
        ((176188578, 176188801), EFK::UTR5),
        ((176188801, 176188804), EFK::StartCodon { frame: Some(0) }),
        ((176188801, 176189453), EFK::CDS { frame: Some(0) })];
    for (eidx, feat) in trx1.exons()[0].features().iter().enumerate() {
        assert_eq!(feat.start(), (trx1_exon1[eidx].0).0);
        assert_eq!(feat.end(), (trx1_exon1[eidx].0).1);
        assert_eq!(feat.kind(), &trx1_exon1[eidx].1)
    }
    let trx1_exon2 = [
        ((176189807, 176190139), EFK::CDS { frame: Some(2) }),
        ((176190139, 176190142), EFK::StopCodon { frame: Some(0) }),
        ((176190142, 176190907), EFK::UTR3)];
    for (eidx, feat) in trx1.exons()[1].features().iter().enumerate() {
        assert_eq!(feat.start(), (trx1_exon2[eidx].0).0);
        assert_eq!(feat.end(), (trx1_exon2[eidx].0).1);
        assert_eq!(feat.kind(), &trx1_exon2[eidx].1)
    }

    let trx2 = next_trx(&mut transcripts);
    assert_eq!(trx2.seq_name(), "chr2");
    assert_eq!(trx2.start(), 176188842);
    assert_eq!(trx2.end(), 176188901);
    assert_eq!(trx2.strand(), &Forward);
    assert_eq!(trx2.exons().len(), 1);
    assert_eq!(trx2.exons()[0].features().len(), 0);

    assert!(transcripts.next().is_none());
}
