extern crate bio;
extern crate linked_hash_map;
extern crate genomic_fx;

use std::io;

use linked_hash_map::LinkedHashMap;

use genomic_fx::{RefFlatReader, RefFlatWriter, RefFlatRecord,
                 RefFlatRecords, RefFlatTranscripts, RefFlatGenes,
                 Transcript, TBuilder, Gene, GBuilder, Strand};


static SINGLE_ROW_NO_CDS: &'static str = include_str!("data/single_row_no_cds.refFlat");
static MULT_ROWS_NO_CDS: &'static str = include_str!("data/mult_rows_no_cds.refFlat");
static MULT_ROWS_MULT_GENES_WITH_CDS: &'static str =
    include_str!("data/mult_rows_mult_genes_with_cds.refFlat");


fn next_rec<'a, R>(rr: &mut RefFlatRecords<'a, R>) -> RefFlatRecord where R: io::Read {
    rr.next().expect("a refflat record result").expect("a refflat record")
}

fn next_trx<'a, R>(rt: &mut RefFlatTranscripts<'a, R>) -> Transcript where R: io::Read {
    rt.next().expect("a transcript result").expect("a transcript")
}

fn next_gx<'a, R>(rg: &mut RefFlatGenes<'a, R>) -> Gene where R: io::Read {
    rg.next().expect("a gene result").expect("a gene")
}

#[test]
fn refflat_reader_records_single_row_no_cds() {
    let mut reader = RefFlatReader::from_reader(SINGLE_ROW_NO_CDS.as_bytes());
    let mut records = reader.records();

    let rec1 = next_rec(&mut records);
    assert_eq!(rec1.gene_id, "DDX11L1".to_owned());

    assert!(records.next().is_none());
}

#[test]
fn refflat_reader_transcripts_single_row_no_cds() {
    let mut reader = RefFlatReader::from_reader(SINGLE_ROW_NO_CDS.as_bytes());
    let mut transcripts = reader.transcripts();

    let trx1 = next_trx(&mut transcripts);
    assert_eq!(trx1.id(), Some("NR_046018"));

    assert!(transcripts.next().is_none());
}

#[test]
fn refflat_reader_genes_single_row_no_cds() {
    let mut reader = RefFlatReader::from_reader(SINGLE_ROW_NO_CDS.as_bytes());
    let mut genes = reader.genes();

    let gx1 = next_gx(&mut genes);
    assert_eq!(gx1.id(), Some("DDX11L1"));

    assert!(genes.next().is_none());
}

#[test]
fn refflat_reader_records_mult_rows_no_cds() {
    let mut reader = RefFlatReader::from_reader(MULT_ROWS_NO_CDS.as_bytes());
    let mut records = reader.records();

    let rec1 = next_rec(&mut records);
    assert_eq!(rec1.gene_id, "DDX11L1".to_owned());

    let rec2 = next_rec(&mut records);
    assert_eq!(rec2.gene_id, "MIR570".to_owned());

    assert!(records.next().is_none());
}

#[test]
fn refflat_reader_transcripts_mult_rows_no_cds() {
    let mut reader = RefFlatReader::from_reader(MULT_ROWS_NO_CDS.as_bytes());
    let mut transcripts = reader.transcripts();

    let trx1 = next_trx(&mut transcripts);
    assert_eq!(trx1.id(), Some("NR_046018"));

    let trx2 = next_trx(&mut transcripts);
    assert_eq!(trx2.id(), Some("NR_030296"));

    assert!(transcripts.next().is_none());
}

#[test]
fn refflat_reader_genes_mult_rows_no_cds() {
    let mut reader = RefFlatReader::from_reader(MULT_ROWS_NO_CDS.as_bytes());
    let mut genes = reader.genes();

    let gx1 = next_gx(&mut genes);
    assert_eq!(gx1.id(), Some("DDX11L1"));

    let gx2 = next_gx(&mut genes);
    assert_eq!(gx2.id(), Some("MIR570"));

    assert!(genes.next().is_none());
}

#[test]
fn refflat_reader_records_mult_rows_mult_genes_with_cds() {
    let mut reader = RefFlatReader::from_reader(MULT_ROWS_MULT_GENES_WITH_CDS.as_bytes());
    let mut records = reader.records();

    let rec1 = next_rec(&mut records);
    assert_eq!(rec1.transcript_name, "NM_001297605".to_owned());

    let _rec2 = next_rec(&mut records);
    let _rec3 = next_rec(&mut records);
    let _rec4 = next_rec(&mut records);

    let rec5 = next_rec(&mut records);
    assert_eq!(rec5.transcript_name, "NM_138428".to_owned());

    assert!(records.next().is_none());
}

#[test]
fn refflat_reader_transcripts_mult_rows_mult_genes_with_cds() {
    let mut reader = RefFlatReader::from_reader(MULT_ROWS_MULT_GENES_WITH_CDS.as_bytes());
    let mut transcripts = reader.transcripts();

    let trx1 = next_trx(&mut transcripts);
    assert_eq!(trx1.id(), Some("NM_001297605"));

    let _trx2 = next_trx(&mut transcripts);
    let _trx3 = next_trx(&mut transcripts);
    let _trx4 = next_trx(&mut transcripts);

    let trx5 = next_trx(&mut transcripts);
    assert_eq!(trx5.id(), Some("NM_138428"));

    assert!(transcripts.next().is_none());
}

#[test]
fn refflat_reader_genes_mult_rows_mult_genes_with_cds() {
    let mut reader = RefFlatReader::from_reader(MULT_ROWS_MULT_GENES_WITH_CDS.as_bytes());
    let mut genes = reader.genes();

    let gx1 = next_gx(&mut genes);
    assert_eq!(gx1.id(), Some("TNFRSF14"));

    let gx2 = next_gx(&mut genes);
    assert_eq!(gx2.id(), Some("SMIM12"));

    assert!(genes.next().is_none());
}

#[test]
fn refflat_writer_rows_single_row_no_cds() {
    let row =
        ("DDX11L1".to_owned(), "NR_046018".to_owned(), "chr1".to_owned(),
         '+', 11873, 14409, 14409, 14409, 3,
         "11873,12612,13220,".to_owned(), "12227,12721,14409,".to_owned());

    let mut writer = RefFlatWriter::from_memory();
    writer.write(&row).expect("a successful write");
    assert_eq!(writer.as_string(), SINGLE_ROW_NO_CDS);
}

#[test]
fn refflat_writer_records_single_row_no_cds() {
    let rec = RefFlatRecord {
        gene_id: "DDX11L1".to_owned(),
        transcript_name: "NR_046018".to_owned(),
        seq_name: "chr1".to_owned(),
        strand: '+',
        trx_start: 11873,
        trx_end: 14409,
        coding_start: 14409,
        coding_end: 14409,
        num_exons: 3,
        exon_starts: "11873,12612,13220,".to_owned(),
        exon_ends: "12227,12721,14409,".to_owned()
    };

    let mut writer = RefFlatWriter::from_memory();
    writer.write_record(&rec).expect("a successful write");
    assert_eq!(writer.as_string(), SINGLE_ROW_NO_CDS);
}

#[test]
fn refflat_writer_transcripts_single_row_no_cds() {
    let trx = TBuilder::new("chr1", 11873,  14409)
        .strand_char('+')
        .id("NR_046018")
        .gene_id("DDX11L1")
        .coords(vec![(11873, 12227), (12612, 12721), (13220, 14409)], None)
        .coding_incl_stop(true)
        .build()
        .expect("a transcript");

    let mut writer = RefFlatWriter::from_memory();
    writer.write_transcript(&trx).expect("a successful write");
    assert_eq!(writer.as_string(), SINGLE_ROW_NO_CDS);
}

#[test]
fn refflat_writer_genes_single_row_no_cds() {
    let mut cs = LinkedHashMap::new();
    cs.insert("NR_046018".to_owned(),
              ((11873, 14409), vec![(11873, 12227), (12612, 12721), (13220, 14409)], None));

    let gx = GBuilder::new("chr1", 11873, 14409)
        .strand_char('+')
        .id("DDX11L1")
        .transcript_coords(cs)
        .transcript_coding_incl_stop(true)
        .build()
        .expect("a gene");

    let mut writer = RefFlatWriter::from_writer(vec![]);
    let res = writer.write_gene(&gx);
    assert!(res.is_ok(), "{:?}");
    assert_eq!(writer.as_string(), SINGLE_ROW_NO_CDS);
}

#[test]
fn refflat_writer_records_mult_rows_mult_genes_with_cds() {
    let recs = [
        RefFlatRecord {
            gene_id: "TNFRSF14".to_owned(),
            transcript_name: "NM_001297605".to_owned(),
            seq_name: "chr1".to_owned(),
            strand: '+',
            trx_start: 2556364,
            trx_end: 2565622,
            coding_start: 2556664,
            coding_end: 2562868,
            num_exons: 7,
            exon_starts: "2556364,2557725,2558342,2559822,2560623,2562864,2563147,".to_owned(),
            exon_ends: "2556733,2557834,2558468,2559978,2560714,2562896,2565622,".to_owned(),
        },
        RefFlatRecord {
            gene_id: "TNFRSF14".to_owned(),
            transcript_name: "NM_003820".to_owned(),
            seq_name: "chr1".to_owned(),
            strand: '+',
            trx_start: 2556364,
            trx_end: 2565622,
            coding_start: 2556664,
            coding_end: 2563273,
            num_exons: 8,
            exon_starts:
                "2556364,2557725,2558342,2559822,2560623,2561672,2562864,2563147,".to_owned(),
            exon_ends:
                "2556733,2557834,2558468,2559978,2560714,2561815,2562896,2565622,".to_owned(),
        },
        RefFlatRecord {
            gene_id: "SMIM12".to_owned(),
            transcript_name: "NM_001164824".to_owned(),
            seq_name: "chr1".to_owned(),
            strand: '-',
            trx_start: 34850361,
            trx_end: 34859045,
            coding_start: 34855698,
            coding_end: 34855977,
            num_exons: 3,
            exon_starts: "34850361,34856555,34858839,".to_owned(),
            exon_ends: "34855982,34856739,34859045,".to_owned(),
        },
        RefFlatRecord {
            gene_id: "SMIM12".to_owned(),
            transcript_name: "NM_001164825".to_owned(),
            seq_name: "chr1".to_owned(),
            strand: '-',
            trx_start: 34850361,
            trx_end: 34859737,
            coding_start: 34855698,
            coding_end: 34855977,
            num_exons: 2,
            exon_starts: "34850361,34859454,".to_owned(),
            exon_ends: "34855982,34859737,".to_owned(),
        },
        RefFlatRecord {
            gene_id: "SMIM12".to_owned(),
            transcript_name: "NM_138428".to_owned(),
            seq_name: "chr1".to_owned(),
            strand: '-',
            trx_start: 34850361,
            trx_end: 34859816,
            coding_start: 34855698,
            coding_end: 34855977,
            num_exons: 2,
            exon_starts: "34850361,34859676,".to_owned(),
            exon_ends: "34855982,34859816,".to_owned(),
        },
    ];

    let mut writer = RefFlatWriter::from_writer(vec![]);
    for rec in recs.iter() {
        writer.write_record(rec).expect("a successful write");
    }
    assert_eq!(writer.as_string(), MULT_ROWS_MULT_GENES_WITH_CDS);
}

#[test]
fn refflat_writer_transcripts_mult_rows_mult_genes_with_cds() {
    let trxs = [
        TBuilder::new("chr1", 2556364, 2565622)
            .strand(Strand::Forward)
            .coords(vec![
                (2556364, 2556733), (2557725, 2557834), (2558342, 2558468), (2559822, 2559978),
                (2560623, 2560714), (2562864, 2562896), (2563147, 2565622)],
                Some((2556664, 2562868)))
            .gene_id("TNFRSF14").id("NM_001297605")
            .coding_incl_stop(true)
            .build().expect("a transcript"),
        TBuilder::new("chr1", 2556364, 2565622)
            .strand(Strand::Forward)
            .coords(vec![
                (2556364, 2556733), (2557725, 2557834), (2558342, 2558468), (2559822, 2559978),
                (2560623, 2560714), (2561672, 2561815), (2562864, 2562896),
                (2563147, 2565622)], Some((2556664, 2563273)))
            .gene_id("TNFRSF14").id("NM_003820")
            .coding_incl_stop(true)
            .build().expect("a transcript"),
        TBuilder::new("chr1", 34850361, 34859045)
            .strand(Strand::Reverse)
            .coords(vec![
                (34850361, 34855982), (34856555, 34856739), (34858839, 34859045)],
                Some((34855698, 34855977)))
            .gene_id("SMIM12").id("NM_001164824")
            .coding_incl_stop(true)
            .build().expect("a transcript"),
        TBuilder::new("chr1", 34850361, 34859737)
            .strand(Strand::Reverse)
            .coords(vec![
                (34850361, 34855982), (34859454, 34859737)], Some((34855698, 34855977)))
            .gene_id("SMIM12").id("NM_001164825")
            .coding_incl_stop(true)
            .build().expect("a transcript"),
        TBuilder::new("chr1", 34850361, 34859816)
            .strand(Strand::Reverse)
            .coords(vec![
                (34850361, 34855982), (34859676, 34859816)], Some((34855698, 34855977)))
            .gene_id("SMIM12").id("NM_138428")
            .coding_incl_stop(true)
            .build().expect("a transcript"),
    ];

    let mut writer = RefFlatWriter::from_writer(vec![]);
    for trx in trxs.iter() {
        writer.write_transcript(trx).expect("a successful write");
    }
    assert_eq!(writer.as_string(), MULT_ROWS_MULT_GENES_WITH_CDS);
}

#[test]
fn refflat_writer_genes_mult_rows_mult_genes_with_cds() {
    let mut trxs1 = LinkedHashMap::new();
    trxs1.insert(
        "NM_001297605".to_owned(),
        ((2556364, 2565622),
            vec![(2556364, 2556733), (2557725, 2557834), (2558342, 2558468), (2559822, 2559978),
                (2560623, 2560714), (2562864, 2562896), (2563147, 2565622)],
            Some((2556664, 2562868))));

    trxs1.insert(
        "NM_003820".to_owned(),
        ((2556364, 2565622),
        vec![(2556364, 2556733), (2557725, 2557834), (2558342, 2558468), (2559822, 2559978),
                (2560623, 2560714), (2561672, 2561815), (2562864, 2562896), (2563147, 2565622)],
        Some((2556664, 2563273))));

    let mut trxs2 = LinkedHashMap::new();
    trxs2.insert(
        "NM_001164824".to_owned(),
        ((34850361, 34859045),
        vec![(34850361, 34855982), (34856555, 34856739), (34858839, 34859045)],
        Some((34855698, 34855977))));

    trxs2.insert(
        "NM_001164825".to_owned(),
        ((34850361, 34859737), vec![(34850361, 34855982), (34859454, 34859737)],
        Some((34855698, 34855977))));

    trxs2.insert(
        "NM_138428".to_owned(),
        ((34850361, 34859816), vec![(34850361, 34855982), (34859676, 34859816)],
        Some((34855698, 34855977))));

    let gxs = [
        GBuilder::new("chr1", 2556364, 2565622)
            .strand(Strand::Forward)
            .id("TNFRSF14")
            .transcript_coords(trxs1)
            .transcript_coding_incl_stop(true)
            .build().expect("a gene"),
        GBuilder::new("chr1", 34850361, 34859816)
            .strand(Strand::Reverse)
            .id("SMIM12")
            .transcript_coords(trxs2)
            .transcript_coding_incl_stop(true)
            .build().expect("a gene"),
    ];

    let mut writer = RefFlatWriter::from_writer(vec![]);
    for gx in gxs.iter() {
        writer.write_gene(gx).expect("a successful write");
    }
    assert_eq!(writer.as_string(), MULT_ROWS_MULT_GENES_WITH_CDS);
}
