use std::cmp::{max, min};
use std::convert::AsRef;
use std::io;
use std::fs;
use std::path::Path;

use bio::io::gff::{self, GffType};
use itertools::{GroupBy, Group, Itertools};
use linked_hash_map::LinkedHashMap;

// use feature::FeatureError;
use {Coord, Gene, GBuilder, Strand, Error};


pub struct Reader<R: io::Read> {
    inner: gff::Reader<R>,
}

impl<R: io::Read> Reader<R> {

    pub fn from_reader(in_reader: R, gff_type: GffType) -> Reader<R> {
        Reader {
            inner: gff::Reader::new(in_reader, gff_type),
        }
    }

    pub fn genes<'a>(&'a mut self) -> GffGenes<'a, R> {
        GffGenes {
            inner: self.records().group_by(GffGenes::<'a, R>::gene_groupf),
        }
    }

    fn records<'a>(&'a mut self) -> GffRecords<'a, R> {
        GffRecords {
            inner: self.inner.records()
        }
    }
}

impl Reader<fs::File> {
    pub fn from_file<P: AsRef<Path>>(path: P, gff_type: GffType) -> io::Result<Self> {
        fs::File::open(path).map(|file| Reader::from_reader(file, gff_type))
    }
}

struct GffRecords<'a, R: 'a> where R: io::Read {
    inner: gff::Records<'a, R>,
}

impl<'a, R> Iterator for GffRecords<'a, R> where R: io::Read {

    type Item = Result<gff::Record, Error>;

    fn next(&mut self) -> Option<Result<gff::Record, Error>> {
        self.inner.next()
            .map(|row| {
                row.or_else(|err| Err(Error::from(err)))
            })
    }
}

pub struct GffGenes<'a, R: 'a> where R: io::Read, {
    inner: GroupBy<GxGroupKey, GffRecords<'a, R>, GxGroupFunc>,
}

type GxGroupKey = Option<(String, String, Strand)>;

type GxGroupFunc = fn(&Result<gff::Record, Error>) -> GxGroupKey;

type GxGroupedRecs<'a, 'b, R> = Group<'b, GxGroupKey, GffRecords<'a, R>, GxGroupFunc>;

type RawTrxCoord = (Coord<u64>, Vec<Coord<u64>>, Option<Coord<u64>>);

const INIT_COORD: (u64, u64) = (::std::u64::MAX, ::std::u64::MIN);

#[inline(always)]
fn adjust_coord(cur_coord: Coord<u64>, record: &gff::Record) -> Coord<u64> {
    (min(cur_coord.0, *record.start() - 1),
     max(cur_coord.1, *record.end()))
}

impl<'a, R> GffGenes<'a, R> where R: io::Read {

    fn gene_groupf(result: &Result<gff::Record, Error>) -> GxGroupKey {
        result.as_ref().ok()
            .map(|res| {
                let gene_id = res.attributes().get("gene_id")
                    .map(|gid| gid.as_str())
                    .unwrap_or("");
                let seq_name = String::from(res.seqname());
                let strand = res.strand().unwrap_or(Strand::Unknown);

                (String::from(gene_id), seq_name, strand)
            })
    }

    fn group_to_gene<'b>(group: (GxGroupKey, GxGroupedRecs<'a, 'b, R>)) -> Result<Gene, Error> {
        let (group_key, records) = group;
        match group_key {

            None => Err(records.filter_map(|x| x.err()).next().unwrap()),

            Some((gid, seq_name, strand)) => {

                let mut gene_coord = INIT_COORD;
                let mut trx_coords: LinkedHashMap<String, RawTrxCoord> = LinkedHashMap::new();

                for record in records {

                    let mut rec = record?;
                    rec.attributes_mut().remove("gene_id");

                    let rtrx_entry = rec.attributes_mut()
                        .remove("transcript_id")
                        .ok_or(Error::Gff("required 'transcript_id' attribute not found"))
                        .map(|id| trx_coords.entry(id).or_insert((INIT_COORD, vec![], None)));

                    match rec.feature_type() {
                        "gene" => {
                            gene_coord = adjust_coord(gene_coord, &rec);
                        },
                        "transcript" => {
                            let trx_entry = rtrx_entry?;
                            trx_entry.0 = adjust_coord(trx_entry.0, &rec);
                        },
                        "exon" => {
                            let trx_entry = rtrx_entry?;
                            (trx_entry.1).push((*rec.start() - 1, *rec.end()));
                        },
                        "cds" | "CDS" => {
                            let trx_entry = rtrx_entry?;
                            trx_entry.2 = (trx_entry.2).or(Some(INIT_COORD))
                                .map(|coord| adjust_coord(coord, &rec));
                        },
                        // TODO: How to best handle nested features? They may be out of order ...
                        _ => {},
                    }
                }

                GBuilder::new(seq_name, gene_coord.0, gene_coord.1)
                    .id(gid)
                    .strand(strand)
                    .transcript_coords(trx_coords)
                    .build()
            },
        }
    }
}

impl<'a, R> Iterator for GffGenes<'a, R> where R: io::Read {

    type Item = Result<Gene, Error>;

    fn next(&mut self) -> Option<Result<Gene, Error>> {
        self.inner.into_iter().map(Self::group_to_gene).next()
    }
}
