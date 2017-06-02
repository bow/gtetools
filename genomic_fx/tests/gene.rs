extern crate bio;
extern crate genomic_fx;
extern crate linked_hash_map;

use std::collections::HashMap;

use linked_hash_map::LinkedHashMap;

use genomic_fx::{Strand, GBuilder};
use Strand::*;

#[test]
fn gbuilder_basic() {
    let mut attribs = HashMap::new();
    attribs.insert("key1".to_owned(), "value1".to_owned());
    attribs.insert("key2".to_owned(), "value2".to_owned());

    let mut coords = LinkedHashMap::new();
    coords.insert("trx01".to_owned(),
                  ((100, 1000), vec![(100, 300), (400, 500), (700, 1000)], Some((200, 800))));
    coords.insert("trx02".to_owned(),
                  ((100, 1000), vec![(100, 300), (400, 550), (700, 1000)], Some((150, 900))));

    let bgx = GBuilder::new("chrT", 100, 1000)
        .strand(Forward)
        .id("gene-1")
        .transcript_coords(coords)
        .attributes(attribs)
        .build();
    assert!(bgx.is_ok(), "{:?}", bgx);
    let gx = bgx.unwrap();
    assert_eq!(gx.span(), 900);
    assert_eq!(gx.seq_name(), "chrT");
    assert_eq!(gx.strand(), &Forward);
    assert_eq!(gx.id(), Some("gene-1"));
    assert_eq!(gx.attributes().get("key1"), Some(&"value1".to_owned()));
    assert_eq!(gx.attributes().get("key2"), Some(&"value2".to_owned()));
    assert_eq!(gx.attributes().len(), 2);
    assert_eq!(gx.transcripts().len(), 2);
}
