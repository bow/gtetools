use std::io::{self, Write};

use clap::{App, Arg, ArgMatches, SubCommand};
use gte::{self, GffReader, RefFlatWriter};

use tools::TEMPLATE_SUBCMD;
use utils;

pub const NAME: &'static str = "gff-to-refflat";


pub fn build_cli<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name(NAME)
        .about("Converts from the GTF or GFF3 format to refFlat")
        .template(TEMPLATE_SUBCMD)
        .arg(Arg::with_name("input")
                .required(true)
                .takes_value(true)
                .help("Path to input annotation file or '-' for stdin"))
        .arg(Arg::with_name("output")
                .required(true)
                .takes_value(true)
                .help("Path to output annotation file or '-' for stdout"))
        .arg(Arg::with_name("gff_type")
                .short("-t")
                .long("--gff-type")
                .required(true)
                .takes_value(true)
                .value_name("TYPE")
                .possible_values(&["gff3", "gtf"])
                .display_order(1)
                .help("Input GFF variant"))
        .arg(Arg::with_name("seq_prefix")
                .long("--seq-prefix")
                .value_name("VALUE")
                .takes_value(true)
                .display_order(2)
                .help("String to prepend to all sequence names"))
        .arg(Arg::with_name("seq_lstrip")
                .long("--seq-lstrip")
                .value_name("VALUE")
                .takes_value(true)
                .display_order(3)
                .help("Left-most string to remove from all sequence names"))
        .arg(Arg::with_name("gene_id_attr")
                .long("--gid")
                .value_name("KEY")
                .default_value("gene_id")
                .takes_value(true)
                .display_order(4)
                .help("Key of GFF record attribute to use as gene identifier"))
        .arg(Arg::with_name("transcript_id_attr")
                .long("--tid")
                .value_name("KEY")
                .default_value("transcript_id")
                .display_order(5)
                .takes_value(true)
                .help("Key of GFF record attribute to use as transcript identifier"))
        .arg(Arg::with_name("loose_codons")
                .long("--loose-codons")
                .display_order(6)
                .takes_value(false)
                .long_help(
                    "If not specified, only GFF transcripts with start and stop codons will be \
                     created. If specified, GFF transcripts without start and/or stop codons \
                     will be created using the min/max coordinates of all their CDS."))
}

pub fn run(args: &ArgMatches) -> ::Result<()> {

    let gff_type = utils::resolve_gff_type(args.value_of("gff_type").unwrap())?;

    let mut reader = utils::resolve_reader(args.value_of("input").unwrap())
        .map(|r| GffReader::from_reader(r, gff_type))?;

    let mut writer = utils::resolve_writer(args.value_of("output").unwrap())
        .map(|w| RefFlatWriter::from_writer(w))?;

    let rtrxs = reader.transcripts(
        args.value_of("gene_id_attr"), args.value_of("transcript_id_attr"),
        args.value_of("seq_prefix"), args.value_of("seq_lstrip"),
        args.is_present("loose_codons"))?;

    for result in rtrxs {
        let wresult = result
            .and_then(|ref trx| writer.write_transcript(trx));
        if let Err(e) = wresult {
            if let gte::Error::Gff(gffe) = e {
                let _ = writeln!(io::stderr(), "skipping: {}", gffe);
            } else {
                return Err(::Error::from(e));
            }
        }
    }

    Ok(())
}
