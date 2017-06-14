use clap::{App, Arg, ArgMatches, SubCommand};

use tools::TEMPLATE_SUBCMD;

pub const NAME: &'static str = "gff-to-refflat";


pub fn build_cli<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name(NAME)
        .about("Converts from the GTF or GFF3 format to refFlat")
        .template(TEMPLATE_SUBCMD)
        .arg(Arg::with_name("input")
                .value_name("input")
                .help("Path to input annotation file or '-' for stdin")
                .takes_value(true)
                .required(true))
        .arg(Arg::with_name("output")
                .value_name("output")
                .help("Path to output annotation file or '-' for stdout")
                .takes_value(true)
                .required(true))
}

pub fn run(args: &ArgMatches) -> Result<(), &'static str> {
    unimplemented!()
}
