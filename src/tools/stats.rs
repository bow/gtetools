use clap::{App, Arg, ArgMatches, SubCommand};

use tools::TEMPLATE_SUBCMD;

pub const NAME: &'static str = "stats";


pub fn build_cli<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name(NAME)
        .about("Gathers various statistics")
        .template(TEMPLATE_SUBCMD)
        .arg(Arg::with_name("input")
                    .value_name("input")
                    .help("Path to input annotation file or '-' for stdin")
                    .takes_value(true)
                    .required(true))
}

pub fn run(args: &ArgMatches) -> Result<(), &'static str> {
    unimplemented!()
}
