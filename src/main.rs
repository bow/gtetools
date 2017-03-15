#[macro_use]
extern crate clap;

use clap::{Arg, App, AppSettings, SubCommand};


fn main() {

    let matches = App::new("gfxtools")
        .version(crate_version!())
        .author(crate_authors!())
        .about("Toolset for working with gene annotation file formats.")
        .max_term_width(80)
        .settings(&[
            AppSettings::SubcommandRequired, AppSettings::GlobalVersion,
            AppSettings::DisableHelpSubcommand])
        .subcommand(
            SubCommand::with_name("stats")
                .about("Gathers various statistics")
                .arg(Arg::with_name("input")
                        .value_name("input")
                        .help("path to input annotation file or '-' for stdin")
                        .takes_value(true)
                        .required(true)))
        .subcommand(
            SubCommand::with_name("gtf2refflat")
                .about("Converts from the GTF format to refFlat")
                .arg(Arg::with_name("input")
                        .value_name("input")
                        .help("path to input annotation file or '-' for stdin")
                        .takes_value(true)
                        .required(true))
                .arg(Arg::with_name("output")
                        .value_name("output")
                        .help("path to output annotation file or '-' for stdin")
                        .takes_value(true)
                        .required(true)))
        .get_matches();

    // TODO
    match matches.subcommand() {
        ("stats", Some(_)) => {}
        ("gtf2refflat", Some(_)) => {}
        _   => {},
    }
}
