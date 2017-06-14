#[macro_use]
extern crate clap;

use std::process;
use std::io::{self, Write};

use clap::{App, AppSettings, ArgMatches};

mod tools;


const TEMPLATE: &'static str = "
{bin} {version}
{about}


USAGE:
    {usage}

SUBCOMMANDS:
{subcommands}

OPTIONS:
{unified}";

const ABOUT: &'static str = "
gtetools is a collection of tools for working with various gene annotation
file formats. Submit bug reports, feature requests, or view the source code
at https://github.com/bow/gtetools.";


/// Constructs a new `clap::App` for argument parsing.
pub fn build_cli<'a, 'b>() -> App<'a, 'b> {
    App::new("gtetools")
        .version(crate_version!())
        .author(crate_authors!())
        .about(ABOUT)
        .template(TEMPLATE)
        .max_term_width(80)
        .settings(&[AppSettings::GlobalVersion,
                    AppSettings::SubcommandRequiredElseHelp,
                    AppSettings::DisableHelpSubcommand,
                    AppSettings::VersionlessSubcommands])
        .subcommand(tools::gff_to_refflat::build_cli::<'a, 'b>())
        .subcommand(tools::stats::build_cli::<'a, 'b>())
}

/// Runs the appropriate tool given the subcommand argument matches.
pub fn run(matches: ArgMatches) -> Result<(), &'static str> {
    match matches.subcommand() {
        (tools::stats::NAME, Some(m)) => tools::stats::run(m),
        (tools::gff_to_refflat::NAME, Some(m)) => tools::gff_to_refflat::run(m),
        // We should not reach this point since we already require
        // that subcommands must be present in the app settings.
        _ => Err("Unexpected subcommand parsing error"),
    }
}

/// Main entry point.
fn main() {
    let matches = build_cli::<'static, 'static>().get_matches();
    if let Err(err) = run(matches) {
        let _ = writeln!(io::stderr(), "error: {}", err);
        process::exit(1);
    }
    process::exit(0);
}
