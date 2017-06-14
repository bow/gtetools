#[macro_use]
extern crate clap;
#[macro_use]
extern crate quick_error;

use std::io::{self, Write};
use std::process;
use std::result;

use clap::{App, AppSettings, ArgMatches};

mod tools;
mod utils;


const TEMPLATE: &'static str = "
{bin} {version}
{about}


Usage: {usage}

Subcommands:
{subcommands}

Options:
{unified}";

const ABOUT: &'static str = "
gtetools is a collection of tools for working with various gene annotation
file formats. Submit bug reports, feature requests, or view the source code
at https://github.com/bow/gtetools.";


quick_error! {
    #[derive(Debug)]
    pub enum Error {
        Gte(err: gte::Error) {
            description(err.description())
            from()
            cause(err)
        }
        Io(err: io::Error) {
            description(err.description())
            from()
            cause(err)
        }
        Other(msg: &'static str) {
            description(msg)
        }
    }
}

type Result<T> = result::Result<T, Error>;


/// Constructs a new `clap::App` for argument parsing.
fn build_cli<'a, 'b>() -> App<'a, 'b> {
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
fn run(matches: ArgMatches) -> ::Result<()> {
    match matches.subcommand() {
        (tools::stats::NAME, Some(m)) => tools::stats::run(m),
        (tools::gff_to_refflat::NAME, Some(m)) => tools::gff_to_refflat::run(m),
        // We should not reach this point since we already require
        // that subcommands must be present in the app settings.
        _ => Err(Error::Other("unexpected command line parsing error")),
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
