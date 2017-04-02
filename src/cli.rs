use clap::{App, AppSettings, Arg, ArgMatches, SubCommand};
use tools;

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
gnftools is a collection of tools for working with various gene annotation
file formats. Submit bug reports, feature requests, or view the source code
at https://github.com/bow/gnftools.";

const TEMPLATE_SUBCMD: &'static str = "
USAGE:
    {usage}

ARGS:
{positionals}

OPTIONS:
{unified}";

/// Constructs a new `clap::App` for argument parsing.
pub fn build_cli() -> App<'static, 'static> {
    App::new("gnftools")
        .version(crate_version!())
        .author(crate_authors!())
        .about(ABOUT)
        .template(TEMPLATE)
        .max_term_width(80)
        .settings(&[AppSettings::GlobalVersion,
                    AppSettings::SubcommandRequiredElseHelp,
                    AppSettings::DisableHelpSubcommand,
                    AppSettings::VersionlessSubcommands])
        .subcommand(SubCommand::with_name("stats")
                        .about("Gathers various statistics")
                        .template(TEMPLATE_SUBCMD)
                        .arg(Arg::with_name("input")
                                 .value_name("input")
                                 .help("Path to input annotation file or '-' for stdin")
                                 .takes_value(true)
                                 .required(true)))
        .subcommand(SubCommand::with_name("gtf2refflat")
                        .about("Converts from the GTF format to refFlat")
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
                                 .required(true)))
}

/// Runs the appropriate tool given the subcommand argument matches.
pub fn run(matches: ArgMatches) -> Result<(), &'static str> {
    match matches.subcommand() {
        ("stats", Some(m)) => tools::stats::run(m),
        ("gtf2refflat", Some(m)) => tools::convert_gtf2refflat::run(m),
        // We should not reach this point since we already require
        // that subcommands must be present in the app settings.
        _ => Err("Unexpected subcommand parsing error"),
    }
}
