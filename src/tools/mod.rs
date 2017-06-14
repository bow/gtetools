//! Functions invoked by the subcommands.

pub mod stats;
pub mod gff_to_refflat;

const TEMPLATE_SUBCMD: &'static str = "
Usage: {usage}

Arguments:
{positionals}

Options:
{unified}";
