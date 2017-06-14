//! Functions invoked by the subcommands.

pub mod stats;
pub mod gff_to_refflat;

const TEMPLATE_SUBCMD: &'static str = "
USAGE:
    {usage}

ARGS:
{positionals}

OPTIONS:
{unified}";
