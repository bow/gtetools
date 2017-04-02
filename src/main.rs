#[macro_use]
extern crate clap;

use std::process;
use std::io::{self, Write};

mod cli;
mod tools;


/// Main entry point.
fn main() {
    let matches = cli::build_cli().get_matches();
    if let Err(err) = cli::run(matches) {
        let _ = writeln!(io::stderr(), "error: {}", err);
        process::exit(1);
    }
    process::exit(0);
}
