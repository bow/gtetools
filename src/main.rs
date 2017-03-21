#[macro_use]
extern crate clap;
extern crate bio;

use std::process;
use std::io::{self, Write};

mod cli;
mod feature;
mod tools;


fn main() {
    let matches = cli::build_cli().get_matches();
    if let Err(err) = cli::run(matches) {
        let _ = writeln!(io::stderr(), "error: {}", err);
        process::exit(1);
    }
    process::exit(0);
}
