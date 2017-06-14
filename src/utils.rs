use std::fs;
use std::io::{self, BufReader, BufWriter, Read, Write};


const STREAM_ARG: &'static str = "-";


pub fn resolve_reader(raw_arg: &str) -> ::Result<Box<Read>>
{
    match raw_arg {
        STREAM_ARG => Ok(Box::new(io::stdin())),
        path => fs::File::open(path)
            .map_err(::Error::from)
            .map(|file| Box::new(BufReader::new(file)) as Box<Read>)
    }
}

pub fn resolve_writer(raw_arg: &str) -> ::Result<Box<Write>>
{
    match raw_arg {
        STREAM_ARG => Ok(Box::new(io::stdout())),
        path => fs::File::create(path)
            .map_err(::Error::from)
            .map(|file| Box::new(BufWriter::new(file)) as Box<Write>)
    }
}
