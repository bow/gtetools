use std::fs;
use std::io::{self, BufReader, BufWriter, Read, Write};

use gte::GffType;

use Error;


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

pub fn resolve_gff_type(raw_arg: &str) -> ::Result<GffType> {

    match raw_arg.to_owned().to_lowercase().as_str() {
        "gtf" | "gtf2" => Ok(GffType::GTF2),
        "gff2" => Ok(GffType::GFF2),
        "gff3" | "gff" => Ok(GffType::GFF3),
        _ => Err(Error::Other("invalid gff type")),
    }
}
