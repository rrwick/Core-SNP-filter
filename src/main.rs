// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/REPO-NAME

// This file is part of REPO-NAME. REPO-NAME is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. REPO-NAME
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with REPO-NAME. If not, see <http://www.gnu.org/licenses/>.

mod misc;

use std::path::PathBuf;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use clap::{Parser, crate_version, crate_description};
use seq_io::fasta::{Reader,Record};

use flate2::read::GzDecoder;


#[derive(Parser)]
#[clap(name = "dropcolumns",
       version = concat!("v", crate_version!()),
       about = crate_description!())]
struct Cli {
    /// Restrict to core genome (0.0 to 1.0, default = 0.0)
    #[arg(short = 'c', long = "core", default_value = "0.0")]
    core: f64,

    /// Exclude invariant sites
    #[arg(short = 'e', long = "exclude_invariant")]
    exclude_invariant: bool,

    /// Input alignment
    input: PathBuf,
}


fn main() {
    let cli = Cli::parse();

    let mut fasta_reader = open_fasta_file(&cli.input);
    let mut alignment_length = 0;
    let mut first_record = true;

    while let Some(record) = fasta_reader.next() {
        let record = record.expect("Error reading record");
        let name = record.id().unwrap();
        let seq = record.full_seq();

        if first_record {
            alignment_length = seq.len();
            first_record = false;
        } else if alignment_length != seq.len() {
            misc::quit_with_error("all sequences must be equal length");
        }







        eprintln!("{} {}", name, seq.len());
    }
}


/// Returns an iterator over a FASTA file - works with either uncompressed or gzipped FASTAs.
fn open_fasta_file(filename: &PathBuf) -> Reader<Box<dyn std::io::Read>> {
    misc::check_if_file_exists(filename);
    let file = match File::open(filename) {
        Ok(file) => file,
        Err(error) => panic!("There was a problem opening the file: {:?}", error),
    };
    let reader: Box<dyn Read> = match misc::is_file_gzipped(filename) {
        true => Box::new(GzDecoder::new(file)),
        _ => Box::new(file),
    };
    Reader::new(reader)
}
