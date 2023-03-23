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
use std::fs::File;
use std::io::prelude::*;
use clap::{Parser, crate_version, crate_description};
use seq_io::fasta::{Reader,Record};
use flate2::read::GzDecoder;
use bit_vec::BitVec;


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

    let alignment_length = get_first_seq_length(&cli.input);
    eprintln!("alignment length: {}", alignment_length);
    let (a, c, g, t) = base_bitvectors(&cli.input, alignment_length);

}


/// Returns a bitvector for each of the four canonical bases for each position of the alignment.
fn base_bitvectors(filename: &PathBuf, alignment_length: usize) -> (BitVec, BitVec, BitVec, BitVec){
    let mut a = BitVec::from_elem(alignment_length, false);
    let mut c = BitVec::from_elem(alignment_length, false);
    let mut g = BitVec::from_elem(alignment_length, false);
    let mut t = BitVec::from_elem(alignment_length, false);

    let mut fasta_reader = open_fasta_file(filename);
    while let Some(record) = fasta_reader.next() {
        let record = record.expect("Error reading record");
        let name = record.id().unwrap();
        let seq = record.full_seq();
        if alignment_length != seq.len() {
            misc::quit_with_error("all sequences must be equal length");
        }
        for i in 0..alignment_length {
            match seq[i] {
                65 =>  a.set(i, true),
                97 =>  a.set(i, true),
                67 =>  c.set(i, true),
                99 =>  c.set(i, true),
                71 =>  g.set(i, true),
                103 => g.set(i, true),
                64 =>  t.set(i, true),
                116 => t.set(i, true),
                _ => (),
            }
        }
        eprintln!("{} {}", name, seq.len());
    }
    (a, c, g, t)
}


fn get_first_seq_length(filename: &PathBuf) -> usize {
    let mut fasta_reader = open_fasta_file(filename);
    while let Some(record) = fasta_reader.next() {
        let record = record.expect("Error reading record");
        return record.full_seq().len();
    }
    misc::quit_with_error("no sequences in input file");
    return 0;
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
