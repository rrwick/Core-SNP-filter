// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Drop-Columns

// This file is part of Drop-Columns. Drop-Columns is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Drop-Columns
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Drop-Columns. If not, see <http://www.gnu.org/licenses/>.

mod misc;

use std::path::PathBuf;
use clap::{Parser, crate_version, crate_description};
use seq_io::fasta::{Record, RefRecord};
use bitvec::prelude::*;


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

    /// Verbose output
    #[arg(long = "verbose")]
    verbose: bool,

    /// Input alignment
    input: PathBuf,
}


fn main() {
    let cli = Cli::parse();

    let alignment_length = misc::get_first_fasta_seq_length(&cli.input);
    eprintln!("alignment length:    {}", alignment_length);
    let (a, c, g, t, seq_count, acgt_counts) = bitvectors_and_counts(&cli.input, alignment_length);
    eprintln!("number of sequences: {}", seq_count);

    let mut keep = bitvec![1; alignment_length];
    if cli.verbose {
        print_verbose_header();
    }
    for i in 0..alignment_length {
        let variation = has_variation(a[i], c[i], g[i], t[i]);
        let frac = acgt_counts[i] as f64 / seq_count as f64;
        keep.set(i, keep_column(variation, frac, cli.exclude_invariant, cli.core));
        if cli.verbose {
            print_verbose_line(i, a[i], c[i], g[i], t[i], acgt_counts[i], variation, frac, keep[i]);
        }
    }
    let output_size = keep.iter().filter(|n| *n == true).count();
    eprintln!("columns to keep:     {}", output_size);

    let mut fasta_reader = misc::open_fasta_file(&cli.input);
    while let Some(record) = fasta_reader.next() {
        let record = record.expect("Error reading record");
        output_sequence(&record, &keep, output_size);
    }
}

fn output_sequence(record: &RefRecord, keep: &BitVec, output_size: usize) {
    let header = get_fasta_header(&record);
    let seq = remove_columns(&record, &keep, output_size);
    println!(">{}\n{}", header, seq);
}


fn remove_columns(record: &RefRecord, keep: &BitVec, output_size: usize) -> String {
    let full_seq = record.full_seq();
    let mut kept_seq = String::with_capacity(output_size);
    for i in 0..full_seq.len() {
        if keep[i] {
            kept_seq.push(full_seq[i] as char)
        }
    }
    assert!(kept_seq.len() == output_size);
    kept_seq
}


fn get_fasta_header(record: &RefRecord) -> String {
    let mut header = String::new();
    header += record.id().unwrap();
    match record.desc() {
        Some(x) => header += &format!(" {}", x.unwrap())[..],
        _       => (),
    }
    return header
}


fn keep_column(variation: bool, frac: f64, exclude_invariant: bool, core: f64) -> bool {
    if !variation && exclude_invariant {
        return false;
    }
    frac >= core
}


fn has_variation(a: bool, c: bool, g: bool, t: bool) -> bool {
    let total = a as i32 + c as i32 + g as i32 + t as i32;
    total > 1
}


fn print_verbose_header() {
    eprintln!();
    eprintln!("pos\ta\tc\tg\tt\tcount\tvar\tfrac\tkeep");
}


fn print_verbose_line(i: usize, a: bool, c: bool, g: bool, t: bool, acgt_counts: usize,
                      variation: bool, frac: f64, keep: bool) {
    eprintln!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}", i, a as i32, c as i32, g as i32, t as i32,
              acgt_counts, variation as i32, frac, keep as i32);}


/// Returns:
/// * a bitvector for each of the four canonical bases for each position of the alignment
/// * the number of sequences in the alignment
/// * how many of the sequences have a canonical base for each position of the alignment
fn bitvectors_and_counts(filename: &PathBuf, alignment_length: usize) -> (BitVec, BitVec, BitVec, BitVec, usize, Vec<usize>){
    let mut a = bitvec![0; alignment_length];
    let mut c = bitvec![0; alignment_length];
    let mut g = bitvec![0; alignment_length];
    let mut t = bitvec![0; alignment_length];
    let mut seq_count = 0;
    let mut acgt_counts = vec![0; alignment_length];

    let mut fasta_reader = misc::open_fasta_file(filename);
    while let Some(record) = fasta_reader.next() {
        let record = record.expect("Error reading record");
        let seq = record.full_seq();
        if alignment_length != seq.len() {
            misc::quit_with_error("all sequences must be equal length");
        }
        seq_count += 1;
        for i in 0..alignment_length {
            match seq[i] {
                65 | 97 =>  {a.set(i, true); acgt_counts[i] += 1;},
                67 | 99 =>  {c.set(i, true); acgt_counts[i] += 1;},
                71 | 103 => {g.set(i, true); acgt_counts[i] += 1;},
                84 | 116 => {t.set(i, true); acgt_counts[i] += 1;},
                _ => (),
            }
        }
    }
    (a, c, g, t, seq_count, acgt_counts)
}


#[cfg(test)]
mod tests {
    use tempfile::tempdir;
    use super::*;

    #[test]
    fn test_keep_column() {
        assert_eq!(keep_column(true, 0.9, false, 0.0), true);
        assert_eq!(keep_column(true, 0.9, true, 0.0), true);
        assert_eq!(keep_column(true, 0.9, true, 0.8), true);
        assert_eq!(keep_column(true, 0.9, true, 0.95), false);
        assert_eq!(keep_column(false, 0.9, true, 0.0), false);
        assert_eq!(keep_column(false, 0.9, true, 0.95), false);
    }

    #[test]
    fn test_has_variation() {
        assert_eq!(has_variation(false, false, false, false), false);
        assert_eq!(has_variation(true, false, false, false), false);
        assert_eq!(has_variation(false, true, false, false), false);
        assert_eq!(has_variation(false, false, true, false), false);
        assert_eq!(has_variation(false, false, false, true), false);
        assert_eq!(has_variation(true, true, false, false), true);
        assert_eq!(has_variation(false, false, true, true), true);
        assert_eq!(has_variation(true, false, true, false), true);
        assert_eq!(has_variation(false, true, false, true), true);
        assert_eq!(has_variation(false, true, true, true), true);
        assert_eq!(has_variation(true, false, true, true), true);
        assert_eq!(has_variation(true, true, false, true), true);
        assert_eq!(has_variation(true, true, true, false), true);
        assert_eq!(has_variation(true, true, true, true), true);
    }

    #[test]
    fn test_bitvectors_and_counts_1() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.fasta");
        let mut file = File::create(&file_path).unwrap();
        writeln!(file, ">seq_1\nACGAT").unwrap();
        writeln!(file, ">seq_2\nGGT-A").unwrap();
        drop(file);
        let alignment_length = get_first_seq_length(&file_path);
        assert_eq!(alignment_length, 5);
        let (a, c, g, t, seq_count, acgt_counts) = bitvectors_and_counts(&file_path, alignment_length);
        println!("{:?} {:?} {:?} {:?}", a, c, g, t);
        assert_eq!(a, bitvec![1, 0, 0, 1, 1]);
        assert_eq!(c, bitvec![0, 1, 0, 0, 0]);
        assert_eq!(g, bitvec![1, 1, 1, 0, 0]);
        assert_eq!(t, bitvec![0, 0, 1, 0, 1]);
        assert_eq!(seq_count, 2);
        assert_eq!(acgt_counts, vec![2, 2, 2, 1, 2]);
    }

    #[test]
    fn test_bitvectors_and_counts_2() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.fasta");
        let mut file = File::create(&file_path).unwrap();
        writeln!(file, ">seq_1\naacgacta").unwrap();
        writeln!(file, ">seq_2\nAGCNACGA").unwrap();
        writeln!(file, ">seq_3\nacgGCTca").unwrap();
        drop(file);
        let alignment_length = get_first_seq_length(&file_path);
        assert_eq!(alignment_length, 8);
        let (a, c, g, t, seq_count, acgt_counts) = bitvectors_and_counts(&file_path, alignment_length);
        println!("{:?} {:?} {:?} {:?}", a, c, g, t);
        assert_eq!(a, bitvec![1, 1, 0, 0, 1, 0, 0, 1]);
        assert_eq!(c, bitvec![0, 1, 1, 0, 1, 1, 1, 0]);
        assert_eq!(g, bitvec![0, 1, 1, 1, 0, 0, 1, 0]);
        assert_eq!(t, bitvec![0, 0, 0, 0, 0, 1, 1, 0]);
        assert_eq!(seq_count, 3);
        assert_eq!(acgt_counts, vec![3, 3, 3, 2, 3, 3, 3, 3]);

    }

}
