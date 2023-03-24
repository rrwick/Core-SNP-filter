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

use bitvec::prelude::*;
use clap::{Parser, crate_version, crate_description};
use seq_io::fasta::{Record, RefRecord};
use std::io;
use std::path::PathBuf;


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
    check_arguments(cli.core);
    drop_columns(&cli.input, cli.exclude_invariant, cli.core, cli.verbose, &mut io::stdout());
}


/// This is the primary function of the program. For easier testing, I factored it out of the main
/// function and use the stdout argument to allow for capturing the output.
fn drop_columns(filename: &PathBuf, exclude_invariant: bool, core: f64, verbose: bool,
                stdout: &mut dyn io::Write) {
    let alignment_length = misc::get_first_fasta_seq_length(filename);
    eprintln!("alignment length:    {}", alignment_length);
    let (a, c, g, t, seq_count, acgt_counts) = bitvectors_and_counts(filename, alignment_length);
    eprintln!("number of sequences: {}", seq_count);

    let mut keep = bitvec![1; alignment_length];
    if verbose {
        print_verbose_header();
    }
    for i in 0..alignment_length {
        let variation = has_variation(a[i], c[i], g[i], t[i]);
        let frac = acgt_counts[i] as f64 / seq_count as f64;
        keep.set(i, keep_column(variation, frac, exclude_invariant, core));
        if verbose {
            print_verbose_line(i, a[i], c[i], g[i], t[i], acgt_counts[i], variation, frac, keep[i]);
        }
    }
    let output_size = keep.iter().filter(|n| *n == true).count();
    eprintln!("columns to keep:     {}", output_size);

    let mut fasta_reader = misc::open_fasta_file(filename);
    while let Some(record) = fasta_reader.next() {
        let record = record.expect("Error reading record");
        output_sequence(&record, &keep, output_size, stdout);
    }
}


fn check_arguments(core: f64) {
    if core < 0.0 || core > 1.0 {
        panic!("--core must be between 0 and 1 (inclusive)")
    }
}


fn output_sequence(record: &RefRecord, keep: &BitVec, output_size: usize,
                   stdout: &mut dyn io::Write) {
    let header = get_fasta_header(&record);
    let seq = remove_columns(&record, &keep, output_size);
    writeln!(stdout, ">{}\n{}", header, seq).unwrap();
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
            panic!("all sequences must be equal length");
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
    use std::fs::File;
    use std::io::Write;
    use std::str::from_utf8;
    use tempfile::{TempDir,tempdir};
    use super::*;

    fn make_test_file(contents: &str) -> (PathBuf, TempDir) {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.fasta");
        let mut file = File::create(&file_path).unwrap();
        write!(file, "{}", contents).unwrap();
        (file_path, dir)
    }

    #[test]
    fn test_check_arguments_1() {
        check_arguments(0.0);
        check_arguments(0.5);
        check_arguments(1.0);
    }

    #[test]
    #[should_panic]
    fn test_check_arguments_2() {
        check_arguments(-0.1);
    }

    #[test]
    #[should_panic]
    fn test_check_arguments_3() {
        check_arguments(1.1);
    }

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
        let (path, _dir) = make_test_file(">seq_1\nACGAT\n\
                                           >seq_2\nGGT-A\n");
        let (a, c, g, t, seq_count, acgt_counts) = bitvectors_and_counts(&path, 5);
        assert_eq!(a, bitvec![1, 0, 0, 1, 1]);
        assert_eq!(c, bitvec![0, 1, 0, 0, 0]);
        assert_eq!(g, bitvec![1, 1, 1, 0, 0]);
        assert_eq!(t, bitvec![0, 0, 1, 0, 1]);
        assert_eq!(seq_count, 2);
        assert_eq!(acgt_counts, vec![2, 2, 2, 1, 2]);
    }

    #[test]
    fn test_bitvectors_and_counts_2() {
        let (path, _dir) = make_test_file(">seq_1\naacgacta\n\
                                           >seq_2\nAGCNACGA\n\
                                           >seq_3\nacgGCTca\n");
        let (a, c, g, t, seq_count, acgt_counts) = bitvectors_and_counts(&path, 8);
        assert_eq!(a, bitvec![1, 1, 0, 0, 1, 0, 0, 1]);
        assert_eq!(c, bitvec![0, 1, 1, 0, 1, 1, 1, 0]);
        assert_eq!(g, bitvec![0, 1, 1, 1, 0, 0, 1, 0]);
        assert_eq!(t, bitvec![0, 0, 0, 0, 0, 1, 1, 0]);
        assert_eq!(seq_count, 3);
        assert_eq!(acgt_counts, vec![3, 3, 3, 2, 3, 3, 3, 3]);
    }

    #[test]
    fn test_drop_columns_1() {
        // No filtering - input is the same as the output.
        let (path, _dir) =       make_test_file(">seq_1\nACGATCAG\n\
                                                 >seq_2\nACCATTAG\n\
                                                 >seq_3\nACGATCAG\n");
        let mut stdout = Vec::new();
        drop_columns(&path, false, 0.0, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\nACGATCAG\n\
                                                 >seq_2\nACCATTAG\n\
                                                 >seq_3\nACGATCAG\n");
    }

    #[test]
    fn test_drop_columns_2() {
        // Dropping invariant sites.
        let (path, _dir) =       make_test_file(">seq_1\nACGATCAG\n\
                                                 >seq_2\nACCATTAG\n\
                                                 >seq_3\nACGATCAG\n");
        let mut stdout = Vec::new();
        drop_columns(&path, true, 0.0, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\nGC\n\
                                                 >seq_2\nCT\n\
                                                 >seq_3\nGC\n");
    }

    #[test]
    fn test_drop_columns_3() {
        // At 60% core, 2 out of 3 sequences is enough.
        let (path, _dir) =       make_test_file(">seq_1\nACGATCAG\n\
                                                 >seq_2\nAC----CG\n\
                                                 >seq_3\nAGGATCAG\n");
        let mut stdout = Vec::new();
        drop_columns(&path, false, 0.6, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\nACGATCAG\n\
                                                 >seq_2\nAC----CG\n\
                                                 >seq_3\nAGGATCAG\n");
    }

    #[test]
    fn test_drop_columns_4() {
        // At 70% core, 2 out of 3 sequences is not enough.
        let (path, _dir) =       make_test_file(">seq_1\nACGATCAG\n\
                                                 >seq_2\nAC----CG\n\
                                                 >seq_3\nAGGATCAG\n");
        let mut stdout = Vec::new();
        drop_columns(&path, false, 0.7, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\nACAG\n\
                                                 >seq_2\nACCG\n\
                                                 >seq_3\nAGAG\n");
    }

    #[test]
    fn test_drop_columns_5() {
        // Same as previous but dropping invariant sites and with verbose output.
        let (path, _dir) =       make_test_file(">seq_1\nACGATCAG\n\
                                                 >seq_2\nAC----CG\n\
                                                 >seq_3\nAGGATCAG\n");
        let mut stdout = Vec::new();
        drop_columns(&path, true, 0.7, true, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\nCA\n\
                                                 >seq_2\nCC\n\
                                                 >seq_3\nGA\n");
    }

    #[test]
    fn test_drop_columns_6() {
        // Same as previous but with some descriptions in the FASTA headers.
        let (path, _dir) =       make_test_file(">seq_1 info\nACGATCAG\n\
                                                 >seq_2\nAC----CG\n\
                                                 >seq_3 lots of stuff\nAGGATCAG\n");
        let mut stdout = Vec::new();
        drop_columns(&path, true, 0.7, true, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1 info\nCA\n\
                                                 >seq_2\nCC\n\
                                                 >seq_3 lots of stuff\nGA\n");
    }

    #[test]
    #[should_panic]
    fn test_drop_columns_7() {
        // Invalid input with different sequence lengths.
        let (path, _dir) =       make_test_file(">seq_1\nACGATCAG\n\
                                                 >seq_2\nAC----CGA\n\
                                                 >seq_3\nAGGATCAG\n");
        let mut stdout = Vec::new();
        drop_columns(&path, true, 0.7, true, &mut stdout);
    }

    #[test]
    fn test_drop_columns_8() {
        // Every column is dropped.
        let (path, _dir) =       make_test_file(">seq_1\nACGATCA-\n\
                                                 >seq_2\nAC----AC\n\
                                                 >seq_3\nACGATCAG\n");
        let mut stdout = Vec::new();
        drop_columns(&path, true, 0.7, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\n\n\
                                                 >seq_2\n\n\
                                                 >seq_3\n\n");
    }

    #[test]
    fn test_drop_columns_9() {
        // Using a mixture of uppercase and lowercase - no columns dropped.
        let (path, _dir) =       make_test_file(">seq_1\nACGAtCaGcAaT\n\
                                                 >seq_2\nAcGaGCaGcAcT\n\
                                                 >seq_3\nACGatTAgCaCT\n");
        let mut stdout = Vec::new();
        drop_columns(&path, false, 0.5, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\nACGAtCaGcAaT\n\
                                                 >seq_2\nAcGaGCaGcAcT\n\
                                                 >seq_3\nACGatTAgCaCT\n");
    }

    #[test]
    fn test_drop_columns_10() {
        // Using a mixture of uppercase and lowercase - invariant columns dropped.
        let (path, _dir) =       make_test_file(">seq_1\nACGAtCaGcAaT\n\
                                                 >seq_2\nAcGaGCaGcAcT\n\
                                                 >seq_3\nACGatTAgCaCT\n");
        let mut stdout = Vec::new();
        drop_columns(&path, true, 0.5, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\ntCa\n\
                                                 >seq_2\nGCc\n\
                                                 >seq_3\ntTC\n");
    }

    #[test]
    fn test_drop_columns_11() {
        // Using a mixture of uppercase and lowercase - non-core columns dropped.
        let (path, _dir) =       make_test_file(">seq_1\nACG--CaGcAaT\n\
                                                 >seq_2\nAcGaGCa--AcT\n\
                                                 >seq_3\nACGa----CaCT\n");
        let mut stdout = Vec::new();
        drop_columns(&path, false, 0.5, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\nACG-CacAaT\n\
                                                 >seq_2\nAcGaCa-AcT\n\
                                                 >seq_3\nACGa--CaCT\n");
    }

    #[test]
    fn test_drop_columns_12() {
        // Using a mixture of uppercase and lowercase - invariant and non-core columns dropped.
        let (path, _dir) =       make_test_file(">seq_1\nACG--CaGcAaT\n\
                                                 >seq_2\nAcGaGCa--AcT\n\
                                                 >seq_3\nACGa----CaCT\n");
        let mut stdout = Vec::new();
        drop_columns(&path, true, 0.5, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\na\n\
                                                 >seq_2\nc\n\
                                                 >seq_3\nC\n");
    }

    #[test]
    fn test_drop_columns_13() {
        // Testing an input with line breaks in the FASTA sequences.
        let (path, _dir) = make_test_file(">seq_1\nACG--\nCaGcA\naT\n\
                                                 >seq_2\nAcGaG\nCa--A\ncT\n\
                                                 >seq_3\nACGa-\n---Ca\nCT\n");
        let mut stdout = Vec::new();
        drop_columns(&path, false, 0.5, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\nACG-CacAaT\n\
                                                 >seq_2\nAcGaCa-AcT\n\
                                                 >seq_3\nACGa--CaCT\n");
    }
}
