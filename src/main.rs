// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Core-SNP-filter

// This file is part of Core-SNP-filter. Core-SNP-filter is free software: you can redistribute it
// and/or modify it under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any later version.
// Core-SNP-filter is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
// the GNU General Public License for more details. You should have received a copy of the GNU
// General Public License along with Core-SNP-filter. If not, see <http://www.gnu.org/licenses/>.

mod misc;

use bitvec::prelude::*;
use clap::{Parser, crate_version, crate_description};
use seq_io::fasta::{Record, RefRecord};
use std::io;
use std::path::{Path, PathBuf};


#[derive(Parser)]
#[clap(name = "Core-SNP-filter",
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
fn drop_columns(filename: &Path, exclude_invariant: bool, core: f64, verbose: bool,
                stdout: &mut dyn io::Write) {
    let alignment_length = misc::get_first_fasta_seq_length(filename);
    let max_width = alignment_length.to_string().len();
    let (a, c, g, t, seq_count, acgt_counts) = bitvectors_and_counts(filename, alignment_length);
    if !verbose {
        stderr_display_1(filename, max_width, seq_count, alignment_length);
    }

    let mut keep = bitvec![1; alignment_length];
    let (mut inv_a, mut inv_c, mut inv_g, mut inv_t, mut inv_other) = (0, 0, 0, 0, 0);
    let mut non_core = 0;
    if verbose {
        print_verbose_header();
    }
    for i in 0..alignment_length {
        let variation = has_variation(a[i], c[i], g[i], t[i]);
        let frac = acgt_counts[i] as f64 / seq_count as f64;
        if exclude_invariant && !variation {
            keep.set(i, false);
            if a[i] { inv_a += 1; }
            else if c[i] { inv_c += 1; }
            else if g[i] { inv_g += 1; }
            else if t[i] { inv_t += 1; }
            else { inv_other += 1; }
        }
        if keep[i] && frac < core {
            keep.set(i, false);
            non_core += 1;
        }
        if verbose {
            print_verbose_line(i, a[i], c[i], g[i], t[i], acgt_counts[i], variation, frac, keep[i]);
        }
    }
    let output_size = keep.iter().filter(|n| *n == true).count();
    let inv_total = inv_a + inv_c + inv_g + inv_t + inv_other;
    let removed_total = inv_total + non_core;
    assert!(alignment_length == output_size + removed_total);
    if !verbose {
        stderr_display_2(max_width, output_size, removed_total, non_core, inv_total,
                         inv_a, inv_c, inv_g, inv_t, inv_other);
    }

    let mut fasta_reader = misc::open_fasta_file(filename);
    while let Some(record) = fasta_reader.next() {
        let record = record.expect("Error reading record");
        output_sequence(&record, &keep, output_size, stdout);
    }
}


fn check_arguments(core: f64) {
    if !(0.0..=1.0).contains(&core) {
        misc::quit_with_error("--core must be between 0 and 1 (inclusive)");
    }
}


fn output_sequence(record: &RefRecord, keep: &BitVec, output_size: usize,
                   stdout: &mut dyn io::Write) {
    let header = get_fasta_header(record);
    let seq = remove_columns(record, keep, output_size);
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
    if let Some(x) = record.desc() {
        header += &format!(" {}", x.unwrap());
    }
    header
}


fn has_variation(a: bool, c: bool, g: bool, t: bool) -> bool {
    let total = a as i32 + c as i32 + g as i32 + t as i32;
    total > 1
}


fn stderr_display_1(filename: &Path, max_width: usize, seq_count: usize, alignment_length: usize) {
    eprintln!();
    eprintln!("Core-SNP-filter");
    eprintln!("{}", "─".repeat(max_width+37));
    eprintln!("input file: {:>w$}", filename.display(), w = max_width+25);
    eprintln!("number of sequences:                 {:>w$}", seq_count, w = max_width);
    eprintln!("input sequence length:               {:>w$}", alignment_length, w = max_width);
}


fn stderr_display_2(max_width: usize, output_size: usize, removed_total: usize, non_core: usize,
                    inv_total: usize, inv_a: usize, inv_c: usize, inv_g: usize, inv_t: usize,
                    inv_other: usize) {
    eprintln!("├ output sequence length:            {:>w$}", output_size, w = max_width);
    eprintln!("└ total sites removed:               {:>w$}", removed_total, w = max_width);
    eprintln!("  ├ non-core sites removed:          {:>w$}", non_core, w = max_width);
    eprintln!("  └ invariant sites removed:         {:>w$}", inv_total, w = max_width);
    eprintln!("    ├ invariant-A sites removed:     {:>w$}", inv_a, w = max_width);
    eprintln!("    ├ invariant-C sites removed:     {:>w$}", inv_c, w = max_width);
    eprintln!("    ├ invariant-G sites removed:     {:>w$}", inv_g, w = max_width);
    eprintln!("    ├ invariant-T sites removed:     {:>w$}", inv_t, w = max_width);
    eprintln!("    └ other invariant sites removed: {:>w$}", inv_other, w = max_width);
    eprintln!();
}


fn print_verbose_header() {
    eprintln!("pos\ta\tc\tg\tt\tcount\tfrac\tvar\tkeep");
}


fn print_verbose_line(i: usize, a: bool, c: bool, g: bool, t: bool, acgt_counts: usize,
                      variation: bool, frac: f64, keep: bool) {
    eprintln!("{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}", i+1, a as i32, c as i32, g as i32, t as i32,
              acgt_counts, frac, variation as i32, keep as i32);
}


/// Returns:
/// * a bitvector for each of the four canonical bases for each position of the alignment
/// * the number of sequences in the alignment
/// * how many of the sequences have a canonical base for each position of the alignment
fn bitvectors_and_counts(filename: &Path, alignment_length: usize)
        -> (BitVec, BitVec, BitVec, BitVec, usize, Vec<usize>){
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

    #[test]
    fn test_drop_columns_14() {
        // Testing lots of non-base characters.
        let (path, _dir) =       make_test_file(">seq_1\nAC---CGG\n\
                                                 >seq_2\nCCCNNNNG\n\
                                                 >seq_3\nACXQVPAG\n");
        let mut stdout = Vec::new();
        drop_columns(&path, true, 0.0, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\nAG\n\
                                                 >seq_2\nCN\n\
                                                 >seq_3\nAA\n");
    }

    #[test]
    fn test_drop_columns_15() {
        // Testing lots of non-base characters.
        let (path, _dir) =       make_test_file(">seq_1\nAC---CGG\n\
                                                 >seq_2\nCCCNNNNG\n\
                                                 >seq_3\nACXQVPAG\n");
        let mut stdout = Vec::new();
        drop_columns(&path, true, 1.0, false, &mut stdout);
        assert_eq!(from_utf8(&stdout).unwrap(), ">seq_1\nA\n\
                                                 >seq_2\nC\n\
                                                 >seq_3\nA\n");
    }
}
