// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Core-SNP-filter

// This file is part of Core-SNP-filter. Core-SNP-filter is free software: you can redistribute it
// and/or modify it under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any later version.
// Core-SNP-filter is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
// the GNU General Public License for more details. You should have received a copy of the GNU
// General Public License along with Core-SNP-filter. If not, see <http://www.gnu.org/licenses/>.

use std::fs::{File, metadata};
use std::io::{prelude::*, BufReader};
use std::path::Path;
use seq_io::fasta::Reader;
use flate2::read::GzDecoder;


#[cfg(not(test))]
/// For friendly error messages, this function normally just prints the error and quits. But when
/// running unit tests, it instead panics so the test can catch it.
pub fn quit_with_error(text: &str) -> ! {
    eprintln!();
    eprintln!("Error: {}", text);
    std::process::exit(1);
}
#[cfg(test)]
pub fn quit_with_error(text: &str) -> ! {
    panic!("{}", text);
}


pub fn check_if_file_is_empty(filename: &Path) {
    if let Ok(metadata) = metadata(filename) {
        if metadata.len() == 0 {
            quit_with_error(&format!("{} is empty", filename.display()));
        }
    } else {
        quit_with_error(&format!("could not access {}", filename.display()));
    }
}


pub fn check_if_file_exists(filename: &Path) {
    if !Path::new(filename).exists() {
        quit_with_error(&format!("{} does not exist", filename.display()));
    }
}


/// This function returns true if the file appears to be gzipped (based on the first two bytes) and
/// false if not. If it can't open the file or read the first two bytes, it will quit with an error
/// message.
pub fn is_file_gzipped(filename: &Path) -> bool {
    let open_result = File::open(filename);
    match open_result {
        Ok(_)  => (),
        Err(e) => quit_with_error(&format!("unable to open {}\n{}", filename.display(), e)),
    }
    let file = open_result.unwrap();

    let mut reader = BufReader::new(file);
    let mut buf = vec![0u8; 2];

    let read_result = reader.read_exact(&mut buf);
    match read_result {
        Ok(_)  => (),
        Err(e) => quit_with_error(&format!("{} is too small\n{}", filename.display(), e)),
    }

    buf[0] == 31 && buf[1] == 139
}


/// Returns an iterator over a FASTA file - works with either uncompressed or gzipped FASTAs.
pub fn open_fasta_file(filename: &Path) -> Reader<Box<dyn std::io::Read>> {
    check_if_file_exists(filename);
    check_if_file_is_empty(filename);
    let file = match File::open(filename) {
        Ok(file) => file,
        Err(e) => quit_with_error(&format!("There was a problem opening {}:\n{}",
                                           filename.display(), e)),
    };
    let reader: Box<dyn Read> = match is_file_gzipped(filename) {
        true => Box::new(GzDecoder::new(file)),
        _ => Box::new(file),
    };
    Reader::new(reader)
}


pub fn get_first_fasta_seq_length(filename: &Path) -> usize {
    let mut fasta_reader = open_fasta_file(filename);
    if let Some(record) = fasta_reader.next() {
        let record = record.expect("Error reading record");
        return record.full_seq().len();
    }
    quit_with_error(&format!("{} contains no sequences", filename.display()));
}


#[cfg(test)]
mod tests {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::fs::File;
    use std::io::Write;
    use std::path::PathBuf;
    use tempfile::{TempDir,tempdir};
    use super::*;

    fn make_test_file(contents: &str) -> (PathBuf, TempDir) {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.fasta");
        let mut file = File::create(&file_path).unwrap();
        write!(file, "{}", contents).unwrap();
        (file_path, dir)
    }

    fn make_gzipped_test_file(contents: &str) -> (PathBuf, TempDir) {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.fasta.gz");
        let mut file = File::create(&file_path).unwrap();
        let mut e = GzEncoder::new(Vec::new(), Compression::default());
        e.write_all(contents.as_bytes()).unwrap();
        let _ = file.write_all(&e.finish().unwrap());
        (file_path, dir)
    }

    #[test]
    fn test_check_if_file_is_empty_1() {
        let (path, _dir) = make_test_file(">seq_1\nACGAT\n");
        check_if_file_is_empty(&path);
    }

    #[test]
    #[should_panic]
    fn test_check_if_file_is_empty_2() {
        let (path, _dir) = make_test_file("");
        check_if_file_is_empty(&path);
    }

    #[test]
    #[should_panic]
    fn test_check_if_file_exists() {
        check_if_file_exists(&PathBuf::from("not_a_real_file"));
    }

    #[test]
    fn test_is_file_gzipped_1() {
        let (path, _dir) = make_test_file(">seq_1\nACGAT\n");
        assert!(!is_file_gzipped(&path));
    }

    #[test]
    fn test_is_file_gzipped_2() {
        let (path, _dir) = make_gzipped_test_file(">seq_1\nACGAT\n");
        assert!(is_file_gzipped(&path));
    }

    #[test]
    #[should_panic]
    fn test_is_file_gzipped_3() {
        let (path, _dir) = make_test_file("");
        is_file_gzipped(&path);
    }

    #[test]
    #[should_panic]
    fn test_is_file_gzipped_4() {
        is_file_gzipped(&PathBuf::from("not_a_real_file"));
    }

    #[test]
    fn test_get_first_fasta_seq_length_1() {
        let (path, _dir) = make_test_file(">seq_1\nACGAT\n\
                                           >seq_2\nGGTA\n\
                                           >seq_3\nCTCGCATCAG\n");
        let first_seq_len = get_first_fasta_seq_length(&path);
        assert_eq!(first_seq_len, 5);
    }

    #[test]
    #[should_panic]
    fn test_get_first_fasta_seq_length_2() {
        let (path, _dir) = make_test_file("");
        get_first_fasta_seq_length(&path);
    }

    #[test]
    #[should_panic]
    fn test_get_first_fasta_seq_length_3() {
        let (path, _dir) = make_gzipped_test_file("");
        get_first_fasta_seq_length(&path);
    }
}
