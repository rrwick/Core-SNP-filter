// Copyright 2023 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/REPO-NAME

// This file is part of REPO-NAME. REPO-NAME is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. REPO-NAME
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with REPO-NAME. If not, see <http://www.gnu.org/licenses/>.

use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::{Path, PathBuf};


pub fn check_if_file_exists(filename: &PathBuf) {
    if !Path::new(filename).exists() {
        let error_message = format!("{:?} file does not exist", filename);
        quit_with_error(&error_message);
    }
}


pub fn quit_with_error(text: &str) {
    eprintln!();
    eprintln!("Error: {}", text);
    std::process::exit(1);
}


/// This function returns true if the file appears to be gzipped (based on the first two bytes) and
/// false if not. If it can't open the file or read the first two bytes, it will quit with an error
/// message.
pub fn is_file_gzipped(filename: &PathBuf) -> bool {
    let open_result = File::open(&filename);
    match open_result {
        Ok(_)  => (),
        Err(_) => quit_with_error(&format!("unable to open {:?}", filename)),
    }
    let file = open_result.unwrap();

    let mut reader = BufReader::new(file);
    let mut buf = vec![0u8; 2];

    let read_result = reader.read_exact(&mut buf);
    match read_result {
        Ok(_)  => (),
        Err(_) => quit_with_error(&format!("{:?} is too small", filename)),
    }

    buf[0] == 31 && buf[1] == 139
}
