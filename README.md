# Core-SNP-filter

This is a tool to filter sites (i.e. columns) in a FASTA-format whole-genome pseudo-alignment based on:
* Whether the site contains variation or not.
* How conserved the site is, i.e. has a base in a sufficient fraction of the sequences.

I wrote this tool because I was using [Snippy](https://github.com/tseemann/snippy), and the `snippy-core` command produces a `core.full.aln` file (contains all sites regardless of variation and conservation) and `core.aln` (only contains invariant sites with 100% conservation). I wanted a tool that could produce a core SNP alignment, but with more flexibility, e.g. including sites with ≥95% conservation.



### Usage

The executable named `coresnpfilter` takes a FASTA file as input. This must be an _aligned_ FASTA file, i.e. all sequences must be the same length. The characters in the FASTA sequences can be bases (e.g. `A` or `c`), gaps (`-`) or any other ASCII character (e.g. `N` for ambiguous bases or `X` for masked bases). The input FASTA can be gzipped, and line breaks (multiple lines per sequence) is okay.

There are two main options:
* `-e`/`--exclude_invariant`: if used, all invariant sites in the alignment are removed. A site counts as invariant if the number of unique unambiguous bases (`A`, `C`, `G` or `T`) at that site is one or zero. For example, a site with only `A` is invariant, but a site with both `A` and `C` is not invariant. Gaps and other characters do not count, e.g. a site with only `A`, `N` and `-` is invariant. Case does not matter, e.g. a site with only `A` and `a` is invariant. 
* `-c`/`--core`: at least this fraction of the sequences must contain an unambiguous base (`A`, `C`, `G` or `T`) at a site for the site to be included. The default is `0.0`, i.e. sites are not filtered based on core fraction. If `1.0` is given, all sites with gaps or other characters will be removed, leaving an alignment containing only unambiguous bases. A more relaxed value of `0.95` will ensure that each site contains mostly unambiguous bases, but up to 5% of the sequences can be gaps or other characters.

`coresnpfilter` outputs a FASTA alignment to stdout. The output will have the same number of sequences as the input, but (depending on the options used) the length of the sequences will likely be shorter. The header lines (names and descriptions) of the output will be the same as the input, and there will be no line breaks in the sequences (each sequence gets one line).

Some example commands:
```bash
# Exclude invariant sites:
coresnpfilter -e core.full.aln > filtered.aln

# With a strict core threshold (same as Snippy's core.aln):
coresnpfilter -e -c 1.0 core.full.aln > filtered.aln

# With a slightly more relaxed core threshold:
coresnpfilter -e -c 0.95 core.full.aln > filtered.aln

# Use gzipped files to save disk space:
coresnpfilter -e -c 0.95 core.full.aln.gz | gzip > filtered.aln.gz

# Running without any options will work, but the output will be the same as the input:
coresnpfilter core.full.aln > filtered.aln
```

Full help text:
```
Core SNP filter

Usage: coresnpfilter [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Input alignment

Options:
  -c, --core <CORE>        Restrict to core genome (0.0 to 1.0, default = 0.0) [default: 0.0]
  -e, --exclude_invariant  Exclude invariant sites
      --verbose            Verbose output
  -h, --help               Print help
  -V, --version            Print version
```



### Installation from pre-built binaries

Core-SNP-filter compiles to a single executable binary (`coresnpfilter`), which makes installation easy!

You can find pre-built binaries for common OSs/CPUs on the [releases page](https://github.com/rrwick/Core-SNP-filter/releases). If you use one of these OSs/CPUs, download the appropriate binary for your system and put the `coresnpfilter` file in a directory that's in your `PATH` variable, e.g. `/usr/local/bin/` or `~/.local/bin/`.

Alternatively, you don't need to install Core-SNP-filter at all. Instead, just run it from wherever the `coresnpfilter` executable happens to be, like this: `/some/path/to/coresnpfilter --help`.



### Installation from source

If you are using incompatible hardware or a different OS, then you'll have to build Core-SNP-filter from source. [Install Rust](https://www.rust-lang.org/tools/install) if you don't already have it. Then clone and build Core-SNP-filter like this:
```
git clone https://github.com/rrwick/Core-SNP-filter.git
cd Core-SNP-filter
cargo build --release
```

You'll find the freshly built executable in `target/release/coresnpfilter`, which you can then move to an appropriate location that's in your `PATH` variable.



### License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
