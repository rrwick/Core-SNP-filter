# Core-SNP-filter

[![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

This is a tool to filter sites (i.e. columns) in a FASTA-format whole-genome pseudo-alignment based on:
* Whether the site contains variation or not.
* How conserved the site is, i.e. contains an unambiguous base in a sufficient fraction of the sequences.

I wrote Core-SNP-filter because I was using [Snippy](https://github.com/tseemann/snippy), and the `snippy-core` command produces a `core.full.aln` file (contains all sites regardless of variation and conservation) and `core.aln` (only contains invariant sites with 100% conservation). I wanted a tool that could produce a core SNP alignment, but with more flexibility, e.g. including sites with ≥95% conservation.

Core-SNP-filter is efficient. On a small input alignment (2 Mbp in length, 100 sequences), it runs in seconds. On a large input alignment (5 Mbp in length, 5000 sequences, 25 GB file size), it takes less than 10 minutes and only uses ~50 MB of RAM. See [this benchmark](https://github.com/mtaouk/Core-SNP-filter-methods/tree/main/benchmarking) for more details.

Important caveat: Core-SNP-filter is only appropriate for DNA alignments, not protein alignments.



## Usage

The executable named `coresnpfilter` takes a FASTA file as input. This must be an _aligned_ FASTA file, i.e. all sequences must be the same length. The characters in the FASTA sequences can be bases (e.g. `A` or `c`), gaps (`-`) or any other ASCII character (e.g. `N` for ambiguous bases or `X` for masked bases). The input FASTA can be gzipped, and line breaks (multiple lines per sequence) are okay.

There are two main options:
* `-e`/`--exclude_invariant`: if used, all invariant sites in the alignment are removed. A site counts as invariant if the number of unique unambiguous bases (`A`, `C`, `G` or `T`) at that site is one or zero. For example, a site with only `A` is invariant, but a site with both `A` and `C` is not invariant. Gaps and other characters do not count, e.g. a site with only `A`, `N` and `-` is invariant. Case does not matter, e.g. a site with only `A` and `a` is invariant.
* `-c`/`--core`: at least this fraction of the sequences must contain an unambiguous base (`A`, `C`, `G` or `T`) at a site for the site to be included. The default is `0.0`, i.e. sites are not filtered based on core fraction. If `1.0` is given, all sites with gaps or other characters will be removed, leaving an alignment containing only unambiguous bases. A more relaxed value of `0.95` will ensure that each site contains mostly unambiguous bases, but up to 5% of the bases can be gaps or other characters.

Core-SNP-filter outputs a FASTA alignment to stdout. The output will have the same number of sequences as the input, but (depending on the options used) the length of the sequences will likely be shorter. The header lines (names and descriptions) of the output will be the same as the input, and there will be no line breaks in the sequences (each sequence gets one line). Some basic information (input file, input sequence length, number of sequences and output sequence length) is printed to stderr.

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
Core-SNP-filter

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



## Installation from pre-built binaries

Core-SNP-filter compiles to a single executable binary (`coresnpfilter`), which makes installation easy!

You can find pre-built binaries for common OSs/CPUs on the [releases page](https://github.com/rrwick/Core-SNP-filter/releases). If you use one of these OSs/CPUs, download the appropriate binary for your system and put the `coresnpfilter` file in a directory that's in your `PATH` variable, e.g. `/usr/local/bin/` or `~/.local/bin/`.

Alternatively, you don't need to install Core-SNP-filter at all. Instead, just run it from wherever the `coresnpfilter` executable happens to be, like this: `/some/path/to/coresnpfilter --help`.



## Installation from source

If you are using incompatible hardware or a different OS, then you'll have to build Core-SNP-filter from source. [Install Rust](https://www.rust-lang.org/tools/install) if you don't already have it. Then clone and build Core-SNP-filter like this:
```
git clone https://github.com/rrwick/Core-SNP-filter.git
cd Core-SNP-filter
cargo build --release
```

You'll find the freshly built executable in `target/release/coresnpfilter`, which you can then move to an appropriate location that's in your `PATH` variable.



## Demo dataset

This repo's [`demo.fasta.gz`](https://raw.githubusercontent.com/rrwick/Core-SNP-filter/main/demo.fasta.gz) file is a pseudo-alignment made from 40 _Klebsiella_ samples (the original was ~5 Mbp long but I subsetted it down to 10 kbp to save space). It has many gaps, invariant sites and Ns, so Core-SNP-filter can make it a lot smaller.

For example, you can use Core-SNP-filter to create an invariant-free 95%-core alignment:
```bash
coresnpfilter -e -c 0.95 demo.fasta.gz > demo_core.fasta
```

The stderr output will look like this:
```
Core-SNP-filter
  input file:                    demo.fasta.gz
  number of sequences:           40
  input sequence length:         10000
  invariant-A sites removed:     1394
  invariant-C sites removed:     1763
  invariant-G sites removed:     1849
  invariant-T sites removed:     1378
  other invariant sites removed: 322
  non-core sites removed:        2143
  output sequence length:        1151
```

You can then build a tree with a program such as [IQ-TREE](http://www.iqtree.org):
```bash
iqtree2 -s demo_core.fasta -T 4
```



## Verbose output

Using the `--verbose` option will make Core-SNP-filter print a per-site table to stderr. Save it to file like this:
```bash
coresnpfilter -e -c 0.95 --verbose core.full.aln 1> filtered.aln 2> core_snp_table.tsv
```

This is mainly for debugging purposes, so you probably don't want to use it. But if you do, the columns are:
1. `pos`: 1-based index of the input alignment site
2. `a`: whether any sequence at this site contains `A` or `a`
3. `c`: whether any sequence at this site contains `C` or `c`
4. `g`: whether any sequence at this site contains `G` or `g`
5. `t`: whether any sequence at this site contains `T` or `t`
6. `count`: the number of sequences at this site which contain an unambiguous base
7. `frac`: the fraction of sequences at this site which contain an unambiguous base
8. `var`: whether there is any variation at this site (i.e. two or more of the `a`/`c`/`g`/`t` columns are true)
9. `keep`: whether the site passed the filter and is included in the output

Boolean columns use `0` for false and `1` for true.

For example, you can use this table to see which sites in your input alignment are included in the output alignment:
```bash
awk '{if ($9==1) print $1;}' core_snp_table.tsv
```



## Citation

[**Taouk ML, Featherstone L, Taiaroa G, Ingle DJ, Stinear TP, Wick RR. Exploring SNP Filtering Strategies: The Influence of Strict vs Soft Core. bioRxiv. 2024. doi:10.1101/2024.08.26.609800.**](https://doi.org/10.1101/2024.08.26.609800)



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
