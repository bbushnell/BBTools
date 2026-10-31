# jasper package

The `jasper` package contains command-line k-mer position counters. Both tools
load a reference FASTA into an in-memory k-mer set, scan FASTQ reads, and write a
tab-delimited table with per-position match counts and percentages for first and
second mates:

```text
#pos	read1_count	read1_perc	read2_count	read2_perc
```

The reported position is the zero-based read offset where a k-mer starts.
Percentages use each position's k-mer opportunity count as the denominator:
all complete String windows for `KmerPosition`, and only fully defined windows
for `KmerPosition3`. A position with no opportunities is reported as zero.
The output writer prints integral percentages as integers, so `100` and `50` are
expected for exact whole percentages rather than `100.000` and `50.000`.

## KmerPosition

`KmerPosition` is the original String-based implementation. It stores literal
forward reference k-mers in a `HashSet<String>` and matches only forward read
k-mers against that set.

Accepted arguments:

- `in=<reads>`: required first FASTQ input.
- `in2=<reads>`: optional second FASTQ input for paired reads.
- `ref=<fasta>`: required reference FASTA input.
- `out=<file>`: optional tabular output path.
- `k=<bases>`: k-mer length; it must be at least 1 and defaults to 6.
- Shared parser flags such as `reads=<n>` are also accepted.

## KmerPosition3

`KmerPosition3` is the 2-bit implementation. It encodes A/C/G/T as two-bit
symbols in a `long`, skips windows containing degenerate bases, and optionally
adds reverse-complemented reference k-mers before scanning reads. The live
counting behavior remains forward-read k-mers against the chosen reference set.

Accepted arguments:

- `in=<reads>`: required first FASTQ input.
- `in2=<reads>`: optional second FASTQ input for paired reads.
- `ref=<fasta>`: required reference FASTA input.
- `out=<file>`: optional tabular output path.
- `k=<bases>`: k-mer length; it must be 1-32 for the packed `long`
  representation and defaults to 19.
- `rcomp=<t|f>`: when true, also insert each reference k-mer's reverse
  complement into the reference set. Defaults to true.
- Shared parser flags such as `reads=<n>` are also accepted.

## Notes

Both tools use BBTools' standard `flag=value` syntax. Flag names are
case-insensitive, while filenames are not modified, and values may contain `=`.
The shared `reads=` limit applies separately to reference fragments loaded and
input fragments scanned; displayed read totals count both mates. Reported stream
errors cause command failure even when JVM assertions are disabled.

Neither tool opens a network listener or writes outside the requested output
file and ordinary BBTools temporary locations.
