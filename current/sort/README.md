# Sorting and shuffling reads

The `sort` package contains read comparators and four command-line tools. Read
pairs move together: sorting or shuffling a primary read also moves its mate.

| Launcher | Class | Operation |
|---|---|---|
| `sortbyname.sh` | `SortByName` | Sorts reads by a selected comparator, spilling batches and merging them when needed. |
| `mergesorted.sh` | `MergeSorted` | Merges inputs that are already sorted by the selected comparator. |
| `shuffle.sh` | `Shuffle` | Loads reads into memory, then shuffles or sorts them. |
| `shuffle2.sh` | `Shuffle2` | Shuffles reads, using temporary files when the configured memory threshold is reached. |

Run a launcher with `--help` for its arguments. Sorting and shuffling accept
`in2=` and `out2=` for separate mate files, or an interleaved stream for pairs in
one file. `MergeSorted` instead takes a list of already-sorted inputs; each paired
input must contain interleaved mates (`int=t`), and `out2=` can split the output.
Specify `int=f` when single-read identifiers could be mistaken for mate identifiers.

```bash
sortbyname.sh in=reads.fq out=sorted.fq
mergesorted.sh in=sorted-a.fq,sorted-b.fq out=merged.fq delete=f
shuffle.sh in=reads.fq out=shuffled.fq
shuffle2.sh in=reads-1.fq in2=reads-2.fq out=random-1.fq out2=random-2.fq seed=17
```

`MergeSorted` requires compatible input ordering; it does not sort each input
from scratch. Its `delete=t` option removes successfully consumed inputs, so use
`delete=f` to keep the original files. Input/output path aliases are rejected.

The `random=t`/`shuffle=t` sort modes are deprecated. Random keys assigned to
reads are not preserved in temporary files, so they cannot support a correct
external merge. `SortByName` permits this legacy mode only with `allowtemp=f`;
otherwise it rejects the configuration before opening streams. `MergeSorted`
rejects random mode, including its underlying file-merge API. Use `shuffle2.sh`
for shuffling that may need temporary files.

`Shuffle` must hold all retained reads in memory. Its modes are `shuffle`,
`name`, `sequence`, `coordinate`, and `id`; either a bare mode flag or
`mode=<name>` selects one. Its `ShuffleThread` helper limits concurrent launches
and records worker failures. Callers must wait for completion and inspect each
worker's `errorState()`; completion alone does not establish success.

`Shuffle2` always shuffles. `seed=<nonnegative number>` makes repeated runs
reproducible with the same input, heap size and settings; `deterministic=t` is
shorthand for `seed=1`. Changing the memory/spill layout can change the order.
Temporary files are created in the working directory. `allowtemp=f` keeps the
operation in memory; `memmult=0.35` sets the default fraction of available memory
per input block, and `maxfiles=16` bounds merge fan-in. The accepted ranges are
`0<memmult<=1` and `maxfiles>=2`. Successful merges delete their temporary inputs
unless `delete=f` is set. Failed writes and merges are reported as errors.

All four tools reject unsupported separate `qfout`/`qfout2` quality files
before opening input or output. Write FASTQ to retain qualities;
the result can then be converted to separate sequence and quality files.

## Comparator API

`ReadComparator` extends `Comparator<Read>` with `getComparator(boolean)`,
`ascendingMult()` and `name()`. The requested direction may select a singleton
or another configured comparator; callers should not assume a fresh instance.
Each concrete comparator defines its own keys and supported input metadata.

The package includes comparators for names, numeric IDs, lengths, quality/errors,
mapping/position, flowcell coordinates, taxonomy, user-supplied lists, random
keys, clump keys, CRISPR metadata and contig length. See each class's contract
before using it outside its existing tool. Several comparators depend on shared
configuration or precomputed fields.

`ReadComparatorPosition` places reads without SAM metadata after reads with
SAM metadata in ascending order. Among SAM-backed reads, unmapped reads follow
mapped reads. Identifiers break ties within equal-position or missing-metadata
groups. Descending order reverses both the position and identifier comparisons.

`ReadComparatorTopological` compares primary bases and length before optional
mate keys, then quality and identifiers. Its public byte-vector helpers compare
only the shared prefix; equality from those helpers does not imply equal lengths.
`ReadComparatorTopological5Bit` first uses a 12-symbol key stored in `numericID`;
call `genKmer` before comparing. The key preparation overwrites the original
numeric ID. Its remaining sequence and mate comparisons also resolve primary
length before considering mate keys.

## Maintenance tests

Regression harnesses live in `current/test/sort` (package `test.sort`), separate
from production classes. They cover comparator ordering, paired merges, path
alias rejection, seeded shuffling, record and mate preservation, and worker
completion/error reporting. Test coverage is scoped to the cases each harness
constructs; it is not a guarantee for arbitrary inputs or concurrent tool instances.
