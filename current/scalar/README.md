# Scalar sequence measurements and plotting

This package computes dimer-composition measurements, plots them, and partitions
reads/contigs by scalar values. The command wrappers are the supported entrypoints:

| Wrapper | Implementation | Purpose |
|---|---|---|
| `scalars.sh` | `Scalars` | Per-file metrics, or means and population standard deviations over rolling windows |
| `scalarintervals.sh` | `ScalarIntervals` | Sampled metric rows, optionally labeled by external taxonomy services |
| `cloudplot.sh` | `CloudPlot` | PNG scatter plots with metric-based axes, rotation, size and color |
| `partition.sh` | `PartitionReads3` | Current partitioning implementation |

`PartitionReads` and `PartitionReads2` remain separate historical implementations.
Do not treat similarity as permission to merge their algorithms. The partition
classes do not consume `ScalarData` tables; they compute their own metrics.
`partition.sh` writes raw, gzip, BGZF, ZIP, SAM and native BAM partition outputs,
including paired outputs and separate FASTA+QUAL files. It bounds per-partition
writer threads with `multiwritermode=auto` by default; `multiwritermode=st` and
`multiwritermode=zt` are available for explicit tuning. In this bounded path,
`.gz` uses Java gzip, `.bgz`/`.bgzip` use native BGZF, and global multithreaded
or external-compressor preferences do not create extra shard workers. Output
compression that requires external per-file compressors, such as bzip2, xz or
dsrc, is rejected for these partition shards.

`Scalars`, `ScalarIntervals` and `CloudPlot` reject `in2=` rather than silently
ignoring it. Interleaved input can carry mates through one reader. Multiple input
files accepted by Scalars/ScalarIntervals are independent files, not an `in`/`in2`
mate pairing. Output paths that refer to an input file are rejected before writing.

## Rows and windows

`ScalarData` stores parallel GC, HH, CAGA, depth, length, position and taxonomy
columns. Optional name storage preserves one slot per row, including missing names.
`ScalarInterval` is a reusable row view populated by copying primitive values and
borrowing the immutable name string. Mutable tables, row views and parsers require
external ownership; they are not thread-safe containers.

Current interval-table output is:

```
#Name  Length  GC  HH  CAGA  Depth  Start  TaxID  TaxID2  End
```

The separators are tabs. `Name` appears only with `printname=t`; `Start` and `End`
appear only with `printpos=t`, including unlimited (`window=0`) measurements.
`End` is appended after the old columns to preserve their positions. The reader uses header names, so
these optional columns do not shift the remaining fields. Keep the header when
passing modern interval tables between tools. Headerless input retains the legacy
layout: `GC HH CAGA TaxID`, `GC HH CAGA Depth TaxID`, or that last form plus `Name`.
A header can describe a legacy named table without a depth column.
If the API requests a Name column after name storage was disabled, that column
remains present with empty data fields.

GC is derived from counted dimer terminal bases. HH is a normalized homopolymer
dimer metric, not Shannon entropy. A `window` in KmerTracker retains valid dimers;
ambiguity resets dimer construction but does not flush existing window entries.
Emission spacing in ScalarIntervals counts valid input bases. Scalars' `break=t`
resets rolling state between individual sequences; the default keeps that state.
Some metrics retain
NaN for undefined denominators; Scalars' histogram conversion places those in bin
zero. These conventions are distinct from a fixed-width contiguous base window.

Sequence-generated `Start` and `End` are zero-based, half-open bounds of the
retained valid dimers: the first dimer's first base through the last dimer's last
base. They include any interior ambiguity gaps; those gaps do not contribute
dimers. Terminal ambiguities are excluded. A partial window uses its actual
dimers, and `window=0` bounds all counted dimers in the contig. No counted dimers,
or `break=f` accumulation that can mix contigs, gives `Start=End=-1`.

`Length` remains the count of newly represented valid bases since the previous
sample. It is a weight, not `End-Start`; overlapping windows must not inflate
CloudPlot's base totals. For `ACGTACGT`, `window=4 interval=2 minlen=1 break=t`
emits spans `[0,5)`, `[2,7)`, `[3,8)` with lengths `5,2,1`. Metrics, sampling
cadence and these weights are unchanged. Loaded coordinates are preserved;
legacy tables without `End` use `-1`, never an end inferred from `Length`.

Per-file Mean/STDev summaries written to stderr use a separate, shorter schema:
label, GC, HH, CAGA, depth, primary TaxID, name. The interval table goes to `out=`.

## Coverage and taxonomy

Coverage can come from pileup/covmaker tables, alignment input or supported contig
headers. ByteFile's normal input backend decodes BAM to SAM; the alignment loader
counts M, = and X CIGAR operations and divides by reference length. It includes
mapped alignments without deduplicating or filtering secondary/supplementary records.
The pileup and alignment loaders close their readers and reject reported read errors
instead of returning partial depth maps. Tadpole coverage comes from the `cov=`
field in headers produced by `Contig`; SPAdes and generic `_cov_` headers are also supported.
Depth-map/configuration state is process-wide: configure it before launching workers,
and do not run independent differently configured scalar jobs concurrently in one JVM.

Sketch/clade labeling is opt-in and can contact remote services. Routine tests use
`sketch=f clade=f`; they do not certify live service availability. Public methods
borrowing an already-started reader leave closure to the caller. Filename/descriptor
overloads own their input; owned table writers wait for completion and report errors.

## Plotting

Use three or four metrics in `order=x,y,rotation[,size]`. Taxonomy is supported for
rotation and color, not X/Y/size. `colorby=` selects the color metric independently.
Constant numeric dimensions use a central normalized value so points remain visible.
Color bounds are computed once per plot, including when color is not an axis metric.

Existing clipping percentiles use FloatList's **sum-weighted** convention rather
than a positional fraction of points. Depth/length percentiles are evaluated after
log transformation; changing that convention would be a separate behavior decision.
The legacy `colorbyname` flag currently enables name retention without implementing
a distinct name-coloring function; prefer explicit supported `colorby=` metrics.

## Validation scope

Finite regression fixtures live under `test/scalar` and `test/repeat`, separate
from production packages. Current review/build/performance evidence is maintained
in the polish project's canonical records. API8 compilation is distinct from
execution on an actual Java8 runtime. Formatting/read coverage does not certify
unexamined concurrency or external-service behavior.
