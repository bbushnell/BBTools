# Flow-cell statistics and read-header tools

The `hiseq` package groups sequencing reads by lane, tile and spatial coordinates,
collects quality statistics, and provides header-conversion and report utilities.
Its plotting tools write tables for subsequent plotting; they do not render
heatmaps or charts themselves. Binary CBCL conversion lives in `illumina`.

## Command-line tools

Run the BBTools shell launchers, which configure the classpath, heap and assertions.
Run a launcher without arguments for its usage text.

| Launcher | Java class | Purpose |
| --- | --- | --- |
| `filterbytile.sh` | `AnalyzeFlowCell` | Collect spatial quality statistics and filter reads using microtile decisions. |
| `plotflowcell.sh` | `PlotFlowCell` | Collect flow-cell statistics and write a tile dump. |
| `tiledump.sh` | `TileDump` | Load, combine spatial bins and analyze an existing tile dump. |
| `plothist.sh` | `PlotHist` | Write histograms for columns in a tile dump or another headered numeric matrix. |
| `plotreadposition.sh` | `PlotReadPosition` | Write read coordinates and barcode Hamming-distance information. |
| `cg2illumina.sh` | `BGI2Illumina` | Convert supported BGI/Complete Genomics read headers to Illumina-style headers. |

For example, the following commands show the input/output relationships. Replace
filenames and choose filtering parameters for the experiment; these examples are
not a calibration or a recommendation to use the defaults for every run.

```bash
filterbytile.sh in=reads.fq.gz out=filtered.fq.gz dump=tiles.tsv
plotflowcell.sh in=reads.fq.gz dump=tiles.tsv
plothist.sh in=tiles.tsv bins=1000
plotreadposition.sh in=reads.fq.gz out=positions.tsv expected=barcodes.txt
cg2illumina.sh in=bgi.fq.gz out=illumina.fq.gz
```

Use separate output paths when comparing commands. FASTQ reads must retain the
header fields needed to identify their lane, tile and coordinates. Quality-based
flow-cell analysis requires per-base quality scores for nonempty primary reads.
Missing qualities produce a clear fatal diagnostic before quality processing,
including when only the second mate lacks them. The separate alignment-only
`sam=` input to `filterbytile.sh` may omit qualities; it supplies alignment
observations rather than primary quality measurements.

`AnalyzeFlowCell` combines spatial quality measurements with optional sampled
Bloom-filter counts, barcode statistics and alignment-derived measurements.
A Bloom count is approximate: collisions can raise a count. Its sampled counts
are distinct from complete read or base totals. `PlotFlowCell` is a separate
implementation with its own options and counting behavior; do not assume it is
an interchangeable alias for `AnalyzeFlowCell`.

## Spatial model

The hierarchy is `FlowCell` → `Lane` → `Tile` → `MicroTile`. A sequencing tile
contains multiple rectangular microtiles. `Tile.xSize` and `Tile.ySize` configure
the grid dimensions; they are shared settings and must remain stable while a
populated grid is used. Loading a tile dump can change these settings.

- `FlowCell` owns lanes, aggregate statistics and operations that widen spatial
  bins to collect more observations.
- `Lane` owns sequencing tiles and lane-level counters.
- `Tile` indexes microtiles by coordinates and provides iteration over populated
  entries.
- `MicroTile` accumulates read, base, quality, alignment, barcode and k-mer
  statistics. These counters have different denominators; consult each method
  or dump column rather than treating all counts as read counts.
- `CycleTracker` collects A/C/G/T/N counts and numeric quality sums by cycle.
  Call `process()` after adding observations before reading derived averages.
  `clear()` retains allocated cycle storage while clearing observations/results.
- `FlowcellCoordinate` is mutable coordinate/parser scratch. Its comparison
  order is lane, tile, y, x. `getFC()` returns a reusable instance local to the
  calling thread; it does not make an instance safe to share between threads.

Model mutation is not generally thread-safe. Some merge paths use locks and some
counter arrays are atomic, but container growth, iteration, multi-field snapshots
and shared configuration still require ownership or coordination by the caller.
A lookup with `create=false` can still create containing lane/tile objects;
it suppresses creation of the requested microtile, not every enclosing container.

## Header parsers

`ReadHeaderParser` defines the common accessors and barcode/coordinate helpers.
The concrete variants retain their input and mutable parsing state; use one
instance per worker, or provide external synchronization.

| Class | Mechanism and scope |
| --- | --- |
| `IlluminaHeaderParser1` | Legacy backward coordinate scan, with separately configured coordinate/comment parsing. Malformed enabled fields use its fatal error path. |
| `IlluminaHeaderParser2` | Indexed delimiter parsing of the seven coordinate fields and optional comment terms. Numeric values are parsed by accessors. |
| `BGIHeaderParser` | Delimiter-based parsing of its supported BGI header layout. |
| `BGIHeaderParser2` | Reverse delimiter parsing used by the BGI conversion tool. |

The variants are retained intentionally. Header support, optional fields and
configuration effects differ; the abstract interface is not a guarantee that
every variant accepts every platform's header. For Illumina comments, `Y` means
chastity-filter failure and `N` means pass. Comment accessors require a valid
comment boundary. `looksValid()` is a layout heuristic for shrinking headers,
not complete validation of all fields or numeric values.

`IlluminaHeaderParser2.encodeCoordinates()` packs lane/tile/x/y into 7/17/20/20
bits. Callers must supply nonnegative values fitting those widths. The method
does not validate that precondition.

## Tile dumps and derived values

`TileDump` handles versioned tabular dumps and applies spatial filtering criteria.
Version 3 remains the default output format. It and earlier formats store some
rounded derived metrics instead of all original raw counters. Reloading reconstructs
synthetic counters and can change ratios through rounding and truncation.

`dumpversion=4` selects an extension that preserves scalar microtile counters:
hits, misses, summed depth and the four quality/probability sums that otherwise
lose precision. It retains the legacy columns in their original order, fills
the seven merge columns even for unmerged tiles, then appends seven `raw_` columns.
Integer counters remain exact; finite floating-point sums round-trip through
Java's double representation. Readers automatically recognize versions 1–4.

`filterbytile.sh`, `plotflowcell.sh` and `tiledump.sh` accept `dumpversion=3|4`.
Use `dumpversion=3` for existing external visualizers that require the legacy
layout; omitting the option produces the same format. `plothist.sh` accepts both
and produces additional histograms for the appended version 4 columns.

Version 4 is an exact scalar-counter format, not a complete flow-cell checkpoint:
legacy summary headers and fitted coefficients retain their existing representation,
and cycle arrays are not serialized. Writing version 4 from an older dump cannot
recover precision already lost by that older file. Keep the input reads when
full recomputation is required.

Read/base counts, alignment counts, error probabilities, percentages and Phred
scores are different quantities. In particular, methods with historical names
containing “percent” need their actual return contract checked before unit
conversion. Some optional report columns depend on collected observations;
consumers should check the dump's version and header rather than assume one layout
across versions. The built-in readers expect the version-specific column order;
header names do not remap columns, so do not reorder them.

The classes' Javadocs describe their individual mutation, ownership and input
contracts. This package overview makes no general speed or cross-platform
compatibility claim; performance and supported input layouts depend on the
selected tool, options and runtime.
