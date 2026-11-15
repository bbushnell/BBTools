# PacBio subread processing

This package contains the PacBio well-based filtering and simulation tools below.
Run them through the BBTools shell launchers, which configure the JVM and enable
assertions. Arguments use `key=value` syntax; input/output paths use `in=` and
`out=` or the tool's named output categories.

| Launcher | Role |
| --- | --- |
| `icecreamfinder.sh` | Classifies inverted-repeat and adapter evidence, with optional trimming, entropy filtering, good/ambiguous/bad outputs, and junction output. |
| `reformatpb.sh` | Rewrites, filters and samples PacBio subreads by well; also has a consensus-building option. |
| `icecreammaker.sh` | Simulates subreads, base-call errors, missing/hidden adapters and optional reference inverted repeats. Requires a positive `zmws=` count. |
| `icecreamgrader.sh` | Counts good/bad reads and bases from the simulator's truth headers. It does not classify contamination from sequence. |

Keep subreads from each well consecutive. `ZMWStreamer` supplies groups to worker
threads; it does not sort arbitrarily shuffled reads. `ZMW` retains the original
mutable `Read` objects, and its trimming/discard helpers can change them. Input
read limits and emitted-well limits have different meanings. The optional SAM
backend drains input after reaching its output-well cap; that cap does not promise
that physical input reading stops early.

## Simulation and truth metadata

`ReadBuilder` writes the positional tab-separated truth fields consumed by the
grader. Ordinary PacBio headers are not a substitute for that format. Header
error rates and identity histograms describe requested simulation parameters,
not a recount of realized substitutions and indels. The simulator does not
compute base qualities.

The simulator's `ccs=t` option selects one internal median-length simulated
subread. It does not perform consensus calling, and it still uses the configured
error model. This differs from ReformatPacBio's consensus option.

## Helpers and alignment

`PBHeader` parses PacBio movie/well/coordinate headers. Coordinate intervals are
half-open; parsing does not establish agreement with a read's sequence length.
Its coordinate parser and numeric-well-only parser accept different suffixes.
`PolymerTrimmer` returns end spans using an exact-byte running-score heuristic;
its fraction parameter is not a separate measured-fraction test. Configure its
shared static settings before starting workers.

`IceCreamAligner.makeAligner` selects JNI when `Shared.USE_JNI` is enabled;
otherwise it requires `bits=32` and returns the Java implementation. Both layouts
compare bytes with fixed integer scoring. The full layout uses reference-width
working arrays and a row-score early exit. The short layout uses query-width
arrays and does not use `minScore`; “short” describes that layout, not necessarily
16-bit arithmetic. JNI chooses its arithmetic width internally and requires the
native library. It does not silently fall back to Java after loading failure.

Alignment endpoint and iteration-count conventions differ between layouts and
backends; use the method contracts rather than assuming they are interchangeable
counters. Each worker owns an aligner with mutable counters; alignment working
arrays are allocated per call.

The CLI constructors change shared I/O/settings state. These tools are intended
as separate invocations, not independent concurrent instances in one JVM. Source
TODOs identify remaining legacy limitations; this package map is not a claim that
every input, backend or failure schedule has been certified.
