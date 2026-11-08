# Bloom filters, counting arrays and read correction

This package contains fixed-capacity, packed counting structures and tools that
use their counts for reference matching, depth filtering and error correction.
Hashed counts are approximate: different keys can share cells, and finite-width
cells saturate. An apparent hit is not proof of an exact sequence match.

The numbered classes preserve different implementations and experiments. Their
numbers are version names, not k-mer lengths: `KmerCount3` does not mean `k=3`.

## Command-line tools

Run the launchers from the BBTools root or put that directory on your path.
They select the classpath, assertions and JVM options. Arguments use `name=value`;
run a launcher without arguments for its detailed help.

| Launcher | Main responsibility |
|---|---|
| `bloomfilter.sh` | Route reads by consecutive reference-k-mer hits. |
| `bbcms.sh` | Correct reads and/or filter by k-mer depth; optionally run two passes. |
| `polyfilter.sh` | Trim and filter suspicious homopolymers using sequence, depth and quality criteria. |
| `bloomfilterparser.sh` | Summarize the historical `BloomFilterWrapper` text-log format. |

For example, reference matching uses explicit output names:

```bash
bloomfilter.sh in=reads.fq ref=reference.fa outm=matched.fq outu=unmatched.fq k=31 minhits=3
```

`outm` receives matches; `outu` receives nonmatches. **`out` is an alias for
`outu`, not `outm`.** The corresponding second-mate outputs are `outm2` and
`outu2` (`out2` aliases `outu2`). By default a pair matches when either mate
matches; `requireboth=t` requires both. `minhits` is a consecutive-hit criterion,
so a read can share some reference k-mers without qualifying as a match.
Reference loading can also exclude windows according to its quality settings.

Two-pass correction on separate mate files can be requested with:

```bash
bbcms.sh in=reads_1.fq in2=reads_2.fq out=corrected_1.fq out2=corrected_2.fq passes=2 k=31 bits=4
```

The first pass builds a cheap filter and removes low-support/junk reads; the
second builds its filter from the survivors and applies its own settings. This
changes the count pool and is not equivalent to a single correction pass.
The intermediate file contains both mates when the input is paired. The wrapper
carries the actual first reader's pairing state into the second pass instead of
reusing flags that described the original separate input files.

`bbcms.sh` also exposes a separate `fixindels` algorithm. It has stricter input
and option constraints than ordinary substitution correction; use its launcher
help rather than assuming the modes are interchangeable. `BloomFilterCorrector1`
handles primitive packed keys; `BloomFilterCorrector2` handles multiword keys.
The launcher's `packed` setting affects the representation of long k-mers.

For a homopolymer-filtering example:

```bash
polyfilter.sh in=reads.fq out=clean.fq outb=bad.fq k=31 polymers=G
```

Here `outb` means rejected reads. The filtering/trimming thresholds and optional
count-pool inputs are specific to PolyFilter; they are not the reference-matching
output conventions of `bloomfilter.sh`.

## Source map

| Classes | Role and limits |
|---|---|
| [BloomFilter](BloomFilter.java) | Loads a counting table and implements query/count-profile helpers. Supports reference files and serialized filters. |
| [ReadCounter](ReadCounter.java) | Builds k-mer count arrays from read input, with worker-local counting state. It is not a general read-statistics reporting class. |
| [KmerCountAbstract](KmerCountAbstract.java) | Shared counting configuration, statistics and histogram helpers; settings are process-wide. |
| [KCountArray](KCountArray.java) | Common counting-array API and factory. The current factory selects `KCountArray7MTA`. |
| [KCountArray7MTA](KCountArray7MTA.java) | Atomic packed storage, optional prefiltering and one-/two-hash-key operations. Compound counting and accounting contracts depend on the selected mode. |
| [KCountArray3](KCountArray3.java) | Direct-address packed array implementing the common `KCountArray` API. |
| `KCountArray4`, `4MT`, `5MT`, `6MT`, `7MT`, `8MT` | Alternative/historical hashed implementations. Queued variants have their own startup, shutdown and query rules. |
| [KCountArray2](KCountArray2.java) | Legacy direct-address packed array with a different API; not the factory's current hashed implementation. |
| [IndexCounter](IndexCounter.java) | Builds a presence table from the indexed reference, rather than a read file. |
| [BloomFilterWrapper](BloomFilterWrapper.java) | Reference-match routing, optional count output and overlap handling. |
| [BloomFilterCorrectorWrapper](BloomFilterCorrectorWrapper.java) | Input/output and worker orchestration for correction and filtering. |
| [BloomFilterCorrectorMultipass](BloomFilterCorrectorMultipass.java) | Coordinates one or two wrapper invocations in the same JVM. |
| `BloomFilterCorrector`, `BloomFilterCorrector1`, `BloomFilterCorrector2` | Shared correction logic plus short-/long-key representations and scratch state. |
| [ErrorCorrect](ErrorCorrect.java) | Older correction and trusted-window helpers, also used by legacy counting paths. Not the BBCMS launcher entry point. |
| `KmerCount3`, `4`, `5`, `6` | Legacy counting programs/helpers with positional interfaces and implementation-specific quality, strand and split-window behavior. |
| `LargeKmerCount`, `LargeKmerCount2` | Historical rolling-hash/collision-estimation experiments. These are not the modern long-k BBCMS implementation. |
| [PolyFilter](PolyFilter.java) | Homopolymer trimming and filtering. |
| [ParseBloomFilter](ParseBloomFilter.java) | Extracts selected log prefixes into a compact table; not a general log parser. |

## Contracts for Java callers

Use the method signatures and Javadocs of the concrete class. There is no
`BloomFilter.contains(String)` convenience API. Short-key methods accept packed
numeric keys; long-key paths use `ukmer.Kmer` and the corresponding hash mode.
Do not pass a long-key hash to an API expecting a packed nucleotide window.

Several `BloomFilter` rolling-query methods require `k<=31`. A separate `kbig`
profile can combine shorter-window counts by taking their minimum; that is a
count estimate for the larger window, not an exact long-key table lookup.
`fillCounts` and `fillCountsBig` reuse caller-owned lists. Output positions include
undefined windows with zero counts, while the returned valid-window count excludes
them. A sequence too short for the requested window leaves an empty list.

`passes` returns true for a read that does **not** meet the configured matching
criterion. `matches` inverts that result. Overloads taking reads and a `LongList`
recompute keys into scratch storage; they do not consume caller-precomputed keys.
The `passes(LongList, int)` overload instead consumes precomputed keys without
changing the list; negative keys separate runs of defined windows.

`hasHighCountFraction` excludes undefined windows from its denominator and rejects
reads with no defined windows, even at fraction zero. It retains the rounded count
threshold: at least `round(fraction * definedWindows)` must meet the depth threshold.
`averageCount` suppresses isolated depth spikes by clipping interior values to the
larger neighboring depth and endpoints to their sole neighbor. A singleton keeps
its depth. Defined windows on opposite sides of an undefined run still share the
existing smoothing state.

`ReadCounter.makeKca` applies its read limit to the primary input. Extra files are
read completely in every counting pass, including intermediate trusted-table passes.

Counts have representation-specific limits. Common packed widths of 1, 2, 4 and
8 bits hold at most 1, 3, 15 and 255, respectively. Saturated counts cannot report
higher depth. A hashed table is not an exact key/value map, and an atomic backing
array does not make every multi-cell operation or derived statistic exact.
The four-neighbor sums returned by `KCountArray7MTA.readLeft` and `readRight`
are capped at `Integer.MAX_VALUE` to preserve their `int` return contract.

In the modern dual-hash path, callers must use the same primary/secondary key
construction for insertion and lookup. This is different from the legacy
`LargeKmerCount` experiment's secondary probe in the same table.

Read the concrete concurrency contract before sharing a counter. Some legacy
MT implementations queue updates and require shutdown before ordinary reads;
`KCountArray7MTA` has a different lifecycle. Many options, seeds and counters are
static. Separate instances do not isolate those settings when several invocations
share a JVM. Corrector scratch buffers belong to the caller/thread and must not
be used concurrently by another worker.
Driver failure contracts also differ: some close streams and restore shared
settings only on normal completion. `IndexCounter` propagates worker failures
to the caller; its partially updated table and counters must then be discarded.

## Reproducibility and maintenance

Record the source revision, input selection, k/representation, bit width, hash
mode, cell count, seed, thread count and relevant quality/correction settings.
Automatic memory sizing can choose different capacities between environments;
use an explicit `cells` value for comparisons that need a fixed table size.
Do not infer a speedup from source shape or compare runs with different outputs
as if they performed the same work.

Regression helpers, including the relocated legacy harnesses, live under
`current/test/bloom`. Existing launchable checks
include `dev/kcountdualhashtest.sh` and `dev/readcounterfailuretest.sh`.
These checks exercise particular contracts, not a proof of the whole package.
New failures should retain a finite input, independent expected result and both
before/after outputs. Documentation-only changes should preserve executable
behavior, including field-initializer order.
