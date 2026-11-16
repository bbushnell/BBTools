# Consensus package

This package collects aligned nucleotide observations, calls reference-based
consensus sequences, and uses paired-read evidence to resize gaps or join contigs.
It also contains the standalone contig-tip polishing tools.

## Entry points

| Launcher | Class | Purpose |
|---|---|---|
| `consensus.sh` | `ConsensusMaker` | Update a reference from aligned reads; optionally write or score graph models. |
| `fixgaps.sh` | `FixScaffoldGaps` | Estimate supported internal N-gap lengths from paired-read insert distributions. |
| `lilypad.sh` | `Lilypad` | Join contigs using cross-contig paired-read evidence and strand/weight thresholds. |
| `tippolish.sh` | `TipPolish` | Run conservative, read-supported contig-tip polishing through child tools. |
| `dev/tipcontextpanel.sh` | `TipContextPanel` | Build and combine candidate tip contexts. |
| `dev/tipcontextguard.sh` | `TipContextGuard` | Check proposed substitutions against read/context and original-kmer evidence. |

Use each launcher's help for its options. The tip tools' evidence rules and
limitations are described in [TipContextGuide.txt](../../docs/guides/TipContextGuide.txt).
The tip algorithms and their existing scientific thresholds are retained.

## Graph data and ownership

`BaseGraph` stores a reference node and deletion node at every padded reference
position. Reference nodes may have insertion chains. The graph constructor adds
one reference observation, including an ambiguous reference symbol; standalone
`BaseNode` construction adds none. Deletion nodes have no A/C/G/T allele arrays.

`BaseGraphPart` supplies the stable REF/INS/DEL type and serialization support.
`BaseNode` stores integer counts and weights, optional probability caches, and
edge references. These integer accumulators do not saturate. Observation weights
may include base quality, mapping quality or the optional identity adjustment;
they are not necessarily plain read counts.

Concurrent `BaseGraph.add` calls synchronize node updates and aggregate counters.
Finish all accumulation before traversing, scoring, initializing probability
caches, or serializing. `BaseNode` and `BaseGraphHelper` do not synchronize their
own operations. `ConsensusObject` policies are mutable process-wide settings;
configure them before processing and coordinate temporary overrides.

With zero padding, `BaseGraph.original` borrows the supplied reference array.
Positive padding creates a new array. Adding or scoring a read may expand its
match string. Orientation search during alignment may reverse-complement the
read and change its strand; callers must own the mutable read and aligner.
The existing terminal-deletion rejection remains part of graph accumulation's
supported match-string boundary.

## Consensus and model scores

`BaseNode.consensus` selects a base using weights, count tie-breaking and the
configured depth/frequency thresholds. It returns the selected base, also writes
that base to `scratch[0]`, and writes numeric confidence to `scratch[1]`.
The confidence is not an ASCII FASTQ character. Called-allele confidence uses selected-versus-
competing weighted odds, clamped to 2–41; protected/low-support/N fallbacks remain
20/2/0. These heuristic scores are not independently calibrated error probabilities.

The default `BaseGraph.qualityByMS=true` derives output qualities from
`(largest count - half the second-largest count) / concrete-base total`.
The alternate mode also returns numeric qualities and retains the existing
insertion-versus-reference/deletion weight denominator. FASTQ writers perform
ASCII encoding. Default-mode output is not intentionally changed by that repair.

Traversal resets reference/substitution/insertion/deletion counters before each
call. Those counters describe the untrimmed traversal. Enabling depth trimming
or N trimming retains the legacy minimum depth of one; uppercase N removal is
additionally controlled by `trimNs`. Depth trimming alone does not imply N removal.

For scoring, initialize probability caches after loading counts, normally through
`BaseGraphHelper.initForScoring`. The helper also imports/exports the five counts
A/C/G/T/deletion for nucleotide model codecs. It does not update insertion chains,
graph-level read totals or an already initialized probability cache.

`score` and `scoreOld` retain distinct indel bounds. Relative scoring uses
`(score-minimum)/(maximum-minimum)`; an undefined or non-finite normalization throws
`ArithmeticException` instead of returning NaN. Raw scoring remains available when
normalization is not requested. Serialization identifiers and field order remain
unchanged.

## Gap and scaffolding behavior

`ConsensusMaker`, `FixScaffoldGaps` and `Lilypad` load all reference records even
when `reads=` limits alignment input. A limited input record may name any contig
in the reference.

Gap resizing uses the global insert distribution and local spanning-pair depth.
Short, terminal, edge-proximal and insufficiently supported N runs are preserved.
`mindepth` must be positive; zero is rejected during argument parsing. Reported changes
use the emitted gap length after the minimum-length clamp, not the requested
pre-clamp adjustment. Empty accepted input has finite summary statistics.

The two scaffolding tools use `Tools.makeHistogram`. It samples `buckets+1`
interior quantiles of a stable, nonnegative count distribution, returning actual
observed values. Rank arithmetic is exact above double's integer precision limit.
An empty distribution retains the zero-filled convention; negative counts and
unrepresentable totals are rejected. The atomic array must not change during
this two-pass operation.

Long-insert reads are preferentially represented among gap-spanning pairs. The
gap estimator uses nearby coverage to choose an insert-distribution percentile;
that adjustment and Lilypad's strand/weight thresholds remain heuristic rules.

## Regression tools

Regression classes live under `test.consensus`; the old tip-test launcher names
are retained after moving their Java classes out of the production package.
The tip tests access existing package-private algorithms through a test-local
reflection facade rather than widening production APIs.

- `dev/consensusgraphboundarytest.sh`: graph text, finite normalization, trimming,
  numeric qualities, model counts and large combined weights.
- `dev/consensuscliboundarytest.sh`: literal CLI inputs and outputs for gap,
  reference-loading, empty-input and failure behavior.
- `dev/consensushistogramtest.sh`: sparse/mixed distributions, exact large-count
  ranks, invalid counts and overflow.
- `dev/consensusreferencelimittest.sh`: independent reference loading under limited
  and unlimited alignment quotas for all three tools.
- `dev/lilypadmixedpairedgetest.sh`: exact edge counts, weights and distances with
  sequential and four-thread accumulation.
- `dev/basegraphdeletionqualitytest.sh`: deletion-flank weights and the supported
  match-string geometry.
- `dev/tipcontextguardtest.sh`, `dev/tipcontextpaneltest.sh`,
  `dev/tipcontextpanelboundstest.sh`, `dev/tippolishtest.sh`: retained tip tests.

Use `phase=fixed` for graph, CLI, histogram and reference-limit fixtures; CLI and
reference-limit fixtures also require `outdir=<new directory>`. `phase=baseline`
is only for reproducing their recorded historical defects. Synthetic score,
traversal and mixed-pair benchmarks measure their named in-memory work, not overall
mapper or assembly throughput. The mixed-pair timing includes a method-handle call
and measures uncontended updates. See the maintenance records for exact source
revisions and results.
