# Assemble Package - Genome Assembly and Scaffolding Tools

## Overview
The `assemble` package contains tools for genome assembly, scaffolding, and assembly improvement.

## Core Components

### AbstractBuildThread.java — Build-worker state

This base class owns each worker's scratch count arrays, sequence builder, result
list and counters. Concrete Tadpole and KmerCompressor workers implement traversal.
The stream array is borrowed; it may be null for table-seeded assembly.

The constructor installs an uncaught-exception handler that stores a worker failure.
Coordinators must start and join the workers, inspect that failure and propagate it
before accepting results. Capturing an exception does not close streams, cancel
other workers or synchronize graph/table updates.

`readsInT` counts individual reads, including mates. `lowqBasesT` counts all bases
in reads marked discarded, rather than only bases whose quality is low.
`insertSizes` is a histogram; bin zero includes unresolved estimates. The base
class makes no guarantee of nonblocking I/O or linear scaling with thread count.

### AbstractRemoveThread.java — Count-range pruning workers

`process()` prunes a caller-owned table set to the inclusive count range
`[min,max]`. It selects the standard worker for an exact `KmerTableSet` instance
and the unlimited-table worker otherwise; supported inputs must match those
concrete types. A positive worker count and exclusive access to the table set
are caller obligations.

Workers claim shards through one atomic counter, zero out-of-range values in
both hash arrays and victim trees, clear ownership and regenerate each processed
shard. Victim-tree processing visits the whole tree; it is not logarithmic work.
The coordinator joins each worker, retrying interrupted joins, before reading
that worker's counters and failure state. A worker failure is propagated after all joins,
without reporting partial removal as success. Counts or shards already modified
are not rolled back. With no failure, the method returns total removed k-mers
and optionally prints timing/statistics.

### ErrorTracker.java — Per-read correction accounting

A correction worker reuses this mutable object for successive reads and clears it
before reuse. It owns no read or sequence data and provides no synchronization.
Callers update the public counters and rollback diagnostics.

- `detected()` sums pincer, tail, reassembly and substitute detections. It excludes
  the dormant brute-force counter.
- `corrected()` includes all correction counters, including brute-force.
- `correctedReassemble()` sums inner and outer reassembly corrections.
- These totals count algorithm events; they do not deduplicate base positions.
- `clearDetected()` and `clearCorrected()` reset their respective counters.
  `clear()` also resets marking and rollback diagnostics.
- `toString()` reports the active detection/correction and marking counters,
  omitting brute-force counters and rollback diagnostics.

### KmerCompressor.java — Short-kmer sequence packing

`KmerCompressor` loads a short-word count table and packs retained DNA k-mers into
sequence paths for compact FASTA storage. It is not a consensus genome assembler:
extension chooses the first eligible unowned neighbor in A/C/G/T order, rather
than selecting the highest-supported branch or stopping at every junction.
The short-word table restricts effective K to 1–31.

Count limits restrict eligible seeds and extensions. Optional same-count mode
also requires equal counts along an extension. Traversal proceeds in both
sequence orientations, using reverse-only keys for the second direction when
strand canonicalization is disabled. Invalid seeds, unavailable neighbors,
ownership conflicts and the extension length cap can stop a path.

Workers keep their own sequence buffers and results while sharing table ownership.
Seeds and appended neighbors require an unowned k-mer; a higher worker ID cannot
replace another worker's claim. Primary cells use compare-and-set and overflow
nodes use their own monitors. An extension may confirm its already-owned starting
seed, but may not append an already-owned neighbor. These checks require stable
tables and do not guarantee arbitrary concurrent access or deterministic scheduling.
Captured worker failures abort before output; shared configuration requires serial
program instances.

Completed paths are sorted longest-first and can be packed together with configured
N padding before FASTA output. Such packing is a storage operation, not evidence
that the paths are adjacent in a genome. Writer errors are retained for the outer
process error check. Use Tadpole for assembly policies rather than interpreting
these greedy storage paths as resolved biological contigs.

### Reassemble.java — Separate genome assemblies in one output

`reassemble.sh` runs Tadpole separately for each input genome and combines the
contigs into one output. It accepts comma-separated input paths or one level of
files from an input directory. Positive taxonomy labels parsed from filenames or
headers are passed through; they are not checked against a taxonomy tree.

Default append mode clears a replaceable output once before processing genomes.
With `tempdir=DIR`, per-genome outputs are concatenated after assembly instead.
Input aliases and overwrite protection are checked before output initialization.
Failed genomes produce a final error even with `failfast=f`; output may be partial.
Requested temporary-file deletion occurs only after all assemblies and final output
completion succeed.

The next genome's numeric ID offset reserves the preceding retained contig-ID
range, including contigs removed by final filtering. IDs can therefore contain
gaps. Reported contig/base counts describe emitted records, not that reserved
range. Empty assemblies preserve the incoming offset. Trials run serially in one
JVM and share parser/Tadpole settings; concurrent instances are unsupported.

### TadpoleWrapper.java — Serial K search and final-output selection

`tadwrapper.sh` runs Tadpole trials at requested K values, rounding them through
`Kmer.getKbig()`, removing duplicates and trying the initial values in ascending
order. An omitted K list defaults to 31. Optional boundary expansion and neighbor
bisection add trials; the search is heuristic and does not guarantee a global
optimum or biologically correct assembly.

Each trial uses Tadpole's normal policy/rescue dispatch. The wrapper computes
statistics from that trial's final named assembly, even when `showstats=f` hides
reporting. Assemblies with contigs rank ahead of those with none, even when the
contigs are shorter than 64 bases. Selection rounds L50, L90 and maximum length
down to multiples of 64 bases, then compares those values in that order, preferring
larger values. Ties prefer fewer contigs, then smaller K. Raw length differences
within a bin do not break ties. Fixed bins give a consistent ordering; a one-base
difference across a bin boundary can still matter. This ranking is not a guarantee
of assembly quality.

`out=` supplies the trial-name pattern containing `%`; `outfinal=` optionally
publishes the selected assembly. Publication checks input aliases and replacement
policy and retains the selected source if publication fails. With `delete=t`,
nonselected intermediates are removed on a best-effort basis after selection and
any requested publication. Other options are forwarded to Tadpole.

Trials share one JVM and static configuration; the wrapper is not reentrant.
Its explicit garbage-collection requests do not establish a memory-use bound.
Use the current shell usage for options and the source contracts for failure and
cleanup details.

### Contig.java — Mutable sequence and graph state

A contig combines a retained mutable base array with coverage, endpoint state,
edge lists and graph-quality annotations. Constructors wrap the supplied array
without copying it. Callers editing bases must keep cached GC/HH/CAGA values
coherent and coordinate exclusive access to graph mutations.

The mean assembly-kmer depth is a float; minimum and maximum depths are integers.
Generated headers format the mean to one decimal and may include sequence scalars,
taxonomy and graph classes. An explicit name bypasses generated metrics but may
still gain graph-class text. Header detail depends on `SHORT_NAMES` and
`VERY_SHORT_NAMES`; the default short form omits endpoint codes and edge lists.
Endpoint codes originate in ShaveObject and are inherited by Tadpole.

`toFasta()` appends the header and optionally wrapped sequence to a supplied
ByteBuilder. The mutable contig ID connects its edge references to the surrounding
graph. Used/associate flags are ordinary booleans. Their setters remove local and
supplied inbound edges; callers must provide the complete inbound list and
coordinate flips, merges and renumbering.

### Edge.java — Directed contig connection

An edge stores source/destination IDs, endpoint orientations, supporting depth,
optional emitted sequence and overlap/trim metadata. Orientation bit 0 selects
the source's right side and bit 1 the destination's right side. Values 0–3
therefore represent left/left, right/left, left/right and right/right endpoints.

`length` means emitted bases for dense paths, k-mer steps for graph-only paths,
or zero for cross-k overlap links. The separate `overlap` field stores an explicit
overlap; zero means use the active graph k. `sourceTrim` and `destTrim` describe
bases replaced or omitted at cross-k joins.

The optional payload array is retained by reference. `flipSource()` complements
and reverses it in place and toggles the source bit; `flipDest()` toggles only the
destination bit. `traversalBase()` reads left-source payloads backward with
complementation. Callers must coordinate mutations of edges and shared payloads.

`merge()` requires matching endpoint IDs and orientation. It adds supporting
depths and replaces payload, length, overlap and trims only when the incoming
depth exceeds this edge's previous depth. It retains the chosen payload without
copying it and leaves path-minimum/maximum/mean depth fields unchanged.
Debug text and Graphviz DOT output describe the stored connection; they do not
validate its biological correctness.

### AbstractExploreThread.java — Exploration workers

Each instance owns a backing thread, two distinct k-mer workspaces, reusable
extension-count arrays and a path builder. The concrete Shaver1/Shaver2 workers
claim primary table partitions through `nextTable`, then overflow-forest
partitions through `nextVictims`. Both scans examine candidate k-mers; the
victim scan is not a second validation pass over confirmed dead ends.

Workers collect local selection counters and matrices. The reported dead-end
counter can include selected bubble paths, and selection counts do not by
themselves measure distinct k-mers ultimately deleted. Historical comments about
nondeterministic counts remain investigation notes, not statistical guarantees.

Uncaught failures are saved for the coordinating Shaver, which joins every worker
before checking failures and aggregating successful phase results. Private
scratch buffers reduce shared state; safe graph mutation still depends on the
table ownership protocol and the coordinator's phase boundary.

### AbstractShaveThread.java — Removal workers

This base thread repeatedly calls the concrete `processNextTable()` until no
partition remains. It does not poll a shared work queue. Shaver1/Shaver2 use an
atomic partition counter so one removal worker handles a table at a time.

For each claimed table, the concrete worker zeros counts marked `STATUS_REMOVE`
in both primary cells and overflow nodes, clears ownership, and calls
`regenerate(0)`. Its removal counter records the regeneration result.
Uncaught failures are retained for Shaver to propagate after all phase workers
join; already-applied table mutations are not rolled back.

### BubblePopper.java — Serial contig-graph simplification

BubblePopper retains the live contig list, destination-edge map and an optional
error classifier. It mutates those graph structures through direct linear merges,
indirect error-arm pruning and biological-bubble unzipping. Callers coordinate
serial access; ordinary used/associate booleans are not synchronization.

`expand()` works rightward and can flip the current contig to process its other
side. `popDirect` defaults true, `popIndirect` false and `unzipBubbles` true.
Indirect pruning and unzipping also require an error classifier; a null classifier
leaves direct merges available. A separate optional debranching pass handles
short terminal branches.

The representative-arm heuristic uses `minLen=2*kbig-1`: keep the first arm
until an eligible-length arm is found, then prefer the highest entry depth among
eligible arms, breaking ties by contig length. If every arm is short, the first
remains selected and the caller rejects it. The source retains a warning that
this preference can favor a long, low-depth error arm over a shorter real arm;
it is not a global optimality rule.

Merge and retirement helpers update edge references as sequences are rewritten.
Used/associate setters remove local edges and every supplied inbound reference;
the caller must supply the complete inbound lists. Optional graph validation
checks structural relationships, rather than proving assembly accuracy.

`expansions` counts accepted merge, prune and unzip operations.
`contigsAbsorbed` counts contigs retired or associated with a surviving path;
`trueBubblesUnzipped` counts unzipped biological bubbles, and `branchesRemoved`
counts removed dead-end edges. These operation counters are distinct from final
assembly statistics.

### AbstractProcessContigThread.java — Contig-end workers

Workers borrow a stable contig list and share an atomic next-index counter.
Each claimed contig goes through its left-end hook and then its right-end hook.
The list must not be resized or have entries replaced while these workers run.
Atomic index assignment distributes work; it does not synchronize graph mutations
inside the concrete hooks.

Count arrays, the path builder, traversal results and statistics are worker-local.
The base does not clear scratch before each hook. The concrete implementation
establishes the values it needs, including the meaning and validity of traversal
results after an unsuccessful search.

The constructor installs an uncaught-failure handler. Tadpole's coordinator joins
workers, checks failures and aggregates their counters. Calling
`processContigs()` directly instead propagates an exception to that caller.
This class does not own resource closure or sibling cancellation.

### Shaver.java — Exploration/removal coordinator

`makeShaver()` selects Shaver1 for an exact `KmerTableSet` instance and Shaver2
for an exact `KmerTableSetU` instance; other runtime classes, including subclasses,
are rejected. The short form defaults to count bounds/seed/extension minimum of
one, branch multiplier three, exploration and removal limits of 100, and both
hair and bubble removal enabled.

A shaving pass initializes ownership and runs two sequential worker phases:
exploration marks candidates, then removal deletes marked k-mers. Each phase
distributes table partitions with atomic counters, joins every worker, checks
captured failures, and only then aggregates that phase's statistics. An
exploration failure prevents the removal phase from starting. A failure in either
phase is rethrown with its cause after the joins; mutations already made by other
workers are not rolled back.

Calls on the same coordinator and tables must be serialized, with the table's
ownership-initialization preconditions satisfied. Atomic partition counters do
not make the entire operation lock-free or safe for unrelated concurrent writers.

`shave(minCount, maxCount)` changes those two thresholds before running a pass.
Assertions require `minCount <= minSeed <= maxCount` and at least one removal
mode. The return value is the current pass's removed-k-mer count; public totals
accumulate completed phase results across calls. Candidate-path counts and
optional event matrices are diagnostics, distinct from the actual removal count.

`startFromHighCounts` and `shaveFast` default true. `shaveVFast` is a compile-time
false option, not an enabled acceleration path. Search/removal thresholds control
which paths can be considered; these heuristics do not guarantee better contiguity
or biological correctness for every input.

### Tadpole1.java — Single-word count and traversal implementation

`Tadpole1` supplies the single-long k-mer implementation used by normal Tadpole
dispatch for effective K up to 31. It owns a `KmerTableSet` and inherits argument
processing, pipeline lifecycle and shared correction policies from `Tadpole`.
The constructor passes a per-entry sizing estimate to the table loader; that
estimate is not a measured whole-program memory requirement.

Workers distribute table and victim-chain work, maintain their own traversal
scratch and counters, and share count/ownership state. Ownership supports the
assembler's coordinated traversal; it does not permit arbitrary concurrent
mutation of tables or reads. Constructing independent pipelines concurrently is
unsupported because configuration also changes process-wide settings.

This implementation supplies count queries, contig extension, graph-end analysis,
insert-distance estimation and legacy correction kernels. Count lookups preserve
the underlying negative absent-key sentinel. Shared reassembly logic calls the
concrete kernel in each sequence orientation; qualities must use the same base
coordinates. An owned extension that reaches its length bound returns `TOO_LONG`
rather than an ownership failure. Reassembly compares the current base's numeric
index with the best-supported base index, rather than comparing an index with
support depth. Other source TODOs retain their separate investigation scope.

### Postfilter.java — Mapping followed by coverage filtering

`postfilter.sh` maps reads to an assembly with BBMap, writes coverage statistics,
unloads the mapping reference through `Data.unloadAll()`, then calls
`FilterByCoverage` on the original contigs. Coverage thresholds select supported
contigs; passing them does not establish that an assembly is biologically correct.

The wrapper defaults to `mincov=2`, `minreads=6`, `minpercent=95` and
`minlength=400`. A contig without a coverage entry is rejected. Optional `trim=N`
trims both ends before the length check; `trim=t` means 100 bases per end and
`trim=f` means zero. Coverage statistics come from mapping before trimming.

Required inputs are `in=` and `ref=`; `in2=` supplies a separate mate file.
When `in2` is absent and a literal `in` path containing `#` does not exist,
the wrapper expands that placeholder to `1` and `2`. Defaults are
`out=filtered.fa`, `covstats=covstats.txt`, `overwrite=t` and `append=f`.
`outdirty=` optionally retains rejected contigs. Path-collision and output
overwrite checks run before mapping, including the rejected-output path.

Unrecognized arguments are forwarded only to BBMap. The wrapper then appends its
own mapping options, including `ambig=all`, `nodisk`, `maxindel=0`, `minhits=2`,
`bw=20`, `tipsearch=0` and `rescue=f`. It directly parses overrides for
`maxindel`, `minhits` and `rescue`; `bw` and `tipsearch` are fixed wrapper fields.
`nodisk` concerns the mapping index, not the coverage and sequence output files.
Filtering receives the wrapper's selected thresholds and paths, rather than all
forwarded mapping arguments. These in-process tool calls use shared global state;
the wrapper does not provide an isolated concurrent pipeline instance.

With `showstats=t` (the default), a non-stdio output is reread by AssemblyStats2
after filtering. Set output paths and overwrite policy explicitly when rerunning:

```bash
postfilter.sh in=reads.fq ref=contigs.fa out=filtered.fa covstats=coverage.txt ow=f
postfilter.sh in=reads#.fq ref=contigs.fa out=clean.fa outdirty=rejected.fa mincov=5 minr=10 minp=90 minl=500 ow=f
```

### TadPipe.java — Serial assembly pipeline with isolated child JVMs

TadPipe expects pretrimmed, filtered, quality-calibrated paired Illumina reads.
It runs each selected stage through its shell launcher in a fresh JVM using a
`ProcessBuilder` argument list. It does not invoke tool `main()` methods in one
shared JVM or perform the old adapter-trimming workflow described by earlier
versions of this README.

The current order is:

| Stage | Default work |
|---|---|
| Pair preparation | Clumpify deduplication without correction, or Reformat when `dedupe=f`. |
| Optional read correction | BBMerge ECCO, four-pass Clumpify correction, then Tadpole correction at K62. |
| Merging | Initial overlap merge, then configured REM passes (default K124,145,93), followed by concatenation of merged subsets. |
| Optional quality trimming | BBDuk trims the remaining unmerged pairs. |
| Optional extension | Merged and unmerged streams are extended separately, each using the original unextended evidence streams. |
| Assembly | Tadpole receives a configured multi-K schedule; defaults include `assemblek=124` and `k=124,300,96,64,32`. |
| Optional polishing | TipPolish uses the retained uncorrected first-stage reads; disabled by default. |
| Final scaffolding | HollowBridge runs after optional polishing; enabled by default. |

Stage switches and defaults are configured in `Config`. Generic phase overrides
are applied before more specific groups; for example, `merge_`, then
`merge124_`, then the applicable `remN_` group. Final assembly aliases such as
`k=` configure only the assembly child. The coordinator reserves managed paths,
pairing and overwrite options; polishing and scaffolding have restricted override
sets. These checks prevent arbitrary stage output redirection, not every possible
invalid tool-option combination.

Inputs must be readable named files. A single input is treated as interleaved
pairs; explicit `interleaved=f` requires `in2=`. Paths are resolved before
children change their working directory; pipeline paths cannot contain commas or
newlines. The selected launchers and enabled model files are checked before any
stage runs. Custom fusion models require an explicit cutoff.

Each invocation creates a unique work directory containing `commands.sh` and
one combined stdout/stderr log per stage. `dryrun=t` still validates prerequisites
and creates that directory/log, but launches no processing stages. Each child must
exit successfully and produce its declared regular output files before the next
stage starts.

Publication requires a regular final FASTA containing sequence text; the check
does not validate its alphabet or biological correctness. The final file is moved
to `out=`, respecting the default `overwrite=f`. Requested managed-read cleanup
runs only afterward. Failures before publication preserve intermediates; cleanup
itself can fail after publication or after deleting some intermediates. Logs and
optional unpolished/unscaffolded assemblies are retained.

Resource settings are passed explicitly to sequential children: default
`childheap=14g`, initial heap 32MiB, the requested thread count and compression
level. There is no claim that this recipe is optimal for every dataset.
Use `tadpipe.sh` for the current options and `dryrun=t` to inspect the concrete
commands for a particular configuration.

### ShaveObject.java — Shared mode, traversal and branch constants

`ShaveObject` holds shared integer codes and mutable diagnostic settings. It does
not implement an assembly pipeline, logging service or synchronization policy.
Its mode constants cover contig assembly, extension, correction, insert-distance
estimation and discarding; Tadpole's additional graph mode is defined elsewhere.

Traversal outcomes include continuation, dead end, length/depth limits, loop,
success and ownership/seed failures. Branch codes combine `BRANCH_BIT=16` with
forward (`F_BRANCH=17`), backward (`B_BRANCH=18`) or both directions
(`D_BRANCH=19`). `isBranchCode()` tests that bit; it does not validate that an
arbitrary integer is a recognized result code.

The separate `STATUS_*` values represent unexplored, explored, remove and keep
states. `codeStrings` supplies display labels, including placeholders for unused
numeric slots. `MAX_CODE` is its length and therefore an exclusive index bound.
The final array reference does not make its elements immutable. Diagnostic flags
and the stream reference are shared mutable fields rather than per-run state.

### Shaver2.java — Multiword k-mer shaving

This implementation binds Shaver's exploration/removal phases to
`KmerTableSetU`, mutable `Kmer` workspaces and its primary/overflow tables.
It provides multiword k-mer handling, not an unlimited-memory guarantee.

Candidate exploration walks rightward, reverse-complements the accumulated path,
then explores the other direction. It stops at configured distance/count bounds
or graph termination conditions. A path exceeding `maxLengthToDiscard` is not
selected for deletion. With hair removal enabled, a dead-end paired with a
dead-end or backward branch can be selected; with bubble removal enabled, two
backward-branch outcomes can be selected. Selection marks k-mers for the later
removal phase; it does not merge two contig paths.

The protected-rejection variant distinguishes `STATUS_KEEP` from ordinary
explored ownership and checks ownership after a removal claim. Callers must
preserve the table ownership and phase contracts described under Shaver.
Worker-local count arrays/builders/Kmer scratch do not make unrelated concurrent
table mutation safe. Primary and overflow scans use the configured low/high-count
handlers; `shaveVFast` remains compiled false.

### Tadpole2.java — Multiword count and traversal implementation

`Tadpole2` supplies the `KmerTableSetU` implementation selected for larger effective
K values or an explicit force option. Its mutable `Kmer` objects represent the
configured multiword layout; it does not automatically fall back to a short-word
assembler during traversal. Common lifecycle and policy remain in `Tadpole`.

The table set owns counts and shared owner IDs. Workers keep separate Kmer objects,
buffers and counters; instances and process-wide configuration are not reentrant.
The concrete kernels provide graph traversal, read extension, insert-distance
estimation and correction. Quality arrays are read at the currently oriented
base index, so shared callers must keep bases and qualities aligned.

Explicit Hollow options select additional experimental paths with compatibility
checks. Requested original-read seed-protection counts use a separate table and
restore the primary pairing policy after loading. These options are not an
unconditional step in every assembly. See the source contracts and shell usage
for supported combinations, and retain the documented outstanding correction
findings when interpreting results.

### Rollback.java — One-use read-field snapshot

`Rollback` copies a read's bases, optional qualities and optional count profile;
it also saves the string ID and flags. It does not snapshot numeric ID, mate links,
alignment data or every other Read field. The caller must own the read and count
list exclusively during saving and restoration.

Restoration reuses equal-length base/quality arrays when possible and otherwise
transfers the saved arrays to the read. Saved null qualities clear any qualities
added later. A supplied count list regains its saved live size and values, growing
its storage if needed; unused capacity is not cleared. Requesting counts when none
were saved throws before changing the read.

The snapshot is not deeply immutable: final array references do not prevent array
mutation, and `counts0` is publicly accessible for correction comparisons. Do not
modify saved counts before restoration, and discard the snapshot after restoring
because the read may now share its saved arrays. This helper is not a general
transaction over the read or the assembly graph.

### Shaver1.java — Single-word k-mer shaving

This implementation binds Shaver's phases to `KmerTableSet` and single-long
k-mer encoding. Exploration and removal visit both primary table cells and
overflow-forest nodes, using the shared partition counters and ownership states.

`exploreAndMark()` explores rightward, reverse-complements the accumulated path,
then explores in the opposite direction. Distance, count, branch and loop
conditions bound the walk. Paths longer than `maxLengthToDiscard` beyond the
initial k-mer are retained. Enabled hair removal selects a dead-end paired with
a dead-end or backward branch; enabled bubble removal selects two backward-branch
outcomes. These paths are claimed as `STATUS_REMOVE`, not collapsed into one
contig. The later removal workers zero marked counts and regenerate their tables.

The class retains distinct traversal methods and low/high-count seed handlers.
Those alternatives have different implementation histories; their presence alone
does not establish that one is faster or equivalent. See ShaveObject for outcome
codes and Shaver for worker joins, failure propagation and serialized-call
requirements. Neither an atomic partition counter nor thread-local scratch is a
general race-freedom guarantee.

### Tadpole.java — Shared assembly and read-processing framework

`Tadpole` owns common configuration, lifecycle, graph processing and correction
logic. `makeTadpole()` chooses `Tadpole1` for effective K up to 31 and `Tadpole2`
for larger K or an explicit force option. The concrete classes supply their table
representation, count queries and traversal kernels. The main entry point expands
configuration before routing to the assembly policy driver or requested same-K
rescue path.

The supported modes are:

| Mode | Work |
|---|---|
| `contig` | Build contigs from loaded k-mer counts and apply configured graph processing. |
| `extend` | Extend input reads through the count graph. |
| `correct` | Correct input reads using the selected correction path. |
| `discard` | Filter reads using count-based criteria. |
| `insert` | Estimate paired-read insert distances and collect their histogram. |
| `graph` | Explore connections between supplied `contigs=` and write `dot=` or `gfa=`; this mode does not write reassembled `out=` sequences. |

A run owns its table data and mutable contig population. Workers borrow streams,
share counts and keep separate scratch buffers and counters. Construction changes
process-wide parsing, I/O and k-mer settings: independent instances are not
isolated or generally reentrant. Atomic work assignment and k-mer ownership do
not make arbitrary concurrent table or graph mutation safe.

Legacy correction combines configurable pincer, tail and reassembly passes.
Reassembly compares forward and reverse proposals and may revisit unresolved
positions for up to six passes. Its window filter suppresses a proposal suffix
when a change-count or summed-quality limit is exceeded; it is not a low-quality
input-base filter. Correction counters count algorithm events, not necessarily
unique positions. Read rollback belongs to the concrete caller and does not make
all shared processing transactional.

`fixindels` selects a separate general substitution/one-base-indel path; `hpindel`
selects specialized homopolymer repair. The two are mutually exclusive. Local
batch, pair-witness and neural-gate options have additional compatibility checks;
see the current shell usage and source contracts before combining them. `pacbio`
supplies bundled model/default choices, with explicit options taking precedence.

Graph output can describe an earlier contig population than final sequence output,
which may undergo further selection or path extraction. Use the configured output
and its corresponding statistics rather than assuming every report describes the
same population. Existing source TODOs identify unresolved behavioral limitations;
this overview is not a blanket correctness, thread-safety or performance guarantee.
