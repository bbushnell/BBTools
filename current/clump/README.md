# Clump package

Clumpify groups reads by selected k-mer pivots. Related sequences can then be
stored near one another, deduplicated, corrected, or condensed into consensus
reads. The command-line entry point is `clumpify.sh` in the BBTools root;
its help describes the supported options.

The implementation uses ConcurrentReadInputStream and
ConcurrentReadOutputStream. The launcher supports
FASTA and FASTQ input, including paired reads, rather than general SAM/BAM input.

## Processing flow

- **One group:** Clumpify delegates to KmerSort1, which holds the reads for
  sorting and clump processing.
- **Multiple groups:** KmerSplit writes the initial partitions. Clumpify invokes
  the selected KmerSort implementation on them; version 3 is the default.
  Subsequent version 2/3 passes reuse partition files instead of starting over
  from the original input.
- **Repartitioning:** choose a pivot using the next pass's comparator, then
  distribute by `hash(pivot)%groups`. The next comparator increments the seed
  and reduces the border, with a lower bound of zero.
- **Final ordering:** requested name sorting or pair repair adds a final ordering
  step. Ordinary clump output is not promised to retain input order.
  With multiple groups, `unpair=t`, and two mate-output files, mates are restored
  after the global name sort. Names must match the usual FASTQ mate convention
  (identical names are also accepted). An unmatched read causes a clear error;
  there is no singleton-output destination in this path. Original inputs are
  retained on that failure. This mode uses an additional temporary-file pass.
- **Memory estimates:** sampled input and heap availability guide group counts.
  Uneven partitions can still exceed the estimate. KmerSort3 has
  an overflow fallback that can stream a group without the in-memory clump
  operation; its diagnostics describe when that occurs.

Keep KmerSort1, KmerSort2, and KmerSort3 distinct. Their loading, grouping and
pass orchestration differ; they are not interchangeable duplicate declarations.
In-memory sorting and clump processing retain the loaded reads.

## Components

| Class | Responsibility |
|---|---|
| Clumpify | Parses orchestration settings, chooses groups, runs passes, and manages temporary files and optional input deletion. |
| KmerSplit | Applies first-pass transformations and partitions reads by hashed pivot. |
| KmerSort | Shared sorter operations, hashing/splitting workers and statistics. |
| KmerSort1 | Original sorter, including the one-group path and its own pass loop. |
| KmerSort2 | Processes previously split groups and writes the next partition set. |
| KmerSort3 | Coordinates grouped fetching, memory estimates and overflow streaming across externally orchestrated passes. |
| Clump | Holds reads sharing a pivot and implements duplicate removal, consensus and correction operations. |
| ClumpList | Builds clumps from pivot-sorted reads, reorders clumps, and dispatches consuming per-clump work. |
| Splitter | Separates clumps using supported sequence variants and optional correlations. |
| KmerComparator | Selects pivots and orders attached read keys, with sequence/quality/name tie breaking. |
| KmerComparator2 | Common optical-comparator tie breaking for reads with attached keys. |
| KmerComparatorX / KmerComparatorY | Preserve the pivot/strand/position prefix and add the selected optical-coordinate tail. |
| ReadKey | Mutable pivot, position, orientation and optional flowcell/UMI metadata attached through Read.obj. |
| PivotSet | Counts selected pivots, including linked mates, in an approximate count array. |
| ClumpTools | Caches pivot counts and checks expanded input/output names before stream startup. |
| Hasher | Whole-sequence hash hints and exact byte comparisons for duplicate filtering. |
| KmerReduce | Emits a k-base representation of each root read's selected pivot; linked mates are not emitted separately. |
| StreamToOutput | Consumes an input stream into borrowed outputs, optionally through a temporary name-sorted FASTQ file. |

## Pivots and metadata

The packed nucleotide representation supports k from 1 through 31. A normal
canonical pivot uses the **larger** packed value of the forward k-mer and its
reverse complement. Selection also depends on comparator seed, border, quality
settings and any count table. A pivot is a grouping key, not proof that two
complete sequences are identical.

With a positive `minprob` and qualities, selection prefers a qualifying k-mer
over every nonqualifying one, then chooses the largest hash among qualifying
candidates. If none qualifies, it retains the largest-hash fallback. The default
`minprob=0` does not filter by quality. Count-table selection has its own threshold.

Short reads and reads with no usable full k-mer use KmerComparator's fallback
encoding. PivotSet counts those fallback values; it does not simply discard
every read shorter than k.

ReadKey.position denotes the zero-based rightmost pivot base. Its comparison
can distinguish keys that its equality operation treats as equal, particularly
under early optical-coordinate sorting. ReadKey.clear() resets pivot metadata,
not every field in the object.

Hasher normalizes ASCII letter case when hashing, while its equality checks
compare bytes exactly. Hash matches therefore require an equality check.
The fixed seed gives repeatability for a particular random-generator
implementation, not an immutable cross-version hash format.

## Ownership and reuse

In the default I/O pipeline, callers supply already-started output streams to
StreamToOutput and retain responsibility for closing those outputs. The bridge
owns consumption and closure of the input during processing.

Clump processing mutates reads and their attached keys. ClumpList.process()
consumes its clumps and advances its work pointer; it is not a repeatable query.
The reordering methods retain added marks, and the cached connectivity map must
not be used after unrelated inherited-list mutations. See the method contracts
before reusing these objects.

Parsing changes shared settings in Clump, KmerComparator, FASTQ and ReadWrite.
In-JVM callers must serialize configuration and processing; the entry points do
not universally restore every setting on every exit. ClumpTools' file cache
keys contain filenames and selected counting parameters, not file-content
fingerprints. Replacing an input file in place requires clearing that cache.

Input/output name checks operate on expanded descriptor names using Tools'
comparison rules. They do not resolve every possible filesystem alias.
Temporary names and memory estimates are conveniences, not guarantees of
exclusive path reservation or bounded heap consumption.

## Optional phase timings

For grouped V3 processing, `phasetimers=t` prints tab-separated diagnostic
rows: `PhaseTimer`, stage class, metric name, seconds. The default is false.
Each sorter invocation prints its own block, so multiple passes remain separate.
Coordinator times accumulate across all groups. Worker `elapsed_sum` values
can overlap and are neither CPU time nor additions to coordinator wall time.
The splitter's hash time is included in its transform/hash/route total.

Measurements themselves have overhead; compare overall CPU and wall time with
profiling disabled before claiming a speed change.

## Regression fixtures

Retained Java fixtures live in [test/clump](../test/clump), outside the production
package. They distinguish saved baseline behavior from repaired expectations
and include finite parameter, name, comparator, filename, stream-lifecycle and
overflow-partition cases, mate preservation, UMI handling and grouped pivot recounts.
Small fixtures check behavior; their runtimes are not throughput benchmarks.
Performance comparisons need representative inputs, consistent JVM/thread settings
and sufficiently long runs.
