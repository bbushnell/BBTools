# Repeat detection and palindrome support

This package contains two related facilities: genome repeat detection through
`findrepeats.sh`, and a Java API for finding reverse-complement stems used by the
CRISPR finder. Repeat detection uses kmer counts, sequence entropy or exact short
tandems. It does not align copies of a repeat.

## Genome repeat detection

Run `findrepeats.sh` from a BBTools installation. Its Java entrypoint is
`repeat.RepeatFinder`. The intended input is unpaired genome sequence, typically
FASTA. Paired input is unsupported: a nonnull `in2`, forced interleaving, or
interleaving detected from the input is rejected before processing streams start.
This prevents silent mate loss; it does not implement paired repeat detection.
Use `interleaved=f` when each sequence should be scanned as a separate contig.

`in=` selects the input, `out=` the tab-separated repeat report, and `outs=` the
optional sequence output. The default report destination is standard output;
`out=null` suppresses it. Sequence output is disabled until `outs=` is supplied.
Sequence output defaults to soft masking (`mask=t`, `softmask=t`, `print=f`).
`print=t` selects extracted repeat sequences and disables masking. `mask=hard`
writes N at covered positions, while `mask=soft`
lowercases them. If both masking and extraction flags are supplied, masking takes
precedence when it remains enabled at the end of parsing.

The detectors can be enabled independently:

| Detector | Controls | Meaning |
|---|---|---|
| Depth (on by default) | `depth=t`, `k=31`, `mindepth=2`, `maxdepth=-1`, `gap=0`, `minrepeat=0` | Counts kmers first, then tracks regions meeting each depth threshold. |
| Low entropy (off by default) | `entropy=t` or a numeric cutoff, `ew=80`, `ek=5` | Retains low-complexity windows; the default entropy cutoff is 0.70. |
| Short tandem (off by default) | `tandem=t`, `stmink=2`, `stmaxk=15`, `stmincount=4`, `stminlen=32` | Finds repeated exact nucleotide units within those size/count bounds. |

For depth detection, `gap` counts consecutive missing qualifying kmer windows,
not substitutions. One substituted base can disturb k consecutive windows.
`minrepeat` is a span in bases. `mindepth` is clamped to at least two; a positive
`maxdepth` caps the thresholds tracked, while observed maximum counts can exceed
that cap. `qhdist` searches substitution neighbors and becomes expensive rapidly
as it increases. The launcher documents depth k from 1 through 31; internal
assertions accepting 32 are not a promise that every dependent path supports it.

Tandem units use a two-bit nucleotide encoding and require k below 16. This is
not an amino-acid tandem detector. Entropy settings must satisfy the tracker's
alphabet/encoded-key limits, positive k no greater than 15, and a window strictly
larger than k. The default 5/80 settings are suitable for ordinary longer regions.
Reported regions containing at most one complete entropy kmer have entropy zero;
they do not construct a tracker with an invalid window equal to k. Longer regions
use the configured k and a window bounded by their length and `ew`.

All repeat coordinates are **zero-based inclusive**: a region from a through b
has length b-a+1. Each region retains its detector label (`R` depth, `E` entropy,
`T` tandem), fixed depth threshold, observed maximum depth, gap statistics,
entropy, GC fraction and distance to the nearer contig tip. Entropy/tandem depth
fields are bookkeeping values, not independently measured copy counts.

The report ranks by depth threshold, length, fewer gaps, higher entropy, lower
GC, then earlier contig/start. It can contain overlapping entries at different
thresholds. Strong subsumption additionally compares gap totals; `weak=t`
ignores that gap comparison. Merged ranges are used to count the union of repeat
bases and to mask sequence, so summing report row lengths is not the union count.
GC excludes undefined bases from its denominator. Sequence previews can be
shortened with `preview=`; cached previews can contain an ellipsis and are not
full sequence data.

The observed maximum includes the first qualifying depth window. Entropy counts
include a qualifying terminal region. Tandem coverage is the union of overlapping
qualifying units, including a second unit ending at the final base. For example,
with two-base units and a minimum of two copies, ATATA has qualifying intervals
0..3 and 1..4 and therefore covers five bases. Output shutdown failures propagate
to the tool's error state, and filenames may contain `=` after the option separator.
Disabling previews removes the preview column from both the header and data rows.

## Palindrome API

`PalindromeFinder` searches a byte array, or an inclusive subregion supplied to
`longestPalindrome(bases, minPos, maxPos)`. Construct it with minimum arm length,
minimum/maximum loop length, minimum matches, maximum mismatches, minimum/maximum
tail length and maximum difference between the two tail lengths. Minimum matches
is clamped to one; minimum arm length is at least that match count.

The returned `Palindrome` is an independent mutable scalar copy, or null when
the scan finds no qualifying candidate. Its fields `a` and `b` are inclusive
bounds; `matches` and `mismatches` count paired positions in **one arm**.
`plen()` is matches+mismatches, and `loop()` is (b-a+1)-2*plen(). Tails are the
uncovered left/right portions of the requested search region. The ranking favors
more matches, then fewer mismatches, then longer span, then a leftward position;
the method name does not mean span alone determines the result.

Use uppercase nucleotide input. Complement lookup returns uppercase complements
and the left side is compared literally, so arbitrary lowercase input does not
have a general case-insensitive matching guarantee. This is not an IUPAC ambiguity
matcher. Input arrays are borrowed and not modified. Finder instances reuse
scratch objects and must be confined to one caller at a time; each concurrent
worker needs its own finder. `tracker` and `trackerFull` are populated by callers,
not automatically by `longestPalindrome`.

`jgi.CrisprFinder` owns one finder per worker, assigns results to CRISPR repeat
pairs, and updates `tracker.PalindromeTracker`. `structures.Crispr` formats those
results, including the one-arm length, loop and tails. `PalindromeFinder.main`
is a separate in-memory diagnostic accepting literal sequences or files; it is
not the streaming genome-repeat launcher.

The scan includes minimal even palindromes and starts outside the smallest loop
allowed for each center's parity. Permitted interior mismatches do not themselves
enlarge the loop. A center is skipped only when its available paired positions
cannot meet the minimum arm length or match the best retained match count.
Equal match counts are still searched so the remaining ranking rules apply.

## Internal ownership and checks

`Repeat` borrows its contig and base array. Its clone is shallow; sequence
extraction explicitly copies bases. `RepeatSet` owns worker-local depth slots,
accepted-repeat lists and mutable entropy/range caches. Its range buffers are
borrowed and overwritten on reuse. Retained repeats keep contigs alive, including
after worker lists are handed to the parent for final reporting. Constructor and
processing code also change shared I/O settings, so independent concurrent
`RepeatFinder` instances are not isolated jobs.

Tests live under `current/test/repeat`. `RepeatBoundaryTest` checks finite count
and coordinate boundaries with independent controls. `PalindromeOracleTest`
independently enumerates intervals and arm lengths for every A/C/G/T sequence of
length two through six, both whole and padded, with zero or one permitted mismatch.
Its 349,184 comparisons cover eight loop ranges and two minimum arm lengths;
they do not certify arbitrary tail constraints, longer inputs or other alphabets.
`RepeatIoBoundaryTest` checks filenames, short statistics, preview columns and
writer failures; `RepeatPairedInputTest` checks paired rejection, unpaired output
and reader closure. `RepeatTimingTest` provides
deterministic entropy, tandem, perfect-palindrome and imperfect-palindrome workloads
with exact output checks, configurable warmup and run length. It measures in-memory
detector calls, including result-list clearing, rather than full CLI throughput. Use matched builds and JVM
settings, and at least ten seconds of compute for a performance comparison.
