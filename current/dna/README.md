# DNA sequences, references and legacy annotations

The `dna` package supplies nucleotide/amino-acid conversions, reference storage,
scaffold metadata, genetic codes and older gene/motif utilities. The classes have
different ownership and coordinate contracts; a shared package name does not make
their mutable state thread-safe.

## Sequence and reference code

| Class | Responsibility and important contract |
|---|---|
| [AminoAcid](AminoAcid.java) | Static alphabet tables, codon translation, reverse complements and packed k-mers. Most lookup helpers require table-supported ASCII inputs. Read each overload's treatment of ambiguity, case and U. |
| [GeneticCode](GeneticCode.java) | Built-in translation tables 4, 11 and 25, plus validated custom TSV tables. Separates initiation from elongation without changing global AminoAcid tables. |
| [ChromosomeArray](ChromosomeArray.java) | Mutable base storage with inclusive logical bounds and spare capacity. Setters normalize/copy input; padded accessors return N outside logical bounds while retaining stored bytes inside them. |
| [Data](Data.java) | Process-wide reference and annotation caches, paths, scaffold lookup and external-tool probes. Initialize before publishing data to workers; coordinate resets and mutation of borrowed objects. |
| [FastaToChromArrays2](FastaToChromArrays2.java) | Byte-oriented FASTA conversion used by reference indexing. Produces packed chromosomes and metadata with configurable padding and scaffold limits. |
| [ChromArrayMaker](ChromArrayMaker.java) | Separate legacy FASTA/FASTQ converter using CrisWrapper. It has different boundary and error-handling behavior; it is not an interchangeable implementation of FastaToChromArrays2. |
| [ChromToFasta](ChromToFasta.java) | Converts serialized chromosome arrays to whole-chromosome or N-run-split FASTA. Helpers borrowing a writer leave it open. Output line widths must be positive. |
| [ScafLoc](ScafLoc.java) | Mutable scaffold name, chromosome number and start coordinate, also returned from Data's shared lookup map. |
| [Scaffold](Scaffold.java) | Scaffold metadata and coverage counters. Counter units matter: `fraghits` is a doubled fragment accumulator for its CoveragePileup caller. |

`AminoAcid.toNTs` chooses representative codons; it cannot reconstruct the original
DNA from a protein. Unresolved amino acids produce NNN. The normal alphabet uses
a packed codon table, with a fallback that preserves substitutions in the public
nucleotide alphabet. Public lookup arrays remain mutable shared state: do not
modify them concurrently with readers.
`splitBase` expands supported nonempty four-base masks and rejects unsupported
symbols or indices instead of returning partially filled character arrays.

Chromosome positions are zero-based and inclusive. `ChromosomeArray.getBytes(a,b)`
returns independent storage. Empty sequence writes leave bounds unchanged, and
resizing rejects a capacity that would discard the base at `maxIndex`. Writes can
leave zero-filled gaps inside logical bounds; those are distinct from N padding
outside the bounds. Numeric packing helpers read storage directly, so they do not
inherit padded-accessor behavior. Two bits per base fit at most 16 bases in an int;
legacy oversized-input limitations are noted at the affected methods.

Data's scaffold lookup selects a scaffold using its own padding convention.
`scaffoldOverlapLength` counts overlap with one selected scaffold, not the sum over
all scaffolds touched by a query. The volatile scaffold-name map publishes a fully
constructed map; it does not protect later mutations or make every other cache safe.
`unloadAll` does not return every static cache to its initial state; remaining
legacy cache/reset limitations are documented in the source.
The genome completion guard is assigned only after normal initialization return;
failed attempts can leave partial state, and retries using in-memory metadata must
supply fresh complete lists. `setPath(null)` follows the literal `"null"` reset,
while an empty string remains a custom `/` prefix. Existing Windows selectors
choose relative reference paths and `C:/workspace/tempdir/`.
For nearest-gene lookup, interval queries select all overlapping padded merged
ranges, or the closer bounding range in a gap (left on ties). The legacy point
overload retains its containing-or-preceding selection instead.

Reference conversion checks retained input and metadata-writer errors before
publishing a new completion summary. This is not a transactional build: files
written earlier can remain after failure. Exceptional cleanup and the deprecated
`writeInfo` path have separately documented limitations. Conversion settings are
shared static state and must not be changed by concurrent conversions.

## Genes and motifs

| Class | Contract |
|---|---|
| [Gene](Gene.java) | Legacy transcript annotation with inclusive coordinates and borrowed exons in ascending genomic order. Strand-aware offsets count exonic bases. The empty constructor creates a sentinel, not a valid annotation for every operation. |
| [Exon](Exon.java) | Inclusive genomic interval and annotation flags. Shared interning requires coordinated access. |
| [GeneSet](GeneSet.java) | Groups transcript records and derives aggregate bounds/flags. Member objects are shared; constructor and comparison preconditions are documented separately. |
| [Matrix](Matrix.java) | Loads named float matrices and copies requested subgrids. Name lookup supports stored original/lowercase aliases, not arbitrary case-insensitive matching. |
| [Motif](Motif.java) | Common motif geometry and scans. Candidate positions are motif centers; scan bounds account for both center and length. |
| [MotifSimple](MotifSimple.java) | Exact and extended alphabet matching under its documented symbol/center preconditions. |
| [MotifMulti](MotifMulti.java) | Combines component motifs with compatible geometry. |
| [MotifProbsN](MotifProbsN.java) | Scores position-specific nucleotide tuples. N is tuple width (1, 2 or 3), not the ambiguous nucleotide symbol. Mutates and retains the supplied probability matrix and reuses scratch state. |

Gene's `toBuild` removes `build`, `b` or `hg` prefixes and an optional equals sign,
returning the numeric suffix. For example, `hg18` returns 18. A caller such as
Translator2 separately maps UCSC build labels to its supported NCBI builds.

Gene's legacy `translationStart/Stop` methods use the outer exon bounds when exons
are present, which can include UTR sequence. Equal coding endpoints use a legacy
zero-length convention; the `untranslated` flag remains caller-supplied. Structural
identity and `equals` are different relations. Whole UTR exons count both inclusive
endpoints; biological 5-prime/3-prime totals follow the strand and remain zero when
`untranslated` is true. Merging preserves a missing receiver source instead of
adopting the donor's source. Stored serialized annotations retain their saved
derived lengths; regenerate them to adopt corrected constructor calculations.
Other legacy annotation-range questions remain explicitly annotated.

Motif scoring is not a general thread-safe service: some implementations retain
mutable matrices or reusable scratch arrays. MotifProbsN rejects negative/nonfinite
probabilities, wholly uninformative matrices and score ranges that cannot be
normalized with finite float constants. Individual background-only rows and zero
probabilities remain valid when the whole matrix has a usable range. Its direct
normalization methods still require suitable input strengths; the logarithmic
variant requires a positive minimum. Construction can modify supplied rows before
rejecting them, and later mutation does not refresh cached extrema.
Matrix loading publishes its registry only after
all default files load and close successfully; a catchable failure can be retried.
This does not make the mutable registry thread-safe. These older
utilities should be used with the documented preconditions, not inferred guarantees
from method names.

## Developer checks

Test classes live in [test/dna](../test/dna/). Build both `dna` and `test.dna` with
their normal dependencies before running the developer wrappers:

- [dev/genetic_code_test.sh](../../dev/genetic_code_test.sh) checks built-in/custom
  translation, initiation, invalid inputs and isolation of shared tables.
- [dev/resourcepathtest.sh](../../dev/resourcepathtest.sh) checks resource lookup.
  It requires exclusive access to the actual runtime and the absence of its named
  audit fixtures before execution. It creates small files under that runtime and
  leaves them for the executing controller to preserve or remove. See
  [DataResourcePathTest](../test/dna/DataResourcePathTest.java) for the full contract.

Other finite harnesses cover alphabet conversion, chromosome bounds, scaffold
coordinates, reference/output errors, subprocess exit-status handling, build
prefixes, matrix/genome initialization retries, path selectors and gene annotations.
Their checks remain active with assertions disabled as well as enabled.
Performance measurements must use matched source, input, JVM settings and work
counts; isolated kernel timings do not establish whole-program speed or heap needs.
