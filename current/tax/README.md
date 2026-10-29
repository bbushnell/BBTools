# Tax Package

This package provides comprehensive taxonomic infrastructure for biological sequence classification, identifier mapping, and phylogenetic analysis. It includes tools for NCBI taxonomy integration, accession number processing, and hierarchical tree management.

## AccessionToTaxid (AccessionToTaxid.java)
Loads NCBI accession-to-taxid rows into resident maps/tables or a disk-backed
`DiskAccessionTable`, then resolves accession strings to TaxIDs for local tools
and `TaxServer`. Four-column accession files, two-column compressed rows,
optional distributed-server filtering, accessions with or without version
suffixes, and `AnalyzeAccession` pattern-map digitization are all kept in one
static lookup service. IMG genome IDs are resolved through loaded IMG metadata
with a small server-query cache for missing local entries.

## AnalyzeAccession (AnalyzeAccession.java)
Counts accession shape patterns across one or more accession tables with worker
threads. Each accession is remapped to an `L`/`D`/delimiter pattern, counted,
and written with its theoretical combination count and bit estimate. Compact
pattern codes are assigned later by `loadCodeMap` when that summary is loaded.
The resulting static pattern table is the encoder used by
`BuildDiskTables`, `AccessionToTaxid` disk lookups, and `TaxServer` to convert
accession text into 64-bit disk-table keys.

## AnalyzeAccession_ST (AnalyzeAccession_ST.java)
Single-threaded predecessor of `AnalyzeAccession`. It reads accession tables,
counts remapped accession patterns, writes pattern/code/count rows, and exposes
the same remapping, combination-counting, pattern loading, and digitization
helpers for older callers that do not need the grouped worker pipeline.

## AsyncDiskAccessor (AsyncDiskAccessor.java)
`DiskAccessor` implementation that opens a disk accession table with an
`AsynchronousFileChannel`. A batch submits each key's first partition probe
concurrently, then resolves any linear-probe collisions through the same channel;
nonpositive keys and missing slots return `-1` in the result position that
matches the input key.

## BuildDiskTables (BuildDiskTables.java)
Builds binary disk tables for low-memory taxonomy lookup. It can load serialized
or raw GI mappings and write a flat `DiskGiTable`, count accession rows after
loading the `AnalyzeAccession` pattern map, and stream accession files through
worker threads into a partitioned `DiskAccessionTable.Builder` using bounded
per-partition buffers. Tree preloading is disabled by default and does not filter
table entries. A filename selects a tree; `tree=t`, `tree=auto`, or bare `tree`
selects the default tree. Use `tree=`, `tree=null`, `tree=f`, or `tree=none` to disable
loading.

## DiskAccessionTable (DiskAccessionTable.java)
Memory-mapped, partitioned open-addressing table from digitized accession keys
to TaxIDs. The file stores a header, a partition offset/capacity table, and
12-byte key/value slots; lookups hash the nonzero key to one partition and linearly
probe until a matching key or an empty slot. The nested `Builder` accumulates
per-partition buffers, sorts each partition by initial slot, and writes the final
mapped table.

## DiskAccessor (DiskAccessor.java)
Small batch-lookup abstraction for disk accession tables. Implementations accept
an array of digitized accession keys, return one TaxID or `-1` per input position,
provide a display name for comparison drivers, share the MurmurHash3 finalizer
used for partition selection, and own any file handles they open.

## DiskGiTable (DiskGiTable.java)
Disk-backed GI-to-TaxID table. It writes a direct-indexed flat binary image from
`GiToTaxid`'s sparse in-memory arrays, then maps that image read-only in segments;
lookups return the TaxID at the GI index, `-1` for negative GI numbers, and `-2`
for GI numbers beyond the stored maximum.

## CanonicalLineage (CanonicalLineage.java)
Collects the nearest ancestor at each of seven exact ranks: kingdom, phylum,
class, order, family, genus and species. Domain and intermediate ranks are omitted.
The view borrows nodes from the supplied tree rather than copying them. Rendering
uses semicolon-separated `k__` through `s__` fields; missing ranks become `NA`,
except that missing species becomes `tid_` followed by the original TaxID.
Null trees and nonpositive IDs produce an entirely missing lineage.

## ExplodeTree (ExplodeTree.java)
Reads one FASTA input, resolves taxonomic headers through a `TaxTree`, and writes
recognized records to gzip FASTA files beneath the taxon's tree-derived output
directory. Per-taxon files are named from the numeric TaxID plus the optional
`prefix=`. `makedirectories=t` precreates one directory and `.name` marker for
each loaded taxon. An optional `results=` table reports emitted base counts by
TaxID, rank and name. Directory/name-write failures are reported, and processing
failures trigger reader/writer cleanup while preserving the primary exception.

## FakeTreeLeafTaxonomy (FakeTreeLeafTaxonomy.java)
Assigns proposed real parents for permanent FakeTree leaf nodes. It joins a
required FakeTree leaf map with FindSSU and QuickClade top-hit tables, loads a
taxonomic tree, and writes a proposal table plus a FakeTree node table suitable
for durable species/subspecies leaves named from `collection__source_id`.
Optional inspection and summary tables expose high-level conflicts, fallback
decisions, evidence-source counts and proposed node-rank counts. Eligible SSU
evidence retains priority when its resolved parent is genus or species. Broader
SSU parents, including the Life fallback for an unresolved hit taxID, yield to
usable QuickClade evidence. Without usable QuickClade evidence, existing SSU
eligibility and fallback rules apply. ANI thresholds and the diagnostic columns
for both sources are unchanged by this preference.

## FilterByTaxa (FilterByTaxa.java)
Filters single- or paired-end reads according to the taxonomy parsed from each
primary read name. It builds a `TaxFilter` from `tree=`, optional GI/accession
tables, explicit names or TaxIDs, regex/substring criteria, include/exclude mode,
promotion and required-rank settings. Accepted pairs are written to the normal
read output stream, while `results=` records the distinct parsed taxonomic nodes
seen in passing primary headers. `besteffort=t` can widen the configured rank by
scanning the first input for a level shared by requested and observed taxa. The
Streamer/writer processing path attempts input and output cleanup independently
after failures while preserving the primary exception.

## FindAncestor (FindAncestor.java)
Reads lines whose first tab-delimited field is a record name and whose remaining
line body is a comma-delimited GI list. It translates those GIs through a local
GI table, writes each valid line's lowest common ancestor, majority taxon and
traversals, and optionally copies invalid lines to `invalid=`. `lines=` limits
processed records. The majority calculation chooses the lowest-ranked taxon
containing a strict majority of translated IDs after count propagation through
the loaded tree.

## GiToTaxid (GiToTaxid.java)
GI and explicit-TaxID resolver backed by either a sparse resident `int[][]` or a
`DiskGiTable`. It can build a sparse table from raw GI mapping files, serialize
or load that table, publish disk-backed low-memory lookups, parse string or byte
headers containing `gi`, `ncbi`, `tid` and related markers, and unload the active
backend while closing any disk table.

## GiToTaxidInt (GiToTaxidInt.java)
Legacy dense GI resolver that stores one `int` slot for every GI from zero
through the largest observed ID. It can build or load a serialized one-dimensional
table, parse GI headers in string or byte form, and return either explicit TaxIDs
or the stored dense-table value; large modern GI ranges make `GiToTaxid` and
`DiskGiTable` the lower-memory path.

## IDNode (IDNode.java)
Mutable cluster node used by `IDTree`. Leaves borrow their lower-triangular
identity rows; merged nodes copy and combine rows by elementwise maximum and
track member leaves in a bitset. A priority queue selects merges by maximum
similarity, breaking ties by node number. Newick output sanitizes labels and
writes four-decimal branch lengths. Zero-similarity components are connected by
selecting another active cluster when a row maximum points into its own membership.
This is not a general-purpose taxonomy node.

## IDTree (IDTree.java)
Builds a greedy maximum-identity clustering tree without loading a taxonomy
database. Each tab-delimited input row starts with a label and contains identities
to preceding rows; values beyond that lower triangle are ignored. Lines beginning
with `#` are skipped. The optional output is a Newick tree ending in a semicolon.
The row limit counts accepted data rows. Input is closed before clustering,
including after parsing failures; cleanup errors propagate or attach to an
earlier exception.

## ImgRecord (ImgRecord.java)
Parses legacy IMG tab-delimited rows whose first columns are IMG ID, genome
name, NCBI TaxID, public flag and high-quality flag. `toArray` skips the header
and can retain only high-quality rows; `toMap` indexes records by IMG ID, with
the last duplicate winning. Empty TaxID fields become `-1`, and the static
`storeName` flag controls whether parsed objects retain genome names.

## ImgRecord2 (ImgRecord2.java)
Parses tab-delimited IMG genome IDs, NCBI TaxIDs and optional names into records.
`toArray` keeps input order; `toMap` indexes by IMG ID, with the last record winning
when an ID occurs more than once. `storeName=true` requires a third field; trailing
empty fields are discarded by the parser. `parseImgId` extracts an IMG ID from a
sequence header and is used by SketchMaker. Present overflowing IMG tokens throw
`NumberFormatException`; missing delimiters, missing tokens, and zero IDs still
return `-1` without assertions. `toArray` closes its owned reader on normal and
exceptional exits and rejects reported close errors. A parsing failure remains
primary if cleanup also fails.

## Lineage (Lineage.java)
Mutable rank-indexed view using `TaxTree`'s already-loaded shared instance.
TaxID construction walks the parent chain and stores canonical nodes. String
construction resolves semicolon-separated `rank__name` fields independently,
ignoring absent or ambiguous names; it does not validate mutual ancestry.
Repeated resolved ranks replace earlier entries. Nodes are borrowed, and the
view is not safe for concurrent mutation.

## NIODiskAccessor (NIODiskAccessor.java)
`DiskAccessor` implementation built on positioned `FileChannel` reads. Small
batches probe the disk accession table on the caller thread; larger batches split
input ranges across worker threads that share the channel, join all workers before
returning, restore caller interruption, and propagate the first slice failure on
the caller thread.

## PrintTaxonomy (PrintTaxonomy.java)
Prints taxonomy rows or semicolon lineage strings for command-line names, a text
file, or read headers. Inputs can be literal names/IDs, `name=`/`id=` lists,
tabular rows with one selected `column=`, or a read stream from `in=`.
GI/accession tables and a required `TaxTree` resolve headers. Unresolved
command-line names or read headers print a not-found line; unresolved selected
columns keep the original token with a `_***NOT_FOUND***` suffix. `level=`,
`minlevel=`, `maxlevel=`, `simple=`, `reverse=` and `printname=` control emitted
ranks and ordering; when a lineage lacks the exact requested `level=`,
`printTaxLevel` emits the next available ancestor at or above that rank.
This is intentional best-effort handling of irregular taxonomy: traversal only
ascends, retains a starting node already above the requested rank, and stops at
the root if the requested rank cannot be reached.
`counts=` can write raw visit counts for loaded-tree nodes.
The static tab-delimited formatter used by legacy BBSketch output always considers
its starting node, even above its maximum rank; that limit only bounds further
ascent. Existing noncanonical/cellular-organisms filters still apply. This preserves
the known classification without attempting to infer a more specific descendant.

## Query (Query.java)
Legacy blocking URL client for GI/accession lookup and repeated-request timing.
It joins identifiers literally into a URL and parses comma-separated integer
responses. Callers must supply URL-compatible values. Unlike `TaxClient`, it
does not encode payloads, retry transport failures or validate response counts.
Its transport prints handled I/O errors and returns null; the lookup helper does
not guard that null result. No timeout or response-size limit is configured here.

## RenameGiToTaxid (RenameGiToTaxid.java)
Rewrites FASTA headers or GFF seqids from GI/accession/name-based identifiers to
taxonomy IDs. It can use local GI, accession and pattern tables or remote server
queries selected by `server=` and `mode=`, writes renamed records to `out=`, and
can write rejected records to `invalid=` and diagnostic rejected headers or seqids
to `badheaders=`. In GFF mode it renames both feature seqids and matching
`##sequence-region` pragmas so a GFF's internal feature and pragma seqids stay
consistent; GFF input rejects `prefix=f` because it would drop feature columns.

## RenameIMG (RenameIMG.java)
Concatenates legacy IMG FASTA files and rewrites their sequence headers. The `in=`
IMG dump selects records and their FASTA paths; an optional `img=` dump supplies
the IMG-to-TaxID mapping. Headers receive `>tid|N|img|IMGID ` when a nonnegative
TaxID is available, or `>img|IMGID ` otherwise, followed by the original header
text. Nonempty sequence lines are copied unchanged. A combined `out=` destination
is required; input files themselves are not renamed. `lines=` caps the total
number of nonempty FASTA lines copied across inputs, with zero meaning unlimited.

## SerialDiskAccessor (SerialDiskAccessor.java)
Caller-thread `DiskAccessor` wrapper around an already opened
`DiskAccessionTable`. It performs each batch lookup in input order through the
table's memory-mapped `get(long)` method and deliberately leaves close ownership
with the caller that supplied the table.

## ShrinkAccession (ShrinkAccession.java)
Converts raw NCBI accession-to-taxid dumps into the smaller table formats used by
BBTools. Sequence mode retains the first accession column, leaves a blank second
column, and keeps TaxID and optionally GI; the separate versioned accession column
of four-column input is discarded. Assembly mode writes accession, a blank second
column, TaxID and selected assembly metrics, including sizes and sequence/gene
counts. Its data loop skips comment lines; downstream accession loaders can ignore
the trailing metrics.

## SplitByTaxa (SplitByTaxa.java)
Splits reads by the taxonomy parsed from each primary read name. The selected
TaxID is climbed to the configured `level=`, then the primary read and its mate
are written to a per-taxon output file made by replacing the first `%` in the
`out=` pattern with a sanitized taxon name. Headers can resolve through a
`TaxTree`, an optional GI table or tree name lookup; unresolved names go to an
`UNKNOWN` bucket. Output-path collision checks reject distinct raw names that
would open the same sanitized path. Replacement treats `$` and backslashes in
taxon names literally. Cleanup after CROS startup failure remains a known
limitation when the output consumer never starts.

## TaxApp (TaxApp.java)
Unfinished Swing demonstration: a button opens an input dialog and discards its
result. It performs no taxonomy lookup and is not an entry point for the working
taxonomy command-line tools.

## TaxClient (TaxClient.java)
Blocking HTTP client for resolving names, GI numbers, accessions and sequence headers through the configured taxonomy server. Static helpers support individual identifiers and comma-separated batches. Messages shorter than 2000 characters use an encoded GET path; longer messages use POST. Null or empty responses trigger at most twelve attempts with backoff. Numeric response fields are parsed as integer TaxIDs; malformed fields throw. The CLI uses individual requests by default, or one batch per identifier type with `fast=t`, and prints input/TaxID pairs to its diagnostic stream.

## TaxFilter (TaxFilter.java)
Reusable include/exclude predicate for read-name and TaxID filtering. It loads
optional tree, GI and accession tables, accepts names, TaxIDs, regex and
case-insensitive substring criteria, can promote selected nodes to a requested
rank, and can require that passing nodes have one or more ancestor ranks present
in their lineage. String filtering resolves parsed GI/NCBI/accession headers
before falling back to tree name lookup; integer filtering works directly on
loaded `TaxNode`s and intentionally passes empty include-mode filters so callers
can perform later name-only checks.

## TaxNode (TaxNode.java)
Mutable taxonomy node with a durable TaxID, parent TaxID, scientific name, current
normal/extended ranks, and saved original-rank/canonical flags. Parent IDs are
resolved through the owning tree; they need not be array indices. Direct and
cumulative abundance counters are mutable and unsynchronized. Name/rank helpers
support lineage matching and rendering. `equals(TaxNode)` compares selected
structural fields while ignoring counters; it is an overload, so equality through
`Object` still uses object identity.

## TaxServer (TaxServer.java)
HTTP service for taxonomy names, TaxIDs, GI numbers, accessions, IMG IDs, sequence
headers and SendSketch routes. Startup loads a `TaxTree`, optional in-memory or
disk GI/accession tables, optional IMG and size metadata, and optional sketch
indexes; request handlers can return plain TaxIDs, JSON path/name objects, or
semicolon-delimited lineages through `gi`, `name`, `taxid`, `header`, `accession`,
`img`, `silvaheader` and prefixed `pt_`/`sc_` route types. Distributed accession
servers route missing local accessions to the shard selected by the same
`Math.floorMod(hashCode, serverCount)` formula used when loading.

## TaxSize (TaxSize.java)
Scans FASTA records whose headers encode taxonomy IDs and counts direct bases
and sequences for each resolved taxon. With `out=`, it emits a taxsize TSV
containing direct base/read counts plus cumulative base, read and descendant-node
counts propagated from every loaded node up to the taxonomy root. `nodesC`
includes the emitted node itself plus all of its descendants. Records whose
headers do not resolve are counted in input totals but omitted from assigned
per-taxon maps.

## TaxTree (TaxTree.java)
Loads and indexes taxonomy trees for ID/header/name lookup, rank traversal,
ancestry and common-ancestor queries. Native TaxIDs normally address node-array
slots; auxiliary IDs use a separate durable-ID-to-slot map. Merged-ID mappings
resolve obsolete IDs. Abundance propagation and node gathering mutate node counts,
so callers sharing a tree must coordinate those operations. Optional name/size
indexes and the process-wide shared-tree cache have their own loading lifecycle.
`TaxTreeText` handles the portable text representation; `TaxTreeAux` attaches
auxiliary taxa.

## TaxTreeAux (TaxTreeAux.java)
Loads a tab-delimited table with named `fake_taxid`, `parent_taxid`, `rank` and
`name` columns. Auxiliary IDs start at 2000000000 and remain external identifiers
while nodes occupy appended packed slots. A nonempty overlay copies base nodes,
shares the merged-ID map, validates new IDs/parents/cycles, and rebuilds rank
indexes. A header-only table returns the base tree unchanged. Successive overlays
preserve earlier auxiliary IDs and private slots, allow new nodes to reference
earlier auxiliary parents, and propagate rank bounds through the combined hierarchy.
Internal packed slots are rejected as external parent IDs.

## TaxTreeText (TaxTreeText.java)
Reads and writes versioned `TaxTree2` text with A48-encoded numeric fields. The
format stores native topology, rank metadata, flags and merged IDs. Writing
rejects auxiliary mappings, optional runtime indexes and nonzero abundance
counters that the format cannot preserve. The CLI converts a saved tree, reloads
the output and compares its persisted and reconstructed fields.
