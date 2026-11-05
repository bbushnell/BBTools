# Short-key k-mer tables

This package supplies primitive-key tables, loading and traversal helpers used by
counting, filtering and assembly programs. Keys are encoded in a Java `long`.
For keys represented by several words, see [ukmer](../ukmer/README.md).

## Main components

| Classes | Role |
|---|---|
| [AbstractKmerTable](AbstractKmerTable.java) | Common table API, formatting and allocation helpers. Concrete backends define supported operations. |
| [KmerTableSet](KmerTableSet.java) | Owns shards, loads reads and exposes count, extension and ownership operations. |
| [AbstractKmerTableSet](AbstractKmerTableSet.java) | Shared prefilter, input, reporting and histogram lifecycle for short and long table sets. |
| [HashArray](HashArray.java), [HashArray1D](HashArray1D.java) | Array probing with forest overflow storage; the 1D implementation stores scalar integer values. |
| [HashArray2D](HashArray2D.java), [HashArrayHybrid](HashArrayHybrid.java), [HashArrayHybridFast](HashArrayHybridFast.java) | Multiple-value and compact scalar/list storage families. Append, replacement and duplicate handling depend on the concrete method. |
| [HashForest](HashForest.java), [KmerNode](KmerNode.java), [KmerNode1D](KmerNode1D.java), [KmerNode2D](KmerNode2D.java) | Bucket trees and their scalar or multiple-value nodes, also used for array overflow. |
| [KmerTable](KmerTable.java), [KmerLink](KmerLink.java) | Linked-bucket table and chain nodes. Several legacy operations are explicitly unsupported. |
| [HashBuffer](HashBuffer.java), [KmerBuffer](KmerBuffer.java) | Buffered update coordination and its primitive staging storage. The storage object does not itself flush a table. |
| [TableLoaderLockFree](TableLoaderLockFree.java), [TableReader](TableReader.java) | Reference-loading and read-matching/masking helpers with their own configuration and ownership requirements. |
| [HistogramMaker](HistogramMaker.java), [DumpThread](DumpThread.java), [OwnershipThread](OwnershipThread.java) | Helpers that distribute table traversal or ownership initialization/clearing among workers. |
| [ScheduleMaker](ScheduleMaker.java) | Plans table capacities and growth within memory budgets. |
| [Walker](Walker.java), [KmerTableInterface](KmerTableInterface.java) | Cursor interface and the small GC-count accumulation interface. |
| [SimpleKmerTable](SimpleKmerTable.java) | Retained placeholder; its constructor allocates no table. |

## Keys, values and storage ownership

A scalar count table and a table of associated integer values serve different
contracts. The `1D` and `2D` names distinguish scalar versus multiple-value storage.
Read the concrete method before treating an increment, set, returned array or
histogram value as interchangeable across these families.

The scalar table absence sentinel is `AbstractKmerTable.NOT_PRESENT` (`-1`). A
stored zero can still be a present entry. Higher-level count helpers may clamp
absence to zero for their particular algorithm. Use the API's documented return
convention rather than inferring presence from every zero-valued result.

Returned arrays, bucket roots and nodes may be borrowed live storage. Resizing or
rebalancing can invalidate structural references. Copy data that must outlive a
cursor or mutation, and follow each backend's copy/retain contract when supplying
value arrays. A walker exposes its current entry only after a successful `next()`.

## Reference loading and query sampling

`TableLoaderLockFree.processData` adds keys to borrowed tables and reports input
statistics for that invocation. Without retained names, calls may append to existing
tables. With `keepNames=true`, every destination must be empty, including overflow
storage, because scaffold IDs restart at one. A rejected nonempty call preserves
the old keys and metadata. Use fresh tables or a backend's supported clear operation
before loading a new named reference set. The accumulated `errorState` is not reset.
Requested paired-read overlap correction is independent of name retention.

Configure `TableReader` and the loader with matching key lengths, orientation and
middle masking. Their speed selector uses the encoded key's remainder modulo 17:
nonpositive settings disable sampling, settings 1 through 16 select remainders at
least that value, and settings at least 17 select no exact-path reference or query
keys. Mutation-expanded reference keys are retained without this filter; queries
still apply it. Query-position skipping is a separate selection step.

## Concurrency and lifecycle

Ordinary table lookup and mutation are not universally synchronized. Callers must
follow the backend's locking or shard-ownership protocol and exclude structural
mutation during traversal. An ownership marker coordinates claims on a key; it
does not by itself make table resizing, insertion or clearing safe concurrently.

Table sets have an allocation/loading lifecycle. `clear()` releases backing-table
references; use a fresh table-set instance for a new loading lifecycle. Loading
also uses process-wide read settings, so conflicting configurations must not run
concurrently in the same JVM. Constructor memory values and capacity schedules
are budgets and estimates, not measurements of live heap consumption.

Histogram and ownership helper entry points join their workers before reporting
worker failures. Callers still own the surrounding tables and resource lifecycle.
The concrete backend determines which values and entries its traversal includes.
For scalar count tables, histogram bin `c` counts distinct keys with count `c`,
with larger counts accumulated in the final bin.

## Output and error conventions

The byte-oriented dump methods produce formatted text, including FASTA or numeric
representations selected by the formatting settings. Their method names describe
the output buffer representation.

The table-set methods `dumpKmersAsBytes` and `dumpKmersAsBytes_MT` return **true
when an error is reported**. They own completion of the writer they create. The
lower-level `DumpThread.dump` helper instead returns whether every worker returned
normally; its caller must also finish and check the writer. Individual backend
traversal returns do not establish asynchronous output success. Keep these layers'
contracts explicit when adding a caller.

Contract regression harnesses live under [test/kmer](../test/kmer/). They cover
lookup semantics, state transitions, failure propagation and finite output oracles.
Existing variants and unsupported operations are retained deliberately; consult
the method Javadocs and source TODOs before extending a backend's supported use.
