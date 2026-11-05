# Multiword k-mer tables

This package represents DNA k-mers in one or more `long` words and supplies table,
loading and traversal implementations for those keys. The [kmer](../kmer/README.md)
package provides related short-key tables and the shared table-set lifecycle.

## Main components

| Classes | Role |
|---|---|
| [Kmer](Kmer.java) | Mutable sequence-key representation, reverse complements, layout geometry and cached lookup state. |
| [AbstractKmerTableU](AbstractKmerTableU.java) | Common multiword table API and formatting/allocation helpers. |
| [KmerTableSetU](KmerTableSetU.java) | Configures and owns shards, loads inputs and exposes higher-level count and assembly operations. |
| [HashArrayU](HashArrayU.java), [HashArrayU1D](HashArrayU1D.java) | Array probing and forest overflow; U1D stores scalar integer values. |
| [HashArrayU2D](HashArrayU2D.java), [HashArrayUHybrid](HashArrayUHybrid.java) | Multiple-value and compact scalar/list storage. Value-array ownership and supported updates are backend-specific. |
| [HashArrayH1D](HashArrayH1D.java), [HollowAnchorTable](HollowAnchorTable.java), [HollowKmerLookup](HollowKmerLookup.java) | Specialized hash/anchor storage and lookup support with distinct configuration requirements. |
| [HashForestU](HashForestU.java), [KmerNodeU](KmerNodeU.java), [KmerNodeU1D](KmerNodeU1D.java), [KmerNodeU2D](KmerNodeU2D.java) | Bucket trees and scalar or multiple-value tree nodes. |
| [HashBufferU](HashBufferU.java), [KmerBufferU](KmerBufferU.java) | Buffered table-update coordination and its staging storage. |
| [HistogramMakerU](HistogramMakerU.java), [DumpThreadU](DumpThreadU.java), [OwnershipThread](OwnershipThread.java) | Worker coordination for table traversal, output and ownership initialization/clearing. |
| [WalkerU](WalkerU.java) | Sequential cursor exposing the current key and value after `next()` succeeds. |

## Key layout and reuse

`Kmer` is mutable scratch state. A query or walker may reuse its key object, so
copy keys that must be retained across subsequent mutations or cursor advances.
Use the documented mutation methods and preserve their cache-invalidation rules.

Packed and symmetric layouts have different word geometry. The effective full
length can differ from a requested length in layouts that round it. Use `Kmer`'s
layout helpers and the table set's effective length instead of reconstructing
geometry from the number of words. Layout and canonical-orientation settings such
as `Kmer.PACKED` and `Kmer.rcomp` are process-wide configuration; keep them stable
while the associated keys, buffers and tables are in use.

Hash-only and hollow modes have additional constructor constraints and different
storage capabilities. They should not be substituted for full-key storage merely
because both offer count lookup methods. Use the configured table-set API and
check the concrete backend before relying on sequence reconstruction or mutation.

## Values, ownership and concurrency

`1D` denotes scalar integer values; `2D` denotes multiple integer values per key.
These names do not establish identical increment, duplicate-removal, regeneration
or array-copy behavior. Some retained operations are unsupported. The concrete
method contract governs what is copied, retained, replaced or returned live.

Ordinary table operations do not all acquire locks. Follow the backend's shard
ownership or explicit locking protocol, and keep tables stable during traversal.
Exclude concurrent access while resizing or rebalancing. Ownership markers
coordinate claims on keys separately from structural synchronization. Borrowed arrays and node references must not be
mutated independently of their owning table.

The shared loading lifecycle is implemented by
[AbstractKmerTableSet](../kmer/AbstractKmerTableSet.java). Its memory helpers return
planned budgets and capacity estimates. `clear()` releases backing references;
it is not a complete reset of all allocation and loading state.

## Traversal and output

A `WalkerU` supports `next()`, `kmer()` and `value()`. Its current `Kmer` may be
reused on the next advance. Histogram helpers accumulate each backend's scalar
view; for count-valued tables, the last bin also receives higher counts.

Dumping writes formatted text through a caller-owned or table-set-owned writer,
as documented by the entry point. `DumpThreadU.dump` joins its table workers and
reports their success, leaving writer completion to its caller. Table-set dump
methods finish their writer and return **true for a reported error**. A successful
traversal return alone does not certify asynchronous output completion.

Contract regression harnesses are under [test/ukmer](../test/ukmer/) and
[test/kmer](../test/kmer/). Dump-worker regressions include both short and long
table sets. Method Javadocs and source TODOs identify backend-specific
preconditions and remaining limitations.
