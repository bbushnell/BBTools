# Primitive maps, sets and sequence indexes

This package provides primitive-key collections, object-key collections, packed
bit sets and sequence-overlap candidate indexes. The numbered variants preserve
different hashing, sizing and storage strategies; they are not interchangeable
solely because their names are similar.

## Contracts to check when choosing a collection

* Check the method's return value. Some `put` methods return the previous value,
  others report whether a new key was inserted, and `LongHashMap.put` leaves an
  existing mapping unchanged. List-valued methods can append, replace or copy:
  check the overload, especially `IntListHashMap.put(key, list)` (replacement)
  versus `put(key, int)` (append).
* Use `contains` to distinguish absence from a stored sentinel-valued result.
  Numeric maps commonly return `-1` for a missing key but can store `-1` too.
  Object-valued maps can store `null`. Empty-slot values are not an occupancy test.
* Primitive-key maps use an internal sentinel that can change when that value is
  inserted as a real key. When reading `keys()`, compare against that instance's
  current `invalid()` value. Object-key maps reserve `null` for empty slots.
* `keys()`, `values()` and list-valued `get()` results expose backing storage.
  Resizing can replace arrays; removing or clearing entries can clear or reuse
  lists. Copy data that must survive subsequent mutation. `toArray()` creates a
  separate key array; its order is not a sorting guarantee.
* Maps and ordinary sets require external synchronization for concurrent mutation.
  Internal synchronization in a constructor or resize helper does not make other
  operations thread-safe. `AtomicBitSet` provides atomic cell updates, but bulk
  operations are not whole-set snapshots and capacity changes need exclusive access.
* Numeric variants differ in overflow behavior. Read the specific `increment`
  contract; do not assume all additions saturate or that negative underflow clamps.
  Legacy counting-only helpers can have stricter requirements than general storage.
  Numeric `setToMax` copies absent keys exactly, including negative values; it
  takes the maximum only when the destination already contains that key.

Pre-size a collection when a useful estimate is available. Reuse collections and
primitive buffers where appropriate, while respecting their ownership contracts.
Do not substitute a boxed map or a different primitive variant without checking
absence, update, iteration and overflow semantics.

## Bit sets and sequence indexes

`RawBitSet` stores one bit per element; `MultiBitSet` stores saturating two-bit
counts. `AbstractBitSet.add` merges into the destination and clears the distinct
source. Its current `setToMax` overloads differ: the RawBitSet source is cleared,
whereas the MultiBitSet source is preserved. That compatibility question remains
documented in the source; callers must not assume uniform source ownership.

`SeqMap` indexes borrowed sequence arrays. Build and mutate it with exclusive
access, then publish it for concurrent queries. Query results are candidate
placements for subsequent validation, not completed alignments. Each thread
reuses its result list, so retain a copy before another query on that map/thread.
`SeqPos` and `SeqPosM` also borrow sequence arrays; their legacy equality methods
expect compatible, non-null sequence-position objects.

## Tests and known limits

Finite regression fixtures and fixed-work benchmarks live in
[`test/map`](../test/map/). They cover relocation return values, overflow,
packed-bit boundaries, sequence metadata and candidate retrieval, and map reuse.
The legacy `main` methods inside production classes often run large hardcoded
benchmarks; they are not quick smoke tests.

Source TODOs retain unresolved contracts and extreme-capacity limitations.
Passing the finite regressions does not establish correctness for every collection operation
or workload. Measure performance on the actual workload when changing a hot path.
