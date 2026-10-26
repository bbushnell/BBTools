# Map regression fixtures and benchmarks

These Java entry points exercise the collections in `map`. Build the complete
`map` and `test.map` packages with the normal BBTools dependencies. Run correctness
fixtures with assertions enabled; the regression oracles also throw explicitly
under `-da`, which supports checking both assertion settings.

## Finite checks

| Entry point in `test.map` | Arguments for repaired code | Coverage |
|---|---|---|
| `LongIntMapIncrementTest` | `expect=fixed` | Return values across growth, signed values and mixed operations |
| `LongMapIncrementTest` | `expect=fixed` | Long-key counter growth and overflow controls |
| `IntMapIncrementTest` | `expect=fixed` | Four int-key counter variants and mixed-operation oracles |
| `ObjectIntMapIncrementTest` | `expect=fixed` | Absent counts, growth, merge and collision controls |
| `LongCounterOverflowTest` | `expect=fixed` | Unit and arbitrary increments, upper saturation and retained lower wrap; 64 boundary and 1,024 random pairs against BigInteger in three maps |
| `MapSetToMaxTest` | `expect=fixed` | Negative absent values, overlapping keys, growth, source preservation, self-merge and sentinel keys across ten maps |
| `IntHashMapLookupTest` | `expect=fixed` | Signed keys, sentinel changes, collision/wrap/removal and 20,000 mixed operations against HashMap |
| `IntHashSetListClearTest` | `expect=fixed` | Sparse/dense clearing and the external-removal restriction |
| `IntListHashMapClearTest` | `expect=fixed` | Collisions, wraparound, removal, sentinel insertion and clear/reuse |
| `LongArrayListHashMapTest` | None | Boxed-Long values, copied input lists, removal and clearing |
| `ListMapResizeTest` | `expect=fixed` | Copying insertion, no-copy adoption/replacement, growth-retained aliases, null/empty inputs, collision rehashing and recycled-list ownership in both long-key list maps |
| `LongLongHashMapHybridTest` | `expect=fixed` | Inline duplicates, side-list appends and clearing retained lists/holes |
| `BitSetCapacityTest` | `type=raw expect=fixed`, `type=multi expect=fixed`, or `type=atomic expect=fixed` | Same-capacity retention at 2^31 logical elements |
| `MapEdgeTest` | `case=NAME expect=fixed` in that order | The cases listed below |

`MapEdgeTest` cases are `double-count`, `render`, `bit-overflow`, `nan-order`,
`double-fetch`, `double-fetch-gc`, `gc-copy`, `null-equals`, and `terminal-capacity`.
Run each case separately. The terminal-capacity fixture uses a small real table
and reflected same-size resize requests to exercise the 95% no-growth policy
without allocating a table at the physical array ceiling.
The large bitset-capacity fixtures need more heap than the other finite checks;
the review used a 2 GiB heap. These are bounded in-memory tests, not file-based
bioinformatics workloads.

Where supported, `expect=baseline` records a particular historical defect. It is
not an alternative valid behavior for the repaired code. Baseline expectations
belong to the source revision that reproduces them.

## Performance measurements

The `*Benchmark` entry points use fixed work, warmup and consumed checksums.
Their `size`, `rounds`, `mode`, and map/type selectors control the workload; read
each entry point's argument parser before changing a recipe. Default rounds are
useful for sizing a run, not a promise of statistically useful timing duration.

Match source, host, JVM, heap, assertions and workload between revisions. Use
multiple repetitions and enough rounds for at least ten seconds of timed compute.
Record both wall and CPU time, and check correctness before interpreting speed.
The process includes setup and warmup; the reported internal timed region does not.
Fresh/growing-map modes intentionally include map allocation; existing-key modes
reuse their populated maps.

`ListMapResizeBenchmark type=array|long rounds=N keys=N values=N` repeatedly builds
maps from scalar insertions, exercising resize without copying caller lists at
insertion time. Defaults are 65,536 keys and eight values per key. Every stored
value and its order are checked. Stdout reports type, rounds, timed wall/CPU seconds
and a checksum; stderr reports main-thread allocated bytes for the timed region.
Allocation volume measures churn, not live heap or minimum required heap size.

`IntHashMapLookupBenchmark pattern=sequential` is the default: its existing keys
occupy their home cells. `pattern=mixed` generates signed odd existing keys and
disjoint even missing keys, checks that some existing keys are displaced from
their home cells, and prints that count to stderr. Mixed result labels start
with `mixed-`. Compare identical patterns and round counts between revisions;
a first-slot benchmark alone does not establish collision-path performance.

`IntListHashMapClearBenchmark mode=nocollide` measures populate/clear cycles where
both old and repaired versions do the required work. `mode=collision` is a repaired
correctness/stress control: the old clear can leave entries behind, so its shorter
time is not a valid performance baseline. The finite wraparound test prints the
occupied count after intentionally reinserting one key; it checks zero occupied
cells immediately after clearing, before that reinsertion.
