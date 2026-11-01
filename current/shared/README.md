# BBTools `shared` package

Process-wide configuration and common utilities used by BBTools programs. This
package contains the classes below; parsing, concrete random generators and SIMD
kernels have their own packages.

| Class | Responsibility |
|---|---|
| [Shared](Shared.java) | Runtime configuration, thread/buffer defaults, environment information and sorting helpers. |
| [Tools](Tools.java) | Array, text, numeric, file-validation and read-processing helpers. |
| [TrimRead](TrimRead.java) | Quality/position trimming and retained state for a later untrim operation. |
| [ErrorRunMarker](ErrorRunMarker.java) | Marks supported low-depth k-mer runs as `N`, using the caller's error predicate. |
| [NameMapper](NameMapper.java) | Exact names and unambiguous first-whitespace-token aliases. |
| [FieldCopier](FieldCopier.java) | Reflective copying of compatible, same-named fields. |
| [Primes](Primes.java) | Sparse capacity-table lookup and a utility for filtering prime lists. |
| [Random](Random.java) | Java 8-compatible random-generator interface with default operations. |
| [Resources](Resources.java) | Resource lookup and download instructions for missing resources. |
| [KillSwitch](KillSwitch.java) | Fail-fast allocation/termination helpers and an optional process watchdog. |
| [SyncHeart](SyncHeart.java) | Locked access to selected `Shared` settings. |
| [MetadataWriter](MetadataWriter.java) | Asynchronous execution-metadata output in JSON or TSV. |
| [Timer](Timer.java) | Stored elapsed-time measurements and optional reporting. |
| [Colors](Colors.java) | ANSI escape strings and display palettes. |

## Ownership and concurrency

`Shared` settings are process-wide. Configure them before starting workers unless
the particular setting supports coordinated changes. `SyncHeart` protects only
accesses made through its own lock; direct reads and writes of `Shared` fields do
not participate. A thread setting also does not automatically limit every Java
library thread pool.

`Tools` methods have individual mutation and ownership contracts. An in-place
operation may replace array contents, while another helper returns a new array or
list. Operations on atomic-array elements do not automatically make an entire
multi-element operation atomic. Callers must coordinate concurrent mutations.

Trimming and marking require exclusive access to the affected read/bases/qualities.
`TrimRead`'s reversible form retains original arrays by reference and uses `Read.obj`
for its undo state. Preserve those original arrays until untrimming. Untrimming is
a one-shot operation: repeating it can extend alignment coordinates and clips again.
Static trimming paths do not create undo state. See individual methods for supported
SAM/match representations and alignment-state limitations.

`ErrorRunMarker` changes bases and available quality positions, while treating the
supplied depth counts as read-only. It does not reconstruct an error-corrected
sequence. Its predicate and arrays are owned by the caller.

`NameMapper` retains values by reference. Exact names take precedence over aliases;
ambiguous fallback matches throw. A lookup can record the first expanded query,
so concurrent lookups also require confinement or external synchronization.

`FieldCopier` copies object references shallowly when object copying is enabled.
Its default excludes static/final destination fields. Reflection/access failures
can produce a partial copy; the return value counts successful assignments. See
its superclass and shadowed-name rules before using it with inherited state.

## Resource and lifecycle behavior

`Primes` loads `?primes.txt.gz` once during class initialization. The table is sparse:
lookup does not promise the nearest mathematical prime. Its legacy beyond-table
product fallback is not a general primality algorithm. The public table must not
be mutated. Resource-read failures are reported instead of accepting partial data;
`makePrimes` closes its input and waits for a successfully started output writer,
preserving an earlier failure when cleanup also fails.

`Resources` delegates path discovery to `dna.Data`. Default lookup overloads may
print instructions and exit the process when a resource is missing. Optional
lookup can avoid exiting; a single-name miss still prints diagnostics, whereas
a candidate-list miss with `exit=false` is silent. `downloadURL` only resolves a
known URL and does not download a file.

`MetadataWriter.write` schedules output and returns without a per-write completion
handle. Returning from it does not mean that the metadata reached its destination.
Its reconstructed shell command is approximate and does not preserve quoting.

`KillSwitch` includes failure paths that halt the VM. Exercise termination behavior
in isolated child processes. `Timer` reports its stored elapsed value, updated by
`start`/`stop`; reading it is not a live clock. `Colors` captures its initial platform
settings when constructing escape strings; later changes do not rebuild palettes.

## Related packages

- [`parse`](../parse/README.md): `Parse`, `Parser`, `PreParser` and the `LineParser` families.
- [`rand`](../rand/README.md): concrete random generators and their seed/threading contracts.
- [`simd`](../simd/README.md): `Vector` dispatch and width-specific SIMD implementations.

The `shared.Random` interface does not promise cryptographic security or thread
safety. Those properties, and seed interpretation, belong to the implementation.

## Tests

Regression and timing entry points live in [`test.shared`](../test/shared/), apart
from production classes. They cover utility boundaries, trimming, reflection,
exceptional lock release, random defaults, prime-file errors, name resolution and
error-run marking. The separately retained `ErrorRunMarkerReference` is a test
oracle for the pre-gate algorithm.

File-based fixtures require a fresh `work=` directory. `PrimesFailureTest` isolates
each selected case in a child JVM, and uses Linux `/proc/self/fd` to check its own
input descriptors. Its injected writer case exercises borrowed stdout flushing,
not every possible filesystem failure.

Timing fixtures report a checked checksum and measured-loop nanoseconds. Compare
matching inputs, warmup, source dependencies, JVM/heap settings and host. Keep
process wall/CPU time as well: `RandomDefaultTimingTest` performs an independent
reference calculation outside the measured loop that still consumes process time.
