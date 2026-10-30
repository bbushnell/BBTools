# SIMD helpers

This package supplies optional vector kernels behind `simd.Vector`, the public
facade used by BBTools callers. Prefer the facade over directly linking a backend:
it selects scalar implementations when vector support is unavailable or an input
is below the relevant length threshold. Individual method contracts still apply.

## Components

| Component | Role |
| --- | --- |
| `Vector` | Public dispatch, scalar implementations and capability probes. |
| `SIMD` | Fixed 256-bit numeric, neural-network and integer-table kernels. |
| `SIMDByte256` | Live byte transforms and scans, using 256-bit vectors and smaller tails. |
| `SIMDByte256_only` | Historical single-width reference; not the live byte backend. Retained caveats are documented in its source. |
| `SIMDLogLog` | Comparisons of packed DDL sketch buckets, including unsigned 16-bit values. |
| `SIMDAminoGlocal`, `SIMDMSA*` | Amino-acid glocal and packed-state nucleotide alignment helpers; callers must honor each helper's encoding and layout. |
| `SIMDAlign`, `SIMDAlignByte` | Alignment kernels and historical experiments. Some prototype methods have documented limitations; they are not interchangeable general-purpose aligners. |
| `QuantumSIMD`, `QuantumDenseTopSIMD` | Independent diagonal/up candidates for sparse and dense QuantumRanker alignment rows. Dependent frontier/deletion work stays in the caller. |
| `VarHandler` | Optional read-only access to private String storage for byte-builder appends. |

## Runtime and switches

The facade avoids direct newer-JDK API imports. Vector backends require the
incubating Vector API; use the BBTools shell launchers so the module and assertion
settings are supplied. `VarHandler` additionally requires a compatible String
layout and permission for the private lookup. An unavailable lookup leaves its
capability flag false.

`Shared.SIMD` is the master runtime switch. Common parser options are:

- `simd=auto`: select according to the capability probe.
- `simd=f`: disable vector dispatch.
- `simdfma=` and `simdff=`: control per-neuron dot products and batched feed-forward
  separately. Disabling batched feed-forward alone does not force every dot product
  to be scalar.
- `simdsparse=`: set both sparse switches; `simdmultsparse=` and `simdfmasparse=`
  set them separately.
- `simdcopy=`: control the optional vector copy path.

Forcing `simd=t` does not provide missing runtime support. The implemented general
backend dispatch uses 256-bit capability; the `simd128` and `simd64` probe results
do not select alternate full backends. Smaller loads and tails within a kernel
are separate implementation details.

## Numeric and ownership contracts

Floating-point vector reductions can group operations differently from scalar
loops. Fused multiplication can also change rounding. Do not require bit-identical
scalar/SIMD floating results unless the particular API promises it. SIMD float
sums retain float-vector accumulation for ordinary results and repeat in double
after a nonfinite result; their ordinary precision is still float, while scalar
sums accumulate in double. Finite inputs cannot overflow the returned double.
Integer sums widen before intermediate int overflow. Long sums retain Java's
wrapping arithmetic.

Generic indexed sparse operations follow their explicit index arrays; equal source
and destination lengths do not imply identity indexing. Specialized neural-network
callers may use dense operations when their validated representation establishes
that identity. Sparse block layout and index validity remain caller obligations.

Transforms mutate the arrays named in their method contracts. Scans and reductions
borrow their inputs. The caller must coordinate concurrent mutation and avoid
sharing mutable output buffers without synchronization. Runtime switches are
process-wide; change them outside concurrent work, not per record or per thread.

`iupacToN` masks ambiguity letters B/D/H/K/M/N/R/S/V/W/X/Y in either case to N.
It preserves A/C/G/T/U and their case, unrelated letters and nonletters.
Uppercasing, U-to-T conversion and junk repair are separate operations; enabling
IUPAC masking alone does not enable them. Read validation applies the requested
flags and its explicit junk policy independently.

## Verification

Regression and benchmark classes live in `test.simd`, under `current/test/simd`.
They exercise dispatch, numeric boundaries, byte handling, alignment boundaries,
sparse mappings and dependent model loading. Performance fixtures retain results
so timed work cannot simply disappear; their reported workload and copy costs
matter when interpreting results. A microbenchmark is not a whole-application
throughput claim. Compile the complete affected packages and check downstream
callers when changing a kernel.
