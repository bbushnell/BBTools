# rand

The `rand` package contains `shared.Random` implementations used for deterministic
or high-throughput non-cryptographic pseudorandom numbers, plus an AES/CTR
generator for callers that need random data derived from AES keystream blocks.

## Generators

| Class | Core generator | Notes |
|---|---|---|
| `FastRandom` | xoroshiro128+ | Small scalar generator with deterministic positive seeds. |
| `FastRandomAES` | AES/CTR | Buffers AES keystream blocks as longs. |
| `FastRandomSIMD` | Four vectorized xoroshiro128+ lanes | Requires the Java Vector API and is opt-in. |
| `FastRandomXoshiro` | Xoshiro256+ | Small scalar generator with a faster power-of-two `nextInt(int)` path. |
| `JavaRandom` | `java.util.Random` | Preserves Java's signed-seed behavior while exposing `shared.Random`. |

## Contract notes

- Custom generators treat negative seeds as requests for a system-time seed,
  matching the `shared.Random` convention.
- `JavaRandom` intentionally delegates seeding to `java.util.Random`, so negative
  seeds are deterministic there.
- The custom generators are mutable and not thread-safe. Give each thread its own
  generator or externally synchronize all access, including Gaussian cache access.
- `JavaRandom` inherits java.util.Random's thread-safe core seed updates, but its
  added `nextLong(long)` and `nextInt(int,int)` methods draw repeatedly and are
  not atomic as a whole.
- Range methods return values in `[origin, bound)` and reject raw values when
  the positive width overflows.
- `nextGaussian()` caches the second Box-Muller sample and reseeding clears that
  cache.
