package fun;

import java.util.Arrays;
import java.util.PriorityQueue;

/**
 * Finds all primes up to a given number using rolling wheel bitmasks and a heap.
 *
 * Commander's algorithm:
 * - Odd-only: position i = number 2*i + 3
 * - For each prime p, a "wheel" of p longs. The pattern of p's multiples in
 *   64-position chunks repeats every p longs. Apply via AND — one op per long.
 * - Small primes combined into product wheels (3*5*7*11 = 1155, 13*17 = 221)
 * - Medium primes get individual wheels (19..61, p longs each)
 * - Large primes tracked via min-heap: entry = (next_odd_composite, prime)
 *   Only relevant starting at prime^2. Checked per-candidate in surviving bits.
 * - Primes under sqrt(max) NOT stored in the sieve data structure.
 */
public class FindPrimes {

    public static void main(String[] args) {
        if (args.length == 0) {
            System.err.println("Usage: FindPrimes <max>");
            return;
        }
        long max = Long.parseLong(args[0]);
        if (max > 1000000) sieve(1000000); // JIT warmup
        long t = System.nanoTime();
        long count = sieve(max);
        double sec = (System.nanoTime() - t) / 1e9;
        System.out.printf("Found %,d primes up to %,d in %.3f seconds.%n", count, max, sec);
    }

    /** Build a wheel mask array for a single prime p (odd-only representation).
     *  Returns p longs. Wheel[k] has bit j CLEAR if position (64*k + j) is a multiple of p. */
    static long[] buildWheel(int p) {
        long[] w = new long[p];
        Arrays.fill(w, ~0L);
        int firstPos = (p - 3) / 2; // position of p itself (which IS a multiple of p)
        // Mark all positions ≡ firstPos (mod p)
        for (int k = 0; k < p; k++) {
            for (int j = 0; j < 64; j++) {
                if (((64 * k + j - firstPos) % p + p) % p == 0) {
                    w[k] &= ~(1L << j);
                }
            }
        }
        return w;
    }

    /** Build a combined wheel for multiple primes. Period = product of all primes. */
    static long[] buildCombinedWheel(int... primes) {
        int period = 1;
        for (int p : primes) period *= p;
        long[] w = new long[period];
        Arrays.fill(w, ~0L);
        for (int p : primes) {
            int firstPos = (p - 3) / 2;
            for (int k = 0; k < period; k++) {
                for (int j = 0; j < 64; j++) {
                    if (((64 * k + j - firstPos) % p + p) % p == 0) {
                        w[k] &= ~(1L << j);
                    }
                }
            }
        }
        return w;
    }

    static long sieve(long max) {
        if (max < 2) return 0;
        if (max < 3) return 1;

        // Wheels — the rolling cookie presses
        long[] wheelA = buildCombinedWheel(3, 5, 7, 11);  // period 1155 longs, ~9KB
        long[] wheelB = buildCombinedWheel(13, 17);         // period 221 longs, ~1.7KB
        // Individual wheels for primes 19..61
        int[] indivPrimes = {19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};
        long[][] indivWheels = new long[indivPrimes.length][];
        for (int i = 0; i < indivPrimes.length; i++) {
            indivWheels[i] = buildWheel(indivPrimes[i]);
        }

        // All wheel primes — counted separately, not in sieve
        int[] allWheelPrimes = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};

        // Small sieve to find primes 67..sqrt(max) for direct marking
        int sqrtMax = (int) Math.sqrt((double) max);
        while ((long)(sqrtMax + 1) * (sqrtMax + 1) <= max) sqrtMax++;

        int sLimit = Math.max((sqrtMax - 3) / 2, 0);
        byte[] sc = new byte[sLimit + 1];
        for (int i = 0; i <= sLimit; i++) {
            if (sc[i] == 0) {
                long p = 2L * i + 3;
                long s = (p * p - 3) / 2;
                if (s > sLimit) break;
                for (long j = s; j <= sLimit; j += p) sc[(int)j] = 1;
            }
        }
        // Collect primes > 61 for direct sieving
        int np = 0;
        for (int i = 0; i <= sLimit; i++) {
            if (sc[i] == 0 && (2 * i + 3) > 61) np++;
        }
        int[] bigPrimes = new int[np];
        for (int i = 0, k = 0; i <= sLimit; i++) {
            if (sc[i] == 0 && (2 * i + 3) > 61) bigPrimes[k++] = 2 * i + 3;
        }

        // Segmented sieve
        int BATCH = 1 << 12; // 4096 longs = 32KB, fits L1
        long[] seg = new long[BATCH];

        long maxOdd = (max % 2 == 0) ? max - 1 : max;
        long maxIdx = (maxOdd < 3) ? -1 : (maxOdd - 3) / 2;

        long count = 1; // prime 2
        for (int p : allWheelPrimes) if (p <= max) count++;

        // Track starting positions for big primes across segments
        long[] nextHit = new long[bigPrimes.length];
        for (int i = 0; i < bigPrimes.length; i++) {
            nextHit[i] = ((long)bigPrimes[i] * bigPrimes[i] - 3) / 2;
        }

        long baseLong = 0; // which global long index we're at

        for (long base = 0; base <= maxIdx; base += (long)BATCH * 64) {
            long top = Math.min(base + (long)BATCH * 64 - 1, maxIdx);
            int len = (int)(top - base + 1); // positions in this segment
            int nL = (len + 63) >>> 6;

            // Roll the wheels — AND per long, no bit manipulation
            int offA = (int)(baseLong % 1155);
            int offB = (int)(baseLong % 221);
            int[] offI = new int[indivPrimes.length];
            for (int w = 0; w < indivPrimes.length; w++) {
                offI[w] = (int)(baseLong % indivPrimes[w]);
            }

            for (int i = 0; i < nL; i++) {
                long val = wheelA[offA] & wheelB[offB];
                for (int w = 0; w < indivPrimes.length; w++) {
                    val &= indivWheels[w][offI[w]];
                    if (++offI[w] >= indivPrimes[w]) offI[w] = 0;
                }
                seg[i] = val;
                if (++offA >= 1155) offA = 0;
                if (++offB >= 221) offB = 0;
            }

            // Clear excess bits in last long
            int rem = len & 63;
            if (rem > 0) seg[nL - 1] &= (1L << rem) - 1;

            // Direct marking for primes > 61
            for (int pi = 0; pi < bigPrimes.length; pi++) {
                long idx = nextHit[pi];
                if (idx > top) continue;
                int p = bigPrimes[pi];
                int local;
                if (idx >= base) {
                    local = (int)(idx - base);
                } else {
                    long diff = base - idx;
                    long steps = (diff + p - 1) / p;
                    local = (int)(idx + steps * (long)p - base);
                }
                while (local < len) {
                    seg[local >>> 6] &= ~(1L << (local & 63));
                    local += p;
                }
                nextHit[pi] = base + local;
            }

            // Count survivors
            for (int j = 0; j < nL; j++) count += Long.bitCount(seg[j]);

            baseLong += nL;
        }

        return count;
    }
}
