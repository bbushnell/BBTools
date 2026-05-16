package ddl;

import rand.FastRandomXoshiro;
import shared.Tools;

/**
 * Empirical collision rate comparison: DDL vs SetSketch register encoding.
 * Both use 2048 buckets and 16-bit registers from the same random hashes.
 * SetSketch encodes register values as floor(-log_b(u)) where u is uniform.
 * DDL encodes as (NLZ << 10) | (~mantissa & 0x3FF).
 */
public class SetSketchCollisionTest {

	static final int BUCKETS = 2048;
	static final int BUCKET_BITS = 11;
	static final long BUCKET_MASK = BUCKETS - 1;
	static final int VALUE_BITS = 64 - BUCKET_BITS; // 53
	static final double VALUE_SCALE = 1.0 / (1L << VALUE_BITS);
	static final int MAX_REGISTER = 65534;

	static char[] buildDDL(int size, long seed) {
		char[] reg = new char[BUCKETS];
		final FastRandomXoshiro rng = new FastRandomXoshiro(seed);
		for (int i = 0; i < size; i++) {
			long hash = Tools.hash64shift(rng.nextLong());
			int bucket = (int)(hash & BUCKET_MASK);
			long remaining = hash >>> BUCKET_BITS;
			int score;
			if (remaining == 0) {
				score = (VALUE_BITS << 10);
			} else {
				int nlz = Long.numberOfLeadingZeros(remaining) - BUCKET_BITS;
				int shift = Long.SIZE - 1 - Long.numberOfLeadingZeros(remaining) - 10;
				int mantissa;
				if (shift >= 0) {
					mantissa = (int)((remaining >>> shift) & 0x3FF);
				} else {
					mantissa = (int)((remaining << -shift) & 0x3FF);
				}
				score = (nlz << 10) | (~mantissa & 0x3FF);
			}
			score = Math.max(score, 1);
			if (score > reg[bucket]) {
				reg[bucket] = (char)score;
			}
		}
		return reg;
	}

	static char[] buildSetSketch(int size, long seed, double logB) {
		char[] reg = new char[BUCKETS];
		final FastRandomXoshiro rng = new FastRandomXoshiro(seed);
		for (int i = 0; i < size; i++) {
			long hash = Tools.hash64shift(rng.nextLong());
			int bucket = (int)(hash & BUCKET_MASK);
			long remaining = hash >>> BUCKET_BITS;
			if (remaining <= 0) remaining = 1;
			double u = remaining * VALUE_SCALE;
			int score = (int)(-Math.log(u) / logB);
			score = Math.min(score, MAX_REGISTER);
			if (score > reg[bucket]) {
				reg[bucket] = (char)score;
			}
		}
		return reg;
	}

	static int countMatches(char[] a, char[] b) {
		int matches = 0;
		for (int i = 0; i < a.length; i++) {
			if (a[i] == b[i]) matches++;
		}
		return matches;
	}

	static void histogram(char[][] sketches, String label) {
		int[] counts = new int[MAX_REGISTER + 1];
		for (char[] s : sketches) {
			for (char v : s) {
				counts[v]++;
			}
		}
		int nonzero = 0;
		for (int c : counts) if (c > 0) nonzero++;
		System.out.println(label + ": " + nonzero + " distinct register values used");
	}

	public static void main(String[] args) {
		int n = 30;
		int size = 2_000_000;
		double b = 1.001;

		for (String arg : args) {
			String[] kv = arg.split("=");
			if (kv[0].equalsIgnoreCase("n")) n = Integer.parseInt(kv[1]);
			else if (kv[0].equalsIgnoreCase("size")) size = Integer.parseInt(kv[1]);
			else if (kv[0].equalsIgnoreCase("b")) b = Double.parseDouble(kv[1]);
		}

		double logB = Math.log(b);
		int pairs = n * (n - 1) / 2;

		System.out.println("SetSketch Collision Test");
		System.out.println("Sketches: " + n + ", Elements: " + size + ", Buckets: " + BUCKETS);
		System.out.println("SetSketch base b=" + b + ", log(b)=" + logB);
		System.out.println();

		System.out.print("Building sketches... ");
		long t0 = System.nanoTime();
		char[][] ddl = new char[n][];
		char[][] ss = new char[n][];
		for (int i = 0; i < n; i++) {
			long seed = (long)i * size * 3;
			ddl[i] = buildDDL(size, seed);
			ss[i] = buildSetSketch(size, seed, logB);
		}
		long t1 = System.nanoTime();
		System.out.printf("%.1f seconds\n", (t1 - t0) / 1e9);

		histogram(ddl, "DDL");
		histogram(ss, "SetSketch");

		System.out.println();
		System.out.println("Comparing " + pairs + " pairs...");

		long ddlTotal = 0, ssTotal = 0;
		int ddlMax = 0, ssMax = 0;

		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				int dm = countMatches(ddl[i], ddl[j]);
				int sm = countMatches(ss[i], ss[j]);
				ddlTotal += dm;
				ssTotal += sm;
				ddlMax = Math.max(ddlMax, dm);
				ssMax = Math.max(ssMax, sm);
			}
		}

		System.out.printf("\nDDL:       avg %.2f false matches/pair, max %d\n",
				(double)ddlTotal / pairs, ddlMax);
		System.out.printf("SetSketch: avg %.2f false matches/pair, max %d\n",
				(double)ssTotal / pairs, ssMax);
		System.out.printf("\nRatio (SS/DDL): %.2fx\n", (double)ssTotal / Math.max(ddlTotal, 1));
	}
}
