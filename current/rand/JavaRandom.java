package rand;

import shared.Random;

/**
 * Java RNG implementing Randy.
 *
 * @author Brian Bushnell
 * @date January 1, 2026
 */
public final class JavaRandom extends java.util.Random implements Random {

	private static final long serialVersionUID = 1L;

	//TODO: Uses Java 9+ library: java.util.Random.nextLong(long) and java.util.Random.nextInt(int,int) - both added in Java 17
	//(via the RandomGenerator interface). Violates the Java-8 bytecode-compliance target (no compiler warning, but the bytecode
	//needs Java 17 at runtime). All FastRandom* siblings implement these bound methods by hand and are Java-8-safe. Address post-compaction.
	public long nextLong(long x) {return super.nextLong(x);}
	public int nextInt(int origin, int bound) {return super.nextInt(origin, bound);}
}