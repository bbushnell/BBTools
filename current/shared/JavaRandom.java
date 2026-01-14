package shared;

/**
 * Java RNG implementing Randy.
 *
 * @author Brian Bushnell
 * @date January 1, 2026
 */
public final class JavaRandom extends java.util.Random implements Random {

	private static final long serialVersionUID = 1L;

	public long nextLong(long x) {return super.nextLong(x);}
	public int nextInt(int origin, int bound) {return super.nextInt(origin, bound);}
}