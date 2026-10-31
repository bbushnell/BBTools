package rand;

import shared.Random;

/**
 * Adapter exposing java.util.Random through the shared.Random interface.
 * Preserves java.util.Random signed-seed behavior, unlike shared.Random
 * implementations that treat negative seeds as system-time requests.
 * Inherited java.util.Random methods provide thread-safe seed updates, but the
 * added multi-draw range methods are not atomic as a whole.
 *
 * @author Brian Bushnell
 * @date January 1, 2026
 */
public final class JavaRandom extends java.util.Random implements Random{

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Version ID. */
	private static final long serialVersionUID=1L;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a JavaRandom with java.util.Random's default seed behavior. */
	public JavaRandom(){}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns a pseudorandom long in [0, bound).
	 * @param bound Positive upper bound, exclusive
	 * @return Pseudorandom value in range */
	@Override
	public long nextLong(long bound){
		if(bound<=0){throw new IllegalArgumentException("bound must be positive");}

		if((bound&(bound-1))==0){return nextLong()&(bound-1);}

		long bits, val;
		do{
			bits=nextLong()>>>1;
			val=bits%bound;
		}while(bits-val+(bound-1)<0);
		return val;
	}

	/** Returns a pseudorandom int in [origin, bound).
	 * @param origin Lower bound, inclusive
	 * @param bound Upper bound, exclusive and greater than origin
	 * @return Pseudorandom value in range */
	@Override
	public int nextInt(int origin, int bound){
		if(origin>=bound){throw new IllegalArgumentException("origin must be less than bound");}
		final int n=bound-origin;
		final int m=n-1;
		int r=nextInt();
		if((n&m)==0){return (r&m)+origin;}
		if(n>0){
			int u=r>>>1;
			while(u+m-(r=u%n)<0){u=nextInt()>>>1;}
			return r+origin;
		}

		while(r<origin || r>=bound){
			r=nextInt();
		}
		return r;
	}

}
