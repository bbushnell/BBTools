package shared;

import dna.AminoAcid;

/**
 * Java-8-compatible pseudorandom generator interface.
 * Defaults derive values from {@link #nextLong()}; concrete implementations may
 * override them and determine seed semantics, distribution and thread safety.
 * This interface makes no cryptographic-security or shared-instance concurrency guarantee.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date January 6, 2026
 */
public interface Random{//Can't extend RandomGenerator, not added until v17.

	/** Returns the next pseudorandom long using the implementation's generator.
	 * @return Next random long */
	public long nextLong();

	/** Returns the low 32 bits of the next long. */
	public default int nextInt(){return (int)nextLong();}

	/** Returns a value in [0,bound) for a positive bound.
	 * Uses low bits for powers of two and rejection sampling otherwise.
	 * @param bound Positive exclusive upper bound
	 * @return Pseudorandom nonnegative value below bound
	 * @throws IllegalArgumentException If bound is not positive */
	public default int nextInt(int bound){
		// FIXED shared/Random#001: zero used to enter the power-of-two path with mask -1.
		// Reject invalid bounds before consuming a random word, including when assertions are disabled.
		if(bound<=0){throw new IllegalArgumentException("bound must be positive: "+bound);}

		// Fast path for powers of 2
		if((bound&(bound-1))==0){
			//return (int)((bound * (nextLong()>>>33))>>>31);//This looks dumb to me
			return ((int)nextLong())&(bound-1);
		}

		// General case for any bound
		int bits, val;
		do{
			bits=(int)(nextLong()>>>33);
			val=bits%bound;
		}while(bits-val+(bound-1)<0); // Reject to avoid modulo bias

		return val;
	}

	/** Returns one low bit, in [0,2). */
	public default int nextInt2(){return ((int)nextLong())&1;}

	/** Returns a masked low-word remainder in [0,3); modulo bias is retained. */
	public default int nextInt3(){return (((int)nextLong())&Integer.MAX_VALUE)%3;}

	/** Returns two low bits, in [0,4). */
	public default int nextInt4(){return ((int)nextLong())&3;}

	/** Returns a masked low-word remainder in [0,5); modulo bias is retained. */
	public default int nextInt5(){return (((int)nextLong())&Integer.MAX_VALUE)%5;}

	/** Returns the nucleotide encoded by the next value from {@link #nextInt4()}. */
	public default byte nextBase(){return AminoAcid.numberToBaseExtended[nextInt4()];}

	/** Returns whether the next long's low bit is set. */
	public default boolean nextBoolean(){return (nextLong()&1)!=0;}

	/** Uses the high 24 bits of a long to produce a value in [0,1). */
	public default float nextFloat(){return (nextLong()>>>40)*0x1.0p-24f;}

	/** Uses the high 53 bits of a long to produce a value in [0,1). */
	public default double nextDouble(){return (nextLong()>>>11)*0x1.0p-53d;}

	/** Fills the caller's non-null array in place, consuming one long per eight bytes.
	 * Bytes are emitted low-byte first; a partial final group discards unused high bits.
	 * An empty array consumes no random values. */
	public default void nextBytes(byte[] bytes){
		int i=0;
		int len=bytes.length;

		// Process 8 bytes at a time for efficiency
		while(i<len-7){
			long rnd=nextLong();
			bytes[i++]=(byte)rnd;
			bytes[i++]=(byte)(rnd>>8);
			bytes[i++]=(byte)(rnd>>16);
			bytes[i++]=(byte)(rnd>>24);
			bytes[i++]=(byte)(rnd>>32);
			bytes[i++]=(byte)(rnd>>40);
			bytes[i++]=(byte)(rnd>>48);
			bytes[i++]=(byte)(rnd>>56);
		}

		// Handle remaining bytes
		if(i<len){
			long rnd=nextLong();
			do{
				bytes[i++]=(byte)rnd;
				rnd>>=8;
			}while(i<len);
		}
	}

	/** Reinitializes RNG state.
	 * @param seed New seed; interpretation of negative values is implementation-specific */
	public void setSeed(long seed);

	/** Returns a pseudorandom value in [0,l); implementations define invalid-bound handling.
	 * @param l Positive exclusive upper bound */
	public long nextLong(long l);

	/** Returns an implementation-generated standard Gaussian variate. */
	public double nextGaussian();

	/** Returns a value in [origin,bound), preserving the bounded path for narrow intervals.
	 * Wider intervals reject full signed-int draws outside the requested range;
	 * their width cannot be passed to the positive-int bounded method.
	 * @param origin Inclusive lower bound
	 * @param bound Exclusive upper bound, greater than origin
	 * @throws IllegalArgumentException If origin is not below bound */
	public default int nextInt(int origin, int bound){
		if(origin>=bound){
			throw new IllegalArgumentException("origin must be less than bound");
		}
		final int width=bound-origin;
		if(width>0){return origin+nextInt(width);}
		// A nonpositive int width here represents a valid interval wider than Integer.MAX_VALUE.
		int value=nextInt();
		while(value<origin || value>=bound){value=nextInt();}
		return value;
	}

}
