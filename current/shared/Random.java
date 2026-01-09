package shared;

import java.util.random.RandomGenerator;

import dna.AminoAcid;

/**
 * PRNG interface.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date January 6, 2026
 */
public interface Random extends RandomGenerator{

	/** Returns the next pseudorandom long using XorShift128+ core.
	 * @return Next random long */
	public long nextLong();
	
	public default int nextInt() {
		return (int)nextLong();
	}
	
	public default int nextInt(int bound) {
		assert(bound>=0) : "bound must be positive: "+bound;
		//if(bound <= 0) {throw new IllegalArgumentException("bound must be positive");}//Slow

		// Fast path for powers of 2
		if((bound&(bound-1))==0) {
			//return (int)((bound * (nextLong()>>>33))>>>31);//This looks dumb to me
			return ((int)nextLong())&(bound-1);
		}

		// General case for any bound
		int bits, val;
		do{
			bits=(int)(nextLong()>>>33);
			val=bits % bound;
		}while(bits-val+(bound-1)<0); // Reject to avoid modulo bias

		return val;
	}
	
	public default int nextInt2() {return ((int)nextLong())&1;}
	
	public default int nextInt3() {return (((int)nextLong())&Integer.MAX_VALUE)%3;}
	
	public default int nextInt4() {return ((int)nextLong())&3;}
	
	public default int nextInt5() {return (((int)nextLong())&Integer.MAX_VALUE)%5;}
	
	public default byte nextBase() {
		return AminoAcid.numberToBaseExtended[nextInt4()];
	}
	
	public default boolean nextBoolean() {
		return (nextLong()&1)!=0;
	}
	
	public default float nextFloat() {
		return (nextLong()>>>40)*0x1.0p-24f;
	}
	
	public default double nextDouble() {
		return (nextLong()>>>11)*0x1.0p-53d;
	}
	
	public default void nextBytes(byte[] bytes) {
		int i=0;
		int len=bytes.length;

		// Process 8 bytes at a time for efficiency
		while(i<len-7) {
			long rnd=nextLong();
			bytes[i++]=(byte)rnd;
			bytes[i++]=(byte)(rnd >> 8);
			bytes[i++]=(byte)(rnd >> 16);
			bytes[i++]=(byte)(rnd >> 24);
			bytes[i++]=(byte)(rnd >> 32);
			bytes[i++]=(byte)(rnd >> 40);
			bytes[i++]=(byte)(rnd >> 48);
			bytes[i++]=(byte)(rnd >> 56);
		}

		// Handle remaining bytes
		if(i<len) {
			long rnd=nextLong();
			do {
				bytes[i++]=(byte)rnd;
				rnd >>= 8;
			} while(i<len);
		}
	}

	/** Reinitializes RNG state.
	 * @param seed New seed (negative for a random seed) */
	public void setSeed(long seed);

}
