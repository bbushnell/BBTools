package shared;

import java.util.Random;

/**
 * A fast, seedable random number generator based on xoshiro256+ algorithm.
 * Designed by David Blackman and Sebastiano Vigna as an improved successor to XorShift128+.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 2025
 */
public final class FastRandomXoshiro extends Random {

	private static final long serialVersionUID=1L;

	// State variables for xoshiro256+
	private long s0, s1, s2, s3;

	/**
	 * Creates a new FastRandomXoshiro with a random seed derived from system time.
	 */
	public FastRandomXoshiro() {
		this(System.nanoTime());
	}

	/**
	 * Creates a new FastRandomXoshiro with the specified seed.
	 * @param seed The initial seed
	 */
	public FastRandomXoshiro(long seed) {
		setSeed(seed>=0 ? seed : System.nanoTime());
	}

	/**
	 * Sets the seed of this random number generator.
	 */
	@Override
	public void setSeed(long seed) {
		// Use SplitMix64 to initialize the state (as recommended by the authors)
		s0=seed;
		s1=mixSeed(s0);
		s2=mixSeed(s1);
		s3=mixSeed(s2);

		// Ensure we don't have all zeros
		if(s0==0 && s1==0 && s2==0 && s3==0) {
			s0=0x5DEECE66DL;
			s1=0xBL;
			s2=0xCCAL;
			s3=0xF00L;
		}

		// Warm up the generator
		for(int i=0; i<4; i++) {
			nextLong();
		}
	}

	/**
	 * Mixes a seed value using SplitMix64 algorithm.
	 */
	private static long mixSeed(long x) {
		x+=0x9E3779B97F4A7C15L;
		x=(x^(x>>>30))*0xBF58476D1CE4E5B9L;
		x=(x^(x>>>27))*0x94D049BB133111EBL;
		return x^(x>>>31);
	}

	@Override
	protected int next(int bits) {
		return (int)(nextLong()>>>(64-bits));
	}

	/**
	 * Returns the next pseudorandom long value.
	 * This is the core generation function using xoshiro256+.
	 */
	@Override
	public long nextLong() {
		long result=s0+s3;

		long t=s1 << 17;

		s2^=s0;
		s3^=s1;
		s1^=s2;
		s0^=s3;

		s2^=t;
		s3=Long.rotateLeft(s3, 45);

		return result;
	}

	/**
	 * Returns a pseudorandom int value.
	 */
	@Override
	public int nextInt() {
		return (int)nextLong();
	}

	/**
	 * Returns a pseudorandom int value between 0 (inclusive) and bound (exclusive).
	 */
	@Override
	public int nextInt(int bound) {
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

	/**
	 * Returns a pseudorandom long value between 0 (inclusive) and bound (exclusive).
	 */
	@Override
	public long nextLong(long bound) {
		if(bound <= 0) {
			throw new IllegalArgumentException("bound must be positive");
		}

		// Fast path for powers of 2
		if((bound&(bound-1))==0) {
			return nextLong()&(bound-1);
		}

		// General case for any bound
		long bits, val;
		do {
			bits=nextLong()>>>1;
		val=bits % bound;
		} while(bits-val+(bound-1)<0); // Reject to avoid modulo bias

		return val;
	}

	/**
	 * Returns a pseudorandom boolean value.
	 */
	@Override
	public boolean nextBoolean() {
		return (nextLong()&1)!=0;
	}

	/**
	 * Returns a pseudorandom float value between 0.0 (inclusive) and 1.0 (exclusive).
	 */
	@Override
	public float nextFloat() {
		return (nextLong()>>>40)*0x1.0p-24f;
	}

	/**
	 * Returns a pseudorandom double value between 0.0 (inclusive) and 1.0 (exclusive).
	 */
	@Override
	public double nextDouble() {
		return (nextLong()>>>11)*0x1.0p-53d;
	}

	/**
	 * Fills the given array with random bytes.
	 */
	@Override
	public void nextBytes(byte[] bytes) {
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

	/**
	 * Main method for benchmarking against other PRNGs.
	 */
	public static void main(String[] args) {
		int iterations=args.length > 0 ? Integer.parseInt(args[0]) : 100_000_000;

		// Test FastRandom
		long startTime=System.nanoTime();
		Random fastRandom=new FastRandom();
		float sum=0;
		for(int i=0; i<iterations; i++) {
			sum+=fastRandom.nextFloat();
		}
		long endTime=System.nanoTime();
		System.out.println("FastRandom time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);

		// Test FastRandomXoshiro
		startTime=System.nanoTime();
		fastRandom=new FastRandomXoshiro();
		sum=0;
		for(int i=0; i<iterations; i++) {
			sum+=fastRandom.nextFloat();
		}
		endTime=System.nanoTime();
		System.out.println("FastRandomXoshiro time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);

		// Test java.util.Random
		startTime=System.nanoTime();
		java.util.Random random=new java.util.Random();
		sum=0;
		for(int i=0; i<iterations; i++) {
			sum+=random.nextFloat();
		}
		endTime=System.nanoTime();
		System.out.println("Random time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);

		// Test ThreadLocalRandom
		startTime=System.nanoTime();
		sum=0;
		Random randy=java.util.concurrent.ThreadLocalRandom.current();
		for(int i=0; i<iterations; i++) {
			sum+=randy.nextFloat();
		}
		endTime=System.nanoTime();
		System.out.println("ThreadLocalRandom time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);
	}
}

