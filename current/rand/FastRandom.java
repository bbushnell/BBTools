package rand;

import shared.Random;

/**
 * Fast, seedable XorShift128+ RNG for non-cryptographic use.
 * Provides faster alternatives to shared.Random with deterministic seeding.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 14, 2025
 */
public final class FastRandom implements shared.Random {

	private static final long serialVersionUID = 1L;

	protected int next(int bits) {
		return (int)(nextLong() >>> (64 - bits));
	}

	private long seed0, seed1;

	/** Creates a FastRandom seeded from system time. */
	public FastRandom() {
		this(System.nanoTime());
	}

	/** Creates a FastRandom with the specified seed (negative seeds use system time).
	 * @param seed Initial seed */
	public FastRandom(long seed) {
		if(seed<0) {seed=System.nanoTime();}

		// Initialize state using SplitMix64 algorithm
		seed0 = seed;
		seed1 = mixSeed(seed0);

		// Ensure we don't have all zeros (would lead to all zeros output)
		if(seed0==0 && seed1==0) {
			seed0 = 0x5DEECE66DL;
			seed1 = 0L;
		}

		// Warm up the generator to avoid initial patterns
		for(int i=0; i<4; i++) {
			nextLong();
		}
	}
	


	// Default to NaN to indicate "empty"
	private double nextNextGaussian=Double.NaN;

	@Override
	public double nextGaussian() {
		// 1. Check if the cache is valid (not NaN)
		if(!Double.isNaN(nextNextGaussian)){
			double temp=nextNextGaussian;
			nextNextGaussian=Double.NaN; // Mark as empty
			return temp;
		}

		// 2. Generate two new numbers
		double v1, v2, s;
		do{
			v1=2*nextDouble()-1; 
			v2=2*nextDouble()-1; 
			s=v1*v1+v2*v2;
		}while(s>=1 || s==0);

		double multiplier=Math.sqrt(-2*Math.log(s)/s);

		// 3. Store the second one
		nextNextGaussian=v2*multiplier;

		// 4. Return the first one
		return v1*multiplier;
	}

	/**
	 * Mixes a seed via SplitMix64 to initialize RNG state.
	 * @param x Seed to mix
	 * @return Mixed seed value
	 */
	private static long mixSeed(long x) {
		x += 0x9E3779B97F4A7C15L;
		x = (x ^ (x >>> 30)) * 0xBF58476D1CE4E5B9L;
		x = (x ^ (x >>> 27)) * 0x94D049BB133111EBL;
		return x ^ (x >>> 31);
	}

	/** Returns the next pseudorandom long using XorShift128+ core.
	 * @return Next random long */
	@Override
	public long nextLong() {
		long s0 = seed0;
		long s1 = seed1;
		long result = s0 + s1;

		s1 ^= s0;
		seed0 = Long.rotateLeft(s0, 24) ^ s1 ^ (s1 << 16);
		seed1 = Long.rotateLeft(s1, 37);

		return result;
	}

	@Override
	public int nextInt() {
		return (int)nextLong();
	}

	/**
	 * Returns a pseudorandom int in [0, bound), rejecting to avoid bias for non-powers of two.
	 * @param bound Upper bound (exclusive, positive)
	 * @return Random int in range
	 */
	@Override
	public int nextInt(int bound) {
		if(bound<=0) {
			throw new IllegalArgumentException("bound must be positive");
		}

		// Fast path for powers of 2
		if((bound & (bound-1))==0) {
			return (int)((bound * (nextLong() >>> 33)) >>> 31);
		}

		// General case for any bound
		int bits, val;
		do {
			bits = (int)(nextLong() >>> 33);
			val = bits % bound;
		} while(bits-val+(bound-1)<0); // Reject to avoid modulo bias

		return val;
	}

	/**
	 * Returns a pseudorandom int in [origin, bound).
	 * @param origin Lower bound (inclusive)
	 * @param bound Upper bound (exclusive)
	 * @return Random int in range
	 */
	@Override
	public int nextInt(int origin, int bound) {
		if(origin>=bound) {
			throw new IllegalArgumentException("origin must be less than bound");
		}
		return origin + nextInt(bound-origin);
	}

	/**
	 * Returns a pseudorandom long in [0, bound), with bias-free rejection for non-powers of two.
	 * @param bound Upper bound (exclusive, positive)
	 * @return Random long in range
	 */
	@Override
	public long nextLong(long bound) {
		if(bound<=0) {
			throw new IllegalArgumentException("bound must be positive");
		}

		// Fast path for powers of 2
		if((bound & (bound-1))==0) {
			return nextLong() & (bound-1);
		}

		// General case for any bound
		long bits, val;
		do {
			bits = nextLong() >>> 1;
			val = bits % bound;
		} while(bits-val+(bound-1)<0); // Reject to avoid modulo bias

		return val;
	}

	@Override
	public boolean nextBoolean() {
		return (nextLong() & 1)!=0;
	}

	@Override
	public float nextFloat() {
		return (nextLong() >>> 40) * 0x1.0p-24f;
	}

	//    @Override
	//    public float nextFloat() {//Not any faster
	//        return Float.intBitsToFloat((int)(0x3f800000 | (nextLong() & 0x7fffff))) - 1.0f;
	//    }

	@Override
	public double nextDouble() {
		return (nextLong() >>> 11) * 0x1.0p-53d;
	}

	@Override
	public void nextBytes(byte[] bytes) {
		int i=0;
		int len=bytes.length;

		// Process 8 bytes at a time for efficiency
		while(i<len-7) {
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
		if(i<len) {
			long rnd=nextLong();
			do {
				bytes[i++]=(byte)rnd;
				rnd>>=8;
			} while(i<len);
		}
	}

	/** Reinitializes RNG state using SplitMix64 seeds and a short warmup.
	 * @param seed New seed (negative uses system time) */
	@Override
	public void setSeed(long seed) {
		seed0=seed;
		seed1=mixSeed(seed0);

		// Ensure we don't have all zeros
		if(seed0==0 && seed1==0) {
			seed0=0x5DEECE66DL;
			seed1=0;
		}

		// Warm up the generator
		for(int i=0; i<4; i++) {
			nextLong();
		}
	}

	/** Benchmarks FastRandom vs shared.Random and ThreadLocalRandom for nextFloat().
	 * @param args Optional iteration count (default 100,000,000) */
	public static void main(String[] args) {
		int iterations=args.length>0 ? Integer.parseInt(args[0]) : 100_000_000;

		// Test FastRandom
		long startTime=System.nanoTime();
		Random fastRandom=new FastRandom();
		float sum=0;
		for(int i=0; i<iterations; i++) {
			sum+=fastRandom.nextFloat();
		}
		long endTime=System.nanoTime();
		System.out.println("FastRandom time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);

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
		java.util.concurrent.ThreadLocalRandom randy=java.util.concurrent.ThreadLocalRandom.current();
		for(int i=0; i<iterations; i++) {
			sum+=randy.nextFloat();
		}
		endTime=System.nanoTime();
		System.out.println("ThreadLocalRandom time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);
	}
}