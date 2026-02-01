package rand;

import shared.Random;

import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * SIMD-accelerated pseudorandom number generator extending shared.Random.
 * Uses vectorized operations to generate multiple random values in parallel,
 * though benchmarks show it performs slower than non-SIMD implementations.
 * May be useful for filling large arrays where parallelism can be leveraged.
 * Based on xoroshiro128+ algorithm with vector operations.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 21, 2025
 */
public final class FastRandomSIMD implements shared.Random {

	private static final long serialVersionUID = 1L;
	private static final VectorSpecies<Long> SPECIES = LongVector.SPECIES_256;
	private static final int VLEN = 4;//SPECIES.length(); // 4 longs per vector

	protected int next(int bits) {
		return (int)(nextLong() >>> (64 - bits));
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

	// State array: first 4 elements are seed0, last 4 are seed1
	private final long[] seedState = new long[8]; 
	// Output buffer
	private final long[] buffer = new long[64];
	private int bufferPos = 64; // Start empty

	private void refillBuffer() {
		// Load seed state
		LongVector s0Vec = LongVector.fromArray(SPECIES, seedState, 0); // First 4 elements
		LongVector s1Vec = LongVector.fromArray(SPECIES, seedState, VLEN); // Last 4 elements

		// Process in chunks of 4 (vector width)
		for(int bufIdx = 0; bufIdx < buffer.length; bufIdx += VLEN) {

			// Calculate results (s0 + s1)
			LongVector resultVec = s0Vec.add(s1Vec);

			// Update state
			s1Vec = s1Vec.lanewise(VectorOperators.XOR, s0Vec);

			// s0' = rotateLeft(s0, 24) ^ s1' ^ (s1' << 16)
			s0Vec = s0Vec.lanewise(VectorOperators.ROL, 24);
			LongVector s1Shifted = s1Vec.lanewise(VectorOperators.LSHL, 16);
			s0Vec = s0Vec.lanewise(VectorOperators.XOR, s1Vec)
				.lanewise(VectorOperators.XOR, s1Shifted);

			// s1' = rotateLeft(s1', 37)
			s1Vec = s1Vec.lanewise(VectorOperators.ROL, 37);

			// Store results to buffer and updated state back to seedState
			resultVec.intoArray(buffer, bufIdx);
		}
		s0Vec.intoArray(seedState, 0);
		s1Vec.intoArray(seedState, VLEN);
		bufferPos = 0;
	}

	public FastRandomSIMD() {
		this(System.nanoTime());
	}

	public FastRandomSIMD(long seed) {
		setSeed(seed);
	}

	private static long mixSeed(long x) {
		x += 0x9E3779B97F4A7C15L;
		x = (x ^ (x >>> 30)) * 0xBF58476D1CE4E5B9L;
		x = (x ^ (x >>> 27)) * 0x94D049BB133111EBL;
		return x ^ (x >>> 31);
	}

	@Override
	public long nextLong() {
		if(bufferPos>=buffer.length) {refillBuffer();}
		return buffer[bufferPos++];
	}

	@Override
	public int nextInt() {
		return (int)nextLong();
	}

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

	@Override
	public int nextInt(int origin, int bound) {
		if(origin>=bound) {
			throw new IllegalArgumentException("origin must be less than bound");
		}
		return origin + nextInt(bound-origin);
	}

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

	/**
	 * Sets the seed of this random number generator.
	 * Uses SplitMix64 to generate 8 distinct seed values for parallel PRNG streams.
	 * Performs warm-up generation to ensure proper state initialization.
	 * @param seed The initial seed value
	 */
	@Override
	public void setSeed(long seed) {
		if(seedState==null) {return;}
		// Use SplitMix64 to generate distinct seeds
		long currentSeed = seed;

		// Initialize first 4 elements (seed0)
		for(int i = 0; i < 4; i++) {
			seedState[i] = currentSeed;
			currentSeed = mixSeed(currentSeed);
		}

		// Initialize last 4 elements (seed1)
		for(int i = 4; i < 8; i++) {
			seedState[i] = currentSeed;
			currentSeed = mixSeed(currentSeed);
		}

		// Warm up
		refillBuffer();
		bufferPos = 16; // Discard first set
		refillBuffer();
	}

	/**
	 * Main method for benchmarking FastRandomSIMD against other PRNG implementations.
	 * Compares performance with FastRandom, shared.Random, and ThreadLocalRandom.
	 * @param args Command-line arguments; first argument sets iteration count (default 100M)
	 */
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

		// Test FastRandomSIMD
		startTime=System.nanoTime();
		fastRandom=new FastRandomSIMD();
		sum=0;
		for(int i=0; i<iterations; i++) {
			sum+=fastRandom.nextFloat();
		}
		endTime=System.nanoTime();
		System.out.println("FastRandomSIMD time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);
		//        for(int i=0; i<32; i++) {System.err.println(i+": "+fastRandom.nextFloat());}

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