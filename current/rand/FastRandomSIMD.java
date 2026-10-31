package rand;

import shared.Random;

import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * SIMD-accelerated pseudorandom number generator implementing shared.Random.
 * Uses vectorized operations to generate multiple random values in parallel,
 * though benchmarks show it performs slower than non-SIMD implementations.
 * May be useful for filling large arrays where parallelism can be leveraged.
 * Based on xoroshiro128+ algorithm with vector operations.
 * Mutable and not thread-safe; callers should give each thread its own instance
 * or externally synchronize all access, including output-buffer and Gaussian
 * cache access.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 21, 2025
 */
public final class FastRandomSIMD implements shared.Random{

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained version marker; this class does not implement Serializable. */
	private static final long serialVersionUID=1L;

	/** Four-lane 256-bit long vector shape. */
	private static final VectorSpecies<Long> SPECIES=LongVector.SPECIES_256;

	/** Number of independent long lanes in SPECIES. */
	private static final int VLEN=4;//SPECIES.length(); // 4 longs per vector

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Eight state words: four seed0 lanes followed by four seed1 lanes. */
	private final long[] seedState=new long[8];

	/** Reused output buffer populated from the four vector lanes. */
	private final long[] buffer=new long[64];

	/** Next unread index in buffer; buffer.length means empty. */
	private int bufferPos=64; // Start empty

	/** Cached second Box-Muller Gaussian sample, or NaN when empty. */
	private double nextNextGaussian=Double.NaN;

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the requested high bits from the next long.
	 * @param bits Number of high bits to return, from 1 through 32
	 * @return Requested pseudorandom bits */
	protected int next(int bits){
		return (int)(nextLong()>>>(64-bits));
	}

	/** Returns a normally distributed value, caching the paired sample until
	 * the next call or a reseed.
	 * @return Pseudorandom Gaussian value with mean 0 and standard deviation 1 */
	@Override
	public double nextGaussian(){
		if(!Double.isNaN(nextNextGaussian)){
			final double temp=nextNextGaussian;
			nextNextGaussian=Double.NaN;
			return temp;
		}

		double v1, v2, s;
		do{
			v1=2*nextDouble()-1;
			v2=2*nextDouble()-1;
			s=v1*v1+v2*v2;
		}while(s>=1 || s==0);

		final double multiplier=Math.sqrt(-2*Math.log(s)/s);
		nextNextGaussian=v2*multiplier;
		return v1*multiplier;
	}

	/** Refills the output buffer from four separately advanced xoroshiro128+ lanes. */
	private void refillBuffer(){
		// Load seed state
		LongVector s0Vec=LongVector.fromArray(SPECIES, seedState, 0); // First 4 elements
		LongVector s1Vec=LongVector.fromArray(SPECIES, seedState, VLEN); // Last 4 elements

		// Process in chunks of 4 (vector width)
		for(int bufIdx=0; bufIdx<buffer.length; bufIdx+=VLEN){

			// Calculate results (s0 + s1)
			LongVector resultVec=s0Vec.add(s1Vec);

			// Update four independent xoroshiro128+ lanes; the scalar
			// transition is result=s0+s1, s1'=s1^s0,
			// s0'=rotl(s0,24)^s1'^(s1'<<16), s1''=rotl(s1',37).
			s1Vec=s1Vec.lanewise(VectorOperators.XOR, s0Vec);

			// s0' = rotateLeft(s0, 24) ^ s1' ^ (s1' << 16)
			s0Vec=s0Vec.lanewise(VectorOperators.ROL, 24);
			final LongVector s1Shifted=s1Vec.lanewise(VectorOperators.LSHL, 16);
			s0Vec=s0Vec.lanewise(VectorOperators.XOR, s1Vec)
					.lanewise(VectorOperators.XOR, s1Shifted);

			// s1' = rotateLeft(s1', 37)
			s1Vec=s1Vec.lanewise(VectorOperators.ROL, 37);

			// Store results to buffer and updated state back to seedState
			resultVec.intoArray(buffer, bufIdx);
		}
		s0Vec.intoArray(seedState, 0);
		s1Vec.intoArray(seedState, VLEN);
		bufferPos=0;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a FastRandomSIMD seeded from system time. */
	public FastRandomSIMD(){
		this(System.nanoTime());
	}

	/** Creates a FastRandomSIMD with the specified seed.
	 * @param seed Initial seed; negative seeds use system time */
	public FastRandomSIMD(long seed){
		setSeed(seed);
	}

	/** Mixes a seed via SplitMix64.
	 * @param x Seed to mix
	 * @return Mixed seed value */
	private static long mixSeed(long x){
		x+=0x9E3779B97F4A7C15L;
		x=(x^(x>>>30))*0xBF58476D1CE4E5B9L;
		x=(x^(x>>>27))*0x94D049BB133111EBL;
		return x^(x>>>31);
	}

	/** Returns the next buffered word, interleaving four generator lanes.
	 * @return Signed pseudorandom long */
	@Override
	public long nextLong(){
		if(bufferPos>=buffer.length){refillBuffer();}
		return buffer[bufferPos++];
	}

	/** Returns the low 32 bits of the next generated long.
	 * @return Signed pseudorandom int */
	@Override
	public int nextInt(){
		return (int)nextLong();
	}

	/** Draws an int from [0, bound), rejecting remainder bias.
	 * @param bound Positive exclusive upper bound
	 * @return Pseudorandom int in range
	 * @throws IllegalArgumentException If bound is nonpositive */
	@Override
	public int nextInt(int bound){
		if(bound<=0){
			throw new IllegalArgumentException("bound must be positive");
		}

		if((bound&(bound-1))==0){return (int)((bound*(nextLong()>>>33))>>>31);}

		int bits, val;
		do{
			bits=(int)(nextLong()>>>33);
			val=bits%bound;
		}while(bits-val+(bound-1)<0);

		return val;
	}

	/** Draws an int from [origin, bound), including widths larger than Integer.MAX_VALUE.
	 * @param origin Inclusive lower bound
	 * @param bound Exclusive upper bound, greater than origin
	 * @return Pseudorandom int in range
	 * @throws IllegalArgumentException If origin is not below bound */
	@Override
	public int nextInt(int origin, int bound){
		if(origin>=bound){throw new IllegalArgumentException("origin must be less than bound");}
		final int n=bound-origin;
		if(n>0){return origin+nextInt(n);}

		int r=nextInt();
		while(r<origin || r>=bound){
			r=nextInt();
		}
		return r;
	}

	/** Draws a long from [0, bound), rejecting remainder bias.
	 * @param bound Positive exclusive upper bound
	 * @return Pseudorandom long in range
	 * @throws IllegalArgumentException If bound is nonpositive */
	@Override
	public long nextLong(long bound){
		if(bound<=0){
			throw new IllegalArgumentException("bound must be positive");
		}

		if((bound&(bound-1))==0){return nextLong()&(bound-1);}

		long bits, val;
		do{
			bits=nextLong()>>>1;
			val=bits%bound;
		}while(bits-val+(bound-1)<0);

		return val;
	}

	/** Tests the low bit of the next generated long.
	 * @return Pseudorandom boolean */
	@Override
	public boolean nextBoolean(){
		return (nextLong()&1)!=0;
	}

	/** Scales the high 24 bits of the next long into [0, 1).
	 * @return Pseudorandom float less than one */
	@Override
	public float nextFloat(){
		return (nextLong()>>>40)*0x1.0p-24f;
	}

	//    @Override
	//    public float nextFloat(){//Not any faster
	//        return Float.intBitsToFloat((int)(0x3f800000 | (nextLong() & 0x7fffff))) - 1.0f;
	//    }

	/** Scales the high 53 bits of the next long into [0, 1).
	 * @return Pseudorandom double less than one */
	@Override
	public double nextDouble(){
		return (nextLong()>>>11)*0x1.0p-53d;
	}

	/** Fills bytes in little-endian chunks from successive nextLong() values.
	 * @param bytes Destination array */
	@Override
	public void nextBytes(byte[] bytes){
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

	/*--------------------------------------------------------------*/
	/*----------------          Reseeding           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Sets the seed of this random number generator.
	 * Uses successive SplitMix64-style mixes to initialize eight state words.
	 * Performs two full-buffer refills so the second refill overwrites the first
	 * 64 generated values before callers read from the buffer.
	 * @param seed Initial seed; negative requests the current system nanotime
	 */
	@Override
	public void setSeed(long seed){
		if(seedState==null){return;}
		if(seed<0){seed=System.nanoTime();}// [rand negative-seed] shared.Random uses negative seeds as system-time requests.
		nextNextGaussian=Double.NaN;
		// Use SplitMix64 to generate distinct seeds
		long currentSeed=seed;

		// Initialize first 4 elements (seed0)
		for(int i=0; i<4; i++){
			seedState[i]=currentSeed;
			currentSeed=mixSeed(currentSeed);
		}

		// Initialize last 4 elements (seed1)
		for(int i=4; i<8; i++){
			seedState[i]=currentSeed;
			currentSeed=mixSeed(currentSeed);
		}

		// Warm up by discarding the whole first 64-value buffer.
		refillBuffer();
		bufferPos=16; // Overwritten by refillBuffer below.
		refillBuffer();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Benchmark           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main method for benchmarking FastRandomSIMD against other PRNG implementations.
	 * Compares performance with FastRandom, shared.Random, and ThreadLocalRandom.
	 * @param args Command-line arguments; first argument sets iteration count (default 100M)
	 */
	public static void main(String[] args){
		int iterations=args.length>0 ? Integer.parseInt(args[0]) : 100_000_000;

		// Test FastRandom
		long startTime=System.nanoTime();
		Random fastRandom=new FastRandom();
		float sum=0;
		for(int i=0; i<iterations; i++){
			sum+=fastRandom.nextFloat();
		}
		long endTime=System.nanoTime();
		System.out.println("FastRandom time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);

		// Test FastRandomSIMD
		startTime=System.nanoTime();
		fastRandom=new FastRandomSIMD();
		sum=0;
		for(int i=0; i<iterations; i++){
			sum+=fastRandom.nextFloat();
		}
		endTime=System.nanoTime();
		System.out.println("FastRandomSIMD time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);
		//        for(int i=0; i<32; i++){System.err.println(i+": "+fastRandom.nextFloat());}

		// Test java.util.Random
		startTime=System.nanoTime();
		java.util.Random random=new java.util.Random();
		sum=0;
		for(int i=0; i<iterations; i++){
			sum+=random.nextFloat();
		}
		endTime=System.nanoTime();
		System.out.println("Random time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);

		// Test ThreadLocalRandom
		startTime=System.nanoTime();
		sum=0;
		java.util.concurrent.ThreadLocalRandom randy=java.util.concurrent.ThreadLocalRandom.current();
		for(int i=0; i<iterations; i++){
			sum+=randy.nextFloat();
		}
		endTime=System.nanoTime();
		System.out.println("ThreadLocalRandom time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);
	}
}
