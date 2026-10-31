package rand;

import shared.Random;

/**
 * Fast, seedable xoroshiro128+ RNG for non-cryptographic use.
 * Implements shared.Random with two state words and deterministic nonnegative seeds.
 * Mutable and not thread-safe; callers should give each thread its own instance
 * or externally synchronize all access, including nextGaussian() cache access.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 14, 2025
 */
public final class FastRandom implements shared.Random{

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained version marker; this class does not implement Serializable. */
	private static final long serialVersionUID=1L;

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** xoroshiro128+ state. */
	private long seed0, seed1;

	/** Cached second Box-Muller Gaussian sample, or NaN when empty. */
	private double nextNextGaussian=Double.NaN;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a FastRandom seeded from system time. */
	public FastRandom(){
		this(System.nanoTime());
	}

	/** Creates a FastRandom with the specified seed (negative seeds use system time).
	 * @param seed Initial seed */
	public FastRandom(long seed){
		if(seed<0){seed=System.nanoTime();}

		// Initialize state using SplitMix64 algorithm
		seed0=seed;
		seed1=mixSeed(seed0);

		// Ensure we don't have all zeros (would lead to all zeros output)
		if(seed0==0 && seed1==0){
			seed0=0x5DEECE66DL;
			seed1=0L;
		}

		// Warm up the generator to avoid initial patterns
		for(int i=0; i<4; i++){
			nextLong();
		}
	}

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

	/** Returns the next pseudorandom long using xoroshiro128+ core.
	 * @return Next random long */
	@Override
	public long nextLong(){
		final long s0=seed0;
		long s1=seed1;
		final long result=s0+s1;

		s1^=s0;
		seed0=Long.rotateLeft(s0, 24)^s1^(s1<<16);
		seed1=Long.rotateLeft(s1, 37);

		return result;
	}

	/** Returns the low 32 bits of the next generated long.
	 * @return Signed pseudorandom int */
	@Override
	public int nextInt(){
		return (int)nextLong();
	}

	/**
	 * Returns a pseudorandom int in [0, bound), rejecting to avoid bias for non-powers of two.
	 * @param bound Upper bound (exclusive, positive)
	 * @return Random int in range
	 */
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

	/**
	 * Returns a pseudorandom int in [origin, bound).
	 * @param origin Lower bound (inclusive)
	 * @param bound Upper bound (exclusive), greater than origin
	 * @return Random int in range
	 */
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

	/**
	 * Returns a pseudorandom long in [0, bound), with bias-free rejection for non-powers of two.
	 * @param bound Upper bound (exclusive, positive)
	 * @return Random long in range
	 */
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

	/** Reinitializes RNG state using SplitMix64 seeds and a short warmup.
	 * @param seed New seed (negative uses system time) */
	@Override
	public void setSeed(long seed){
		if(seed<0){seed=System.nanoTime();}// [rand negative-seed] shared.Random uses negative seeds as system-time requests.
		nextNextGaussian=Double.NaN;
		seed0=seed;
		seed1=mixSeed(seed0);

		// Ensure we don't have all zeros
		if(seed0==0 && seed1==0){
			seed0=0x5DEECE66DL;
			seed1=0;
		}

		// Warm up the generator
		for(int i=0; i<4; i++){
			nextLong();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Mixes a seed via SplitMix64 to initialize RNG state.
	 * @param x Seed to mix
	 * @return Mixed seed value */
	private static long mixSeed(long x){
		x+=0x9E3779B97F4A7C15L;
		x=(x^(x>>>30))*0xBF58476D1CE4E5B9L;
		x=(x^(x>>>27))*0x94D049BB133111EBL;
		return x^(x>>>31);
	}

	/** Benchmarks FastRandom vs shared.Random and ThreadLocalRandom for nextFloat().
	 * @param args Optional iteration count (default 100,000,000) */
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
