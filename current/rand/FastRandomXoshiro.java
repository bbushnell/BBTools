package rand;

import shared.Random;

/**
 * Fast, seedable Xoshiro256+ RNG for non-cryptographic use.
 * Implements the shared.Random interface with deterministic positive seeds.
 * Mutable and not thread-safe; callers should give each thread its own instance
 * or externally synchronize all access, including nextGaussian() cache access.
 *
 * @author Brian Bushnell
 */
public final class FastRandomXoshiro implements Random{

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained version marker; this class does not implement Serializable. */
	private static final long serialVersionUID=1L;

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** State variables for xoshiro256+. */
	private long s0, s1, s2, s3;

	/** Cached second Box-Muller Gaussian sample, or NaN when empty. */
	private double nextNextGaussian=Double.NaN;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a FastRandomXoshiro seeded from system time. */
	public FastRandomXoshiro(){
		this(System.nanoTime());
	}

	/** Creates a FastRandomXoshiro with the specified seed.
	 * @param seed Initial seed; negative seeds use system time */
	public FastRandomXoshiro(long seed){
		setSeed(seed>=0 ? seed : System.nanoTime());
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

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

	/** Replaces all four state words, clears the Gaussian cache and advances four draws.
	 * @param seed Seed; negative requests the current system nanotime */
	@Override
	public void setSeed(long seed){
		if(seed<0){seed=System.nanoTime();}// [rand negative-seed] shared.Random uses negative seeds as system-time requests.
		nextNextGaussian=Double.NaN;
		// Use SplitMix64 to initialize the state (as recommended by the authors)
		s0=seed;
		s1=mixSeed(s0);
		s2=mixSeed(s1);
		s3=mixSeed(s2);

		// Ensure we don't have all zeros
		if(s0==0 && s1==0 && s2==0 && s3==0){
			s0=0x5DEECE66DL;
			s1=0xBL;
			s2=0xCCAL;
			s3=0xF00L;
		}

		// Warm up the generator
		for(int i=0; i<4; i++){
			nextLong();
		}
	}

	/**
	 * Mixes a seed value using SplitMix64 algorithm.
	 * Applies avalanche function to improve distribution quality.
	 * @param x Seed value to mix
	 * @return Mixed seed value with better distribution properties
	 */
	private static long mixSeed(long x){
		x+=0x9E3779B97F4A7C15L;
		x=(x^(x>>>30))*0xBF58476D1CE4E5B9L;
		x=(x^(x>>>27))*0x94D049BB133111EBL;
		return x^(x>>>31);
	}

	/** Returns the requested high bits from the next long.
	 * @param bits Number of high bits to return, from 1 through 32
	 * @return Requested pseudorandom bits */
	protected int next(int bits){
		return (int)(nextLong()>>>(64-bits));
	}

	/** Advances the xoshiro256+ state and returns its sum output.
	 * @return Signed pseudorandom long */
	@Override
	public long nextLong(){
		final long result=s0+s3;

		final long t=s1<<17;

		s2^=s0;
		s3^=s1;
		s1^=s2;
		s0^=s3;

		s2^=t;
		s3=Long.rotateLeft(s3, 45);

		return result;
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
	 * @throws AssertionError If bound is nonpositive and assertions are enabled */
	@Override
	public int nextInt(int bound){
		// [rand/FastRandomXoshiro#001] Zero would mask with -1 and return an unrestricted int.
		assert(bound>0) : "bound must be positive: "+bound;
		//if(bound<=0){throw new IllegalArgumentException("bound must be positive");}//Slow

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

		// Fast path for powers of 2
		if((bound&(bound-1))==0){
			return nextLong()&(bound-1);
		}

		// General case for any bound
		long bits, val;
		do{
			bits=nextLong()>>>1;
			val=bits%bound;
		}while(bits-val+(bound-1)<0); // Reject to avoid modulo bias

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
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Benchmarks FastRandomXoshiro against other random number generators.
	 * Tests performance of FastRandom, FastRandomXoshiro, shared.Random,
	 * and ThreadLocalRandom by generating floats and measuring execution time.
	 * @param args Optional number of iterations (default: 100,000,000)
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

		// Test FastRandomXoshiro
		startTime=System.nanoTime();
		fastRandom=new FastRandomXoshiro();
		sum=0;
		for(int i=0; i<iterations; i++){
			sum+=fastRandom.nextFloat();
		}
		endTime=System.nanoTime();
		System.out.println("FastRandomXoshiro time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);

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

