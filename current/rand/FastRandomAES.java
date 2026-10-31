package rand;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.security.NoSuchAlgorithmException;
import shared.Random;

import javax.crypto.Cipher;
import javax.crypto.NoSuchPaddingException;
import javax.crypto.spec.IvParameterSpec;
import javax.crypto.spec.SecretKeySpec;

/**
 * AES/CTR-based deterministic pseudorandom generator implementing shared.Random.
 * May be useful for filling large arrays from AES keystream bytes, but the
 * time-derived default seed and explicit 64-bit seed do not make it suitable
 * for secret generation.
 * Numeric sequences depend on ByteOrder.nativeOrder() when AES output bytes are
 * copied into longs.
 * Mutable and not thread-safe; callers should give each thread its own instance
 * or externally synchronize all access, including buffered values and Gaussian
 * cache access.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 21, 2025
 */
public final class FastRandomAES implements Random{

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained version marker; this class does not implement Serializable. */
	private static final long serialVersionUID=1L;

	/** Buffered longs per AES refill. */
	private static final int BUFFER_SIZE=64;

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Native-endian view of outputBlock for efficient long access. */
	private final ByteBuffer directBuffer;

	/** Reused copy of native-endian longs read from the current AES output block. */
	private final long[] longBuffer=new long[BUFFER_SIZE];

	/** Next unread index in longBuffer; BUFFER_SIZE means empty. */
	private int bufferPos=BUFFER_SIZE; // Start empty

	/** Cached second Box-Muller Gaussian sample, or NaN when empty. */
	private double nextNextGaussian=Double.NaN;

	/** AES/CTR cipher holding the internal counter state. */
	private final Cipher cipher;

	/** Zero plaintext encrypted to extract raw AES/CTR keystream bytes. */
	private final byte[] counterBlock;

	/** Reused AES/CTR keystream block. */
	private final byte[] outputBlock;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a FastRandomAES seeded from system time. */
	public FastRandomAES(){
		this(System.nanoTime());
	}

	/** Creates a FastRandomAES with the specified seed.
	 * @param seed Initial seed; negative seeds use system time */
	public FastRandomAES(long seed){
		// Initialize counter and output blocks
		counterBlock=new byte[BUFFER_SIZE*8];
		outputBlock=new byte[BUFFER_SIZE*8];

		// Create a ByteBuffer view for efficient access
		directBuffer=ByteBuffer.wrap(outputBlock).order(ByteOrder.nativeOrder());

		// Set up AES in counter mode
		try{
			cipher=Cipher.getInstance("AES/CTR/NoPadding");
		}catch(NoSuchAlgorithmException | NoSuchPaddingException e){
			throw new RuntimeException(e);
		}
		setSeed(seed);
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

	/** Refills the buffered random longs from the next AES/CTR keystream block. */
	private void refillBuffer(){
		try{
			// AES/CTR advances its internal counter on each update(); encrypting
			// zero bytes returns successive raw keystream blocks.
			cipher.update(counterBlock, 0, counterBlock.length, outputBlock);

			// Convert bytes to longs efficiently using ByteBuffer
			directBuffer.clear();
			for(int i=0; i<BUFFER_SIZE; i++){
				longBuffer[i]=directBuffer.getLong();
			}

			bufferPos=0;
		}catch(Exception e){
			throw new RuntimeException("AES encryption failed", e);
		}
	}

	/** Returns the requested high bits from the next long.
	 * @param bits Number of high bits to return, from 1 through 32
	 * @return Requested pseudorandom bits */
	protected int next(int bits){
		return (int)(nextLong()>>>(64-bits));
	}

	/** Returns the next native-endian buffered AES keystream word.
	 * @return Signed pseudorandom long */
	@Override
	public long nextLong(){
		if(bufferPos>=BUFFER_SIZE){
			refillBuffer();
		}
		return longBuffer[bufferPos++];
	}

	//    /**
	//     * Creates a new FastRandom with the specified seed.
	//     * @param seed The initial seed
	//     */
	//    public FastRandomAES(long seed){
	//        setSeed(seed);
	//    }

	/** Mixes a seed via SplitMix64.
	 * @param x Seed to mix
	 * @return Mixed seed value */
	private static long mixSeed(long x){
		x+=0x9E3779B97F4A7C15L;
		x=(x^(x>>>30))*0xBF58476D1CE4E5B9L;
		x=(x^(x>>>27))*0x94D049BB133111EBL;
		return x^(x>>>31);
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
		if(bound<=0){throw new IllegalArgumentException("bound must be positive");}

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
		if(bound<=0){throw new IllegalArgumentException("bound must be positive");}

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

	/** Reinitializes the AES key and discards the first keystream block.
	 * @param seed New seed; negative seeds use system time */
	@Override
	public void setSeed(long seed){
		if(cipher==null){return;}
		if(seed<0){seed=System.nanoTime();}// [rand negative-seed] shared.Random uses negative seeds as system-time requests.
		nextNextGaussian=Double.NaN;
		try{
			// Use seed to derive key
			byte[] keyBytes=new byte[16];
			ByteBuffer.wrap(keyBytes).order(ByteOrder.LITTLE_ENDIAN)
					.putLong(0, seed).putLong(8, ~seed);
			SecretKeySpec key=new SecretKeySpec(keyBytes, "AES");

			cipher.init(Cipher.ENCRYPT_MODE, key, new IvParameterSpec(new byte[16]));

			// Warm up
			refillBuffer();
			bufferPos=BUFFER_SIZE; // Discard first set
			refillBuffer();

		}catch(Exception e){
			throw new RuntimeException("Failed to initialize AES PRNG", e);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Benchmarks FastRandomAES against other random number generators.
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

		// Test FastRandomAES
		startTime=System.nanoTime();
		fastRandom=new FastRandomAES();
		sum=0;
		for(int i=0; i<iterations; i++){
			sum+=fastRandom.nextFloat();
		}
		endTime=System.nanoTime();
		System.out.println("FastRandomAES time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);
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
