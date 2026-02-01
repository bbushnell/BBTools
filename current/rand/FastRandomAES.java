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
 * AES-based cryptographically secure pseudorandom number generator using counter mode.
 * Uses SIMD operations but performance is slower than other PRNGs.
 * May be useful for filling large arrays with cryptographically secure random data.
 * Extends Random and provides all standard PRNG methods using AES encryption as the source.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 21, 2025
 */
public final class FastRandomAES implements Random {
	private static final long serialVersionUID = 1L;

	// Buffer size in longs (512 bytes = 64 longs)
	private static final int BUFFER_SIZE = 64;

	// Direct ByteBuffer for efficient native access
	private final ByteBuffer directBuffer;
	private final long[] longBuffer = new long[BUFFER_SIZE];
	private int bufferPos = BUFFER_SIZE; // Start empty

	// AES internals
	private final Cipher cipher;
	private final byte[] counterBlock;
	private final byte[] outputBlock;

	public FastRandomAES() {
		this(System.nanoTime());
	}

	public FastRandomAES(long seed) {
		// Initialize counter and output blocks
		counterBlock = new byte[BUFFER_SIZE * 8];
		outputBlock = new byte[BUFFER_SIZE * 8];

		// Create direct buffer for efficient access
		directBuffer = ByteBuffer.wrap(outputBlock).order(ByteOrder.nativeOrder());

		// Set up AES in counter mode
		try {
			cipher = Cipher.getInstance("AES/CTR/NoPadding");
		} catch (NoSuchAlgorithmException | NoSuchPaddingException e) {
			throw new RuntimeException(e);
		}
		setSeed(seed);
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

	private void refillBuffer() {
		try {
			// Encrypt counter blocks to generate random bytes
			cipher.update(counterBlock, 0, counterBlock.length, outputBlock);

			// Convert bytes to longs efficiently using ByteBuffer
			directBuffer.clear();
			for (int i = 0; i < BUFFER_SIZE; i++) {
				longBuffer[i] = directBuffer.getLong();
			}

			bufferPos = 0;
		} catch (Exception e) {
			throw new RuntimeException("AES encryption failed", e);
		}
	}

	protected int next(int bits) {
		return (int)(nextLong() >>> (64 - bits));
	}

	@Override
	public long nextLong() {
		if (bufferPos >= BUFFER_SIZE) {
			refillBuffer();
		}
		return longBuffer[bufferPos++];
	}

	//    /**
	//     * Creates a new FastRandom with the specified seed.
	//     * @param seed The initial seed
	//     */
	//    public FastRandomAES(long seed) {
	//        setSeed(seed);
	//    }

	private static long mixSeed(long x) {
		x += 0x9E3779B97F4A7C15L;
		x = (x ^ (x >>> 30)) * 0xBF58476D1CE4E5B9L;
		x = (x ^ (x >>> 27)) * 0x94D049BB133111EBL;
		return x ^ (x >>> 31);
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

	@Override
	public void setSeed(long seed) {
		if(cipher==null) {return;}
		try {
			// Use seed to derive key
			byte[] keyBytes = new byte[16];
			ByteBuffer.wrap(keyBytes).order(ByteOrder.LITTLE_ENDIAN)
			.putLong(0, seed).putLong(8, ~seed);
			SecretKeySpec key = new SecretKeySpec(keyBytes, "AES");

			cipher.init(Cipher.ENCRYPT_MODE, key, new IvParameterSpec(new byte[16]));

			// Warm up
			refillBuffer();
			bufferPos = BUFFER_SIZE; // Discard first set
			refillBuffer();

		} catch (Exception e) {
			throw new RuntimeException("Failed to initialize AES PRNG", e);
		}
	}

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
		fastRandom=new FastRandomAES();
		sum=0;
		for(int i=0; i<iterations; i++) {
			sum+=fastRandom.nextFloat();
		}
		endTime=System.nanoTime();
		System.out.println("FastRandomAES time: "+(endTime-startTime)/1_000_000+" ms, sum: "+sum);
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