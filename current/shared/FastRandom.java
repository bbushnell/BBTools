package shared;

import java.util.Random;

/**
 * A fast, seedable random number generator for non-cryptographic purposes.
 * Based on XorShift128+ algorithm (Sebastiano Vigna and David Blackman) for speed and quality.
 * Safe for use in constructors and thread initialization.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 14, 2025
 */
public final class FastRandom extends java.util.Random {
	
    private static final long serialVersionUID = 1L;
    
    @Override
    protected int next(int bits) {
        return (int)(nextLong() >>> (64 - bits));
    }
    
    /** State variables for XorShift128+ */
    private long seed0, seed1;
    
    /**
     * Creates a new FastRandom with a random seed derived from system time.
     */
    public FastRandom() {
        this(System.nanoTime());
    }
    
    /**
     * Creates a new FastRandom with the specified seed.
     * @param seed The initial seed
     */
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
    
    /**
     * Mixes a seed value using SplitMix64 algorithm.
     */
    private static long mixSeed(long x) {
        x += 0x9E3779B97F4A7C15L;
        x = (x ^ (x >>> 30)) * 0xBF58476D1CE4E5B9L;
        x = (x ^ (x >>> 27)) * 0x94D049BB133111EBL;
        return x ^ (x >>> 31);
    }
    
    /**
     * Returns the next pseudorandom long value.
     * This is the core generation function using XorShift128+.
     */
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
     * Returns a pseudorandom int value between origin (inclusive) and bound (exclusive).
     */
    @Override
    public int nextInt(int origin, int bound) {
        if(origin>=bound) {
            throw new IllegalArgumentException("origin must be less than bound");
        }
        return origin + nextInt(bound-origin);
    }
    
    /**
     * Returns a pseudorandom long value between 0 (inclusive) and bound (exclusive).
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
    
    /**
     * Returns a pseudorandom boolean value.
     */
    @Override
    public boolean nextBoolean() {
        return (nextLong() & 1)!=0;
    }
    
    /**
     * Returns a pseudorandom float value between 0.0 (inclusive) and 1.0 (exclusive).
     */
    @Override
    public float nextFloat() {
        return (nextLong() >>> 40) * 0x1.0p-24f;
    }

//    @Override
//    public float nextFloat() {//Not any faster
//        return Float.intBitsToFloat((int)(0x3f800000 | (nextLong() & 0x7fffff))) - 1.0f;
//    }
    
    /**
     * Returns a pseudorandom double value between 0.0 (inclusive) and 1.0 (exclusive).
     */
    @Override
    public double nextDouble() {
        return (nextLong() >>> 11) * 0x1.0p-53d;
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
     */
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
    
    /**
     * Main method for benchmarking against other PRNGs.
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