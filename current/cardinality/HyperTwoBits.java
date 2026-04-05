package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * HyperTwoBits: 2-bit-per-register cardinality estimator.
 * <p>
 * Port of the algorithm by Janson, Lumbroso, and Sedgewick (TCS 2025,
 * "Bit-array-based alternatives to HyperLogLog"). Uses 2M total bits
 * (2 bits per substream) with standard error ~1.46/sqrt(M).
 * <p>
 * Core idea: instead of recording the maximum NLZ per substream (like HLL/DLL),
 * record only whether a threshold T has been exceeded, using a 2-bit counter.
 * The three sketches of HyperBitBitBit (thresholds T, T+4, T+8) are encoded
 * compactly as a single 2-bit value per register:
 * <ul>
 *   <li>0 = none of the thresholds exceeded
 *   <li>1 = NLZ >= T seen
 *   <li>2 = NLZ >= T+4 seen
 *   <li>3 = NLZ >= T+8 seen
 * </ul>
 * When the fraction of nonzero entries exceeds FILL_THRESHOLD (98.8%),
 * T increments by 4 and all nonzero counters decrement by 1.
 * <p>
 * Estimation uses the LinearCounting-like formula:
 * {@code N = U * M * ln(1/beta)} where U = 2^T, beta = fraction of zeros,
 * and the correction factor c(beta) = sqrt(1/beta - 1) / ln(1/beta) gives
 * the relative standard error as c(beta)/sqrt(M).
 * <p>
 * Memory: 2 bits per register, stored as two long[] arrays (s1=MSB, s0=LSB).
 * For M=2048: 2*2048 = 4096 bits = 512 bytes + 6 bits for T.
 * <p>
 * <b>Limitations:</b>
 * <ul>
 *   <li>No LinearCounting fallback for small cardinalities — accuracy degrades
 *       when beta approaches 0 or 1 (c(beta) diverges).
 *   <li>No microIndex support (unlike DLL4).
 *   <li>Order-independent for distinct elements, but the threshold advancement
 *       mechanism means different T values at estimation time depending on
 *       arrival order when elements cause fills at different rates.
 * </ul>
 *
 * @author Eru (port), Janson/Lumbroso/Sedgewick (algorithm)
 * @date April 2026
 * @see <a href="https://doi.org/10.1016/j.tcs.2025.115450">TCS 2025</a>
 */
public final class HyperTwoBits extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	HyperTwoBits(){
		this(2048, 31, -1, 0);
	}

	HyperTwoBits(Parser p){
		super(p);
		s0=new long[wordsForBuckets(buckets)];
		s1=new long[wordsForBuckets(buckets)];
	}

	HyperTwoBits(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		s0=new long[wordsForBuckets(buckets)];
		s1=new long[wordsForBuckets(buckets)];
	}

	@Override
	public HyperTwoBits copy(){return new HyperTwoBits(buckets, k, -1, minProb);}

	private static int wordsForBuckets(int b){return (b+63)>>>6;}

	/*--------------------------------------------------------------*/
	/*----------------        2-bit Accessors       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reads the 2-bit value for register k from packed arrays s1 (MSB) and s0 (LSB).
	 * Value = 2*s1[k] + s0[k], range [0,3].
	 */
	private int get(int k){
		final int wordIdx=k>>>6;
		final int bitIdx=k&63;
		return (int)(((s1[wordIdx]>>>bitIdx)&1L)*2 + ((s0[wordIdx]>>>bitIdx)&1L));
	}

	/**
	 * Writes a 2-bit value for register k into packed arrays.
	 * v must be in [0,3].
	 */
	private void set(int k, int v){
		final int wordIdx=k>>>6;
		final long mask=1L<<(k&63);
		// Clear bit positions
		s1[wordIdx]=(s1[wordIdx]&~mask) | (((v>>>1)&1L)*mask);
		s0[wordIdx]=(s0[wordIdx]&~mask) | ((v&1L)*mask);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		final int bucket=(int)(key&bucketMask);
		final int nlz=Long.numberOfLeadingZeros(key);

		// Check thresholds T, T+4, T+8 and update 2-bit counter
		// Value 1 = nlz >= T, value 2 = nlz >= T+4, value 3 = nlz >= T+8
		if(nlz>=threshold){
			final int oldVal=get(bucket);
			int newVal;
			if(nlz>=threshold+8){
				newVal=3;
			}else if(nlz>=threshold+4){
				newVal=2;
			}else{
				newVal=1;
			}
			if(newVal>oldVal){
				if(oldVal==0){nonzeroCount++;}
				set(bucket, newVal);
				lastCardinality=-1;
			}
		}

		// Check fill threshold: when 98.8% of registers are nonzero, advance T
		if(nonzeroCount>=fillLimit){
			advanceThreshold();
		}
	}

	/**
	 * Advances threshold by 4 and decrements all nonzero counters.
	 * Uses bitwise operations on the packed long[] arrays for efficiency.
	 * <p>
	 * Decrement truth table (from paper Section 5):
	 * <pre>
	 *   before(s1,s0) → after(t1,t0)
	 *   0(0,0) → 0(0,0)   no change
	 *   1(0,1) → 0(0,0)   drops to zero
	 *   2(1,0) → 1(0,1)
	 *   3(1,1) → 2(1,0)
	 * </pre>
	 * Derived: t1 = s1 AND s0, t0 = s1 AND NOT s0.
	 * Then: s0 = t0 = s1 AND NOT s0, s1 = t1 = s1 AND s0.
	 * Order matters: compute new s0 first from old s1 and old s0.
	 */
	private void advanceThreshold(){
		threshold+=4;
		nonzeroCount=0;
		final int words=s0.length;
		for(int j=0; j<words; j++){
			final long oldS1=s1[j];
			final long oldS0=s0[j];
			s0[j]=oldS1&(~oldS0);  // t0 = s1 AND NOT s0
			s1[j]=oldS1&oldS0;     // t1 = s1 AND s0
		}
		// Recount nonzero entries
		for(int j=0; j<words; j++){
			nonzeroCount+=Long.bitCount(s0[j]|s1[j]);
		}
		lastCardinality=-1;
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}

		// beta = fraction of zeros in the sketch
		final double beta=1.0-(double)nonzeroCount/buckets;

		if(beta<=0 || beta>=1){
			// Degenerate: sketch fully occupied or fully empty
			// Fully empty → 0; fully occupied → estimate from last valid state
			if(beta<=0){
				// All registers nonzero. Best guess: scale by U*M.
				// This is outside the valid estimation range.
				lastCardinality=(long)((double)(1L<<threshold)*buckets);
			}else{
				lastCardinality=0;
			}
			return lastCardinality;
		}

		// N_hat = U * M * ln(1/beta) where U = 2^T
		final double U=Math.pow(2.0, threshold);
		final double estimate=U*buckets*Math.log(1.0/beta);
		long card=Math.max(0, (long)estimate);
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((HyperTwoBits)log);
	}

	/**
	 * Merges another HyperTwoBits into this one.
	 * <p>
	 * Merge strategy (from paper Section 7):
	 * <ol>
	 *   <li>If thresholds match: per-register max of 2-bit values.
	 *   <li>If thresholds differ by >= 8: use the one with higher T.
	 *   <li>If thresholds differ by 4: align by decrementing the higher-T
	 *       sketch's registers, then take per-register max.
	 * </ol>
	 * We generalize: bring both to the higher threshold by decrementing
	 * the lower-T sketch (d/4) times, then per-register max.
	 */
	public void add(HyperTwoBits log){
		added+=log.added;
		lastCardinality=-1;

		// Align thresholds by advancing the lower one
		// Each advance = decrement all nonzero by 1, threshold += 4
		// For simplicity, copy the other's arrays if we need to decrement them
		long[] otherS0, otherS1;
		int otherT=log.threshold;

		if(threshold==otherT){
			otherS0=log.s0;
			otherS1=log.s1;
		}else if(threshold>otherT){
			// Other is behind; decrement other's values (threshold-otherT)/4 times
			int steps=(threshold-otherT)/4;
			otherS0=log.s0.clone();
			otherS1=log.s1.clone();
			for(int s=0; s<steps; s++){
				for(int j=0; j<otherS0.length; j++){
					final long oS1=otherS1[j];
					final long oS0=otherS0[j];
					otherS0[j]=oS1&(~oS0);
					otherS1[j]=oS1&oS0;
				}
			}
		}else{
			// We are behind; advance ourselves
			int steps=(otherT-threshold)/4;
			for(int s=0; s<steps; s++){
				for(int j=0; j<s0.length; j++){
					final long oS1=s1[j];
					final long oS0=s0[j];
					s0[j]=oS1&(~oS0);
					s1[j]=oS1&oS0;
				}
			}
			threshold=otherT;
			otherS0=log.s0;
			otherS1=log.s1;
		}

		// Per-register max of 2-bit values: max(a,b) for 2-bit packed values.
		// For each bit position: result = (a > b) ? a : b.
		// With 2-bit encoding (s1=MSB, s0=LSB):
		// a > b when a1>b1, or (a1==b1 and a0>b0)
		// Simplification: just OR the arrays. This gives max(a,b) >= actual max,
		// but OR can turn (2,1) into 3 which is wrong.
		// Correct approach: per-register max via bit manipulation.
		// a_gt_b = a1 & ~b1  |  (a1 ~^ b1) & a0 & ~b0
		//        = a1 & ~b1  |  ~(a1^b1) & a0 & ~b0
		// result = a_gt_b ? a : b  (mux per bit)
		final int words=s0.length;
		for(int j=0; j<words; j++){
			final long a1=s1[j], a0=s0[j];
			final long b1=otherS1[j], b0=otherS0[j];
			// a > b when MSB of a is set and MSB of b isn't,
			// OR MSBs are equal and LSB of a is set and LSB of b isn't
			final long aGtB=(a1&~b1) | (~(a1^b1)&a0&~b0);
			// Mux: for each bit position, select a if aGtB, else b
			s1[j]=(a1&aGtB)|(b1&~aGtB);
			s0[j]=(a0&aGtB)|(b0&~aGtB);
		}

		// Recount
		nonzeroCount=0;
		for(int j=0; j<words; j++){
			nonzeroCount+=Long.bitCount(s0[j]|s1[j]);
		}

		// May need to advance if merged result exceeds fill
		while(nonzeroCount>=fillLimit && threshold<MAX_THRESHOLD){
			advanceThreshold();
		}
	}

	public int filledBuckets(){return nonzeroCount;}
	public double occupancy(){return (double)nonzeroCount/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public double[] rawEstimates(){
		final double card=(double)cardinality();
		// Return a minimal array compatible with the calibration driver
		return new double[]{card, card, card, card, card, card, card, card,
			card, card, card, card, card, card, card, card, card};
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Packed 2-bit counter arrays: s1[i]=MSB, s0[i]=LSB. M registers packed into ceil(M/64) longs. */
	private final long[] s0, s1;
	/** Current NLZ threshold. Starts at 1, increments by 4 when sketch fills. */
	private int threshold=1;
	/** Number of registers with nonzero value. */
	private int nonzeroCount=0;
	/** Fill limit: advance threshold when nonzeroCount >= this. 98.8% of buckets. */
	private final int fillLimit=(int)(0.988*buckets);

	public int getThreshold(){return threshold;}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum threshold before we stop advancing (NLZ can't exceed 64). */
	private static final int MAX_THRESHOLD=60;

}
