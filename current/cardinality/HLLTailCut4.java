package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * HLL-TailCut+ with 4-bit offset registers (K=16).
 * <p>
 * Based on Xiao, Zhou, Chen "Better with Fewer Bits" (IEEE INFOCOM 2017).
 * Uses 4-bit offset registers (K=16) plus a shared base register B.
 * The "TailCut" truncates the right tail of the register distribution
 * by capping offsets at K-1.  Estimation uses harmonic mean with
 * LinearCounting blend for small cardinalities.
 * <p>
 * Per XZC17, the harmonic mean estimator is correct (unbiased) with K=16.
 * The 3-bit variant (K=8) exhibits -5.2% bias with the same estimator.
 * This 4-bit version trades 33% more memory per register for correct
 * estimation without bias compensation.
 * <p>
 * <b>Online component (Algorithm 1 from paper):</b>
 * <ol>
 *   <li>Hash element, extract bucket (low bits) and rho (NLZ+1 of key).</li>
 *   <li>If rho - B &ge; K (overflow): advance B by min(all offsets), subtract from all.</li>
 *   <li>Store: offset[j] = max(offset[j], min(rho - B, K - 1)).</li>
 * </ol>
 * <p>
 * <b>Key invariants:</b>
 * <ul>
 *   <li>B = 0 and offset = 0 means empty (never received an element).</li>
 *   <li>B &gt; 0 implies ALL registers have been filled at least once
 *       (because B only advances when min(all offsets) &gt; 0).</li>
 *   <li>Absolute rho = B + offset for filled registers; NLZ = rho - 1.</li>
 * </ul>
 *
 * @author Chloe
 * @date April 2026
 */
public final class HLLTailCut4 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	HLLTailCut4(){
		this(2048, 31, -1, 0);
	}

	HLLTailCut4(Parser p){
		super(p);
		offsets=new byte[buckets];
	}

	HLLTailCut4(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		offsets=new byte[buckets];
	}

	@Override
	public HLLTailCut4 copy(){return new HLLTailCut4(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		final int bucket=(int)(key&bucketMask);
		final int nlz=Long.numberOfLeadingZeros(key);
		final int rho=nlz+1; // 1-indexed rank (paper's rho(x'))

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);

		// Algorithm 1, line 4: overflow detection
		if(rho-B>=K){
			// Lines 5-8: find minimum offset and advance base
			int deltaB=K; // sentinel; max offset is K-1
			for(int j=0; j<buckets; j++){
				final int off=offsets[j]&0xFF;
				if(off<deltaB){deltaB=off;}
			}
			if(deltaB>0){
				B+=deltaB;
				for(int j=0; j<buckets; j++){
					offsets[j]=(byte)((offsets[j]&0xFF)-deltaB);
				}
			}
		}

		// Algorithm 1, line 9: store
		final int rawOffset=rho-B;
		if(rawOffset<=0){return;} // at or below base
		final int newOffset=Math.min(rawOffset, K-1);
		final int oldOffset=offsets[bucket]&0xFF;
		if(newOffset<=oldOffset){return;}

		lastCardinality=-1;
		offsets[bucket]=(byte)newOffset;
	}

	/**
	 * Estimates cardinality using HLL harmonic mean formula adapted for
	 * base register + offset encoding (paper Eq. 9).
	 * With K=16 (4-bit registers), this estimator is correct per XZC17.
	 * Switches to LinearCounting for small cardinalities when empty registers exist.
	 */
	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}

		int V=0; // empty register count
		double harmonicSum=0;
		for(int j=0; j<buckets; j++){
			final int offset=offsets[j]&0xFF;
			if(B==0 && offset==0){
				V++;
				harmonicSum+=1.0; // 2^{-0} = 1
			}else{
				harmonicSum+=Math.pow(2.0, -(B+offset));
			}
		}

		// Harmonic mean estimate (Eq. 9): n = alpha * m^2 / sum(2^{-(B+offset_j)})
		final double alpha=(buckets>=128) ?
			0.7213/(1.0+1.079/buckets) :
			(buckets==64) ? 0.709 :
			(buckets==32) ? 0.697 : 0.673;
		double est=alpha*buckets*buckets/harmonicSum;

		// Small range correction: LinearCounting when estimate is small and empties exist
		if(est<=2.5*buckets && V>0){
			est=buckets*Math.log((double)buckets/V);
		}

		long card=(long)est;
		final int microCard=Long.bitCount(microIndex);
		card=Math.max(card, microCard);
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);

		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((HLLTailCut4)log);
	}

	/**
	 * Merges another HLLTailCut4 into this one.
	 * Reconstructs absolute rho values, takes per-register max,
	 * and re-encodes relative to the higher base.
	 */
	public void add(HLLTailCut4 other){
		added+=other.added;
		lastCardinality=-1;
		microIndex|=other.microIndex;

		if(offsets==other.offsets){return;}

		final int newB=Math.max(B, other.B);
		for(int j=0; j<buckets; j++){
			final int offA=offsets[j]&0xFF;
			final int offB=other.offsets[j]&0xFF;
			// Reconstruct absolute rho (0 = empty)
			final int rhoA=(B==0 && offA==0) ? 0 : B+offA;
			final int rhoB=(other.B==0 && offB==0) ? 0 : other.B+offB;
			final int maxRho=Math.max(rhoA, rhoB);

			if(maxRho==0){
				offsets[j]=0;
			}else{
				final int newOff=maxRho-newB;
				offsets[j]=(byte)Math.max(0, Math.min(newOff, K-1));
			}
		}
		B=newB;
	}

	public int filledBuckets(){
		int count=0;
		for(int j=0; j<buckets; j++){if(B>0 || (offsets[j]&0xFF)>0){count++;}}
		return count;
	}
	public double occupancy(){return (double)filledBuckets()/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public double[] rawEstimates(){
		final double est=(double)cardinality();
		final int total=17+AbstractCardStats.NUM_DLC_TIERS+AbstractCardStats.NUM_EXTRA;
		final double[] r=new double[total];
		java.util.Arrays.fill(r, est);
		return r;
	}

	/**
	 * Builds an NLZ histogram compatible with CardStats.
	 * Index 0 = empty buckets; index k (k&ge;1) = count with absNlz = k-1.
	 * @return NLZ histogram array of length 66
	 */
	public int[] buildNlzHistogram(){
		int[] nlzCounts=new int[66];
		for(int j=0; j<buckets; j++){
			final int offset=offsets[j]&0xFF;
			if(B==0 && offset==0){
				nlzCounts[0]++;
			}else{
				final int absNlz=B+offset-1; // rho = B+offset, nlz = rho-1
				if(absNlz>=0 && absNlz<64){
					nlzCounts[absNlz+1]++;
				}
			}
		}
		return nlzCounts;
	}

	/** Returns the current shared base register value. */
	public int getBaseRegister(){return B;}

	/** Returns the number of registers that are still empty (only meaningful when B==0). */
	public int emptyCount(){
		if(B>0){return 0;}
		int count=0;
		for(int j=0; j<buckets; j++){
			if(offsets[j]==0){count++;}
		}
		return count;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Offset registers: one byte per bucket, values in [0, K-1]. */
	private final byte[] offsets;

	/** Shared base register (paper's B). */
	private int B=0;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of bits per offset register. */
	public static final int REGISTER_BITS=4;

	/** Cutoff bound: K = 2^REGISTER_BITS.  Offset values in [0, K-1]. */
	public static final int K=1<<REGISTER_BITS; // 16

	/** Mask for a single register: 0xF for 4-bit registers. */
	public static final int REGISTER_MASK=0xF;

}
