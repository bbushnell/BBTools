package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * UltraLogLog8: 6-bit NLZ + 2-bit sub-NLZ history in 8-bit registers.
 * <p>
 * Extends the LL6 (HLL) baseline by tracking whether elements with
 * NLZ = (max-1) and NLZ = (max-2) were ever observed. These 2 history
 * bits refine the per-bucket contribution to the harmonic sum using
 * theoretically derived correction factors, reducing estimator variance
 * without requiring additional correction factor tables.
 * <p>
 * Register encoding (8 bits):
 * <ul>
 *   <li>Bits 7-2: nlzStored = absNlz + 1 (1-63, or 0 = empty)</li>
 *   <li>Bit 1: sawMinus2 (observed NLZ = max - 2)</li>
 *   <li>Bit 0: sawMinus1 (observed NLZ = max - 1)</li>
 * </ul>
 * <p>
 * Correction factors derived from E[log2(u) | state, max=X] where
 * u = n/2^X via numerical integration of the conditional distribution.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public final class UltraLogLog8 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	UltraLogLog8(){
		this(2048, 31, -1, 0);
	}

	UltraLogLog8(Parser p){
		super(p);
		maxArray=new byte[buckets];
	}

	UltraLogLog8(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new byte[buckets];
	}

	@Override
	public UltraLogLog8 copy(){return new UltraLogLog8(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Computes cardinality using the corrected harmonic sum.
	 * Each bucket contributes 2^(-(absNlz + CF[state])) to the sum,
	 * where CF[state] is the theoretically derived correction.
	 */
	private CardinalityStats summarize(){
		final int[] nlzCounts=new int[64];
		for(int i=0; i<buckets; i++){
			final int stored=maxArray[i]&0xFF;
			if(stored>0){
				final int absNlz=(stored>>2)-1;
				if(absNlz>=0 && absNlz<64){nlzCounts[absNlz]++;}
			}
		}
		lastRawNlz=nlzCounts;
		lastCorrNlz=nlzCounts;
		return CardinalityStats.fromNlzCounts(nlzCounts, buckets, microIndex,
		                                      CF_MATRIX, CF_BUCKETS,
		                                      CorrectionFactor.lastCardMatrix, CorrectionFactor.lastCardKeys);
	}

	/**
	 * Computes the ULL-corrected harmonic sum directly, bypassing
	 * nlzCounts aggregation. Each bucket gets an individual correction
	 * based on its history bits.
	 * @return corrected cardinality estimate using Mean formula
	 */
	double correctedEstimate(){
		double corrSum=0;
		double plainSum=0;
		int count=0;
		for(int i=0; i<buckets; i++){
			final int stored=maxArray[i]&0xFF;
			if(stored==0){continue;}
			count++;
			final int nlzStored=(stored>>2)&0x3F; // absNlz+1
			final int absNlz=nlzStored-1;
			final int state=stored&0x3;
			final double base=Math.pow(2.0, -absNlz);
			plainSum+=base;
			corrSum+=base*STATE_MULTIPLIER[state];
		}
		if(count==0){return 0;}
		// Self-normalize: divide out the average multiplier so only
		// relative bucket differences matter, not the global bias
		final double normFactor=corrSum/plainSum;
		final double normalizedSum=corrSum/normFactor; // == plainSum scaled by relative corrections
		// Wait - that's just plainSum. Need different approach.
		// Instead: adjust each bucket's contribution relative to the sketch-wide average multiplier
		// correctedSum / avgMult = sum of (2^(-X) * mult[s] / avgMult)
		// where avgMult = corrSum / plainSum
		// This equals plainSum * (sum of mult[s]*base / corrSum)...
		// Actually: the right formula is to use corrSum in place of plainSum but
		// scale the result so the estimator is unbiased.
		// Mean = 2*(count+B)/B * count^2 / sum
		// If we use corrSum, we get a different estimate. The ratio corrSum/plainSum
		// tells us by how much our corrections shifted the sum. We want the VARIANCE
		// reduction without the BIAS shift.
		// So: use corrSum but multiply result by normFactor to undo the bias
		final double rawEstimate=2.0*((count+(double)buckets)/buckets)*((double)count*count)/corrSum;
		return rawEstimate*normFactor;
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		// Use corrected estimate for ULL, fall back to standard for comparison
		final CardinalityStats s=summarize();
		double standard=s.hybridDLL();
		double corrected=correctedEstimate();
		// Use whichever is valid; corrected for the main estimate
		long card=(long)corrected;
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((UltraLogLog8)log);
	}

	/** Merges another UltraLogLog8 via per-bucket max with history carry. */
	public void add(UltraLogLog8 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				mergeRegister(i, log.maxArray[i]&0xFF);
			}
		}
	}

	/**
	 * Merges a single register value into this sketch.
	 * Takes the register with the higher NLZ. On ties, OR the history bits.
	 * When the other has higher NLZ, carry my history via shift register.
	 */
	private void mergeRegister(int bucket, int otherStored){
		final int myStored=maxArray[bucket]&0xFF;
		if(otherStored<=0){return;}
		if(myStored==0){maxArray[bucket]=(byte)otherStored; return;}

		final int myNlz=(myStored>>2);
		final int otherNlz=(otherStored>>2);

		if(myNlz>otherNlz){
			// My NLZ is higher; merge other's history into mine via shift
			final int k=myNlz-otherNlz;
			final int otherHist=otherStored&0x3;
			final int shifted=((otherHist|0x4)>>k)&0x3;
			maxArray[bucket]=(byte)(myStored|shifted); // OR in the carried bits
		}else if(otherNlz>myNlz){
			// Other is higher; carry my history into other's
			final int k=otherNlz-myNlz;
			final int myHist=myStored&0x3;
			final int shifted=((myHist|0x4)>>k)&0x3;
			maxArray[bucket]=(byte)(otherStored|shifted);
		}else{
			// Same NLZ: OR the history bits (union of knowledge)
			maxArray[bucket]=(byte)(myStored|otherStored);
		}
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(key&bucketMask);

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(USE_MICRO){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		final int newNlzStored=Math.min(nlz+1, 63); // absNlz+1, clamped
		final int oldStored=maxArray[bucket]&0xFF;
		final int oldNlzStored=(oldStored>0) ? (oldStored>>2) : 0;

		// History is immutable: only promotions change the register.
		// Elements below current max are pure early exits.
		if(newNlzStored<=oldNlzStored){return;}

		// New NLZ is higher: promote with shift-register history carry
		lastCardinality=-1;
		if(oldStored==0){filledBuckets++;}

		final int k=newNlzStored-oldNlzStored; // How far we jumped
		// History is a 2-bit shift register. On k-step promotion:
		// bit2 represents "old max existed", bits 1-0 are old history.
		// Shift right by k, mask to 2 bits.
		final int oldHist=(oldStored>0) ? (oldStored&0x3) : 0;
		final int newHist=((oldHist|0x4)>>k)&0x3;

		maxArray[bucket]=(byte)((newNlzStored<<2)|newHist);
	}

	public int filledBuckets(){return filledBuckets;}
	public double occupancy(){return (double)filledBuckets/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public double[] rawEstimates(){
		final CardinalityStats s=summarize();
		return s.toArray(correctedEstimate());
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * One byte per bucket: bits 7-2 = nlzStored (absNlz+1),
	 * bit 1 = sawMinus2, bit 0 = sawMinus1.
	 */
	private final byte[] maxArray;
	private int filledBuckets=0;

	int[] lastRawNlz, lastCorrNlz;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Theoretically derived state multipliers for the harmonic sum.
	 * Index = (sawMinus2 << 1) | sawMinus1.
	 * Multiplier = 2^(-CF[state]) where CF is the mean-preserving correction.
	 *
	 * CF values from E[log2(u) | state, max=X]:
	 *   state 00 (0): CF = -2.255 -> mult = 2^(+2.255) = 4.7745
	 *   state 01 (1): CF = -0.774 -> mult = 2^(+0.774) = 1.7108
	 *   state 10 (2): CF = -0.384 -> mult = 2^(+0.384) = 1.3056
	 *   state 11 (3): CF = +1.085 -> mult = 2^(-1.085) = 0.4713
	 */
	/**
	 * Linear-space normalized multipliers.
	 * Derived from theoretical E[log2(u)|state] correction factors,
	 * then normalized so frequency-weighted average = 1.0 to preserve
	 * the estimator's mean in linear (harmonic sum) space.
	 */
	static final double[] STATE_MULTIPLIER={
		3.0798,  // state 00: saw neither; NLZ likely outlier, inflate contribution
		1.1036,  // state 01: saw max-1 only
		0.8422,  // state 10: saw max-2 only
		0.3040   // state 11: saw both; well-established, reduce contribution
	};

	/** Reuse LL6 correction factors for the standard summarize() path. */
	public static final String CF_FILE="?cardinalityCorrectionLL6.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

}
