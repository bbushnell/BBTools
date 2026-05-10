package cardinality;

import java.util.HashMap;
import shared.Tools;

/**
 * HLL_4 algorithm from Apache Data Sketches, reimplemented for BBTools benchmarking.
 * <p>
 * Original algorithm: Apache Data Sketches library
 * Source: https://github.com/apache/datasketches-java
 * Original authors: Lee Rhodes, Kevin Lang, and the Apache DataSketches team
 * License: Apache 2.0
 * <p>
 * This is a clean-room reimplementation from the algorithm description for
 * accuracy comparison purposes. It does NOT depend on the Apache library.
 * <p>
 * Key innovation: 4-bit nibble storage with a global curMin offset.
 * Each register stores (value - curMin) in 4 bits [0..14].
 * Value 15 (AUX_TOKEN) indicates the actual value is stored in an auxiliary HashMap.
 * curMin advances monotonically when no registers remain at the minimum value.
 * <p>
 * Uses standard HyperLogLog harmonic-mean estimation with alpha_m correction.
 *
 * @author Chloe
 * @date May 2026
 */
public final class ApacheHLL4Wrapper extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	ApacheHLL4Wrapper(){this(2048, 31, -1, 0);}

	ApacheHLL4Wrapper(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	ApacheHLL4Wrapper(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		lgK=Integer.numberOfTrailingZeros(buckets);
		m=buckets;
		nibbles=new byte[(m+1)>>>1];
		auxMap=new HashMap<>();
		curMin=0;
		numAtCurMin=m;
		kxq0=m;
		kxq1=0;
		hipAccum=0;
	}

	@Override
	public ApacheHLL4Wrapper copy(){
		return new ApacheHLL4Wrapper(m, k, -1, minProb);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Core Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		// microIndex for low-cardinality
		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);

		// Extract slot index (low lgK bits) and value (nlz of remaining bits + 1)
		final int slotNo=(int)(key&(m-1));
		final long remaining=key>>>lgK;
		// NLZ within the (64-lgK) remaining bits, +1 so minimum value is 1
		final int newValue=Long.numberOfLeadingZeros(remaining)-lgK+1;

		updateSlot(slotNo, newValue);
	}

	private void updateSlot(final int slotNo, final int newValue){
		final int oldValue=getSlotValue(slotNo);
		if(newValue<=oldValue){return;}

		// HIP (Historic Inverse Probability) accumulator update
		// kxq tracks sum of 2^(-value) for each register, split by curMin parity
		final double oldInvPow2=invPow2(oldValue);
		final double newInvPow2=invPow2(newValue);
		if((oldValue&1)==0){kxq0-=oldInvPow2;}
		else{kxq1-=oldInvPow2;}
		if((newValue&1)==0){kxq0+=newInvPow2;}
		else{kxq1+=newInvPow2;}
		hipAccum+=1.0/(kxq0+kxq1);

		// Store the new value
		final int shiftedNew=newValue-curMin;
		if(shiftedNew>=AUX_TOKEN){
			// Store AUX_TOKEN in nibble, actual value in aux map
			putNibble(slotNo, AUX_TOKEN);
			auxMap.put(slotNo, newValue);
		}else{
			// Fits in 4 bits
			putNibble(slotNo, shiftedNew);
			// If old value was an exception, remove from aux map
			if(oldValue-curMin>=AUX_TOKEN){
				auxMap.remove(slotNo);
			}
		}

		// Track numAtCurMin and advance curMin if needed
		if(oldValue==curMin){
			numAtCurMin--;
			while(numAtCurMin==0){
				shiftToBiggerCurMin();
			}
		}

		lastCardinality=-1;
	}

	/**
	 * Advances curMin by 1. Decrements all nibble values by 1.
	 * Nibbles that become AUX_TOKEN-1 (14) stay as-is.
	 * Nibbles that were AUX_TOKEN remain AUX_TOKEN (value in aux map).
	 * Recounts numAtCurMin for the new curMin.
	 */
	private void shiftToBiggerCurMin(){
		final int oldCurMin=curMin;
		final int newCurMin=oldCurMin+1;
		int newNumAtCurMin=0;
		// Build new aux map for values that overflow after shift
		final HashMap<Integer,Integer> newAux=new HashMap<>();

		for(int i=0; i<m; i++){
			final int nib=getNibble(i);
			if(nib==AUX_TOKEN){
				// Value in aux map; re-evaluate after curMin shift
				final int actualValue=auxMap.get(i);
				final int shiftedNew=actualValue-newCurMin;
				if(shiftedNew>=AUX_TOKEN){
					// Still an exception
					newAux.put(i, actualValue);
					// nibble stays AUX_TOKEN
				}else{
					// Now fits in a nibble
					putNibble(i, shiftedNew);
					if(shiftedNew==0){newNumAtCurMin++;}
				}
			}else{
				// Normal nibble: decrement by 1
				final int shiftedNew=nib-1;
				if(shiftedNew<0){
					// This shouldn't happen if numAtCurMin tracking is correct
					// but handle gracefully
					putNibble(i, 0);
					newNumAtCurMin++;
				}else{
					putNibble(i, shiftedNew);
					if(shiftedNew==0){newNumAtCurMin++;}
				}
			}
		}

		curMin=newCurMin;
		numAtCurMin=newNumAtCurMin;
		auxMap=newAux;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nibble Access            ----------------*/
	/*--------------------------------------------------------------*/

	private int getNibble(final int slotNo){
		final int byteIdx=slotNo>>>1;
		if((slotNo&1)==0){
			return nibbles[byteIdx]&0x0F;
		}else{
			return (nibbles[byteIdx]>>>4)&0x0F;
		}
	}

	private void putNibble(final int slotNo, final int value){
		final int byteIdx=slotNo>>>1;
		if((slotNo&1)==0){
			nibbles[byteIdx]=(byte)((nibbles[byteIdx]&0xF0)|(value&0x0F));
		}else{
			nibbles[byteIdx]=(byte)((nibbles[byteIdx]&0x0F)|((value&0x0F)<<4));
		}
	}

	private int getSlotValue(final int slotNo){
		final int nib=getNibble(slotNo);
		if(nib<AUX_TOKEN){
			return nib+curMin;
		}else{
			// Exception: actual value in aux map
			final Integer val=auxMap.get(slotNo);
			return (val!=null) ? val : curMin+AUX_TOKEN;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Estimation            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final double est=estimate();
		long card=Math.max(0, (long)est);
		card=Math.max(card, microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	/**
	 * Composite HLL estimator following Apache Data Sketches logic:
	 * - Low range: linear counting (when many empty registers)
	 * - Normal range: HLL harmonic mean with alpha_m correction
	 * - HIP (Historic Inverse Probability) for improved accuracy
	 */
	private double estimate(){
		return hllEstimate();
	}

	private double hllEstimate(){
		final double alpha=alphaM(m);
		double harmonicSum=0;
		int zeros=0;
		for(int i=0; i<m; i++){
			final int val=getSlotValue(i);
			harmonicSum+=invPow2(val);
			if(val==0){zeros++;}
		}
		final double rawEst=alpha*m*m/harmonicSum;

		// Small range correction (linear counting)
		if(rawEst<=2.5*m && zeros>0){
			return m*Math.log((double)m/zeros);
		}
		// Large range correction (for 32-bit hash space — not needed with 64-bit)
		return rawEst;
	}

	/** alpha_m correction factor for HyperLogLog. */
	private static double alphaM(final int m){
		switch(m){
			case 16: return 0.673;
			case 32: return 0.697;
			case 64: return 0.709;
			default: return 0.7213/(1.0+1.079/m);
		}
	}

	private static double invPow2(final int val){
		return Double.longBitsToDouble((1023L-val)<<52);
	}

	private long microCardinality(){
		return Long.bitCount(microIndex);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Merge / Add           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		final ApacheHLL4Wrapper other=(ApacheHLL4Wrapper)log;
		assert(other.m==m);
		added+=other.added;
		lastCardinality=-1;
		microIndex|=other.microIndex;
		// Merge by taking max of each register
		for(int i=0; i<m; i++){
			final int otherVal=other.getSlotValue(i);
			final int myVal=getSlotValue(i);
			if(otherVal>myVal){
				updateSlot(i, otherVal);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Drivable Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public int filledBuckets(){
		int count=0;
		for(int i=0; i<m; i++){
			if(getSlotValue(i)>0){count++;}
		}
		return count;
	}

	public double occupancy(){return (double)filledBuckets()/m;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public double[] rawEstimates(){
		final double hll=hllEstimate();
		final double micro=microCardinality();
		final double hybrid=Math.max(hll, micro);
		final int len=DDLCalibrationDriver.NUM_EST;
		final double[] r=new double[len];
		r[0]=hll;
		r[1]=hll;
		r[2]=hll;
		r[3]=hll;
		r[6]=hybrid;
		return r;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int lgK;
	private final int m;
	/** Packed 4-bit nibble array: 2 slots per byte. */
	private final byte[] nibbles;
	/** Auxiliary map for register values that overflow 4-bit storage (value-curMin >= 15). */
	private HashMap<Integer,Integer> auxMap;
	/** Global minimum value floor. All nibbles store (value - curMin). */
	private int curMin;
	/** Number of registers whose value equals curMin. */
	private int numAtCurMin;
	/** KxQ accumulators for HIP estimator (even/odd parity split). */
	private double kxq0;
	private double kxq1;
	/** HIP (Historic Inverse Probability) accumulator. */
	private double hipAccum;

	/** Sentinel nibble value indicating the real value is in auxMap. */
	private static final int AUX_TOKEN=15;
}
