package cardinality;

import structures.IntHashMap;
import shared.Tools;

/**
 * HyperLogLogLog (HLLL) cardinality estimator, ported for BBTools benchmarking.
 * <p>
 * Original paper: "HyperLogLogLog: Cardinality Estimation With One Log More"
 * Authors: Matti Karppa, Rasmus Pagh (KDD 2022 / ESA 2022)
 * Paper: https://dl.acm.org/doi/abs/10.1145/3534678.3539246
 * Source: https://github.com/mkarppa/hyperlogloglog
 * License: MIT
 * <p>
 * Ported from the C++ reference implementation (HyperLogLogLog.hpp) to Java.
 * Exception map uses BBTools IntHashMap for performance.
 * <p>
 * Core idea: m registers of mBits width (default 3) store offsets relative to
 * a global base B. Values that exceed the offset range are stored in a sparse
 * exception map S. Periodic rebasing minimizes |S|.
 * <p>
 * Estimation uses the standard HLL harmonic mean with small-range correction.
 *
 * @author Chloe (port), Matti Karppa (original algorithm and C++ implementation)
 * @date May 2026
 */
public final class HLLLWrapper extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	HLLLWrapper(){this(2048, 31, -1, 0);}

	HLLLWrapper(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	HLLLWrapper(int buckets_, int k_, long seed, float minProb_){
		this(buckets_, k_, seed, minProb_, DEFAULT_MBITS);
	}

	HLLLWrapper(int buckets_, int k_, long seed, float minProb_, int mBits_){
		super(buckets_, k_, seed, minProb_);
		this.mBits=mBits_;
		this.maxOffset=(1<<mBits)-1;
		this.logM=Integer.numberOfTrailingZeros(buckets);
		this.M=new byte[buckets];
		this.S=new IntHashMap(64);
		this.B=0;
		this.lowerBound=0;
	}

	@Override
	public HLLLWrapper copy(){return new HLLLWrapper(buckets, k, -1, minProb, mBits);}

	/*--------------------------------------------------------------*/
	/*----------------         Core Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		final long micro=(key>>>58)&0x3FL;
		microIndex|=(1L<<micro);

		final int j=(int)(key&bucketMask);
		final long w=key>>>logM;
		final int r=(w==0) ? (64-logM+1) : (Long.numberOfLeadingZeros(w)-logM+1);

		addJr(j, r);
	}

	private void addJr(int j, int r){
		if(r<=lowerBound){return;}

		final int sVal=S.get(j);
		final int r0=(sVal>=0) ? sVal : (M[j]&0xFF)+B;

		if(r0<r){
			if(B<=r && r<=B+maxOffset){
				if(sVal>=0){S.remove(j);}
				M[j]=(byte)(r-B);
			}else{
				final int oldSize=S.size();
				S.put(j, r);
				if(S.size()>oldSize){compressFull();}
			}
			lastCardinality=-1;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Compression            ----------------*/
	/*--------------------------------------------------------------*/

	private void compressFull(){
		final int m=buckets;
		final int[] regs=allRegisters();
		int bestNs=S.size();
		int bestBase=B;

		int minVal=Integer.MAX_VALUE;
		for(int i=0; i<m; i++){
			if(regs[i]<minVal){minVal=regs[i];}
		}
		lowerBound=minVal;

		int potentialBase=minVal;
		int nBelowB=0;
		while(nBelowB<bestNs && potentialBase<=64){
			int ns=0;
			int nextPotentialBase=Integer.MAX_VALUE;
			for(int i=0; i<m; i++){
				final int v=regs[i];
				if(v<potentialBase || v>potentialBase+maxOffset){
					ns++;
				}
				if(v==potentialBase){nBelowB++;}
				if(v>potentialBase && v<nextPotentialBase){
					nextPotentialBase=v;
				}
			}
			if(ns<bestNs){
				bestNs=ns;
				bestBase=potentialBase;
			}
			if(nextPotentialBase==Integer.MAX_VALUE){break;}
			potentialBase=nextPotentialBase;
		}

		if(bestBase!=B){
			rebase(regs, bestBase);
		}
	}

	private int[] allRegisters(){
		final int m=buckets;
		if(regsBuf==null){regsBuf=new int[m];}
		for(int i=0; i<m; i++){
			final int sVal=S.get(i);
			regsBuf[i]=(sVal>=0) ? sVal : (M[i]&0xFF)+B;
		}
		return regsBuf;
	}

	private void rebase(int[] regs, int newB){
		final int m=buckets;
		S.clear();
		for(int i=0; i<m; i++){
			final int r=regs[i];
			if(newB<=r && r<=newB+maxOffset){
				M[i]=(byte)(r-newB);
			}else{
				S.put(i, r);
			}
		}
		B=newB;
	}

	private int getRegister(int i){
		final int sVal=S.get(i);
		return (sVal>=0) ? sVal : (M[i]&0xFF)+B;
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

	private double estimate(){
		final int m=buckets;
		double harmonicSum=0;
		int V=0;
		for(int i=0; i<m; i++){
			final int r=getRegister(i);
			V+=(r==0) ? 1 : 0;
			harmonicSum+=Math.pow(2.0, -r);
		}
		double E=alpha(m)*((double)m)*((double)m)/harmonicSum;

		if(E<=2.5*m && V!=0){
			return m*Math.log((double)m/V);
		}
		return E;
	}

	private static double alpha(int m){
		switch(m){
			case 16: return 0.673;
			case 32: return 0.697;
			case 64: return 0.709;
			default: return 0.7213/(1.0+1.079/m);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Merge / Add            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		final HLLLWrapper other=(HLLLWrapper)log;
		assert(other.mBits==mBits && other.buckets==buckets);
		added+=other.added;
		lastCardinality=-1;
		microIndex|=other.microIndex;

		final int m=buckets;
		for(int i=0; i<m; i++){
			final int r1=getRegister(i);
			final int r2=other.getRegister(i);
			final int r=Math.max(r1, r2);
			if(B<=r && r<=B+maxOffset){
				M[i]=(byte)(r-B);
				S.remove(i);
			}else{
				S.put(i, r);
			}
		}
		compressFull();
	}

	/*--------------------------------------------------------------*/
	/*----------------      Drivable Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public int filledBuckets(){
		final int m=buckets;
		int count=0;
		for(int i=0; i<m; i++){
			if(getRegister(i)!=0){count++;}
		}
		return count;
	}

	public double occupancy(){return (double)filledBuckets()/buckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}
	@Override public int bitsPerWord(){return mBits+1;}
	@Override public int bucketsPerWord(){return 1;}

	@Override
	public double[] rawEstimates(){
		final double est=estimate();
		final double micro=microCardinality();
		final double hybrid=Math.max(est, micro);
		final int len=DDLCalibrationDriver.NUM_EST;
		final double[] r=new double[len];
		r[0]=est;
		r[1]=est;
		r[2]=est;
		r[3]=est;
		r[6]=hybrid;
		return r;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Utility              ----------------*/
	/*--------------------------------------------------------------*/

	private long microCardinality(){
		return Long.bitCount(microIndex);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final byte[] M;
	private final IntHashMap S;
	private final int mBits;
	private final int maxOffset;
	private final int logM;
	private int B;
	private int lowerBound;
	private int[] regsBuf;

	static final int DEFAULT_MBITS=3;
}
