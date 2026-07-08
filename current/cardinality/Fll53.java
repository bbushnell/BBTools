package cardinality;

import shared.Tools;

/**
 * Fll53: Brian's 21-bit organism — 4-bit localExp + 5 x 3-bit antennae + 2-bit
 * future-rehit counter, packed 3 per 64-bit long (1 spare bit).  The 3-bit
 * antennae cover tiers {floor, +1, +2}: wider aperture = ~4x less overflow and
 * much less recovered-dropout duplicate distortion, at ~25% fewer organisms
 * per byte.  512 bytes = 64 longs = 192 organisms = 960 tails.
 * Promotion when all 5 floor bits set: tails shift down one tier, counter
 * clears, localExp++.  IOT: delta outside [0,2] is dropped.
 * Counter increments on re-hits at FUTURE levels (delta 1 or 2), saturating.
 *
 * @author Amber (design: Brian)
 * @date July 2026
 */
public final class Fll53 extends CardinalityTracker {

	static final int BUCKETS_PER_WORD=5;
	static final int ORG_BITS=21;
	static final int ORGS_PER_LONG=3;
	/** Within an organism: [20:17] exp, [16:2] antennae, [1:0] counter. */
	static final int EXP_SHIFT=17;
	static final int FIELD_SHIFT=2;
	static final int FIELD_MASK=0x7FFF;      // 15 antenna bits
	/** Floor bits within the 15-bit field: bits 0,3,6,9,12. */
	static final int LSB_MASK=0x1249;
	/** Promotion keep-mask: per-tail bits {1,2} -> after >>>1 land on {0,1}. */
	static final int SHIFT_KEEP=0x36DB;
	static final int CTR_MASK=0x3;
	static final int ORG_MASK=(1<<ORG_BITS)-1;

	/** Counter-thinning policy for the 2-bit future-rehit counter (sat 3),
	 * ported from Fll83.CTR_MODE with thresholds scaled to the 2-bit range:
	 *   0 = off (full bump, old behavior; default)
	 *   1 = bit-13 thinning: ctr += (key>>>13)&1 (halves effective count,
	 *       doubles range to ~6; the 83 real-stream winner)
	 *   2 = piecewise: full bump below 2, bit-13 at 2+
	 * Other values fall through to full bump.  Bit 13 is free for the same
	 * reason as in Fll83.  Tables/nets are calibrated per counter semantic.
	 * Set per-JVM via -Dfll53.ctrMode=N; announced on stderr when nonzero. */
	public static int CTR_MODE=Integer.getInteger("fll53.ctrMode", 0);
	static{
		if(CTR_MODE!=0){System.err.println("Fll53.CTR_MODE="+CTR_MODE
			+(CTR_MODE==1 ? " (bit-13 thinning)"
			: CTR_MODE==2 ? " (piecewise: full<2, bit-13 at 2+)"
			: " (unknown -> full bump)"));}
	}

	Fll53(int buckets_, int k_, long seed, float minProb_){
		super(nextPow2(roundUpToOrgs(buckets_)*BUCKETS_PER_WORD), k_, seed, minProb_);
		numOrgs=roundUpToOrgs(buckets_);
		modBuckets=numOrgs*BUCKETS_PER_WORD;
		packed=new long[(numOrgs+ORGS_PER_LONG-1)/ORGS_PER_LONG];
		numLocalZeros=numOrgs;
	}

	private static int roundUpToOrgs(int b){return (b+BUCKETS_PER_WORD-1)/BUCKETS_PER_WORD;}
	private static int nextPow2(int n){return Integer.highestOneBit(Math.max(1, n-1))<<1;}

	@Override
	public Fll53 copy(){return new Fll53(modBuckets, k, -1, minProb);}

	public int getOrg(int i){
		return (int)((packed[i/ORGS_PER_LONG]>>>((i%ORGS_PER_LONG)*ORG_BITS))&ORG_MASK);
	}

	private void setOrg(int i, int v){
		final int idx=i/ORGS_PER_LONG, sh=(i%ORGS_PER_LONG)*ORG_BITS;
		packed[idx]=(packed[idx]&~((long)ORG_MASK<<sh))|((long)(v&ORG_MASK)<<sh);
	}

	@Override
	public final void hashAndStore(final long number){
		final long key=Tools.hash64shift(number^hashXor);
		if(Long.compareUnsigned(key, eeMask)>0){return;}

		final int hashNLZ=Long.numberOfLeadingZeros(key);
		final int bucket=(int)Long.remainderUnsigned(key, modBuckets);
		final int orgIdx=bucket/BUCKETS_PER_WORD;
		final int register=bucket%BUCKETS_PER_WORD;

		int org=getOrg(orgIdx);
		final int localExp=(org>>>EXP_SHIFT)&0xF;
		final int delta=hashNLZ-(globalExp+localExp);
		if(delta<0 || delta>2){return;}   // IOT, 3-tier aperture

		final int bitToSet=1<<(FIELD_SHIFT+register*3+delta);
		if((org&bitToSet)!=0){
			// Future-level re-hit bumps the saturating 2-bit counter,
			// subject to the CTR_MODE thinning policy (see field javadoc)
			if(delta>=1 && (org&CTR_MASK)<3){
				final int c=org&CTR_MASK;
				final boolean bump;
				switch(CTR_MODE){
					case 1: bump=((key>>>13)&1)!=0; break;
					case 2: bump=c<2 || ((key>>>13)&1)!=0; break;
					default: bump=true;
				}
				if(bump){setOrg(orgIdx, org+1);}
			}
			return;
		}

		lastCardinality=-1;
		org|=bitToSet;
		final int field=(org>>>FIELD_SHIFT)&FIELD_MASK;
		if((field&LSB_MASK)==LSB_MASK){
			final boolean wasZeroExp=(localExp==0);
			org=promote(org);
			setOrg(orgIdx, org);
			if(wasZeroExp){
				numLocalZeros--;
				if(numLocalZeros==0){advanceGlobal();}
			}
		}else{
			setOrg(orgIdx, org);
		}
	}

	/** Tails shift down one tier; counter clears; localExp++.  Cascades. */
	static int promote(int org){
		while(true){
			final int field=(org>>>FIELD_SHIFT)&FIELD_MASK;
			if((field&LSB_MASK)!=LSB_MASK){break;}
			final int localExp=(org>>>EXP_SHIFT)&0xF;
			if(localExp>=15){break;}
			final int newField=(field>>>1)&SHIFT_KEEP;
			org=((localExp+1)<<EXP_SHIFT)|(newField<<FIELD_SHIFT);
		}
		return org;
	}

	private void advanceGlobal(){
		while(numLocalZeros==0){
			globalExp++;
			numLocalZeros=0;
			for(int i=0; i<numOrgs; i++){
				final int org=getOrg(i);
				final int le=((org>>>EXP_SHIFT)&0xF)-1;
				assert le>=0 : "localExp=0 during advanceGlobal, org "+i;
				setOrg(i, (le<<EXP_SHIFT)|(org&((1<<EXP_SHIFT)-1)));
				if(le==0){numLocalZeros++;}
			}
			eeMask=(-1L)>>>globalExp;
		}
	}

	@Override
	public long cardinality(){
		double meanExp=0;
		for(int i=0; i<numOrgs; i++){meanExp+=(getOrg(i)>>>EXP_SHIFT)&0xF;}
		meanExp=globalExp+meanExp/numOrgs;
		return (long)(modBuckets*Math.pow(2.0, meanExp));
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public void add(CardinalityTracker log){throw new UnsupportedOperationException();}

	public int getNumOrgs(){return numOrgs;}
	public int getGlobalExp(){return globalExp;}

	private final long[] packed;
	private final int numOrgs;
	private final int modBuckets;
	private int globalExp=0;
	private int numLocalZeros;
	private long eeMask=-1L;
	private long lastCardinality=-1;
}
