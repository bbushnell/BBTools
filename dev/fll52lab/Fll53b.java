package cardinality;

import shared.Tools;

/**
 * Fll53b: 20-bit organism — 4-bit localExp + 4 x 3-bit antennae + TWO 2-bit
 * counters at adjacent future levels (ctrA watches delta=1, ctrB delta=2),
 * packed 3 per long.  On promotion, Brian's FLOW mode shifts ctrB into ctrA
 * (the tier it watched just became ctrA's tier) instead of resetting both —
 * a two-stage shift register of duplicate pressure aging with the antennae.
 * Toggle Fll53b.FLOW to compare flow vs reset.  512B = 192 organisms = 768 tails.
 *
 * @author Amber (design: Brian)
 * @date July 2026
 */
public final class Fll53b extends CardinalityTracker {

	static final int BUCKETS_PER_WORD=4;
	static final int ORG_BITS=20;
	static final int ORGS_PER_LONG=3;
	/** [19:16] exp, [15:4] 4x3-bit antennae, [3:2] ctrB (level-2), [1:0] ctrA (level-1). */
	static final int EXP_SHIFT=16;
	static final int FIELD_SHIFT=4;
	static final int FIELD_MASK=0xFFF;       // 12 antenna bits
	/** Floor bits within the 12-bit field: bits 0,3,6,9. */
	static final int LSB_MASK=0x249;
	static final int SHIFT_KEEP=0x6DB;
	static final int CTRA_SHIFT=0, CTRB_SHIFT=2;
	static final int ORG_MASK=(1<<ORG_BITS)-1;
	/** Promotion counter behavior: true = ctrB flows into ctrA; false = both reset. */
	public static boolean FLOW=true;

	Fll53b(int buckets_, int k_, long seed, float minProb_){
		super(nextPow2(roundUpToOrgs(buckets_)*BUCKETS_PER_WORD), k_, seed, minProb_);
		numOrgs=roundUpToOrgs(buckets_);
		modBuckets=numOrgs*BUCKETS_PER_WORD;
		packed=new long[(numOrgs+ORGS_PER_LONG-1)/ORGS_PER_LONG];
		numLocalZeros=numOrgs;
	}

	private static int roundUpToOrgs(int b){return (b+BUCKETS_PER_WORD-1)/BUCKETS_PER_WORD;}
	private static int nextPow2(int n){return Integer.highestOneBit(Math.max(1, n-1))<<1;}

	@Override
	public Fll53b copy(){return new Fll53b(modBuckets, k, -1, minProb);}

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
			// Per-level re-hit counters: A watches future-1, B watches future-2
			if(delta==1){
				final int c=(org>>>CTRA_SHIFT)&3;
				if(c<3){setOrg(orgIdx, org+(1<<CTRA_SHIFT));}
			}else if(delta==2){
				final int c=(org>>>CTRB_SHIFT)&3;
				if(c<3){setOrg(orgIdx, org+(1<<CTRB_SHIFT));}
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

	/** Tails shift down one tier; counters flow or reset; localExp++.  Cascades. */
	static int promote(int org){
		while(true){
			final int field=(org>>>FIELD_SHIFT)&FIELD_MASK;
			if((field&LSB_MASK)!=LSB_MASK){break;}
			final int localExp=(org>>>EXP_SHIFT)&0xF;
			if(localExp>=15){break;}
			final int newField=(field>>>1)&SHIFT_KEEP;
			// Brian's flow-on-promotion: the tier ctrB watched becomes the tier
			// ctrA watches, so the count moves house with it.  Else amnesia.
			final int ctrs=FLOW ? ((org>>>CTRB_SHIFT)&3)<<CTRA_SHIFT : 0;
			org=((localExp+1)<<EXP_SHIFT)|(newField<<FIELD_SHIFT)|ctrs;
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
	public int ctrA(int i){return (getOrg(i)>>>CTRA_SHIFT)&3;}
	public int ctrB(int i){return (getOrg(i)>>>CTRB_SHIFT)&3;}
	public int getGlobalExp(){return globalExp;}

	private final long[] packed;
	private final int numOrgs;
	private final int modBuckets;
	private int globalExp=0;
	private int numLocalZeros;
	private long eeMask=-1L;
	private long lastCardinality=-1;
}
