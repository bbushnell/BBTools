package cardinality;

import parse.Parser;
import shared.Tools;

/**
 * LogLog30: 6-bit LogLog with 5 registers packed per 32-bit word.
 * Uses modulo bucket selection; bucket count is always a multiple of 5.
 * 30 bits used per word (5 × 6), 2 bits wasted per word.
 * At 2KB: 512 ints × 5 = 2560 registers (vs LL6's 2048 at 2KB).
 *
 * @author Nahida
 * @date May 2026
 */
public final class LogLog30 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	LogLog30(){this(2048, 31, -1, 0);}

	LogLog30(Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	LogLog30(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		modBuckets=roundTo5(buckets_);
		words=new int[(modBuckets+4)/5];
	}

	@Override
	public LogLog30 copy(){return new LogLog30(modBuckets, k, -1, minProb);}

	@Override public int actualBuckets(){return modBuckets;}
	@Override public int bitsPerWord(){return 32;}
	@Override public int bucketsPerWord(){return 5;}

	private static int roundTo5(int x){
		return Math.max(5, ((x+4)/5)*5);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private int getReg(int idx){
		int word=words[idx/5];
		int shift=(idx%5)*6;
		return (word>>>shift)&0x3F;
	}

	private void setReg(int idx, int val){
		int wi=idx/5;
		int shift=(idx%5)*6;
		words[wi]=(words[wi]&~(0x3F<<shift))|(val<<shift);
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);

		final int nlz=Long.numberOfLeadingZeros(key);
		final int bucket=(int)(Long.remainderUnsigned(key, modBuckets));

		final long micro=(key>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);
		if(LAZY_ALLOCATE){
			if(Long.bitCount(microIndex)<MICRO_CUTOFF_BITS){return;}
		}

		final int newStored=Math.min(nlz+1, 63);
		final int oldStored=getReg(bucket);

		if(newStored<=oldStored){return;}
		lastCardinality=-1;
		if(oldStored==0){filledCount++;}
		setReg(bucket, newStored);
	}

	private CardStats summarize(){
		nlzCounts=new int[66];
		int filled=0;
		for(int i=0; i<modBuckets; i++){
			final int stored=getReg(i);
			if(stored>0){
				final int absNlz=stored-1;
				if(absNlz<64){nlzCounts[absNlz+1]++;}
				filled++;
			}
		}
		nlzCounts[0]=modBuckets-filled;
		lastRawNlz=nlzCounts;
		lastCorrNlz=nlzCounts;
		return new CardStats(null, nlzCounts, 0, 0, 0, 0,
				modBuckets, 0, added, CF_MATRIX, CF_BUCKETS, 0,
				terminalMeanCF(), terminalMeanPlusCF());
	}

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final CardStats s=summarize();
		final double rawHyb=s.hybridDLL();
		long card=(long)(rawHyb);
		card=Math.max(card, s.microCardinality());
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((LogLog30)log);
	}

	public void add(LogLog30 log){
		added+=log.added;
		lastCardinality=-1;
		microIndex|=log.microIndex;
		assert(modBuckets==log.modBuckets);
		for(int i=0; i<modBuckets; i++){
			final int a=getReg(i);
			final int b=log.getReg(i);
			if(b>a){setReg(i, b);}
		}
	}

	@Override public int filledBuckets(){return filledCount;}
	@Override public double occupancy(){return (double)filledCount/modBuckets;}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	public int[] bucketValues(){
		int[] out=new int[modBuckets];
		for(int i=0; i<modBuckets; i++){out[i]=getReg(i);}
		return out;
	}

	@Override
	public double[] rawEstimates(){
		final CardStats s=summarize();
		final double hybridEst=s.hybridDLL();
		return AbstractCardStats.buildLegacyArray(s, hybridEst);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] words;
	private final int modBuckets;
	private int filledCount=0;
	int[] lastRawNlz, lastCorrNlz;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	public static final String CF_FILE="?cardinalityCorrectionLL6.tsv.gz";
	private static int CF_BUCKETS=2048;
	private static float[][] CF_MATRIX=initializeCF(CF_BUCKETS);
	public static float[][] initializeCF(int buckets){
		CF_BUCKETS=buckets;
		return CF_MATRIX=CorrectionFactor.loadFile(CF_FILE, buckets);
	}

	@Override public float terminalMeanCF(){return 0.720969f;}

}
