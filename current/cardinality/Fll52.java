package cardinality;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Tools;

/**
 * Fll52: the FLL4 candidate species — 16-bit organisms with 5 antennae and an
 * in-word duplicate sense, NOTHING auxiliary (honest memory accounting).
 * Word layout: [15:12] localExp, [11:2] 5 x 2-bit antennae (LSB=floor hit,
 * MSB=future hit), [1:0] 2-bit saturating counter of FUTURE-BIT RE-HITS
 * (the shape-separation winner), reset on promotion (recent-window semantics).
 * Promotion fires when all 5 floor bits set; MSBs shift to LSB positions,
 * counter clears, localExp++.  Global advance as in FLL2.
 * 512 bytes = 256 words = 1280 tails, counter included (Brian's number).
 *
 * @author Amber (design: Brian)
 * @date July 2026
 */
public final class Fll52 extends CardinalityTracker {

	static final int BUCKETS_PER_WORD=5;
	/** Bits 2,4,6,8,10 — floor bits. */
	static final int LSB_MASK=0x554;
	/** Bits 3,5,7,9,11 — future bits (antennae tips). */
	static final int MSB_MASK=0xAA8;
	static final int CTR_MASK=0x3;

	Fll52(){this(2560, 31, -1, 0);}   // 512 words = 1KB default

	Fll52(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	Fll52(int buckets_, int k_, long seed, float minProb_){
		super(nextPow2(roundUpToWords(buckets_)*BUCKETS_PER_WORD), k_, seed, minProb_);
		numWords=roundUpToWords(buckets_);
		modBuckets=numWords*BUCKETS_PER_WORD;
		words=new short[numWords];
		numLocalZeros=numWords;
	}

	private static int roundUpToWords(int b){return (b+BUCKETS_PER_WORD-1)/BUCKETS_PER_WORD;}
	private static int nextPow2(int n){return Integer.highestOneBit(Math.max(1, n-1))<<1;}

	@Override
	public Fll52 copy(){return new Fll52(modBuckets, k, -1, minProb);}

	@Override
	public final void hashAndStore(final long number){
		final long key=Tools.hash64shift(number^hashXor);
		if(Long.compareUnsigned(key, eeMask)>0){return;}

		final int hashNLZ=Long.numberOfLeadingZeros(key);
		final int bucket=(int)Long.remainderUnsigned(key, modBuckets);
		final int wordIdx=bucket/BUCKETS_PER_WORD;
		final int register=bucket%BUCKETS_PER_WORD;

		int word=words[wordIdx]&0xFFFF;
		final int localExp=(word>>>12)&0xF;
		final int delta=hashNLZ-(globalExp+localExp);
		if(delta<0 || delta>1){return;}   // IOT

		final int bitToSet=1<<(2+delta+register*2);
		if((word&bitToSet)!=0){
			// Future-bit re-hit: bump the in-word saturating counter
			if(delta==1 && (word&CTR_MASK)<3){words[wordIdx]=(short)(word+1);}
			return;
		}

		lastCardinality=-1;
		word|=bitToSet;
		if((word&LSB_MASK)==LSB_MASK){
			final boolean wasZeroExp=(localExp==0);
			word=promote(word);
			words[wordIdx]=(short)word;
			if(wasZeroExp){
				numLocalZeros--;
				if(numLocalZeros==0){advanceGlobal();}
			}
		}else{
			words[wordIdx]=(short)word;
		}
	}

	/** Promotion: future bits become floor bits, counter clears, localExp++. */
	static int promote(int word){
		while((word&LSB_MASK)==LSB_MASK){
			final int localExp=(word>>>12)&0xF;
			if(localExp>=15){break;}
			word=((word&MSB_MASK)>>>1)|((localExp+1)<<12);
		}
		return word;
	}

	private void advanceGlobal(){
		while(numLocalZeros==0){
			globalExp++;
			numLocalZeros=0;
			for(int i=0; i<numWords; i++){
				final int w=words[i]&0xFFFF;
				final int le=((w>>>12)&0xF)-1;
				assert le>=0 : "localExp=0 during advanceGlobal, word "+i;
				words[i]=(short)((le<<12)|(w&0xFFF));
				if(le==0){numLocalZeros++;}
			}
			eeMask=(-1L)>>>globalExp;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------   Sufficient statistic (56)   ----------------*/
	/*--------------------------------------------------------------*/

	static final int NCLASS=56;
	static final int NC=8*NCLASS;
	/** Prefix sums over a=count(00 buckets) for 5 antennae. */
	static final int[] OFFSET_A={0, 21, 36, 46, 52, 55};

	/** 10-bit antenna field (word bits 2..11 shifted down) -> 56-class index. */
	static int idx56(int hist){
		final int h1=hist>>>1;
		final int LSB5=0x155;
		final int d=Integer.bitCount(hist&h1&LSB5);
		final int a=5-Integer.bitCount((hist|h1)&LSB5);
		final int c=Integer.bitCount(~hist&h1&LSB5);
		final int b=5-a-c-d;
		final int r=5-a;
		return OFFSET_A[a]+b*(r+1)-b*(b-1)/2+c;
	}

	static int classOf(int word){
		final int le=Math.min((word>>>12)&0xF, 7);
		return le*NCLASS+idx56((word>>>2)&0x3FF);
	}

	/*--------------------------------------------------------------*/
	/*----------------     Complexity features       ----------------*/
	/*--------------------------------------------------------------*/

	/** Fixed window for the modular complexity net: 256 words = 512 bytes. */
	public static final int WIN=256;

	/** 16 features over words [w0, w0+WIN), for the modular window net. */
	static double[] windowFeatures(Fll52 f, long adds, int w0){
		final int nw=f.numWords;
		final double winAdds=Math.max(1.0, adds*(double)WIN/nw);
		long lsb=0, msb=0;
		double expSum=0, expSat=0;
		int satCnt=0, zeroCnt=0;
		final int[] cc=new int[4];
		for(int w=w0; w<w0+WIN; w++){
			final int word=f.words[w]&0xFFFF;
			final int le=(word>>>12)&0xF;
			expSum+=le;
			lsb+=Integer.bitCount(word&LSB_MASK);
			msb+=Integer.bitCount(word&MSB_MASK);
			final int c=word&CTR_MASK;
			cc[c]++;
			if(c==3){expSat+=le; satCnt++;}
			if((word&0xFFC)==0){zeroCnt++;}
		}
		final double meanExp=f.globalExp+expSum/WIN;
		final double lf=(double)lsb/(WIN*5), mf=(double)msb/(WIN*5);
		final double f1=(double)cc[1]/WIN, f2=(double)cc[2]/WIN, f3=(double)cc[3]/WIN;
		final double r=Math.log(winAdds)*INV_LN2-meanExp;
		final double dExpSat=(satCnt>0) ? expSat/satCnt-expSum/WIN : 0;
		return new double[]{r, lf, mf, dExpSat, f1, f2, f3,
			(f1+2*f2+3*f3)/3, f2+f3, lf*mf, r*f3, mf*f3, lf*lf, mf*mf, f3*f3,
			(double)zeroCnt/WIN};
	}

	private static final double INV_LN2=1.0/Math.log(2.0);

	/*--------------------------------------------------------------*/
	/*----------------      Window net (16-H-1)      ----------------*/
	/*--------------------------------------------------------------*/

	/** Tiny tanh MLP with linear output; text format from the training farm. */
	static final class WindowNet {
		final int H, D=16;
		final double[][] w1;
		final double[] b1, w2, mean, sd;
		final double b2;

		private WindowNet(java.util.ArrayList<String> lines){
			H=Integer.parseInt(lines.get(0).trim());
			mean=new double[D]; sd=new double[D]; b1=new double[H]; w2=new double[H];
			w1=new double[H][D];
			final double[][] heads={mean, sd, b1, w2};
			for(int i=0; i<4; i++){
				final String[] f=lines.get(1+i).trim().split("\t");
				for(int j=0; j<heads[i].length; j++){heads[i][j]=Double.parseDouble(f[j]);}
			}
			for(int h=0; h<H; h++){
				final String[] f=lines.get(5+h).trim().split("\t");
				for(int j=0; j<D; j++){w1[h][j]=Double.parseDouble(f[j]);}
			}
			b2=Double.parseDouble(lines.get(5+H).trim());
		}

		double predict(double[] raw){
			double out=b2;
			for(int h=0; h<H; h++){
				double z=b1[h];
				for(int d=0; d<D; d++){z+=w1[h][d]*(raw[d]-mean[d])/sd[d];}
				out+=w2[h]*Math.tanh(z);
			}
			return out;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Resource loading         ----------------*/
	/*--------------------------------------------------------------*/

	static final String NET_FILE="fll52net.tsv.gz";
	private static WindowNet net=null;
	private static boolean netLoadTried=false;

	static synchronized WindowNet net(){
		if(net!=null || netLoadTried){return net;}
		netLoadTried=true;
		final String path=Data.findPath("?"+NET_FILE);
		if(path==null){
			System.err.println("WARNING: Fll52 window net not found: "+NET_FILE
				+" (complexity blending disabled; duplicate-heavy streams will overcount)");
			return null;
		}
		try{
			final FileFormat ff=FileFormat.testInput(path, FileFormat.TEXT, null, false, true);
			final ByteFile bf=ByteFile.makeByteFile(ff, 1);
			final java.util.ArrayList<String> lines=new java.util.ArrayList<>();
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lines.add(new String(line));
			}
			bf.close();
			net=new WindowNet(lines);
		}catch(Exception e){
			System.err.println("WARNING: Fll52 window net failed to load: "+e);
		}
		return net;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Estimation             ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		final double mleEst;
		final Fll52MLE mle=Fll52MLE.forWords(numWords);
		mleEst=(mle!=null) ? mle.estimate(this) : crudeEstimate();
		double est=mleEst;
		// Complexity blend: modular window net, k predictions averaged in log2.
		final WindowNet wn=net();
		if(wn!=null && numWords>=WIN && mleEst>=32.0*modBuckets){
			double ySum=0;
			int k=0;
			for(int w0=0; w0+WIN<=numWords; w0+=WIN){
				ySum+=wn.predict(windowFeatures(this, added, w0));
				k++;
			}
			final double cHat=Math.pow(2.0, ySum/k);
			final double w=Math.max(0, Math.min(1, (cHat-0.85)/0.10));
			est=w*mleEst+(1-w)*added*cHat;
		}
		long card=Math.max(0, (long)est);
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	/** Moment fallback when no likelihood table matches (rough; unbiased-ish). */
	private long crudeEstimate(){
		double meanExp=0;
		for(int i=0; i<numWords; i++){meanExp+=(words[i]>>>12)&0xF;}
		meanExp=globalExp+meanExp/numWords;
		return (long)(modBuckets*Math.pow(2.0, meanExp));
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public void add(CardinalityTracker log){throw new UnsupportedOperationException();}

	public int getWord(int i){return words[i]&0xFFFF;}
	public int getNumWords(){return numWords;}
	public int getGlobalExp(){return globalExp;}

	private final short[] words;
	private final int numWords;
	private final int modBuckets;
	private int globalExp=0;
	private int numLocalZeros;
	private long eeMask=-1L;
	private long lastCardinality=-1;
}
