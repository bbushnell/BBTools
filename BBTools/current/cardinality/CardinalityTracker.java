package cardinality;

import java.util.Arrays;
import java.util.Random;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import jgi.Dedupe;
import shared.Parser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.LongList;
import structures.SuperLongList;
import ukmer.Kmer;

/**
 * Abstract superclass for cardinality-tracking structures like LogLog.
 * @author Brian Bushnell
 * @date Feb 20, 2020
 *
 */
public abstract class CardinalityTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Factory method; creates a tracker using default settings.
	 * Subclass is determined by static Parser.loglogType field.
	 */
	public static CardinalityTracker makeTracker(){
		if(trackCounts || "BBLog".equalsIgnoreCase(Parser.loglogType)){
			return new BBLog();//Fastest, most accurate
		}else if("LogLog".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog();//Least accurate
		}else if("LogLog2".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog2();//Slowest, uses mantissa
		}else if("LogLog16".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog16();//Uses 10-bit mantissa
		}else if("LogLog8".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog8();//Lowest memory
		}
		assert(false) : "TODO: "+Parser.loglogType;
		throw new RuntimeException(Parser.loglogType);
	}
	
	/** 
	 * Factory method; creates a tracker using parsed settings.
	 * Subclass is determined by static Parser.loglogType field.
	 */
	public static CardinalityTracker makeTracker(Parser p){
		if(trackCounts || "BBLog".equalsIgnoreCase(Parser.loglogType)){
			return new BBLog(p);
		}else if("LogLog".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog(p);
		}else if("LogLog2".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog2(p);
		}else if("LogLog16".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog16(p);
		}else if("LogLog8".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog8(p);
		}
		assert(false) : "TODO: "+Parser.loglogType;
		throw new RuntimeException(Parser.loglogType);
	}
	
	/** 
	 * Factory method; creates a tracker using specified settings.
	 * Subclass is determined by static Parser.loglogType field.
	 */
	public static CardinalityTracker makeTracker(int buckets_, int k_, long seed, float minProb_){
		if(trackCounts || "BBLog".equalsIgnoreCase(Parser.loglogType)){
			return new BBLog(buckets_, k_, seed, minProb_);
		}else if("LogLog".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog(buckets_, k_, seed, minProb_);
		}else if("LogLog2".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog2(buckets_, k_, seed, minProb_);
		}else if("LogLog16".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog16(buckets_, k_, seed, minProb_);
		}else if("LogLog8".equalsIgnoreCase(Parser.loglogType)){
			return new LogLog8(buckets_, k_, seed, minProb_);
		}
		assert(false) : "TODO: "+Parser.loglogType;
		throw new RuntimeException(Parser.loglogType);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create a tracker with parsed parameters. */
	public CardinalityTracker(Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}
	
	/**
	 * Create a tracker with specified parameters.
	 * @param buckets_ Number of buckets (counters)
	 * @param k_ Kmer length
	 * @param seed Random number generator seed; -1 for a random seed
	 * @param minProb_ Ignore kmers with under this probability of being correct
	 */
	public CardinalityTracker(int buckets_, int k_, long seed, float minProb_){
//		if((buckets_&1)==0){buckets_=(int)Primes.primeAtLeast(buckets_);} //Legacy code, needed modulo operation
		buckets=powerOf2AtLeast(buckets_);
		assert(buckets>0 && Integer.bitCount(buckets)==1) : "Buckets must be a power of 2: "+buckets;
		bucketMask=buckets-1;
		k=Kmer.getKbig(k_);
		minProb=minProb_;
		
		//For old hash function
//		tables=new long[numTables][][];
//		for(int i=0; i<numTables; i++){
//			tables[i]=makeCodes(steps, bits, (seed<0 ? -1 : seed+i));
//		}
		
		Random randy=Shared.threadLocalRandom(seed<0 ? -1 : seed);
		hashXor=randy.nextLong();
	}
	
	/** 
	 * Return the lowest power of 2 that is >= target. 
	 * Provided because buckets must currently be a power of 2. 
	 */
	public static final int powerOf2AtLeast(int target){
		if(target<1){return 1;}
		int ret=1, limit=Tools.min(target, 0x40000000);
		while(ret<limit){ret<<=1;}
		return ret;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Deprecated; use LogLogWrapper */
	public static final void main(String[] args){
		LogLogWrapper llw=new LogLogWrapper(args);
		
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		llw.process();
		
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
	}
	
	/** 
	 * Old table-based hash function.
	 * Slower than the new function.
	 */
	public final long hash(final long value0, final long[][] table){
		long value=value0, code=0;
		long mask=(bits>63 ? -1L : ~((-1L)<<bits));

		for(int i=0; i<steps; i++){//I could also do while value!=0
			int x=(int)(value&mask);
			value>>=bits;
			code=code^table[i][x];
		}
		return code;
	}
	
	/** Hash and add a number to this tracker */
	public final void add(long number){
		hashAndStore(number);
	}
	
	/** 
	 * Hash and track the Read 
	 * */
	public final void hash(Read r){
		if(r==null){return;}
		if(r.length()>=k){hash(r.bases, r.quality);}
		if(r.mateLength()>=k){hash(r.mate.bases, r.mate.quality);}
	}
	
	/** 
	 * Hash and track the sequence 
	 * */
	public final void hash(byte[] bases, byte[] quals){
		if(k<32){hashSmall(bases, quals);}
		else{hashBig(bases, quals);}
	}
	
	/** 
	 * Hash and track the sequence using short kmers 0<k<32 
	 * */
	public final void hashSmall(byte[] bases, byte[] quals){
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		int len=0;
		
		long kmer=0, rkmer=0;
		
		if(minProb>0 && quals!=null){//Debranched loop
			assert(quals.length==bases.length) : quals.length+", "+bases.length;
			float prob=1;
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				
				{//Update probability
					byte q=quals[i];
					prob=prob*PROB_CORRECT[q];
					if(len>k){
						byte oldq=quals[i-k];
						prob=prob*PROB_CORRECT_INVERSE[oldq];
					}
				}
				if(x>=0){
					len++;
				}else{
					len=0;
					kmer=rkmer=0;
					prob=1;
				}
				if(len>=k && prob>=minProb){
					add(Tools.max(kmer, rkmer));
				}
			}
		}else{

			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				
				if(x>=0){
					len++;
				}else{
					len=0;
					kmer=rkmer=0;
				}
				if(len>=k){
					add(Tools.max(kmer, rkmer));
				}
			}
		}
	}
	
	/** 
	 * Hash and track the sequence using long kmers >31 
	 * */
	public final void hashBig(byte[] bases, byte[] quals){
		
		Kmer kmer=getLocalKmer();
		int len=0;
		float prob=1;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=Dedupe.baseToNumber[b];
			kmer.addRightNumeric(x);
			if(minProb>0 && quals!=null){//Update probability
				prob=prob*PROB_CORRECT[quals[i]];
				if(len>k){
					byte oldq=quals[i-k];
					prob=prob*PROB_CORRECT_INVERSE[oldq];
				}
			}
			if(AminoAcid.isFullyDefined(b)){
				len++;
			}else{
				len=0;
				prob=1;
			}
			if(len>=k && prob>=minProb){
				add(kmer.xor());
			}
		}
	}
	
	/** 
	 * Make a table of random bitmasks for hashing. 
	 * Superceded by new hash method. 
	 * */
	private static final long[][] makeCodes(int length, int bits, long seed){
		if(true) {return null;}//Short circuit
		Random randy=Shared.threadLocalRandom(seed);
		int modes=1<<bits;
		long[][] r=new long[length][modes];
		for(int i=0; i<length; i++){
			for(int j=0; j<modes; j++){
				long x=randy.nextLong();
				while(Long.bitCount(x)>33){
					x&=(~(1L<<randy.nextInt(64)));
				}
				while(Long.bitCount(x)<31){
					x|=(1L<<randy.nextInt(64));
				}
				r[i][j]=x;
				
			}
		}
		return r;
	}
	
	public final float compensationFactorBuckets(){
		assert(Integer.bitCount(buckets)==1) : buckets;
		int zeros=Integer.numberOfTrailingZeros(buckets);
		return compensationFactorLogBuckets(zeros);
	}
	
	/** 
	 * Multiplier to compensate for overestimating genome size
	 * when the number of buckets is too small,
	 * as a function of log2(buckets)
	 * @param logBuckets
	 * @return Multiplier for final estimate
	 */
	public final float compensationFactorLogBuckets(int logBuckets){
		float[] array=compensationFactorLogBucketsArray();
		return (array!=null && logBuckets<array.length) ? array[logBuckets] : 1/(1+(1<<logBuckets));
	}
	
	public SuperLongList toFrequency(){
		SuperLongList list=new SuperLongList(1000);
		int[] counts=getCounts();
		for(int x : counts){
			if(x>0){list.add(x);}
		}
		list.sort();
		return list;
	}
	
	/** 
	 * Print a kmer frequency histogram.
	 * @param path File to which to print
	 * @param overwrite
	 * @param append
	 * @param supersample Adjust counts for the effect of subsampling
	 */
	public void printKhist(String path, boolean overwrite, boolean append, boolean supersample, int decimals){
		SuperLongList sll=toFrequency();
		ByteStreamWriter bsw=new ByteStreamWriter(path, overwrite, append, false);
		bsw.start();
		bsw.print("#Depth\tCount\n");
		final double mult=Tools.max(1.0, (supersample ? cardinality()/(double)buckets : 1));
		final long[] array=sll.array();
		final LongList list=sll.list();
		
		for(int depth=0; depth<array.length; depth++){
			long count=array[depth];
			if(count>0){
				bsw.print(depth).tab();
				if(supersample){
					if(decimals>0){
						bsw.print(count*mult, decimals).nl();
					}else{
						bsw.print(Tools.max(1, Math.round(count*mult))).nl();
					}
				}else{
					bsw.print(count).nl();
				}
			}
		}
		int count=0;
		long prevDepth=-1;
		for(int i=0; i<list.size; i++){
			long depth=list.get(i);
			if(depth!=prevDepth && count>0){
				assert(depth>prevDepth);
				bsw.print(prevDepth).tab();
				if(supersample){
					if(decimals>0){
						bsw.print(count*mult, decimals).nl();
					}else{
						bsw.print(Tools.max(1, Math.round(count*mult))).nl();
					}
				}else{
					bsw.print(count).nl();
				}
				count=0;
			}else{
				count++;
			}
			prevDepth=depth;
		}
		if(count>0){
			bsw.print(prevDepth).tab();
			if(supersample){
				if(decimals>0){
					bsw.print(count*mult, decimals).nl();
				}else{
					bsw.print(Tools.max(1, Math.round(count*mult))).nl();
				}
			}else{
				bsw.print(count).nl();
			}
		}
		bsw.poisonAndWait();
	}
	
	/** Sum of counts array, if present */
	public final long countSum(){
		int[] counts=getCounts();
		return counts==null ? 0 : shared.Vector.sum(counts);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Calculate cardinality estimate from this tracker */
	public abstract long cardinality();

	/** Counts array, if present.
	 * Should be overridden for classes that track counts. */
	public int[] getCounts(){
		return null;
	}
	
	/** Add another tracker to this one */
	public abstract void add(CardinalityTracker log);
	
	/** Generate a 64-bit hashcode from a number, and add it to this tracker */
	public abstract void hashAndStore(final long number);
	
	/** TODO: Deprecate?
	 * Designed to compensate for overestimate with small numbers of buckets.
	 * Appears to be handled by using the harmonic mean in the case of multiple trials. */
	public abstract float[] compensationFactorLogBucketsArray();
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Kmer length */
	public final int k;
	
	/** Ignore kmers under this probability of correctness */
	public final float minProb;
	
	/** 
	 * Number of buckets for tracking.
	 * Must be a power of 2.
	 * Larger is slower and uses more memory, but more accurate.
	 */
	public final int buckets;
	
	/** Mask for lower bits of hashcode to yield bucket; bucketMask==buckets-1. */
	final int bucketMask;
	
	/** Reusable kmer for hashing. */
	private final ThreadLocal<Kmer> localKmer=new ThreadLocal<Kmer>();
	
	/** The hash function will yield a different mapping for each mask here. */
	private final long hashXor;
	
	/** For long kmer mode */
	protected Kmer getLocalKmer(){
		Kmer kmer=localKmer.get();
		if(kmer==null){
			localKmer.set(new Kmer(k));
			kmer=localKmer.get();
		}
		kmer.clearFast();
		return kmer;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Deprecated Table Fields   ----------------*/
	/*--------------------------------------------------------------*/
	
	static final int numTables=4;
	static final int numTablesMask=numTables-1;
	/** Bits hashed per cycle; no longer needed */
	private static final int bits=8;
	private static final int steps=(63+bits)/bits;;
//	final long[][][] tables;
	
	/*--------------------------------------------------------------*/
	/*----------------            Statics           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Converts quality score to probability of correctness */
	public static final float[] PROB_CORRECT=Arrays.copyOf(align2.QualityTools.PROB_CORRECT, 128);
	/** Inverse of PROB_CORRECT */
	public static final float[] PROB_CORRECT_INVERSE=Arrays.copyOf(align2.QualityTools.PROB_CORRECT_INVERSE, 128);

	/** Atomic mode allows less memory when multithreaded, but lower peak concurrency due to contention. */
	public static final boolean atomic=false;//non-atomic is faster.
	/** Track number of times the value in each bucket was seen. */
	public static boolean trackCounts=false;
	/** Ignores kmers such that x%SKIPMOD!=0, to skip expensive hash and store functions */ 
	static final long SKIPMOD=1;//No longer used; requires a modulo operation
	/** Records the last cardinality tracked, for use in static contexts (RQCFilter) */
	public static long lastCardinality=-1;
//	/** Ignore hashed values above this, to skip expensive read and store functions. */
//	static final long maxHashedValue=((-1L)>>>3);//No longer used
	
	/** These determine how to combine multiple buckets to yield a value. 
	 * Mean is best even though mwa is pretty close.  Median is much worse. */
	public static boolean USE_MEAN=true;//Arithmetic mean
	public static boolean USE_MEDIAN=false;
	public static boolean USE_MWA=false;//Median-weighted-average
	public static boolean USE_HMEAN=false;//Harmonic mean
	public static boolean USE_GMEAN=false;//Geometric mean
	
}
