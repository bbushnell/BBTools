package cardinality;

import java.util.Arrays;
import shared.Random;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import jgi.Dedupe;
import parse.Parser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.SamLine;
import structures.LongList;
import structures.SuperLongList;
import ukmer.Kmer;

/**
 * Abstract superclass for cardinality-tracking structures like LogLog.
 * Provides probabilistic estimation of unique k-mer counts in sequences using various LogLog-based algorithms.
 * Supports different implementations optimized for accuracy, speed, or memory usage.
 * @author Brian Bushnell
 * @date Feb 20, 2020
 */
public abstract class CardinalityTracker implements Drivable {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final String pickType(String type) {
		if(trackCounts){
			if("BBLog".equalsIgnoreCase(type)) {return type;}
			else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type) || 
				"DDL".equalsIgnoreCase(type) || "DynamicDemiLog".equalsIgnoreCase(type)){
				return type;
			}else {return "LogLog16";}
		}
		return type;
	}
	
	/**
	 * Factory method that creates a tracker using default settings.
	 * Subclass is determined by static Parser.loglogType field.
	 * BBLog is preferred when trackCounts is enabled for optimal accuracy and speed.
	 * @return New CardinalityTracker instance of the configured type
	 */
	public static CardinalityTracker makeTracker(String type){
		type=pickType(Parser.loglogType);
		if("BBLog".equalsIgnoreCase(type)){
			return new BBLog();//Fastest, most accurate
		}else if("LogLog".equalsIgnoreCase(type)){
			return new LogLog();//Least accurate
		}else if("LogLog2".equalsIgnoreCase(type)){
			return new LogLog2();//Slowest, uses mantissa
		}else if("LogLog16".equalsIgnoreCase(type) || "LL16".equalsIgnoreCase(type)){
			return new LogLog16();//Uses 10-bit mantissa
//		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type)){
//			return new DynamicLogLog();//Uses 10-bit mantissa
		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type) || 
			"DDL".equalsIgnoreCase(type) || "DynamicDemiLog".equalsIgnoreCase(type)){
			return new DynamicDemiLog();//Uses 10-bit mantissa
		}else if("LogLog8".equalsIgnoreCase(type)){
			return new LogLog8();//Lowest memory
		}else if("InvertedLogLog".equalsIgnoreCase(type) || "ILL".equalsIgnoreCase(type)){
			return new InvertedLogLog();//Experimental: stores max low-bits per NLZ bucket
		}else if("DDL2".equalsIgnoreCase(type) || "DynamicDemiLog2".equalsIgnoreCase(type)){
			return new DynamicDemiLog2();
		}else if("DDL8".equalsIgnoreCase(type) || "DynamicDemiLog8".equalsIgnoreCase(type)){
			return new DynamicDemiLog8();
		}else if("DLL4".equalsIgnoreCase(type) || "DynamicLogLog4".equalsIgnoreCase(type) ||
			"DDL4".equalsIgnoreCase(type) || "DynamicDemiLog4".equalsIgnoreCase(type)){
			return new DynamicLogLog4();
		}else if("DLL3".equalsIgnoreCase(type) || "DynamicLogLog3".equalsIgnoreCase(type) ||
			"DDL3".equalsIgnoreCase(type) || "DynamicDemiLog3".equalsIgnoreCase(type)){
			return new DynamicLogLog3();
		}else if("DLL3v2".equalsIgnoreCase(type) || "DynamicLogLog3v2".equalsIgnoreCase(type)){
			return new DynamicLogLog3v2();
		}else if("DLL2".equalsIgnoreCase(type) || "DynamicLogLog2".equalsIgnoreCase(type)){
			return new DynamicLogLog2();
		}else if("LL6".equalsIgnoreCase(type) || "LogLog6".equalsIgnoreCase(type)){
			return new LogLog6();
		}
		assert(false) : "TODO: "+type;
		throw new RuntimeException(type);
	}

	/**
	 * Factory method that creates a tracker using default settings.
	 * Subclass is determined by static Parser.loglogType field.
	 * BBLog is preferred when trackCounts is enabled for optimal accuracy and speed.
	 * @return New CardinalityTracker instance of the configured type
	 */
	public static CardinalityTracker makeTracker(){return makeTracker(Parser.loglogType);}
	
	/**
	 * Factory method that creates a tracker using parsed settings.
	 * Subclass is determined by static type field.
	 * Parameters are extracted from the Parser object.
	 * @param p Parser containing configuration parameters
	 * @return New CardinalityTracker instance configured from parser
	 */
	public static CardinalityTracker makeTracker(Parser p){
		final String type=pickType(Parser.loglogType);
		if("BBLog".equalsIgnoreCase(type)){
			return new BBLog(p);
		}else if("LogLog".equalsIgnoreCase(type)){
			return new LogLog(p);
		}else if("LogLog2".equalsIgnoreCase(type)){
			return new LogLog2(p);
		}else if("LogLog16".equalsIgnoreCase(type) || "LL16".equalsIgnoreCase(type)){
			return new LogLog16(p);
//		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type)){
//			return new DynamicLogLog(p);
		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type) || 
			"DDL".equalsIgnoreCase(type) || "DynamicDemiLog".equalsIgnoreCase(type)){
			return new DynamicDemiLog(p);
		}else if("LogLog8".equalsIgnoreCase(type)){
			return new LogLog8(p);
		}else if("InvertedLogLog".equalsIgnoreCase(type) || "ILL".equalsIgnoreCase(type)){
			return new InvertedLogLog(p);//Experimental: stores max low-bits per NLZ bucket
		}else if("DDL2".equalsIgnoreCase(type) || "DynamicDemiLog2".equalsIgnoreCase(type)){
			return new DynamicDemiLog2(p);
		}else if("DDL8".equalsIgnoreCase(type) || "DynamicDemiLog8".equalsIgnoreCase(type)){
			return new DynamicDemiLog8(p);
		}else if("DLL4".equalsIgnoreCase(type) || "DynamicLogLog4".equalsIgnoreCase(type) ||
			"DDL4".equalsIgnoreCase(type) || "DynamicDemiLog4".equalsIgnoreCase(type)){
			return new DynamicLogLog4(p);
		}else if("DLL3".equalsIgnoreCase(type) || "DynamicLogLog3".equalsIgnoreCase(type) ||
			"DDL3".equalsIgnoreCase(type) || "DynamicDemiLog3".equalsIgnoreCase(type)){
			return new DynamicLogLog3(p);
		}else if("DLL3v2".equalsIgnoreCase(type) || "DynamicLogLog3v2".equalsIgnoreCase(type)){
			return new DynamicLogLog3v2(p);
		}else if("DLL2".equalsIgnoreCase(type) || "DynamicLogLog2".equalsIgnoreCase(type)){
			return new DynamicLogLog2(p);
		}else if("LL6".equalsIgnoreCase(type) || "LogLog6".equalsIgnoreCase(type)){
			return new LogLog6(p);
		}
		assert(false) : "TODO: "+type;
		throw new RuntimeException(type);
	}

	/**
	 * Factory method that creates a tracker with specified settings.
	 * Subclass is determined by static type field.
	 * Allows direct specification of all key parameters.
	 * @param buckets_ Number of buckets (will be rounded to next power of 2)
	 * @param k_ K-mer length for hashing
	 * @param seed Random number generator seed (-1 for random seed)
	 * @param minProb_ Ignore k-mers with correctness probability below this threshold
	 * @return New CardinalityTracker instance with specified configuration
	 */
	public static CardinalityTracker makeTracker(int buckets_, int k_, long seed, float minProb_){
		final String type=pickType(Parser.loglogType);
		if("BBLog".equalsIgnoreCase(type)){
			return new BBLog(buckets_, k_, seed, minProb_);
		}else if("LogLog".equalsIgnoreCase(type)){
			return new LogLog(buckets_, k_, seed, minProb_);
		}else if("LogLog2".equalsIgnoreCase(type)){
			return new LogLog2(buckets_, k_, seed, minProb_);
		}else if("LogLog16".equalsIgnoreCase(type) || "LL16".equalsIgnoreCase(type)){
			return new LogLog16(buckets_, k_, seed, minProb_);
//		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type)) {
//			return new DynamicLogLog(buckets_, k_, seed, minProb_);
		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type) || 
			"DDL".equalsIgnoreCase(type) || "DynamicDemiLog".equalsIgnoreCase(type)){
			return new DynamicDemiLog(buckets_, k_, seed, minProb_);
		}else if("LogLog8".equalsIgnoreCase(type)){
			return new LogLog8(buckets_, k_, seed, minProb_);
		}else if("InvertedLogLog".equalsIgnoreCase(type) || "ILL".equalsIgnoreCase(type)){
			return new InvertedLogLog(buckets_, k_, seed, minProb_);//Experimental: stores max low-bits per NLZ bucket
		}else if("DDL2".equalsIgnoreCase(type) || "DynamicDemiLog2".equalsIgnoreCase(type)){
			return new DynamicDemiLog2(buckets_, k_, seed, minProb_);
		}else if("DDL8".equalsIgnoreCase(type) || "DynamicDemiLog8".equalsIgnoreCase(type)){
			return new DynamicDemiLog8(buckets_, k_, seed, minProb_);
		}else if("DLL4".equalsIgnoreCase(type) || "DynamicLogLog4".equalsIgnoreCase(type) ||
			"DDL4".equalsIgnoreCase(type) || "DynamicDemiLog4".equalsIgnoreCase(type)){
			return new DynamicLogLog4(buckets_, k_, seed, minProb_);
		}else if("DLL3".equalsIgnoreCase(type) || "DynamicLogLog3".equalsIgnoreCase(type) ||
			"DDL3".equalsIgnoreCase(type) || "DynamicDemiLog3".equalsIgnoreCase(type)){
			return new DynamicLogLog3(buckets_, k_, seed, minProb_);
		}else if("DLL3v2".equalsIgnoreCase(type) || "DynamicLogLog3v2".equalsIgnoreCase(type)){
			return new DynamicLogLog3v2(buckets_, k_, seed, minProb_);
		}else if("DLL2".equalsIgnoreCase(type) || "DynamicLogLog2".equalsIgnoreCase(type)){
			return new DynamicLogLog2(buckets_, k_, seed, minProb_);
		}else if("LL6".equalsIgnoreCase(type) || "LogLog6".equalsIgnoreCase(type)){
			return new LogLog6(buckets_, k_, seed, minProb_);
		}
		assert(false) : "TODO: "+type;
		throw new RuntimeException(type);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a tracker with parameters extracted from a Parser.
	 * @param p Parser containing loglog configuration values */
	public CardinalityTracker(Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}
	
	/**
	 * Creates a tracker with specified parameters.
	 * Buckets will be rounded up to the next power of 2 for efficient bit masking.
	 * Initializes hash function with random XOR value from the specified seed.
	 * @param buckets_ Number of buckets (counters) - will be made power of 2
	 * @param k_ K-mer length for sequence hashing
	 * @param seed Random number generator seed; -1 for a random seed
	 * @param minProb_ Ignore k-mers with under this probability of being correct
	 */
	public CardinalityTracker(int buckets_, int k_, long seed, float minProb_){
		buckets=powerOf2AtLeast(buckets_);
		assert(buckets>0 && Integer.bitCount(buckets)==1) : "Buckets must be a power of 2: "+buckets;
		bucketMask=buckets-1;
		bucketBits=Long.bitCount(bucketMask);
		k=Kmer.getKbig(k_);
		minProb=minProb_;
		hashXor=(seed==-1 ? defaultSeed : seed);
//		assert(false) : seed+", "+randy+", "+hashXor;
		//Xor is not used anymore, but could be, as long as each copy gets the same one.
	}

	public abstract CardinalityTracker copy();
	
	/**
	 * Returns the lowest power of 2 that is greater than or equal to target.
	 * Required because buckets must be a power of 2 for efficient bit masking.
	 * @param target The minimum value needed
	 * @return Smallest power of 2 >= target, capped at 0x40000000
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
	
	/**
	 * Program entry point (deprecated).
	 * Redirects to LogLogWrapper for actual processing.
	 * @param args Command-line arguments
	 */
	public static final void main(String[] args){
		LogLogWrapper llw=new LogLogWrapper(args);
		
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		llw.process();
		
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
	}
	
	/**
	 * Hashes and adds a number to this tracker.
	 * This is the ONLY correct way to add elements externally.
	 * Handles added-count tracking for clampToAdded support.
	 * @param number The value to hash and track
	 */
	public final void add(long number){
		hashAndStore(number);
		added++;
	}
	
	/**
	 * Hashes and tracks all k-mers from a Read and its mate.
	 * Processes both forward and mate sequences if they meet minimum length requirements.
	 * @param r The Read to process (may be null)
	 */
	public final void hash(Read r){
		if(r==null){return;}
		if(r.length()>=k){hash(r.bases, r.quality);}
		if(r.mateLength()>=k){hash(r.mate.bases, r.mate.quality);}
	}
	
	/**
	 * Hashes and tracks all k-mers from a Read and its mate.
	 * Processes both forward and mate sequences if they meet minimum length requirements.
	 * @param r The Read to process (may be null)
	 */
	public final void hash(SamLine r){
		if(r==null || r.seq==null){return;}
		if(r.length()>=k){hash(r.seq, r.qual);}
	}
	
	/**
	 * Hashes and tracks all k-mers from a sequence with quality scores.
	 * Routes to appropriate method based on k-mer size (small vs big).
	 * @param bases Sequence bases as byte array
	 * @param quals Quality scores (may be null)
	 */
	public final void hash(byte[] bases, byte[] quals){
		if(k<32){hashSmall(bases, quals);}
		else{hashBig(bases, quals);}
	}
	
	/**
	 * Hashes and tracks sequence using short k-mers (k < 32).
	 * Uses bit-shifting operations for efficient k-mer rolling.
	 * Applies quality-based probability filtering when quality scores provided.
	 * Uses canonical k-mers (maximum of forward and reverse complement).
	 * @param bases Sequence bases as byte array
	 * @param quals Quality scores for probability calculation (may be null)
	 */
	public final void hashSmall(byte[] bases, byte[] quals){
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		final int lenmask=(-1)>>>k;
		final int lenmask2=lenmask>>1;
		
		long kmer=0, rkmer=0;
		
		if(minProb>0 && quals!=null){//Debranched loop
			assert(quals.length==bases.length) : quals.length+", "+bases.length;
			float prob=1;
			int len=0;
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				
//				long x=AminoAcid.baseToNumber[b];
//				final long x2=x^3;
				
				final long x=((b>>1)&3)|((b&8)<<28);//Different encoding
				final long x2=x^2;
				
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

				//Update probability
				prob=prob*PROB_CORRECT[quals[i]];
				if(len>=k){prob=prob*PROB_CORRECT_INVERSE[quals[i-k]];}
				
				if(x>=0){
					len++;
				}else{
					len=0;
					kmer=rkmer=0;
					prob=1;
				}
				
				if(len>=k && prob>=minProb){add(Math.max(kmer, rkmer));}
			}
		}else{

			//Old, slower logic
//			for(int i=0; i<bases.length; i++){
//				byte b=bases[i];
//				long x=AminoAcid.baseToNumber[b];
//				long x2=AminoAcid.baseToComplementNumber[b];
//				kmer=((kmer<<2)|x)&mask;
//				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
//				
//				if(x>=0){
//					len++;
//				}else{
//					len=0;
//					kmer=rkmer=0;
//				}
//				if(len>=k){
//					add(Math.max(kmer, rkmer));
//				}
//			}
			
			//Does not handle IUPAC
			int len=-1;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
//				final long x=AminoAcid.baseToNumber[b];//Normal 2-bit encoding
//				final long x2=x^3;
				final long x=((b>>1)&3)|((b&8)<<28);//Alternate 2-bit encoding
				final long x2=x^2;
				
				// shift register, tracks valid bases in upper zeros
				len|=x;
				len>>>=1;
				
				// Let the kmers roll. Garbage flushes out naturally after k shifts!
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				
				// Countdown timer reached the threshold
				if(len<=lenmask){
					add(Math.max(kmer, rkmer));
				}
			}
		}
	}
	
	/**
	 * Hashes and tracks sequence using long k-mers (k >= 32).
	 * Uses Kmer objects for handling k-mers longer than 31 bases.
	 * Applies quality-based probability filtering when quality scores provided.
	 * @param bases Sequence bases as byte array
	 * @param quals Quality scores for probability calculation (may be null)
	 */
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
				if(len>=k){
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
	 * Makes a table of random bitmasks for hashing (deprecated).
	 * Creates random bit patterns with controlled bit density for XOR operations.
	 * Short-circuited to return null as this method is superseded.
	 * @param length Number of hash rounds
	 * @param bits Bits per hash step
	 * @param seed Random seed for table generation
	 * @return Random bitmask table (currently returns null)
	 */
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
	
	public final float compensationFactorLogBuckets(int logBuckets){
		float[] array=compensationFactorLogBucketsArray();
		return (array!=null && logBuckets<array.length) ? array[logBuckets] : 1/(1+(1<<logBuckets));
	}
	
	public SuperLongList toFrequency(){
		SuperLongList list=new SuperLongList(1000);
		if(getClass()==LogLog16.class) {
			char[] counts=counts16();
			for(char x : counts){
				if(x>0){list.add(x);}
			}
		}else {
			int[] counts=getCounts();
			for(int x : counts){
				if(x>0){list.add(x);}
			}
		}
		list.sort();
		return list;
	}

	/**
	 * Prints a k-mer frequency histogram to file.
	 * Outputs depth-count pairs with optional supersampling adjustment.
	 * Handles both array-based and list-based frequency data.
	 * @param path File path for output
	 * @param overwrite Whether to overwrite existing file
	 * @param append Whether to append to existing file
	 * @param supersample Adjust counts for the effect of subsampling
	 * @param decimals Number of decimal places for supersampled counts
	 */
	public void printKhist(String path, boolean overwrite, boolean append, boolean supersample, int decimals){
		if(this.getClass()==LogLog16.class) {
			System.err.println("a");
			printKhist32(path, overwrite, append, supersample, decimals);
		}else {
			System.err.println("b");
			printKhist32(path, overwrite, append, supersample, decimals);
		}
	}
	
	public void printKhist32(String path, boolean overwrite, boolean append, boolean supersample, int decimals){
		SuperLongList sll=toFrequency();
		ByteStreamWriter bsw=new ByteStreamWriter(path, overwrite, append, false);
		bsw.start();
		bsw.print("#Depth\tCount\n");
		final double mult=Math.max(1.0, (supersample ? cardinality()/(double)buckets : 1));
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
						bsw.print(Math.max(1, Math.round(count*mult))).nl();
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
						bsw.print(Math.max(1, Math.round(count*mult))).nl();
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
					bsw.print(Math.max(1, Math.round(count*mult))).nl();
				}
			}else{
				bsw.print(count).nl();
			}
		}
		bsw.poisonAndWait();
	}
	
//	public void printKhist16(String path, boolean overwrite, boolean append, boolean supersample, int decimals){
//		char[] array=counts16();
//		ByteStreamWriter bsw=new ByteStreamWriter(path, overwrite, append, false);
//		bsw.start();
//		bsw.print("#Depth\tCount\n");
//		final double mult=Math.max(1.0, (supersample ? cardinality()/(double)buckets : 1));
//		
//		for(int depth=0; depth<array.length; depth++){
//			long count=array[depth];
//			if(count>0){
//				bsw.print(depth).tab();
//				if(supersample){
//					if(decimals>0){
//						bsw.print(count*mult, decimals).nl();
//					}else{
//						bsw.print(Math.max(1, Math.round(count*mult))).nl();
//					}
//				}else{
//					bsw.print(count).nl();
//				}
//			}
//		}
//		bsw.poisonAndWait();
//	}
	
	public final long countSum(){
		int[] counts=getCounts();
		return counts==null ? 0 : simd.Vector.sum(counts);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Calculates cardinality estimate from this tracker.
	 * Implementation varies by subclass algorithm.
	 * @return Estimated number of unique k-mers observed
	 */
	public abstract long cardinality();

	/**
	 * Returns the counts array if present.
	 * Should be overridden for classes that track counts.
	 * @return Array of bucket counts, or null if not tracking counts
	 */
	public int[] getCounts(){
		return null;
	}
	
	public char[] counts16(){return null;}
	
	/**
	 * Merges another tracker into this one.
	 * Combines cardinality estimates from both trackers.
	 * @param log The tracker to add to this one
	 */
	public abstract void add(CardinalityTracker log);
	
	/**
	 * Generates a 64-bit hashcode from a number and adds it to this tracker.
	 * <p>
	 * <b>DO NOT CALL DIRECTLY.</b> Use {@link #add(long)} instead, which calls
	 * this method and also increments the {@code added} counter needed for
	 * {@code clampToAdded} support. Calling hashAndStore directly will leave
	 * {@code added=0}, causing {@code cardinality()} to return 0 when clamping
	 * is enabled (the default).
	 * <p>
	 * This method exists as the subclass extension point — subclasses override
	 * this to implement their specific hashing and storage logic.
	 *
	 * @param number The value to hash and store
	 */
	public abstract void hashAndStore(final long number);
	
	/**
	 * Returns array of compensation factors indexed by log2(buckets).
	 * Designed to compensate for overestimate with small numbers of buckets.
	 * May be deprecated in favor of harmonic mean handling multiple trials.
	 * @return Array of compensation factors, or null if not implemented
	 */
	public abstract float[] compensationFactorLogBucketsArray();

	/*--------------------------------------------------------------*/
	/*----------------       Drivable Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Default: throws UnsupportedOperationException. Override in calibratable subclasses. */
	@Override
	public double[] rawEstimates(){throw new UnsupportedOperationException(getClass().getSimpleName()+" does not support rawEstimates()");}

	/** Default: throws UnsupportedOperationException. Override in calibratable subclasses. */
	@Override
	public int filledBuckets(){throw new UnsupportedOperationException(getClass().getSimpleName()+" does not support filledBuckets()");}

	/** Default: throws UnsupportedOperationException. Override in calibratable subclasses. */
	@Override
	public double occupancy(){throw new UnsupportedOperationException(getClass().getSimpleName()+" does not support occupancy()");}

	@Override
	public long getLastCardinality(){return lastCardinality;}

	@Override
	public void setLastCardinality(long val){lastCardinality=val;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** K-mer length for sequence hashing */
	public final int k;
	
	/** Minimum probability threshold for k-mer correctness filtering */
	public final float minProb;
	
	/** Number of buckets for tracking; must be a power of 2 for efficiency */
	public final int buckets;
	
	public final int bucketBits;
	
	/** Bit mask for extracting bucket index from hashcode; equals buckets-1 */
	final int bucketMask;
	
	/** Thread-local storage for Kmer objects used in long k-mer hashing */
	private final ThreadLocal<Kmer> localKmer=new ThreadLocal<Kmer>();
	
	/** Value XORed with hash inputs to vary hash function behavior */
	protected final long hashXor;
	
	/**
	 * Gets a thread-local Kmer object for long k-mer mode processing.
	 * Creates new Kmer if none exists for current thread.
	 * Clears the Kmer before returning for reuse.
	 * @return Thread-local Kmer object ready for use
	 */
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
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	long added=0;
	long microIndex=0;
	/** Absolute NLZ histogram: nlzCounts[k] = number of buckets with absoluteNlz==k.
	 * Lazy-allocated on first summarize() call; cleared and reused on subsequent calls.
	 * Used by DynamicLC (DLC) estimator in CardinalityStats. Null for non-DLL/DDL classes. */
	protected int[] nlzCounts;
	/** Sorted dif-value buffer for median, MWA, and trimmed-mean estimators.
	 * Lazy-allocated on first summarize() when USE_SORTBUF is true; null otherwise.
	 * Gated by USE_SORTBUF so HotSpot eliminates the dead code when disabled. */
	protected LongList sortBuf;
	/** Cached cardinality estimate; -1 means stale. Shared by all calibratable subclasses. */
	public long lastCardinality=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------            Statics           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Array converting quality scores to probability of base correctness */
	public static final float[] PROB_CORRECT=Arrays.copyOf(align2.QualityTools.PROB_CORRECT, 128);
	public static final float[] PROB_CORRECT_INVERSE=Arrays.copyOf(align2.QualityTools.PROB_CORRECT_INVERSE, 128);
	
	/** Whether to track occurrence counts for each bucket value */
	public static boolean trackCounts=false;
	public static boolean clampToAdded=true;
	/** Records the most recent cardinality estimate for static contexts */
	public static long lastCardinalityStatic=-1;
//	/** Ignore hashed values above this, to skip expensive read and store functions. */
//	static final long maxHashedValue=((-1L)>>>3);//No longer used
	
	public static long defaultSeed=0;
	
	/** Whether to use arithmetic mean for combining multiple bucket estimates */
	public static boolean USE_MEAN=true;//Arithmetic mean of inverted differences
	public static boolean USE_MEDIAN=false;
	public static boolean USE_MWA=false;//Median-weighted-average
	public static boolean USE_HMEAN=false;//Harmonic mean
	public static boolean USE_HMEANM=false;//Harmonic mean including mantissa
	public static boolean USE_GMEAN=false;//Geometric mean
	public static boolean USE_HLL=false;//HLL formula
	public static boolean USE_HYBRID=false;//Hybrid of LC and Mean
	public static final boolean USE_MICRO=false;
	/** When true, allocate and fill sortBuf for median/MWA/mean99 estimators.
	 * When false (default), those estimators return fallback values and no LongList is allocated. */
	public static final boolean USE_SORTBUF=false;
	/** When true, use cardinality-indexed CF table in addition to occupancy-indexed table.
	 * Only applies when CorrectionFactor.USE_CORRECTION is also true; cardinality CF
	 * cannot be active without the occupancy CF. Default true. */
	public static boolean USE_CARD_CF=true;
	public static final int MICRO_CUTOFF_BITS=56;//Higher is less accurate, max is 64
	
}
