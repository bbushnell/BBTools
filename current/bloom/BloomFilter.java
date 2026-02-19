package bloom;

import java.io.File;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;
import stream.Read;
import structures.IntList;
import structures.LongList;

/**
 * Wraps a KCountArray and provides multithreaded reference loading.
 * 
 * @author Brian Bushnell
 * @date April 23, 2018
 *
 */
public class BloomFilter implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -3987955563503838492L;
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		BloomFilter x=new BloomFilter(args);
		
		System.err.println(x.filter.toShortString());
		t.stop("Time: \t");
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public BloomFilter(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		KmerCountAbstract.CANONICAL=true;
		
		//Create a parser object
		Parser parser=new Parser();

		int k_=31;
		int kbig_=31;
		int bits_=1;
		int hashes_=2;
		int minConsecutiveMatches_=3;
		float memFraction=1;
		boolean rcomp_=true;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("k") || a.equals("ksmall")){
				k_=Integer.parseInt(b);
				assert(k_<=31 && k_>=1);
			}else if(a.equals("kbig")){
				kbig_=Integer.parseInt(b);
				assert(kbig_>=1);
			}else if(a.equals("hashes")){
				hashes_=Integer.parseInt(b);
				assert(hashes_<=10000 && hashes_>=1);
			}else if(a.equals("minhits")){
				minConsecutiveMatches_=Integer.parseInt(b);
				assert(minConsecutiveMatches_>=1);
			}else if(a.equals("bits")){
				bits_=Integer.parseInt(b);
			}else if(a.equals("memfraction")){
				memFraction=Float.parseFloat(b);
			}else if(a.equals("extra")){
				if(b==null){extra.clear();}
				else{
					for(String s : b.split(",")){extra.add(s);}
				}
			}else if(a.equals("rcomp")){
				rcomp_=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		k=k_;
		kbig=Tools.max(k_, kbig_);
		smallPerBig=kbig-k+1;
		bits=bits_;
		hashes=hashes_;
		minConsecutiveMatches=minConsecutiveMatches_;
		rcomp=rcomp_;
		

		assert(bits==1 || bits==2 || bits==4 || bits==8 || bits==16 || bits==32) : "Bits must be a power of 2.";
		
		{//Process parser fields
			in1=parser.in1;
			in2=parser.in2;
		}
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		filterMemory=setMemory(memFraction);
		
		shift=bitsPerBase*k;
		shift2=shift-bitsPerBase;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		filter=load();
	}
	
	/**
	 * Constructor with direct parameter specification.
	 * Creates a BloomFilter with specified k-mer parameters and memory allocation.
	 *
	 * @param k_ Small k-mer length
	 * @param kbig_ Large k-mer length
	 * @param bits_ Bits per k-mer count
	 * @param hashes_ Number of hash functions
	 * @param minConsecutiveMatches_ Minimum consecutive k-mer matches required
	 * @param rcomp_ Use reverse complement matching
	 * @param ecco_ Enable error correction
	 * @param merge_ Merge overlapping reads
	 * @param memFraction Fraction of available memory to use
	 */
	public BloomFilter(int k_, int kbig_, int bits_, int hashes_,
			int minConsecutiveMatches_, boolean rcomp_, boolean ecco_, boolean merge_, float memFraction){
		this(null, null, null, k_, kbig_, bits_, hashes_, minConsecutiveMatches_, rcomp_, ecco_, merge_, memFraction);
	}
	
	/**
	 * Constructor with input files and parameter specification.
	 * Creates a BloomFilter from specified input files with given parameters.
	 *
	 * @param in1_ Primary input file path
	 * @param in2_ Secondary input file path
	 * @param extra_ Additional input file paths
	 * @param k_ Small k-mer length
	 * @param kbig_ Large k-mer length
	 * @param bits_ Bits per k-mer count
	 * @param hashes_ Number of hash functions
	 * @param minConsecutiveMatches_ Minimum consecutive k-mer matches required
	 * @param rcomp_ Use reverse complement matching
	 * @param ecco_ Enable error correction
	 * @param merge_ Merge overlapping reads
	 * @param memFraction Fraction of available memory to use
	 */
	public BloomFilter(String in1_, String in2_, ArrayList<String> extra_, int k_, int kbig_, int bits_, int hashes_,
			int minConsecutiveMatches_, boolean rcomp_, boolean ecco_, boolean merge_, float memFraction){
		if(extra_!=null){
			for(String s : extra_){extra.add(s);}
		}
		filterMemory=setMemory(memFraction);
		
		in1=in1_;
		in2=in2_;
		k=k_;
		kbig=Tools.max(k_, kbig_);
		smallPerBig=kbig-k+1;
		bits=bits_;
		hashes=hashes_;
		minConsecutiveMatches=minConsecutiveMatches_;
		rcomp=rcomp_;
		ecco=ecco_;
		merge=merge_;

		shift=bitsPerBase*k;
		shift2=shift-bitsPerBase;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		filter=load();
	}

	/**
	 * Constructor for loading from BBMap index.
	 * Creates a BloomFilter by loading k-mer data from an existing BBMap index.
	 *
	 * @param bbmapIndex_ Must be true to indicate BBMap index loading
	 * @param k_ Small k-mer length
	 * @param kbig_ Large k-mer length
	 * @param bits_ Bits per k-mer count
	 * @param hashes_ Number of hash functions
	 * @param minConsecutiveMatches_ Minimum consecutive k-mer matches required
	 * @param rcomp_ Use reverse complement matching
	 */
	public BloomFilter(boolean bbmapIndex_, int k_, int kbig_, int bits_, int hashes_, int minConsecutiveMatches_, boolean rcomp_) {
		assert(bbmapIndex_);
		filterMemory=setMemory(0.75);
		
		in1=null;
		in2=null;
		k=k_;
		kbig=Tools.max(k_, kbig_);
		smallPerBig=kbig-k+1;
		bits=bits_;
		hashes=hashes_;
		minConsecutiveMatches=minConsecutiveMatches_;
		rcomp=rcomp_;

		shift=bitsPerBase*k;
		shift2=shift-bitsPerBase;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		filter=loadFromIndex();
	}
	
	/**
	 * Calculates available memory for filter allocation.
	 * Determines usable memory based on JVM settings and applies multiplier.
	 * @param mult Fraction of available memory to allocate
	 * @return Memory in bytes available for filter
	 */
	private static long setMemory(double mult){
		if(printMem) {Shared.printMemory();}
		
		Runtime rt=Runtime.getRuntime();
		final long mmemory=rt.maxMemory();
		final long tmemory=rt.totalMemory();
		final long fmemory=rt.freeMemory();
		final long umemory=tmemory-fmemory;
		
		double xmsRatio=Shared.xmsRatio();
		double usableMemory=Tools.max(((mmemory-96000000)*(xmsRatio>0.97 ? 0.82 : 0.72)), mmemory*0.45);
		double availableMemory=usableMemory-umemory;
		double filterMemory=availableMemory*mult;
		
//		System.err.println((long)(usableMemory/1000000)+", "+(long)(availableMemory/1000000)+", "+(long)(filterMemory/1000000));
		
		return (long)filterMemory;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads k-mer data from input files into a KCountArray.
	 * Creates and populates the count array from reference sequences.
	 * @return Populated KCountArray7MTA with k-mer counts
	 */
	private KCountArray7MTA load(){
		final int cbits=bits;
		final long totalBits=8*filterMemory;
		final long cells=(OVERRIDE_CELLS>0 ? OVERRIDE_CELLS : totalBits/cbits);
		//		System.err.println("filterMemory="+filterMemory+", cells="+cells);
		KCountArray7MTA kca;
		if(in1!=null || (extra!=null && extra.size()>0)) {
			ReadCounter rc=new ReadCounter(k, rcomp, ecco, merge, Shared.AMINO_IN);
			kca=(KCountArray7MTA)rc.makeKca(in1, in2, extra==null || extra.isEmpty() ? null : extra,
					cbits, cells, hashes, minq,
					maxReads, 1, 1, 1, null, 0);
		}else {
			kca=(KCountArray7MTA) KCountArray.makeNew(cells, cbits, hashes, null, 0);
		}
		return kca;
	}

	/**
	 * Loads k-mer data from an existing BBMap index.
	 * Creates a KCountArray from pre-built index rather than raw sequences.
	 * @return KCountArray7MTA populated from index data
	 */
	private KCountArray7MTA loadFromIndex(){
		KmerCountAbstract.CANONICAL=true;
		final int cbits=bits;
		final long totalBits=8*filterMemory;
		final long cells=totalBits/cbits;
		outstream.println("Filter Memory = "+Tools.format("%.2f GB", filterMemory/(double)(1024*1024*1024)));
		IndexCounter ic=new IndexCounter(k, rcomp);
		KCountArray7MTA kca=(KCountArray7MTA)ic.makeKcaFromIndex(cells, cbits, hashes);
		return kca;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Determines if both reads pass the filter threshold.
	 * Both reads must individually pass for the pair to pass.
	 *
	 * @param r1 First read
	 * @param r2 Second read
	 * @param thresh Minimum k-mer count threshold
	 * @return true if both reads pass, false otherwise
	 */
	public boolean passes(Read r1, Read r2, final int thresh) {
		boolean pass=passes(r1, thresh);
		return pass && passes(r2, thresh);
	}
	
	/**
	 * Calculates average k-mer count across a sequence.
	 * Uses smoothing to reduce impact of hash collision spikes.
	 * @param bases Sequence bases to analyze
	 * @return Average k-mer count, or 0 if sequence too short
	 */
	public float averageCount(final byte[] bases) {
		if(bases==null || bases.length<k-1){return 0;}

		long kmer=0;
		long rkmer=0;
		long sum=0;
		int len=0;
		int counted=0;
		
		int prev2=0, prev=0, count=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
				prev2=prev;
				prev=count;
				count=getCount(kmer, rkmer);
				//This should get rid of collision spikes.
				prev=Tools.min(prev, Tools.max(prev2, count));
				sum+=prev;
				counted++;
			}
		}
		sum+=count; //Last kmer did not get counted; a 0 was used instead.
		return sum/Tools.max(counted, 1f);
	}
	
	/**
	 * Finds the minimum k-mer count in a read.
	 * Useful for detecting contamination or low-coverage regions.
	 * @param r Read to analyze
	 * @return Minimum k-mer count, or -1 if read too short
	 */
	public int minCount(Read r) {
		if(r==null || r.length()<k-1){return -1;}
		final byte[] bases=r.bases;

		long kmer=0;
		long rkmer=0;
		int len=0;
		int min=Integer.MAX_VALUE;
		int counted=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
				counted++;
				int count=getCount(kmer, rkmer);
				min=Tools.min(min, count);
				if(count==0){
//					assert(false) : counted+", "+min;
					break;
				}
			}
		}
//		assert(false) : counted+", "+min;
		return counted>0 ? min : -1;
	}
	
	/**
	 * Tests if a read has sufficient high-count k-mers.
	 * Determines if at least the specified fraction of k-mers exceed threshold.
	 *
	 * @param r Read to evaluate
	 * @param thresh Count threshold for "high" k-mers
	 * @param fraction Minimum fraction of k-mers that must be high
	 * @return true if read meets high-count fraction requirement
	 */
	public boolean hasHighCountFraction(Read r, final int thresh, final float fraction) {
		if(r==null || r.length()<k-1){return false;}
		final byte[] bases=r.bases;
		final int kmers=r.length()-k+1;
		
		final int minHigh=Math.round(fraction*kmers);
		final int maxLow=kmers-minHigh;

		long kmer=0;
		long rkmer=0;
		int len=0;
		int counted=0;
		int low=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
				counted++;
				int count=getCount(kmer, rkmer);
				if(count<thresh){
					low++;
					if(low>maxLow){return false;}
				}
			}
		}
		return true;
	}
	
	/**
	 * Returns fraction of k-mers with count below threshold.
	 * Complement of highCountFraction.
	 *
	 * @param r Read to examine
	 * @param thresh Count threshold
	 * @param smooth Apply count smoothing
	 * @return Fraction of low-count k-mers
	 */
	public float lowCountFraction(final Read r, final int thresh, final boolean smooth) {
		return 1-highCountFraction(r, thresh, smooth);
	}
	
	/**
	 * Return the fraction of kmers with count at least thresh.
	 * @param r Read to examine.
	 * @param thresh Minimum count to be 'high'.
	 * @param smooth Reduce each count to min(count, prevCount).
	 * @return High count fraction, or zero if no valid kmers.
	 */
	public float highCountFraction(final Read r, final int thresh, final boolean smooth) {
		return r==null ? 0 : highCountFraction(r.bases, thresh, smooth);
	}
	
	/**
	 * Return the fraction of kmers with count at least thresh.
	 * @param bases Bases to examine.
	 * @param thresh Minimum count to be 'high'.
	 * @param smooth Reduce each count to min(count, prevCount).
	 * @return High count fraction, or zero if no valid kmers.
	 */
	public float highCountFraction(final byte[] bases, final int thresh, boolean smooth) {
		if(bases==null || bases.length<k-1){return 0;}

		long kmer=0;
		long rkmer=0;
		int len=0;
		int counted=0;
		int highCount=0;
		int prevHigh=1;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
				counted++;
				final int count=getCount(kmer, rkmer);
				final int high=count>=thresh ? 1 : 0;
				highCount+=(high&prevHigh);
				prevHigh=(smooth ? high : prevHigh);
			}
		}
		return counted<1 ? 0 : highCount/(float)counted;
	}
	
	/**
	 * Determines if a read pair represents junk/contamination.
	 * Tests k-mer counts at read ends within specified range.
	 *
	 * @param r1 First read
	 * @param r2 Second read
	 * @param range Number of bases from ends to examine
	 * @return true if reads appear to be junk/contamination
	 */
	public boolean isJunk(Read r1, Read r2, int range){
		assert(bits>1);
		if(r2==null || r2.length()<k){return isJunk(r1, range);}
		if(r1.length()<k){return isJunk(r2, range);}
		if(getLeftCount(r1.bases, range)>1 || getLeftCount(r2.bases, range)>1){return false;}
//		return getRightCount(r1.bases, range)<3 && getRightCount(r2.bases, range)<3; //&& is more correct; || allows for a fuller filter. 
		return getRightCount(r1.bases, range)<3 || getRightCount(r2.bases, range)<3;
	}
	
	/**
	 * Determines if a single read represents junk/contamination.
	 * Tests k-mer counts at both ends within specified range.
	 *
	 * @param r Read to evaluate
	 * @param range Number of bases from ends to examine
	 * @return true if read appears to be junk/contamination
	 */
	public boolean isJunk(Read r, int range){
		assert(bits>1);
		if(r.length()<k){return true;}
		return getLeftCount(r.bases, range)<2 && getRightCount(r.bases, range)<2;
	}
	
	/**
	 * Gets minimum k-mer count from left end of sequence.
	 * Examines k-mers within specified range from start.
	 *
	 * @param bases Sequence bases
	 * @param range Number of bases from start to examine
	 * @return Minimum k-mer count in left region, or -1 if too short
	 */
	private int getLeftCount(byte[] bases, int range){
		assert(range>0) : range;
		if(bases.length<k){return -1;}
		long kmer=0, rkmer=0;
		int len=0;
		final int stop=Tools.min(bases.length, k+range-1);
		int min=Integer.MAX_VALUE;
		int counted=0;
		for(int i=0; i<stop; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(x<0){len=0; rkmer=0;}
			else{
				len++;
				if(len>=k){
					counted++;
					int count=getCount(kmer, rkmer);
					min=Tools.min(min, count);
				}
			}
		}
		return counted>0 ? min : -1;
	}
	
	/**
	 * Gets minimum k-mer count from right end of sequence.
	 * Examines k-mers within specified range from end.
	 *
	 * @param bases Sequence bases
	 * @param range Number of bases from end to examine
	 * @return Minimum k-mer count in right region, or -1 if too short
	 */
	private int getRightCount(byte[] bases, int range){
		assert(range>0) : range;
		if(bases.length<k){return -1;}
		long kmer=0, rkmer=0;
		int len=0;
		final int start=Tools.max(0, bases.length-k-range+1);
		int min=Integer.MAX_VALUE;
		int counted=0;
		for(int i=start; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(x<0){len=0; rkmer=0;}
			else{
				len++;
				if(len>=k){
					counted++;
					int count=getCount(kmer, rkmer);
					min=Tools.min(min, count);
				}
			}
		}
		return counted>0 ? min : -1;
	}
	
	/**
	 * Determines if a read passes the consecutive k-mer match filter.
	 * Rejects reads with too many consecutive high-count k-mers.
	 *
	 * @param r Read to evaluate
	 * @param thresh Count threshold for k-mer matches
	 * @return true if read passes (not contaminated), false if rejected
	 */
	public boolean passes(Read r, final int thresh) {
		if(r==null || r.length()<k+minConsecutiveMatches-1){return true;}
		final byte[] bases=r.bases;

		long kmer=0;
		long rkmer=0;
		int len=0;
		int streak=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
				boolean found=contains(kmer, rkmer, thresh);
				if(found){
					streak++;
					if(streak>=minConsecutiveMatches){return false;}
				}else{streak=0;}
			}
		}
		return true;
	}
	
	/**
	 * Tests if a read matches filter criteria using precomputed keys.
	 * Inverse of passes() method.
	 *
	 * @param r Read to test
	 * @param keys Precomputed k-mer keys
	 * @param thresh Count threshold
	 * @return true if read matches (should be filtered)
	 */
	public boolean matches(Read r, LongList keys, final int thresh) {
		return !passes(r, keys, thresh);
	}
	
	/**
	 * Tests if either read in a pair matches filter criteria.
	 * Returns true if at least one read matches.
	 *
	 * @param r1 First read
	 * @param r2 Second read
	 * @param keys Precomputed k-mer keys
	 * @param thresh Count threshold
	 * @return true if either read matches filter
	 */
	public boolean matchesEither(Read r1, Read r2, LongList keys, final int thresh) {
		boolean match=!passes(r1, keys, thresh);
		return match || !passes(r2, keys, thresh);
	}
	
	/**
	 * Tests if both reads in a pair pass filter using precomputed keys.
	 * Both reads must pass for the pair to pass.
	 *
	 * @param r1 First read
	 * @param r2 Second read
	 * @param keys Precomputed k-mer keys
	 * @param thresh Count threshold
	 * @return true if both reads pass filter
	 */
	public boolean passes(Read r1, Read r2, LongList keys, final int thresh) {
		boolean pass=passes(r1, keys, thresh);
		return pass && passes(r2, keys, thresh);
	}
	
	/**
	 * Tests if a read passes filter using precomputed k-mer keys.
	 * More efficient than recomputing k-mers for repeated testing.
	 *
	 * @param r Read to test
	 * @param keys List to store/reuse k-mer keys
	 * @param thresh Count threshold
	 * @return true if read passes consecutive match filter
	 */
	public boolean passes(Read r, LongList keys, final int thresh) {
		if(r==null || r.length()<k+minConsecutiveMatches-1){return true;}
		if(minConsecutiveMatches<2){return passes(r, thresh);}
		keys.clear();
		final byte[] bases=r.bases;

		long kmer=0;
		long rkmer=0;
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(x<0){
				if(len>=k){keys.add(-1);}
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
				if(len>=k){keys.add(toKey(kmer, rkmer));}
			}
		}
		return passes(keys, thresh);
	}
	
	/**
	 * Tests if precomputed k-mer keys pass consecutive match filter.
	 * Core filtering logic using key array for efficiency.
	 *
	 * @param keys List of k-mer keys to evaluate
	 * @param thresh Count threshold for matches
	 * @return true if keys pass consecutive match criteria
	 */
	public boolean passes(final LongList keys, final int thresh) {
		assert(minConsecutiveMatches>1);
		final long[] array=keys.array;
		final int len=keys.size;
		
		for(int i=minConsecutiveMatches-1; i<len; i+=minConsecutiveMatches){
			final boolean found;
			{
				final long key=array[i];
				found=(key<0 ? false : filter.read(key)>=thresh);
			}
			if(found){
				int streak=1;
				for(int j=1; j<minConsecutiveMatches; j++){
					final long key=array[i-j];
					if(key<0 || filter.read(key)<thresh){break;}
					else{streak++;}
				}
				if(streak>=minConsecutiveMatches){return false;}
				for(int j=1; j<minConsecutiveMatches && j+i<len; j++){
					final long key=array[i+j];
					if(key<0 || filter.read(key)<thresh){break;}
					else{streak++;}
					if(streak>=minConsecutiveMatches){return false;}
				}
			}
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Returns number of counts */
	public int fillCounts(byte[] bases, IntList counts){
		final int blen=bases.length;
		if(blen<k){return 0;}
		final int min=k-1;
		long kmer=0, rkmer=0;
		int len=0;
		int valid=0;

		counts.clear();

		/* Loop through the bases, maintaining a forward kmer via bitshifts */
		for(int i=0; i<blen; i++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			final long x2=AminoAcid.baseToComplementNumber[base];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
			}
			
			if(i>=min){
				if(len>=k){
					int count=getCount(kmer, rkmer);
					counts.add(count);
					valid++;
				}else{
					counts.add(0);
				}
			}
		}
		return valid;
	}
	
	/**
	 * Gets count for a k-mer/reverse-complement pair.
	 * Uses canonical k-mer representation if reverse complement enabled.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @return Count value from filter
	 */
	public int getCount(final long kmer, final long rkmer){
		final long key=toKey(kmer, rkmer);
		return filter.read(key);
	}
	
	/**
	 * Gets count for a canonical k-mer key.
	 * Direct lookup in the count array.
	 * @param key Canonical k-mer key
	 * @return Count value from filter
	 */
	public int getCount(final long key){
//		assert(key==toKey(key, AminoAcid.reverseComplementBinaryFast(key, k))); //slow
		return filter.read(key);
	}
	
	/**
	 * Tests if k-mer pair has count meeting threshold.
	 * Convenience method combining getCount with threshold test.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param thresh Count threshold
	 * @return true if k-mer count >= threshold
	 */
	public boolean contains(final long kmer, final long rkmer, final int thresh){
		final long key=toKey(kmer, rkmer);
		return filter.read(key)>=thresh;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Returns number of counts */
	public int fillCountsBig(byte[] bases, IntList counts){
		assert(smallPerBig>1) : smallPerBig;
		final int valid0=fillCounts(bases, counts);
		if(valid0<smallPerBig){return 0;}
//		System.err.println(counts.size);
		for(int i=0, lim=counts.size()-smallPerBig+1; i<lim; i++){
			int count=smallToBig(counts, i);
			counts.set(i, count);
		}
		counts.size-=(smallPerBig-1);
//		System.err.println(counts.size+", "+k+", "+kbig+", "+smallPerBig);
		return valid0-smallPerBig+1; //Normally correct
	}
	
	
	/**
	 * Fills counts for precomputed big k-mers.
	 * More efficient when k-mers are already available.
	 * @param kmers List of big k-mers
	 * @param counts List to fill with corresponding counts
	 */
	public void fillCountsBig(LongList kmers, IntList counts){
		assert(smallPerBig>1) : smallPerBig;
		counts.clear();
		for(int i=0; i<kmers.size; i++){
			long kmer=kmers.get(i);
			int count=getCountBig(kmer);
			counts.add(count);
		}
//		assert(false) : counts;
	}
	
	/**
	 * Converts small k-mer counts to big k-mer count.
	 * Takes minimum count across constituent small k-mers.
	 *
	 * @param counts Array of small k-mer counts
	 * @param start Starting position in counts array
	 * @return Big k-mer count (minimum of small k-mer counts)
	 */
	private int smallToBig(IntList counts, final int start){
		assert(smallPerBig>1) : smallPerBig;
		final int[] array=counts.array;
		int min=array[start];
		for(int i=start+1; i<start+smallPerBig && min>0; i++){
			min=Tools.min(min, array[i]);
		}
		return min;
	}
	
	/**
	 * Gets big k-mer count for k-mer/reverse-complement pair.
	 * Delegates to single-parameter version.
	 *
	 * @param kmer Forward k-mer (unused)
	 * @param rkmer Reverse complement k-mer (unused)
	 * @return Big k-mer count
	 */
	@SuppressWarnings("unused")
	public int getCountBig(final long kmer, final long rkmer){
		return getCountBig(kmer);
	}
	
	/**
	 * Gets big k-mer count by breaking into constituent small k-mers.
	 * Returns minimum count across all small k-mers in the big k-mer.
	 * @param kmer Big k-mer to analyze
	 * @return Minimum count among constituent small k-mers
	 */
	public int getCountBig(long kmer){
		int min=Integer.MAX_VALUE;
		for(int i=0; i<smallPerBig && min>0; i++){
			long small=kmer&mask;
			long key=toKey(small);
			int count=getCount(key);
			min=Tools.min(min, count);
			kmer>>=bitsPerBase;
		}
		return min;
	}
	
	/**
	 * Tests if big k-mer meets count threshold.
	 * Uses same logic as regular contains but for big k-mers.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param thresh Count threshold
	 * @return true if big k-mer count >= threshold
	 */
	public boolean containsBig(final long kmer, final long rkmer, final int thresh){
		final long key=toKey(kmer, rkmer);
		return filter.read(key)>=thresh;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Converts k-mer to canonical key representation.
	 * Uses reverse complement if enabled to get canonical form.
	 * @param kmer K-mer to convert
	 * @return Canonical key for hash table lookup
	 */
	public long toKey(final long kmer){
		return (rcomp ? toKey(kmer, AminoAcid.reverseComplementBinaryFast(kmer, k)) : kmer);
	}
	
	/**
	 * Converts k-mer pair to canonical key.
	 * Uses maximum of forward/reverse if reverse complement enabled.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @return Canonical key (max of kmer/rkmer if rcomp, else kmer)
	 */
	public long toKey(final long kmer, final long rkmer){
		return (rcomp ? Tools.max(kmer, rkmer) : kmer);
	}
	
	/**
	 * Extracts k-mers from a read with quality filtering.
	 * Only includes k-mers meeting quality and probability thresholds.
	 *
	 * @param r Read to process
	 * @param list List to fill with k-mer keys
	 * @param k K-mer length
	 * @param minQuality Minimum base quality required
	 * @param minProb Minimum k-mer probability required
	 * @param rcomp Use reverse complement canonical form
	 */
	public static final void toKmers(Read r, final LongList list, int k, final int minQuality, final float minProb, final boolean rcomp){
		assert(k<=32);
		assert(list!=null);
		
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		
		if(bases==null || bases.length<k){return;}
		
		final int shift=bitsPerBase*k;
		final int shift2=shift-bitsPerBase;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		
		long kmer=0;
		long rkmer=0;
		int len=0;
		float prob=1;
		
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

			final byte q;
			if(quals==null){
				q=50;
			}else{
				q=quals[i];
				prob=prob*align2.QualityTools.PROB_CORRECT[q];
				if(len>=k){
					byte oldq=quals[i-k];
					prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
				}
			}

			if(x<0 || q<minQuality){
				len=0;
				kmer=rkmer=0;
				prob=1;
			}else{
				len++;
				if(len>=k && prob>=minProb){
					long key=(rcomp ? Tools.max(kmer, rkmer) : kmer);
					list.add(key);
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Maximum number of reads to process for filter construction */
	long maxReads=-1;
	/** Enable error correction during k-mer counting */
	boolean ecco=false;
	/** Merge overlapping read pairs before k-mer extraction */
	boolean merge=false;
	/** Minimum quality score for k-mer inclusion */
	byte minq=0;

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	/** Additional input file paths beyond primary and secondary */
	private ArrayList<String> extra=new ArrayList<String>();
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** The underlying k-mer count array for filtering operations */
	public final KCountArray7MTA filter;
	/** Small k-mer length for filtering */
	final int k;
	/** Large k-mer length for enhanced specificity */
	final int kbig;
	/** Number of small k-mers per big k-mer (kbig-k+1) */
	final int smallPerBig;
	/** Bits allocated per k-mer count in the array */
	final int bits;
	/** Number of hash functions used in the count array */
	final int hashes;
	/** Minimum consecutive k-mer matches required to trigger filtering */
	final int minConsecutiveMatches;//Note this is similar to smallPerBig
	
	/** Bit shift amount for k-mer encoding (bitsPerBase * k) */
	final int shift;
	/** Secondary bit shift for reverse complement operations */
	final int shift2;
	/** Bit mask for extracting k-mer values during rolling hash */
	final long mask;
	/** Whether to use reverse complement canonical k-mer representation */
	final boolean rcomp;

//	private final long usableMemory;
	/** Memory allocated for the k-mer count array in bytes */
	private final long filterMemory;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Override value for number of hash table cells (for testing) */
	public static long OVERRIDE_CELLS=-1;
	/** Number of bits used to encode each DNA base */
	static final int bitsPerBase=2;
	/** Whether to print memory usage information during initialization */
	public static boolean printMem=true;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private transient PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
}
