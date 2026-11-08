package bloom;

import java.io.File;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;
import stream.Read;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/**
 * Owns a counting Bloom table loaded from reads or an existing BBMap index.
 * Counts saturate at the table's cell limit and can be inflated by hash collisions;
 * a positive count is evidence of support, not proof of an exact reference match.
 *
 * Primitive sequence queries use packed small-k windows. The kbig helpers can
 * combine their depths by taking a minimum across constituent windows; this is
 * distinct from querying a multiword ukmer.Kmer when the stored k itself exceeds
 * thirty-one. Callers must use the query representation used during construction.
 *
 * Completed tables can be queried concurrently when callers supply separate
 * scratch lists. The public table remains mutable; this class does not isolate
 * construction or mutation from concurrent queries. Construction also uses
 * process-wide counting, encoding and I/O settings.
 *
 * @author Brian Bushnell
 * @date April 23, 2018
 *
 */
public class BloomFilter implements Serializable{

	/** Serialization identity retained for saved reference filters. */
	private static final long serialVersionUID=-3987955563503838492L;

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
	 * Parses command-line settings and loads the reference table before returning.
	 * Expands paired input templates and checks the primary inputs. Standard parser
	 * input paths are copied here; other parser-owned fields are not automatically
	 * forwarded to the reference loader.
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
		boolean dualHash_=false;

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Preserve '=' within values and the historical null value for empty options.
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=(split.length>1 && !split[1].isEmpty() ? split[1] : null);

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
			}else if(a.equals("dualhash") || a.equals("hash2")){
				dualHash_=Parse.parseBoolean(b);
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
		dualHash=dualHash_ && hashes>1;


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
	 * Allocates an empty table with the requested geometry, without independent
	 * secondary hashes. No reference reads are loaded by this overload.
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
		this(null, null, null, k_, kbig_, bits_, hashes_, minConsecutiveMatches_, rcomp_, ecco_, merge_, memFraction, false);
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
		this(in1_, in2_, extra_, k_, kbig_, bits_, hashes_, minConsecutiveMatches_, rcomp_, ecco_, merge_, memFraction, false);
	}

	/**
	 * Loads the supplied inputs, or creates an empty table if no primary or extra
	 * input is supplied. Copies extra paths into an owned list. The memory fraction
	 * scales a JVM-based estimate; a positive OVERRIDE_CELLS overrides that estimate
	 * for the table size. Dual hashing is effective only with more than one lane.
	 * @param in1_ Primary input path, or null
	 * @param in2_ Mate input path, or null
	 * @param extra_ Additional reference paths, or null
	 * @param k_ Stored k-mer length
	 * @param kbig_ Combined-window length, raised to k_ when smaller
	 * @param bits_ Counter width in bits
	 * @param hashes_ Number of table lanes
	 * @param minConsecutiveMatches_ Required run of supported windows for matching
	 * @param rcomp_ Select reverse-complement canonical keys when enabled
	 * @param ecco_ Correct overlapping mates during reference counting
	 * @param merge_ Merge overlapping mates during reference counting
	 * @param memFraction Multiplier on the estimated available table memory
	 * @param dualHash_ Request independent secondary hashes for lane selection
	 */
	public BloomFilter(String in1_, String in2_, ArrayList<String> extra_, int k_, int kbig_, int bits_, int hashes_,
			int minConsecutiveMatches_, boolean rcomp_, boolean ecco_, boolean merge_, float memFraction, boolean dualHash_){
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
		dualHash=dualHash_ && hashes>1;
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
	public BloomFilter(boolean bbmapIndex_, int k_, int kbig_, int bits_, int hashes_, int minConsecutiveMatches_, boolean rcomp_){
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
		dualHash=false;

		shift=bitsPerBase*k;
		shift2=shift-bitsPerBase;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		filter=loadFromIndex();
	}

	/**
	 * Estimates a table-memory budget from maximum heap, current used heap and
	 * the initial-to-maximum heap ratio, then applies the requested multiplier.
	 * This heuristic does not measure the minimum heap needed by the workload.
	 * @param mult Fraction of available memory to allocate
	 * @return Memory in bytes available for filter
	 */
	private static long setMemory(double mult){
		if(printMem){Shared.printMemory();}

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
	 * Loads reference counts through ReadCounter, or allocates an empty table.
	 * Positive OVERRIDE_CELLS takes precedence over the memory-derived cell count.
	 * @return Populated KCountArray7MTA with k-mer counts
	 */
	private KCountArray7MTA load(){
		final int cbits=bits;
		final long totalBits=8*filterMemory;
		final long cells=(OVERRIDE_CELLS>0 ? OVERRIDE_CELLS : totalBits/cbits);
		//		System.err.println("filterMemory="+filterMemory+", cells="+cells);
		KCountArray7MTA kca;
		if(in1!=null || (extra!=null && extra.size()>0)){
			ReadCounter rc=new ReadCounter(k, rcomp, ecco, merge, Shared.AMINO_IN, dualHash);
			kca=(KCountArray7MTA)rc.makeKca(in1, in2, extra==null || extra.isEmpty() ? null : extra,
					cbits, cells, hashes, minq,
					maxReads, 1, 1, 1, null, 0);
		}else{
			kca=(KCountArray7MTA) KCountArray.makeNew(cells, cbits, hashes, null, 0);
		}
		return kca;
	}

	/**
	 * Loads k-mer data from an existing BBMap index.
	 * Creates a KCountArray from pre-built index rather than raw sequences.
	 * Uses memory-derived capacity; OVERRIDE_CELLS does not apply to this path.
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
	public boolean passes(Read r1, Read r2, final int thresh){
		boolean pass=passes(r1, thresh);
		return pass && passes(r2, thresh);
	}

	/**
	 * Averages smoothed depths over fully defined windows, for stored k at most 31.
	 * Clips interior collision spikes to the maximum of the neighboring depths and
	 * each endpoint to its sole neighbor. A singleton retains its original depth.
	 * Smoothing state spans undefined runs; this is not a moving average.
	 * @param bases Sequence bases to analyze
	 * @return Smoothed depth average, or zero for null input or no defined windows
	 */
	public float averageCount(final byte[] bases){
		assert(k<=31) : KillSwitch.assertDie("BloomFilter.averageCount is k<=31 only (primitive-long rolling); k="+k+" needs a ukmer.Kmer path.");
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
			//Only rkmer is zeroed here (kmer is not, and several sibling loops zero neither).
			//Safe: the len<k gate blocks all use of kmer/rkmer until k more valid bases have
			//shifted any garbage past the 2k-bit mask, so both self-heal before len>=k again.
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
		//The delayed final depth needs the same one-neighbor clipping as the first.
		sum+=(counted>1 ? Tools.min(count, prev) : count);
		return sum/Tools.max(counted, 1f);
	}

	/**
	 * Finds the minimum k-mer count in a read.
	 * Ignores windows containing undefined bases and stops early at depth zero.
	 * @param r Read to analyze
	 * @return Minimum defined-window depth, or -1 for null input or no defined windows
	 */
	public int minCount(Read r){
		assert(k<=31) : KillSwitch.assertDie("BloomFilter.minCount is k<=31 only (primitive-long rolling); k="+k+" needs a ukmer.Kmer path.");
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
	 * Requires at least round(fraction*definedWindows) high-count windows.
	 * Undefined windows are excluded from both counts. An input without a fully
	 * defined window returns false even when the requested fraction is zero.
	 *
	 * @param r Read to evaluate
	 * @param thresh Count threshold for "high" k-mers
	 * @param fraction Minimum fraction of defined k-mers that must be high, from zero to one
	 * @return true if a defined window exists and the rounded required count is met
	 */
	public boolean hasHighCountFraction(Read r, final int thresh, final float fraction){
		if(r==null || r.length()<k){return false;}
		//k>31 primitive-long rolling would silently truncate the kmer to its last 32
		//bases and query a key the filter was never built with (the filter is built via
		//ukmer.Kmer.xor(); see ReadCounter.addReadBig). Roll a real ukmer.Kmer instead.
		if(k>31){return hasHighCountFractionBig(r, thresh, fraction);}
		final byte[] bases=r.bases;
		final int kmers=r.length()-k+1;
		//Large float-rounded budgets need the exact final count, without early rejection.
		if(kmers>0x1000000 || !(fraction>=0 && fraction<=1)){
			final int maxLow=(kmers>0x1000000 ? Integer.MAX_VALUE : kmers-requiredHighCount(fraction, kmers));
			return finishHighCountFraction(bases, thresh, fraction, maxLow, 0, 0, 0);
		}
		final int maxLow=kmers-requiredHighCount(fraction, kmers);
		long kmer=0;
		long rkmer=0;
		int low=0;
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			if(x<0){
				//The prefix has i defined bases. Continue after this separator without rereading it.
				return finishHighCountFraction(bases, thresh, fraction, maxLow, i+1, low, Tools.max(0, i-k+1));
			}
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(i>=k-1 && getCount(kmer, rkmer)<thresh){
				if(++low>maxLow){return false;}
			}
		}
		return true;
	}

	/** Continues after an undefined base, or scans from zero when early bounds are unsafe.
	 * Prefix counts are carried in; no base or table query from the prefix is repeated. */
	private boolean finishHighCountFraction(final byte[] bases, final int thresh, final float fraction,
			final int maxLow, final int start, int low, int counted){
		assert(start>=0 && start<=bases.length && counted>=low && low>=0) :
				"Fraction continuation needs a valid separator boundary and prefix counts: start="+start+
				" counted="+counted+" low="+low;
		long kmer=0, rkmer=0;
		int len=0;
		for(int i=start; i<bases.length; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			if(x<0){len=0; kmer=rkmer=0; continue;}
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(++len>=k){
				counted++;
				if(getCount(kmer, rkmer)<thresh && ++low>maxLow){return false;}
			}
		}
		return counted>0 && low<=counted-requiredHighCount(fraction, counted);
	}

	/** Rounds the float product exactly as {@link Math#round(float)}.
	 * The product must stay float; the double half-offset avoids premature float rounding.
	 * Negative values and NaNs retain the standard-library path. */
	private static int requiredHighCount(final float fraction, final int windows){
		final float product=fraction*windows;
		if(product>=0){return (int)(product+0.5d);}
		return Math.round(product);
	}

	/**
	 * Returns fraction of k-mers with count below threshold.
	 * Complement of highCountFraction, including a result of one when there are
	 * no defined windows or the read is null.
	 *
	 * @param r Read to examine
	 * @param thresh Count threshold
	 * @param smooth Apply count smoothing
	 * @return Fraction of low-count k-mers
	 */
	public float lowCountFraction(final Read r, final int thresh, final boolean smooth){
		return 1-highCountFraction(r, thresh, smooth);
	}

	/**
	 * Return the fraction of kmers with count at least thresh.
	 * @param r Read to examine.
	 * @param thresh Minimum count to be 'high'.
	 * @param smooth Require the previous defined window to be high too; the first has no such restriction
	 * @return High count fraction, or zero if no valid kmers.
	 */
	public float highCountFraction(final Read r, final int thresh, final boolean smooth){
		return r==null ? 0 : highCountFraction(r.bases, thresh, smooth);
	}

	/**
	 * Returns the fraction of fully defined windows with depth at least thresh.
	 * With smoothing, a high window following a low defined window contributes
	 * zero; the first defined window is counted normally. This previous-window
	 * state is retained across undefined-base runs.
	 * @param bases Bases to examine.
	 * @param thresh Minimum count to be 'high'.
	 * @param smooth Require both current and previous defined windows to be high
	 * @return High count fraction, or zero if no valid kmers.
	 */
	public float highCountFraction(final byte[] bases, final int thresh, boolean smooth){
		assert(k<=31) : KillSwitch.assertDie("BloomFilter.highCountFraction is k<=31 only (primitive-long rolling); k="+k+" needs a ukmer.Kmer path.");
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
	 * Applies a weak-end-support heuristic to a pair. Both left ends must have
	 * minimum depth at most one, and either right end below three, to reject a
	 * full-length pair. If a mate is absent or short, uses the single-read rule.
	 * This depth heuristic does not establish contamination provenance.
	 *
	 * @param r1 First read
	 * @param r2 Second read
	 * @param range Maximum number of window positions to inspect at each end; positive
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
	 * Rejects a short read, or a read whose left and right minimum depths are both
	 * below two. An end with no defined window has depth -1 for this comparison.
	 *
	 * @param r Read to evaluate
	 * @param range Maximum number of window positions to inspect at each end; positive
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
	 * @param range Maximum window positions, spanning at most k+range-1 bases
	 * @return Minimum defined-window depth in the left region, or -1 if none
	 */
	private int getLeftCount(byte[] bases, int range){
		assert(range>0) : range;
		if(bases.length<k){return -1;}
		if(k>31){return getLeftCountBig(bases, range);}//See hasHighCountFraction for why.
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
	 * @param range Maximum window positions, spanning at most k+range-1 bases
	 * @return Minimum defined-window depth in the right region, or -1 if none
	 */
	private int getRightCount(byte[] bases, int range){
		assert(range>0) : range;
		if(bases.length<k){return -1;}
		if(k>31){return getRightCountBig(bases, range);}//See hasHighCountFraction for why.
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

	/*--------------------------------------------------------------*/
	/*----------------      k>31 depth queries      ----------------*/
	/*--------------------------------------------------------------*/
	//These mirror getLeftCount/getRightCount/hasHighCountFraction, but roll a real
	//ukmer.Kmer and look up filter.read(kmer.xor()) -- the exact key the filter was
	//built with (ReadCounter.addReadBig) -- instead of a bit-packed long, which cannot
	//represent a k>31 kmer and would silently query the last-32-base truncation. Used
	//only when this BloomFilter's k>31 (the BBCMS k>31 path, where ksmall==k); the
	//small-k/kbig min-combine mechanism keeps its own k<=31 primitive path unchanged.
	//A fresh Kmer is allocated per call to keep query state private without introducing
	//a shared mutable cursor or transient ThreadLocal state. Its allocation cost must
	//be measured on the relevant workload rather than assumed negligible. new Kmer(k) picks
	//up the same process-global PACKED/FULL_MIX/MASK_CORE statics the build used.

	/** k>31 counterpart of {@link #getLeftCount}. */
	private int getLeftCountBig(byte[] bases, int range){
		final Kmer kmer=new Kmer(k);
		final int stop=Tools.min(bases.length, k+range-1);
		int min=Integer.MAX_VALUE;
		int counted=0;
		for(int i=0; i<stop; i++){
			kmer.addRight(bases[i]);
			if(kmer.len>=k){
				counted++;
				int count=readKmer(kmer);
				min=Tools.min(min, count);
			}
		}
		return counted>0 ? min : -1;
	}

	/** k>31 counterpart of {@link #getRightCount}. */
	private int getRightCountBig(byte[] bases, int range){
		final Kmer kmer=new Kmer(k);
		final int start=Tools.max(0, bases.length-k-range+1);
		int min=Integer.MAX_VALUE;
		int counted=0;
		for(int i=start; i<bases.length; i++){
			kmer.addRight(bases[i]);
			if(kmer.len>=k){
				counted++;
				int count=readKmer(kmer);
				min=Tools.min(min, count);
			}
		}
		return counted>0 ? min : -1;
	}

	/** k>31 counterpart of {@link #hasHighCountFraction}. */
	private boolean hasHighCountFractionBig(Read r, final int thresh, final float fraction){
		final byte[] bases=r.bases;
		final int kmers=r.length()-k+1;
		final int minHigh=Math.round(fraction*kmers);
		//Same conservative early-rejection bound as the primitive-window path.
		final int maxLow=(kmers<=0x1000000 ? kmers-minHigh : Integer.MAX_VALUE);
		final Kmer kmer=new Kmer(k);
		int low=0, counted=0;
		for(int i=0; i<bases.length; i++){
			kmer.addRight(bases[i]);
			if(kmer.len>=k){
				counted++;
				int count=readKmer(kmer);
				if(count<thresh){
					low++;
					if(low>maxLow){return false;}
				}
			}
		}
		return counted>0 && low<=counted-Math.round(fraction*counted);
	}

	/**
	 * Determines if a read passes the consecutive k-mer match filter.
	 * Rejects reads with too many consecutive high-count k-mers.
	 *
	 * @param r Read to evaluate
	 * @param thresh Count threshold for k-mer matches
	 * @return true if no required run is found, including null or too-short reads; false for a match
	 */
	public boolean passes(Read r, final int thresh){
		assert(k<=31) : KillSwitch.assertDie("BloomFilter.passes is k<=31 only (primitive-long rolling); k="+k+" needs a ukmer.Kmer path.");
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
			if(x<0){len=0; rkmer=0; streak=0;}else{len++;}
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
	 * Tests whether a read has the required consecutive supported windows.
	 * Inverts passes(r, keys, thresh), which recomputes keys into the scratch list.
	 *
	 * @param r Read to test
	 * @param keys Reusable scratch list; contents may be replaced
	 * @param thresh Count threshold
	 * @return true if read matches (should be filtered)
	 */
	public boolean matches(Read r, LongList keys, final int thresh){
		return !passes(r, keys, thresh);
	}

	/**
	 * Tests if either read in a pair matches filter criteria.
	 * Returns true if at least one read matches.
	 *
	 * @param r1 First read
	 * @param r2 Second read
	 * @param keys Caller-owned scratch list, reused for each mate actually queried
	 * @param thresh Count threshold
	 * @return true if either read matches filter
	 */
	public boolean matchesEither(Read r1, Read r2, LongList keys, final int thresh){
		boolean match=!passes(r1, keys, thresh);
		return match || !passes(r2, keys, thresh);
	}

	/**
	 * Tests if both reads in a pair pass, recomputing their keys in shared scratch.
	 * Both reads must pass for the pair to pass.
	 *
	 * @param r1 First read
	 * @param r2 Second read
	 * @param keys Caller-owned scratch list, reused for each mate actually queried
	 * @param thresh Count threshold
	 * @return true if both reads pass filter
	 */
	public boolean passes(Read r1, Read r2, LongList keys, final int thresh){
		boolean pass=passes(r1, keys, thresh);
		return pass && passes(r2, keys, thresh);
	}

	/**
	 * Tests a read using a reusable list of encoded k-mer keys.
	 * Recomputes valid keys and run separators in caller-owned scratch storage.
	 * A minimum streak below two delegates to the direct query without changing scratch.
	 * Null and too-short reads return true without clearing existing scratch contents.
	 *
	 * @param r Read to test
	 * @param keys Non-null scratch when a sufficiently long read needs a streak of two or more
	 * @param thresh Count threshold
	 * @return true if read passes consecutive match filter
	 */
	public boolean passes(Read r, LongList keys, final int thresh){
		assert(k<=31) : KillSwitch.assertDie("BloomFilter.passes is k<=31 only (primitive-long rolling); k="+k+" needs a ukmer.Kmer path.");
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
	 * Negative keys break runs. Requires minConsecutiveMatches greater than one;
	 * the input list is not modified.
	 *
	 * @param keys List of k-mer keys to evaluate
	 * @param thresh Count threshold for matches
	 * @return true if keys pass consecutive match criteria
	 */
	public boolean passes(final LongList keys, final int thresh){
		assert(minConsecutiveMatches>1);
		final long[] array=keys.array;
		final int len=keys.size;

		//Stride by minConsecutiveMatches, probing only positions i≡(mcm-1) mod mcm.
		//Any run of mcm consecutive matches spans mcm consecutive indices, which contain
		//exactly one such position, so this stride can never skip a real run.
		for(int i=minConsecutiveMatches-1; i<len; i+=minConsecutiveMatches){
			final boolean found;
			{
				final long key=array[i];
				found=(key<0 ? false : readKey(key)>=thresh);
			}
			if(found){
				int streak=1;
				//array[i-j] needs no bounds check: i>=mcm-1 (loop start) and j<mcm, so i-j>=0.
				for(int j=1; j<minConsecutiveMatches; j++){
					final long key=array[i-j];
					if(key<0 || readKey(key)<thresh){break;}
					else{streak++;}
				}
				if(streak>=minConsecutiveMatches){return false;}
				for(int j=1; j<minConsecutiveMatches && j+i<len; j++){
					final long key=array[i+j];
					if(key<0 || readKey(key)<thresh){break;}
					else{streak++;}
					if(streak>=minConsecutiveMatches){return false;}
				}
			}
		}
		return true;
	}

	/*--------------------------------------------------------------*/

	/** Replaces counts with one entry per small-k window; undefined windows have zero.
	 * @param bases Non-null sequence, unchanged by this method
	 * @param counts Reusable non-null destination, cleared even for short sequences
	 * @return Number of fully defined windows, not the number of output entries */
	public int fillCounts(byte[] bases, IntList counts){
		assert(k<=31) : KillSwitch.assertDie("BloomFilter.fillCounts is k<=31 only (primitive-long rolling); k="+k+" needs a ukmer.Kmer path. (k>31 counts go through fillCountsFromKmers in the corrector.)");
		counts.clear();
		final int blen=bases.length;
		if(blen<k){return 0;}
		final int min=k-1;
		long kmer=0, rkmer=0;
		int len=0;
		int valid=0;

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
		return readKey(key);
	}

	/**
	 * Looks up an already selected key without canonicalizing it again.
	 * @param key Canonical key when rcomp is enabled, forward key otherwise
	 * @return Count value from filter
	 */
	public int getCount(final long key){
//		assert(key==toKey(key, AminoAcid.reverseComplementBinaryFast(key, k))); //slow
		return readKey(key);
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
		return readKey(key)>=thresh;
	}

	/*--------------------------------------------------------------*/

	/** Replaces counts with the minimum small-window count for each kbig window.
	 * Undefined windows produce zero; sequences shorter than kbig clear the list.
	 * @param bases Non-null sequence, unchanged by this method
	 * @param counts Reusable non-null destination
	 * @return Number of fully defined kbig windows, which may be less than counts.size */
	public int fillCountsBig(byte[] bases, IntList counts){
		assert(smallPerBig>1) : smallPerBig;
		if(bases.length<kbig){counts.clear(); return 0;}
		final int valid0=fillCounts(bases, counts);
		final int bigWindows=bases.length-kbig+1;
		final boolean allDefined=(valid0==counts.size);
		if(valid0<smallPerBig){
			//No full big window can exist, but the output still has one slot per position.
			Arrays.fill(counts.array, 0, bigWindows, 0);
			counts.size=bigWindows;
			return 0;
		}
		for(int i=0; i<bigWindows; i++){
			int count=smallToBig(counts, i);
			counts.set(i, count);
		}
		counts.size=bigWindows;
		if(allDefined){return bigWindows;}
		//Only ambiguous input needs this scan; total valid small windows cannot
		//identify how many big windows fit wholly inside separate defined runs.
		int valid=0, len=0;
		for(byte base : bases){
			if(AminoAcid.baseToNumber[base]<0){len=0;}
			else if(++len>=kbig){valid++;}
		}
		return valid;
	}


	/**
	 * Replaces counts with the minimum constituent small-window depth for each
	 * packed kbig-mer. This overload consumes packed bases, not ukmer hash values;
	 * it does not recognize negative values as undefined-window sentinels.
	 * @param kmers Non-null list of fully defined, packed kbig-mers; unchanged
	 * @param counts Non-null reusable destination, cleared before appending results
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
	 * Takes the minimum of smallPerBig adjacent counts, stopping early at zero.
	 * The caller supplies a complete constituent window within the logical list.
	 *
	 * @param counts List containing nonnegative small-window depths
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
	 * Queries the constituent windows of the forward packed kbig-mer.
	 * The supplied reverse complement is unused: each small window is separately
	 * canonicalized according to this filter's strand setting.
	 *
	 * @param kmer Forward packed kbig-mer
	 * @param rkmer Reverse complement k-mer (unused)
	 * @return Big k-mer count
	 */
	@SuppressWarnings("unused")
	public int getCountBig(final long kmer, final long rkmer){
		return getCountBig(kmer);
	}

	/**
	 * Returns the minimum depth of the smallPerBig constituent windows in a
	 * packed kbig-mer, canonicalizing each small window when rcomp is enabled.
	 * Hash collisions can inflate these depths; the minimum is not an exact
	 * observed count of the complete kbig-mer.
	 * @param kmer Fully defined kbig-mer encoded in one long, not a ukmer hash
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
	 * Tests the selected key's raw table depth, exactly as contains does.
	 * Despite its name, this method does not split a packed big k-mer into small
	 * windows or call getCountBig. Callers needing that minimum must query it directly.
	 *
	 * @param kmer Forward key
	 * @param rkmer Reverse-complement key
	 * @param thresh Count threshold
	 * @return true if the selected key's table count is at least thresh
	 */
	public boolean containsBig(final long kmer, final long rkmer, final int thresh){
		final long key=toKey(kmer, rkmer);
		return readKey(key)>=thresh;
	}

	/** Reads an already selected key, deriving a secondary lane hash only in dual-hash mode. */
	private int readKey(final long key){
		return dualHash ? filter.read(key, KCountArray.secondaryHash(key)) : filter.read(key);
	}

	/** Reads the supplied Kmer's primary and, when enabled, independent secondary hash. */
	private int readKmer(final Kmer kmer){
		return dualHash ? filter.read(kmer.xor(), kmer.xor2()) : filter.read(kmer.xor());
	}

	/**
	 * Gets the Bloom-filter depth for a ukmer.Kmer using this filter's exact
	 * configured hash mode.
	 * @param kmer K-mer to query
	 * @return Raw Bloom-filter depth
	 */
	public int getCount(final Kmer kmer){
		return readKmer(kmer);
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
	 * Appends qualifying packed windows without clearing the destination.
	 * Undefined or low-quality bases break the current run. With qualities present,
	 * the probability threshold applies to the product of per-base correctness
	 * probabilities in each window. Missing qualities use quality fifty for the
	 * minimum-quality check and probability one for the probability check.
	 *
	 * @param r Non-null read; null bases or fewer than k bases append nothing
	 * @param list Non-null destination; existing entries are retained
	 * @param k Packed DNA window length, from one through thirty-two
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

	/** Input-fragment limit passed to ReadCounter during construction; negative means unlimited. */
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
	/** Stored window length; primitive query methods require at most thirty-one bases. */
	final int k;
	/** Combined-window length, at least k; constituent-minimum queries use this geometry. */
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
	/** Whether kmers use independent primary and secondary Bloom-lane hashes. */
	public final boolean dualHash;

//	private final long usableMemory;
	/** Estimated table-memory budget in bytes; OVERRIDE_CELLS can supersede it during read loading. */
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
	/** Retained status field; this class does not currently assign it after initialization. */
	public boolean errorState=false;
}
