package kmer;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicLong;

import bloom.KCountArray;
import bloom.KmerCountAbstract;
import bloom.ReadCounter;
import fileIO.ByteStreamWriter;
import jgi.CallPeaks;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/**
 * Loads and holds kmers for Tadpole
 * @author Brian Bushnell
 * @date Jun 22, 2015
 *
 */
public abstract class AbstractKmerTableSet {
	
	/**
	 * Checks if a command-line argument is valid for k-mer table configuration.
	 * Validates arguments related to input files, table sizing, threading, filtering,
	 * and processing parameters.
	 *
	 * @param a The argument string to validate
	 * @return true if the argument is recognized and valid, false otherwise
	 */
	public static final boolean isValidArgument(String a){
			if(a.equals("in") || a.equals("in1")){
			}else if(a.equals("in2")){
			}else if(a.equals("append") || a.equals("app")){
			}else if(a.equals("overwrite") || a.equals("ow")){
			}else if(a.equals("initialsize")){
			}else if(a.equals("showstats") || a.equals("stats")){
			}else if(a.equals("ways")){
			}else if(a.equals("buflen") || a.equals("bufflen") || a.equals("bufferlength")){
			}else if(a.equals("k")){
			}else if(a.equals("threads") || a.equals("t")){
			}else if(a.equals("showspeed") || a.equals("ss")){
			}else if(a.equals("ecco")){
			}else if(a.equals("merge")){
			}else if(a.equals("verbose")){
			}else if(a.equals("verbose2")){
			}else if(a.equals("minprob")){
			}else if(a.equals("reads") || a.startsWith("maxreads")){
			}else if(a.equals("prealloc") || a.equals("preallocate")){
			}else if(a.equals("prefilter")){
			}else if(a.equals("prefiltersize") || a.equals("prefilterfraction") || a.equals("pff")){
			}else if(a.equals("minprobprefilter") || a.equals("mpp")){
			}else if(a.equals("minprobmain") || a.equals("mpm")){
			}else if(a.equals("prefilterpasses") || a.equals("prepasses")){
			}else if(a.equals("prehashes") || a.equals("hashes")){
			}else if(a.equals("onepass")){
			}else if(a.equals("passes")){
			}else if(a.equals("rcomp")){
			}else if(a.equals("maskmiddle")){
			}else if(a.equals("filtermemory") || a.equals("prefiltermemory") || a.equals("filtermem")){
			}else{
				return false;
			}
			return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Main processing method that executes the complete k-mer loading pipeline.
	 * Counts k-mers from input files, tracks timing statistics, and handles error states.
	 * @param t Timer for tracking execution time
	 * @throws RuntimeException if processing encounters errors
	 */
	public final void process(Timer t){
		
		/* Count kmers */
		long added=processInput();
		
		/* Stop timer and calculate speed statistics */
		t.stop();
		
		showStats(t, added);
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getSimpleName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	
	/** Clears all k-mer data from the table set, freeing memory for reuse */
	public abstract void clear();
	
	
	/**
	 * Processes all input files to build k-mer tables with optional prefiltering.
	 * Creates prefilter arrays for memory efficiency, allocates main tables,
	 * and loads k-mers from all specified input sources.
	 * @return Number of k-mers successfully added to tables
	 */
	public final long processInput(){
		
		/* Start phase timer */
		Timer t=new Timer();

//		if(DISPLAY_PROGRESS){
//			outstream.println("Before loading:");
//			Shared.printMemory();
//			outstream.println();
//		}
		
		prefilterArray=makePrefilter(new KCountArray[1], null);
		if(prefilterArray!=null){
			prefilterArray.purgeFilter();
			filterMax2=Tools.min(filterMax, prefilterArray.maxValue-1);
			
			/* This is already getting printed in makePrefilter */
//			if(DISPLAY_PROGRESS){
//				outstream.println("After prefilter:");
//				Shared.printMemory();
//				outstream.println();
//			}
		}
//		assert(false) : prefilterArray.cellBits+", "+prefilterArray.maxValue+", "+filterMax+", "+filterMax2;
		
		if(DISPLAY_STATS){System.err.println("Estimated kmer capacity: \t"+estimatedKmerCapacity());}
		
		assert(!allocated);
		allocateTables();
		allocated=true;
		
		if(DISPLAY_PROGRESS){
			outstream.println("After table allocation:");
			Shared.printMemory();
			outstream.println();
		}
		
		/* Fill tables with kmers */
		long added=loadKmers();
		
		/* Clear prefilter; no longer needed */
		prefilterArray=null;
		
//		long removed=0;
//		if(prefilter && filterMax>0){
//			removed=removeKmersAtMost(filterMax);
//			System.err.println("Removed "+removed+" low-depth kmers.");
//		}
		
		return added;
	}
	
	
	/**
	 * Creates a count-min sketch prefilter for memory-efficient k-mer processing.
	 * Uses multiple passes to estimate k-mer frequencies before main table allocation.
	 * Automatically adjusts cell bits and hash count based on available memory and target accuracy.
	 *
	 * @param filter Array to store the created prefilter (modified in place)
	 * @param ht Timer for tracking prefilter creation time (may be null)
	 * @return The created KCountArray prefilter, or null if prefiltering is disabled
	 */
	public final KCountArray makePrefilter(final KCountArray[] filter, Timer ht){
//		assert(false) : lastFilter+", "+prefilter+", "+filterMax()+", "+currentPass+", "+filterMemory(currentPass);
		if(!prefilter){return null;}
		assert(!MASK_MIDDLE) : "TODO: MaskMiddle not yet supported by prefilter.";
		
		if(filter[0]!=null){
			filter[0].purgeFilter();
			assert(filter[0].prefilter()==null);
		}
		
		KmerCountAbstract.CANONICAL=true;

		long precells=-1;
		int cbits=1;
		if(onePass){
			while(filterMax>=(1<<cbits)){cbits*=2;}
		}else{
			while(filterMax+1>=(1<<cbits)){cbits*=2;}
		}
		if(prepasses>2 && currentPass==prepasses-1){cbits=1;}
		
		byte minq=0;
		if(precells<1){
			long prebits=(filterMemory(currentPass)-10)*8;
			
//			System.err.println("prebits="+prebits+", currentPass="+currentPass);
			
			precells=prebits/cbits;
			if(precells<100000){ //Not enough memory - no point.
				prefilter=false;
				return null;
			}
		}
		if(prehashes<1){prehashes=2;}

		ReadCounter rc=new ReadCounter(kbig(), true, ecco(), false, Shared.AMINO_IN);
		if(onePass){
			assert(filter==null || filter.length==1) : "Multiple filtering passes are not allowed in onepass mode.\n"+filter.length+","+prepasses+", "+onePass+", "+prefilter;
			filter[0]=rc.makeKca(null, null, null, cbits, precells, prehashes, minq,
					maxReads, 1, 1, 1, null, 0);
		}else{
			if(ht==null){ht=new Timer();}
			ht.start();
			filter[0]=rc.makeKca_als(in1, in2, extra, cbits, precells, prehashes, minq,
					maxReads, 1, 1, 1, filter[0], filterMax);
			assert(filterMax<filter[0].maxValue || (currentPass>0 && currentPass==prepasses-1));
			outstream.println("Made prefilter:   \t"+filter[0].toShortString(prehashes));
			double uf=filter[0].usedFraction();
//			System.err.println("cellsUsed: "+filter[0].cellsUsed(1)+" //123"); //123
			if(uf>0.5){
				outstream.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.95 ? "incredibly" : uf>0.9 ? "extremely" : uf>0.8 ? "very" :
					uf>0.7 ? "rather" : uf>0.6 ? "fairly" : "somewhat")+" full.  Ideal load is under 50% used." +
						"\nFor better accuracy, run on a node with more memory; quality-trim or error-correct reads; or increase prefiltersize.");
			}
			ht.stop();
			currentPass++;
			
			final double kmers=filter[0].estimateUniqueKmers(prehashes, Tools.min(filterMax+1, filter[0].maxValue));
			outstream.println("Estimated valid kmers: \t\t"+(long)kmers);
			
//			outstream.println("Estimated valid kmers 1+: "+(long)filter[0].estimateUniqueKmers(prehashes, 1));
//			outstream.println("Estimated valid kmers 2+: "+(long)filter[0].estimateUniqueKmers(prehashes, 2));
//			outstream.println("Estimated valid kmers 3+: "+(long)filter[0].estimateUniqueKmers(prehashes, 3));
//			outstream.println("Estimated valid kmers 4+: "+(long)filter[0].estimateUniqueKmers(prehashes, 4));
			
			if(prepasses<0){//auto
				if((currentPass&1)==0){
					return makePrefilter(filter, ht);
				}else if(currentPass<5){
					if(kmers>estimatedKmerCapacity()){
						return makePrefilter(filter, ht);
					}
				}
			}else if(currentPass<prepasses){
				return makePrefilter(filter, ht);
			}
			
			if(DISPLAY_PROGRESS){
				outstream.println("Prefilter time:\t"+ht);
				outstream.println("After prefilter:");
				Shared.printMemory();
				outstream.println();
			}
		}
		
		return filter[0];
	}
	
	
	/**
	 * Displays processing statistics including read counts, k-mer counts, and timing.
	 * Shows memory usage information and throughput statistics if enabled.
	 * @param t Timer containing elapsed processing time
	 * @param added Number of unique k-mers loaded into tables
	 */
	public final void showStats(Timer t, long added){
		
		if(!DISPLAY_STATS){return;}
		
		if(DISPLAY_PROGRESS){
			outstream.println("After loading:");
			Shared.printMemory();
			outstream.println();
		}
		
		t.stop();
		outstream.println("Input:                      \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		outstream.println("Unique Kmers:               \t"+added);
		outstream.println("Load Time:                  \t"+t);
		
		if(showSpeed){
			outstream.println();
			outstream.println(Tools.readsBasesProcessed(t.elapsed, readsIn, basesIn, 8));
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads k-mers from all configured input files including paired-end handling.
	 * Processes both primary input files (in1/in2) and extra reference files.
	 * Handles automatic paired-end file detection using '#' placeholder syntax.
	 * @return Total number of k-mers loaded across all files
	 */
	private final long loadKmers(){
		//allocateTables();
		assert(allocated);
		kmersLoaded=0;
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=false;
//		Timer t=new Timer();
		for(int i=0; i<in1.size(); i++){
			String a=in1.get(i);
			String b=in2.size()>i ? in2.get(i) : null;
			int idx=a.indexOf('#');
			if(idx>=0 && b==null && !new File(a).exists()){
				b=a.replaceFirst("#", "2");
				a=a.replaceFirst("#", "1");
			}
			kmersLoaded+=loadKmers(a, b);
		}
		for(int i=0; i<extra.size(); i++){
			String a=extra.get(i);
			String b=null;
			int idx=a.indexOf('#');
			if(idx>=0 && b==null && !new File(a).exists()){
				b=a.replaceFirst("#", "2");
				a=a.replaceFirst("#", "1");
			}
			kmersLoaded+=loadKmers(a, b);
		}
//		t.stop("loadKmers:");
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		return kmersLoaded;
	}
	
	/**
	 * Load reads into tables, using multiple LoadThread.
	 */
	public abstract long loadKmers(String fname1, String fname2);
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Regenerates k-mer counts up to specified depth limit.
	 * Used for error correction and quality improvement operations.
	 * @param limit Maximum count value to regenerate
	 * @return Number of k-mers regenerated
	 */
	public abstract long regenerate(int limit);
	
	/**
	 * Retrieves a specific table from the table set by index.
	 * @param tnum Table number/index to retrieve
	 * @return The table object at the specified index
	 */
	public abstract Object getTable(int tnum);
	
	/**
	 * Creates a histogram of k-mer count frequencies up to specified maximum.
	 * @param histMax Maximum count value to include in histogram
	 * @return Array where index represents count and value represents frequency
	 */
	public abstract long[] fillHistogram(int histMax);
	
	/**
	 * Creates and fills an array with GC content counts for k-mers.
	 * @param histMax Maximum count value to process
	 * @return Array of GC counts indexed by k-mer frequency
	 */
	public final long[] fillGcCounts(int histMax){
		long[] gcCounts=new long[histMax+1];
		countGC(gcCounts, histMax);
		return gcCounts;
	}
	
	/**
	 * Counts GC bases in k-mers up to specified maximum count.
	 * @param gcCounts Array to store GC count results (modified in place)
	 * @param max Maximum k-mer count to consider
	 */
	public final void countGC(long[] gcCounts, int max) {
		// Determine thread count, capped by the number of tables to avoid empty threads
		final int threads=Tools.mid(1, Shared.threads(), tables().length);
		
		// If single threaded, just run it directly to save setup overhead
		if(threads<2){
			for(KmerTableInterface set : tables()){
				set.countGC(gcCounts, max);
			}
			return;
		}

		// Setup threads
		ArrayList<CountGCThread> al=new ArrayList<CountGCThread>(threads);
		
		for(int i=0; i<threads; i++){
			CountGCThread t=new CountGCThread(i, threads, max);
			al.add(t);
			t.start();
		}
		
		// Wait for completion and reduce results
		for(CountGCThread t : al){
			while(t.getState()!=Thread.State.TERMINATED){
				try{
					t.join();
				}catch(InterruptedException e){
					// Ideally handle interrupt, but standard BBTools practice often ignores it here
				}
			}
			// Reduction step: Sum thread-local counts into the main array
			for(int i=0; i<gcCounts.length; i++){
				gcCounts[i]+=t.myGcCounts[i];
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	private class CountGCThread extends Thread {
		CountGCThread(int tid_, int threads_, int max_){
			tid=tid_;
			threads=threads_;
			// Each thread gets a private array to avoid atomic contention
			myGcCounts=new long[max_+1];
			max=max_;
		}
		
		@Override
		public void run(){
			KmerTableInterface[] tables=tables();
			for(int i=tid; i<tables.length; i+=threads){
				// This calls the existing single-threaded method in AbstractKmerTable
				tables()[i].countGC(myGcCounts, max);
			}
		}
		
		final int tid, threads;
		final int max;
		final long[] myGcCounts;
	}
	
	/**
	 * Creates a GC content histogram as fractions from count and GC data.
	 * Calculates GC percentage for each frequency bin based on k-mer length.
	 *
	 * @param counts K-mer frequency counts per bin
	 * @param gcCounts GC base counts per bin
	 * @return Array of GC fractions (0.0 to 1.0) indexed by frequency
	 */
	public final float[] makeGcHistogram(long[] counts, long[] gcCounts){
		float[] gcHist=new float[counts.length];
		final long k=kbig();
		for(int i=0; i<counts.length; i++){
			long gc=gcCounts[i];
			double bases=Tools.max(counts[i], 1)*k;
			gcHist[i]=(float)(gc/bases);
		}
		return gcHist;
	}
	
	/** Initializes ownership tracking for thread-safe k-mer table access */
	public abstract void initializeOwnership();
	
	/** Clears all ownership claims on k-mer table entries */
	public abstract void clearOwnership();
	
	/** Returns the number of ways (hash table segments) in the table set.
	 * @return Number of parallel table segments for load balancing */
	public abstract int ways();
	
	/**
	 * Fills count values for all k-mers in a sequence.
	 * Convenience method that calls fillSpecificCounts without position filtering.
	 *
	 * @param bases Sequence bases to process
	 * @param counts List to store k-mer counts (modified in place)
	 * @param kmer Reusable Kmer object for efficiency
	 * @return Number of k-mers processed
	 */
	public final int fillCounts(byte[] bases, IntList counts, Kmer kmer){
		return fillSpecificCounts(bases, counts, null, kmer);
	}
	
	/**
	 * Fills count values for k-mers at specific positions in a sequence.
	 *
	 * @param bases Sequence bases to process
	 * @param counts List to store k-mer counts (modified in place)
	 * @param positions BitSet indicating which positions to process (null for all)
	 * @param kmer Reusable Kmer object for efficiency
	 * @return Number of k-mers processed
	 */
	public abstract int fillSpecificCounts(byte[] bases, IntList counts, BitSet positions, Kmer kmer);
	
	/**
	 * Regenerates k-mer counts for a sequence, tracking which positions changed.
	 *
	 * @param bases Sequence bases to process
	 * @param counts List to store updated k-mer counts (modified in place)
	 * @param kmer Reusable Kmer object for efficiency
	 * @param changed BitSet tracking which positions were modified (modified in place)
	 * @return Number of k-mers regenerated
	 */
	public abstract int regenerateCounts(byte[] bases, IntList counts, Kmer kmer, BitSet changed);
	
	/*--------------------------------------------------------------*/
	/*----------------       Printing Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Dumps k-mers to file in binary format within specified count range.
	 * Single-threaded implementation.
	 *
	 * @param fname Output file path
	 * @param mincount Minimum k-mer count to include
	 * @param maxcount Maximum k-mer count to include
	 * @param printTime Whether to display timing information
	 * @param remaining Counter for remaining k-mers to process
	 * @return true if dump completed successfully, false on error
	 */
	public abstract boolean dumpKmersAsBytes(String fname, int mincount, int maxcount, boolean printTime, AtomicLong remaining);
	/**
	 * Dumps k-mers to file in binary format within specified count range.
	 * Multi-threaded implementation for improved performance on large datasets.
	 *
	 * @param fname Output file path
	 * @param mincount Minimum k-mer count to include
	 * @param maxcount Maximum k-mer count to include
	 * @param printTime Whether to display timing information
	 * @param remaining Counter for remaining k-mers to process
	 * @return true if dump completed successfully, false on error
	 */
	public abstract boolean dumpKmersAsBytes_MT(String fname, int mincount, int maxcount, boolean printTime, AtomicLong remaining);
	
	/**
	 * Creates and optionally writes a k-mer count histogram with various analysis options.
	 * Supports smoothing, log scaling, GC content analysis, and flexible output formatting.
	 *
	 * @param fname Output file path (null to skip file output)
	 * @param cols Number of columns in output format
	 * @param max Maximum count value to include in histogram
	 * @param printHeader Whether to include column headers in output
	 * @param printZeros Whether to include zero-count entries
	 * @param printTime Whether to display timing information
	 * @param smooth Whether to apply progressive smoothing to histogram
	 * @param calcGC Whether to calculate GC content statistics
	 * @param doLogScale Whether to apply logarithmic scaling
	 * @param logWidth Width parameter for log scaling
	 * @param logPasses Number of passes for log scaling
	 * @param smoothRadius Radius for smoothing operations
	 * @return 2D array containing [counts, gc_counts] histograms
	 */
	public final long[][] makeKhist(String fname, int cols, int max, boolean printHeader, boolean printZeros, boolean printTime, 
			boolean smooth, boolean calcGC, boolean doLogScale, double logWidth, int logPasses, int smoothRadius){
		Timer t=new Timer(), t2=new Timer();
		boolean printSteps=true;
		long[] ca=fillHistogram(max);
		if(printSteps) {t2.stopAndStart("fillHistogram:");}
		float[] gcHist=null;
		if(calcGC){
//			assert(false) : max+", "+ca.length;
			long[] gc=(calcGC ? fillGcCounts(max) : null);
			if(printSteps) {t2.stopAndStart("fillGcCounts:");}
//			assert(false) : max+", "+ca.length+", "+gc.length;
			gcHist=makeGcHistogram(ca, gc);
			if(printSteps) {t2.stopAndStart("makeGcHistogram:");}
		}
		
		long[] logScale=null;
		
		if(smooth){
			ca=CallPeaks.smoothProgressive(ca, smoothRadius);
			if(printSteps) {t2.stopAndStart("smoothProgressive:");}
		}
		if(doLogScale){
			logScale=CallPeaks.logScale(ca, logWidth, 1, logPasses);
			if(printSteps) {t2.stopAndStart("logScale:");}
		}
		
		long[][] ret=new long[2][];
		ret[0]=ca;
		if(gcHist!=null){
			final int k=kbig();
			ret[1]=new long[ca.length];
			for(int i=1; i<ca.length; i++){
				ret[1][i]=Math.round(ca[i]*gcHist[i]*k);
			}
			if(printSteps) {t2.stopAndStart("gcHist:");}
		}
		
		if(fname==null){return ret;}
		
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, false, true);
		bsw.start();
		if(printHeader){
			bsw.print("#Depth\t"+(cols==3 ? "RawCount\t" : "")+"Count"+(doLogScale ? "\tlogScale" : "")+(calcGC ? "\tGC%\n" : "\n"));
		}
		
		for(int i=1; i<ca.length; i++){
			long count=ca[i];
			if(printZeros || count>0){
				bsw.print(i);
				bsw.print('\t');
				if(cols==3){
					bsw.print(i*count);
					bsw.print('\t');
				}
				bsw.print(count);
				if(doLogScale){
					bsw.print('\t').print(logScale[i]);
				}
				if(gcHist!=null){
					bsw.print(Tools.format("\t%.2f", 100f*gcHist[i]));
				}
				bsw.print('\n');
			}
		}
		if(printSteps) {t2.stopAndStart("fileWrite:");}
		bsw.poisonAndWait();
		t.stop();
		if(printTime){outstream.println("Histogram Write Time:       \t"+t);}
		return ret;
	}
	
	public abstract <T extends KmerTableInterface> T[] tables();
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Whether to display processing statistics */
	public boolean showStats=true;
	
//	public boolean silent=false;
	
	/** Has this class encountered errors while processing? */
	public boolean errorState=false;
	
	/** Use a count-min prefilter for low-depth kmers */
	public boolean prefilter=false;
	/** Fill the prefilter at the same time as the main table */
	public boolean onePass=false;
	
	/** Whether to process amino acid sequences instead of nucleotides */
	public boolean amino=false;
	/** Number of hashes used by prefilter */
	public int prehashes=2;
	/** Fraction of memory used by prefilter */
	public double prefilterFraction=0.2;
	
	/** Initial size of data structures */
	public int initialSize=-1;
	/** Fraction of available memory preallocated to arrays */
	public double preallocFraction=1.0;
	
	/** The active prefilter array for k-mer frequency estimation */
	public KCountArray prefilterArray=null;
	
	/** Whether to use minimum probability filtering during prefilter phase */
	public boolean minProbPrefilter=true;
	/** Whether to use minimum probability filtering during main processing */
	public boolean minProbMain=true;

	/** Input reads for kmers */
	public ArrayList<String> in1=new ArrayList<String>(), in2=new ArrayList<String>();
	
	/** Extra files for use as kmers */
	public ArrayList<String> extra=new ArrayList<String>();
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	public long maxReads=-1;
	
	/** Buffer length for I/O operations */
	public int buflen=1000;
	
	/** Filter kmers up to this level; don't store them in primary data structure */
	protected int filterMax=0;
	/** Secondary filter maximum after prefilter adjustment */
	protected int filterMax2=0;
	
	/** Total number of reads processed */
	public long readsIn=0;
	/** Total number of bases processed */
	public long basesIn=0;
	/** Number of reads with low quality scores */
	public long lowqReads=0;
	/** Number of bases with low quality scores */
	public long lowqBases=0;
	/** Number of reads that were trimmed during processing */
	public long readsTrimmed=0;
	/** Number of bases removed during trimming operations */
	public long basesTrimmed=0;
	
	/** Total number of k-mers encountered in input */
	public long kmersIn=0;
	/** Number of k-mers successfully loaded into tables */
	public long kmersLoaded=0;
	
	/** Current preprocessing pass number for multi-pass filtering */
	private int currentPass=0;
	/** Total number of preprocessing passes to perform */
	protected int prepasses=1;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns the k-mer length used by this table set.
	 * @return Length of k-mers in bases */
	public abstract int kbig();
	/**
	 * Calculates memory allocated for filtering operations in specified pass.
	 * @param pass Pass number to calculate memory for
	 * @return Memory allocation in bytes
	 */
	public abstract long filterMemory(int pass);
	/** Returns total memory used by the main k-mer tables.
	 * @return Memory usage in bytes */
	public abstract long tableMemory();
	/** Estimates the maximum number of k-mers this table set can hold.
	 * @return Estimated k-mer capacity */
	public abstract long estimatedKmerCapacity();
	/** Returns whether error correction and counting overlap is enabled.
	 * @return true if ECCO mode is active, false otherwise */
	public abstract boolean ecco();
	/** Returns whether quality trimming from left end is enabled.
	 * @return true if left-end quality trimming is active */
	public abstract boolean qtrimLeft();
	/** Returns whether quality trimming from right end is enabled.
	 * @return true if right-end quality trimming is active */
	public abstract boolean qtrimRight();
	/** Returns the minimum average quality score required for sequences.
	 * @return Minimum average quality threshold */
	public abstract float minAvgQuality();
	/** Returns the maximum k-mer count for prefiltering */
	public final int filterMax(){return filterMax;}
	/** Returns whether reverse complement processing is enabled.
	 * @return true if reverse complement k-mers are processed */
	public abstract boolean rcomp();
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Allocates memory for the main k-mer tables based on current configuration */
	protected abstract void allocateTables();
	/** Whether the main k-mer tables have been allocated */
	protected boolean allocated=false;

	/** Print messages to this stream */
	public static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Print speed statistics upon completion */
	public static boolean showSpeed=true;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Display kmer loading information */
	public static boolean DISPLAY_STATS=true;
	/** Verbose messages */
	public static boolean verbose=false;
	/** Debugging verbose messages */
	public static boolean verbose2=false;
	
	/** Maximum number of N bases allowed in k-mers */
	public static int maxNs=Integer.MAX_VALUE;
	/** Minimum sequence length to process */
	public static int minLen=0;
	
	/** Increment owner by this much to indicate claim is final. */
	public static final int CLAIM_OFFSET=100000;
	
	/** Default initial table size */
	public static final int initialSizeDefault=128000;
	
	/** Probability lookup table for base quality scores */
	public static final float[] PROB_CORRECT=Arrays.copyOf(align2.QualityTools.PROB_CORRECT, 127);
	/** Inverse probability lookup table for base quality scores */
	public static final float[] PROB_CORRECT_INVERSE=Arrays.copyOf(align2.QualityTools.PROB_CORRECT_INVERSE, 127);
	
	/** Whether to ignore unrecognized command-line arguments */
	public static boolean IGNORE_UNKNOWN_ARGS=true;

	/** Constant indicating hash collision occurred during lookup */
	/** Constant indicating k-mer is not present in table */
	public static final int NOT_PRESENT=AbstractKmerTable.NOT_PRESENT, HASH_COLLISION=AbstractKmerTable.HASH_COLLISION;
	/** Constant indicating k-mer has no thread ownership */
	public static final int NO_OWNER=AbstractKmerTable.NO_OWNER;
	
	/** Default minimum probability threshold for k-mer inclusion */
	public static double defaultMinprob=0;
	
	/** IIRC this has to do with allowing 32-mers */
	public static boolean MASK_CORE=false;
	/** Whether to mask the middle region of k-mers */
	public static boolean MASK_MIDDLE=false;
	/** Whether to use fast filling algorithm for table population */
	public static boolean FAST_FILL=true;
	
}
