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
import stream.FASTQ;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/**
 * Shared loading, prefiltering and reporting lifecycle for short- and long-k-mer tables.
 * Concrete backends allocate and own the shards. Loading mutates instance state and
 * process-wide read settings; callers must not run conflicting configurations concurrently.
 * Traversal requires stable tables. Recorded worker and I/O failure limitations remain
 * explicit at the affected methods rather than implying every partial result is detected.
 * @author Brian Bushnell
 * @date Jun 22, 2015
 *
 */
public abstract class AbstractKmerTableSet{

	/**
	 * Recognizes argument names forwarded to table-set constructors.
	 * This is a name whitelist, not validation of values or supported combinations.
	 *
	 * @param a Lowercase argument name without its value; nonnull
	 * @return Whether the name is recognized, including names beginning with maxreads
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
		}else if(a.equals("hashonly") || a.equals("hashedkmers")){
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
	 * Loads inputs, stops the supplied timer, reports enabled statistics and checks errorState.
	 * A loader exception propagates before reporting; an errorState failure is thrown
	 * after reporting. Neither path rolls back tables or counters.
	 * @param t Running timer whose elapsed interval includes loading
	 * @throws RuntimeException If loading throws or errorState is set after reporting
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

	/** Releases backing-table references; does not promise to reset the allocation lifecycle. */
	public abstract void clear();

	/**
	 * Builds an optional prefilter, allocates main tables once and loads configured inputs.
	 * Drops the prefilter reference after normal loading. This entry point does not
	 * itself inspect errorState; process(Timer) checks it after reporting statistics.
	 * @return Sum of backend loader returns, normally the number of newly inserted keys
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

		if(DISPLAY_STATS && showStats){System.err.println("Estimated kmer capacity: \t"+estimatedKmerCapacity());}

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
	 * Allocates or fills the count-min prefilter using the configured memory and count threshold.
	 * One-pass mode allocates an empty sketch for the main loader to populate. Other
	 * modes read inputs now, chaining prior sketches and advancing currentPass.
	 * A negative prepasses enables the bounded automatic pass schedule; a positive
	 * value is the target total pass count. Fewer than 100000 available cells disables
	 * prefiltering and returns null. CANONICAL is set globally and is not restored here.
	 *
	 * @param filter Nonnull one-element holder; a prior sketch may be supplied in slot zero
	 * @param ht Optional timer restarted for each non-one-pass input pass, not an accumulated total
	 * @return Final sketch, or null when disabled initially or by insufficient planned memory
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
		//Use long arithmetic so selecting 32-bit cells cannot wrap the shift or threshold.
		if(onePass){
			while(filterMax>=(1L<<cbits)){cbits*=2;}
		}else{
			while(filterMax+1L>=(1L<<cbits)){cbits*=2;}
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
			//Exact saturation is supported: processInput caps its subsequent gate at maxValue-1.
			assert(filterMax<=filter[0].maxValue || (currentPass>0 && currentPass==prepasses-1)) :
				"Counter must represent the threshold except on the final one-bit pass: threshold="+filterMax+
				", maximum="+filter[0].maxValue+", pass="+currentPass+", passes="+prepasses;
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

			final int minDepth=Tools.min(filterMax, filter[0].maxValue-1)+1;
			final double kmers=filter[0].estimateUniqueKmers(prehashes, minDepth);
			outstream.println("Estimated valid kmers: \t\t"+(long)kmers);

//			outstream.println("Estimated valid kmers 1+: "+(long)filter[0].estimateUniqueKmers(prehashes, 1));
//			outstream.println("Estimated valid kmers 2+: "+(long)filter[0].estimateUniqueKmers(prehashes, 2));
//			outstream.println("Estimated valid kmers 3+: "+(long)filter[0].estimateUniqueKmers(prehashes, 3));
//			outstream.println("Estimated valid kmers 4+: "+(long)filter[0].estimateUniqueKmers(prehashes, 4));

			//Comprehension (verified termination): currentPass is incremented once per call (line ~230) and every
			//recursive branch is gated on it, so recursion is bounded and terminates. auto mode (prepasses<0): recurse
			//while currentPass is even, else only while currentPass<5 AND the estimate exceeds capacity -> hard cap at 5 passes.
			//Fixed mode: recurse while currentPass<prepasses. onePass short-circuits before here (returns in the if-branch).
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
	 * Both DISPLAY_STATS and the instance showStats flag must be enabled.
	 * DISPLAY_PROGRESS gates memory output and showSpeed gates throughput within this summary.
	 * @param t Timer stopped again here when statistics are enabled
	 * @param added Loader return total, printed as Unique Kmers
	 */
	public final void showStats(Timer t, long added){

		if(!DISPLAY_STATS || !showStats){return;}

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
	 * Expands the first '#' into mate numbers only when no explicit mate is supplied
	 * and the literal first path does not exist. Each backend load gets the configured limit.
	 * Restores temporary validation/interleaving settings on failure as well as success;
	 * earlier table mutations and loaded-count contributions are not rolled back.
	 * @return Sum of the backend loader returns across primary and extra inputs
	 */
	private final long loadKmers(){
		assert(allocated) : "processInput must allocate tables before invoking the input loader";
		kmersLoaded=0;
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=false;
		try{
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
			//KMR109: retain the post-primary policy and restore it even when an extra load fails.
			final boolean oldForceInterleaved=FASTQ.FORCE_INTERLEAVED;
			FASTQ.FORCE_INTERLEAVED=false;
			try{
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
			}finally{
				FASTQ.FORCE_INTERLEAVED=oldForceInterleaved;
			}
			return kmersLoaded;
		}finally{
			Read.VALIDATE_IN_CONSTRUCTOR=vic;
		}
	}

	/**
	 * Loads one primary/mate path pair into allocated tables using the concrete backend.
	 * Counter, pairing and worker behavior follow that backend's contract.
	 * @param fname1 First input path
	 * @param fname2 Optional mate path, or null
	 * @return Backend loader contribution to the processInput total
	 */
	public abstract long loadKmers(String fname1, String fname2);

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Rebuilds supported backend tables after discarding entries at or below the limit.
	 * Backend ownership and regeneration preconditions still apply.
	 * @param limit Inclusive count threshold for removal
	 * @return Number of removed entries, summed across tables
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
	 * @param histMax Final bin index; backends may fold higher counts into this bin
	 * @return Array where index represents count and value represents frequency
	 */
	public abstract long[] fillHistogram(int histMax);

	/**
	 * Allocates GC totals grouped by the backend's scalar count bins.
	 * @param histMax Final bin index, also receiving higher counts
	 * @return New array of length histMax+1; each key contributes its G/C base count once
	 */
	public final long[] fillGcCounts(int histMax){
		long[] gcCounts=new long[histMax+1];
		countGC(gcCounts, histMax);
		return gcCounts;
	}

	/**
	 * Adds per-key GC totals to count bins without clearing existing destination values.
	 * Tables must remain stable. Joins all parallel workers and propagates failure before
	 * reducing their accumulators. Direct traversal propagates failures without rollback.
	 * @param gcCounts Accumulator of length max+1
	 * @param max Final bin index, also receiving contributions from higher counts
	 */
	public final void countGC(long[] gcCounts, int max){
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

		//KMR108: join every worker before reporting failure or reducing private accumulators.
		Throwable failure=null;
		for(CountGCThread t : al){
			while(t.getState()!=Thread.State.TERMINATED){
				try{
					t.join();
				}catch(InterruptedException e){
					// Ideally handle interrupt, but standard BBTools practice often ignores it here
				}
			}
			if(failure==null){failure=t.failure;}
		}
		rethrowGCFailure(failure);
		for(CountGCThread t : al){
			// Reduction step: Sum thread-local counts into the main array
			for(int i=0; i<gcCounts.length; i++){
				gcCounts[i]+=t.myGcCounts[i];
			}
		}
	}

	/** Preserves worker failure types after all GC workers have terminated. */
	private static void rethrowGCFailure(final Throwable failure){
		if(failure==null){return;}
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		if(failure instanceof Error){throw (Error)failure;}
		throw new RuntimeException("K-mer GC counting worker failed.", failure);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Accumulates GC totals for a strided subset of shards into a private array. */
	private class CountGCThread extends Thread{
		/** Creates worker tid_ of threads_; max_ is the inclusive last histogram bin. */
		CountGCThread(int tid_, int threads_, int max_){
			tid=tid_;
			threads=threads_;
			// Each thread gets a private array to avoid atomic contention
			myGcCounts=new long[max_+1];
			max=max_;
		}

		/** Traverses assigned shards and retains any failure for the joining coordinator. */
		@Override
		public void run(){
			try{
				KmerTableInterface[] tables=tables();
				for(int i=tid; i<tables.length; i+=threads){
					// This calls the existing single-threaded method in AbstractKmerTable
					tables()[i].countGC(myGcCounts, max);
				}
			}catch(Throwable t){
				failure=t;
			}
		}

		/** Zero-based worker index and stride between assigned table indices. */
		final int tid, threads;
		/** Inclusive last histogram bin. */
		final int max;
		/** Private accumulator, merged only after every worker finishes successfully. */
		final long[] myGcCounts;
		/** Captured failure read by the coordinator after termination. */
		volatile Throwable failure;
	}

	/**
	 * Divides per-bin GC totals by max(key count,1) times the full key length.
	 * Inputs must describe matching bins and nonoverflowing base totals. Fractions
	 * are not validated or clamped; a consistently empty bin yields zero.
	 *
	 * @param counts K-mer frequency counts per bin
	 * @param gcCounts GC base totals, with at least counts.length entries
	 * @return New per-bin fraction array; valid consistent inputs yield values from zero to one
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

	/** Initializes backend ownership after loading; does not make arbitrary table access thread-safe. */
	public abstract void initializeOwnership();

	/** Clears backend ownership storage; subsequent owner access may require reinitialization. */
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
	 * @return Number of valid windows reported by the backend, not the count-list length
	 */
	public final int fillCounts(byte[] bases, IntList counts, Kmer kmer){
		return fillSpecificCounts(bases, counts, null, kmer);
	}

	/**
	 * Fills count values for selected k-mer start positions in a sequence.
	 * The concrete backends differ in traversal extent. Missing selected k-mers
	 * are reported as zero.
	 *
	 * @param bases Sequence bases to process
	 * @param counts Output list, cleared and refilled by the backend
	 * @param positions Selected zero-based window starts, or null for all
	 * @param kmer Reusable matching-geometry scratch key for long keys; ignored by short-key implementation
	 * @return Backend's valid-window total, distinct from the number of emitted list entries
	 */
	public abstract int fillSpecificCounts(byte[] bases, IntList counts, BitSet positions, Kmer kmer);

	/**
	 * Recomputes counts across windows overlapping the span of changed bases.
	 * Keeps unaffected entries; invalid recomputed windows receive zero. Backend
	 * lookup results are retained without the fillSpecificCounts missing-value clamp.
	 *
	 * @param bases Sequence bases to process
	 * @param counts Existing count list covering the sequence, updated in place
	 * @param kmer Reusable matching-geometry scratch for long keys; ignored for short keys
	 * @param changed Nonempty mask of changed base indices; read without modification
	 * @return Number of valid recomputed windows, not total list length
	 */
	public abstract int regenerateCounts(byte[] bases, IntList counts, Kmer kmer, BitSet changed);

	/*--------------------------------------------------------------*/
	/*----------------       Printing Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Traverses tables in the caller thread and writes formatted k-mer/count text.
	 * Concrete implementations own and finish the output writer; null filenames skip output.
	 *
	 * @param fname Output file path
	 * @param mincount Inclusive minimum k-mer count to include
	 * @param maxcount Inclusive maximum k-mer count to include
	 * @param printTime Whether to display timing information
	 * @param remaining Shared record allowance, or null for unlimited output; may become negative
	 * @return true when an output error is reported; false on normal completion or skipped output
	 */
	public abstract boolean dumpKmersAsBytes(String fname, int mincount, int maxcount, boolean printTime, AtomicLong remaining);
	/**
	 * Writes formatted k-mer/count text using workers where the concrete backend permits.
	 * May fall back to caller-thread traversal; implementations own writer completion.
	 *
	 * @param fname Output file path
	 * @param mincount Inclusive minimum k-mer count to include
	 * @param maxcount Inclusive maximum k-mer count to include
	 * @param printTime Whether to display timing information
	 * @param remaining Shared record allowance, or null for unlimited output; may become negative
	 * @return true when a worker or output error is reported; false on normal completion or skipped output
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
	 * @param printZeros Whether to emit empty bins; bin zero is always omitted
	 * @param printTime Whether to print final write timing; intermediate step timers always print
	 * @param smooth Whether to apply progressive smoothing to histogram
	 * @param calcGC Whether to calculate GC content statistics
	 * @param doLogScale Whether to apply logarithmic scaling
	 * @param logWidth Width parameter for log scaling
	 * @param logPasses Number of passes for log scaling
	 * @param smoothRadius Radius for smoothing operations
	 * @return Two arrays: counts after optional smoothing, and GC totals reconstructed
	 * from those counts and the original per-bin GC fractions (null when calcGC is false).
	 * Log scaling affects the optional output column, not these returned arrays.
	 * @throws RuntimeException If the output writer reports failure after completion
	 */
	public final long[][] makeKhist(String fname, int cols, int max, boolean printHeader, boolean printZeros, boolean printTime,
			boolean smooth, boolean calcGC, boolean doLogScale, double logWidth, int logPasses, int smoothRadius){
		Timer t=new Timer(), t2=new Timer();
		boolean printSteps=true;
		long[] ca=fillHistogram(max);
		if(printSteps){t2.stopAndStart("fillHistogram:");}
		float[] gcHist=null;
		if(calcGC){
//			assert(false) : max+", "+ca.length;
			long[] gc=(calcGC ? fillGcCounts(max) : null);
			if(printSteps){t2.stopAndStart("fillGcCounts:");}
//			assert(false) : max+", "+ca.length+", "+gc.length;
			gcHist=makeGcHistogram(ca, gc);
			if(printSteps){t2.stopAndStart("makeGcHistogram:");}
		}

		long[] logScale=null;

		if(smooth){
			ca=CallPeaks.smoothProgressive(ca, smoothRadius);
			if(printSteps){t2.stopAndStart("smoothProgressive:");}
		}
		if(doLogScale){
			logScale=CallPeaks.logScale(ca, logWidth, 1, logPasses);
			if(printSteps){t2.stopAndStart("logScale:");}
		}

		long[][] ret=new long[2][];
		ret[0]=ca;
		if(gcHist!=null){
			final int k=kbig();
			ret[1]=new long[ca.length];
			for(int i=1; i<ca.length; i++){
				ret[1][i]=Math.round(ca[i]*gcHist[i]*k);
			}
			if(printSteps){t2.stopAndStart("gcHist:");}
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
		if(printSteps){t2.stopAndStart("fileWrite:");}
		final boolean writerError=bsw.poisonAndWait();
		errorState|=writerError;
		t.stop();
		if(printTime){outstream.println("Histogram Write Time:       \t"+t);}
		if(writerError){throw new RuntimeException("Failed to write k-mer histogram: "+fname);}
		return ret;
	}

	/** Borrows the backing shards; callers must not mutate their membership during traversal. */
	public abstract <T extends KmerTableInterface> T[] tables();

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Enables capacity and load-summary reporting when DISPLAY_STATS is also true.
	 * Does not suppress separate allocation progress or prefilter diagnostics. */
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

	/** Main-loader fragment limit per input; negative means unlimited.
	 * Prefilter counting applies it to the first primary input but treats remaining files as unlimited extras. */
	public long maxReads=-1;

	/** Key capacity of each worker's per-shard insertion buffer. */
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

	/** Number of completed non-one-pass prefilter passes; not reset by processInput. */
	private int currentPass=0;
	/** Target total prefilter passes; negative selects automatic scheduling. */
	protected int prepasses=1;

	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the effective full k-mer length, which packed/layout selection may adjust.
	 * @return Number of sequence symbols represented by a complete key */
	public abstract int kbig();
	/**
	 * Returns the planned prefilter memory budget for a pass, not measured live usage.
	 * @param pass Zero-based pass number used to select the alternating budget
	 * @return Planned budget in bytes
	 */
	public abstract long filterMemory(int pass);
	/** Returns the planned main-table memory budget, not allocated or live heap usage.
	 * @return Planned budget in bytes */
	public abstract long tableMemory();
	/** Returns the constructor's capacity estimate under its layout and load assumptions.
	 * @return Estimated number of keys, not a measured live limit */
	public abstract long estimatedKmerCapacity();
	/** Returns whether error correction from overlapping mate pairs was requested.
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
	/** Returns whether key selection considers reverse complements.
	 * @return Whether forward and reverse-complement keys share canonical matching */
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

	/** Undefined-base limit per read in the short-key loader; MAX_VALUE disables this filter. */
	public static int maxNs=Integer.MAX_VALUE;
	/** Post-trimming minimum read length enforced by the short-key loader. */
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

	/** Distinct backend lookup sentinels for absence and exhausted primary probing. */
	public static final int NOT_PRESENT=AbstractKmerTable.NOT_PRESENT, HASH_COLLISION=AbstractKmerTable.HASH_COLLISION;
	/** Constant indicating k-mer has no thread ownership */
	public static final int NO_OWNER=AbstractKmerTable.NO_OWNER;

	/** Default minimum probability threshold for k-mer inclusion */
	public static double defaultMinprob=0;

	/** Enables end-base-masked short-key routing and compatible fast extension lookups.
	 * Long-key hash geometry additionally uses the separate Kmer.MASK_CORE setting. */
	public static boolean MASK_CORE=false;
	/** Requests middle-base masking in supporting short-key paths; prefiltering rejects it. */
	public static boolean MASK_MIDDLE=false;
	/** Enables eligible fast left/right extension-count lookup paths when core masking is active. */
	public static boolean FAST_FILL=true;

}
