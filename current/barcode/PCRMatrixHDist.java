package barcode;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * Assigns single or dual barcodes by Hamming distance and a second-best distance margin.
 * Expected names and split candidates must be populated before searches. Counts represent
 * observation weights, not distinct-name counts. populateCounts accumulates them;
 * makeAssignmentMap replaces mapping statistics for its current filtered observations.
 * Candidate arrays, query lists and static settings must remain stable during worker execution.
 * A failed worker causes caller failure after joining; partial mutations are not rolled back
 * and must not be used as a completed result. Concurrent operations on one instance are unsupported.
 *
 * @author Brian Bushnell
 * @date March 22, 2024
 */
public class PCRMatrixHDist extends PCRMatrix implements Accumulator<PCRMatrixHDist.PopThread>{

	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Allocates positional counts; expected candidates are populated separately.
	 * @param length1_ First index length
	 * @param length2_ Second index length, or 0 for a single index
	 * @param delimiter_ Separator character, or 0 when absent
	 * @param hdistSum_ Sum dual-index distances rather than taking their maximum
	 */
	public PCRMatrixHDist(int length1_, int length2_, int delimiter_, boolean hdistSum_){
		super(length1_, length2_, delimiter_, hdistSum_);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses shared distance/margin defaults before processing; the inherited placeholder flag is inert. */
	public static boolean parseStatic(String arg, String a, String b){
		if(a.equals("maxhdist") || a.equals("hdist") || a.equals("maxhdist0") || a.equals("hdist0")){
			maxHDist0=Integer.parseInt(b);
		}else if(a.equals("clearzone") || a.equals("cz") || a.equals("clearzone0") || a.equals("cz0")){
			clearzone0=Integer.parseInt(b);
		}else if(a.equals("parse_flag_goes_here")){
			//set something
		}else{
			return false;
		}
		return true;
	}

	/**
	 * Parses instance-specific parameters (currently no instance parameters supported).
	 *
	 * @param arg Complete argument string
	 * @param a Parameter name
	 * @param b Parameter value
	 * @return false (no instance parameters currently supported)
	 */
	@Override
	public boolean parse(String arg, String a, String b){
		return false;
	}

	/** No additional Hamming configuration finalization is required. */
	public static void postParseStatic(){}

	/*--------------------------------------------------------------*/
	/*----------------            HDist             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Finds the closest matching barcode using default distance and clearzone parameters.
	 * @param s Query barcode string
	 * @return Borrowed candidate object, or null when distance/margin checks reject assignment
	 */
	@Override
	public Barcode findClosest(String s){
		return length2<1 ? findClosestSingleHDist(s, maxHDist0, clearzone0) :
			findClosestDualHDist(s, maxHDist0, clearzone0);
	}

	/**
	 * Finds the closest matching barcode with specified distance and clearzone thresholds.
	 * Routes to single or dual barcode search based on configuration.
	 *
	 * @param s Query barcode string
	 * @param maxHDist Maximum Hamming distance allowed for matching
	 * @param clearzone Minimum distance gap required between best and second-best matches
	 * @return Closest matching Barcode object or null if no match within thresholds
	 */
	public Barcode findClosest(String s, int maxHDist, int clearzone){
		return length2<1 ? findClosestSingleHDist(s, maxHDist, clearzone) :
			findClosestDualHDist(s, maxHDist, clearzone);
	}

	/**
	 * Throws RuntimeException as this Hamming distance implementation does not support probability calculations
	 */
	@Override
	public void makeProbs(){
		throw new RuntimeException("This class does not support this method.");
	}

	/** Initializes data structures (no-op for this implementation) */
	@Override
	public void initializeData(){}

	/**
	 * Refines barcode assignments (no-op for this Hamming distance implementation).
	 * @param codeCounts Collection of barcode counts to refine
	 * @param minCount Minimum count threshold for refinement
	 */
	@Override
	public void refine(Collection<Barcode> codeCounts, long minCount){}

	/**
	 * Creates mapping from observed barcodes to assigned reference barcodes.
	 * Uses highpass filtering, which bypasses the threshold for small collections. Resets the
	 * positional counts, totals and all candidate weights, and returns only assignments to
	 * expected references. Chooses serial or worker processing from list size, expected count
	 * and configured thread count.
	 * Serial and worker paths produce the same statistics. Caller owns exclusive mutation;
	 * observed Barcode objects must be distinct from this matrix's owned reference candidates.
	 *
	 * @param codeCounts Collection of observed barcode counts
	 * @param minCount Inclusive count threshold when highpass filtering is enabled
	 * @return New map from assigned observed names to expected reference names; rejected/chimeric names omitted
	 * @throws IllegalStateException if a worker fails; partial state is not a completed result
	 */
	@Override
	public HashMap<String, String> makeAssignmentMap(Collection<Barcode> codeCounts, long minCount){
		Timer t=new Timer();
		assert(expectedList!=null && expectedList.size()>0) : expectedList;
		assert(codeCounts!=null);
		ArrayList<Barcode> list=highpass(codeCounts, minCount);
		HashMap<String, String> map=new HashMap<String, String>(Tools.max(200, list.size()/10));
		resetMappingCounts();
		final long ops=list.size()*(long)expectedList.size();
		if(list.size()<2 || ops<100000 || Shared.threads()<2){//Singlethreaded mode
			if(verbose){System.err.println(list.size());}
			for(Barcode query : list){
				final String s=query.name;
				assert(s.length()==counts.length);
				Barcode ref=findClosest(s);
				final long count=query.count();
				//BAR102: count the same positional/reference statistics as the worker paths.
				add(s, ref, count);
				if(ref!=null && ref.expected==1){map.put(s, ref.name);}
			}
		}else{
			populateCountsMT(list, maxHDist0, clearzone0, map);
		}
		t.stop();
		if(verbose){
			if(verbose){
				System.err.println("Pair Assignment Rate:   \tTotal\tGood\tBad");
			}
			System.err.println(String.format("Final Assignment Rate:  \t%.4f\t%.4f\t%.6f",
					assignedFraction(), expectedFraction(), chimericFraction())+"\t"+t.timeInSeconds(2)+"s");
		}
		return map;
	}

	/** Replaces mapping-owned statistics before workers start; caller owns exclusive mutation. */
	private void resetMappingCounts(){
		assert(allCodes!=null) : "populateSplitCodes must precede mapping so all candidate counts can be replaced.";
		Tools.fill(counts, 0);
		totalCounted=totalAssigned=totalAssignedToExpected=0;
		for(Barcode ref : expectedList){ref.setCount(0);}
		//BAR102: base clearCounts omits synthetic candidates; mapping must reset their weights too.
		for(Barcode ref : allCodes){ref.setCount(0);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Populating          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Populates barcode assignment counts by finding closest matches for each query barcode.
	 * Adds to existing counts rather than clearing them. Automatically selects serial or worker
	 * execution. Callers own the list and must not modify it or its counts during processing.
	 * @param list List of query barcodes to process
	 * @param minCount Must be less than 2 under assertions; no filtering is implemented here
	 * @throws IllegalStateException if a worker fails; discard the partial result
	 */
	@Override
	public void populateCounts(ArrayList<Barcode> list, long minCount){
		assert(minCount<2) : "TODO";
		assert(expectedList!=null && expectedList.size()>0) : expectedList;
		assert(list!=null);
		final long ops=list.size()*(long)expectedList.size();
		if(list.size()<2 || ops<100000 || Shared.threads()<2){
			populateCountsST(list, maxHDist0, clearzone0);
		}else{
			populateCountsMT(list, maxHDist0, clearzone0, null);
		}
	}

	/** Adds every query's weight using the specified search thresholds; caller owns exclusive mutation. */
	private void populateCountsST(ArrayList<Barcode> countList,
			int maxHDist, int clearzone){
		for(Barcode query : countList){
			final String s=query.name;
			assert(s.length()==counts.length);
			Barcode ref=findClosest(s, maxHDist, clearzone);
			add(s, ref, query.count());
		}
	}

	/**
	 * Multi-threaded population of barcode assignment counts using worker threads.
	 * Distributes work by stride, joins workers, accumulates successful local results and merges
	 * expected assignments. Failure sets errorState and throws before map merging; earlier shared
	 * increments or successful-worker accumulation may remain. errorState is sticky across calls.
	 *
	 * @param list List of barcodes to process
	 * @param maxHDist Maximum Hamming distance for matching
	 * @param clearzone Minimum distance gap between matches
	 * @param map Optional assignment map to populate during processing
	 */
	private void populateCountsMT(ArrayList<Barcode> list,
			int maxHDist, int clearzone, HashMap<String, String> map){
		//Do anything necessary prior to processing

		//Determine how many threads may be used
		final int threads=Tools.mid(1, Tools.min(matrixThreads, Shared.threads()), list.size()/8);

		//Fill a list with PopThreads
		ArrayList<PopThread> alpt=new ArrayList<PopThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new PopThread(list, maxHDist, clearzone, map, i, threads));
		}

		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState|=!success;
		if(!success){
			throw new IllegalStateException("Barcode matrix worker failed; partial counts and assignments must not be used.");
		}

		//Do anything necessary after processing
		if(localCounts && map!=null){
			for(PopThread pt : alpt){
				synchronized(pt){map.putAll(pt.map);}
			}
		}
	}

	/**
	 * Unsupported operation: asserts under -ea and otherwise does nothing.
	 */
	@Override
	public void populateUnexpected(){assert(false) : "Not valid.";}

	/**
	 * Throws RuntimeException as probability output is not supported.
	 * @param bb ByteBuilder for output
	 * @return Never returns (always throws exception)
	 */
	@Override
	public ByteBuilder toBytesProb(ByteBuilder bb){
		throw new RuntimeException("This class does not support this method.");
	}

	/** The public Hamming implementation passes the constructor's validity check. */
	protected boolean valid(){return true;}

	/*--------------------------------------------------------------*/

	/** Owns one stride of a borrowed query list. Local mode isolates positional counts/map;
	 * reference Barcode increments remain shared. Shared mode serializes updates on counts. */
	final class PopThread extends Thread{

		/** Borrows input and selects local count/map storage or the shared caller map at construction. */
		public PopThread(ArrayList<Barcode> list_,
				int maxHDist_, int clearzone_, HashMap<String, String> map_, int tid_, int threads_){
			list=list_;
			maxHDist=maxHDist_;
			clearzone=clearzone_;
			tid=tid_;
			threads=threads_;
			map=(map_==null ? null : localCounts ? new HashMap<String, String>() : map_);
			countsT=(localCounts ? new long[length][5][5] : null);
		}

		/**
		 * Processes the assigned stride and publishes success only after every query completes.
		 * An uncaught failure leaves success false for the joining accumulator to reject.
		 */
		@Override
		public void run(){
			for(int i=tid; i<list.size(); i+=threads){
				Barcode query=list.get(i);
				final String s=query.name;
				assert(s.length()==length);
				Barcode ref=findClosest(s, maxHDist, clearzone);
				if(localCounts){
					addT(s, ref, query.count());
					if(map!=null && ref!=null && ref.expected==1){map.put(s, ref.name);}
				}else{
					synchronized(counts){
						add(s, ref, query.count());
						if(map!=null && ref!=null && ref.expected==1){map.put(s, ref.name);}
					}
				}
			}
			success=true;
		}

		/**
		 * Thread-local method to add barcode assignment data to thread-specific counts.
		 * Updates position-specific mismatch matrices and assignment statistics.
		 *
		 * @param query Query barcode string
		 * @param ref Reference barcode (null if no match)
		 * @param count Number of occurrences to add
		 */
		public void addT(String query, Barcode ref, long count){
			assert(ref==null || ref.length()==countsT.length);
			for(int i=0; i<query.length(); i++){
				final int q=query.charAt(i), r=(ref==null ? 'N' : ref.charAt(i));
				final byte xq=baseToNumber[q], xr=baseToNumber[r];
				countsT[i][xq][xr]+=count;
			}
			totalCountedT+=count;
			if(ref!=null){
				ref.incrementSync(count);
				totalAssignedT+=count;
				totalAssignedToExpectedT+=ref.expected*count;
			}
		}

		/** Borrowed stable input; each query is assigned to exactly one worker by index stride. */
		final ArrayList<Barcode> list;
		/** Per-search mismatch limit. */
		final int maxHDist;
		/** Required second-best distance margin. */
		final int clearzone;
		/** First list index owned by this worker. */
		final int tid;
		/** Stride between this worker's list indices. */
		final int threads;
		/** Null, worker-local map, or shared map according to the caller and localCounts setting. */
		final HashMap<String, String> map;

		/** Worker-owned positional counts in local mode; null in shared mode. */
		final long[][][] countsT;
		/** Observation weights accumulated locally, merged only after successful completion. */
		long totalCountedT;
		/** Local weight assigned to any candidate. */
		long totalAssignedT;
		/** Local weight assigned to an expected candidate. */
		long totalAssignedToExpectedT;
		/** Published only after every assigned query completes; thread termination alone is not success. */
		volatile boolean success=false;
	}

	/**
	 * Accumulates thread-local results from a worker thread into global counts.
	 * Called by ThreadWaiter after joining, under the accumulator lock. Failed workers set
	 * errorState and are skipped. Shared-count mode has already applied its updates and is not
	 * merged again. Skipping a failed worker does not roll back shared reference increments.
	 * @param t PopThread containing thread-local results to accumulate
	 */
	@Override
	public final void accumulate(PopThread t){
		if(!t.success){errorState=true; return;}
		//Guarded on localCounts DELIBERATELY: in !localCounts mode run() already merged into the shared
		//`counts` under synchronized(counts), so there is nothing to accumulate here. Removing this guard
		//(always merging t.countsT) would DOUBLE-COUNT every read in shared-counts mode.
		if(localCounts){
			synchronized(t){
				Tools.add(counts, t.countsT);
				totalCounted+=t.totalCountedT;
				totalAssigned+=t.totalAssignedT;
				totalAssignedToExpected+=t.totalAssignedToExpectedT;
			}
		}
	}

	/** Reports the accumulated error flag; this is not a worker-completion or liveness probe.
	 * @return true if errorState is false */
	@Override
	public boolean success(){
		return !errorState;
	}

	/*--------------------------------------------------------------*/

	/** Process-wide default mismatch limit; configure before concurrent work. */
	static int maxHDist0=2;
	/** Process-wide default best-to-second-best distance margin. */
	static int clearzone0=1;
	/** Local diagnostic switch, distinct from the inherited instance verbosity field. */
	private static final boolean verbose=false;
}
