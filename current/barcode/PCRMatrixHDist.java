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
 * Tracks data about bar code mismatches by position.
 * Uses split barcodes instead of contiguous.
 * 
 * @author Brian Bushnell
 * @date March 22, 2024
 *
 */
public class PCRMatrixHDist extends PCRMatrix implements Accumulator<PCRMatrixHDist.PopThread> {

	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a PCRMatrixHDist with specified barcode configuration.
	 *
	 * @param length1_ Length of first barcode segment
	 * @param length2_ Length of second barcode segment (0 for single barcodes)
	 * @param delimiter_ ASCII value of delimiter character between segments
	 * @param hdistSum_ Whether to sum Hamming distances across segments
	 */
	public PCRMatrixHDist(int length1_, int length2_, int delimiter_, boolean hdistSum_) {
		super(length1_, length2_, delimiter_, hdistSum_);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses static configuration parameters for Hamming distance calculations.
	 *
	 * @param arg Complete argument string
	 * @param a Parameter name (e.g., "maxhdist", "clearzone")
	 * @param b Parameter value
	 * @return true if parameter was recognized and parsed, false otherwise
	 */
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
	
	@Override
	public boolean parse(String arg, String a, String b) {
		return false;
	}
	
	/** Performs post-parsing static initialization (currently no operations) */
	public static void postParseStatic(){}
	
	/*--------------------------------------------------------------*/
	/*----------------            HDist             ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public Barcode findClosest(String s) {
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
	public Barcode findClosest(String s, int maxHDist, int clearzone) {
		return length2<1 ? findClosestSingleHDist(s, maxHDist, clearzone) : 
			findClosestDualHDist(s, maxHDist, clearzone);
	}

	@Override
	public void makeProbs() {
		throw new RuntimeException("This class does not support this method.");
	}

	@Override
	public void initializeData() {}
	
	@Override
	public void refine(Collection<Barcode> codeCounts, long minCount) {}
	
	@Override
	public HashMap<String, String> makeAssignmentMap(Collection<Barcode> codeCounts, long minCount) {
		Timer t=new Timer();
		assert(expectedList!=null && expectedList.size()>0) : expectedList;
		assert(codeCounts!=null);
		ArrayList<Barcode> list=highpass(codeCounts, minCount);
		HashMap<String, String> map=new HashMap<String, String>(Tools.max(200, list.size()/10));
		totalCounted=totalAssigned=totalAssignedToExpected=0;
		final long ops=list.size()*(long)expectedList.size();
		if(list.size()<2 || ops<100000 || Shared.threads()<2) {//Singlethreaded mode
			for(Barcode query : list) {
				final String s=query.name;
				assert(s.length()==counts.length);
				Barcode ref=findClosest(s);
				final long count=query.count();
				totalCounted+=count;
				if(ref!=null) {
					totalAssigned+=count;
					if(ref.expected==1) {
						totalAssignedToExpected+=count;
						map.put(s, ref.name);
					}
				}
			}
		}else {
			populateCountsMT(list, maxHDist0, clearzone0, map);
		}
		t.stop();
		if(verbose) {
			if(verbose) {
				System.err.println("Pair Assignment Rate:   \tTotal\tGood\tBad");
			}
			System.err.println(String.format("Final Assignment Rate:  \t%.4f\t%.4f\t%.6f", 
					assignedFraction(), expectedFraction(), chimericFraction())+"\t"+t.timeInSeconds(2)+"s");
		}
		return map;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Populating          ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void populateCounts(ArrayList<Barcode> list, long minCount) {
		assert(minCount<2) : "TODO";
		assert(expectedList!=null && expectedList.size()>0) : expectedList;
		assert(list!=null);
		final long ops=list.size()*(long)expectedList.size();
		if(list.size()<2 || ops<100000 || Shared.threads()<2) {
			populateCountsST(list, maxHDist0, clearzone0);
		}else {
			populateCountsMT(list, maxHDist0, clearzone0, null);
		}
	}

	/**
	 * Single-threaded population of barcode assignment counts.
	 * @param countList List of barcodes to process
	 * @param maxHDist Maximum Hamming distance for matching
	 * @param clearzone Minimum distance gap between best and second-best matches
	 */
	private void populateCountsST(ArrayList<Barcode> countList,
			int maxHDist, int clearzone) {
		for(Barcode query : countList) {
			final String s=query.name;
			assert(s.length()==counts.length);
			Barcode ref=findClosest(s, maxHDist, clearzone);
			add(s, ref, query.count());
		}
	}

	/**
	 * Multi-threaded population of barcode assignment counts using worker threads.
	 * Distributes work across available threads with stride-based partitioning.
	 *
	 * @param list List of barcodes to process
	 * @param maxHDist Maximum Hamming distance for matching
	 * @param clearzone Minimum distance gap between matches
	 * @param map Optional assignment map to populate during processing
	 */
	private void populateCountsMT(ArrayList<Barcode> list,
			int maxHDist, int clearzone, HashMap<String, String> map) {
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
		errorState&=!success;
		
		//Do anything necessary after processing
		if(localCounts && map!=null) {
			for(PopThread pt : alpt) {
				synchronized(pt) {map.putAll(pt.map);}
			}
		}
	}
	
	@Override
	public void populateUnexpected() {assert(false) : "Not valid.";}
	
	@Override
	public ByteBuilder toBytesProb(ByteBuilder bb) {
		throw new RuntimeException("This class does not support this method.");
	}
	
	/** Returns true indicating this implementation is valid for use */
	protected boolean valid() {return true;}
	
	/*--------------------------------------------------------------*/

	/** Worker thread for multi-threaded barcode assignment processing.
	 * Each thread processes a subset of barcodes using stride-based work partitioning. */
	final class PopThread extends Thread {

		/**
		 * Constructs a worker thread for barcode processing.
		 *
		 * @param list_ List of barcodes to process
		 * @param maxHDist_ Maximum Hamming distance threshold
		 * @param clearzone_ Clearzone threshold for disambiguation
		 * @param map_ Optional assignment map to populate
		 * @param tid_ Thread ID for stride-based work partitioning
		 * @param threads_ Total number of worker threads
		 */
		public PopThread(ArrayList<Barcode> list_,
				int maxHDist_, int clearzone_, HashMap<String, String> map_, int tid_, int threads_) {
			list=list_;
			maxHDist=maxHDist_;
			clearzone=clearzone_;
			tid=tid_;
			threads=threads_;
			map=(map_==null ? null : localCounts ? new HashMap<String, String>() : map_);
			countsT=(localCounts ? new long[length][5][5] : null);
		}

		@Override
		public void run() {
			for(int i=tid; i<list.size(); i+=threads) {
				Barcode query=list.get(i);
				final String s=query.name;
				assert(s.length()==length);
				Barcode ref=findClosest(s, maxHDist, clearzone);
				if(localCounts) {
					addT(s, ref, query.count());
					if(map!=null && ref!=null && ref.expected==1) {map.put(s, ref.name);}
				}else {
					synchronized(counts) {
						add(s, ref, query.count());
						if(map!=null && ref!=null && ref.expected==1) {map.put(s, ref.name);}
					}
				}
			}
		}
		
		/**
		 * Thread-local method to add barcode assignment data to thread-specific counts.
		 * Updates position-specific mismatch matrices and assignment statistics.
		 *
		 * @param query Query barcode string
		 * @param ref Reference barcode (null if no match)
		 * @param count Number of occurrences to add
		 */
		public void addT(String query, Barcode ref, long count) {
			assert(ref==null || ref.length()==countsT.length);
			for(int i=0; i<query.length(); i++) {
				final int q=query.charAt(i), r=(ref==null ? 'N' : ref.charAt(i));
				final byte xq=baseToNumber[q], xr=baseToNumber[r];
				countsT[i][xq][xr]+=count;
			}
			totalCountedT+=count;
			if(ref!=null) {
				ref.incrementSync(count);
				totalAssignedT+=count;
				totalAssignedToExpectedT+=ref.expected*count;
			}
		}

		/** List of barcodes assigned to this thread for processing */
		final ArrayList<Barcode> list;
		/** Maximum Hamming distance threshold for barcode matching */
		final int maxHDist;
		/** Minimum distance gap required between best and second-best matches */
		final int clearzone;
		/** Thread ID used for stride-based work partitioning */
		final int tid;
		/** Total number of worker threads in the processing pool */
		final int threads;
		/** Optional assignment map for storing barcode-to-reference mappings */
		final HashMap<String, String> map;
		
		/**
		 * Thread-local 3D array tracking [position][query_base][ref_base] mismatches
		 */
		final long[][][] countsT;
		/** Thread-local count of total barcodes processed by this thread */
		long totalCountedT;
		/** Thread-local count of barcodes successfully assigned by this thread */
		long totalAssignedT;
		/**
		 * Thread-local count of barcodes assigned to expected references by this thread
		 */
		long totalAssignedToExpectedT;
	}

	@Override
	public final void accumulate(PopThread t) {
		if(localCounts) {
			synchronized(t) {
				Tools.add(counts, t.countsT);
				totalCounted+=t.totalCountedT;
				totalAssigned+=t.totalAssignedT;
				totalAssignedToExpected+=t.totalAssignedToExpectedT;
			}
		}
	}

	@Override
	public boolean success() {
		return !errorState;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Default maximum Hamming distance threshold for barcode matching */
	static int maxHDist0=2;
	/**
	 * Default clearzone threshold requiring minimum distance gap between matches
	 */
	static int clearzone0=1;
	
}
