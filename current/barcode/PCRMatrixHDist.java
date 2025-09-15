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
	
	public PCRMatrixHDist(int length1_, int length2_, int delimiter_, boolean hdistSum_) {
		super(length1_, length2_, delimiter_, hdistSum_);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	public static void postParseStatic(){}
	
	/*--------------------------------------------------------------*/
	/*----------------            HDist             ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public Barcode findClosest(String s) {
		return length2<1 ? findClosestSingleHDist(s, maxHDist0, clearzone0) : 
			findClosestDualHDist(s, maxHDist0, clearzone0);
	}
	
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

	private void populateCountsST(ArrayList<Barcode> countList,
			int maxHDist, int clearzone) {
		for(Barcode query : countList) {
			final String s=query.name;
			assert(s.length()==counts.length);
			Barcode ref=findClosest(s, maxHDist, clearzone);
			add(s, ref, query.count());
		}
	}

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
	
	protected boolean valid() {return true;}
	
	/*--------------------------------------------------------------*/

	final class PopThread extends Thread {

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

		final ArrayList<Barcode> list;
		final int maxHDist;
		final int clearzone;
		final int tid;
		final int threads;
		final HashMap<String, String> map;
		
		final long[][][] countsT;
		long totalCountedT;
		long totalAssignedT;
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
	
	static int maxHDist0=2;
	static int clearzone0=1;
	
}
