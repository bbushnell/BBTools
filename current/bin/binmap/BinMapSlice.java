package bin.binmap;

import java.util.ArrayList;
import java.util.Collection;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

import bin.Bin;
import bin.BinObject;
import bin.Cluster;
import bin.Oracle;

/**
 * A slice of the BinMap corresponding to a single quantized GC level.
 * Maintains its own bounds for dimensions 2-5 to minimize search space.
 * 
 * @author Brian Bushnell
 * @Contributor Amber
 * @date January 20, 2026
 */
public class BinMapSlice {
	
	public BinMapSlice(int gcLevel_) {
		gcLevel=gcLevel_;
	}
	
	/**
	 * Adds a cluster to this slice.
	 * Updates local slice bounds for dimensions 2-5.
	 */
	public synchronized void add(Cluster c, Key key) {
		ArrayList<Cluster> list=map.get(key);
		if(list==null) {
			map.putIfAbsent((Key)key.clone(), new ArrayList<Cluster>(8));
			list=map.get(key);
			
			// Update local bounds for this slice
			minDim2=Math.min(minDim2, key.dim2);
			maxDim2=Math.max(maxDim2, key.dim2);
			minDim3=Math.min(minDim3, key.dim3);
			maxDim3=Math.max(maxDim3, key.dim3);
			minDim4=Math.min(minDim4, key.dim4);
			maxDim4=Math.max(maxDim4, key.dim4);
			minDim5=Math.min(minDim5, key.dim5);
			maxDim5=Math.max(maxDim5, key.dim5);
		}
		synchronized(list) {
			list.add(c);
		}
	}
	
	/**
	 * Searches for best match within this specific GC slice.
	 * Uses local bounds to restrict the inner loop search range.
	 */
	public void findBest(Bin a, Key key, Oracle oracle, 
			int minD2, int maxD2, 
			int minD3, int maxD3, 
			int minD4, int maxD4, 
			int minD5, int maxD5, 
			long minSizeToCompare, boolean verbose) {
		
		// Intersect requested search range with actual data bounds for this slice
		// This is the "Pruning" magic
		final int start2=Math.max(minD2, minDim2);
		final int stop2=Math.min(maxD2, maxDim2);
		
		final int start3=Math.max(minD3, minDim3);
		final int stop3=Math.min(maxD3, maxDim3);
		
		final int start4=Math.max(minD4, minDim4);
		final int stop4=Math.min(maxD4, maxDim4);
		
		final int start5=Math.max(minD5, minDim5);
		final int stop5=Math.min(maxD5, maxDim5);
		
		if(start2>stop2 || start3>stop3 || start4>stop4 || start5>stop5) {return;}

		// Inner loops for dimensions 5->2
		long lookups=0, valid=0;
		for(int dim5=start5; dim5<=stop5; dim5++) {
			for(int dim4=start4; dim4<=stop4; dim4++) {
				for(int dim3=start3; dim3<=stop3; dim3++) {
					for(int dim2=start2; dim2<=stop2; dim2++) {
						
						key.setLevel(gcLevel, dim2, dim3, dim4, dim5);
						lookups++;
						
						ArrayList<Cluster> list=map.get(key);
						if(list!=null) {
							valid++;
							findBestBinIndex(a, list, minSizeToCompare, oracle);
						}
					}
				}
			}
		}
		hashLookups.addAndGet(lookups);
		hashLookupsValid.addAndGet(valid);
//		assert(false) : lookups+", "+valid;
	}
	
	boolean isValid() {
		boolean valid=true;
		for(ArrayList<Cluster> list : map.values()) {
			assert(valid=BinObject.isValid(list, false) && valid);
		}
		return valid;
	}
	
	/**
	 * Finds the best matching bin within a specific cluster list.
	 * Iterates through clusters in size order, updating oracle with best match.
	 * Stops early when clusters become too small for comparison.
	 *
	 * @param a Query bin
	 * @param clusters List of candidate clusters
	 * @param minSizeToCompare Minimum size threshold for comparison
	 * @param oracle Similarity calculator that tracks best match
	 * @return Index of best match or -1 if none found
	 */
	private static int findBestBinIndex(Bin a, ArrayList<? extends Bin> clusters, 
			long minSizeToCompare, Oracle oracle) {
		
		int bestIdx=-1;
		for(int i=0; i<clusters.size(); i++) {
			Bin b=clusters.get(i);
			if(b==null || a==b) {continue;}
			if(b.size()<minSizeToCompare) {break;}
			
			float f=oracle.similarity(a, b, 1f);
			assert(f!=0) : b.id();
			if(verbose) {
				System.err.println("Comparing to "+b.id()+"; score="+oracle.score+"/"+oracle.topScore);
			}
			if(f>oracle.topScore) {
				assert(f>0);//Actually, could be -1; or clear should set to -1
				oracle.best=b;
				oracle.bestIdx=bestIdx=i;
				oracle.topScore=f;
			}
		}
		return bestIdx;
	}
	
	public Collection<ArrayList<Cluster>> values(){return map.values();}
	
	public final int gcLevel;
	public final ConcurrentHashMap<Key, ArrayList<Cluster>> map=new ConcurrentHashMap<Key, ArrayList<Cluster>>(1024, 0.75f, 4);
	
	// Local bounds for this specific GC slice
	private int minDim2=9999999;
	private int maxDim2=0;
	private int minDim3=9999999;
	private int maxDim3=0;
	private int minDim4=9999999;
	private int maxDim4=0;
	private int minDim5=9999999;
	private int maxDim5=0;
	
	// Stats
	static AtomicLong hashLookups=new AtomicLong(0);
	static AtomicLong hashLookupsValid=new AtomicLong(0);
	
	private static final boolean verbose=false;
}