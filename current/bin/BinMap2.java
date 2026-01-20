package bin;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

/**
 * 5-Dimensional BinMap optimized with a "Sliced" architecture.
 * Organizes clusters primarily by GC level (Slice), then by 4 Depth dimensions.
 * Reduces search space by maintaining independent bounds for each GC slice.
 * 
 * @author Brian Bushnell
 * @contributor Amber
 * @date January 20, 2026
 */
public class BinMap2 extends BinObject implements Iterable<ArrayList<Cluster>> {
	
	public BinMap2(ArrayList<Contig> contigs) {
		contigList=contigs;
	}

	/**
	 * Adds a bin to the map, creating the appropriate GC slice if needed.
	 */
	public Cluster add(Bin a, Key key) {
		Cluster c=a.toCluster();
		if(key==null) {key=Key.makeKey();}
		key.set(a);
		
		// Update global GC bounds (cheap)
		updateGCBounds(key.gcLevel);
		
		// Get or create the slice for this GC level
		BinMapSlice slice=slices.get(key.gcLevel);
		if(slice==null) {
			slices.putIfAbsent(key.gcLevel, new BinMapSlice(key.gcLevel));
			slice=slices.get(key.gcLevel);
		}
		
		// Delegate to slice to update inner bounds and store cluster
		slice.add(c, key);
		
		return c;
	}
	
	private synchronized void updateGCBounds(int gc) {
		minGridGC=Math.min(minGridGC, gc);
		maxGridGC=Math.max(maxGridGC, gc);
	}

	/**
	 * Finds best matching cluster.
	 * Iterates over GC slices first, skipping entire slices that don't exist.
	 */
	public Cluster findBestCluster(Bin a, long minSizeToCompare, Key key, int matrixRange, Oracle oracle) {
		if(key==null) {key=Key.makeKey();}
		oracle.clear();
		final float gc=a.gc();
		key.set(a); // Sets initial key state
		
		float mult=Binner.sizeAdjustMult(a.size());
		final float maxDepthRatio=1+(oracle.maxDepthRatio0-1)*mult;
		final float maxGCDif=oracle.maxGCDif0*mult;

		// Calculate Global Search Bounds
		final int minGCLevel=key.lowerBoundDim1(gc, matrixRange, minGridGC, maxGCDif);
		final int maxGCLevel=key.upperBoundDim1(gc, matrixRange, maxGridGC, maxGCDif);
		
		final int minDim2=key.lowerBoundDim2(a, matrixRange, 0, maxGCDif, maxDepthRatio);
		final int maxDim2=key.upperBoundDim2(a, matrixRange, 999999, maxGCDif, maxDepthRatio);
		
		final int minDim3=key.lowerBoundDim3(a, matrixRange, 0, maxGCDif, maxDepthRatio);
		final int maxDim3=key.upperBoundDim3(a, matrixRange, 999999, maxGCDif, maxDepthRatio);
		
		final int minDim4=key.lowerBoundDim4(a, matrixRange, 0, maxGCDif, maxDepthRatio);
		final int maxDim4=key.upperBoundDim4(a, matrixRange, 999999, maxGCDif, maxDepthRatio);
		
		final int minDim5=key.lowerBoundDim5(a, matrixRange, 0, maxGCDif, maxDepthRatio);
		final int maxDim5=key.upperBoundDim5(a, matrixRange, 999999, maxGCDif, maxDepthRatio);

		// Outer Loop: GC Level (The Slices)
		for(int gcLevel=minGCLevel; gcLevel<=maxGCLevel; gcLevel++) {
			BinMapSlice slice=slices.get(gcLevel);
			
			// OPTIMIZATION: If no data exists at this GC level, skip entirely.
			if(slice==null) {continue;}
			
			// Delegate to slice, passing the GLOBAL requested bounds.
			// The slice will intersect these with its LOCAL bounds.
			slice.findBest(a, key, oracle, 
					minDim2, maxDim2, 
					minDim3, maxDim3, 
					minDim4, maxDim4, 
					minDim5, maxDim5, 
					minSizeToCompare, verbose);
		}
		
		return (Cluster)oracle.best;
	}
	
	/** Static helper to process a single list of clusters */
	public static int findBestBinIndex(Bin a, ArrayList<Cluster> list, long minSizeToCompare, Oracle oracle, boolean verbose) {
		int bestIdx=-1;
		// Iterate backwards or forwards? Original was forwards.
		for(int i=0; i<list.size(); i++) {
			Cluster b=list.get(i);
			if(b==null || a==b) {continue;}
			if(b.size()<minSizeToCompare) {break;} // Assuming sorted descending
			
			float f=oracle.similarity(a, b, 1f);
			if(f>oracle.topScore) {
				oracle.best=b;
				oracle.bestIdx=bestIdx=i;
				oracle.topScore=f;
			}
		}
		return bestIdx;
	}

	/**
	 * Iterator that walks through all ArrayLists in all Slices.
	 */
	@Override
	public Iterator<ArrayList<Cluster>> iterator() {
		return new Iterator<ArrayList<Cluster>>() {
			final Iterator<BinMapSlice> sliceIter = slices.values().iterator();
			Iterator<ArrayList<Cluster>> listIter = null;

			@Override
			public boolean hasNext() {
				while(listIter == null || !listIter.hasNext()) {
					if(!sliceIter.hasNext()) return false;
					listIter = sliceIter.next().values().iterator();
				}
				return true;
			}

			@Override
			public ArrayList<Cluster> next() {
				if(!hasNext()) throw new java.util.NoSuchElementException();
				return listIter.next();
			}
		};
	}
	
	public void clear() {
		slices.clear();
		minGridGC=999999;
		maxGridGC=0;
	}

	public ConcurrentHashMap<Integer, BinMapSlice> slices = new ConcurrentHashMap<>();
	public ArrayList<Contig> contigList;
	private int minGridGC=999999;
	private int maxGridGC=0;
	
	// Stats
	public AtomicLong hashLookups=new AtomicLong(0); // You can wire these up if you want
}
