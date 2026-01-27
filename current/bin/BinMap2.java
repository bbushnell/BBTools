package bin;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicLong;

/**
 * 5-Dimensional BinMap optimized with a "Sliced" architecture and Array-based lookup.
 * Organizes clusters primarily by GC level (Slice), then by 4 Depth dimensions.
 * Uses a fixed array for GC slices to eliminate boxing/hashing overhead.
 * 
 * @author Brian Bushnell
 * @contributor Amber
 * @date January 20, 2026
 */
public class BinMap2 extends BinObject implements Iterable<ArrayList<Cluster>> {
	
	public BinMap2(ArrayList<Contig> contigs) {
		contigList=contigs;
		// Pre-allocate slices. Max GC is 1.0. 
		// If width is 0.02, max index is 50. Size 60 is plenty safe.
		// We could calculate (int)(1/Key.gcLevelWidth)+5, but 100 covers any weird float edge cases.
		int arraySize=(int)(1.0f/Key.gcLevelWidth)+1;
		slices=new BinMapSlice[arraySize];
	}
	
	/** Higher is less stringent; 1.0 is neutral, 0 is exact match */
	public Cluster addOrMerge(Bin a, int minSizeToCompare, int minSizeToMerge, 
			int minSizeToAdd, Oracle oracle, Key key, int matrixRange) {
		if(minSizeToMerge>=0 && a.size()>=minSizeToMerge) {
			Cluster best=findBestCluster(a, minSizeToCompare, key, matrixRange, oracle);
			if(best!=null) {
				best.add(a);
				return best;
			}
		}
		if(minSizeToAdd>=0 && a.size()>=minSizeToAdd) {return add(a, key.set(a));}
		residual.add(a);
		return null;
	}

	/**
	 * Adds a bin to the map, creating the appropriate GC slice if needed (lazy init).
	 * Note: If map build is single-threaded, no synchronization needed here.
	 */
	public Cluster add(Bin a, Key key) {
		Cluster c=a.toCluster();
		if(key==null) {key=Key.makeKey();}
		key.set(a);
		
		int gc=key.gcLevel;
		
		// Update global GC bounds
		minGridGC=Math.min(minGridGC, gc);
		maxGridGC=Math.max(maxGridGC, gc);
		
		// Array access - O(1) and no boxing
		if(slices[gc]==null) {
			synchronized(slices) {
				if(slices[gc]==null) {slices[gc]=new BinMapSlice(gc);}
			}
		}
		
		// Delegate to slice
		slices[gc].add(c, key);
		
		return c;
	}
	
	/**
	 * Adds all bins from a collection to the map.
	 * Bins below minimum size are added to residual collection.
	 */
	public void addAll(Collection<? extends Bin> bins, int minSize) {
		Key key=Key.makeKey();
		for(Bin b : bins) {
			if(b.size()<minSize) {
				residual.add(b);
			}else {
				add(b, key);
			}
		}
	}
	
	/**
	 * Finds best matching cluster.
	 * Iterates over GC slices using array indices.
	 */
	public Cluster findBestCluster(Bin a, long minSizeToCompare, Key key, int matrixRange, Oracle oracle) {
		if(key==null) {key=Key.makeKey();}
		oracle.clear();
		final float gc=a.gc();
		key.set(a);
		
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

		// Prune invalid GC ranges
		if(minGCLevel > maxGCLevel) return null;
		
		// Array Bounds Check to prevent IOOBE
		int start=Math.max(0, minGCLevel);
		int stop=Math.min(slices.length-1, maxGCLevel);

		// Outer Loop: GC Level (Array Iteration)
		for(int gcLevel=start; gcLevel<=stop; gcLevel++) {
			BinMapSlice slice=slices[gcLevel];
			
			// NULL CHECK: The massive optimization
			if(slice==null) {continue;}
			
			// Delegate to slice
			slice.findBest(a, key, oracle, 
					minDim2, maxDim2, 
					minDim3, maxDim3, 
					minDim4, maxDim4, 
					minDim5, maxDim5, 
					minSizeToCompare, verbose);
		}
		
		return (Cluster)oracle.best;
	}
	
	/**
	 * Converts the hash map contents to a flat list of clusters.
	 * Optionally includes residual bins as single-contig clusters.
	 */
	public ArrayList<Cluster> toList(boolean addResidue){
		ArrayList<Cluster> list=new ArrayList<Cluster>();
		
		// Iterate array
		for(int i=0; i<slices.length; i++) {
			if(slices[i]!=null) {
				for(ArrayList<Cluster> subList : slices[i].values()) {
					list.addAll(subList);
				}
			}
		}
		
		if(addResidue) {
			for(Bin b : residual) {
				list.add(b.cluster()==null ? b.toCluster() : b.cluster());
			}
		}
		return list;
	}
	
	/**
	 * Iterator that walks through all ArrayLists in all non-null Slices.
	 */
	@Override
	public Iterator<ArrayList<Cluster>> iterator() {
		return new Iterator<ArrayList<Cluster>>() {
			int sliceIndex=0;
			Iterator<ArrayList<Cluster>> listIter=null;

			@Override
			public boolean hasNext() {
				while(listIter==null || !listIter.hasNext()) {
					if(sliceIndex>=slices.length) return false;
					BinMapSlice slice=slices[sliceIndex++];
					if(slice!=null) {
						listIter=slice.values().iterator();
					}
				}
				return true;
			}

			@Override
			public ArrayList<Cluster> next() {
				if(!hasNext()) {throw new java.util.NoSuchElementException();}
				return listIter.next();
			}
		};
	}
	
	/** Counts the total number of clusters stored in the map.
	 * @return Total cluster count across all hash buckets */
	public int countClusters() {
		int sum=0;
		for(ArrayList<Cluster> list : this) {sum+=list.size();}
		return sum;
	}
	
	public void clear(boolean clearResidual) {
		for(int i=0; i<slices.length; i++) {
			slices[i]=null; // Nuke the slices
		}
		if(clearResidual) {residual.clear();}
		minGridGC=999999;
		maxGridGC=0;
	}
	
	/**
	 * Validates the internal consistency of the map structure.
	 * Checks contig list, residual collection, and all stored clusters.
	 * @return true if all validations pass
	 */
	public boolean isValid() {
		assert(isValid(contigList, true));
		assert(isValid(residual, false));
		for(BinMapSlice slice : slices) {
			assert(slice==null || slice.isValid());
		}
		return true;
	}
	
	public int size() {
		int size=0;
		for(BinMapSlice s : slices) {
			if(s!=null) {size+=s.map.size();}
		}
		return size;
	}
	
	long hashLookups() {
		return BinMapSlice.hashLookups.get();
//		long x=0;
//		for(BinMapSlice slice : slices) {
//			if(slice!=null) {x+=slice.hashLookups.get();}
//		}
//		return x;
	}

	long hashLookupsValid() {
		return BinMapSlice.hashLookupsValid.get();
//		long x=0;
//		for(BinMapSlice slice : slices) {
//			if(slice!=null) {x+=slice.hashLookupsValid.get();}
//		}
//		return x;
	}
	
	public final BinMapSlice[] slices;
	public ArrayList<Bin> residual=new ArrayList<Bin>();
	public ArrayList<Contig> contigList;
	
	private int minGridGC=999999;
	private int maxGridGC=0;
}