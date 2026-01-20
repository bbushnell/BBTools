package bin;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

import shared.Tools;

/**
 * Multi-dimensional hash map for organizing genomic bins by GC content and coverage depth.
 * Uses quantized GC and depth levels as keys to efficiently store similar genomic bins.
 * Enables fast nearest-neighbor searches within specified similarity tolerances.
 *
 * @author Brian Bushnell
 * @date Feb 4, 2025
 */
public class BinMap extends BinObject implements Iterable<Cluster> {
	
	/** Constructs a BinMap with the specified list of contigs.
	 * @param contigs List of contigs for reference */
	public BinMap(ArrayList<Contig> contigs) {contigList=contigs;}
	
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
	 * Adds all bins from a collection to the map.
	 * Bins below minimum size are added to residual collection.
	 * @param bins Collection of bins to add
	 * @param minSize Minimum size threshold for map inclusion
	 */
	void addAll(Collection<? extends Bin> bins, int minSize) {
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
	 * Converts a bin to a cluster and adds it to the appropriate hash bucket.
	 * Thread-safe operation using synchronized access to bucket lists.
	 *
	 * @param a Bin to convert and add
	 * @param key Hash key for bucket selection (created if null)
	 * @return The newly created cluster
	 */
	Cluster add(Bin a, Key key) {
		Cluster c=a.toCluster();
		if(key==null) {key=Key.makeKey();}
		key.set(a);
		ArrayList<Cluster> list=getOrMakeList(key);
		synchronized(list) {list.add(c);}
		return c;
	}
	
	/**
	 * Retrieves or creates the cluster list for the specified key.
	 * Updates global grid bounds for GC and depth levels when creating new buckets.
	 * Thread-safe using putIfAbsent for concurrent access.
	 *
	 * @param key Hash key identifying the bucket
	 * @return List of clusters for this key
	 */
	public ArrayList<Cluster> getOrMakeList(Key key){
		ArrayList<Cluster> list=map.get(key);
		if(list==null) {
			map.putIfAbsent((Key)(key.clone()), new ArrayList<Cluster>(8));
			list=map.get(key);
			minGridGC=Tools.min(minGridGC, key.gcLevel);
			maxGridGC=Tools.max(maxGridGC, key.gcLevel);
			minGridDim2=Math.min(minGridDim2, key.dim2);
			maxGridDim2=Math.max(maxGridDim2, key.dim2);
			minGridDim3=Math.min(minGridDim3, key.dim3);
			maxGridDim3=Math.max(maxGridDim3, key.dim3);
			minGridDim4=Math.min(minGridDim4, key.dim4);
			maxGridDim4=Math.max(maxGridDim4, key.dim4);
			minGridDim5=Math.min(minGridDim5, key.dim5);
			maxGridDim5=Math.max(maxGridDim5, key.dim5);
			assert(minGridGC<=maxGridGC) : minGridGC+", "+maxGridGC+", "+key;
			assert(minGridDim2<=maxGridDim2);
			assert(minGridDim3<=maxGridDim3);
			assert(minGridDim4<=maxGridDim4);
			assert(minGridDim5<=maxGridDim5);
		}
		return list;
	}
	
	/**
	 * Searches for the best matching cluster within specified GC/depth tolerances.
	 * Uses oracle similarity scoring with size-adjusted thresholds.
	 * Explores a 3D grid centered on the bin's quantized GC/depth coordinates.
	 *
	 * @param a Bin seeking a matching cluster
	 * @param minSizeToCompare Minimum cluster size for comparison
	 * @param key Working key for hash operations (created if null)
	 * @param matrixRange Search radius in quantized space
	 * @param oracle Similarity calculator that tracks best match
	 * @return Best matching cluster or null if none found
	 */
	public Cluster findBestCluster(Bin a, int minSizeToCompare, Key key, int matrixRange, Oracle oracle) {
		if(key==null) {key=Key.makeKey();}
		oracle.clear();
		final float gc=a.gc();
		key.set(a);
		
		float mult=Binner.sizeAdjustMult(a.size());
		final float maxDepthRatio=1+(oracle.maxDepthRatio0-1)*mult;
		final float maxGCDif=oracle.maxGCDif0*mult;
//		final int minGCLevel=Tools.max(minGridGC, key.gcLevel-matrixRange, Key.quantizeGC(gc-maxGCDif));
//		final int maxGCLevel=Tools.min(maxGridGC, key.gcLevel+matrixRange, Key.quantizeGC(gc+maxGCDif));

		final int minGCLevel=key.lowerBoundDim1(gc, matrixRange, minGridGC, maxGCDif);
		final int maxGCLevel=key.upperBoundDim1(gc, matrixRange, maxGridGC, maxGCDif);
		
		assert(maxDepthRatio>=1) : maxDepthRatio;
		assert(minGCLevel<=key.gcLevel || minGCLevel==999999 || minGridGC>key.gcLevel) : 
			key.gcLevel+", "+minGCLevel+", "+minGridGC+", "+
			matrixRange+", "+gc+", "+maxGCDif+", "+Key.quantizeGC(gc-maxGCDif);
		assert(maxGCLevel>=key.gcLevel || minGCLevel==999999 || maxGridGC<key.gcLevel) : 
			key.gcLevel+", "+maxGCLevel+", "+maxGridGC+", "+
			matrixRange+", "+gc+", "+maxGCDif+", "+Key.quantizeGC(gc+maxGCDif);
		
		final int minDim2=key.lowerBoundDim2(a, matrixRange, minGridDim2, maxGCDif, maxDepthRatio);
		final int maxDim2=key.upperBoundDim2(a, matrixRange, maxGridDim2, maxGCDif, maxDepthRatio);
		
		final int minDim3=key.lowerBoundDim3(a, matrixRange, minGridDim3, maxGCDif, maxDepthRatio);
		final int maxDim3=key.upperBoundDim3(a, matrixRange, maxGridDim3, maxGCDif, maxDepthRatio);
		
		final int minDim4=key.lowerBoundDim4(a, matrixRange, minGridDim4, maxGCDif, maxDepthRatio);
		final int maxDim4=key.upperBoundDim4(a, matrixRange, maxGridDim4, maxGCDif, maxDepthRatio);
		
		final int minDim5=key.lowerBoundDim5(a, matrixRange, minGridDim5, maxGCDif, maxDepthRatio);
		final int maxDim5=key.upperBoundDim5(a, matrixRange, maxGridDim5, maxGCDif, maxDepthRatio);

		assert(minGCLevel<=maxGCLevel || maxGCLevel==0 || minGridGC>key.gcLevel || maxGridGC<key.gcLevel) : 
			minGCLevel+", "+maxGCLevel+", "+key.gcLevel+", "+matrixRange;
		assert(minDim2<=maxDim2 || maxDim2==0 || minGridDim2>key.dim2 || maxGridDim2<key.dim2) : 
			key.getClass()+"\n"+
			"r="+minDim2+"-"+maxDim2+", key="+key+", mr="+matrixRange+", gr="+minGridDim2+"-"+maxGridDim2+
			"\nlbd2="+key.lowerBoundDim2(a, matrixRange, minGridDim2, maxGCDif, maxDepthRatio)+
			"\nubd2="+key.upperBoundDim2(a, matrixRange, maxGridDim2, maxGCDif, maxDepthRatio);
		assert(minDim3<=maxDim3 || maxDim3==0 || minGridDim3>key.dim3 || maxGridDim3<key.dim3) : 
			minDim3+", "+maxDim3;
		assert(minDim4<=maxDim4 || maxDim4==0 || minGridDim4>key.dim4 || maxGridDim4<key.dim4) : 
			minDim4+", "+maxDim4;
		assert(minDim5<=maxDim5 || maxDim5==0 || minGridDim5>key.dim5 || maxGridDim5<key.dim5) : 
			minDim5+", "+maxDim5;
		
		//Bad assertion, it fails sometimes and is basically wrong, just for testing.
//		assert((minDim3==maxDim3)==(minDim4==maxDim4) || maxDim4==0 || 
//			minGridDim3>key.dim3 || maxGridDim3<key.dim3 || minGridDim4>key.dim4 || maxGridDim4<key.dim4) : 
//			minDim3+"-"+maxDim3+", "+minDim4+"-"+maxDim4+"\n"+
//			minGridDim3+"-"+maxGridDim3+", "+minGridDim4+"-"+maxGridDim4+"\n"+a;
		
		if(verbose) {
			System.err.println("Using params minSizeToCompare="+minSizeToCompare+
					", maxKmerDif="+oracle.max4merDif0+", maxDepthRatio="+maxDepthRatio+
					", maxProduct="+oracle.maxProduct0+", maxGCDif="+maxGCDif+
					", maxCovariance="+oracle.maxCovariance0+", matrixRange="+matrixRange);
		}

		long lookups=0, valid=0;
		for(int dim5=minDim5; dim5<=maxDim5; dim5++) {
			for(int dim4=minDim4; dim4<=maxDim4; dim4++) {
				for(int dim3=minDim3; dim3<=maxDim3; dim3++) {
					for(int dim2=minDim2; dim2<=maxDim2; dim2++) {
						for(int gcLevel=minGCLevel; gcLevel<=maxGCLevel; gcLevel++) {
							if(verbose) {System.err.println("Looking at depth "+dim2+", gc "+gcLevel);}
							key.setLevel(gcLevel, dim2, dim3, dim4, dim5);
							lookups++;
							ArrayList<Cluster> list=map.get(key);
							if(list!=null) {
								valid++;
								int idx=findBestBinIndex(a, list, minSizeToCompare, oracle);
								if(verbose && idx>=0) {System.err.println("***Set best to "+oracle.best.id());}
							}
						}
					}
				}
			}
		}
		hashLookups.addAndGet(lookups);
		hashLookupsValid.addAndGet(valid);
		return (Cluster)oracle.best;
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
	private int findBestBinIndex(Bin a, ArrayList<? extends Bin> clusters, 
			int minSizeToCompare, Oracle oracle) {
		
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
	
	/**
	 * Converts the hash map contents to a flat list of clusters.
	 * Optionally includes residual bins as single-contig clusters.
	 * @param addResidue Whether to include residual bins in the output
	 * @return List containing all clusters from the map
	 */
	ArrayList<Cluster> toList(boolean addResidue){
		ArrayList<Cluster> list=new ArrayList<Cluster>();
		for(Entry<Key, ArrayList<Cluster>> e : map.entrySet()) {
			list.addAll(e.getValue());
		}
		if(addResidue) {
			for(Bin b : residual) {
				list.add(b.cluster()==null ? b.toCluster() : b.cluster());
			}
		}
		return list;
	}
	
	@Override
	public Iterator<Cluster> iterator() {
		return toList(false).iterator();
	}
	
	/**
	 * Clears the hash map and optionally the residual collection.
	 * Resets grid bounds to initial values.
	 * @param clearResidual Whether to also clear the residual bin collection
	 */
	public void clear(boolean clearResidual) {
		map.clear();
		if(clearResidual) {residual.clear();}	
		minGridGC=999999;
		maxGridGC=0;
		minGridDim2=999999;
		maxGridDim2=0;
		minGridDim3=999999;
		maxGridDim3=0;
		minGridDim4=999999;
		maxGridDim4=0;
		minGridDim5=999999;
		maxGridDim5=0;
	}
	
	/**
	 * Validates the internal consistency of the map structure.
	 * Checks contig list, residual collection, and all stored clusters.
	 * @return true if all validations pass
	 */
	public boolean isValid() {
		assert(isValid(contigList, true));
		assert(isValid(residual, false));
		for(ArrayList<Cluster> list : map.values()) {
			assert(isValid(list, false));
		}
		return true;
	}
	
	/** Counts the total number of clusters stored in the map.
	 * @return Total cluster count across all hash buckets */
	public int countClusters() {
		int sum=0;
		for(ArrayList<Cluster> list : map.values()) {sum+=list.size();}
		return sum;
	}
	
	/**
	 * Primary hash map storing clusters indexed by quantized GC and depth levels
	 */
	public ConcurrentHashMap<Key, ArrayList<Cluster>> map=
			new ConcurrentHashMap<Key, ArrayList<Cluster>>(2000, 0.6f, 32);
	/** Collection of bins that don't meet size thresholds for clustering */
	public ArrayList<Bin> residual=new ArrayList<Bin>();//Should really be contigs
	/** Reference list of all contigs being processed */
	public ArrayList<Contig> contigList;
	/** Minimum GC level observed in the current grid */
	private int minGridGC=999999;
	/** Maximum GC level observed in the current grid */
	private int maxGridGC=0;
	private int minGridDim2=9999999;//Generally depth1 or HH
	private int maxGridDim2=0;
	private int minGridDim3=9999999;
	private int maxGridDim3=0;
	private int minGridDim4=9999999;
	private int maxGridDim4=0;
	private int minGridDim5=9999999;
	private int maxGridDim5=0;

	AtomicLong hashLookups=new AtomicLong(0);
	AtomicLong hashLookupsValid=new AtomicLong(0);
	
}
