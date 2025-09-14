package bin;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

import shared.Tools;

public class BinMap extends BinObject implements Iterable<Cluster> {
	
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
	
	void addAll(Collection<? extends Bin> bins, int minSize) {
		Key key=new Key();
		for(Bin b : bins) {
			if(b.size()<minSize) {
				residual.add(b);
			}else {
				add(b, key);
			}
		}
	}
	
	Cluster add(Bin a, Key key) {
		Cluster c=a.toCluster();
		if(key==null) {key=new Key();}
		key.set(a);
		ArrayList<Cluster> list=getOrMakeList(key);
		synchronized(list) {list.add(c);}
		return c;
	}
	
	public ArrayList<Cluster> getOrMakeList(Key key){
		ArrayList<Cluster> list=map.get(key);
		if(list==null) {
			map.putIfAbsent((Key)(key.clone()), new ArrayList<Cluster>(8));
			list=map.get(key);
			minGridGC=Tools.min(minGridGC, key.gcLevel);
			maxGridGC=Tools.max(maxGridGC, key.gcLevel);
			minGridDepth=Tools.min(minGridDepth, key.covLevel, key.covLevel2);
			maxGridDepth=Tools.max(maxGridDepth, key.covLevel, key.covLevel2);
			assert(minGridGC<=maxGridGC) : minGridGC+", "+maxGridGC+", "+key;
			assert(minGridDepth<=maxGridDepth);
		}
		return list;
	}
	
	public Cluster findBestCluster(Bin a, int minSizeToCompare, Key key, int matrixRange, Oracle oracle) {
		if(key==null) {key=new Key();}
		oracle.clear();
		final float gc=a.gc();
		final float depth=a.depth(0);
		final float depth2=a.depth(1);
		key.set(a);
		
		float mult=Binner.sizeAdjustMult(a.size());
		final float maxDepthRatio=1+(oracle.maxDepthRatio0-1)*mult;
		final float maxGCDif=oracle.maxGCDif0*mult;
		assert(maxDepthRatio>=1) : maxDepthRatio;
		final int minDepthLevel=Tools.max(minGridDepth, 
				key.covLevel-matrixRange, Key.quantizeDepth(depth/maxDepthRatio));
		final int maxDepthLevel=Tools.min(maxGridDepth,
				key.covLevel+matrixRange, Key.quantizeDepth(depth*maxDepthRatio));
		final int minGCLevel=Tools.max(minGridGC, key.gcLevel-matrixRange, Key.quantizeGC(gc-maxGCDif));
		final int maxGCLevel=Tools.min(maxGridGC, key.gcLevel+matrixRange, Key.quantizeGC(gc+maxGCDif));
		
		int minDepthLevel2=0;
		int maxDepthLevel2=0;
		if(a.numDepths()>1) {
			minDepthLevel2=Tools.max(minGridDepth, key.covLevel2-matrixRange, Key.quantizeDepth(depth2/maxDepthRatio));
			maxDepthLevel2=Tools.min(maxGridDepth, key.covLevel2+matrixRange, Key.quantizeDepth(depth2*maxDepthRatio));
		}
		
		if(verbose) {
			System.err.println("Using params minSizeToCompare="+minSizeToCompare+
					", maxKmerDif="+oracle.max4merDif0+", maxDepthRatio="+maxDepthRatio+
					", maxProduct="+oracle.maxProduct0+", maxGCDif="+maxGCDif+
					", maxCovariance="+oracle.maxCovariance0+", matrixRange="+matrixRange);
		}
		
		for(int depthLevel2=minDepthLevel2; depthLevel2<=maxDepthLevel2; depthLevel2++) {
			for(int depthLevel=minDepthLevel; depthLevel<=maxDepthLevel; depthLevel++) {
				for(int gcLevel=minGCLevel; gcLevel<=maxGCLevel; gcLevel++) {
					if(verbose) {System.err.println("Looking at depth "+depthLevel+", gc "+gcLevel);}
					key.setLevel(gcLevel, depthLevel, depthLevel2);
					ArrayList<Cluster> list=map.get(key);
					if(list!=null) {
						int idx=findBestBinIndex(a, list, minSizeToCompare, oracle);
						if(verbose && idx>=0) {System.err.println("***Set best to "+oracle.best.id());}
					}
				}
			}
		}
		return (Cluster)oracle.best;
	}
	
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
	
	public void clear(boolean clearResidual) {
		map.clear();
		if(clearResidual) {residual.clear();}	
		minGridGC=999999;
		maxGridGC=0;
		minGridDepth=999999;
		maxGridDepth=0;
	}
	
	public boolean isValid() {
		assert(isValid(contigList, true));
		assert(isValid(residual, false));
		for(ArrayList<Cluster> list : map.values()) {
			assert(isValid(list, false));
		}
		return true;
	}
	
	public int countClusters() {
		int sum=0;
		for(ArrayList<Cluster> list : map.values()) {sum+=list.size();}
		return sum;
	}
	
	public ConcurrentHashMap<Key, ArrayList<Cluster>> map=
			new ConcurrentHashMap<Key, ArrayList<Cluster>>(2000, 0.6f, 32);
	public ArrayList<Bin> residual=new ArrayList<Bin>();//Should really be contigs
	public ArrayList<Contig> contigList;
	private int minGridGC=999999;
	private int maxGridGC=0;
	private int minGridDepth=999999;
	private int maxGridDepth=0;
	
//	public AtomicLong fastComparisons=new AtomicLong(0);
//	public AtomicLong slowComparisons=new AtomicLong(0);
	
}
