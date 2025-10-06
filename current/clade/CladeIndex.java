package clade;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.concurrent.ConcurrentHashMap;

import prok.GeneCaller;
import shared.Parse;
import shared.Shared;

/**
 * Indexes Clade objects by GC content for efficient similarity searching.
 * Implements a GC-bucket based search strategy to quickly find the best matching reference 
 * for a query Clade using a two-stage comparison process.
 * 
 * @author Brian Bushnell
 * @date April 19, 2024
 */
public class CladeIndex implements Cloneable {

	/**
	 * Creates a new CladeIndex containing the specified Clades.
	 * Each Clade will be indexed by its GC content for fast retrieval.
	 * 
	 * @param coll Collection of Clades to index
	 */
	public CladeIndex(Collection<Clade> coll) {
		for(Clade c : coll) {add(c);}
		sort();
	}
	
	private void sort() {
		for(ArrayList<Clade> list : gcDex) {
			if(list!=null) {Shared.sort(list);}
		}
	}

	/**
	 * Creates a clone of this CladeIndex.
	 * The clone shares the same indexed Clades but has zeroed comparison counters.
	 * 
	 * @return A clone of this CladeIndex
	 */
	public CladeIndex clone() {
		CladeIndex ci=null;
		try {
			ci = (CladeIndex) super.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ci.comparisons=ci.slowComparisons=0;
		return ci;
	}

	/**
	 * Loads and indexes Clades from the specified reference files.
	 * 
	 * @param ref Array of reference file paths
	 * @return A new CladeIndex containing Clades from the references
	 */
	public static CladeIndex loadIndex(String... ref) {
		return loadIndex(Arrays.asList(ref));
	}

	/**
	 * Loads and indexes Clades from the specified reference files.
	 * 
	 * @param ref Collection of reference file paths
	 * @return A new CladeIndex containing Clades from the references
	 */
	public static CladeIndex loadIndex(Collection<String> ref) {
		CladeLoader loader=new CladeLoader();
		ConcurrentHashMap<Integer, Clade> map=loader.load(ref, null);
		CladeIndex index=new CladeIndex(map.values());
		return index;
	}

	/**
	 * Parses a command-line parameter and sets the corresponding configuration.
	 * Supports a wide range of tuning parameters for the search algorithm and comparison methods.
	 * 
	 * @param arg Original argument string (unused)
	 * @param a Parameter name
	 * @param b Parameter value
	 * @return true if the parameter was recognized and processed, false otherwise
	 */
	public static boolean parse(String arg, String a, String b) {
		if(a.equals("steps") || a.equals("maxsteps")){
			maxSteps=Integer.parseInt(b);
		}else if(a.equals("maxk") || a.equals("kmax")){
			Comparison.maxK=Clade.MAXK=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("callSSU")){
			Clade.callSSU=Parse.parseBoolean(b);
		}else if(a.equals("aligner") || a.equals("idaligner")){
			GeneCaller.useIDAligner=(b==null || !("f".equals(b) || "false".equals(b)));
			if(GeneCaller.useIDAligner) {aligner.Factory.setType(b);}
		}else if(a.equals("heapsize") || a.equals("heap")){
			heapSize=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("comparisonCutoffMult") || a.equals("ccm")){
			Comparison.setComparisonCutoffMult(Float.parseFloat(b));
		}else if(a.equalsIgnoreCase("comparisonCutoffMult2") || a.equals("ccm2")){
			Comparison.setComparisonCutoffMult2(Float.parseFloat(b));
		}else if(a.equalsIgnoreCase("mink5bases") || a.equalsIgnoreCase("k5len")){
			Comparison.minK5Bases=Parse.parseKMG(b);
		}else if(a.equalsIgnoreCase("mink4bases") || a.equalsIgnoreCase("k4len")){
			Comparison.minK4Bases=Parse.parseKMG(b);
		}else if(a.equalsIgnoreCase("k4Mult")){
			Comparison.k4Mult=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("k3Mult")){
			Comparison.k3Mult=Float.parseFloat(b);
		}else if(a.equals("gcdelta") || a.equals("gcdif")){
			gcDelta=Float.parseFloat(b);
		}else if(a.equals("strdelta") || a.equals("strdif")){
			strDelta=Float.parseFloat(b);
		}else if(a.equals("hhdelta") || a.equals("hhdif")){
			hhDelta=Float.parseFloat(b);
		}else if(a.equals("cagadelta") || a.equals("cagadif")){
			cagaDelta=Float.parseFloat(b);
		}else if(a.equals("gcmult")){
			gcMult=Float.parseFloat(b);
		}else if(a.equals("strmult")){
			strMult=Float.parseFloat(b);
		}else if(a.equals("hhmult")){
			hhMult=Float.parseFloat(b);
		}else if(a.equals("cagamult")){
			cagaMult=Float.parseFloat(b);
		}else if(a.equals("abs") || a.equals("absdif")){
			Comparison.method=Comparison.ABS;
		}else if(a.equals("abscomp")){
			Comparison.method=Comparison.ABSCOMP;
		}else if(a.equals("cos") || a.equals("cosine")){
			Comparison.method=Comparison.COS;
		}else if(a.equals("hel") || a.equals("hellinger")){
			Comparison.method=Comparison.HEL;
		}else if(a.equals("euc") || a.equals("euclidian")){
			Comparison.method=Comparison.EUC;
		}else if(a.equalsIgnoreCase("earlyExit") || a.equals("ee")){
			Comparison.earlyExit=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("calcCladeEntropy") || a.equals("entropy")){
			CladeObject.calcCladeEntropy=Parse.parseBoolean(b);
		}else if(a.equals("banself")){
			banSelf=Parse.parseBoolean(b);
		}else if(a.equals("includeself")){
			banSelf=!Parse.parseBoolean(b);
		}else if(a.equals("linear")){
			LINEAR_SEARCH=Parse.parseBoolean(b);
		}else if(a.equals("binary")){
			LINEAR_SEARCH=!Parse.parseBoolean(b);
		}else {
			return false;
		}
		return true;
	}

	/**
	 * Adds a Clade to the index.
	 * The Clade is indexed by its GC content, rounded to the nearest percentage.
	 * 
	 * @param c The Clade to add
	 */
	public void add(Clade c) {
		int idx=Math.round(c.gc*100);
		ArrayList<Clade> list=gcDex[idx];
		if(list==null) {gcDex[idx]=list=new ArrayList<Clade>(16);}
		list.add(c);
		cladesLoaded++;
	}
	
	public void setFromBest(final Clade clade) {
		Comparison comp=findSingleBest(clade);
		if(comp==null) {return;}
		clade.name=comp.ref.name;
		clade.taxID=comp.ref.taxID;
		clade.lineage=comp.ref.lineage;
	}
	
	public Comparison findSingleBest(final Clade c) {
		ArrayList<Comparison> list=findBest(c, 1);
		return list==null || list.isEmpty() ? null : list.get(0);
	}
	
	/**
	 * Finds the best match for a given Clade in the index.
	 * Uses a GC-bucket based search strategy that first checks Clades with similar GC content,
	 * then expands outward to neighboring GC buckets if needed.
	 *
	 * @param c The query Clade to find matches for
	 * @param maxHits The number of results to return
	 * @return An ordered list of Comparison objects containing the best matches found.
	 */
	public ArrayList<Comparison> findBest(final Clade c, final int maxHits) {
		assert(c.finished());
		final int center=Math.round(c.gc*100);
		final Comparison temp=new Comparison();
		final ComparisonHeap heap=new ComparisonHeap(maxHits);
		{
			Comparison best=new Comparison();
			best.query=c;
			heap.offer(best);
		}
		synchronized(heap) {
			synchronized(c) {//Probably unnecessary...
				findBestBinary(c, gcDex[center], heap, temp, center);
				for(int i=1, lim=maxSteps*2+2; i<=maxSteps || (heap.worst().ref==null && i<lim); i++) {
					int low=center-i, high=center+i;
					if(low>=0) {findBestBinary(c, gcDex[low], heap, temp, low);}
					if(high<gcDex.length) {findBestBinary(c, gcDex[high], heap, temp, high);}
				}
			}
		}
//		if(heap.worst().ref==null) {//Nothing was found
//			assert(heap.size()==1) : heap;
//			return new ArrayList<Comparison>(1);//Could return null but this can be used as a placeholder
//		}
		return heap.toList();
	}

	/**
	 * Searches a specific GC bucket for the best match to the given Clade.
	 * Uses a two-stage comparison process with early filtering to improve performance.
	 * Updates the best match if a better one is found.
	 * 
	 * @param a The query Clade to match
	 * @param list List of Clades in this GC bucket
	 * @param best Current best matches (will be updated if better matches are found)
	 * @param temp Temporary comparison object for reuse
	 * @param gcLevel The GC percentage (0-100) of this bucket
	 */
	private void findBestLinear(Clade a, ArrayList<Clade> list, ComparisonHeap heap,
		Comparison temp, int gcLevel) {
		if(list==null || list.isEmpty()) {return;}
		//		System.err.println("\nSearching a list of size "+list.size());
		Comparison worst=heap.worst();
		float k5Limit=worst.k5dif;
		float gcLimit=worst.gcdif+Math.min(gcDelta, k5Limit*gcMult);
		float strLimit=worst.strdif+Math.min(strDelta, k5Limit*strMult);
		float hhLimit=worst.hhdif+Math.min(hhDelta, k5Limit*hhMult);
		float cagaLimit=worst.cagadif+Math.min(cagaDelta, k5Limit*cagaMult);
		//Early exit because GC won't match this list
		if(Math.abs((gcLevel*0.01f)-a.gc)>gcLimit+0.005f) {return;}
		for(Clade b : list) {//TODO: binary search using hh
			if(b==a || (b.taxID==a.taxID && banSelf)) {continue;}
			//			System.err.println("Comparing to "+b);
			comparisons++;
			if(!temp.quickCompare(a, b, gcLimit, strLimit, cagaLimit)) {continue;}
			slowComparisons++;
			float ret=temp.slowCompare(a, b, k5Limit);//ret is not currently used
			//			System.err.println("Comparison: "+temp);
			boolean added=heap.offer(temp);
			if(added) {
				//				System.err.println("***New worst!");
				worst=heap.worst();
				k5Limit=worst.k5dif;		
				gcLimit=worst.gcdif+Math.min(gcDelta, k5Limit*gcMult);
				strLimit=worst.strdif+Math.min(strDelta, k5Limit*strMult);
				hhLimit=worst.hhdif+Math.min(hhDelta, k5Limit*hhMult);
				cagaLimit=worst.cagadif+Math.min(cagaDelta, k5Limit*cagaMult);
			}
		}
	}

	/**
	 * Searches a specific GC bucket for the best match to the given Clade.
	 * Uses a two-stage comparison process with early filtering to improve performance.
	 * Updates the best match if a better one is found.
	 * Assumes list is sorted by HH, and expands from the best position using binary search.
	 * 
	 * @param a The query Clade to match
	 * @param list List of Clades in this GC bucket
	 * @param best Current best matches (will be updated if better matches are found)
	 * @param temp Temporary comparison object for reuse
	 * @param gcLevel The GC percentage (0-100) of this bucket
	 */
	private void findBestBinary(Clade a, ArrayList<Clade> list, ComparisonHeap heap,
		Comparison temp, int gcLevel) {
		if(list==null || list.isEmpty()) {return;}
		if(LINEAR_SEARCH) {findBestLinear(a, list, heap, temp, gcLevel);return;}
		//		System.err.println("\nSearching a list of size "+list.size());
		Comparison worst=heap.worst();
		float k5Limit=worst.k5dif;
		float gcLimit=worst.gcdif+Math.min(gcDelta, k5Limit*gcMult);
		float strLimit=worst.strdif+Math.min(strDelta, k5Limit*strMult);
		float hhLimit=worst.hhdif+Math.min(hhDelta, k5Limit*hhMult);
		float cagaLimit=worst.cagadif+Math.min(cagaDelta, k5Limit*cagaMult);
		//Early exit because GC won't match this list
		if(Math.abs((gcLevel*0.01f)-a.gc)>gcLimit+0.005f) {return;}
		final int center=binarySearchHH(list, a.hh);
		for(int i=center; i>=0; i--) {
			Clade b=list.get(i);
			if(b==a || (b.taxID==a.taxID && banSelf)) {continue;}
			comparisons++;
			boolean pass=temp.quickCompare(a, b, gcLimit, strLimit, cagaLimit);
			if(temp.hhdif>hhLimit && i<center) {break;}//Never break at center
			if(!pass) {continue;}
			slowComparisons++;
			float ret=temp.slowCompare(a, b, k5Limit);//ret is not currently used
			boolean added=heap.offer(temp);
			if(added) {
				worst=heap.worst();
				k5Limit=worst.k5dif;		
				gcLimit=worst.gcdif+Math.min(gcDelta, k5Limit*gcMult);
				strLimit=worst.strdif+Math.min(strDelta, k5Limit*strMult);
				hhLimit=worst.hhdif+Math.min(hhDelta, k5Limit*hhMult);
				cagaLimit=worst.cagadif+Math.min(cagaDelta, k5Limit*cagaMult);
//				System.err.println("B gcLimit="+gcLimit+", strLimit="+strLimit+", hhLimit="+hhLimit);
			}
		}
		for(int i=center+1; i<list.size(); i++) {
			Clade b=list.get(i);
			if(b==a || (b.taxID==a.taxID && banSelf)) {continue;}
			comparisons++;
			boolean pass=temp.quickCompare(a, b, gcLimit, strLimit, cagaLimit);
			if(temp.hhdif>hhLimit) {break;}
			if(!pass) {continue;}
			slowComparisons++;
			float ret=temp.slowCompare(a, b, k5Limit);//ret is not currently used
			boolean added=heap.offer(temp);
			if(added) {
				worst=heap.worst();
				k5Limit=worst.k5dif;		
				gcLimit=worst.gcdif+Math.min(gcDelta, k5Limit*gcMult);
				strLimit=worst.strdif+Math.min(strDelta, k5Limit*strMult);
				hhLimit=worst.hhdif+Math.min(hhDelta, k5Limit*hhMult);
				cagaLimit=worst.cagadif+Math.min(cagaDelta, k5Limit*cagaMult);
//				System.err.println("C gcLimit="+gcLimit+", strLimit="+strLimit+", hhLimit="+hhLimit);
			}
		}
	}

	/** Returns index of the closest element */	
	public static final int binarySearchHH(ArrayList<Clade> list, final float key) {
		final int length=(list==null ? 0 : list.size());
		if(length<2) {return 0;}
		else if(length==2) {
			return Math.abs(key-list.get(0).hh)<=Math.abs(key-list.get(1).hh) ? 0 : 1;
		}
		int a=0, b=length-1;
		while(b>a){
			final int mid=(a+b)/2;
			final float f=list.get(mid).hh;
			if(key<f){b=mid;}
			else if(key>f){a=mid+1;}
			else{return mid;}
		}
		assert(a==b) : a+", "+b;
		if(a==0 || a==length-1) {return a;}
		float dif1=Math.abs(key-list.get(a-1).hh);
		float dif2=Math.abs(key-list.get(a).hh);
		float dif3=Math.abs(key-list.get(a+1).hh);
		if(dif1<dif2) {return a-1;}
		if(dif3<dif2) {return a+1;}
		return a;
	}

	/**
	 * Gets the number of Clades indexed.
	 * 
	 * @return Number of Clades in the index
	 */
	public int size() {
		return cladesLoaded;
	}

	/** Array of Clade lists indexed by GC percentage (0-100) */
	@SuppressWarnings("unchecked")
	final ArrayList<Clade>[] gcDex=new ArrayList[101];

//	/** Array of Clade lists indexed by GC percentage (0-100) */
//	@SuppressWarnings("unchecked")
//	final ArrayList<Clade>[][] gcDex2=new ArrayList[101][101];

	/** Number of Clades loaded in this index */
	int cladesLoaded=0;
	/** Total number of quick comparisons performed */
	long comparisons=0;
	/** Number of detailed comparisons that passed the quick filter */
	long slowComparisons=0;

	/** Number of intermediate comparisons to retain */
	static int heapSize=1;
	/** Whether to exclude Clades with the same taxonomic ID from matches */
	static boolean banSelf=false;
	/** Maximum number of GC buckets to search in each direction */
	static int maxSteps=6;
	/** Maximum allowed GC difference (absolute) */
	static float gcDelta=0.05f;
	/** Maximum allowed strandedness difference (absolute) */
	static float strDelta=0.12f;
	/** Maximum allowed hh difference (absolute) */
	static float hhDelta=0.025f;
	/** Maximum allowed caga difference (absolute) */
	static float cagaDelta=0.017f;
	/** Multiplier for dynamic GC difference threshold based on k-mer similarity */
	static float gcMult=0.5f; //These are optimized for ABS; higher is safer
	/** Multiplier for dynamic strandedness difference threshold based on k-mer similarity */
	static float strMult=1.2f;
	/** Multiplier for dynamic hh difference threshold based on k-mer similarity */
	static float hhMult=0.5f;
	/** Multiplier for dynamic caga difference threshold based on k-mer similarity */
	static float cagaMult=0.8f;
	
	static boolean LINEAR_SEARCH=false;

}