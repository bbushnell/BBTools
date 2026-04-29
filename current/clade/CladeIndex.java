package clade;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import ddl.DDLIndex;
import ddl.DDLRecord;
import parse.Parse;
import prok.GeneCaller;
import shared.Shared;
import tax.TaxTree;

public class CladeIndex implements Cloneable {

	public CladeIndex(Collection<Clade> coll) {
		for(Clade c : coll) {add(c);}
		sort();
	}
	
	private void sort() {
		for(ArrayList<Clade> list : gcDex) {
			if(list!=null) {Shared.sort(list);}
		}
	}

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

	public static CladeIndex loadIndex(String... ref) {
		return loadIndex(Arrays.asList(ref));
	}

	public static CladeIndex loadIndex(Collection<String> ref) {
		CladeLoader loader=new CladeLoader();
		ConcurrentHashMap<Integer, Clade> map=loader.load(ref, null);
		CladeIndex index=new CladeIndex(map.values());
		if(USE_SKETCHES){
			index.cladeMap=map;
		}
		return index;
	}

	public static boolean parse(String arg, String a, String b) {
		if(a.equals("steps") || a.equals("maxsteps")){
			maxSteps=Integer.parseInt(b);
		}else if(a.equals("maxk") || a.equals("kmax")){
			Comparison.maxK=Clade.MAXK=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("callSSU")){
			Clade.callSSU=Parse.parseBoolean(b);
		}else if(a.equals("aligner") || a.equals("idaligner")){
			GeneCaller.useIDAligner=(b==null || !("f".equals(b) || "false".equals(b)));
			if(GeneCaller.useIDAligner) {idaligner.Factory.setType(b);}
		}else if(a.equals("heapsize") || a.equals("heap") || a.equals("buffer")){
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
		}else if(a.equalsIgnoreCase("makeddls") || a.equals("ddls")){
			Clade.MAKE_DDLS=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("ddlk")){
			Clade.DDL_K=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("ddlbuckets")){
			Clade.DDL_BUCKETS=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("sketch") || a.equalsIgnoreCase("usesketches")
				|| a.equals("ddl")){
			USE_SKETCHES=Parse.parseBoolean(b);
			if(USE_SKETCHES){Clade.MAKE_DDLS=true;}
		}else if(a.equalsIgnoreCase("sketchfile") || a.equalsIgnoreCase("ddlfile")
				|| a.equalsIgnoreCase("ddlref") || a.equalsIgnoreCase("sketchref")){
			sketchFile=b;
			USE_SKETCHES=true;
			Clade.MAKE_DDLS=true;
		}else if(a.equalsIgnoreCase("sketchindex") || a.equalsIgnoreCase("sketchidx")
				|| a.equalsIgnoreCase("ddlindex") || a.equals("index")){
			USE_SKETCH_INDEX=Parse.parseBoolean(b);
			if(USE_SKETCH_INDEX){USE_SKETCHES=true; Clade.MAKE_DDLS=true;}
		}else if(a.equalsIgnoreCase("sketchhits") || a.equalsIgnoreCase("maxsketchhits")){
			maxSketchHits=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minsketchhits") || a.equalsIgnoreCase("minsketches")){
			minSketchMatches=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("loadthreads")){
			CladeLoaderMT.loadThreads=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("ddlloadthreads") || a.equalsIgnoreCase("sketchloadthreads")){
			ddlLoadThreads=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("indexthreads")){
			indexThreads=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("bundlesize")){
			assert(false) : "recompile to change "+arg;
			//DDLLoaderMT.RECORDS_PER_BUNDLE=Integer.parseInt(b);
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
		ArrayList<Comparison> results=heap.toList();
		if(Clade.MAKE_DDLS){
			for(Comparison comp : results){comp.compareDDL();}
		}
		if(ddlIndex!=null){
			addSketchInfo(results, c);
		}
		return results;
	}

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

	public int size() {
		return cladesLoaded;
	}

	@SuppressWarnings("unchecked")
	final ArrayList<Clade>[] gcDex=new ArrayList[101];

//	/** Array of Clade lists indexed by GC percentage (0-100) */
//	@SuppressWarnings("unchecked")
//	final ArrayList<Clade>[][] gcDex2=new ArrayList[101][101];

	int cladesLoaded=0;
	/** Total number of quick comparisons performed during search operations */
	long comparisons=0;
	/** Number of detailed comparisons that passed the initial quick filter */
	long slowComparisons=0;

	/** Number of intermediate comparisons to retain in the result heap */
	static int heapSize=1;
	/** Whether to exclude Clades with the same taxonomic ID from match results */
	static boolean banSelf=false;
	/**
	 * Maximum number of GC buckets to search in each direction from the query GC
	 */
	static int maxSteps=6;
	/** Maximum allowed GC content difference for initial filtering */
	static float gcDelta=0.05f;
	/** Maximum allowed strandedness difference for initial filtering */
	static float strDelta=0.12f;
	static float hhDelta=0.025f;
	static float cagaDelta=0.017f;
	static float gcMult=0.5f; //These are optimized for ABS; higher is safer
	static float strMult=1.2f;
	static float hhMult=0.5f;
	static float cagaMult=0.8f;
	
	static boolean LINEAR_SEARCH=false;

	static boolean USE_SKETCHES=false;
	static boolean USE_SKETCH_INDEX=false;
	static int maxSketchHits=5;
	static int minSketchMatches=3;
	static int ddlLoadThreads=0;
	static int indexThreads=0;
	static String sketchFile=null;
	static final String DEFAULT_SKETCH_FILE="refseqSketchDDL.tsv.gz";

	DDLIndex ddlIndex;
	ArrayList<DDLRecord> sketchRecords;
	ConcurrentHashMap<Integer, Clade> cladeMap;

	public void loadSketches(String resourceDir){
		sketchRecords=new ArrayList<DDLRecord>();
		File dir=new File(resourceDir);
		if(!dir.isDirectory()){
			System.err.println("WARNING: Sketch directory not found: "+resourceDir);
			return;
		}
		File[] files=dir.listFiles((d,name)->(name.endsWith(".ddl.gz") || name.endsWith(".ddl")
				|| (name.contains("_ddl") && name.endsWith(".tsv.gz"))));
		if(files==null || files.length==0){
			System.err.println("WARNING: No .ddl.gz files found in "+resourceDir);
			return;
		}
		Arrays.sort(files, (a, b)->Long.compare(b.length(), a.length()));

		final int totalThreads=ddlLoadThreads>0 ? ddlLoadThreads : Math.min(Shared.threads(), 16);
		long t0=System.nanoTime();

		if(files.length<2 || totalThreads<4){
			for(File f : files){
				ArrayList<DDLRecord> records=ddl.DDLLoaderMT.loadFile(f.getAbsolutePath(), Clade.DDL_K);
				sketchRecords.addAll(records);
			}
		}else{
			final int filesParallel=Math.min(files.length, Math.min(4, totalThreads/4));
			final int parseThreadsPerFile=Math.min(16, totalThreads);
			System.err.println("Loading "+files.length+" DDL files: "+filesParallel+" parallel, "+parseThreadsPerFile+" parse threads each.");
			final ExecutorService pool=Executors.newFixedThreadPool(filesParallel);
			@SuppressWarnings("unchecked")
			final Future<ArrayList<DDLRecord>>[] futures=new Future[files.length];
			for(int i=0; i<files.length; i++){
				final File f=files[i];
				futures[i]=pool.submit(()->{
					return ddl.DDLLoaderMT.loadFile(f.getAbsolutePath(), Clade.DDL_K, parseThreadsPerFile, 1);
				});
			}
			pool.shutdown();
			for(int i=0; i<futures.length; i++){
				try{
					sketchRecords.addAll(futures[i].get());
				}catch(Exception e){
					System.err.println("Error loading "+files[i].getName());
					e.printStackTrace();
				}
			}
		}
		long elapsed=System.nanoTime()-t0;
		System.err.println("Loaded "+sketchRecords.size()+" sketches from "+files.length+" files in "+String.format("%.3f", elapsed/1e9)+" seconds.");
	}

	public void finishSketches(ArrayList<DDLRecord> externalRecords){
		if(externalRecords!=null && !externalRecords.isEmpty()){
			sketchRecords=externalRecords;
			attachSketchesToClades(cladeMap);
			if(USE_SKETCH_INDEX){
				buildSketchIndex();
			}
		}
	}

	public void finishSketches(){
		finishSketches(sketchRecords);
	}

	public void attachSketchesToClades(ConcurrentHashMap<Integer, Clade> cladeMap){
		if(sketchRecords==null){return;}
		int attached=0;
		for(DDLRecord rec : sketchRecords){
			Clade c=cladeMap.get(rec.taxID);
			if(c!=null){
				boolean fresh=(c.ddl==null || c.ddl.filledBuckets()<1);
				if(fresh || rec.cardinality>=c.ddl.cardinality()){
					if(fresh){attached++;}
					c.ddl=rec.ddl;
				}
			}
		}
		System.err.println("Attached "+attached+" sketches to clades.");
	}

	public void buildSketchIndex(){
		if(sketchRecords==null || sketchRecords.isEmpty()){return;}
		long t0=System.nanoTime();
		ddlIndex=new DDLIndex();
		int it=indexThreads>0 ? indexThreads : Math.min(Shared.threads(), 32);
		ddlIndex.addAll(sketchRecords, it);
		long elapsed=System.nanoTime()-t0;
		System.err.println("Built sketch index with "+sketchRecords.size()+" entries, "+
				ddlIndex.populatedCells()+" populated cells in "+
				String.format("%.3f", elapsed/1e9)+" seconds.");
	}

	public void addSketchInfo(ArrayList<Comparison> results, Clade query){
		if(ddlIndex==null || query.ddl==null){return;}
		int[][] topSketch=ddlIndex.topHits(query.ddl, maxSketchHits);
		if(topSketch==null || topSketch.length==0){return;}

		int bestIdx=topSketch[0][0];
		int bestMatches=topSketch[0][1];
		if(bestMatches<minSketchMatches){return;}

		DDLRecord bestRec=sketchRecords.get(bestIdx);
		TaxTree tree=TaxTree.getTree();

		Clade bestSketchClade=(cladeMap!=null && bestRec.taxID>0 ? cladeMap.get(bestRec.taxID) : null);

		for(Comparison comp : results){
			if(comp.ref==null){continue;}
			comp.sketchTaxID=bestRec.taxID;
			comp.sketchName=bestRec.name;
			comp.sketchMatches=bestMatches;
			if(tree!=null && bestRec.taxID>0 && comp.ref.taxID>0){
				comp.sketchLCA=tree.commonAncestorLevel(comp.ref.taxID, bestRec.taxID);
			}else if(bestSketchClade!=null){
				comp.sketchLCA=lineageLCA(comp.ref.lineage(), bestSketchClade.lineage());
			}
		}

		if(cladeMap!=null){
			for(int s=0; s<topSketch.length; s++){
				int idx=topSketch[s][0];
				int matches=topSketch[s][1];
				if(matches<minSketchMatches){break;}
				DDLRecord rec=sketchRecords.get(idx);
				if(rec.taxID<=0){continue;}

				boolean found=false;
				for(Comparison comp : results){
					if(comp.ref!=null && comp.ref.taxID==rec.taxID){found=true; break;}
				}

				if(!found){
					Clade refClade=cladeMap.get(rec.taxID);
					if(refClade!=null){
						Comparison sketchComp=new Comparison(query, refClade);
						sketchComp.sketchTaxID=rec.taxID;
						sketchComp.sketchName=rec.name;
						sketchComp.sketchMatches=matches;
						if(tree!=null){
							sketchComp.sketchLCA=tree.commonAncestorLevel(refClade.taxID, rec.taxID);
						}
						if(Clade.MAKE_DDLS){sketchComp.compareDDL();}
						sketchComp.isSketchHit=true;
						results.add(sketchComp);
					}
				}
			}
		}
	}

	private static final String[] LCA_PREFIXES={"s__","g__","f__","o__","c__","p__","k__","sk__"};
	private static final int[] LCA_LEVELS={
		TaxTree.SPECIES, TaxTree.GENUS, TaxTree.FAMILY, TaxTree.ORDER,
		TaxTree.CLASS, TaxTree.PHYLUM, TaxTree.KINGDOM, TaxTree.SUPERKINGDOM
	};

	static int lineageLCA(CharSequence lineageA, CharSequence lineageB){
		if(lineageA==null || lineageB==null){return -1;}
		String a=lineageA.toString(), b=lineageB.toString();
		for(int i=0; i<LCA_PREFIXES.length; i++){
			String prefix=LCA_PREFIXES[i];
			int posA=a.indexOf(prefix);
			int posB=b.indexOf(prefix);
			if(posA<0 || posB<0){continue;}
			int endA=a.indexOf(';', posA);
			int endB=b.indexOf(';', posB);
			String nameA=a.substring(posA+prefix.length(), endA<0 ? a.length() : endA);
			String nameB=b.substring(posB+prefix.length(), endB<0 ? b.length() : endB);
			if(nameA.equals(nameB)){return LCA_LEVELS[i];}
		}
		return -1;
	}

}