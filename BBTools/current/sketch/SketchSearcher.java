package sketch;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import shared.Parse;
import shared.Shared;
import shared.Tools;
import structures.AbstractBitSet;
import structures.ByteBuilder;
import structures.Heap;
import structures.IntHashMap;
import tax.TaxNode;
import tax.TaxTree;

public class SketchSearcher extends SketchObject {
	
	public SketchSearcher(){
		
	}

	public boolean parse(String arg, String a, String b, boolean addFileIfNotFound){
		
//		System.err.println("Parsing "+arg+"; ref="+refFiles); //123
		
		if(parseSketchFlags(arg, a, b)){
			//Do nothing
		}else if(defaultParams.parse(arg, a, b)){
			//Do nothing
		}else if(a.equals("verbose")){
			verbose=Parse.parseBoolean(b);
		}else if(a.equals("ref")){
			addRefFiles(b);
		}else if(arg.equalsIgnoreCase("nt") || arg.equalsIgnoreCase("RefSeq") || arg.equalsIgnoreCase("refseqbig") || arg.equalsIgnoreCase("nr")
				|| arg.equalsIgnoreCase("img") || arg.equalsIgnoreCase("silva") || arg.equalsIgnoreCase("ribo")
				 || arg.equalsIgnoreCase("mito") || arg.equalsIgnoreCase("fungi") 
				 || arg.equalsIgnoreCase("prokprot") || arg.equalsIgnoreCase("prokprotbig") || arg.equalsIgnoreCase("protein") || 
				 arg.equalsIgnoreCase("protien") || a.equalsIgnoreCase("prot")){
			addRefFiles(arg);
		}else if(a.equals("threads") || a.equals("sketchthreads") || a.equals("t")){
			threads=Integer.parseInt(b);
		}
		
		else if(a.equalsIgnoreCase("minLevelExtended") || a.equalsIgnoreCase("minLevel")){
			minLevelExtended=TaxTree.parseLevelExtended(b);
		}else if(a.equals("index") || a.equals("makeindex")){
			if(b!=null && "auto".equalsIgnoreCase(b)){
				autoIndex=true;
				makeIndex=true;
			}else{
				autoIndex=false;
				makeIndex=Parse.parseBoolean(b);
			}
		}else if(a.equals("indexsize") || a.equals("indexlimit")){
			SketchIndex.indexLimit=Integer.parseInt(b);
		}
		
		else if(b==null && arg.indexOf('=')<0 && addFileIfNotFound && (arg.indexOf(',')>=0 || new File(arg).exists())){
			addRefFiles(arg);
		}else{
			return false;
		}
//		System.err.println("Parsed "+arg+"; ref="+refFiles); //123
		return true;
	}

	public boolean compare(ArrayList<Sketch> querySketches, ByteBuilder sb, DisplayParams params, int maxThreads){
		assert(params.postParsed);
		final boolean json=params.json();
		ConcurrentHashMap<Integer, Comparison> map=new ConcurrentHashMap<Integer, Comparison>();
		
		SketchResults[] alca=new SketchResults[querySketches.size()];

		if(verbose2){System.err.println("At compare.");}
		
		boolean success=true;
		final CompareBuffer buffer=new CompareBuffer(false);
		AtomicInteger fakeID=new AtomicInteger(minFakeID);
		for(int i=0; i<querySketches.size(); i++){
			fakeID.set(minFakeID);
			Sketch a=querySketches.get(i);
			
			SketchResults results=processSketch(a, buffer, fakeID, map, params, maxThreads);
			a.clearRefHitCounts();
			alca[i]=results;
//			System.out.println(a.present);
		}
		
		if(verbose2){System.err.println("Made results.");}
		
		for(int i=0; i<alca.length; i++){
//			Sketch s=sketches.get(i);
			SketchResults results=alca[i];

			if(json && alca.length>1 && i==0){
				sb.append('[');
			}

			sb.append(results.toText(params));

			if(json && alca.length>1){
				if(i<alca.length-1){
					sb.append(',');
				}else{
					sb.append(']');
				}
			}
		}
		return success;
	}
	
	private class CompareThread extends Thread {
		
		CompareThread(Sketch a_, ArrayList<Sketch> localRefSketches_, int pid_, int incr_,
				AtomicInteger fakeID_, ConcurrentHashMap<Integer, Comparison> map_, DisplayParams params_){
			a=a_;
			pid=pid_;
			incr=incr_;
			fakeID=fakeID_;
			map=map_;
			params=params_;
			localRefSketches=localRefSketches_;
			buffer=new CompareBuffer(params.needContamCounts());
			if(buffer.cbs!=null){buffer.cbs.setCapacity(a.length(), 0);}
		}
		
		@Override
		public void run(){
			if(a.length()<1 || a.length()<params.minHits || (params.requireSSU && !a.hasSSU())){return;}//TODO: Change to 'require16S'
			assert(a.compareBitSet()==null || buffer.cbs!=null) : (a.compareBitSet()==null)+", "+(buffer.cbs==null); //Unsafe to use a.cbs multithreaded unless atomic
			final AbstractBitSet cbs=(buffer.cbs==null ? a.compareBitSet() : buffer.cbs);
			for(int i=pid; i<localRefSketches.size(); i+=incr){
				Sketch b=localRefSketches.get(i);
				if(params.passesFilter(b)){
					processPair(a, b, buffer, cbs, fakeID, map, params);
					localComparisons++;
				}
			}
			comparisons.getAndAdd(localComparisons);
		}
		
		final AtomicInteger fakeID;
		final ConcurrentHashMap<Integer, Comparison> map;
		final CompareBuffer buffer;
		final int incr;
		final int pid;
		final Sketch a;
		final DisplayParams params;
		final ArrayList<Sketch> localRefSketches;
		long localComparisons=0;
		
	}
	
	public SketchResults processSketch(Sketch a, CompareBuffer buffer, AtomicInteger fakeID, 
			ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, int maxThreads){
		if(a.length()<1 || a.length()<params.minHits || (params.requireSSU && !a.hasSSU())){return new SketchResults(a);}
		//		Timer t=new Timer();
		//		t.start("Began query.");
		assert(a.compareBitSet()==null);
		assert(a.indexBitSet()==null);
		
		if(verbose2){System.err.println("At processSketch 1");} //123
		
		a.makeBitSets(params.needContamCounts(), index!=null);
		
		final SketchResults sr;
		if(index!=null){
			sr=index.getSketches(a, params);
		}else{
			sr=new SketchResults(a, refSketches, null);
		}
		
		if(verbose2){System.err.println("At processSketch 2");} //123
		
		if(sr==null || sr.refSketchList==null || sr.refSketchList.isEmpty()){
			if(verbose2){System.err.println("At processSketch 2.0");} //123
			return sr;
		}
		
		if(verbose2){System.err.println("At processSketch 2.1");} //123
		
		if(verbose2){System.err.println("At processSketch 2.2");} //123
		
		if(maxThreads>1 && Shared.threads()>1 && sr.refSketchList.size()>31){
			if(verbose2){System.err.println("At processSketch 2.3");} //123
			assert((buffer.cbs==null)==(params.needContamCounts()));
			spawnThreads(a, sr.refSketchList, fakeID, map, params, maxThreads);
			if(verbose2){System.err.println("At processSketch 2.4");} //123
		}else{
			if(verbose2){System.err.println("At processSketch 2.5");} //123
			assert(buffer.cbs==null);
			long comp=0;
			for(Sketch b : sr.refSketchList){
				if(params.passesFilter(b)){
					comp++;
					processPair(a, b, buffer, a.compareBitSet(), /*sr.taxHits,*/ fakeID, map, params);
				}
			}
			comparisons.getAndAdd(comp);
			if(verbose2){System.err.println("At processSketch 2.6");} //123
		}
		if(verbose2){System.err.println("At processSketch 3");} //123
		
		sr.addMap(map, params, buffer);
		
		fakeID.set(minFakeID);
		map.clear();
		if(verbose2){System.err.println("At processSketch 4");} //123
		a.clearRefHitCounts();
		
		return sr;
	}
	
	//For remote homology
	boolean passesTax(Sketch q, Sketch ref){
		assert(minLevelExtended>=0);
		final int qid=q.taxID;
		if(qid<0 || qid>=minFakeID){return false;}
		TaxNode qtn=taxtree.getNode(qid);
		if(qtn==null){return false;}
		if(qtn.levelExtended>minLevelExtended){return false;}
		final int rid=(ref==null ? -1 : ref.taxID);
		if(rid>=0 && rid<minFakeID){
			TaxNode rtn=taxtree.getNode(rid);
			if(rtn!=null && rtn.levelExtended<=minLevelExtended){
				TaxNode ancestor=taxtree.commonAncestor(qtn, rtn);
				if(ancestor!=null && ancestor.levelExtended>=minLevelExtended){
					return true;
				}
			}
		}
		return false;
	}

	private void spawnThreads(Sketch a, ArrayList<Sketch> refs, AtomicInteger fakeID,
			ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, int maxThreads){
		final int toSpawn=Tools.max(1, Tools.min((refs.size()+7)/8, threads, maxThreads, Shared.threads()));
		ArrayList<CompareThread> alct=new ArrayList<CompareThread>(toSpawn);
		if(verbose2){System.err.println("At spawnThreads");} //123
		for(int t=0; t<toSpawn; t++){
			alct.add(new CompareThread(a, refs, t, toSpawn, fakeID, map, params));
		}
		for(CompareThread ct : alct){ct.start();}
		for(CompareThread ct : alct){

			//Wait until this thread has terminated
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					ct.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
		if(params.needContamCounts()){
			for(CompareThread ct : alct){
				if(ct.buffer.cbs==null){
					assert((AUTOSIZE || AUTOSIZE_LINEAR) && index!=null);//Not really what this does
					break;
				}
				a.addToBitSet(ct.buffer.cbs);
			}
		}
		a.clearRefHitCounts();
		alct=null;
	}
	
//	private void writeResults(ArrayList<Comparison> al, Sketch s, StringBuilder sb){
//		sb.append("\nResults for "+s.name()+":\n\n");
//
//		ArrayList<TaxNode> tnl=new ArrayList<TaxNode>();
//		for(Comparison c : al){
//			formatComparison(c, format, sb, printTax);
//		}
//	}
	
	boolean processPair(Sketch a, Sketch b, CompareBuffer buffer, AbstractBitSet abs,
			AtomicInteger fakeID, ConcurrentHashMap<Integer, Comparison> map, DisplayParams params){
//		System.err.println("Comparing "+a.name()+" and "+b.name());
		assert(!params.printRefHits || a.refHitCounts()!=null || !SketchObject.makeIndex);
		
		
		if(b.genomeSizeBases<params.minBases){return false;}
		if(minLevelExtended>-1 && !passesTax(a, b)){return false;}
		if(params.minSizeRatio>0){
			long sea=a.genomeSizeEstimate();
			long seb=b.genomeSizeEstimate();
			if(Tools.min(sea, seb)<params.minSizeRatio*Tools.max(sea, seb)){return false;}
		}
		Comparison c=compareOneToOne(a, b, buffer, abs, /*taxHits, params.contamLevel(),*/ params.minHits, params.minWKID, params.minANI, params.requireSSU, null);
		if(c==null){return false;}
		if(c.taxID()<1){c.taxID=fakeID.getAndIncrement();}
		
//		System.err.println("TID: "+c.taxID()+", "+fakeID);
		
		TaxNode tn=(taxtree==null ? null : taxtree.getNode(b.taxID));
		if(tn!=null){
			c.taxName=tn.name;
			if(tn.level<params.taxLevel){
				TaxNode tn2=taxtree.getNodeAtLevel(b.taxID, params.taxLevel);
				tn=tn2;
			}
		}
		Integer key=(tn==null ? c.taxID : tn.id);
		
		Comparison old=map.get(key);
//		System.err.println("A. Old: "+(old==null ? 0 : old.hits)+", new: "+c.hits);
		if(old!=null && params.compare(old, c)>0){return false;}
		
		old=map.put(key, c);
		while(old!=null && params.compare(old, c)>0){
//			System.err.println("B. Old: "+(old==null ? 0 : old.hits)+", new: "+c.hits);
			c=old;
			old=map.put(key, c);
		}
		return true;
	}
	
//	//TODO:  Interestingly, the heap never seems to be created by anything...  not sure what it's for.
//	private static Comparison compareOneToOne(final Sketch a, final Sketch b, CompareBuffer buffer, AbstractBitSet abs,
//			int minHits, float minWKID, float minANI, boolean aniFromWKID, Heap<Comparison> heap){
////		assert(heap!=null); //Optional, for testing.
//		if(a==b && !compareSelf){return null;}
//		final int matches=a.countMatches(b, buffer, abs, true/*!makeIndex || !AUTOSIZE*/, null, -1);
//		assert(matches==buffer.hits());
//		if(matches<minHits){return null;}
////		asdf //TODO: handle k1 and k2 WKIDs here.
//		{
////			final int div=aniFromWKID ? buffer.minDivisor() : buffer.maxDivisor();
////			final float xkid=matches/(float)div;//This could be kid or wkid at this point...
////			if(xkid<minWKID){return null;}
//			
//			final int div=aniFromWKID ? buffer.minDivisor() : buffer.maxDivisor();
//			final float xkid=matches/(float)div;//This could be kid or wkid at this point...
//			if(xkid<minWKID){return null;}
//			
//			//TODO (?)  This is only necessary because of the order of setting minwkid and minani.
//			//minWKID can be deterministically determined from minANI so if it is set correctly this can be skipped.
//			if(minANI>0){
//				final float ani=wkidToAni(xkid, a.k1Fraction());
//				if(ani<minANI){return null;}
//			}
//		}
//		
//		if(heap!=null && !heap.hasRoom() && heap.peek().hits()>matches){return null;} //TODO:  Should be based on score
//		
////		System.err.print("*");
//		Comparison c=new Comparison(buffer, a, b);
//		if(heap==null || heap.add(c)){return c;}
//		return null;
//	}
	
	//TODO:  Interestingly, the heap never seems to be created by anything...  not sure what it's for.
	private static Comparison compareOneToOne(final Sketch a, final Sketch b, CompareBuffer buffer, AbstractBitSet abs,
			int minHits, float minWKID, float minANI, boolean requireSSU, Heap<Comparison> heap){
//		assert(heap!=null); //Optional, for testing.
//		assert(a.refHitCounts!=null);
		if(a==b && !compareSelf){return null;}
		if(requireSSU && !a.sharesSSU(b)){return null;}
		final int matches=a.countMatches(b, buffer, abs, true/*!makeIndex || !AUTOSIZE*/, null, -1);
		assert(matches==buffer.hits());
		if(matches<minHits){return null;}
		
		{
			final float wkid=buffer.wkid();
			if(wkid<minWKID){return null;}
			
			if(minANI>0){
				final float ani=buffer.ani();
				if(ani<minANI){return null;}
			}
		}
		
		if(heap!=null && !heap.hasRoom() && heap.peek().hits()>matches){return null;} //TODO:  Should be based on score
		
//		System.err.print("*");
		Comparison c=new Comparison(buffer, a, b);
		if(heap==null || heap.add(c)){return c;}
		return null;
	}
	
	public void addRefFiles(String a){
		if(a.equalsIgnoreCase("nr")){
			addRefFiles(NR_PATH());
			if(blacklist==null){blacklist=Blacklist.nrBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="nr";}
			if(!setK){k=defaultKAmino; k2=defaultK2Amino;}
		}else if(a.equalsIgnoreCase("nt")){
			addRefFiles(NT_PATH());
			if(blacklist==null){blacklist=Blacklist.ntBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="nt";}
			if(!setK){k=defaultK; k2=defaultK2;}
		}else if(a.equalsIgnoreCase("refseq")){
			addRefFiles(REFSEQ_PATH());
			if(blacklist==null){blacklist=Blacklist.refseqBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="RefSeq";}
			if(!setK){k=defaultK; k2=defaultK2;}
			if(!SET_AUTOSIZE_FACTOR){AUTOSIZE_FACTOR=2.0f;}
		}else if(a.equalsIgnoreCase("refseqbig")){
			addRefFiles(REFSEQ_PATH_BIG());
			if(blacklist==null){blacklist=Blacklist.refseqBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="RefSeq";}
			if(!setK){k=defaultK; k2=defaultK2;}
			if(!SET_AUTOSIZE_FACTOR){AUTOSIZE_FACTOR=4.5f;}
		}else if(a.equalsIgnoreCase("silva")){
//			TaxTree.SILVA_MODE=Parse.parseBoolean(b);
			addRefFiles(SILVA_PATH());
			if(blacklist==null){blacklist=Blacklist.silvaBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="Silva";}
			if(!setK){k=defaultK; k2=defaultK2;}
		}else if(a.equalsIgnoreCase("img")){
			addRefFiles(IMG_PATH());
			if(blacklist==null){blacklist=Blacklist.imgBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="IMG";}
			if(!setK){k=defaultK; k2=defaultK2;}
		}else if(a.equalsIgnoreCase("prokprot") || a.equalsIgnoreCase("protein")){
			addRefFiles(PROKPROT_PATH());
			if(blacklist==null){blacklist=Blacklist.prokProtBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="ProkProt";}
			if(!setK){k=defaultKAmino; k2=defaultK2Amino;}
			if(!amino && !translate) {
				translate=true;
				System.err.println("Setting translate to true because a protein dataset is being used.");
			}
			if(!SET_AUTOSIZE_FACTOR){AUTOSIZE_FACTOR=3.0f;}
		}else if(a.equalsIgnoreCase("prokprotbig") || a.equalsIgnoreCase("proteinbig")){
			addRefFiles(PROKPROT_PATH_BIG());
			if(blacklist==null){blacklist=Blacklist.prokProtBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="ProkProt";}
			if(!setK){k=defaultKAmino; k2=defaultK2Amino;}
			if(!amino && !translate) {
				translate=true;
				System.err.println("Setting translate to true because a protein dataset is being used.");
			}
			if(!SET_AUTOSIZE_FACTOR){AUTOSIZE_FACTOR=7.5f;}
		}else if(a.equalsIgnoreCase("mito") || a.equalsIgnoreCase("refseqmito")){
			addRefFiles(MITO_PATH());
			if(blacklist==null){blacklist=Blacklist.mitoBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="RefSeqMito";}
			if(!setK){k=defaultK; k2=defaultK2;}
		}else if(a.equalsIgnoreCase("fungi") || a.equalsIgnoreCase("refseqfungi")){
			addRefFiles(FUNGI_PATH());
			if(blacklist==null){blacklist=Blacklist.fungiBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="RefSeqFungi";}
			if(!setK){k=defaultK; k2=defaultK2;}
		}else{
			addFiles(a, refFiles);
		}
	}
	
	static void addFiles(String a, Set<String> list){
		if(a==null){return;}
		File f=new File(a);
		assert(!list.contains(a)) : "Duplicate file "+a;
		
		if(f.exists()){
			list.add(a);
		}else if(a.indexOf(',')>0){
			for(String s : a.split(",")){addFiles(s, list);}
		}else if(a.indexOf('#')>=0 && new File(a.replaceFirst("#", "0")).exists()){
			for(int i=0; true; i++){
				String temp=a.replaceFirst("#", ""+i);
				if(!new File(temp).exists()){break;}
				list.add(temp);
			}
		}else{
			list.add(a);
		}
	}
	
	public void makeIndex(){
		assert(index==null);
		index=new SketchIndex(refSketches);
		index.load();
	}
	
	public void loadReferences(int mode_, DisplayParams params){
		loadReferences(mode_, params.minKeyOccuranceCount, params.minEntropy, params.minProb, params.minQual);
	}
	
	public void loadReferences(int mode_, int minKeyOccuranceCount, float minEntropy, float minProb, byte minQual) {
		makeTool(minKeyOccuranceCount, false, false);
		refSketches=tool.loadSketches_MT(mode_, 1f, -1, minEntropy, minProb, minQual, refFiles);
		assert(refSketches!=null) : refFiles;
		if(mode_==PER_FILE){
			Collections.sort(refSketches, SketchIdComparator.comparator);
		}
		taxIDToSketchIDMap=new IntHashMap(Tools.max(3, (int)(refSketches.size()*1.2f)));
		for(int i=0; i<refSketches.size(); i++){
			Sketch sk=refSketches.get(i);
			if(sk!=null && sk.taxID>0){
				taxIDToSketchIDMap.set(sk.taxID, i);
			}
		}
//		System.err.println("Sketches: "+refSketches.get(0).name());
		if(makeIndex){
			makeIndex();
		}
	}
	
	public void makeTool(int minKeyOccuranceCount, boolean trackCounts, boolean mergePairs){
		if(tool==null){
			tool=new SketchTool(targetSketchSize, minKeyOccuranceCount, trackCounts, mergePairs, rcomp);
		}
	}
	
	public ArrayList<Sketch> loadSketchesFromString(String sketchString){
		return tool.loadSketchesFromString(sketchString);
	}
	
	public int refFileCount(){return refFiles==null ? 0 : refFiles.size();}
	public int refSketchCount(){return refSketches==null ? 0 : refSketches.size();}
	
	public Sketch findReferenceSketch(int taxID){
		if(taxID<1){return null;}
		int skid=taxIDToSketchIDMap.get(taxID);
		return skid<0 ? null : refSketches.get(skid);
	}
	
	/*--------------------------------------------------------------*/
	
	public SketchIndex index=null;
	public boolean autoIndex=true;
	
	public SketchTool tool=null;
	public ArrayList<Sketch> refSketches;
	LinkedHashSet<String> refFiles=new LinkedHashSet<String>();
	/** For ref sketch lookups by TaxID */
	private IntHashMap taxIDToSketchIDMap;
	public int threads=Shared.threads();
	boolean verbose;
	boolean errorState=false;
	AtomicLong comparisons=new AtomicLong(0);
	
	int minLevelExtended=-1;
	
}
