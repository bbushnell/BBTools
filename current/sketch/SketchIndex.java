package sketch;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import kmer.AbstractKmerTable;
import kmer.HashBuffer;
import kmer.KmerTableSet;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.AbstractBitSet;
import structures.IntHashMap;
import structures.IntHashSetList;
import structures.IntList;

public class SketchIndex extends SketchObject {
	
	public SketchIndex(ArrayList<Sketch> refs){
		refSketches=refs;
		tables=new KmerTableSet(new String[] {"ways="+WAYS, "tabletype="+AbstractKmerTable.ARRAYHF, "prealloc="+(prealloc>0 ? ""+prealloc : "f")}, 
				20+(defaultParams.trackCounts() ? 4 : 0)+18);//An extra 18 bytes per kmer because a lot of them occur multiple times
		tables.allocateTables();
		tableArray=tables.tables();
	}
	
	public void load(){
		spawnIndexThreads();
		if(useWhitelist){
			assert(!Whitelist.exists());
			Whitelist.initialize(tableArray);
		}
	}
	
	/** Spawn index threads */
	private void spawnIndexThreads(){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		ArrayList<IndexThread> alht=new ArrayList<IndexThread>(threads);
		AtomicInteger ai=new AtomicInteger(0);
		AtomicLong totalKeys=new AtomicLong(0);
		AtomicLong uniqueKeys=new AtomicLong(0);
		for(int i=0; i<threads; i++){
			alht.add(new IndexThread(ai, totalKeys, uniqueKeys));
		}
		
		//Start the threads
		for(IndexThread pt : alht){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		long codesProcessed=0;
		for(IndexThread pt : alht){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
					synchronized(pt){
						codesProcessed+=pt.codesProcessedT;
					}
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		System.err.println("Indexed "+uniqueKeys+" unique and "+totalKeys+" total hashcodes."); //For some reason codesProcessed is nondeterministic.
		
		//Do anything necessary after processing
//		System.gc();
	}
	
	/*--------------------------------------------------------------*/
	
	public SketchResults getSketches(Sketch a, DisplayParams params){
		if(useIntMap){
			return getSketchesMap(a, params);
		}else{
			return getSketchesList(a, params);
		}
	}
	
	/** Return true if added. */
	private boolean addToTaxSet(int sketchID, IntHashSetList taxSet, int taxLevelExtended){
		Sketch sk=refSketches.get(sketchID);
		int taxID=sk.taxID;
		if(taxID<0 || taxID>=minFakeID){return false;}
		taxID=taxtree.getIdAtLevelExtended(taxID, taxLevelExtended);
		return taxSet.add(taxID);
	}
	
	public SketchResults getSketchesList(final Sketch a, DisplayParams params){
		final int minHits=params.minHits, contamLevel=params.contamLevel();
		final boolean countContamHits=params.needContamCounts();//, metaFilter=params.hasMetaFilters(), taxFilter=params.hasTaxFilters();
		final Timer t=(printTime ? new Timer() : null);
		
		final int[] singleton=new int[1];
		final IntList idList=new IntList(Tools.min(targetSketchSize, indexLimit, 1000));
		AbstractBitSet abs=a.indexBitSet();
		assert((abs==null)!=countContamHits);
		
		final IntHashSetList taxSet;
		final int[][] taxHits;
		if(contamLevel>=0){
			taxSet=new IntHashSetList(31);
			taxHits=new int[a.length()][];
			assert(taxtree!=null) : "A TaxTree is required for this operation.";
		}else{
			taxSet=null;
			taxHits=null;
		}
		
		for(int i=0; i<a.keys.length; i++){
			long key=a.keys[i];
			AbstractKmerTable set=tableArray[(int)(key%WAYS)];
//			System.err.println(set.getValue(key));
			final int[] ids=set.getValues(key, singleton);
//			System.err.println(Arrays.toString(ids));
			if(ids!=null && ids[0]>=0){
				int incr=0;
				for(int id : ids){
					if(id>=0){
						final int trueID=id-1;//Minimum id is 1, indicating sketch 0.
						idList.add(trueID);//Minimum id is 1, indicating sketch 0.
						incr++;
						if(taxSet!=null){addToTaxSet(trueID, taxSet, contamLevel);}
					}
				}
				if(countContamHits && incr>0){abs.increment(i, incr);}
				if(taxSet!=null && taxSet.size()>0){
					taxHits[i]=taxSet.toArray();
					taxSet.clear();
				}
			}
			
		}
		
//		assert(abs!=null);
//		assert(false) : abs.cardinality();
		
		if(printTime){
			t.stop("\nTime for searching index: \t");
			t.start();
		}
		
//		System.err.println("idList.size:"+idList.size);
		if(idList.size<minHits){return new SketchResults(a);}//null breaks some things
		idList.sort();
		
		if(printTime){
			t.stop("Time for sorting "+idList.size()+" hits:\t");
			t.start();
		}
		
		ArrayList<Sketch> list=new ArrayList<Sketch>(Tools.min(8, idList.size));

		int last=-1;
		int hits=0;
		for(int i=0; i<idList.size; i++){
			int id=idList.get(i);
			if(id==last){
//				System.err.println("A: "+last+", "+id+", "+count+", "+minHits);
				hits++;
			}else{
//				System.err.println("B: "+last+", "+id+", "+count+", "+minHits);
				if(last>-1 && (hits>=minHits)){
					final Sketch ref=refSketches.get(last);
					list.add(ref);
				}
				last=id;
				hits=0;
			}
		}
		if(last>-1 && (hits>=minHits)){
			final Sketch ref=refSketches.get(last);
//			if((!metaFilter || ref.passesMeta(params)) && (!taxFilter || params.passesFilter(ref))){list.add(ref);}
			list.add(ref);
		}
		if(printTime){
			t.stop("Time for fetching sketches: \t");
		}
		return list.isEmpty() ? new SketchResults(a) : new SketchResults(a, list, taxHits);
	}
	
//	static ThreadLocal<IntHashMap> intMapHolder=new ThreadLocal<IntHashMap>();
	
	public final int[] getSketchIdsMap(long key, int[] singleton){
		AbstractKmerTable set=tableArray[(int)(key%WAYS)];
		final int[] ids=set.getValues(key, singleton);
		return ids;
	}
	
	public SketchResults getSketchesMap(final Sketch a, DisplayParams params){
		final int minHits=params.minHits, contamLevel=params.contamLevel();
		final boolean countContamHits=params.needContamCounts();//, metaFilter=params.hasMetaFilters(), taxFilter=params.hasTaxFilters();
		final Timer t=(printTime ? new Timer() : null);
		final int[] singleton=new int[1];
		
		final IntHashMap idMap=new IntHashMap(Tools.min(targetSketchSize, indexLimit, intMapSize), 0.7f);
		
		AbstractBitSet abs=a.indexBitSet();
		assert((abs==null)!=countContamHits);
		
		final IntHashSetList taxSet;
		final int[][] taxHits;
		if(contamLevel>=0){
			taxSet=new IntHashSetList(31);
			taxHits=new int[a.length()][];
			assert(taxtree!=null) : "A TaxTree is required for this operation.";
		}else{
			taxSet=null;
			taxHits=null;
		}
//		assert(false) : (taxHits==null)+", "+contamLevel;
		
		if(printTime){
			t.stop("\nTime for allocation:      \t");
			t.start();
		}
		
		int[] refHitCounts;
		if(params.printRefHits){
			refHitCounts=new int[a.keys.length];
			a.setRefHitCounts(refHitCounts);
		}else{refHitCounts=null;}
		
		for(int i=0; i<a.keys.length; i++){
			long key=a.keys[i];
			final int[] ids=getSketchIdsMap(key, singleton);
			
			if(ids!=null && ids[0]>=0){
				int incr=0;
				for(int id : ids){
					if(id>=0){
						final int trueID=id-1;//Minimum id is 1, indicating sketch 0.
						idMap.increment(trueID);
						if(!allToAll || compareSelf){incr++;}
						if(taxSet!=null){addToTaxSet(trueID, taxSet, contamLevel);}
						if(refHitCounts!=null && trueID!=a.sketchID){refHitCounts[i]++;}
					}
				}
				if(countContamHits && incr>0){abs.increment(i, incr);}
				if(taxSet!=null && taxSet.size()>0){
					taxHits[i]=taxSet.toArray();
					taxSet.clear();
				}
			}
		}
		
		if(printTime){
			t.stop("Time for searching index: \t");
			System.err.println("Size:                   \t"+idMap.size());
			t.start();
		}
		
//		System.err.println("idList.size:"+idList.size);
		final int size=idMap.size();
		if(size==0){return new SketchResults(a);}//null breaks some things
		
		ArrayList<Sketch> list=new ArrayList<Sketch>(Tools.min(8, size));

		final int[] keys=idMap.keys();
		final int[] values=idMap.values();
		
		for(int i=0; i<keys.length; i++){
			int value=values[i];
			if(value>=minHits){
				int id=keys[i];
				final Sketch ref=refSketches.get(id);
//				if((!metaFilter || ref.passesMeta(params)) && (!taxFilter || params.passesFilter(ref))){list.add(ref);}
				list.add(ref);
			}
		}
		if(printTime){
			t.stop("Time for fetching sketches: \t");
		}
		return list.isEmpty() ? new SketchResults(a) : new SketchResults(a, list, taxHits);
	}
	
	/*--------------------------------------------------------------*/
	
	public class IndexThread extends Thread {
		
		public IndexThread(AtomicInteger nextIndex_, AtomicLong keyCount_, AtomicLong uniqueKeyCount_){
			buffer=new HashBuffer(tableArray, 1000, 31, true, false);
			nextIndex=nextIndex_;
			keyCount=keyCount_;
			uniqueKeyCount=uniqueKeyCount_;
		}
		
		@Override
		public void run(){
//			System.err.println("Thread running.");
			int id=nextIndex.getAndIncrement();
			final int numSketches=refSketches.size();
			final int limit0=Tools.min((AUTOSIZE || AUTOSIZE_LINEAR ? Integer.MAX_VALUE : targetSketchSize), indexLimit);
//			System.err.println("numSketches="+numSketches);
			while(id<numSketches){
				final Sketch sk=refSketches.get(id);
				final long[] array=sk.keys;
				final int limit=Tools.min(array.length, limit0);
//				System.err.println("limit="+limit);
				for(int i=0; i<limit; i++){
					long key=array[i];
					buffer.set(key, id+1);//Must be greater than zero
					codesProcessedT++;
				}
				id=nextIndex.getAndIncrement();
			}
			long temp=buffer.flush();
			
			synchronized(this){
				codesProcessedT+=0;
				success=true;
				keyCount.getAndAdd(codesProcessedT);
				uniqueKeyCount.getAndAdd(buffer.uniqueAdded);
//				if(codesProcessedT>0){System.err.println(codesProcessedT);}
			}
		}
		
		AtomicInteger nextIndex;
		AtomicLong keyCount;
		AtomicLong uniqueKeyCount;
		long codesProcessedT=0;
		HashBuffer buffer;
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	
	public final KmerTableSet tables;
	public final AbstractKmerTable[] tableArray;
	public final ArrayList<Sketch> refSketches;
	
	public boolean errorState=false;

	private static final boolean printTime=false;
	public static boolean useIntMap=true;
//	public static boolean useIntMapBinary=false;
	public static int intMapSize=1000;
	public static int indexLimit=Integer.MAX_VALUE;
	public static final int WAYS=31;
	
}
