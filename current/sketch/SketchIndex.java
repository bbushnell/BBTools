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

/**
 * Index structure for fast k-mer lookup across multiple reference sketches.
 * Provides efficient search capabilities using hash tables distributed across multiple ways.
 * Supports both list-based and map-based search algorithms for different performance needs.
 * Used for rapid similarity searches in large sketch databases.
 *
 * @author Brian Bushnell
 */
public class SketchIndex extends SketchObject {
	
	/**
	 * Constructs a SketchIndex with reference sketches and allocates hash tables.
	 * Creates KmerTableSet with specified ways and preallocated space for efficient indexing.
	 * @param refs ArrayList of reference sketches to index
	 */
	public SketchIndex(ArrayList<Sketch> refs){
		refSketches=refs;
		tables=new KmerTableSet(new String[] {"ways="+WAYS, "tabletype="+AbstractKmerTable.ARRAYHF, "prealloc="+(prealloc>0 ? ""+prealloc : "f")}, 
				20+(defaultParams.trackCounts() ? 4 : 0)+18);//An extra 18 bytes per kmer because a lot of them occur multiple times
		tables.allocateTables();
		tableArray=tables.tables();
	}
	
	/** Loads the index by spawning indexing threads and optionally initializing whitelist.
	 * Populates hash tables with k-mers from all reference sketches using multithreading. */
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
	
	/**
	 * Main entry point for sketch similarity search.
	 * Delegates to either map-based or list-based search algorithm based on useIntMap setting.
	 *
	 * @param a Query sketch to search for
	 * @param params Display and filtering parameters
	 * @return SketchResults containing matching reference sketches
	 */
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
	
	/**
	 * List-based sketch search algorithm using IntList for hit accumulation.
	 * Searches index for k-mer matches, accumulates hits per sketch, and filters by minimum hits.
	 * Supports contamination level tracking and taxonomic filtering when enabled.
	 *
	 * @param a Query sketch to search for
	 * @param params Display parameters including minimum hits threshold
	 * @return SketchResults with matching sketches above hit threshold
	 */
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
	
	/**
	 * Retrieves sketch IDs that contain the specified k-mer key.
	 * Uses hash table lookup with key modulo WAYS for distribution.
	 *
	 * @param key K-mer hash code to search for
	 * @param singleton Reusable array for single value results
	 * @return Array of sketch IDs containing the key, or null if not found
	 */
	public final int[] getSketchIdsMap(long key, int[] singleton){
		AbstractKmerTable set=tableArray[(int)(key%WAYS)];
		final int[] ids=set.getValues(key, singleton);
		return ids;
	}
	
	/**
	 * Map-based sketch search algorithm using IntHashMap for hit counting.
	 * More memory efficient than list-based approach for large result sets.
	 * Tracks hit counts per sketch and supports contamination analysis.
	 *
	 * @param a Query sketch to search for
	 * @param params Display parameters including minimum hits and contamination level
	 * @return SketchResults with matching sketches above hit threshold
	 */
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
	
	/**
	 * Worker thread for parallel indexing of reference sketches.
	 * Processes assigned sketches by adding their k-mers to hash tables with sketch IDs.
	 * Uses HashBuffer for efficient batch insertions and tracks processing statistics.
	 */
	public class IndexThread extends Thread {
		
		/**
		 * Constructs IndexThread with shared counters for work distribution.
		 * Initializes HashBuffer for efficient k-mer insertion into hash tables.
		 *
		 * @param nextIndex_ Atomic counter for next sketch to process
		 * @param keyCount_ Atomic counter for total k-mers processed
		 * @param uniqueKeyCount_ Atomic counter for unique k-mers added
		 */
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
		
		/** Shared counter for thread work distribution across reference sketches */
		AtomicInteger nextIndex;
		/** Total number of k-mer keys processed by all threads */
		AtomicLong keyCount;
		/** Number of unique k-mer keys added to the index */
		AtomicLong uniqueKeyCount;
		/** Thread-local counter of k-mer codes processed by this thread */
		long codesProcessedT=0;
		/** Buffer for efficient batch insertion of k-mers into hash tables */
		HashBuffer buffer;
		/** Indicates successful completion of this thread's indexing work */
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** Set of hash tables for distributed k-mer storage */
	public final KmerTableSet tables;
	/** Array of individual hash tables from the table set */
	public final AbstractKmerTable[] tableArray;
	/** Reference sketches being indexed for similarity searches */
	public final ArrayList<Sketch> refSketches;
	
	/** Indicates if an error occurred during indexing operations */
	public boolean errorState=false;

	/** Debug flag to enable timing output for performance analysis */
	private static final boolean printTime=false;
	/** Controls whether to use IntHashMap or IntList for hit accumulation */
	public static boolean useIntMap=true;
//	public static boolean useIntMapBinary=false;
	/** Initial size for IntHashMap when using map-based search algorithm */
	public static int intMapSize=1000;
	/** Maximum number of k-mers to index per sketch */
	public static int indexLimit=Integer.MAX_VALUE;
	/** Number of hash table ways for distributing k-mers across tables */
	public static final int WAYS=31;
	
}
