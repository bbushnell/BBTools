package var2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import ml.CellNet;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
* Thread-safe container for managing large collections of genomic variants.
* Uses sharded ConcurrentHashMap architecture for optimal concurrent performance
* with minimal lock contention. Provides efficient storage, retrieval, and
* processing of variants across multiple threads.
* 
* Key features:
* - Sharded storage using position-based hashing for load distribution
* - Multithreaded variant processing with statistical accumulation
* - Custom iterator for seamless traversal across all shards
* - Nearby variant analysis for artifact detection
* - Comprehensive filtering and quality assessment pipeline
* 
* @author Brian Bushnell
* @author Isla
* @date December 2024
*/
public class VarMap implements Iterable<Var> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Construction          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates VarMap with scaffold mapping but default processing parameters.
	 * @param scafMap_ Scaffold mapping for variant coordinate resolution
	 */
	VarMap(ScafMap scafMap_){
		this(scafMap_, -1, -1, -1, -1, -1);
	}

	/**
	 * Creates VarMap with full initialization of processing parameters.
	 * Initializes sharded storage with WAYS concurrent hash maps for optimal performance.
	 * 
	 * @param scafMap_ Scaffold mapping for coordinate resolution
	 * @param ploidy_ Expected organism ploidy level
	 * @param pairingRate_ Proper pair rate for dataset normalization
	 * @param totalQualityAvg_ Average base quality across dataset
	 * @param mapqAvg_ Average mapping quality across dataset
	 * @param readLengthAvg_ Average read length for bias corrections
	 */
	@SuppressWarnings("unchecked")
	/**
	* Creates VarMap with full initialization of processing parameters.
	* Initializes sharded storage with WAYS concurrent hash maps for optimal performance.
	* 
	* @param scafMap_ Scaffold mapping for coordinate resolution
	* @param ploidy_ Expected organism ploidy level
	* @param pairingRate_ Proper pair rate for dataset normalization
	* @param totalQualityAvg_ Average base quality across dataset
	* @param mapqAvg_ Average mapping quality across dataset
	* @param readLengthAvg_ Average read length for bias corrections
	*/
	VarMap(ScafMap scafMap_, int ploidy_, double pairingRate_, double totalQualityAvg_,
			double mapqAvg_, double readLengthAvg_){
		scafMap=scafMap_;
		ploidy=ploidy_;
		properPairRate=pairingRate_;
		totalQualityAvg=totalQualityAvg_;
		totalMapqAvg=mapqAvg_;
		readLengthAvg=readLengthAvg_;
		
		//Initialize sharded storage for concurrent access
		maps=new ConcurrentHashMap[WAYS];
		for(int i=0; i<WAYS; i++){
			maps[i]=new ConcurrentHashMap<Var, Var>();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Counts nearby variants using default filter parameters.
	 * Convenience wrapper for nearby variant analysis.
	 * @param varFilter Filter containing nearby variant parameters
	 * @return Number of variants exceeding nearby variant thresholds
	 */
	public int countNearbyVars(VarFilter varFilter){
		return countNearbyVars(varFilter, varFilter.maxNearbyCount, varFilter.nearbyDist, 
				varFilter.nearbyGap, varFilter.flagNearby);
	}
	
	/**
	* Analyzes clustering of variants to detect potential artifacts.
	* Scans for variants within specified distance and gap thresholds
	* to identify regions with suspiciously high variant density.
	* 
	* @param varFilter Filter for quality assessment during scanning
	* @param maxCount0 Maximum allowed nearby variants
	* @param maxDist Maximum distance to scan for nearby variants
	* @param maxGap Maximum gap between consecutive variants in cluster
	* @param flag Whether to mark variants exceeding thresholds
	* @return Number of variants exceeding nearby count threshold
	*/
	public int countNearbyVars(VarFilter varFilter, final int maxCount0, final int maxDist, final int maxGap, final boolean flag){
		final int maxCount=maxCount0<0 ? 19 : Tools.mid(maxCount0, 8, 19); // Clamp to reasonable range
		final Var[] array=toArray(true); // Get sorted array for positional scanning
		int failed=0;
		
		for(int vloc=0; vloc<array.length; vloc++){
			int x=countNearbyVars(varFilter, array, vloc, maxCount, maxDist, maxGap, flag);
			if(x>maxCount){failed++;}
		}
		return failed;
	}
	
	/**
	 * Tests if a variant passes quality filters independently.
	 * Used during nearby variant analysis to determine which variants count toward clustering.
	 * Only determines whether the nearby variants pass.
	 * 
	 * @param v Variant to test
	 * @param varFilter Quality filter to apply
	 * @return True if variant passes quality thresholds
	 */
	private boolean passesSolo(Var v, VarFilter varFilter){
		assert(varFilter!=null);
		if(varFilter==null){return true;}
		
		boolean pass=varFilter.passesFast(v); // Quick filter checks first
		if(pass){
			v.calcCoverage(scafMap); // Calculate coverage if needed
			pass=v.forced() || varFilter.passesFilter(v, properPairRate, totalQualityAvg,
					totalMapqAvg, readLengthAvg, ploidy, scafMap, null, false);
		}
		return pass;
	}
	
	/**
	 * Counts nearby variants for a specific variant position.
	 * Scans left and right from target position to identify clustering.
	 * Only counts variants that pass quality filters.
	 * 
	 * @param varFilter Quality filter for determining which variants to count
	 * @param array Sorted array of all variants
	 * @param vloc0 Index of target variant in array
	 * @param maxCount Maximum nearby variants before flagging
	 * @param maxDist Maximum distance to scan from target
	 * @param maxGap Maximum gap between consecutive variants
	 * @param flag Whether to flag variants exceeding threshold
	 * @return Number of nearby variants found
	 */
	public int countNearbyVars(VarFilter varFilter, final Var[] array, final int vloc0, final int maxCount, 
			final int maxDist, final int maxGap, final boolean flag){
		final Var v0=array[vloc0];
		assert(v0.nearbyVarCount==-1) : "Nearby vars were already counted?";
		int nearby=0;
		
		//Scan leftward from target position
		{
			Var prev=v0;
			for(int i=vloc0-1; i>=0 && nearby<=maxCount; i--){
				final Var v=array[i];
				//Stop if gap too large or distance too far
				if(prev.start-v.stop>maxGap || v0.start-v.stop>maxDist){break;}
				
				if(!v.forced() || passesSolo(v, varFilter)){
					nearby++;
					prev=v; // Update for gap calculation
				}
			}
		}
		
		//Scan rightward from target position
		{
			Var prev=v0;
			for(int i=vloc0+1; i<array.length && nearby<=maxCount; i++){
				final Var v=array[i];
				//Stop if gap too large or distance too far
				if(v.start-prev.stop>maxGap || v.start-v0.stop>maxDist){break;}
				
				if(!v.forced() || passesSolo(v, varFilter)){
					nearby++;
					prev=v; // Update for gap calculation
				}
			}
		}
		
		v0.nearbyVarCount=nearby;
		if(flag && nearby>varFilter.maxNearbyCount){
			v0.setFlagged(true); // Mark for special handling
		}
		return nearby;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Tests if VarMap contains a specific variant.
	 * @param v Variant to search for
	 * @return True if variant is present
	 */
	public boolean containsKey(Var v){
		return get(v)!=null;
	}
	
	/**
	 * Retrieves variant from appropriate shard.
	 * Uses position-based hashing to determine correct map.
	 * @param v Variant key to search for
	 * @return Stored variant or null if not found
	 */
	Var get(final Var v){
		final int way=v.start&MASK; // Hash position to determine shard
		return maps[way].get(v);
	}
	
	/**
	 * Returns total number of variants across all shards.
	 * @return Total variant count
	 */
	public long size(){
		long size=0;
		for(int i=0; i<maps.length; i++){size+=maps[i].size();}
		return size;
	}
	
	/**
	 * Alternative size calculation using iterator (slow).
	 * Used for debugging and validation only.
	 * @return Total variant count via iteration
	 */
	public long size2(){
		assert(false) : "Slow";
		int i=0;
		for(Var v : this){i++;}
		return i;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Adders            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds variant to appropriate shard with thread safety.
	 * Merges with existing variant if already present.
	 * @param v Variant to add
	 * @return 1 if new variant added, 0 if merged with existing
	 */
	private int add(final Var v){
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		synchronized(map){
			Var old=map.get(v);
			if(old==null){
				map.put(v, v); // Add new variant
				return 1;
			}
			else{
				synchronized(old){
					old.add(v); // Merge statistics
				}
			}
		}
		return 0;
	}
	
	/**
	 * Adds variant without synchronization (for single-threaded use).
	 * @param v Variant to add
	 * @return 1 if new variant added, 0 if merged
	 */
	int addUnsynchronized(final Var v){
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		Var old=map.get(v);
		if(old==null){
			map.put(v, v);
			return 1;
		}
		old.add(v);
		return 0;
	}
	
	/**
	 * Removes variant without synchronization.
	 * @param v Variant to remove
	 * @return 1 if variant was present and removed, 0 otherwise
	 */
	int removeUnsynchronized(Var v){
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		return map.remove(v)==null ? 0 : 1;
	}
	
	/**
	* Efficiently merges thread-local variant collections into main storage.
	* Groups variants by shard for batch processing with minimal locking.
	* 
	* @param mapT Thread-local variant map to merge
	* @return Number of new variants added
	*/
	int dumpVars(HashMap<Var, Var> mapT){
		int added=0;
		@SuppressWarnings("unchecked")
		ArrayList<Var>[] absent=new ArrayList[WAYS];
		
		//Initialize per-shard collections
		for(int i=0; i<WAYS; i++){
			absent[i]=new ArrayList<Var>();
		}
		
		//Sort variants by target shard and attempt unlocked merges
		for(Entry<Var, Var> e : mapT.entrySet()){
			Var v=e.getValue();
			final int way=v.start&MASK;
			ConcurrentHashMap<Var, Var> map=maps[way];
			Var old=map.get(v);
			if(old==null){
				absent[way].add(v); // Defer for locked insertion
			}
			else{
				synchronized(old){
					old.add(v); // Safe to merge immediately
				}
			}
		}
		
		mapT.clear(); // Free thread-local memory
		
		//Process deferred insertions with proper locking
		for(int way=0; way<WAYS; way++){
			ConcurrentHashMap<Var, Var> map=maps[way];
			ArrayList<Var> list=absent[way];
			synchronized(map){
				for(Var v : list){
					Var old=get(v);
					if(old==null){
						map.put(v, v);
						added++;
					}
					else{
						synchronized(old){
							old.add(v); // Race condition resolved
						}
					}
				}
			}
		}
		return added;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Other             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Single-threaded variant processing for debugging and validation.
	 * @param filter Quality filter to apply
	 * @param net CellNet for prediction
	 * @param scoreArray Score histogram arrays to populate
	 * @param ploidyArray Ploidy distribution array
	 * @param avgQualityArray Average quality histograms
	 * @param maxQualityArray Maximum quality histogram
	 * @param ADArray Allele depth arrays
	 * @param AFArray Allele frequency arrays
	 * @return Array of variant counts by type
	 */
	public long[] processVariantsST(VarFilter filter, CellNet net, long[][] scoreArray, long[] ploidyArray, long[][] avgQualityArray,
			long[] maxQualityArray, long[][] ADArray, double[] AFArray){
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		long[] types=new long[Var.VAR_TYPES];
		for(ConcurrentHashMap<Var, Var> map : maps){
			//Three-pass processing: fast filter, insertion bias correction, full statistics
			long[] types2=processVariants(map, filter, null, null, null, null, null, null, null, false, false);
			types2=processVariants(map, filter, null, null, null, null, null, null, null, true, false);
			types2=processVariants(map, filter, net, scoreArray, ploidyArray, avgQualityArray, maxQualityArray, ADArray, AFArray, false, false);
			Tools.add(types, types2);
		}
		return types;
	}
	
	/**
	 * Multithreaded variant processing with statistical accumulation.
	 * Processes all shards in parallel for optimal performance.
	 * 
	 * @param filter Quality filter for variant assessment
	 * @param net CellNet for prediction
	 * @param scoreArray Score histogram collection
	 * @param ploidyArray Ploidy distribution tracking
	 * @param avgQualityArray Quality histogram collection  
	 * @param maxQualityArray Maximum quality distribution
	 * @param ADArray Allele depth statistics
	 * @param AFArray Allele frequency statistics
	 * @return Variant type counts across all processed variants
	 */
	public long[] processVariantsMT(VarFilter filter, CellNet net, long[][] scoreArray, long[] ploidyArray, 
			long[][] avgQualityArray, long[] maxQualityArray, long[][] ADArray, double[] AFArray){
		//Three-pass processing for statistical analysis
		processVariantsMT_inner(filter, null, null, null, null, null, null, null, false); //Initial filtering
		processVariantsMT_inner(filter, null, null, null, null, null, null, null, true); //Insertion bias correction
		return processVariantsMT_inner(filter, net, scoreArray, ploidyArray, avgQualityArray, maxQualityArray, 
				ADArray, AFArray, false);
	}
	
	/**
	 * Internal multithreaded processing implementation.
	 * Creates one processing thread per shard for parallel execution.
	 * 
	 * @param filter Quality filter to apply
	 * @param net CellNet for prediction
	 * @param scoreArray Optional score histogram arrays
	 * @param ploidyArray Optional ploidy distribution array
	 * @param avgQualityArray Optional quality histogram arrays
	 * @param maxQualityArray Optional maximum quality histogram
	 * @param ADArray Optional allele depth arrays
	 * @param AFArray Optional allele frequency arrays
	 * @param processInsertions Whether to process insertion bias corrections
	 * @return Accumulated variant type counts
	 */
	private long[] processVariantsMT_inner(VarFilter filter, CellNet net, long[][] scoreArray, long[] ploidyArray, 
			long[][] avgQualityArray, long[] maxQualityArray, long[][] ADArray, double[] AFArray, boolean processInsertions){
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		//Create processing threads for each shard
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(WAYS);
		for(int i=0; i<WAYS; i++){
			ProcessThread pt=new ProcessThread(maps[i], filter, net, scoreArray!=null, ploidyArray!=null, 
					avgQualityArray!=null, ADArray!=null, processInsertions);
			alpt.add(pt);
			pt.start();
		}
		
		//Collect results from all threads
		long[] types=new long[Var.VAR_TYPES];
		boolean success=true;
		for(ProcessThread pt : alpt){
			//Wait for thread completion
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e){
					e.printStackTrace();
				}
			}
			
			//Accumulate statistics from each thread
			if(pt.types!=null){
				Tools.add(types, pt.types);
			}
			if(scoreArray!=null){Tools.add(scoreArray, pt.scoreArray);}
			if(ploidyArray!=null){Tools.add(ploidyArray, pt.ploidyArray);}
			if(avgQualityArray!=null){Tools.add(avgQualityArray, pt.avgQualityArray);}
			if(maxQualityArray!=null){Tools.add(maxQualityArray, pt.maxQualityArray);}
			if(ADArray!=null){Tools.add(ADArray, pt.ADArray);}
			if(ADArray!=null){Tools.add(AFArray, pt.AFArray);} // Note: triggered on ADArray check
			success&=pt.success;
		}
		
		return types;
	}
	
	/**
	 * Worker thread for parallel variant processing.
	 * Handles one shard of the variant collection with independent statistics.
	 */
	private class ProcessThread extends Thread {
		
		/**
		 * Creates processing thread for single shard.
		 * @param map_ Shard to process
		 * @param filter_ Quality filter to apply
		 * @param trackScores Whether to collect score statistics
		 * @param trackPloidy Whether to collect ploidy statistics
		 * @param trackQuality Whether to collect quality statistics
		 * @param trackAD Whether to collect allele depth statistics
		 * @param processInsertions_ Whether to process insertion bias corrections
		 */
		ProcessThread(Map<Var, Var> map_, VarFilter filter_, CellNet net0_, boolean trackScores, boolean trackPloidy, 
				boolean trackQuality, boolean trackAD, boolean processInsertions_){
			map=map_;
			filter=filter_;
			net0=net0_;
			scoreArray=(trackScores ? new long[8][200] : null);
			ploidyArray=(trackPloidy ? new long[ploidy+1] : null);
			avgQualityArray=(trackQuality ? new long[8][100] : null);
			maxQualityArray=(trackQuality ? new long[100] : null);
			ADArray=(trackAD ? new long[2][7] : null);
			AFArray=(trackAD ? new double[7] : null);
			processInsertions=processInsertions_;
		}
		
		/**
		 * Executes variant processing for assigned shard.
		 */
		@Override
		public void run(){
			net=(net0==null ? null : net0.copy(false));
			types=processVariants(map, filter, net, scoreArray, ploidyArray, avgQualityArray, maxQualityArray, 
					ADArray, AFArray, processInsertions, false);
			success=true;
		}
		
		/**Quality filter for variant assessment*/
		final VarFilter filter;
		/**Shard to process*/
		final Map<Var, Var> map;
		/**CellNet instance for machine learning predictions*/
		final CellNet net0;
		/**Thread-local copy of CellNet*/
		CellNet net;
		/**Variant type counts*/
		long[] types;
		/**Score histogram arrays*/
		final long[][] scoreArray;
		/**Ploidy distribution array*/
		final long[] ploidyArray;
		/**Average quality histogram arrays*/
		final long[][] avgQualityArray;
		/**Maximum quality histogram*/
		final long[] maxQualityArray;
		/**Allele depth arrays*/
		final long[][] ADArray;
		/**Allele frequency arrays*/
		final double[] AFArray;
		/**Whether to process insertion bias corrections*/
		boolean processInsertions;
		/**Processing success flag*/
		boolean success=false;
	}

	/**
	* Core variant processing logic for filtering and statistics collection.
	* Handles both insertion bias correction and comprehensive quality assessment.
	* 
	* @param map Single shard to process
	* @param filter Quality filter for variant assessment
	* @param scoreArray Optional score histogram collection
	* @param ploidyArray Optional ploidy distribution tracking
	* @param avgQualityArray Optional quality histogram collection
	* @param maxQualityArray Optional maximum quality tracking
	* @param ADArray Optional allele depth collection
	* @param AFArray Optional allele frequency collection
	* @param processInsertions Whether to handle insertion bias corrections
	* @param considerNearby Whether to consider nearby variant counts in filtering
	* @return Variant type counts for processed variants
	*/
	private long[] processVariants(Map<Var, Var> map, VarFilter filter, CellNet net, long[][] scoreArray, long[] ploidyArray, 
			long[][] avgQualityArray, long[] maxQualityArray, long[][] ADArray, double[] AFArray, boolean processInsertions, boolean considerNearby){
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		Iterator<Entry<Var, Var>> iterator=map.entrySet().iterator();
		long[] types=new long[Var.VAR_TYPES];
		
		while(iterator.hasNext()){
			Entry<Var, Var> entry=iterator.next();
			final Var v=entry.getValue();
			
			if(processInsertions){
				//Handle insertion bias correction pass
				assert(readLengthAvg>0);
				if(v.type()==Var.INS){
					synchronized(v){
						v.reviseAlleleFraction(readLengthAvg, scafMap.getScaffold(v.scafnum), this);
					}
				}
			}else{
				//Handle filtering and statistics collection pass
				boolean pass=filter.passesFast(v);
				if(pass){
					v.calcCoverage(scafMap);
					pass=v.forced() || filter.passesFilter(v, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap, net, considerNearby);
				}
				
				if(pass){
					types[v.type()]++;
					
					//Collect score statistics if requested
					//TODO: phredScore is calculated twice
					if(scoreArray!=null){
						int score=(int)v.phredScore(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, 
								filter.rarity, ploidy, scafMap, net);
						scoreArray[0][score]++; //Overall scores
						scoreArray[v.type()+1][score]++; //Type-specific scores
					}
					
					//Collect ploidy statistics
					if(ploidyArray!=null){ploidyArray[v.calcCopies(ploidy)]++;}
					
					//Collect quality statistics
					if(avgQualityArray!=null){
						int q=(int)v.baseQAvg();
						avgQualityArray[0][q]++; // Overall quality
						avgQualityArray[v.type()+1][q]++; // Type-specific quality
					}
					if(maxQualityArray!=null){maxQualityArray[(int)v.baseQMax]++;}
					
					//Collect depth and frequency statistics
					if(ADArray!=null){
						ADArray[0][v.type()]+=v.alleleCount(); // Allele depth by type
						ADArray[1][v.type()]+=v.coverage(); // Reference depth by type
					}
					if(AFArray!=null){AFArray[v.type()]+=v.alleleFraction();}
				}else{
					iterator.remove(); // Remove variants that fail filtering
				}
			}
		}
		return types;
	}

	/**
	 * Adds shared variants from another map (used in multi-sample processing).
	 * @param map Target shard to modify
	 * @param sharedMap Source of shared variants
	 * @return Variant type counts for added variants
	 */
	private long[] addSharedVariants(Map<Var, Var> map, Map<Var, Var> sharedMap){
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		//Add missing shared variants
		for(Var v : sharedMap.keySet()){
			if(!map.containsKey(v)){
				Var v2=new Var(v); // Create copy for this sample
				map.put(v2, v2);
			}
		}
		
		//Count added variants by type
		long[] types=new long[Var.VAR_TYPES];
		for(Var v : sharedMap.keySet()){
			v.calcCoverage(scafMap);
			types[v.type()]++;
		}
		return types;
	}
	
	/**
	 * Creates sorted array of all variants for positional analysis.
	 * @param sort Whether to sort array by genomic position
	 * @return Array containing all variants
	 */
	public Var[] toArray(boolean sort){
		Var[] array=new Var[(int)size()];
		int i=0;
		
		for(Var v : this){
			assert(i<array.length);
			array[i]=v;
			i++;
		}
		if(sort){Shared.sort(array);} // Sort by position for analysis
		return array;
	}
	
	/**
	 * Validates internal data structure consistency (slow debugging method).
	 * Checks that all key-value pairs are properly mapped and no variants
	 * appear in multiple shards.
	 * @param quiet Whether to suppress detailed output
	 * @return True if all validation checks pass
	 */
	private boolean mappedToSelf(boolean quiet){
		assert(false) : "Slow";
		
		//Validate each shard independently
		for(ConcurrentHashMap<Var, Var> map : maps){
			for(Var key : map.keySet()){
				Var value=map.get(key);
				assert(value!=null);
				assert(value.equals(key));
				assert(value==key); // Should be same reference
				assert(map.get(value).equals(key));
			}
			
			//Validate entry set consistency
			for(Entry<Var, Var> e : map.entrySet()){
				Var key=e.getKey();
				Var value=e.getValue();
				assert(value!=null);
				assert(value.equals(key));
				assert(value==key);
			}
			
			//Ensure no cross-shard contamination
			for(ConcurrentHashMap<Var, Var> map2 : maps){
				if(map2!=map){
					for(Var key : map.keySet()){
						assert(!map2.containsKey(key));
					}
				}
			}
		}
		
		//Validate iterator consistency
		int i=0;
		for(Var v : this){
			if(!quiet){
				System.err.println(i+"\t"+v.start+"\t"+v.stop+"\t"+v.toKey()+"\t"+v.hashcode+"\t"+v.hashCode()+"\t"+new String(v.allele)+"\t"+((Object)v).hashCode());
			}
			Var v2=get(v);
			assert(v==v2);
			assert(get(v2)==v);
			assert(get(v)==v) : "\n"+i+"\t"+v2.start+"\t"+v2.stop+"\t"+v2.toKey()+"\t"+v2.hashcode+"\t"+v2.hashCode()+"\t"+new String(v2.allele)+"\t"+((Object)v2).hashCode();
			i++;
		}
		assert(i==size()) : i+", "+size()+", "+size2();
		return true;
	}
	
	/**
	 * Calculates coverage for all variants and returns type distribution.
	 * @param scafMap Scaffold mapping for coverage calculation
	 * @return Array of variant counts by type
	 */
	public long[] calcCoverage(ScafMap scafMap){
		long[] types=new long[Var.VAR_TYPES];
		for(Var v : this){
			v.calcCoverage(scafMap);
			types[v.type()]++;
		}
		return types;
	}
	
	/**
	 * Counts variants by type without additional processing.
	 * @return Array of variant counts by type
	 */
	public long[] countTypes(){
		long[] types=new long[Var.VAR_TYPES];
		for(Var v : this){
			types[v.type()]++;
		}
		return types;
	}
	
	/**
	 * Clears all variants and resets processing parameters.
	 * Reinitializes all shards for fresh usage.
	 */
	public void clear(){
		properPairRate=-1;
		pairedInSequencingRate=-1;
		totalQualityAvg=-1;
		totalMapqAvg=-1;
		readLengthAvg=-1;
		
		//Reinitialize all shards
		for(int i=0; i<maps.length; i++){
			maps[i]=new ConcurrentHashMap<Var, Var>();
		}
	}
	
	/**
	 * Creates string representation of all variants.
	 * @return Formatted string containing all variant data
	 */
	@Override
	/**
	* Creates string representation of all variants.
	* @return Formatted string containing all variant data
	*/
	public String toString(){
		ByteBuilder sb=new ByteBuilder();
		for(ConcurrentHashMap<Var, Var> map : maps){
			for(Var v : map.keySet()){
				v.toTextQuick(sb);
				sb.nl();
			}
		}
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Iteration           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates iterator for traversing all variants across shards.
	 * @return Custom iterator that handles sharded storage
	 */
	@Override
	/**
	* Creates iterator for traversing all variants across shards.
	* @return Custom iterator that seamlessly handles sharded storage
	*/
	public VarMapIterator iterator(){
		return new VarMapIterator();
	}

	/**
	* Custom iterator for traversing all variants across multiple hash map shards.
	* Seamlessly moves through each shard's entries without exposing the internal
	* partitioning structure to the caller.
	*/
	private class VarMapIterator implements Iterator<Var> {

		/**Creates iterator and prepares first available shard for iteration.*/
		VarMapIterator(){
			makeReady(); // Initialize to first non-empty shard
		}

		/**
		 * Checks if more variants are available across any remaining shards.
		 * @return True if more variants exist
		 */
		@Override
		/**
		* Checks if more variants are available across any remaining shards.
		* @return True if more variants exist
		*/
		public boolean hasNext(){
			return iter.hasNext();
		}

		/**
		 * Returns next variant and advances iterator position.
		 * Automatically transitions to next shard when current shard is exhausted.
		 * @return Next Var object in iteration sequence
		 */
		@Override
		/**
		* Returns next variant and advances iterator position.
		* Automatically transitions to next shard when current shard is exhausted.
		* @return Next Var object in iteration sequence
		*/
		public Var next(){
			Entry<Var, Var> e=iter.next();
			if(!iter.hasNext()){makeReady();} // Prepare next shard if current is exhausted
			Var v=e==null ? null : e.getValue(); // Extract variant from map entry
			return v;
		}

		/**
		 * Advances to next available shard with variants.
		 * Skips empty shards and prepares iterator for next non-empty collection.
		 */
		private void makeReady(){
			while((iter==null || !iter.hasNext()) && nextMap<maps.length){
				iter=maps[nextMap].entrySet().iterator(); // Get iterator for current shard
				nextMap++; // Advance to next shard for subsequent calls
			}
		}

		/**Index of next shard to examine*/
		private int nextMap=0;
		/**Current shard's entry iterator*/
		private Iterator<Entry<Var, Var>> iter=null;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/**Expected organism ploidy level for variant calling*/
	public int ploidy=-1;
	/**Fraction of reads that mapped as proper pairs*/
	public double properPairRate=-1;
	/**Fraction of reads that were paired in sequencing*/
	public double pairedInSequencingRate=-1;
	/**Average base quality across all processed reads*/
	public double totalQualityAvg=-1;
	/**Average mapping quality across all processed reads*/
	public double totalMapqAvg=-1;
	/**Average read length across all processed reads*/
	public double readLengthAvg=-1;
	/**Scaffold mapping for coordinate resolution and reference access*/
	public final ScafMap scafMap;
	/**Array of concurrent hash maps for sharded variant storage*/
	final ConcurrentHashMap<Var, Var>[] maps;

	/*--------------------------------------------------------------*/
	/*----------------        Static fields         ----------------*/
	/*--------------------------------------------------------------*/

	/**Serialization version identifier*/
	private static final long serialVersionUID=1L;
	/** Number of hash map shards (must be power of 2 for efficient masking) */
	private static final int WAYS=8;
	/**Bit mask for shard selection (WAYS-1)*/
	public static final int MASK=WAYS-1;

}
