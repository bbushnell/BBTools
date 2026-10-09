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
 * Eight position-sharded maps holding mutable variants as both keys and values.
 * New insertions retain the supplied Var; equivalent variants merge read evidence.
 * ConcurrentHashMap protects individual map operations, not all compound operations
 * or mutable Var fields. dumpVars coordinates concurrent merges; unsynchronized
 * methods require exclusive access or an external lock shared by their callers.
 *
 * Finish accumulation before counting, processing, snapshotting or clearing. Callers
 * coordinate mutable configuration, scaffold coverage and the lifetimes of borrowed
 * variants. clear replaces shards; it is not a concurrent reset. Processing removes
 * rejected variants and mutates coverage/revised-AF state. Iteration is not a snapshot.
 * @author Brian Bushnell
 * @author Isla
 * @date December 2024
 */
public class VarMap implements Iterable<Var>{

	/*--------------------------------------------------------------*/
	/*----------------        Construction          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Retains the scaffold map and initializes processing parameters to -1 sentinels.
	 * Callers populate the required parameters before scoring or processing.
	 * @param scafMap_ Borrowed scaffold map for coverage/reference access
	 */
	VarMap(ScafMap scafMap_){
		this(scafMap_, -1, -1, -1, -1, -1);
	}

	/**
	 * Retains configuration and allocates eight initially empty shards.
	 * The scaffold map is borrowed, and processing parameters are stored without validation.
	 * @param scafMap_ Scaffold map for coverage/reference access
	 * @param ploidy_ Sample ploidy
	 * @param pairingRate_ Dataset proper-pair fraction
	 * @param totalQualityAvg_ Dataset mean base quality
	 * @param mapqAvg_ Dataset mean mapping quality
	 * @param readLengthAvg_ Dataset mean read length in bases
	 */
	@SuppressWarnings("unchecked")
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
	 * Counts and optionally flags neighbors using this filter's limits.
	 * Requires stable contents and unset nearby counts; mutates the stored variants.
	 * @param varFilter Non-null filter supplying count, distance, gap and flag settings
	 * @return Count above maxNearbyCount, or above19 for a negative count limit
	 */
	public int countNearbyVars(VarFilter varFilter){
		return countNearbyVars(varFilter, varFilter.maxNearbyCount, varFilter.nearbyDist,
				varFilter.nearbyGap, varFilter.flagNearby);
	}

	/**
	 * Counts neighbors in a sorted reference snapshot and reports counts above maxCount0.
	 * A negative maxCount0 retains threshold19. The scan cap is at least8 and at least
	 * the reporting threshold; it may store cap+1 as the first excess. Counts therefore
	 * need not be complete cluster sizes. Requires unset nearby counts and stable contents.
	 * Flagging uses varFilter.maxNearbyCount separately; existing flags are not cleared.
	 * @param varFilter Non-null contributor/flag configuration
	 * @param maxCount0 Reporting threshold; negative selects the historical19 threshold
	 * @param maxDist Maximum signed coordinate gap from target, in bases
	 * @param maxGap Maximum signed gap from last accepted neighbor, in bases
	 * @param flag Whether to mark variants above varFilter.maxNearbyCount
	 * @return Number above the reporting threshold, including forced variants
	 */
	public int countNearbyVars(VarFilter varFilter, final int maxCount0, final int maxDist, final int maxGap, final boolean flag){
		//VAR2-025: callers gate rejection on failed, so report against their actual limit.
		//Keep the historical minimum count cap8, but never truncate below a larger limit.
		final int threshold=(maxCount0<0 ? 19 : maxCount0);
		final int maxCount=Tools.max(threshold, 8);
		final Var[] array=toArray(true); // Get sorted array for positional scanning
		int failed=0;

		for(int vloc=0; vloc<array.length; vloc++){
			int x=countNearbyVars(varFilter, array, vloc, maxCount, maxDist, maxGap, flag);
			if(x>threshold){failed++;}
		}
		return failed;
	}

	/**
	 * Runs fast checks, calculates missing coverage, then applies emission filters
	 * without a network or nearby rejection. Forced variants bypass both filter stages,
	 * but still reach coverage calculation. Null filter asserts, or passes with assertions off.
	 * @param v Non-null variant whose coverage cache may change
	 * @param varFilter Required emission filter
	 * @return Whether the implemented checks accept the variant
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
	 * Tests the independent contributor depth/AF thresholds used for nearby counts.
	 * Uses cached revised AF when available, otherwise raw AF; does not compute a revision.
	 * Does not exempt forced variants or apply the other emission thresholds.
	 * @param v Non-null variant with usable counts/coverage when AF is requested
	 * @param varFilter Non-null contributor configuration
	 * @return Whether this variant may contribute to another variant's nearby count
	 */
	private static boolean countsTowardNVC(Var v, VarFilter varFilter){
		if(v.alleleCount()<varFilter.nvcMinCount){return false;}
		if(varFilter.nvcMaf>0){
			final double af=v.revisedAlleleFraction==-1 ? v.alleleFraction() : v.revisedAlleleFraction;
			if(af<varFilter.nvcMaf){return false;}
		}
		return true;
	}

	/**
	 * Counts left then right on the target scaffold, stopping at the distance/gap
	 * limits or after the first count above maxCount. Gap is measured from the last
	 * accepted neighbor, while distance is measured from the target. Ordinary variants
	 * skip passesSolo; callers normally prefilter them. Forced variants call passesSolo,
	 * which itself bypasses emission thresholds. All contributors receive the NVC gate.
	 * Writes the target's nearbyVarCount and may set its flag; an existing flag is not cleared.
	 * @param varFilter Non-null filter supplying contributor and flagging thresholds
	 * @param array Stable array sorted by Var.compareTo, containing borrowed variants
	 * @param vloc0 Target index; its nearby count must be -1 under assertions
	 * @param maxCount Scan cap; a result may be maxCount+1
	 * @param maxDist Maximum signed coordinate gap from target to neighbor, in bases
	 * @param maxGap Maximum signed gap from last accepted neighbor, in bases
	 * @param flag Whether to set the target flag above varFilter.maxNearbyCount
	 * @return Stored nearby count, potentially truncated at the first excess
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
				//VAR2-024: coordinates are local to each scaffold; sorted neighbors stop here.
				if(v.scafnum!=v0.scafnum){break;}
				//Stop if gap too large or distance too far
				if(prev.start-v.stop>maxGap || v0.start-v.stop>maxDist){break;}

				if((!v.forced() || passesSolo(v, varFilter)) && countsTowardNVC(v, varFilter)){
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
				//VAR2-024: coordinates are local to each scaffold; sorted neighbors stop here.
				if(v.scafnum!=v0.scafnum){break;}
				//Stop if gap too large or distance too far
				if(v.start-prev.stop>maxGap || v.start-v0.stop>maxDist){break;}

				if((!v.forced() || passesSolo(v, varFilter)) && countsTowardNVC(v, varFilter)){
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
	 * Looks up an equivalent variant by identity and start-derived shard.
	 * @param v Non-null lookup key
	 * @return Whether a stored value is found
	 */
	public boolean containsKey(Var v){
		return get(v)!=null;
	}

	/**
	 * Returns the borrowed mutable value from the shard selected by start's low bits.
	 * Caller must preserve key identity/allele contents and coordinate value mutation.
	 * @param v Non-null lookup key
	 * @return Stored variant, or null
	 */
	Var get(final Var v){
		final int way=v.start&MASK; // Hash position to determine shard
		return maps[way].get(v);
	}

	/**
	 * Sums shard sizes. Concurrent modification can make this an inexact observation;
	 * no whole-map snapshot or cross-shard lock is taken.
	 * @return Observed total size
	 */
	public long size(){
		long size=0;
		for(int i=0; i<maps.length; i++){size+=maps[i].size();}
		return size;
	}

	/**
	 * Legacy iterator count with an unconditional failing assertion on entry.
	 * With assertions disabled it uses an int accumulator; this is not a production validator.
	 * @return Iterator count when the entry guard is disabled
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
	 * Inserts the borrowed object or merges Var.add read evidence into an equal value.
	 * Locks the shard for the compound lookup/insertion and the existing value for merging.
	 * These locks do not protect concurrent resets or callers using unsynchronized methods.
	 * @param v Non-null variant; a new insertion transfers its mutation ownership to the map
	 * @return 1 for an insertion, 0 for a merge
	 */
	private int add(final Var v){
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		synchronized(map){
			Var old=map.get(v);
			if(old==null){
				map.put(v, v); // Add new variant
				return 1;
			}else{
				synchronized(old){
					old.add(v); // Merge statistics
				}
			}
		}
		return 0;
	}

	/**
	 * Inserts a borrowed variant or merges read evidence without compound-operation locks.
	 * Requires exclusive access or an external lock shared by all relevant writers.
	 * @param v Non-null variant; newly stored objects must not be independently mutated
	 * @return 1 for an insertion, 0 for a merge
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
	 * Removes the matching entry from its shard. Coordinate with accumulation/processing
	 * so an ongoing merge cannot update a value that has just been detached.
	 * @param v Non-null identity key
	 * @return 1 if an entry was removed, otherwise 0
	 */
	int removeUnsynchronized(Var v){
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		return map.remove(v)==null ? 0 : 1;
	}

	/**
	 * Consumes a caller-owned local map, merging evidence and installing absent variants.
	 * Existing values are locked for merging; absent values are grouped by shard and
	 * rechecked under that shard's lock. Clears mapT before deferred insertions finish.
	 * The operation is not transactional if a later merge/insertion fails. Newly stored
	 * objects are retained, so callers must stop mutating them after transfer. Var.add
	 * merges read evidence, not coverage, forced flags or all cached annotations.
	 * @param mapT Non-null local map, exclusively owned during this call
	 * @return Number of newly installed identities, excluding merges
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
			}else{
				synchronized(old){
					old.add(v); // Safe to merge immediately
				}
			}
		}

		mapT.clear(); // Release local-map references; deferred variants are retained above

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
	 * Processes one shard at a time: composite filtering, insertion-bias revision,
	 * then final filtering/statistics. Completes all three stages on a shard before
	 * moving to the next; unlike MT, this is not a whole-map barrier between stages.
	 * Requires exclusive processing and initialized dataset/scaffold state. Mutates or
	 * removes stored variants; provided statistics arrays accumulate rather than reset.
	 * @param filter Non-null filter; forced variants bypass emission rejection
	 * @param net Optional network used in the final stage without making a local copy
	 * @param scoreArray Optional overall/type score histograms; the final bin is overflow
	 * @param ploidyArray Optional copy-count histogram
	 * @param avgQualityArray Optional overall/type mean-quality histograms
	 * @param maxQualityArray Optional maximum-quality histogram
	 * @param ADArray Optional two-row allele/total-coverage sums by type
	 * @param AFArray Optional allele-fraction sums by type
	 * @return Type counts from the final stage of each shard
	 */
	public long[] processVariantsST(VarFilter filter, CellNet net, long[][] scoreArray, long[] ploidyArray, long[][] avgQualityArray,
			long[] maxQualityArray, long[][] ADArray, double[] AFArray){
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);

		long[] types=new long[Var.VAR_TYPES];
		for(ConcurrentHashMap<Var, Var> map : maps){
			//Per-shard passes: composite filtering, insertion revision, final filtering/statistics
			long[] types2=processVariants(map, filter, null, null, null, null, null, null, null, false, false);
			types2=processVariants(map, filter, null, null, null, null, null, null, null, true, false);
			types2=processVariants(map, filter, net, scoreArray, ploidyArray, avgQualityArray, maxQualityArray, ADArray, AFArray, false, false);
			Tools.add(types, types2);
		}
		return types;
	}

	/**
	 * Runs three whole-map passes: composite filtering, insertion-bias revision, then
	 * final filtering/statistics with the optional network. Each pass starts eight workers
	 * and waits for completion before the next. Callers must finish accumulation first.
	 * Histogram parameters use processVariantsMT_inner's optional-buffer shape contract.
	 * @param filter Non-null emission filter
	 * @param net Optional final-stage network, copied per worker
	 * @param scoreArray Optional overall/type score histograms
	 * @param ploidyArray Optional copy-count histogram
	 * @param avgQualityArray Optional mean-quality histograms
	 * @param maxQualityArray Optional maximum-quality histogram
	 * @param ADArray Optional allele/coverage sums by type
	 * @param AFArray Optional AF sums by type
	 * @return Final-pass counts by variant type
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
	 * Starts one worker per shard and adds its results to the supplied output arrays.
	 * Joins retry after printing InterruptedException; interruption is not restored.
	 * For valid output buffers, all workers are joined and unsuccessful workers cause
	 * a RuntimeException. Variants and arrays can already be partially modified on failure.
	 *
	 * Workers mirror requested destination histogram widths locally and merge after
	 * the shard pass. Depth/AF scratch buffers are allocated when either ADArray or
	 * AFArray is requested. Score overflow lands in the final requested score bin;
	 * quality indices must fit their requested ranges. Values are added to destinations
	 * without clearing them. Configuration must stay stable during the pass.
	 * @param filter Shared read-only filter during processing
	 * @param net Optional network copied separately by each worker
	 * @param scoreArray Optional overall/type score counts
	 * @param ploidyArray Optional copy-count histogram
	 * @param avgQualityArray Optional mean-quality counts
	 * @param maxQualityArray Optional maximum-quality counts
	 * @param ADArray Optional allele/coverage sums by type
	 * @param AFArray Optional AF sums by type
	 * @param processInsertions Select insertion-bias revision instead of filtering/statistics
	 * @return Sum of worker type counts; insertion passes return zero counts
	 */
	private long[] processVariantsMT_inner(VarFilter filter, CellNet net, long[][] scoreArray, long[] ploidyArray,
			long[][] avgQualityArray, long[] maxQualityArray, long[][] ADArray, double[] AFArray, boolean processInsertions){
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);

		//Create processing threads for each shard
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(WAYS);
		for(int i=0; i<WAYS; i++){
			ProcessThread pt=new ProcessThread(maps[i], filter, net,
					scoreArray==null ? 0 : scoreArray[0].length,
					ploidyArray==null ? 0 : ploidyArray.length,
					avgQualityArray==null ? 0 : avgQualityArray[0].length,
					maxQualityArray==null ? 0 : maxQualityArray.length,
					ADArray!=null || AFArray!=null, processInsertions);
			alpt.add(pt);
			pt.start();
		}

		//Collect results from all threads
		long[] types=new long[Var.VAR_TYPES];
		boolean success=true;
		for(ProcessThread pt : alpt){
			//Wait for thread completion
			while(pt.getState()!=Thread.State.TERMINATED){
				try{
					pt.join();
				}catch(InterruptedException e){
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
			if(AFArray!=null){Tools.add(AFArray, pt.AFArray);}
			success&=pt.success;
		}
		if(!success){
			Throwable failure=null;
			for(ProcessThread pt : alpt){
				if(pt.failure!=null){
					failure=pt.failure;
					break;
				}
			}
			throw new RuntimeException("VarMap.processVariantsMT_inner: a ProcessThread failed", failure);
		}

		return types;
	}

	/**
	 * Worker for one shard, with local histograms and an optional network evaluation copy.
	 * The shard, filter and enclosing map configuration are borrowed and must stay stable
	 * except for mutations performed by the coordinated processing passes.
	 */
	private class ProcessThread extends Thread{

		/**
		 * Stores borrowed inputs and allocates the enabled local histograms.
		 * Histogram widths mirror requested destinations and merge after the pass.
		 * @param map_ Shard exclusively assigned to this worker for the pass
		 * @param filter_ Filter shared without concurrent mutation
		 * @param net0_ Optional network to copy inside run
		 * @param scoreBins Score histogram width, or zero to disable score histograms
		 * @param ploidyBins Copy-count histogram width, or zero to disable ploidy histograms
		 * @param avgQualityBins Average-quality histogram width, or zero to disable it
		 * @param maxQualityBins Maximum-quality histogram width, or zero to disable it
		 * @param trackAD Allocate both depth and allele-fraction arrays
		 * @param processInsertions_ Select insertion-bias revision
		 */
		ProcessThread(Map<Var, Var> map_, VarFilter filter_, CellNet net0_, int scoreBins, int ploidyBins,
				int avgQualityBins, int maxQualityBins, boolean trackAD, boolean processInsertions_){
			map=map_;
			filter=filter_;
			net0=net0_;
			scoreArray=(scoreBins>0 ? new long[Var.VAR_TYPES+1][scoreBins] : null);
			ploidyArray=(ploidyBins>0 ? new long[ploidyBins] : null);
			avgQualityArray=(avgQualityBins>0 ? new long[Var.VAR_TYPES+1][avgQualityBins] : null);
			maxQualityArray=(maxQualityBins>0 ? new long[maxQualityBins] : null);
			ADArray=(trackAD ? new long[2][Var.VAR_TYPES] : null);
			AFArray=(trackAD ? new double[Var.VAR_TYPES] : null);
			processInsertions=processInsertions_;
		}

		/**
		 * Copies net0 with copy(false), processes the shard, then marks success.
		 * An exception or error leaves success false; there is no rollback of earlier mutations.
		 */
		@Override
		public void run(){
			try{
				net=(net0==null ? null : net0.copy(false));
				types=processVariants(map, filter, net, scoreArray, ploidyArray, avgQualityArray, maxQualityArray,
						ADArray, AFArray, processInsertions, false);
				success=true;
			}catch(Throwable t){
				failure=t;
			}
		}

		/**Quality filter for variant assessment*/
		final VarFilter filter;
		/**Shard to process*/
		final Map<Var, Var> map;
		/**CellNet instance for machine learning predictions*/
		final CellNet net0;
		/**Thread-local copy of CellNet*/
		CellNet net;
		/**Failure thrown by this shard, if any*/
		Throwable failure;
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
	 * Runs one final filtering/statistics pass over the current survivors.
	 * Used after composite filtering and nearby counting; does not perform either step
	 * itself. Nearby counts are available as NN features, but considerNearby is false,
	 * so hard nearby rejection remains caller-owned. Uses the same worker/histogram
	 * optional-buffer contract as processVariantsMT. A null net selects ordinary
	 * composite scoring.
	 * @param filter Non-null emission filter
	 * @param net Optional network copied per worker
	 * @param scoreArray Optional overall/type score histograms
	 * @param ploidyArray Optional copy-count histogram
	 * @param avgQualityArray Optional mean-quality histograms
	 * @param maxQualityArray Optional maximum-quality histogram
	 * @param ADArray Optional allele/coverage sums by type
	 * @param AFArray Optional AF sums by type
	 * @return Counts of surviving variants by type
	 */
	public long[] rescoreWithNetMT(VarFilter filter, CellNet net, long[][] scoreArray, long[] ploidyArray,
			long[][] avgQualityArray, long[] maxQualityArray, long[][] ADArray, double[] AFArray){
		return processVariantsMT_inner(filter, net, scoreArray, ploidyArray, avgQualityArray, maxQualityArray,
				ADArray, AFArray, false);
	}

	/**
	 * Processes one shard exclusively. In insertion mode, revises insertion AF and
	 * may update nearby substitutions in other shards through Var.reviseAlleleFraction.
	 * Otherwise applies fast checks, fills coverage, bypasses full filtering for forced
	 * variants, removes failures, and accumulates statistics for survivors. Score may
	 * be calculated once for filtering and again for its histogram. Scores overflow
	 * into the final score bin; quality stats with support must fit their bounded
	 * histograms. Zero-support survivors have undefined average base quality and are
	 * counted in the Q0 average-quality bin for legacy histogram compatibility. All
	 * public processing entry points currently pass considerNearby=false.
	 * @param map Shard with no concurrent accumulation/removal
	 * @param filter Non-null emission filter
	 * @param net Optional scoring network, owned by this invocation/worker
	 * @param scoreArray Optional score histograms; the final bin is overflow
	 * @param ploidyArray Optional copy-count histogram
	 * @param avgQualityArray Optional bounded mean-quality histograms
	 * @param maxQualityArray Optional bounded maximum-quality histogram
	 * @param ADArray Optional allele/total-coverage sums by type
	 * @param AFArray Optional AF sums by type
	 * @param processInsertions Select bias revision; no statistics or removals in this mode
	 * @param considerNearby Whether full filtering may enforce the nearby threshold
	 * @return Type counts for survivors, or zeros in insertion mode
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
						final int score=overflowHistogramIndex(v.phredScore(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg,
								filter.rarity, ploidy, scafMap, net), scoreArray[0]);
						scoreArray[0][score]++; //Overall scores
						scoreArray[v.type()+1][score]++; //Type-specific scores
					}

					//Collect ploidy statistics
					if(ploidyArray!=null){ploidyArray[v.calcCopies(ploidy)]++;}

					//Collect quality statistics
					if(avgQualityArray!=null){
						final int q=boundedHistogramIndex(v.baseQAvg(), avgQualityArray[0], "average base quality", v.alleleCount()==0);
						avgQualityArray[0][q]++; // Overall quality
						avgQualityArray[v.type()+1][q]++; // Type-specific quality
					}
					if(maxQualityArray!=null){
						maxQualityArray[boundedHistogramIndex(v.baseQMax, maxQualityArray, "maximum base quality", false)]++;
					}

					//Collect depth and frequency statistics
					if(ADArray!=null){
						ADArray[0][v.type()]+=v.alleleCount(); // Allele depth by type
						ADArray[1][v.type()]+=v.coverage(); // Total coverage by type
					}
					if(AFArray!=null){AFArray[v.type()]+=v.alleleFraction();}
				}else{
					iterator.remove(); // Remove variants that fail filtering
				}
			}
		}
		return types;
	}

	/** Returns the raw score bin, with the final slot reserved for overflow. */
	private static int overflowHistogramIndex(final double value, final long[] array){
		return Tools.mid(0, (int)value, array.length-1);
	}

	/**
	 * Returns a direct histogram index for values with a finite valid domain.
	 * A zero-support average has an undefined NaN mean and is counted in bin 0.
	 */
	private static int boundedHistogramIndex(final double value, final long[] array, final String label, final boolean allowNaNZero){
		if(allowNaNZero && Double.isNaN(value)){return 0;}
		if(value>=0 && value<array.length){return (int)value;}
		throw new IllegalArgumentException(label+" histogram value "+value
				+" is outside [0,"+(array.length-1)+"]");
	}

	/**
	 * Unused helper that installs key-only Var copies for absent shared identities.
	 * The Var copy constructor shares allele bytes and does not copy accumulated evidence.
	 * Then calculates coverage on the source sharedMap variants and counts all source
	 * variants, including identities already present in the target.
	 * @param map Target shard, exclusively owned
	 * @param sharedMap Source identities; their coverage caches may be updated
	 * @return Type counts for all sharedMap variants, not merely the new insertions
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

		//Count all shared source variants by type
		long[] types=new long[Var.VAR_TYPES];
		for(Var v : sharedMap.keySet()){
			v.calcCoverage(scafMap);
			types[v.type()]++;
		}
		return types;
	}

	/**
	 * Allocates a reference snapshot using size(), then optionally sorts by Var.compareTo.
	 * Variants are borrowed, not copied. Requires stable contents and a size fitting an
	 * int array; concurrent growth can overrun allocation and shrinkage can leave nulls.
	 * @param sort Whether to sort by scaffold, adjusted position and remaining identity keys
	 * @return New array of references to stored variants
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
	 * Legacy assertion-based diagnostic skeleton with an unconditional failing entry
	 * assertion. With assertions disabled, its invariant checks are disabled too; it is
	 * not a callable production validator. Optional printing still walks variant fields.
	 * @param quiet Suppress per-variant debug output
	 * @return True if execution reaches the end
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
	 * Fills missing coverage caches through Var.calcCoverage and counts variants by type.
	 * Already cached coverage is retained even if a different scaffold map is supplied.
	 * Requires coordinated access to mutable variants and scaffold coverage.
	 * @param scafMap Borrowed map for uncached coverage
	 * @return New type-count array
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
	 * Counts current variants by type without filtering or filling coverage.
	 * Concurrent mutation prevents a whole-map snapshot guarantee.
	 * @return New type-count array
	 */
	public long[] countTypes(){
		long[] types=new long[Var.VAR_TYPES];
		for(Var v : this){
			types[v.type()]++;
		}
		return types;
	}

	/**
	 * Replaces all shards and resets dataset rates/means to -1. Retains ploidy and scafMap.
	 * Requires exclusive access: existing iterators/borrowed values are not invalidated
	 * or cleared in place, and outstanding users can still reference the old shards.
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
	 * Formats every variant through toTextQuick, in unspecified shard/map order,
	 * with a newline after each and no header. Uses that helper's fixed scoring defaults,
	 * not this map's dataset settings. The delegate enables global Var.useIdentity and
	 * can update variant caches; this is not a side-effect-free snapshot.
	 * @return Newly allocated diagnostic text
	 */
	@Override
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
	 * Creates a traversal over borrowed values, visiting shards in index order.
	 * Each shard uses its map iterator; there is no whole-map snapshot, ordering by
	 * variant position, or coordination with clear. Require stable shards for a complete traversal.
	 * @return Iterator over stored mutable variants
	 */
	@Override
	public VarMapIterator iterator(){
		return new VarMapIterator();
	}

	/**
	 * Traverses shard entry iterators and returns their borrowed values.
	 * Shards are entered lazily, so concurrent writes/reset can change what is observed.
	 * No remove override is provided.
	 */
	private class VarMapIterator implements Iterator<Var>{

		/**
		 * Advances to the first nonempty shard, or the last empty iterator.
		 */
		VarMapIterator(){
			makeReady(); // Initialize to first non-empty shard
		}

		/**
		 * Tests the current shard iterator. Shard advancement occurs during construction
		 * and after next consumes the current shard's last observed entry.
		 * @return Whether the current iterator has another entry
		 */
		@Override
		public boolean hasNext(){
			return iter.hasNext();
		}

		/**
		 * Consumes the current entry and prepares the next nonempty shard if needed.
		 * Delegates exhaustion behavior to the current entry iterator.
		 * @return Borrowed variant from the consumed entry
		 */
		@Override
		public Var next(){
			Entry<Var, Var> e=iter.next();
			if(!iter.hasNext()){makeReady();} // Prepare next shard if current is exhausted
			Var v=e==null ? null : e.getValue(); // Extract variant from map entry
			return v;
		}

		/**
		 * Skips exhausted/empty shards, retaining the final iterator when all are exhausted.
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
	/**Fixed shard-array object; clear replaces its mutable map elements.*/
	final ConcurrentHashMap<Var, Var>[] maps;

	/*--------------------------------------------------------------*/
	/*----------------        Static fields         ----------------*/
	/*--------------------------------------------------------------*/

	/**Legacy unused identifier; this class does not implement Serializable.*/
	private static final long serialVersionUID=1L;
	/** Number of hash map shards (must be power of 2 for efficient masking) */
	private static final int WAYS=8;
	/**Bit mask for shard selection (WAYS-1)*/
	public static final int MASK=WAYS-1;

}
