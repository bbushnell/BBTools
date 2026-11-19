package prot.cluster;

import java.math.BigDecimal;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;

import ddl.FamilySketchBuilder;
import ddl.MutableFamilyIndex;
import ddl.SnapshotFamilyIndex;
import fileIO.ByteStreamWriter;
import prot.ProteinSequence;
import prot.HbmPositionModel;
import template.Accumulator;
import template.ThreadWaiter;
import structures.IntList;
import structures.ByteBuilder;

/** Serial reference or explicit asynchronous greedy admission, followed by frozen passes.
 * Input ordinals and family IDs are stable within a run. Concurrent ranking uses an
 * immutable prefix plus a current tail; model versions can advance during a query.
 * Screens precede alignment, then the best passing shortlist candidate wins.
 * @author Keqing
 */
public final class GreedyClusterer {

	public GreedyClusterer(ProteinSequence[] sequences, int minimumLength, int shortlist,
			int minimumCluster, int passes, int topN, long seed, FamilySketchBuilder sketches,
			TailAlignmentGate gate, ConsensusUpdater updater, ClusteringMetrics metrics){
		this(sequences, minimumLength, shortlist, minimumCluster, passes, topN, seed, sketches,
			gate, updater, metrics, false);
	}

	/** Indexed ranking is exact; false retains the exhaustive reference implementation. */
	public GreedyClusterer(ProteinSequence[] sequences, int minimumLength, int shortlist,
			int minimumCluster, int passes, int topN, long seed, FamilySketchBuilder sketches,
			TailAlignmentGate gate, ConsensusUpdater updater, ClusteringMetrics metrics, boolean indexed){
		this(sequences,minimumLength,shortlist,minimumCluster,passes,topN,seed,sketches,gate,updater,metrics,indexed,1);
	}

	/** Greedy remains serial; only frozen queries and independent family rebuilds use workers. */
	public GreedyClusterer(ProteinSequence[] sequences, int minimumLength, int shortlist,
			int minimumCluster, int passes, int topN, long seed, FamilySketchBuilder sketches,
			TailAlignmentGate gate, ConsensusUpdater updater, ClusteringMetrics metrics, boolean indexed,int threads){
		this(sequences,minimumLength,shortlist,minimumCluster,passes,topN,seed,sketches,gate,updater,metrics,indexed,threads,null);
	}

	/** Optional histograms observe decisions before mutation; they never influence selection. */
	public GreedyClusterer(ProteinSequence[] sequences, int minimumLength, int shortlist,
			int minimumCluster, int passes, int topN, long seed, FamilySketchBuilder sketches,
			TailAlignmentGate gate, ConsensusUpdater updater, ClusteringMetrics metrics, boolean indexed,int threads,
			AssignmentHistograms histograms){
		this(sequences,minimumLength,shortlist,minimumCluster,passes,topN,seed,sketches,gate,updater,metrics,indexed,threads,histograms,CandidateFilter.disabled());
	}

	/** Screening changes alignment eligibility only; the ranked shortlist stays unchanged. */
	public GreedyClusterer(ProteinSequence[] sequences, int minimumLength, int shortlist,
			int minimumCluster, int passes, int topN, long seed, FamilySketchBuilder sketches,
			TailAlignmentGate gate, ConsensusUpdater updater, ClusteringMetrics metrics, boolean indexed,int threads,
			AssignmentHistograms histograms,CandidateFilter filter){
		this(sequences,minimumLength,shortlist,minimumCluster,passes,topN,seed,sketches,
			gate,updater,metrics,indexed,threads,histograms,filter,null);
	}

	/** An optional serial audit records the final frozen pass without changing decisions. */
	public GreedyClusterer(ProteinSequence[] sequences, int minimumLength, int shortlist,
			int minimumCluster, int passes, int topN, long seed, FamilySketchBuilder sketches,
			TailAlignmentGate gate, ConsensusUpdater updater, ClusteringMetrics metrics, boolean indexed,int threads,
			AssignmentHistograms histograms,CandidateFilter filter,AssignmentDecisionAudit audit){
		this(sequences,minimumLength,shortlist,minimumCluster,passes,topN,seed,sketches,
			gate,updater,metrics,indexed,threads,histograms,filter,audit,1);
	}

	/** Optional parallel immutable construction for fixed seeds; ordinary constructors retain mutable indexing. */
	public GreedyClusterer(ProteinSequence[] sequences, int minimumLength, int shortlist,
			int minimumCluster, int passes, int topN, long seed, FamilySketchBuilder sketches,
			TailAlignmentGate gate, ConsensusUpdater updater, ClusteringMetrics metrics, boolean indexed,int threads,
			AssignmentHistograms histograms,CandidateFilter filter,AssignmentDecisionAudit audit,int frozenIndexThreads){
		if(sequences==null || sketches==null || gate==null || updater==null || metrics==null ||
			filter==null || minimumLength<gate.minSpan || shortlist<1 || minimumCluster<1 || passes<0 || topN<0 || threads<1 ||
			frozenIndexThreads<1 || frozenIndexThreads>threads || frozenIndexThreads>1 && !indexed){
			throw new IllegalArgumentException("Invalid clustering inputs, components or bounds");
		}
		if(histograms!=null){histograms.validateSchema(passes,sketches.width(),shortlist);}
		if(audit!=null){audit.validateSchema(passes,shortlist,sequences.length,threads);}
		filter.validateWidth(sketches.width());
		this.sequences=sequences; this.minimumLength=minimumLength; this.shortlist=shortlist;
		this.minimumCluster=minimumCluster; this.passes=passes; this.topN=topN; this.seed=seed;
		this.sketches=sketches; this.gate=gate; this.updater=updater; this.metrics=metrics;
		scorer=new HbmPositionModel.ConsensusScorer(gate.matrix);
		this.threads=threads; this.histograms=histograms; this.filter=filter; this.audit=audit;
		this.frozenIndexThreads=frozenIndexThreads;
		greedyThreads=Math.min(threads,16);
		index=new MutableFamilyIndex(sketches.width(), indexed && frozenIndexThreads==1);
		metrics.attachIndex(index);
		querySketch=new char[sketches.width()]; candidateIds=new int[shortlist];
		assignments=new int[sequences.length]; Arrays.fill(assignments, -1);
		int count=0;
		for(int i=0; i<sequences.length; i++){
			if(sequences[i]==null){throw new IllegalArgumentException("Missing sequence at input ordinal "+i);}
			if(sequences[i].length()>=minimumLength){count++;}
		}
		order=new int[count];
		for(int i=0, j=0; i<sequences.length; i++){
			if(sequences[i].length()>=minimumLength){order[j++]=i;}
		}
		final Random random=new Random(seed);
		for(int i=order.length-1; i>0; i--){
			final int j=random.nextInt(i+1), old=order[i]; order[i]=order[j]; order[j]=old;
		}
	}

	/** Enables the released member pivot before seed import or greedy publication.
	 * Existing constructors retain legacy ranking. The cache must belong to this
	 * exact ordered population; seed descriptors must be attached to codec entries.
	 * Consensus-input merging needs a separate original-member binding and is rejected.
	 */
	public GreedyClusterer memberRanking(final MemberRankInputs inputs, final int layer){
		if(started || initialSeeds!=null || memberInputs!=null || inputs==null || !inputs.belongsTo(sequences) ||
			layer<0 || merge || mergeInputWeights!=null){
			throw new IllegalArgumentException("Attach member ranking once before seeds/execution to its original population; consensus-input merging is not bound");
		}
		memberInputs=inputs; memberLayer=layer;
		memberQuery=new MemberRanking(this){
			@Override Family family(final int id){return families.get(id);}
			@Override public double score(final int id){return cosine.score(id);}
		};
		return this;
	}

	/** Selects the explicit asynchronous mode before execution; legacy constructors stay serial. */
	public GreedyClusterer concurrentGreedy(double fraction){
		if(started || initialSeeds!=null || !index.indexed() || audit!=null || greedyObserver!=null || !(updater instanceof ConsensusRefiner) ||
			gate.mode==TailAlignmentGate.Mode.FULL_QUERY_UNPINNED || !Double.isFinite(fraction) || fraction<=0){
			throw new IllegalArgumentException("Concurrent greedy requires an unstarted indexed right-pinned engine, native refiner, positive snapshot fraction and no serial audit");
		}
		concurrentGreedy=true; snapshotFraction=fraction;
		return this;
	}

	/** Attaches a serial decision-time observer; caller owns and closes its writers. */
	public GreedyClusterer greedyAudit(GreedyDecisionAudit observer){
		if(initialSeeds!=null){throw new IllegalArgumentException("Fixed seeds have no initial greedy audit phase");}
		if(merge){throw new IllegalArgumentException("Interpass merging does not support decision audits");}
		if(scavenge){throw new IllegalArgumentException("Scavenging does not support greedyaudit");}
		if(started || observer==null || greedyObserver!=null || filter.greedyHits>4){
			throw new IllegalStateException("Attach one greedy audit before execution with a supported hit floor");
		}
		observer.validateSerial(concurrentGreedy);
		greedyObserver=observer;
		return this;
	}

	/** Retries initial singletons before pruning; unsuccessful tails remain unassigned.
	 * Pair histograms and the initial greedy observer do not describe this extra phase.
	 */
	public GreedyClusterer scavenge(final boolean enabled){
		if(started){throw new IllegalStateException("Select scavenging before execution");}
		if(enabled && (initialSeeds!=null || histograms!=null || greedyObserver!=null)){
			throw new IllegalArgumentException("Scavenging does not support pairhist or greedyaudit");
		}
		scavenge=enabled;
		return this;
	}

	/** Optional progress-dependent greedy floor; no tracker exists on the default path. */
	public GreedyClusterer greedyRamp(int targetFamilies,int minimumMembers){
		if(started || initialSeeds!=null || merge || ramp!=null || filter.greedyHits>=0){
			throw new IllegalStateException("Select the greedy ramp once before execution, with an occupancy-based filter");
		}
		ramp=new GreedyHitRamp(order.length,targetFamilies,minimumMembers);
		return this;
	}

	/** Starts from exported models, without greedy admission or invented original memberships.
	 * Only ordinary all-input reassignment/prune/rebuild passes follow seed loading.
	 */
	public GreedyClusterer fixedSeeds(ClusterSketchCodec.Loaded loaded){
		return fixedSeeds(loaded,0);
	}

	/** Restricts imported references before any assignment or index insertion.
	 * Exported member counts rank seeds; they never invent current input memberships.
	 * Zero retains the original order and behavior; ties use ascending stable ID.
	 */
	public GreedyClusterer fixedSeeds(ClusterSketchCodec.Loaded loaded,final int limit){
		return fixedSeeds(loaded,limit,null);
	}

	/** Explicit exclusions precede source-size ranking and require a fully populated positive limit.
	 * Null/empty exclusions preserve the original selector, including its order and short-library behavior.
	 */
	public GreedyClusterer fixedSeeds(ClusterSketchCodec.Loaded loaded,final int limit,final int[] excludedIds){
		return fixedSeeds(loaded,limit,excludedIds,false);
	}

	/** Explicit full-roster continuation admits empty serialized templates before import.
	 * The old overloads keep their positive-source-count contract.
	 */
	public GreedyClusterer fixedSeeds(ClusterSketchCodec.Loaded loaded,final int limit,final int[] excludedIds,
			final boolean preserveAll){
		if(limit<0){throw new IllegalArgumentException("seedlimit must be nonnegative");}
		if(preserveAll && (limit!=0 || excludedIds!=null && excludedIds.length>0)){
			throw new IllegalArgumentException("Full-roster continuation cannot limit or exclude imported seeds");
		}
		if(started || initialSeeds!=null || loaded==null || passes<1 || concurrentGreedy || scavenge ||
			ramp!=null || histograms!=null || audit!=null || greedyObserver!=null){
			throw new IllegalArgumentException("Fixed seeds require passes>=1 and no greedy mode, scavenging, ramp or decision/pair audit");
		}
		final ClusterSketchCodec.Header h=loaded.header;
		if(!sketches.matchesGeometry(h.width,h.geometry,h.sketchSeed,h.positionMode,h.partitions,h.overlap,h.blacklistSha80)){
			throw new IllegalArgumentException("Fixed-seed sketch geometry or blacklist differs from query sketches");
		}
		final java.util.HashSet<String> names=new java.util.HashSet<String>();
		final structures.IntHashMap ids=new structures.IntHashMap();
		for(int i=0; i<loaded.families.length; i++){
			final ClusterSketchCodec.Entry entry=loaded.families[i];
			if(entry==null || entry.familyId==Integer.MAX_VALUE || entry.consensus.length<gate.minSpan ||
				ids.get(entry.familyId)>=0 || !names.add(entry.name)){
				throw new IllegalArgumentException("Fixed seeds need unique bounded IDs and consensuses at least minspan long");
			}
			ids.put(entry.familyId,1);
		}
		final boolean exclude=excludedIds!=null && excludedIds.length>0;
		final int count;
		final ClusterSketchCodec.Entry[] selected;
		if(exclude){
			if(limit<=0){throw new IllegalArgumentException("seedexclude requires positive seedlimit");}
			final structures.IntHashMap excluded=new structures.IntHashMap(excludedIds.length);
			for(int id:excludedIds){
				if(id<0 || id==Integer.MAX_VALUE || excluded.get(id)>=0){
					throw new IllegalArgumentException("seedexclude needs distinct bounded nonnegative IDs");
				}
				if(ids.get(id)<0){throw new IllegalArgumentException("Unknown seedexclude ID: "+id);}
				excluded.put(id,1);
			}
			final int survivors=loaded.families.length-excludedIds.length;
			if(survivors<limit){throw new IllegalArgumentException("seedexclude leaves fewer seeds than seedlimit");}
			selected=new ClusterSketchCodec.Entry[survivors];
			int next=0;
			for(ClusterSketchCodec.Entry entry:loaded.families){
				if(excluded.get(entry.familyId)<0){selected[next++]=entry;}
			}
			assert next==survivors : "Every distinct excluded ID was present in the validated source roster";
			Arrays.sort(selected,SEED_SIZE_ORDER);
			count=limit;
		}else{
			count=limit==0 ? loaded.families.length : Math.min(limit,loaded.families.length);
			if(count<loaded.families.length){
				selected=loaded.families.clone();
				Arrays.sort(selected,SEED_SIZE_ORDER);
			}else{selected=loaded.families;}
		}
		final Family[] models=new Family[count];
		for(int i=0; i<count; i++){
			models[i]=Family.fromSeed(selected[i],preserveAll);
			if(memberInputs!=null){models[i].enableMemberRanking(memberInputs,memberLayer,selected[i].memberRanking);}
		}
		initialSeeds=models;
		if(preserveAll){preserveAllSeeds(true);}
		return this;
	}

	/** Records the actual selected models before assignment consumes source counts or prunes empty references.
	 * Rows retain selector order; member counts belong to the immutable source codec, not this run.
	 */
	public void writeSelectedSeeds(final String path,final boolean overwrite){
		if(started || initialSeeds==null || path==null || path.isEmpty()){
			throw new IllegalStateException("Selected-seed roster requires fixed seeds before execution and an output path");
		}
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,true);
		writer.start();
		try{
			writer.println("family_id\tsource_members");
			final ByteBuilder row=new ByteBuilder(64);
			for(Family family:initialSeeds){
				assert(family.size()==0) : "Selected-seed roster precedes assignment; source counts must not be replaced by current memberships";
				row.clear().append(family.familyId).tab().append(family.referenceSize()); writer.println(row);
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Selected-seed roster output failed: "+path);}}
	}

	/** One-pass evaluation against fixed models can omit terminal model work, including after seed selection. */
	public GreedyClusterer skipFinalRebuild(final boolean enabled){
		if(started){throw new IllegalStateException("Select final-rebuild behavior before execution");}
		if(enabled){validateFinalRebuildSkip();}
		skipFinalRebuild=enabled; return this;
	}

	private void validateFinalRebuildSkip(){
		if(initialSeeds==null || passes!=1 || minimumCluster!=1 || topN!=0 || merge ||
			consensusQuality!=null || !(updater instanceof ConsensusRefiner)){
			throw new IllegalArgumentException("Final-rebuild skip requires fixed seeds, one pass, mincluster1, topn0, no merge/quality and native refinement");
		}
	}

	/** Opt-in observational accounting; never changes pass count or pruning policy. */
	public GreedyClusterer observePassChanges(){
		if(started || initialSeeds==null || merge || topN!=0 || passChanges!=null || !(updater instanceof ConsensusRefiner)){
			throw new IllegalArgumentException("passchanges requires fixed seeds, merge=f, topn=0 and a native refiner before execution");
		}
		previousAssignments=new int[assignments.length];
		passChanges=new ArrayList<long[]>();
		return this;
	}

	/** Stops after completed refinement pass 2+ when fewer than this fraction of all input tails changed winners. */
	public GreedyClusterer stopChangeFraction(final double fraction){
		if(started || !Double.isFinite(fraction) || fraction<0 || fraction>=1){
			throw new IllegalArgumentException("stopchangefraction must be finite and in [0,1) before execution");
		}
		stopChangeFraction=fraction;
		if(fraction>0 && previousAssignments==null){previousAssignments=new int[assignments.length];}
		return this;
	}

	/** Exact decimal strictness prevents a binary product from turning equality into convergence. */
	public static boolean belowStopFraction(final long changed,final int inputs,final double fraction){
		if(inputs<0 || changed<0 || changed>inputs || !Double.isFinite(fraction) || fraction<0 || fraction>=1){
			throw new IllegalArgumentException("Stop comparison requires valid original-input counts and a finite fraction in [0,1)");
		}
		return BigDecimal.valueOf(changed).compareTo(BigDecimal.valueOf(fraction).multiply(BigDecimal.valueOf(inputs)))<0;
	}

	/** Emits only completed-run observations, separate from legacy checkpoint/stat schemas.
	 * Pass1 compares against all-unassigned inputs, not unavailable seed member IDs.
	 * Consensus changes compare exact start/end bytes for surviving stable IDs;
	 * removals are separate. Settled is advisory and requires a previous completed pass.
	 */
	public void writePassChanges(final String path,final boolean overwrite){
		if(!completedRun || passChanges==null || passChanges.size()!=completedPasses){
			throw new IllegalStateException("Pass changes require every completed pass of a completed observed run");
		}
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,true);
		writer.start();
		try{
			writer.println("pass\teligible\tassigned_after_reassign\tassigned_after_prune\tunassigned_after_prune\tchanged_assignments_reassign\tchanged_assignments_end\tchanged_consensuses\tremoved_references\tpruned_assignments\tsettled");
			final ByteBuilder row=new ByteBuilder();
			for(long[] values:passChanges){
				row.clear();
				for(int i=0; i<values.length; i++){if(i>0){row.tab();} row.append(values[i]);}
				writer.println(row);
			}
		}finally{
			if(writer.poisonAndWait()){throw new IllegalStateException("Pass-change report failed: "+path);}
		}
	}

	private void beginPassChanges(){
		System.arraycopy(assignments,0,previousAssignments,0,assignments.length);
		if(passChanges==null){return;}
		if(previousConsensuses==null){previousConsensuses=new byte[families.size()][];}
		for(int id=0; id<families.size(); id++){
			final Family family=families.get(id);
			// Native Family.setConsensus replaces/copies arrays; old payloads are immutable.
			previousConsensuses[id]=family==null ? null : family.consensusUnsafe();
		}
	}

	private long changedAssignments(){
		long changed=0;
		for(int input:order){if(assignments[input]!=previousAssignments[input]){changed++;}}
		return changed;
	}

	private void endPassChanges(final int pass,final long afterReassign,final long changedReassign,final long changedEnd){
		long changedConsensuses=0,removed=0;
		for(int id=0; id<previousConsensuses.length; id++){
			final byte[] before=previousConsensuses[id];
			if(before==null){continue;}
			final Family family=families.get(id);
			if(family==null){removed++;}
			else if(!Arrays.equals(before,family.consensusUnsafe())){changedConsensuses++;}
			previousConsensuses[id]=null;
		}
		final long pruned=afterReassign-assigned;
		final long settled=pass>1 && changedReassign==0 && changedEnd==0 && changedConsensuses==0 && removed==0 && pruned==0 ? 1 : 0;
		passChanges.add(new long[]{pass,order.length,afterReassign,assigned,order.length-assigned,
			changedReassign,changedEnd,changedConsensuses,removed,pruned,settled});
	}

	private boolean trackPassChanges(){return passChanges!=null || stopChangeFraction>0;}

	/** Defers removal of fixed-seed references, including empty models, until final reassignment. */
	public GreedyClusterer pruneFinal(final boolean enabled){
		if(started){throw new IllegalStateException("Select final-only pruning before execution");}
		pruneFinal=enabled;
		return this;
	}

	/** Keeps every imported reference through one reassignment and optional rebuild, including empty templates.
	 * Configure after fixedSeeds; empty references keep their last model for the next live layer.
	 */
	public GreedyClusterer preserveAllSeeds(final boolean enabled){
		if(started){throw new IllegalStateException("Select seed preservation before execution");}
		if(enabled){validateSeedPreservation();}
		preserveAllSeeds=enabled;
		return this;
	}

	private void validateSeedPreservation(){
		if(initialSeeds==null || passes!=1 || minimumCluster!=1 || topN!=0 || merge ||
			pruneFinal || consensusQuality!=null){
			throw new IllegalArgumentException("Seed preservation requires fixed seeds, one pass, mincluster1, topn0, and no merge, pruning or quality report");
		}
	}

	/** Limits asynchronous greedy workers independently of the deterministic phase worker count. */
	public GreedyClusterer greedyThreads(final int count){
		if(started || count<1 || count>threads){
			throw new IllegalArgumentException("Greedy workers must be selected before execution and lie in [1,t]");
		}
		greedyThreads=count;
		return this;
	}

	/** Stably orders the seeded eligible-input permutation before any worker starts. */
	public GreedyClusterer inputOrder(ClusteringOrder.Mode mode){
		if(started || orderConfigured){throw new IllegalStateException("Select input ordering once before execution");}
		ClusteringOrder.apply(order,sequences,mode); orderConfigured=true;
		return this;
	}

	/** Attaches one final-pass observer; only the driver records completed winner slots. */
	public GreedyClusterer assignmentMetrics(AssignmentMetricSummary observer){
		//TODO: Probable bug - stopchangefraction can finish before the configured final pass,
		// but assignment-metric capture still uses pass==passes in serial/parallel queries.
		// AssignmentMetricSummary.write then rejects retained inputs without recorded metrics.
		// The cake driver uses final-consensus quality instead; this separate observer needs
		// per-pass reset/capture support before combining it with early stopping.
		if(started || merge || assignmentMetrics!=null || observer==null || passes<1 || observer.inputs()!=sequences.length){
			throw new IllegalArgumentException("Assignment metrics require one observer matching all inputs before a run with passes>=1");
		}
		assignmentMetrics=observer;
		return this;
	}

	/** Captures final-consensus cohesion during existing native member alignments only. */
	public GreedyClusterer consensusQuality(final ConsensusQualityReport observer){
		if(started || observer==null || consensusQuality!=null || !(updater instanceof ConsensusRefiner)){
			throw new IllegalArgumentException("Consensus quality requires one observer and a native refiner before execution");
		}
		consensusQuality=observer;
		return this;
	}

	/** Merges similar consensus families after each ordinary pass, never during greedy admission. */
	public GreedyClusterer merge(final boolean enabled){
		return merge(enabled,gate);
	}

	/** Uses an optional merge-only threshold while member assignment and union validation keep the ordinary gate. */
	public GreedyClusterer merge(final boolean enabled,final TailAlignmentGate nestedGate){
		if(started){throw new IllegalStateException("Select interpass merging before execution");}
		if(enabled && memberInputs!=null){
			throw new IllegalArgumentException("Member ranking requires original-member descriptor binding before consensus-input interpass merging");
		}
		if(enabled && (passes<1 || !(updater instanceof ConsensusRefiner) ||
			gate.mode==TailAlignmentGate.Mode.FULL_QUERY_UNPINNED || histograms!=null || audit!=null ||
			greedyObserver!=null || assignmentMetrics!=null || ramp!=null)){
			throw new IllegalArgumentException("Interpass merging requires passes>=1, a native pinned refiner and no pairhist, audits, assignmetrics or ramp");
		}
		if(enabled){checkMergeGate(nestedGate);}
		merge=enabled;
		mergeGate=nestedGate;
		if(enabled && mergeSummaries==null){mergeSummaries=new ArrayList<MergeSummary>();}
		return this;
	}

	/** Reports after final pruning/top-N selection without exposing mutable assignment storage. */
	public void writeAssignmentMetrics(String path,boolean overwrite){
		if(!completedRun || assignmentMetrics==null){
			throw new IllegalStateException("Assignment metrics require a completed observed run");
		}
		assignmentMetrics.write(path,overwrite,assignments);
	}

	/** Executes once. The returned result owns the completed memberships and assignment array. */
	public Result run(){
		if(started){throw new IllegalStateException("A reference clusterer cannot run twice");}
		if(stopChangeFraction>0 && (initialSeeds==null || merge || topN!=0 || pruneFinal || skipFinalRebuild || !(updater instanceof ConsensusRefiner))){
			throw new IllegalArgumentException("stopchangefraction requires fixed seeds, ordinary rebuilds, merge=false, topn=0 and prunefinal=false");
		}
		if(frozenIndexThreads>1 && (initialSeeds==null || merge)){
			throw new IllegalArgumentException("Parallel frozen indexing requires fixed seeds and merge=false; greedy indexing stays unchanged");
		}
		if(passChanges!=null && (initialSeeds==null || merge || topN!=0)){
			throw new IllegalArgumentException("passchanges requires fixed seeds, merge=f and topn=0");
		}
		if(skipFinalRebuild){validateFinalRebuildSkip();}
		if(preserveAllSeeds){validateSeedPreservation();}
		validateWeightedMerge();
		if(pruneFinal && (initialSeeds==null || passes<1 || merge)){
			throw new IllegalArgumentException("prunefinal requires fixed seeds, passes>=1 and merge=false");
		}
		started=true;
		if(initialSeeds!=null){loadSeeds();}
		else{
			if(concurrentGreedy){greedyConcurrent();}else{greedy();}
			if(scavenge){scavengeSingletons();}
			pruneSmall();
			rebuild(0);
		}
		for(int pass=1; pass<=passes; pass++){
			if(trackPassChanges()){beginPassChanges();}
			reassign(pass);
			final long afterReassign=assigned;
			final long changedReassign=trackPassChanges() ? changedAssignments() : 0;
			// Rebuilding changes only models, so final pruning here has the same
			// membership outcome as after the final rebuild and avoids wasted work.
			if(!preserveAllSeeds && (!pruneFinal || pass==passes)){pruneSmall();}
			if(!skipFinalRebuild){rebuild(pass);}
			if(merge){mergePass(pass);}
			final long changedEnd=trackPassChanges() ? changedAssignments() : 0;
			completedPasses=pass;
			if(passChanges!=null){endPassChanges(pass,afterReassign,changedReassign,changedEnd);}
			if(stopChangeFraction>0 && pass>1 && belowStopFraction(changedEnd,sequences.length,stopChangeFraction)){break;}
		}
		final ArrayList<Family> retained=liveFamilies();
		if(preserveAllSeeds && retained.size()!=initialSeeds.length){
			throw new IllegalStateException("Seed-preserving assignment/refinement must retain every imported model, including empty references");
		}
		Collections.sort(retained, SIZE_ORDER);
		if(topN>0 && retained.size()>topN){
			for(int i=topN; i<retained.size(); i++){remove(retained.get(i));}
			retained.subList(topN, retained.size()).clear();
		}
		verifyAssignments(retained);
		metrics.begin("final", completedPasses, order.length);
		metrics.end(order.length, assigned, order.length-assigned, index.size(), createdFamilies());
		completedRun=true;
		return new Result(retained.toArray(new Family[retained.size()]), assignments,
			order.length, assigned, createdFamilies(), completedPasses);
	}

	/** Model counts describe loaded seeds; no input is assigned until the first frozen pass. */
	private void loadSeeds(){
		metrics.begin("seeds",0,order.length);
		for(Family family:initialSeeds){
			while(families.size()<=family.familyId){families.add(null);}
			families.set(family.familyId,family); index.upsert(family.familyId,family.sketchUnsafe());
		}
		assert(assigned==0 && index.size()==initialSeeds.length) : "Loaded seeds are models, not assignments to the current inputs";
		metrics.end(order.length,0,order.length,index.size(),createdFamilies());
	}

	private int createdFamilies(){return initialSeeds==null ? families.size() : initialSeeds.length;}

	/** Restricts greedy consensus unions using original-tail weights, not consensus-record counts. */
	public GreedyClusterer mergeSizes(int[] originalMembers,int maximumSize){
		if(memberInputs!=null){throw new IllegalArgumentException("Weighted consensus inputs cannot substitute for original-member ranking descriptors");}
		if(started || mergeInputWeights!=null || originalMembers==null || originalMembers.length!=sequences.length || maximumSize<1){
			throw new IllegalArgumentException("Configure one positive merge size limit and one original-tail weight per input before execution");
		}
		for(int i=0; i<originalMembers.length; i++){
			if(originalMembers[i]<1 || sequences[i].length()<minimumLength){
				throw new IllegalArgumentException("Weighted consensus merging requires positive weights and no length-filtered inputs: "+i);
			}
		}
		mergeInputWeights=originalMembers.clone(); mergeMaximumSize=maximumSize;
		validateWeightedMerge();
		return this;
	}

	private void validateWeightedMerge(){
		if(mergeInputWeights!=null && (passes!=0 || minimumCluster!=1 || topN!=0 || initialSeeds!=null || merge ||
			scavenge || ramp!=null || histograms!=null || audit!=null || greedyObserver!=null || assignmentMetrics!=null)){
			throw new IllegalArgumentException("Weighted consensus merging requires passes=0, mincluster=1, topn=0 and no seeds, interpass merge, scavenging, ramp or decision observers");
		}
	}

	/** Snapshot screening saves alignments; admission must still recheck the live count atomically. */
	private boolean mergeEligible(int input,Family family){
		return mergeInputWeights==null || mergeInputWeights[input]<mergeMaximumSize || family.mergeMemberCount()<mergeMaximumSize;
	}

	private Family newGreedyFamily(int id,int input,byte[] query){
		assert(id>=0 && input>=0 && input<sequences.length) : "Greedy families retain stable nonnegative family and input ordinals";
		final Family family=mergeInputWeights==null ? new Family(id,"family_"+id,input,query) :
			Family.mergeSeed(id,"family_"+id,input,query,mergeInputWeights[input]);
		if(memberInputs!=null){family.enableMemberRanking(memberInputs,memberLayer,null);}
		return family;
	}

	private boolean admitGreedy(Family family,int input){
		if(mergeInputWeights!=null){return family.tryAddMergeMember(input,mergeInputWeights[input],mergeMaximumSize);}
		family.addMember(input); return true;
	}

	/** Size filtering precedes the reported shortlist; all actual DP work remains counted once. */
	public String[] mergeSizeStatistics(){
		if(!completedRun || mergeInputWeights==null){throw new IllegalStateException("Merge size statistics require a completed weighted merge run");}
		return new String[]{"merge_size_filtered\t"+mergeSizeFiltered,
			"merge_seed_size_filtered\t"+mergeSeedSizeFiltered,
			"merge_reservation_failures\t"+mergeReservationFailures,
			"merge_shortlist_scope\tsize_eligible_before_alignment"};
	}

	/** Every successful join is visible to the next input, including scheduled model updates. */
	private void greedy(){
		metrics.begin("greedy", 0, order.length);
		for(int processed=0; processed<order.length; processed++){
			final int input=order[processed];
			Family family=bestFamily(input,"greedy",0);
			if(family==null){
				metrics.alignment(true);
				scorer.resetConsensus(sequences[input].enc);
				final TailAlignmentDecision self=gate.decide(sequences[input].enc,sequences[input].enc,scorer);
				if(!self.accepted()){
					throw new IllegalArgumentException("Input "+input+" ("+sequences[input].id+
						") cannot seed a valid family under the configured gate: "+self);
				}
				final int id=families.size();
				family=newGreedyFamily(id,input,sequences[input].enc);
				families.add(family);
			}else if(!admitGreedy(family,input)){
				throw new IllegalStateException("Serial weighted admission changed eligibility without a concurrent writer");
			}
			if(ramp!=null){ramp.admitted(family.familyId,family.size());}
			if(greedyObserver!=null){greedyObserver.admit(input,family.familyId);}
			assignments[input]=family.familyId; assigned++;
			final boolean consensusDue=family.shouldRebuildConsensus();
			if(consensusDue){
				updater.update(family, sequences, gate, metrics, true);
				family.markConsensusRebuilt();
			}
			if(consensusDue || family.shouldRebuildSketch()){refresh(family);}
			metrics.checkpoint(processed+1, assigned, order.length-assigned, index.size(), createdFamilies());
		}
		metrics.end(order.length, assigned, order.length-assigned, index.size(), createdFamilies());
	}

	/** Detaches every initial singleton before any query so no singleton can be a target.
	 * The original input permutation defines replay order even after concurrent greedy.
	 * Targets only grow; no seed is created and stable family IDs are never reused.
	 */
	private void scavengeSingletons(){
		final long wallStart=System.nanoTime();
		final ThreadMXBean cpu=ManagementFactory.getThreadMXBean();
		final long cpuStart=cpu.isCurrentThreadCpuTimeSupported() && cpu.isThreadCpuTimeEnabled() ?
			cpu.getCurrentThreadCpuTime() : -1;
		metrics.begin("scavenge",0,order.length);
		final IntList singletons=new IntList();
		for(int input : order){
			final int id=assignments[input];
			if(id<0){throw new IllegalStateException("Greedy must assign every eligible tail before scavenging");}
			final Family family=families.get(id);
			if(family==null){throw new IllegalStateException("Scavenging found a removed greedy family: "+id);}
			if(family.size()==1){singletons.add(input); remove(family);}
		}
		scavengeInputs=singletons.size;
		final long retainedAssignments=assigned;
		final int retainedFamilies=index.size();
		long processed=order.length-singletons.size;
		try(ConsensusRefiner replay=concurrentGreedy ? ((ConsensusRefiner)updater).fork() : null){
			for(int i=0; i<singletons.size; i++){
				final int input=singletons.get(i);
				assert(assignments[input]==-1) : "All initial singleton memberships must be removed before replay";
				Family family=bestFamily(input,"scavenge",0);
				if(family!=null){
					assert(family.size()>=2) : "Scavenging may target only surviving nonsingleton families";
					family.addMember(input); assignments[input]=family.familyId; assigned++; scavengeRescued++;
					if(ramp!=null){ramp.admitted(family.familyId,family.size());}
					final boolean consensusDue=family.shouldRebuildConsensus();
					if(consensusDue){
						if(replay==null){updater.update(family,sequences,gate,metrics,true);}
						else{
							// Late concurrent admissions need inclusive old-model placement until
							// ordinary reassignment. Proposed models still admit every member.
							final Family.MembershipView view=family.membershipView();
							family=new Family(family);
							final ConsensusRefiner.Update update=replay.updateOnline(family,view,sequences,gate,metrics);
							families.set(family.familyId,family);
							((ConsensusRefiner)updater).merge(update);
							metrics.consensusUpdated(update.accepted());
						}
						family.markConsensusRebuilt();
					}
					if(consensusDue || family.shouldRebuildSketch()){refresh(family);}
				}
				metrics.checkpoint(++processed,assigned,order.length-assigned,index.size(), createdFamilies());
			}
		}
		if(index.size()!=retainedFamilies || assigned!=retainedAssignments+scavengeRescued){
			throw new IllegalStateException("Scavenging must conserve target families and rescued memberships");
		}
		metrics.end(processed,assigned,order.length-assigned,index.size(), createdFamilies());
		scavengeWallNanos=System.nanoTime()-wallStart;
		final long cpuEnd=cpuStart>=0 && cpu.isThreadCpuTimeEnabled() ? cpu.getCurrentThreadCpuTime() : -1;
		scavengeCpuAvailable=cpuStart>=0 && cpuEnd>=cpuStart;
		scavengeCpuNanos=scavengeCpuAvailable ? cpuEnd-cpuStart : 0;
	}

	/** Stage-local rescues precede size pruning and later passes; CPU is the replay driver only. */
	public String[] scavengeStatistics(){
		if(!completedRun || !scavenge){throw new IllegalStateException("Scavenge statistics require a completed enabled run");}
		return new String[]{"scavenge\ttrue","scavenge_inputs\t"+scavengeInputs,
			"scavenge_rescued\t"+scavengeRescued,"scavenge_unassigned\t"+(scavengeInputs-scavengeRescued),
			"scavenge_elapsed_seconds\t"+Double.toString(scavengeWallNanos*1e-9),
			"scavenge_driver_cpu_seconds\t"+Double.toString(scavengeCpuNanos*1e-9),
			"scavenge_cpu_available\t"+(scavengeCpuAvailable ? 1 : 0)};
	}

	/** Workers publish immutable versions; the ordinary mutable index is populated only after joining. */
	private void greedyConcurrent(){
		metrics.begin("greedy",0,order.length);
		if(order.length>0){
			final ConcurrentPass work=new ConcurrentPass(this);
			work.run();
			for(Family family : work.slots.completed()){
				families.add(family); index.upsert(family.familyId,family.sketchUnsafe());
			}
			assigned=order.length;
			verifyAssignments(liveFamilies());
		}
		metrics.end(order.length,assigned,order.length-assigned,index.size(), createdFamilies());
	}

	/** Additional concurrent-only diagnostics; old serial configuration and output stay unchanged. */
	public String[] concurrentStatistics(){
		if(!concurrentGreedy || !completedRun){throw new IllegalStateException("Concurrent statistics require a completed concurrent run");}
		return new String[]{"concurrentgreedy\ttrue","snapshotfraction\t"+snapshotFraction,
			"greedy_threads\t"+greedyThreads,"greedy_workers_started\t"+greedyWorkersStarted,
			"seed_collisions\t"+seedCollisions,"version_collisions\t"+versionCollisions,
			"discarded_version_proposals\t"+discardedVersionProposals,
			"discarded_version_passing_proposals\t"+discardedVersionPassingProposals,
			"late_member_checks\t"+lateChecks,"late_member_failures\t"+lateFailures,
			"snapshot_build_nanos\t"+snapshotBuildNanos,"snapshot_rebuilds\t"+snapshotRebuilds,
			"seed_candidates\t"+metrics.seedCandidates(),"seed_alignments\t"+metrics.seedAlignments(),
			"seed_hit_filtered\t"+metrics.seedHitFiltered(),"seed_hydro_filtered\t"+metrics.seedHydroFiltered(),
			"seed_matches\t"+metrics.seedMatches(),
			"incompatible_snapshots\t"+((ConsensusRefiner)updater).incompatibleSnapshots()};
	}

	/** One asynchronous admission phase. No worker mutates the serial family table or index. */
	private static final class ConcurrentPass implements Accumulator<GreedyWorker>{
		ConcurrentPass(GreedyClusterer owner){
			this.owner=owner; slots=new FamilySlots(owner.order.length);
			index=new SnapshotFamilyIndex(owner.sketches.width(),owner.snapshotFraction,slots);
		}
		void run(){
			final ArrayList<GreedyWorker> workers=new ArrayList<GreedyWorker>();
			for(int i=0; i<Math.min(owner.greedyThreads,owner.order.length); i++){workers.add(new GreedyWorker(this));}
			owner.greedyWorkersStarted=workers.size();
			try{
				index.start();
				ThreadWaiter.startThreads(workers);
				long lastReport=System.nanoTime();
				for(GreedyWorker worker : workers){
					while(worker.isAlive()){
						worker.join(1000);
						final long now=System.nanoTime();
						if(now-lastReport>=5000000000L){
							long done=0;
							for(GreedyWorker running : workers){done=Math.addExact(done,running.completed);}
							owner.metrics.concurrentProgress(done,slots.size()); lastReport=now;
						}
					}
				}
				for(GreedyWorker worker : workers){accumulate(worker);}
				if(!success()){fail(new IllegalStateException("Incomplete concurrent greedy pass"));}
			}catch(Throwable problem){
				fail(problem);
				final ArrayList<GreedyWorker> started=new ArrayList<GreedyWorker>();
				for(GreedyWorker worker : workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}finally{
				try{index.close();}catch(Throwable problem){fail(problem);}
			}
			if(failure!=null){throw new IllegalStateException("Concurrent greedy failed before frozen-pass handoff",failure);}
			if(completed!=owner.order.length){throw new IllegalStateException("Concurrent greedy omitted input ordinals");}
			for(GreedyWorker worker : workers){
				owner.metrics.mergeCompletedWorker(worker.local);
				if(owner.histograms!=null){owner.histograms.merge(worker.histograms);}
				for(ConsensusRefiner.Update update : worker.updates){((ConsensusRefiner)owner.updater).merge(update);}
				owner.seedCollisions=Math.addExact(owner.seedCollisions,worker.seedCollisions);
				owner.versionCollisions=Math.addExact(owner.versionCollisions,worker.versionCollisions);
				owner.discardedVersionProposals=Math.addExact(owner.discardedVersionProposals,worker.discardedVersionProposals);
				owner.discardedVersionPassingProposals=Math.addExact(owner.discardedVersionPassingProposals,worker.discardedVersionPassingProposals);
				owner.lateChecks=Math.addExact(owner.lateChecks,worker.lateChecks);
				owner.lateFailures=Math.addExact(owner.lateFailures,worker.lateFailures);
				owner.mergeSizeFiltered=Math.addExact(owner.mergeSizeFiltered,worker.sizeFiltered);
				owner.mergeSeedSizeFiltered=Math.addExact(owner.mergeSeedSizeFiltered,worker.seedSizeFiltered);
				owner.mergeReservationFailures=Math.addExact(owner.mergeReservationFailures,worker.reservationFailures);
			}
			owner.snapshotBuildNanos=index.buildNanos(); owner.snapshotRebuilds=index.rebuilds();
		}
		int claim(){final long value=next.getAndAdd(200); return (int)Math.min(value,owner.order.length);}
		synchronized void fail(Throwable problem){if(failure==null){failure=problem;}}
		@Override public void accumulate(GreedyWorker worker){
			if(!worker.success){fail(new IllegalStateException("Concurrent greedy worker did not complete normally"));}
			completed=Math.addExact(completed,worker.completed);
		}
		@Override public boolean success(){return failure==null && completed==owner.order.length;}
		@Override public ReadWriteLock rwlock(){return null;}// Driver-only accumulation after join.
		final GreedyClusterer owner;
		final FamilySlots slots;
		final SnapshotFamilyIndex index;
		final AtomicLong next=new AtomicLong();
		volatile Throwable failure;
		long completed;
	}

	/** Each worker owns every scorer, sketch builder, query buffer and diagnostic accumulator. */
	private static final class GreedyWorker extends Thread implements MutableFamilyIndex.Cosine{
		GreedyWorker(ConcurrentPass work){
			this.work=work; owner=work.owner;
			scorer=new HbmPositionModel.ConsensusScorer(owner.gate.matrix);
			memberQuery=owner.memberInputs==null ? null : new MemberRanking(owner){
				@Override Family family(final int id){return work.slots.get(id);}
				@Override public double score(final int id){return GreedyWorker.this.score(id);}
			};
		}
		@Override public void run(){
			try(ClusteringMetrics metrics=new ClusteringMetrics(null,null,false,owner.shortlist,Long.MAX_VALUE);
				ConsensusRefiner refiner=((ConsensusRefiner)owner.updater).fork()){
				local=metrics; this.refiner=refiner; local.begin("greedy",0,owner.order.length);
				sketches=owner.sketches.fork(); querySketch=new char[sketches.width()]; ids=new int[owner.shortlist];
				if(owner.mergeInputWeights!=null){
					acceptedFamilies=new Family[owner.shortlist]; acceptedDecisions=new TailAlignmentDecision[owner.shortlist];
					acceptedRanks=new int[owner.shortlist];
				}
				histograms=owner.histograms==null ? null : new AssignmentHistograms(owner.passes,sketches.width(),owner.shortlist);
				for(int start=work.claim(); start<owner.order.length && work.failure==null; start=work.claim()){
					final int stop=(int)Math.min((long)start+200,owner.order.length);
					for(int i=start; i<stop && work.failure==null; i++){
						work.index.checkFailure(); process(owner.order[i]); completed++;
					}
				}
				success=work.failure==null;
			}catch(Throwable problem){success=false; work.fail(problem);}
		}
		private void process(int input){
			if(acceptedFamilies!=null){Arrays.fill(acceptedFamilies,null); Arrays.fill(acceptedDecisions,null);}
			final boolean earlyGreedy=owner.ramp!=null && owner.ramp.early(work.slots.size());
			final byte[] query=owner.sequences[input].enc;
			queryFilter.reset(query); sketches.sketchEncoded(query,querySketch);
			Arrays.fill(monomers,0);
			for(byte aa : query){if(aa>=0 && aa<20){monomers[aa]++;}}
			norm=0; for(double value : monomers){norm+=value*value;}
			if(histograms!=null){histograms.beginQuery("greedy",0,query);}
			if(memberQuery!=null){memberQuery.prepare(query.length,queryFilter,true,earlyGreedy,null);}
			final int candidates=work.index.rank(querySketch,owner.shortlist,ids,null,memberQuery==null ? this : memberQuery,workspace);
			Family best=null; TailAlignmentDecision bestDecision=null;
			int rank=0,passing=0,hitsFiltered=0,hydroFiltered=0,eligibleCandidates=0,bestSlot=-1;
			for(int i=0; i<candidates; i++){
				final Family family=work.slots.get(ids[i]);
				if(!owner.mergeEligible(input,family)){sizeFiltered++; continue;}
				eligibleCandidates++;
				final int hits=workspace.scoreFor(family.familyId);
				final int screen=owner.filter.rejection(true,hits,occupied(family),queryFilter,family,earlyGreedy);
				if(screen==CandidateFilter.LOW_HITS){hitsFiltered++; continue;}
				if(screen==CandidateFilter.HYDRO){hydroFiltered++; continue;}
				local.alignment(false);
				final TailAlignmentDecision decision=decide(query,family);
				if(decision.accepted()){
					passing++;
					if(acceptedFamilies!=null){acceptedFamilies[i]=family; acceptedDecisions[i]=decision; acceptedRanks[i]=eligibleCandidates;}
					if(best==null || TailAlignmentDecision.compareAccepted(decision,family.familyId,bestDecision,best.familyId)<0){
						best=family; bestDecision=decision; rank=eligibleCandidates; bestSlot=i;
					}
				}else if(histograms!=null){histograms.rejected(hits,cosine(family),family.scalarsUnsafe(),family.scalarMask());}
			}
			boolean joined=false;
			if(acceptedFamilies!=null){
				// A lost reservation does not repeat search or alignment: try the
				// remaining scored candidates, then the ordinary seed competition.
				while(best!=null && !owner.admitGreedy(best,input)){
					reservationFailures++; acceptedFamilies[bestSlot]=null;
					best=null; bestDecision=null; bestSlot=-1; rank=0;
					for(int i=0; i<candidates; i++){
						final Family candidate=acceptedFamilies[i];
						if(candidate==null){continue;}
						if(best==null || TailAlignmentDecision.compareAccepted(acceptedDecisions[i],candidate.familyId,bestDecision,best.familyId)<0){
							best=candidate; bestDecision=acceptedDecisions[i]; bestSlot=i; rank=acceptedRanks[i];
						}
					}
				}
				joined=best!=null;
			}
			local.filtered(hitsFiltered,hydroFiltered);
			if(histograms!=null){
				if(best!=null){histograms.chosen(workspace.scoreFor(best.familyId),cosine(best),best.scalarsUnsafe(),best.scalarMask(),rank,best.size());}
				histograms.finishQuery(candidates,passing,best!=null,hitsFiltered,hydroFiltered);
			}
			boolean seeded=false;
			if(best==null){
				boolean selfChecked=false;
				for(int id=workspace.rankedFamilies(); best==null; id++){
					if(work.failure!=null){return;}
					Family current=work.slots.get(id);
					if(current==null){
						if(!selfChecked){
							local.alignment(true); scorer.resetConsensus(query);
							if(!owner.gate.decide(query,query,scorer).accepted()){throw new IllegalArgumentException("Input cannot seed a valid family: "+input);}
							selfChecked=true;
						}
						final Family candidate=owner.newGreedyFamily(id,input,query);
						candidate.refreshStats(owner.sequences,sketches); candidate.markSketchRebuilt(); local.sketchUpdated();
						if(work.slots.seed(id,candidate)){best=candidate; seeded=true; break;}
						seedCollisions++; current=work.slots.get(id);
					}
					if(!owner.mergeEligible(input,current)){seedSizeFiltered++; continue;}
					final int screen=owner.filter.rejection(true,matches(current.sketchUnsafe()),occupied(current),queryFilter,current,earlyGreedy);
					local.seedCandidate(screen);
					if(screen==CandidateFilter.PASS && decide(query,current).accepted()){
						if(acceptedFamilies!=null){
							if(!owner.admitGreedy(current,input)){reservationFailures++; continue;}
							joined=true;
						}
						best=current; local.seedMatch();
					}
				}
			}
			if(!seeded && !joined && !owner.admitGreedy(best,input)){
				throw new IllegalStateException("Unweighted admission cannot lose a size reservation");
			}
			local.query(eligibleCandidates,rank);
			if(owner.ramp!=null){owner.ramp.admitted(best.familyId,best.size());}
			owner.assignments[input]=best.familyId;
			final Family current=work.slots.get(best.familyId);
			if(!seeded && current!=best){checkLate(query,current);}
			rebuild(current);
		}
		private void rebuild(Family current){
			final boolean consensusDue=current.shouldRebuildConsensus();
			if(!consensusDue && !current.shouldRebuildSketch()){return;}
			final Family.MembershipView view=current.membershipView();
			final Family candidate=new Family(current);
			ConsensusRefiner.Update update=null;
			if(consensusDue){
				update=refiner.updateOnline(candidate,view,owner.sequences,owner.gate,local);
				candidate.markConsensusRebuilt(view);
			}
			candidate.refreshStats(owner.sequences,sketches,view); candidate.markSketchRebuilt(view); local.sketchUpdated();
			if(!work.slots.replace(current,candidate)){
				versionCollisions++;
				if(update!=null){
					discardedVersionProposals+=update.proposals();
					discardedVersionPassingProposals+=update.accepted();
				}
				return;
			}
			if(update!=null){updates.add(update); local.consensusUpdated(update.accepted());}
			if(update!=null && update.accepted()>0){
				final Family.MembershipView after=candidate.membershipView();
				for(int i=view.size(); i<after.size(); i++){checkLate(owner.sequences[after.memberIndex(i)].enc,candidate);}
			}
		}
		private void checkLate(byte[] query,Family family){
			lateChecks++; local.alignment(true);
			if(!decide(query,family).accepted()){lateFailures++;}
		}
		private TailAlignmentDecision decide(byte[] query,Family family){
			final byte[] reference=family.consensusUnsafe(); scorer.resetConsensus(reference);
			return owner.gate.decide(query,reference,scorer);
		}
		private int occupied(Family family){
			return owner.filter.greedyHits>=0 ? 0 : family.occupiedBins();
		}
		private int matches(char[] row){
			int hits=0; for(int i=0; i<row.length; i++){if(querySketch[i]!=0 && querySketch[i]==row[i]){hits++;}} return hits;
		}
		@Override public double score(int id){return cosine(work.slots.get(id));}
		private double cosine(Family family){
			double dot=0,targetNorm=0; final float[] target=family.monomersUnsafe();
			for(int i=0; i<20; i++){dot+=monomers[i]*target[i]; targetNorm+=target[i]*(double)target[i];}
			return norm==0 || targetNorm==0 ? 0 : Math.min(1,Math.max(0,dot/Math.sqrt(norm*targetNorm)));
		}
		final ConcurrentPass work;
		final GreedyClusterer owner;
		final MemberRanking memberQuery;
		final MutableFamilyIndex.Workspace workspace=new MutableFamilyIndex.Workspace();
		final CandidateFilter.Context queryFilter=new CandidateFilter.Context();
		final HbmPositionModel.ConsensusScorer scorer;
		final double[] monomers=new double[20];
		final ArrayList<ConsensusRefiner.Update> updates=new ArrayList<ConsensusRefiner.Update>();
		FamilySketchBuilder sketches;
		ConsensusRefiner refiner;
		ClusteringMetrics local;
		AssignmentHistograms histograms;
		char[] querySketch;
		int[] ids;
		Family[] acceptedFamilies;
		TailAlignmentDecision[] acceptedDecisions;
		int[] acceptedRanks;
		double norm;
		volatile long completed;
		long seedCollisions,versionCollisions,lateChecks,lateFailures,discardedVersionProposals,discardedVersionPassingProposals;
		long sizeFiltered,seedSizeFiltered,reservationFailures;
		boolean success;
	}

	/** Membership changes, but reference sequences and all ranking statistics stay frozen until the end. */
	private void reassign(int pass){
		metrics.begin("reassign", pass, order.length);
		if(memberInputs!=null){
			frozenMemberRanks=new MemberRankStatistics.Snapshot[families.size()];
			for(int id=0; id<families.size(); id++){
				final Family family=families.get(id);
				if(family!=null){frozenMemberRanks[id]=family.memberRankingSnapshot();}
			}
		}
		if(referenceSizes==null || referenceSizes.length!=families.size()){referenceSizes=new int[families.size()];}
		for(int id=0; id<families.size(); id++){
			final Family family=families.get(id);
			referenceSizes[id]=family==null ? 0 : family.referenceSize();
			assert(family==null || referenceSizes[id]>0 ||
				(preserveAllSeeds || pruneFinal && pass>1) && family.isEmpty()) :
				"Zero-member references require explicit seed preservation or a deferred-pruning intermediate pass";
		}
		if(audit!=null && pass==passes){
			audit.beginPass(pass,index.size(),order.length);
			for(Family family : families){if(family!=null){audit.reference(family,referenceSizes[family.familyId]);}}
		}
		if(threads>1 && order.length>0){
			reassignParallel(pass);
			return;
		}
		for(Family family : families){if(family!=null){family.clearMembers();}}
		Arrays.fill(assignments, -1); assigned=0;
		for(int processed=0; processed<order.length; processed++){
			final int input=order[processed];
			final Family family=bestFamily(input,"reassign",pass);
			if(family!=null){
				family.addMember(input); assignments[input]=family.familyId; assigned++;
			}
			metrics.checkpoint(processed+1, assigned, order.length-assigned, index.size(), createdFamilies());
		}
		metrics.end(order.length, assigned, order.length-assigned, index.size(), createdFamilies());
		if(audit!=null && pass==passes){audit.endPass();}
	}

	/** Workers read immutable references; only the driver changes memberships after all joins. */
	private void reassignParallel(int pass){
		final MutableFamilyIndex.FrozenView view=frozenIndexThreads==1 ? index.freeze() : index.freezeIndexed(frozenIndexThreads);
		final ParallelPass work=new ParallelPass(this,pass,null,view);
		work.run();
		for(Family family : families){if(family!=null){family.clearMembers();}}
		Arrays.fill(assignments,-1); assigned=0;
		for(int processed=0; processed<order.length; processed++){
			final int input=order[processed], id=work.winners[processed];
			if(work.candidates[processed]<0){throw new IllegalStateException("Worker omitted input ordinal "+input);}
			metrics.alignments(false,work.aligned[processed]);
			metrics.filtered(work.hitFiltered[processed],work.hydroFiltered[processed]);
			metrics.query(work.candidates[processed],work.ranks[processed]);
			if(id>=0){
				final Family family=families.get(id);
				if(family==null){throw new IllegalStateException("Frozen query selected a removed family "+id);}
				if(work.identities!=null){assignmentMetrics.record(input,id,work.identities[processed],work.scoreRatios[processed]);}
				family.addMember(input); assignments[input]=id; assigned++;
			}
			metrics.checkpoint(processed+1,assigned,order.length-assigned,index.size(), createdFamilies());
		}
		metrics.end(order.length,assigned,order.length-assigned,index.size(), createdFamilies());
	}

	/** Preserves ranking, screens before alignment, then selects the best admitted candidate. */
	private Family bestFamily(int input,String stage,int pass){
		final boolean earlyGreedy=pass==0 && ramp!=null && ramp.early(families.size());
		final boolean captureGreedy=pass==0 && greedyObserver!=null;
		final byte[] query=sequences[input].enc;
		queryFilter.reset(query);
		if(histograms!=null){histograms.beginQuery(stage,pass,query);}
		sketches.sketchEncoded(query, querySketch);
		Arrays.fill(queryMonomers, 0);
		for(byte aa : query){if(aa>=0 && aa<20){queryMonomers[aa]++;}}
		queryNorm=0;
		for(double x : queryMonomers){queryNorm+=x*x;}
		if(memberQuery!=null){memberQuery.prepare(query.length,queryFilter,pass==0,earlyGreedy,pass==0 ? null : frozenMemberRanks);}
		final MutableFamilyIndex.Cosine ranking=memberQuery==null ? cosine : memberQuery;
		final int candidates=index.indexed() ? index.rankIndexed(querySketch, shortlist, candidateIds, null, ranking, workspace) :
			index.rankFullScan(querySketch, shortlist, candidateIds, null, ranking, workspace);
		if(captureGreedy){
			greedyObserver.beginQuery(input,sequences[input].id,candidates,queryFilter);
			final int target=greedyObserver.historicalFamily();
			if(target>=0){
				final int hits=workspace.scoreFor(target), occupied=index.occupiedBins(target);
				final int required=filter.requiredHits(true,occupied,earlyGreedy);
				// Observe the current ranking before any family admission or model update.
				// A capped rank of 501 means greater than 500, never an exact rank.
				int observedRank=0;
				if(hits>=required){
					// A target already in the shortlist has an exact rank for free.
					for(int i=0; i<candidates; i++){
						if(candidateIds[i]==target){observedRank=i+1; break;}
					}
					if(observedRank==0){observedRank=index.rankOf(target,ranking,workspace,500);}
				}
				greedyObserver.historicalState(hits,required,occupied,observedRank,index.size());
			}
		}
		final boolean capture=audit!=null && pass==passes;
		if(capture){
			final int target=audit.targetFor(input);
			final Family family=target>=0 && target<families.size() ? families.get(target) : null;
			audit.beginQuery(input,sequences[input].id,candidates,family,
				family==null ? 0 : workspace.scoreFor(target),queryFilter);
		}
		Family best=null; TailAlignmentDecision bestDecision=null; int rank=0,passing=0,hitFiltered=0,hydroFiltered=0,eligibleCandidates=0;
		for(int i=0; i<candidates; i++){
			final Family family=families.get(candidateIds[i]);
			assert(family!=null) : "The live sketch index must never refer to a pruned family";
			if(!mergeEligible(input,family)){mergeSizeFiltered++; continue;}
			eligibleCandidates++;
			final int hits=workspace.scoreFor(family.familyId);
			final int occupied=pass==0 ? (filter.greedyHits<0 || captureGreedy ? index.occupiedBins(family.familyId) : 0) :
				(filter.passHitsFlex ? family.occupiedBins() : 0);
			final int screened=filter.rejection(pass==0,hits,occupied,queryFilter,family,earlyGreedy);
			if(capture && screened!=CandidateFilter.PASS){
				audit.candidate(i+1,family,referenceSizes[family.familyId],hits,queryFilter,screened,null);
			}
			if(captureGreedy && screened!=CandidateFilter.PASS){
				greedyObserver.candidate(i+1,family,hits,filter.requiredHits(true,occupied,earlyGreedy),
					occupied,queryFilter,screened,null);
			}
			if(screened==CandidateFilter.LOW_HITS){hitFiltered++; continue;}
			if(screened==CandidateFilter.HYDRO){hydroFiltered++; continue;}
			metrics.alignment(false);
			final byte[] reference=family.consensusUnsafe();
			scorer.resetConsensus(reference);
			final TailAlignmentDecision decision=gate.decide(query,reference,scorer);
			if(capture){audit.candidate(i+1,family,referenceSizes[family.familyId],hits,queryFilter,screened,decision);}
			if(captureGreedy){greedyObserver.candidate(i+1,family,hits,filter.requiredHits(true,occupied,earlyGreedy),
				occupied,queryFilter,screened,decision);}
			if(histograms!=null){
				if(decision.accepted()){passing++;}
				else{histograms.rejected(workspace.scoreFor(family.familyId),cosine.score(family.familyId),
					family.scalarsUnsafe(),family.scalarMask());}
			}
			if(decision.accepted() && (best==null ||
				TailAlignmentDecision.compareAccepted(decision, family.familyId, bestDecision, best.familyId)<0)){
				best=family; bestDecision=decision; rank=eligibleCandidates;
			}
		}
		if(histograms!=null){
			if(best!=null){histograms.chosen(workspace.scoreFor(best.familyId),cosine.score(best.familyId),
				best.scalarsUnsafe(),best.scalarMask(),rank,pass==0 ? best.size() : referenceSizes[best.familyId]);}
			histograms.finishQuery(candidates,passing,best!=null,hitFiltered,hydroFiltered);
		}
		metrics.filtered(hitFiltered,hydroFiltered);
		metrics.query(eligibleCandidates, rank);
		if(capture){audit.finishQuery(best==null ? -1 : best.familyId,rank);}
		if(captureGreedy){greedyObserver.finishQuery(best==null ? -1 : best.familyId,rank);}
		if(assignmentMetrics!=null && pass==passes && best!=null){
			assignmentMetrics.record(input,best.familyId,bestDecision.identity,bestDecision.scoreRatio);
		}
		return best;
	}

	/** Rebuilds every retained family's consensus and ranking statistics between frozen passes. */
	private void rebuild(int pass){
		metrics.begin("consensus", pass, order.length);
		final boolean capture=captureQuality(pass);
		long processed=order.length-assigned;
		// Custom updater implementations have no worker-factory contract; retain their serial path.
		if((threads>1 || concurrentGreedy && pass==0) && updater instanceof ConsensusRefiner){
			final ArrayList<Family> rebuildFamilies=liveFamilies();
			if(preserveAllSeeds || pruneFinal && pass<passes){
				// Empty references stay queryable with their last published model;
				// they cannot contribute a new consensus or modal sketch yet.
				int kept=0;
				for(int i=0; i<rebuildFamilies.size(); i++){
					final Family family=rebuildFamilies.get(i);
					if(!family.isEmpty()){rebuildFamilies.set(kept++,family);}
				}
				rebuildFamilies.subList(kept,rebuildFamilies.size()).clear();
			}
			final Family[] retained=rebuildFamilies.toArray(new Family[rebuildFamilies.size()]);
			final ParallelPass work=new ParallelPass(this,pass,retained,null);
			work.run();
			for(int i=0; i<retained.length; i++){
				final Family family=retained[i];
				final ConsensusRefiner.Update update=work.updates[i];
				if(update==null){throw new IllegalStateException("Worker omitted family "+family.id);}
				((ConsensusRefiner)updater).merge(update);
				if(capture){publishQuality(update.quality());}
				metrics.alignments(true,update.alignments());
				metrics.consensusUpdated(update.accepted());
				family.markSketchRebuilt();
				index.upsert(family.familyId,family.sketchUnsafe());
				metrics.sketchUpdated();
				processed+=family.size();
				metrics.checkpoint(processed,assigned,order.length-assigned,index.size(), createdFamilies());
			}
		}else{
			for(Family family : families){
				if(family==null){continue;}
				if((preserveAllSeeds || pruneFinal && pass<passes) && family.isEmpty()){continue;}
				if(capture){publishQuality(((ConsensusRefiner)updater).updateObserved(family,sequences,gate,metrics,false));}
				else{updater.update(family, sequences, gate, metrics, false);}
				refresh(family);
				processed+=family.size();
				metrics.checkpoint(processed, assigned, order.length-assigned, index.size(), createdFamilies());
			}
		}
		// Pruning leaves stale postings; compact only when another frozen query pass follows.
		if(index.indexed() && pass<passes && !merge){index.compact();}
		metrics.end(processed, assigned, order.length-assigned, index.size(), createdFamilies());
	}

	/** Builds replacements privately, then changes ownership only after every worker succeeds.
	 * A rejected union remains provisional until the next all-input reassignment;
	 * it must not enter the ordinary strict rebuild again in this pass.
	 */
	private void mergePass(final int pass){
		metrics.begin("merge",pass,order.length);
		final MergeSummary summary=new MergeSummary(pass,index.size(),assigned);
		final Family[] before=families.toArray(new Family[families.size()]);
		final DdlPassUnionBuilder.ModelInput input=DdlPassUnionBuilder.collectModels(before);
		DdlPassUnionBuilder.Group[] groups=new DdlPassUnionBuilder.Group[0];
		if(input.size()>1){
			try(ConsensusRefiner refiner=((ConsensusRefiner)updater).freshDriver();
				ClusteringMetrics local=new ClusteringMetrics(null,null,false,shortlist,Long.MAX_VALUE)){
				final GreedyClusterer nested=new GreedyClusterer(input.sequences,gate.minSpan,shortlist,1,
					0,0,seed,sketches.fork(),mergeGate,refiner,local,index.indexed(),threads,null,filter);
				if(concurrentGreedy){nested.concurrentGreedy(snapshotFraction).greedyThreads(greedyThreads);}
				groups=DdlPassUnionBuilder.mergeGroups(input,nested.run());
				summary.modelAlignments=local.totalAlignments();
			}
		}
		final MergeWork work=new MergeWork(this,pass,before,groups);
		work.run();
		final boolean[] touched=new boolean[before.length];
		for(DdlPassUnionBuilder.Replacement replacement:work.replacements){
			if(replacement==null){throw new IllegalStateException("Union worker omitted a replacement");}
			long sourceMembers=0;
			for(int id:replacement.sourceFamilyIds){
				if(id<0 || id>=before.length || before[id]==null || touched[id] || families.get(id)!=before[id]){
					throw new IllegalStateException("Merge groups must own disjoint live source families: "+id);
				}
				touched[id]=true; sourceMembers=Math.addExact(sourceMembers,before[id].size());
			}
			final Family union=replacement.family;
			if(union.familyId!=replacement.sourceFamilyIds[0] || union.size()!=sourceMembers){
				throw new IllegalStateException("Union must preserve source membership and its smallest family ID");
			}
			for(int i=0; i<union.size(); i++){
				final int member=union.memberIndex(i);
				if(Arrays.binarySearch(replacement.sourceFamilyIds,assignments[member])<0){
					throw new IllegalStateException("Union member belongs to a different source family: "+member);
				}
			}
			if(consensusQuality!=null && pass==passes && replacement.quality==null){
				throw new IllegalStateException("Final union omitted its requested cohesion snapshot");
			}
		}
		for(DdlPassUnionBuilder.Replacement replacement:work.replacements){
			final Family union=replacement.family;
			for(int id:replacement.sourceFamilyIds){index.remove(id); families.set(id,null);}
			families.set(union.familyId,union); index.upsert(union.familyId,union.sketchUnsafe());
			for(int i=0; i<union.size(); i++){assignments[union.memberIndex(i)]=union.familyId;}
			((ConsensusRefiner)updater).merge(replacement.update);
			if(consensusQuality!=null && pass==passes){publishQuality(replacement.quality);}
			metrics.sketchUpdated(); metrics.consensusUpdated(replacement.update.accepted());
			summary.add(replacement);
		}
		if(assigned!=summary.assignedBefore || index.size()!=summary.familiesBefore-summary.sourceFamilies+groups.length){
			throw new IllegalStateException("Pass merging must conserve assignments and remove only grouped source families");
		}
		verifyAssignments(liveFamilies());
		summary.familiesAfter=index.size(); summary.assignedAfter=assigned;
		metrics.alignments(true,Math.addExact(summary.modelAlignments,summary.unionAlignments));
		if(index.indexed() && pass<passes){index.compact();}
		metrics.end(order.length,assigned,order.length-assigned,index.size(),createdFamilies());
		mergeSummaries.add(summary);
	}

	/** One record per enabled pass; model/union work is distinct from original-query alignment. */
	public String[] mergeStatistics(){
		if(!completedRun || !merge || mergeSummaries.size()!=completedPasses){
			throw new IllegalStateException("Merge statistics require every completed enabled pass to complete");
		}
		final ArrayList<String> rows=new ArrayList<String>();
		rows.add("merge\ttrue"); rows.add("merge_passes\t"+passes);
		for(MergeSummary summary:mergeSummaries){summary.append(rows);}
		return rows.toArray(new String[rows.size()]);
	}

	private void checkMergeGate(final TailAlignmentGate nestedGate){
		if(nestedGate==null){throw new IllegalArgumentException("Interpass merging needs a consensus-family gate");}
		if(nestedGate.mode!=gate.mode || nestedGate.denominator!=gate.denominator ||
				nestedGate.maxLength!=gate.maxLength || nestedGate.minSpan!=gate.minSpan || nestedGate.minPaired!=gate.minPaired || !nestedGate.matrix.same(gate.matrix)){
			throw new IllegalArgumentException("Interpass merging may override thresholds but not gate geometry or substitution matrix");
		}
	}

	private static final class MergeSummary{
		MergeSummary(int pass,int before,long assigned){this.pass=pass; familiesBefore=before; assignedBefore=assigned;}
		void add(DdlPassUnionBuilder.Replacement replacement){
			final int count=replacement.sourceFamilyIds.length;
			if(count==2){groups2++;}else if(count==3){groups3++;}else{groups4plus++;}
			sourceFamilies=Math.addExact(sourceFamilies,count);
			mergedMembers=Math.addExact(mergedMembers,replacement.family.size());
			accepted=Math.addExact(accepted,replacement.update.accepted());
			unchanged=Math.addExact(unchanged,replacement.update.unchanged());
			rejected=Math.addExact(rejected,replacement.update.rejected());
			incompatible=Math.addExact(incompatible,replacement.update.incompatibleSnapshots());
			unionAlignments=Math.addExact(unionAlignments,replacement.update.alignments());
		}
		void append(ArrayList<String> rows){
			final String prefix="merge_pass_"+pass+"_";
			final String[] names={"families_before","families_after","groups_2","groups_3","groups_4plus",
				"source_families","merged_members","assigned_before","assigned_after","accepted","unchanged",
				"rejected","incompatible_snapshots","model_alignments","union_alignments"};
			final long[] values={familiesBefore,familiesAfter,groups2,groups3,groups4plus,sourceFamilies,mergedMembers,
				assignedBefore,assignedAfter,accepted,unchanged,rejected,incompatible,modelAlignments,unionAlignments};
			assert(names.length==values.length) : "Per-pass merge report columns must retain their matching counts";
			for(int i=0; i<names.length; i++){rows.add(prefix+names[i]+"\t"+values[i]);}
		}
		final int pass,familiesBefore;
		final long assignedBefore;
		int familiesAfter;
		long assignedAfter,groups2,groups3,groups4plus,sourceFamilies,mergedMembers;
		long accepted,unchanged,rejected,incompatible,modelAlignments,unionAlignments;
	}

	/** Independent unions read an immutable family table and publish into disjoint slots. */
	private static final class MergeWork implements Accumulator<MergeWorker>{
		MergeWork(GreedyClusterer owner,int pass,Family[] before,DdlPassUnionBuilder.Group[] groups){
			this.owner=owner; this.pass=pass; this.before=before; this.groups=groups;
			replacements=new DdlPassUnionBuilder.Replacement[groups.length];
		}
		void run(){
			final ArrayList<MergeWorker> workers=new ArrayList<MergeWorker>();
			for(int i=0; i<Math.min(owner.threads,groups.length); i++){workers.add(new MergeWorker(this));}
			try{
				if(!ThreadWaiter.startAndWait(workers,this)){throw new IllegalStateException("Incomplete union worker pass");}
			}catch(Throwable cause){
				fail(cause);
				final ArrayList<MergeWorker> started=new ArrayList<MergeWorker>();
				for(MergeWorker worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null){throw new IllegalStateException("Pass unions failed before engine publication",failure);}
			if(completed!=groups.length){throw new IllegalStateException("Union workers omitted groups");}
		}
		int claim(){return (int)Math.min(next.getAndIncrement(),groups.length);}
		synchronized void fail(Throwable cause){if(failure==null){failure=cause;}}
		@Override public void accumulate(MergeWorker worker){
			if(!worker.success){fail(new IllegalStateException("Union worker did not complete normally"));}
			completed=Math.addExact(completed,worker.completed);
		}
		@Override public boolean success(){return failure==null && completed==groups.length;}
		@Override public ReadWriteLock rwlock(){return null;}// Driver only, after joining all workers.
		final GreedyClusterer owner;
		final int pass;
		final Family[] before;
		final DdlPassUnionBuilder.Group[] groups;
		final DdlPassUnionBuilder.Replacement[] replacements;
		final AtomicLong next=new AtomicLong();
		volatile Throwable failure;
		long completed;
	}

	private static final class MergeWorker extends Thread{
		MergeWorker(MergeWork work){this.work=work;}
		@Override public void run(){
			final GreedyClusterer owner=work.owner;
			try(ConsensusRefiner refiner=((ConsensusRefiner)owner.updater).fork(1);
				ClusteringMetrics local=new ClusteringMetrics(null,null,false,owner.shortlist,Long.MAX_VALUE)){
				local.begin("consensus",work.pass,owner.order.length);
				final DdlPassUnionBuilder builder=new DdlPassUnionBuilder(owner.sequences,owner.sketches.fork(),owner.gate,refiner,local);
				for(int i=work.claim(); i<work.groups.length && work.failure==null; i=work.claim()){
					work.replacements[i]=builder.build(work.groups[i],work.before,owner.consensusQuality!=null && work.pass==owner.passes);
					completed++;
				}
				success=work.failure==null;
			}catch(Throwable cause){work.fail(cause);}
		}
		final MergeWork work;
		long completed;
		boolean success;
	}

	/** Capture requested for a live family must never silently omit its snapshot. */
	private void publishQuality(final ConsensusQualityReport.Snapshot snapshot){
		if(snapshot==null){throw new IllegalStateException("Final consensus rebuild omitted a requested quality snapshot");}
		consensusQuality.publish(snapshot);
	}

	/** run() may stop after any completed pass 2+; retain its quality while those
	 * member alignments already exist. Later passes replace the same family snapshots.
	 * Fixed-length runs retain their previous final-pass-only behavior.
	 */
	private boolean captureQuality(final int pass){
		assert(pass>=0 && pass<=passes) : "Quality capture must describe a rebuild in this run's configured pass range";
		return consensusQuality!=null && (pass==passes || stopChangeFraction>0 && pass>1);
	}

	private void refresh(Family family){
		assert(!family.isEmpty()) : "Empty families must be pruned before profile construction";
		family.refreshStats(sequences, sketches);
		family.markSketchRebuilt();
		index.upsert(family.familyId, family.sketchUnsafe());
		metrics.sketchUpdated();
	}

	private void pruneSmall(){
		for(Family family : families){if(family!=null && family.size()<minimumCluster){remove(family);}}
	}

	/** Removes membership during singleton detachment, size pruning or final top-N selection. */
	private void remove(Family family){
		assert(families.get(family.familyId)==family) : "Stable IDs address the original family slot";
		for(int i=0; i<family.size(); i++){
			final int input=family.memberIndex(i);
			if(assignments[input]!=family.familyId){throw new IllegalStateException("Family removal found inconsistent membership at "+input);}
			assignments[input]=-1; assigned--;
		}
		index.remove(family.familyId); families.set(family.familyId, null);
	}

	private ArrayList<Family> liveFamilies(){
		final ArrayList<Family> live=new ArrayList<Family>(index.size());
		for(Family family : families){if(family!=null){live.add(family);}}
		assert(live.size()==index.size()) : "Membership table and sketch index must have identical live families";
		return live;
	}

	/** Checks unique ownership without allocating a second per-input census array. */
	private void verifyAssignments(ArrayList<Family> retained){
		long members=0, labels=0;
		for(Family family : retained){
			int previous=-1;
			long originalMembers=0;
			// Sorting a detached compact list makes duplicate member IDs visible without mutating order.
			final int[] ids=family.memberIndexes(); Arrays.sort(ids);
			for(int input : ids){
				if(input<=previous || assignments[input]!=family.familyId){
					throw new IllegalStateException("Duplicate or inconsistent family member "+input+" in "+family.id);
				}
				previous=input; members++;
				if(mergeInputWeights!=null){originalMembers=Math.addExact(originalMembers,mergeInputWeights[input]);}
			}
			if(mergeInputWeights!=null && originalMembers!=family.mergeMemberCount()){
				throw new IllegalStateException("Weighted merge family does not conserve original-tail counts: "+family.id);
			}
		}
		for(int input=0; input<assignments.length; input++){
			final int id=assignments[input];
			if(id>=0){
				if(id>=families.size() || families.get(id)==null || sequences[input].length()<minimumLength){
					throw new IllegalStateException("Invalid final assignment at input "+input);
				}
				labels++;
			}
		}
		if(labels!=members || labels!=assigned){throw new IllegalStateException("Final membership does not conserve assigned input ordinals");}
	}

	/** A model update never owns membership deletion; pruning and reassignment are explicit engine steps. */
	public interface ConsensusUpdater {
		void update(Family family, ProteinSequence[] sequences, TailAlignmentGate gate,
			ClusteringMetrics metrics, boolean online);
	}

	/** One phase owns its workers and disjoint result slots. No worker waits on another worker. */
	private static final class ParallelPass implements Accumulator<ProcessThread>{
		ParallelPass(GreedyClusterer owner,int pass,Family[] retained,MutableFamilyIndex.FrozenView view){
			this.owner=owner; this.pass=pass; this.retained=retained; this.view=view;
			limit=retained==null ? owner.order.length : retained.length;
			final boolean capture=retained==null && pass==owner.passes && owner.assignmentMetrics!=null;
			identities=capture ? new double[limit] : null;
			scoreRatios=capture ? new double[limit] : null;
			if(retained==null){
				winners=new int[limit]; candidates=new int[limit]; ranks=new int[limit]; updates=null;
				aligned=new int[limit]; hitFiltered=new int[limit]; hydroFiltered=new int[limit];
				Arrays.fill(winners,-1); Arrays.fill(candidates,-1);
			}else{
				winners=candidates=ranks=aligned=hitFiltered=hydroFiltered=null; updates=new ConsensusRefiner.Update[limit];
			}
		}
		void run(){
			final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>();
			for(int i=0; i<Math.min(owner.threads,limit); i++){workers.add(new ProcessThread(this));}
			try{
				if(!ThreadWaiter.startAndWait(workers,this)){throw new IllegalStateException("Incomplete clustering worker pass");}
			}catch(Throwable cause){
				fail(cause);
				// Thread.start may fail after starting only a prefix; never join unstarted NEW threads.
				final ArrayList<ProcessThread> started=new ArrayList<ProcessThread>();
				for(ProcessThread worker : workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null){throw new IllegalStateException("Clustering "+(retained==null ? "reassign" : "consensus")+
				" pass "+pass+" failed before ordered result publication",failure);}
			if(completed!=limit){throw new IllegalStateException("Workers completed "+completed+" of "+limit+" items");}
			// Merge only after all workers succeeded; partially collected passes never publish diagnostics.
			if(owner.histograms!=null && retained==null){
				for(ProcessThread worker : workers){owner.histograms.merge(worker.histograms);}
			}
		}
		synchronized void fail(Throwable cause){if(failure==null){failure=cause;}}
		int claim(int count){
			final long start=next.getAndAdd(count);
			return start>=limit ? limit : (int)start;
		}
		@Override public void accumulate(ProcessThread worker){
			if(!worker.success){fail(new IllegalStateException("Worker did not complete normally"));}
			completed=Math.addExact(completed,worker.completed);
		}
		@Override public boolean success(){return failure==null && completed==limit;}
		@Override public ReadWriteLock rwlock(){return null;}// Joining establishes visibility; accumulation is driver-only.
		final GreedyClusterer owner;
		final int pass,limit;
		final Family[] retained;
		final MutableFamilyIndex.FrozenView view;
		final int[] winners,candidates,ranks,aligned,hitFiltered,hydroFiltered;
		final double[] identities,scoreRatios;
		final ConsensusRefiner.Update[] updates;
		final AtomicLong next=new AtomicLong();
		volatile Throwable failure;
		long completed;
	}

	/** Template-shaped worker: claim a coarse batch, use private scratch, report success after completion. */
	private static final class ProcessThread extends Thread{
		ProcessThread(ParallelPass work){this.work=work;}
		@Override public void run(){
			try{
				if(work.retained==null){queries();}else{families();}
				success=work.failure==null;
			}catch(Throwable cause){work.fail(cause);}
		}
		private void queries(){
			final QueryContext query=new QueryContext(work.owner,work.view,work.pass);
			histograms=query.histograms;
			for(int start=work.claim(200); start<work.limit && work.failure==null; start=work.claim(200)){
				final int stop=(int)Math.min((long)start+200,work.limit);
				for(int i=start; i<stop; i++){
					query.classify(work.owner.order[i]);
					work.winners[i]=query.winner; work.candidates[i]=query.candidates; work.ranks[i]=query.rank;
					work.aligned[i]=query.aligned; work.hitFiltered[i]=query.hitFiltered; work.hydroFiltered[i]=query.hydroFiltered;
					if(work.identities!=null && query.winner>=0){
						work.identities[i]=query.identity; work.scoreRatios[i]=query.scoreRatio;
					}
					completed++;
				}
			}
		}
		private void families(){
			final GreedyClusterer owner=work.owner;
			final boolean capture=owner.captureQuality(work.pass);
			final FamilySketchBuilder builder=owner.sketches.fork();
			try(ConsensusRefiner refiner=((ConsensusRefiner)owner.updater).fork();
				ClusteringMetrics local=new ClusteringMetrics(null,null,false,owner.shortlist,Long.MAX_VALUE)){
				for(int i=work.claim(1); i<work.limit && work.failure==null; i=work.claim(1)){
					final Family family=work.retained[i];
					local.begin("consensus",work.pass,family.size());
					work.updates[i]=owner.concurrentGreedy && work.pass==0 ?
						refiner.updateAfterConcurrentGreedy(family,owner.sequences,owner.gate,local,capture) :
						refiner.updateFamily(family,owner.sequences,owner.gate,local,capture);
					family.refreshStats(owner.sequences,builder);
					local.end(family.size(),family.size(),0,1,1);
					completed++;
				}
			}
		}
		final ParallelPass work;
		AssignmentHistograms histograms;
		long completed;
		boolean success;
	}

	/** Private query buffers and scorer; family statistics and index contents stay frozen. */
	private static final class QueryContext implements MutableFamilyIndex.Cosine{
		QueryContext(GreedyClusterer owner,MutableFamilyIndex.FrozenView view,int pass){
			this.owner=owner; this.view=view; sketches=owner.sketches.fork();
			scorer=new HbmPositionModel.ConsensusScorer(owner.gate.matrix);
			this.pass=pass;
			histograms=owner.histograms==null ? null : new AssignmentHistograms(owner.passes,sketches.width(),owner.shortlist);
			querySketch=new char[sketches.width()]; candidateIds=new int[owner.shortlist];
			memberQuery=owner.memberInputs==null ? null : new MemberRanking(owner){
				@Override Family family(final int id){return owner.families.get(id);}
				@Override public double score(final int id){return QueryContext.this.score(id);}
			};
		}
		void classify(int input){
			final byte[] query=owner.sequences[input].enc;
			queryFilter.reset(query);
			if(histograms!=null){histograms.beginQuery("reassign",pass,query);}
			sketches.sketchEncoded(query,querySketch);
			Arrays.fill(monomers,0);
			for(byte aa : query){if(aa>=0 && aa<20){monomers[aa]++;}}
			norm=0; for(double x : monomers){norm+=x*x;}
			if(memberQuery!=null){memberQuery.prepare(query.length,queryFilter,false,false,owner.frozenMemberRanks);}
			final MutableFamilyIndex.Cosine ranking=memberQuery==null ? this : memberQuery;
			candidates=view.indexed() ? view.rankIndexed(querySketch,owner.shortlist,candidateIds,null,ranking,workspace) :
				view.rankFullScan(querySketch,owner.shortlist,candidateIds,null,ranking,workspace);
			winner=-1; rank=aligned=hitFiltered=hydroFiltered=0; TailAlignmentDecision best=null; int passing=0;
			for(int i=0; i<candidates; i++){
				final Family family=owner.families.get(candidateIds[i]);
				if(family==null){throw new IllegalStateException("Frozen index selected a missing family");}
				final int screened=owner.filter.rejection(false,workspace.scoreFor(family.familyId),
					owner.filter.passHitsFlex ? family.occupiedBins() : 0,queryFilter,family);
				if(screened==CandidateFilter.LOW_HITS){hitFiltered++; continue;}
				if(screened==CandidateFilter.HYDRO){hydroFiltered++; continue;}
				final byte[] reference=family.consensusUnsafe();
				aligned++;
				scorer.resetConsensus(reference);
				final TailAlignmentDecision decision=owner.gate.decide(query,reference,scorer);
				if(histograms!=null){
					if(decision.accepted()){passing++;}
					else{histograms.rejected(workspace.scoreFor(family.familyId),score(family.familyId),
						family.scalarsUnsafe(),family.scalarMask());}
				}
				if(decision.accepted() && (best==null || TailAlignmentDecision.compareAccepted(decision,family.familyId,best,winner)<0)){
					winner=family.familyId; best=decision; rank=i+1;
				}
			}
			if(histograms!=null){
				if(winner>=0){
					final Family family=owner.families.get(winner);
					histograms.chosen(workspace.scoreFor(winner),score(winner),family.scalarsUnsafe(),family.scalarMask(),rank,owner.referenceSizes[winner]);
				}
				histograms.finishQuery(candidates,passing,winner>=0,hitFiltered,hydroFiltered);
			}
			if(owner.assignmentMetrics!=null && pass==owner.passes && best!=null){
				identity=best.identity; scoreRatio=best.scoreRatio;
			}
		}
		@Override public double score(int id){
			final float[] target=owner.families.get(id).monomersUnsafe();
			double dot=0,targetNorm=0;
			for(int i=0; i<20; i++){dot+=monomers[i]*target[i]; targetNorm+=target[i]*(double)target[i];}
			return norm==0 || targetNorm==0 ? 0 : Math.min(1,Math.max(0,dot/Math.sqrt(norm*targetNorm)));
		}
		final GreedyClusterer owner;
		final MutableFamilyIndex.FrozenView view;
		final MemberRanking memberQuery;
		final int pass;
		final AssignmentHistograms histograms;
		final FamilySketchBuilder sketches;
		final char[] querySketch;
		final int[] candidateIds;
		final double[] monomers=new double[20];
		final MutableFamilyIndex.Workspace workspace=new MutableFamilyIndex.Workspace();
		final HbmPositionModel.ConsensusScorer scorer;
		final CandidateFilter.Context queryFilter=new CandidateFilter.Context();
		double norm,identity,scoreRatio;
		int winner,candidates,rank,aligned,hitFiltered,hydroFiltered;
	}

	/** One query worker owns this callback and its scratch. Reassignment uses only
	 * pre-clear snapshots; greedy reads a coherent current member prefix. Raw-hit
	 * admission precedes pivot evaluation/capping through the native index contract.
	 */
	private abstract static class MemberRanking implements MutableFamilyIndex.Ranking{
		MemberRanking(final GreedyClusterer owner_){owner=owner_;}
		void prepare(final int length, final CandidateFilter.Context query, final boolean greedy_,
				final boolean early_, final MemberRankStatistics.Snapshot[] frozen_){
			if(length<1 || query==null || !query.defined() || (!greedy_ && frozen_==null)){
				throw new IllegalArgumentException("Member ranking requires measured query KD and pre-clear descriptors for frozen passes");
			}
			queryLength=length; queryKd=query.hydrophobicity(); greedy=greedy_; early=early_; frozen=frozen_;
		}
		abstract Family family(int id);
		@Override public boolean accept(final int id, final int hits){
			return hits>=owner.filter.requiredHits(greedy,family(id).occupiedBins(),early);
		}
		@Override public double score(final int id, final int hits){
			final double length,kd;
			if(frozen==null){
				family(id).readMemberRanking(scratch); length=scratch.medianLength; kd=scratch.meanKd;
			}else{
				final MemberRankStatistics.Snapshot value=frozen[id];
				if(value==null){throw new IllegalStateException("Frozen index candidate has no retained-member descriptor: "+id);}
				length=value.medianLength; kd=value.meanKd;
			}
			final double ratio=Math.min(queryLength,length)/Math.max(queryLength,length);
			final double cosine=score(id), similarity=1-Math.abs(queryKd-kd)/9;
			assert(hits>=0 && ratio>0 && ratio<=1 && cosine>=0 && cosine<=1 && similarity>=0 && similarity<=1) :
				"Member pivot bounds require native cosine/KD ranges and positive measured lengths";
			return hits*Math.sqrt(Math.sqrt(ratio))+1.5*cosine+owner.memberInputs.hydroWeight*similarity;
		}
		@Override public double lowerBound(final int hits){return 0;}
		@Override public double upperBound(final int hits){return hits+1.5+owner.memberInputs.hydroWeight;}
		final GreedyClusterer owner;
		final MemberRankStatistics.View scratch=new MemberRankStatistics.View();
		MemberRankStatistics.Snapshot[] frozen;
		int queryLength;
		double queryKd;
		boolean greedy,early;
	}

	public static final class Result {
		private Result(Family[] families, int[] assignments, int eligible, long assigned, long created,
				int completedPasses){
			this.families=families; this.assignments=assignments; this.eligible=eligible;
			this.assigned=assigned; this.created=created; this.completedPasses=completedPasses;
		}
		public int familyFor(int input){return assignments[input];}
		public int inputs(){return assignments.length;}
		public final Family[] families;
		public final int eligible;
		public final int completedPasses;
		public final long assigned, created;
		private final int[] assignments;
	}

	private final MutableFamilyIndex.Cosine cosine=new MutableFamilyIndex.Cosine(){
		@Override
		public double score(int id){
			final float[] target=families.get(id).monomersUnsafe();
			double dot=0, norm=0;
			for(int i=0; i<20; i++){dot+=queryMonomers[i]*target[i]; norm+=target[i]*(double)target[i];}
			return queryNorm==0 || norm==0 ? 0 : Math.min(1, Math.max(0, dot/Math.sqrt(queryNorm*norm)));
		}
	};
	private static final Comparator<ClusterSketchCodec.Entry> SEED_SIZE_ORDER=new Comparator<ClusterSketchCodec.Entry>(){
		@Override public int compare(ClusterSketchCodec.Entry a,ClusterSketchCodec.Entry b){
			final int size=Integer.compare(b.members,a.members);
			return size!=0 ? size : Integer.compare(a.familyId,b.familyId);
		}
	};
	private static final Comparator<Family> SIZE_ORDER=new Comparator<Family>(){
		@Override
		public int compare(Family a, Family b){
			final int size=Integer.compare(b.size(), a.size());
			return size!=0 ? size : Integer.compare(a.familyId, b.familyId);
		}
	};
	private final ProteinSequence[] sequences;
	private final int minimumLength, shortlist, minimumCluster, passes, topN,threads,frozenIndexThreads;
	private int greedyThreads,greedyWorkersStarted;
	private boolean pruneFinal,skipFinalRebuild,preserveAllSeeds;
	private final long seed;
	private final FamilySketchBuilder sketches;
	private final TailAlignmentGate gate;
	private final ConsensusUpdater updater;
	private final ClusteringMetrics metrics;
	private final AssignmentHistograms histograms;
	private final CandidateFilter filter;
	private GreedyHitRamp ramp;
	private final AssignmentDecisionAudit audit;
	private GreedyDecisionAudit greedyObserver;
	private AssignmentMetricSummary assignmentMetrics;
	private ConsensusQualityReport consensusQuality;
	private Family[] initialSeeds;
	private MemberRankInputs memberInputs;
	private int memberLayer;
	private MemberRanking memberQuery;
	private MemberRankStatistics.Snapshot[] frozenMemberRanks;
	private final CandidateFilter.Context queryFilter=new CandidateFilter.Context();
	private final HbmPositionModel.ConsensusScorer scorer;
	private final MutableFamilyIndex index;
	private final MutableFamilyIndex.Workspace workspace=new MutableFamilyIndex.Workspace();
	private final ArrayList<Family> families=new ArrayList<Family>();
	private final int[] order, assignments, candidateIds;
	private int[] referenceSizes;
	private int[] previousAssignments;
	private byte[][] previousConsensuses;
	private ArrayList<long[]> passChanges;
	private int[] mergeInputWeights;
	private int mergeMaximumSize;
	private int completedPasses;
	private long mergeSizeFiltered,mergeSeedSizeFiltered,mergeReservationFailures;
	private final char[] querySketch;
	private final double[] queryMonomers=new double[20];
	private double queryNorm;
	private long assigned;
	private boolean started,completedRun,concurrentGreedy,orderConfigured;
	private boolean scavenge,scavengeCpuAvailable;
	private boolean merge;
	private TailAlignmentGate mergeGate;
	private ArrayList<MergeSummary> mergeSummaries;
	private long scavengeInputs,scavengeRescued,scavengeWallNanos,scavengeCpuNanos;
	private double snapshotFraction=0.02,stopChangeFraction;
	private long seedCollisions,versionCollisions,lateChecks,lateFailures,snapshotBuildNanos,snapshotRebuilds;
	private long discardedVersionProposals,discardedVersionPassingProposals;
}
