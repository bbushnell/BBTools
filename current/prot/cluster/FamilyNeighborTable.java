package prot.cluster;

import java.lang.management.ManagementFactory;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;

import com.sun.management.OperatingSystemMXBean;

import ddl.FamilySketchBuilder;
import ddl.MutableFamilyIndex;
import prot.HbmPositionModel;
import template.Accumulator;
import template.ThreadWaiter;

/** Cached directed consensus comparisons for one immutable selected family roster.
 * Each row contains the exact native DDL shortlist, not the global closest
 * alignments. Consensus query sketches/monomers rank against family-modal target
 * sketches/monomers, as in MutualNearestMerge. Alignment metrics are separate
 * from the configured rank and are never assumed symmetric. Explicit member-ranked
 * entrypoints use frozen retained-member length/KD with the same model cosine;
 * the original entrypoints retain shared-hit-plus-cosine ordering.
 * Model payloads are borrowed only until build returns. The completed table
 * retains primitive results and copied IDs/counts, with no model or index owner.
 * @author Brian Bushnell, Fischl, Keqing
 */
public final class FamilyNeighborTable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds each selected row once, including models with all-zero sketches.
	 * maxSizeRatio=0 means explicitly unfiltered evidence, excluding only self.
	 * Otherwise the inclusive size ratio is applied BEFORE shortlist admission;
	 * zero-current-member models have empty rows and cannot be targets, but remain
	 * in the roster. Counts are current post-assignment counts, not historical sizes.
	 * Callers keep model payloads immutable through this joined construction.
	 * Geometry retains mode, denominator and span limits; acceptance cutoffs are
	 * replaced by zero so weak valid neighbors remain available for later policy.
	 */
	public static FamilyNeighborTable build(final Family[] models, final int[] currentMembers,
			final FamilySketchBuilder builder, final TailAlignmentGate geometry, final int limit,
			final double maxSizeRatio, final int indexThreads, final int queryThreads){
		final Construction work=new Construction(models, currentMembers, builder, geometry,
			limit, maxSizeRatio, indexThreads, queryThreads);
		work.run();
		return new FamilyNeighborTable(work);
	}

	/** Freezes descriptors from one original cache, including last-nonempty empty templates.
	 * Models stay immutable until workers join; alignment and eligibility are unchanged.
	 */
	public static FamilyNeighborTable buildMemberRanked(final Family[] models, final int[] currentMembers,
			final FamilySketchBuilder builder, final TailAlignmentGate geometry, final int limit,
			final double maxSizeRatio, final int indexThreads, final int queryThreads, final MemberRankInputs inputs){
		require(inputs!=null,"Member-ranked neighbor tables require the original ordinal cache");
		final Construction work=new Construction(models,currentMembers,builder,geometry,
			limit,maxSizeRatio,indexThreads,queryThreads,null,inputs);
		work.run();
		return new FamilyNeighborTable(work);
	}

	/** Diagnostic-only serial queries against the COMPLETE target roster and index.
	 * queryFamilyIds must be explicit, unique stable IDs; output keeps their input
	 * order. The result contains only these measured rows, never an incomplete
	 * FamilyNeighborTable. Load timing belongs to the caller. Preparation, index,
	 * scratch and query costs are separate. First calls are reported separately
	 * from repeated batches whose summed query wall AND process CPU reach the
	 * requested duration. minSeconds=0 runs functional parity only, not a timing
	 * benchmark; positive windows must be at least ten seconds. Result comparison
	 * and retained-row copying occur outside query timing intervals.
	 */
	public static Measurement measureRows(final Family[] models, final int[] currentMembers,
			final FamilySketchBuilder builder, final TailAlignmentGate geometry, final int limit,
			final double maxSizeRatio, final int indexThreads, final int[] queryFamilyIds,
			final double minSeconds){
		return measureRows(models,currentMembers,builder,geometry,limit,maxSizeRatio,indexThreads,queryFamilyIds,minSeconds,null);
	}

	/** Same full-roster diagnostic with explicit retained-member ranking; zero seconds is functional parity only. */
	public static Measurement measureMemberRankedRows(final Family[] models, final int[] currentMembers,
			final FamilySketchBuilder builder, final TailAlignmentGate geometry, final int limit,
			final double maxSizeRatio, final int indexThreads, final int[] queryFamilyIds,
			final double minSeconds, final MemberRankInputs inputs){
		require(inputs!=null,"Member-ranked row measurements require the original ordinal cache");
		return measureRows(models,currentMembers,builder,geometry,limit,maxSizeRatio,indexThreads,queryFamilyIds,minSeconds,inputs);
	}

	/** Shares construction and measurement accounting while keeping the legacy entrypoint explicit. */
	private static Measurement measureRows(final Family[] models, final int[] currentMembers,
			final FamilySketchBuilder builder, final TailAlignmentGate geometry, final int limit,
			final double maxSizeRatio, final int indexThreads, final int[] queryFamilyIds,
			final double minSeconds, final MemberRankInputs inputs){
		require(Double.isFinite(minSeconds) && (minSeconds==0 || minSeconds>=10) && minSeconds<=3600,
			"Diagnostic timing needs 10..3600 seconds, or explicit zero for functional checks only");
		final CpuClock clock=new CpuClock();
		long cpuStart=clock.now(), wallStart=System.nanoTime();
		final int[] queries=resolveQueries(models, queryFamilyIds);
		final Construction work=new Construction(models, currentMembers, builder, geometry,
			limit, maxSizeRatio, indexThreads, 1, clock, inputs);
		final long setupWall=elapsed(wallStart), setupCpu=clock.since(cpuStart);
		final long prepareWall=setupWall-work.indexWall, prepareCpu=setupCpu-work.indexCpu;
		require(prepareWall>=0 && prepareCpu>=0, "Nested index timing must fit within complete preparation timing");
		boolean searchable=false;
		int positives=0;
		for(int count:work.members){if(count>0){positives++;}}
		for(int row:queries){searchable|=work.view!=null && (maxSizeRatio==0 || work.members[row]>0);}
		require(minSeconds==0 || searchable, "A positive timing window requires a searchable query, not only empty-row skips");
		cpuStart=clock.now(); wallStart=System.nanoTime();
		final Worker worker=new Worker(work);
		final long scratchWall=elapsed(wallStart), scratchCpu=clock.since(cpuStart);
		final RowSample[] samples=new RowSample[queries.length];
		long coldWall=0, coldCpu=0;
		for(int i=0; i<queries.length; i++){
			final long attempted=worker.attempted, valid=worker.valid, visited=worker.visited, cosines=worker.cosines;
			cpuStart=clock.now(); wallStart=System.nanoTime();
			worker.process(queries[i]);
			final long wall=elapsed(wallStart), cpu=clock.since(cpuStart);
			coldWall=Math.addExact(coldWall, wall); coldCpu=Math.addExact(coldCpu, cpu);
			samples[i]=new RowSample(work, queries[i], worker.attempted-attempted, worker.valid-valid,
				worker.visited-visited, worker.cosines-cosines, wall, cpu);
		}
		long repeatedWall=0, repeatedCpu=0, repetitions=0;
		final long minimum=(long)Math.ceil(minSeconds*1e9);
		do{
			for(RowSample sample:samples){
				final long attempted=worker.attempted, valid=worker.valid, visited=worker.visited, cosines=worker.cosines;
				cpuStart=clock.now(); wallStart=System.nanoTime();
				worker.process(sample.sourceRow);
				final long wall=elapsed(wallStart), cpu=clock.since(cpuStart);
				sample.repeatedWall=Math.addExact(sample.repeatedWall, wall);
				sample.repeatedCpu=Math.addExact(sample.repeatedCpu, cpu);
				repeatedWall=Math.addExact(repeatedWall, wall); repeatedCpu=Math.addExact(repeatedCpu, cpu);
				sample.check(work, worker.attempted-attempted, worker.valid-valid, worker.visited-visited, worker.cosines-cosines);
			}
			repetitions=Math.addExact(repetitions, 1);
		}while(repeatedWall<minimum || repeatedCpu<minimum);
		return new Measurement(work, samples, positives, prepareWall, prepareCpu, scratchWall, scratchCpu,
			coldWall, coldCpu, repeatedWall, repeatedCpu, repetitions, minSeconds);
	}

	/** Maps requested IDs before building the expensive index, without changing caller arrays. */
	private static int[] resolveQueries(final Family[] models, final int[] ids){
		require(models!=null && ids!=null && ids.length>0, "Diagnostic queries require a roster and nonempty explicit stable-ID list");
		final int[] sortedIds=new int[models.length], requested=ids.clone();
		for(int i=0; i<models.length; i++){
			require(models[i]!=null, "Diagnostic target roster cannot contain null models");
			sortedIds[i]=models[i].familyId;
		}
		Arrays.sort(sortedIds); Arrays.sort(requested);
		for(int i=1; i<sortedIds.length; i++){require(sortedIds[i]!=sortedIds[i-1], "Diagnostic target IDs must be unique");}
		for(int i=1; i<requested.length; i++){require(requested[i]!=requested[i-1], "Diagnostic query IDs must be unique");}
		final int[] rows=new int[ids.length];
		for(int i=0; i<ids.length; i++){
			rows[i]=Arrays.binarySearch(sortedIds, ids[i]);
			require(rows[i]>=0, "A diagnostic query ID is absent from the full selected roster");
		}
		return rows;
	}

	/** Publishes arrays only after all workers joined and all rows completed successfully. */
	private FamilyNeighborTable(final Construction work){
		assert(work.failure==null && work.done==work.familyIds.length) :
			"A neighbor table may publish only complete rows from one frozen population";
		familyIds=work.familyIds; members=work.members; rowSizes=work.rowSizes;
		nearest=work.nearest; targets=work.targets; shared=work.shared; reasons=work.reasons;
		ranks=work.ranks; identities=work.identities; ratios=work.ratios;
		scores=work.scores; querySpans=work.querySpans; referenceSpans=work.referenceSpans;
		stride=work.stride; limit=work.limit; maxSizeRatio=work.maxSizeRatio;
		width=work.builder.width(); geometry=work.gate;
		attemptedAlignments=work.attempted; validAlignments=work.valid;
		postingsVisited=work.visited; cosineCalls=work.cosines;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Read-only API        ----------------*/
	/*--------------------------------------------------------------*/

	public int size(){return familyIds.length;}
	public int width(){return width;}
	public int shortlistLimit(){return limit;}
	public double maxSizeRatio(){return maxSizeRatio;}
	/** Immutable geometry used for measurement; its identity and ratio cutoffs are zero. */
	public TailAlignmentGate geometry(){return geometry;}
	/** Returns the stable-ID-ordered row, or -1 when the family is absent. */
	public int rowForFamilyId(final int id){
		final int row=Arrays.binarySearch(familyIds, id);
		return row<0 ? -1 : row;
	}
	public int familyId(final int row){checkRow(row); return familyIds[row];}
	public int currentMembers(final int row){checkRow(row); return members[row];}
	/** Only [start(row),end(row)) contains edges; padding between rows is not evidence. */
	public int start(final int row){checkRow(row); return row*stride;}
	public int end(final int row){checkRow(row); return row*stride+rowSizes[row];}
	public int edgeCount(final int row){checkRow(row); return rowSizes[row];}
	/** Nearest valid alignment within this shortlist: identity, raw ratio, then stable ID.
	 * Returns -1 for an empty row or a row with no span-valid alignment. Inspect
	 * edgeCount/reasons to distinguish these cases. No merge threshold or fallback
	 * is applied here, and a missing reverse edge has not been measured.
	 */
	public int nearestEdge(final int row){checkRow(row); return nearest[row];}
	public int targetRow(final int edge){checkEdge(edge); return targets[edge];}
	public int targetFamilyId(final int edge){return familyIds[targetRow(edge)];}
	public int sharedHits(final int edge){checkEdge(edge); return shared[edge];}
	/** Exact configured rank: legacy sharedHits+cosine or the retained-member pivot; never a distance. */
	public double rankScore(final int edge){checkEdge(edge); return ranks[edge];}
	public TailAlignmentDecision.Reason reason(final int edge){checkEdge(edge); return REASONS[reasons[edge]];}
	public double identity(final int edge){checkEdge(edge); return identities[edge];}
	/** Raw query-normalized ratio, deliberately not clamped to one. */
	public double scoreRatio(final int edge){checkEdge(edge); return ratios[edge];}
	public int rawScore(final int edge){checkEdge(edge); return scores[edge];}
	public int querySpan(final int edge){checkEdge(edge); return querySpans[edge];}
	public int referenceSpan(final int edge){checkEdge(edge); return referenceSpans[edge];}
	/** Allocated primitive payload, including unused row capacity; not a JVM heap requirement. */
	public long primitiveBytes(){return 16L*familyIds.length+45L*targets.length;}

	private void checkRow(final int row){
		if(row<0 || row>=familyIds.length){throw new IllegalArgumentException("Unknown neighbor-table row: "+row);}
	}
	private void checkEdge(final int edge){
		if(edge<0 || edge>=targets.length || edge%stride>=rowSizes[edge/stride]){
			throw new IllegalArgumentException("Unpopulated neighbor-table edge: "+edge);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Construction         ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows frozen payloads, owns the shared index and gathers joined worker counters. */
	private static final class Construction implements Accumulator<Worker>{
		/** Validates the frozen roster and reserves fixed primitive row capacity before starting workers. */
		Construction(final Family[] input, final int[] counts, final FamilySketchBuilder builder_,
				final TailAlignmentGate geometry, final int limit_, final double ratio,
				final int indexThreads, final int queryThreads){
			this(input, counts, builder_, geometry, limit_, ratio, indexThreads, queryThreads, null, null);
		}

		/** Optional clock measures only snapshot construction; ordinary builds do not read clocks. */
		Construction(final Family[] input, final int[] counts, final FamilySketchBuilder builder_,
				final TailAlignmentGate geometry, final int limit_, final double ratio,
				final int indexThreads, final int queryThreads, final CpuClock clock, final MemberRankInputs inputs){
			require(input!=null && counts!=null && input.length==counts.length && builder_!=null && geometry!=null,
				"Neighbor table needs models, current counts, sketch builder and pinned geometry");
			require(limit_>0 && limit_<Integer.MAX_VALUE && indexThreads>0 && queryThreads>0 &&
				Double.isFinite(ratio) && (ratio==0 || ratio>=1) &&
				geometry.mode!=TailAlignmentGate.Mode.FULL_QUERY_UNPINNED,
				"Neighbor table needs bounded positive workers/limit, pinned geometry and zero or >=1 size ratio");
			builder=builder_; limit=limit_; maxSizeRatio=ratio;
			hydroWeight=inputs==null ? MemberRankInputs.DEFAULT_HYDRO_WEIGHT : inputs.hydroWeight;
			gate=new TailAlignmentGate(geometry.mode, geometry.denominator, geometry.maxLength, geometry.minSpan, 0, 0, null, geometry.matrix, geometry.minPaired);
			final int n=input.length;
			final long[] order=new long[n];
			int positives=0;
			for(int i=0; i<n; i++){
				final Family model=input[i];
				require(model!=null && counts[i]>=0 && model.consensusUnsafe()!=null &&
					model.consensusUnsafe().length>=gate.minSpan && model.sketchUnsafe()!=null &&
					model.sketchUnsafe().length==builder.width(),
					"Every selected model needs nonnegative current count, valid consensus span and exact sketch width");
				order[i]=((long)model.familyId<<32)|(i&0xffffffffL);
				if(counts[i]>0){positives++;}
			}
			Arrays.sort(order);
			models=new Family[n]; familyIds=new int[n]; members=new int[n];
			memberRanks=inputs==null ? null : new MemberRankStatistics.Snapshot[n];
			for(int row=0; row<n; row++){
				final int original=(int)order[row];
				models[row]=input[original]; familyIds[row]=models[row].familyId; members[row]=counts[original];
				require(row==0 || familyIds[row]!=familyIds[row-1], "Neighbor-table family IDs must be unique");
				if(memberRanks!=null){
					require(models[row].memberRankingUses(inputs),
						"Neighbor descriptors must use the same original ordinal cache for family "+familyIds[row]);
					memberRanks[row]=models[row].memberRankingSnapshot();
				}
			}
			stride=Math.min(limit, Math.max(0, n-1));
			final long capacity=(long)n*stride;
			require(capacity<=Integer.MAX_VALUE-8, "Neighbor table exceeds primitive edge-array capacity");
			final int cells=(int)capacity;
			rowSizes=new int[n]; nearest=new int[n]; Arrays.fill(nearest, -1);
			targets=new int[cells]; shared=new int[cells]; reasons=new byte[cells];
			ranks=new double[cells]; identities=new double[cells]; ratios=new double[cells];
			scores=new int[cells]; querySpans=new int[cells]; referenceSpans=new int[cells];
			final boolean searchable=stride>0 && (ratio==0 || positives>=2);
			final long indexCpuStart=clock==null ? 0 : clock.now();
			final long indexWallStart=clock==null ? 0 : System.nanoTime();
			view=searchable ? MutableFamilyIndex.snapshot(builder.width(), n, new MutableFamilyIndex.SketchSource(){
				@Override public int size(){return models.length;}
				@Override public char[] sketch(final int row){return models[row].sketchUnsafe();}
			}, indexThreads) : null;
			indexWall=clock==null ? 0 : elapsed(indexWallStart);
			indexCpu=clock==null ? 0 : clock.since(indexCpuStart);
			workerCount=Math.min(n, searchable ? queryThreads : 1);
		}

		/** One joined batch pass fills independent rows; any failure prevents publication. */
		void run(){
			final ArrayList<Worker> workers=new ArrayList<Worker>(workerCount);
			for(int i=0; i<workerCount; i++){workers.add(new Worker(this));}
			try{
				if(!ThreadWaiter.startAndWait(workers, this)){fail(new IllegalStateException("Neighbor-table workers did not finish"));}
			}catch(Throwable cause){
				fail(cause);
				final ArrayList<Worker> started=new ArrayList<Worker>();
				for(Worker worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null){throw new IllegalStateException("Neighbor-table construction failed before publication", failure);}
			if(done!=models.length){throw new IllegalStateException("Neighbor table omitted a selected model row");}
		}
		/** Claims disjoint row batches without synchronization in the per-alignment path. */
		int claim(){return (int)Math.min(next.getAndAdd(BATCH), models.length);}
		synchronized void fail(final Throwable cause){if(failure==null){failure=cause;}}
		@Override public void accumulate(final Worker worker){
			if(!worker.success){fail(new IllegalStateException("Neighbor-table worker failed"));}
			done=Math.addExact(done, worker.done);
			attempted=Math.addExact(attempted, worker.attempted); valid=Math.addExact(valid, worker.valid);
			visited=Math.addExact(visited, worker.visited); cosines=Math.addExact(cosines, worker.cosines);
		}
		@Override public boolean success(){return failure==null && done==models.length;}
		@Override public ReadWriteLock rwlock(){return null;}// Joins precede caller-owned aggregation.

		final Family[] models;
		final MemberRankStatistics.Snapshot[] memberRanks;
		final FamilySketchBuilder builder;
		final TailAlignmentGate gate;
		final MutableFamilyIndex.FrozenView view;
		final int limit, stride, workerCount;
		final double maxSizeRatio,hydroWeight;
		final long indexWall, indexCpu;
		final int[] familyIds, members, rowSizes, nearest, targets, shared, scores, querySpans, referenceSpans;
		final double[] ranks, identities, ratios;
		final byte[] reasons;
		final AtomicLong next=new AtomicLong();
		volatile Throwable failure;
		int done;
		long attempted, valid, visited, cosines;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Diagnostic results     ----------------*/
	/*--------------------------------------------------------------*/

	/** Completed diagnostic, exposing only measured rows and the full target census.
	 * All durations are nonnegative integer nanoseconds. Query durations are sums
	 * of timed native row calls, excluding comparisons and report materialization.
	 * Loader cost is not included. No unmeasured row can be queried through this API.
	 */
	public static final class Measurement{
		private Measurement(final Construction work, final RowSample[] samples, final int positives,
				final long prepareWall, final long prepareCpu, final long scratchWall, final long scratchCpu,
				final long coldWall, final long coldCpu, final long repeatedWall, final long repeatedCpu,
				final long repetitions_, final double requestedSeconds){
			rows=new MeasuredRow[samples.length];
			for(int i=0; i<rows.length; i++){rows[i]=new MeasuredRow(samples[i]);}
			targetCount=work.models.length; positiveTargetCount=positives;
			width=work.builder.width(); shortlistLimit=work.limit; maxSizeRatio=work.maxSizeRatio;
			geometry=work.gate; minimumQuerySeconds=requestedSeconds; timingValid=requestedSeconds>0;
			indexBuilt=work.view!=null;
			prepareWallNanos=prepareWall; prepareCpuNanos=prepareCpu;
			indexWallNanos=work.indexWall; indexCpuNanos=work.indexCpu;
			scratchWallNanos=scratchWall; scratchCpuNanos=scratchCpu;
			coldWallNanos=coldWall; coldCpuNanos=coldCpu;
			repeatedWallNanos=repeatedWall; repeatedCpuNanos=repeatedCpu; repetitions=repetitions_;
			assert(repetitions>=1 && (!timingValid || repeatedWallNanos>=(long)Math.ceil(requestedSeconds*1e9) &&
				repeatedCpuNanos>=(long)Math.ceil(requestedSeconds*1e9))) : "A timing-valid diagnostic must meet both repeated query-clock windows";
		}
		public int measuredQueries(){return rows.length;}
		/** index refers to the explicit query list's order, not the full roster's row number. */
		public MeasuredRow row(final int index){
			if(index<0 || index>=rows.length){throw new IllegalArgumentException("Query was not measured: "+index);}
			return rows[index];
		}
		private final MeasuredRow[] rows;
		public final int targetCount, positiveTargetCount, width, shortlistLimit;
		public final double maxSizeRatio, minimumQuerySeconds;
		public final TailAlignmentGate geometry;
		public final boolean timingValid, indexBuilt;
		public final long prepareWallNanos, prepareCpuNanos, indexWallNanos, indexCpuNanos;
		public final long scratchWallNanos, scratchCpuNanos, coldWallNanos, coldCpuNanos;
		public final long repeatedWallNanos, repeatedCpuNanos, repetitions;
	}

	/** One explicitly queried model, with native first-call evidence and repeated-call costs.
	 * Metrics and candidates describe ONE query, not the sum over repetitions.
	 * An empty measured result is explicit; absent query IDs have no row at all.
	 */
	public static final class MeasuredRow{
		private MeasuredRow(final RowSample sample){
			familyId=sample.familyId; currentMembers=sample.currentMembers;
			consensusLength=sample.consensusLength; modelOccupiedBins=sample.modelOccupiedBins;
			zeroCountSkipped=sample.zeroCountSkipped; nearestFamilyId=sample.nearestFamilyId;
			targetIds=sample.targetIds; shared=sample.shared; reasons=sample.reasons;
			ranks=sample.ranks; identities=sample.identities; ratios=sample.ratios;
			scores=sample.scores; querySpans=sample.querySpans; referenceSpans=sample.referenceSpans;
			attemptedAlignments=sample.attempted; validAlignments=sample.valid;
			postingsVisited=sample.visited; cosineCalls=sample.cosines;
			firstWallNanos=sample.firstWall; firstCpuNanos=sample.firstCpu;
			repeatedWallNanos=sample.repeatedWall; repeatedCpuNanos=sample.repeatedCpu;
		}
		public int candidates(){return targetIds.length;}
		public int targetFamilyId(final int rank){checkRank(rank); return targetIds[rank];}
		public int sharedHits(final int rank){checkRank(rank); return shared[rank];}
		public double rankScore(final int rank){checkRank(rank); return ranks[rank];}
		public TailAlignmentDecision.Reason reason(final int rank){checkRank(rank); return REASONS[reasons[rank]];}
		public double identity(final int rank){checkRank(rank); return identities[rank];}
		public double scoreRatio(final int rank){checkRank(rank); return ratios[rank];}
		public int rawScore(final int rank){checkRank(rank); return scores[rank];}
		public int querySpan(final int rank){checkRank(rank); return querySpans[rank];}
		public int referenceSpan(final int rank){checkRank(rank); return referenceSpans[rank];}
		private void checkRank(final int rank){
			if(rank<0 || rank>=targetIds.length){throw new IllegalArgumentException("Unmeasured candidate rank: "+rank);}
		}
		public final int familyId, currentMembers, consensusLength, modelOccupiedBins, nearestFamilyId;
		public final boolean zeroCountSkipped;
		public final long attemptedAlignments, validAlignments, postingsVisited, cosineCalls;
		public final long firstWallNanos, firstCpuNanos, repeatedWallNanos, repeatedCpuNanos;
		private final int[] targetIds, shared, scores, querySpans, referenceSpans;
		private final double[] ranks, identities, ratios;
		private final byte[] reasons;
	}

	/** Private accumulation; evidence arrays are exclusively transferred to the immutable measured row. */
	private static final class RowSample{
		/** Copies only this measured row after its timers stop; full-roster padding never escapes. */
		RowSample(final Construction work, final int row, final long attempted_, final long valid_,
				final long visited_, final long cosines_, final long wall, final long cpu){
			sourceRow=row; familyId=work.familyIds[row]; currentMembers=work.members[row];
			consensusLength=work.models[row].consensusUnsafe().length;
			modelOccupiedBins=work.models[row].occupiedBins();
			zeroCountSkipped=work.maxSizeRatio>0 && currentMembers==0;
			final int count=work.rowSizes[row], first=row*work.stride, last=first+count;
			assert(attempted_==count && valid_>=0 && valid_<=attempted_ && visited_>=0 && cosines_>=0) :
				"Each retained diagnostic candidate must have one native alignment and nonnegative query-work deltas";
			final int nearest=work.nearest[row];
			nearestFamilyId=nearest<0 ? -1 : work.familyIds[work.targets[nearest]];
			targetIds=new int[count];
			for(int i=0; i<count; i++){targetIds[i]=work.familyIds[work.targets[first+i]];}
			shared=Arrays.copyOfRange(work.shared, first, last); reasons=Arrays.copyOfRange(work.reasons, first, last);
			ranks=Arrays.copyOfRange(work.ranks, first, last); identities=Arrays.copyOfRange(work.identities, first, last);
			ratios=Arrays.copyOfRange(work.ratios, first, last); scores=Arrays.copyOfRange(work.scores, first, last);
			querySpans=Arrays.copyOfRange(work.querySpans, first, last);
			referenceSpans=Arrays.copyOfRange(work.referenceSpans, first, last);
			attempted=attempted_; valid=valid_; visited=visited_; cosines=cosines_;
			firstWall=wall; firstCpu=cpu;
		}

		/** Exact comparison outside timed spans prevents a faster but changed query from passing. */
		void check(final Construction work, final long attempted_, final long valid_, final long visited_, final long cosines_){
			final int best=work.nearest[sourceRow];
			if(work.rowSizes[sourceRow]!=targetIds.length || attempted!=attempted_ || valid!=valid_ ||
				visited!=visited_ || cosines!=cosines_ || nearestFamilyId!=(best<0 ? -1 : work.familyIds[work.targets[best]])){
				throw new IllegalStateException("Repeated neighbor query changed its outcome or work counts: "+familyId);
			}
			for(int i=0; i<targetIds.length; i++){
				final int edge=sourceRow*work.stride+i;
				if(targetIds[i]!=work.familyIds[work.targets[edge]] || shared[i]!=work.shared[edge] || reasons[i]!=work.reasons[edge] ||
					Double.doubleToLongBits(ranks[i])!=Double.doubleToLongBits(work.ranks[edge]) ||
					Double.doubleToLongBits(identities[i])!=Double.doubleToLongBits(work.identities[edge]) ||
					Double.doubleToLongBits(ratios[i])!=Double.doubleToLongBits(work.ratios[edge]) ||
					scores[i]!=work.scores[edge] || querySpans[i]!=work.querySpans[edge] || referenceSpans[i]!=work.referenceSpans[edge]){
					throw new IllegalStateException("Repeated neighbor query changed candidate evidence: "+familyId+" rank="+i);
				}
			}
		}
		final int sourceRow, familyId, currentMembers, consensusLength, modelOccupiedBins, nearestFamilyId;
		final boolean zeroCountSkipped;
		final int[] targetIds, shared, scores, querySpans, referenceSpans;
		final double[] ranks, identities, ratios;
		final byte[] reasons;
		final long attempted, valid, visited, cosines, firstWall, firstCpu;
		long repeatedWall, repeatedCpu;
	}

	/** Same process-CPU source used by the existing real member-graph cost diagnostic. */
	private static final class CpuClock{
		CpuClock(){
			final java.lang.management.OperatingSystemMXBean nativeBean=ManagementFactory.getOperatingSystemMXBean();
			require(nativeBean instanceof OperatingSystemMXBean, "Process CPU timing requires the native JVM operating-system bean");
			bean=(OperatingSystemMXBean)nativeBean;
			now();
		}
		long now(){
			final long value=bean.getProcessCpuTime();
			if(value<0){throw new IllegalStateException("Process CPU clock is unavailable");}
			return value;
		}
		long since(final long start){
			final long delta=now()-start;
			if(delta<0){throw new IllegalStateException("Process CPU clock moved backward");}
			return delta;
		}
		final OperatingSystemMXBean bean;
	}
	private static long elapsed(final long start){
		final long delta=System.nanoTime()-start;
		if(delta<0){throw new IllegalStateException("Monotonic wall clock moved backward");}
		return delta;
	}

	/** Owns all ranking and native alignment scratch; writes only claimed table rows. */
	private static final class Worker extends Thread implements MutableFamilyIndex.Cosine, MutableFamilyIndex.Eligibility{
		/** Forks configuration only; ranking and scorer buffers belong to this worker. */
		Worker(final Construction owner_){
			super("ddl-family-neighbor");
			owner=owner_; builder=owner.builder.fork();
			scorer=new HbmPositionModel.ConsensusScorer(owner.gate.matrix);
			querySketch=new char[builder.width()]; candidateIds=new int[owner.stride];
			candidateRanks=new double[owner.stride];
			ranking=owner.memberRanks==null ? this : new MutableFamilyIndex.Ranking(){
				@Override public double score(final int target){return Worker.this.score(target);}
				@Override public double score(final int target,final int hits){
					final MemberRankStatistics.Snapshot query=owner.memberRanks[queryRow],other=owner.memberRanks[target];
					final double ratio=Math.min(query.medianLength,other.medianLength)/Math.max(query.medianLength,other.medianLength);
					final double cosine=Worker.this.score(target),similarity=1-Math.abs(query.meanKd-other.meanKd)/9;
					assert(hits>=0 && ratio>0 && ratio<=1 && cosine>=0 && cosine<=1 && similarity>=0 && similarity<=1) :
						"Frozen member descriptors and cosine must satisfy the complete-score shortlist bounds";
					return hits*Math.sqrt(Math.sqrt(ratio))+1.5*cosine+owner.hydroWeight*similarity;
				}
				@Override public double lowerBound(final int hits){return 0;}
				@Override public double upperBound(final int hits){return hits+1.5+owner.hydroWeight;}
			};
		}
		@Override public void run(){
			try{
				for(int first=owner.claim(); first<owner.models.length && owner.failure==null; first=owner.claim()){
					final int last=(int)Math.min((long)first+BATCH, owner.models.length);
					for(int row=first; row<last && owner.failure==null; row++){process(row); done++;}
				}
				success=owner.failure==null;
			}catch(Throwable cause){owner.fail(cause);}
		}

		/** Uses the configured shortlist pivot and unchanged geometry-only nearest rule. */
		private void process(final int row){
			queryRow=row;
			if(owner.view==null || owner.maxSizeRatio>0 && owner.members[row]==0){return;}
			final byte[] query=owner.models[row].consensusUnsafe();
			builder.sketchEncoded(query, querySketch);
			Arrays.fill(monomers, 0); norm=0;
			for(byte aa:query){if(aa>=0 && aa<20){monomers[aa]++;}}
			for(double count:monomers){norm+=count*count;}
			final int found=owner.view.rankIndexed(querySketch, owner.stride, candidateIds,
				candidateRanks, ranking, scratch, this);
			visited+=scratch.postingsVisited(); cosines+=scratch.cosineCalls();
			owner.rowSizes[row]=found;
			int best=-1;
			for(int rank=0; rank<found; rank++){
				final int target=candidateIds[rank], edge=row*owner.stride+rank;
				final byte[] reference=owner.models[target].consensusUnsafe();
				scorer.resetConsensus(reference);
				final TailAlignmentDecision decision=owner.gate.decide(query, reference, scorer);
				attempted++;
				if(!Double.isFinite(decision.identity) || !Double.isFinite(decision.scoreRatio) ||
					decision.identity<0 || decision.identity>1){
					throw new IllegalStateException("Native neighbor metrics must be finite with identity in [0,1]");
				}
				owner.targets[edge]=target; owner.shared[edge]=scratch.scoreFor(target);
				owner.ranks[edge]=candidateRanks[rank]; owner.reasons[edge]=(byte)decision.reason.ordinal();
				owner.identities[edge]=decision.identity; owner.ratios[edge]=decision.scoreRatio;
				owner.scores[edge]=decision.score; owner.querySpans[edge]=decision.querySpan;
				owner.referenceSpans[edge]=decision.referenceSpan;
				if(decision.accepted()){
					valid++;
					if(best<0 || decision.identity>owner.identities[best] || decision.identity==owner.identities[best] &&
						(decision.scoreRatio>owner.ratios[best] || decision.scoreRatio==owner.ratios[best] &&
						owner.familyIds[target]<owner.familyIds[owner.targets[best]])){best=edge;}
				}
			}
			owner.nearest[row]=best;
			assert(found<=owner.stride && (best<0 || best>=row*owner.stride && best<row*owner.stride+found)) :
				"Nearest alignment must belong to this row's complete ranked shortlist";
		}

		/** Self and count/ratio restrictions apply to positive and zero-hit admission alike. */
		@Override public boolean accept(final int target){
			assert(queryRow>=0) : "Eligibility is defined only within an active row query";
			if(target==queryRow){return false;}
			if(owner.maxSizeRatio==0){return true;}
			final int a=owner.members[queryRow], b=owner.members[target];
			return a>0 && b>0 && Math.max(a, b)/(double)Math.min(a, b)<=owner.maxSizeRatio;
		}

		/** Native query-consensus monomer counts versus saved family monomer frequencies. */
		@Override public double score(final int target){
			final float[] values=owner.models[target].monomersUnsafe();
			double dot=0, targetNorm=0;
			for(int i=0; i<20; i++){dot+=monomers[i]*values[i]; targetNorm+=values[i]*(double)values[i];}
			return norm==0 || targetNorm==0 ? 0 : Math.min(1, Math.max(0, dot/Math.sqrt(norm*targetNorm)));
		}

		final Construction owner;
		final MutableFamilyIndex.Cosine ranking;
		final FamilySketchBuilder builder;
		final char[] querySketch;
		final int[] candidateIds;
		final double[] candidateRanks, monomers=new double[20];
		final MutableFamilyIndex.Workspace scratch=new MutableFamilyIndex.Workspace();
		final HbmPositionModel.ConsensusScorer scorer;
		int queryRow=-1, done;
		double norm;
		long attempted, valid, visited, cosines;
		boolean success;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private static void require(final boolean valid, final String message){if(!valid){throw new IllegalArgumentException(message);}}
	private final int[] familyIds, members, rowSizes, nearest, targets, shared, scores, querySpans, referenceSpans;
	private final double[] ranks, identities, ratios;
	private final byte[] reasons;
	private final int stride, limit, width;
	private final double maxSizeRatio;
	private final TailAlignmentGate geometry;
	public final long attemptedAlignments, validAlignments, postingsVisited, cosineCalls;
	private static final TailAlignmentDecision.Reason[] REASONS=TailAlignmentDecision.Reason.values();
	private static final int BATCH=16;
}
