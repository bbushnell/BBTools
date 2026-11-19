package ddl;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;

import shared.Shared;
import simd.SIMDFamilyProfiles;
import structures.IntHashMap;
import template.Accumulator;
import template.ThreadWaiter;

/** Single-threaded mutable family sketches with exact exhaustive and indexed ranking.
 * Each published row is an owned copy. Stable nonnegative IDs may have holes.
 * Frozen CSR postings are combined with versioned delta postings; replacement
 * invalidates every previous posting for that ID, including unchanged keys.
 * Neither cosine callbacks nor another thread may mutate/query this owner while ranking.
 * Parallel readers use a FrozenView with worker-local scratch and frozen cosine inputs.
 * @author Brian Bushnell, Fischl
 */
public final class MutableFamilyIndex{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Exhaustive mode maintains rows only; indexed mode also maintains postings. */
	public MutableFamilyIndex(final int width, final boolean indexed){
		this(width, indexed, 1<<20);
	}

	/** The delta is compacted after max(minDeltaRecords, frozenPostingCount/4) records. */
	public MutableFamilyIndex(final int width_, final boolean indexed_, final int minDeltaRecords_){
		require(width_>0 && width_<=DDLFamilyProfiles.MAX_BUCKETS && Integer.bitCount(width_)==1,
			"Width must be a power of two representable by the (bucket<<16)|value key");
		require(minDeltaRecords_>0 && minDeltaRecords_<=MAX_ARRAY, "Positive bounded delta threshold required");
		width=width_; indexed=indexed_; minDeltaRecords=minDeltaRecords_;
		deltaHeads=indexed ? new IntHashMap(16) : null;
		Arrays.fill(activePosition, -1);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/

	/** Publishes an owned copy; identical active rows require no posting changes. */
	public boolean upsert(final int id, final char[] sketch){
		require(!busy, "A ranking callback must not mutate its family index");
		require(id>=0 && id<MAX_ARRAY && sketch!=null && sketch.length==width,
			"A family needs a nonnegative stable ID and one uint16 value per bucket");
		ensureFamilies(id+1);
		if(Arrays.equals(rows[id], sketch)){return false;}
		require(versions[id]<Integer.MAX_VALUE, "Family version exhausted; wrapping would revive stale postings");
		require(mutationEpoch<Long.MAX_VALUE, "Score freshness epoch must not wrap and validate stale query scores");
		final char[] copy=sketch.clone();
		int occupied=0;
		for(char value : copy){if(value!=0){occupied++;}}
		if(indexed){ensureDelta((long)deltaSize+occupied);}
		if(rows[id]==null){activePosition[id]=size; activeIds[size++]=id;}
		livePostings+=occupied-occupancies[id];
		occupancies[id]=occupied; rows[id]=copy; versions[id]++; mutationEpoch++;
		if(indexed){
			for(int bucket=0; bucket<width; bucket++){
				final int value=copy[bucket];
				if(value==0){continue;}
				final int key=(bucket<<16)|value;
				deltaFamilies[deltaSize]=id;
				deltaVersions[deltaSize]=versions[id];
				deltaNext[deltaSize]=deltaHeads.get(key);
				deltaHeads.put(key, deltaSize++);
			}
			if(deltaSize>=Math.max(minDeltaRecords, frozenPostings.length/4)){compact();}
		}
		assert(livePostings>=0 && size>0) : "Owned live rows determine nonnegative posting and family counts";
		return true;
	}

	/** Removes an active family; historical postings remain invalid until compaction. */
	public boolean remove(final int id){
		require(!busy && id>=0, "Removal requires a nonnegative ID outside an active ranking call");
		if(id>=rows.length || rows[id]==null){return false;}
		require(mutationEpoch<Long.MAX_VALUE, "Score freshness epoch must not wrap and validate stale query scores");
		final int position=activePosition[id], last=activeIds[--size];
		activeIds[position]=last; activePosition[last]=position; activePosition[id]=-1;
		livePostings-=occupancies[id]; occupancies[id]=0; rows[id]=null; mutationEpoch++;
		// Reusing this ID increments its version in upsert; null rows are rejected meanwhile.
		assert(size>=0 && livePostings>=0) : "Removing one active row cannot remove another family or its postings";
		return true;
	}

	/** Rebuilds the frozen CSR from current rows and discards all historical delta entries.
	 * Exhaustive mode has no postings to compact. Peak build storage includes a dense cursor.
	 */
	public void compact(){
		require(!busy, "Compaction cannot change posting arrays during ranking");
		if(!indexed){return;}
		final PostingTable table=buildPostings(width, rows, activeIds, size, livePostings);
		frozenOffsets=table.offsets; frozenPostings=table.postings; frozenVersions=versions.clone();
		deltaHeads.clear(); deltaSize=0; compactions++;
		// Reuse bounded delta buffers across compactions instead of allocating per update.
	}

	/** Builds current-row CSR directly; snapshot construction never creates a mutable delta. */
	private static PostingTable buildPostings(final int width, final char[][] rows, final int[] activeIds,
			final int size, final long livePostings){
		require(livePostings>=0 && livePostings<=MAX_ARRAY, "Live CSR postings exceed Java array capacity");
		final int[] offsets=size==0 ? new int[0] : new int[(width<<16)+1];
		final int[] postings=new int[(int)livePostings];
		for(int i=0; i<size; i++){
			final char[] row=rows[activeIds[i]];
			for(int bucket=0; bucket<width; bucket++){
				if(row[bucket]!=0){offsets[((bucket<<16)|row[bucket])+1]++;}
			}
		}
		for(int i=1; i<offsets.length; i++){offsets[i]+=offsets[i-1];}
		if(offsets.length>0){
			assert(offsets[offsets.length-1]==postings.length) : "CSR prefix sum must cover each current nonzero register once";
			final int[] cursor=Arrays.copyOf(offsets, offsets.length-1);
			for(int i=0; i<size; i++){
				final int id=activeIds[i]; final char[] row=rows[id];
				for(int bucket=0; bucket<width; bucket++){
					if(row[bucket]!=0){postings[cursor[(bucket<<16)|row[bucket]]++]=id;}
				}
			}
		}
		return new PostingTable(offsets, postings);
	}

	/** Bin owners count and fill disjoint CSR ranges; no atomics or locks occur per posting. */
	private static PostingTable buildPostingsParallel(final int width,final char[][] rows,
			final int[] activeIds,final int size,final int threads){
		require(threads>0, "Frozen index construction needs positive worker count");
		if(size==0){return new PostingTable(new int[0],new int[0]);}
		final PostingBuild work=new PostingBuild(width,rows,activeIds,size,threads);
		work.runPhase(false);
		long total=0;
		for(PostingSlice slice:work.slices){
			slice.base=(int)total;
			total+=slice.count;
			require(total<=MAX_ARRAY, "Live CSR postings exceed Java array capacity");
		}
		work.postings=new int[(int)total];
		work.cursor=new int[work.offsets.length-1];
		work.runPhase(true);
		assert(work.offsets[work.offsets.length-1]==work.postings.length) :
			"Disjoint bin prefixes must cover every counted nonzero register exactly once";
		return new PostingTable(work.offsets,work.postings);
	}

	/** Coordinates two joined phases; partial arrays never escape a failed worker. */
	private static final class PostingBuild implements Accumulator<PostingWorker>{
		PostingBuild(final int width,final char[][] rows_,final int[] ids_,final int size_,final int threads){
			assert(width>0 && rows_!=null && ids_!=null && size_>0 && size_<=ids_.length && threads>0) :
				"Validated snapshot rows and a nonempty active roster must precede bin partitioning";
			rows=rows_; ids=ids_; size=size_;
			offsets=new int[(width<<16)+1];
			final int workers=Math.min(threads,width);
			slices=new PostingSlice[workers];
			for(int i=0; i<workers; i++){
				slices[i]=new PostingSlice(this,(int)((long)width*i/workers),(int)((long)width*(i+1)/workers));
			}
		}

		/** Thread.start failures and worker failures both join every thread that actually started. */
		void runPhase(final boolean fill){
			assert(failure==null) : "A failed counting phase cannot publish or begin filling a partial CSR";
			done=0;
			final ArrayList<PostingWorker> workers=new ArrayList<PostingWorker>(slices.length);
			for(PostingSlice slice:slices){workers.add(new PostingWorker(this,slice,fill));}
			try{
				if(!ThreadWaiter.startAndWait(workers,this)){fail(new IllegalStateException("Incomplete frozen index worker phase"));}
			}catch(Throwable cause){
				fail(cause);
				final ArrayList<PostingWorker> started=new ArrayList<PostingWorker>();
				for(PostingWorker worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null){throw new IllegalStateException("Parallel frozen index failed before snapshot publication",failure);}
			if(done!=slices.length){throw new IllegalStateException("Frozen index construction omitted a bin range");}
		}
		synchronized void fail(final Throwable cause){if(failure==null){failure=cause;}}
		@Override public void accumulate(final PostingWorker worker){
			if(!worker.success){fail(new IllegalStateException("Frozen index worker failed"));}
			done++;
		}
		@Override public boolean success(){return failure==null && done==slices.length;}
		@Override public ReadWriteLock rwlock(){return null;}// Joining precedes driver-only accumulation.
		final char[][] rows;
		final int[] ids,offsets;
		final int size;
		final PostingSlice[] slices;
		int[] postings,cursor;
		int done;
		volatile Throwable failure;
	}

	/** Owns offsets in (firstKey,lastKey] and cursor entries in [firstKey,lastKey). */
	private static final class PostingSlice{
		PostingSlice(final PostingBuild work_,final int first_,final int last_){
			assert(first_>=0 && last_>first_) : "Every frozen-index worker must own a nonempty bin range";
			work=work_; first=first_; last=last_;
		}
		/** Counts rows in the same active-ID order as the serial constructor, then computes local prefixes. */
		void count(){
			for(int i=0; i<work.size; i++){
				final char[] row=work.rows[work.ids[i]];
				for(int bucket=first; bucket<last; bucket++){
					if(row[bucket]!=0){work.offsets[((bucket<<16)|row[bucket])+1]++; count++;}
				}
			}
			require(count<=MAX_ARRAY, "A frozen CSR bin range exceeds Java array capacity");
			int sum=0;
			for(int key=(first<<16)+1,end=last<<16; key<=end; key++){sum+=work.offsets[key]; work.offsets[key]=sum;}
			assert(sum==count) : "Each local prefix must conserve its own nonzero-register count";
		}
		/** Adjacent owners never read each other's boundary offsets while global bases are being added. */
		void fill(){
			final int firstKey=first<<16,lastKey=last<<16;
			work.cursor[firstKey]=base;
			for(int key=firstKey+1; key<lastKey; key++){work.cursor[key]=base+work.offsets[key];}
			for(int key=firstKey+1; key<=lastKey; key++){work.offsets[key]+=base;}
			final int end=(int)(base+count);
			for(int i=0; i<work.size; i++){
				final int id=work.ids[i]; final char[] row=work.rows[id];
				for(int bucket=first; bucket<last; bucket++){
					if(row[bucket]!=0){
						final int position=work.cursor[(bucket<<16)|row[bucket]]++;
						assert(position>=base && position<end) : "A bin owner must fill only its counted posting slice";
						work.postings[position]=id;
					}
				}
			}
		}
		final PostingBuild work;
		final int first,last;
		long count;
		int base;
	}

	/** Worker failure reaches the coordinator before a FrozenView can reference the output arrays. */
	private static final class PostingWorker extends Thread{
		PostingWorker(final PostingBuild work_,final PostingSlice slice_,final boolean fill_){
			super("ddl-frozen-build-"+(fill_ ? "fill-" : "count-")+slice_.first);
			work=work_; slice=slice_; fill=fill_;
		}
		@Override public void run(){
			try{
				if(fill){slice.fill();}else{slice.count();}
				success=true;
			}catch(Throwable cause){work.fail(cause);}
		}
		final PostingBuild work;
		final PostingSlice slice;
		final boolean fill;
		boolean success;
	}

	/** Grows all ID-addressed state together; no family ID changes during growth. */
	private void ensureFamilies(final int needed){
		if(needed<=rows.length){return;}
		final int old=rows.length, capacity=grow(old, needed);
		rows=Arrays.copyOf(rows, capacity); versions=Arrays.copyOf(versions, capacity);
		occupancies=Arrays.copyOf(occupancies, capacity); activeIds=Arrays.copyOf(activeIds, capacity);
		activePosition=Arrays.copyOf(activePosition, capacity);
		Arrays.fill(activePosition, old, capacity, -1);
		assert(rows.length>=needed && versions.length==rows.length) : "Posting IDs address the same generation and row arrays";
	}

	/** Grows primitive delta columns once per capacity change, not once per posting. */
	private void ensureDelta(final long needed){
		require(needed<=MAX_ARRAY, "Mutable posting count exceeds Java array capacity");
		if(needed<=deltaFamilies.length){return;}
		final int capacity=grow(deltaFamilies.length, (int)needed);
		deltaFamilies=Arrays.copyOf(deltaFamilies, capacity);
		deltaVersions=Arrays.copyOf(deltaVersions, capacity);
		deltaNext=Arrays.copyOf(deltaNext, capacity);
		assert(deltaVersions.length==deltaNext.length) : "Every delta posting needs a version and linked-list successor";
	}

	/*--------------------------------------------------------------*/
	/*----------------       Scoring / Ranking      ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact matching-register counts by comparing every active row. Zero/zero never matches. */
	public void scoreFullScan(final char[] query, final Workspace work){
		validateQuery(query, work);
		scoreRows(query, work);
		work.scoredOwner=this; work.scoredEpoch=mutationEpoch;
	}

	/** Exact matching-register counts from valid frozen and delta postings. */
	public void scoreIndexed(final char[] query, final Workspace work){
		validateQuery(query, work);
		require(indexed, "This instance was created for exhaustive scoring only");
		scorePostings(query, work);
		work.scoredOwner=this; work.scoredEpoch=mutationEpoch;
	}

	/** Reference path: compares all sketches and evaluates every active family cosine. */
	public int rankFullScan(final char[] query, final int limit, final int[] outIds,
			final double[] outScores, final Cosine cosine, final Workspace work){
		return rank(query, limit, outIds, outScores, cosine, work, false);
	}

	/** Sparse exact path, including zero-hit completion whenever its bound can still win. */
	public int rankIndexed(final char[] query, final int limit, final int[] outIds,
			final double[] outScores, final Cosine cosine, final Workspace work){
		require(indexed, "This instance was created for exhaustive ranking only");
		return rank(query, limit, outIds, outScores, cosine, work, true);
	}

	/** Returns all requested winners, including zero-hit families, in score/ID order. */
	private int rank(final char[] query, final int limit, final int[] outIds,
			final double[] outScores, final Cosine cosine, final Workspace work, final boolean sparse){
		validateQuery(query, work);
		validateRanking(limit, outIds, outScores, cosine);
		busy=true; work.inUse=true;
		try{
			if(sparse){scorePostings(query, work);}else{scoreRows(query, work);}
			final int found=rankScores(activeIds, size, limit, outIds, outScores, cosine, work, sparse);
			work.scoredOwner=this; work.scoredEpoch=mutationEpoch;
			return found;
		}finally{work.inUse=false; busy=false;}
	}

	/** Read-only diagnostic rank from this owner's last completed query scores.
	 * The caller supplies the SAME query scoring state used for the shortlist and
	 * must not change callback-visible family state in between. Returns an exact
	 * one-based rank up to maxRank, otherwise maxRank+1 (censored, not exact).
	 * Diagnostic callbacks do not enter workspace ranking counters. Scores, heap,
	 * index rows/postings and the callers shortlist remain unchanged.
	 */
	public int rankOf(final int familyId, final Cosine cosine, final Workspace work, final int maxRank){
		require(!busy && work!=null && !work.inUse && cosine!=null && maxRank>0 && maxRank<Integer.MAX_VALUE,
			"Rank lookup needs exclusive completed scores, cosine and a positive nonoverflowing cap");
		require(familyId>=0 && familyId<rows.length && rows[familyId]!=null,
			"Rank lookup requires an active target family: "+familyId);
		require(work.scoredOwner==this && work.scoredEpoch==mutationEpoch,
			"Rank lookup scores must belong to this unchanged index after a completed query");
		busy=true; work.inUse=true;
		try{
			require(eligible(familyId, cosine, work, null),
				"Diagnostic rank target must pass this query's admission filter: "+familyId);
			final double target=rankingScore(familyId, work.scores[familyId], cosine);
			int rank=1;
			for(int i=0; i<size; i++){
				final int id=activeIds[i];
				if(id==familyId || !eligible(id, cosine, work, null)){continue;}
				final int hits=work.scores[id];
				// Bounds describe the complete score, including optional hit discounts.
				if(!better(upperBound(hits, cosine), id, target, familyId)){continue;}
				if(better(lowerBound(hits, cosine), id, target, familyId) ||
					better(rankingScore(id, hits, cosine), id, target, familyId)){
					if(++rank>maxRank){return maxRank+1;}
				}
			}
			assert(rank>=1 && rank<=size) : "A live target ranks within the unchanged active roster";
			return rank;
		}finally{work.inUse=false; busy=false;}
	}

	/** Uses the normal cosine contract without incrementing scientific ranking counters. */
	private static double diagnosticCosine(final int id, final Cosine cosine){
		final double value=cosine.score(id);
		if(!Double.isFinite(value) || value<0 || value>1){
			throw new IllegalArgumentException("Cosine must be finite in [0,1] for family "+id+": "+value);
		}
		return value;
	}

	/** Full-score callbacks own their arithmetic; ordinary cosine retains its exact addition. */
	private static double rankingScore(final int id, final int hits, final Cosine cosine){
		if(!(cosine instanceof Ranking)){return hits+diagnosticCosine(id, cosine);}
		final Ranking ranking=(Ranking)cosine;
		final double high=upperBound(hits, cosine), low=ranking.lowerBound(hits);
		final double value=ranking.score(id, hits);
		if(!Double.isFinite(value) || !Double.isFinite(low) || low>high || value<low || value>high){
			throw new IllegalArgumentException("Full family score must lie within its finite pruning bounds; family="+id+
				", hits="+hits+", score="+value+", lower="+low+", upper="+high);
		}
		return value;
	}

	/** Global bounds must cover every eligible family, including those never evaluated. */
	private static double upperBound(final int hits, final Cosine cosine){
		if(!(cosine instanceof Ranking)){return hits+1.0;}
		final Ranking ranking=(Ranking)cosine;
		final double bound=ranking.upperBound(hits), low=ranking.lowerBound(hits);
		if(!Double.isFinite(bound) || !Double.isFinite(low) || low>bound){
			throw new IllegalArgumentException("Ranking bounds must be finite and ordered before pruning; hits="+hits+
				", lower="+low+", upper="+bound);
		}
		return bound;
	}

	/** Length-discounted scores need not be at least their raw matching-register count. */
	private static double lowerBound(final int hits, final Cosine cosine){
		if(!(cosine instanceof Ranking)){return hits;}
		final double bound=((Ranking)cosine).lowerBound(hits);
		if(!Double.isFinite(bound) || bound>upperBound(hits, cosine)){
			throw new IllegalArgumentException("Ranking lower bound must be finite and no greater than its upper bound; hits="+hits);
		}
		return bound;
	}

	/** Applies the callback's complete-score bound with the final ascending-ID tie rule. */
	private static void offer(final int id, final int wanted, final Cosine cosine, final Workspace work, final boolean bounded){
		assert(id>=0 && id<work.scores.length) : "Only IDs from the current active roster/postings may enter the heap";
		if(bounded && work.heapSize==wanted && !better(upperBound(work.scores[id], cosine), id, work.heapScores[0], work.heapIds[0])){
			work.boundSkips++; return;
		}
		final double value=rankingScore(id, work.scores[id], cosine);
		work.cosineCalls++;
		work.offer(id, value, wanted);
	}

	/** Validates before changing query scratch; callbacks cannot reenter this index. */
	private void validateQuery(final char[] query, final Workspace work){
		require(!busy && query!=null && query.length==width && work!=null && !work.inUse,
			"A non-reentrant query needs the index width and a workspace");
	}

	/** Scalar all-row oracle path; operates on the same immutable row copies as compaction. */
	private void scoreRows(final char[] query, final Workspace work){
		scoreRows(query, work, rows, activeIds, size, width);
	}

	/** Shared scalar kernel for the mutable owner and immutable query snapshots. */
	private static void scoreRows(final char[] query, final Workspace work, final char[][] rows,
			final int[] activeIds, final int size, final int width){
		work.begin(rows.length);
		for(int i=0; i<size; i++){
			final int id=activeIds[i]; final char[] row=rows[id];
			int count=0;
			for(int bucket=0; bucket<width; bucket++){
				if(query[bucket]!=0 && query[bucket]==row[bucket]){count++;}
			}
			if(count>0){work.scores[id]=count; work.touched[work.touchedSize++]=id;}
		}
		work.comparedRegisters=(long)size*width;
		assert(work.touchedSize<=size) : "A scalar row contributes at most one touched family";
	}

	/** Counts only the current version of each active row; stale postings cannot duplicate hits. */
	private void scorePostings(final char[] query, final Workspace work){
		work.begin(rows.length);
		for(int bucket=0; bucket<width; bucket++){
			if(query[bucket]==0){continue;}
			final int key=(bucket<<16)|query[bucket];
			if(frozenOffsets.length>0){
				final int end=frozenOffsets[key+1];
				for(int p=frozenOffsets[key]; p<end; p++){
					final int id=frozenPostings[p]; work.postingsVisited++;
					if(rows[id]!=null && versions[id]==frozenVersions[id]){work.hit(id);}
				}
			}
			for(int p=deltaHeads.get(key); p>=0; p=deltaNext[p]){
				final int id=deltaFamilies[p]; work.postingsVisited++;
				if(rows[id]!=null && versions[id]==deltaVersions[p]){work.hit(id);}
			}
		}
		assert(work.touchedSize<=size) : "Version-filtered postings count only live family rows";
	}

	/** Publishes a read-only query snapshot before worker startup. Only row references and
	 * the active-ID vector are copied; immutable owned rows and CSR arrays are shared.
	 * Pending delta/stale postings are compacted first. A clean CSR is not rebuilt again.
	 * Drop the view and worker references after join so future updates can reclaim old arrays.
	 */
	public FrozenView freeze(){
		require(!busy, "Only the serial owner may freeze an index outside a ranking callback");
		if(indexed && (deltaSize>0 || frozenPostings.length!=livePostings)){compact();}
		return new FrozenView(width, indexed, rows.clone(), Arrays.copyOf(activeIds, size), frozenOffsets, frozenPostings);
	}

	/** Builds an indexed view from current rows even when the owner maintains no mutable postings.
	 * This does not change owner storage, deltas or ordinary freeze/ranking behavior.
	 */
	public FrozenView freezeIndexed(final int buildThreads){
		require(!busy && buildThreads>0, "Freeze indexed rows outside a ranking callback with positive build threads");
		final char[][] captured=rows.clone();
		final int[] ids=Arrays.copyOf(activeIds,size);
		final PostingTable table=buildThreads==1 ? buildPostings(width,captured,ids,size,livePostings) :
			buildPostingsParallel(width,captured,ids,size,buildThreads);
		return new FrozenView(width,true,captured,ids,table.offsets,table.postings);
	}

	/** Captures a dense prefix of immutable published sketches and builds its CSR once.
	 * A source may publish replacement rows while this runs, but may never mutate a
	 * returned array or remove a published ID. Each row is read once, then retained.
	 */
	public static FrozenView snapshot(final int width, final int count, final SketchSource source){
		require(width>0 && width<=DDLFamilyProfiles.MAX_BUCKETS && Integer.bitCount(width)==1 &&
			count>=0 && count<=MAX_ARRAY && source!=null, "Snapshot requires bounded width/count and immutable sketch source");
		final char[][] rows=new char[count][];
		final int[] ids=new int[count];
		long postings=0;
		for(int id=0; id<count; id++){
			final char[] row=source.sketch(id);
			require(row!=null && row.length==width, "Published prefix must contain an immutable exact-width sketch in every slot");
			rows[id]=row; ids[id]=id;
			for(char value : row){if(value!=0){postings++;}}
		}
		final PostingTable table=buildPostings(width, rows, ids, count, postings);
		return new FrozenView(width, true, rows, ids, table.offsets, table.postings);
	}

	/** Optional bin-partitioned bulk construction; the ordinary three-argument path stays serial.
	 * Row capture is serial and read-once; register counting/prefix/fill are parallel.
	 */
	public static FrozenView snapshot(final int width,final int count,final SketchSource source,final int buildThreads){
		require(buildThreads>0, "Snapshot construction needs positive build threads");
		if(buildThreads==1){return snapshot(width,count,source);}
		require(width>0 && width<=DDLFamilyProfiles.MAX_BUCKETS && Integer.bitCount(width)==1 &&
			count>=0 && count<=MAX_ARRAY && source!=null, "Snapshot requires bounded width/count and immutable sketch source");
		final char[][] rows=new char[count][];
		final int[] ids=new int[count];
		for(int id=0; id<count; id++){
			final char[] row=source.sketch(id);
			require(row!=null && row.length==width, "Published prefix must contain an immutable exact-width sketch in every slot");
			rows[id]=row; ids[id]=id;
		}
		final PostingTable table=buildPostingsParallel(width,rows,ids,count,buildThreads);
		return new FrozenView(width,true,rows,ids,table.offsets,table.postings);
	}

	/** Acquire-visible dense family prefix; row arrays are immutable after publication.
	 * size() is monotonic for the lifetime of a greedy pass. No slot below size() is null.
	 */
	public interface SketchSource{
		int size();
		char[] sketch(int familyId);
	}

	/** One construction result; arrays become immutable when published in a view. */
	private static final class PostingTable{
		private PostingTable(final int[] offsets_, final int[] postings_){
			offsets=offsets_; postings=postings_;
		}
		private final int[] offsets, postings;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Observability / API    ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of currently active families, including all-zero sketches. */
	public int size(){return size;}
	public int width(){return width;}
	public boolean indexed(){return indexed;}
	public long compactions(){return compactions;}
	public long livePostings(){return livePostings;}
	public int deltaRecords(){return deltaSize;}
	public int frozenPostingCount(){return frozenPostings.length;}

	/** Returns the maintained nonzero-register count for an active family in constant time.
	 * Inactive IDs are rejected so a removed family cannot masquerade as an empty sketch.
	 */
	public int occupiedBins(final int familyId){
		if(familyId<0 || familyId>=rows.length || rows[familyId]==null){
			throw new IllegalArgumentException("Occupied-bin lookup requires an active family ID: "+familyId);
		}
		final int count=occupancies[familyId];
		assert(count>=0 && count<=width) : "Upsert counts at most one occupied value per sketch register";
		return count;
	}

	/** Returns a detached row, or null for an absent ID. */
	public char[] profileCopy(final int id){
		require(id>=0, "Family IDs cannot be negative");
		return id>=rows.length || rows[id]==null ? null : rows[id].clone();
	}

	/** Allocated primitive-array payload only, not heap usage or a minimum-Xmx estimate.
	 * Excludes object headers, reference arrays, callback state and Family-owned copies.
	 */
	public long primitivePayloadBytes(){
		long bytes=4L*((long)versions.length+occupancies.length+activeIds.length+activePosition.length+
			frozenOffsets.length+frozenPostings.length+(long)frozenVersions.length+
			deltaFamilies.length+deltaVersions.length+(long)deltaNext.length)+2L*size*width;
		if(deltaHeads!=null){bytes+=4L*(deltaHeads.keys().length+(long)deltaHeads.values().length);}
		assert(bytes>=0) : "Array payload accounting uses long arithmetic before summing capacities";
		return bytes;
	}

	/** Supplies the current family cosine for this query; clamp floating roundoff to 1.
	 * Every active family must define a finite [0,1] value; zero norms return zero.
	 * The indexed path may omit callbacks whose upper bound cannot enter the shortlist.
	 */
	public interface Cosine{
		double score(int familyId);
	}

	/** Optional complete ranking score in raw-hit units, with query-stable global bounds.
	 * score(id,hits) includes all terms: the index never adds hits a second time.
	 * Bounds must enclose every active family's possible score at that hit count,
	 * even if its callback is skipped. Bad bounds can hide candidates; callers must
	 * establish them analytically and compare indexed ranking to an exhaustive oracle.
	 * Raw Workspace.scoreFor values remain unchanged for admission screening.
	 * Extending Cosine keeps every existing rank/frozen/tail signature compatible;
	 * implementations still provide score(id) for their ordinary cosine component.
	 */
	public interface Ranking extends Cosine{
		/** Raw-hit admission precedes score evaluation and shortlist capacity. */
		default boolean accept(int familyId, int hits){return true;}
		double score(int familyId, int hits);
		double lowerBound(int hits);
		double upperBound(int hits);
	}

	/** Query-local eligibility over snapshot IDs, applied before shortlist admission.
	 * The predicate must remain stable for a query and must not reenter its workspace.
	 * It may reject self or incompatible family sizes without consuming a shortlist slot.
	 */
	public interface Eligibility{
		boolean accept(int familyId);
	}

	/** Immutable frozen-pass index. Each worker supplies its OWN Workspace, cosine state,
	 * query and output arrays. The caller keeps all callback-visible family statistics frozen
	 * until workers join. This object has no owner reference, mutable delta or shared query flag.
	 */
	public static final class FrozenView{

		/** Captures exclusively owned outer references and already immutable row/posting payloads. */
		private FrozenView(final int width_, final boolean indexed_, final char[][] rows_,
				final int[] activeIds_, final int[] offsets_, final int[] postings_){
			width=width_; indexed=indexed_; rows=rows_; activeIds=activeIds_; offsets=offsets_; postings=postings_;
			boolean dense=true;
			for(int i=0; i<activeIds.length && dense; i++){dense=activeIds[i]==i;}
			contiguous=dense;
			assert(activeIds.length<=rows.length && (offsets.length==0 || offsets[offsets.length-1]==postings.length)) :
				"Freeze publishes only a complete current-row CSR and a bounded active roster";
		}

		/** Same exhaustive ordering as the serial owner, using this snapshot only. */
		public int rankFullScan(final char[] query, final int limit, final int[] outIds,
				final double[] outScores, final Cosine cosine, final Workspace work){
			return rank(query, limit, outIds, outScores, cosine, work, false);
		}

		/** Same exact sparse ordering as the serial owner, with worker-local query state. */
		public int rankIndexed(final char[] query, final int limit, final int[] outIds,
				final double[] outScores, final Cosine cosine, final Workspace work){
			require(indexed, "This snapshot contains exhaustive rows only");
			return rank(query, limit, outIds, outScores, cosine, work, true);
		}

		/** Ranks only eligible IDs; null retains the ordinary unfiltered ordering.
		 * Eligibility precedes heap admission and zero-hit completion. Fewer than limit
		 * eligible families return fewer winners, with unused outputs cleared as usual.
		 */
		public int rankIndexed(final char[] query, final int limit, final int[] outIds,
				final double[] outScores, final Cosine cosine, final Workspace work, final Eligibility eligible){
			require(indexed, "This snapshot contains exhaustive rows only");
			return rank(query, limit, outIds, outScores, cosine, work, true, eligible);
		}

		/** Exhaustive reference for the same eligibility and exact score/ID ordering. */
		public int rankFullScan(final char[] query, final int limit, final int[] outIds,
				final double[] outScores, final Cosine cosine, final Workspace work, final Eligibility eligible){
			return rank(query, limit, outIds, outScores, cosine, work, false, eligible);
		}

		/** Ranks one captured prefix plus every newer family up to the captured count.
		 * Prefix rows intentionally remain stale; tail rows are acquired once per ID.
		 * Both populations feed the same heap and exact zero-hit completion rule.
		 */
		public int rankWithTail(final char[] query, final int count, final SketchSource source,
				final int limit, final int[] outIds, final double[] outScores, final Cosine cosine, final Workspace work){
			require(indexed && contiguous && count>=activeIds.length && count<=MAX_ARRAY && source!=null,
				"Tail ranking requires a dense indexed prefix and a monotonic published family count");
			require(query!=null && query.length==width && work!=null && !work.inUse,
				"Tail ranking needs matching width and exclusive non-reentrant scratch");
			validateRanking(limit, outIds, outScores, cosine);
			work.inUse=true;
			try{
				work.begin(Math.max(rows.length, count));
				addPostings(query, work);
				for(int id=activeIds.length; id<count; id++){
					final char[] row=source.sketch(id);
					require(row!=null && row.length==width, "Published tail must contain an immutable exact-width sketch in every slot");
					final int hits=Shared.SIMD ? SIMDFamilyProfiles.matches(query, row) : scalarMatches(query, row);
					if(hits>0){work.scores[id]=hits; work.touched[work.touchedSize++]=id;}
				}
				work.comparedRegisters=(long)(count-activeIds.length)*width;
				return rankScores(null, count, limit, outIds, outScores, cosine, work, true);
			}finally{work.inUse=false;}
		}

		/** Workspace reentrancy protection is thread-local by contract, never a snapshot-wide flag. */
		private int rank(final char[] query, final int limit, final int[] outIds,
				final double[] outScores, final Cosine cosine, final Workspace work, final boolean sparse){
			return rank(query, limit, outIds, outScores, cosine, work, sparse, null);
		}

		/** Predicate state belongs to the caller; exception recovery releases the same query guard. */
		private int rank(final char[] query, final int limit, final int[] outIds,
				final double[] outScores, final Cosine cosine, final Workspace work, final boolean sparse,
				final Eligibility eligible){
			require(query!=null && query.length==width && work!=null && !work.inUse,
				"Frozen ranking needs matching width and an exclusively owned, non-reentrant workspace");
			validateRanking(limit, outIds, outScores, cosine);
			work.inUse=true;
			try{
				if(sparse){scorePostings(query, work);}
				else{MutableFamilyIndex.scoreRows(query, work, rows, activeIds, activeIds.length, width);}
				return rankScores(activeIds, activeIds.length, limit, outIds, outScores, cosine, work, sparse, eligible);
			}finally{work.inUse=false;}
		}

		/** Frozen postings contain no stale versions; their arrays are never mutated after publication. */
		private void scorePostings(final char[] query, final Workspace work){
			work.begin(rows.length);
			addPostings(query, work);
		}

		/** Adds prefix hits to already cleared scratch, leaving room for a disjoint live tail. */
		private void addPostings(final char[] query, final Workspace work){
			if(offsets.length==0){return;}
			for(int bucket=0; bucket<width; bucket++){
				if(query[bucket]==0){continue;}
				final int key=(bucket<<16)|query[bucket], end=offsets[key+1];
				for(int p=offsets[key]; p<end; p++){work.postingsVisited++; work.hit(postings[p]);}
			}
			assert(work.touchedSize<=activeIds.length) : "A frozen posting can name only a snapshot-active family";
		}

		public int size(){return activeIds.length;}
		public int width(){return width;}
		public boolean indexed(){return indexed;}

		private final int width;
		private final boolean indexed, contiguous;
		private final char[][] rows;
		private final int[] activeIds, offsets, postings;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Workspace           ----------------*/
	/*--------------------------------------------------------------*/

	/** Reused per-query primitive state, exclusively owned by one worker.
	 * A workspace can follow index growth and removals; inUse rejects callback reentry,
	 * not concurrent access. Sharing a workspace between workers is unsupported.
	 */
	public static final class Workspace{

		/** Clears touched scores only; stale scores from removed families are cleared too. */
		private void begin(final int capacity){
			scoredOwner=null;
			for(int i=0; i<touchedSize; i++){scores[touched[i]]=0;}
			if(scores.length<capacity){
				// Live-tail watermarks can grow by one; avoid copying the whole workspace per new family.
				final int grown=grow(scores.length, capacity);
				scores=Arrays.copyOf(scores, grown); touched=Arrays.copyOf(touched, grown);
			}
			touchedSize=0; heapSize=0; rankedFamilies=0;
			postingsVisited=comparedRegisters=cosineCalls=boundSkips=zeroHitVisited=0;
			assert(scores.length==touched.length) : "Every nonzero family score needs room in the touched list";
		}

		/** Adds exactly one matching nonzero bucket for an active current row. */
		private void hit(final int id){
			assert(id>=0 && id<scores.length) : "Posting IDs must fit the workspace grown from the index";
			if(scores[id]++==0){touched[touchedSize++]=id;}
		}

		/** Resizes once when a caller raises its shortlist capacity. */
		private void prepareHeap(final int wanted){
			assert(wanted>=0) : "Shortlist capacity is clipped to the nonnegative live family count";
			if(heapIds.length<wanted){heapIds=new int[wanted]; heapScores=new double[wanted];}
			heapSize=0;
		}

		/** Maintains a worst-at-root heap with descending-score/ascending-ID final ordering. */
		private void offer(final int id, final double value, final int limit){
			assert(limit>0 && heapSize<=limit) : "A ranking offer requires a bounded positive-capacity heap";
			if(heapSize<limit){
				int child=heapSize++;
				while(child>0){
					final int parent=(child-1)>>>1;
					if(!better(heapScores[parent], heapIds[parent], value, id)){break;}
					heapScores[child]=heapScores[parent]; heapIds[child]=heapIds[parent]; child=parent;
				}
				heapScores[child]=value; heapIds[child]=id;
			}else if(better(value, id, heapScores[0], heapIds[0])){
				heapScores[0]=value; heapIds[0]=id; siftDown();
			}
		}

		/** Repairs the root after replacement or removal, choosing the worse child at ties. */
		private void siftDown(){
			assert(heapSize>0) : "Sift-down requires a remaining root";
			int parent=0; final int id=heapIds[0]; final double value=heapScores[0];
			while(parent<(heapSize>>>1)){
				int child=(parent<<1)+1;
				if(child+1<heapSize && better(heapScores[child], heapIds[child], heapScores[child+1], heapIds[child+1])){child++;}
				if(!better(value, id, heapScores[child], heapIds[child])){break;}
				heapScores[parent]=heapScores[child]; heapIds[parent]=heapIds[child]; parent=child;
			}
			heapScores[parent]=value; heapIds[parent]=id;
		}

		/** Drains worst-first into reverse output positions so the result is best-first. */
		private void drain(final int[] outIds, final double[] outScores){
			assert(outIds.length>=heapSize) : "Caller validated room for the complete shortlist";
			while(heapSize>0){
				final int target=--heapSize;
				outIds[target]=heapIds[0];
				if(outScores!=null){outScores[target]=heapScores[0];}
				if(heapSize>0){heapIds[0]=heapIds[heapSize]; heapScores[0]=heapScores[heapSize]; siftDown();}
			}
		}

		/** Exact hit count from the last completed scoring/ranking call. */
		public int scoreFor(final int familyId){
			require(familyId>=0, "Family IDs cannot be negative");
			return familyId<scores.length ? scores[familyId] : 0;
		}
		public int touchedFamilies(){return touchedSize;}
		public long postingsVisited(){return postingsVisited;}
		public long comparedRegisters(){return comparedRegisters;}
		/** Number of evaluated cosine or optional full-score callbacks, excluding diagnostic ranks. */
		public long cosineCalls(){return cosineCalls;}
		public long boundSkips(){return boundSkips;}
		public long zeroHitVisited(){return zeroHitVisited;}
		/** Captured population of the last completed rank, not its truncated shortlist size. */
		public int rankedFamilies(){return rankedFamilies;}

		private int[] scores=new int[0], touched=new int[0], heapIds=new int[0];
		private double[] heapScores=new double[0];
		private int touchedSize, heapSize, rankedFamilies;
		private long postingsVisited, comparedRegisters, cosineCalls, boundSkips, zeroHitVisited;
		private boolean inUse;
		private MutableFamilyIndex scoredOwner;
		private long scoredEpoch;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Shared argument contract for mutable and frozen rankers. */
	private static void validateRanking(final int limit, final int[] outIds, final double[] outScores, final Cosine cosine){
		require(limit>0 && outIds!=null && outIds.length>=limit && cosine!=null &&
			(outScores==null || outScores.length>=limit), "Ranking needs positive capacity, outputs and a cosine callback");
	}

	/** One exact heap/bound/zero-hit algorithm serves serial and frozen concurrent queries. */
	private static int rankScores(final int[] activeIds, final int size, final int limit, final int[] outIds,
			final double[] outScores, final Cosine cosine, final Workspace work, final boolean sparse){
		return rankScores(activeIds, size, limit, outIds, outScores, cosine, work, sparse, null);
	}

	/** Filters candidates before they can consume heap capacity or affect its pruning bound. */
	private static int rankScores(final int[] activeIds, final int size, final int limit, final int[] outIds,
			final double[] outScores, final Cosine cosine, final Workspace work, final boolean sparse,
			final Eligibility eligible){
		work.rankedFamilies=size;
		Arrays.fill(outIds, 0, limit, -1);
		if(outScores!=null){Arrays.fill(outScores, 0, limit, 0);}
		final int wanted=Math.min(limit, size);
		work.prepareHeap(wanted);
		if(wanted==0){return 0;}
		if(!sparse){
			for(int i=0; i<size; i++){
				final int id=activeIds==null ? i : activeIds[i];
				if(eligible(id, cosine, work, eligible)){offer(id, wanted, cosine, work, false);}
			}
		}else{
			for(int i=0; i<work.touchedSize; i++){
				final int id=work.touched[i];
				if(eligible(id, cosine, work, eligible)){offer(id, wanted, cosine, work, true);}
			}
			// Equality still requires canonical-ID completion, even for zero-hit scores.
			if(work.heapSize<wanted || work.heapScores[0]<=upperBound(0, cosine)){
				for(int i=0; i<size; i++){
					final int id=activeIds==null ? i : activeIds[i];
					if(work.scores[id]==0 && eligible(id, cosine, work, eligible)){
						work.zeroHitVisited++; offer(id, wanted, cosine, work, true);
					}
				}
			}
		}
		assert(work.heapSize<=wanted && (eligible!=null || cosine instanceof Ranking || work.heapSize==wanted)) :
			"Unfiltered completion fills the quota; filtered completion may have fewer eligible families";
		final int found=work.heapSize;
		work.drain(outIds, outScores);
		return found;
	}

	/** Optional raw-hit admission never computes the pivot for a rejected family. */
	private static boolean eligible(final int id, final Cosine cosine, final Workspace work, final Eligibility filter){
		return (filter==null || filter.accept(id)) &&
			(!(cosine instanceof Ranking) || ((Ranking)cosine).accept(id, work.scores[id]));
	}

	/** Integer oracle/fallback: empty/empty agreement contributes no evidence. */
	private static int scalarMatches(final char[] query, final char[] reference){
		assert(query.length==reference.length) : "A register comparison needs identical sketch geometry";
		int count=0;
		for(int b=0; b<query.length; b++){if(query[b]!=0 && query[b]==reference[b]){count++;}}
		return count;
	}

	/** Exact total order shared by the bound and final heap, without epsilon ties. */
	private static boolean better(final double a, final int idA, final double b, final int idB){
		return a!=b ? a>b : idA<idB;
	}

	/** Geometric growth without signed-int overflow; Java arrays reserve header headroom. */
	private static int grow(final int old, final int needed){
		require(needed>old && needed<=MAX_ARRAY, "Invalid primitive-array growth request");
		return (int)Math.min(MAX_ARRAY, Math.max((long)needed, Math.max(16L, old*2L)));
	}

	private static void require(final boolean ok, final String message){
		if(!ok){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private final int width, minDeltaRecords;
	private final boolean indexed;
	private final IntHashMap deltaHeads;
	private char[][] rows=new char[16][];
	private int[] versions=new int[16], occupancies=new int[16], activeIds=new int[16], activePosition=new int[16];
	private int[] frozenOffsets=new int[0], frozenPostings=new int[0], frozenVersions=new int[0];
	private int[] deltaFamilies=new int[0], deltaVersions=new int[0], deltaNext=new int[0];
	private int size, deltaSize;
	private long livePostings, compactions, mutationEpoch;
	private boolean busy;
	private static final int MAX_ARRAY=Integer.MAX_VALUE-8;
}
