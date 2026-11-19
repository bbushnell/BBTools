package prot;

import java.util.Arrays;

import ddl.DDLFamilySearch;

/** Selects raw-DDL hard negative pairs in both directions from one pooled query view.
 * The same raw positive stream feeds query-local heaps and family-owner heaps;
 * serial shards are exact parity oracles, while production workers can stripe
 * the private family heap by family rank to avoid multiplying top-k storage.
 * @author Yelan
 */
public final class FamilyRawPairCollector{

	private FamilyRawPairCollector(){}

	/** Scans every query once and returns the sorted union of query-top and family-top pairs.
	 * Query order defines the global negative ordinal used for family-side tie breaks.
	 */
	public static Result collect(final DDLFamilySearch.Index index, final Query[] queries,
			final int queryTop, final int familyTop){
		return collect(index, queries, 0, queries==null ? 0 : queries.length, 0,
			queries==null ? 0 : queries.length, queryTop, familyTop, 1);
	}

	/** Scans one contiguous query shard. firstQueryOrdinal is the global ordinal
	 * of queries[from] so family-side top-k ties can merge exactly.
	 */
	public static Result collect(final DDLFamilySearch.Index index, final Query[] queries,
			final int from, final int to, final int firstQueryOrdinal, final int totalQueries,
			final int queryTop, final int familyTop){
		return collect(index, queries, from, to, firstQueryOrdinal, totalQueries,
			queryTop, familyTop, 1);
	}

	/** Scans one query tile and routes raw positives to exclusive family-rank
	 * owners. The owner stripes retain families*familyTop cells total.
	 */
	public static Result collect(final DDLFamilySearch.Index index, final Query[] queries,
			final int from, final int to, final int firstQueryOrdinal, final int totalQueries,
			final int queryTop, final int familyTop, final int familyOwners){
		if(index==null || queries==null || queryTop<1 || familyTop<1){
			throw new IllegalArgumentException("Raw pair collection requires an index, queries and positive top-k limits");
		}
		final int families=index.families();
		if(families<1){throw new IllegalArgumentException("Raw pair collection requires at least one family");}
		if(!index.supportsInversePostingRarity()){
			throw new IllegalArgumentException("Raw pair collection requires an unweighted indexed uint16 CSR");
		}
		if(from<0 || to<from || to>queries.length || firstQueryOrdinal<0 || totalQueries<1 || familyOwners<1 ||
				to==from || (long)firstQueryOrdinal+to-from>totalQueries){
			throw new IllegalArgumentException("Invalid pooled query shard bounds");
		}
		final int queryCount=to-from;
		final int effectiveQueryTop=Math.min(queryTop, families);
		final int effectiveFamilyTop=Math.min(familyTop, queryCount);
		final long maxQueryPairs=(long)queryCount*effectiveQueryTop;
		final long maxFamilyPairs=(long)families*effectiveFamilyTop;
		requireCapacity(maxQueryPairs, "query-side pair capacity");
		requireCapacity(maxFamilyPairs, "family-side pair capacity");
		requireCapacity(maxQueryPairs+maxFamilyPairs, "combined pair capacity");
		final long[] queryKeys=new long[queryCount];
		final int[] panels=new int[queryCount], localOrdinals=new int[queryCount];
		for(int q=from; q<to; q++){
			final Query query=queries[q];
			if(query==null || query.sketch==null || query.panel<0 || query.localOrdinal<0){
				throw new IllegalArgumentException("Invalid raw query "+q);
			}
			final int slot=q-from;
			panels[slot]=query.panel;
			localOrdinals[slot]=query.localOrdinal;
			queryKeys[slot]=panelLocal(query.panel, query.localOrdinal);
		}
		Arrays.sort(queryKeys);
		for(int i=1; i<queryKeys.length; i++){
			if(queryKeys[i]==queryKeys[i-1]){throw new IllegalArgumentException("Duplicate pooled query identity");}
		}
		final FamilyStripe[] familyHeaps=makeFamilyStripes(families, effectiveFamilyTop, familyOwners);
		final long[] queryPairs=new long[(int)maxQueryPairs];
		final QueryBatch batch=collectQueryFrontier(index, queries, from, to, firstQueryOrdinal,
			queryTop, new RawPositiveObserver(){
				@Override public void observe(int family, int queryOrdinal, int hits){
					familyHeaps[owner(family, families, familyHeaps.length)].offer(family, queryOrdinal, hits);
				}
			}, queryPairs, 0);
		final int queryPairCount=batch.nextOffset;
		final int familyPairCount=size(familyHeaps);
		final long[] familyPairs=new long[familyPairCount];
		final int[] familyScores=new int[familyPairCount];
		drain(familyHeaps, familyPairs, familyScores);
		final long[] pairs=union(queryPairs, queryPairCount, familyPairs, familyPairCount);
		return new Result(pairs, panels, localOrdinals, firstQueryOrdinal, totalQueries, families,
			effectiveQueryTop, effectiveFamilyTop, batch.positivePairs, queryPairCount, familyPairCount,
			Arrays.copyOf(queryPairs, queryPairCount), familyPairs, familyScores);
	}

	/** Scans a caller-owned query slice. firstQueryOrdinal is the global ordinal
	 * of queries[from]. This writes each query-side top-k pair into queryPairs
	 * starting at offset, and synchronously emits every raw positive so the
	 * caller can feed persistent family-owner heaps or queues.
	 */
	public static QueryBatch collectQueryFrontier(final DDLFamilySearch.Index index, final Query[] queries,
			final int from, final int to, final int firstQueryOrdinal, final int queryTop,
			final RawPositiveObserver observer, final long[] queryPairs, int offset){
		final int families=validateQueryFrontier(index, queries, from, to, firstQueryOrdinal,
			queryTop, queryPairs, offset);
		return collectQueryFrontier(index, queries, from, to, firstQueryOrdinal, queryTop,
			new DDLFamilySearch.Workspace(families), observer, queryPairs, offset);
	}

	/** Scans a query slice with caller-owned score scratch. firstQueryOrdinal is
	 * the global ordinal of queries[from].
	 */
	public static QueryBatch collectQueryFrontier(final DDLFamilySearch.Index index, final Query[] queries,
			final int from, final int to, final int firstQueryOrdinal, final int queryTop,
			final DDLFamilySearch.Workspace work, final RawPositiveObserver observer,
			final long[] queryPairs, int offset){
		final int families=validateQueryFrontier(index, queries, from, to, firstQueryOrdinal,
			queryTop, queryPairs, offset);
		if(work==null || work.scores.length!=families){
			throw new IllegalArgumentException("Raw query collection requires workspace coverage for all families");
		}
		final int effectiveQueryTop=Math.min(queryTop, families);
		final SingleTopK queryHeap=new SingleTopK(effectiveQueryTop);
		final int startOffset=offset;
		long positivePairs=0;
		for(int q=from; q<to; q++){
			final Query query=queries[q];
			if(query==null || query.sketch==null || query.panel<0 || query.localOrdinal<0){
				throw new IllegalArgumentException("Invalid raw query "+q);
			}
			final int ordinal=firstQueryOrdinal+q-from;
			index.score(query.sketch, work);
			queryHeap.clear();
			positivePairs+=work.size;
			for(int t=0; t<work.size; t++){
				final int family=work.touched[t];
				final int hits=work.scores[family];
				if(hits<=0){throw new IllegalStateException("Touched family with nonpositive raw hits: "+family);}
				queryHeap.offer(family, hits);
				if(observer!=null){observer.observe(family, ordinal, hits);}
			}
			offset=queryHeap.drain(ordinal, queryPairs, offset);
		}
		return new QueryBatch(startOffset, offset, positivePairs);
	}

	private static int validateQueryFrontier(final DDLFamilySearch.Index index, final Query[] queries,
			final int from, final int to, final int firstQueryOrdinal, final int queryTop,
			final long[] queryPairs, final int offset){
		if(index==null || queries==null || queryTop<1 || queryPairs==null){
			throw new IllegalArgumentException("Raw query collection requires an index, queries, top-k limit and output");
		}
		final int families=index.families();
		if(families<1){throw new IllegalArgumentException("Raw query collection requires at least one family");}
		if(!index.supportsInversePostingRarity()){
			throw new IllegalArgumentException("Raw query collection requires an unweighted indexed uint16 CSR");
		}
		if(from<0 || to<from || to>queries.length || firstQueryOrdinal<0 || offset<0){
			throw new IllegalArgumentException("Invalid raw query shard bounds");
		}
		final int queryCount=to-from;
		final int effectiveQueryTop=Math.min(queryTop, families);
		final long maxQueryPairs=(long)queryCount*effectiveQueryTop;
		requireCapacity(maxQueryPairs, "query-side pair capacity");
		if((long)offset+maxQueryPairs>queryPairs.length){
			throw new IllegalArgumentException("Raw query output is too small");
		}
		return families;
	}

	public static Result collectStriped(final DDLFamilySearch.Index index, final Query[] queries,
			final int queryTop, final int familyTop, final int familyOwners){
		return collect(index, queries, 0, queries==null ? 0 : queries.length, 0,
			queries==null ? 0 : queries.length, queryTop, familyTop, familyOwners);
	}

	/** Drains persistent family-owner heaps once and returns the sorted unique
	 * union of their family-frontier pairs and the caller-retained query pairs.
	 */
	public static Result finishStreaming(final int families, final int queryTop, final int familyTop,
			final long positivePairs, final long[] queryPairs, final int queryPairCount,
			final int[] panels, final int[] localOrdinals, final FamilyStripe[] familyHeaps){
		if(families<1 || queryTop<1 || familyTop<1 || queryPairs==null || panels==null ||
				localOrdinals==null || panels.length!=localOrdinals.length ||
				panels.length<1 || queryPairCount<0 || queryPairCount>queryPairs.length ||
				positivePairs<0 || familyHeaps==null || familyHeaps.length<1){
			throw new IllegalArgumentException("Invalid streamed raw pair inputs");
		}
		for(int q=0; q<panels.length; q++){
			if(panels[q]<0 || localOrdinals[q]<0){
				throw new IllegalArgumentException("Invalid pooled query identity at "+q);
			}
		}
		for(int i=0; i<queryPairCount; i++){
			final long pair=queryPairs[i];
			if(family(pair)<0 || family(pair)>=families || query(pair)<0 || query(pair)>=panels.length){
				throw new IllegalArgumentException("Invalid streamed query pair at "+i);
			}
		}
		final int effectiveFamilyTop=validateFamilyHeaps(families, familyTop, panels.length, familyHeaps);
		final int familyPairCount=size(familyHeaps);
		requireCapacity((long)queryPairCount+familyPairCount, "streamed combined pair capacity");
		final long[] familyPairs=new long[familyPairCount];
		final int[] familyScores=new int[familyPairCount];
		drain(familyHeaps, familyPairs, familyScores);
		for(int i=0; i<familyPairCount; i++){
			if(query(familyPairs[i])<0 || query(familyPairs[i])>=panels.length){
				throw new IllegalArgumentException("Invalid streamed family pair at "+i);
			}
		}
		final long[] pairs=union(queryPairs, queryPairCount, familyPairs, familyPairCount);
		final int effectiveQueryTop=Math.min(queryTop, families);
		return new Result(pairs, panels.clone(), localOrdinals.clone(), 0, panels.length,
			families, effectiveQueryTop, effectiveFamilyTop, positivePairs, queryPairCount,
			familyPairCount, Arrays.copyOf(queryPairs, queryPairCount), familyPairs, familyScores);
	}

	/** Reduces exact per-shard frontiers. The global family top-k is contained in
	 * the union of each shard's local family top-k, so raw hits need not be rescanned.
	 */
	public static Result merge(final Result[] shards, final int families, final int totalQueries,
			final int familyTop){
		if(shards==null || shards.length<1 || families<1 || totalQueries<1 || familyTop<1){
			throw new IllegalArgumentException("Raw pair merge requires shards, families and positive dimensions");
		}
		final int[] panels=new int[totalQueries], localOrdinals=new int[totalQueries];
		Arrays.fill(panels, -1); Arrays.fill(localOrdinals, -1);
		long queryCapacity=0, positivePairs=0;
		final int effectiveFamilyTop=Math.min(familyTop, totalQueries);
		final FamilyTopK familyHeaps=new FamilyTopK(families, effectiveFamilyTop);
		final long[] queryKeys=new long[totalQueries];
		for(Result shard:shards){
			if(shard==null || shard.families!=families || shard.totalQueries!=totalQueries){
				throw new IllegalArgumentException("Incompatible raw pair shard");
			}
			if(shard.queryLimit!=shards[0].queryLimit){
				throw new IllegalArgumentException("Raw pair shards use incompatible query frontiers");
			}
			if(shard.familyLimit<Math.min(effectiveFamilyTop, shard.panels.length)){
				throw new IllegalArgumentException("Raw pair shard family frontier was retained too shallowly");
			}
			queryCapacity+=shard.queryPairs;
			requireCapacity(queryCapacity, "merged query-side pair capacity");
			positivePairs+=shard.positivePairs;
			for(int i=0; i<shard.panels.length; i++){
				final int ordinal=shard.firstQuery+i;
				if(ordinal<0 || ordinal>=totalQueries || panels[ordinal]>=0){
					throw new IllegalArgumentException("Overlapping or out-of-range raw query shard");
				}
				panels[ordinal]=shard.panels[i];
				localOrdinals[ordinal]=shard.localOrdinals[i];
				queryKeys[ordinal]=panelLocal(shard.panels[i], shard.localOrdinals[i]);
			}
			for(int i=0; i<shard.familyFrontier.length; i++){
				familyHeaps.offer(family(shard.familyFrontier[i]), query(shard.familyFrontier[i]), shard.familyScores[i]);
			}
		}
		for(int i=0; i<panels.length; i++){
			if(panels[i]<0){throw new IllegalArgumentException("Missing raw query shard for global ordinal "+i);}
		}
		Arrays.sort(queryKeys);
		for(int i=1; i<queryKeys.length; i++){
			if(queryKeys[i]==queryKeys[i-1]){throw new IllegalArgumentException("Duplicate merged pooled query identity");}
		}
		final int familyPairCount=familyHeaps.size();
		requireCapacity(queryCapacity+familyPairCount, "merged combined pair capacity");
		final long[] queryPairs=new long[(int)queryCapacity];
		int offset=0;
		for(Result shard:shards){
			System.arraycopy(shard.queryFrontier, 0, queryPairs, offset, shard.queryFrontier.length);
			offset+=shard.queryFrontier.length;
		}
		final long[] familyPairs=new long[familyPairCount];
		final int[] familyScores=new int[familyPairCount];
		final int familyOffset=familyHeaps.drain(familyPairs, familyScores, 0);
		if(familyOffset!=familyPairs.length){throw new IllegalStateException("Merged family heap did not fill result array");}
		final long[] pairs=union(queryPairs, queryPairs.length, familyPairs, familyPairCount);
		return new Result(pairs, panels, localOrdinals, 0, totalQueries, families,
			shards[0].queryLimit, effectiveFamilyTop, positivePairs, (int)queryCapacity,
			familyPairCount, queryPairs, familyPairs, familyScores);
	}

	public static FamilyStripe[] makeFamilyStripes(final int families, final int familyTop,
			final int owners){
		if(families<1 || familyTop<1 || owners<1){throw new IllegalArgumentException("Invalid family owner dimensions");}
		final FamilyStripe[] heaps=new FamilyStripe[Math.min(owners, families)];
		for(int i=0; i<heaps.length; i++){
			final int first=(int)(((long)families*i+heaps.length-1)/heaps.length);
			final int last=(int)(((long)families*(i+1)+heaps.length-1)/heaps.length);
			heaps[i]=new FamilyStripe(first, last, familyTop);
		}
		return heaps;
	}

	public static int owner(final int family, final int families, final int owners){
		if(family<0 || families<1 || family>=families || owners<1){
			throw new IllegalArgumentException("Invalid family owner lookup");
		}
		final int ownerCount=Math.min(owners, families);
		return (int)((long)family*ownerCount/families);
	}

	private static int size(final FamilyStripe[] heaps){
		int sum=0;
		for(FamilyStripe heap:heaps){sum=Math.addExact(sum, heap.size());}
		return sum;
	}

	private static void drain(final FamilyStripe[] heaps, final long[] pairs, final int[] scores){
		int offset=0;
		for(FamilyStripe heap:heaps){offset=heap.drain(pairs, scores, offset);}
		if(offset!=pairs.length){throw new IllegalStateException("Family owners did not fill result array");}
	}

	private static long[] union(final long[] queryPairs, final int queryPairCount,
			final long[] familyPairs, final int familyPairCount){
		final long[] pairs=new long[queryPairCount+familyPairCount];
		System.arraycopy(queryPairs, 0, pairs, 0, queryPairCount);
		System.arraycopy(familyPairs, 0, pairs, queryPairCount, familyPairCount);
		Arrays.sort(pairs);
		int unique=0;
		for(int i=0; i<pairs.length; i++){
			if(i==0 || pairs[i]!=pairs[i-1]){pairs[unique++]=pairs[i];}
		}
		return Arrays.copyOf(pairs, unique);
	}

	public static long pack(final int family, final int query){
		if(family<0 || query<0){throw new IllegalArgumentException("Negative pair coordinate");}
		return ((long)family<<32)|(query&0xffffffffL);
	}

	public static int family(final long pair){return (int)(pair>>>32);}

	public static int query(final long pair){return (int)pair;}

	public interface RawPositiveObserver{
		void observe(int family, int queryOrdinal, int hits);
	}

	/** Streaming completion drains supplied heaps directly; stripe depth must
	 * exactly match the requested effective family frontier to preserve Result
	 * familyLimit semantics.
	 */
	private static int validateFamilyHeaps(final int families, final int familyTop,
			final int totalQueries, final FamilyStripe[] heaps){
		final int effectiveFamilyTop=Math.min(familyTop, totalQueries);
		int expectedFirst=0;
		for(int i=0; i<heaps.length; i++){
			final FamilyStripe heap=heaps[i];
			if(heap==null || heap.firstFamily!=expectedFirst || heap.lastFamily<=heap.firstFamily ||
					heap.lastFamily>families || Math.min(heap.limit(), totalQueries)!=effectiveFamilyTop){
				throw new IllegalArgumentException("Invalid streamed family stripe at "+i);
			}
			expectedFirst=heap.lastFamily;
		}
		if(expectedFirst!=families){throw new IllegalArgumentException("Streamed family stripes do not cover every family");}
		return effectiveFamilyTop;
	}

	private static long panelLocal(final int panel, final int local){
		return ((long)panel<<32)|(local&0xffffffffL);
	}

	private static void requireCapacity(final long cells, final String label){
		if(cells>Integer.MAX_VALUE){
			throw new IllegalArgumentException(label+" exceeds Java array bounds: "+cells);
		}
	}

	/** Borrowed immutable query row plus its pooled panel identity. */
	public static final class Query{
		public Query(final int panel_, final int localOrdinal_, final char[] sketch_){
			panel=panel_; localOrdinal=localOrdinal_; sketch=sketch_;
		}
		public final int panel, localOrdinal;
		public final char[] sketch;
	}

	/** Sorted unique union of the two directional hard-negative frontiers. */
	public static final class Result{
		private Result(final long[] pairs_, final int[] panels_, final int[] localOrdinals_,
				final int firstQuery_, final int totalQueries_, final int families_,
				final int queryLimit_, final int familyLimit_, final long positivePairs_,
				final int queryPairs_, final int familyPairs_, final long[] queryFrontier_,
				final long[] familyFrontier_, final int[] familyScores_){
			pairs=pairs_; panels=panels_; localOrdinals=localOrdinals_;
			firstQuery=firstQuery_; totalQueries=totalQueries_; families=families_;
			queryLimit=queryLimit_; familyLimit=familyLimit_; positivePairs=positivePairs_;
			queryPairs=queryPairs_; familyPairs=familyPairs_;
			queryFrontier=queryFrontier_==null ? new long[0] : queryFrontier_;
			familyFrontier=familyFrontier_==null ? new long[0] : familyFrontier_;
			familyScores=familyScores_==null ? new int[0] : familyScores_;
		}
		public int size(){return pairs.length;}
		public long pairAt(final int index){return pairs[index];}
		public int familyAt(final int index){return family(pairs[index]);}
		public int queryAt(final int index){return query(pairs[index]);}
		public int panelAt(final int index){return panels[queryAt(index)-firstQuery];}
		public int localOrdinalAt(final int index){return localOrdinals[queryAt(index)-firstQuery];}
		public int queryFrontierSize(){return queryFrontier.length;}
		public long queryFrontierAt(final int index){return queryFrontier[index];}
		public int familyFrontierSize(){return familyFrontier.length;}
		public long familyFrontierAt(final int index){return familyFrontier[index];}
		public final int families;
		public final long positivePairs;
		public final int queryPairs, familyPairs;
		private final long[] pairs;
		private final long[] queryFrontier, familyFrontier;
		private final int[] familyScores;
		private final int[] panels, localOrdinals;
		private final int firstQuery, totalQueries, queryLimit, familyLimit;
	}

	/** Frontier emitted by one streamed query slice. */
	public static final class QueryBatch{
		private QueryBatch(final int firstOffset_, final int nextOffset_, final long positivePairs_){
			firstOffset=firstOffset_; nextOffset=nextOffset_;
			pairs=nextOffset-firstOffset; positivePairs=positivePairs_;
		}
		public final int firstOffset, nextOffset, pairs;
		public final long positivePairs;
	}

	/** Exclusive owner for one contiguous family-rank stripe. Query workers can
	 * route raw positives by owner(family,totalFamilies,stripes.length) and feed
	 * each stripe through a bounded buffer owned by that stripe's thread.
	 */
	public static final class FamilyStripe{
		private FamilyStripe(final int firstFamily_, final int lastFamily_, final int limit_){
			heap=new FamilyTopK(firstFamily_, lastFamily_, limit_);
			firstFamily=firstFamily_; lastFamily=lastFamily_;
		}
		public boolean owns(final int family){return family>=firstFamily && family<lastFamily;}
		public void offer(final int family, final int ordinal, final int score){heap.offer(family, ordinal, score);}
		public int size(){return heap.size();}
		public int drain(final long[] out, final int[] outScores, final int offset){
			final int next=offset+heap.size();
			heap.drain(out, outScores, offset);
			return next;
		}
		public final int firstFamily, lastFamily;
		private final FamilyTopK heap;
		private int limit(){return heap.limit;}
	}

	private static final class SingleTopK{
		SingleTopK(final int limit_){
			limit=limit_; families=new int[limit]; scores=new int[limit];
		}
		void clear(){size=0;}
		void offer(final int family, final int score){
			if(size<limit){
				families[size]=family; scores[size]=score; siftUp(size++);
			}else if(better(score, family, scores[0], families[0])){
				families[0]=family; scores[0]=score; siftDown(0);
			}
		}
		int drain(final int query, final long[] out, final int offset){
			final int n=size;
			for(int i=0; i<n; i++){out[offset+i]=pack(families[i], query);}
			size=0;
			return offset+n;
		}
		private void siftUp(int child){
			while(child>0){
				final int parent=(child-1)>>>1;
				if(!worse(child, parent)){break;}
				swap(child, parent); child=parent;
			}
		}
		private void siftDown(int parent){
			while(true){
				final int left=parent*2+1;
				if(left>=size){break;}
				final int right=left+1;
				int worst=left;
				if(right<size && worse(right, left)){worst=right;}
				if(!worse(worst, parent)){break;}
				swap(parent, worst); parent=worst;
			}
		}
		private boolean worse(final int a, final int b){
			return scores[a]!=scores[b] ? scores[a]<scores[b] : families[a]>families[b];
		}
		private void swap(final int a, final int b){
			final int tf=families[a], ts=scores[a];
			families[a]=families[b]; scores[a]=scores[b];
			families[b]=tf; scores[b]=ts;
		}
		private final int limit;
		private final int[] families, scores;
		private int size;
	}

	private static final class FamilyTopK{
		FamilyTopK(final int families_, final int limit_){
			this(0, families_, limit_);
		}
		FamilyTopK(final int firstFamily_, final int lastFamily_, final int limit_){
			firstFamily=firstFamily_; lastFamily=lastFamily_; limit=limit_;
			if(firstFamily<0 || lastFamily<firstFamily || limit<1){
				throw new IllegalArgumentException("Invalid family heap stripe");
			}
			final int familyCount=lastFamily-firstFamily;
			final long cellCount=(long)familyCount*limit;
			requireCapacity(cellCount, "family heap stripe capacity");
			final int cells=(int)cellCount;
			sizes=new int[familyCount]; ordinals=new int[cells]; scores=new int[cells];
		}
		void offer(final int family, final int ordinal, final int score){
			final int index=family-firstFamily;
			if(index<0 || family>=lastFamily){throw new IllegalArgumentException("Family outside heap stripe");}
			final int offset=index*limit;
			int size=sizes[index];
			if(size<limit){
				final int cell=offset+size;
				ordinals[cell]=ordinal; scores[cell]=score; sizes[index]=++size;
				siftUp(offset, size-1);
			}else if(better(score, ordinal, scores[offset], ordinals[offset])){
				ordinals[offset]=ordinal; scores[offset]=score;
				siftDown(offset, size, 0);
			}
		}
		int size(){
			int sum=0;
			for(int size:sizes){sum=Math.addExact(sum, size);}
			return sum;
		}
		int drain(final long[] out, final int[] outScores, int offset){
			for(int index=0; index<sizes.length; index++){
				final int size=sizes[index];
				final int base=index*limit, family=firstFamily+index;
				for(int i=0; i<size; i++){
					out[offset]=pack(family, ordinals[base+i]);
					if(outScores!=null){outScores[offset]=scores[base+i];}
					offset++;
				}
				sizes[index]=0;
			}
			return offset;
		}
		private void siftUp(final int offset, int child){
			while(child>0){
				final int parent=(child-1)>>>1;
				if(!worse(offset+child, offset+parent)){break;}
				swap(offset+child, offset+parent); child=parent;
			}
		}
		private void siftDown(final int offset, final int size, int parent){
			while(true){
				final int left=parent*2+1;
				if(left>=size){break;}
				final int right=left+1;
				int worst=left;
				if(right<size && worse(offset+right, offset+left)){worst=right;}
				if(!worse(offset+worst, offset+parent)){break;}
				swap(offset+parent, offset+worst); parent=worst;
			}
		}
		private boolean worse(final int a, final int b){
			return scores[a]!=scores[b] ? scores[a]<scores[b] : ordinals[a]>ordinals[b];
		}
		private void swap(final int a, final int b){
			final int to=ordinals[a], ts=scores[a];
			ordinals[a]=ordinals[b]; scores[a]=scores[b];
			ordinals[b]=to; scores[b]=ts;
		}
		private final int firstFamily, lastFamily, limit;
		private final int[] sizes, ordinals, scores;
	}

	private static boolean better(final int scoreA, final int ordinalA,
			final int scoreB, final int ordinalB){
		return scoreA!=scoreB ? scoreA>scoreB : ordinalA<ordinalB;
	}
}
