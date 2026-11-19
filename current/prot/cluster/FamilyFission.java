package prot.cluster;

import ddl.FamilySketchBuilder;
import ddl.MemberSketchIndex;
import prot.ProteinSequence;
import structures.IntList;
import structures.LongList;

/** DDL-guided single-family fission proposal.
 * Builds a bounded member graph for one frozen family, splits obvious disconnected
 * DDL neighborhoods into two balanced children, then rebuilds child consensuses
 * from exact original input ordinals. Recursion, stable new-ID allocation,
 * objective scoring and publication belong to the caller.
 * @author Brian Bushnell, Yelan
 */
public final class FamilyFission{

	public FamilyFission(final int minChildSize_, final int neighborLimit_,
			final int minHits_, final int buildThreads_){
		if(minChildSize_<1 || neighborLimit_<1 || minHits_<1 || buildThreads_<1){
			throw new IllegalArgumentException(
				"Fission needs positive child, neighbor, hit and build-thread limits");
		}
		minChildSize=minChildSize_; neighborLimit=neighborLimit_;
		minHits=minHits_; buildThreads=buildThreads_;
	}

	/** Proposes one binary split; inputs and the parent family are never mutated. */
	public Result split(final Family parent, final ProteinSequence[] sequences,
			final FamilySketchBuilder sketches, final TailAlignmentGate gate,
			final ConsensusRefiner refiner, final ClusteringMetrics metrics,
			final int rightFamilyId){
		return split(parent,sequences,sketches,gate,refiner,metrics,rightFamilyId,false);
	}

	/** Proposes one binary split and optionally captures per-child rebuild summaries. */
	public Result split(final Family parent, final ProteinSequence[] sequences,
			final FamilySketchBuilder sketches, final TailAlignmentGate gate,
			final ConsensusRefiner refiner, final ClusteringMetrics metrics,
			final int rightFamilyId, final boolean summarizeRebuild){
		return split(parent, sequences, sketches, gate, refiner, metrics,
			rightFamilyId, summarizeRebuild, false);
	}

	/** Captures rejected replacement models only for a request-local diagnostic.
	 * Ordinary status, children and skipped-right accounting retain their meaning.
	 * Audit-only right construction uses separate metrics and is never published.
	 */
	public Result split(final Family parent, final ProteinSequence[] sequences,
			final FamilySketchBuilder sketches, final TailAlignmentGate gate,
			final ConsensusRefiner refiner, final ClusteringMetrics metrics,
			final int rightFamilyId, final boolean summarizeRebuild, final boolean auditRejected){
		return split(parent, sequences, sketches, gate, refiner, metrics, rightFamilyId,
			summarizeRebuild, auditRejected, 0);
	}

	/** Opt-in central seeds and bounded competitive reassignment; zero preserves the historical path. */
	public Result split(final Family parent, final ProteinSequence[] sequences,
			final FamilySketchBuilder sketches, final TailAlignmentGate gate,
			final ConsensusRefiner refiner, final ClusteringMetrics metrics,
			final int rightFamilyId, final boolean summarizeRebuild, final boolean auditRejected,
			final int maxRounds){
		validate(parent, sequences, sketches, gate, refiner, metrics, rightFamilyId);
		if(maxRounds<0 || maxRounds>2){throw new IllegalArgumentException("Fission refinement rounds must be 0, 1 or 2");}
		final int size=parent.size();
		if(size<2L*minChildSize){
			return Result.rejected(Status.TOO_SMALL, 0, 0, 0, 0).withoutCandidates(auditRejected);
		}
		final int[] members=parent.memberIndexes();
		final char[][] rows=new char[size][];
		final FamilySketchBuilder localSketches=sketches.fork();
		for(int i=0; i<size; i++){
			final int member=members[i];
			if(member<0 || member>=sequences.length || sequences[member]==null){
				throw new IllegalArgumentException("Parent contains missing input ordinal "+member);
			}
			rows[i]=new char[localSketches.width()];
			localSketches.sketchEncoded(sequences[member].enc, rows[i]);
		}
		final MemberSketchIndex index=new MemberSketchIndex(localSketches.width(),
			members, rows, buildThreads);
		final Partition partition=partition(index, maxRounds>0);
		if(partition.left.size<1 || partition.right.size<1){
			return Result.rejected(Status.CONNECTED, partition.components,
				partition.retainedEdges, partition.unionEdges, partition.postingsVisited).withoutCandidates(auditRejected);
		}
		if(partition.left.size<minChildSize || partition.right.size<minChildSize){
			return Result.rejected(Status.UNBALANCED, partition.components,
				partition.retainedEdges, partition.unionEdges, partition.postingsVisited).withoutCandidates(auditRejected);
		}
		int leftSeed=partition.left.array[0], rightSeed=partition.right.array[0];
		if(maxRounds>0){
			leftSeed=center(partition.left, members, partition.weights, index, -1);
			rightSeed=center(partition.right, members, partition.weights, index, leftSeed);
			leftSeed=center(partition.left, members, partition.weights, index, rightSeed);
		}
		final Family left=child(parent.familyId, parent.id, partition.left, members, sequences, leftSeed);
		final Family right=child(rightFamilyId, parent.id+".split"+rightFamilyId,
			partition.right, members, sequences, rightSeed);
		if(maxRounds>0){
			return refineChildren(parent, left, right, sequences, localSketches, gate, refiner,
				metrics, partition, maxRounds, auditRejected).withRepresentatives(members[leftSeed], members[rightSeed]);
		}
		final ConsensusRefiner.Update leftUpdate=refiner.updateAfterConcurrentGreedy(left,
			sequences, gate, metrics, true, summarizeRebuild || auditRejected, auditRejected);
		if(leftUpdate.accepted()+leftUpdate.unchanged()<1){
			final Result result=Result.rejected(Status.LEFT_REJECTED, leftUpdate, null,
				partition.components, partition.retainedEdges,
				partition.unionEdges, partition.postingsVisited);
			return auditRejected ? result.withAudit(audit(left, right, leftUpdate, null,
				sequences, localSketches, gate, refiner)) : result;
		}
		final ConsensusRefiner.Update rightUpdate=refiner.updateAfterConcurrentGreedy(right,
			sequences, gate, metrics, true, summarizeRebuild || auditRejected, auditRejected);
		if(rightUpdate.accepted()+rightUpdate.unchanged()<1){
			final Result result=Result.rejected(Status.RIGHT_REJECTED, leftUpdate, rightUpdate,
				partition.components, partition.retainedEdges,
				partition.unionEdges, partition.postingsVisited);
			return auditRejected ? result.withAudit(audit(left, right, leftUpdate, rightUpdate,
				sequences, localSketches, gate, refiner)) : result;
		}
		left.refreshStats(sequences, localSketches);
		right.refreshStats(sequences, localSketches);
		left.markSketchRebuilt();
		right.markSketchRebuilt();
		return Result.accepted(left, right, leftUpdate, rightUpdate,
			partition.retainedEdges, partition.unionEdges, partition.postingsVisited,
			partition.components);
	}

	/** Builds both initial proposals, then refines actual candidate bytes before any final feasibility decision. */
	private Result refineChildren(final Family parent, final Family left, final Family right,
			final ProteinSequence[] sequences, final FamilySketchBuilder sketches, final TailAlignmentGate gate,
			final ConsensusRefiner refiner, final ClusteringMetrics metrics, final Partition partition,
			final int maxRounds, final boolean auditRejected){
		final ConsensusRefiner.Update initialLeft=refiner.updateAfterConcurrentGreedy(left, sequences, gate, metrics, true, true, true);
		final ConsensusRefiner.Update initialRight=refiner.updateAfterConcurrentGreedy(right, sequences, gate, metrics, true, true, true);
		final Family first=auditCandidate(left, initialLeft, sequences, sketches, gate.minSpan);
		final Family second=auditCandidate(right, initialRight, sequences, sketches, gate.minSpan);
		if(first==null || second==null){
			final Result missing=Result.rejected(first==null ? Status.LEFT_REJECTED : Status.RIGHT_REJECTED,
				initialLeft, initialRight, partition.components, partition.retainedEdges, partition.unionEdges, partition.postingsVisited);
			return auditRejected ? missing.withAudit(new Audit(null, null, initialLeft.rebuildSummary(), initialRight.rebuildSummary(),
				first==null ? "initial_left_candidate_unavailable" : "initial_right_candidate_unavailable", 0, 0, 0, 0)) : missing;
		}
		final CompetitiveRefinement.Result refined=CompetitiveRefinement.refine(first, second,
			parent.memberIndexes(), sequences, sketches, gate, refiner, minChildSize, maxRounds);
		if(refined.strictAccepted()){
			return Result.accepted(refined.first, refined.second, refined.firstUpdate, refined.secondUpdate,
				partition.retainedEdges, partition.unionEdges, partition.postingsVisited, partition.components)
				.withRefinement(refined.summary);
		}
		final Status status=!refined.available() ? Status.REFINEMENT_UNAVAILABLE :
			refined.firstUpdate.accepted()+refined.firstUpdate.unchanged()==0 ? Status.LEFT_REJECTED : Status.RIGHT_REJECTED;
		Result result=Result.rejected(status, refined.firstUpdate, refined.secondUpdate, partition.components,
			partition.retainedEdges, partition.unionEdges, partition.postingsVisited).withRefinement(refined.summary);
		if(auditRejected){
			result=result.withAudit(new Audit(refined.first, refined.second,
				refined.firstUpdate==null ? null : refined.firstUpdate.rebuildSummary(),
				refined.secondUpdate==null ? null : refined.secondUpdate.rebuildSummary(),
				refined.unavailableReason, 0, 0, 0, 0));
		}
		return result;
	}

	/** Completes only a missing sibling of a usable left proposal, without ordinary metrics. */
	private static Audit audit(final Family left, final Family right,
			final ConsensusRefiner.Update leftUpdate, final ConsensusRefiner.Update rightUpdate,
			final ProteinSequence[] sequences, final FamilySketchBuilder sketches,
			final TailAlignmentGate gate, final ConsensusRefiner refiner){
		assert(leftUpdate!=null) : "Fission audit starts only after the ordinary left rebuild returned";
		final Family auditLeft=auditCandidate(left, leftUpdate, sequences, sketches, gate.minSpan);
		if(auditLeft==null){
			return new Audit(null, null, leftUpdate.rebuildSummary(),
				rightUpdate==null ? null : rightUpdate.rebuildSummary(),
				"left_candidate_unavailable", 0, 0, 0, 0);
		}
		ConsensusRefiner.Update completedRight=rightUpdate;
		long rebuilds=0, accepted=0, rejected=0, alignments=0;
		if(completedRight==null){
			try(ClusteringMetrics extra=new ClusteringMetrics(null, null, true, 1, Long.MAX_VALUE)){
				extra.begin("consensus", 0, right.size());
				completedRight=refiner.updateAfterConcurrentGreedy(right, sequences, gate,
					extra, true, true, true);
				extra.end(right.size(), right.size(), 0, 1, 1);
				rebuilds=completedRight.proposals(); accepted=completedRight.accepted();
				rejected=completedRight.rejected(); alignments=completedRight.alignments();
				assert(alignments==extra.constructionAlignments()) :
					"Audit-only sibling work must remain wholly outside ordinary fission metrics";
			}
		}
		final Family auditRight=auditCandidate(right, completedRight, sequences, sketches, gate.minSpan);
		return new Audit(auditLeft, auditRight, leftUpdate.rebuildSummary(),
			completedRight.rebuildSummary(), auditRight==null ? "right_candidate_unavailable" : null,
			rebuilds, accepted, rejected, alignments);
	}

	/** Copies exact current or vetoed consensus bytes and owns a fresh membership list.
	 * Family's copy constructor shares membership, so it cannot provide audit isolation.
	 * A missing rejected proposal never falls back to the child's original seed.
	 */
	private static Family auditCandidate(final Family child, final ConsensusRefiner.Update update,
			final ProteinSequence[] sequences, final FamilySketchBuilder sketches, final int minSpan){
		assert(update!=null && child.size()>0) : "Only rebuilt nonempty child partitions may become audit candidates";
		final byte[] consensus=update.accepted()+update.unchanged()>0 ? child.consensus() : update.rejectedProposal();
		if(consensus==null || consensus.length<minSpan){return null;}
		final Family detached=new Family(child.familyId, child.id, child.memberIndex(0), child.representative());
		for(int i=1; i<child.size(); i++){detached.addMember(child.memberIndex(i));}
		detached.setConsensus(consensus);
		detached.refreshStats(sequences, sketches);
		detached.markSketchRebuilt();
		assert(detached.size()==child.size()) : "Diagnostic replacements must charge every original child member";
		return detached;
	}

	/** Runs the exact DDL graph partitioner without rebuilding child consensuses. */
	public GraphStats analyzeGraph(final MemberSketchIndex index){
		if(index==null || index.size()<1){
			throw new IllegalArgumentException(
				"Fission graph analysis needs a nonempty member index");
		}
		final Partition partition=partition(index);
		return new GraphStats(partition.retainedEdges, partition.unionEdges,
			partition.components, partition.left.size, partition.right.size,
			partition.postingsVisited);
	}

	private void validate(final Family parent, final ProteinSequence[] sequences,
			final FamilySketchBuilder sketches, final TailAlignmentGate gate,
			final ConsensusRefiner refiner, final ClusteringMetrics metrics,
			final int rightFamilyId){
		if(parent==null || parent.size()<1 || parent.sketchUnsafe()==null || sequences==null ||
			sketches==null || gate==null || refiner==null || metrics==null){
			throw new IllegalArgumentException("Fission needs a frozen parent and all rebuild components");
		}
		if(parent.familyId==rightFamilyId || rightFamilyId<0){
			throw new IllegalArgumentException("Right fission child needs a new nonnegative family ID");
		}
		if(refiner.rounds()!=1){
			throw new IllegalArgumentException("Fission rebuilds need a one-round worker refiner");
		}
		if(gate.mode==TailAlignmentGate.Mode.FULL_QUERY_UNPINNED){
			throw new IllegalArgumentException("Fission rebuilds require right-pinned validation");
		}
	}

	/** Splits disconnected components, falling back to a weak-bridge tree cut. */
	private Partition partition(final MemberSketchIndex index){
		return partition(index, false);
	}

	/** Optional central-seed metadata reuses retained graph edges without changing the cut. */
	private Partition partition(final MemberSketchIndex index, final boolean centers){
		final int size=index.size();
		final MemberSketchIndex.Workspace work=new MemberSketchIndex.Workspace();
		final int limit=Math.min(neighborLimit, Math.max(0, size-1));
		final int[] neighbors=new int[limit];
		final int[] hits=new int[neighbors.length];
		final LongList edgeList=new LongList(size);
		long postings=0;
		for(int i=0; i<size; i++){
			final int found=limit<1 ? 0 : index.topNeighbors(i, limit, minHits, neighbors, hits, work);
			postings+=work.postingsVisited();
			for(int j=0; j<found; j++){
				final int a=Math.min(i, neighbors[j]), b=Math.max(i, neighbors[j]);
				edgeList.add((((long)a)<<32) | (b&0xffffffffL));
			}
		}
		if(edgeList.size>0){
			edgeList.sort();
			edgeList.condense();
		}
		final int retainedEdges=edgeList.size;
		final int[] edgeA=new int[retainedEdges], edgeB=new int[retainedEdges];
		final int[] edgeHits=new int[retainedEdges];
		final UnionFind uf=new UnionFind(size);
		int unionEdges=0;
		for(int i=0; i<retainedEdges; i++){
			final long edge=edgeList.array[i];
			final int a=(int)(edge>>>32), b=(int)edge;
			edgeA[i]=a; edgeB[i]=b; edgeHits[i]=index.sharedHits(a, b);
			if(uf.union(a, b)){unionEdges++;}
		}
		final int[] roots=new int[size], counts=new int[size];
		int components=0;
		for(int i=0; i<size; i++){
			final int root=uf.find(i);
			if(counts[root]++==0){roots[components++]=root;}
		}
		if(components<2){
			final Partition cut=balancedCut(index, edgeA, edgeB, edgeHits, retainedEdges);
			if(cut!=null){
				final Partition result=new Partition(cut.left, cut.right, 1, retainedEdges, unionEdges, postings);
				return centers ? weighted(result, size, edgeA, edgeB, edgeHits) : result;
			}
			return new Partition(new IntList(size), new IntList(0), 1, retainedEdges,
				unionEdges, postings);
		}
		sortRootsBySize(roots, counts, components, size);
		final boolean[] leftByRoot=new boolean[size];
		int leftSize=0, rightSize=0;
		for(int i=0; i<components; i++){
			final int root=roots[i];
			if(leftSize<=rightSize){leftByRoot[root]=true; leftSize+=counts[root];}
			else{rightSize+=counts[root];}
		}
		final IntList left=new IntList(leftSize), right=new IntList(rightSize);
		for(int i=0; i<size; i++){
			(leftByRoot[uf.find(i)] ? left : right).add(i);
		}
		final Partition result=new Partition(left, right, components, retainedEdges, unionEdges, postings);
		return centers ? weighted(result, size, edgeA, edgeB, edgeHits) : result;
	}

	/** Sum each retained edge once at each endpoint, excluding edges crossing the selected cut. */
	private static Partition weighted(final Partition partition, final int size,
			final int[] edgeA, final int[] edgeB, final int[] hits){
		final boolean[] right=new boolean[size];
		for(int i=0; i<partition.right.size; i++){right[partition.right.array[i]]=true;}
		final long[] weights=new long[size];
		for(int i=0; i<edgeA.length; i++){
			final int a=edgeA[i], b=edgeB[i];
			if(right[a]==right[b]){weights[a]+=hits[i]; weights[b]+=hits[i];}
		}
		return new Partition(partition.left, partition.right, partition.components, partition.retainedEdges,
			partition.unionEdges, partition.postingsVisited, weights);
	}

	/** Chooses a maximum-weight center, then least shared-hit similarity to the other center and stable ordinal. */
	private static int center(final IntList side, final int[] members, final long[] weights,
			final MemberSketchIndex index, final int other){
		assert(side.size>0 && weights!=null) : "Representative selection needs a nonempty cut and its within-half edge weights";
		int best=-1, bestShared=Integer.MAX_VALUE;
		long bestWeight=-1;
		for(int i=0; i<side.size; i++){
			final int local=side.array[i];
			if(weights[local]<bestWeight){continue;}
			final int shared=other<0 ? 0 : index.sharedHits(local, other);
			if(best<0 || weights[local]>bestWeight || shared<bestShared ||
				shared==bestShared && members[local]<members[best]){
				best=local; bestShared=shared; bestWeight=weights[local];
			}
		}
		return best;
	}

	private Partition balancedCut(final MemberSketchIndex index, final int[] edgeA,
			final int[] edgeB, final int[] edgeHits, final int edges){
		final int size=index.size();
		final LongList order=new LongList(edges);
		for(int i=0; i<edges; i++){
			order.add((((long)(index.width()-edgeHits[i]))<<32) | (i&0xffffffffL));
		}
		order.sort();
		final UnionFind uf=new UnionFind(size);
		final int[] treeA=new int[size-1], treeB=new int[size-1], treeHits=new int[size-1];
		int treeEdges=0, bestHit=0;
		for(int i=0; i<order.size; i++){
			final int edge=(int)order.array[i];
			if(uf.union(edgeA[edge], edgeB[edge])){
				treeA[treeEdges]=edgeA[edge]; treeB[treeEdges]=edgeB[edge];
				treeHits[treeEdges]=edgeHits[edge]; treeEdges++;
				bestHit=Math.max(bestHit, edgeHits[edge]);
			}
		}
		if(treeEdges<size-1){return null;}
		final int[] head=new int[size], next=new int[treeEdges*2], to=new int[treeEdges*2];
		final int[] edgeId=new int[treeEdges*2];
		for(int i=0; i<size; i++){head[i]=-1;}
		for(int i=0, arc=0; i<treeEdges; i++){
			arc=addArc(treeA[i], treeB[i], i, arc, head, next, to, edgeId);
			arc=addArc(treeB[i], treeA[i], i, arc, head, next, to, edgeId);
		}
		final int[] parent=new int[size], parentEdge=new int[size], orderNodes=new int[size];
		for(int i=0; i<size; i++){parent[i]=-2;}
		parent[0]=-1; parentEdge[0]=-1;
		int nodes=0;
		orderNodes[nodes++]=0;
		for(int cursor=0; cursor<nodes; cursor++){
			final int node=orderNodes[cursor];
			for(int arc=head[node]; arc>=0; arc=next[arc]){
				final int child=to[arc];
				if(parent[child]==-2){
					parent[child]=node; parentEdge[child]=edgeId[arc];
					orderNodes[nodes++]=child;
				}
			}
		}
		if(nodes<size){return null;}
		final int[] below=new int[size];
		for(int i=size-1; i>=0; i--){
			final int node=orderNodes[i];
			below[node]++;
			if(parent[node]>=0){below[parent[node]]+=below[node];}
		}
		int bestEdge=-1, bestBalance=-1, bestWeakness=0;
		for(int i=1; i<size; i++){
			final int node=orderNodes[i], cut=parentEdge[node];
			final int childSize=below[node], balance=Math.min(childSize, size-childSize);
			if(balance<minChildSize || treeHits[cut]>=bestHit){continue;}
			final int weakness=bestHit-treeHits[cut];
			if(weakness>bestWeakness || weakness==bestWeakness && balance>bestBalance ||
					(weakness==bestWeakness && balance==bestBalance && cut<bestEdge)){
				bestEdge=cut; bestBalance=balance; bestWeakness=weakness;
			}
		}
		if(bestEdge<0){return null;}
		final int leftStart=treeA[bestEdge];
		final int rightStart=treeB[bestEdge];
		final boolean[] right=markSide(rightStart, leftStart, head, next, to);
		final IntList leftList=new IntList(size-bestBalance), rightList=new IntList(bestBalance);
		for(int i=0; i<size; i++){(right[i] ? rightList : leftList).add(i);}
		return new Partition(leftList, rightList, 1, edges, treeEdges, 0);
	}

	private static int addArc(final int from, final int dest, final int id, final int arc,
			final int[] head, final int[] next, final int[] to, final int[] edgeId){
		to[arc]=dest; edgeId[arc]=id; next[arc]=head[from]; head[from]=arc;
		return arc+1;
	}

	private static boolean[] markSide(final int start, final int blocked,
			final int[] head, final int[] next, final int[] to){
		final boolean[] side=new boolean[head.length];
		final int[] stack=new int[head.length];
		int size=0;
		side[start]=true; stack[size++]=start;
		for(int cursor=0; cursor<size; cursor++){
			final int node=stack[cursor];
			for(int arc=head[node]; arc>=0; arc=next[arc]){
				final int child=to[arc];
				if(child!=blocked && !side[child]){
					side[child]=true; stack[size++]=child;
				}
			}
		}
		return side;
	}

	private static void sortRootsBySize(final int[] roots, final int[] counts,
			final int length, final int size){
		final LongList order=new LongList(length);
		for(int i=0; i<length; i++){
			final int root=roots[i];
			order.add((((long)(size-counts[root]))<<32) | (root&0xffffffffL));
		}
		order.sort();
		for(int i=0; i<length; i++){roots[i]=(int)order.array[i];}
	}

	private static Family child(final int familyId, final String id, final IntList locals,
			final int[] members, final ProteinSequence[] sequences, final int representative){
		final int first=members[representative];
		final Family child=new Family(familyId, id, first, sequences[first].enc);
		for(int i=0; i<locals.size; i++){
			if(locals.array[i]!=representative){child.addMember(members[locals.array[i]]);}
		}
		return child;
	}

	private static final class Partition{
		Partition(final IntList left_, final IntList right_, final int components_,
				final int retainedEdges_, final int unionEdges_, final long postingsVisited_){
			this(left_, right_, components_, retainedEdges_, unionEdges_, postingsVisited_, null);
		}
		Partition(final IntList left_, final IntList right_, final int components_,
				final int retainedEdges_, final int unionEdges_, final long postingsVisited_, final long[] weights_){
			left=left_; right=right_; components=components_;
			retainedEdges=retainedEdges_; unionEdges=unionEdges_; postingsVisited=postingsVisited_;
			weights=weights_;
		}
		final IntList left, right;
		final int components, retainedEdges, unionEdges;
		final long postingsVisited;
		final long[] weights;
	}

	private static final class UnionFind{
		UnionFind(final int size){
			parent=new int[size]; rank=new byte[size];
			for(int i=0; i<size; i++){parent[i]=i;}
		}
		int find(final int x){
			int root=x;
			while(parent[root]!=root){root=parent[root];}
			for(int y=x; parent[y]!=y; ){
				final int next=parent[y]; parent[y]=root; y=next;
			}
			return root;
		}
		boolean union(final int a, final int b){
			int ra=find(a), rb=find(b);
			if(ra==rb){return false;}
			if(rank[ra]<rank[rb]){final int temp=ra; ra=rb; rb=temp;}
			parent[rb]=ra;
			if(rank[ra]==rank[rb]){rank[ra]++;}
			return true;
		}
		final int[] parent;
		final byte[] rank;
	}

	public enum Status{
		ACCEPTED, TOO_SMALL, CONNECTED, UNBALANCED, LEFT_REJECTED, RIGHT_REJECTED, REFINEMENT_UNAVAILABLE
	}

	public static final class GraphStats{
		GraphStats(final int retainedEdges_, final int unionEdges_,
				final int components_, final int leftSize_, final int rightSize_,
				final long postingsVisited_){
			retainedEdges=retainedEdges_; unionEdges=unionEdges_;
			components=components_; leftSize=leftSize_; rightSize=rightSize_;
			postingsVisited=postingsVisited_;
		}
		public final int retainedEdges, unionEdges, components, leftSize, rightSize;
		public final long postingsVisited;
	}

	public static final class Result{
		/** Retains graph/status data while recording that no pair was structurally proposed. */
		private Result withoutCandidates(final boolean auditRejected){
			return auditRejected ? withAudit(new Audit(null, null, null, null,
				status.name(), 0, 0, 0, 0)) : this;
		}
		/** Adds request-local diagnostics without making a rejected split accepted. */
		private Result withAudit(final Audit diagnostic){
			assert(!accepted() && diagnostic!=null) : "Accepted fissions reuse ordinary children rather than extra audit models";
			return new Result(status, left, right, leftUpdate, rightUpdate, retainedEdges,
				unionEdges, postingsVisited, components, diagnostic, refinement, leftRepresentative, rightRepresentative);
		}
		private Result withRefinement(final CompetitiveRefinement.Summary summary){
			return new Result(status, left, right, leftUpdate, rightUpdate, retainedEdges,
				unionEdges, postingsVisited, components, audit, summary, leftRepresentative, rightRepresentative);
		}
		private Result withRepresentatives(final int first, final int second){
			return new Result(status, left, right, leftUpdate, rightUpdate, retainedEdges,
				unionEdges, postingsVisited, components, audit, refinement, first, second);
		}
		private static Result accepted(final Family left_, final Family right_,
				final ConsensusRefiner.Update leftUpdate_,
				final ConsensusRefiner.Update rightUpdate_,
				final int retainedEdges_, final int unionEdges_,
				final long postingsVisited_, final int components_){
			return new Result(Status.ACCEPTED, left_, right_, leftUpdate_, rightUpdate_,
				retainedEdges_, unionEdges_, postingsVisited_, components_);
		}
		private static Result rejected(final Status status_, final int components_,
				final int retainedEdges_, final int unionEdges_, final long postingsVisited_){
			return rejected(status_, null, null, components_, retainedEdges_, unionEdges_, postingsVisited_);
		}
		private static Result rejected(final Status status_,
				final ConsensusRefiner.Update leftUpdate_,
				final ConsensusRefiner.Update rightUpdate_, final int components_,
				final int retainedEdges_, final int unionEdges_, final long postingsVisited_){
			if(status_==Status.ACCEPTED){
				throw new IllegalArgumentException("Accepted results need children");
			}
			return new Result(status_, null, null, leftUpdate_, rightUpdate_, retainedEdges_,
				unionEdges_, postingsVisited_, components_);
		}
		private Result(final Status status_, final Family left_, final Family right_,
				final ConsensusRefiner.Update leftUpdate_,
				final ConsensusRefiner.Update rightUpdate_,
				final int retainedEdges_, final int unionEdges_,
				final long postingsVisited_, final int components_){
			this(status_, left_, right_, leftUpdate_, rightUpdate_, retainedEdges_,
				unionEdges_, postingsVisited_, components_, null, null, -1, -1);
		}
		private Result(final Status status_, final Family left_, final Family right_,
				final ConsensusRefiner.Update leftUpdate_,
				final ConsensusRefiner.Update rightUpdate_,
				final int retainedEdges_, final int unionEdges_,
				final long postingsVisited_, final int components_, final Audit audit_,
				final CompetitiveRefinement.Summary refinement_, final int firstRepresentative, final int secondRepresentative){
			status=status_; left=left_; right=right_;
			leftUpdate=leftUpdate_; rightUpdate=rightUpdate_;
			retainedEdges=retainedEdges_; unionEdges=unionEdges_;
			postingsVisited=postingsVisited_; components=components_;
			audit=audit_;
			refinement=refinement_; leftRepresentative=firstRepresentative; rightRepresentative=secondRepresentative;
		}
		public boolean accepted(){return status==Status.ACCEPTED;}
		public final Status status;
		public final Family left, right;
		public final ConsensusRefiner.Update leftUpdate, rightUpdate;
		public final int retainedEdges, unionEdges, components;
		public final long postingsVisited;
		public final Audit audit;
		public final CompetitiveRefinement.Summary refinement;
		/** Original ordinals of the initial central representatives; -1 when central seeding was not requested. */
		public final int leftRepresentative, rightRepresentative;
	}

	/** Request-local rejected candidates; never retain these models in batch summaries.
	 * Side summaries include an audit-only right rebuild, whereas Result.rightUpdate
	 * remains null after LEFT_REJECTED. Extra counts keep native Update meanings:
	 * unchanged rebuilds contribute to extraRebuilds but not extraAccepted.
	 */
	public static final class Audit{
		private Audit(final Family left_, final Family right_,
				final ConsensusRefiner.RebuildSummary leftRebuild_,
				final ConsensusRefiner.RebuildSummary rightRebuild_, final String unavailableReason_,
				final long extraRebuilds_, final long extraAccepted_, final long extraRejected_,
				final long extraAlignments_){
			if((unavailableReason_==null)!=(left_!=null && right_!=null) ||
				extraRebuilds_<0 || extraAccepted_<0 || extraRejected_<0 || extraAlignments_<0 ||
				extraAccepted_>extraRebuilds_ || extraRejected_>extraRebuilds_-extraAccepted_){
				throw new IllegalArgumentException("Fission audit availability and native extra-work counters must agree");
			}
			left=left_; right=right_; leftRebuild=leftRebuild_; rightRebuild=rightRebuild_;
			unavailableReason=unavailableReason_; extraRebuilds=extraRebuilds_;
			extraAccepted=extraAccepted_; extraRejected=extraRejected_; extraAlignments=extraAlignments_;
		}
		public boolean available(){return unavailableReason==null;}
		public final Family left, right;
		public final ConsensusRefiner.RebuildSummary leftRebuild, rightRebuild;
		public final String unavailableReason;
		public final long extraRebuilds, extraAccepted, extraRejected, extraAlignments;
	}

	private final int minChildSize, neighborLimit, minHits, buildThreads;
}
