package prok;

import java.util.Arrays;
import structures.IntList;

/** Invocation-local observation state allocated only for enabled window tracing.
 * It borrows the seed stream read-only; every emitted group owns fresh member arrays.
 * Coordinates are inclusive and oriented. Unknown integral values use Long.MIN_VALUE;
 * negative raw coordinates are valid predictions, not missing-value sentinels.
 * @author Keqing
 */
final class NcrnaVoteDiagnostics {

	NcrnaVoteDiagnostics(NcrnaVoteDiagSink sink_, String family_, String contig_, int strand_,
			int length_, int pass_, int slack_, int pad_, int[] centers_, long[] keys_, SeedOffsetTable table_){
		if(sink_==null || centers_==null || (keys_!=null && keys_.length!=centers_.length)
			|| (strand_!=0 && strand_!=1) || (pass_!=1 && pass_!=2)){
			throw new IllegalArgumentException("Vote provenance needs a sink, aligned seed streams, strand and actual scavenger pass");
		}
		sink=sink_;family=family_;contig=contig_;strand=strand_;length=length_;pass=pass_;
		slack=slack_;pad=pad_;centers=centers_;keys=keys_;table=table_;
		parents=new int[centers.length];Arrays.fill(parents, -1);resetVote();
	}

	/** Clears only observations; never clears or feeds the caller's voting scratch. */
	void resetVote(){
		voters=0;rawLeft=rawRight=leftSum=rightSum=leftSd=rightSd=Double.NaN;
		leftMin=rightMin=Double.POSITIVE_INFINITY;leftMax=rightMax=Double.NEGATIVE_INFINITY;
	}
	boolean hasWeightedVote(){return !Double.isNaN(rawLeft) && !Double.isNaN(rawRight);}
	void contribution(double left, double right){
		assert(Double.isFinite(left) && Double.isFinite(right)) : "A contributed seed prediction must be finite before it can explain a weighted vote";
		voters++;leftMin=Math.min(leftMin, left);leftMax=Math.max(leftMax, left);
		rightMin=Math.min(rightMin, right);rightMax=Math.max(rightMax, right);
	}
	void weighted(double left, double right, double wl, double wr, double sdLeft, double sdRight){
		assert(voters>0 && wl>0 && wr>0) : "Weighted diagnostic values must come from actual positive-weight seed contributions";
		rawLeft=left;rawRight=right;leftSum=wl;rightSum=wr;leftSd=sdLeft;rightSd=sdRight;
	}

	/** Records exactly the legacy pre-vote interval membership, preserving input order. */
	void legacy(int from, int to, int start, int stop, String disposition, String reason){
		int count=0;for(int center:centers){if(center>=from && center<=to){count++;}}
		final int[] members=new int[count];int j=0;
		for(int i=0; i<centers.length; i++){if(centers[i]>=from && centers[i]<=to){members[j++]=i;}}
		emit("legacy", disposition, reason, from, to, start, stop, members, false);
	}

	/** Uses the builder's actual sorted cluster and already-computed means/sums. */
	void core(int from, int to, long[] order, int[] original, double[] left, double[] right,
			double l, double r, double wl, double wr, long start, long stop, boolean retained){
		assert(to>from) : "Each core proposal must correspond to an executed nonempty builder group";
		resetVote();final int[] members=new int[to-from];
		for(int j=from; j<to; j++){
			final int i=(int)order[j];members[j-from]=original[i];contribution(left[i], right[i]);
		}
		weighted(l, r, wl, wr, Double.NaN, Double.NaN);
		emit("core", retained ? "retained_vote" : "fallback", retained ? "VOTE_RETAINED" : "SHORT_VOTE",
			UNKNOWN, UNKNOWN, start, stop, members, !retained);
	}

	/** The fallback subrange is selected by the builder, never a proximity search. */
	void fallback(IntList fallback, int from, int to, int start, int stop, boolean retained){
		assert(to>from) : "Each fallback proposal owns a nonempty range of the builder's sorted fallback stream";
		resetVote();final int[] members=new int[to-from];boolean untrained=false, shortVote=false;
		for(int j=from; j<to; j++){
			final int i=Arrays.binarySearch(centers, fallback.get(j));
			if(i<0){throw new IllegalStateException("Fallback seed is missing from its original invocation: "+fallback.get(j));}
			members[j-from]=i;if(parents[i]>=0){shortVote=true;}else{untrained=true;}
		}
		final String reason=!retained ? "SHORT_FALLBACK" : untrained && shortVote ? "MIXED_FALLBACK"
			: shortVote ? "SHORT_VOTE_SEEDS" : "UNTRAINED_SEEDS";
		emit("fallback", retained ? "retained_fallback" : "dropped", reason, UNKNOWN, UNKNOWN, start, stop, members, false);
	}

	private void emit(String policy, String disposition, String reason, long inputStart, long inputEnd,
			long start, long stop, int[] members, boolean assignParent){
		final int id=nextId++, n=members.length;
		final int[] c=new int[n], trained=new int[n], parent=new int[n];final long[] k=new long[n];
		for(int j=0; j<n; j++){
			final int i=members[j];
			assert(i>=0 && i<centers.length) : "Member indices must come from this invocation's actual seed stream";
			c[j]=centers[i];k[j]=keys==null ? UNKNOWN : keys[i];parent[j]=parents[i];
			final boolean known=table!=null && keys!=null && keys[i]!=UNKNOWN;
			final SeedOffsetTable.KmerInfo info=known ? table.get(keys[i]) : null;
			trained[j]=!known ? -1 : info!=null && info.trained() ? 1 : 0;
			if(assignParent){parents[i]=id;}
		}
		final Group group=new Group(family, contig, strand, length, pass, id, policy, disposition, reason,
			inputStart, inputEnd, rawLeft, rawRight, Double.isNaN(rawLeft) ? UNKNOWN : Math.round(rawLeft),
			Double.isNaN(rawRight) ? UNKNOWN : Math.round(rawRight), start, stop, slack, pad, voters,
			c, k, trained, parent, leftSum, rightSum, voters==0 ? Double.NaN : leftMax-leftMin,
			voters==0 ? Double.NaN : rightMax-rightMin, leftSd, rightSd);
		sink.group(group);
	}

	/** Immutable scalar evidence with observer-owned mutable arrays, detached from the caller. */
	static final class Group {
		Group(String family_, String contig_, int strand_, int length_, int pass_, int id_, String policy_,
				String disposition_, String reason_, long inputStart_, long inputEnd_, double rawStart_, double rawEnd_,
				long roundedStart_, long roundedEnd_, long start_, long stop_, int slack_, int pad_, int voters_,
				int[] centers_, long[] keys_, int[] trained_, int[] parents_, double leftWeight_, double rightWeight_,
				double leftRange_, double rightRange_, double leftSd_, double rightSd_){
			assert(centers_.length==keys_.length && centers_.length==trained_.length && centers_.length==parents_.length)
				: "Aligned diagnostic lists preserve exact seed identity and ancestry for every member";
			family=family_;contig=contig_;strand=strand_;length=length_;pass=pass_;id=id_;policy=policy_;
			disposition=disposition_;reason=reason_;inputStart=inputStart_;inputEnd=inputEnd_;rawStart=rawStart_;rawEnd=rawEnd_;
			roundedStart=roundedStart_;roundedEnd=roundedEnd_;start=start_;stop=stop_;slack=slack_;pad=pad_;voters=voters_;
			centers=centers_;keys=keys_;trained=trained_;parents=parents_;leftWeight=leftWeight_;rightWeight=rightWeight_;
			leftRange=leftRange_;rightRange=rightRange_;leftSd=leftSd_;rightSd=rightSd_;
		}
		final String family, contig, policy, disposition, reason;
		final int strand, length, pass, id, slack, pad, voters;
		final long inputStart, inputEnd, roundedStart, roundedEnd, start, stop;
		final double rawStart, rawEnd, leftWeight, rightWeight, leftRange, rightRange, leftSd, rightSd;
		final int[] centers, trained, parents;
		final long[] keys;
	}

	static final long UNKNOWN=Long.MIN_VALUE;
	private final NcrnaVoteDiagSink sink;
	private final String family, contig;
	private final int strand, length, pass, slack, pad;
	private final int[] centers, parents;
	private final long[] keys;
	private final SeedOffsetTable table;
	private int nextId, voters;
	private double rawLeft, rawRight, leftSum, rightSum, leftSd, rightSd, leftMin, leftMax, rightMin, rightMax;
}
