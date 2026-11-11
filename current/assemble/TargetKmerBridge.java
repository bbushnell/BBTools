package assemble;

import java.util.Arrays;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import ukmer.Kmer;

/**
 * Finds a bounded, unbranched kmer path between inward-facing paired-read tips.
 * Both input reads must already be oriented in the same genomic direction.
 * Each instance owns its scratch space and is intended for one worker thread.
 *
 * @author Brian, Noelle
 */
public final class TargetKmerBridge{

	/**
	 * Allocates worker-local scratch for one Tadpole count graph.
	 * @param tadpole_ graph owner used for kmer counts, extension counts and thresholds
	 * @param maxDistance_ maximum number of extension bases to search; must be in [0, MAX_BRIDGE_DISTANCE]
	 * @param reciprocal_ whether to require the reverse walk to reproduce the forward insert
	 */
	public TargetKmerBridge(final Tadpole tadpole_, final int maxDistance_, final boolean reciprocal_){
		if(maxDistance_<0 || maxDistance_>MAX_BRIDGE_DISTANCE){
			throw new IllegalArgumentException("maxDistance must be in [0, "+MAX_BRIDGE_DISTANCE+"], not "+maxDistance_);
		}
		tadpole=tadpole_;
		k=tadpole.k();
		maxDistance=maxDistance_;
		reciprocal=reciprocal_;
		forwardSource=new Kmer(k);
		forwardTarget=new Kmer(k);
		reverseSource=new Kmer(k);
		reverseTarget=new Kmer(k);
		walkKmer=new Kmer(k);
		int visitedCapacity=2;
		while(visitedCapacity<(maxDistance+2)*2){visitedCapacity<<=1;}
		visitedHashes=KillSwitch.allocLong1D(visitedCapacity);
		visitedEpochs=KillSwitch.allocInt1D(visitedCapacity);
		visitedMask=visitedCapacity-1;
	}

	/** Returns inferred insert size, or -1.  Does not mutate either read. */
	public int find(final Read left, final Read right){
		status=NO_PATH;
		forwardDistance=forwardHits=0;
		if(left==null || right==null || left.length()<k || right.length()<k){return -1;}
		if(!loadKmer(forwardSource, left.bases, left.length()-k)){return -1;}
		if(!loadKmer(forwardTarget, right.bases, 0)){status=NO_TARGET; return -1;}

		//Both endpoints are cheap to reject and must be viable before either walk begins.
		if(tadpole.bridgeCount(forwardSource)<tadpole.minCountSeed){return -1;}
		if(tadpole.bridgeCount(forwardTarget)<tadpole.minCountSeed){status=NO_TARGET; return -1;}

		forwardExtension.clear();
		final int forward=trace(forwardSource, forwardTarget, left.length(), right.length(), forwardExtension);
		if(forward<0){return -1;}
		forwardDistance=traceDistance;
		forwardHits=1;

		if(reciprocal){
			reverseSource.setFrom(forwardTarget).rcomp();
			reverseTarget.setFrom(forwardSource).rcomp();
			final int reverse=trace(reverseSource, reverseTarget, right.length(), left.length(), null);
			if(reverse<0 || reverse!=forward){
				status=NONRECIPROCAL;
				return -1;
			}
		}
		status=SUCCESS;
		return forward;
	}

	/** Applies the successful forward graph path to the left read. */
	public int apply(final Read left){
		assert(status==SUCCESS);
		final int oldLength=left.length();
		final int newLength=oldLength+forwardDistance;
		assert(forwardDistance<=forwardExtension.length()) : forwardDistance+", "+forwardExtension.length();
		left.bases=KillSwitch.copyOf(left.bases, newLength);
		System.arraycopy(forwardExtension.array, 0, left.bases, oldLength, forwardDistance);
		if(left.quality!=null){
			left.quality=KillSwitch.copyOf(left.quality, newLength);
			for(int i=oldLength; i<newLength; i++){left.quality[i]=Shared.FAKE_QUAL;}
		}
		return forwardDistance;
	}

	/** Walks the strongest rightward path from {@code source} toward {@code target}. */
	private int trace(final Kmer source, final Kmer target, final int seedLength,
			final int targetLength, final ByteBuilder extension){
		traceDistance=0;
		walkKmer.setFrom(source);
		final long targetXor=target.xor();
		clearVisited();
		addVisited(walkKmer.xor());

		for(int distance=1; distance<=maxDistance; distance++){
			final int maxPos=tadpole.bridgeRightCounts(walkKmer, rightCounts);
			final int max=rightCounts[maxPos];
			final int second=rightCounts[Tools.secondHighestPosition(rightCounts)];
			if(max<tadpole.minCountExtend || tadpole.isJunction(max, second)){status=NO_PATH; return -1;}

			walkKmer.addRightNumeric(maxPos);
			if(extension!=null){extension.append(AminoAcid.numberToBase[maxPos]);}
			final long xor=walkKmer.xor();
			if(!addVisited(xor)){status=CYCLE; return -1;}
			if(xor==targetXor && walkKmer.sameOrientation(target)){
				traceDistance=distance;
				return seedLength+targetLength-k+distance;
			}
		}
		status=NO_TARGET;
		return -1;
	}

	/** Loads a concrete base window into {@code kmer}; ambiguous bases abort the candidate. */
	private boolean loadKmer(final Kmer kmer, final byte[] bases, final int start){
		kmer.clear();
		for(int i=start, lim=start+k; i<lim; i++){
			final int x=AminoAcid.baseToNumber[bases[i]];
			if(x<0){return false;}
			kmer.addRightNumeric(x);
		}
		return kmer.len()>=k;
	}

	/** Advances the visited-set epoch, clearing the epoch table on integer wraparound. */
	private void clearVisited(){
		visitedEpoch++;
		if(visitedEpoch==0){
			Arrays.fill(visitedEpochs, 0);
			visitedEpoch=1;
		}
	}

	/**
	 * Adds one kmer xor to the trace-local visited set.
	 * @return false on a repeated xor hash, including any 64-bit collision
	 */
	private boolean addVisited(final long hash){
		//A 64-bit collision can only reject a valid path as cyclic; it cannot create a false bridge.
		int cell=((int)hash)&visitedMask;
		while(visitedEpochs[cell]==visitedEpoch){
			if(visitedHashes[cell]==hash){return false;}
			cell=(cell+1)&visitedMask;
		}
		visitedEpochs[cell]=visitedEpoch;
		visitedHashes[cell]=hash;
		return true;
	}

	/** Returns the terminal status from the most recent {@link #find(Read, Read)} call. */
	public int status(){return status;}
	/** Returns the most recent forward hit distance, even if reciprocal checking later failed. */
	public int distance(){return forwardDistance;}
	/** Returns 1 when the latest search reached its forward target, otherwise 0. */
	public int hits(){return forwardHits;}

	/** A bounded, unbranched path reached the target; {@link #apply(Read)} is valid. */
	public static final int SUCCESS=0;
	/** The source is unusable or the forward walk reaches a dead end or junction. */
	public static final int NO_PATH=1;
	/** The target kmer is absent, ambiguous, or outside the bounded walk. */
	public static final int NO_TARGET=2;
	/** Reserved status for ambiguous bridge candidates. */
	public static final int AMBIGUOUS=3;
	/** The walk repeated a kmer xor hash before reaching the target. */
	public static final int CYCLE=4;
	/** The reverse path failed or produced a different insert size. */
	public static final int NONRECIPROCAL=5;
	/** Largest distance whose visited-table capacity can grow to a positive power of two. */
	public static final int MAX_BRIDGE_DISTANCE=(1<<29)-2;

	/** Count graph and threshold owner. */
	private final Tadpole tadpole;
	/** Active kmer length and maximum extension distance in bases. */
	private final int k, maxDistance;
	/** Whether a successful forward bridge needs reverse confirmation. */
	private final boolean reciprocal;
	/** Reused right-extension counts for the current walk kmer. */
	private final int[] rightCounts=new int[4];
	/** Reused Kmer workspaces for endpoints and the active walk. */
	private final Kmer forwardSource, forwardTarget, reverseSource, reverseTarget, walkKmer;
	/** Open-addressed xor table for cycle rejection within one trace. */
	private final long[] visitedHashes;
	/** Epoch tags parallel to {@link #visitedHashes}. */
	private final int[] visitedEpochs;
	/** Power-of-two mask for the visited table. */
	private final int visitedMask;
	/** Bases that should be appended to the left read after a successful forward walk. */
	private final ByteBuilder forwardExtension=new ByteBuilder();
	/** Terminal status from the latest candidate. */
	private int status=NO_PATH;
	/** Extension bases walked in the active trace. */
	private int traceDistance=0;
	/** Successful forward extension distance and hit count. */
	private int forwardDistance=0, forwardHits=0;
	/** Visited-set epoch; zero at construction and nonzero while tracing paths. */
	private int visitedEpoch=0;
}
