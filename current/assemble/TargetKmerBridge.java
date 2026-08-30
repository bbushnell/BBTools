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
public final class TargetKmerBridge {

	public TargetKmerBridge(final Tadpole tadpole_, final int maxDistance_, final boolean reciprocal_){
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

	private boolean loadKmer(final Kmer kmer, final byte[] bases, final int start){
		kmer.clear();
		for(int i=start, lim=start+k; i<lim; i++){
			final int x=AminoAcid.baseToNumber[bases[i]];
			if(x<0){return false;}
			kmer.addRightNumeric(x);
		}
		return kmer.len()>=k;
	}

	private void clearVisited(){
		visitedEpoch++;
		if(visitedEpoch==0){
			Arrays.fill(visitedEpochs, 0);
			visitedEpoch=1;
		}
	}

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

	public int status(){return status;}
	public int distance(){return forwardDistance;}
	public int hits(){return forwardHits;}

	public static final int SUCCESS=0;
	public static final int NO_PATH=1;
	public static final int NO_TARGET=2;
	public static final int AMBIGUOUS=3;
	public static final int CYCLE=4;
	public static final int NONRECIPROCAL=5;

	private final Tadpole tadpole;
	private final int k, maxDistance;
	private final boolean reciprocal;
	private final int[] rightCounts=new int[4];
	private final Kmer forwardSource, forwardTarget, reverseSource, reverseTarget, walkKmer;
	private final long[] visitedHashes;
	private final int[] visitedEpochs;
	private final int visitedMask;
	private final ByteBuilder forwardExtension=new ByteBuilder();
	private int status=NO_PATH;
	private int traceDistance=0;
	private int forwardDistance=0, forwardHits=0;
	private int visitedEpoch=0;
}
