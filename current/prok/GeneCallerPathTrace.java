package prok;

import java.util.ArrayList;

/** Contig-local observation state for GeneCaller's actual bounded path search.
 * Only primitive snapshots reach the observer; candidate objects remain caller-owned.
 * This records decisions, not an exhaustive optimum or a counterfactual bonus threshold.
 * @author Brian, Keqing
 */
public final class GeneCallerPathTrace{

	/** Callbacks run synchronously on the caller thread. Do not reenter the caller
	 * or change its scoring policy. IDs pack frame-list number and index; -1 is absent.
	 * candidate selects which nodes' incoming/outgoing comparisons need edge/window
	 * detail. scored still reports every node so the selected predecessor chain can
	 * be reconstructed in memory. A false selection never changes the search.
	 */
	public interface Observer{
		void begin(String contig, int contigLength, int candidates);
		boolean candidate(long id, int start, int stop, int strand, int frame, int type,
			int stopLeft, float score, float bonus, boolean nativeMeasured, float nativeScore);
		void edge(long node, long predecessor, int predecessorStrand, int overlap, float predecessorScore,
			int predecessorLength, float overlapScore, float beforeTransitions, float proposedScore,
			int proposedLength, boolean replaces);
		void window(long node, int predecessorList, int lastPosition, int bestPosition, int nominalMinimum, int lastVisited);
		void scored(long node, float plusScore, float minusScore, int plusLength, int minusLength,
			long plusPrev, long minusPrev);
		void end(long bestNode, float bestScore, int selectedNodes);
	}

	GeneCallerPathTrace(Observer observer_, ArrayList<Orf>[] frames_, int length, boolean cascade){
		require(observer_!=null && frames_!=null && frames_.length==6 && length>=0,
			"Path observations require six retained frame lists and a real contig length");
		observer=observer_; frames=frames_; int count=0; String name=null;
		for(int f=0; f<6; f++){
			require(frames[f]!=null, "Retained path lists must be explicit even when empty");
			count=Math.addExact(count, frames[f].size()); wanted[f]=new boolean[frames[f].size()];
			Orf previous=null;
			for(Orf orf:frames[f]){
				require(orf!=null && orf.prevPlus==null && orf.prevMinus==null,
					"Observe a fresh candidate graph; stale predecessor links cannot describe this path search");
				if(name==null){name=orf.scafName;}
				require(name!=null && name.equals(orf.scafName) && orf.scaflen==length &&
					orf.strand>=0 && orf.strand<2 && orf.frame>=0 && orf.frame<3 && 3*orf.strand+orf.frame==f,
					"Path node identity requires one contig and the caller's original frame-list membership");
				require(previous==null || previous.compareTo(orf)<=0,
					"Path search advances each original frame list in ascending endpoint order");
				previous=orf;
			}
		}
		observer.begin(name, length, count);// A completely empty graph has no candidate-derived name.
		for(int f=0; f<6; f++){
			for(int i=0; i<frames[f].size(); i++){
				final Orf orf=frames[f].get(i);
				final boolean nativeMeasured=cascade && orf.type==ProkObject.CDS;
				final int stopLeft=orf.type==ProkObject.CDS && !orf.stopTruncated ? GeneCaller.tailStopLeft(orf) : -1;
				require(Float.isFinite(orf.orfScore) && Float.isFinite(orf.proteinBonus) &&
					(!nativeMeasured || Float.isFinite(orf.tailCascadeNativeScore)),
					"Path observations cannot publish nonfinite candidate scores");
				wanted[f][i]=observer.candidate(id(f, i), orf.start, orf.stop, orf.strand, orf.frame, orf.type,
					stopLeft, orf.orfScore, orf.proteinBonus, nativeMeasured, nativeMeasured ? orf.tailCascadeNativeScore : 0);
				wantedFrame[f]|=wanted[f][i];
			}
		}
	}

	void start(Orf orf, int index){
		currentList=3*orf.strand+orf.frame; currentIndex=index;
		require(index>=0 && index<frames[currentList].size() && frames[currentList].get(index)==orf,
			"The observed current node must be the next actual scored entry in its frame list");
		plusPrev=minusPrev=-1;
	}

	void edge(int list, int index, int overlap, float predecessorScore, int predecessorLength,
			float overlapScore, float proposedScore, int proposedLength, boolean replaces){
		final long previous=id(list, index); final int strand=list/3;
		if(replaces){if(strand==0){plusPrev=previous;}else{minusPrev=previous;}}
		if(!wanted[currentList][currentIndex] && !wanted[list][index]){return;}
		require(overlap>=0 && Float.isFinite(predecessorScore) && Float.isFinite(overlapScore) &&
			Float.isFinite(predecessorScore+overlapScore) && Float.isFinite(proposedScore),
			"An observed evaluated edge must retain finite actual path scores and a nonnegative overlap");
		observer.edge(id(currentList, currentIndex), previous, strand, overlap, predecessorScore, predecessorLength,
			overlapScore, predecessorScore+overlapScore, proposedScore, proposedLength, replaces);
	}

	void window(int list, int last, int best, int minimum, int visited){
		if(wanted[currentList][currentIndex] || wantedFrame[list]){
			observer.window(id(currentList, currentIndex), list, last, best, minimum, visited);
		}
	}

	void scored(Orf orf){
		require(frames[currentList].get(currentIndex)==orf && node(plusPrev)==orf.prevPlus && node(minusPrev)==orf.prevMinus,
			"Recorded predecessor IDs must equal the caller's final predecessor objects for this node");
		require(Float.isFinite(orf.pathScorePlus) && Float.isFinite(orf.pathScoreMinus),
			"Observed predecessor-state scores must be finite");
		observer.scored(id(currentList, currentIndex), orf.pathScorePlus, orf.pathScoreMinus,
			orf.pathLengthPlus, orf.pathLengthMinus, plusPrev, minusPrev);
	}

	void finish(Orf best, int selected){
		if(best==null){require(selected==0, "An empty path has no selected nodes"); observer.end(-1, 0, 0); return;}
		for(int f=0; f<6; f++){
			for(int i=0; i<frames[f].size(); i++){
				if(frames[f].get(i)==best){observer.end(id(f, i), best.pathScore(), selected); return;}
			}
		}
		throw new IllegalStateException("Selected path endpoint is absent from the retained graph");
	}

	private Orf node(long id){return id<0 ? null : frames[(int)(id>>>32)].get((int)id);}
	private static long id(int list, int index){return ((long)list<<32)|(index&0xffffffffL);}
	private static void require(boolean valid, String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final Observer observer;
	private final ArrayList<Orf>[] frames;
	private final boolean[][] wanted=new boolean[6][];
	private final boolean[] wantedFrame=new boolean[6];
	private int currentList, currentIndex;
	private long plusPrev=-1, minusPrev=-1;
}
