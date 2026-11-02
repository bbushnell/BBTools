package repeat;

import java.util.ArrayList;
import java.util.Collections;

import repeat.Repeat.PosComparator2;
import shared.Tools;
import stream.Read;
import structures.CRange;
import tracker.EntropyTracker;

/** Worker-confined depth slots, closed repeat snapshots and reusable range/statistics buffers.
 * Lists and Repeat objects are mutable; returned buffers are borrowed unless a method
 * explicitly creates a new list. Contigs remain aliased through retained repeats.
 */
public class RepeatSet{

	/** Creates empty tracking state and a lazy entropy cache for the configured maximum window. */
	RepeatSet(int k_, int minDepth_, int maxDepth_, int minRepeat_, int maxGap_, boolean weakSubsumes_, boolean amino_, int ek_, int ew_){
		k=k_;
		minDepth=minDepth_;
		maxDepth=maxDepth_;
		minRepeat=minRepeat_;
		maxGap=maxGap_;
		weakSubsumes=weakSubsumes_;
		amino=amino_;
		ek=ek_;
		ew=ew_;
		eta=new EntropyTracker[ew+1];
	}

	/** Returns a cached mutable tracker, or null when the region contains at most one full kmer.
	 * That degenerate region has zero entropy and must not construct a window equal
	 * to k. Otherwise the window is the median of k, maximum window and length.
	 * The caller resets the tracker before use; the cache is not safe for concurrent access. */
	EntropyTracker getET(int len){
		//REPEAT-005: retain the observed zero result without an invalid tracker or a changed normalization window.
		if(len<=ek){return null;}
		int window=Tools.mid(ek, ew, len);
		if(eta[window]==null){eta[window]=new EntropyTracker(ek, window, amino, 0, true);}
		return eta[window];
	}

	//Replaced by collectResidual
//	void flushOpen() {
//		for(Repeat r : openRepeats) {
//			if(r.contigNum>=0 && r.depth>1 && r.length()>=minRepeat){
//				assert(false);
//				oldRepeats.add(r.clone());
//			}
//		}
//	}

	/** Moves closed references into retained history and empties the closed list, without cloning. */
	void retireClosed(){
		oldRepeats.addAll(closedRepeats);
		closedRepeats.clear();
	}

	/** Retains r unless the last accepted repeat subsumes it; calculates statistics before retention.
	 * The same reference enters closed/recent lists, so callers must stop mutating it. */
	void addRepeat(Repeat r){
		if(lastRepeat!=null && lastRepeat.subsumes(r, weakSubsumes)){
			//for more efficiency I could prevent r from being created in the first place
			//but that would require adding complexity to Repeat
			return;
		}
		//Actually this can happen in rare cases due to lazy collection.
//		assert(lastRepeat==null || !r.subsumes(lastRepeat)) : "\n"+lastRepeat+"\n"+r;

		r.calcStats(getET(r.length()));
		closedRepeats.add(r);
		recent.add(r);
		lastRepeat=r;
	}

	/** Closes qualifying residual depths through maxDepthSeen, cloning and clearing those slots.
	 * Short active slots are left for increment's next region/contig transition. */
	void collectResidual(int maxDepthSeen){
		maxDepthSeen=Tools.min(maxDepthSeen, openRepeats.size()-1); //In case maxDepth is set.
		assert(maxDepthSeen<openRepeats.size()) : maxDepthSeen+", "+openRepeats.size();
		for(int i=maxDepthSeen; i>=minDepth; i--){
			Repeat r=openRepeats.get(i);
			if(r.contigNum>=0 && r.depth>1 && r.length()>=minRepeat){
				addRepeat(r.clone());
				r.clear();
			}
		}
	}

	/**
	 * Increments repeat tracking at the specified position and depth.
	 * Expands open repeats array as needed and processes all depth levels
	 * from min(depth,maxDepth) down to minDepth. Each slot receives the original
	 * observed depth; the cap limits tracked thresholds, not the observed argument.
	 *
	 * @param contig The contig being processed
	 * @param pos Position in the contig (must be >= k-1)
	 * @param depth Coverage depth at this position
	 */
	void increment(Read contig, int pos, final int depth){
		assert(pos>=k-1) : pos;
		int depth2=Tools.min(depth, maxDepth);
		while(depth2>=openRepeats.size()){
			openRepeats.add(new Repeat(null, -2, openRepeats.size(), k, maxGap, minRepeat, 'R'));
		}
		for(int i=depth2; i>=minDepth; i--){//Potentially quadratic
			Repeat r=openRepeats.get(i);
			Repeat old=r.increment(contig, pos, depth);
			if(old!=null){addRepeat(old);}
		}
	}

	/** Sorts and greedily removes subsumed closed snapshots in place. */
	int subsumeClosed(boolean weak){return subsume(closedRepeats, weak);}

	/** Returns the shared range buffer for closed repeats; the next range call overwrites it. */
	ArrayList<CRange> closedToRanges(boolean merge){
		return toRanges(closedRepeats, merge, rangeBuffer);
	}

	/** Returns the same shared range buffer for recently accepted repeats. */
	ArrayList<CRange> recentToRanges(boolean merge){
		return toRanges(recent, merge, rangeBuffer);
	}

	/** Creates reads from closed regions after positional containment filtering. */
	ArrayList<Read> fetchRepeatSequence(){return fetchRepeatSequence(closedRepeats, maxGap, k);}

	/** Copies the list structure, filters contained regions, and creates independent base arrays.
	 * Repeat and contig objects themselves remain borrowed while copying. */
	static ArrayList<Read> fetchRepeatSequence(ArrayList<Repeat> repeats0, int maxGap, int k){
		ArrayList<Read> reads=new ArrayList<Read>();
		if(repeats0.isEmpty()){return reads;}
		ArrayList<Repeat> repeats=(ArrayList<Repeat>) repeats0.clone();
		removeFullyContained(repeats, maxGap, k);
		for(Repeat pete : repeats){
			Read r=pete.toRead();
			reads.add(r);
		}
		return reads;
	}

	/** Sorts by position and removes following regions subsumed by the last retained region.
	 * Mutates/compacts the supplied list and returns the number removed. */
	public static int subsume(ArrayList<Repeat> list, boolean weak){
		if(list.isEmpty()){return 0;}
		list.sort(Repeat.PosComparator.comparator);

		int removed=0;
		Repeat current=list.get(0);
		for(int i=1; i<list.size(); i++){
			Repeat r=list.get(i);
			assert(current!=r);
			if(current.subsumes(r, weak)){
//				if(printSubsumes) {
//					System.err.println(current+"\n"+r);
//				}
				list.set(i, null);
				removed++;
			}else{
				current=r;
			}
		}
		if(removed>0){Tools.condenseStrict(list);}
		return removed;
	}

	/**
	 * Removes repeats that are fully contained within other repeats.
	 * Sorts repeats by position and removes those spanned by earlier repeats.
	 * Includes assertion checking for unexpected overlaps.
	 *
	 * @param repeats Repeat list to process (will be modified)
	 * @param maxGap Maximum gap size (used for overlap assertions)
	 * @param k K-mer length (used for overlap assertions)
	 * @return Number of repeats removed
	 */
	public static int removeFullyContained(final ArrayList<Repeat> repeats, int maxGap, int k){//maxGap and k are just for an assertion
		if(repeats.size()<2){return 0;}
		repeats.sort(PosComparator2.comparator);
		Repeat current=repeats.get(0);
		int removed=0;
		for(int i=1; i<repeats.size(); i++){
			final Repeat r=repeats.get(i);
			if(current.spans(r)){
				repeats.set(i, null);
				removed++;
			}else{
				assert(!current.overlaps(r) || maxGap<=k) : "\n"+current+"\n"+r;
				current=r;
			}
		}
		if(removed>0){
			Tools.condenseStrict(repeats);
		}
		return removed;
	}

	/** Replaces ranges with sorted inclusive coordinates, optionally using CRange's merge policy.
	 * Allocates a list only when ranges is null. New range objects retain the contig
	 * reference; the source repeat list is not changed. */
	public static ArrayList<CRange> toRanges(ArrayList<Repeat> repeats, boolean merge, ArrayList<CRange> ranges){
		if(ranges==null){ranges=new ArrayList<CRange>();}
		ranges.clear();
		if(repeats.size()==1){ranges.add(repeats.get(0).toRange());}
		if(repeats.size()<2){return ranges;}
		{
			Repeat current=repeats.get(repeats.size()-1);
			for(int i=repeats.size()-2; i>=0; i--){
				Repeat r=repeats.get(i);
				if(current.spans(r)){
					//do nothing
				}else{
					ranges.add(current.toRange());
					current=r;
				}
			}
			ranges.add(current.toRange());
			Collections.sort(ranges);
		}
		if(merge){CRange.mergeList(ranges, false);}
		return ranges;
	}

	/** Depth-window length. */
	final int k;
	/** Lowest tracked depth threshold. */
	final int minDepth;
	/** Highest tracked threshold, independent of observed counts passed to each slot. */
	final int maxDepth;
	/** Minimum span for retaining a closed depth/entropy region. */
	final int minRepeat;
	/** Maximum missing-window gap within a depth region. */
	final int maxGap;
	/** Whether last-repeat subsumption ignores gap totals. */
	final boolean weakSubsumes;
	/** Entropy alphabet selection for lazily created trackers. */
	boolean amino;
	/** Entropy k and maximum cached window size. */
	final int ek, ew;

	/** Reusable open slots indexed by depth threshold. */
	final ArrayList<Repeat> openRepeats=new ArrayList<Repeat>();
	/** Accepted snapshots awaiting batch output/masking. */
	final ArrayList<Repeat> closedRepeats=new ArrayList<Repeat>();
	/** Accepted snapshots for the current detector's base-count calculation. */
	final ArrayList<Repeat> recent=new ArrayList<Repeat>();
	/** Retained snapshots accumulated across retired batches. */
	final ArrayList<Repeat> oldRepeats=new ArrayList<Repeat>();
	/** Borrowed output buffer shared by both range-conversion methods. */
	private final ArrayList<CRange> rangeBuffer=new ArrayList<CRange>();
	/** Last accepted snapshot, possibly already retired into oldRepeats. */
	Repeat lastRepeat=null;
	/** Lazy worker-local trackers indexed by window length. */
	final EntropyTracker[] eta;
//	final EntropyTracker et;

}
