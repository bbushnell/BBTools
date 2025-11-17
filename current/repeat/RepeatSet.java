package repeat;

import java.util.ArrayList;
import java.util.Collections;

import repeat.Repeat.PosComparator2;
import shared.Tools;
import stream.Read;
import structures.CRange;
import tracker.EntropyTracker;

/**
 * Manages collections of sequence repeats during analysis and processing.
 * Tracks open repeats being extended, closed completed repeats, and provides
 * utilities for filtering, merging, and converting repeats to various formats.
 * Includes entropy tracking and subsumption logic for repeat deduplication.
 *
 * @author Brian Bushnell
 * @date 2025
 */
public class RepeatSet {
	
	/**
	 * Constructs a RepeatSet with specified parameters for repeat detection.
	 * Initializes entropy tracker array for sequence complexity analysis.
	 *
	 * @param k_ K-mer length for repeat detection
	 * @param minDepth_ Minimum depth threshold for valid repeats
	 * @param maxDepth_ Maximum depth to track in open repeats
	 * @param minRepeat_ Minimum repeat length to retain
	 * @param maxGap_ Maximum gap size allowed within repeats
	 * @param weakSubsumes_ Whether to use weak subsumption logic
	 * @param amino_ Whether processing amino acid sequences
	 * @param ek_ Entropy tracker k-mer length
	 * @param ew_ Entropy tracker maximum window size
	 */
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
	
	/**
	 * Gets or creates an EntropyTracker for the specified sequence length.
	 * Uses window size between ek and ew based on sequence length.
	 * @param len Sequence length for entropy calculation
	 * @return EntropyTracker instance for the appropriate window size
	 */
	EntropyTracker getET(int len) {
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
	
	/** Moves all closed repeats to the old repeats collection.
	 * Clears the closed repeats list for new additions. */
	void retireClosed() {
		oldRepeats.addAll(closedRepeats);
		closedRepeats.clear();
	}
	
	/**
	 * Adds a repeat to the closed collection if not subsumed by last repeat.
	 * Calculates repeat statistics using entropy tracking and updates recent list.
	 * @param r The repeat to add
	 */
	void addRepeat(Repeat r){
		if(lastRepeat!=null && lastRepeat.subsumes(r, weakSubsumes)) {
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
	
	/**
	 * Collects residual repeats from open collection that meet criteria.
	 * Processes repeats from maxDepthSeen down to minDepth, adding valid ones
	 * to closed collection and clearing their open slots.
	 * @param maxDepthSeen Maximum depth level to process residual repeats from
	 */
	void collectResidual(int maxDepthSeen) {
		maxDepthSeen=Tools.min(maxDepthSeen, openRepeats.size()-1); //In case maxDepth is set.
		assert(maxDepthSeen<openRepeats.size()) : maxDepthSeen+", "+openRepeats.size();
		for(int i=maxDepthSeen; i>=minDepth; i--) {
			Repeat r=openRepeats.get(i);
			if(r.contigNum>=0 && r.depth>1 && r.length()>=minRepeat){
				addRepeat(r.clone());
				r.clear();
			}
		}
	}
	
	/** TODO: This can go from O(N^2) to O(N) if it follows contours instead of incrementing all.
	 * Also, make sure min depth is implemented. */
	void increment(Read contig, int pos, final int depth) {
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
	
	/**
	 * Removes subsumed repeats from the closed repeats collection.
	 * @param weak Whether to use weak subsumption logic
	 * @return Number of repeats removed
	 */
	int subsumeClosed(boolean weak) {return subsume(closedRepeats, weak);}
	
	/**
	 * Converts closed repeats to coordinate ranges.
	 * @param merge Whether to merge overlapping ranges
	 * @return List of coordinate ranges representing closed repeats
	 */
	ArrayList<CRange> closedToRanges(boolean merge){
		return toRanges(closedRepeats, merge, rangeBuffer);
	}
	
	/**
	 * Converts recent repeats to coordinate ranges.
	 * @param merge Whether to merge overlapping ranges
	 * @return List of coordinate ranges representing recent repeats
	 */
	ArrayList<CRange> recentToRanges(boolean merge){
		return toRanges(recent, merge, rangeBuffer);
	}
	
	/**
	 * Extracts repeat sequences from closed repeats as Read objects.
	 * Removes fully contained repeats before extraction.
	 * @return List of Read objects containing repeat sequences
	 */
	ArrayList<Read> fetchRepeatSequence(){return fetchRepeatSequence(closedRepeats, maxGap, k);}
	
	/**
	 * Extracts repeat sequences from specified repeat collection.
	 * Clones the input list, removes fully contained repeats, then converts
	 * remaining repeats to Read objects.
	 *
	 * @param repeats0 Input repeat collection
	 * @param maxGap Maximum gap size for overlap detection
	 * @param k K-mer length for overlap detection
	 * @return List of Read objects containing repeat sequences
	 */
	static ArrayList<Read> fetchRepeatSequence(ArrayList<Repeat> repeats0, int maxGap, int k){
		ArrayList<Read> reads=new ArrayList<Read>();
		if(repeats0.isEmpty()){return reads;}
		ArrayList<Repeat> repeats=(ArrayList<Repeat>) repeats0.clone();
		removeFullyContained(repeats, maxGap, k);
		for(Repeat pete : repeats) {
			Read r=pete.toRead();
			reads.add(r);
		}
		return reads;
	}
	
	/**
	 * Removes subsumed repeats from a sorted repeat list.
	 * Sorts the list by position, then removes repeats that are subsumed
	 * by earlier repeats using the specified subsumption logic.
	 *
	 * @param list Repeat list to process (will be modified)
	 * @param weak Whether to use weak subsumption logic
	 * @return Number of repeats removed
	 */
	public static int subsume(ArrayList<Repeat> list, boolean weak) {
		if(list.isEmpty()) {return 0;}
		list.sort(Repeat.PosComparator.comparator);
		
		int removed=0;
		Repeat current=list.get(0);
		for(int i=1; i<list.size(); i++) {
			Repeat r=list.get(i);
			assert(current!=r);
			if(current.subsumes(r, weak)) {
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
	
	public static int removeFullyContained(final ArrayList<Repeat> repeats, int maxGap, int k) {//maxGap and k are just for an assertion
		if(repeats.size()<2) {return 0;}
		repeats.sort(PosComparator2.comparator);
		Repeat current=repeats.get(0);
		int removed=0;
		for(int i=1; i<repeats.size(); i++) {
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
	
	/**
	 * Converts repeat collection to coordinate ranges.
	 * Processes repeats in reverse order to handle spanning relationships,
	 * then sorts ranges and optionally merges overlapping ones.
	 *
	 * @param repeats Input repeat collection
	 * @param merge Whether to merge overlapping ranges
	 * @param ranges Output range collection (may be null)
	 * @return List of coordinate ranges, optionally merged
	 */
	public static ArrayList<CRange> toRanges(ArrayList<Repeat> repeats, boolean merge, ArrayList<CRange> ranges){
		if(ranges==null){ranges=new ArrayList<CRange>();}
		ranges.clear();
		if(repeats.size()==1) {ranges.add(repeats.get(0).toRange());}
		if(repeats.size()<2){return ranges;}
		{
			Repeat current=repeats.get(repeats.size()-1);
			for(int i=repeats.size()-2; i>=0; i--) {
				Repeat r=repeats.get(i);
				if(current.spans(r)) {
					//do nothing
				}else {
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
	
	/** K-mer length for repeat detection */
	final int k;
	/** Minimum depth threshold for valid repeats */
	final int minDepth;
	/** Maximum depth to track in open repeats */
	final int maxDepth;
	/** Minimum repeat length to retain */
	final int minRepeat;
	/** Maximum gap size allowed within repeats */
	final int maxGap;
	/** Whether to use weak subsumption logic for repeat comparison */
	final boolean weakSubsumes;
	/** Whether processing amino acid sequences */
	boolean amino;
	final int ek, ew;
	
	/** Collection of repeats currently being extended */
	final ArrayList<Repeat> openRepeats=new ArrayList<Repeat>();
	/** Collection of completed repeats ready for processing */
	final ArrayList<Repeat> closedRepeats=new ArrayList<Repeat>();
	/** Collection of recently added repeats */
	final ArrayList<Repeat> recent=new ArrayList<Repeat>();
	/** Collection of previously processed repeats */
	final ArrayList<Repeat> oldRepeats=new ArrayList<Repeat>();
	/** Reusable buffer for coordinate range operations */
	private final ArrayList<CRange> rangeBuffer=new ArrayList<CRange>();
	/** Most recently added repeat for subsumption checking */
	Repeat lastRepeat=null;
	/** Array of entropy trackers indexed by window size */
	final EntropyTracker[] eta;
//	final EntropyTracker et;
	
}
