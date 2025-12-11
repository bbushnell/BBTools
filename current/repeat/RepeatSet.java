package repeat;

import java.util.ArrayList;
import java.util.Collections;

import repeat.Repeat.PosComparator2;
import shared.Tools;
import stream.Read;
import structures.CRange;
import tracker.EntropyTracker;

public class RepeatSet {
	
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
	
	void retireClosed() {
		oldRepeats.addAll(closedRepeats);
		closedRepeats.clear();
	}
	
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
	
	/**
	 * Increments repeat tracking at the specified position and depth.
	 * Expands open repeats array as needed and processes all depth levels
	 * from current depth down to minDepth. May add completed repeats.
	 *
	 * @param contig The contig being processed
	 * @param pos Position in the contig (must be >= k-1)
	 * @param depth Coverage depth at this position
	 */
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
	
	int subsumeClosed(boolean weak) {return subsume(closedRepeats, weak);}
	
	ArrayList<CRange> closedToRanges(boolean merge){
		return toRanges(closedRepeats, merge, rangeBuffer);
	}
	
	ArrayList<CRange> recentToRanges(boolean merge){
		return toRanges(recent, merge, rangeBuffer);
	}
	
	ArrayList<Read> fetchRepeatSequence(){return fetchRepeatSequence(closedRepeats, maxGap, k);}
	
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
	
	final int k;
	final int minDepth;
	final int maxDepth;
	final int minRepeat;
	final int maxGap;
	final boolean weakSubsumes;
	boolean amino;
	final int ek, ew;
	
	final ArrayList<Repeat> openRepeats=new ArrayList<Repeat>();
	final ArrayList<Repeat> closedRepeats=new ArrayList<Repeat>();
	final ArrayList<Repeat> recent=new ArrayList<Repeat>();
	final ArrayList<Repeat> oldRepeats=new ArrayList<Repeat>();
	private final ArrayList<CRange> rangeBuffer=new ArrayList<CRange>();
	Repeat lastRepeat=null;
	final EntropyTracker[] eta;
//	final EntropyTracker et;
	
}
