package sort;

import hiseq.FlowcellCoordinate;
import stream.Read;


/**
 * Comparator for reads based on flowcell coordinates from their identifiers.
 * Extracts HiSeq-style coordinate information from read IDs and sorts by position.
 * Uses thread-local FlowcellCoordinate objects for parsing efficiency.
 *
 * @author Brian Bushnell
 * @date Oct 27, 2014
 */
public final class ReadComparatorFlowcell extends ReadComparator {

	/** Private constructor; use the ascending/descending singletons. */
	private ReadComparatorFlowcell(int mult_){mult=mult_;}

	/**
	 * Compares two reads based on flowcell coordinates with ascending/descending control.
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return Negative if r1 < r2, positive if r1 > r2, zero if equal
	 */
	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return mult*x;
	}

	/** Compares two reads by flowcell coordinates parsed from their ids (lane, tile, x, y),
	 * with pair number as a tiebreaker. Uses thread-local parse buffers. */
	public int compareInner(Read r1, Read r2) {
		if(r1.id==null && r2.id==null){return r1.pairnum()-r2.pairnum();}
		if(r1.id==null){return -1;}
		if(r2.id==null){return 1;}

		FlowcellCoordinate fc1=tlc1.get(), fc2=tlc2.get();
		if(fc1==null){
			fc1=new FlowcellCoordinate();
			fc2=new FlowcellCoordinate();
			tlc1.set(fc1);
			tlc2.set(fc2);
		}
		fc1.setFrom(r1.id);
		fc2.setFrom(r2.id);

		int x=fc1.compareTo(fc2);
		if(x==0){return r1.pairnum()-r2.pairnum();}
		return x;
	}

	@Override
	public ReadComparator getComparator(boolean asc){return asc ? ascending : descending;}

	@Override
	public final int ascendingMult() {return mult;}
	@Override
	public final String name() {return "Flowcell Coordinate";}

	/** Sort direction multiplier: +1 for ascending, -1 for descending (immutable). */
	private final int mult;

	/** Thread-local scratch buffer for parsing the first read's flowcell coordinate. */
	public ThreadLocal<FlowcellCoordinate> tlc1=new ThreadLocal<FlowcellCoordinate>();
	/** Thread-local scratch buffer for parsing the second read's flowcell coordinate. */
	public ThreadLocal<FlowcellCoordinate> tlc2=new ThreadLocal<FlowcellCoordinate>();

	/** Ascending-order singleton. */
	public static final ReadComparatorFlowcell ascending=new ReadComparatorFlowcell(1);
	/** Descending-order singleton. */
	public static final ReadComparatorFlowcell descending=new ReadComparatorFlowcell(-1);
	/** Default (ascending) singleton; alias retained for selection/identity call sites. */
	public static final ReadComparatorFlowcell comparator=ascending;

}
