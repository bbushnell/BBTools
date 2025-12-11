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
	
	private ReadComparatorFlowcell(){}
	
	/**
	 * Compares two reads based on flowcell coordinates with ascending/descending control.
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return Negative if r1 < r2, positive if r1 > r2, zero if equal
	 */
	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return ascending*x;
	}
	
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
	
	private int ascending=1;
	
	/** Sets the sort direction for comparisons.
	 * @param asc true for ascending order, false for descending */
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}

	public ThreadLocal<FlowcellCoordinate> tlc1=new ThreadLocal<FlowcellCoordinate>();
	public ThreadLocal<FlowcellCoordinate> tlc2=new ThreadLocal<FlowcellCoordinate>();
	
	public static final ReadComparatorFlowcell comparator=new ReadComparatorFlowcell();
	
}
