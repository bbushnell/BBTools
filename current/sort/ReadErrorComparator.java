package sort;

import stream.Read;

/**
 * Comparator for sorting Read objects by error rates and quality metrics.
 * Orders reads from lowest to highest error count, with secondary sorting
 * by length (longest first), expected errors, numeric ID, and string ID.
 * Handles both single reads and paired-end reads with mates.
 *
 * @author Brian Bushnell
 * @date May 30, 2013
 */
public final class ReadErrorComparator extends ReadComparator{
	
	/**
	 * Compares two reads using hierarchical sorting criteria.
	 * Primary sort: total error count (read + mate) ascending
	 * Secondary sort: total length (read + mate) descending
	 * Tertiary sort: expected errors (read + mate) ascending
	 * Quaternary sort: numeric ID ascending
	 * Final sort: string ID lexicographically
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return negative if r1 < r2, positive if r1 > r2, zero if equal
	 */
	@Override
	public int compare(Read r1, Read r2) {
		return compareInner(r1, r2)*ascending;
	}
	
	public int compareInner(Read r1, Read r2) {

		int a=(r1.errors+(r1.mate==null ? 0 : r1.mate.errors));
		int b=(r2.errors+(r2.mate==null ? 0 : r2.mate.errors));
		if(a!=b){return a-b;}
		
		a=(r1.length()+(r1.mate==null ? 0 : r1.mateLength()));
		b=(r2.length()+(r2.mate==null ? 0 : r2.mateLength()));
		if(a!=b){return b-a;}
		
		float a2=(r1.expectedErrors(true, 0)+(r1.mate==null ? 0 : r1.mate.expectedErrors(true, 0)));
		float b2=(r2.expectedErrors(true, 0)+(r2.mate==null ? 0 : r2.mate.expectedErrors(true, 0)));
		if(a2!=b2){return a2>b2 ? 1 : -1;}
		
		if(r1.numericID<r2.numericID){return -1;}
		else if(r1.numericID>r2.numericID){return 1;}
		
		if(!r1.id.equals(r2.id)){return r1.id.compareTo(r2.id);}
		return 0;
	}
	
	public static final ReadErrorComparator comparator=new ReadErrorComparator();

	@Override
	public void setAscending(boolean asc){ascending=asc ? 1 : -1;}
	
	@Override
	public final int ascendingMult() {return ascending;}
	@Override
	public final String name() {return "Error";}
	
	int ascending=1;
	
}
