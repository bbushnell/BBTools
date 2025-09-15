package sort;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */

public final class ReadComparatorName extends ReadComparator {
	
	private ReadComparatorName(){}
	
	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return ascending*x;
	}
	
	public static int compareInner(Read r1, Read r2) {
		
		if(r1.id==null && r2.id==null){return r1.pairnum()-r2.pairnum();}
		if(r1.id==null){return -1;}
		if(r2.id==null){return 1;}
		int x=r1.id.compareTo(r2.id);
		if(x==0){return r1.pairnum()-r2.pairnum();}
		return x;
	}
	
	private int ascending=1;
	
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}

	public static final ReadComparatorName comparator=new ReadComparatorName();
	
}
