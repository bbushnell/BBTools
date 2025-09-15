package sort;

import stream.Read;

/**
 * Sorts longest reads first
 * @author Brian Bushnell
 * @date Jul 19, 2013
 *
 */
public final class ReadLengthComparator extends ReadComparator {
	
	private ReadLengthComparator(){}
	
	@Override
	public int compare(Read a, Read b) {
		int x=compareInner(a, b);
		if(x==0){x=compareInner(a.mate, b.mate);}
		if(x==0){x=a.id.compareTo(b.id);}
		if(x==0){x=a.numericID>b.numericID ? 1 : a.numericID<b.numericID ? -1 : 0;}
		return ascending*x;
	}

	private static int compareInner(Read a, Read b) {
		if(a==b){return 0;}
		if(a==null){return 1;}
		if(b==null){return -1;}
		int x=a.length()-b.length();
		return x;
	}
	
	public static final ReadLengthComparator comparator=new ReadLengthComparator();
	
	private int ascending=-1;
	
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}
	
}
