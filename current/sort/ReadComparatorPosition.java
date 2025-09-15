package sort;

import stream.Read;
import stream.SamLine;
import var2.ScafMap;

/**
 * @author Brian Bushnell
 * @date November 20, 2016
 *
 */

public final class ReadComparatorPosition extends ReadComparator {
	
	private ReadComparatorPosition(){}
	
	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return ascending*x;
	}
	
	public static int compareInner(Read r1, Read r2) {
		int x=compareInner(r1.samline, r2.samline);
		if(x!=0){return x;}
		if(r1.id==null && r2.id==null){return 0;}
		if(r1.id==null){return -1;}
		if(r2.id==null){return 1;}
		return r1.id.compareTo(r2.id);
	}
	
	public static int compareInner(SamLine a, SamLine b) {
		if(a.scafnum<0){a.setScafnum(scafMap);}
		if(b.scafnum<0){b.setScafnum(scafMap);}
		if(a.scafnum!=b.scafnum){return a.scafnum-b.scafnum;}
		if(a.pos!=b.pos){return a.pos-b.pos;}
		if(a.strand()!=b.strand()){return a.strand()-b.strand();}
		if(a.pnext!=b.pnext){return a.pnext-b.pnext;}
		if(a.pairnum()!=b.pairnum()){return a.pairnum()-b.pairnum();}
		return 0;
	}
	
	private int ascending=1;
	
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}

	public static final ReadComparatorPosition comparator=new ReadComparatorPosition();
	public static ScafMap scafMap=null;
	
}
