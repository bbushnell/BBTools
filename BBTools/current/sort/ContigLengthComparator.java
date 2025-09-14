package sort;

import java.util.Comparator;

import assemble.Contig;
import shared.Tools;

/**
 * Sorts longest contig first
 * @author Brian Bushnell
 * @date Jul 12, 2018
 *
 */
public final class ContigLengthComparator implements Comparator<Contig> {
	
	private ContigLengthComparator(){}
	
	@Override
	public int compare(Contig a, Contig b) {
		int x=compareInner(a, b);
		if(x==0){x=a.coverage>b.coverage ? 1 : a.coverage<b.coverage ? -1 : 0;}
		if(x==0){x=compareVectors(a.bases, b.bases);}
		if(x==0){x=a.id>b.id ? 1 : a.id<b.id ? -1 : 0;}
		return ascending*x;
	}
	
	private static int compareInner(Contig a, Contig b) {
		if(a==b){return 0;}
		if(a==null){return 1;}
		if(b==null){return -1;}
		int x=a.length()-b.length();
		return x;
	}
	
	private static int compareVectors(final byte[] a, final byte[] b){
		if(a==null || b==null){
			if(a==null && b!=null){return 1;}
			if(a!=null && b==null){return -1;}
			return 0;
		}
		final int lim=Tools.min(a.length, b.length);
		for(int i=0; i<lim; i++){
			if(a[i]<b[i]){return -1;}
			if(a[i]>b[i]){return 1;}
		}
		return 0;
	}
	
	public static final ContigLengthComparator comparator=new ContigLengthComparator();
	
	private int ascending=-1;
	
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}
	
}
