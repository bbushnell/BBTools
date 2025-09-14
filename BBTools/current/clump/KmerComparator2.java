package clump;

import java.util.Comparator;

import shared.Tools;
import stream.Read;

/** A minimal KmerComparator without the heavyweight auto-hashing methods of the original */
public abstract class KmerComparator2 implements Comparator<Read>{

	@Override
	public final int compare(Read a, Read b) {
		final ReadKey keyA=(ReadKey)a.obj;
		final ReadKey keyB=(ReadKey)b.obj;
		
		int x=compare(keyA, keyB);
		if(x==0){
			x=compareSequence(a, b, 0);
		}
		return x==0 ? a.id.compareTo(b.id) : x;
	}
	
	//This gets overriden
	public abstract int compare(ReadKey a, ReadKey b);
	
	public static final int compareSequence(Read a, Read b, int depth){
		int x=compareSequence(a.bases, b.bases);
		if(x!=0){return x;}
		if(a.mate!=null){x=compareSequence(a.mate.bases, b.mate.bases);}
		if(x!=0){return x;}
		return compareQuality(a, b);
	}
	
	public static final int compareSequence(final byte[] a, final byte[] b){
		if(a==null || b==null){
			if(a==null && b!=null){return 1;}
			if(a!=null && b==null){return -1;}
			return 0;
		}
		if(a.length!=b.length){
			return b.length-a.length;
		}
		for(int i=0, lim=a.length; i<lim; i++){
			int x=a[i]-b[i];
			if(x!=0){return x;}
		}
		return 0;
	}
	
	//Not optimal, but fast.  This function is probably not very important.
	public static final int compareQuality(Read a, Read b){
		if(a.quality==null){return 0;}
		int qa=Tools.sumInt(a.quality);
		int qb=Tools.sumInt(b.quality);
		return qb-qa;
	}
	
}
