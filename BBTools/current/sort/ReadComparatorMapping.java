package sort;

import java.util.Comparator;

import shared.Shared;
import shared.Tools;
import stream.Read;

/**
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */

public class ReadComparatorMapping implements Comparator<Read> {

	@Override
	public int compare(Read a, Read b) {
		
		if(a.mate==null){
			int x=compare2(a, b);
			if(x!=0){return x;}
			return compare3(a, b);
		}else{
			
			if(a.mapped() && b.mapped()){
				int x=compare2(a, b);
				if(x!=0){return x;}
				
				if(a.paired() && b.paired()){
					x=compare2(a.mate, b.mate);
					if(x!=0){return x;}
					x=compare3(a, b);
					if(x!=0){return x;}
					x=compare3(a.mate, b.mate);
					return x;
				}else{
					assert(!a.paired() && !b.paired());
					return compare3(a, b);
				}
			}
			
			if(!a.mapped() && !b.mapped()){
				int x=compare2(a.mate, b.mate);
				if(x!=0){return x;}
				return compare3(a.mate, b.mate);
			}else if(a.mapped()){
				if(a.paired()){
					int x=compare2(a.mate, b.mate);
					if(x!=0){return x;}
					return -1;
				}else{
					int x=compareCross(a, b.mate);
					if(x!=0){return x;}
					return -1;
				}
			}else if(b.mapped()){
				if(b.paired()){
					int x=compare2(a.mate, b.mate);
					if(x!=0){return x;}
					return 1;
				}else{
					int x=compareCross(b, a.mate);
					if(x!=0){return 0-x;}
					return 1;
				}
			}else{
				assert(false) : a.mapped()+", "+a.paired()+", "+b.mapped()+", "+b.paired()+", "+a.mate.mapped()+", "+b.mate.mapped();
			}
			
			//I think this is unreachable...
			return compare3(a, b);
		}
	}
	
	public int compare2(Read a, Read b) {//TODO: This could all be packed in a long (the ID).
		if(a.mapped() && !b.mapped()){return -1;}
		if(b.mapped() && !a.mapped()){return 1;}
		if(a.chrom!=b.chrom){return a.chrom-b.chrom;}
		if(a.strand()!=b.strand()){return a.strand()-b.strand();}
		
		assert(!SAME_STRAND_PAIRS) : "TODO";
		if(a.strand()==Shared.PLUS){
			if(a.start!=b.start){return a.start-b.start;}
		}else{
			if(a.stop!=b.stop){return a.stop-b.stop;}
		}
		
		if(a.paired()!=b.paired()){return a.paired() ? -1 : 1;}
		return 0;
	}
	
	public int compareCross(Read a, Read b) {
		if(a.mapped() && !b.mapped()){return -1;}
		if(b.mapped() && !a.mapped()){return 1;}
		if(a.chrom!=b.chrom){return a.chrom-b.chrom;}
		if(SAME_STRAND_PAIRS){
			if(a.strand()!=b.strand()){
				return a.strand()-b.strand();
			}
		}else{
			if(a.strand()==b.strand()){
				return a.strand()==0 ? -1 : 1;
			}
		}
		if(a.start!=b.start){return a.start-b.start;}
		if(a.paired()!=b.paired()){return a.paired() ? -1 : 1;}
		return 0;
	}
	
	public int compare3(Read a, Read b){
		if(a.length()!=b.length()){
			return b.length()-a.length(); //Preferentially puts longer reads first
		}
		if(a.perfect() != b.perfect()){return a.perfect() ? -1 : 1;}
		int x;
		
		if(a.match!=null && b.match!=null){
			x=compareMatchStrings(a.match, b.match);
			if(x!=0){return x;}
		}
		
		assert(!SAME_STRAND_PAIRS) : "TODO";
		if(a.strand()==Shared.PLUS){
			if(a.start!=b.start){return a.start-b.start;} //This line should be dead code
			if(a.stop!=b.stop){return a.stop-b.stop;}
		}else{
			if(a.stop!=b.stop){return a.stop-b.stop;} //This line should be dead code
			if(a.start!=b.start){return a.start-b.start;}
		}
		
		x=compareVectors(a.quality, b.quality);
		if(x!=0){return 0-x;}
//		if(a.stop!=b.stop){return a.stop-b.stop;}
		if(a.numericID!=b.numericID){return a.numericID>b.numericID ? 1 : -1;}
		return a.id.compareTo(b.id);
	}
	
	public int compareVectors(final byte[] a, final byte[] b){
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
	
	public int compareMatchStrings(final byte[] a, final byte[] b){
		if(a==null || b==null){
			if(a==null && b!=null){return 1;}
			if(a!=null && b==null){return -1;}
			return 0;
		}
		final int lim=Tools.min(a.length, b.length);
		for(int i=0; i<lim; i++){
			if(a[i]!=b[i]){
				boolean ad=(a[i]=='D');
				boolean bd=(b[i]=='D');
				boolean ai=(a[i]=='I' || a[i]=='X' || a[i]=='Y');
				boolean bi=(b[i]=='I' || b[i]=='X' || b[i]=='Y');
				if(ai!=bi){return ai ? 1 : -1;}
				if(ad!=bd){return ad ? 1 : -1;}
			}
		}
		return 0;
	}

	public static boolean SAME_STRAND_PAIRS=false;
	
}
