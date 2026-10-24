package sort;

import java.util.Comparator;

import shared.Shared;
import shared.Tools;
import stream.Read;


/**
 * Comparator for sorting mapped reads by genomic position.
 * Linked roots use one stable coordinate anchor: mapped primary, mapped mate, or root
 * when neither end is mapped. Equal anchors are ordered by root shape, then by the
 * same primary-or-mate position and detail fields used inside each shape.
 *
 * @author Brian Bushnell
 * @date Oct 27, 2014
 */
public class ReadComparatorMapping implements Comparator<Read>{

	/**
	 * Compares two root reads by stable anchor coordinate, then root mapping shape:
	 * both ends mapped and properly paired, both ends mapped but not paired, mapped
	 * primary with an unmapped mate, mapped unlinked primary, unmapped primary with
	 * mapped mate, linked unmapped, then unlinked unmapped. Same-shape ties preserve
	 * the prior primary/mate detail path: fully paired roots compare both ends, primary
	 * shapes compare the primary only, and mate-only or fully unmapped linked shapes
	 * compare the mate only.
	 *
	 * @param a First read to compare
	 * @param b Second read to compare
	 * @return Negative if {@code a < b}, positive if {@code a > b}, zero if equal
	 */
	@Override
	public int compare(Read a, Read b){
		int x=compare2(anchor(a), anchor(b), false);
		if(x!=0){return x;}

		final int shapeA=shape(a), shapeB=shape(b);
		x=shapeA-shapeB;
		if(x!=0){return x;}

		if(shapeA==0){
			x=compare2(a, b);
			if(x!=0){return x;}
			x=compare2(a.mate, b.mate);
			if(x!=0){return x;}
			x=compare3(a, b);
			if(x!=0){return x;}
			return compare3(a.mate, b.mate);
		}else if(shapeA<4 || shapeA==6){
			x=compare2(a, b);
			if(x!=0){return x;}
			return compare3(a, b);
		}else{
			x=compare2(a.mate, b.mate);
			if(x!=0){return x;}
			return compare3(a.mate, b.mate);
		}
	}
	
	/** Position-level comparison: mapped before unmapped, then chromosome, strand, and coordinate
	 * (start on the plus strand, stop on the minus), then paired before unpaired. */
	public int compare2(Read a, Read b){//TODO: This could all be packed in a long (the ID).
		return compare2(a, b, true);
	}

	/** Position-level comparison with optional paired-flag comparison. */
	private int compare2(Read a, Read b, boolean comparePaired){
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
		
		if(comparePaired && a.paired()!=b.paired()){return a.paired() ? -1 : 1;}
		return 0;
	}
	
	/** Position-level comparison between a read and its mate's counterpart, used when only one
	 * end of a pair is mapped. */
	public int compareCross(Read a, Read b){
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
	
	/** Secondary comparison for reads at the same position: longer first, then perfect-match,
	 * then match string, then coordinate, then quality (higher first), numericID, and string id. */
	public int compare3(Read a, Read b){
		if(a.length()!=b.length()){
			return b.length()-a.length(); //Preferentially puts longer reads first
		}
		if(a.perfect()!=b.perfect()){return a.perfect() ? -1 : 1;}
		int x;
		
		x=compareMatchStrings(a.match, b.match);
		if(x!=0){return x;}
		
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

	/** Stable coordinate anchor: a mapped primary, else a mapped mate when present, else the root. */
	private static Read anchor(Read r){
		return r.mapped() || r.mate==null || !r.mate.mapped() ? r : r.mate;
	}

	/** Root shape after the stable anchor coordinate is tied. */
	private static int shape(Read r){
		if(r.mapped()){
			if(r.mate!=null && r.mate.mapped()){return r.paired() ? 0 : 1;}
			return r.mate==null ? 3 : 2;
		}
		if(r.mate!=null && r.mate.mapped()){return 4;}
		return r.mate==null ? 6 : 5;
	}
	
	/** Lexicographically compares two byte arrays to the shorter length; null sorts after non-null. */
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
	
	/** Compares two alignment match strings, ordering by the first position whose insertion/deletion
	 * status differs (insertions, then deletions, sort later), then by match length. */
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
		return a.length-b.length;
	}

	/** When true, paired reads are expected on the same strand (alters strand handling in the position
	 * comparisons); currently unsupported and guarded by asserts. */
	public static boolean SAME_STRAND_PAIRS=false;
	
}
