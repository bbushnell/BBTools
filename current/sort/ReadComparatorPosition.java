package sort;

import stream.Read;
import stream.SamLine;
import var2.ScafMap;


/**
 * Comparator for sorting genomic reads based on their position and metadata.
 * Provides multi-level comparison using scaffold, position, strand, and pair information.
 * Orders reads systematically for alignment processing workflows.
 *
 * @author Brian Bushnell
 * @date November 20, 2016
 */
public final class ReadComparatorPosition extends ReadComparator {

	/** Private constructor; use the ascending/descending singletons. */
	private ReadComparatorPosition(int mult_){mult=mult_;}

	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return mult*x;
	}

	/** Compares two reads by SAM alignment position when both have SamLines (a read with a
	 * SamLine sorts before one without); falls back to string id. Caller applies sort direction. */
	public static int compareInner(Read r1, Read r2) {
		final SamLine sl1=r1.samline, sl2=r2.samline;
		if(sl1!=null && sl2!=null) {
			int x=compareInner(r1.samline, r2.samline);
			if(x!=0){return x;}
		}else if(sl1==null && sl2==null) {
			//fall through
		}else if(sl1==null) {
			//TODO: Possible bug [sort/ReadComparatorPosition#001] - a read WITHOUT a SamLine sorts FIRST here (sl1==null ->
			//-1 => r1 before r2; sl2==null -> +1 => r1 after r2). This contradicts BOTH (a) this method's javadoc ("a read
			//with a SamLine sorts before one without") AND (b) the SamLine-level convention below (line ~54) where an
			//UNMAPPED read (scafnum<0, also positionless) sorts LAST. So no-samline and unmapped-samline reads - both
			//positionless - land at OPPOSITE ends. Still a valid total order (no crash), but the -1/+1 look swapped: to make
			//positionless reads sort last (matching unmapped + the doc), return +1 here and -1 below. MASKED when reads
			//homogeneously have/lack SamLines (typical: a SAM sort has them all; a FASTQ sort has none). QUESTION for Brian.
			return -1;
		}else {
			return 1;
		}
		if(r1.id==null && r2.id==null){return 0;}
		if(r1.id==null){return -1;}
		if(r2.id==null){return 1;}
		return r1.id.compareTo(r2.id);
	}

	/** Compares two SAM lines by scaffold, then position, strand, mate position, and pair number.
	 * Unmapped reads (negative scaffold number) sort after all mapped reads. */
	public static int compareInner(SamLine a, SamLine b) {
		if(a.scafnum<0){a.setScafnum(scafMap);}
		if(b.scafnum<0){b.setScafnum(scafMap);}
		boolean aUnmapped=(a.scafnum<0);
		boolean bUnmapped=(b.scafnum<0);
		if(aUnmapped!=bUnmapped){return aUnmapped ? 1 : -1;}
		if(a.scafnum!=b.scafnum){return a.scafnum-b.scafnum;}
		if(a.pos!=b.pos){return a.pos-b.pos;}
		if(a.strand()!=b.strand()){return a.strand()-b.strand();}
		if(a.pnext!=b.pnext){return a.pnext-b.pnext;}
		if(a.pairnum()!=b.pairnum()){return a.pairnum()-b.pairnum();}
		return 0;
	}

	@Override
	public ReadComparator getComparator(boolean asc){return asc ? ascending : descending;}

	@Override
	public final int ascendingMult() {return mult;}
	@Override
	public final String name() {return "Position";}

	/** Sort direction multiplier: +1 for ascending, -1 for descending (immutable). */
	private final int mult;

	/** Ascending-order singleton. */
	public static final ReadComparatorPosition ascending=new ReadComparatorPosition(1);
	/** Descending-order singleton. */
	public static final ReadComparatorPosition descending=new ReadComparatorPosition(-1);
	/** Default (ascending) singleton; alias retained for selection/identity call sites. */
	public static final ReadComparatorPosition comparator=ascending;
	/** Maps scaffold names to numeric indices so SamLines can be ordered by scaffold; set before sorting. */
	public static ScafMap scafMap=null;

}
