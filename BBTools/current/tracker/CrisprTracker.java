package tracker;

import shared.Tools;
import structures.ByteBuilder;
import structures.Crispr;
import structures.LongList;
import structures.Range;

/**
 * Tracks crispr stats.
 * 
 * @author Brian Bushnell
 * @date Sept 5, 2023
 *
 */
public class CrisprTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void add(final Crispr p, final byte[] s) {
		Range r=(p.b.length()>p.a.length() ? p.b : p.a);
		int sstart=p.a.b+1, sstop=p.b.a-1;
		float rgc=Tools.calcGC(s, r.a, r.b);
		float sgc=Tools.calcGC(s, sstart, sstop);
		int rgci=(int)(Math.round(rgc*gcMult));
		int sgci=(int)(Math.round(sgc*gcMult));
		int rlen=r.length();
		int slen=sstop-sstart+1;
		int ulen=rlen+slen;
		int tlen=p.b.b-p.a.a+1;

		rlenList.increment(rlen);
		slenList.increment(slen);
		ulenList.increment(ulen);
		tlenList.increment(tlen);
		rgcList.increment(rgci);
		sgcList.increment(sgci);
		matchList.increment(p.matches);
		mismatchList.increment(p.mismatches);
		crisprsFound++;
//		partialTipRepeats+=((p.a.length()==p.b.length()) ? 0 : 1);
		partialTipRepeats+=((p.a.length()==p.b.length() && p.a.a>0 && p.b.b<s.length-1) ? 0 : 1);
		//Matches, mismatches, and copies need to be incremented manually
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public CrisprTracker add(CrisprTracker p) {
		for(int i=0; i<lists.length; i++) {
			lists[i].incrementBy(p.lists[i]);
		}
		crisprsFound+=p.crisprsFound;
		readsWithCrisprs+=p.readsWithCrisprs;
		trimmedByConsensus+=p.trimmedByConsensus;
		partialTipRepeats+=p.partialTipRepeats;
		modifiedByRef+=p.modifiedByRef;
		alignedToRef+=p.alignedToRef;
		failedAlignment+=p.failedAlignment;
		alignments+=p.alignments;
		alignmentRequested+=p.alignmentRequested;
		return this;
	}

	
	public ByteBuilder appendTo(ByteBuilder bb) {
		return PalindromeTracker.append(bb, 
			"#Value\trepeat\tspacer\tperiod\ttotal\trGC\tsGC\tmatch\tmismtch\tcopies\trefmm\trefmmv", 
			lists, histmax);
	}
	
	@Override
	public String toString() {
		return appendTo(new ByteBuilder()).toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] acgtn=new int[5];

	public long crisprsFound=0;
	public long readsWithCrisprs=0;
	public long trimmedByConsensus=0;
	public long partialTipRepeats=0;
	public long alignedToRef=0;
	public long modifiedByRef=0;
	public long failedAlignment=0;
	public long alignments=0;
	public long alignmentRequested=0;
	
	/** Repeat length */
	public LongList rlenList=new LongList();
	/** Spacer length */
	public LongList slenList=new LongList();
	/** Length of a spacer+repeat */
	public LongList ulenList=new LongList();
	/** Total length */
	public LongList tlenList=new LongList();
	/** Number of matches */
	public LongList matchList=new LongList();
	/** Number of mismatches between repeats */
	public LongList mismatchList=new LongList();
	/** Copy count of this repeat */
	public LongList copyList=new LongList();
	/** Repeat gc, in 2% increments */
	public LongList rgcList=new LongList(51);
	/** Spacer gc, in 2% increments */
	public LongList sgcList=new LongList(51);
	/** Number of mismatches with reference */
	public LongList refMismatchList=new LongList(51);
	/** Number of mismatches with reference for valid alignments */
	public LongList refMismatchListValid=new LongList(51);
	
	public final LongList[] lists={rlenList, slenList, ulenList, tlenList, 
			rgcList, sgcList, matchList, mismatchList, 
			copyList, refMismatchList, refMismatchListValid};
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	public static int histmax=150;
	public static int gcMult=100;
	
}
