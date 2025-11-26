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
	
	/**
	 * Adds CRISPR statistics from a single CRISPR structure to the tracker.
	 * Calculates GC content, lengths, and increments appropriate counters and histograms.
	 * @param p The CRISPR structure containing repeat ranges and alignment data
	 * @param s The sequence bytes for GC content calculation
	 */
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
	
	/**
	 * Merges statistics from another CrisprTracker into this one.
	 * Combines all counters, histograms, and metrics from both trackers.
	 * @param p The CrisprTracker to merge into this one
	 * @return This CrisprTracker after merging
	 */
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

	
	/**
	 * Appends formatted CRISPR statistics to a ByteBuilder.
	 * Uses PalindromeTracker formatting for consistency with other trackers.
	 * @param bb The ByteBuilder to append statistics to
	 * @return The ByteBuilder with appended statistics
	 */
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

	/** Array for tracking nucleotide composition (A, C, G, T, N) */
	private final int[] acgtn=new int[5];

	/** Total number of CRISPR sequences found */
	public long crisprsFound=0;
	/** Number of reads containing CRISPR sequences */
	public long readsWithCrisprs=0;
	/** Number of sequences trimmed based on consensus analysis */
	public long trimmedByConsensus=0;
	/** Number of repeats that are partial or at sequence tips */
	public long partialTipRepeats=0;
	/** Number of sequences successfully aligned to reference */
	public long alignedToRef=0;
	/** Number of sequences modified based on reference alignment */
	public long modifiedByRef=0;
	/** Number of sequences that failed alignment to reference */
	public long failedAlignment=0;
	/** Total number of alignments performed */
	public long alignments=0;
	/** Number of alignments that were requested */
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
	
	/** Array containing all histogram lists for batch operations */
	public final LongList[] lists={rlenList, slenList, ulenList, tlenList, 
			rgcList, sgcList, matchList, mismatchList, 
			copyList, refMismatchList, refMismatchListValid};
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum value for histogram binning */
	public static int histmax=150;
	/**
	 * Multiplier for converting GC percentages to histogram bins (100 = 1% resolution)
	 */
	public static int gcMult=100;
	
}
