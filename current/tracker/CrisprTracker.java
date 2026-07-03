package tracker;

import shared.Tools;
import structures.ByteBuilder;
import structures.Crispr;
import structures.LongList;
import structures.Range;

/**
 * Tracks statistics for CRISPR repeats and spacers, including lengths, GC content,
 * matches, mismatches, and alignment outcomes. Maintains histograms and counters for
 * downstream reporting and merging across reads or contigs.
 *
 * @author Brian Bushnell
 * @date Sept 5, 2023
 */
public class CrisprTracker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Accumulates metrics for a single CRISPR structure.
	 * Computes GC content for repeats and spacers, updates length histograms, and
	 * increments match/mismatch counters and discovery tallies.
	 *
	 * @param p CRISPR structure containing repeat ranges and alignment info
	 * @param s Sequence bytes used to compute GC content
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
		//n [tracker/CrisprTracker] comprehension: all increment() args below must be >=0 or LongList.increment AIOOBEs (a
		//n negative loc is <array.length so it skips the auto-grow resize and indexes array[neg] directly, verified
		//n structures/LongList#93/#83). rlen=r.length()>=0. slen=sstop-sstart+1=(p.b.a-1)-(p.a.b+1)+1=p.b.a-p.a.b-1>=1
		//n REQUIRES the Crispr invariant p.a.b<p.b.a (first repeat ends before second begins, i.e. a real spacer exists) —
		//n guaranteed by the repeat-spacer-repeat detector that builds Crispr; overlapping repeats would give slen<=0. ulen,
		//n tlen>=1 follow. rgci/sgci=round(gc*gcMult) with gc in [0,1] -> [0,100], never negative; NaN gc (all-N region,
		//n calcGC 0/0) -> Math.round(NaN)=0 -> bin 0, in-bounds (same safe degenerate as scalar/Scalars#001).
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
		//n studied praise: this line is a real refinement over the commented-out predecessor above. The old test called a
		//n repeat-pair "complete" whenever lengths matched; the new one ALSO demands the pair not be flush against either
		//n sequence end (p.a.a>0 && p.b.b<s.length-1). A repeat touching position 0 or s.length-1 is likely truncated by the
		//n contig boundary ("tip"), so equal observed lengths there are coincidental, not evidence of completeness. Verified
		//n against the field doc "partial or tip repeats that may be incomplete" — the code now matches the stated intent.
		partialTipRepeats+=((p.a.length()==p.b.length() && p.a.a>0 && p.b.b<s.length-1) ? 0 : 1);
		//Matches, mismatches, and copies need to be incremented manually
		//n comprehension: matchList/mismatchList ARE filled above (this add()); this note's "manually" refers to copyList,
		//n refMismatchList, refMismatchListValid — the 3 lists never touched here, incremented by the alignment caller.
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Merges statistics from another tracker into this one.
	 * Adds histogram counts and aggregates all counters.
	 * @param p Tracker to merge from
	 * @return This tracker after merging
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
	 * Appends a tab-delimited summary of CRISPR statistics to the provided buffer.
	 * Uses PalindromeTracker formatting for consistency with other trackers.
	 * @param bb Buffer to append results to
	 * @return The same buffer with statistics appended
	 */
	public ByteBuilder appendTo(ByteBuilder bb) {
		return PalindromeTracker.append(bb, 
			"#Value\trepeat\tspacer\tperiod\ttotal\trGC\tsGC\tmatch\tmismtch\tcopies\trefmm\trefmmv", 
			lists, histmax);
	}
	
	/** Returns a string representation of the tracked CRISPR statistics. */
	@Override
	public String toString() {
		return appendTo(new ByteBuilder()).toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts of A/C/G/T/N bases observed across processed CRISPR regions. */
	//n [tracker/CrisprTracker#002] LOW/dead: acgtn is never read or written anywhere in this class (no callers set it —
	//n it is private). Dead field; harmless. Leaving in place (may be a stub for a planned per-base tally). Flag, don't remove.
	private final int[] acgtn=new int[5];

	/** Total number of CRISPR structures discovered. */
	public long crisprsFound=0;
	/** Number of reads that contained at least one CRISPR structure. */
	public long readsWithCrisprs=0;
	/** Number of reads trimmed based on CRISPR consensus analysis. */
	public long trimmedByConsensus=0;
	/** Count of CRISPRs with partial or tip repeats that may be incomplete. */
	public long partialTipRepeats=0;
	/** Number of CRISPRs successfully aligned to a reference. */
	public long alignedToRef=0;
	/** Number of CRISPRs modified after reference alignment. */
	public long modifiedByRef=0;
	/** Number of CRISPR alignments that failed. */
	public long failedAlignment=0;
	/** Total number of alignment attempts performed. */
	public long alignments=0;
	/** Number of times an alignment was requested. */
	public long alignmentRequested=0;
	
	/** Histogram of repeat lengths. */
	public LongList rlenList=new LongList();
	/** Histogram of spacer lengths. */
	public LongList slenList=new LongList();
	/** Histogram of combined spacer+repeat lengths. */
	public LongList ulenList=new LongList();
	/** Histogram of total CRISPR structure lengths. */
	public LongList tlenList=new LongList();
	/** Histogram of match counts between repeats. */
	public LongList matchList=new LongList();
	/** Histogram of mismatch counts between repeats. */
	public LongList mismatchList=new LongList();
	/** Histogram of repeat copy counts. */
	public LongList copyList=new LongList();
	//n [tracker/CrisprTracker#001] LOW/DOC: "2% increments" + capacity 51 are STALE. gcMult=100 (see field below, doc says
	//n "100 = 1% steps") so rgci/sgci=round(gc*100) span bins 0..100 (101 values), i.e. 1% increments, not 2%. Not a crash:
	//n LongList(51) is only a growth hint and increment() auto-grows (resize) at the first bin>50, which fires on any GC>50%
	//n (extremely common). Effect is a one-time array regrow + a misleading comment. Fix = doc "1% increments" and, if desired,
	//n size these new LongList(101). Verified 1% via gcMult doc + Math.round; no correctness impact.
	/** Histogram of repeat GC content in 1% increments (bins 0..100). */
	public LongList rgcList=new LongList(51);
	/** Histogram of spacer GC content in 1% increments (bins 0..100). */
	public LongList sgcList=new LongList(51);
	/** Histogram of reference mismatch counts for all alignments. */
	public LongList refMismatchList=new LongList(51);
	/** Histogram of reference mismatches for valid alignments only. */
	public LongList refMismatchListValid=new LongList(51);
	
	/** Array of all histogram lists for batch operations and merging. */
	public final LongList[] lists={rlenList, slenList, ulenList, tlenList, 
			rgcList, sgcList, matchList, mismatchList, 
			copyList, refMismatchList, refMismatchListValid};
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum value for histogram binning. */
	public static int histmax=150;
	/** Multiplier for converting GC fraction to histogram bins (100 = 1% steps). */
	public static int gcMult=100;
	
}
