package assemble;

import structures.IntList;

/** One-pass geometry for bounded, supported-flank low-depth regions.
 * Returns candidate intervals, never an error verdict or an accepted edit.
 * Counts are indexed by original kmer START. Bases and the depth list are borrowed
 * without copying and must remain stable until scanning ends or reset replaces them.
 * Mutable scan/result state belongs to one worker; concurrent calls are unsupported.
 * @author Fischl */
public final class LocalEditTroughLocator{
	/** Fixes K>=5 and the positive depth threshold separating low windows from flanks. */
	public LocalEditTroughLocator(final int k_, final int minCount_){
		if(k_<5 || minCount_<1){throw new IllegalArgumentException("Locator requires K>=5 and positive supported-depth threshold.");}
		k=k_; minCount=minCount_;
	}
	/** Begins a dense scan; -1 denotes absent depth, and sparse sentinels are not allowed. */
	public void reset(final byte[] bases_, final IntList counts_){
		reset(bases_, counts_, false);
	}
	/**
	 * Borrows exactly max(0,readLength-K+1) depths and resets scan/counters on success.
	 * Invalid arguments discard the previous scan/result but do not reset its counters.
	 * Sparse profiles distinguish UNKNOWN from INVALID_N: neither is low evidence.
	 * A discovered low region must have measured flanks or N/edge barriers; an UNKNOWN
	 * immediate flank is a producer-contract failure, not evidence to silently skip.
	 */
	public void reset(final byte[] bases_, final IntList counts_, final boolean sparse_){
		bases=null; counts=null; clearResult();
		if(bases_==null || counts_==null || counts_.size!=Math.max(0, bases_.length-k+1)){
			throw new IllegalArgumentException("Locator needs exactly max(0,readLength-K+1) original-window depths.");
		}
		bases=bases_; counts=counts_; sparse=sparse_; cursor=checkedBases=0; lastUndefined=-1;
		lowRegions=skippedEdge=skippedWide=skippedUndefined=0; clearResult();
	}
	/**
	 * Advances to the next interior low-depth run of at most K windows with defined
	 * uppercase A/C/G/T from the left supporting window through the right one.
	 * Skips edge, wide and undefined regions; never rescans earlier positions.
	 * A true result publishes half-open window, base and insertion-boundary ranges.
	 * A false result means exhaustion and clears all ranges to -1. Read result fields
	 * before the next call/reset; enumeration alone does not validate an edit.
	 */
	public boolean next(){
		if(counts==null){throw new IllegalStateException("Reset the locator before scanning.");}
		clearResult();
		while(cursor<counts.size){
			if(!low(cursor)){cursor++; continue;}
			final int a=cursor++;
			while(cursor<counts.size && low(cursor)){cursor++;}
			final int end=cursor;
			lowRegions++;
			if(a==0 || end==counts.size){skippedEdge++; continue;}
			if(sparse){
				if(counts.get(a-1)==LocalEditDepthProfile.UNKNOWN || counts.get(end)==LocalEditDepthProfile.UNKNOWN){
					throw new IllegalStateException("Sparse low regions require measured flanks or N/edge barriers.");
				}
				if(counts.get(a-1)==LocalEditDepthProfile.INVALID_N || counts.get(end)==LocalEditDepthProfile.INVALID_N){skippedUndefined++; continue;}
			}
			if(end-a>k){skippedWide++; continue;}
			// Include both supporting flank windows. Advance a shared validity
			// cursor rather than rescanning their union for every trough.
			final int from=a-1, to=end+k;
			assert(from>=0 && to<=bases.length) : "Interior trough must have a complete high-depth window on each side.";
			while(checkedBases<to){
				final byte b=bases[checkedBases];
				if(b!='A' && b!='C' && b!='G' && b!='T'){lastUndefined=checkedBases;}
				checkedBases++;
			}
			if(lastUndefined>=from){skippedUndefined++; continue;}
			depthStart=a; depthEnd=end;
			// Intersect [j,j+k) over all low-window starts j in [a,end).
			baseStart=end-1; baseEnd=a+k;
			// Insertion boundary g crosses a window iff j<g<j+k.
			gapStart=end; gapEnd=a+k;
			assert(baseStart<baseEnd && gapStart<=gapEnd && baseStart>=0 && baseEnd<=bases.length) :
				"A trough of width<=K yields a nonempty base intersection and possibly empty boundary intersection.";
			return true;
		}
		return false;
	}
	/** Classifies one visited depth below minCount, including absent=-1. Rejects other
	 * negative values except the two recognized sparse sentinels, which are barriers. */
	private boolean low(final int position){
		final int depth=counts.get(position);
		if(sparse && (depth==LocalEditDepthProfile.UNKNOWN || depth==LocalEditDepthProfile.INVALID_N)){return false;}
		if(depth<-1){throw new IllegalArgumentException("Depth must be nonnegative or absent=-1: "+depth);}
		return depth<minCount;
	}
	/** Invalidates published ranges without rewinding the scan or clearing counters. */
	private void clearResult(){depthStart=depthEnd=baseStart=baseEnd=gapStart=gapEnd=-1;}
	/** Half-open low-window starts and their common base/insertion-boundary intersections.
	 * Gaps are integer boundaries BEFORE a base; a width-K run has an empty gap range. */
	public int depthStart, depthEnd, baseStart, baseEnd, gapStart, gapEnd;
	/** Low runs encountered, including rejected runs, and skips by their first rejection reason. */
	public int lowRegions, skippedEdge, skippedWide, skippedUndefined;
	/** Borrowed original sequence and its start-indexed depths; null after a failed reset. */
	private byte[] bases;
	private IntList counts;
	/** Next depth index, next unchecked base, and most recent undefined base index. */
	private int cursor, checkedBases, lastUndefined;
	/** Enables the depth-profile sentinel contract; does not fill missing measurements. */
	private boolean sparse;
	/** Fixed word length and supported-flank depth threshold. */
	private final int k, minCount;
}
