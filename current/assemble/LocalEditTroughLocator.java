package assemble;

import structures.IntList;

/** One-pass geometry for bounded, supported-flank low-depth regions.
 * Returns candidate intervals, never an error verdict or an accepted edit.
 * Counts are indexed by original kmer START; read/profile remain immutable.
 * @author Fischl */
final class LocalEditTroughLocator {
	LocalEditTroughLocator(final int k_,final int minCount_){
		if(k_<5 || minCount_<1){throw new IllegalArgumentException("Locator requires K>=5 and positive supported-depth threshold.");}
		k=k_;minCount=minCount_;
	}
	void reset(final byte[] bases_,final IntList counts_){
		reset(bases_,counts_,false);
	}
	/** Sparse profiles distinguish unmeasured windows from N-invalid windows. */
	void reset(final byte[] bases_,final IntList counts_,final boolean sparse_){
		bases=null;counts=null;clearResult();
		if(bases_==null || counts_==null || counts_.size!=Math.max(0,bases_.length-k+1)){
			throw new IllegalArgumentException("Locator needs exactly max(0,readLength-K+1) original-window depths.");
		}
		bases=bases_;counts=counts_;sparse=sparse_;cursor=checkedBases=0;lastUndefined=-1;
		lowRegions=skippedEdge=skippedWide=skippedUndefined=0;clearResult();
	}
	boolean next(){
		if(counts==null){throw new IllegalStateException("Reset the locator before scanning.");}
		clearResult();
		while(cursor<counts.size){
			if(!low(cursor)){cursor++;continue;}
			final int a=cursor++;
			while(cursor<counts.size && low(cursor)){cursor++;}
			final int end=cursor;
			lowRegions++;
			if(a==0 || end==counts.size){skippedEdge++;continue;}
			if(sparse){
				if(counts.get(a-1)==LocalEditDepthProfile.UNKNOWN || counts.get(end)==LocalEditDepthProfile.UNKNOWN){
					throw new IllegalStateException("Sparse low regions require measured flanks or N/edge barriers.");
				}
				if(counts.get(a-1)==LocalEditDepthProfile.INVALID_N || counts.get(end)==LocalEditDepthProfile.INVALID_N){skippedUndefined++;continue;}
			}
			if(end-a>k){skippedWide++;continue;}
			// Include both supporting flank windows. Advance a shared validity
			// cursor rather than rescanning their union for every trough.
			final int from=a-1,to=end+k;
			assert(from>=0 && to<=bases.length) : "Interior trough must have a complete high-depth window on each side.";
			while(checkedBases<to){
				final byte b=bases[checkedBases];
				if(b!='A' && b!='C' && b!='G' && b!='T'){lastUndefined=checkedBases;}
				checkedBases++;
			}
			if(lastUndefined>=from){skippedUndefined++;continue;}
			depthStart=a;depthEnd=end;
			// Intersect [j,j+k) over all low-window starts j in [a,end).
			baseStart=end-1;baseEnd=a+k;
			// Insertion boundary g crosses a window iff j<g<j+k.
			gapStart=end;gapEnd=a+k;
			assert(baseStart<baseEnd && gapStart<=gapEnd && baseStart>=0 && baseEnd<=bases.length) :
				"A trough of width<=K yields a nonempty base intersection and possibly empty boundary intersection.";
			return true;
		}
		return false;
	}
	private boolean low(final int position){
		final int depth=counts.get(position);
		if(sparse && (depth==LocalEditDepthProfile.UNKNOWN || depth==LocalEditDepthProfile.INVALID_N)){return false;}
		if(depth<-1){throw new IllegalArgumentException("Depth must be nonnegative or absent=-1: "+depth);}
		return depth<minCount;
	}
	private void clearResult(){depthStart=depthEnd=baseStart=baseEnd=gapStart=gapEnd=-1;}
	/** All intervals are half-open; gaps are integer boundaries BEFORE a base. */
	int depthStart,depthEnd,baseStart,baseEnd,gapStart,gapEnd;
	int lowRegions,skippedEdge,skippedWide,skippedUndefined;
	private byte[] bases;
	private IntList counts;
	private int cursor,checkedBases,lastUndefined;
	private boolean sparse;
	private final int k,minCount;
}
