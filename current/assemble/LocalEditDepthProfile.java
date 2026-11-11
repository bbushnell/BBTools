package assemble;

import structures.IntList;
import ukmer.Kmer;

/**
 * Worker-local dense or sampled depths indexed by original kmer-start position.
 * Sparse low-depth components found by sampling are expanded to measured depth>=3,
 * N barriers or edges; a low region missed by every sample need not be discovered.
 * Borrows a stable count provider and reuses mutable keys/profile storage. Calls on
 * one instance must be serial. Input bases and table counts are not modified.
 * @author Fischl */
public final class LocalEditDepthProfile{
	/** Borrows counts for an exactly representable K>=5 under the current Kmer layout. */
	public LocalEditDepthProfile(final int k_, final HomopolymerIndelProposal.CountLookup lookup_){
		if(k_<5 || lookup_==null){throw new IllegalArgumentException("Profile requires K>=5 and immutable counts.");}
		k=k_; lookup=lookup_; key=new Kmer(k);
		if(key.kbig!=k){throw new IllegalArgumentException("Profile K must equal the table K.");}
	}
	/**
	 * Replaces depths and resets the lookup-attempt counter for uppercase A/C/G/T/N.
	 * Short reads have no windows. Dense scans map N-containing windows to zero;
	 * sparse scans retain separate INVALID_N and unmeasured UNKNOWN states.
	 * Window start 0 and the last window are sampled, with requestedStride between samples;
	 * K or window counts no larger than that stride force a dense scan instead.
	 * Sampled depths 0–2 trigger refinement, matching the corrector's depth-3 locator.
	 * The reused depth list remains valid until the next fill; an exception may leave
	 * a partial profile, which must not be consumed as a completed scan.
	 */
	public void fill(final byte[] bases, final int requestedStride){
		if(bases==null || requestedStride<1){throw new IllegalArgumentException("Profile needs bases and a positive stride.");}
		queries=0; depths.clear(); key.clearFast(); keyStart=-1;
		final int windows=Math.max(0, bases.length-k+1);
		// Small K/short domains use exact scanning rather than dropping their
		// entire error signature between samples. Dense consumers see no sentinels.
		stride=k<=requestedStride || windows<=requestedStride ? 1 : requestedStride;
		sparse=stride>1;
		long nextSample=0;
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			if(b=='N'){key.clearFast();}else{
				if(b!='A' && b!='C' && b!='G' && b!='T'){throw new IllegalArgumentException("Profile accepts only A/C/G/T/N.");}
				key.addRight(b);
			}
			if(i<k-1){continue;}
			final int start=i-k+1; keyStart=start;
			final boolean sample=start==nextSample || start==windows-1;
			if(start==nextSample){nextSample+=stride;}
			depths.add(key.len()<k ? (sparse ? INVALID_N : 0) : sample ? count() : UNKNOWN);
		}
		assert(depths.size==windows) : "Each depth entry addresses one original kmer-start position.";
		if(!sparse){return;}
		for(int p=0; p<windows; p++){
			final int depth=depths.get(p);
			if(depth<0 || depth>=3){continue;}
			// Finish this sampled component once. depthAt reuses all previously
			// sampled/refined values; only UNKNOWN entries cause new queries.
			for(int left=p-1; left>=0; left--){
				final int d=depthAt(bases, left); if(d==INVALID_N || d>=3){break;}
			}
			int end=p+1;
			for(; end<windows; end++){
				final int d=depthAt(bases, end); if(d==INVALID_N || d>=3){break;}
			}
			p=end-1;
		}
	}
	/** Returns a cached state or measures one valid UNKNOWN window. Adjacent requests
	 * roll the scratch key left/right; other requests rebuild it from the same bases. */
	private int depthAt(final byte[] bases, final int start){
		assert(start>=0 && start<depths.size) : "Refinement stays within the original kmer-start domain.";
		final int existing=depths.get(start); if(existing!=UNKNOWN){return existing;}
		if(start==keyStart+1 && key.len()>=k){key.addRight(bases[start+k-1]);}else if(start==keyStart-1 && key.len()>=k){key.addLeft(bases[start]);}else{
			key.clearFast(); for(int i=start; i<start+k; i++){key.addRight(bases[i]);}
		}
		keyStart=start;
		assert(key.len()>=k) : "Only valid, unsampled windows are UNKNOWN; N-containing windows were marked INVALID_N during the initial pass.";
		final int depth=count(); depths.set(start, depth); return depth;
	}
	/** Counts one complete scratch key, normalizing the provider's absent=-1 to zero.
	 * The counter includes a lookup that throws or returns a forbidden depth below-1. */
	private int count(){
		assert(key.len()>=k) : "Count-table queries require a complete valid kmer.";
		queries++; final int depth=lookup.count(key);
		if(depth<-1){throw new IllegalStateException("Invalid count depth: "+depth);}
		return Math.max(0, depth);
	}
	/** Sparse states: valid but unmeasured, and invalid because a window contains N. */
	public static final int UNKNOWN=-2, INVALID_N=-3;
	/** Reused window-depth output, addressed by original start rather than end position. */
	public final IntList depths=new IntList();
	/** Count-provider invocations attempted during the current fill. */
	public long queries;
	/** Whether the effective stride permits unmeasured entries. */
	public boolean sparse;
	/** Effective sample spacing; may be one despite a larger requested stride. */
	public int stride;
	/** Start associated with the scratch key, possibly incomplete after N; -1 before scanning. */
	private int keyStart;
	/** Exact window length, shared with the immutable count table. */
	private final int k;
	/** Mutable scan/refinement key; never retained by the borrowed lookup. */
	private final Kmer key;
	/** Borrowed count provider, kept stable by the caller for the scan. */
	private final HomopolymerIndelProposal.CountLookup lookup;
}
