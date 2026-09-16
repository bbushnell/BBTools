package assemble;

import structures.IntList;
import ukmer.Kmer;

/** Worker-local dense or sampled depths with explicit unknown/invalid states.
 * Sparse components are expanded to measured high flanks, N barriers or edges.
 * Never changes reads or the immutable count table.
 * @author Fischl */
final class LocalEditDepthProfile {
	LocalEditDepthProfile(final int k_,final HomopolymerIndelProposal.CountLookup lookup_){
		if(k_<5 || lookup_==null){throw new IllegalArgumentException("Profile requires K>=5 and immutable counts.");}
		k=k_;lookup=lookup_;key=new Kmer(k);
		if(key.kbig!=k){throw new IllegalArgumentException("Profile K must equal the table K.");}
	}
	void fill(final byte[] bases,final int requestedStride){
		if(bases==null || requestedStride<1){throw new IllegalArgumentException("Profile needs bases and a positive stride.");}
		queries=0;depths.clear();key.clearFast();keyStart=-1;
		final int windows=Math.max(0,bases.length-k+1);
		// Small K/short domains use exact scanning rather than dropping their
		// entire error signature between samples. Dense consumers see no sentinels.
		stride=k<=requestedStride || windows<=requestedStride ? 1 : requestedStride;
		sparse=stride>1;
		long nextSample=0;
		for(int i=0;i<bases.length;i++){
			final byte b=bases[i];
			if(b=='N'){key.clearFast();}
			else{
				if(b!='A' && b!='C' && b!='G' && b!='T'){throw new IllegalArgumentException("Profile accepts only A/C/G/T/N.");}
				key.addRight(b);
			}
			if(i<k-1){continue;}
			final int start=i-k+1;keyStart=start;
			final boolean sample=start==nextSample || start==windows-1;
			if(start==nextSample){nextSample+=stride;}
			depths.add(key.len()<k ? (sparse ? INVALID_N : 0) : sample ? count() : UNKNOWN);
		}
		assert(depths.size==windows) : "Each depth entry addresses one original kmer-start position.";
		if(!sparse){return;}
		for(int p=0;p<windows;p++){
			final int depth=depths.get(p);
			if(depth<0 || depth>=3){continue;}
			// Finish this sampled component once. depthAt reuses all previously
			// sampled/refined values; only UNKNOWN entries cause new queries.
			for(int left=p-1;left>=0;left--){
				final int d=depthAt(bases,left);if(d==INVALID_N || d>=3){break;}
			}
			int end=p+1;
			for(;end<windows;end++){
				final int d=depthAt(bases,end);if(d==INVALID_N || d>=3){break;}
			}
			p=end-1;
		}
	}
	private int depthAt(final byte[] bases,final int start){
		assert(start>=0 && start<depths.size) : "Refinement stays within the original kmer-start domain.";
		final int existing=depths.get(start);if(existing!=UNKNOWN){return existing;}
		if(start==keyStart+1 && key.len()>=k){key.addRight(bases[start+k-1]);}
		else if(start==keyStart-1 && key.len()>=k){key.addLeft(bases[start]);}
		else{
			key.clearFast();for(int i=start;i<start+k;i++){key.addRight(bases[i]);}
		}
		keyStart=start;
		assert(key.len()>=k) : "Only valid, unsampled windows are UNKNOWN; N-containing windows were marked INVALID_N during the initial pass.";
		final int depth=count();depths.set(start,depth);return depth;
	}
	private int count(){
		assert(key.len()>=k) : "Count-table queries require a complete valid kmer.";
		queries++;final int depth=lookup.count(key);
		if(depth<-1){throw new IllegalStateException("Invalid count depth: "+depth);}
		return Math.max(0,depth);
	}
	static final int UNKNOWN=-2,INVALID_N=-3;
	final IntList depths=new IntList();
	long queries;
	boolean sparse;
	int stride;
	private int keyStart;
	private final int k;
	private final Kmer key;
	private final HomopolymerIndelProposal.CountLookup lookup;
}
