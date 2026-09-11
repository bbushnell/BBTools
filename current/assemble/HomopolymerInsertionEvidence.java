package assemble;

import structures.IntList;
import ukmer.Kmer;

/** Read-only virtual single-base insertion profiles; not an acceptance rule.
 * One instance and its scratch belong to one worker.
 * @author Fischl
 */
final class HomopolymerInsertionEvidence {
	HomopolymerInsertionEvidence(final int k_, final HomopolymerIndelProposal.CountLookup lookup_){
		if(k_<5 || lookup_==null){throw new IllegalArgumentException("Insertion evidence needs an exact k>=5 and immutable count lookup.");}
		k=k_;lookup=lookup_;key=new Kmer(k);
		if(key.kbig!=k){throw new IllegalArgumentException("Requested k is rounded by the current kmer layout: "+k+" -> "+key.kbig);}
	}

	/** Every kmer of original [from,to) with base inserted at original gap.
	 * Counts are normalized like HomopolymerIndelProposal: absent -1 becomes0.
	 * Output order follows the oriented candidate, with no per-window allocation. */
	void profile(final byte[] bases, final int from, final int to, final int gap,
			final byte base, final IntList depths){
		if(bases==null || depths==null || from<0 || to<from || to>bases.length ||
				gap<from || gap>to || to-from<k || !defined(base)){
			throw new IllegalArgumentException("Insertion profile requires a defined base, complete context and an in-context gap.");
		}
		for(int i=from;i<to;i++){
			if(!defined(bases[i])){throw new IllegalArgumentException("Undefined insertion context base at "+i);}
		}
		depths.clear();key.clearFast();minimum=Integer.MAX_VALUE;maximum=0;
		final int relative=gap-from,size=to-from+1;
		for(int i=0;i<size;i++){
			final byte b=i==relative ? base : bases[from+i-(i>relative ? 1 : 0)];
			key.addRight(b);if(i<k-1){continue;}
			final int raw=lookup.count(key);
			if(raw< -1){throw new IllegalStateException("Count lookup returned invalid depth "+raw+"; only -1 means absent.");}
			final int depth=Math.max(0,raw);depths.add(depth);
			minimum=Math.min(minimum,depth);maximum=Math.max(maximum,depth);
		}
		assert(depths.size==size-k+1 && depths.size>0) : "Every complete virtual candidate kmer must contribute exactly one count.";
	}
	private static boolean defined(final byte b){return b=='A' || b=='C' || b=='G' || b=='T';}
	int minimum,maximum;
	private final int k;
	private final Kmer key;
	private final HomopolymerIndelProposal.CountLookup lookup;
}
