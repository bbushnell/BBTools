package assemble;

import structures.IntList;
import ukmer.Kmer;

/** Read-only complete profiles of a virtual single-base deletion.
 * One instance/scratch per worker; no correction or truth classification.
 * @author Fischl */
final class HomopolymerDeletionEvidence {

	HomopolymerDeletionEvidence(final int k_,final HomopolymerIndelProposal.CountLookup lookup_){
		if(k_<5 || lookup_==null){throw new IllegalArgumentException("Deletion evidence needs exact k>=5 and an immutable lookup.");}
		k=k_;lookup=lookup_;key=new Kmer(k);
		if(key.kbig!=k){throw new IllegalArgumentException("Deletion profile k rounded from "+k+" to "+key.kbig);}
	}
	/** Materializes no sequence: skip exactly deletedPos within original[from,to).
	 * Context endpoints stay fixed even for an immediately adjacent run flank. */
	void profile(final byte[] bases,final int from,final int to,final int deletedPos,final IntList depths){
		if(bases==null || depths==null || from<0 || to> bases.length || to-from-1<k || deletedPos<from || deletedPos>=to){
			throw new IllegalArgumentException("Deletion profile needs a complete in-bounds context and an included deletion position.");
		}
		for(int i=from;i<to;i++){if(!defined(bases[i])){throw new IllegalArgumentException("Undefined deletion-context base at "+i);}}
		depths.clear();key.clearFast();minimum=Integer.MAX_VALUE;maximum=0;
		for(int i=from;i<to;i++){
			if(i==deletedPos){continue;}key.addRight(bases[i]);
			if(key.len()<k){continue;}
			final int raw=lookup.count(key);
			if(raw< -1){throw new IllegalStateException("Invalid deletion count "+raw+"; only -1 means absent.");}
			final int depth=Math.max(0,raw);depths.add(depth);minimum=Math.min(minimum,depth);maximum=Math.max(maximum,depth);
		}
		assert(depths.size==to-from-k && depths.size>0) : "One deleted base leaves exactly originalLength-k complete kmer windows.";
	}
	private static boolean defined(final byte b){return b=='A' || b=='C' || b=='G' || b=='T';}
	int minimum,maximum;
	private final int k;
	private final Kmer key;
	private final HomopolymerIndelProposal.CountLookup lookup;
}
