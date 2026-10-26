package assemble;

import map.LongHashMap;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/** Reusable multiplicities of unfiltered Hollow keys wholly inside one original seed.
 * Uses the same two fingerprints as the global pair table, including strand canonicalization.
 * Primary-hash collisions retain separate secondary fingerprints. Never edits a shared table.
 * @author Brian Bushnell, Fischl */
public final class HollowSeedCounts {

	/** Binds the ledger to the same full span as the immutable global table. */
	public HollowSeedCounts(final int k_){
		if(k_<=64 || !Kmer.PACKED){throw new IllegalArgumentException("Hollow seed counts require packed K>64.");}
		k=k_;
		word=new Kmer(k, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Counting           ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces the ledger with all defined windows inside the specified original prefix.
	 * Qualities and loader transformations are deliberately outside this diagnostic contract. */
	public void reset(final byte[] bases, final int length){
		if(bases==null || length<0 || length>bases.length){throw new IllegalArgumentException("Invalid original seed prefix length: "+length);}
		assert(word.hollow && word.kbig==k) : "Local and global pair hashes must encode the same Hollow span.";
		clear();
		word.clear();
		for(int i=0; i<length; i++){
			word.addRight(bases[i]);
			if(word.len>=k){add(word.xor(), word.xor2());}
		}
	}

	/** Clears only local state; buffers are retained for the next seed. */
	public void clear(){
		heads.clear();
		fingerprints.clear();
		next.clear();
		counts.clear();
	}

	/** Adds one pair-table key; exposed for deterministic primary-collision fixtures. */
	public void add(final long primary, final long secondary){
		final int head=heads.get(primary);
		int previous=-1;
		for(int index=head<0 ? -1 : head-1; index>=0; index=next.get(index)){
			if(fingerprints.get(index)==secondary){
				assert(counts.get(index)<Integer.MAX_VALUE) : "A seed multiplicity must not overflow and masquerade as external support.";
				counts.increment(index);
				return;
			}
			previous=index;
		}
		final int index=fingerprints.size;
		fingerprints.add(secondary);
		next.add(-1);
		counts.add(1);
		//LongHashMap.put does not overwrite. Append collisions to the existing chain.
		if(previous<0){heads.put(primary, index+1);}
		else{next.set(previous, index);}
		assert(fingerprints.size==next.size && next.size==counts.size) : "Every secondary fingerprint needs exactly one chain link and count.";
	}

	/** Reads a complete nominal or alternate-span query using its ordinary Hollow hashes. */
	public int count(final Kmer query){
		assert(query!=null && query.hollow && query.kbig==k && query.len>=k) :
				"Subtracting a different geometry or incomplete query would invent an external count.";
		return count(query.xor(), query.xor2());
	}

	/** Looks up both fingerprints; sharing only a placement hash is not self support. */
	public int count(final long primary, final long secondary){
		final int head=heads.get(primary);
		for(int index=head<0 ? -1 : head-1; index>=0; index=next.get(index)){
			assert(index<fingerprints.size) : "Collision-chain links must refer to this seed's active ledger, not an earlier seed.";
			if(fingerprints.get(index)==secondary){return counts.get(index);}
		}
		return 0;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int k;
	private final Kmer word;
	private final LongHashMap heads=new LongHashMap(128);
	private final LongList fingerprints=new LongList(128);
	private final IntList next=new IntList(128), counts=new IntList(128);
}
