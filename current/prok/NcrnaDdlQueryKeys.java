package prok;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import ddl.DDLIndexSmall;
import kmer.AbstractKmerTable;
import map.IntHashSet;
import shared.Tools;
import structures.IntList;

/** Raw-query keys for the member-voted ncRNA ranking/caller experiment.
 * CardinalityTracker shares the DNA loop with member construction. This prok-local
 * collector uses Tools.hash64shift and DynamicDemiLog.encodeKey, fixed at
 * exponent4/seed12345, before native floor/bucket-max/blacklist suppression.
 * The insertion API cannot return displaced raw keys; encoding is shared with
 * that insertion path rather than copied here.
 *
 * A worker owns and reuses one instance for one k/width. Each build replaces its
 * previous query. All nonzero raw keys survive, including keys displaced by a
 * higher key in the same bucket; repeated occurrences need no extra storage.
 * The optional native sketch is borrowed after its owner has cleared it; this
 * class does not clear or retain that sketch. Raw-only use neither reads nor
 * changes global DDL exponent/blacklist state. A supplied native sketch still
 * requires the matching exponent4/no-blacklist construction environment.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlQueryKeys implements CardinalityTracker.KmerConsumer, NcrnaDdlLayoutVotes.KeyPrefix{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlQueryKeys(final int buckets_, final int k_){
		if(buckets_<32 || buckets_>2048 || Integer.bitCount(buckets_)!=1 || k_<1 || k_>31){
			throw new IllegalArgumentException("Raw-key experiment requires a power-of-two width32..2048 and k1..31");
		}
		buckets=buckets_;
		k=k_;
		bucketMask=buckets-1;
		keys=new IntHashSet(256);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Query Keys          ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows an entire captured window, with no slicing or base conversion.
	 * Preserves the native loop's ambiguity handling; that loop does not support
	 * general IUPAC symbols. No qualities, probability filtering or blacklist.
	 */
	public void build(final byte[] bases){
		build(bases, null, true, false);
	}

	/** Replace raw keys and optionally populate an already-cleared native sketch
	 * from the same canonical/plain or forward/plain/masked words. */
	public void build(final byte[] bases, final DynamicDemiLog sketch, final boolean canonical, final boolean middle){
		if(bases==null || middle && (canonical || k<3)){
			throw new IllegalArgumentException("Raw keys need captured bytes; middle masking requires forward words of length at least three");
		}
		if(sketch!=null && sketch.maxArray().length!=buckets){
			throw new IllegalArgumentException("Raw keys and native query must retain the same bucket width");
		}
		if(sketch!=null && (DynamicDemiLog.exponentBits()!=4 || DynamicDemiLog.blacklistExists())){
			throw new IllegalStateException("A native sketch destination requires exponent4 and no process-global DDL blacklist");
		}
		keys.clear();
		packed.clear();
		wordMask=AbstractKmerTable.makeMiddleMask(k, false, middle ? 1 : 0);
		destination=sketch;
		try{CardinalityTracker.forEachKmer(bases, k, k, canonical, false, this);}
		finally{destination=null;}
		assert(keys.size()<=Math.max(0, bases.length-k+1))
			: "Each complete native k-mer adds at most one distinct packed bucket/key";
	}

	/** Internal primitive sink; the native destination is borrowed only during build. */
	@Override
	public void accept(final long number, final int length){
		assert(length==k) : "Raw query collection must keep its fixed native word length";
		final long word=number&wordMask;
		if(destination!=null){destination.add(word);}
		add(word);
	}

	/** Counts profile buckets hit by any raw query key, at most once per bucket. */
	public int equalBuckets(final char[] profile){
		if(profile==null || profile.length!=buckets){
			throw new IllegalArgumentException("Query keys and member profile must use the same width");
		}
		int matches=0;
		for(int bucket=0; bucket<buckets; bucket++){
			final int key=profile[bucket];
			if(key!=0 && keys.contains((bucket<<16)|key)){matches++;}
		}
		assert(matches<=buckets) : "Raw-key matching counts each profile bucket at most once";
		return matches;
	}

	/** Independent diagnostic snapshot; the ranking loop does not call this. */
	public int[] packedKeys(){return keys.toArray();}
	/** Borrows read-only storage for keyCount() live keys until the next build.
	 * Unused capacity is not part of the query; callers must not modify it. */
	public int[] packedArray(){return packed.array;}

	/** Replaces worker counters through the shared index, preserving codec width. */
	public int[] count(final DDLIndexSmall index, final int[] counts){
		if(index==null || index.buckets()!=buckets){
			throw new IllegalArgumentException("Raw query and shared index must retain the same bucket width");
		}
		return index.query(packed.array, packed.size, counts);
	}
	/** Primitive read-only iteration for the inverted profile scorer. */
	public int keyCount(){return packed.size;}
	public int buckets(){return buckets;}
	/** Read-only membership for construction-time consensus preference. */
	public boolean contains(final int packedKey){return keys.contains(packedKey);}
	public int packedKey(final int index){
		assert(index>=0 && index<packed.size) : "Inverted lookup may visit only the current query's unique live keys";
		return packed.array[index];
	}

	/** Absolute native key, deliberately before bucket-max/floor suppression. */
	void add(final long number){
		final long hash=Tools.hash64shift(number^12345L);
		final int bucket=(int)(hash&bucketMask);
		final int value=(bucket<<16)|encodeKey(hash);
		if(keys.add(value)){packed.add(value);}
	}

	/** Fixed exponent4 key from the shared native encoder, before insertion. */
	public static int encodeKey(final long hash){
		return DynamicDemiLog.encodeKey(hash, 4, 12);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private final int buckets, k, bucketMask;
	private final IntHashSet keys;
	private final IntList packed=new IntList();
	private DynamicDemiLog destination;
	private long wordMask=-1L;
}
