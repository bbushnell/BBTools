package prok;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import ddl.DDLIndexSmall;
import kmer.AbstractKmerTable;
import map.IntHashSet;
import shared.Tools;
import structures.IntList;

/** Two forward word lengths in one bucket/key space, enumerated in one pass.
 * Short words are plain; long words mask the same single middle base as BBDuk.
 * Competing lengths use disjoint raw-word tags before hashing, then compete for
 * every bucket. These tags do not prevent collisions in the compressed keys.
 * Partitioned lengths retain the untagged hash and occupy separate bucket halves;
 * the long half therefore equals the existing long-only codec at half the width.
 * One owner reuses this object; returned member arrays are borrowed until reuse.
 * @author Brian, Raiden
 */
public final class NcrnaDdlDualKeys implements CardinalityTracker.KmerConsumer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlDualKeys(final int shortK_, final int longK_, final int buckets_, final boolean partitioned_){
		if((shortK_<7 || shortK_>9) || (longK_!=10 && longK_!=11)
				|| (buckets_!=2048 && buckets_!=4096)){
			throw new IllegalArgumentException("Dual-length experiment requires short k7/8/9, masked long k10/11 and 2048/4096 total buckets");
		}
		shortK=shortK_;
		longK=longK_;
		buckets=buckets_;
		partitioned=partitioned_;
		shortMask=(1L<<(2*shortK))-1;
		longMask=(1L<<(2*longK))-1;
		middleMask=AbstractKmerTable.makeMiddleMask(longK, false, 1);
		half=buckets/2;
		ends=partitioned ? new int[]{half, buckets} : null;
		bucketMask=(partitioned ? half : buckets)-1;
		member=new char[buckets];
		assert((shortMask&SHORT_TAG)==0 && (longMask&LONG_TAG)==0)
			: "Competing length tags must lie above every supported raw word bit";
	}

	/*--------------------------------------------------------------*/
	/*----------------          Query Keys          ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces the live query without slices, temporary words or per-key objects.
	 * The long rolling word's low bits also contain the current short word.
	 * Ambiguity state follows the corrected CardinalityTracker.hashSmall loop.
	 * Shared strict enumeration retains its A/C/G/T bit codes and sends every
	 * undefined symbol through the same sign-extended N state, including IUPAC.
	 */
	public void build(final byte[] bases){
		if(bases==null){throw new IllegalArgumentException("Dual-length keys require actual captured window bytes");}
		keys.clear();
		packed.clear();
		lastBases=bases;
		memberReady=false;
		CardinalityTracker.forEachKmer(bases, shortK, longK, false, true, this);
		assert(packed.size<=Math.max(0, bases.length-shortK+1)+Math.max(0, bases.length-longK+1))
			: "Each valid span contributes at most one raw bucket/key per length";
	}

	/** Internal primitive sink; masking follows complete-span ambiguity checks. */
	@Override
	public void accept(final long word, final int k){
		assert(k==shortK || k==longK) : "Dual collection must receive only its two declared word lengths";
		add(k==shortK ? word : word&middleMask, k==shortK);
	}

	private void add(final long word, final boolean shortWord){
		assert((word&~(shortWord ? shortMask : longMask))==0)
			: "Only the declared raw word bits may precede a competing length tag";
		final long tag=partitioned ? 0L : shortWord ? SHORT_TAG : LONG_TAG;
		if(buildingMember){
			(partitioned ? shortWord ? nativeShort : nativeLong : nativeMixed).hashAndStore(word^tag);
			return;
		}
		final long hash=Tools.hash64shift(word^tag^12345L);
		final int offset=partitioned && !shortWord ? half : 0;
		final int bucket=(int)(hash&bucketMask)+offset;
		final int value=(bucket<<16)|NcrnaDdlQueryKeys.encodeKey(hash);
		if(keys.add(value)){packed.add(value);}
	}

	public int buckets(){return buckets;}
	public boolean partitioned(){return partitioned;}
	public int keyCount(){return packed.size;}
	/** Borrows read-only storage for keyCount() live keys until the next build.
	 * Unused capacity is not part of the query; callers must not modify it. */
	public int[] packedArray(){return packed.array;}

	/** Replaces worker counters through the shared index, preserving codec width. */
	public int[] count(final DDLIndexSmall index, final int[] counts){
		checkIndex(index);
		return index.query(packed.array, packed.size, counts);
	}

	/** Counts the two intrinsic bin halves with distinct worker-owned rows. */
	public int[][] countHalves(final DDLIndexSmall index, final int[][] counts){
		checkIndex(index);
		if(!partitioned){throw new IllegalArgumentException("Separate half counts require a partitioned dual query");}
		return index.query(packed.array, packed.size, ends, counts);
	}

	private void checkIndex(final DDLIndexSmall index){
		if(index==null || index.buckets()!=buckets){
			throw new IllegalArgumentException("Dual query and shared index must retain the same bucket width");
		}
	}
	public int packedKey(final int index){
		assert(index>=0 && index<packed.size) : "Inverted lookup visits only unique live query keys";
		return packed.array[index];
	}

	/** One winning key per member/bucket before plurality voting. Called only at
	 * profile construction; normal queries never clear or scan the dense array.
	 * Taking the maximum absolute key matches native DDL insertion, including
	 * competition BETWEEN lengths before a member contributes its one vote.
	 */
	public char[] memberKeys(){
		assert(member.length==buckets) : "Member votes and query keys must use the same combined bucket space";
		if(memberReady){return member;}
		if(lastBases==null || DynamicDemiLog.exponentBits()!=4){
			throw new IllegalStateException("Build the sequence and configure native exponent4 before constructing its member DDL");
		}
		if(partitioned){
			if(nativeShort==null){
				nativeShort=DynamicDemiLog.create(half, shortK, 12345, 0, false, false);
				nativeLong=DynamicDemiLog.create(half, longK, 12345, 0, false, false);
			}
			nativeShort.clear();
			nativeLong.clear();
		}else{
			if(nativeMixed==null){nativeMixed=DynamicDemiLog.create(buckets, longK, 12345, 0, false, false);}
			nativeMixed.clear();
		}
		buildingMember=true;
		try{CardinalityTracker.forEachKmer(lastBases, shortK, longK, false, true, this);}
		finally{buildingMember=false;}
		if(partitioned){
			System.arraycopy(nativeShort.maxArray(), 0, member, 0, half);
			System.arraycopy(nativeLong.maxArray(), 0, member, half, half);
		}else{System.arraycopy(nativeMixed.maxArray(), 0, member, 0, buckets);}
		memberReady=true;
		return member;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private final int shortK, longK, buckets, half, bucketMask;
	private final boolean partitioned;
	private final long shortMask, longMask, middleMask;
	private final IntHashSet keys=new IntHashSet(256);
	private final IntList packed=new IntList();
	private final char[] member;
	private final int[] ends;
	private byte[] lastBases;
	private DynamicDemiLog nativeShort, nativeLong, nativeMixed;
	private boolean buildingMember, memberReady;
	private static final long SHORT_TAG=1L<<62, LONG_TAG=1L<<63;
}
