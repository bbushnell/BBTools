package prok;

import java.util.Arrays;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import ddl.DDLIndexSmall;
import kmer.AbstractKmerTable;
import map.IntHashSet;
import shared.Tools;
import structures.IntList;

/** Compact member/consensus partitions in one fixed-width bucket/key space.
 * Short k0 selects three-quarter masked10 members plus one-quarter consensus.
 * Short k7/8 selects one-quarter short, two-quarter masked10 members, and one
 * shared consensus quarter containing both tagged lengths. Members never enter
 * the consensus quarter; consensus never enters a member area. Queries probe
 * both applicable areas and retain every distinct final bucket/key pair.
 *
 * One owner reuses rolling state, primitive query keys and construction maxima.
 * Construction returns a borrowed array consumed before the next build; normal
 * queries never clear or scan that dense array. Hash seed12345/exponent4 and the
 * native forward DNA bit codes match NcrnaDdlQueryKeys and DynamicDemiLog.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlPartitionKeys implements CardinalityTracker.KmerConsumer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlPartitionKeys(final int shortK_, final int buckets_){
		this(shortK_, buckets_, shortK_!=0);
	}

	/** Member areas stay fixed; this option selects the lengths admitted to the
	 * consensus quarter and probed there by queries. Existing geometry uses both. */
	public NcrnaDdlPartitionKeys(final int shortK_, final int buckets_, final boolean consensusShort_){
		if((shortK_!=0 && shortK_!=7 && shortK_!=8)
				|| (buckets_!=512 && buckets_!=1024 && buckets_!=2048) || shortK_==0 && consensusShort_){
			throw new IllegalArgumentException("Compact partitions require short k0/7/8 and total width512/1024/2048");
		}
		shortK=shortK_;
		consensusShort=consensusShort_;
		buckets=buckets_;
		quarter=buckets/4;
		memberBins=3*quarter;
		ends=new int[]{quarter, memberBins, buckets};
		maxima=new char[buckets];
	}

	/*--------------------------------------------------------------*/
	/*----------------          Enumeration         ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces raw query keys without a window slice or per-word allocation. */
	public void build(final byte[] bases){build(bases, false);}

	/** Preserves the legacy bin/hash geometry while selecting query form only. */
	public void build(final byte[] bases, final boolean queryDdl){
		if(queryDdl){
			if(shortK==0 || DynamicDemiLog.exponentBits()!=4){
				throw new IllegalArgumentException("A native query DDL requires power-of-two short/long parts and native exponent4");
			}
			if(nativeShort==null){
				nativeShort=DynamicDemiLog.create(quarter, shortK, 12345, 0, false, false);
				nativeLong=DynamicDemiLog.create(2*quarter, 10, 12345, 0, false, false);
				nativeConsensus=DynamicDemiLog.create(quarter, 10, 12345, 0, false, false);
			}
			nativeShort.clear();
			nativeLong.clear();
			nativeConsensus.clear();
			nativeQuery=true;
			try{enumerate(bases, QUERY);}
			finally{nativeQuery=false;}
			System.arraycopy(nativeShort.maxArray(), 0, maxima, 0, quarter);
			System.arraycopy(nativeLong.maxArray(), 0, maxima, quarter, 2*quarter);
			System.arraycopy(nativeConsensus.maxArray(), 0, maxima, memberBins, quarter);
			for(int bin=0; bin<buckets; bin++){
				if(maxima[bin]!=0){packed.add((bin<<16)|maxima[bin]);}
			}
			assert(packed.size<=buckets) : "A legacy-geometry query DDL emits only nonempty greatest keys";
		}else{enumerate(bases, QUERY);}
	}

	/** Native maximum per member/bin, before that member supplies one vote. */
	public char[] memberKeys(final byte[] bases){
		enumerate(bases, MEMBER);
		return maxima;
	}

	/** Independent maximum-key sketch in the consensus quarter only. */
	public char[] consensusKeys(final byte[] bases){
		enumerate(bases, CONSENSUS);
		return maxima;
	}

	private void enumerate(final byte[] bases, final int mode){
		if(bases==null){throw new IllegalArgumentException("Partition keys require the actual sequence bytes");}
		keys.clear();
		packed.clear();
		if(mode!=QUERY){Arrays.fill(maxima, (char)0);}
		buildMode=mode;
		CardinalityTracker.forEachKmer(bases, shortK==0 ? 10 : shortK, 10, false, true, this);
		assert(mode!=QUERY || packed.size<=2L*Math.max(0, bases.length-9)
			+2L*(shortK==0 ? 0 : Math.max(0, bases.length-shortK+1)))
			: "Each valid span contributes at most one member and one consensus query key";
	}

	/** Internal primitive sink; complete-span validation precedes middle masking. */
	@Override
	public void accept(final long word, final int k){
		assert(k==10 || shortK!=0 && k==shortK) : "Partition queries must use only their declared member/consensus spans";
		add(k==10 ? word&MIDDLE_MASK : word, k!=10, buildMode);
	}

	private void add(final long word, final boolean shortWord, final int mode){
		if(nativeQuery){
			(shortWord ? nativeShort : nativeLong).hashAndStore(word);
			if(!shortWord || consensusShort){nativeConsensus.hashAndStore(word^(shortWord ? SHORT_TAG : LONG_TAG));}
			return;
		}
		final long memberHash=mode==CONSENSUS && shortK!=0 ? 0L : Tools.hash64shift(word^12345L);
		if(mode!=CONSENSUS){
			final int bucket=shortK==0 ? familyBucket(memberHash, memberBins)
				: shortWord ? (int)(memberHash&(quarter-1)) : quarter+(int)(memberHash&(2*quarter-1));
			store(bucket, memberHash, mode);
		}
		if(mode!=MEMBER && (!shortWord || consensusShort)){
			final long hash=shortK==0 ? memberHash
				: Tools.hash64shift(word^(shortWord ? SHORT_TAG : LONG_TAG)^12345L);
			store(memberBins+(int)(hash&(quarter-1)), hash, mode);
		}
	}

	private void store(final int bucket, final long hash, final int mode){
		assert(bucket>=0 && bucket<buckets) : "Partition mapping must stay within the one stored profile";
		final char key=(char)NcrnaDdlQueryKeys.encodeKey(hash);
		if(mode==QUERY){
			final int value=(bucket<<16)|key;
			if(keys.add(value)){packed.add(value);}
		}else{
			DynamicDemiLog.addKey(maxima, bucket, key);
		}
	}

	/** Unsigned low32 range reduction; the product must remain a long so the
	 * three-quarter area includes every bin without overflowing or division. */
	public static int familyBucket(final long hash, final int bins){
		assert(bins>0 && bins<=1536) : "Three-quarter compact member areas have at most1536bins";
		return (int)(((hash&0xffffffffL)*bins)>>>32);
	}

	public int buckets(){return buckets;}
	/** Geometry guard for separate short/long/consensus score accounting. */
	public int shortK(){return shortK;}
	public int memberBins(){return memberBins;}
	public int keyCount(){return packed.size;}
	/** Borrows read-only storage for keyCount() live keys until the next build or
	 * member/consensus construction. Unused capacity must not be inspected. */
	public int[] packedArray(){return packed.array;}

	/** Replaces worker counters through the shared index, preserving codec width. */
	public int[] count(final DDLIndexSmall index, final int[] counts){
		checkIndex(index);
		return index.query(packed.array, packed.size, counts);
	}

	/** Counts the original short7/member10/consensus areas with distinct rows. */
	public int[][] countParts(final DDLIndexSmall index, final int[][] counts){
		checkIndex(index);
		if(shortK!=7){throw new IllegalArgumentException("Three-part compact counts require the original short7 geometry");}
		return index.query(packed.array, packed.size, ends, counts);
	}

	private void checkIndex(final DDLIndexSmall index){
		if(index==null || index.buckets()!=buckets){
			throw new IllegalArgumentException("Partition query and shared index must retain the same bucket width");
		}
	}
	public int packedKey(final int index){
		assert(index>=0 && index<packed.size) : "Inverted scoring visits only live unique query keys";
		return packed.array[index];
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private final int shortK, buckets, quarter, memberBins;
	private final boolean consensusShort;
	private int buildMode;
	private boolean nativeQuery;
	private DynamicDemiLog nativeShort, nativeLong, nativeConsensus;
	private final char[] maxima;
	private final IntHashSet keys=new IntHashSet(256);
	private final IntList packed=new IntList();
	private final int[] ends;
	private static final int QUERY=0, MEMBER=1, CONSENSUS=2;
	private static final long MIDDLE_MASK=AbstractKmerTable.makeMiddleMask(10, false, 1);
	private static final long SHORT_TAG=1L<<62, LONG_TAG=1L<<63;
}
