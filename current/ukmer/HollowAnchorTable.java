package ukmer;

import java.util.ConcurrentModificationException;

import dna.AminoAcid;
import shared.Tools;

/**
 * Experimental exact outer32bp pair counts with lookup from either endpoint.
 * Each reverse-complement-equivalent pair is stored once at a fixed span.
 * The primary probe sequence depends only on the first arm, allowing all of
 * its partners to be enumerated. A primitive reverse index exposes the other
 * direction without duplicating counts or retaining the masked interior.
 * Loading is single-threaded; freeze and safely publish before parallel reads.
 * This is not installed as a Tadpole count-table backend.
 * @author Brian Bushnell, Fischl
 */
public final class HollowAnchorTable {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a table whose literal arms are separated by span-32 start positions. */
	public HollowAnchorTable(final int span_, final int initialCapacity){
		if(span_<=64){throw new IllegalArgumentException("Nonoverlapping Hollow arms require span>64: "+span_);}
		if(initialCapacity<1 || initialCapacity>MAX_CAPACITY){
			throw new IllegalArgumentException("Invalid Hollow anchor capacity: "+initialCapacity);
		}
		span=span_;
		int capacity=4;
		while(capacity<initialCapacity){capacity<<=1;}
		allocate(capacity);
		assert(first.length>=initialCapacity) : "Power-of-two rounding must retain requested capacity.";
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds one observation, sharing its count with the reverse-complement pair. */
	public int increment(final long a, final long b){return add(a, b, 1);}

	/** Adds positive support and returns the saturated count, including after growth. */
	public int add(long a, long b, final int amount){
		if(frozen){throw new IllegalStateException("A frozen Hollow anchor table cannot be modified.");}
		if(amount<=0){throw new IllegalArgumentException("Support must be positive: "+amount);}
		final long ra=reverseComplement(a), rb=reverseComplement(b);
		if(preferReverse(a, b, ra, rb)){a=rb; b=ra;}
		int cell=findCell(a, b);
		if(counts[cell]>0){
			counts[cell]=(int)Math.min(Integer.MAX_VALUE, counts[cell]+(long)amount);
		}else{
			if(size>=sizeLimit){grow(); cell=findCell(a, b);}
			insert(cell, a, b, amount);
		}
		modifications++;
		assert(counts[cell]>0) : "Zero is reserved for an empty slot; saturation must preserve occupied entries.";
		return counts[cell];
	}

	/** Returns exact pair support; neither key word reserves a sentinel bit pattern. */
	public int get(long a, long b){
		final long ra=reverseComplement(a), rb=reverseComplement(b);
		if(preferReverse(a, b, ra, rb)){a=rb; b=ra;}
		final int cell=findCell(a, b);
		assert(cell>=0 && cell<counts.length) : "Pair probing must terminate at an occupied match or empty slot.";
		return counts[cell];
	}

	/**
	 * Visits every oriented right partner of the supplied literal32bp arm once.
	 * Counts are per pair, not sums across different partners. Callback order is
	 * unspecified and can change after resizing. The callback must not mutate
	 * this table; mutation is detected and throws rather than silently skipping keys.
	 */
	public int forEachPartner(final long anchor, final PartnerConsumer consumer){
		if(consumer==null){throw new IllegalArgumentException("Partner enumeration requires a consumer.");}
		final long version=modifications;
		int visited=0;
		for(int cell=bucket(anchor); counts[cell]>0; cell=(cell+1)&mask){
			if(first[cell]==anchor){
				consumer.accept(second[cell], counts[cell]);
				checkUnmodified(version);
				visited++;
			}
		}
		for(int link=reverseHeads[bucket(anchor)]; link!=0; ){
			final int cell=link-1;
			assert(counts[cell]>0) : "Reverse links must reference live pair entries (insert/grow).";
			link=reverseNext[cell];
			if(reverseComplement(second[cell])!=anchor){continue;}
			final long partner=reverseComplement(first[cell]);
			// A self-reverse-complement pair was already emitted by the primary scan.
			if(first[cell]==anchor && second[cell]==partner){continue;}
			consumer.accept(partner, counts[cell]);
			checkUnmodified(version);
			visited++;
		}
		return visited;
	}

	/** Prevents subsequent writes; callers still provide safe publication between threads. */
	public void freeze(){frozen=true;}

	public int size(){return size;}
	public int capacity(){return counts.length;}
	public boolean frozen(){return frozen;}

	/** Exact primitive-array payload only, not a measurement of JVM heap requirements. */
	public long arrayBytes(){return 28L*counts.length;}

	/** Callback receives the other literal32bp arm in the queried direction. */
	public interface PartnerConsumer {
		void accept(long partner, int count);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Internal methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses unsigned lexicographic max, matching Hollow Kmer strand selection. */
	private static boolean preferReverse(final long a, final long b, final long ra, final long rb){
		final int cmp=Long.compareUnsigned(a, rb);
		return cmp<0 || (cmp==0 && Long.compareUnsigned(b, ra)<0);
	}

	/** Complements and reverses all32 bases; there is no unused high-bit sentinel. */
	public static long reverseComplement(final long word){
		return AminoAcid.reverseComplementBinaryFast(word, 32);
	}

	/** Hashes only the first arm; both words remain necessary for equality. */
	private int bucket(final long arm){return (int)Tools.hash64shift(arm)&mask;}

	/** Finds either the exact pair or its insertion slot in a non-full table. */
	private int findCell(final long a, final long b){
		assert(size<counts.length) : "The load limit must leave an empty probe terminator.";
		int cell=bucket(a);
		while(counts[cell]!=0 && (first[cell]!=a || second[cell]!=b)){cell=(cell+1)&mask;}
		return cell;
	}

	/** Inserts a new canonical pair and links its reverse direction without changing its count. */
	private void insert(final int cell, final long a, final long b, final int count){
		assert(count>0 && counts[cell]==0) : "Pair insertion needs an empty slot and positive occupancy count.";
		first[cell]=a;
		second[cell]=b;
		counts[cell]=count;
		final int bucket=bucket(reverseComplement(b));
		reverseNext[cell]=reverseHeads[bucket];
		reverseHeads[bucket]=cell+1;
		size++;
	}

	/** Rebuilds both indices after slot positions change; old reverse links are never reused. */
	private void grow(){
		if(counts.length>=MAX_CAPACITY){throw new IllegalStateException("Hollow anchor table reached its maximum capacity.");}
		final long[] oldFirst=first, oldSecond=second;
		final int[] oldCounts=counts;
		final int expected=size;
		allocate(counts.length<<1);
		size=0;
		for(int i=0; i<oldCounts.length; i++){
			if(oldCounts[i]>0){insert(findCell(oldFirst[i], oldSecond[i]), oldFirst[i], oldSecond[i], oldCounts[i]);}
		}
		assert(size==expected) : "Resizing must conserve all distinct pairs and rebuild reverse links: "+size+" != "+expected;
	}

	/** Allocates zeroed arrays; occupancy and one-based reverse links share zero-as-empty. */
	private void allocate(final int capacity){
		assert(capacity>=4 && Integer.bitCount(capacity)==1) : "Mask-based probing requires a power-of-two capacity.";
		first=new long[capacity];
		second=new long[capacity];
		counts=new int[capacity];
		reverseHeads=new int[capacity];
		reverseNext=new int[capacity];
		mask=capacity-1;
		sizeLimit=(capacity*2)/3;
	}

	/** Prevents a callback from invalidating the arrays or counts currently being enumerated. */
	private void checkUnmodified(final long expected){
		if(modifications!=expected){throw new ConcurrentModificationException("Partner callback modified its Hollow anchor table.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final int span;
	private long[] first, second;
	private int[] counts, reverseHeads, reverseNext;
	private int size, sizeLimit, mask;
	private long modifications;
	private boolean frozen;
	private static final int MAX_CAPACITY=1<<29;
}
