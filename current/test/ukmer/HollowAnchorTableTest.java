package test.ukmer;

import java.util.ConcurrentModificationException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import ukmer.HollowAnchorTable;

/** Deterministic exact-key, strand, collision and growth fixtures for the experimental endpoint index.
 * @author Brian Bushnell, Fischl */
public final class HollowAnchorTableTest {

	public static void main(final String[] args){
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Hollow anchor fixtures require -ea.");}
		testSentinelsAndPalindromes();
		testSharedFirstArm();
		testRandomOracle();
		testMutationContract();
		System.err.println("HOLLOW_ANCHOR_TABLE_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	/** All bits are biological sequence; reverse-complement palindromes must be emitted once. */
	private static void testSentinelsAndPalindromes(){
		final HollowAnchorTable table=new HollowAnchorTable(224, 1);
		check(table.get(0, 0)==0, "Absent all-zero pair must have zero count.");
		check(table.increment(0, 0)==1, "All-zero pair is a valid key.");
		check(table.increment(-1L, -1L)==2, "All-one pair shares the reverse-complement count.");
		check(table.size()==1, "Reverse-complement insertion must not duplicate the pair.");
		check(table.add(-1L, 0, Integer.MAX_VALUE)==Integer.MAX_VALUE, "Saturation fixture setup failed.");
		check(table.increment(-1L, 0)==Integer.MAX_VALUE, "Support must saturate, never wrap negative.");
		final HashMap<Long, Integer> found=partners(table, -1L);
		check(found.size()==2 && found.get(0L)==Integer.MAX_VALUE && found.get(-1L)==2,
				"Self-reverse-complement pair must occur once alongside the all-one partner.");
		check(table.arrayBytes()==28L*table.capacity(), "Array payload accounting excludes object/JVM overhead.");
		boolean rejected=false;
		try{new HollowAnchorTable(64, 4);}catch(IllegalArgumentException e){rejected=true;}
		check(rejected, "Overlapping/non-Hollow span must be rejected.");
	}

	/** A deliberately shared first arm forces long probe chains and repeated growth. */
	private static void testSharedFirstArm(){
		final HollowAnchorTable table=new HollowAnchorTable(128, 4);
		for(int i=1; i<=256; i++){
			check(table.add(-1L, i, i)==i, "Insertion must return the count after resizing.");
		}
		check(table.size()==256 && table.capacity()>4, "Fixture must grow while retaining all partners.");
		final HashMap<Long, Integer> found=partners(table, -1L);
		check(found.size()==256, "Lookup must enumerate all partners, not only the first hash hit.");
		for(int i=1; i<=256; i++){
			check(table.get(-1L, i)==i && found.get((long)i)==i, "Shared-arm count changed during growth.");
			final HashMap<Long, Integer> reversed=partners(table, reverse(i));
			check(reversed.size()==1 && reversed.get(0L)==i, "Reverse index must survive slot relocation.");
		}
	}

	/** A small boxed oracle is test-only; its base-by-base RC does not reuse the table helper. */
	private static void testRandomOracle(){
		final Random random=new Random(517294);
		final HollowAnchorTable table=new HollowAnchorTable(95, 4);
		final HashMap<Pair, Integer> oracle=new HashMap<Pair, Integer>();
		final long[] anchors=new long[48];
		for(int i=0; i<anchors.length; i++){anchors[i]=random.nextLong();}
		for(int i=0; i<12000; i++){
			final long a=anchors[random.nextInt(anchors.length)];
			final long b=(i%3==0 ? anchors[random.nextInt(anchors.length)] : random.nextLong());
			check(reverse(b)==HollowAnchorTable.reverseComplement(b), "Fast RC must match literal32-base reversal.");
			final Pair key=canonical(a, b);
			final Integer old=oracle.get(key);
			final int count=old==null ? 1 : old+1;
			oracle.put(key, count);
			check(table.increment(a, b)==count, "Exact pair count disagrees with independent oracle.");
			check(table.get(reverse(b), reverse(a))==count, "RC lookup must return the same support.");
		}
		check(table.size()==oracle.size(), "Unique pair cardinality disagrees with oracle.");
		for(long anchor : anchors){
			checkPartners(table, oracle, anchor);
			checkPartners(table, oracle, reverse(anchor));
		}
		table.freeze();
		check(table.frozen(), "Freeze state must be observable by the owner.");
		checkPartners(table, oracle, anchors[0]);
		boolean rejected=false;
		try{table.increment(0, 0);}catch(IllegalStateException e){rejected=true;}
		check(rejected && table.size()==oracle.size(), "Frozen writes must fail without changing the table.");
	}

	/** Callback mutation must fail visibly instead of silently changing iteration coverage. */
	private static void testMutationContract(){
		final HollowAnchorTable table=new HollowAnchorTable(96, 4);
		table.increment(-1L, 17);
		boolean rejected=false;
		try{
			table.forEachPartner(-1L, new HollowAnchorTable.PartnerConsumer(){
				@Override
				public void accept(final long partner, final int count){table.increment(-1L, partner);}
			});
		}catch(ConcurrentModificationException e){rejected=true;}
		check(rejected, "Mutating callback must be rejected rather than producing a partial partner list.");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Oracle helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Compares all partners/counts in the requested orientation, independent of iteration order. */
	private static void checkPartners(final HollowAnchorTable table, final Map<Pair, Integer> oracle, final long anchor){
		final HashMap<Long, Integer> expected=new HashMap<Long, Integer>();
		for(Map.Entry<Pair, Integer> entry : oracle.entrySet()){
			final Pair pair=entry.getKey();
			if(pair.a==anchor){expected.put(pair.b, entry.getValue());}
			if(reverse(pair.b)==anchor){expected.put(reverse(pair.a), entry.getValue());}
		}
		check(expected.equals(partners(table, anchor)), "Partner set/count differs from literal pair oracle.");
	}

	/** Rejects duplicate callbacks even when they would overwrite an equal map value. */
	private static HashMap<Long, Integer> partners(final HollowAnchorTable table, final long anchor){
		final HashMap<Long, Integer> found=new HashMap<Long, Integer>();
		final int visited=table.forEachPartner(anchor, new HollowAnchorTable.PartnerConsumer(){
			@Override
			public void accept(final long partner, final int count){
				check(found.put(partner, count)==null, "A directed pair was emitted twice.");
			}
		});
		check(visited==found.size(), "Returned visitation count must equal delivered unique pairs.");
		return found;
	}

	/** Slow literal RC oracle, intentionally independent of the production bit-swap routine. */
	private static long reverse(long word){
		long reverse=0;
		for(int i=0; i<32; i++){reverse=(reverse<<2)|(3-(word&3)); word>>>=2;}
		return reverse;
	}

	/** Picks the unsigned larger of the two equivalent oriented endpoint pairs. */
	private static Pair canonical(final long a, final long b){
		final long ra=reverse(a), rb=reverse(b);
		final int cmp=Long.compareUnsigned(a, rb);
		return cmp<0 || (cmp==0 && Long.compareUnsigned(b, ra)<0) ? new Pair(rb, ra) : new Pair(a, b);
	}

	private static void check(final boolean condition, final String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}

	/** Immutable test-oracle key; not used by the primitive production table. */
	private static final class Pair {
		Pair(final long a_, final long b_){a=a_; b=b_;}
		@Override
		public int hashCode(){return Long.hashCode(a)^Integer.rotateLeft(Long.hashCode(b), 13);}
		@Override
		public boolean equals(final Object other){
			if(!(other instanceof Pair)){return false;}
			final Pair pair=(Pair)other;
			return a==pair.a && b==pair.b;
		}
		final long a, b;
	}

	private static int checks;
}
