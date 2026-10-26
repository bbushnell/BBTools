package test.assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import assemble.Contig;
import assemble.HollowTipAnchors;
import ukmer.HollowAnchorTable;

/** Assembly-wide uniqueness and oriented endpoint-geometry fixtures.
 * @author Brian Bushnell, Fischl */
public final class HollowTipAnchorsTest {

	public static void main(final String[] args){
		boolean assertions=false;
		assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Tip-anchor fixtures require -ea.");}
		testUniqueGeometry();
		testInternalDuplicates();
		testAmbiguousAndShort();
		System.err.println("HOLLOW_TIP_ANCHORS_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	/** Overlapping tip bands count physical positions once and preserve oriented gap geometry. */
	private static void testUniqueGeometry(){
		final ArrayList<Contig> contigs=new ArrayList<Contig>();
		contigs.add(new Contig(sequence(400, 113), 91));
		contigs.add(new Contig(sequence(300, 117), 52));
		contigs.add(new Contig(sequence(40, 119), 17));
		final HollowTipAnchors index=new HollowTipAnchors(contigs, 64);
		check(index.size()==141 && index.uniqueCount()==141, "Expected66+66+9 unique candidate positions.");
		for(int i=0; i<index.size(); i++){
			check(index.isUnique(i), "Random deterministic fixture must contain no duplicated32bp anchors.");
			final long word=index.forwardWord(i);
			check(index.lookup(word)==i && index.lookup(HollowAnchorTable.reverseComplement(word))==i,
					"Both strands must identify the same assembly position.");
		}
		final long a=encode(contigs.get(0).bases, 350), b=encode(contigs.get(1).bases, 12);
		final int left=index.lookup(a), right=index.lookup(b);
		check(left>=0 && right>=0, "Geometry fixture anchors must be present and unique.");
		check(index.contigIndex(left)==0 && index.start(left)==350 && index.contigLength(left)==400,
				"Owners use list indices, never the unrelated mutable Contig.id.");
		check(!index.reverse(left, a) && index.reverse(left, HollowAnchorTable.reverseComplement(a)),
				"Observed arm strand must orient its contig.");
		check(index.exitDistance(left, false)==18 && index.entryDistance(right, false)==12,
				"Distances exclude the32bp arm itself.");
		check(index.canExit(left, false) && !index.canEnter(left, false), "Right tip must not be admitted as a forward entry.");
		check(index.canEnter(left, true) && !index.canExit(left, true), "Reverse strand must swap the eligible end.");
		check(index.gap(left, false, right, false, 128)==34, "Gap must be128-64-18-12=34.");
		check(index.gap(right, true, left, true, 128)==34, "Reverse connection must preserve signed gap.");
		check(index.gap(left, false, right, false, 65)==-29, "Negative gap represents overlap, not an unsigned distance.");
		boolean rejected=false;
		try{index.reverse(left, a^1);}catch(IllegalArgumentException e){rejected=true;}
		check(rejected, "An unrelated observed word must not silently choose a strand.");
	}

	/** An internal copy, including a reverse-complement copy, invalidates a tip candidate. */
	private static void testInternalDuplicates(){
		for(int mode=0; mode<3; mode++){
			final byte[] first=sequence(400, 221), second=sequence(400, 223);
			final long word=encode(first, 11);
			if(mode==0){System.arraycopy(first, 11, second, 180, 32);}
			else if(mode==1){
				for(int i=0; i<32; i++){second[180+i]=complement(first[42-i]);}
			}else{System.arraycopy(first, 11, first, 180, 32);}
			final ArrayList<Contig> contigs=new ArrayList<Contig>();
			contigs.add(new Contig(first));
			contigs.add(new Contig(second));
			final HollowTipAnchors index=new HollowTipAnchors(contigs, 64);
			check(index.lookup(word)<0 && index.lookup(HollowAnchorTable.reverseComplement(word))<0,
					"Global uniqueness must reject internal, reverse, and same-contig duplicates; mode="+mode);
			check(index.uniqueCount()<index.size(), "The repeated candidate must remain visibly nonunique.");
		}
	}

	/** Short/ambiguous words do not enter the index; all-zero words remain valid biological keys. */
	private static void testAmbiguousAndShort(){
		final ArrayList<Contig> contigs=new ArrayList<Contig>();
		contigs.add(new Contig(sequence(31, 18)));
		final byte[] ambiguous=new byte[80];
		Arrays.fill(ambiguous, (byte)'N');
		contigs.add(new Contig(ambiguous));
		final byte[] palindrome=new byte[32];
		final byte[] motif={'A','C','G','T'};
		for(int i=0; i<palindrome.length; i++){palindrome[i]=motif[i&3];}
		check(encode(palindrome, 0)==HollowAnchorTable.reverseComplement(encode(palindrome, 0)),
				"Palindrome fixture must be self-reverse-complementary.");
		contigs.add(new Contig(palindrome));
		check(new HollowTipAnchors(contigs, 64).size()==0, "Undefined, short and strand-ambiguous words must be excluded.");
		final byte[] adenines=new byte[32], thymines=new byte[32];
		Arrays.fill(adenines, (byte)'A');
		Arrays.fill(thymines, (byte)'T');
		contigs.add(new Contig(adenines));
		final HollowTipAnchors unique=new HollowTipAnchors(contigs, 64);
		check(unique.uniqueCount()==1 && unique.lookup(0)>=0 && unique.lookup(-1L)>=0,
				"All-zero/all-one keys must find the sole homopolymer32bp position.");
		contigs.add(new Contig(thymines));
		check(new HollowTipAnchors(contigs, 64).uniqueCount()==0, "A reverse-complement contig copy must invalidate uniqueness.");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Test helpers         ----------------*/
	/*--------------------------------------------------------------*/

	/** Produces tiny deterministic in-memory sequence fixtures, not a benchmark read corpus. */
	private static byte[] sequence(final int length, final int seed){
		final Random random=new Random(seed);
		final byte[] bases=new byte[length], alphabet={'A','C','G','T'};
		for(int i=0; i<length; i++){bases[i]=alphabet[random.nextInt(4)];}
		return bases;
	}

	/** Literal32bp encoder independent of production rolling-anchor discovery. */
	private static long encode(final byte[] bases, final int start){
		long word=0;
		for(int i=start; i<start+32; i++){
			final int code=bases[i]=='A' ? 0 : bases[i]=='C' ? 1 : bases[i]=='G' ? 2 : 3;
			word=(word<<2)|code;
		}
		return word;
	}

	private static byte complement(final byte base){
		return (byte)(base=='A' ? 'T' : base=='C' ? 'G' : base=='G' ? 'C' : 'A');
	}

	private static void check(final boolean condition, final String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}

	private static int checks;
}
