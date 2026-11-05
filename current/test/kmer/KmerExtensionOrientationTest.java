package test.kmer;

import java.util.Arrays;

import kmer.AbstractKmerTableSet;
import kmer.KmerTableSet;
import shared.Shared;

/**
 * Checks DNA extension orientation using literal short-key counts.
 * Run each case in a fresh JVM under both assertion modes. Core/middle masking
 * are disabled so dispatchers use the direct safe paths. No reads are loaded.
 * @author Sucrose
 */
public final class KmerExtensionOrientationTest{

	/** Accepts case=left-forward|right-forward|left-canonical|right-canonical. */
	public static void main(String[] args){
		String mode=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(!("left-forward".equals(mode) || "right-forward".equals(mode) ||
				"left-canonical".equals(mode) || "right-canonical".equals(mode))){
			throw new IllegalArgumentException("Expected one of the four orientation case names");
		}
		final boolean left=mode.startsWith("left"), canonical=mode.endsWith("canonical");
		final int[] expected;
		final int expectedMax;
		if(left){
			expected=canonical ? new int[] {19, 11, 17, 13} : new int[] {2, 7, 5, 3};
			expectedMax=canonical ? 0 : 1;
		}else{
			expected=canonical ? new int[] {19, 37, 31, 41} : new int[] {2, 23, 43, 29};
			expectedMax=canonical ? 3 : 2;
		}

		Shared.AMINO_IN=false;
		AbstractKmerTableSet.MASK_CORE=false;
		AbstractKmerTableSet.MASK_MIDDLE=false;
		AbstractKmerTableSet.FAST_FILL=true;
		final KmerTableSet tables=new KmerTableSet(new String[] {
			"k=3", "ways=3", "threads=1", "prefilter=0", "prealloc=f", "minprob=0", "rcomp="+canonical}, 12);
		check(tables.k==3 && tables.ways==3 && tables.rcomp()==canonical &&
			tables.coreMask==-1L && tables.middleMask==-1L && tables.cmMask==-1L &&
			!Shared.AMINO_IN && !AbstractKmerTableSet.MASK_CORE,
			"Fixture requires three-base DNA and unmasked safe-path dispatch");
		tables.allocateTables();
		if(!canonical){populate(tables, FORWARD_KEYS, FORWARD_COUNTS);}
		populate(tables, REVERSE_KEYS, REVERSE_COUNTS);

		final int[] direct={-101, -102, -103, -104}, dispatched={-201, -202, -203, -204};
		final int directMax, dispatchedMax;
		//AAA=0 and TTT=63 are literal reverse complements; 63 masks three two-bit symbols.
		if(left){
			directMax=tables.fillLeftCounts_safe(0L, 63L, direct, 63L, 4);
			dispatchedMax=tables.fillLeftCounts(0L, 63L, dispatched, 63L, 4);
		}else{
			directMax=tables.fillRightCounts_safe(0L, 63L, direct, 63L, 4);
			dispatchedMax=tables.fillRightCounts(0L, 63L, dispatched, 63L, 4);
		}
		//Evaluate both paths before checking either, so a target failure cannot skip the dispatcher.
		final String diagnostic="case="+mode+", expected="+Arrays.toString(expected)+", expectedMax="+expectedMax+
			", direct="+Arrays.toString(direct)+", directMax="+directMax+
			", dispatched="+Arrays.toString(dispatched)+", dispatchedMax="+dispatchedMax;
		System.out.println(diagnostic);
		check(Arrays.equals(direct, expected) && directMax==expectedMax &&
			Arrays.equals(dispatched, expected) && dispatchedMax==expectedMax,
			"KMR023: extension counts must use the configured orientation; "+diagnostic);
		System.out.println("PASS case="+mode);
	}

	/** Inserts literal raw keys without using the orientation selector to construct the oracle. */
	private static void populate(final KmerTableSet tables, final long[] keys, final int[] counts){
		for(int i=0; i<keys.length; i++){
			check(tables.getTableForKey(keys[i]).set(keys[i], counts[i])==1,
				"Literal fixture keys must be distinct; key="+keys[i]);
			check(tables.getCount(keys[i])==counts[i], "Raw fixture count must round-trip; key="+keys[i]);
		}
	}

	/** Preserves the oracle under -da and reports the failed contract. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** AAA, CAA, GAA, TAA, AAC, AAG, AAT in two-bit encoding. */
	private static final long[] FORWARD_KEYS={0, 16, 32, 48, 1, 2, 3};
	/** Distinct counts for the literal forward keys. */
	private static final int[] FORWARD_COUNTS={2, 7, 5, 3, 23, 43, 29};
	/** TTT, TTG, TTC, TTA, GTT, CTT, ATT; each exceeds its forward partner. */
	private static final long[] REVERSE_KEYS={63, 62, 61, 60, 47, 31, 15};
	/** Counts used for reverse keys and for the seven canonical stored keys. */
	private static final int[] REVERSE_COUNTS={19, 11, 17, 13, 37, 31, 41};
}
