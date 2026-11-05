package test.ukmer;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.util.BitSet;

import kmer.AbstractKmerTableSet;
import structures.IntList;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayU1D;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/**
 * Baselines selected-position count filling on actual long-kmer table geometry.
 * Each case has one failing target so baseline runners can capture every known
 * selected-window mask or missing-count error in a fresh JVM.
 * @author Jean
 */
public final class KmerTableSetUFillSpecificCountsTest{

	/** Exercises one masked, unmasked or invalid-window count-filling case. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected case="+CASES);}
		if(!args[0].startsWith("case=")){throw new IllegalArgumentException("Expected case= argument: "+args[0]);}
		final String mode=args[0].substring(5);

		final Field orientation=Kmer.class.getDeclaredField("rcomp");
		orientation.setAccessible(true);
		final boolean oldPacked=Kmer.PACKED, oldRcomp=orientation.getBoolean(null);
		final boolean oldMaskCore=AbstractKmerTableSet.MASK_CORE;
		final boolean oldMaskMiddle=AbstractKmerTableSet.MASK_MIDDLE;
		try{
			Kmer.PACKED=false;
			AbstractKmerTableSet.MASK_CORE=false;
			AbstractKmerTableSet.MASK_MIDDLE=false;

			final FixtureSet set=new FixtureSet("in=/dev/null", "k=62", "ways=1", "threads=1",
					"prealloc=f", "prefilter=0", "onepass=f", "rcomp=f", "qin=33", "interleaved=f");
			set.installForTest(new HashArrayU1D(new int[] {101}, set.k, set.kbig));
			check(set.kbig()==62, mode, "The test expects the requested 62-base long-kmer geometry.");
			check(set.tables().length==1, mode, "The fixture should install exactly one shard.");
			check(set.tables()[0].arrayLength()==161, mode, "The fixture must keep a tiny 101-prime primary table plus probe tail.");

			final byte[] bases="ACGTTGCAACGAGTCCGATGATCGTACGTTAGCATCGATGCTAGCATGACCGTTAACGTTGCATT".
					getBytes(StandardCharsets.US_ASCII);
			check(bases.length==set.kbig()+3, mode, "The fixture must expose exactly four candidate starts.");
			if(!isMissingFirst(mode)){put(set, bases, 0, 7, mode);}
			put(set, bases, 1, 11, mode);
			put(set, bases, 2, 13, mode);

			runCase(set, bases, mode);
		}finally{
			AbstractKmerTableSet.MASK_MIDDLE=oldMaskMiddle;
			AbstractKmerTableSet.MASK_CORE=oldMaskCore;
			orientation.setBoolean(null, oldRcomp);
			Kmer.PACKED=oldPacked;
		}
		System.out.println("PASS "+bugId(mode)+" case="+mode);
	}

	/** Dispatches the selected one-target fixture case. */
	private static void runCase(final KmerTableSetU set, final byte[] bases, final String mode){
		if(mode.equals("first")){
			expect(set, bases, bits(0), mode, 1, 7, 0, 0, 0);
		}else if(mode.equals("second")){
			expect(set, bases, bits(1), mode, 1, 0, 11, 0, 0);
		}else if(mode.equals("empty")){
			expect(set, bases, new BitSet(), mode, 0, 0, 0, 0, 0);
		}else if(mode.equals("all")){
			expect(set, bases, null, mode, 4, 7, 11, 13, 0);
		}else if(mode.equals("invalid")){
			bases[0]='N';
			expect(set, bases, bits(1), mode, 1, 0, 11, 0, 0);
		}else if(mode.equals("missing-first-null")){
			expect(set, bases, null, mode, 4, 0, 11, 13, 0);
		}else if(mode.equals("missing-first-selected")){
			expect(set, bases, bits(0), mode, 1, 0, 0, 0, 0);
		}else{
			throw new IllegalArgumentException("Expected case="+CASES+", not case="+mode);
		}
	}

	/** Returns true for UKM111 cases with an absent first-window key. */
	private static boolean isMissingFirst(final String mode){
		return mode.startsWith("missing-first-");
	}

	/** Returns the finding family exercised by the selected case. */
	private static String bugId(final String mode){
		return isMissingFirst(mode) ? "UKM111" : "UKM110";
	}

	/** Installs a bounded table shard without invoking ScheduleMaker allocation. */
	private static final class FixtureSet extends KmerTableSetU{
		FixtureSet(String... args){super(args, 0);}

		void installForTest(final AbstractKmerTableU table) throws Exception{
			TABLES.set(this, new AbstractKmerTableU[] {table});
		}

		private static final Field TABLES;
		static{
			try{
				TABLES=KmerTableSetU.class.getDeclaredField("tables");
				TABLES.setAccessible(true);
			}catch(Exception e){throw new ExceptionInInitializerError(e);}
		}
	}

	/** Inserts a literal count for the real k-mer at one candidate start. */
	private static void put(final KmerTableSetU set, final byte[] bases, final int start,
			final int count, final String mode){
		final Kmer kmer=kmerAt(set, bases, start, mode);
		final AbstractKmerTableU table=set.tables()[set.kmerToWay(kmer)];
		check(table.set(kmer, count)==1, mode, "Fixture start "+start+" must add a distinct key.");
	}

	/** Checks the returned valid count and every reported start count. */
	private static void expect(final KmerTableSetU set, final byte[] bases, final BitSet positions,
			final String mode, final int valid, final int... expected){
		final IntList counts=new IntList(expected.length);
		final int observed=set.fillSpecificCounts(bases, counts, positions, new Kmer(set.kbig()));
		check(observed==valid, mode, "Wrong valid count: expected "+valid+", got "+observed);
		check(counts.size==expected.length, mode,
				"Wrong result length: expected "+expected.length+", got "+counts.size);
		for(int i=0; i<expected.length; i++){
			check(counts.get(i)==expected[i], mode,
					"Wrong count at start "+i+": expected "+expected[i]+", got "+counts.get(i));
		}
	}

	/** Builds one complete k-mer through the production rolling API. */
	private static Kmer kmerAt(final KmerTableSetU set, final byte[] bases, final int start,
			final String mode){
		final Kmer kmer=new Kmer(set.kbig());
		for(int i=start; i<start+set.kbig(); i++){kmer.addRight(bases[i]);}
		check(kmer.len()==set.kbig(), mode, "Fixture start "+start+" must build one full k-mer.");
		return kmer;
	}

	/** Returns a mask with one selected start. */
	private static BitSet bits(final int bit){
		final BitSet bits=new BitSet();
		bits.set(bit);
		return bits;
	}

	/** Keeps result comparisons active when Java assertions are disabled. */
	private static void check(final boolean condition, final String mode, final String message){
		if(!condition){throw new AssertionError(bugId(mode)+" case="+mode+": "+message);}
	}

	private static final String CASES=
			"first|second|empty|all|invalid|missing-first-null|missing-first-selected";
}
