package test.clump2;

import java.nio.charset.StandardCharsets;

import bloom.KCountArray;
import clump2.KmerComparator;
import stream.Read;

/** Finite oracle for KmerComparator minprob fallback replacement.
 * @author Sucrose
 */
public final class KmerComparatorMinProbTest{

	/** Runs literal k=1 pivot cases for minprob and count-path fallback replacement. */
	public static void main(String[] args){
		check(args.length==1, "Expected baseline or fixed");
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown expectation: "+args[0]);

		final float oldMinProb=KmerComparator.minProb;
		final boolean oldVerbose=KmerComparator.verbose;
		try{
			KmerComparator.verbose=false;
			System.out.println("case\texpected\tobserved\tpass");
			checkProb("lower-hash-qualifying", "TG", new byte[] {0, 40}, 0.99f,
					expected(fixed, 3, 2));
			checkProb("higher-hash-fallback", "GT", new byte[] {0, 0}, 0.99f,
					expected(fixed, 2, 3));
			checkProb("higher-hash-qualifying", "GT", new byte[] {0, 40}, 0.99f, 3);
			checkProb("all-qualifying", "TG", new byte[] {40, 40}, 0.99f, 3);
			checkProb("no-qualifying", "TG", new byte[] {0, 0}, 0.99f, 3);
			checkProb("default-minprob0", "TG", new byte[] {0, 40}, 0f, 3);
			checkProb("missing-quality", "TG", null, 0.99f, 3);
			checkCount("count-lower-qualifying", counted(2, 1), 2);
			checkCount("count-no-qualifying", counted(0, 0), 3);
		}finally{
			KmerComparator.minProb=oldMinProb;
			KmerComparator.verbose=oldVerbose;
		}
	}

	/** Checks one probability-filtered read against its literal expected canonical pivot. */
	private static void checkProb(final String name, final String bases,
			final byte[] quals, final float minProb, final long expected){
		KmerComparator.minProb=minProb;
		final KmerComparator comparator=comparator();
		print(name+"-fillMax", expected,
				comparator.hash(read(name+"-fillMax", bases, quals), null, 0, false));
		print(name+"-example3", expected, comparator.example3(read(name+"-example3", bases, quals)));
	}

	/** Checks the table-filtered path against its unchanged count-threshold pivot. */
	private static void checkCount(final String name, final KCountArray table,
			final long expected){
		KmerComparator.minProb=0.99f;
		final long observed=comparator().hash(read(name, "TG", new byte[] {0, 40}),
				table, 2, false);
		print(name, expected, observed);
	}

	/** Selects a literal baseline or fixed expectation for policy-changing rows. */
	private static long expected(final boolean fixed, final long baseline,
			final long repaired){
		return fixed ? repaired : baseline;
	}

	/** Builds a tiny table for the canonical k=1 C/G and A/T buckets. */
	private static KCountArray counted(final int cg, final int at){
		final KCountArray table=KCountArray.makeNew(65536, 2, 1);
		for(int i=0; i<cg; i++){table.increment(2);}
		for(int i=0; i<at; i++){table.increment(3);}
		check(table.read(2)==cg, "Expected count "+cg+" for key 2, got "+
				table.read(2));
		check(table.read(3)==at, "Expected count "+at+" for key 3, got "+
				table.read(3));
		table.shutdown();
		return table;
	}

	/** Creates a deterministic identity-hash k=1 comparator. */
	private static KmerComparator comparator(){
		return new KmerComparator(1, 1, 0, 0, false, false);
	}

	/** Creates a read from literal bases and quality scores. */
	private static Read read(final String id, final String bases, final byte[] quals){
		return new Read(bases.getBytes(StandardCharsets.US_ASCII), quals, id, 0);
	}

	/** Emits and checks a literal pivot expectation. */
	private static void print(final String name, final long expected, final long observed){
		check(observed==expected, "Unexpected pivot for "+name+": expected "+
				expected+", got "+observed);
		System.out.println(name+'\t'+expected+'\t'+observed+"\ttrue");
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private KmerComparatorMinProbTest(){}
}
