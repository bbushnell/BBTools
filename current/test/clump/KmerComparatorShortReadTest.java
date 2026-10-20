package test.clump;

import java.nio.charset.StandardCharsets;

import clump.KmerComparator;
import stream.Read;

/**
 * Finite oracle for KmerComparator teaching examples on fillShort fallbacks.
 * @author Sucrose
 */
public final class KmerComparatorShortReadTest{

	/** Runs baseline or fixed expectations for short and all-N literal reads. */
	public static void main(String[] args){
		check(args.length==1, "Expected baseline or fixed");
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown expectation: "+args[0]);
		final KmerComparator comparator=new KmerComparator(5, false, false);

		System.out.println("case\texample\tobserved\texpected");
		checkFallback(comparator, fixed, "short", "ACG", ACG_FILL_SHORT);
		checkFallback(comparator, fixed, "all-n", "NNNNNNNN", ALL_N_FILL_SHORT);
		checkFullLength(comparator);
	}

	/** Checks all public teaching examples against literal fallback pivots. */
	private static void checkFallback(final KmerComparator comparator, final boolean fixed,
			final String name, final String bases, final long expected){
		for(int example=1; example<=3; example++){
			long observed=Long.MIN_VALUE;
			boolean threw=false;
			try{
				observed=callExample(comparator, example, read(name+"-"+example, bases));
			}catch(NullPointerException e){
				threw=true;
			}
			check(threw!=fixed, "Unexpected fillShort fallback for "+name+
					" example"+example);
			check(threw || observed==expected, "Unexpected fallback kmer for "+
					name+" example"+example+": expected "+expected+", got "+observed);
			print(name, example, threw, observed, expected);
		}
	}

	/** Checks full-length valid reads that never enter the fillShort fallback. */
	private static void checkFullLength(final KmerComparator comparator){
		checkExample(comparator, "full", "AAAAA", 1, 0);
		checkExample(comparator, "full", "AAAAA", 2, 1023);
		checkExample(comparator, "full", "AAAAA", 3, 1023);
	}

	/** Checks one teaching example against an expected literal pivot. */
	private static void checkExample(final KmerComparator comparator, final String name,
			final String bases, final int example, final long expected){
		final long observed=callExample(comparator, example, read(name+"-"+example, bases));
		check(observed==expected, "Unexpected kmer for "+name+" example"+example+
				": expected "+expected+", got "+observed);
		print(name, example, false, observed, expected);
	}

	/** Calls one of the three public teaching examples. */
	private static long callExample(final KmerComparator comparator, final int example,
			final Read read){
		if(example==1){return comparator.example1(read);}
		if(example==2){return comparator.example2(read);}
		return comparator.example3(read);
	}

	/** Creates a Read from a literal base string. */
	private static Read read(final String id, final String bases){
		return new Read(bases.getBytes(StandardCharsets.US_ASCII), null, id, 0);
	}

	/** Prints one case result. */
	private static void print(final String name, final int example, final boolean threw,
			final long observed, final long expected){
		System.out.println(name+"\texample"+example+'\t'+
				(threw ? "NullPointerException" : Long.toString(observed))+
				'\t'+expected);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private KmerComparatorShortReadTest(){}

	/** fillShort fallback: ACG encodes as max(6, 432). */
	private static final long ACG_FILL_SHORT=432;

	/** fillShort fallback: Ns encode as zero for the short-read pivot. */
	private static final long ALL_N_FILL_SHORT=0;
}
