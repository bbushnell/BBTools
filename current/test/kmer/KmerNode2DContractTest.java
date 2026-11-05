package test.kmer;

import java.util.Arrays;

import kmer.AbstractKmerTable;
import kmer.HashForest;
import kmer.KmerNode2D;

/**
 * Checks Node2D logical-prefix merging and scalar sentinel initialization.
 * Runs one literal-oracle case per fresh JVM, with assertions enabled or disabled.
 * Exercises public node and forest APIs; does not establish affected Seal workloads
 * or change the separate bounded-duplicate policy.
 * @author Sucrose
 */
public final class KmerNode2DContractTest{

	/** Accepts case=node-|forest- followed by merge-prefix, merge-full, sentinel-init or positive-init. */
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("case="), "Expected exactly one case= argument");
		final String name=args[0].substring(5);
		final boolean forest=name.startsWith("forest-");
		check(forest || name.startsWith("node-"), "Expected node- or forest- case prefix: "+name);
		final String mode=name.substring(forest ? 7 : 5);
		final boolean prefix=mode.equals("merge-prefix"), full=mode.equals("merge-full");
		final boolean sentinel=mode.equals("sentinel-init"), positive=mode.equals("positive-init");
		check(prefix || full || sentinel || positive, "Unknown case: "+name);
		final int[] expected=prefix ? new int[] {5, 7} : full ? new int[] {5, 7, 91} : new int[] {5, 7, 11};
		System.out.println("case="+name+", expected="+Arrays.toString(expected));
		final AbstractKmerTable table=create(forest, sentinel ? -1 : 5);
		if(prefix || full){
			final int vlen=prefix ? 1 : 2;
			System.out.println("merge=[7, 91, -1], vlen="+vlen);
			check(table.set(KEY, new int[] {7, 91, -1}, vlen)==0,
				"Merging values into the existing key must not report a newly created key");
		}else{
			if(sentinel){add(table, 5);}
			add(table, 7);
			add(table, 11);
		}
		checkValues(table, expected, name);
		System.out.println("PASS case="+name);
	}

	/** Initializes the same key through a direct node or a nonresizing multi-value forest. */
	private static AbstractKmerTable create(final boolean forest, final int initial){
		if(!forest){return new KmerNode2D(KEY, initial);}
		final HashForest table=new HashForest(1, false, true);
		check(table.set(KEY, initial)==1, "The first insertion into an empty forest must create one key");
		return table;
	}

	/** Announces each insertion before it runs so a baseline exception identifies its triggering value. */
	private static void add(final AbstractKmerTable table, final int value){
		System.out.println("add="+value+", before="+Arrays.toString(table.getValues(KEY, null)));
		check(table.set(KEY, value)==0, "Adding a value to the existing key must not create another key");
	}

	/** Checks exact values and negative padding without prescribing the backing array capacity. */
	private static void checkValues(final AbstractKmerTable table, final int[] expected, final String name){
		final int[] actual=table.getValues(KEY, null);
		final String diagnostic="case="+name+", expected="+Arrays.toString(expected)+", actual="+Arrays.toString(actual);
		System.out.println(diagnostic);
		check(actual!=null && actual.length>=expected.length, "Missing literal values; "+diagnostic);
		for(int i=0; i<expected.length; i++){
			check(actual[i]==expected[i], "Literal prefix differs at index "+i+"; "+diagnostic);
		}
		for(int i=expected.length; i<actual.length; i++){
			check(actual[i]<0, "Unexpected value beyond the logical prefix at index "+i+"; "+diagnostic);
		}
	}

	/** Keeps the oracle active under -da and reports the violated contract. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** One nonnegative literal key, shared by all cases; no collisions or duplicate values are needed. */
	private static final long KEY=42;
}
