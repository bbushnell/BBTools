package test.kmer;

import kmer.HashArray1D;

/**
 * Checks saturating counts in primary cells and their overflow trees.
 * Each case uses 61 distinct keys hashing to the same cell, with exactly one
 * victim. Run cases separately with -ea and -da; checks do not rely on assertions.
 * @author Sucrose
 */
public final class HashCountOverflowTest {

	/** Selects a primary/victim case and either count-returning or creation-returning update. */
	public static void main(String[] args){
		if(args.length!=1){throw new IllegalArgumentException("Expected case=<main|victim|ordinary>-<count|created>");}
		final String mode=args[0];
		final boolean victim, ordinary, created;
		if(mode.equals("case=main-count")){victim=false; ordinary=false; created=false;}
		else if(mode.equals("case=main-created")){victim=false; ordinary=false; created=true;}
		else if(mode.equals("case=victim-count")){victim=true; ordinary=false; created=false;}
		else if(mode.equals("case=victim-created")){victim=true; ordinary=false; created=true;}
		else if(mode.equals("case=ordinary-count")){victim=true; ordinary=true; created=false;}
		else if(mode.equals("case=ordinary-created")){victim=true; ordinary=true; created=true;}
		else{throw new IllegalArgumentException("Unknown case: "+mode);}
		run(victim, ordinary, created);
		System.out.println("PASS "+mode+" keys=61 expectedCount="+(ordinary ? 7 : Integer.MAX_VALUE));
	}

	/** Uses a literal 2-count increment across the signed-int boundary, or a 5+2 control. */
	private static void run(final boolean victim, final boolean ordinary, final boolean created){
		final HashArray1D table=new HashArray1D(new int[] {127}, Long.MAX_VALUE);
		for(int i=0; i<61; i++){
			check(table.set(127L*i, 1)==1, "Fixture must create each collision key once: "+i);
		}
		check(table.size()==60 && table.victims().size()==1, "Fixture must have 60 primary keys and one victim");
		final long key=victim ? 7620L : 0L;
		check(table.victims().contains(key)==victim, "Target is in the wrong storage tier");
		final int initial=ordinary ? 5 : Integer.MAX_VALUE-1;
		final int expected=ordinary ? 7 : Integer.MAX_VALUE;
		check(table.set(key, initial)==0 && table.getValue(key)==initial, "Could not seed the existing target count");
		final int result=created ? table.incrementAndReturnNumCreated(key, 2) : table.increment(key, 2);
		check(result==(created ? 0 : expected), "Wrong update return: "+result+"; expected "+(created ? 0 : expected));
		check(table.getValue(key)==expected, "Count must saturate rather than wrap: "+table.getValue(key));
		check(table.size()+table.victims().size()==61, "Updating a count must preserve all keys");
		for(int i=0; i<61; i++){
			final long other=127L*i;
			check(table.getValue(other)==(other==key ? expected : 1), "Changed count for key "+other);
		}
	}

	/** Keeps the independent oracle enabled even when production assertions are disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
