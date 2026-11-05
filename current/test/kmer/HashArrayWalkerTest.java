package test.kmer;

import kmer.HashArray1D;
import kmer.Walker;

/**
 * Checks that a short-key walker emits exactly the two live key/count pairs.
 * The regenerated case removes a low-count colliding key before walking; the
 * fresh control inserts only the same two survivors. No backing arrays are edited
 * by the fixture, and expected records are literal rather than derived from lookup.
 * Run one case per JVM, with either assertion mode. No large allocation or I/O.
 * @author Sucrose
 */
public final class HashArrayWalkerTest{

	/** Builds one finite table and checks every emitted record without filtering it.
	 * @param args Exactly case=regenerated or case=fresh-control */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected case=regenerated|fresh-control.");
		}
		final String mode=args[0].substring(5);
		final boolean regenerated=mode.equals("regenerated");
		if(!regenerated && !mode.equals("fresh-control")){
			throw new IllegalArgumentException("Unknown walker case: "+mode);
		}

		final HashArray1D table=new HashArray1D(new int[] {17}, -1L);
		check(table.arrayLength()==77, mode, "Fixture requires 17 primary cells and the 60-cell probe tail.");
		check(table.kmerToCell(1)==1 && table.kmerToCell(18)==1 && table.kmerToCell(35)==1,
			mode, "Fixture keys must share the same initial probe cell.");
		if(regenerated){check(table.set(1, 1)==1, mode, "Fixture low-count key must be new.");}
		check(table.set(18, 3)==1, mode, "Fixture first survivor must be new.");
		check(table.set(35, 5)==1, mode, "Fixture second survivor must be new.");
		if(regenerated){
			check(table.regenerate(1)==1, mode, "Fixture regeneration must remove exactly the count-1 key.");
		}
		check(table.arrayLength()==77 && table.size()==2 && table.victims().size()==0,
			mode, "Fixture must retain exactly two primary entries without growth or victims.");
		check(!table.contains(1) && table.getValue(18)==3 && table.getValue(35)==5,
			mode, "Fixture lookup state must contain only the two expected survivors.");

		verifyWalk(table.walk(), mode);
		System.out.println("PASS KMR004 case="+mode+" records=2");
	}

	/** Requires the literal unordered pairs (18,3) and (35,5), exactly once each.
	 * Nonpositive or unknown entries fail immediately; the test does not apply a
	 * consumer-side filter that would hide the post-regeneration empty-cell bug.
	 * @param walker Fresh cursor over the stable fixture table
	 * @param mode Case name included in failure diagnostics */
	private static void verifyWalk(final Walker walker, final String mode){
		int seen=0, records=0;
		while(walker.next()){
			final long key=walker.kmer();
			final int value=walker.value();
			check(key>=0 && value>0, mode, "Walker emitted an empty cell: key="+key+", value="+value);
			final int bit=key==18 ? 1 : key==35 ? 2 : 0;
			check(bit!=0, mode, "Walker emitted an unexpected key: "+key);
			check((seen&bit)==0, mode, "Walker emitted a duplicate key: "+key);
			check(value==(key==18 ? 3 : 5), mode, "Walker emitted the wrong count for key "+key+": "+value);
			seen|=bit;
			records++;
		}
		check(seen==3 && records==2, mode, "Walker lost a survivor: seen="+seen+", records="+records);
	}

	/** Keeps all oracle checks active under both -ea and -da. */
	private static void check(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("KMR004 case="+mode+": "+message);}
	}
}
