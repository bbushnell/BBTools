package test.kmer;

import kmer.HashArray1D;
import kmer.KmerNode;

/**
 * Checks calcMem's existing approximate byte model on tiny, verified table layouts.
 * Literal totals are controls for arithmetic widening, not measured heap requirements
 * or a large-table overflow reproduction. Run each case with -ea and -da.
 * @author Sucrose
 */
public final class HashArrayMemoryTest{

	/** Accepts exactly one case=empty|empty-owned|victim|victim-owned argument. */
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("case="), "Expected exactly one case= argument");
		final String name=args[0].substring(5);
		final boolean victim, owned;
		final long expected;
		switch(name){
			case "empty": victim=false; owned=false; expected=820; break;
			case "empty-owned": victim=false; owned=true; expected=1064; break;
			case "victim": victim=true; owned=false; expected=868; break;
			case "victim-owned": victim=true; owned=true; expected=1112; break;
			default: throw new IllegalArgumentException("Unknown case: "+name);
		}
		final HashArray1D table=new HashArray1D(1, -1L, false);
		if(victim){
			for(int key=0; key<61; key++){
				check(table.set(key, 1)==1, "Each distinct fixture key must be newly inserted: "+key);
			}
		}
		if(owned){table.initializeOwnership();}
		verifyLayout(table, victim, owned);
		final long actual=table.calcMem();
		final String diagnostic="case="+name+", expectedBytes="+expected+", actualBytes="+actual;
		System.out.println(diagnostic);
		check(actual==expected, "Verified tiny layout must retain its literal byte estimate; "+diagnostic);
		System.out.println("PASS case="+name);
	}

	/** Establishes actual storage lengths and all inserted keys before trusting a byte total. */
	private static void verifyLayout(final HashArray1D table, final boolean victim, final boolean owned){
		final long[] keys=table.array();
		final int[] values=table.values();
		check(table.arrayLength()==61 && keys.length==61 && values.length==61,
			"Fixture requires one hash cell plus 60 probe-tail slots in both primary arrays");
		check(table.victims().arrayLength()==11,
			"Fixture requires the 11-bucket victim forest selected for the minimum requested size");
		check(table.size()==(victim ? 60 : 0) && table.victims().size()==(victim ? 1 : 0),
			"Fixture must have either no keys or exactly 60 primary keys and one victim");
		if(owned){
			check(table.owners()!=null && table.owners().length()==61,
				"Ownership storage must cover every allocated primary slot");
		}else{
			check(table.owners()==null, "Unowned controls must not allocate ownership storage");
		}
		for(int i=0; i<61; i++){
			final long expectedKey=victim && i<60 ? i : -1;
			check(keys[i]==expectedKey, "Unexpected primary key at slot "+i+": "+keys[i]);
			if(victim && i<60){check(values[i]==1, "Unexpected primary count at slot "+i);}
		}
		int nodes=0;
		for(KmerNode node : table.victims()){
			nodes++;
			check(victim && node.getValue(60)==1, "The sole victim must contain key 60 with count 1");
		}
		check(nodes==(victim ? 1 : 0), "Victim traversal must agree with the literal node count");
		if(victim){
			for(int key=0; key<61; key++){
				check(table.getValue(key)==1, "Inserted key must remain retrievable: "+key);
			}
		}
		check(table.getValue(61)==-1, "Key beyond the fixture range must be absent");
	}

	/** Keeps layout and literal-result checks active when Java assertions are disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
