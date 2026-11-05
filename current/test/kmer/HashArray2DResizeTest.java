package test.kmer;

import java.util.Arrays;

import kmer.HashArray2D;

/**
 * Finite value-preservation checks for multi-value table resizing.
 * Run each case in a fresh JVM with assertions enabled and disabled.
 * Uses public APIs only; no production fields or resize thresholds are changed.
 * @author Sucrose
 */
public final class HashArray2DResizeTest {

	/** Selects one independently runnable case; every check remains active under -da. */
	public static void main(String[] args){
		if(args.length!=1){throw new IllegalArgumentException("Expected case=resize-victims, resize-distinct, or fixed-victims");}
		final String mode=args[0];
		final boolean resize, collide;
		if(mode.equals("case=resize-victims")){resize=true; collide=true;}
		else if(mode.equals("case=resize-distinct")){resize=true; collide=false;}
		else if(mode.equals("case=fixed-victims")){resize=false; collide=true;}
		else{throw new IllegalArgumentException("Unknown case: "+mode);}
		run(resize, collide);
		System.out.println("PASS "+mode+" keys=112 valuesPerKey=3");
	}

	/**
	 * 257*i has distinct residues modulo 127 for i=0..111, but residue zero
	 * modulo 257. Inserting key 112 crosses floor(127*0.88)=111 and resizes.
	 * With the current 60-cell probe window, 52 keys must move to victim trees.
	 * The fixed-size control creates the same victim layout without migration;
	 * the distinct-key control resizes without needing victim trees.
	 */
	private static void run(final boolean resize, final boolean collide){
		final int[] schedule=resize ? new int[] {127, 257} : new int[] {257};
		final HashArray2D table=new HashArray2D(schedule, Long.MAX_VALUE);
		final int stride=collide ? 257 : 1;
		for(int i=0; i<111; i++){insert(table, stride, i);}
		check(table.arrayLength()==(resize ? 187 : 317), "Unexpected initial capacity");
		check(table.size()+table.victims().size()==111, "Initial key count must be 111");
		check(table.victims().size()==(!resize && collide ? 51 : 0), "Unexpected initial overflow count");
		for(int i=0; i<111; i++){checkValues(table, (long)stride*i, i, false);}

		insert(table, stride, 111);
		check(table.arrayLength()==317, "Resize must finish at modulus 257 plus 60 probe cells");
		check(table.size()+table.victims().size()==112, "Resize must preserve all 112 keys");
		check(table.victims().size()==(collide ? 52 : 0), "Fixture did not exercise the expected victim layout");
		for(int i=0; i<112; i++){checkValues(table, (long)stride*i, i, false);}

		// Append to every key, including every migrated victim. A readable borrowed
		// array alone does not prove that a node retained its logical value count.
		for(int i=0; i<112; i++){
			final long key=(long)stride*i;
			check(table.set(key, 2000+i)==0, "Append must not create another key: "+key);
			check(table.set(key, 1000+i)==0, "Duplicate must not create another key: "+key);
			check(table.set(key, 2000+i)==0, "Repeated append must retain the same key: "+key);
			checkValues(table, key, i, true);
		}
		check(table.size()+table.victims().size()==112, "Updates must preserve the key count");
		check(!table.contains(1000000L), "Absent key unexpectedly present");
	}

	/** Inserts a fresh packed two-value list with unused capacity and a known terminator. */
	private static void insert(final HashArray2D table, final int stride, final int i){
		final long key=(long)stride*i;
		check(table.set(key, new int[] {10+i, 1000+i, -1, -1}, 2)==1,
			"Fresh key must be counted exactly once: "+key);
	}

	/** Checks the literal input values and appended value, allowing only negative padding afterward. */
	private static void checkValues(final HashArray2D table, final long key, final int i, final boolean appended){
		check(table.contains(key), "Lost key "+key);
		final int[] actual=table.getValues(key, new int[1]);
		final int count=appended ? 3 : 2;
		check(actual!=null && actual.length>=count, "Missing values for key "+key);
		check(actual[0]==10+i && actual[1]==1000+i,
			"Changed input values for key "+key+": "+Arrays.toString(actual));
		if(appended){check(actual[2]==2000+i, "Lost appended value for key "+key+": "+Arrays.toString(actual));}
		for(int j=count; j<actual.length; j++){
			check(actual[j]<0, "Unexpected extra value for key "+key+": "+Arrays.toString(actual));
		}
	}

	/** Throws in both assertion modes so disabled Java assertions cannot disable the oracle. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
