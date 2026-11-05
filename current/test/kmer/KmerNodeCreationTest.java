package test.kmer;

import kmer.KmerNode;
import kmer.KmerNode1D;

/**
 * Checks direct-node creation reports independently of the resulting count.
 * Each case uses at most three scalar nodes and literal key/count expectations.
 * HashForest has its own creation path; these cases exercise the node API directly.
 * Run one case per JVM with either assertion mode. No large allocation or I/O.
 * @author Sucrose
 */
public final class KmerNodeCreationTest{

	/** Selects one affected case or nearby unaffected control.
	 * @param args Exactly one case=name argument from the cases below */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected one case=name argument.");
		}
		final String mode=args[0].substring(5);
		switch(mode){
			case "new-left": child(mode, 20, 10, 2); break;
			case "new-right": child(mode, 60, 80, 2); break;
			case "empty-root": root(mode, -1, 0, 40, 2, 1, 2); break;
			case "existing-zero": root(mode, 40, 0, 40, 1, 0, 1); break;
			case "new-one-control": child(mode, 20, 10, 1); break;
			case "existing-positive-control": root(mode, 40, 9, 40, 2, 0, 11); break;
			default: throw new IllegalArgumentException("Unknown creation case: "+mode);
		}
		System.out.println("PASS KMR008 case="+mode);
	}

	/** Inserts below an existing child, checking recursive insertion and preserved ancestors.
	 * @param mode Diagnostic case name
	 * @param branchKey Existing child key, either 20 or 60, with count 5
	 * @param key Absent descendant key, either 10 or 80
	 * @param incr Positive initial count of the new descendant */
	private static void child(final String mode, final long branchKey, final long key, final int incr){
		final KmerNode1D tree=new KmerNode1D(40, 9);
		check(tree.set(branchKey, 5)==1 && !tree.contains(key), mode, "Fixture requires an existing branch and absent descendant.");
		final int created=tree.incrementAndReturnNumCreated(key, incr);
		check(tree.pivot()==40 && tree.getValue(40)==9 && tree.getValue(branchKey)==5
			&& tree.getValue(key)==incr && tree.size()==3, mode, "Insertion must preserve both ancestors and store the new count.");
		final boolean left=branchKey<40;
		final KmerNode branch=left ? tree.left() : tree.right();
		check(branch!=null && branch.pivot()==branchKey && (left ? tree.right() : tree.left())==null,
			mode, "Only the expected root branch may exist.");
		final KmerNode leaf=left ? branch.left() : branch.right();
		check(leaf!=null && leaf.pivot()==key && leaf.left()==null && leaf.right()==null
			&& (left ? branch.right() : branch.left())==null, mode, "The new key must be the sole descendant leaf.");
		check(created==1, mode, "Expected created=1, got="+created);
	}

	/** Initializes an empty pivot or updates a present key without creating child links.
	 * @param mode Diagnostic case name
	 * @param pivot Initial pivot, -1 for an empty root or 40 for a present key
	 * @param initial Initial scalar count
	 * @param key Key to initialize or update
	 * @param incr Positive increment
	 * @param expectedCreated Literal creation oracle, independent of the final count
	 * @param expectedCount Literal final-count oracle */
	private static void root(final String mode, final long pivot, final int initial, final long key,
			final int incr, final int expectedCreated, final int expectedCount){
		final KmerNode1D tree=new KmerNode1D(pivot, initial);
		check(tree.contains(key)==(pivot==key), mode, "Fixture presence must depend on the key, not its count.");
		final int created=tree.incrementAndReturnNumCreated(key, incr);
		check(tree.pivot()==key && tree.getValue(key)==expectedCount && tree.size()==1
			&& tree.left()==null && tree.right()==null, mode, "Only the expected root key and count may remain.");
		check(created==expectedCreated, mode, "Expected created="+expectedCreated+", got="+created);
	}

	/** Keeps fixture and oracle checks active under both -ea and -da. */
	private static void check(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("KMR008 case="+mode+": "+message);}
	}
}
