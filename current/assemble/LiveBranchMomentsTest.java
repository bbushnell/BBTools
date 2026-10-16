package assemble;

import java.util.Random;

import ukmer.HashArrayU1D;
import ukmer.Kmer;
import ukmer.WalkerU;

/** Exact count-only/walker parity, including overflowing hashes and large squared counts. @author Fischl */
public final class LiveBranchMomentsTest {

	/** Exercises word boundaries without relying on the model's narrower training domain. */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("This deterministic fixture takes no arguments.");}
		for(int k : new int[]{32, 33, 40, 62, 64, 65, 95, 96, 124, 127, 128, 300}){test(k);}
		System.err.println("LIVE_BRANCH_MOMENTS_TEST_PASS checks="+checks);
	}

	/** A fixed small table forces victims; counts include INT_MAX to expose summation-order changes. */
	private static void test(final int k){
		final Kmer key=new Kmer(k);
		final HashArrayU1D table=new HashArrayU1D(new int[]{31}, key.k, key.kbig);
		compare(table, 0);
		final Random random=new Random(7030+k);
		for(int i=0; i<500; i++){
			key.clear();
			for(int j=0; j<k; j++){key.addRightNumeric(random.nextInt(4));}
			table.set(key, i%3==0 ? Integer.MAX_VALUE : 1+random.nextInt(1000000000));
		}
		check(table.victims().size()>0, "The fixed-size fixture must exercise collision victims.");
		compare(table, 500);
		final int[] values=table.values();
		int cell=0;
		while(cell<values.length && values[cell]<=0){cell++;}
		check(cell<values.length, "The fixture must contain a main-array entry.");
		final int saved=values[cell];
		// Count values are authoritative; empty/deleted cells must not leak into sums.
		values[cell]=0; compare(table, 499);
		values[cell]=-1; compare(table, 499);
		values[cell]=saved; compare(table, 500);
	}

	/** Compares raw double bits, not only rounded model features or an approximate tolerance. */
	private static void compare(final HashArrayU1D table, final long expected){
		final SameKJoinFeatures.Moments direct=new SameKJoinFeatures.Moments();
		SameKJoinFeatures.addLongCounts(direct, table);
		final SameKJoinFeatures.Moments walked=new SameKJoinFeatures.Moments();
		final WalkerU walker=table.walk();
		while(walker.next()){if(walker.value()>0){walked.add(walker.value());}}
		check(direct.unique==expected && direct.unique==walked.unique, "Count-only population differs from the key-aware walker.");
		check(Double.doubleToRawLongBits(direct.observations)==Double.doubleToRawLongBits(walked.observations),
				"Observation summation order changed.");
		check(Double.doubleToRawLongBits(direct.squares)==Double.doubleToRawLongBits(walked.squares),
				"Squared-count summation order changed.");
	}

	/** Keeps fixture checks active even if assertions are accidentally disabled. */
	private static void check(final boolean pass, final String message){
		checks++;
		if(!pass){throw new AssertionError(message);}
	}

	private static int checks;
}
