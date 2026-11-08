package test.bloom;

import java.util.Arrays;

import bloom.BloomFilter;
import bloom.BloomFilterCorrector1;
import bloom.KCountArray7MTA;
import shared.Shared;
import structures.LongList;
import ukmer.Kmer;

/** Finite public-API check for saturated-count overflow in hasErrorsFast (BLM213).
 * Both near-maximum counts followed by depth ten are abrupt drops; a flat ten/ten
 * profile is the unaffected control. k=2 exercises the sampled-loop comparison,
 * while k=31 reaches the final-window comparison. No live production caller is claimed.
 * @author Brian Bushnell, Sucrose */
public final class BloomCorrectorFastOverflowTest{

	/** Runs six real-table cases without replacing any correction/counting implementation. */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("This finite test accepts no arguments.");}
		final int oldThreads=Shared.threads();
		final long oldCells=BloomFilter.OVERRIDE_CELLS;
		final boolean oldMask=Kmer.MASK_CORE;
		int failures=0;
		try{
			Shared.setThreads(1);
			BloomFilter.OVERRIDE_CELLS=2000003;
			Kmer.MASK_CORE=false;
			for(int k : new int[]{2, 31}){
				failures+=exercise(k, Integer.MAX_VALUE, true);
				failures+=exercise(k, Integer.MAX_VALUE-1, true);
				failures+=exercise(k, 10, false);
			}
		}finally{
			Shared.setThreads(oldThreads);
			BloomFilter.OVERRIDE_CELLS=oldCells;
			Kmer.MASK_CORE=oldMask;
		}
		if(failures>0){throw new AssertionError("Fast saturated-count detection: "+failures+" failures of 6 cases.");}
		System.out.println("BLOOM_CORRECTOR_FAST_OVERFLOW_TEST_PASS cases=6");
	}

	/** Loads two literal window depths and checks the result and read-only contract. */
	private static int exercise(final int k, final int high, final boolean expected){
		KCountArray7MTA.setSeed(20261105);
		final BloomFilter filter=new BloomFilter(null, null, null, k, k, 32, 3, 1, true, false, false, 0.01f, false);
		check(filter.filter.maxValue==Integer.MAX_VALUE, "Fixture needs supported 32-bit counter saturation.");
		final BloomFilterCorrector1 corrector=new BloomFilterCorrector1(filter, k, k);
		final byte[] bases=new byte[k+1];
		Arrays.fill(bases, (byte)'C');
		bases[k]='A';
		final Kmer rolling=new Kmer(k);
		for(int i=0; i<k; i++){rolling.addRight(bases[i]);}
		final long first=filter.toKey(rolling.array1()[0], rolling.array2()[0]);
		rolling.addRight(bases[k]);
		final long last=filter.toKey(rolling.array1()[0], rolling.array2()[0]);
		check(first!=last, "Literal windows must have distinct canonical keys.");
		filter.filter.increment(first, high);
		filter.filter.increment(last, 10);
		final LongList kmers=new LongList();
		check(corrector.fillKmers(bases, kmers)==2 && kmers.size()==2, "Input has exactly two valid windows.");
		final long storedFirst=kmers.get(0), storedLast=kmers.get(1);
		check(corrector.getCount2(storedFirst)==high && corrector.getCount2(storedLast)==10,
				"Collision or saturation mismatch invalidates the literal profile.");
		final boolean observed=corrector.hasErrorsFast(kmers);
		final boolean unchanged=kmers.size()==2 && kmers.get(0)==storedFirst && kmers.get(1)==storedLast &&
				corrector.getCount2(storedFirst)==high && corrector.getCount2(storedLast)==10;
		final boolean pass=observed==expected && unchanged;
		System.out.println("k="+k+" high="+high+" low=10 expected="+expected+" actual="+observed+" unchanged="+unchanged+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Keeps fixture prerequisites active with assertions enabled or disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
