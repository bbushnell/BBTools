package test.bloom;

import java.util.Arrays;

import bloom.BloomFilter;
import bloom.BloomFilterCorrector;
import bloom.BloomFilterCorrector1;
import bloom.BloomFilterCorrector2;
import bloom.KCountArray;
import bloom.KCountArray7MTA;
import stream.Read;
import structures.LongList;
import ukmer.Kmer;

/** Checks merge boundaries against finite real count tables and uniform support.
 * Run in a fresh JVM; this test sets the count-array seed for reproducibility.
 * @author Brian Bushnell, Sucrose */
public final class BloomCorrectorMergeBoundaryTest{

	/** Tests short and long representations, with either hash mode and boundary controls. */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("This finite test accepts no arguments.");}
		final long oldCells=BloomFilter.OVERRIDE_CELLS;
		final boolean oldMask=Kmer.MASK_CORE;
		int failures=0;
		try{
			check(BloomFilterCorrector2.FORCE_MULT_SPLIT==null,
					"The test requires natural long-kmer geometry, without a diagnostic split override.");
			Kmer.MASK_CORE=false;
			BloomFilter.OVERRIDE_CELLS=2000003;
			for(int k : new int[]{31, 62, 63}){
				for(boolean dual : new boolean[]{false, true}){failures+=exercise(k, dual);}
			}
		}finally{
			BloomFilter.OVERRIDE_CELLS=oldCells;
			Kmer.MASK_CORE=oldMask;
		}
		if(failures>0){
			throw new AssertionError("Uniformly supported merges must not depend on an unrelated key: "+
					failures+" failed cases of 48.");
		}
		System.out.println("BLOOM_CORRECTOR_MERGE_BOUNDARY_TEST_PASS cases=48");
	}

	/** Adds one supported sequence key, then a high-count key absent from that sequence. */
	private static int exercise(final int k, final boolean dual){
		assert(k>=31) : "The fixture requires room for two overlapping mates and four-window boundary checks.";
		KCountArray7MTA.setSeed(20261105);
		final BloomFilter filter=new BloomFilter(null, null, null, k, k, 16, 3, 1, true, false, false, 0.01f, dual);
		final Kmer repeated=new Kmer(k);
		for(int i=0; i<k; i++){repeated.addRight((byte)'C');}
		final long key=(k<=31 ? filter.toKey(repeated.array1()[0], repeated.array2()[0]) : repeated.xor());
		final long key2=(k<=31 ? KCountArray.secondaryHash(key) : repeated.xor2());
		//A short stored zero decodes as all A and canonicalizes to all T. Long
		//stored zeros are already hash values, including the second hash lane.
		final long unrelated=(k<=31 ? (1L<<(2*k))-1 : 0);
		final long unrelated2=(k<=31 ? KCountArray.secondaryHash(unrelated) : 0);
		check(key>=0 && key!=unrelated, "The sequence key must be valid and distinct from the unrelated key.");
		add(filter, key, key2, 5);
		final BloomFilterCorrector corrector=(k<=31 ? new BloomFilterCorrector1(filter, k, k) : new ReadyLongCorrector(filter, k));
		final int windows=2*k+18, length=windows+k-1, width=4;
		final byte[] bases=new byte[length];
		Arrays.fill(bases, (byte)'C');
		final Read merged=new Read(bases, null, "uniform-support", 0, false);
		int failures=0;
		for(int background : new int[]{0, 1000}){
			if(background>0){add(filter, unrelated, unrelated2, background);}
			check(count(filter, key, key2)==5 && count(filter, unrelated, unrelated2)==background,
					"A hash collision or count error would invalidate the independently specified table.");
			//The first case ends exactly at the last valid window. The other
			//boundaries extend beyond it; every existing window still has count 5.
			for(int firstLength : new int[]{windows-width-1, windows-width, windows, windows+2}){
				final int secondLength=length-10;
				final boolean accepted=corrector.mergeOK(merged, firstLength, secondLength, new LongList(), width, 3, 80);
				System.out.println("k="+k+" dual="+dual+" mergedLength="+length+" firstLength="+firstLength+
						" secondLength="+secondLength+" windows="+windows+" unrelatedCount="+background+
						" expected=true actual="+accepted);
				if(!accepted){failures++;}
			}
		}
		return failures;
	}

	/** Inserts through the same real table mode that the corrector will query. */
	private static void add(final BloomFilter filter, final long key, final long key2, final int amount){
		assert(amount>0) : "Only the two explicit positive fixture counts are inserted.";
		if(filter.dualHash){filter.filter.increment(key, key2, amount);}
		else{filter.filter.increment(key, amount);}
	}

	/** Reads back fixture counts so collisions cannot masquerade as boundary failures. */
	private static int count(final BloomFilter filter, final long key, final long key2){
		return filter.dualHash ? filter.filter.read(key, key2) : filter.filter.read(key);
	}

	/** Fails under either assertion mode when the test's prerequisites do not hold. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Exposes only the normal thread-local initialization; all merge logic stays real. */
	private static final class ReadyLongCorrector extends BloomFilterCorrector2{
		ReadyLongCorrector(final BloomFilter filter, final int k){
			super(filter, k, k);
			initializeThreadLocals();
		}
	}
}
