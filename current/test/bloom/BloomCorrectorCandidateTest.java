package test.bloom;

import java.util.Arrays;

import bloom.BloomFilter;
import bloom.BloomFilterCorrector;
import bloom.BloomFilterCorrector1;
import bloom.BloomFilterCorrector2;
import bloom.KCountArray7MTA;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/** Distinguishes candidate base codes from their supporting counts in directional correction.
 * Real tables support only the all-C window at depth two or three. G/T terminal
 * substitutions must become C; an already-C terminal base must remain unchanged.
 * This is a finite helper-contract test, not end-to-end correction accuracy.
 * @author Brian Bushnell, Sucrose */
public final class BloomCorrectorCandidateTest{

	/** Exercises both code/count coincidences and unaffected controls at short and long k. */
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
			check(BloomFilterCorrector2.FORCE_MULT_SPLIT==null, "Natural long-window geometry is required for this fixture.");
			for(int k : new int[]{31, 62, 63}){
				for(int depth : new int[]{2, 3}){failures+=exercise(k, depth);}
			}
		}finally{
			Shared.setThreads(oldThreads);
			BloomFilter.OVERRIDE_CELLS=oldCells;
			Kmer.MASK_CORE=oldMask;
		}
		if(failures>0){throw new AssertionError("Supported C corrections must not be skipped because G/T codes equal counts: "+failures+" failures of 18 cases.");}
		System.out.println("BLOOM_CORRECTOR_CANDIDATE_TEST_PASS cases=18");
	}

	/** Loads the sole supported window, verifies all candidates, then checks full output and recounts. */
	private static int exercise(final int k, final int depth){
		KCountArray7MTA.setSeed(20261105);
		final BloomFilter filter=new BloomFilter(null, null, null, k, k, 4, 3, 1, true, false, false, 0.01f, false);
		final Kmer repeated=new Kmer(k);
		for(int i=0; i<k; i++){repeated.addRight((byte)'C');}
		final long key=(k<=31 ? filter.toKey(repeated.array1()[0], repeated.array2()[0]) : repeated.xor());
		filter.filter.increment(key, depth);
		final BloomFilterCorrector corrector=(k<=31 ? new BloomFilterCorrector1(filter, k, k) : new ReadyLongCorrector(filter, k));
		final byte[] bases=new byte[k+1];
		Arrays.fill(bases, (byte)'C');
		final LongList kmers=new LongList();
		//The only positive extension must be C. Verify all four real table queries
		//before invoking the candidate selector, so collisions cannot fake the defect.
		for(byte candidate : new byte[]{'A', 'C', 'G', 'T'}){
			bases[k]=candidate;
			check(corrector.fillKmers(bases, kmers)==2 && kmers.size()==2, "Each input must have exactly two valid windows.");
			check(corrector.getCount2(kmers.get(0))==depth && corrector.getCount2(kmers.get(1))==(candidate=='C' ? depth : 0),
					"A collision or encoding mismatch invalidates the literal candidate table.");
		}
		int failures=0;
		for(byte candidate : new byte[]{'C', 'G', 'T'}){
			bases[k]=candidate;
			final ByteBuilder bb=new ByteBuilder(k+1);
			bb.append(bases);
			check(bb.array.length==bb.length() && bb.length()==k+1,
					"Exact backing capacity isolates candidate selection from the separate recount capacity contract.");
			final IntList counts=new IntList();
			counts.add(depth);
			counts.add(candidate=='C' ? depth : 0);
			final int expectedCorrections=(candidate=='C' ? 0 : 1);
			final int corrected=corrector.reassemble_inner(bb, null, new int[4], counts, 0);
			boolean pass=corrected==expectedCorrections && bb.length()==k+1 && counts.size==2 &&
					counts.get(0)==depth && counts.get(1)==depth;
			for(int i=0; i<bb.length(); i++){pass&=bb.array[i]=='C';}
			System.out.println("k="+k+" depth="+depth+" originalLast="+(char)candidate+
					" expectedLast=C actualLast="+(char)bb.array[k]+" expectedCorrections="+expectedCorrections+
					" actualCorrections="+corrected+" counts="+counts+" pass="+pass);
			if(!pass){failures++;}
		}
		return failures;
	}

	/** Keeps fixture prerequisites active under either assertion mode. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Initializes normal long-corrector scratch without replacing any production algorithm. */
	private static final class ReadyLongCorrector extends BloomFilterCorrector2{
		ReadyLongCorrector(final BloomFilter filter, final int k){super(filter, k, k); initializeThreadLocals();}
	}
}
