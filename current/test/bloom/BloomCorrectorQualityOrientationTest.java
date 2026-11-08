package test.bloom;

import java.util.Arrays;
import java.util.BitSet;

import assemble.ErrorTracker;
import bloom.BloomFilter;
import bloom.BloomFilterCorrector;
import bloom.BloomFilterCorrector1;
import bloom.BloomFilterCorrector2;
import bloom.KCountArray7MTA;
import dna.AminoAcid;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/** Tests that reverse proposal scanning uses qualities aligned to its bases (BLM212).
 * An all-C window has depth 100 and a G-substituted endpoint window has depth 6.
 * Default thresholds accept the substitution at quality 0, reject it at quality 40,
 * and accept absent qualities (the scan uses 20). Uniform and opposite-end controls
 * separate quality orientation from candidate selection. This is a finite helper test.
 * @author Brian Bushnell, Sucrose */
public final class BloomCorrectorQualityOrientationTest{

	/** Runs ten cases per geometry with prerequisites active in either assertion mode. */
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
			check(BloomFilterCorrector2.FORCE_MULT_SPLIT==null, "Natural k-mer geometry is required.");
			for(int k : new int[]{31, 62, 63}){
				for(boolean first : new boolean[]{true, false}){
					failures+=exercise(k, first, 0, 40);
					failures+=exercise(k, first, 40, 0);
					failures+=exercise(k, first, 0, 0);
					failures+=exercise(k, first, 40, 40);
					failures+=exercise(k, first, -1, -1);
				}
			}
		}finally{
			Shared.setThreads(oldThreads);
			BloomFilter.OVERRIDE_CELLS=oldCells;
			Kmer.MASK_CORE=oldMask;
		}
		if(failures>0){throw new AssertionError("Directional quality alignment: "+failures+" failures of 30 cases.");}
		System.out.println("BLOOM_CORRECTOR_QUALITY_ORIENTATION_TEST_PASS cases=30");
	}

	/** Builds a literal table, validates all endpoint alternatives, then checks one complete pass. */
	private static int exercise(final int k, final boolean first, final int errorQuality, final int oppositeQuality){
		KCountArray7MTA.setSeed(20261105);
		final BloomFilter filter=new BloomFilter(null, null, null, k, k, 8, 3, 1, true, false, false, 0.01f, false);
		final BloomFilterCorrector corrector=(k<=31 ? new BloomFilterCorrector1(filter, k, k) : new ReadyLongCorrector(filter, k));
		final byte[] bases=new byte[k+1];
		Arrays.fill(bases, (byte)'C');
		final int pos=(first ? 0 : k), alteredWindow=(first ? 0 : 1);
		filter.filter.increment(key(filter, bases, 0, k), 100);
		bases[pos]='G';
		filter.filter.increment(key(filter, bases, alteredWindow, k), 6);
		final LongList kmers=new LongList();
		for(byte candidate : new byte[]{'A', 'C', 'G', 'T'}){
			bases[pos]=candidate;
			check(corrector.fillKmers(bases, kmers)==2 && kmers.size()==2, "Exactly two valid windows required.");
			final int expected=(candidate=='C' ? 100 : candidate=='G' ? 6 : 0);
			check(corrector.getCount2(kmers.get(alteredWindow))==expected &&
					corrector.getCount2(kmers.get(1-alteredWindow))==100, "Collision or incorrect forward table encoding.");
			final byte[] reverse=bases.clone();
			AminoAcid.reverseComplementBasesInPlace(reverse);
			check(corrector.fillKmers(reverse, kmers)==2 && kmers.size()==2, "Exactly two reverse windows required.");
			check(corrector.getCount2(kmers.get(1-alteredWindow))==expected &&
					corrector.getCount2(kmers.get(alteredWindow))==100, "Collision or incorrect reverse table encoding.");
		}
		bases[pos]='G';
		final byte[] quals=(errorQuality<0 ? null : new byte[k+1]);
		if(quals!=null){
			Arrays.fill(quals, (byte)40);
			quals[pos]=(byte)errorQuality;
			quals[k-pos]=(byte)oppositeQuality;
		}
		final int expectedCorrections=(errorQuality==40 ? 0 : 1);
		final byte[] expectedBases=bases.clone(), expectedQuals=(quals==null ? null : quals.clone());
		if(expectedCorrections==1){
			expectedBases[pos]='C';
			if(expectedQuals!=null){expectedQuals[pos]=20;}
		}
		final IntList counts=new IntList();
		counts.add(first ? 6 : 100);
		counts.add(first ? 100 : 6);
		final ErrorTracker tracker=new ErrorTracker();
		final ByteBuilder left=new ByteBuilder(k+1), right=new ByteBuilder(k+1);
		check(left.array.length==k+1 && right.array.length==k+1, "Exact proposal capacity isolates orientation from recount bounds.");
		final int corrected=corrector.reassemble_pass(bases, quals, left, right, new int[4], counts,
				new IntList(), tracker, 0, null, null, new BitSet());
		final boolean pass=corrected==expectedCorrections && Arrays.equals(bases, expectedBases) &&
				Arrays.equals(quals, expectedQuals) && counts.size()==2 && counts.get(1-alteredWindow)==100 &&
				counts.get(alteredWindow)==(expectedCorrections==1 ? 100 : 6) &&
				tracker.corrected()==expectedCorrections && tracker.correctedReassembleOuter==expectedCorrections &&
				tracker.correctedReassembleInner==0 && tracker.detectedReassemble==expectedCorrections && !tracker.rollback;
		System.out.println("k="+k+" end="+(first ? "first" : "last")+" errorQuality="+errorQuality+
				" oppositeQuality="+oppositeQuality+" expectedCorrections="+expectedCorrections+
				" actualCorrections="+corrected+" actualBase="+(char)bases[pos]+" counts="+counts+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Encodes one complete window using the filter's canonical single-hash representation. */
	private static long key(final BloomFilter filter, final byte[] bases, final int start, final int k){
		check(start>=0 && start+k<=bases.length, "Table fixture window must fit its literal sequence.");
		final Kmer kmer=new Kmer(k);
		for(int i=start; i<start+k; i++){kmer.addRight(bases[i]);}
		return k<=31 ? filter.toKey(kmer.array1()[0], kmer.array2()[0]) : kmer.xor();
	}

	/** Retains fixture prerequisites under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Initializes existing long-corrector scratch without replacing its algorithms. */
	private static final class ReadyLongCorrector extends BloomFilterCorrector2{
		ReadyLongCorrector(final BloomFilter filter, final int k){super(filter, k, k); initializeThreadLocals();}
	}
}
