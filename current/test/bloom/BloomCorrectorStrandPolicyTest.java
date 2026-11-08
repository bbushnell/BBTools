package test.bloom;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import bloom.BloomFilter;
import bloom.BloomFilterCorrector1;
import bloom.KCountArray;
import bloom.KCountArray7MTA;
import structures.ByteBuilder;
import structures.LongList;

/** Compares short-corrector lookup APIs against explicit strand-policy counts.
 * Run in a fresh JVM because this test selects a deterministic count-array seed.
 * @author Brian Bushnell, Sucrose */
public final class BloomCorrectorStrandPolicyTest{

	/** Covers canonical/forward-only tables, both orientations, both hash modes and k boundaries. */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("This finite test accepts no arguments.");}
		final long oldCells=BloomFilter.OVERRIDE_CELLS;
		int failures=0;
		try{
			BloomFilter.OVERRIDE_CELLS=2000003;
			final BloomFilterCorrector1 unattached=new BloomFilterCorrector1(null, 3, 3);
			check(unattached.toValue(6, 27)==27 && unattached.toValue(27, 6)==27,
					"Key selection without an attached filter must retain its existing canonical behavior.");
			for(int k : new int[]{3, 31}){
				for(boolean canonical : new boolean[]{false, true}){
					for(boolean dual : new boolean[]{false, true}){
						for(boolean storeReverse : new boolean[]{false, true}){
							failures+=exercise(k, canonical, dual, storeReverse);
						}
					}
				}
			}
			for(boolean canonical : new boolean[]{false, true}){
				for(boolean dual : new boolean[]{false, true}){failures+=leftExtension(canonical, dual);}
			}
		}finally{BloomFilter.OVERRIDE_CELLS=oldCells;}
		if(failures>0){
			throw new AssertionError("Both lookup APIs must honor the table's strand policy: "+
					failures+" failed cases of 36.");
		}
		System.out.println("BLOOM_CORRECTOR_STRAND_POLICY_TEST_PASS strandCases=32 leftExtensionCases=4");
	}

	/** Loads five copies of one literal strand and queries it and its reverse complement. */
	private static int exercise(final int k, final boolean canonical, final boolean dual, final boolean storeReverse){
		assert(k==3 || k==31) : "The literal fixtures cover a small k and the largest short packed k.";
		//Forward sequence is A repeated k-2 times, then CG: packed value 6.
		//Its reverse complement is CG followed by k-2 T bases, encoded directly
		//here without calling the production reverse-complement/canonical helpers.
		final long forward=6;
		long reverse=6;
		for(int i=2; i<k; i++){reverse=(reverse<<2)|3;}
		check(reverse>forward, "The two literal strands must have distinct ordered packed keys.");
		final long populated=(canonical || storeReverse ? reverse : forward);
		KCountArray7MTA.setSeed(20261105);
		final BloomFilter filter=new BloomFilter(null, null, null, k, k, 8, 3, 1, canonical, false, false, 0.01f, dual);
		if(dual){filter.filter.increment(populated, KCountArray.secondaryHash(populated), 5);}
		else{filter.filter.increment(populated, 5);}
		check(rawCount(filter, populated)==5 && rawCount(filter, populated==forward ? reverse : forward)==0,
				"A hash collision or count error would invalidate the literal two-key fixture.");
		final BloomFilterCorrector1 corrector=new BloomFilterCorrector1(filter, k, k);
		int failures=0;
		for(boolean queryReverse : new boolean[]{false, true}){
			final long key=(queryReverse ? reverse : forward), mateKey=(queryReverse ? forward : reverse);
			final int expected=(canonical || queryReverse==storeReverse ? 5 : 0);
			final int pair=corrector.getCount(key, mateKey);
			final int stored=corrector.getCount2(key);
			final long selected=corrector.toValue(key, mateKey), expectedKey=(canonical ? reverse : key);
			final boolean pass=pair==expected && stored==expected && selected==expectedKey;
			System.out.println("k="+k+" canonical="+canonical+" dual="+dual+" storeReverse="+storeReverse+
					" queryReverse="+queryReverse+" expected="+expected+" pair="+pair+" stored="+stored+
					" expectedKey="+expectedKey+" selectedKey="+selected+" pass="+pass);
			if(!pass){failures++;}
		}
		return failures;
	}

	/** Checks both left-lookup callers with literal ACG anchor, AAC predecessor and CGA successor. */
	private static int leftExtension(final boolean canonical, final boolean dual){
		KCountArray7MTA.setSeed(20261105);
		final BloomFilter filter=new BloomFilter(null, null, null, 3, 3, 8, 3, 1, canonical, false, false, 0.01f, dual);
		//Packed forward/reverse pairs: ACG/CGT=6/27, AAC/GTT=1/47, CGA/TCG=24/54.
		for(long key : canonical ? new long[]{27, 47, 54} : new long[]{6, 1, 24}){
			if(dual){filter.filter.increment(key, KCountArray.secondaryHash(key), 5);}
			else{filter.filter.increment(key, 5);}
		}
		//A/C/G/T prepended to ACG produce AAC/CAC/GAC/TAC, whose literal
		//canonical keys are 47/46/45/49; only AAC is supported in this fixture.
		final long[] leftKeys=(canonical ? new long[]{47, 46, 45, 49} : new long[]{1, 17, 33, 49});
		for(int i=0; i<4; i++){
			check(rawCount(filter, leftKeys[i])==(i==0 ? 5 : 0), "Left-extension fixture has a collision or wrong count.");
		}
		final LeftProbe corrector=new LeftProbe(filter);
		final byte[] bases="ACG".getBytes(StandardCharsets.US_ASCII);
		final ByteBuilder bb=new ByteBuilder();
		bb.append(bases);
		final int[] left=new int[4], right=new int[4];
		final int maximum=corrector.leftCount(bases);
		final int appended=corrector.extendToRight2(bb, left, right, 0, false);
		final boolean pass=maximum==5 && Arrays.equals(left, new int[]{5, 0, 0, 0}) && appended==0 &&
				bb.length()==3 && bb.array[0]=='A' && bb.array[1]=='C' && bb.array[2]=='G';
		System.out.println("leftExtension canonical="+canonical+" dual="+dual+" expectedMaximum=5 maximum="+maximum+
				" expectedLeft=[5, 0, 0, 0] left="+Arrays.toString(left)+" appended="+appended+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Exposes the inherited left-edge query without overriding its lookup implementation. */
	private static final class LeftProbe extends BloomFilterCorrector1{
		LeftProbe(final BloomFilter filter){super(filter, 3, 3);}
		int leftCount(final byte[] bases){
			final LongList kmers=new LongList();
			kmers.add(6);
			return maxLeftCount(bases, kmers, 0);
		}
	}

	/** Queries the raw table without applying a strand policy. */
	private static int rawCount(final BloomFilter filter, final long key){
		return filter.dualHash ? filter.filter.read(key, KCountArray.secondaryHash(key)) : filter.filter.read(key);
	}

	/** Rejects invalid fixture setup under either assertion mode. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
