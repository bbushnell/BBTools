package test.bloom;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.BitSet;

import bloom.BloomFilterCorrector1;
import structures.IntList;

/** Exhaustive ten-window support patterns checked by an independent base-centric union.
 * Calls only the shared marking helper; the supplied counts do not encode DNA keys.
 * @author Brian Bushnell, Sucrose */
public final class BloomMarkingUnionTest{

	/** Checks overlap chains, adjacency, gaps and final endpoints for seven interval widths. */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("This finite test accepts no arguments.");}
		int cases=0, failures=0;
		for(int k : new int[]{1, 2, 3, 4, 7, 15, 31}){
			final BloomFilterCorrector1 corrector=new BloomFilterCorrector1(null, k, k);
			final BitSet covered=new BitSet();
			for(int pattern=0; pattern<1024; pattern++){
				final IntList counts=new IntList();
				for(int start=0; start<10; start++){counts.add((pattern&(1<<start))==0 ? 0 : 2);}
				for(boolean quality : new boolean[]{false, true}){
					final int length=k+9;
					final byte[] bases=new byte[length], quals=(quality ? new byte[length] : null);
					final BitSet expected=new BitSet();
					int expectedMarked=0;
					for(int base=0; base<length; base++){
						bases[base]=DNA[base&3];
						if(quals!=null){quals[base]=30;}
						boolean supported=false;
						for(int start=0; start<10; start++){
							if((pattern&(1<<start))!=0 && start<=base && base<start+k){supported=true;}
						}
						if(supported){expected.set(base);}else{expectedMarked++;}
					}
					covered.set(0, length+64);//The helper must clear stale scratch, including beyond this read.
					final int marked=corrector.markBadBases(bases, quals, counts, covered, 2, false, (byte)0);
					boolean pass=marked==expectedMarked && covered.equals(expected);
					for(int base=0; base<length; base++){
						pass&=bases[base]==(expected.get(base) ? DNA[base&3] : 'N');
						if(quals!=null){pass&=quals[base]==(expected.get(base) ? 30 : 0);}
					}
					cases++;
					if(!pass){
						failures++;
						if(failures<=16){System.out.println("FAIL k="+k+" pattern="+pattern+" qualities="+quality+
								" expectedMarked="+expectedMarked+" actualMarked="+marked+
								" expectedCoverage="+expected+" actualCoverage="+covered);}
					}
				}
			}
		}
		failures+=policy("middle", 4, "ACGTACG", new byte[]{30, 30, 30, 30, 30, 30, 30}, true,
				new byte[]{30, 10, 30, 30, 30, 10, 30}, 2);
		failures+=policy("tail", 1, "ACGTACG", new byte[]{30, 30, 30, 30, 30, 30, 30}, true,
				new byte[]{30, 30, 30, 10, 30, 30, 10}, 2);
		failures+=policy("interrupted-tail", 1, "ACGTNCG", new byte[]{30, 30, 30, 30, 0, 5, 30}, true,
				new byte[]{30, 30, 30, 10, 0, 5, 30}, 1);
		failures+=policy("all-low", 0, "ANCGTNA", new byte[]{30, 0, 5, 30, 30, 30, 30}, false,
				new byte[]{10, 0, 5, 10, 10, 0, 10}, 6);
		failures+=policy("no-quality", 0, "ANCGTNA", null, false, null, 5);
		failures+=policy("fully-supported", 31, "ANCGTNA", new byte[]{30, 0, 5, 30, 30, 30, 30}, true,
				new byte[]{30, 0, 5, 30, 30, 30, 30}, 0);
		cases+=6;
		if(failures>0){throw new AssertionError("Marking union/policy failed "+failures+" of "+cases+" cases.");}
		System.out.println("BLOOM_MARKING_UNION_TEST_PASS cases="+cases);
	}

	/** Uses literal positive-quality outcomes, including existing N/low-quality reset boundaries. */
	private static int policy(final String name, final int pattern, final String sequence, final byte[] quals,
			final boolean deltaOnly, final byte[] expectedQuals, final int expectedMarked){
		final byte[] bases=sequence.getBytes(StandardCharsets.US_ASCII);
		final IntList counts=new IntList();
		for(int start=0; start<5; start++){counts.add((pattern&(1<<start))==0 ? 0 : 2);}
		final BitSet expected=new BitSet(), actual=new BitSet();
		for(int base=0; base<7; base++){
			for(int start=0; start<5; start++){
				if((pattern&(1<<start))!=0 && start<=base && base<start+3){expected.set(base);}
			}
		}
		final int marked=new BloomFilterCorrector1(null, 3, 3).markBadBases(bases, quals, counts, actual, 2, deltaOnly, (byte)10);
		final boolean pass=marked==expectedMarked && expected.equals(actual) &&
				Arrays.equals(bases, sequence.getBytes(StandardCharsets.US_ASCII)) && Arrays.equals(quals, expectedQuals);
		System.out.println("policy="+name+" expectedMarked="+expectedMarked+" actualMarked="+marked+
				" expectedQualities="+Arrays.toString(expectedQuals)+" actualQualities="+Arrays.toString(quals)+" pass="+pass);
		return pass ? 0 : 1;
	}

	private static final byte[] DNA={'A', 'C', 'G', 'T'};
}
