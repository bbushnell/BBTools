package test.bloom;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.BitSet;

import bloom.BloomFilterCorrector1;
import structures.IntList;

/** Finite union-of-supported-windows oracle for the base-marking helper.
 * Tests the supplied-count contract, not Bloom-filter accuracy or CLI reachability.
 * @author Brian Bushnell, Sucrose */
public final class BloomCorrectorMarkingTest{

	/** Exercises every support pattern for three overlapping length-three windows. */
	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("This finite test accepts no arguments.");}
		final BloomFilterCorrector1 corrector=new BloomFilterCorrector1(null, 3, 3);
		int cases=0, failures=0;
		for(int pattern=0; pattern<8; pattern++){
			for(int qualityMode=0; qualityMode<2; qualityMode++){
				final byte[] bases="ACGTA".getBytes(StandardCharsets.US_ASCII);
				final byte[] quals=(qualityMode==0 ? null : new byte[]{20, 21, 22, 23, 24});
				final byte[] expected=bases.clone();
				final byte[] expectedQuals=(quals==null ? null : quals.clone());
				final IntList counts=new IntList();
				for(int start=0; start<3; start++){counts.add((pattern&(1<<start))==0 ? 0 : 2);}
				int expectedMarked=0;
				//Independent base-centric oracle: test every covering window, without
				//skipping windows or sharing the implementation's BitSet traversal.
				for(int base=0; base<expected.length; base++){
					boolean supported=false;
					for(int start=0; start<3; start++){
						if((pattern&(1<<start))!=0 && start<=base && base<start+3){supported=true;}
					}
					if(!supported){
						expected[base]='N';
						if(expectedQuals!=null){expectedQuals[base]=0;}
						expectedMarked++;
					}
				}
				final int marked=corrector.markBadBases(bases, quals, counts, new BitSet(), 2, false, (byte)0);
				final boolean pass=marked==expectedMarked && Arrays.equals(bases, expected)
						&& Arrays.equals(quals, expectedQuals);
				System.out.println("pattern="+pattern+" qualities="+(quals!=null)+" expectedMarked="+expectedMarked+
						" marked="+marked+" expected="+new String(expected, StandardCharsets.US_ASCII)+
						" actual="+new String(bases, StandardCharsets.US_ASCII)+" pass="+pass);
				cases++;
				if(!pass){failures++;}
			}
		}
		if(failures>0){
			throw new AssertionError("markBadBases must preserve every base covered by any supported window: "+
					failures+" failed cases of "+cases);
		}
		System.out.println("BLOOM_CORRECTOR_MARKING_TEST_PASS cases="+cases);
	}
}
