package test.bloom;

import java.util.Random;

import bloom.KCountArray;
import bloom.KCountArray7MTA;
import dna.AminoAcid;
import ukmer.Kmer;

/** Regression checks for independent secondary hashes in Bloom-filter lanes. */
public class KCountArrayDualHashTest{

	/** Runs all finite dual-hash lane checks. */
	public static void main(String[] args){
		testZeroSecondaryParity();
		testIndependentLanes();
		testPrimitiveKeys();
		testPrefilter();
		testReverseComplements();
		testDistributionAndFalsePositives();
		System.out.println("KCountArrayDualHashTest PASS: "+checks+" checks");
	}

	/** Verifies that xor2=0 preserves the legacy single-key lane sequence. */
	private static void testZeroSecondaryParity(){
		final KCountArray7MTA array=new KCountArray7MTA(10007, 8, 4, null, 0);
		final long[] keys={0, 1, -1L, Long.MIN_VALUE, Long.MAX_VALUE,
				0x5555555555555555L, 0xAAAAAAAAAAAAAAAAL, 0x0123456789ABCDEFL};
		for(int round=0; round<3; round++){
			for(long key : keys){
				array.incrementDual(key, 0);
			}
		}
		for(long key : keys){
			check(array.read(key)==array.read(key, 0), "xor2=0 count parity failed for "+key);
			long oldLane=array.diagnosticHash(key, 0);
			check(oldLane==array.hashForLane(key, 0, 0), "xor2=0 lane parity failed for "+key+", lane 0");
			for(int lane=1; lane<4; lane++){
				oldLane=array.diagnosticHash(Long.rotateRight(oldLane, 6), lane);
				check(oldLane==array.hashForLane(key, 0, lane),
						"xor2=0 lane parity failed for "+key+", lane "+lane);
			}
		}
	}

	/** Verifies that primary and secondary keys control disjoint lane subsets. */
	private static void testIndependentLanes(){
		final KCountArray7MTA array=new KCountArray7MTA(1000003, 8, 4, null, 0);
		final long xor1=0x123456789ABCDEFL;
		final long xor2a=0x0FEDCBA987654321L;
		array.incrementDual(xor1, xor2a);
		check(array.read(xor1, xor2a)==1, "Inserted dual-hash key was absent");

		long xor2b=1;
		while(xor2b==xor2a || array.read(xor1, xor2b)>0){xor2b++;}
		check(array.hashForLane(xor1, xor2a, 0)==array.hashForLane(xor1, xor2b, 0),
				"Equal primary hashes did not share lane zero");
		boolean separated=false;
		for(int lane=1; lane<4; lane++){
			separated|=array.hashForLane(xor1, xor2a, lane)!=array.hashForLane(xor1, xor2b, lane);
		}
		check(separated, "Different secondary hashes did not separate later lanes");
		check(array.read(xor1, xor2b)==0, "Primary collision leaked through every secondary lane");

		final long xor1b=xor1+1;
		check(array.hashForLane(xor1, xor2a, 0)!=array.hashForLane(xor1b, xor2a, 0),
				"Different primary hashes did not separate lane zero");
		check(array.read(xor1b, xor2a)==0, "Secondary collision overrode the primary hash");
	}

	/** Verifies that dual-key updates must satisfy both prefilter hashes. */
	private static void testPrefilter(){
		final long key=17, key2=29;
		final KCountArray7MTA pre=new KCountArray7MTA(1009, 8, 3, null, 0);
		pre.increment(key, key2, 2);
		final KCountArray7MTA main=new KCountArray7MTA(1009, 8, 3, pre, 2);
		main.incrementDual(key, key2);
		check(main.read(key, key2)==1, "Dual hashes did not pass consistently through the prefilter");
		check(main.read(key, key2+1)==0, "Prefilter admitted a primary-only collision");
	}

	/** Verifies deterministic secondary hashes for primitive long keys. */
	private static void testPrimitiveKeys(){
		final KCountArray7MTA array=new KCountArray7MTA(10007, 8, 3, null, 0);
		final long[] keys={0, 1, 3, 0x3FFFFFFFL, 0x5555555555555555L, Long.MAX_VALUE};
		for(long key : keys){
			final long key2=KCountArray.secondaryHash(key);
			check(key2==KCountArray.secondaryHash(key), "Primitive secondary hash was nondeterministic");
			array.incrementDual(key, key2);
		}
		for(long key : keys){
			check(array.read(key, KCountArray.secondaryHash(key))>0, "Primitive dual-hash key was absent: "+key);
		}
	}

	/** Verifies that ukmer secondary hashes are canonical across reverse complements. */
	private static void testReverseComplements(){
		final Random randy=new Random(1);
		for(int k : new int[] {32, 64, 95, 128, 300}){
			for(int sample=0; sample<100; sample++){
				final byte[] bases=new byte[k];
				for(int i=0; i<k; i++){bases[i]=AminoAcid.numberToBase[randy.nextInt(4)];}
				final Kmer a=makeKmer(bases);
				final Kmer b=makeKmer(AminoAcid.reverseComplementBases(bases));
				check(a.xor()==b.xor(), "Primary hash was not canonical at k="+k);
				check(a.xor2()==b.xor2(), "Secondary hash was not canonical at k="+k);
			}
		}
	}

	/** Samples dual-hash distribution and finite-table false-positive behavior. */
	private static void testDistributionAndFalsePositives(){
		final KCountArray7MTA array=new KCountArray7MTA(20011, 1, 3, null, 0);
		final Random randy=new Random(2);
		long ones=0, bits=0;
		for(int i=0; i<2000; i++){
			final long key=randy.nextLong()&Long.MAX_VALUE;
			final long key2=randy.nextLong()&Long.MAX_VALUE;
			array.incrementDual(key, key2);
			for(int lane=0; lane<3; lane++){
				ones+=Long.bitCount(array.hashForLane(key, key2, lane));
				bits+=64;
			}
		}
		final double oneFraction=ones/(double)bits;
		check(oneFraction>0.48 && oneFraction<0.52, "Lane hashes were imbalanced: "+oneFraction);

		int falsePositives=0;
		final int queries=20000;
		for(int i=0; i<queries; i++){
			final long key=randy.nextLong()&Long.MAX_VALUE;
			final long key2=randy.nextLong()&Long.MAX_VALUE;
			if(array.read(key, key2)>0){falsePositives++;}
		}
		final double rate=falsePositives/(double)queries;
		check(rate<0.05, "Dual-hash false-positive rate was excessive: "+rate);
	}

	/** Builds a ukmer.Kmer from a concrete base array. */
	private static Kmer makeKmer(byte[] bases){
		final Kmer kmer=new Kmer(bases.length);
		for(byte b : bases){kmer.addRight(b);}
		return kmer;
	}

	/** Counts and enforces a required fixture invariant. */
	private static void check(boolean condition, String message){
		checks++;
		if(!condition){throw new RuntimeException(message);}
	}

	private static long checks=0;
}
