package prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.AminoAcid;
import map.LongHashSet;
import structures.IntList;

/** Correctness gates for the shared conserved-RNA 17-mer front end. */
public class ConservedRnaSeedIndexTest {

	public static void main(String[] args){
		testCollisionsAllBitsAndStrands();
		testCapacityGuard();
		System.out.println("PASS ConservedRnaSeedIndexTest");
	}

	private static void testCollisionsAllBitsAndStrands(){
		final LongHashSet[] sets=new LongHashSet[32];
		final long a=encode(repeat('A', 17));
		final long c=encode(repeat('C', 17));
		final long g=encode(repeat('G', 17));
		for(int i=0; i<sets.length; i++){
			sets[i]=new LongHashSet(4);
			sets[i].add(a); // One seed occupies all 32 bits, producing mask -1.
		}
		sets[1].add(c);
		sets[2].add(g);
		final ConservedRnaSeedIndex index=new ConservedRnaSeedIndex(sets);
		if(index.slotCount()!=32){throw new AssertionError("Expected 32 slots, got "+index.slotCount());}

		final byte[] forward=bytes("NN"+repeat('A', 17)+"N"+repeat('C', 17)+"N");
		final byte[] reverse=AminoAcid.reverseComplementBases(forward);
		assertMatchesLegacy("forward", index.scan(forward), forward, sets);
		assertMatchesLegacy("reverse", index.scan(reverse), reverse, sets);
		if(!Arrays.equals(index.scan(forward).hits(0), new int[]{10})){
			throw new AssertionError("Unexpected all-family collision positions: "+Arrays.toString(index.scan(forward).hits(0)));
		}
		if(!Arrays.equals(index.scan(forward).hits(1), new int[]{10, 28})){
			throw new AssertionError("Unexpected multi-seed slot positions: "+Arrays.toString(index.scan(forward).hits(1)));
		}
	}

	private static void testCapacityGuard(){
		final LongHashSet[] sets=new LongHashSet[33];
		for(int i=0; i<sets.length; i++){
			sets[i]=new LongHashSet(2);
			sets[i].add(i);
		}
		boolean threw=false;
		try{new ConservedRnaSeedIndex(sets);}
		catch(IllegalArgumentException e){threw=e.getMessage().contains("capacity is 32");}
		if(!threw){throw new AssertionError("Expected loud failure above 32 family/type slots");}
	}

	private static void assertMatchesLegacy(String label, ConservedRnaSeedIndex.ScanResult result,
			byte[] bases, LongHashSet[] sets){
		for(int slot=0; slot<sets.length; slot++){
			final int[] expected=legacyScan(bases, sets[slot]);
			final int[] observed=result.hits(slot);
			if(!Arrays.equals(expected, observed)){
				throw new AssertionError(label+" slot "+slot+" mismatch: expected="
					+Arrays.toString(expected)+", observed="+Arrays.toString(observed));
			}
		}
	}

	/** Reference implementation copied from the pre-unification per-family scanners. */
	private static int[] legacyScan(byte[] bases, LongHashSet set){
		final int k=ConservedRnaSeedIndex.K;
		final long mask=~((-1L)<<(2*k));
		final IntList hits=new IntList();
		long kmer=0;
		int len=0;
		for(int i=0; i<bases.length; i++){
			final int x=AminoAcid.baseToNumber[bases[i]];
			if(x<0){len=0; kmer=0;}
			else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				if(len>=k && set.contains(kmer)){hits.add(i-k/2);}
			}
		}
		return hits.toArray();
	}

	private static long encode(String s){
		long kmer=0;
		for(byte b : bytes(s)){
			final int x=AminoAcid.baseToNumber[b];
			if(x<0){throw new IllegalArgumentException("Non-ACGT test seed: "+s);}
			kmer=(kmer<<2)|x;
		}
		return kmer;
	}

	private static String repeat(char c, int n){
		final char[] array=new char[n];
		Arrays.fill(array, c);
		return new String(array);
	}

	private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
}
