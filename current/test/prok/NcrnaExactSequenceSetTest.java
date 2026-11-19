package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import prok.NcrnaExactSequenceSet;

/** Exact-gene lookup checks used by the saved-window overlap fixture.
 * Includes two distinct DNA sequences with a deliberate polynomial-key collision.
 * @author Raiden
 */
public final class NcrnaExactSequenceSetTest{

	public static void run(){
		final NcrnaExactSequenceSet set=new NcrnaExactSequenceSet();
		check(set.add(dna("acgtu")), "The first complete member is retained");
		check(!set.add(dna("ACGTT")) && set.size()==1, "Case and U/T normalization collapse presence without adding a vote");
		check(set.contains(dna("AcGtT")), "Lookup applies the same DNA normalization as insertion");
		check(!set.contains(dna("AACGT")), "Reverse-complement identity is not sense-oriented gene identity");
		check(!set.contains(dna("ACGTA")), "A full-length single-base mismatch is not an exact member");
		check(!set.contains(dna("ACGT")) && !set.contains(dna("TACGTT")), "Neither a prefix nor a containing sequence is exact identity");

		final byte[] a=new byte[61], b=new byte[61];
		Arrays.fill(a, (byte)'A');
		Arrays.fill(b, (byte)'A');
		// Difference polynomial (x^4-1)(x^8-1)(x^16-1)(x^32-1).
		// At x=31 it is divisible by 2^34, so both 32-bit lookup keys collide.
		// ASCII A/C and normalized A/C differ by a fixed factor; both collide.
		for(int mask=0; mask<16; mask++){
			final int position=60-4*mask;
			if((Integer.bitCount(mask)&1)==0){a[position]='C';}
			else{b[position]='C';}
		}
		check(!Arrays.equals(a, b) && Arrays.hashCode(a)==Arrays.hashCode(b), "The literal collision fixture must contain distinct equal-length DNA arrays with equal polynomial keys");
		check(set.add(a) && !set.contains(b), "A colliding lookup key cannot establish sequence identity");
		check(set.add(b) && set.contains(a) && set.contains(b) && set.size()==3, "Both colliding sequences remain independently reachable");
		reject(set, null, false);
		reject(set, new byte[0], true);
		reject(set, dna("ACNT"), false);
		reject(set, dna("AC-T"), true);
	}

	private static void reject(final NcrnaExactSequenceSet set, final byte[] bases, final boolean insert){
		boolean rejected=false;
		try{
			if(insert){set.add(bases);}
			else{set.contains(bases);}
		}catch(final IllegalArgumentException expected){rejected=true;}
		check(rejected, "Absent, partial or undefined input cannot silently become a negative exact-gene result");
	}

	private static byte[] dna(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean condition, final String why){if(!condition){throw new AssertionError(why);}}
}
