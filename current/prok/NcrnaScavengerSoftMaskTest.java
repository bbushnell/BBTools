package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import map.LongHashSet;
import shared.Shared;

/** Regression checks for case-insensitive ncRNA alignment on soft-masked input. */
public class NcrnaScavengerSoftMaskTest {

	public static void main(String[] args){
		final byte[] consensus=bytes("ACGTTGCAAGTCGATCGTACGATGC");
		final byte[] upper=bytes("TTTTTTTTTTACGTTGCAAGTCGATCGTACGATGCGGGGGGGGGG");
		final byte[] masked=upper.clone();
		for(int i=10; i<10+consensus.length; i++){
			if(masked[i]>='A' && masked[i]<='Z'){masked[i]=(byte)(masked[i]+('a'-'A'));}
		}
		final byte[] upperBefore=upper.clone();
		final byte[] maskedBefore=masked.clone();

		final NcrnaScavenger scavenger=new NcrnaScavenger(
			new byte[][]{consensus}, null, null, null, 17, 10, 10);
		final Orf upperOrf=new Orf("upper", 8, 36, Shared.PLUS, 0, upper, false, ProkObject.RNA);
		final Orf maskedOrf=new Orf("masked", 8, 36, Shared.PLUS, 0, masked, false, ProkObject.RNA);
		scavenger.trimToAlignmentExtent(upperOrf, upper, 0, 8, 36);
		scavenger.trimToAlignmentExtent(maskedOrf, masked, 0, 8, 36);

		if(upperOrf.start!=maskedOrf.start || upperOrf.stop!=maskedOrf.stop){
			throw new AssertionError("Soft masking changed trimmed coordinates: upper="
				+upperOrf.start+"-"+upperOrf.stop+", masked="+maskedOrf.start+"-"+maskedOrf.stop);
		}
		if(!Arrays.equals(upper, upperBefore) || !Arrays.equals(masked, maskedBefore)){
			throw new AssertionError("ncRNA refinement modified the shared genome bases");
		}
		testFullScavengeCaseEquivalence(consensus, upper, masked);
		System.out.println("PASS NcrnaScavengerSoftMaskTest "+upperOrf.start+"-"+upperOrf.stop);
	}

	/** Exercises the public seed -> candidate -> alignWindow path, including the
	 * primary copyRegionUpper call site. */
	private static void testFullScavengeCaseEquivalence(byte[] consensus, byte[] upper, byte[] masked){
		final LongHashSet kmers=kmers(consensus, 17);
		final NcrnaScavenger scavenger=new NcrnaScavenger(new byte[][]{consensus}, null, null, kmers, 17, 20, 20,
			7, 10, false, 0f, 0f, 0f, 1, 0f, 1f, 0.90f, 0.80f);
		final byte[] upperBefore=upper.clone(), maskedBefore=masked.clone();
		final ArrayList<Orf> upperCalls=scavenger.scavenge("upper", upper, Shared.PLUS, new ArrayList<int[]>());
		final ArrayList<Orf> maskedCalls=scavenger.scavenge("masked", masked, Shared.PLUS, new ArrayList<int[]>());
		if(upperCalls.size()!=1 || maskedCalls.size()!=1){
			throw new AssertionError("Expected one call from each case variant, got upper="+upperCalls.size()+", masked="+maskedCalls.size());
		}
		final Orf a=upperCalls.get(0), b=maskedCalls.get(0);
		if(a.start!=b.start || a.stop!=b.stop || a.strand!=b.strand){
			throw new AssertionError("Soft masking changed full-scavenge call: upper="+a.start+"-"+a.stop+", masked="+b.start+"-"+b.stop);
		}
		if(!Arrays.equals(upper, upperBefore) || !Arrays.equals(masked, maskedBefore)){
			throw new AssertionError("Full scavenging modified shared genome bases");
		}
	}

	private static LongHashSet kmers(byte[] bases, int k){
		final LongHashSet set=new LongHashSet(bases.length);
		final long mask=~((-1L)<<(2*k));
		long kmer=0; int len=0;
		for(byte base : bases){
			final int x=AminoAcid.baseToNumber[base];
			if(x<0){len=0; kmer=0;}
			else{
				kmer=((kmer<<2)|x)&mask; len++;
				if(len>=k){set.add(kmer);}
			}
		}
		return set;
	}

	private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
}
