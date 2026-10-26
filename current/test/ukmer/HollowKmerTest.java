package test.ukmer;

import java.util.Arrays;
import java.util.Random;

import dna.AminoAcid;
import shared.Tools;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayH1D;
import ukmer.HashBufferU;
import ukmer.Kmer;

/** Checks hollow fingerprints against independently encoded endpoint sequences.
 * @author Brian Bushnell, Fischl */
public final class HollowKmerTest {

	/** Tests word boundaries, canonical ties and buffered reconstruction without external data. */
	public static void main(final String[] args){
		Kmer.PACKED=true;
		final Random random=new Random(42019);
		for(int span : new int[]{65, 72, 96, 124, 128, 160, 192, 224}){
			for(int trial=0; trial<100; trial++){
				final byte[] bases=randomBases(span, random);
				exercise(bases);
				if(trial==0){buffered(bases);}
			}
			final byte[] palindromeArms=randomBases(span, random);
			final byte[] right=reverse(Arrays.copyOf(palindromeArms, 32));
			System.arraycopy(right, 0, palindromeArms, span-32, 32);
			exercise(palindromeArms);
		}
		boolean rejected=false;
		try{new Kmer(64, true);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "K64 has no hollow interior.");
		Kmer.PACKED=false;
		rejected=false;
		try{new Kmer(128, true);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Legacy symmetric storage must not masquerade as two full32 arms.");
		Kmer.PACKED=true;
		System.out.println("HOLLOW_KMER_TEST_PASS checks="+checks);
	}

	/** A slow base-level oracle must agree with the rolling packed implementation. */
	private static void exercise(final byte[] bases){
		final Kmer k=make(bases, true);
		final long[] oracle=oracle(bases);
		check(k.xor()==oracle[0] && k.xor2()==oracle[1], "Independent endpoint encoding differs.");
		final Kmer rc=make(reverse(bases), true);
		equalHashes(k, rc, "Reverse complements must address the same count.");
		equalHashes(k, k.clone(), "Clone must preserve hollow mode.");
		final Kmer copied=new Kmer(bases.length, true);
		copied.setFrom(k);
		equalHashes(k, copied, "setFrom must preserve hashes.");
		copied.clear();
		copied.setFrom(k.array1());
		equalHashes(k, copied, "Raw sequence reconstruction must retain hollow mode.");
		final Kmer solid=make(bases, false);
		final long solidHash=solid.xor2();
		for(int pos=32; pos<bases.length-32; pos++){
			final long old=k.substituteBase(pos, (AminoAcid.baseToNumber[bases[pos]]+1)&3);
			check(k.xor()==oracle[0] && k.xor2()==oracle[1], "Middle substitution leaked into a fingerprint.");
			k.substituteBase(pos, old);
		}
		for(int pos : new int[]{0, 31, bases.length-32, bases.length-1}){
			final long old=k.substituteBase(pos, (AminoAcid.baseToNumber[bases[pos]]+1)&3);
			check(k.xor()!=oracle[0] || k.xor2()!=oracle[1], "Anchor substitution must change the count key.");
			k.substituteBase(pos, old);
		}
		final long evicted=k.addRightNumeric(2);
		final byte[] shifted=Arrays.copyOfRange(bases, 1, bases.length+1);
		shifted[shifted.length-1]='G';
		equalHashes(k, make(shifted, true), "Rolling right must match fresh encoding.");
		k.addLeftNumeric(evicted);
		check(k.xor()==oracle[0] && k.xor2()==oracle[1], "Roll/restore changed the original key.");
		final long last=k.addLeftNumeric(1);
		final byte[] left=new byte[bases.length];
		left[0]='C';
		System.arraycopy(bases, 0, left, 1, bases.length-1);
		equalHashes(k, make(left, true), "Rolling left must match fresh encoding.");
		k.addRightNumeric(last);
		check(Arrays.equals(bases, k.toString().getBytes(java.nio.charset.StandardCharsets.US_ASCII)),
				"Hollow hashing must not change full sequence storage.");
		check(solid.xor2()==solidHash && !solid.hollow, "Hollow instances must not modify solid instances.");
	}

	/** Buffered flushing must collapse interior variants while separating changed anchors. */
	private static void buffered(final byte[] bases){
		final HashArrayH1D table=new HashArrayH1D(new int[]{101});
		final HashBufferU buffer=new HashBufferU(new AbstractKmerTableU[]{table}, 4, bases.length, false, true);
		final Kmer original=make(bases, true), variant=original.clone();
		variant.substituteBase(32, (AminoAcid.baseToNumber[bases[32]]+1)&3);
		buffer.incrementAndReturnNumCreated(original);
		buffer.incrementAndReturnNumCreated(variant);
		buffer.incrementAndReturnNumCreated(make(reverse(bases), true));
		buffer.flush();
		check(table.getValue(original)==3 && table.getValue(variant)==3, "Buffered interiors must share one count.");
		variant.substituteBase(0, (AminoAcid.baseToNumber[bases[0]]+1)&3);
		buffer.incrementAndReturnNumCreated(variant);
		buffer.flush();
		check(table.getValue(original)==3 && table.getValue(variant)==1, "Distinct anchors must remain separate.");
	}

	/** Encodes each arm base by base, independently of packed-word boundaries. */
	private static long[] oracle(final byte[] bases){
		assert(bases.length>64) : "Oracle arms must not overlap.";
		final byte[] rc=reverse(bases);
		final long f1=arm(bases, 0), f2=arm(bases, bases.length-32);
		final long r1=arm(rc, 0), r2=arm(rc, rc.length-32);
		final int first=Long.compareUnsigned(f1, r1);
		final boolean forward=first>0 || (first==0 && Long.compareUnsigned(f2, r2)>=0);
		final long a=forward ? f1 : r1, b=forward ? f2 : r2;
		return new long[]{
			(Long.rotateLeft(Tools.hash64shift(a), 25)^Tools.hash64shift(b)^Tools.hash64shift(bases.length))&Long.MAX_VALUE,
			Tools.splitMix64(Long.rotateRight(Tools.splitMix64(a), 23)^b^Tools.splitMix64(bases.length))&Long.MAX_VALUE};
	}

	/** Packs exactly one independently selected arm into a long. */
	private static long arm(final byte[] bases, final int start){
		assert(start>=0 && start+32<=bases.length) : "Both oracle arms must fit the sequence.";
		long word=0;
		for(int i=start; i<start+32; i++){word=(word<<2)|AminoAcid.baseToNumber[bases[i]];}
		return word;
	}

	/** Makes a complete rolling object; all input bases are deterministic ACGT. */
	private static Kmer make(final byte[] bases, final boolean hollow){
		final Kmer k=new Kmer(bases.length, hollow);
		for(byte b : bases){k.addRight(b);}
		return k;
	}

	/** Copies before reverse-complementing so test inputs remain immutable. */
	private static byte[] reverse(final byte[] bases){
		final byte[] copy=bases.clone();
		AminoAcid.reverseComplementBasesInPlace(copy);
		return copy;
	}

	/** Produces reproducible fixtures without a data dependency. */
	private static byte[] randomBases(final int length, final Random random){
		final byte[] bases=new byte[length], alphabet={'A', 'C', 'G', 'T'};
		for(int i=0; i<length; i++){bases[i]=alphabet[random.nextInt(4)];}
		return bases;
	}

	/** Both fingerprints, rather than only placement, define count equivalence. */
	private static void equalHashes(final Kmer a, final Kmer b, final String message){
		check(a.xor()==b.xor() && a.xor2()==b.xor2(), message);
	}

	/** Test checks remain fatal even when a caller disables Java assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
}
