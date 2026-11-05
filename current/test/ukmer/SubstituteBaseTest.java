package test.ukmer;

import ukmer.Kmer;

import dna.AminoAcid;

/**
 * Checks every position and all four replacement bases on deterministic sequences.
 * Rebuilds each expected Kmer through addRight from independently edited bases,
 * then compares restoration with the untouched original, including lazy caches.
 * Covers packed lengths32,33,63,64,95,127 and symmetric lengths62,93; it does not
 * exhaust every possible input sequence or layout setting. Same-base edits are included.
 * This is the committed Item3c/Fix3 regression invoked by dev/substitutebasetest.sh;
 * the earlier scratch-only description was incorrect, as Noelle noted.
 * Use a fresh JVM: counters are cumulative and normal completion sets PACKED=true
 * rather than restoring its incoming value. Failure throws after printing details.
 */
public class SubstituteBaseTest{

	/** Cumulative edit cases and detected mismatches; not reset by main. */
	static int total=0, failures=0;
	/** Details of all mismatches retained until the final report. */
	static StringBuilder failLog=new StringBuilder();

	/** Runs the fixed packed/symmetric grids and reports all collected mismatches.
	 * @param args Unused command-line arguments
	 * @throws RuntimeException If any explicit comparison failed */
	public static void main(String[] args){
		String[] seqs={
			"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
			"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",//poly-T
			"GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA",//RC-pair-prone (palindromic 4-mers, canonical-flip-prone)
			"ATCGGGCATTAGCGATCCATGGCTAGCTAGGCATCGATCGGATCCGATCGTAGCTAGCATGCTAGCTGATCGATCGGCATCGATCGTAGCATGCATCGGATCGCATGCTAGCATCGATGCATGCATG",
		};

		System.out.println("=== PACKED=true ===");
		Kmer.PACKED=true;
		int[] ksPacked={32, 33, 63, 64, 95, 127};
		for(int k : ksPacked){
			for(String fullSeq : seqs){
				runGrid(k, fullSeq.substring(0, k));
			}
		}

		System.out.println("=== PACKED=false (legacy symmetric, representative subset) ===");
		Kmer.PACKED=false;
		int[] ksLegacy={62, 93};
		for(int kbig : ksLegacy){
			for(String fullSeq : seqs){
				runGrid(kbig, fullSeq.substring(0, kbig));
			}
		}
		Kmer.PACKED=true;//restore default

		System.out.println("Total checks: "+total);
		System.out.println("Failures: "+failures);
		if(failures>0){
			System.out.println("--- failure log ---");
			System.out.println(failLog.toString());
			throw new RuntimeException("SubstituteBaseTest FAILED: "+failures+"/"+total+" checks failed -- see log above.");
		}
		System.out.println("ALL PASS");
	}

	/** Checks each edit and restoration against separately built sequence keys.
	 * Forces lazy key/hash/palindrome regeneration before comparing cached state.
	 * @param kbig Requested length, equal to the effective length for these cases
	 * @param seq Defined DNA sequence containing exactly kbig bases */
	private static void runGrid(int kbig, String seq){
		byte[] bases=seq.getBytes();

		Kmer baseline=buildKmer(kbig, bases);
		long[] baseArray1=baseline.array1().clone();
		long[] baseArray2=baseline.array2().clone();
		long[] baseKey=baseline.key().clone();//FULL key array, not just key[0]
		long baseXor=baseline.xor();
		int baseLen=baseline.len();
		boolean baseCorePalindrome=baseline.corePalindrome();
		String baseStr=baseline.toString();
		boolean baseVerify=baseline.verify(true);
		if(!baseVerify){
			failures++;
			failLog.append("BASELINE-VERIFY-FAIL kbig=").append(kbig).append(" seq=").append(seq).append("\n");
		}

		for(int pos=0; pos<kbig; pos++){
			long origBase=AminoAcid.baseToNumber[bases[pos]];
			for(long alt=0; alt<4; alt++){
				total++;
				Kmer scratch=buildKmer(kbig, bases);

				long returned=scratch.substituteBase(pos, alt);
				if(returned!=origBase){
					failures++;
					failLog.append("RETURN-MISMATCH kbig=").append(kbig).append(" pos=").append(pos)
						.append(" alt=").append(alt).append(" expected=").append(origBase)
						.append(" got=").append(returned).append("\n");
				}

				byte[] editedBases=bases.clone();
				editedBases[pos]=AminoAcid.numberToBase[(int)alt];
				Kmer expected=buildKmer(kbig, editedBases);//oracle: independent rebuild, never derived from scratch

				//Force lazy regen on BOTH sides before comparing key/xor/corePalindrome (Amber's point).
				long[] scratchKey=scratch.key(); long scratchXor=scratch.xor(); boolean scratchCP=scratch.corePalindrome();
				long[] expectedKey=expected.key(); long expectedXor=expected.xor(); boolean expectedCP=expected.corePalindrome();

				boolean arraysMatch=java.util.Arrays.equals(scratch.array1(), expected.array1())
					&& java.util.Arrays.equals(scratch.array2(), expected.array2());
				boolean keyMatch=java.util.Arrays.equals(scratchKey, expectedKey);//full array, not just [0]
				boolean xorMatch=(scratchXor==expectedXor);
				boolean lenMatch=(scratch.len()==expected.len());
				boolean cpMatch=(scratchCP==expectedCP);
				boolean strMatch=scratch.toString().equals(expected.toString());
				boolean verifyOk=scratch.verify(true);

				if(!(arraysMatch && keyMatch && xorMatch && lenMatch && cpMatch && strMatch && verifyOk)){
					failures++;
					failLog.append("EDIT-MISMATCH kbig=").append(kbig).append(" pos=").append(pos)
						.append(" alt=").append(alt)
						.append(" arrays=").append(arraysMatch)
						.append(" key=").append(keyMatch)
						.append(" xor=").append(xorMatch)
						.append(" len=").append(lenMatch)
						.append(" corePalindrome=").append(cpMatch)
						.append(" str=").append(strMatch)
						.append(" verify=").append(verifyOk)
						.append(" scratchStr=").append(scratch.toString())
						.append(" expectedStr=").append(expected.toString())
						.append("\n");
				}

				final long restored=scratch.substituteBase(pos, returned);//restore
				if(restored!=alt){
					failures++;
					failLog.append("RESTORE-RETURN-MISMATCH kbig=").append(kbig).append(" pos=").append(pos)
						.append(" alt=").append(alt).append(" expected=").append(alt)
						.append(" got=").append(restored).append("\n");
				}
				long[] restoredKey=scratch.key(); long restoredXor=scratch.xor(); boolean restoredCP=scratch.corePalindrome();

				boolean restoreArraysMatch=java.util.Arrays.equals(scratch.array1(), baseArray1)
					&& java.util.Arrays.equals(scratch.array2(), baseArray2);
				boolean restoreKeyMatch=java.util.Arrays.equals(restoredKey, baseKey);//full array
				boolean restoreXorMatch=(restoredXor==baseXor);
				boolean restoreLenMatch=(scratch.len()==baseLen);
				boolean restoreCpMatch=(restoredCP==baseCorePalindrome);
				boolean restoreStrMatch=scratch.toString().equals(baseStr);
				boolean restoreVerifyOk=scratch.verify(true);

				if(!(restoreArraysMatch && restoreKeyMatch && restoreXorMatch && restoreLenMatch
						&& restoreCpMatch && restoreStrMatch && restoreVerifyOk)){
					failures++;
					failLog.append("RESTORE-MISMATCH kbig=").append(kbig).append(" pos=").append(pos)
						.append(" alt=").append(alt)
						.append(" arrays=").append(restoreArraysMatch)
						.append(" key=").append(restoreKeyMatch)
						.append(" xor=").append(restoreXorMatch)
						.append(" len=").append(restoreLenMatch)
						.append(" corePalindrome=").append(restoreCpMatch)
						.append(" str=").append(restoreStrMatch)
						.append(" verify=").append(restoreVerifyOk)
						.append("\n");
				}
			}
		}
	}

	/** Rebuilds one key from the first kbig bases using the current packing settings.
	 * @param kbig Requested length, equal to the effective length for these cases
	 * @param bases Defined DNA bases with at least kbig entries
	 * @return Newly allocated, independently populated key */
	private static Kmer buildKmer(int kbig, byte[] bases){
		Kmer kmer=new Kmer(kbig);
		kmer.clear();
		for(int i=0; i<kbig; i++){
			kmer.addRight(bases[i]);
		}
		return kmer;
	}
}
