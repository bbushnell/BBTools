package test.jgi;

import java.util.Arrays;

import jgi.KmerSample;

/**
 * Tiny production-linked fixture for KmerSample's public bit-vector helpers.
 * The FASTA input keeps K=1, K=2 and K=3 finite while crossing a sequence line.
 * @author Jean
 */
public final class KmerSampleTest{

	/** Runs the literal FASTA fixture and prints every case result. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1 || !args[0].regionMatches(true, 0, "in=", 0, 3)){
			throw new IllegalArgumentException("Expected one in=<FASTA> argument");
		}
		final String fastaName=args[0].substring(3);

		int failures=0;
		int cases=0;
		for(final Case c : CASES){
			cases++;
			failures+=check(c, fastaName);
		}

		System.out.println("KmerSampleTest cases="+cases+" failures="+failures);
		if(failures!=0){
			throw new AssertionError("KmerSample failures="+failures);
		}
	}

	/** Executes one K value and verifies both public methods against its literal bitset. */
	private static int check(final Case c, final String fastaName){
		int failures=0;
		try{
			final int[] observed=KmerSample.makeKmerSet(c.k, fastaName);
			if(Arrays.equals(c.expected, observed)){
				System.out.println("PASS makeKmerSet K="+c.k+" array="+Arrays.toString(observed));
			}else{
				failures++;
				System.err.println("FAIL makeKmerSet K="+c.k+": expected="+
						Arrays.toString(c.expected)+", observed="+Arrays.toString(observed));
			}
		}catch(final RuntimeException t){
			failures++;
			System.err.println("FAIL makeKmerSet K="+c.k+": "+t);
		}

		for(long kmer=0; kmer<c.domain; kmer++){
			final boolean expected=expectedContains(kmer, c.presentKmers);
			try{
				final boolean observed=KmerSample.containsKmer(kmer, c.expected);
				if(observed!=expected){
					failures++;
					System.err.println("FAIL containsKmer K="+c.k+" kmer="+kmer+
							": expected="+expected+", observed="+observed);
				}
			}catch(final RuntimeException t){
				failures++;
				System.err.println("FAIL containsKmer K="+c.k+" kmer="+kmer+
						": expected="+expected+", exception="+t);
			}
		}
		if(failures==0){
			System.out.println("PASS containsKmer K="+c.k+" domain="+c.domain);
		}
		return failures;
	}

	/** Checks a literal list of kmers that should have set bits. */
	private static boolean expectedContains(final long kmer, final int[] presentKmers){
		for(final int present : presentKmers){
			if(present==kmer){return true;}
		}
		return false;
	}

	/** One K value, all possible encoded kmers for that K, and its expected bitset. */
	private static final class Case{

		Case(final int k_, final long domain_, final int[] expected_, final int[] presentKmers_){
			k=k_;
			domain=domain_;
			expected=expected_;
			presentKmers=presentKmers_;
		}

		final int k;
		final long domain;
		final int[] expected;
		final int[] presentKmers;
	}

	private static final Case[] CASES={
			new Case(1, 4, new int[]{0xF}, new int[]{0, 1, 2, 3}),
			new Case(2, 16, new int[]{0x8842}, new int[]{1, 6, 11, 15}),
			new Case(3, 64, new int[]{0x08000040, 0x80000000}, new int[]{6, 27, 63})
	};
}
