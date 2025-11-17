package fun;

import java.util.HashSet;
import java.util.Random;

import dna.AminoAcid;

/**
 * Probabilistic simulation tool for calculating k-mer overlap between random DNA sequences.
 * Uses Monte Carlo simulation to estimate the probability of finding shared k-mers
 * between two randomly generated sequences of specified lengths.
 * Employs efficient bit manipulation for k-mer encoding and HashSet for unique k-mer tracking.
 *
 * @author Brian Bushnell
 */
public class ProbShared2 {

	/**
	 * Program entry point for k-mer overlap probability simulation.
	 * @param args Command-line arguments: k-mer length, sequence1 length,
	 * sequence2 length, simulation rounds
	 */
	public static void main(String args[]){
		int k=Integer.parseInt(args[0]);
		int len1=Integer.parseInt(args[1]);
		int len2=Integer.parseInt(args[2]);
		int rounds=Integer.parseInt(args[3]);
		
		System.out.println("Probability:   "+simulate(k, len1, len2, rounds));
	}
	
	/**
	 * Performs Monte Carlo simulation to estimate k-mer overlap probability.
	 * Runs multiple simulation rounds and calculates the proportion of successful
	 * k-mer matches between randomly generated sequence pairs.
	 *
	 * @param k K-mer length for sequence analysis
	 * @param len1 Length of first sequence
	 * @param len2 Length of second sequence
	 * @param rounds Number of simulation rounds to perform
	 * @return Probability of k-mer overlap as a double between 0.0 and 1.0
	 */
	static double simulate(int k, int len1, int len2, int rounds){
		int successes=0;
		final HashSet<Long> set=new HashSet<Long>();
		for(int i=0; i<rounds; i++){
			successes+=simulateOnePair(k, len1, len2, set);
		}
		return successes/(double)rounds;
	}
	
	/**
	 * Simulates k-mer overlap between a single pair of random sequences.
	 * Generates first sequence, extracts k-mers into HashSet, then generates
	 * second sequence and checks for k-mer matches.
	 *
	 * @param k K-mer length for analysis
	 * @param len1 Length of first sequence
	 * @param len2 Length of second sequence
	 * @param set Reusable HashSet for k-mer storage
	 * @return 1 if sequences share at least one k-mer, 0 otherwise
	 */
	static int simulateOnePair(int k, int len1, int len2, HashSet<Long> set){
		set.clear();
		
		final int shift=2*k;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		int len=0;
		
		byte[] bases=randomSequence(len2);
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){len=0;}else{len++;}
			if(len>=k){
				set.add(kmer);
			}
		}
		
		bases=randomSequence(len1);
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){len=0;}else{len++;}
			if(len>=k){
				if(set.contains(kmer)){return 1;}
			}
		}
		return 0;
	}
	
	/**
	 * Generates a random DNA sequence of specified length.
	 * Uses uniform random distribution to select bases (A, T, G, C)
	 * from the four-letter DNA alphabet.
	 *
	 * @param len Desired sequence length in bases
	 * @return Byte array representing random DNA sequence
	 */
	static byte[] randomSequence(int len){
		byte[] array=new byte[len];
		for(int i=0; i<len; i++){
			int number=randy.nextInt(4);
			array[i]=AminoAcid.numberToBase[number];
		}
		return array;
	}
	
	/** Random number generator for sequence generation */
	static final Random randy=new Random();
	/** Lookup table for converting numeric values to DNA bases */
	static final byte[] numberToBase=AminoAcid.numberToBase;
	/** Lookup table for converting DNA bases to numeric values */
	static final byte[] baseToNumber=AminoAcid.baseToNumber;
	
}
