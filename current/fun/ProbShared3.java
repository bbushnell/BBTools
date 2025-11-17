package fun;

import java.util.HashSet;
import java.util.Random;

/**
 * Probabilistic simulation tool for estimating k-mer occurrence probability
 * between random sequences using Monte Carlo methods.
 * Calculates the likelihood of finding shared k-mers between two sequences
 * of specified lengths through statistical sampling.
 *
 * @author Brian Bushnell
 */
public class ProbShared3 {

	/**
	 * Program entry point for k-mer probability simulation.
	 * Runs Monte Carlo simulation with specified parameters and outputs
	 * probability estimate.
	 *
	 * @param args Command-line arguments: [k-mer_size] [sequence1_length]
	 * [sequence2_length] [simulation_rounds]
	 */
	public static void main(String args[]){
		int k=Integer.parseInt(args[0]);
		int len1=Integer.parseInt(args[1]);
		int len2=Integer.parseInt(args[2]);
		int rounds=Integer.parseInt(args[3]);
		
		System.out.println("Probability:   "+simulate(k, len1, len2, rounds));
	}
	
	/**
	 * Performs Monte Carlo simulation to estimate probability of shared k-mers.
	 * Runs multiple simulation rounds and calculates success rate as
	 * proportion of rounds containing at least one matching k-mer.
	 *
	 * @param k K-mer size (bases per k-mer)
	 * @param len1 Length of first sequence
	 * @param len2 Length of second sequence
	 * @param rounds Number of simulation rounds to perform
	 * @return Probability estimate as double between 0.0 and 1.0
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
	 * Simulates one pair of sequences to test for k-mer intersection.
	 * Generates random k-mers for second sequence in set, then tests if
	 * any k-mer from first sequence matches.
	 *
	 * @param k K-mer size in bases
	 * @param len1 Length of first sequence (tested against set)
	 * @param len2 Length of second sequence (used to fill set)
	 * @param set HashSet for storing k-mers from second sequence
	 * @return 1 if intersection found, 0 otherwise
	 */
	static int simulateOnePair(int k, int len1, int len2, HashSet<Long> set){
		fillRandomSet(k, len2, set);
		final long space=(long)Math.pow(4, k);
		final int kmers=len1-k+1;
		for(int i=0; i<kmers; i++){
			long kmer=(randy.nextLong()&Long.MAX_VALUE)%space;
			if(set.contains(kmer)){return 1;}
		}
		return 0;
	}
	
	/**
	 * Fills HashSet with random k-mers representing a sequence of given length.
	 * Generates (len - k + 1) random k-mers within the k-mer space and
	 * adds them to the provided set.
	 *
	 * @param k K-mer size in bases
	 * @param len Sequence length to simulate
	 * @param set HashSet to populate with random k-mers
	 */
	static void fillRandomSet(int k, int len, HashSet<Long> set){
		set.clear();
		final long space=(long)Math.pow(4, k);
		final int kmers=len-k+1;
		for(int i=0; i<kmers; i++){
			set.add((randy.nextLong()&Long.MAX_VALUE)%space);
		}
	}
	
	/** Random number generator for k-mer simulation */
	static final Random randy=new Random();
	
}
