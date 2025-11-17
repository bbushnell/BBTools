package fun;

/**
 * Calculates probabilistic k-mer intersection and cardinality for genomic
 * sequence analysis. Computes statistical probabilities of k-mer overlaps
 * between two DNA/RNA sequences of different lengths using combinatorial
 * probability theory.
 *
 * @author Brian Bushnell
 */
public class ProbShared {

	/**
	 * Program entry point that demonstrates k-mer cardinality and intersection
	 * probability calculations.
	 * @param args Command-line arguments: [k-mer_length] [sequence1_length] [sequence2_length]
	 */
	public static void main(String args[]){
		int k=Integer.parseInt(args[0]);
		int len1=Integer.parseInt(args[1]);
		int len2=Integer.parseInt(args[2]);

		System.out.println("Cardinality 1: "+cardinality(k, len1));
		System.out.println("Cardinality 2: "+cardinality(k, len2));
		System.out.println("Probability:   "+probIntersect(k, len1, len2));
		
	}
	
	/**
	 * Estimates the expected number of unique k-mers in a sequence of given length.
	 * Uses probability theory to account for potential k-mer collisions in the
	 * 4^k k-mer space. Iterates through each k-mer position, calculating the
	 * probability that it represents a new unique k-mer.
	 *
	 * @param k Length of k-mers to analyze
	 * @param seqLength Total length of the sequence
	 * @return Estimated count of unique k-mers in the sequence
	 */
	static int cardinality(int k, int seqLength){
		double space=Math.pow(4, k);
		int kmers=seqLength-k+1;
		double unique=0;
		for(int i=0; i<kmers; i++){
			double prob=(space-unique)/space;
			unique+=prob;
		}
		return (int)Math.round(unique);
	}

	/**
	 * Calculates the probability that two sequences share at least one common k-mer.
	 * Uses the cardinality estimates of both sequences and iteratively computes
	 * the cumulative probability of finding shared k-mers in the 4^k k-mer space.
	 *
	 * @param k Length of k-mers to analyze
	 * @param len1 Length of the first sequence
	 * @param len2 Length of the second sequence
	 * @return Probability (0.0-1.0) that the sequences share at least one k-mer
	 */
	static double probIntersect(int k, int len1, int len2){
		int card1=cardinality(k, len1);
		int card2=cardinality(k, len2);
		double space=Math.pow(4, k);
		double cumulativeProbUnshared=1;
		for(int i=0; i<card1; i++){
			double probShared=card2/space;
			double probUnshared=1-probShared;
			space-=probUnshared;
			cumulativeProbUnshared*=probUnshared;
		}
		return 1-cumulativeProbUnshared;
	}
	
}
