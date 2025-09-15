package var2;

import shared.Tools;

/**
 * Statistical probability calculations for variant calling quality assessment.
 * Provides binomial probability analysis to evaluate the likelihood of observed
 * variant patterns occurring by chance, helping distinguish real variants from
 * sequencing artifacts or random noise.
 * 
 * Uses precomputed probability matrices for efficient calculation of binomial
 * event probabilities with bias adjustments for real-world sequencing data.
 * 
 * @author Brian Bushnell
 * @author Isla Winglet
 * @date June 25, 2025
 */
public class VarProb {

	/*--------------------------------------------------------------*/
	/*----------------    Probability Calculation   ----------------*/
	/*--------------------------------------------------------------*/

	/** 
	 * Calculates adjusted probability of a binomial event being at least this lopsided.
	 * Accounts for expected sequencing bias and provides realistic probability estimates
	 * for variant calling decisions. Used to assess whether observed allele ratios
	 * are likely due to chance or represent genuine biological variation.
	 * 
	 * @param a Count of first outcome (e.g., reference alleles)
	 * @param b Count of second outcome (e.g., alternative alleles) 
	 * @return Probability that this outcome or more extreme could occur by chance
	 */
	public static double eventProb(int a, int b){
		double allowedBias=0.75;    // Expected sequencing bias tolerance
		double slopMult=0.95;       // Bias adjustment multiplier

		double n=a+b;               // Total observations
		double k=Tools.min(a, b);   // Minority count

		// Apply bias tolerance adjustments
		double slop=n*(allowedBias*0.5);
		double dif=n-k*2;           // Difference from expected 50/50
		dif=dif-(Tools.min(slop, dif)*slopMult);
		n=k*2+dif;                  // Adjusted total

		assert(k<=n*0.5) : a+", "+b+", "+n+", "+k+", "+slop+", "+dif;

		// Scale down for large values to fit precomputed matrix
		if(n>PROBLEN){
			double mult=PROBLEN/n;
			n=PROBLEN;
			k=(int)(k*mult);
		}

		int n2=(int)Math.round(n);
		int k2=Tools.min(n2/2, (int)(k+1));

		double result=prob[n2][k2];
		
		// Handle edge cases for nearly balanced outcomes
		if(result<1 || a==b || a+1==b || a==b+1){return result;}

		// Apply slope correction for moderate imbalances
		double slope=Tools.min(a, b)/(double)Tools.max(a, b);
		return (0.998+slope*0.002);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Matrix Initialization     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates factorial lookup table for combinatorial calculations.
	 * Precomputes n! for values 0 to len-1 to avoid repeated calculation.
	 * 
	 * @param len Maximum factorial to compute
	 * @return Array where factorial[n] = n!
	 */
	private static double[] makeFactorialArray(int len) {
		double[] x=new double[len];
		x[0]=1;                     // 0! = 1 by definition
		for(int i=1; i<len; i++){
			x[i]=x[i-1]*i;          // n! = (n-1)! * n
		}
		return x;
	}

	/**
	 * Creates binomial coefficient matrix for "n choose k" calculations.
	 * Precomputes combinations for efficient probability calculation.
	 * Only stores values for k <= n/2 due to symmetry: C(n,k) = C(n,n-k)
	 * 
	 * @param len Maximum n value to compute
	 * @return Matrix where binomial[n][k] = C(n,k) = n!/(k!(n-k)!)
	 */
	private static double[][] makeBinomialMatrix(int len) {
		double[][] matrix=new double[len][];
		for(int n=0; n<len; n++){
			final int kmax=n/2;         // Only need half due to symmetry
			final double nf=factorial[n];
			matrix[n]=new double[kmax+1];
			for(int k=0; k<=kmax; k++){
				final double kf=factorial[k];
				final double nmkf=factorial[n-k];
				double combinations=nf/kf;  // n!/k!
				combinations=combinations/nmkf; // (n!/k!)/(n-k)!
				matrix[n][k]=combinations;
			}
		}
		return matrix;
	}

	/** 
	 * Calculates binomial coefficients for large values using iterative approach.
	 * Used when values exceed precomputed matrix bounds.
	 * Carefully manages numerical precision to avoid overflow.
	 * 
	 * @param n Total items
	 * @param k Items to choose
	 * @return Binomial coefficient C(n,k)
	 */
	private static double bigBinomial(int n, int k){
		double combinations=1;
		
		// Iteratively build combination value while managing precision
		for(int a=-1, b=-1, c=-1; a<=n || b<=k || c<=n-k; ) {
			double mult=1;
			if(combinations<1000000000 && a<=n){    // Multiply phase
				a++;
				mult=Tools.max(1, a);
			}else if(b<=k){                         // Divide by k! phase
				b++;
				mult=1.0/Tools.max(1, b);
			}else if(c<=n-k){                       // Divide by (n-k)! phase
				c++;
				mult=1.0/Tools.max(1, c);
			}else{
				assert(false) : a+", "+b+", "+c+", "+n+", "+k+", "+(n-k)+", "+combinations;
			}
			combinations*=mult;
			assert(combinations<Double.MAX_VALUE) : a+", "+b+", "+c+", "+n+", "+k+", "+(n-k)+", "+combinations;
		}

		return combinations;
	}

	/**
	 * Creates cumulative probability matrix for binomial events.
	 * Calculates probability of observing k or fewer minority outcomes in n trials,
	 * accounting for the decreasing probability as n increases (multiplying by 0.5^n).
	 * 
	 * @param len Maximum n value to compute
	 * @return Matrix where prob[n][k] = P(X <= k) for binomial(n, 0.5)
	 */
	private static double[][] makeProbMatrix(int len) {
		double[][] matrix=new double[len][];
		double mult=2;              // Starting multiplier (2^0)
		
		for(int n=0; n<len; n++){
			final int kmax=n/2;
			final double[] array=matrix[n]=new double[kmax+1];
			
			// Calculate probability for each k value
			for(int k=0; k<=kmax; k++){
				final double combinations=binomial[n][k];
				array[k]=combinations*mult; // Scale by 2^(-n+1)
			}
			
			// Convert to cumulative probabilities
			for(int k=0; k<=kmax; k++){
				array[k]=Tools.min(1, (k==0 ? 0 : array[k-1])+array[k]);
			}
			
			mult*=0.5;              // Multiply by 2^(-1) for next n
		}
		return matrix;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Data           ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum n value for precomputed probability matrices */
	final static int PROBLEN=100;

	/** Precomputed factorial values: factorial[n] = n! */
	private static final double[] factorial=makeFactorialArray(PROBLEN+1);
	
	/** Precomputed binomial coefficients: binomial[n][k] = C(n,k) */
	private static final double[][] binomial=makeBinomialMatrix(PROBLEN+1);
	
	/** Precomputed cumulative probabilities: prob[n][k] = P(X ≤ k) for binomial(n,0.5) */
	static final double[][] prob=makeProbMatrix(PROBLEN+1);
}