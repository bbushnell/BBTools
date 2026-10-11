package var2;

import shared.Tools;

/**
 * Binomial-derived balance score used by variant strand/read-bias calculations
 * and feature vectors. The lookup table contains capped two-sided fair-coin tails;
 * eventProb adjusts the counts and sometimes the result, so its return is not the
 * unadjusted probability of the original observed counts.
 * Tables are initialized once, in dependency order, through a maximum total of 100.
 *
 * @author Brian Bushnell
 * @author Isla
 * @contributor Sucrose
 * @date June 25, 2025
 */
public class VarProb{

	/*--------------------------------------------------------------*/
	/*----------------    Probability Calculation   ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Calculates a symmetric balance score, with larger values indicating less bias.
	 * Reduces the count imbalance using a 0.75 allowance and 0.95 multiplier, then
	 * rescales an adjusted total above PROBLEN to fit the table. The scaled minority
	 * count is truncated; the total is rounded, and the minority lookup index gets
	 * one extra count capped at half the rounded total.
	 * If the table returns one for counts differing by more than one, returns
	 * {@code 0.998+0.002*min(a,b)/max(a,b)} instead. These adjustments make this
	 * a score rather than a literal binomial test on the original counts.
	 * Inputs are nonnegative counts; their sum is computed as int before conversion
	 * to double, so totals beyond the int range are not supported by this arithmetic.
	 *
	 * @param a Count of the first group, such as plus-strand or read-1 observations
	 * @param b Count of the second group, such as minus-strand or read-2 observations
	 * @return Adjusted balance score; zero/zero and equal counts return one
	 */
	public static double eventProb(int a, int b){
		double allowedBias=0.75;// Allowance used to reduce the observed imbalance
		double slopMult=0.95;// Bias adjustment multiplier

		double n=a+b;// Total observations
		double k=Tools.min(a, b);// Minority count

		// Apply bias tolerance adjustments
		double slop=n*(allowedBias*0.5);
		double dif=n-k*2;// Difference from expected 50/50
		dif=dif-(Tools.min(slop, dif)*slopMult);
		n=k*2+dif;// Adjusted total

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
	 * Precomputes floating-point n! for indices zero through len-1.
	 * No overflow check is made; the current caller supplies PROBLEN+1.
	 * @param len Positive number of entries
	 * @return New array of factorial values
	 */
	private static double[] makeFactorialArray(int len){
		double[] x=new double[len];
		x[0]=1;// 0! = 1 by definition
		for(int i=1; i<len; i++){
			x[i]=x[i-1]*i;// n! = (n-1)! * n
		}
		return x;
	}

	/**
	 * Creates binomial coefficient matrix for "n choose k" calculations.
	 * Precomputes combinations for efficient probability calculation.
	 * Only stores {@code 0<=k<=floor(n/2)}, using factorial symmetry.
	 * Requires the factorial table to have been initialized with at least len entries.
	 *
	 * @param len Number of rows; the largest n is len-1
	 * @return New ragged matrix of floating-point binomial coefficients
	 */
	private static double[][] makeBinomialMatrix(int len){
		double[][] matrix=new double[len][];
		for(int n=0; n<len; n++){
			final int kmax=n/2;// Only need half due to symmetry
			final double nf=factorial[n];
			matrix[n]=new double[kmax+1];
			for(int k=0; k<=kmax; k++){
				final double kf=factorial[k];
				final double nmkf=factorial[n-k];
				double combinations=nf/kf;// n!/k!
				combinations=combinations/nmkf;// (n!/k!)/(n-k)!
				matrix[n][k]=combinations;
			}
		}
		return matrix;
	}

	/**
	 * Unused iterative combinatorial helper, not called by eventProb or table setup.
	 * Alternates multiplication and division around a magnitude threshold.
	 * The loop bounds need review before use: counters advance before factors are
	 * consumed, including n+1, k+1 and n-k+1. This is not a validated n-choose-k routine.
	 *
	 * @param n Intended total count
	 * @param k Intended selected count
	 * @return Product accumulated by the legacy loop
	 */
	private static double bigBinomial(int n, int k){
		//TODO: Probable bug - preincrement with <= limits includes extra numerator/denominator factors (VAR2-009); no active caller.
		double combinations=1;

		// Iteratively build combination value while managing precision
		for(int a=-1, b=-1, c=-1; a<=n || b<=k || c<=n-k; ){
			double mult=1;
			if(combinations<1000000000 && a<=n){// Multiply phase
				a++;
				mult=Tools.max(1, a);
			}else if(b<=k){// Divide by k! phase
				b++;
				mult=1.0/Tools.max(1, b);
			}else if(c<=n-k){// Divide by (n-k)! phase
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
	 * Builds capped two-sided fair-coin tails for the stored minority-count range.
	 * Each row accumulates twice the binomial mass through k, capping the sum at one:
	 * {@code prob[n][k]=min(1, 2*P(X<=k))}, with X binomial(n, 0.5).
	 * Requires the binomial table to be initialized with at least len rows.
	 *
	 * @param len Number of rows; the largest n is len-1
	 * @return New ragged matrix, storing k through floor(n/2)
	 */
	private static double[][] makeProbMatrix(int len){
		double[][] matrix=new double[len][];
		double mult=2;// Row multiplier 2^(1-n), initially 2 for n=0

		for(int n=0; n<len; n++){
			final int kmax=n/2;
			final double[] array=matrix[n]=new double[kmax+1];

			// Calculate probability for each k value
			for(int k=0; k<=kmax; k++){
				final double combinations=binomial[n][k];
				array[k]=combinations*mult;// Scale by 2^(-n+1)
			}

			// Convert to cumulative probabilities
			for(int k=0; k<=kmax; k++){
				array[k]=Tools.min(1, (k==0 ? 0 : array[k-1])+array[k]);
			}

			mult*=0.5;// Multiply by 2^(-1) for next n
		}
		return matrix;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Data           ----------------*/
	/*--------------------------------------------------------------*/

	/** Largest total represented by a table row. */
	static final int PROBLEN=100;

	/** Factorials for indices 0..PROBLEN; initialized before binomial. */
	private static final double[] factorial=makeFactorialArray(PROBLEN+1);

	/** Half-rows of binomial coefficients; initialized before prob. */
	private static final double[][] binomial=makeBinomialMatrix(PROBLEN+1);

	/** Capped two-sided tails; callers must not mutate the shared rows. */
	static final double[][] prob=makeProbMatrix(PROBLEN+1);
}
