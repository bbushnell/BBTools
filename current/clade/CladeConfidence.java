package clade;

import tax.TaxTree;

/**
 * Estimates the probability that a QuickClade hit is taxonomically correct
 * at a given level, based on sequence length and k-mer frequency distances.
 *
 * Model: p = floor + (1-floor) * sigmoid(a0 + a1*log2(len/2500) + a2*k3dif + a3*kdif)
 * Uses K5 parameters when k5dif is available (< 1.0), K4 otherwise.
 *
 * Trained on 559M hit observations from 7062 prokaryotic genomes shredded at
 * both fixed and variable lengths, with top-20 to top-100 hits per query.
 *
 * @author Brian Bushnell, Ady
 */
public class CladeConfidence {

	/**
	 * Probability that a hit is correct at the given taxonomic level.
	 *
	 * @param length Query sequence length in bases
	 * @param k3dif Sum of absolute 3-mer frequency differences (0-1)
	 * @param k4dif Sum of absolute 4-mer frequency differences (0-1)
	 * @param k5dif Sum of absolute 5-mer frequency differences (0-1); 1.0 = not computed
	 * @param taxLevel TaxTree level constant (e.g. TaxTree.SPECIES, TaxTree.GENUS)
	 * @return Estimated probability of correctness (0-1)
	 */
	public static float probCorrect(int length, float k3dif, float k4dif, float k5dif, int taxLevel) {
		int idx = levelToIndex(taxLevel);
		if (idx < 0) { return -1; }
		double[] params = (k5dif < 1.0f) ? PROB_K5[idx] : PROB_K4[idx];
		float kdif = (k5dif < 1.0f) ? k5dif : k4dif;
		return predict(params, length, k3dif, kdif);
	}

	/**
	 * Probability that the hit organism is the same species as the query.
	 */
	public static float probSameSpecies(int length, float k3dif, float k4dif, float k5dif) {
		return probCorrect(length, k3dif, k4dif, k5dif, TaxTree.SPECIES);
	}

	/**
	 * Probability that the hit organism is in the same genus as the query.
	 */
	public static float probSameGenus(int length, float k3dif, float k4dif, float k5dif) {
		return probCorrect(length, k3dif, k4dif, k5dif, TaxTree.GENUS);
	}

	/**
	 * Probability that the hit organism is in the same family as the query.
	 */
	public static float probSameFamily(int length, float k3dif, float k4dif, float k5dif) {
		return probCorrect(length, k3dif, k4dif, k5dif, TaxTree.FAMILY);
	}

	/**
	 * Probability that the hit organism is in the same order as the query.
	 */
	public static float probSameOrder(int length, float k3dif, float k4dif, float k5dif) {
		return probCorrect(length, k3dif, k4dif, k5dif, TaxTree.ORDER);
	}

	/**
	 * Probability that the hit organism is in the same class as the query.
	 */
	public static float probSameClass(int length, float k3dif, float k4dif, float k5dif) {
		return probCorrect(length, k3dif, k4dif, k5dif, TaxTree.CLASS);
	}

	/**
	 * Probability that the hit organism is in the same phylum as the query.
	 */
	public static float probSamePhylum(int length, float k3dif, float k4dif, float k5dif) {
		return probCorrect(length, k3dif, k4dif, k5dif, TaxTree.PHYLUM);
	}

	/**
	 * Probability that the hit organism is in the same kingdom as the query.
	 */
	public static float probSameKingdom(int length, float k3dif, float k4dif, float k5dif) {
		return probCorrect(length, k3dif, k4dif, k5dif, TaxTree.KINGDOM);
	}

	/**
	 * Probability that the hit organism is in the same domain as the query.
	 */
	public static float probSameDomain(int length, float k3dif, float k4dif, float k5dif) {
		return probCorrect(length, k3dif, k4dif, k5dif, TaxTree.DOMAIN);
	}

	/*---------- internals ----------*/

	private static float predict(double[] params, int length, float k3dif, float kdif) {
		double floor = params[0];
		double log2len = Math.log(Math.max(length, 1)) * INV_LN2 - LOG2_REF;
		double z = params[1] + params[2] * log2len + params[3] * k3dif + params[4] * kdif;
		double sigmoid = 1.0 / (1.0 + Math.exp(z));
		return (float) (floor + (1.0 - floor) * sigmoid);
	}

	private static int levelToIndex(int taxLevel) {
		switch (taxLevel) {
			case TaxTree.SPECIES: return 0;
			case TaxTree.GENUS:   return 1;
			case TaxTree.FAMILY:  return 2;
			case TaxTree.ORDER:   return 3;
			case TaxTree.CLASS:   return 4;
			case TaxTree.PHYLUM:  return 5;
			case TaxTree.KINGDOM: return 6;
			case TaxTree.DOMAIN:  return 7;
			default: return -1;
		}
	}

	private static final double INV_LN2 = 1.0 / Math.log(2);
	private static final double LOG2_REF = Math.log(2500) * INV_LN2;

	// K4 model: p = floor + (1-floor) * sigmoid(a0 + a1*log2(len/2500) + a2*k3dif + a3*k4dif)
	// Rows: species, genus, family, order, class, phylum, kingdom, domain
	// Cols: floor, a0, a1_log2len, a2_k3dif, a3_k4dif
	static final double[][] PROB_K4 = {
		{0.0113446, -3.5630508, 0.6403488, 51.9226674, 13.6112786},
		{0.0795275, -5.0429731, 0.4555063, 59.6982223, 10.8967730},
		{0.2127537, -4.4862228, 0.2292117, 48.4028932,  7.7315173},
		{0.3256019, -5.2138624, 0.1925119, 49.9547608,  7.6156285},
		{0.5084922, -5.4920396, 0.0667606, 63.4318036,  0.5094376},
		{0.6071749, -6.5750378, 0.0950648, 57.9167172,  6.3128603},
		{0.6801771, -7.1822199, 0.1554127, 43.4878516, 15.1889807},
		{0.8182326, -8.3517371, 0.1448564, 32.9901037, 21.5209495},
	};

	// K5 model: p = floor + (1-floor) * sigmoid(a0 + a1*log2(len/2500) + a2*k3dif + a3*k5dif)
	static final double[][] PROB_K5 = {
		{0.0058409, -2.5413033, 0.4705238, 70.8753575, -0.7640568},
		{0.0556486, -3.6911274, 0.3030928, 71.3175058, -2.1744591},
		{0.1795990, -3.3621772, 0.1272263, 57.1529654, -2.5656541},
		{0.2788433, -4.2430087, 0.1167844, 57.7251810, -2.0374220},
		{0.4376369, -4.2088439,-0.0017810, 63.6425910, -5.5906412},
		{0.5352429, -5.0486971, 0.0025206, 64.2620784, -4.7840435},
		{0.6134413, -5.4231686, 0.0501707, 55.7766015, -2.2099322},
		{0.7916569, -7.1943760, 0.0995638, 49.9422788,  2.2346173},
	};
}
