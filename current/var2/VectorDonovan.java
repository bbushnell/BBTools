package var2;

import java.io.*;
import java.util.*;

import dna.Data;

public class VectorDonovan {
	/**
	 * Java implementation of sklearn's QuantileTransformer with uniform output distribution.
	 * Replicates the exact transformation used in Python for neural network feature preprocessing.
	 */
	public static class QuantileTransformer {
		private final double[] quantiles;
		private final String featureName;

		private QuantileTransformer(double[] quantiles, String featureName) {
			this.quantiles = quantiles.clone();
			this.featureName = featureName;
		}

		/**
		 * Constructor for Total_Depth transformer
		 */
		public static QuantileTransformer forTotalDepth() {
			String path=Data.findPath("?total_depth_quantiles.txt");
			double[] quantiles = loadQuantiles(path);
			return new QuantileTransformer(quantiles, "Total_Depth");
		}

		/**
		 * Constructor for End_Distance_Average transformer
		 */
		public static QuantileTransformer forEndDistanceAverage() {
			String path=Data.findPath("?end_distance_quantiles.txt");
			double[] quantiles = loadQuantiles(path);
			return new QuantileTransformer(quantiles, "End_Distance_Average");
		}

		/**
		 * Loads quantile values from a text file
		 */
		private static double[] loadQuantiles(String filename) {
			try {
				// Try to load from resources first, then from file system
				InputStream is = QuantileTransformer.class.getResourceAsStream("/" + filename);
				if (is == null) {
					is = new FileInputStream(filename);
				}

				List<Double> quantileList = new ArrayList<>();
				try (BufferedReader reader = new BufferedReader(new InputStreamReader(is))) {
					String line;
					while ((line = reader.readLine()) != null) {
						line = line.trim();
						if (!line.isEmpty() && !line.startsWith("#")) {
							quantileList.add(Double.parseDouble(line));
						}
					}
				}

				double[] result = new double[quantileList.size()];
				for(int i = 0; i < quantileList.size(); i++){
					result[i] = quantileList.get(i).doubleValue();
				}
				return result;

			} catch (IOException | NumberFormatException e) {
				throw new RuntimeException("Failed to load quantiles from " + filename, e);
			}
		}

		/**
		 * Transforms input value using quantile mapping to uniform distribution [0,1]
		 */
		public double transform(double value) {
			// Handle edge cases
			if (value <= quantiles[0]) return 0.0;
			if (value >= quantiles[quantiles.length - 1]) return 1.0;

			// Find the quantile interval containing this value
			int leftIndex = findQuantileIndex(value);

			if (leftIndex == quantiles.length - 1) {
				return 1.0; // Value equals the maximum quantile
			}

			// Linear interpolation between quantiles
			double leftQuantile = quantiles[leftIndex];
			double rightQuantile = quantiles[leftIndex + 1];
			double leftRank = leftIndex / (double)(quantiles.length - 1);
			double rightRank = (leftIndex + 1) / (double)(quantiles.length - 1);

			// Interpolate the rank
			double fraction = (value - leftQuantile) / (rightQuantile - leftQuantile);
			return leftRank + fraction * (rightRank - leftRank);
		}

		/**
		 * Binary search to find the largest quantile index where quantiles[index] <= value
		 */
		private int findQuantileIndex(double value) {
			int left = 0;
			int right = quantiles.length - 1;

			while (left < right) {
				int mid = left + (right - left + 1) / 2;
				if (quantiles[mid] <= value) {
					left = mid;
				} else {
					right = mid - 1;
				}
			}
			return left;
		}
	}

	/**
	 * Donovan's feature vector for quality score prediction and true/false classification.
	 * Uses 16 features with normalization pipeline including quantile transformations.
	 * 
	 * @param v Variant to convert
	 * @param pairingRate Overall proper pairing rate (unused in Donovan's model)
	 * @param totalQualityAvg Average base quality from dataset (unused in Donovan's model)
	 * @param totalMapqAvg Average mapping quality from dataset (unused in Donovan's model)
	 * @param readLengthAvg Average read length
	 * @param ploidy Sample ploidy (unused in Donovan's model)
	 * @param map Scaffold mapping (unused in Donovan's model)
	 * @return 16-element normalized feature vector for Donovan's neural network
	 */
	public static float[] makeDonovanVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){

		// Extract raw 16 features
		double[] rawFeatures = extractDonovanRawFeatures(v, readLengthAvg, map);

		// Apply complete normalization pipeline
		double[] normalizedFeatures = normalizeDonovanFeatures(rawFeatures);

		// Convert to float array for consistency with existing interface
		float[] vec = new float[16];
		for (int i = 0; i < 16; i++) {
			vec[i] = (float) normalizedFeatures[i];
		}

		return vec;
	}

	/**
	 * Extracts the 16 raw features for Donovan's model.
	 */
	private static double[] extractDonovanRawFeatures(Var v, double readLengthAvg, ScafMap scafMap) {

		final int scafEndDist = !Var.doNscan ? Var.nScan : 
			(scafMap == null ? v.start : v.contigEndDist(scafMap));

		double[] features = new double[16];

		// 1. Total_Depth (DP field in VCF)
		features[0] = v.coverage();

		// 2. Allelic_Depth 
		features[1] = v.alleleCount();

		// 3. Revised_Allele_Fraction
		features[2] = v.revisedAlleleFraction(v.alleleFraction(), readLengthAvg);

		// 4. Nearby_Variant_Count
		features[3] = v.nearbyVarCount;

		// 5. Mapping_Quality_Average
		features[4] = v.mapQAvg();

		// 6. Mapping_Quality_Max
		features[5] = v.mapQMax;

		// 7. Base_Quality_Max
		features[6] = v.baseQMax;

		// 8. Base_Quality_Average
		features[7] = v.baseQAvg();

		// 9. Identity_Average
		features[8] = v.identityAvg();

		// 10. Identity_Max
		features[9] = v.idMax * 0.001; // Convert from per-mille to fraction

		// 11. End_Distance_Max
		features[10] = v.endDistMax;

		// 12. End_Distance_Average
		features[11] = v.edistAvg();

		// 13. Length_Average
		features[12] = v.lengthAvg();

		// 14. Strand_Bias (SB field in VCF)
		features[13] = v.strandBiasScore(scafEndDist);

		// 15. Read1_Plus_Count
		features[14] = v.r1plus;

		// 16. Read1_Minus_Count
		features[15] = v.r1minus;

		return features;
	}

	/**
	 * Cached transformer instances for performance
	 */
	private static final QuantileTransformer TOTAL_DEPTH_TRANSFORMER = QuantileTransformer.forTotalDepth();
	private static final QuantileTransformer END_DISTANCE_TRANSFORMER = QuantileTransformer.forEndDistanceAverage();

	/**
	 * Percentile parameters from training data {p0, p99}
	 */
	private static final double[][] PERCENTILES = {
			{2.000000, 3127.000000}, // 0: Total_Depth
			{2.000000, 481.000000}, // 1: Allelic_Depth  
			{0.050000, 1.000000}, // 2: Revised_Allele_Fraction
			{0.000000, 20.000000}, // 3: Nearby_Variant_Count
			{1.000000, 42.982500}, // 4: Mapping_Quality_Average
			{1.000000, 43.000000}, // 5: Mapping_Quality_Max
			{5.000000, 41.000000}, // 6: Base_Quality_Max
			{5.000000, 39.000000}, // 7: Base_Quality_Average
			{144.000000, 992.000000}, // 8: Identity_Average
			{144.000000, 992.000000}, // 9: Identity_Max
			{5.000000, 73.000000}, // 10: End_Distance_Max
			{5.000000, 63.833300}, // 11: End_Distance_Average
			{3.000000, 138.000000}, // 12: Length_Average
			{0.000000, 1.000000}, // 13: Strand_Bias
			{0.000000, 125.000000}, // 14: Read1_Plus_Count
			{0.000000, 116.000000} // 15: Read1_Minus_Count
	};

	/**
	 * Applies Donovan's complete normalization pipeline.
	 * Replicates the Python normalization sequence used in training.
	 */
	private static double[] normalizeDonovanFeatures(double[] features) {

		double[] working = new double[16];
		System.arraycopy(features, 0, working, 0, 16);

		// Step 1: Log transformations (log1p = log(1+x))
		working[1] = Math.log1p(working[1]); // Allelic_Depth
		working[3] = Math.log1p(working[3]); // Nearby_Variant_Count  
		working[14] = Math.log1p(working[14]); // Read1_Plus_Count
		working[15] = Math.log1p(working[15]); // Read1_Minus_Count

		// Step 2: Power transformations
		working[12] = Math.pow(working[12], 10); // Length_Average^10
		working[8] = Math.pow(working[8], 8); // Identity_Average^8
		working[5] = Math.pow(working[5], 5); // Mapping_Quality_Max^5
		working[10] = Math.pow(working[10], 4); // End_Distance_Max^4
		working[6] = Math.pow(working[6], 5); // Base_Quality_Max^5
		working[4] = Math.pow(working[4], 3); // Mapping_Quality_Average^3
		working[9] = Math.pow(working[9], 7); // Identity_Max^7
		working[13] = Math.pow(working[13], 100); // Strand_Bias^100

		// Step 3: Quantile transformations for specific features
		working[0] = TOTAL_DEPTH_TRANSFORMER.transform(working[0]); // Total_Depth
		working[11] = END_DISTANCE_TRANSFORMER.transform(working[11]); // End_Distance_Average

		// Step 4: Percentile scaling with clipping for remaining features
		for (int i = 0; i < 16; i++) {
			// Skip features that already had quantile transformation
			if (i == 0 || i == 11) {
				continue;
			}

			double min_val = PERCENTILES[i][0];
			double max_val = PERCENTILES[i][1];

			if (max_val - min_val == 0) {
				continue; // No variance, leave unchanged
			}

			// Scale to 0-1 range using percentiles
			working[i] = (working[i] - min_val) / (max_val - min_val);

			// Clip to [-0.5, 1.25] range
			working[i] = Math.max(-0.5, Math.min(1.25, working[i]));
		}

		return working;
	}
}