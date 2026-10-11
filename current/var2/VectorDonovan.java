package var2;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import dna.Data;

/**
 * Legacy 16-feature extractor retained from Donovan's model work.
 * Uses fixed power/log transforms, two file-backed quantile tables and clipped
 * percentile scaling. Historical intent was Python/sklearn preprocessing parity;
 * this class does not establish that equivalence or validate a model schema.
 * Current Java scoring callers use VectorUMP45 directly.
 * @author Brian Bushnell
 */
//TODO: Possible bug [var2/VectorDonovan#001] (LOW, dead) - no callers found in current Java sources (superseded intern
//scaffolding; live vector is VectorUMP45). Also latent: the static QuantileTransformer fields below initialize at
//first active use of the outer class from resource files and would throw ExceptionInInitializerError if those files are absent - a
//trap if this class is ever re-wired. Cleanup candidate.
public class VectorDonovan{
	/**
	 * Interpolates between equally spaced ranks of a stored quantile table.
	 * The table is copied once and then read without mutation. Nonempty sorted,
	 * finite quantiles are a caller/resource obligation, not validated here.
	 * Duplicate entries and boundary ties follow the literal search and endpoint
	 * branches below; no exact sklearn-compatibility guarantee is made.
	 */
	public static class QuantileTransformer{
		/** Private copy of the resource values, expected in ascending order. */
		private final double[] quantiles;
		/** Historical feature label; retained but not used by the transform. */
		private final String featureName;

		/** Copies the supplied table without sorting or validating its entries. */
		private QuantileTransformer(double[] quantiles, String featureName){
			this.quantiles=quantiles.clone();
			this.featureName=featureName;
		}

		/**
		 * Loads the total-depth quantiles through Data.findPath and loadQuantiles.
		 * @return Fresh transformer owning a copy of the loaded values
		 */
		public static QuantileTransformer forTotalDepth(){
			String path=Data.findPath("?total_depth_quantiles.txt");
			double[] quantiles=loadQuantiles(path);
			return new QuantileTransformer(quantiles, "Total_Depth");
		}

		/**
		 * Loads the mean-end-distance quantiles through Data.findPath and loadQuantiles.
		 * @return Fresh transformer owning a copy of the loaded values
		 */
		public static QuantileTransformer forEndDistanceAverage(){
			String path=Data.findPath("?end_distance_quantiles.txt");
			double[] quantiles=loadQuantiles(path);
			return new QuantileTransformer(quantiles, "End_Distance_Average");
		}

		/**
		 * Tries an absolute class resource name, then the supplied filesystem path.
		 * Reads using the default charset and closes the reader/stream after parsing.
		 * Blank and comment lines are ignored; remaining trimmed lines are parsed as
		 * one double each. Empty, unsorted and nonfinite tables are not rejected.
		 * @param filename Resolved non-null path, also tried with a leading slash as a resource
		 * @return Newly allocated values in file order
		 * @throws RuntimeException if opening, reading or numeric parsing fails
		 * @throws NullPointerException if filename is null and no resource is found
		 */
		private static double[] loadQuantiles(String filename){
			try{
				// Try to load from resources first, then from file system
				InputStream is=QuantileTransformer.class.getResourceAsStream("/"+filename);
				if(is==null){is=new FileInputStream(filename);}

				List<Double> quantileList=new ArrayList<>();
				try(BufferedReader reader=new BufferedReader(new InputStreamReader(is))){
					String line;
					while((line=reader.readLine())!=null){
						line=line.trim();
						if(!line.isEmpty() && !line.startsWith("#")){
							quantileList.add(Double.parseDouble(line));
						}
					}
				}

				double[] result=new double[quantileList.size()];
				for(int i=0; i<quantileList.size(); i++){
					result[i]=quantileList.get(i).doubleValue();
				}
				return result;

			}catch(IOException | NumberFormatException e){
				throw new RuntimeException("Failed to load quantiles from "+filename, e);
			}
		}

		/**
		 * Maps values at/below the first quantile to zero and at/above the last to one.
		 * Otherwise interpolates rank between the enclosing quantiles. Requires a
		 * nonempty sorted finite table; NaN input is not rejected and may propagate.
		 * The first boundary test wins when both endpoints are equal.
		 * @param value Raw feature value
		 * @return Interpolated rank; [0,1] assumes ordered data and finite intermediate arithmetic
		 */
		public double transform(double value){
			// Handle edge cases
			if(value<=quantiles[0]){return 0.0;}
			if(value>=quantiles[quantiles.length-1]){return 1.0;}

			// Find the quantile interval containing this value
			int leftIndex=findQuantileIndex(value);

			if(leftIndex==quantiles.length-1){return 1.0;}

			// Linear interpolation between quantiles
			double leftQuantile=quantiles[leftIndex];
			double rightQuantile=quantiles[leftIndex+1];
			double leftRank=leftIndex/(double)(quantiles.length-1);
			double rightRank=(leftIndex+1)/(double)(quantiles.length-1);

			// Interpolate the rank
			double fraction=(value-leftQuantile)/(rightQuantile-leftQuantile);
			return leftRank+fraction*(rightRank-leftRank);
		}

		/**
		 * Binary search to find the largest quantile index where {@code quantiles[index]<=value}.
		 * @param value Value to locate in quantile array
		 * @return Index of quantile interval containing the value
		 */
		private int findQuantileIndex(double value){
			int left=0;
			int right=quantiles.length-1;

			while(left<right){
				int mid=left+(right-left+1)/2;
				if(quantiles[mid]<=value){
					left=mid;
				}else{
					right=mid-1;
				}
			}
			return left;
		}
	}

	/**
	 * Allocates raw, transformed and float arrays for the 16-feature layout.
	 * Revised allele fraction may be read from or written to the variant's cache.
	 * Mean accessors divide by allele count without a zero-count guard; this
	 * extractor does not reject nonfinite results. Quantile entries use table ranks;
	 * other entries are scaled and clipped to [-0.5,1.25] when not NaN.
	 * @param v Non-null variant with accumulated statistics
	 * @param pairingRate Unused compatibility parameter
	 * @param totalQualityAvg Unused compatibility parameter
	 * @param totalMapqAvg Unused compatibility parameter
	 * @param readLengthAvg Mean read length in bases for revised allele fraction
	 * @param ploidy Unused compatibility parameter
	 * @param map Optional scaffold map for end-distance calculation
	 * @return Newly allocated 16-element vector
	 */
	public static float[] makeDonovanVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){
		// Extract raw 16 features
		double[] rawFeatures=extractDonovanRawFeatures(v, readLengthAvg, map);

		// Apply complete normalization pipeline
		double[] normalizedFeatures=normalizeDonovanFeatures(rawFeatures);

		// Convert to float array for consistency with existing interface
		float[] vec=new float[16];
		for(int i=0; i<16; i++){
			vec[i]=(float)normalizedFeatures[i];
		}
		return vec;
	}

	/**
	 * Extracts raw values in the order named at each assignment below.
	 * Mean identity retains per-mille units, while maximum identity is converted
	 * to a fraction. Mean accessors require positive supporting count for finite
	 * results. Revised AF may update the variant's insertion cache. End-distance
	 * context uses Var.nScan when scanning is disabled, v.start without a map,
	 * or v.contigEndDist with a map.
	 * @param v Non-null variant
	 * @param readLengthAvg Mean read length in bases for revised AF
	 * @param scafMap Optional scaffold map
	 * @return Fresh 16-element raw array
	 */
	private static double[] extractDonovanRawFeatures(Var v, double readLengthAvg, ScafMap scafMap){
		final int scafEndDist=!Var.doNscan ? Var.nScan :
			(scafMap==null ? v.start : v.contigEndDist(scafMap));

		double[] features=new double[16];

		// 1. Total_Depth (DP field in VCF)
		features[0]=v.coverage();

		// 2. Allelic_Depth
		features[1]=v.alleleCount();

		// 3. Revised_Allele_Fraction
		features[2]=v.revisedAlleleFraction(v.alleleFraction(), readLengthAvg);

		// 4. Nearby_Variant_Count
		features[3]=v.nearbyVarCount;

		// 5. Mapping_Quality_Average
		features[4]=v.mapQAvg();

		// 6. Mapping_Quality_Max
		features[5]=v.mapQMax;

		// 7. Base_Quality_Max
		features[6]=v.baseQMax;

		// 8. Base_Quality_Average
		features[7]=v.baseQAvg();

		// 9. Identity_Average
		features[8]=v.identityAvg();

		//TODO: Probable bug [VAR2-015] - this fraction is later raised to power 7, then
		//scaled against literal 144..992 endpoints, unlike the per-mille mean above.
		//Check original training preprocessing before changing units; no current caller found.
		// 10. Identity_Max
		features[9]=v.idMax*0.001; // Convert from per-mille to fraction

		// 11. End_Distance_Max
		features[10]=v.endDistMax;

		// 12. End_Distance_Average
		features[11]=v.edistAvg();

		// 13. Length_Average
		features[12]=v.lengthAvg();

		// 14. Strand_Bias (SB field in VCF)
		features[13]=v.strandBiasScore(scafEndDist);

		// 15. Read1_Plus_Count
		features[14]=v.r1plus;

		// 16. Read1_Minus_Count
		features[15]=v.r1minus;
		return features;
	}

	/** Loaded during outer-class initialization; failure prevents active use. */
	private static final QuantileTransformer TOTAL_DEPTH_TRANSFORMER=QuantileTransformer.forTotalDepth();
	/** Loaded after total depth; the same initialization-failure contract applies. */
	private static final QuantileTransformer END_DISTANCE_TRANSFORMER=QuantileTransformer.forEndDistanceAverage();

	/**
	 * Historical training endpoints, originally labeled {p0, p99} per feature.
	 * Applied after nonlinear transforms; rows 0 and 11 are skipped. Compatibility
	 * with the original training pipeline must be established before reuse.
	 */
	private static final double[][] PERCENTILES={
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
	 * Copies the first 16 values, applies log/power and quantile transforms, then
	 * scales the non-quantile entries with literal percentile endpoints and clips
	 * them to [-0.5,1.25]. Equal endpoints leave the transformed entry unchanged.
	 * No finite-value validation occurs; NaN survives the clipping operations.
	 * @param features Raw array with at least 16 elements; not modified
	 * @return Newly allocated transformed array
	 */
	private static double[] normalizeDonovanFeatures(double[] features){
		double[] working=new double[16];
		System.arraycopy(features, 0, working, 0, 16);

		// Step 1: Log transformations (log1p = log(1+x))
		working[1]=Math.log1p(working[1]); // Allelic_Depth
		working[3]=Math.log1p(working[3]); // Nearby_Variant_Count
		working[14]=Math.log1p(working[14]); // Read1_Plus_Count
		working[15]=Math.log1p(working[15]); // Read1_Minus_Count

		// Step 2: Power transformations
		working[12]=Math.pow(working[12], 10); // Length_Average^10
		working[8]=Math.pow(working[8], 8); // Identity_Average^8
		working[5]=Math.pow(working[5], 5); // Mapping_Quality_Max^5
		working[10]=Math.pow(working[10], 4); // End_Distance_Max^4
		working[6]=Math.pow(working[6], 5); // Base_Quality_Max^5
		working[4]=Math.pow(working[4], 3); // Mapping_Quality_Average^3
		working[9]=Math.pow(working[9], 7); // Identity_Max^7
		working[13]=Math.pow(working[13], 100); // Strand_Bias^100

		// Step 3: Quantile transformations for specific features
		working[0]=TOTAL_DEPTH_TRANSFORMER.transform(working[0]); // Total_Depth
		working[11]=END_DISTANCE_TRANSFORMER.transform(working[11]); // End_Distance_Average

		// Step 4: Percentile scaling with clipping for remaining features
		for(int i=0; i<16; i++){
			// Skip features that already had quantile transformation
			if(i==0 || i==11){continue;}

			double min_val=PERCENTILES[i][0];
			double max_val=PERCENTILES[i][1];

			if(max_val-min_val==0){continue;} // Equal endpoints: leave unchanged

			// Scale relative to the literal percentile endpoints; clipping follows
			working[i]=(working[i]-min_val)/(max_val-min_val);

			// Clip to [-0.5, 1.25] range
			working[i]=Math.max(-0.5, Math.min(1.25, working[i]));
		}
		return working;
	}
}
