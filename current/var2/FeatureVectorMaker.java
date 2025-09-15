package var2;

import shared.*;

/**
 * Converts Var objects into feature vectors for neural network analysis.
 * Implements a factory pattern supporting three different neural network models
 * developed by summer interns, each with specialized feature extraction algorithms.
 * 
 * Feature extraction modes:
 * - ELBA: Quality score prediction using 32-dimensional vectors with quality metrics,
 *   strand bias analysis, and context normalization against dataset statistics
 * - LAWRENCE: Genotype calling using 8-dimensional min-max scaled vectors with
 *   revised allele fractions, mapping qualities, and strand bias scoring
 * - DONOVAN: Template framework for additional model development (16-dimensional placeholder)
 * 
 * The feature vectors are designed for specific neural network architectures trained
 * on variant calling datasets, with each intern optimizing for their particular task.
 * 
 * @author Brian Bushnell
 * @author Isla
 * @date July 21, 2025
 */
public class FeatureVectorMaker {
	
	/*--------------------------------------------------------------*/
	/*----------------        Mode Constants        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Feature extraction modes for different neural network models trained by summer interns */
	public static final int ELBA=0, LAWRENCE=1, DONOVAN=2;
	
	/** Current active feature extraction mode, defaults to ELBA's quality prediction model */
	private static int MODE=ELBA;
	
	/** String representations of modes for command-line parsing and configuration */
	private static final String[] MODE_NAMES={"ELBA", "LAWRENCE", "DONOVAN"};
	
	/*--------------------------------------------------------------*/
	/*----------------      Main Entry Point       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Convert a variant to a feature vector using the currently selected extraction method.
	 * Implements factory pattern routing to intern-specific feature extraction algorithms.
	 * Each algorithm produces vectors of different dimensionality and scaling appropriate
	 * for their trained neural network architectures.
	 * 
	 * @param v Variant to convert with coverage, quality, and mapping statistics
	 * @param pairingRate Dataset-wide proper pairing rate for normalization context
	 * @param totalQualityAvg Average base quality across entire sequencing dataset
	 * @param totalMapqAvg Average mapping quality across entire sequencing dataset
	 * @param readLengthAvg Average read length for the sequencing technology used
	 * @param ploidy Expected sample ploidy (1 for haploid, 2 for diploid organisms)
	 * @param map Scaffold mapping for calculating positional features and end distances
	 * @return Float array feature vector ready for neural network input layer
	 */
	public static float[] toVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){
		
		if(MODE==ELBA){
			return makeElbaVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, map);
		}else if(MODE==LAWRENCE){
			return makeLawrenceVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, map);
		}else if(MODE==DONOVAN){
			return makeDonovanVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, map);
		}else{
			throw new RuntimeException("Unknown feature vector mode: "+MODE);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Intern-Specific Methods   ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Generate Elba's 32-dimensional feature vector optimized for quality score prediction.
	 * Implements quality analysis combining raw statistics, averages, and normalization
	 * against dataset-wide metrics. Features include variant properties, quality metrics
	 * (both raw and averaged), strand bias analysis, pairing statistics, and contextual
	 * normalization ratios.
	 * Algorithm combines absolute quality metrics with relative measures normalized
	 * against dataset averages to provide context-aware features for neural network training.
	 * @param v Variant with accumulated quality and coverage statistics
	 * @param pairingRate Dataset proper pairing rate for normalization
	 * @param totalQualityAvg Dataset average base quality for relative scaling
	 * @param totalMapqAvg Dataset average mapping quality for relative scaling
	 * @param readLengthAvg Technology-specific average read length
	 * @param ploidy Expected sample ploidy for genotype context
	 * @param map Scaffold mapping for positional feature calculation
	 * @return 32-element float array with normalized quality prediction features
	 */
	private static float[] makeElbaVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){
		
		//Elba's vector focuses on quality prediction - allocate appropriate size
		float[] vec=new float[32]; //32-element vector matches Elba's network architecture
		int idx=0;
		
		//Basic variant properties
		vec[idx++]=(float)v.type(); //Variant type (SUB, INS, DEL)
		vec[idx++]=(float)v.alleleCount(); //Supporting read count
		vec[idx++]=(float)v.coverage(); //Total coverage depth
		vec[idx++]=(float)v.alleleFraction(); //Allele frequency
		
		//Quality metrics (raw values)
		vec[idx++]=(float)v.baseQMax; //Maximum base quality
		vec[idx++]=(float)v.mapQMax; //Maximum mapping quality
		vec[idx++]=(float)v.idMax*0.001f; //Maximum identity scaled to [0,1]
		vec[idx++]=(float)v.endDistMax; //Maximum end distance
		
		//Quality metrics (averages)
		final int count=v.alleleCount();
		vec[idx++]=(count>0 ? (float)(v.baseQSum/(double)count) : 0f); //Average base quality
		vec[idx++]=(count>0 ? (float)(v.mapQSum/(double)count) : 0f); //Average mapping quality
		vec[idx++]=(count>0 ? (float)(v.idSum/(double)count*0.001) : 0f); //Average identity scaled to [0,1]
		vec[idx++]=(count>0 ? (float)(v.endDistSum/(double)count) : 0f); //Average end distance
		
		//Strand bias and pairing
		vec[idx++]=(float)v.strandRatio(); //Strand balance ratio
		vec[idx++]=(count>0 ? (float)(v.properPairCount/(double)count) : 0f); //Proper pair rate
		
		//Context features relative to dataset averages
		vec[idx++]=(float)(totalQualityAvg>0 ? (v.baseQSum/(double)count)/totalQualityAvg : 1.0); //Quality ratio vs dataset
		vec[idx++]=(float)(totalMapqAvg>0 ? (v.mapQSum/(double)count)/totalMapqAvg : 1.0); //MapQ ratio vs dataset
		vec[idx++]=(float)(pairingRate>0 ? (v.properPairCount/(double)count)/pairingRate : 1.0); //Pairing ratio vs dataset
		
		//Positional features
		vec[idx++]=(float)ploidy; //Sample ploidy
		vec[idx++]=(float)readLengthAvg; //Average read length
		
		//Nearby variant context (if available)
		vec[idx++]=(float)v.nearbyVarCount; //Nearby variant density
		
		//Reserved slots for additional features
		while(idx<vec.length){vec[idx++]=0f;}
		
		return vec;
	}
	
	/**
	 * Generate Lawrence's 8-dimensional feature vector optimized for genotype calling.
	 * Implements min-max scaling algorithm with pre-trained bounds derived from training data.
	 * Features focus on revised allele fractions, quality metrics, and strand bias scoring
	 * specifically tuned for distinguishing true variants from sequencing artifacts.
	 * 
	 * The algorithm applies min-max normalization: (value - min) / (max - min)
	 * with bounds clipped to [0,1] to handle values outside the training range.
	 * Feature bounds were empirically derived from variant calling training datasets.
	 * 
	 * @param v Variant with quality and mapping statistics
	 * @param pairingRate Dataset proper pairing rate (unused in current implementation)
	 * @param totalQualityAvg Dataset average base quality (unused in current implementation)
	 * @param totalMapqAvg Dataset average mapping quality (unused in current implementation)
	 * @param readLengthAvg Average read length for revised allele fraction calculation
	 * @param ploidy Sample ploidy (unused in current implementation)
	 * @param map Scaffold mapping for strand bias end distance calculation
	 * @return 8-element float array with min-max scaled genotype calling features
	 */
	private static float[] makeLawrenceVector(Var v, double pairingRate, double totalQualityAvg, 
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){

		float[] vec=new float[8];

		//1. RAF: Revised Allele Fraction
		vec[0]=(float)v.revisedAlleleFraction(v.alleleFraction(), readLengthAvg);
		//2. AVG_MQS: Average Map Quality Score
		vec[1]=(float)v.mapQAvg();
		//3. MQM: Mapping Quality Max
		vec[2]=(float)v.mapQMax;
		//4. AVG_IDS: Average Identity Sum
		vec[3]=(float)v.identityAvg();
		//5. EDM: End Distance Max
		vec[4]=(float)v.endDistMax;
		//6. BQM: Base Quality Max
		vec[5]=(float)v.baseQMax;
		//7. AVG_LS: Average Length Sum
		vec[6]=(float)v.lengthAvg();
		//8. SB: Strand Bias
		final int scafEndDist=!Var.doNscan ? Var.nScan : (map==null ? v.start : v.contigEndDist(map));
		vec[7]=(float)v.strandBiasScore(scafEndDist);
		
		//Min-max scaling with empirically-derived bounds from training data
		float[] scaledVec=new float[8];

		final float[] FEATURE_MIN_VALUES={
			0.050000f, //RAF
			1.000000f, //AVG_MQS
			1.000000f, //MQM
			144.000000f, //AVG_IDS
			5.000000f, //EDM
			5.000000f, //BQM
			3.000000f, //AVG_LS
			0.000000f //SB
		};

		final float[] FEATURE_MAX_VALUES={
			1.000000f, //RAF
			43.000000f, //AVG_MQS
			43.000000f, //MQM
			992.000000f, //AVG_IDS
			74.000000f, //EDM
			41.000000f, //BQM
			138.000000f, //AVG_LS
			1.000000f //SB
		};

		for(int i=0; i<vec.length; i++){
			float rawValue=vec[i];
			float min=FEATURE_MIN_VALUES[i];
			float max=FEATURE_MAX_VALUES[i];

			//The scaling formula: (value - min) / (max - min)
			float denominator=max-min;

			if(denominator==0){
				//Edge case: If max equals min, the feature is constant
				scaledVec[i]=0.0f;
			}else{
				float scaledValue=(rawValue-min)/denominator;
				//Clip the value to be strictly between 0 and 1
				//This handles cases where a new variant's feature is outside the training range
				scaledVec[i]=Math.max(0.0f, Math.min(1.0f, scaledValue));
			}
		}

		return scaledVec;
	}
	
	/**
	 * Generate Donovan's 16-dimensional feature vector template for future model development.
	 * Currently implements basic variant features as a foundation for custom neural network
	 * architecture development. This method serves as a starting framework that can be
	 * extended with domain-specific features as requirements become defined.
	 * 
	 * Template includes fundamental variant properties (allele count, coverage, frequency,
	 * and type) with remaining vector slots reserved for custom feature engineering.
	 * 
	 * @param v Variant with basic statistics for feature extraction
	 * @param pairingRate Dataset proper pairing rate (reserved for future use)
	 * @param totalQualityAvg Dataset average base quality (reserved for future use)
	 * @param totalMapqAvg Dataset average mapping quality (reserved for future use)
	 * @param readLengthAvg Average read length (reserved for future use)
	 * @param ploidy Sample ploidy (reserved for future use)
	 * @param map Scaffold mapping (reserved for future use)
	 * @return 16-element float array with basic features and reserved slots for expansion
	 */
	private static float[] makeDonovanVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){
		
		//Placeholder for third intern's feature extraction
		float[] vec=new float[16]; //Adjust based on requirements
		int idx=0;
		
		//Basic features as starting point
		vec[idx++]=(float)v.alleleCount();
		vec[idx++]=(float)v.coverage();
		vec[idx++]=(float)v.alleleFraction();
		vec[idx++]=(float)v.type();
		
		//Donovan can customize this method for their specific model
		//Reserved slots for custom features
		while(idx<vec.length){vec[idx++]=0f;}
		
		return vec;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Utility Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Set the active feature extraction mode for all subsequent toVector() calls.
	 * Changes the global MODE variable to route to the specified intern's algorithm.
	 * @param mode Feature extraction mode constant (ELBA=0, LAWRENCE=1, DONOVAN=2)
	 */
	public static void setMode(int mode){
		assert(mode>=0 && mode<=2) : "Invalid mode: "+mode;
		MODE=mode;
	}
	
	/**
	 * Set the active feature extraction mode from string representation.
	 * Parses case-insensitive mode names and converts to internal mode constants.
	 * @param modeStr Mode name string ("ELBA", "LAWRENCE", "DONOVAN", case-insensitive)
	 * @throws RuntimeException if modeStr does not match any known mode name
	 */
	public static void setMode(String modeStr){
		int mode=Tools.find(modeStr.toUpperCase(), MODE_NAMES);
		if(mode<0){throw new RuntimeException("Unknown feature mode: "+modeStr);}
		MODE=mode;
	}
	
	/**
	 * Get the currently active feature extraction mode constant.
	 * @return Current mode constant (ELBA=0, LAWRENCE=1, DONOVAN=2)
	 */
	public static int getMode(){
		return MODE;
	}
	
	/**
	 * Get the expected feature vector length for the currently active mode.
	 * Each intern's neural network architecture requires different input dimensions.
	 * @return Expected vector length (ELBA=32, LAWRENCE=8, DONOVAN=16)
	 */
	public static int getVectorLength(){
		if(MODE==ELBA){return 32;}
		else if(MODE==LAWRENCE){return 8;} //Possible bug: was 24, but makeLawrenceVector returns 8 elements
		else if(MODE==DONOVAN){return 16;}
		else{throw new RuntimeException("Unknown mode: "+MODE);}
	}
	
	/**
	 * Get the string representation of the currently active mode.
	 * @return Current mode name ("ELBA", "LAWRENCE", or "DONOVAN")
	 */
	public static String getModeName(){
		return MODE_NAMES[MODE];
	}
}