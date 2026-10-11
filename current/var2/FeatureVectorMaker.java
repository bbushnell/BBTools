package var2;

import shared.Tools;

/**
 * Legacy dispatcher for four feature layouts, retained from the intern models.
 * ELBA has 32 slots (20 populated), LAWRENCE has eight scaled features, DONOVAN
 * has 16 slots (four populated), and UMP45 delegates to VectorUMP45 (33 slots).
 * The private ELBA and DONOVAN implementations are distinct from the standalone
 * VectorElba and VectorDonovan classes. The current Java scoring path bypasses
 * this dispatcher and calls VectorUMP45 directly; see the historical issue below.
 *
 * Mode is process-global and unsynchronized; configure before use, not during
 * concurrent extraction. This class neither selects a network nor validates its
 * input schema. LAWRENCE and UMP45 may update the variant's revised-AF cache.
 * @author Brian Bushnell
 * @author Isla
 * @date July 21, 2025
 */
public class FeatureVectorMaker{
	/*--------------------------------------------------------------*/
	/*----------------        Mode Constants        ----------------*/
	/*--------------------------------------------------------------*/

	/** Stable indices into MODE_NAMES. */
	public static final int ELBA=0, LAWRENCE=1, DONOVAN=2, UMP45=3;

	/** Legacy selection, initially ELBA; not read by current direct scoring callers. */
	private static int MODE=ELBA;

	/** Names in the same order as the mode constants. */
	private static final String[] MODE_NAMES={"ELBA", "LAWRENCE", "DONOVAN", "UMP45"};

	/*--------------------------------------------------------------*/
	/*----------------      Main Entry Point       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Allocates a vector using the process-global mode and its historical layout.
	 * No finiteness or network-dimension check is added by the dispatcher; UMP45
	 * retains its own checks. Parameter usage depends on the selected branch.
	 * @param v Non-null variant with coverage and supporting-read statistics
	 * @param pairingRate Dataset proper-pair fraction for ELBA or UMP45
	 * @param totalQualityAvg Dataset mean base quality for ELBA or UMP45
	 * @param totalMapqAvg Dataset mean mapping quality for ELBA or UMP45
	 * @param readLengthAvg Dataset mean read length in bases; unused by DONOVAN
	 * @param ploidy Sample ploidy for ELBA or UMP45
	 * @param map Optional scaffold context for LAWRENCE or UMP45
	 * @return Fresh array in the selected layout
	 * @throws RuntimeException if MODE contains an unsupported value
	 */
	//TODO: Possible bug [var2/FeatureVectorMaker#001] (LOW) - no toVector callers found in current
	//Java sources; live scoring calls VectorUMP45.makeVector directly. netmode= parsing sets MODE
	//and useNet, but the chosen layout is not consulted by that scoring path. Rewiring would require
	//compatible model dimensions and preprocessing. Cleanup/architecture decision remains for Brian.
	public static float[] toVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){
		if(MODE==ELBA){
			return makeElbaVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, map);
		}else if(MODE==LAWRENCE){
			return makeLawrenceVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, map);
		}else if(MODE==DONOVAN){
			return makeDonovanVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, map);
		}else if(MODE==UMP45){
			return VectorUMP45.makeVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, map);
		}else{
			throw new RuntimeException("Unknown feature vector mode: "+MODE);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Intern-Specific Methods   ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds the historical 32-slot ELBA layout: 20 features followed by 12 zeros.
	 * Combines raw statistics, guarded means and three dataset-relative ratios.
	 * The ratios guard dataset denominators but not allele count (see issue #002).
	 * This is not the standalone VectorElba 11-feature preprocessing pipeline.
	 * @param v Non-null variant with accumulated statistics
	 * @param pairingRate Dataset proper-pair fraction; nonpositive values use ratio 1
	 * @param totalQualityAvg Dataset mean base quality; nonpositive values use ratio 1
	 * @param totalMapqAvg Dataset mean mapping quality; nonpositive values use ratio 1
	 * @param readLengthAvg Dataset mean read length, emitted directly
	 * @param ploidy Sample ploidy, emitted directly
	 * @param map Unused compatibility parameter
	 * @return Fresh 32-element array; no finiteness or range validation
	 */
	private static float[] makeElbaVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){
		//Elba's vector focuses on quality prediction - allocate appropriate size
		float[] vec=new float[32]; //Historical layout; no model-dimension check here
		int idx=0;

		//Basic variant properties
		vec[idx++]=(float)v.type(); //Variant type (SUB, INS, DEL)
		vec[idx++]=(float)v.alleleCount(); //Supporting read count
		vec[idx++]=(float)v.coverage(); //Total coverage depth
		vec[idx++]=(float)v.alleleFraction(); //Allele frequency

		//Quality metrics (raw values)
		vec[idx++]=(float)v.baseQMax; //Maximum base quality
		vec[idx++]=(float)v.mapQMax; //Maximum mapping quality
		vec[idx++]=(float)v.idMax*0.001f; //Maximum identity converted from per-mille to fraction
		vec[idx++]=(float)v.endDistMax; //Maximum end distance

		//Quality metrics (averages)
		final int count=v.alleleCount();
		vec[idx++]=(count>0 ? (float)(v.baseQSum/(double)count) : 0f); //Average base quality
		vec[idx++]=(count>0 ? (float)(v.mapQSum/(double)count) : 0f); //Average mapping quality
		vec[idx++]=(count>0 ? (float)(v.idSum/(double)count*0.001) : 0f); //Mean identity converted from per-mille to fraction
		vec[idx++]=(count>0 ? (float)(v.endDistSum/(double)count) : 0f); //Average end distance

		//Strand bias and pairing
		vec[idx++]=(float)v.strandRatio(); //Strand balance ratio
		vec[idx++]=(count>0 ? (float)(v.properPairCount/(double)count) : 0f); //Proper pair rate

		//Context features relative to dataset averages
		//TODO: Possible bug [var2/FeatureVectorMaker#002] (LOW, latent - this method is DEAD per #001) - these
		//three ratio dims divide by count with NO count>0 guard (the avg dims above DO guard), so count==0
		//gives baseQSum/count = 0.0/0.0 = NaN into the vector, and there is no validateVector here (unlike
		//VectorUMP45). If toVector is ever revived, add the count>0 guard (NaN->1.0).
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
	 * Builds eight historical features and scales with literal min/max endpoints.
	 * The source originally attributed these bounds to training data. Non-NaN
	 * scaled values are clipped to the inclusive interval [0,1]; NaN is retained.
	 * Revised AF may update the variant's insertion cache. Mean accessors have no
	 * zero-count guard. Scaffold context follows Var.nScan when scanning is off,
	 * v.start without a map, or v.contigEndDist with a map.
	 * @param v Non-null variant with supporting-read statistics
	 * @param pairingRate Unused compatibility parameter
	 * @param totalQualityAvg Unused compatibility parameter
	 * @param totalMapqAvg Unused compatibility parameter
	 * @param readLengthAvg Mean read length in bases for revised allele fraction
	 * @param ploidy Unused compatibility parameter
	 * @param map Optional scaffold map for strand-bias context
	 * @return Fresh eight-element scaled vector
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
				//Clip to the inclusive interval [0,1]; NaN remains NaN
				//This handles cases where a new variant's feature is outside the training range
				scaledVec[i]=Math.max(0.0f, Math.min(1.0f, scaledValue));
			}
		}
		return scaledVec;
	}

	/**
	 * Retains the 16-slot DONOVAN placeholder: count, coverage, AF and type in
	 * slots 0-3, with the remaining slots zero. This does not call VectorDonovan
	 * or load its quantile resources.
	 * @param v Non-null variant
	 * @param pairingRate Unused compatibility parameter
	 * @param totalQualityAvg Unused compatibility parameter
	 * @param totalMapqAvg Unused compatibility parameter
	 * @param readLengthAvg Unused compatibility parameter
	 * @param ploidy Unused compatibility parameter
	 * @param map Unused compatibility parameter
	 * @return Fresh 16-element placeholder vector
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
	 * Sets the global mode for subsequent dispatcher calls, not direct scoring.
	 * Validation is assertion-only; with assertions disabled an invalid value is stored.
	 * @param mode ELBA=0, LAWRENCE=1, DONOVAN=2 or UMP45=3
	 */
	public static void setMode(int mode){
		assert(mode>=0 && mode<=3) : "Invalid mode: "+mode;
		MODE=mode;
	}

	/**
	 * Sets the global mode by name after uppercasing with the default locale.
	 * @param modeStr Non-null ELBA, LAWRENCE, DONOVAN or UMP45 name
	 * @throws RuntimeException if the uppercased name is unknown
	 */
	public static void setMode(String modeStr){
		int mode=Tools.find(modeStr.toUpperCase(), MODE_NAMES);
		if(mode<0){throw new RuntimeException("Unknown feature mode: "+modeStr);}
		MODE=mode;
	}

	/** Returns the legacy global selection. */
	public static int getMode(){return MODE;}

	/**
	 * Returns the selected layout length without checking any network.
	 * @return ELBA=32, LAWRENCE=8, DONOVAN=16 or VectorUMP45.DIMS (currently 33)
	 * @throws RuntimeException if MODE contains an unsupported value
	 */
	public static int getVectorLength(){
		if(MODE==ELBA){return 32;}
		else if(MODE==LAWRENCE){return 8;}
		else if(MODE==DONOVAN){return 16;}
		else if(MODE==UMP45){return VectorUMP45.DIMS;}
		else{throw new RuntimeException("Unknown mode: "+MODE);}
	}

	/** Returns the name indexed by MODE; an invalid stored mode throws. */
	public static String getModeName(){return MODE_NAMES[MODE];}
}
