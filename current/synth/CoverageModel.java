package synth;

import shared.Random;
import shared.Shared;

/** Fixed spatial read-acceptance model used by RandomReadsMG.
 * Combines transformed sine waves with a linear origin-of-replication bias,
 * then maps their sum into the requested probability range. A probability is an
 * acceptance weight, not a guarantee of realized read coverage. Contigs shorter
 * than 2000 bases bypass the model and always have probability 1.
 * Construction consumes the supplied RNG but does not retain it. Probability
 * lookup reads fixed parameters; read decisions consume a caller-supplied RNG,
 * whose sharing and synchronization remain the caller's responsibility.
 * Numeric ranges are not validated here; zero normalization and other invalid
 * floating-point configurations can produce NaN rather than a bounded result.
 * @author Brian Bushnell, Jean
 */
class CoverageModel{

	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses four waves, amplitude 0.7, origin bias up to 0.25 and probability floor 0.1.
	 * @param randy Non-null parameter generator, consumed but not retained
	 */
	public CoverageModel(Random randy){
		this(4, 0.7, 0.25, 0.1, randy);
	}

	/** Uses periods spanning 2000 through 80000 bases.
	 * @param numWaves Number of component waves
	 * @param maxAmplitude_ Total wave-amplitude budget, narrowed to float
	 * @param maxOriBias Upper bound on the randomly chosen origin-bias strength
	 * @param minProb Probability floor for modeled contigs, normally in [0,1]
	 * @param randy Non-null parameter generator, consumed but not retained
	 */
	public CoverageModel(int numWaves, double maxAmplitude_, double maxOriBias, 
			double minProb, Random randy){
		this(numWaves, maxAmplitude_, maxOriBias, minProb, randy, 2000, 80000);
	}

	/** Samples fixed waves from equal-width bands of the linear period interval.
	 * Longer-period bands receive larger relative weights; positive total weight
	 * is normalized to the amplitude budget. A sole wave receives the whole weight.
	 * Draw order is one origin-bias float, then one period float and one phase
	 * double per wave, including when a parameter range has zero width.
	 * Callers supply usable numeric values: periods must be positive and ordered,
	 * and the resulting floating-point normalization must be positive and finite.
	 * These expectations are not checked or converted into fallback behavior.
	 * @param numWaves Number of waves; ordinary wave modeling uses at least one
	 * @param maxAmplitude_ Total wave-amplitude budget, narrowed to float
	 * @param maxOriBias Upper bound on origin bias; zero disables this contribution
	 * @param minProb Probability floor, normally in [0,1]
	 * @param randy Non-null parameter generator, consumed but not retained
	 * @param minPeriod Lower period bound in bases
	 * @param maxPeriod Upper period bound in bases; may equal minPeriod
	 */
	public CoverageModel(int numWaves, double maxAmplitude_, double maxOriBias, 
			double minProb, Random randy, int minPeriod, int maxPeriod){
		amplitudes=new float[numWaves];
		inversePeriods=new float[numWaves];
		offsets=new float[numWaves];
		maxOriBiasStrength=(float)(randy.nextFloat()*maxOriBias);
		minProbability=(float)minProb;
		float maxAmplitude=(float)maxAmplitude_;

		// Define period ranges in LINEAR space
		float periodRange=maxPeriod-minPeriod;
		float rangePerWave=periodRange/numWaves;

		// Generate parameters for each sine wave
		float totalAmplitude=0;
		for(int i=0; i<numWaves; i++){
			// Assign amplitude, with larger amplitudes for longer periods
			// A sole wave receives the full weight; interpolation would divide by zero.
			float relativeSize=numWaves==1 ? 1f : 0.4f+0.6f*i/(numWaves-1); // 0.4-1.0 range
			amplitudes[i]=maxAmplitude*relativeSize/numWaves;
			totalAmplitude+=amplitudes[i];

			// Assign each wave to a specific range of the spectrum (LINEAR)
			float rangeMin=minPeriod+i*rangePerWave;
			float rangeMax=rangeMin+rangePerWave;
			float period=rangeMin+randy.nextFloat()*(rangeMax-rangeMin);

			// Store inverse period with 2*PI baked in for faster calculation
			inversePeriods[i]=(float)(pi2/period);

			// Random phase
			offsets[i]=(float)(randy.nextDouble()*pi2);
		}

		// Normalize amplitudes to fit within maxAmplitude
		if(totalAmplitude>0){
			for(int i=0; i<numWaves; i++){
				amplitudes[i]=amplitudes[i]*maxAmplitude/totalAmplitude;
			}
		}

		// Max possible value for normalization
		maxPossibleSum=maxAmplitude+maxOriBiasStrength;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Computes the acceptance weight without drawing random numbers or changing state.
	 * Each wave contributes one minus the eighth power of its sine mapped into [0,1].
	 * Origin bias decreases linearly with position; the combined contribution is
	 * normalized and scaled above the floor. Finite results are clamped to the
	 * requested range, but the clamp does not repair NaN. Contigs shorter than 2000
	 * return 1 before any wave, position or normalization arithmetic.
	 * @param position Zero-based position within the contig
	 * @param contigLength Contig length in bases
	 * @return Acceptance probability, or NaN for an invalid floating-point model
	 */
	public float getProbability(int position, int contigLength){
		// For very short contigs, use uniform probability
		if(contigLength<2000){return 1;}

		// Normalized position for ORI bias
		float relPos=position/(float)contigLength;

		// Sum transformed sine wave contributions
		float sinSum=0;
		for(int i=0; i<amplitudes.length; i++){
			// Optimized angle calculation using pre-computed inverse period
			float angle=(position*inversePeriods[i])+offsets[i];
			float normalizedSine=(float)((Math.sin(angle)+1f)*0.5f); // 0-1 range

			// Efficient power calculation for sharper valleys
			float squared=normalizedSine*normalizedSine;
			float fourth=squared*squared;
			float eighth=fourth*fourth;

			sinSum+=amplitudes[i]*(1f-eighth);
		}

		// Add ORI bias (linear decrease from start to end)
		float maxOriBiasValue=maxOriBiasStrength*(1f-relPos);

		// Scale to available probability range (above minimum)
		float availableRange=1f-minProbability;
		//TODO: Probable bug - SYN009: waveamp=0 and oribias=0 give 0/0 here, so
		//long-contig probabilities become NaN. RandomReadsMG accepts both options.
		float scaledProbability=minProbability+
				availableRange*(sinSum+maxOriBiasValue)/maxPossibleSum;

		// Ensure bounds
		return Math.max(minProbability, Math.min(1f, scaledProbability));
	}

	/** Tests one double draw against the probability using an inclusive comparison.
	 * Consumes one draw even for a short contig or a NaN probability; NaN rejects.
	 * @param position Zero-based candidate position within the contig
	 * @param contigLength Contig length in bases
	 * @param randy Non-null decision generator owned by the caller
	 * @return True when the draw is no greater than the modeled probability
	 */
	public boolean shouldGenerateReadAt(int position, int contigLength, Random randy){
		float prob=getProbability(position, contigLength);
		return randy.nextDouble()<=prob;
	}

	/** Diagnostic that convolves probabilities with a fixed read length and rescales
	 * the resulting coverage track to the requested mean. Samples model parameters,
	 * not individual reads. Prints 100-base bin means to stdout and parameters to
	 * stderr; the last bin uses its actual base count. Requires positive genome
	 * and read lengths and a usable model; argument values are otherwise unchecked.
	 * @param args Positional numWaves, coverage, readLength and genomeLength values
	 */
	public static void main(String[] args){
		if(args.length<4){
			System.err.println("Usage: java CoverageModel <numWaves> <coverage> <readLength> <genomeLength>");
			System.exit(1);
		}

		// Parse arguments
		int numWaves=Integer.parseInt(args[0]);
		float targetCoverage=Float.parseFloat(args[1]);
		int readLength=Integer.parseInt(args[2]);
		int genomeLength=Integer.parseInt(args[3]);

		// Constants
		float maxAmplitude=0.8f; // Increase for more dramatic effect
		float maxOriBias=0.25f;
		float minProb=0.05f; // Lower minimum for better contrast
		int binSize=100;

		// Initialize the model
		Random randy=shared.Shared.random();
		CoverageModel model=new CoverageModel(numWaves, maxAmplitude, maxOriBias, minProb, randy);

		// Create a perfect coverage track based on probabilities
		float[] rawProbabilities=new float[genomeLength];
		for(int i=0; i<genomeLength; i++){
			rawProbabilities[i]=model.getProbability(i, genomeLength);
		}

		// Convolve the probabilities with read length to get actual coverage
		float[] coverage=new float[genomeLength];
		for(int i=0; i<genomeLength; i++){
			float readProb=rawProbabilities[i];
			for(int j=0; j<readLength && i+j<genomeLength; j++){
				coverage[i+j]+=readProb;
			}
		}

		// Scale the coverage to match target
		float avgCov=0;
		for(int i=0; i<genomeLength; i++){
			avgCov+=coverage[i];
		}
		avgCov/=genomeLength;

		float scaleFactor=targetCoverage/avgCov;
		for(int i=0; i<genomeLength; i++){
			coverage[i]*=scaleFactor;
		}

		// Output binned results
		System.out.println("#pos\tavgCov");
		for(int i=0; i<genomeLength; i+=binSize){
			float binTotal=0;
			int validBases=0;
			for(int j=0; j<binSize && i+j<genomeLength; j++){
				binTotal+=coverage[i+j];
				validBases++;
			}
			System.out.printf("%d\t%.2f%n", i, binTotal/validBases);
		}

		// Report model parameters
		System.err.println("Model parameters:");
		for(int i=0; i<numWaves; i++){
			float period=pi2/model.inversePeriods[i];
			System.err.printf("Wave %d: period=%.1f, amplitude=%.4f, phase=%.4f%n", 
					i, period, model.amplitudes[i], model.offsets[i]);
		}
		System.err.println("ORI bias: "+model.maxOriBiasStrength);
		System.err.println("Min probability: "+model.minProbability);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Main              ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Legacy sampling diagnostic, not invoked by main. Candidate positions advance
	 * by 31 modulo genome length; acceptance uses the model and a random draw.
	 * Stops after a target read count truncated from coverage*genomeLength/readLength,
	 * then prints 100-base coverage bins and generation statistics. This can retry
	 * indefinitely if no visited position accepts. Terminal bucket accounting has
	 * the unresolved SYN010 issue annotated below; this is not an exact coverage oracle.
	 * @param args Positional numWaves, coverage, readLength and genomeLength values
	 */
	public static void main_old(String[] args){
		if(args.length<4){
			System.err.println("Usage: java CoverageModel <numWaves> <coverage> <readLength> <genomeLength>");
			System.exit(1);
		}

		// Parse arguments
		int numWaves=Integer.parseInt(args[0]);
		float targetCoverage=Float.parseFloat(args[1]);
		int readLength=Integer.parseInt(args[2]);
		int genomeLength=Integer.parseInt(args[3]);

		// Constants
		float maxAmplitude=0.6f;
		float maxOriBias=0.3f;
		float minProb=0.2f;
		int binSize=100;

		// Initialize the model
		Random randy=shared.Shared.random();
		CoverageModel model=new CoverageModel(numWaves, maxAmplitude, maxOriBias, minProb, randy);

		// Calculate how many reads to simulate
		long totalReads=(long)(targetCoverage*genomeLength/readLength);

		// Initialize coverage bins
		int numBins=(genomeLength+binSize-1)/binSize; // Ceiling division
		int[] basesPerBucket=new int[numBins]; // Count total bases in each bucket
		int[] bucketSizes=new int[numBins]; // Actual size of each bucket

		// Set up bucket sizes
		for(int i=0; i<numBins-1; i++){
			bucketSizes[i]=binSize;
		}
		// Last bucket might be smaller
		bucketSizes[numBins-1]=genomeLength-(numBins-1)*binSize;

		// Simulate reads
		long readsAttempted=0;
		long readsGenerated=0;
		long ptr=0;
		while(readsGenerated<totalReads){
			int pos;
			do{
				// Advance the candidate position; acceptance, not position, is random.
				//        		pos=randy.nextInt(genomeLength);
				pos=(int)(ptr%genomeLength);
				ptr+=31;//A prime
				readsAttempted++;
			}while(!model.shouldGenerateReadAt(pos, genomeLength, randy));

			readsGenerated++;

			// Add coverage for this read
			int readEnd=Math.min(pos+readLength, genomeLength);
			int startBucket=pos/binSize;
			int endBucket=(readEnd-1)/binSize;

			// First bucket - partial
			//TODO: Probable bug - SYN010: this ignores readEnd when the terminal bucket
			//is partial. genomeLength=150, pos=149, readLength=10 counts 10 instead of 1.
			int basesInFirstBucket=Math.min((startBucket+1)*binSize-pos, readLength);
			basesPerBucket[startBucket]+=basesInFirstBucket;

			// Middle buckets - full
			for(int bucket=startBucket+1; bucket<endBucket; bucket++){
				basesPerBucket[bucket]+=binSize;
			}

			// Last bucket - partial (if not same as first bucket)
			if(startBucket<endBucket){
				int basesInLastBucket=readEnd-endBucket*binSize;
				if(basesInLastBucket>0){
					basesPerBucket[endBucket]+=basesInLastBucket;
				}
			}
		}

		// Output results
		System.out.println("#pos\tavgCov");
		for(int i=0; i<numBins; i++){
			int position=i*binSize;
			float avgCoverage=basesPerBucket[i]/(float)bucketSizes[i];
			System.out.printf("%d\t%.2f%n", position, avgCoverage);
		}

		// Report stats
		//        System.err.println("Period: "+Arrays.toString(model.periods));
		System.err.println("Reads generated: "+readsGenerated);
		System.err.println("Positions attempted: "+readsAttempted);
		System.err.println("Acceptance rate: "+(100.0*readsGenerated/readsAttempted)+"%");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Per-wave weights, normalized when their initial total is positive. */
	private final float[] amplitudes;
	/** Angular increments in radians per base, with 2*PI included. */
	private final float[] inversePeriods;
	/** Sampled phase offsets in radians. */
	private final float[] offsets;
	/** Sampled linear-bias strength at the contig origin. */
	private final float maxOriBiasStrength;
	/** Lower bound requested for finite modeled probabilities. */
	private final float minProbability;
	/** Amplitude budget plus sampled origin strength; probability divisor. */
	private final float maxPossibleSum;
	
	/*--------------------------------------------------------------*/
	
	/** One full turn, rounded to float for angle arithmetic. */
	private static final float pi2=(float)(2*Math.PI);

}
