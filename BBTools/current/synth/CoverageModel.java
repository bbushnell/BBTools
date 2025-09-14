package synth;

import java.util.Random;

/**
 * Models realistic genomic coverage variation using summed sine waves plus ORI bias.
 * Optimized for performance with pre-computed inverse periods and efficient power calculations.
 *
 *@author Brian Bushnell
 *@contributor Isla
 *@date May 18, 2025
 */
class CoverageModel{

	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	public CoverageModel(Random randy){
		this(4, 0.7, 0.25, 0.1, randy);
	}

	public CoverageModel(int numWaves, double maxAmplitude_, double maxOriBias, 
			double minProb, Random randy){
		this(numWaves, maxAmplitude_, maxOriBias, minProb, randy, 2000, 80000);
	}

	/**
	 * Creates a coverage model with multiple sine waves plus ORI bias.
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
			float relativeSize=0.4f+0.6f*i/(numWaves-1); // 0.4-1.0 range
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

	/**
	 * Calculates probability of a read starting at the given position.
	 * Optimized version using pre-computed values and efficient calculations.
	 * Ignores contigs under 2kbp.
	 * Otherwise, contigLength is only used for the ORI linear slope.
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
		float scaledProbability=minProbability+
				availableRange*(sinSum+maxOriBiasValue)/maxPossibleSum;

		// Ensure bounds
		return Math.max(minProbability, Math.min(1f, scaledProbability));
	}

	/**
	 * Determines whether to generate a read at the given position.
	 */
	public boolean shouldGenerateReadAt(int position, int contigLength, Random randy){
		float prob=getProbability(position, contigLength);
		return randy.nextDouble()<=prob;
	}

	/**
	 * Simulates coverage directly from the probability model without random sampling.
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
		Random randy=new Random();
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
	 * Simulates reads and outputs coverage histogram.
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
		Random randy=new Random();
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
			do {
				// Pick a random position
				//        		pos=randy.nextInt(genomeLength);
				pos=(int)(ptr%genomeLength);
				ptr+=31;//A prime
				readsAttempted++;
			} while(!model.shouldGenerateReadAt(pos, genomeLength, randy));

			readsGenerated++;

			// Add coverage for this read
			int readEnd=Math.min(pos+readLength, genomeLength);
			int startBucket=pos/binSize;
			int endBucket=(readEnd-1)/binSize;

			// First bucket - partial
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
	
	private final float[] amplitudes;
	private final float[] inversePeriods; // Store 2*PI/period for performance
	private final float[] offsets;
	private final float maxOriBiasStrength;
	private final float minProbability;
	private final float maxPossibleSum;
	
	/*--------------------------------------------------------------*/
	
	private static final float pi2=(float)(2*Math.PI);

}
