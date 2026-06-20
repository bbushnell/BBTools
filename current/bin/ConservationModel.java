package bin;

import map.IntHashSet;
import shared.Random;

/**
 * Models realistic genomic conservation variation using summed sine waves.
 * Used to determine position-dependent mutation rates during sequence simulation.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 30, 2025
 */
public class ConservationModel {

    /*--------------------------------------------------------------*/
    /*----------------             Init             ----------------*/
    /*--------------------------------------------------------------*/

    /**
     * Creates a conservation model with default parameters.
     * Uses zero base mutation rate and 3 sine waves.
     * @param randy Random number generator for wave parameters
     */
    public ConservationModel(Random randy) {
        this(0.0f, 3, randy);
    }

    /**
     * Creates a conservation model with specified base rate and wave count.
     * Uses default amplitude (0.3) and period range (200-500 bp).
     * @param baseMutationRate Base mutation probability (0-1)
     * @param numWaves Number of sine waves to combine
     * @param randy Random number generator for wave parameters
     */
    public ConservationModel(float baseMutationRate, int numWaves, Random randy) {
        this(baseMutationRate, numWaves, 0.3f, randy, 200, 500);
    }

    /**
     * Creates a conservation model with multiple sine waves for mutation rate variation.
     * @param baseMutationRate Base mutation probability (0-1)
     * @param numWaves Number of sine waves to combine
     * @param maxAmplitude Maximum amplitude for sine wave oscillation
     * @param randy Random number generator
     * @param minPeriod Minimum period length in base pairs
     * @param maxPeriod Maximum period length in base pairs
     */
    public ConservationModel(float baseMutationRate, int numWaves, float maxAmplitude, 
                           Random randy, int minPeriod, int maxPeriod) {
        //TODO: Possible bug [bin/ConservationModel#001] - FIX APPLIED (guard): without this, numWaves > unique periods available in [minPeriod,maxPeriod] makes the unique-period do-while loop spin forever (hang). Reachable via uncapped waves=/sinewaves= user args (MutateGenome:108, TestAlignerSuite:61). Guard converts the hang into a loud AssertionError.
        assert(numWaves <= maxPeriod - minPeriod + 1) : "numWaves=" + numWaves +
            " exceeds the " + (maxPeriod - minPeriod + 1) + " unique periods available in [" +
            minPeriod + ", " + maxPeriod + "]; the unique-period loop would never terminate.";
        this.baseMutationRate = baseMutationRate;
        amplitudes = new float[numWaves];
        inversePeriods = new float[numWaves];
        offsets = new float[numWaves];
        
        // Use IntSet to ensure unique periods
        IntHashSet usedPeriods = new IntHashSet();
        
        // Generate parameters for each sine wave
        for(int i = 0; i < numWaves; i++) {
            // Assign equal amplitudes that sum to maxAmplitude
            amplitudes[i] = maxAmplitude / numWaves;
            
            // Pick unique period in range
            //claim: terminates - the #001 guard ensures numWaves <= available unique periods; period uniform in [minPeriod, maxPeriod] inclusive. Worst case (numWaves==available) is coupon-collector slow but finite.
            int period;
            do {
                period = minPeriod + randy.nextInt(maxPeriod - minPeriod + 1);
            } while(usedPeriods.contains(period));
            usedPeriods.add(period);
            
            // Store inverse period with 2*PI baked in for efficiency
            inversePeriods[i] = (float)(pi2 / period);
            
            // Random phase offset
            offsets[i] = randy.nextFloat() * pi2;
        }
    }

    /*--------------------------------------------------------------*/
    /*----------------           Methods            ----------------*/
    /*--------------------------------------------------------------*/

    /**
     * Calculates mutation probability at the given position.
     * @param position Position in sequence (0-based)
     * @return Mutation probability (0-1)
     */
    public float getMutationProbability(int position) {
        //claim: returns mutation prob clamped to [0,1]. Each wave adds amplitudes[i]*((sin+1)/2) in [0,amplitudes[i]]; amplitudes sum to maxAmplitude, so rate in [base, base+maxAmplitude] pre-clamp.
        float sineSum = 0;
        
        // Sum all sine wave contributions
        for(int i = 0; i < amplitudes.length; i++) {
            float angle = (position * inversePeriods[i]) + offsets[i];
            float normalizedSine = (float)((Math.sin(angle) + 1f) * 0.5f); // 0-1 range
            sineSum += amplitudes[i] * normalizedSine;
        }
        
        // Add base mutation rate
        float totalRate = baseMutationRate + sineSum;
        
        // Ensure bounds
        return Math.max(0f, Math.min(1f, totalRate));
    }

    /**
     * Determines whether to mutate at the given position.
     * @param position Position in sequence (0-based)  
     * @param randy Random number generator
     * @return true if position should be mutated
     */
    public boolean shouldMutatePosition(int position, Random randy) {
        //obvious: Bernoulli draw; nextFloat() in [0,1); '<=' means rate>=1 always mutates, rate 0 effectively never (only if nextFloat()==0).
        float mutationRate = getMutationProbability(position);
        return randy.nextFloat() <= mutationRate;
    }

    /*--------------------------------------------------------------*/
    /*----------------           Fields             ----------------*/
    /*--------------------------------------------------------------*/
    
    /** Base mutation probability applied uniformly across sequence */
    private final float baseMutationRate;
    /** Amplitude values for each sine wave component */
    private final float[] amplitudes;
    /** Pre-computed 2*PI/period values for efficient sine calculation */
    private final float[] inversePeriods; // Store 2*PI/period for performance
    /** Random phase offsets for each sine wave component */
    private final float[] offsets;
    
    /*--------------------------------------------------------------*/
    
    /** Pre-computed 2*PI constant for sine wave calculations */
    private static final float pi2 = (float)(2 * Math.PI);
    
}