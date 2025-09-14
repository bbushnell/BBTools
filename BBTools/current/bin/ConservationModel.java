package bin;

import java.util.Random;

import structures.IntHashSet;

/**
 * Models realistic genomic conservation variation using summed sine waves.
 * Used to determine position-dependent mutation rates during sequence simulation.
 *
 * @author Brian Bushnell
 * @contributor Isla Winglet
 * @date May 30, 2025
 */
public class ConservationModel {

    /*--------------------------------------------------------------*/
    /*----------------             Init             ----------------*/
    /*--------------------------------------------------------------*/

    public ConservationModel(Random randy) {
        this(0.0f, 3, randy);
    }

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
        float mutationRate = getMutationProbability(position);
        return randy.nextFloat() <= mutationRate;
    }

    /*--------------------------------------------------------------*/
    /*----------------           Fields             ----------------*/
    /*--------------------------------------------------------------*/
    
    private final float baseMutationRate;
    private final float[] amplitudes;
    private final float[] inversePeriods; // Store 2*PI/period for performance
    private final float[] offsets;
    
    /*--------------------------------------------------------------*/
    
    private static final float pi2 = (float)(2 * Math.PI);
    
}