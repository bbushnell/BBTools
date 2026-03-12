package cardinality;

import java.util.BitSet;
import rand.FastRandomXoshiro;
import parse.Parse;

/**
 * Low-complexity cardinality calibration driver.
 * <p>
 * Tests estimator accuracy on datasets with bounded cardinality and repeated values.
 * Each estimator draws with replacement from a fixed array of unique values, biased
 * toward lower indices via min(rand(), rand()) to simulate skewed frequency distributions.
 * <p>
 * Estimates are recorded after EVERY add (not just on new-value events), indexed by
 * the current true cardinality. This captures how the estimator behaves while "parked"
 * at a given cardinality as more duplicate values arrive and internal state evolves.
 * <p>
 * Two run modes:
 * <ul>
 *   <li>iterations=0 (default): each estimator stops upon saturation (all unique values seen).</li>
 *   <li>iterations=N (float): each estimator runs for cardinality*N total adds regardless,
 *       continuing to accumulate at trueCard=cardinality after saturation.</li>
 * </ul>
 * <p>
 * Multithreaded; each thread processes a disjoint range of estimators sequentially
 * (for(estimator){for(add)}) for cache efficiency, then merges into shared accumulators
 * under a coarse lock.
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public class LowComplexityCalibrationDriver {

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		int cardinality=5000;
		int numDDLs=128;
		float iterations=0;     // 0 = stop on saturation; >0 = run for cardinality*iterations adds
		int buckets=2048;
		int k=31;
		int threads=4;
		long masterSeed=1;
		String loglogtype="dll4";

		for(String arg : args){
			final String[] split=arg.split("=");
			if(split.length!=2){continue;}
			final String a=split[0].toLowerCase();
			final String b=split[1];
			if(a.equals("cardinality") || a.equals("card")){cardinality=Parse.parseIntKMG(b);}
			else if(a.equals("iterations") || a.equals("iter")){iterations=Float.parseFloat(b);}
			else if(a.equals("ddls") || a.equals("dlls")){numDDLs=Parse.parseIntKMG(b);}
			else if(a.equals("buckets")){buckets=Parse.parseIntKMG(b);}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("t") || a.equals("threads")){threads=Integer.parseInt(b);}
			else if(a.equals("seed")){masterSeed=Long.parseLong(b);}
			else if(a.equals("loglogtype") || a.equals("type")){loglogtype=b.toLowerCase();}
			else if(a.equals("cf") || a.equals("loglogcf")){CorrectionFactor.USE_CORRECTION=Parse.parseBoolean(b);}
			else if(a.equals("cardcf")){CardinalityTracker.USE_CARD_CF=Parse.parseBoolean(b);}
			else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}

		// totalAdds: 0 means "run until saturated" (handled via break in loop)
		final int totalAdds=(iterations>0 ? (int)(cardinality*iterations) : Integer.MAX_VALUE);
		final boolean stopOnSaturation=(iterations<=0);

		// Generate fixed value array: cardinality unique random longs
		final FastRandomXoshiro masterRng=new FastRandomXoshiro(masterSeed);
		final long[] valueArray=new long[cardinality];
		for(int i=0; i<cardinality; i++){valueArray[i]=masterRng.nextLong();}

		// Generate per-estimator seeds
		final long[] seeds=new long[numDDLs];
		for(int i=0; i<numDDLs; i++){seeds[i]=masterRng.nextLong();}

		// Determine number of estimate types from a dummy instance
		final CardinalityTracker dummy=DDLCalibrationDriver.makeInstance(loglogtype, buckets, k, 1, 0);
		final int numTypes=dummy.rawEstimates().length;

		// Shared accumulator arrays (indexed by true cardinality 0..cardinality)
		final double[][] sums  =new double[cardinality+1][numTypes];
		final double[][] sumAbs=new double[cardinality+1][numTypes];
		final double[][] sumSq =new double[cardinality+1][numTypes];
		final long[]     counts=new long[cardinality+1];
		final Object     lock  =new Object();

		// Finals for lambda capture
		final int finalCardinality=cardinality;
		final int finalTotalAdds=totalAdds;
		final boolean finalStop=stopOnSaturation;
		final String finalType=loglogtype;
		final int finalBuckets=buckets;
		final int finalK=k;
		final int finalNumTypes=numTypes;
		final int estPerThread=(numDDLs+threads-1)/threads;

		final Thread[] threadArray=new Thread[threads];
		for(int tid=0; tid<threads; tid++){
			final int estStart=tid*estPerThread;
			final int estEnd=Math.min(estStart+estPerThread, numDDLs);

			// Per-thread local accumulators (reused across all estimators in this thread)
			final double[][] lSums  =new double[finalCardinality+1][finalNumTypes];
			final double[][] lSumAbs=new double[finalCardinality+1][finalNumTypes];
			final double[][] lSumSq =new double[finalCardinality+1][finalNumTypes];
			final long[]     lCounts=new long[finalCardinality+1];

			threadArray[tid]=new Thread(()->{
				final FastRandomXoshiro rng=new FastRandomXoshiro(seeds[estStart]);
				final BitSet seen=new BitSet(finalCardinality);

				for(int estIdx=estStart; estIdx<estEnd; estIdx++){
					// Reset state for this estimator
					rng.setSeed(seeds[estIdx]);
					seen.clear();
					int trueCard=0;
					final CardinalityTracker est=DDLCalibrationDriver.makeInstance(
						finalType, finalBuckets, finalK, seeds[estIdx], 0);

					for(int add=0; add<finalTotalAdds; add++){
						// Biased draw: min(rand, rand) favors lower indices
						final int pos=Math.min(rng.nextInt(finalCardinality), rng.nextInt(finalCardinality));
						est.hashAndStore(valueArray[pos]);
						if(seen!=null && !seen.get(pos)){
							seen.set(pos);
							trueCard++;
						}
						// Record after every add, at current true cardinality
						if(trueCard>0){
							final double[] estimates=est.rawEstimates();
							for(int t=0; t<finalNumTypes; t++){
								final double relErr=(estimates[t]-trueCard)/(double)trueCard;
								lSums[trueCard][t]  +=relErr;
								lSumAbs[trueCard][t]+=Math.abs(relErr);
								lSumSq[trueCard][t] +=relErr*relErr;
							}
							lCounts[trueCard]++;
						}
						if(finalStop && trueCard==finalCardinality){break;}
					}
				}

				// Merge local accumulators into shared arrays under lock
				synchronized(lock){
					for(int c=1; c<=finalCardinality; c++){
						if(lCounts[c]==0){continue;}
						counts[c]+=lCounts[c];
						for(int t=0; t<finalNumTypes; t++){
							sums[c][t]  +=lSums[c][t];
							sumAbs[c][t]+=lSumAbs[c][t];
							sumSq[c][t] +=lSumSq[c][t];
						}
					}
				}
			});
			threadArray[tid].start();
		}

		for(Thread t : threadArray){
			try{t.join();}catch(InterruptedException e){throw new RuntimeException(e);}
		}

		// Print header
		final StringBuilder header=new StringBuilder();
		header.append("TrueCard\tOccupancy");
		for(String name : EST_NAMES){
			header.append('\t').append(name).append("_err");
			header.append('\t').append(name).append("_abs");
			header.append('\t').append(name).append("_std");
		}
		System.out.println(header);

		// Print rows sampled at ~1% cardinality intervals
		long prevBucket=-1;
		for(int c=1; c<=cardinality; c++){
			if(counts[c]==0){continue;}
			final long bucket=c*100L/cardinality;
			if(bucket==prevBucket && c!=cardinality){continue;}
			prevBucket=bucket;
			final double n=counts[c];
			final double occ=(double)c/buckets;
			final StringBuilder sb=new StringBuilder();
			sb.append(c).append('\t').append(String.format("%.6f", occ));
			for(int t=0; t<numTypes && t<EST_NAMES.length; t++){
				final double meanErr =sums[c][t]/n;
				final double meanAbs =sumAbs[c][t]/n;
				final double variance=sumSq[c][t]/n-meanErr*meanErr;
				final double std     =Math.sqrt(Math.max(0, variance));
				sb.append('\t').append(String.format("%.6f", meanErr));
				sb.append('\t').append(String.format("%.6f", meanAbs));
				sb.append('\t').append(String.format("%.6f", std));
			}
			System.out.println(sb);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Estimate type names matching CardinalityStats.toArray() output order. */
	static final String[] EST_NAMES;
	static{
		final String[] base={"Mean","HMean","HMeanM","GMean","HLL","LC","Hybrid","MWA","MedianCorr","Mean99","Micro","DLC","DLC3B","DLCBest","HybDLC"};
		EST_NAMES=new String[base.length+CardinalityStats.NUM_DLC_TIERS];
		System.arraycopy(base, 0, EST_NAMES, 0, base.length);
		for(int i=0; i<CardinalityStats.NUM_DLC_TIERS; i++){EST_NAMES[base.length+i]="DLC"+i;}
	}

}
