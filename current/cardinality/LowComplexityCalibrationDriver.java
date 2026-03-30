package cardinality;

import java.util.ArrayList;
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
 * Estimates are recorded at exponentially-spaced reporting thresholds (reportfrac intervals,
 * matching the high-complexity driver's spacing). rawEstimates() is called on every add
 * while trueCard is parked at a threshold, capturing tier-promotion effects from duplicates.
 * Between thresholds, only add() runs — no rawEstimates() overhead.
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
		double reportFrac=0.01;
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
			else if(a.equals("reportfrac")){reportFrac=Double.parseDouble(b);}
			else if(a.equals("loglogtype") || a.equals("type")){loglogtype=b.toLowerCase();}
			else if(a.equals("cf") || a.equals("loglogcf")){CorrectionFactor.USE_CORRECTION=Parse.parseBoolean(b);}
			else if(a.equals("cardcf")){CardinalityTracker.USE_CARD_CF=Parse.parseBoolean(b);}
			else if(a.equals("frozenhistory") || a.equals("frozen")){UltraLogLog8.FROZEN_HISTORY=Parse.parseBoolean(b);}
			else if(a.equals("printcv") || a.equals("cv")){DDLCalibrationDriver.PRINT_CV=Parse.parseBoolean(b);}
			else if(a.equals("correctoverflow") || a.equals("co")){
				DynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
				BankedDynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
			}else if(a.equals("earlypromote") || a.equals("ep")){
				DynamicLogLog3.EARLY_PROMOTE=Parse.parseBoolean(b);
				DynamicLogLog4.EARLY_PROMOTE=Parse.parseBoolean(b);
				BankedDynamicLogLog3.EARLY_PROMOTE=Parse.parseBoolean(b);
			}
			else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}

		// totalAdds: 0 means "run until saturated" (handled via break in loop)
		final long totalAdds=(iterations>0 ? (long)(cardinality*iterations) : Long.MAX_VALUE);
		final boolean stopOnSaturation=(iterations<=0);

		// Pre-compute exponentially-spaced reporting thresholds (same as DDLCalibrationDriver)
		final int[] thresholds=computeThresholds(cardinality, reportFrac);
		final int numThresholds=thresholds.length;
		System.err.println("Thresholds: "+numThresholds+" points from 1 to "+cardinality);

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

		// Shared accumulator arrays indexed by threshold index (tiny: ~numThresholds rows)
		final double[][] sums  =new double[numThresholds][numTypes];
		final double[][] sumAbs=new double[numThresholds][numTypes];
		final double[][] sumSq =new double[numThresholds][numTypes];
		final long[]     counts=new long[numThresholds];
		final Object     lock  =new Object();

		// Finals for lambda capture
		final int finalCardinality=cardinality;
		final long finalTotalAdds=totalAdds;
		final boolean finalStop=stopOnSaturation;
		final String finalType=loglogtype;
		final int finalBuckets=buckets;
		final int finalK=k;
		final int finalNumTypes=numTypes;
		final int finalNumThresholds=numThresholds;
		final int estPerThread=(numDDLs+threads-1)/threads;

		final long t0=System.nanoTime();
		final Thread[] threadArray=new Thread[threads];
		for(int tid=0; tid<threads; tid++){
			final int estStart=tid*estPerThread;
			final int estEnd=Math.min(estStart+estPerThread, numDDLs);

			// Per-thread local accumulators (indexed by threshold — tiny footprint)
			final double[][] lSums  =new double[finalNumThresholds][finalNumTypes];
			final double[][] lSumAbs=new double[finalNumThresholds][finalNumTypes];
			final double[][] lSumSq =new double[finalNumThresholds][finalNumTypes];
			final long[]     lCounts=new long[finalNumThresholds];

			threadArray[tid]=new Thread(()->{
				final FastRandomXoshiro rng=new FastRandomXoshiro(seeds[estStart]);
				final BitSet seen=new BitSet(finalCardinality);

				for(int estIdx=estStart; estIdx<estEnd; estIdx++){
					// Reset state for this estimator
					rng.setSeed(seeds[estIdx]);
					seen.clear();
					int trueCard=0;
					int ti=0; // threshold index
					final CardinalityTracker est=DDLCalibrationDriver.makeInstance(
						finalType, finalBuckets, finalK, seeds[estIdx], 0);

					for(long add=0; add<finalTotalAdds; add++){
						// Biased draw: min(rand, rand) favors lower indices
						final int pos=Math.min(rng.nextInt(finalCardinality), rng.nextInt(finalCardinality));
						est.add(valueArray[pos]);
						if(!seen.get(pos)){
							seen.set(pos);
							trueCard++;
							// Advance threshold index past any thresholds we've passed
							while(ti<finalNumThresholds && trueCard>thresholds[ti]){ti++;}
						}
						// Record on every add while parked at a threshold cardinality.
						// This captures tier-promotion effects from duplicate adds.
						if(ti<finalNumThresholds && trueCard==thresholds[ti]){
							final double[] estimates=est.rawEstimates();
							for(int t=0; t<finalNumTypes; t++){
								final double relErr=(estimates[t]-trueCard)/(double)trueCard;
								lSums[ti][t]  +=relErr;
								lSumAbs[ti][t]+=Math.abs(relErr);
								lSumSq[ti][t] +=relErr*relErr;
							}
							lCounts[ti]++;
						}
						if(finalStop && trueCard==finalCardinality){break;}
					}
				}

				// Merge local accumulators into shared arrays under lock
				synchronized(lock){
					for(int ti2=0; ti2<finalNumThresholds; ti2++){
						if(lCounts[ti2]==0){continue;}
						counts[ti2]+=lCounts[ti2];
						for(int t=0; t<finalNumTypes; t++){
							sums[ti2][t]  +=lSums[ti2][t];
							sumAbs[ti2][t]+=lSumAbs[ti2][t];
							sumSq[ti2][t] +=lSumSq[ti2][t];
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
			header.append('\t').append(name).append("_cv");
		}
		System.out.println(header);

		// Print one row per threshold
		for(int ti=0; ti<numThresholds; ti++){
			if(counts[ti]==0){continue;}
			final double n=counts[ti];
			final int card=thresholds[ti];
			final double occ=(double)card/buckets;
			final StringBuilder sb=new StringBuilder();
			sb.append(card).append('\t').append(String.format("%.6f", occ));
			for(int t=0; t<numTypes && t<EST_NAMES.length; t++){
				final double meanErr =sums[ti][t]/n;
				final double meanAbs =sumAbs[ti][t]/n;
				final double variance=sumSq[ti][t]/n-meanErr*meanErr;
				final double std     =Math.sqrt(Math.max(0, variance));
				sb.append('\t').append(String.format("%.6f", meanErr));
				sb.append('\t').append(String.format("%.6f", meanAbs));
				sb.append('\t').append(String.format("%.6f", std));
				final double cvDenom=Math.abs(1.0+meanErr);
				sb.append('\t').append(String.format("%.6f", cvDenom>0 ? std/cvDenom : 0));
			}
			System.out.println(sb);
		}

		// Stderr summary: avg and peak abs error per estimator + runtime
		{
			final double elapsed=(System.nanoTime()-t0)*1e-9;
			System.err.println();
			System.err.println("=== Low-Complexity Calibration Summary ===");
			System.err.println("Type: "+loglogtype+"  Buckets: "+buckets+"  DDLs: "+numDDLs
				+"  Card: "+cardinality+"  Iter: "+iterations
				+"  Rows: "+numThresholds+"  Elapsed: "+String.format("%.1f", elapsed)+"s");
			System.err.println("--- Avg and Peak Mean Absolute Error, Avg CV (lower = better) ---");
			final double[] totalAbsErr=new double[numTypes];
			final double[] peakAbsErr=new double[numTypes];
			final double[] totalCV=new double[numTypes];
			int rowsWithData=0;
			for(int ti=0; ti<numThresholds; ti++){
				if(counts[ti]==0){continue;}
				rowsWithData++;
				for(int t=0; t<numTypes; t++){
					final double n2=counts[ti];
					final double meanErr=sums[ti][t]/n2;
					final double meanAbsAtRow=sumAbs[ti][t]/n2;
					totalAbsErr[t]+=meanAbsAtRow;
					if(meanAbsAtRow>peakAbsErr[t]){peakAbsErr[t]=meanAbsAtRow;}
					final double variance=sumSq[ti][t]/n2-meanErr*meanErr;
					final double std=Math.sqrt(Math.max(0, variance));
					final double cvDenom=Math.abs(1.0+meanErr);
					if(cvDenom>0){totalCV[t]+=std/cvDenom;}
				}
			}
			System.err.println(String.format("%-12s %-12s %-12s %s", "", "AvgAbsErr", "PeakAbsErr", "AvgCV"));
			final int[] keyIdx={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
			for(int ki=0; ki<keyIdx.length; ki++){
				final int e=keyIdx[ki];
				if(e>=numTypes || e>=EST_NAMES.length){continue;}
				System.err.println(String.format("%-12s %.8f  %.8f  %.8f",
					EST_NAMES[e], totalAbsErr[e]/rowsWithData, peakAbsErr[e],
					rowsWithData>0 ? totalCV[e]/rowsWithData : 0));
			}
			System.err.println();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	/** Generates exponentially-spaced thresholds from 1 to maxCard at reportFrac intervals. */
	static int[] computeThresholds(final int maxCard, final double reportFrac){
		final ArrayList<Integer> list=new ArrayList<>();
		int next=1;
		while(next<maxCard){
			list.add(next);
			next=Math.max(next+1, (int)(next*(1+reportFrac)));
		}
		list.add(maxCard);
		final int[] arr=new int[list.size()];
		for(int i=0; i<arr.length; i++){arr[i]=list.get(i);}
		return arr;
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
