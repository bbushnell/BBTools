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
 * Output format matches DDLCalibrationDriver2 exactly: same columns, same stderr summary,
 * including LDLC columns for history-capable classes (UDLL6, PLL16c).
 * <p>
 * Multithreaded; each thread processes a disjoint range of estimators sequentially
 * (for(estimator){for(add)}) for cache efficiency, then merges into shared accumulators
 * under a coarse lock.
 *
 * @author Brian Bushnell, Chloe, Eru
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
			else if(a.equals("printcv") || a.equals("cv")){DDLCalibrationDriver.PRINT_CV=Parse.parseBoolean(b);}
			else if(a.equals("printstd") || a.equals("std")){DDLCalibrationDriver.PRINT_STD=Parse.parseBoolean(b);}
			else if(a.equals("correctoverflow") || a.equals("co")){
				DynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
				BankedDynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
			}else if(a.equals("earlypromote") || a.equals("ep")){
				DynamicLogLog3.EARLY_PROMOTE=Parse.parseBoolean(b);
				DynamicLogLog4.EARLY_PROMOTE=Parse.parseBoolean(b);
				BankedDynamicLogLog3.EARLY_PROMOTE=Parse.parseBoolean(b);
			}else if(a.equals("sbsformula") || a.equals("usesbsformula")){
				CorrectionFactor.USE_SBS_FORMULA=Parse.parseBoolean(b);
			}else if(a.equals("formulas") || a.equals("useformulas")){
				CorrectionFactor.USE_FORMULAS=Parse.parseBoolean(b);
			}else if(a.equals("usesbs") || a.equals("sbsinhybrid")){
				CardinalityTracker.USE_SBS_IN_HYBRID=Parse.parseBoolean(b);
			}
			else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}

		// Reuse constants from DDLCalibrationDriver/Driver2
		final int NUM_EST=DDLCalibrationDriver.NUM_EST;
		final int NUM_OUT1=DDLCalibrationDriver.NUM_OUT1;
		final int NUM_LDLC=DDLCalibrationDriver2.NUM_LDLC;
		final String[] LDLC_NAMES=DDLCalibrationDriver2.LDLC_NAMES;

		// Set up per-class CF/formula whitelist (same as Driver2)
		DDLCalibrationDriver.makeInstance(loglogtype, buckets, k, 0L, 0);
		DDLCalibrationDriver.v3ColsForType(loglogtype);

		// totalAdds: 0 means "run until saturated" (handled via break in loop)
		final long totalAdds=(iterations>0 ? (long)(cardinality*iterations) : Long.MAX_VALUE);
		final boolean stopOnSaturation=(iterations<=0);

		// Pre-compute exponentially-spaced reporting thresholds (same as DDLCalibrationDriver)
		final long[] thresholds=DDLCalibrationDriver.computeThresholds(cardinality, reportFrac);
		final int numThresholds=thresholds.length;
		System.err.println("Thresholds: "+numThresholds+" points from 1 to "+cardinality);

		// Generate fixed value array: cardinality unique random longs
		final FastRandomXoshiro masterRng=new FastRandomXoshiro(masterSeed);
		final long[] valueArray=new long[cardinality];
		for(int i=0; i<cardinality; i++){valueArray[i]=masterRng.nextLong();}

		// Generate per-estimator seeds
		final long[] seeds=new long[numDDLs];
		for(int i=0; i<numDDLs; i++){seeds[i]=masterRng.nextLong();}

		// Shared accumulator arrays
		final double[][] sumErr=new double[numThresholds][NUM_OUT1];
		final double[][] sumAbsErr=new double[numThresholds][NUM_OUT1];
		final double[][] sumSqErr=new double[numThresholds][NUM_OUT1];
		final double[][] ldlcSumErr=new double[numThresholds][NUM_LDLC];
		final double[][] ldlcSumAbsErr=new double[numThresholds][NUM_LDLC];
		final double[][] ldlcSumSqErr=new double[numThresholds][NUM_LDLC];
		final double[] occSum=new double[numThresholds];
		final int[] nArr=new int[numThresholds];
		final Object lock=new Object();

		// Finals for lambda capture
		final int finalCardinality=cardinality;
		final long finalTotalAdds=totalAdds;
		final boolean finalStop=stopOnSaturation;
		final String finalType=loglogtype;
		final int finalBuckets=buckets;
		final int finalK=k;
		final int finalNumThresholds=numThresholds;
		final int estPerThread=(numDDLs+threads-1)/threads;
		final int[] ldlcIdx={0, 1, 2, 4, 5, 6, 7}; // Indices into ldlcEstimate() array

		final long t0=System.nanoTime();
		final Thread[] threadArray=new Thread[threads];
		for(int tid=0; tid<threads; tid++){
			final int estStart=tid*estPerThread;
			final int estEnd=Math.min(estStart+estPerThread, numDDLs);

			// Per-thread local accumulators
			final double[][] lSumErr=new double[finalNumThresholds][NUM_OUT1];
			final double[][] lSumAbsErr=new double[finalNumThresholds][NUM_OUT1];
			final double[][] lSumSqErr=new double[finalNumThresholds][NUM_OUT1];
			final double[][] lLdlcErr=new double[finalNumThresholds][NUM_LDLC];
			final double[][] lLdlcAbsErr=new double[finalNumThresholds][NUM_LDLC];
			final double[][] lLdlcSqErr=new double[finalNumThresholds][NUM_LDLC];
			final double[] lOccSum=new double[finalNumThresholds];
			final int[] lN=new int[finalNumThresholds];

			threadArray[tid]=new Thread(()->{
				final FastRandomXoshiro rng=new FastRandomXoshiro(seeds[estStart]);
				final BitSet seen=new BitSet(finalCardinality);

				for(int estIdx=estStart; estIdx<estEnd; estIdx++){
					rng.setSeed(seeds[estIdx]);
					seen.clear();
					int trueCard=0;
					int ti=0;
					final CardinalityTracker est=DDLCalibrationDriver.makeInstance(
						finalType, finalBuckets, finalK, seeds[estIdx], 0);

					for(long add=0; add<finalTotalAdds; add++){
						final int pos=Math.min(rng.nextInt(finalCardinality), rng.nextInt(finalCardinality));
						est.add(valueArray[pos]);
						if(!seen.get(pos)){
							seen.set(pos);
							trueCard++;
							while(ti<finalNumThresholds && trueCard>thresholds[ti]){ti++;}
						}
						// Record on every add while parked at a threshold
						if(ti<finalNumThresholds && trueCard==thresholds[ti]){
							final double occ=est.occupancy();
							final double[] estimates=est.rawEstimates();
							lOccSum[ti]+=occ;
							lN[ti]++;
							for(int e=0; e<Math.min(NUM_EST, estimates.length); e++){
								final double v=estimates[e];
								final double err=(v-trueCard)/(double)trueCard;
								lSumErr[ti][e]+=err;
								lSumAbsErr[ti][e]+=Math.abs(err);
								lSumSqErr[ti][e]+=err*err;
							}
							// LDLC columns (UDLL6 or PLL16c with history)
							final double[] ldlcR;
							if(est instanceof UltraDynamicLogLog6){
								ldlcR=((UltraDynamicLogLog6)est).ldlcEstimate();
							}else if(est instanceof ProtoLogLog16c){
								ldlcR=((ProtoLogLog16c)est).ldlcEstimate();
							}else{ldlcR=null;}
							if(ldlcR!=null){
								for(int e=0; e<NUM_LDLC; e++){
									final double v=ldlcR[ldlcIdx[e]];
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][e]+=lerr;
									lLdlcAbsErr[ti][e]+=Math.abs(lerr);
									lLdlcSqErr[ti][e]+=lerr*lerr;
								}
							}
						}
						if(finalStop && trueCard==finalCardinality){break;}
					}
				}

				// Merge into shared accumulators
				synchronized(lock){
					for(int ti2=0; ti2<finalNumThresholds; ti2++){
						if(lN[ti2]==0){continue;}
						nArr[ti2]+=lN[ti2];
						occSum[ti2]+=lOccSum[ti2];
						for(int e=0; e<NUM_OUT1; e++){
							sumErr[ti2][e]+=lSumErr[ti2][e];
							sumAbsErr[ti2][e]+=lSumAbsErr[ti2][e];
							sumSqErr[ti2][e]+=lSumSqErr[ti2][e];
						}
						for(int e=0; e<NUM_LDLC; e++){
							ldlcSumErr[ti2][e]+=lLdlcErr[ti2][e];
							ldlcSumAbsErr[ti2][e]+=lLdlcAbsErr[ti2][e];
							ldlcSumSqErr[ti2][e]+=lLdlcSqErr[ti2][e];
						}
					}
				}
			});
			threadArray[tid].start();
		}

		for(Thread t : threadArray){
			try{t.join();}catch(InterruptedException e){throw new RuntimeException(e);}
		}

		// Build ReportRows (reuses DDLCalibrationDriver's format logic)
		final ArrayList<DDLCalibrationDriver.ReportRow> mergedRows=new ArrayList<>(numThresholds);
		for(int ti=0; ti<numThresholds; ti++){
			if(nArr[ti]<1){continue;}
			mergedRows.add(new DDLCalibrationDriver.ReportRow(
				thresholds[ti], nArr[ti], occSum[ti], sumErr[ti], sumAbsErr[ti], sumSqErr[ti]));
		}

		// Print stderr summary (matching DDLCalibrationDriver2 format exactly)
		{
			final double elapsed=(System.nanoTime()-t0)*1e-9;
			System.err.println();
			System.err.println("=== DDL2 Calibration Summary ===");
			System.err.println("Type: "+loglogtype+"  Buckets: "+buckets+"  DDLs: "+numDDLs
				+"  MaxCard: "+cardinality+"  Rows: "+mergedRows.size()
				+"  Elapsed: "+String.format("%.1f", elapsed)+"s"
				+"  LowComplexity: iter="+iterations);
			System.err.println("--- Log-Weighted and Linear-Weighted Avg Absolute Error, Peak, Signed (lower = better) ---");
			final double[] totalMeanAbsErr=new double[NUM_EST];
			final double[] linWtAbsErr=new double[NUM_EST];
			final double[] peakMeanAbsErr=new double[NUM_EST];
			final double[] totalMeanSignedErr=new double[NUM_EST];
			final double[] totalCV=new double[NUM_EST];
			double linWeightSum=0;
			int cvRows=0;
			for(DDLCalibrationDriver.ReportRow row : mergedRows){
				final double card=row.trueCard;
				linWeightSum+=card;
				for(int e=0; e<NUM_EST; e++){
					final double meanErr=row.sumErr[e]/row.n;
					final double meanAbsAtRow=row.sumAbsErr[e]/row.n;
					totalMeanAbsErr[e]+=meanAbsAtRow;
					linWtAbsErr[e]+=meanAbsAtRow*card;
					totalMeanSignedErr[e]+=meanErr;
					if(meanAbsAtRow>peakMeanAbsErr[e]){peakMeanAbsErr[e]=meanAbsAtRow;}
					final double variance=row.sumSqErr[e]/row.n-meanErr*meanErr;
					final double stdev=Math.sqrt(Math.max(0, variance));
					final double denom=Math.abs(1.0+meanErr);
					if(denom>0){totalCV[e]+=stdev/denom;}
				}
				cvRows++;
			}
			final int rows=mergedRows.size();
			System.err.println(String.format("%-12s %-12s %-12s %-12s %-12s %s",
				"", "LogWtAbsErr", "LinWtAbsErr", "PeakAbsErr", "AvgSignErr", "AvgCV"));
			final String[] ESTIMATOR_NAMES=DDLCalibrationDriver.ESTIMATOR_NAMES;
			for(int e=0; e<Math.min(17, NUM_EST); e++){
				System.err.println(String.format("%-12s %.8f  %.8f  %.8f  %+.8f  %.8f",
					ESTIMATOR_NAMES[e], totalMeanAbsErr[e]/rows,
					linWeightSum>0 ? linWtAbsErr[e]/linWeightSum : 0,
					peakMeanAbsErr[e],
					rows>0 ? totalMeanSignedErr[e]/rows : 0,
					cvRows>0 ? totalCV[e]/cvRows : 0));
			}
			// LDLC summary rows
			final double[] ldlcTotalAbsErr=new double[NUM_LDLC];
			final double[] ldlcLinWtAbsErr=new double[NUM_LDLC];
			final double[] ldlcPeakAbsErr=new double[NUM_LDLC];
			final double[] ldlcTotalSignedErr=new double[NUM_LDLC];
			double ldlcLinWeightSum=0;
			for(int ti=0; ti<numThresholds; ti++){
				if(nArr[ti]<1){continue;}
				final double card=thresholds[ti];
				ldlcLinWeightSum+=card;
				for(int e=0; e<NUM_LDLC; e++){
					final double meanAbsAtRow=ldlcSumAbsErr[ti][e]/nArr[ti];
					ldlcTotalAbsErr[e]+=meanAbsAtRow;
					ldlcLinWtAbsErr[e]+=meanAbsAtRow*card;
					ldlcTotalSignedErr[e]+=ldlcSumErr[ti][e]/nArr[ti];
					if(meanAbsAtRow>ldlcPeakAbsErr[e]){ldlcPeakAbsErr[e]=meanAbsAtRow;}
				}
			}
			for(int e=0; e<NUM_LDLC; e++){
				System.err.println(String.format("%-12s %.8f  %.8f  %.8f  %+.8f  %.8f",
					LDLC_NAMES[e], ldlcTotalAbsErr[e]/rows,
					ldlcLinWeightSum>0 ? ldlcLinWtAbsErr[e]/ldlcLinWeightSum : 0,
					ldlcPeakAbsErr[e],
					rows>0 ? ldlcTotalSignedErr[e]/rows : 0, 0.0));
			}
			System.err.println();
		}

		// Print stdout: File 1 format with LDLC columns appended
		{
			final StringBuilder hdr=new StringBuilder(DDLCalibrationDriver.header1());
			for(int e=0; e<NUM_LDLC; e++){
				hdr.append('\t').append(LDLC_NAMES[e]).append("_err");
				hdr.append('\t').append(LDLC_NAMES[e]).append("_abs");
			}
			System.out.println(hdr);
		}
		int finalTi=0;
		for(DDLCalibrationDriver.ReportRow row : mergedRows){
			final StringBuilder sb=new StringBuilder(DDLCalibrationDriver.formatRow1(row));
			final int ti=finalTi;
			final int ni=row.n;
			for(int e=0; e<NUM_LDLC; e++){
				sb.append('\t').append(String.format("%.6f", ldlcSumErr[ti][e]/ni));
				sb.append('\t').append(String.format("%.6f", ldlcSumAbsErr[ti][e]/ni));
			}
			System.out.println(sb);
			finalTi++;
		}
	}

}
