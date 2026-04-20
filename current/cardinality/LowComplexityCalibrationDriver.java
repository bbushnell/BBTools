package cardinality;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicLong;
import rand.FastRandomXoshiro;
import parse.Parse;
import parse.Parser;
import shared.Shared;

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
		float iterations=4;     // 0 = stop on saturation; >0 = run for cardinality*iterations adds
		int buckets=2048;
		int k=31;
		int threads=4;
		long masterSeed=1;
		double reportFrac=0.01;
		String loglogtype="dll4";
		boolean minRand=true;   // true=min(rand,rand) low-complexity skew; false=uniform like HC
		String cffile=null;
		Shared.capThreads(4);
		Parser.printSetThreads=false;

		for(String arg : args){
			final String[] split=arg.split("=");
			if(split.length!=2){continue;}
			final String a=split[0].toLowerCase();
			final String b=split[1];
			if(a.equals("cardinality") || a.equals("card") || a.equals("maxcard")){cardinality=Parse.parseIntKMG(b);}
			else if(a.equals("iterations") || a.equals("iter")){iterations=Float.parseFloat(b);}
			else if(a.equals("ddls") || a.equals("dlls")){numDDLs=Parse.parseIntKMG(b);}
			else if(a.equals("buckets")){buckets=Parse.parseIntKMG(b);}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("seed")){masterSeed=Long.parseLong(b);}
			else if(a.equals("reportfrac")){reportFrac=Double.parseDouble(b);}
			else if(a.equals("cffile")){cffile=b; CorrectionFactor.USE_CORRECTION=true;}
			else if(a.startsWith("hsbtable")){System.err.println("Note: hsbtable= is deprecated; HSB values are now hardcoded in StateTable.");}
			else if(CardinalityParser.parse(arg, a, b)){}
			else if(a.equals("minrand") || a.equals("skew")){minRand=Parse.parseBoolean(b);}
			else if(!Parser.parseStatic(arg, a, b)){throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}
		loglogtype=Parser.loglogType;
		threads=Shared.threads();

		// Reuse constants from DDLCalibrationDriver/Driver2
		final int NUM_EST=DDLCalibrationDriver.NUM_EST;
		final int NUM_OUT1=DDLCalibrationDriver.NUM_OUT1;
		final int NUM_LDLC=DDLCalibrationDriver2.NUM_LDLC;
		final String[] LDLC_NAMES=DDLCalibrationDriver2.LDLC_NAMES;

		// Initialize global cardinality state (CF tables, formula coefficients, etc.)
		CardinalityParser.initializeAll(loglogtype, buckets, k, cffile, null, true);

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
		final int[] ldlcIdx={0, 1, 2, 4, 5, 6, 7}; // Indices into ldlcEstimate() array
		final boolean finalMinRand=minRand;
		final AtomicLong nextDDL=new AtomicLong(0);
		final int finalNumDDLs=numDDLs;

		final long t0=System.nanoTime();
		final Thread[] threadArray=new Thread[threads];
		for(int tid=0; tid<threads; tid++){

			threadArray[tid]=new Thread(()->{
				// Allocate accumulators on this thread's NUMA node
				final double[][] lSumErr=new double[finalNumThresholds][NUM_OUT1];
				final double[][] lSumAbsErr=new double[finalNumThresholds][NUM_OUT1];
				final double[][] lSumSqErr=new double[finalNumThresholds][NUM_OUT1];
				final double[][] lLdlcErr=new double[finalNumThresholds][NUM_LDLC];
				final double[][] lLdlcAbsErr=new double[finalNumThresholds][NUM_LDLC];
				final double[][] lLdlcSqErr=new double[finalNumThresholds][NUM_LDLC];
				final double[] lOccSum=new double[finalNumThresholds];
				final int[] lN=new int[finalNumThresholds];
				final FastRandomXoshiro rng=new FastRandomXoshiro(0);
				final BitSet seen=new BitSet(finalCardinality);

				long estIdx;
				while((estIdx=nextDDL.getAndIncrement())<finalNumDDLs){
					rng.setSeed(seeds[(int)estIdx]);
					seen.clear();
					int trueCard=0;
					int ti=0;
					final CardinalityTracker est=DDLCalibrationDriver.makeInstance(
						finalType, finalBuckets, finalK, seeds[(int)estIdx], 0);

					for(long add=0; add<finalTotalAdds; add++){
						final int pos=finalMinRand ? Math.min(rng.nextInt(finalCardinality), rng.nextInt(finalCardinality)) : rng.nextInt(finalCardinality);
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
							// LDLC columns: named getters from CardStats
							if(est.getClass()==UltraDynamicLogLog6.class){
								final UltraDynamicLogLog6 u=(UltraDynamicLogLog6)est;
								final CardStats cs=u.consumeLastSummarized();
								final double[] ldlcVals={cs.ldlc(), cs.dlcSbs(), cs.hc(),
									u.fgraEstimate(), cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2()};
								for(int e=0; e<DDLCalibrationDriver2.NUM_LDLC_BASE; e++){
									final double v=ldlcVals[e];
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][e]+=lerr;
									lLdlcAbsErr[ti][e]+=Math.abs(lerr);
									lLdlcSqErr[ti][e]+=lerr*lerr;
								}
								// HLDLC: 50/50 blend of Hybrid+2 and LDLC
								{
									final float hw=est.hldlcWeight(); final double hldlc=hw*ldlcVals[0]+(1-hw)*ldlcVals[6];
									final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr; lLdlcAbsErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=Math.abs(lerr); lLdlcSqErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr*lerr;
								}
							}else if(est.getClass()==BankedDynamicLogLog5.class){
								final BankedDynamicLogLog5 c=(BankedDynamicLogLog5)est;
								final CardStats cs=c.consumeLastSummarized();
								final double[] ldlcVals={cs.ldlc(), cs.dlcSbs(), cs.hc(),
									0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2()};
								for(int e=0; e<DDLCalibrationDriver2.NUM_LDLC_BASE; e++){
									final double v=ldlcVals[e];
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][e]+=lerr;
									lLdlcAbsErr[ti][e]+=Math.abs(lerr);
									lLdlcSqErr[ti][e]+=lerr*lerr;
								}
								{
									final float hw=est.hldlcWeight(); final double hldlc=hw*ldlcVals[0]+(1-hw)*ldlcVals[6];
									final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr; lLdlcAbsErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=Math.abs(lerr); lLdlcSqErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr*lerr;
								}
							}else if(est.getClass()==CompressedDynamicLogLog4.class){
								final CompressedDynamicLogLog4 c=(CompressedDynamicLogLog4)est;
								final CardStats cs=c.consumeLastSummarized();
								final double[] ldlcVals={cs.ldlc(), cs.dlcSbs(), cs.hc(),
									0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2()};
								for(int e=0; e<DDLCalibrationDriver2.NUM_LDLC_BASE; e++){
									final double v=ldlcVals[e];
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][e]+=lerr;
									lLdlcAbsErr[ti][e]+=Math.abs(lerr);
									lLdlcSqErr[ti][e]+=lerr*lerr;
								}
								{
									final float hw=est.hldlcWeight(); final double hldlc=hw*ldlcVals[0]+(1-hw)*ldlcVals[6];
									final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr; lLdlcAbsErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=Math.abs(lerr); lLdlcSqErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr*lerr;
								}
							}else if(est.getClass()==CompressedDynamicLogLog5.class){
								final CompressedDynamicLogLog5 c=(CompressedDynamicLogLog5)est;
								final CardStats cs=c.consumeLastSummarized();
								final double[] ldlcVals={cs.ldlc(), cs.dlcSbs(), cs.hc(),
									0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2()};
								for(int e=0; e<7; e++){
									final double v=ldlcVals[e];
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][e]+=lerr;
									lLdlcAbsErr[ti][e]+=Math.abs(lerr);
									lLdlcSqErr[ti][e]+=lerr*lerr;
								}
								{
									final float hw=est.hldlcWeight(); final double hldlc=hw*ldlcVals[0]+(1-hw)*ldlcVals[6];
									final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][11]+=lerr; lLdlcAbsErr[ti][11]+=Math.abs(lerr); lLdlcSqErr[ti][11]+=lerr*lerr;
								}
							}else if(est.getClass()==BankedCompressedDynamicLogLog5.class){
								final BankedCompressedDynamicLogLog5 c=(BankedCompressedDynamicLogLog5)est;
								final CardStats cs=c.consumeLastSummarized();
								final double[] ldlcVals={cs.ldlc(), cs.dlcSbs(), cs.hc(),
									0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2()};
								for(int e=0; e<7; e++){
									final double v=ldlcVals[e];
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][e]+=lerr;
									lLdlcAbsErr[ti][e]+=Math.abs(lerr);
									lLdlcSqErr[ti][e]+=lerr*lerr;
								}
								{
									final float hw=est.hldlcWeight(); final double hldlc=hw*ldlcVals[0]+(1-hw)*ldlcVals[6];
									final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][11]+=lerr; lLdlcAbsErr[ti][11]+=Math.abs(lerr); lLdlcSqErr[ti][11]+=lerr*lerr;
								}
							}else if(est.getClass()==ArithmeticCompressedDynamicLogLog5.class){
								final ArithmeticCompressedDynamicLogLog5 c=(ArithmeticCompressedDynamicLogLog5)est;
								final CardStats cs=c.consumeLastSummarized();
								final double[] ldlcVals={cs.ldlc(), cs.dlcSbs(), cs.hc(),
									0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2()};
								for(int e=0; e<7; e++){
									final double v=ldlcVals[e];
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][e]+=lerr;
									lLdlcAbsErr[ti][e]+=Math.abs(lerr);
									lLdlcSqErr[ti][e]+=lerr*lerr;
								}
								{
									final float hw=est.hldlcWeight(); final double hldlc=hw*ldlcVals[0]+(1-hw)*ldlcVals[6];
									final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][11]+=lerr; lLdlcAbsErr[ti][11]+=Math.abs(lerr); lLdlcSqErr[ti][11]+=lerr*lerr;
								}
							}else if(est.getClass()==CompressedDynamicLogLog3.class){
								final CompressedDynamicLogLog3 c=(CompressedDynamicLogLog3)est;
								final CardStats cs=c.consumeLastSummarized();
								final double[] ldlcVals={cs.ldlc(), cs.dlcSbs(), cs.hc(),
									0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2()};
								for(int e=0; e<DDLCalibrationDriver2.NUM_LDLC_BASE; e++){
									final double v=ldlcVals[e];
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][e]+=lerr;
									lLdlcAbsErr[ti][e]+=Math.abs(lerr);
									lLdlcSqErr[ti][e]+=lerr*lerr;
								}
								{
									final float hw=est.hldlcWeight(); final double hldlc=hw*ldlcVals[0]+(1-hw)*ldlcVals[6];
									final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr; lLdlcAbsErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=Math.abs(lerr); lLdlcSqErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr*lerr;
								}
							}else if(est.getClass()==BankedCompressedDynamicLogLog3.class){
								final BankedCompressedDynamicLogLog3 c=(BankedCompressedDynamicLogLog3)est;
								final CardStats cs=c.consumeLastSummarized();
								final double[] ldlcVals={cs.ldlc(), cs.dlcSbs(), cs.hc(),
									0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2()};
								for(int e=0; e<DDLCalibrationDriver2.NUM_LDLC_BASE; e++){
									final double v=ldlcVals[e];
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][e]+=lerr;
									lLdlcAbsErr[ti][e]+=Math.abs(lerr);
									lLdlcSqErr[ti][e]+=lerr*lerr;
								}
								{
									final float hw=est.hldlcWeight(); final double hldlc=hw*ldlcVals[0]+(1-hw)*ldlcVals[6];
									final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr; lLdlcAbsErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=Math.abs(lerr); lLdlcSqErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr*lerr;
								}
							}else if(est.getClass()==ProtoLogLog16c.class){
								final double[] ldlcR=((ProtoLogLog16c)est).ldlcEstimate();
								for(int e=0; e<DDLCalibrationDriver2.NUM_LDLC_BASE; e++){
									final double v=ldlcR[ldlcIdx[e]];
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][e]+=lerr;
									lLdlcAbsErr[ti][e]+=Math.abs(lerr);
									lLdlcSqErr[ti][e]+=lerr*lerr;
								}
							}
							// LC_noMicro and SBS_noMicro from rawEstimates array
							{
								final int[] extraIdx={AbstractCardStats.LC_NOMICRO_IDX, AbstractCardStats.SBS_NOMICRO_IDX};
								for(int e=0; e<extraIdx.length; e++){
									final double v=(extraIdx[e]<estimates.length ? estimates[extraIdx[e]] : 0);
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									lLdlcErr[ti][7+e]+=lerr;
									lLdlcAbsErr[ti][7+e]+=Math.abs(lerr);
									lLdlcSqErr[ti][7+e]+=lerr*lerr;
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
			System.err.println("=== LC Calibration Summary ===");
			final int actualBk=DDLCalibrationDriver.makeInstance(loglogtype, buckets, k, 0, 0).actualBuckets();
			System.err.println("Type: "+loglogtype+"  Buckets: "+actualBk+"  DDLs: "+numDDLs
				+"  MaxCard: "+cardinality+"  Rows: "+mergedRows.size()
				+"  Elapsed: "+String.format("%.1f", elapsed)+"s"
				+"  LowComplexity: iter="+iterations);
			System.err.println("--- Log-Weighted and Width-Weighted Avg Absolute Error, Peak, Signed (lower = better) ---");
			final int rows=mergedRows.size();
			// Bucket widths: trailing difference of trueCard (defined for row 0 as trueCard[0]-0, typically 1).
			final double[] widths=new double[rows];
			widths[0]=Math.max(1, mergedRows.get(0).trueCard);
			for(int i=1; i<rows; i++){
				widths[i]=Math.max(1, mergedRows.get(i).trueCard-mergedRows.get(i-1).trueCard);
			}

			final double[] totalMeanAbsErr=new double[NUM_EST];
			final double[] widthWtAbsErr=new double[NUM_EST];
			final double[] countWtAbsErr=new double[NUM_EST];
			final double[] peakMeanAbsErr=new double[NUM_EST];
			final double[] totalMeanSignedErr=new double[NUM_EST];
			final double[] totalCV=new double[NUM_EST];
			double widthWeightSum=0;
			double countWeightSum=0;
			int cvRows=0;
			int rowIdx=0;
			for(DDLCalibrationDriver.ReportRow row : mergedRows){
				final double w=widths[rowIdx];
				widthWeightSum+=w;
				countWeightSum+=row.n;
				for(int e=0; e<NUM_EST; e++){
					final double meanErr=row.sumErr[e]/row.n;
					final double meanAbsAtRow=row.sumAbsErr[e]/row.n;
					totalMeanAbsErr[e]+=meanAbsAtRow;
					widthWtAbsErr[e]+=meanAbsAtRow*w;
					countWtAbsErr[e]+=meanAbsAtRow*row.n;
					totalMeanSignedErr[e]+=meanErr;
					if(meanAbsAtRow>peakMeanAbsErr[e]){peakMeanAbsErr[e]=meanAbsAtRow;}
					final double variance=row.sumSqErr[e]/row.n-meanErr*meanErr;
					final double stdev=Math.sqrt(Math.max(0, variance));
					final double denom=Math.abs(1.0+meanErr);
					if(denom>0){totalCV[e]+=stdev/denom;}
				}
				cvRows++;
				rowIdx++;
			}
			System.err.println(String.format("%-12s %-13s %-13s %-13s %-13s %-13s %s",
				"", "LogWtAbsErr", "WidthWtAbsErr", "CountWtAbsErr", "PeakAbsErr", "AvgSignErr", "AvgCV"));
			final String[] ESTIMATOR_NAMES=DDLCalibrationDriver.ESTIMATOR_NAMES;
			for(int e=0; e<Math.min(17, NUM_EST); e++){
				System.err.println(String.format("%-12s %.8f   %.8f   %.8f   %.8f   %+.8f   %.8f",
					ESTIMATOR_NAMES[e], totalMeanAbsErr[e]/rows,
					widthWeightSum>0 ? widthWtAbsErr[e]/widthWeightSum : 0,
					countWeightSum>0 ? countWtAbsErr[e]/countWeightSum : 0,
					peakMeanAbsErr[e],
					rows>0 ? totalMeanSignedErr[e]/rows : 0,
					cvRows>0 ? totalCV[e]/cvRows : 0));
			}
			// LDLC summary rows — trailing-difference widths, gap-aware across skipped empty thresholds
			final double[] tWidths=new double[numThresholds];
			{
				long prev=0;
				for(int ti=0; ti<numThresholds; ti++){
					if(nArr[ti]<1){tWidths[ti]=0; continue;}
					tWidths[ti]=Math.max(1, thresholds[ti]-prev);
					prev=thresholds[ti];
				}
			}

			final double[] ldlcTotalAbsErr=new double[NUM_LDLC];
			final double[] ldlcWidthWtAbsErr=new double[NUM_LDLC];
			final double[] ldlcCountWtAbsErr=new double[NUM_LDLC];
			final double[] ldlcPeakAbsErr=new double[NUM_LDLC];
			final double[] ldlcTotalSignedErr=new double[NUM_LDLC];
			double ldlcWidthWeightSum=0;
			double ldlcCountWeightSum=0;
			for(int ti=0; ti<numThresholds; ti++){
				if(nArr[ti]<1){continue;}
				final double w=tWidths[ti];
				ldlcWidthWeightSum+=w;
				ldlcCountWeightSum+=nArr[ti];
				for(int e=0; e<NUM_LDLC; e++){
					final double meanAbsAtRow=ldlcSumAbsErr[ti][e]/nArr[ti];
					ldlcTotalAbsErr[e]+=meanAbsAtRow;
					ldlcWidthWtAbsErr[e]+=meanAbsAtRow*w;
					ldlcCountWtAbsErr[e]+=meanAbsAtRow*nArr[ti];
					ldlcTotalSignedErr[e]+=ldlcSumErr[ti][e]/nArr[ti];
					if(meanAbsAtRow>ldlcPeakAbsErr[e]){ldlcPeakAbsErr[e]=meanAbsAtRow;}
				}
			}
			for(int e=0; e<NUM_LDLC; e++){
				System.err.println(String.format("%-12s %.8f   %.8f   %.8f   %.8f   %+.8f   %.8f",
					LDLC_NAMES[e], ldlcTotalAbsErr[e]/rows,
					ldlcWidthWeightSum>0 ? ldlcWidthWtAbsErr[e]/ldlcWidthWeightSum : 0,
					ldlcCountWeightSum>0 ? ldlcCountWtAbsErr[e]/ldlcCountWeightSum : 0,
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
