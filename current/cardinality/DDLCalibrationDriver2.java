package cardinality;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import parse.Parse;
import rand.FastRandomXoshiro;
import shared.Shared;
import shared.Tools;

/**
 * Cache-friendly multithreaded calibration driver for DynamicDemiLog cardinality estimators.
 * <p>
 * Each thread processes one DDL at a time: create, feed all values sequentially,
 * record measurements at reporting thresholds, discard, repeat. This keeps the DDL's
 * working set in L1/L2 cache throughout its lifetime, eliminating the cache-miss storm
 * of DDLCalibrationDriver where each value touched all N DDLs.
 * <p>
 * After all threads finish, the main thread merges per-threshold accumulators and
 * prints all rows (File 1), then writes File 3 (v4 CF table).
 * <p>
 * Each thread is seeded with (masterSeed + threadId), giving one deterministic
 * FastRandomXoshiro per thread that generates both DDL seeds and hash values.
 * This guarantees bit-for-bit reproducibility across runs with the same parameters.
 * <p>
 * Parameters (key=value format):
 * <ul>
 *   <li>ddls=1000        - total DDL instances across all threads</li>
 *   <li>buckets=2048     - buckets per DDL (must be a multiple of 256)</li>
 *   <li>k=31             - k-mer length (needed by DDL constructor)</li>
 *   <li>maxmult=10       - stop when trueCard reaches buckets*maxmult</li>
 *   <li>reportfrac=0.01  - report every this fraction of current cardinality</li>
 *   <li>seed=12345       - master seed; thread t uses seed+t</li>
 *   <li>threads=1        - number of parallel simulation threads</li>
 *   <li>out1=path        - output file for per-threshold statistics (File 1)</li>
 *   <li>out3=path        - output file for v5 CF table (File 3)</li>
 *   <li>type=ddl         - DDL type: ddl, ddl2, ddl8, dll2, dll3, dll3v2, dll4</li>
 * </ul>
 *
 * @author Chloe
 * @date March 2026
 */
public class DDLCalibrationDriver2 {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of estimators reported by rawEstimates(). */
	static final int NUM_EST=DDLCalibrationDriver.NUM_EST;
	/** Estimator names in rawEstimates() index order. */
	static final String[] ESTIMATOR_NAMES=DDLCalibrationDriver.ESTIMATOR_NAMES;

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final long t0=System.nanoTime();
		int numDDLs=128;
		int buckets=2048;
		int k=31;
		long maxMult=10;
		double reportFrac=0.01;
		long masterSeed=12345L;
		int threads=Shared.threads();
		String out1="stdout.txt";
		String out3=null;
		String out4=null;
		String loglogtype="ddl";
		String notes="";
		String cffile=null;

		for(String arg : args){
			final String[] split=arg.split("=");
			if(split.length!=2){continue;}
			final String a=split[0].toLowerCase();
			final String b=split[1];
			if(a.equals("ddls") || a.equals("dlls")){numDDLs=Parse.parseIntKMG(b);}
			else if(a.equals("buckets")){buckets=Parse.parseIntKMG(b);}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("maxmult") || a.equals("mult")){maxMult=Parse.parseIntKMG(b);}
			else if(a.equals("reportfrac")){reportFrac=Double.parseDouble(b);}
			else if(a.equals("seed")){masterSeed=Long.parseLong(b);}
			else if(a.equals("valseed")){/*ignored: single-rng design uses seed only*/}
			else if(a.equals("threads") || a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("out") || a.equals("out1")){out1=b;}
			else if(a.equals("out3")){out3=b;}
			else if(a.equals("out4")){out4=b;}
			else if(a.equals("loglogtype") || a.equals("type")){loglogtype=b.toLowerCase();}
			else if(a.equals("cf") || a.equals("loglogcf")){CorrectionFactor.USE_CORRECTION=Parse.parseBoolean(b);}
			else if(a.equals("cardcf")){CardinalityTracker.USE_CARD_CF=Parse.parseBoolean(b);}
			else if(a.equals("cffile")){cffile=b; CorrectionFactor.USE_CORRECTION=true;}
			else if(a.equals("dlcalpha") || a.equals("alpha")){CardinalityStats.DLC_ALPHA=Float.parseFloat(b);}
			else if(a.equals("cfiters") || a.equals("cfiterations")){CardinalityStats.DEFAULT_CF_ITERS=Integer.parseInt(b);}
			else if(a.equals("cfdif") || a.equals("cfconvergence")){CardinalityStats.DEFAULT_CF_DIF=Double.parseDouble(b);}
			else if(a.equals("cfmult") || a.equals("minseedmult")){CardinalityStats.MIN_SEED_CF_MULT=Float.parseFloat(b);}
			else if(a.equals("minvfraction") || a.equals("minvk")){
				float x=Float.parseFloat(b);
				if(x<1){CardinalityStats.DLC_MIN_VK_FRACTION=x;}
				else{CardinalityStats.DLC_MIN_VK=(int)x;}
			}else if(a.equals("promotethreshold") || a.equals("pt")){
				DynamicLogLog3.PROMOTE_THRESHOLD=Integer.parseInt(b);
				DynamicLogLog3v2.PROMOTE_THRESHOLD=Integer.parseInt(b);
				DynamicLogLog4.PROMOTE_THRESHOLD=Integer.parseInt(b);
			}else if(a.equals("promotefrac") || a.equals("pf")){
				DynamicLogLog3v2.PROMOTE_FRAC=Float.parseFloat(b);
			}else if(a.equals("resetonpromote") || a.equals("rop")){
				DynamicLogLog3v2.RESET_ON_PROMOTE=Parse.parseBoolean(b);
			}else if(a.equals("earlypromote") || a.equals("ep")){
				DynamicLogLog3.EARLY_PROMOTE=Parse.parseBoolean(b);
				DynamicLogLog4.EARLY_PROMOTE=Parse.parseBoolean(b);
			}else if(a.equals("correctoverflow") || a.equals("co")){DynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
			}else if(a.equals("overflowscale") || a.equals("os")){DynamicLogLog3.OVERFLOW_SCALE=Double.parseDouble(b);
			}else if(a.equals("usestoredoverflow") || a.equals("uso")){DynamicLogLog3.USE_STORED_OVERFLOW=Parse.parseBoolean(b);
		}else if(a.equals("printdlctiers")){DDLCalibrationDriver.PRINT_DLC_TIERS=Parse.parseBoolean(b);
			}else if(a.equals("printstd")){DDLCalibrationDriver.PRINT_STD=Parse.parseBoolean(b);
			}else if(a.equals("out2")){
				System.err.println("Note: out2= is not supported by DDLCalibrationDriver2; ignoring.");
			}else if(a.equals("notes")){notes=b.replace('_',' ');
			}else if(a.equals("frozenhistory") || a.equals("frozen")){
				UltraLogLog8.FROZEN_HISTORY=Parse.parseBoolean(b);
			}else if(a.equals("printcv") || a.equals("cv")){
				DDLCalibrationDriver.PRINT_CV=Parse.parseBoolean(b);
			}else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}

		final long maxTrue=(long)buckets*maxMult;
		final int numThreads=Math.min(threads, numDDLs);

		// Force class initialization of the DDL type BEFORE loading the user CF file.
		// Each DDL class loads its own default (v4) CF on first use, which resets
		// CorrectionFactor.v1Buckets=0 and clobbers any v5 table loaded earlier.
		// Creating a dummy instance here triggers that class init; then we overwrite
		// the default CF with the user-provided file.
		DDLCalibrationDriver.makeInstance(loglogtype, buckets, k, 0L, 0);
		if(cffile!=null){CorrectionFactor.initialize(cffile, buckets);}
		final long[] thresholds=DDLCalibrationDriver.computeThresholds(maxTrue, reportFrac);
		final int numThresholds=thresholds.length;

		// Create and start worker threads.
		// Each thread t is seeded with (masterSeed+t) — one RNG generates both DDL seeds and values.
		final CalibrationThread2[] calThreads=new CalibrationThread2[numThreads];
		for(int t=0; t<numThreads; t++){
			final int threadDDLs=numDDLs/numThreads+(t<numDDLs%numThreads ? 1 : 0);
			calThreads[t]=new CalibrationThread2(
				masterSeed+t, threadDDLs, thresholds, buckets, k, maxTrue, loglogtype, out4);
		}
		for(CalibrationThread2 ct : calThreads){ct.start();}

		// Wait for all threads
		for(CalibrationThread2 ct : calThreads){
			try{ct.join();}catch(InterruptedException e){Thread.currentThread().interrupt();}
			if(!ct.success){System.err.println("Warning: a calibration thread did not complete successfully.");}
		}

		// Merge per-threshold accumulators across threads
		// sumErr[ti][e], sumAbsErr[ti][e], sumSqErr[ti][e], occSum[ti], n[ti]
		final double[][] sumErr=new double[numThresholds][NUM_EST];
		final double[][] sumAbsErr=new double[numThresholds][NUM_EST];
		final double[][] sumSqErr=new double[numThresholds][NUM_EST];
		final double[] occSum=new double[numThresholds];
		final int[] n=new int[numThresholds];
		for(CalibrationThread2 ct : calThreads){
			for(int ti=0; ti<numThresholds; ti++){
				n[ti]+=ct.n[ti];
				occSum[ti]+=ct.occSum[ti];
				for(int e=0; e<NUM_EST; e++){
					sumErr[ti][e]+=ct.sumErr[ti][e];
					sumAbsErr[ti][e]+=ct.sumAbsErr[ti][e];
					sumSqErr[ti][e]+=ct.sumSqErr[ti][e];
				}
			}
		}

		// Build ReportRows for output (reuses DDLCalibrationDriver's row/format logic)
		final ArrayList<DDLCalibrationDriver.ReportRow> mergedRows=new ArrayList<>(numThresholds);
		for(int ti=0; ti<numThresholds; ti++){
			if(n[ti]<1){continue;}
			// ReportRow constructor takes arrays of length NUM_OUT1; pad with zeros for MedianLC slot
			final int NUM_OUT1=DDLCalibrationDriver.NUM_OUT1;
			final double[] rErr=new double[NUM_OUT1];
			final double[] rAbsErr=new double[NUM_OUT1];
			final double[] rSqErr=new double[NUM_OUT1];
			System.arraycopy(sumErr[ti], 0, rErr, 0, NUM_EST);
			System.arraycopy(sumAbsErr[ti], 0, rAbsErr, 0, NUM_EST);
			System.arraycopy(sumSqErr[ti], 0, rSqErr, 0, NUM_EST);
			mergedRows.add(new DDLCalibrationDriver.ReportRow(
				thresholds[ti], n[ti], occSum[ti], rErr, rAbsErr, rSqErr));
		}

		// Print stderr summary
		{
			final double elapsed=(System.nanoTime()-t0)*1e-9;
			System.err.println();
			System.err.println("=== DDL2 Calibration Summary ===");
			System.err.println("Type: "+loglogtype+"  Buckets: "+buckets+"  DDLs: "+numDDLs
				+"  MaxCard: "+maxTrue+"  Rows: "+mergedRows.size()
				+"  Elapsed: "+String.format("%.1f", elapsed)+"s");
			System.err.println("--- Avg and Peak Mean Absolute Error, Avg CV (lower = better) ---");
			final double[] totalMeanAbsErr=new double[NUM_EST];
			final double[] peakMeanAbsErr=new double[NUM_EST];
			final double[] totalCV=new double[NUM_EST];
			int cvRows=0;
			for(DDLCalibrationDriver.ReportRow row : mergedRows){
				for(int e=0; e<NUM_EST; e++){
					final double meanErr=row.sumErr[e]/row.n;
					final double meanAbsAtRow=row.sumAbsErr[e]/row.n;
					totalMeanAbsErr[e]+=meanAbsAtRow;
					if(meanAbsAtRow>peakMeanAbsErr[e]){peakMeanAbsErr[e]=meanAbsAtRow;}
					final double variance=row.sumSqErr[e]/row.n-meanErr*meanErr;
					final double stdev=Math.sqrt(Math.max(0, variance));
					if(meanAbsAtRow>0){totalCV[e]+=stdev/meanAbsAtRow;}
				}
				cvRows++;
			}
			final int rows=mergedRows.size();
			System.err.println(String.format("%-12s %-12s %-12s %s", "", "AvgAbsErr", "PeakAbsErr", "AvgCV"));
			final int[] keyIdx={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
			for(int ki=0; ki<keyIdx.length; ki++){
				final int e=keyIdx[ki];
				if(e>=NUM_EST){continue;}
				System.err.println(String.format("%-12s %.8f  %.8f  %.8f",
					ESTIMATOR_NAMES[e], totalMeanAbsErr[e]/rows, peakMeanAbsErr[e],
					cvRows>0 ? totalCV[e]/cvRows : 0));
			}
			System.err.println();
		}

		// Write File 1
		final ByteStreamWriter bsw1=new ByteStreamWriter(out1, true, false, false);
		bsw1.start();
		bsw1.println(DDLCalibrationDriver.header1());
		for(DDLCalibrationDriver.ReportRow row : mergedRows){
			bsw1.println(DDLCalibrationDriver.formatRow1(row));
		}
		bsw1.poisonAndWait();

		// Write File 3 (v5 CF table)
		if(out3!=null){
			try(PrintStream ps=new PrintStream(new FileOutputStream(out3))){
				DDLCalibrationDriver.printHistogramV5(mergedRows, ps,
						loglogtype, buckets, numDDLs, maxTrue,
						CorrectionFactor.USE_CORRECTION,
						DynamicLogLog3.EARLY_PROMOTE,
						DynamicLogLog3.PROMOTE_THRESHOLD,
						notes);
				System.err.println("V5 CF table written to "+out3);
			}catch(Exception e){
				System.err.println("Error writing CF file: "+e.getMessage());
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Cache-friendly calibration thread.
	 * Processes one DDL at a time: create, feed maxTrue values, record stats at each threshold, discard.
	 * Uses FastRandomXoshiro for value generation — no LongHashSet needed (PRNG uniqueness assumed).
	 * Accumulates per-threshold sums locally; main thread merges after all threads join.
	 */
	static final class CalibrationThread2 extends Thread {

		CalibrationThread2(long seed, int numDDLs, long[] thresholds,
			int buckets, int k, long maxTrue, String loglogtype, String out4){
			this.seed=seed;
			this.numDDLs=numDDLs;
			this.thresholds=thresholds;
			this.buckets=buckets;
			this.k=k;
			this.maxTrue=maxTrue;
			this.loglogtype=loglogtype;
			this.out4=out4;
			final int nt=thresholds.length;
			n=new int[nt];
			occSum=new double[nt];
			sumErr=new double[nt][NUM_EST];
			sumAbsErr=new double[nt][NUM_EST];
			sumSqErr=new double[nt][NUM_EST];
		}

		@Override
		public void run(){
			try{
				runInner();
				success=true;
			}catch(Exception e){
				e.printStackTrace();
			}
		}

		void runInner(){
			// One RNG per thread, seeded with (masterSeed+threadId).
			// It generates the DDL seed for each DDL, then all hash values for that DDL.
			// Fully deterministic: same seed → same output regardless of system state.
			final FastRandomXoshiro rng=new FastRandomXoshiro(seed);
			java.io.PrintWriter out4Pw=null;
			if(out4!=null){try{out4Pw=new java.io.PrintWriter(new java.io.FileOutputStream(out4));}catch(Exception e){e.printStackTrace();}}
			for(int di=0; di<numDDLs; di++){
				final long ddlSeed=rng.nextLong()&Long.MAX_VALUE; // positive: literal seed, not nanoTime
				final CardinalityTracker ddl=DDLCalibrationDriver.makeInstance(
					loglogtype, buckets, k, ddlSeed, 0);

				int ti=0;
				for(long trueCard=1; trueCard<=maxTrue; trueCard++){
					ddl.hashAndStore(rng.nextLong());


					// Check reporting threshold
					if(trueCard>=thresholds[ti]){
						final double occ=ddl.occupancy();
						final double[] est=ddl.rawEstimates();
						occSum[ti]+=occ;
						n[ti]++;
						if(di==0 && out4Pw!=null){
							int[] raw=null, corr=null;
							int mz=0;
							if(ddl instanceof DynamicLogLog3){
								final DynamicLogLog3 d=(DynamicLogLog3)ddl;
								raw=d.lastRawNlz; corr=d.lastCorrNlz; mz=d.getMinZeros();
							}else if(ddl instanceof DynamicLogLog4){
								final DynamicLogLog4 d=(DynamicLogLog4)ddl;
								raw=d.lastRawNlz; corr=d.lastCorrNlz; mz=d.getMinZeros();
							}
							if(raw!=null){
								out4Pw.print(trueCard+"\t"+mz);
								for(int tt=0; tt<20; tt++){out4Pw.print("\t"+raw[tt]);}
								for(int tt=0; tt<20; tt++){out4Pw.print("\t"+corr[tt]);}
								out4Pw.println();
							}
						}
						for(int e=0; e<NUM_EST; e++){
							final double err=(est[e]-trueCard)/(double)trueCard;
							sumErr[ti][e]+=err;
							sumAbsErr[ti][e]+=Math.abs(err);
							sumSqErr[ti][e]+=err*err;
						}
						ti++;
						if(ti>=thresholds.length){break;}
					}
				}
			}
			if(out4Pw!=null){out4Pw.close();}
		}

		/** Seed for this thread's RNG (masterSeed + threadId). */
		final long seed;
		/** Number of DDLs this thread processes. */
		final int numDDLs;
		/** Pre-computed reporting thresholds (shared read-only). */
		final long[] thresholds;
		/** Bucket count. */
		final int buckets;
		/** k-mer length. */
		final int k;
		/** Maximum true cardinality. */
		final long maxTrue;
		/** DDL type string. */
		final String loglogtype;
		/** Output file for per-tier nlzCounts (first DDL only). */
		final String out4;

		/** Number of DDLs that contributed to each threshold row. */
		final int[] n;
		/** Sum of occupancy at each threshold. */
		final double[] occSum;
		/** Sum of signed relative errors per threshold per estimator. */
		final double[][] sumErr;
		/** Sum of absolute relative errors per threshold per estimator. */
		final double[][] sumAbsErr;
		/** Sum of squared relative errors per threshold per estimator. */
		final double[][] sumSqErr;

		/** Set to true after runInner() completes without exception. */
		volatile boolean success=false;
	}

}
