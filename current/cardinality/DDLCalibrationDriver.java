package cardinality;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;

import map.LongHashSet;
import parse.Parse;
import shared.Shared;
import shared.Tools;

/**
 * Multithreaded calibration driver for DynamicDemiLog cardinality estimators.
 * <p>
 * N independent threads each simulate a full DDL run (own seeds, own value streams).
 * Each thread posts per-interval statistics to a dedicated queue; the main thread
 * collects one row per interval from all threads, merges them, and prints live —
 * preserving the single-threaded line-by-line output behavior.
 * <p>
 * At program end, per-thread occupancy histograms are aggregated and written to a
 * second output file (out2=). The histogram records raw estimator values at exact
 * 1/256-occupancy checkpoints, enabling derivation of per-occupancy compensation
 * curves for every estimator.
 * <p>
 * Parameters (key=value format):
 * <ul>
 *   <li>ddls=1000        - total DDL instances across all threads</li>
 *   <li>buckets=2048     - buckets per DDL (must be a multiple of 256)</li>
 *   <li>k=31             - k-mer length (needed by DDL constructor)</li>
 *   <li>maxmult=10       - stop when trueCard reaches buckets*maxmult</li>
 *   <li>reportfrac=0.01  - report every this fraction of current cardinality</li>
 *   <li>seed=12345       - master seed for DDL seed generation</li>
 *   <li>valseed=42       - master seed for per-thread value generation</li>
 *   <li>threads=1        - number of parallel simulation threads</li>
 *   <li>out2=path        - output file for occupancy histogram (File 2)</li>
 * </ul>
 *
 * @author Chloe
 * @date March 2026
 */
public class DDLCalibrationDriver {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of estimators reported by rawEstimates(). */
	static final int NUM_EST=10;
	/** Estimator names in rawEstimates() index order. */
	static final String[] ESTIMATOR_NAMES={
		"Mean","HMean","HMeanM","GMean","HLL","LC","Hybrid","MWA","MedianCorr","Mean99"
	};
	/** Which estimators get a CF column in File 2. LC and Hybrid are excluded:
	 *  LC never uses CF; Hybrid is pre-corrected from its components. */
	static final boolean[] NEEDS_CF={true,true,true,true,true,false,false,true,true,true};
	/** Extra out1-only derived columns appended after NUM_EST. */
	static final int NUM_OUT1=NUM_EST+1; // +1 for MedianLC
	static final String MEDIAN_LC="MedianLC";

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		int numDDLs=128;
		int buckets=2048;
		int k=31;
		long maxMult=10;
		double reportFrac=0.01;
		long masterSeed=12345L;
		long valSeed=42L;
		int threads=Shared.threads();
		int step=1;
		String out2=null;
		String loglogtype="ddl"; // ddl, ddl2, ddl8, dll4

		for(String arg : args){
			final String[] split=arg.split("=");
			if(split.length!=2){continue;}
			final String a=split[0].toLowerCase();
			final String b=split[1];
			if(a.equals("ddls")){numDDLs=Parse.parseIntKMG(b);}
			else if(a.equals("buckets")){buckets=Parse.parseIntKMG(b);}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("maxmult")){maxMult=Parse.parseIntKMG(b);}
			else if(a.equals("reportfrac")){reportFrac=Double.parseDouble(b);}
			else if(a.equals("seed")){masterSeed=Long.parseLong(b);}
			else if(a.equals("valseed")){valSeed=Long.parseLong(b);}
			else if(a.equals("threads") || a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("step") || a.equals("resolution") || a.equals("res")){step=Integer.parseInt(b);}
			else if(a.equals("out2")){out2=b;}
			else if(a.equals("loglogtype") || a.equals("type")){loglogtype=b.toLowerCase();}
			else if(a.equals("cf") || a.equals("loglogcf")){CorrectionFactor.USE_CORRECTION=Parse.parseBoolean(b);}
			else{assert(false) : "Unknown parameter '"+arg+"'";}
		}

		assert(buckets%step==0) : "buckets must be a multiple of step; buckets="+buckets+", step="+step;
		final int numSlots=buckets/step+1;

		final long maxTrue=(long)buckets*maxMult;
		final int numThreads=Math.min(threads, numDDLs);
		final long[] thresholds=computeThresholds(maxTrue, reportFrac);

		// Pre-generate DDL seeds sequentially from masterSeed for reproducibility
		final Random seedRng=new Random(masterSeed);
		final int[] threadSizes=new int[numThreads];
		for(int t=0; t<numThreads; t++){
			threadSizes[t]=numDDLs/numThreads+(t<numDDLs%numThreads ? 1 : 0);
		}
		final long[][] ddlSeeds=new long[numThreads][];
		for(int t=0; t<numThreads; t++){
			ddlSeeds[t]=new long[threadSizes[t]];
			for(int i=0; i<threadSizes[t]; i++){
				ddlSeeds[t][i]=seedRng.nextLong()&Long.MAX_VALUE;
			}
		}

		// Pre-generate independent val seeds for each thread
		final Random valSeedRng=new Random(valSeed);
		final long[] threadValSeeds=new long[numThreads];
		for(int t=0; t<numThreads; t++){threadValSeeds[t]=valSeedRng.nextLong();}

		// Create threads
		final CalibrationThread[] calThreads=new CalibrationThread[numThreads];
		for(int t=0; t<numThreads; t++){
			final CardinalityTracker[] ddls=new CardinalityTracker[threadSizes[t]];
			for(int i=0; i<ddls.length; i++){
				ddls[i]=makeInstance(loglogtype, buckets, k, ddlSeeds[t][i], 0);
			}
			calThreads[t]=new CalibrationThread(ddls, thresholds, threadValSeeds[t], buckets, maxTrue, numSlots, step);
		}

		// Print File 1 header and start threads
		System.out.println(header1());
		for(CalibrationThread ct : calThreads){ct.start();}

		// Main polling loop: collect one merged row per threshold, print live
		for(int ti=0; ti<thresholds.length; ti++){
			final ReportRow[] rows=new ReportRow[numThreads];
			for(int t=0; t<numThreads; t++){
				try{
					rows[t]=calThreads[t].queue.take();
				}catch(InterruptedException e){
					Thread.currentThread().interrupt();
					throw new RuntimeException(e);
				}
			}
			printRow1(mergeRows(rows));
		}

		// Wait for all threads to finish
		for(CalibrationThread ct : calThreads){
			try{ct.join();}catch(InterruptedException e){Thread.currentThread().interrupt();}
			if(!ct.success){System.err.println("Warning: a calibration thread did not complete successfully.");}
		}

		// Aggregate per-thread histograms
		final long[] histCount=new long[numSlots];
		final long[] histTrueCard=new long[numSlots];
		final double[][] histRawEst=new double[numSlots][NUM_EST];
		for(CalibrationThread ct : calThreads){
			for(int s=0; s<numSlots; s++){
				histCount[s]+=ct.histCount[s];
				histTrueCard[s]+=ct.histTrueCard[s];
				for(int e=0; e<NUM_EST; e++){histRawEst[s][e]+=ct.histRawEst[s][e];}
			}
		}

		// Print File 2 (occupancy histogram)
		if(out2!=null){
			try(PrintStream ps=new PrintStream(new FileOutputStream(out2))){
				printHistogram(histCount, histTrueCard, histRawEst, numSlots, step, buckets, ps);
				System.err.println("Histogram written to "+out2);
			}catch(Exception e){
				System.err.println("Error writing histogram file: "+e.getMessage());
			}
		}else{
			System.err.println("Note: out2 not specified; occupancy histogram not written.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Threshold Logic        ----------------*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/*----------------          Factory             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a CardinalityTracker of the specified type.
	 * Recognized types: ddl, ddl2, ddl8, dll4 (and their full class name equivalents).
	 */
	static CardinalityTracker makeInstance(String type, int buckets, int k, long seed, float minProb){
		if("ddl".equals(type) || "dynamicDemiLog".equalsIgnoreCase(type)){
			return new DynamicDemiLog(buckets, k, seed, minProb);
		}else if("ddl2".equals(type) || "dynamicdemilog2".equalsIgnoreCase(type)){
			return new DynamicDemiLog2(buckets, k, seed, minProb);
		}else if("ddl8".equals(type)){
			return new DynamicDemiLog8(buckets, k, seed, minProb);
		}else if("dll4".equals(type) || "dynamicdemilog4".equalsIgnoreCase(type)){
			return new DynamicLogLog4(buckets, k, seed, minProb);
		}
		throw new RuntimeException("Unknown loglogtype: "+type);
	}

	/**
	 * Pre-computes logarithmically-spaced reporting thresholds from 1 to maxTrue.
	 * Each successive threshold is at least 1 greater and at least (1+reportFrac) times
	 * the previous, giving dense reporting at low cardinality and sparse at high.
	 * maxTrue is always the final threshold.
	 */
	static long[] computeThresholds(final long maxTrue, final double reportFrac){
		final ArrayList<Long> list=new ArrayList<>();
		long next=1;
		while(next<maxTrue){
			list.add(next);
			next=Math.max(next+1, (long)(next*(1+reportFrac)));
		}
		list.add(maxTrue);
		final long[] arr=new long[list.size()];
		for(int i=0; i<arr.length; i++){arr[i]=list.get(i);}
		return arr;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Row Output            ----------------*/
	/*--------------------------------------------------------------*/

	/** Merges N ReportRows (one per thread, all at the same threshold) into one combined row. */
	static ReportRow mergeRows(final ReportRow[] rows){
		final double[] sumErr=new double[NUM_OUT1];
		final double[] sumAbsErr=new double[NUM_OUT1];
		final double[] sumSqErr=new double[NUM_OUT1];
		long trueCard=rows[0].trueCard;
		int n=0;
		double occSum=0;
		for(ReportRow row : rows){
			n+=row.n;
			occSum+=row.occSum;
			for(int e=0; e<NUM_OUT1; e++){
				sumErr[e]+=row.sumErr[e];
				sumAbsErr[e]+=row.sumAbsErr[e];
				sumSqErr[e]+=row.sumSqErr[e];
			}
		}
		return new ReportRow(trueCard, n, occSum, sumErr, sumAbsErr, sumSqErr);
	}

	/** Prints one merged row in File 1 (interval) format to stdout. */
	static void printRow1(final ReportRow row){
		final int n=row.n;
		final double avgOcc=row.occSum/n;
		final StringBuilder sb=new StringBuilder();
		sb.append(row.trueCard).append('\t').append(String.format("%.5f", avgOcc));
		for(int e=0; e<NUM_OUT1; e++){
			final double meanErr=row.sumErr[e]/n;
			final double meanAbsErr=row.sumAbsErr[e]/n;
			// Variance via computational formula: E[X^2] - E[X]^2
			final double variance=row.sumSqErr[e]/n-meanErr*meanErr;
			final double stdev=Math.sqrt(Math.max(0, variance));
			sb.append('\t').append(String.format("%.6f", meanErr));
			sb.append('\t').append(String.format("%.6f", meanAbsErr));
			sb.append('\t').append(String.format("%.6f", stdev));
		}
		System.out.println(sb);
	}

	/** Prints the occupancy histogram (File 2) to the given stream. */
	static void printHistogram(final long[] histCount, final long[] histTrueCard,
		final double[][] histRawEst, final int numSlots, final int step,
		final int buckets, final PrintStream ps){
		ps.println(header2());
		for(int s=0; s<numSlots; s++){
			final long cnt=histCount[s];
			if(cnt==0){continue;}
			final StringBuilder sb=new StringBuilder();
			sb.append(s);
			final double avgTrueCard=(s==0 ? 0 : (double)histTrueCard[s]/cnt);
			for(int e=0; e<NUM_EST; e++){
				if(e==6){continue;} // Hybrid: no column
				// Slot 0, LC (e==5), and Hybrid pad all get 1.0
				if(s==0 || e==5){sb.append('\t').append("1.00000000");}
				else{
					final double avgRaw=histRawEst[s][e]/cnt;
					final double cf=computeCF(avgTrueCard, avgRaw);
					sb.append('\t').append(String.format("%.8f", cf));
				}
			}
			ps.println(sb);
		}
	}

	/** Computes correction factor: the multiplier that converts rawEstimate to trueCard.
	 *  Returns 1 for empty slot (both zero), NaN only when avgRaw is negative. */
	static double computeCF(final double avgTrueCard, final double avgRaw){
		if(avgTrueCard==0 && avgRaw==0){return 1.0;}// empty slot
		if(!Double.isFinite(avgRaw) || avgRaw<=0){return Double.NaN;}
		return avgTrueCard/avgRaw;
	}

	static String header1(){
		final StringBuilder sb=new StringBuilder("TrueCard\tOccupancy");
		for(String name : ESTIMATOR_NAMES){
			sb.append('\t').append(name).append("_err");
			sb.append('\t').append(name).append("_abs");
			sb.append('\t').append(name).append("_std");
		}
		sb.append('\t').append(MEDIAN_LC).append("_err");
		sb.append('\t').append(MEDIAN_LC).append("_abs");
		sb.append('\t').append(MEDIAN_LC).append("_std");
		return sb.toString();
	}

	static String header2(){
		// Writes one column per CorrectionFactor type constant (0..8), so matrix[type] works directly.
		// e=5 (LC) gets a Pad_cf column at type index 6 (LINEAR); e=6 (Hybrid) is skipped entirely.
		final StringBuilder sb=new StringBuilder("#Slot");
		for(int e=0; e<NUM_EST; e++){
			if(e==6){continue;} // Hybrid: no type constant, pre-corrected, skip
			if(e==5){sb.append("\tPad_cf");} // placeholder for LINEAR=6; always 1.0
			else{sb.append('\t').append(ESTIMATOR_NAMES[e]).append("_cf");}
		}
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Per-interval statistics row passed from a CalibrationThread to the main thread
	 * via an ArrayBlockingQueue.  Contains sufficient information to compute mean error,
	 * mean absolute error, and standard deviation across all DDLs.
	 */
	static final class ReportRow {
		/** True cardinality at the time of this report. */
		final long trueCard;
		/** Number of DDLs contributing to this row. */
		final int n;
		/** Sum of occupancy values across all DDLs. */
		final double occSum;
		/** Sum of signed relative errors per estimator. */
		final double[] sumErr;
		/** Sum of absolute relative errors per estimator. */
		final double[] sumAbsErr;
		/** Sum of squared relative errors per estimator (for stdev via E[X^2]-E[X]^2). */
		final double[] sumSqErr;

		ReportRow(long trueCard, int n, double occSum,
			double[] sumErr, double[] sumAbsErr, double[] sumSqErr){
			this.trueCard=trueCard;
			this.n=n;
			this.occSum=occSum;
			this.sumErr=sumErr;
			this.sumAbsErr=sumAbsErr;
			this.sumSqErr=sumSqErr;
		}
	}

	/**
	 * Independent simulation thread.  Maintains its own DDLs, value generator, and
	 * true-cardinality tracker.  Posts one ReportRow to its queue per reporting threshold.
	 * Also records raw estimator values into a 257-slot occupancy histogram at exact
	 * 1/256-occupancy checkpoints; the main thread aggregates these after all threads finish.
	 */
	static final class CalibrationThread extends Thread {

		CalibrationThread(CardinalityTracker[] ddls, long[] thresholds,
			long valSeed, int buckets, long maxTrue, int numSlots, int step){
			this.ddls=ddls;
			this.thresholds=thresholds;
			this.valSeed=valSeed;
			this.buckets=buckets;
			this.maxTrue=maxTrue;
			this.numSlots=numSlots;
			this.step=step;
			// Queue capacity = total thresholds: thread can fully finish before main thread consumes
			this.queue=new ArrayBlockingQueue<>(thresholds.length+1);
			histCount=new long[numSlots];
			histTrueCard=new long[numSlots];
			histRawEst=new double[numSlots][NUM_EST];
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

		void runInner() throws InterruptedException {

			// Record slot 0 before any elements: all DDLs are empty, trueCard=0
			for(CardinalityTracker ddl : ddls){
				histCount[0]++;
				// histTrueCard[0] += 0L; (zero, implicit)
				final double[] raw=ddl.rawEstimates();
				for(int e=0; e<NUM_EST; e++){histRawEst[0][e]+=raw[e];}
			}

			final LongHashSet trueSet=new LongHashSet(
				(int)Math.min((long)maxTrue+(long)maxTrue/4+1L, Integer.MAX_VALUE));
			final Random valRng=new Random(valSeed);

			// Per-DDL cached rawEstimates; recomputed only when lastCardinality==-1 (bucket changed).
			final double[][] cachedRaw=new double[ddls.length][];

			int ti=0;
			long prevTrueCard=0;

			while(ti<thresholds.length){
				final long val=valRng.nextLong();
				trueSet.add(val);
				final long trueCard=trueSet.size();
				for(CardinalityTracker ddl : ddls){ddl.hashAndStore(val);}

				// Record every DDL at its current occupancy slot for every unique element.
				// If the DDL's state didn't change (lastCardinality>=0), reuse cached estimates.
				// If it changed (lastCardinality==-1), recompute and mark current (set to 0).
				if(trueCard>prevTrueCard){
					for(int d=0; d<ddls.length; d++){
						if(ddls[d].getLastCardinality()<0 || cachedRaw[d]==null){
							cachedRaw[d]=ddls[d].rawEstimates();
							ddls[d].setLastCardinality(0);
						}
						final int idx=ddls[d].filledBuckets()/step;
						histCount[idx]++;
						histTrueCard[idx]+=trueCard;
						for(int e=0; e<NUM_EST; e++){histRawEst[idx][e]+=cachedRaw[d][e];}
					}
					prevTrueCard=trueCard;
				}

				// Reporting threshold: post merged stats row when reached
				if(trueCard>=thresholds[ti]){
					queue.put(buildReportRow(trueCard));
					ti++;
				}
			}
		}

		/** Computes statistics over all DDLs at the given trueCard and returns a ReportRow. */
		ReportRow buildReportRow(final long trueCard){
			final double[] sumErr=new double[NUM_OUT1];
			final double[] sumAbsErr=new double[NUM_OUT1];
			final double[] sumSqErr=new double[NUM_OUT1];
			double occSum=0;
			for(CardinalityTracker ddl : ddls){
				occSum+=ddl.occupancy();
				final double[] est=ddl.rawEstimates();
				for(int e=0; e<NUM_EST; e++){
					final double err=(est[e]-trueCard)/(double)trueCard;
					sumErr[e]+=err;
					sumAbsErr[e]+=Math.abs(err);
					sumSqErr[e]+=err*err;
				}
				// MedianLC: (LC + MedianCorr) / 2 (indices 5 and 8)
				final double mlc=(est[5]+est[8])*0.5;
				final double mlce=(mlc-trueCard)/(double)trueCard;
				sumErr[NUM_EST]+=mlce;
				sumAbsErr[NUM_EST]+=Math.abs(mlce);
				sumSqErr[NUM_EST]+=mlce*mlce;
			}
			return new ReportRow(trueCard, ddls.length, occSum, sumErr, sumAbsErr, sumSqErr);
		}

		/** CardinalityTracker instances owned by this thread. */
		final CardinalityTracker[] ddls;
		/** Pre-computed reporting thresholds (shared read-only with all threads). */
		final long[] thresholds;
		/** Seed for this thread's value generator. */
		final long valSeed;
		/** Bucket count (for histogram index computation). */
		final int buckets;
		/** Maximum true cardinality (= buckets * maxMult). */
		final long maxTrue;
		/** Number of histogram slots (= buckets/step + 1). */
		final int numSlots;
		/** Filled-bucket stride per histogram slot. */
		final int step;
		/** Queue for passing interval rows to the main thread. */
		final ArrayBlockingQueue<ReportRow> queue;

		/** Occupancy histogram: sample counts per slot. */
		final long[] histCount;
		/** Occupancy histogram: sum of true cardinalities per slot. */
		final long[] histTrueCard;
		/** Occupancy histogram: sum of raw estimator values per slot. [slot][estimator] */
		final double[][] histRawEst;

		/** Set to true after runInner() completes without exception. */
		volatile boolean success=false;
	}

}
