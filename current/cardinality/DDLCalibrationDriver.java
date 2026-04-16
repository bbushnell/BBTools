package cardinality;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteStreamWriter;
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
 * NOTE! This class is largely superceded by DDLCalibrationDriver2 which is more efficient.
 * DDLCalibrationDriver remains for memory bandwidth testing and some parsing duties. 
 *
 * @author Chloe
 * @date March 2026
 */
public class DDLCalibrationDriver {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of estimators reported by rawEstimates(). Includes MeanH and MeanM at end. */
	static final int NUM_EST=11+6+CardinalityStats.NUM_DLC_TIERS+AbstractCardStats.NUM_EXTRA+2;// +2 for WordEst, WordEstCV
	/** Estimator names in rawEstimates() index order. */
	static final String[] ESTIMATOR_NAMES;
	/** Which estimators get a CF column in File 2. LC and Micro excluded (no CF applied).
	 *  Hybrid included at column 10 (after main loop) to enable its own CF curve.
	 *  DLC columns have no CF. */
	static final boolean[] NEEDS_CF;
	static{
		final String[] base={"Mean","HMean","HMeanM","GMean","HLL","LC","Hybrid","HybDLC50","DThHyb","LCmin","DLCPure","DLC","DLC3B","DLCBest","HybDLC","SBS","SBSMult"};
		ESTIMATOR_NAMES=new String[NUM_EST];
		System.arraycopy(base, 0, ESTIMATOR_NAMES, 0, base.length);
		for(int i=0; i<CardinalityStats.NUM_DLC_TIERS; i++){ESTIMATOR_NAMES[base.length+i]="DLC"+i;}
		ESTIMATOR_NAMES[AbstractCardStats.MEANH_IDX]="MeanH";
		ESTIMATOR_NAMES[AbstractCardStats.MEANM_IDX]="MeanM";
		ESTIMATOR_NAMES[AbstractCardStats.LC_NOMICRO_IDX]="LC_noMicro";
		ESTIMATOR_NAMES[AbstractCardStats.SBS_NOMICRO_IDX]="SBS_noMicro";
		NEEDS_CF=new boolean[NUM_EST];
		//                      Mean  HMean HMnM  GMean HLL   LC    Hybr  HD50  DThH  LCmin DPure DLC
		final boolean[] baseCF={true, true, true, true, false,false,true, true, true, false,false,false};
		System.arraycopy(baseCF, 0, NEEDS_CF, 0, baseCF.length);
		// DLC entries (11..NUM_EST-1): all false (no CF)
	}
	/** Extra out1-only derived columns appended after NUM_EST. */
	static final int NUM_OUT1=NUM_EST+1; // +1 for MedianLC
	static final String MEDIAN_LC="MedianLC";
	/** When true, appends an extra RawMean column (avg uncorrected meanEst) to File 1 output. */
	static boolean OUTPUT_RAW_MEAN=false;
	/** CF table output version: 0 = legacy bipartite, 1 = unified DLC3B-indexed. */
	static int cfVersion=0;
	/** When true, asserts DLC estimate within 50% of true cardinality and dumps state on failure. */
	static boolean ASSERT_DLC=false;
	/** When false (default), suppresses per-tier DLC0..DLC63 and MedianLC columns from File 1.
	 *  These are useful for DLC algorithm development but waste space in production runs. */
	static boolean PRINT_DLC_TIERS=false;
	/** When true, runs pure hashAndStore throughput benchmark: no rawEstimates(), no histograms,
	 *  no filledBuckets(). Use for cache-thrashing speed comparisons between estimator types. */
	static boolean BENCHMARK_MODE=false;
	/** When false (default), suppresses the _std (standard deviation) column for each estimator.
	 *  Saves ~1/3 of output width; std is rarely needed when inspecting mean errors. */
	static boolean PRINT_STD=false;
	/** When true (default), prints the _cv (coefficient of variation = std/abs) column.
	 *  This magnitude-independent variance metric shows how tight the estimate is,
	 *  regardless of systematic bias (which CF tables correct). */
	static boolean PRINT_CV=true;

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		if(true){throw new RuntimeException("Use ddlcalibrate.sh.");}
		final long t0=System.nanoTime();
		int numDDLs=128;
		int buckets=2048;
		int k=31;
		long maxMult=10;
		double reportFrac=0.01;
		long masterSeed=12345L;
		long valSeed=42L;
		int threads=Shared.threads();
		int step=1;
		String out1="stdout.txt";
		String out2=null;
		String out3=null;
		String loglogtype="ddl"; // ddl, ddl2, ddl8, dll4

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
			else if(a.equals("valseed")){valSeed=Long.parseLong(b);}
			else if(a.equals("threads") || a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("step") || a.equals("resolution") || a.equals("res")){step=Integer.parseInt(b);}
			else if(a.equals("out") || a.equals("out1")){out1=b;}
			else if(a.equals("out2")){out2=b;}
			else if(a.equals("out3")){out3=b;}
			else if(a.equals("loglogtype") || a.equals("type")){loglogtype=b.toLowerCase();}
			else if(a.equals("cf") || a.equals("loglogcf")){CorrectionFactor.USE_CORRECTION=Parse.parseBoolean(b);}
			else if(a.equals("cardcf")){CardinalityTracker.USE_CARD_CF=Parse.parseBoolean(b);}
			else if(a.equals("rawmean")){OUTPUT_RAW_MEAN=Parse.parseBoolean(b);}
			else if(a.equals("cfversion")){cfVersion=Integer.parseInt(b);}
			else if(a.equals("cffile")){CorrectionFactor.initialize(b, buckets); CorrectionFactor.USE_CORRECTION=true;}
			else if(a.equals("dlcalpha") || a.equals("alpha")){CardinalityStats.DLC_ALPHA=Float.parseFloat(b);}
			else if(a.equals("cfiters") || a.equals("cfiterations")){CardinalityStats.DEFAULT_CF_ITERS=Integer.parseInt(b);}
			else if(a.equals("cfdif") || a.equals("cfconvergence")){CardinalityStats.DEFAULT_CF_DIF=Double.parseDouble(b);}
			else if(a.equals("tracecf")){CorrectionFactor.TRACE_CF=Parse.parseBoolean(b);}
			else if(a.equals("minvfraction") || a.equals("minvk")){
				float x=Float.parseFloat(b);
				if(x<1) {
					CardinalityStats.DLC_MIN_VK_FRACTION=x;
				}else {
					CardinalityStats.DLC_MIN_VK=(int)x;
				}
			}else if(a.equals("promotethreshold") || a.equals("pt")){
				DynamicLogLog3.PROMOTE_THRESHOLD=Integer.parseInt(b);
				DynamicLogLog4.PROMOTE_THRESHOLD=Integer.parseInt(b);
			}else if(a.equals("assertdlc")){
				ASSERT_DLC=Parse.parseBoolean(b);
			}else if(a.equals("benchmark") || a.equals("bench")){
				BENCHMARK_MODE=Parse.parseBoolean(b);
			}else if(a.equals("printcv") || a.equals("cv")){
				PRINT_CV=Parse.parseBoolean(b);
			}else if(a.equals("empiricalmantissa") || a.equals("em")){
				DynamicDemiLog8.USE_EMPIRICAL_MANTISSA=Parse.parseBoolean(b);
			}else if(a.equals("mantissacfoffset") || a.equals("mco")){
				DynamicDemiLog8.MANTISSA_CF_OFFSET=Double.parseDouble(b);
			}
			else if(a.equals("correctoverflow") || a.equals("co")){
				DynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
				BankedDynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
				CompressedDynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
				if(Parse.parseBoolean(b)){DynamicLogLog3.IGNORE_OVERFLOW=false; DynamicLogLog2.IGNORE_OVERFLOW=false; BankedDynamicLogLog3.IGNORE_OVERFLOW=false; CompressedDynamicLogLog3.IGNORE_OVERFLOW=false;}
			}else if(a.equals("ignoreoverflow") || a.equals("io")){
				DynamicLogLog3.IGNORE_OVERFLOW=Parse.parseBoolean(b);
				DynamicLogLog2.IGNORE_OVERFLOW=Parse.parseBoolean(b);
				BankedDynamicLogLog3.IGNORE_OVERFLOW=Parse.parseBoolean(b);
				DynamicLogLog4.IGNORE_OVERFLOW=Parse.parseBoolean(b);
				CompressedDynamicLogLog3.IGNORE_OVERFLOW=Parse.parseBoolean(b);
				if(Parse.parseBoolean(b)){DynamicLogLog3.CORRECT_OVERFLOW=false; BankedDynamicLogLog3.CORRECT_OVERFLOW=false; CompressedDynamicLogLog3.CORRECT_OVERFLOW=false;}
			}else if(a.equals("earlypromote") || a.equals("ep")){
				DynamicLogLog3.EARLY_PROMOTE=Parse.parseBoolean(b);
				DynamicLogLog4.EARLY_PROMOTE=Parse.parseBoolean(b);
				CompressedDynamicLogLog3.EARLY_PROMOTE=Parse.parseBoolean(b);
				DualHashDynamicLogLog4.EARLY_PROMOTE=Parse.parseBoolean(b);
			}else if(a.equals("promotefrac") || a.equals("pf")){
				final float pf=Float.parseFloat(b);
				BankedDynamicLogLog3.PROMOTE_FRAC=pf;
				BankedDynamicLogLog5.PROMOTE_FRAC=pf;
			}
			else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}

		assert(buckets%step==0) : "buckets must be a multiple of step; buckets="+buckets+", step="+step;
		final int numSlots=buckets/step+1;

		final long maxTrue=(long)buckets*maxMult;
		final int numThreads=Math.min(threads, numDDLs);
		final long[] thresholds=computeThresholds(maxTrue, reportFrac);
		v3ColsForType(loglogtype);

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

		// Open File 1 output (ByteStreamWriter: "stdout.txt" → System.out)
		final ByteStreamWriter bsw1=new ByteStreamWriter(out1, true, false, false);
		bsw1.start();
		bsw1.println(header1());

		// Start threads
		for(CalibrationThread ct : calThreads){ct.start();}

		// Benchmark mode: just wait for threads and print elapsed time
		if(BENCHMARK_MODE){
			for(CalibrationThread ct : calThreads){
				try{ct.join();}catch(InterruptedException e){throw new RuntimeException(e);}
			}
			final double elapsed=(System.nanoTime()-t0)*1e-9;
			System.err.println("Benchmark: "+loglogtype+"  Buckets: "+buckets+"  DDLs: "+numDDLs
				+"  MaxCard: "+maxTrue+"  Elapsed: "+String.format("%.3f", elapsed)+"s");
			bsw1.poisonAndWait();
			return;
		}

		// Main polling loop: collect one merged row per threshold, print live
		final double[] totalMeanAbsErr=new double[NUM_OUT1];
		final ArrayList<ReportRow> mergedRows=new ArrayList<>();
		int numRows=0;
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
			final ReportRow merged=mergeRows(rows);
			mergedRows.add(merged);
			bsw1.println(formatRow1(merged));
			for(int e=0; e<NUM_OUT1; e++){totalMeanAbsErr[e]+=merged.sumAbsErr[e]/merged.n;}
			numRows++;
		}
		bsw1.poisonAndWait();

		// Wait for all threads to finish
		for(CalibrationThread ct : calThreads){
			try{ct.join();}catch(InterruptedException e){Thread.currentThread().interrupt();}
			if(!ct.success){System.err.println("Warning: a calibration thread did not complete successfully.");}
		}

		// Print stderr summary
		{
			final double elapsed=(System.nanoTime()-t0)*1e-9;
			System.err.println();
			System.err.println("=== DDL Calibration Summary ===");
			System.err.println("Type: "+loglogtype+"  Buckets: "+buckets+"  DDLs: "+numDDLs
				+"  MaxCard: "+maxTrue+"  Rows: "+numRows+"  Elapsed: "+String.format("%.1f", elapsed)+"s");
			System.err.println("--- Avg Mean Absolute Error (lower = better) ---");
			// Key estimators only: indices 0-6, 11-14, then MedianLC
			final int[] keyIdx={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
			for(int ki=0; ki<keyIdx.length; ki++){
				final int e=keyIdx[ki];
				if(e>=NUM_EST){continue;}
				System.err.println(String.format("%-12s %.8f", ESTIMATOR_NAMES[e], totalMeanAbsErr[e]/numRows));
			}
			System.err.println(String.format("%-12s %.8f", MEDIAN_LC, totalMeanAbsErr[NUM_EST]/numRows));
			// Individual DLC tiers summary (only non-trivial ones)
			System.err.println("--- DLC Tier Detail (non-trivial) ---");
			for(int e=17; e<NUM_EST; e++){
				final double avg=totalMeanAbsErr[e]/numRows;
				if(avg<0.999){
					System.err.println(String.format("%-12s %.8f", ESTIMATOR_NAMES[e], avg));
				}
			}
			System.err.println();
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

		// Aggregate per-thread cardinality histograms
		final int numCardSlots=computeNumCardSlots(maxTrue, buckets);
		final long[] histCardCount=new long[numCardSlots];
		final long[] histCardTrueCard=new long[numCardSlots];
		final double[] histCardMeanEstSum=new double[numCardSlots];
		final double[][] histCardRawEst=new double[numCardSlots][NUM_EST];
		for(CalibrationThread ct : calThreads){
			for(int s=0; s<numCardSlots; s++){
				histCardCount[s]+=ct.histCardCount[s];
				histCardTrueCard[s]+=ct.histCardTrueCard[s];
				histCardMeanEstSum[s]+=ct.histCardMeanEstSum[s];
				for(int e=0; e<NUM_EST; e++){histCardRawEst[s][e]+=ct.histCardRawEst[s][e];}
			}
		}

		// Aggregate per-thread v1 histograms (if cfVersion 1 or 2)
		final int numV1Slots=(cfVersion>=1 && cfVersion<3 ? computeNumV1Slots(maxTrue, buckets) : 0);
		final long[] histV1Count=(numV1Slots>0 ? new long[numV1Slots] : null);
		final long[] histV1TrueCard=(numV1Slots>0 ? new long[numV1Slots] : null);
		final double[] histV1DLC3BSum=(numV1Slots>0 ? new double[numV1Slots] : null);
		final double[][] histV1RawEst=(numV1Slots>0 ? new double[numV1Slots][NUM_EST] : null);
		if(numV1Slots>0){
			for(CalibrationThread ct : calThreads){
				for(int s=0; s<numV1Slots; s++){
					histV1Count[s]+=ct.histV1Count[s];
					histV1TrueCard[s]+=ct.histV1TrueCard[s];
					histV1DLC3BSum[s]+=ct.histV1DLC3BSum[s];
					for(int e=0; e<NUM_EST; e++){histV1RawEst[s][e]+=ct.histV1RawEst[s][e];}
				}
			}
		}

		// Aggregate per-thread v3 histograms (if cfVersion>=3)
		final int numV3Slots=(cfVersion>=3 ? computeNumV3Slots(maxTrue) : 0);
		final long[] histV3Count=(numV3Slots>0 ? new long[numV3Slots] : null);
		final long[] histV3TrueCard=(numV3Slots>0 ? new long[numV3Slots] : null);
		final double[] histV3DLCSum=(numV3Slots>0 ? new double[numV3Slots] : null);
		final double[][] histV3RawEst=(numV3Slots>0 ? new double[numV3Slots][NUM_EST] : null);
		final double[][] histV3CfSum=(numV3Slots>0 ? new double[numV3Slots][NUM_EST] : null);
		if(numV3Slots>0){
			for(CalibrationThread ct : calThreads){
				for(int s=0; s<numV3Slots; s++){
					histV3Count[s]+=ct.histV3Count[s];
					histV3TrueCard[s]+=ct.histV3TrueCard[s];
					histV3DLCSum[s]+=ct.histV3DLCSum[s];
					for(int e=0; e<NUM_EST; e++){
						histV3RawEst[s][e]+=ct.histV3RawEst[s][e];
						histV3CfSum[s][e]+=ct.histV3CfSum[s][e];
					}
				}
			}
		}

		// Print File 2
		if(out2!=null){
			try(PrintStream ps=new PrintStream(new FileOutputStream(out2))){
				if(cfVersion>=3){
					printHistogramV3(mergedRows, ps);
				}else if(cfVersion>=1){
					printHistogramV1(histV1Count, histV1TrueCard, histV1DLC3BSum, histV1RawEst,
						numV1Slots, buckets, ps);
				}else{
					printHistogram(histCount, histTrueCard, histRawEst,
						histCardCount, histCardTrueCard, histCardMeanEstSum, histCardRawEst,
						numSlots, numCardSlots, step, buckets, ps);
				}
				System.err.println("Histogram written to "+out2);
			}catch(Exception e){
				System.err.println("Error writing histogram file: "+e.getMessage());
			}
		}else{
			System.err.println("Note: out2 not specified; occupancy histogram not written.");
		}

		// Print File 3: v3 CF table (always row-based, independent of cfVersion)
		if(out3!=null){
			try(PrintStream ps=new PrintStream(new FileOutputStream(out3))){
				printHistogramV3(mergedRows, ps);
				System.err.println("V3 CF table written to "+out3);
			}catch(Exception e){
				System.err.println("Error writing v3 CF file: "+e.getMessage());
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Threshold Logic        ----------------*/
	/*--------------------------------------------------------------*/

	/** Minimum load factor recorded in the cardinality histogram (slot 0 = this load). */
	static final double MIN_CARD_LOAD=0.2;

	/** Number of log-spaced cardinality slots from MIN_CARD_LOAD to maxTrue/buckets at 1% increments.
	 *  Load-factor indexed so the table is valid for any bucket count. */
	static int computeNumCardSlots(long maxTrue, int buckets){
		return (int)(Math.log((double)maxTrue/(buckets*MIN_CARD_LOAD))/Math.log(1.01))+5;
	}

	/**
	 * Returns the cardinality histogram slot index for rawMeanEst.
	 * Slot 0 corresponds to loadFactor=MIN_CARD_LOAD; each slot is 1% larger.
	 * Load-factor indexed so the table is valid for any bucket count.
	 * Returns -1 if rawMeanEst is below MIN_CARD_LOAD*buckets (out of range).
	 */
	static int cardIdx(double rawMeanEst, int numCardSlots, int buckets){
		if(rawMeanEst<buckets*MIN_CARD_LOAD){return -1;}
		final int idx=(int)(Math.log(rawMeanEst/(buckets*MIN_CARD_LOAD))/Math.log(1.01));
		if(idx<0||idx>=numCardSlots){return -1;}
		return idx;
	}

	/*--------------------------------------------------------------*/
	/*----------------       v1 Slot Logic          ----------------*/
	/*--------------------------------------------------------------*/

	/** DLC3B index for v1 rawEstimates() array. */
	static final int DLC3B_IDX=12;

	/** Number of v1 histogram slots: integer slots 0..buckets, then 1% log-spaced above. */
	static int computeNumV1Slots(long maxTrue, int buckets){
		if(maxTrue<=buckets){return buckets+1;}
		return buckets+1+(int)(Math.log((double)maxTrue/buckets)/Math.log(1.01))+5;
	}

	/**
	 * Returns the v1 histogram slot for a raw DLC3B estimate.
	 * Slots 0..buckets: one per integer DLC3B value (slot = (int)dlc3b, clamped to [0,buckets]).
	 * Slots above buckets: 1% log-spaced (slot = buckets + log(dlc3b/buckets)/log(1.01)).
	 * Returns -1 if dlc3b < 1.
	 */
	static int v1Idx(double dlc3b, int numV1Slots, int buckets){
		if(dlc3b<1){return -1;}
		if(dlc3b<=buckets){return Math.min((int)dlc3b, buckets);}
		final int idx=buckets+(int)(Math.log(dlc3b/buckets)/Math.log(1.01));
		if(idx<0||idx>=numV1Slots){return -1;}
		return idx;
	}

	/** Returns the representative DLC3B value for a v1 slot index. */
	static double v1SlotToKey(int slot, int buckets){
		if(slot<=buckets){return slot;}
		return buckets*Math.pow(1.01, slot-buckets);
	}

	/*--------------------------------------------------------------*/
	/*----------------       v3 Slot Logic          ----------------*/
	/*--------------------------------------------------------------*/

	/** DLC (dlcLogSpace025) index in rawEstimates() array. */
	static final int DLC_IDX=11;

	/** Base v3 CF column definitions: rawEstimates() index → output column name.
	 *  Extended per-class by v3ColsForType() to include MeanH_cf, MeanM_cf as needed. */
	static final int[] V3_BASE_INDICES={0, 1, 2, 3, 4, 6, 11, 13, 12, 8}; // Mean,HMean,HMeanM,GMean,HLL,Hybrid,DLC,DLCBest,DLC3B,DThHyb
	static final String[] V3_BASE_NAMES={"Mean_cf","HMean_cf","HMeanM_cf","GMean_cf","HLL_cf","Hybrid_cf","DLC_cf","DLCBest_cf","DLC3B_cf","DThHyb_cf"};
	/** Active CF column definitions — set per-class by v3ColsForType(). */
	static int[] V3_EST_INDICES=V3_BASE_INDICES;
	static String[] V3_COL_NAMES=V3_BASE_NAMES;

	/** Configures V3_COL_NAMES/V3_EST_INDICES based on the DDL class type.
	 *  Classes with history get MeanH_cf; classes with mantissa get MeanM_cf.
	 *  No wasted columns: only outputs columns the class actually uses. */
	/** Index of WordEst in rawEstimates() when appended by DynamicLogLog4.
	 *  = 17 + NUM_DLC_TIERS + NUM_EXTRA = 85. */
	static final int WORDEST_RAW_IDX=17+AbstractCardStats.NUM_DLC_TIERS+AbstractCardStats.NUM_EXTRA;

	static void v3ColsForType(String type){
		final boolean hasHistory=type.equals("udll6") || type.equals("pll16c") || type.equals("ttll") || type.equals("bdll4") || type.equals("bdll5") || type.equals("cdll4") || type.equals("cdll5");
		final boolean hasMantissa=type.equals("ddl") || type.equals("ddl10") || type.equals("ddl8")
			|| type.equals("ddl8v2") || type.equals("ddl2");
		final boolean hasWordEst=type.equals("dll4") || type.equals("dll4m");
		final boolean hasMean16=type.equals("ttll");
		final int extra=(hasHistory ? 2 : 0)+(hasMantissa ? 1 : 0)+(hasWordEst ? 1 : 0)+(hasMean16 ? 1 : 0);
		if(extra==0){
			V3_EST_INDICES=V3_BASE_INDICES;
			V3_COL_NAMES=V3_BASE_NAMES;
			return;
		}
		final int base=V3_BASE_INDICES.length;
		V3_EST_INDICES=new int[base+extra];
		V3_COL_NAMES=new String[base+extra];
		System.arraycopy(V3_BASE_INDICES, 0, V3_EST_INDICES, 0, base);
		System.arraycopy(V3_BASE_NAMES, 0, V3_COL_NAMES, 0, base);
		int idx=base;
		if(hasHistory){V3_EST_INDICES[idx]=AbstractCardStats.MEANH_IDX; V3_COL_NAMES[idx]="MeanH_cf"; idx++;
			V3_EST_INDICES[idx]=AbstractCardStats.HC_IDX; V3_COL_NAMES[idx]="HC_cf"; idx++;}
		if(hasMantissa){V3_EST_INDICES[idx]=AbstractCardStats.MEANM_IDX; V3_COL_NAMES[idx]="MeanM_cf"; idx++;}
		if(hasWordEst){V3_EST_INDICES[idx]=WORDEST_RAW_IDX; V3_COL_NAMES[idx]="WordEst_cf"; idx++;}
		if(hasMean16){V3_EST_INDICES[idx]=TwinTailLogLog.MEAN16_RAW_IDX; V3_COL_NAMES[idx]="Mean16_cf"; idx++;}
	}

	/** Number of v3 histogram slots: integer slots 1..100, then 1% log-spaced above. */
	static int computeNumV3Slots(long maxTrue){
		if(maxTrue<=100){return (int)maxTrue+1;}
		return 101+(int)(Math.log((double)maxTrue/100)/Math.log(1.01))+5;
	}

	/**
	 * Returns the v3 histogram slot for a raw DLC estimate.
	 * Slots 1..100: one per integer DLC value.
	 * Slots above 100: 1% log-spaced (slot = 100 + log(dlc/100)/log(1.01)).
	 * Returns -1 if dlc < 1.
	 */
	static int v3Idx(double dlc, int numV3Slots){
		if(dlc<1){return -1;}
		if(dlc<=100){return (int)dlc;}
		final int idx=100+(int)(Math.log(dlc/100)/Math.log(1.01));
		if(idx<0||idx>=numV3Slots){return -1;}
		return idx;
	}

	/** Returns the representative DLC value for a v3 slot index. */
	static double v3SlotToKey(int slot){
		if(slot<=100){return slot;}
		return 100*Math.pow(1.01, slot-100);
	}

	/**
	 * Prints v5 CF table: identical data to v4 but with structured metadata header.
	 * Header lines are tab-delimited key-value pairs (#Key\tValue).
	 * Includes bucket count so the loader can scale lookup keys when bucket counts differ.
	 * CF = 1/(1+avgErr) where avgErr = sumErr/n for each estimator.
	 * Keyed by trueCard (absolute cardinality); at lookup time the key is scaled by
	 * (tableBuckets / currentBuckets) to normalize across bucket counts.
	 */
	static void printHistogramV5(final ArrayList<ReportRow> rows, final PrintStream ps,
			String loglogtype, int buckets, int numDdls, long maxTrue,
			boolean cfEnabled, boolean earlyPromote, int promoteThreshold, String notes){
		ps.println("#Version\t5");
		ps.println("#Type\t"+loglogtype);
		ps.println("#Buckets\t"+buckets);
		ps.println("#Estimators\t"+numDdls);
		ps.println("#MaxCardinality\t"+maxTrue);
		ps.println("#CF\t"+cfEnabled);
		ps.println("#clampToAdded\t"+CardinalityTracker.clampToAdded);
		ps.println("#EarlyPromote\t"+earlyPromote);
		ps.println("#PROMOTE_THRESHOLD\t"+promoteThreshold);
		if(notes!=null && !notes.isEmpty()){ps.println("#Notes\t"+notes);}
		ps.println(terminalRow(rows));
		final StringBuilder hdr=new StringBuilder("#TrueCard");
		for(String name : V3_COL_NAMES){hdr.append('\t').append(name);}
		ps.println(hdr);
		printHistogramRows(rows, ps);
	}

	/** Prints v3 CF table directly from per-trueCard report rows.
	 *  CF = 1/(1+avgErr) where avgErr = sumErr/n for each estimator.
	 *  Keyed by trueCard; DLC is used as an approximation at lookup time.
	 *  This avoids all histogram binning bias (Jensen's inequality + selection bias). */
	static void printHistogramV3(final ArrayList<ReportRow> rows, final PrintStream ps){
		ps.println("#VERSION=4");
		ps.println("#ESTIMATOR=DLC");
		final StringBuilder hdr=new StringBuilder("#TrueCard");
		for(String name : V3_COL_NAMES){hdr.append('\t').append(name);}
		ps.println(hdr);
		printHistogramRows(rows, ps);
	}

	/** Builds the "#terminal" header row: per-column average CF over the last octave
	 *  (rows where trueCard >= maxTrueCard/2). Used by downstream code to read the
	 *  asymptotic per-column bias directly from the CF file header. */
	private static String terminalRow(final ArrayList<ReportRow> rows){
		long maxCard=0;
		for(final ReportRow r : rows){if(r.n>=1 && r.trueCard>maxCard){maxCard=r.trueCard;}}
		final long threshold=maxCard/2;
		final int ncols=V3_EST_INDICES.length;
		final double[] sums=new double[ncols];
		final int[] counts=new int[ncols];
		for(final ReportRow r : rows){
			if(r.trueCard<threshold || r.n<1){continue;}
			for(int c=0; c<ncols; c++){
				final double avgErr=r.sumErr[V3_EST_INDICES[c]]/r.n;
				final double cf=(avgErr<=-0.99) ? 1.0 : 1.0/(1.0+avgErr);
				if(Double.isFinite(cf)){sums[c]+=cf; counts[c]++;}
			}
		}
		final StringBuilder sb=new StringBuilder("#terminal");
		for(int c=0; c<ncols; c++){
			final double avg=counts[c]>0 ? sums[c]/counts[c] : Double.NaN;
			sb.append('\t').append(String.format("%.6f", avg));
		}
		return sb.toString();
	}

	/** Shared row-printing logic for v3/v4/v5 CF tables. */
	private static void printHistogramRows(final ArrayList<ReportRow> rows, final PrintStream ps){
		for(final ReportRow row : rows){
			if(row.trueCard<1 || row.n<1){continue;}
			final StringBuilder sb=new StringBuilder();
			sb.append(row.trueCard);
			for(int c=0; c<V3_EST_INDICES.length; c++){
				final int eIdx=V3_EST_INDICES[c];
				final double avgErr=row.sumErr[eIdx]/row.n;
				final double cf=(avgErr<=-0.99) ? 1.0 : 1.0/(1.0+avgErr);
				sb.append('\t').append(Double.isFinite(cf) ?
					String.format("%.8f", cf) : "1.00000000");
			}
			ps.println(sb);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Factory             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a CardinalityTracker of the specified type.
	 * Delegates to CardinalityTracker.makeTracker for the canonical factory.
	 */
	static CardinalityTracker makeInstance(String type, int buckets, int k, long seed, float minProb){
		return CardinalityTracker.makeTracker(type, buckets, k, seed, minProb);
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
		double sumRawMean=0;
		for(ReportRow row : rows){
			n+=row.n;
			occSum+=row.occSum;
			sumRawMean+=row.sumRawMean;
			for(int e=0; e<NUM_OUT1; e++){
				sumErr[e]+=row.sumErr[e];
				sumAbsErr[e]+=row.sumAbsErr[e];
				sumSqErr[e]+=row.sumSqErr[e];
			}
		}
		final ReportRow merged=new ReportRow(trueCard, n, occSum, sumErr, sumAbsErr, sumSqErr);
		merged.sumRawMean=sumRawMean;
		return merged;
	}

	/** Formats one merged row in File 1 (interval) format. */
	static String formatRow1(final ReportRow row){
		final int n=row.n;
		final double avgOcc=row.occSum/n;
		final StringBuilder sb=new StringBuilder();
		sb.append(row.trueCard).append('\t').append(String.format("%.5f", avgOcc)).append('\t').append(n);
		for(int e=0; e<NUM_OUT1; e++){
			if(!PRINT_DLC_TIERS && e>=17){break;}
			final double meanErr=row.sumErr[e]/n;
			final double meanAbsErr=row.sumAbsErr[e]/n;
			sb.append('\t').append(String.format("%.6f", meanErr));
			sb.append('\t').append(String.format("%.6f", meanAbsErr));
			if(PRINT_STD || PRINT_CV){
				final double variance=row.sumSqErr[e]/n-meanErr*meanErr;
				final double stdev=Math.sqrt(Math.max(0, variance));
				if(PRINT_STD){sb.append('\t').append(String.format("%.6f", stdev));}
				if(PRINT_CV){
					// CV = std / |1+meanErr|: coefficient of variation of the raw estimates,
					// independent of CF scaling. Measures relative spread around the
					// (possibly biased) mean, so comparing CV across estimators shows
					// which has tighter estimates regardless of bias correction.
					final double denom=Math.abs(1.0+meanErr);
					sb.append('\t').append(String.format("%.6f", denom>0 ? stdev/denom : 0));
				}
			}
		}
		if(OUTPUT_RAW_MEAN){sb.append('\t').append(String.format("%.4f", row.sumRawMean/n));}
		return sb.toString();
	}

	/** Prints the bipartite occupancy+cardinality histogram (File 2) to the given stream. */
	static void printHistogram(final long[] histCount, final long[] histTrueCard,
		final double[][] histRawEst,
		final long[] histCardCount, final long[] histCardTrueCard,
		final double[] histCardMeanEstSum, final double[][] histCardRawEst,
		final int numSlots, final int numCardSlots, final int step,
		final int buckets, final PrintStream ps){
		ps.println("#SECTION=occupancy");
		ps.println(header2());
		for(int s=0; s<numSlots; s++){
			final long cnt=histCount[s];
			if(cnt==0){continue;}
			final StringBuilder sb=new StringBuilder();
			sb.append(s);
			final double avgTrueCard=(s==0 ? 0 : (double)histTrueCard[s]/cnt);
			for(int e=0; e<NUM_EST; e++){
				if(e==6){continue;} // Hybrid appended after loop at column 10
				// Slot 0, LC (e==5), and Micro (e==10) always get 1.0
				if(s==0 || e==5 || e==10){sb.append('\t').append("1.00000000");}
				else{
					final double avgRaw=histRawEst[s][e]/cnt;
					final double cf=computeCF(avgTrueCard, avgRaw);
					sb.append('\t').append(String.format("%.8f", cf));
				}
			}
			// Hybrid_cf at column 10 (type HYBRID=10)
			if(s==0){sb.append('\t').append("1.00000000");}
			else{
				final double avgRaw=histRawEst[s][6]/cnt;
				final double cf=computeCF(avgTrueCard, avgRaw);
				sb.append('\t').append(String.format("%.8f", cf));
			}
			ps.println(sb);
		}
		// Cardinality-indexed section: log-spaced rows grouped by rawMeanEst.
		if(histCardCount!=null && numCardSlots>0){
			// Truncate at the peak slot. Post-peak slots are biased: they contain only
			// overestimating DDLs (by definition of having meanEst above the modal value),
			// so CF = trueCard/meanEst is systematically low there. getCardCF() clamps to
			// the last table entry for any rawMeanEst beyond the table, so the peak CF
			// naturally extends to arbitrarily high cardinality — which is correct.
			int lastPrintSlot=0;
			long peakCount=0;
			for(int cIdx=0; cIdx<numCardSlots; cIdx++){
				if(histCardCount[cIdx]>peakCount){peakCount=histCardCount[cIdx]; lastPrintSlot=cIdx;}
			}
			lastPrintSlot=Math.max(0, lastPrintSlot-5); // 5-slot safety margin before peak
			ps.println();
			ps.println("#SECTION=cardinality");
			ps.println(header2card());
			for(int cIdx=0; cIdx<=lastPrintSlot; cIdx++){
				final long cnt=histCardCount[cIdx];
				if(cnt<10){continue;} // skip very sparse early rows
				final double avgMeanEst=histCardMeanEstSum[cIdx]/cnt;
				final double avgTrueCard=(double)histCardTrueCard[cIdx]/cnt;
				final StringBuilder sb=new StringBuilder();
				sb.append(String.format("%.8f", avgMeanEst/buckets));
				sb.append('\t').append(cnt); // Samples column
				for(int e=0; e<NUM_EST; e++){
					if(e==6){continue;} // Hybrid appended after
					if(e==5||e==10){sb.append('\t').append("1.00000000");} // LC/Micro: no CF
					else{
						final double avgRaw=histCardRawEst[cIdx][e]/cnt;
						final double cf=computeCF(avgTrueCard, avgRaw);
						sb.append('\t').append(Double.isFinite(cf) ?
							String.format("%.8f", cf) : "1.00000000");
					}
				}
				final double avgRawHybrid=histCardRawEst[cIdx][6]/cnt;
				final double cfHybrid=computeCF(avgTrueCard, avgRawHybrid);
				sb.append('\t').append(Double.isFinite(cfHybrid) ?
					String.format("%.8f", cfHybrid) : "1.00000000");
				ps.println(sb);
			}
		}
	}

	/** v1 CF column definitions: rawEstimates() index → output column name. */
	static final int[] V1_EST_INDICES={0, 1, 2, 3, 12, 4, 6, 14}; // Mean,HMean,HMeanM,GMean,DLC3B,HLL,Hybrid,HybDLC
	static final String[] V1_COL_NAMES={"Mean_cf","HMean_cf","HMeanM_cf","GMean_cf","DLC3B_cf","HLL_cf","Hybrid_cf","HybDLC_cf"};

	/** Prints v1 CF table: unified, cardinality-indexed by raw DLC3B estimate. */
	static void printHistogramV1(final long[] histV1Count, final long[] histV1TrueCard,
		final double[] histV1DLC3BSum, final double[][] histV1RawEst,
		final int numV1Slots, final int buckets, final PrintStream ps){
		ps.println("#VERSION=1");
		ps.println("#ESTIMATOR=DLC3B");
		// Header
		final StringBuilder hdr=new StringBuilder("#DLC3B_est");
		for(String name : V1_COL_NAMES){hdr.append('\t').append(name);}
		ps.println(hdr);
		// Truncate at peak slot (same logic as v0 cardinality section)
		int lastPrintSlot=0;
		long peakCount=0;
		for(int s=0; s<numV1Slots; s++){
			if(histV1Count[s]>peakCount){peakCount=histV1Count[s]; lastPrintSlot=s;}
		}
		lastPrintSlot=Math.max(0, lastPrintSlot-5);
		// Data rows
		for(int s=1; s<=lastPrintSlot; s++){
			final long cnt=histV1Count[s];
			if(cnt<10){continue;}
			final double avgDLC3B=histV1DLC3BSum[s]/cnt;
			final double avgTrueCard=(double)histV1TrueCard[s]/cnt;
			final StringBuilder sb=new StringBuilder();
			sb.append(String.format("%.2f", avgDLC3B));
			for(int c=0; c<V1_EST_INDICES.length; c++){
				final int eIdx=V1_EST_INDICES[c];
				final double avgRaw=histV1RawEst[s][eIdx]/cnt;
				final double cf=computeCF(avgTrueCard, avgRaw);
				sb.append('\t').append(Double.isFinite(cf) ?
					String.format("%.8f", cf) : "1.00000000");
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
		final StringBuilder sb=new StringBuilder("TrueCard\tOccupancy\tCount");
		for(int e=0; e<NUM_EST; e++){
			if(!PRINT_DLC_TIERS && e>=17){break;}
			final String name=ESTIMATOR_NAMES[e];
			sb.append('\t').append(name).append("_err");
			sb.append('\t').append(name).append("_abs");
			if(PRINT_STD){sb.append('\t').append(name).append("_std");}
			if(PRINT_CV){sb.append('\t').append(name).append("_cv");}
		}
		if(PRINT_DLC_TIERS){
			sb.append('\t').append(MEDIAN_LC).append("_err");
			sb.append('\t').append(MEDIAN_LC).append("_abs");
			if(PRINT_STD){sb.append('\t').append(MEDIAN_LC).append("_std");}
			if(PRINT_CV){sb.append('\t').append(MEDIAN_LC).append("_cv");}
		}
		if(OUTPUT_RAW_MEAN){sb.append('\t').append("RawMean");}
		return sb.toString();
	}

	static String header2(){
		// Writes one column per CorrectionFactor type constant (0..9), so matrix[type] works directly.
		// e=5 (LC) gets a Pad_cf column at type index 6 (LINEAR); e=6 (Hybrid) is skipped in loop.
		// Hybrid_cf is appended at the end as column 10 (type HYBRID=10).
		final StringBuilder sb=new StringBuilder("#Slot");
		for(int e=0; e<NUM_EST; e++){
			if(e==6){continue;} // Hybrid appended after loop at column 10
			if(e==5||e==10){sb.append("\tPad_cf");} // LC/Micro: no CF, placeholder
			else{sb.append('\t').append(ESTIMATOR_NAMES[e]).append("_cf");}
		}
		sb.append("\tHybrid_cf"); // column 10, type HYBRID=10
		return sb.toString();
	}

	/** Header line for the #SECTION=cardinality block of File 2. Same columns as header2() but #MeanEst first. */
	static String header2card(){
		final StringBuilder sb=new StringBuilder("#MeanEst\tSamples");
		for(int e=0; e<NUM_EST; e++){
			if(e==6){continue;}
			if(e==5||e==10){sb.append("\tPad_cf");}
			else{sb.append('\t').append(ESTIMATOR_NAMES[e]).append("_cf");}
		}
		sb.append("\tHybrid_cf");
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
		/** Sum of raw (uncorrected) Mean estimates across all DDLs; for OUTPUT_RAW_MEAN. */
		double sumRawMean=0;

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
			numCardSlots=computeNumCardSlots(maxTrue, buckets);
			histCardCount=new long[numCardSlots];
			histCardTrueCard=new long[numCardSlots];
			histCardMeanEstSum=new double[numCardSlots];
			histCardRawEst=new double[numCardSlots][NUM_EST];
			// v1 histogram (only when cfVersion>=1)
			if(cfVersion>=1 && cfVersion<3){
				numV1Slots=computeNumV1Slots(maxTrue, buckets);
				histV1Count=new long[numV1Slots];
				histV1TrueCard=new long[numV1Slots];
				histV1DLC3BSum=new double[numV1Slots];
				histV1RawEst=new double[numV1Slots][NUM_EST];
			}else{
				numV1Slots=0;
				histV1Count=null;
				histV1TrueCard=null;
				histV1DLC3BSum=null;
				histV1RawEst=null;
			}
			// v3 histogram (only when cfVersion>=3)
			if(cfVersion>=3){
				numV3Slots=computeNumV3Slots(maxTrue);
				histV3Count=new long[numV3Slots];
				histV3TrueCard=new long[numV3Slots];
				histV3DLCSum=new double[numV3Slots];
				histV3RawEst=new double[numV3Slots][NUM_EST];
				histV3CfSum=new double[numV3Slots][NUM_EST];
			}else{
				numV3Slots=0;
				histV3Count=null;
				histV3TrueCard=null;
				histV3DLCSum=null;
				histV3RawEst=null;
				histV3CfSum=null;
			}
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

			// Benchmark mode: pure add() throughput, no estimates or histograms.
			if(BENCHMARK_MODE){
				final Random valRng=new Random(valSeed);
				for(long trueCard=1; trueCard<=maxTrue; trueCard++){
					final long val=valRng.nextLong();
					for(CardinalityTracker ddl : ddls){ddl.add(val);}
				}
				return;
			}

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
				for(CardinalityTracker ddl : ddls){ddl.add(val);}

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
						// Also record into cardinality histogram (indexed by raw Mean estimate).
						final double rawMeanEst=cachedRaw[d][0]; // uncorrected when cf=f
						final int cIdx=cardIdx(rawMeanEst, numCardSlots, buckets);
						if(cIdx>=0){
							histCardCount[cIdx]++;
							histCardTrueCard[cIdx]+=trueCard;
							histCardMeanEstSum[cIdx]+=rawMeanEst;
							for(int e=0; e<NUM_EST; e++){histCardRawEst[cIdx][e]+=cachedRaw[d][e];}
						}
						// v1 histogram: indexed by raw DLC3B estimate
						if(histV1Count!=null){
							final double dlc3b=cachedRaw[d][DLC3B_IDX];
							final int vIdx=v1Idx(dlc3b, numV1Slots, buckets);
							if(vIdx>=0){
								histV1Count[vIdx]++;
								histV1TrueCard[vIdx]+=trueCard;
								histV1DLC3BSum[vIdx]+=dlc3b;
								for(int e=0; e<NUM_EST; e++){histV1RawEst[vIdx][e]+=cachedRaw[d][e];}
							}
						}
					}
					// v3 histogram: use per-trueCard AVERAGES to avoid selection bias.
					// Average across all DDLs first, then insert one point into histogram.
					if(histV3Count!=null){
						double avgDLC=0;
						final double[] avgRaw=new double[NUM_EST];
						int validCount=0;
						for(int d=0; d<ddls.length; d++){
							if(cachedRaw[d]!=null){
								avgDLC+=cachedRaw[d][DLC_IDX];
								for(int e=0; e<NUM_EST; e++){avgRaw[e]+=cachedRaw[d][e];}
								validCount++;
							}
						}
						if(validCount>0){
							avgDLC/=validCount;
							for(int e=0; e<NUM_EST; e++){avgRaw[e]/=validCount;}
							final int v3i=v3Idx(avgDLC, numV3Slots);
							if(v3i>=0){
								histV3Count[v3i]++;
								histV3TrueCard[v3i]+=trueCard;
								histV3DLCSum[v3i]+=avgDLC;
								for(int e=0; e<NUM_EST; e++){
									histV3RawEst[v3i][e]+=avgRaw[e];
									if(avgRaw[e]>0){histV3CfSum[v3i][e]+=trueCard/avgRaw[e];}
								}
							}
						}
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
			double sumRawMean=0;
			for(CardinalityTracker ddl : ddls){
				occSum+=ddl.occupancy();
				final double[] est=ddl.rawEstimates();
				// DLC assertion: uncomment to debug DLC divergence.
				// When enabled, dumps CardinalityStats.toString() at the first >50% error.
//				if(ASSERT_DLC && trueCard>0){
//					final double dlcEst=est[DLC_IDX];
//					final double relErr=Math.abs(dlcEst-trueCard)/(double)trueCard;
//					if(relErr>0.5){
//						System.err.println("\n=== DLC ASSERTION FAILURE ===");
//						System.err.println("trueCard="+trueCard+" dlcEst="+String.format("%.2f", dlcEst)
//							+" relErr="+String.format("%.4f", relErr));
//						if(ddl instanceof DynamicLogLog3v2){
//							final DynamicLogLog3v2 d=(DynamicLogLog3v2)ddl;
//							System.err.println("minZeros="+d.minZeros());
//							if(d.lastStats!=null){System.err.println(d.lastStats.toString());}
//						}
//						System.err.println("filledBuckets="+ddl.filledBuckets());
//						assert false : "DLC error "+String.format("%.1f%%", relErr*100)+" at trueCard="+trueCard;
//					}
//				}
				sumRawMean+=est[0];
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
			final ReportRow row=new ReportRow(trueCard, ddls.length, occSum, sumErr, sumAbsErr, sumSqErr);
			row.sumRawMean=sumRawMean;
			return row;
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

		/** Number of log-spaced cardinality histogram slots. */
		final int numCardSlots;
		/** Cardinality histogram: sample counts per slot. */
		final long[] histCardCount;
		/** Cardinality histogram: sum of true cardinalities per slot. */
		final long[] histCardTrueCard;
		/** Cardinality histogram: sum of raw Mean estimates per slot (for avg MeanEst key). */
		final double[] histCardMeanEstSum;
		/** Cardinality histogram: sum of raw estimator values per slot. [slot][estimator] */
		final double[][] histCardRawEst;

		/** v1 histogram: number of DLC3B-indexed slots. 0 when cfVersion<1. */
		final int numV1Slots;
		/** v1 histogram: sample counts per slot. Null when cfVersion<1. */
		final long[] histV1Count;
		/** v1 histogram: sum of true cardinalities per slot. */
		final long[] histV1TrueCard;
		/** v1 histogram: sum of raw DLC3B estimates per slot (for avg key). */
		final double[] histV1DLC3BSum;
		/** v1 histogram: sum of raw estimator values per slot. */
		final double[][] histV1RawEst;

		/** v3 histogram: number of DLC-indexed slots. 0 when cfVersion<3. */
		final int numV3Slots;
		/** v3 histogram: sample counts per slot. Null when cfVersion<3. */
		final long[] histV3Count;
		/** v3 histogram: sum of true cardinalities per slot. */
		final long[] histV3TrueCard;
		/** v3 histogram: sum of raw DLC estimates per slot (for avg key). */
		final double[] histV3DLCSum;
		/** v3 histogram: sum of raw estimator values per slot. */
		final double[][] histV3RawEst;
		/** v3 histogram: sum of per-sample CF (trueCard/rawEst) per slot per estimator.
		 *  Used to compute average-of-ratios instead of ratio-of-averages (Jensen fix). */
		final double[][] histV3CfSum;

		/** Set to true after runInner() completes without exception. */
		volatile boolean success=false;
	}

}
