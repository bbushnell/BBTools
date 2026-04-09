package cardinality;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.PreParser;
import rand.FastRandomXoshiro;
import shared.Shared;
import shared.Tools;

/**
 * Cache-friendly multithreaded calibration driver for DynamicDemiLog cardinality estimators.
 * <p>
 * Each thread processes one DDL at a time: create, feed all values sequentially,
 * record measurements at reporting thresholds, discard, repeat. This keeps the DDL's
 * working set in L1/L2 cache throughout its lifetime.
 * <p>
 * Follows BBTools MT template: constructor wraps global state changes in synchronized(class),
 * threads take defensive local copies of statics at initialization, accumulator merge
 * synchronizes on each thread before reading.
 *
 * @author Chloe, Brian Bushnell
 * @date April 2026
 */
public class DDLCalibrationDriver2 {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static boolean CLAMP_TO_ADDED=false;
	static final int NUM_EST=DDLCalibrationDriver.NUM_EST;
	static final String[] ESTIMATOR_NAMES=DDLCalibrationDriver.ESTIMATOR_NAMES;
	static final int NUM_LDLC=9;
	static final String[] LDLC_NAMES={"LDLC", "DLC_L", "HC", "FGRA", "HLL+H", "Mean+H", "Hybrid+2", "LC_noMicro", "SBS_noMicro"};

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final long t0=System.nanoTime();
		DDLCalibrationDriver2 driver=new DDLCalibrationDriver2(args);
		driver.process(t0);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	DDLCalibrationDriver2(String[] args){
		synchronized(DDLCalibrationDriver2.class){
			{
				PreParser pp=new PreParser(args, null, false);
				args=pp.args;
			}
			parse(args);
			initGlobalState();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	private void parse(String[] args){
		for(String arg : args){
			final String[] split=arg.split("=");
			if(split.length!=2){throw new RuntimeException("Unknown parameter '"+arg+"'");}
			final String a=split[0].toLowerCase();
			final String b=split[1];
			if(a.equals("ddls") || a.equals("dlls")){numDDLs=Parse.parseIntKMG(b);}
			else if(a.equals("buckets")){buckets=Parse.parseIntKMG(b);}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("maxmult") || a.equals("mult")){maxMult=Double.parseDouble(b);}
			else if(a.equals("dupfactor") || a.equals("dupefactor") || a.equals("dupes")){dupFactor=Integer.parseInt(b);}
			else if(a.equals("reportfrac")){reportFrac=Double.parseDouble(b);}
			else if(a.equals("seed")){masterSeed=Long.parseLong(b);}
			else if(a.equals("valseed")){/*ignored*/}
			else if(a.equals("threads") || a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("out") || a.equals("out1")){out1=b;}
			else if(a.equals("out3")){out3=b;}
			else if(a.equals("out4")){out4=b;}
			else if(a.equals("loglogtype") || a.equals("type")){loglogtype=b.toLowerCase();}
			else if(a.equals("cf") || a.equals("loglogcf")){CorrectionFactor.USE_CORRECTION=Parse.parseBoolean(b);}
			else if(a.equals("cardcf")){CardinalityTracker.USE_CARD_CF=Parse.parseBoolean(b);}
			else if(a.equals("cffile")){cffile=b; CorrectionFactor.USE_CORRECTION=true;}
			else if(a.equals("dlcalpha") || a.equals("alpha")){CardinalityStats.DLC_ALPHA=Float.parseFloat(b); AbstractCardStats.DLC_ALPHA=Float.parseFloat(b);}
			else if(a.equals("dlctarget")){CardinalityStats.DLC_TARGET_FRAC=Float.parseFloat(b); AbstractCardStats.DLC_TARGET_FRAC=Float.parseFloat(b);}
			else if(a.equals("dlcblendlo")){CardinalityStats.DLC_BLEND_LO=Float.parseFloat(b); AbstractCardStats.DLC_BLEND_LO=Float.parseFloat(b);}
			else if(a.equals("dlcblendhi")){CardinalityStats.DLC_BLEND_HI=Float.parseFloat(b); AbstractCardStats.DLC_BLEND_HI=Float.parseFloat(b);}
			else if(a.equals("cfiters") || a.equals("cfiterations")){CardinalityStats.DEFAULT_CF_ITERS=Integer.parseInt(b); AbstractCardStats.DEFAULT_CF_ITERS=Integer.parseInt(b);}
			else if(a.equals("cfdif") || a.equals("cfconvergence")){CardinalityStats.DEFAULT_CF_DIF=Double.parseDouble(b); AbstractCardStats.DEFAULT_CF_DIF=Double.parseDouble(b);}
			else if(a.equals("cfmult") || a.equals("minseedmult")){CardinalityStats.MIN_SEED_CF_MULT=Float.parseFloat(b); AbstractCardStats.MIN_SEED_CF_MULT=Float.parseFloat(b);}
			else if(a.equals("minvfraction") || a.equals("minvk")){
				float x=Float.parseFloat(b);
				if(x<1){CardinalityStats.DLC_MIN_VK_FRACTION=x; AbstractCardStats.DLC_MIN_VK_FRACTION=x;}
				else{CardinalityStats.DLC_MIN_VK=(int)x; AbstractCardStats.DLC_MIN_VK=(int)x;}
			}
			else if(a.equals("dlcinfopow") || a.equals("infopow")){AbstractCardStats.DLC_INFO_POWER=Float.parseFloat(b);}
			else if(a.equals("hcinfopow") || a.equals("hcpow")){AbstractCardStats.HC_INFO_POWER=Float.parseFloat(b);}
			else if(a.equals("dlchistblendlo")){AbstractCardStats.DLCSBS_BLEND_LO=Float.parseFloat(b);}
			else if(a.equals("dlchistblendhi")){AbstractCardStats.DLCSBS_BLEND_HI=Float.parseFloat(b);}
			else if(a.equals("ldlcblo")){AbstractCardStats.LDLC_B_LO=Float.parseFloat(b);}
			else if(a.equals("ldlcbhi")){AbstractCardStats.LDLC_B_HI=Float.parseFloat(b);}
			else if(a.equals("dlcinfomode") || a.equals("infomode")){AbstractCardStats.DLC_INFO_MODE=Integer.parseInt(b);}
			else if(a.equals("promotethreshold") || a.equals("pt")){
				DynamicLogLog3.PROMOTE_THRESHOLD=Integer.parseInt(b);
				DynamicLogLog3v2.PROMOTE_THRESHOLD=Integer.parseInt(b);
				DynamicLogLog4.PROMOTE_THRESHOLD=Integer.parseInt(b);
			}else if(a.equals("promotefrac") || a.equals("pf")){DynamicLogLog3v2.PROMOTE_FRAC=Float.parseFloat(b);
			}else if(a.equals("resetonpromote") || a.equals("rop")){DynamicLogLog3v2.RESET_ON_PROMOTE=Parse.parseBoolean(b);
			}else if(a.equals("earlypromote") || a.equals("ep")){
				DynamicLogLog3.EARLY_PROMOTE=Parse.parseBoolean(b);
				DynamicLogLog4.EARLY_PROMOTE=Parse.parseBoolean(b);
				BankedDynamicLogLog3.EARLY_PROMOTE=Parse.parseBoolean(b);
			}else if(a.equals("correctoverflow") || a.equals("co")){
				DynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
				BankedDynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
				if(Parse.parseBoolean(b)){DynamicLogLog3.IGNORE_OVERFLOW=false; DynamicLogLog2.IGNORE_OVERFLOW=false; BankedDynamicLogLog3.IGNORE_OVERFLOW=false;}
			}else if(a.equals("ignoreoverflow") || a.equals("io")){
				DynamicLogLog3.IGNORE_OVERFLOW=Parse.parseBoolean(b);
				DynamicLogLog2.IGNORE_OVERFLOW=Parse.parseBoolean(b);
				BankedDynamicLogLog3.IGNORE_OVERFLOW=Parse.parseBoolean(b);
				if(Parse.parseBoolean(b)){DynamicLogLog3.CORRECT_OVERFLOW=false; BankedDynamicLogLog3.CORRECT_OVERFLOW=false;}
			}else if(a.equals("overflowscale") || a.equals("os")){
				DynamicLogLog3.OVERFLOW_SCALE=Double.parseDouble(b);
				BankedDynamicLogLog3.OVERFLOW_SCALE=Double.parseDouble(b);
			}else if(a.equals("usestoredoverflow") || a.equals("uso")){
				DynamicLogLog3.USE_STORED_OVERFLOW=Parse.parseBoolean(b);
				BankedDynamicLogLog3.USE_STORED_OVERFLOW=Parse.parseBoolean(b);
			}else if(a.equals("calcfgra") || a.equals("fgra")){UltraDynamicLogLog6.CALC_FGRA=Parse.parseBoolean(b);
			}else if(a.equals("printdlctiers")){DDLCalibrationDriver.PRINT_DLC_TIERS=Parse.parseBoolean(b);
			}else if(a.equals("printstd")){DDLCalibrationDriver.PRINT_STD=Parse.parseBoolean(b);
			}else if(a.equals("out2")){System.err.println("Note: out2= is not supported by DDLCalibrationDriver2; ignoring.");
			}else if(a.equals("notes")){notes=b.replace('_',' ');
			}else if(a.equals("frozenhistory") || a.equals("frozen")){UltraLogLog8.FROZEN_HISTORY=Parse.parseBoolean(b);
			}else if(a.equals("statepower") || a.equals("sp")){UltraLogLog8.STATE_POWER=Double.parseDouble(b);
			}else if(a.equals("statecfoffset") || a.equals("sco")){UltraLogLog8.STATE_CF_OFFSET=Double.parseDouble(b);
			}else if(a.equals("pllmode") || a.equals("pmode")){pllmode=b;
			}else if(a.equals("plloffset") || a.equals("pco")){
				ProtoLogLog16b.CF_OFFSET=Double.parseDouble(b);
				ProtoLogLog16.CF_OFFSET=Double.parseDouble(b);
				ProtoLogLog16c.CF_OFFSET=Double.parseDouble(b);
			}else if(a.equals("hbits")){int v=Integer.parseInt(b); ProtoLogLog16.HISTORY_BITS=v; ProtoLogLog16b.HISTORY_BITS=v; ProtoLogLog16c.HISTORY_BITS=v;
			}else if(a.equals("lbits")){int v=Integer.parseInt(b); ProtoLogLog16.LUCK_BITS=v; ProtoLogLog16b.LUCK_BITS=v; ProtoLogLog16c.LUCK_BITS=v;
			}else if(a.equals("mbits")){int v=Integer.parseInt(b); ProtoLogLog16.MANTISSA_BITS=v; ProtoLogLog16b.MANTISSA_BITS=v; ProtoLogLog16c.MANTISSA_BITS=v;
			}else if(a.equals("abits")){int v=Integer.parseInt(b); ProtoLogLog16.ANDTISSA_BITS=v; ProtoLogLog16b.ANDTISSA_BITS=v; ProtoLogLog16c.ANDTISSA_BITS=v;
			}else if(a.equals("nbits")){int v=Integer.parseInt(b); ProtoLogLog16.NLZ2_BITS=v; ProtoLogLog16b.NLZ2_BITS=v; ProtoLogLog16c.NLZ2_BITS=v;
			}else if(a.equals("emode") || a.equals("estimatemode")){ProtoLogLog16b.ESTIMATE_MODE=Integer.parseInt(b);
			}else if(a.equals("nocf")){if(Parse.parseBoolean(b)){ProtoLogLog16b.CORRECTION_TABLE=null; ProtoLogLog16.CORRECTION_TABLE=null;}
			}else if(a.equals("empiricalmantissa") || a.equals("em")){DynamicDemiLog8.USE_EMPIRICAL_MANTISSA=Parse.parseBoolean(b);
			}else if(a.equals("mantissacfoffset") || a.equals("mco")){DynamicDemiLog8.MANTISSA_CF_OFFSET=Double.parseDouble(b);
			}else if(a.equals("printcv") || a.equals("cv")){DDLCalibrationDriver.PRINT_CV=Parse.parseBoolean(b);
			}else if(a.equals("clamp") || a.equals("clamptoadded")){CLAMP_TO_ADDED=Parse.parseBoolean(b);
			}else if(a.equals("fll2mult")){FutureLogLog2.TERMINAL_CORRECTION=Double.parseDouble(b);
			}else if(a.equals("hllhistcf") || a.equals("histcf")){CardinalityStats.HLL_HIST_TERMINAL_CF=Double.parseDouble(b);
			}else if(a.equals("tracecf")){CorrectionFactor.TRACE_CF=Parse.parseBoolean(b);
			}else if(a.equals("saturate") || a.equals("sat")){UltraDynamicLogLog6.SATURATE_ON_OVERFLOW=Parse.parseBoolean(b);
			}else if(a.equals("histlc") || a.equals("historylc")){CardinalityTracker.USE_HISTORY_FOR_LC=Parse.parseBoolean(b);
			}else if(a.equals("lchistfile")){CorrectionFactor.sbsFile=b;
			}else if(a.equals("lchistmultfile")){CorrectionFactor.sbsMultFile=b;
			}else if(a.equals("microlc")){CardinalityStats.USE_MICRO_FOR_LC=Parse.parseBoolean(b); AbstractCardStats.USE_MICRO_FOR_LC=Parse.parseBoolean(b);
			}else if(a.equals("usemicro") || a.equals("usemicroindex")){AbstractCardStats.USE_MICRO_INDEX=Parse.parseBoolean(b);
			}else if(a.equals("lchisthybrid") || a.equals("lchybrid")){CardinalityTracker.USE_SBS_IN_HYBRID=Parse.parseBoolean(b);
			}else if(a.equals("sbsformula") || a.equals("usesbsformula")){CorrectionFactor.USE_SBS_FORMULA=Parse.parseBoolean(b);
			}else if(a.equals("meancfformula") || a.equals("usemeancfformula")){CorrectionFactor.USE_MEAN_CF_FORMULA=Parse.parseBoolean(b);
			}else if(a.equals("hccfformula") || a.equals("usehccfformula")){CorrectionFactor.USE_HC_CF_FORMULA=Parse.parseBoolean(b);
			}else if(a.equals("formulas") || a.equals("useformulas")){CorrectionFactor.USE_FORMULAS=Parse.parseBoolean(b);
			}else if(a.equals("hcweight") || a.equals("ldlcweight")){CardinalityTracker.LDLC_HC_WEIGHT=Double.parseDouble(b); hcWeightExplicit=true;
			}else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}
	}

	/** Configure global state: class init, CF loading, formula coefficients.
	 *  Must be called inside synchronized(DDLCalibrationDriver2.class). */
	private void initGlobalState(){
		if(pllmode!=null){
			if(pllmode.equals("mantissa")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_MANTISSA); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_MANTISSA); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_MANTISSA);}
			else if(pllmode.equals("andtissa")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_ANDTISSA); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_ANDTISSA); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_ANDTISSA);}
			else if(pllmode.equals("nlz2")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_NLZ2); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_NLZ2); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_NLZ2);}
			else if(pllmode.equals("history")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_HISTORY); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_HISTORY); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_HISTORY);}
			else if(pllmode.equals("luck")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_LUCK); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_LUCK); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_LUCK);}
			else if(pllmode.equals("none")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_NONE); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_NONE); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_NONE);}
			else if(pllmode.equals("histmant") || pllmode.equals("historymantissa")){ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_HISTORY|ProtoLogLog16c.MODE_MANTISSA);}
			else{throw new RuntimeException("Unknown pllmode: "+pllmode);}
		}
		if(!hcWeightExplicit){
			final int hb=ProtoLogLog16c.HISTORY_BITS;
			if(hb==1){CardinalityTracker.LDLC_HC_WEIGHT=0.456;}
			else if(hb==2){CardinalityTracker.LDLC_HC_WEIGHT=0.50;}
			else if(hb==3){CardinalityTracker.LDLC_HC_WEIGHT=0.475;}
		}
		maxTrue=Math.max(1, (long)(buckets*maxMult));
		DDLCalibrationDriver.makeInstance(loglogtype, buckets, k, 0L, 0);
		DDLCalibrationDriver.v3ColsForType(loglogtype);

		// Per-class Mean CF formula coefficients
		boolean classHasMeanCf=true, classHasHcCf=false, classHasHmeanmCf=false, classHasMeanhCf=false;
		if(loglogtype.equals("dll4") || loglogtype.equals("dll4m")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4;}
		else if(loglogtype.equals("ll6")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_LL6;}
		else if(loglogtype.equals("dll3") || loglogtype.equals("dll3v2")){
			if(DynamicLogLog3.IGNORE_OVERFLOW){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL3_IOF;}
			else{CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL3_IOT;}
		}else if(loglogtype.equals("dll2")){
			if(DynamicLogLog2.IGNORE_OVERFLOW){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL2_IOF;}
			else{CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL2_IOT;}
		}else if(loglogtype.equals("bdll3")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_BDLL3_COF;}
		else if(loglogtype.equals("udll6")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_UDLL6;
			CorrectionFactor.meanhCfCoeffs=CorrectionFactor.MCF_UDLL6_MEANH;
			classHasHcCf=true; classHasMeanhCf=true;
		}else if(loglogtype.equals("ddl") || loglogtype.equals("ddl10")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DDL; CorrectionFactor.hmeanmCfCoeffs=CorrectionFactor.HMCF_DDL; classHasHmeanmCf=true;}
		else if(loglogtype.equals("ddl8")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DDL8; CorrectionFactor.hmeanmCfCoeffs=CorrectionFactor.HMCF_DDL8; classHasHmeanmCf=true;}
		else{classHasMeanCf=false;}
		if(CorrectionFactor.USE_FORMULAS){
			if(!classHasMeanCf){CorrectionFactor.USE_MEAN_CF_FORMULA=false;}
			if(!classHasHcCf){CorrectionFactor.USE_HC_CF_FORMULA=false;}
			if(!classHasHmeanmCf){CorrectionFactor.USE_HMEANM_CF_FORMULA=false;}
		}
		CorrectionFactor.loadSbsTable();
		CorrectionFactor.loadSbsMultTable();
		if("fll2".equals(loglogtype)){FutureLogLog2.loadCFTable(); FutureLogLog2.loadCardCFTable();}
		if(cffile!=null){
			CorrectionFactor.initialize(cffile, buckets);
			if("pll16b".equals(loglogtype)){ProtoLogLog16b.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("pll16c".equals(loglogtype)){ProtoLogLog16c.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("udll6".equals(loglogtype)){UltraDynamicLogLog6.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("ertl".equals(loglogtype)){ErtlULL.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
		}
		// Publish immutable CF snapshot for worker threads
		CorrectionFactor.publishSnapshot();
		thresholds=DDLCalibrationDriver.computeThresholds(maxTrue, reportFrac);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void process(final long t0){
		final int numThreads=Math.min(threads, numDDLs);
		final int numThresholds=thresholds.length;
		final CalThread[] calThreads=spawnThreads(numThreads);
		for(CalThread ct : calThreads){
			try{ct.join();}catch(InterruptedException e){Thread.currentThread().interrupt();}
			if(!ct.success){System.err.println("Warning: a calibration thread did not complete successfully.");}
		}
		final MergedResults merged=accumulate(calThreads, numThresholds);
		final ArrayList<DDLCalibrationDriver.ReportRow> mergedRows=buildRows(merged, numThresholds);
		printSummary(mergedRows, merged, numThresholds, t0);
		writeFile1(mergedRows, merged, numThresholds);
		if(out3!=null){writeFile3(mergedRows);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	private CalThread[] spawnThreads(int numThreads){
		final CalThread[] calThreads=new CalThread[numThreads];
		for(int t=0; t<numThreads; t++){
			final int threadDDLs=numDDLs/numThreads+(t<numDDLs%numThreads ? 1 : 0);
			calThreads[t]=new CalThread(masterSeed+t, threadDDLs, thresholds,
				buckets, k, maxTrue, loglogtype, out4, dupFactor);
		}
		for(CalThread ct : calThreads){ct.start();}
		return calThreads;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Accumulation          ----------------*/
	/*--------------------------------------------------------------*/

	static final class MergedResults {
		MergedResults(int nt){
			sumErr=new double[nt][NUM_EST]; sumAbsErr=new double[nt][NUM_EST]; sumSqErr=new double[nt][NUM_EST];
			ldlcErr=new double[nt][NUM_LDLC]; ldlcAbsErr=new double[nt][NUM_LDLC]; ldlcSqErr=new double[nt][NUM_LDLC];
			occSum=new double[nt]; n=new int[nt];
		}
		final double[][] sumErr, sumAbsErr, sumSqErr;
		final double[][] ldlcErr, ldlcAbsErr, ldlcSqErr;
		final double[] occSum;
		final int[] n;
	}

	private MergedResults accumulate(CalThread[] calThreads, int numThresholds){
		final MergedResults m=new MergedResults(numThresholds);
		for(CalThread ct : calThreads){
			synchronized(ct){
				for(int ti=0; ti<numThresholds; ti++){
					m.n[ti]+=ct.n[ti]; m.occSum[ti]+=ct.occSum[ti];
					for(int e=0; e<NUM_EST; e++){
						m.sumErr[ti][e]+=ct.sumErr[ti][e];
						m.sumAbsErr[ti][e]+=ct.sumAbsErr[ti][e];
						m.sumSqErr[ti][e]+=ct.sumSqErr[ti][e];
					}
					for(int e=0; e<NUM_LDLC; e++){
						m.ldlcErr[ti][e]+=ct.ldlcSumErr[ti][e];
						m.ldlcAbsErr[ti][e]+=ct.ldlcSumAbsErr[ti][e];
						m.ldlcSqErr[ti][e]+=ct.ldlcSumSqErr[ti][e];
					}
				}
			}
		}
		return m;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Report Building      ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<DDLCalibrationDriver.ReportRow> buildRows(MergedResults m, int numThresholds){
		final int NUM_OUT1=DDLCalibrationDriver.NUM_OUT1;
		final ArrayList<DDLCalibrationDriver.ReportRow> rows=new ArrayList<>(numThresholds);
		for(int ti=0; ti<numThresholds; ti++){
			if(m.n[ti]<1){continue;}
			final double[] rErr=new double[NUM_OUT1], rAbsErr=new double[NUM_OUT1], rSqErr=new double[NUM_OUT1];
			System.arraycopy(m.sumErr[ti], 0, rErr, 0, NUM_EST);
			System.arraycopy(m.sumAbsErr[ti], 0, rAbsErr, 0, NUM_EST);
			System.arraycopy(m.sumSqErr[ti], 0, rSqErr, 0, NUM_EST);
			rows.add(new DDLCalibrationDriver.ReportRow(thresholds[ti], m.n[ti], m.occSum[ti], rErr, rAbsErr, rSqErr));
		}
		return rows;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Output Methods       ----------------*/
	/*--------------------------------------------------------------*/

	private void printSummary(ArrayList<DDLCalibrationDriver.ReportRow> mergedRows,
			MergedResults m, int numThresholds, long t0){
		final double elapsed=(System.nanoTime()-t0)*1e-9;
		System.err.println();
		System.err.println("=== DDL2 Calibration Summary ===");
		System.err.println("Type: "+loglogtype+"  Buckets: "+buckets+"  DDLs: "+numDDLs
			+"  MaxCard: "+maxTrue+"  Rows: "+mergedRows.size()
			+"  Elapsed: "+String.format("%.1f", elapsed)+"s"
			+(dupFactor>1 ? "  DupFactor: "+dupFactor : ""));
		System.err.println("--- Log-Weighted and Linear-Weighted Avg Absolute Error, Peak, Signed (lower = better) ---");
		final int rows=mergedRows.size();
		final double[] totalAbsErr=new double[NUM_EST], linWtAbsErr=new double[NUM_EST];
		final double[] peakAbsErr=new double[NUM_EST], totalSignErr=new double[NUM_EST], totalCV=new double[NUM_EST];
		double linWeightSum=0; int cvRows=0;
		for(DDLCalibrationDriver.ReportRow row : mergedRows){
			final double card=row.trueCard; linWeightSum+=card;
			for(int e=0; e<NUM_EST; e++){
				final double meanErr=row.sumErr[e]/row.n;
				final double absAtRow=row.sumAbsErr[e]/row.n;
				totalAbsErr[e]+=absAtRow; linWtAbsErr[e]+=absAtRow*card; totalSignErr[e]+=meanErr;
				if(absAtRow>peakAbsErr[e]){peakAbsErr[e]=absAtRow;}
				final double var=row.sumSqErr[e]/row.n-meanErr*meanErr;
				final double sd=Math.sqrt(Math.max(0, var));
				final double denom=Math.abs(1.0+meanErr);
				if(denom>0){totalCV[e]+=sd/denom;}
			}
			cvRows++;
		}
		System.err.println(String.format("%-12s %-12s %-12s %-12s %-12s %s",
			"", "LogWtAbsErr", "LinWtAbsErr", "PeakAbsErr", "AvgSignErr", "AvgCV"));
		for(int e=0; e<Math.min(17, NUM_EST); e++){
			System.err.println(String.format("%-12s %.8f  %.8f  %.8f  %+.8f  %.8f",
				ESTIMATOR_NAMES[e], totalAbsErr[e]/rows,
				linWeightSum>0 ? linWtAbsErr[e]/linWeightSum : 0, peakAbsErr[e],
				rows>0 ? totalSignErr[e]/rows : 0, cvRows>0 ? totalCV[e]/cvRows : 0));
		}
		// LDLC summary
		final double[] ldlcTotAbs=new double[NUM_LDLC], ldlcLinWt=new double[NUM_LDLC];
		final double[] ldlcPeak=new double[NUM_LDLC], ldlcTotSign=new double[NUM_LDLC];
		double ldlcLinSum=0;
		for(int ti=0; ti<numThresholds; ti++){
			if(m.n[ti]<1){continue;}
			final double card=thresholds[ti]; ldlcLinSum+=card;
			for(int e=0; e<NUM_LDLC; e++){
				final double absAtRow=m.ldlcAbsErr[ti][e]/m.n[ti];
				ldlcTotAbs[e]+=absAtRow; ldlcLinWt[e]+=absAtRow*card;
				ldlcTotSign[e]+=m.ldlcErr[ti][e]/m.n[ti];
				if(absAtRow>ldlcPeak[e]){ldlcPeak[e]=absAtRow;}
			}
		}
		for(int e=0; e<NUM_LDLC; e++){
			System.err.println(String.format("%-12s %.8f  %.8f  %.8f  %+.8f  %.8f",
				LDLC_NAMES[e], ldlcTotAbs[e]/rows,
				ldlcLinSum>0 ? ldlcLinWt[e]/ldlcLinSum : 0, ldlcPeak[e],
				rows>0 ? ldlcTotSign[e]/rows : 0, 0.0));
		}
		System.err.println();
	}

	private void writeFile1(ArrayList<DDLCalibrationDriver.ReportRow> mergedRows,
			MergedResults m, int numThresholds){
		final ByteStreamWriter bsw=new ByteStreamWriter(out1, true, false, false);
		bsw.start();
		{
			final StringBuilder hdr=new StringBuilder(DDLCalibrationDriver.header1());
			for(int e=0; e<NUM_LDLC; e++){
				hdr.append('\t').append(LDLC_NAMES[e]).append("_err");
				hdr.append('\t').append(LDLC_NAMES[e]).append("_abs");
			}
			bsw.println(hdr.toString());
		}
		int finalTi=0;
		for(DDLCalibrationDriver.ReportRow row : mergedRows){
			final StringBuilder sb=new StringBuilder(DDLCalibrationDriver.formatRow1(row));
			for(int e=0; e<NUM_LDLC; e++){
				final int ni=row.n;
				sb.append('\t').append(String.format("%.6f", m.ldlcErr[finalTi][e]/ni));
				sb.append('\t').append(String.format("%.6f", m.ldlcAbsErr[finalTi][e]/ni));
			}
			bsw.println(sb.toString());
			finalTi++;
		}
		bsw.poisonAndWait();
	}

	private void writeFile3(ArrayList<DDLCalibrationDriver.ReportRow> mergedRows){
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

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	static final class CalThread extends Thread {

		CalThread(long seed, int numDDLs, long[] thresholds,
				int buckets, int k, long maxTrue, String loglogtype, String out4, int dupFactor){
			this.seed=seed; this.numDDLs=numDDLs; this.thresholds=thresholds;
			this.buckets=buckets; this.k=k; this.maxTrue=maxTrue;
			this.loglogtype=loglogtype; this.out4=out4; this.dupFactor=dupFactor;
			final int nt=thresholds.length;
			n=new int[nt]; occSum=new double[nt];
			sumErr=new double[nt][NUM_EST]; sumAbsErr=new double[nt][NUM_EST]; sumSqErr=new double[nt][NUM_EST];
			ldlcSumErr=new double[nt][NUM_LDLC]; ldlcSumAbsErr=new double[nt][NUM_LDLC]; ldlcSumSqErr=new double[nt][NUM_LDLC];
		}

		@Override
		public void run(){
			synchronized(DDLCalibrationDriver2.class){
				// Ensures visibility of all global state set before thread creation on NUMA
			}
			synchronized(this){
				try{runInner(); success=true;}
				catch(Exception e){e.printStackTrace();}
			}
		}

		void runInner(){
			final FastRandomXoshiro rng=new FastRandomXoshiro(seed);
			java.io.PrintWriter out4Pw=null;
			if(out4!=null){try{out4Pw=new java.io.PrintWriter(new java.io.FileOutputStream(out4));}catch(Exception e){e.printStackTrace();}}
			for(int di=0; di<numDDLs; di++){
				final long ddlSeed=rng.nextLong()&Long.MAX_VALUE;
				final CardinalityTracker ddl=DDLCalibrationDriver.makeInstance(loglogtype, buckets, k, ddlSeed, 0);
				int ti=0;
				for(long trueCard=1; trueCard<=maxTrue; trueCard++){
					{final long val=rng.nextLong(); for(int d=0; d<dupFactor; d++){ddl.add(val);}}
					if(trueCard>=thresholds[ti]){
						final double occ=ddl.occupancy();
						final double[] est=ddl.rawEstimates();
						occSum[ti]+=occ; n[ti]++;
						if(di==0 && out4Pw!=null){
							int[] raw=null, corr=null; int mz=0;
							if(ddl.getClass()==DynamicLogLog3.class){
								final DynamicLogLog3 d=(DynamicLogLog3)ddl; raw=d.lastRawNlz; corr=d.lastCorrNlz; mz=d.getMinZeros();
							}else if(ddl.getClass()==DynamicLogLog4.class){
								final DynamicLogLog4 d=(DynamicLogLog4)ddl; raw=d.lastRawNlz; corr=d.lastCorrNlz; mz=d.getMinZeros();
							}
							if(raw!=null){
								out4Pw.print(trueCard+"\t"+mz);
								for(int tt=0; tt<20; tt++){out4Pw.print("\t"+raw[tt]);}
								for(int tt=0; tt<20; tt++){out4Pw.print("\t"+corr[tt]);}
								out4Pw.println();
							}
						}
						for(int e=0; e<Math.min(NUM_EST, est.length); e++){
							final double v=CLAMP_TO_ADDED ? Math.min(est[e], trueCard) : est[e];
							final double err=(v-trueCard)/(double)trueCard;
							sumErr[ti][e]+=err; sumAbsErr[ti][e]+=Math.abs(err); sumSqErr[ti][e]+=err*err;
						}
						if(ddl.getClass()==UltraDynamicLogLog6.class){
							final UltraDynamicLogLog6 u=(UltraDynamicLogLog6)ddl;
							final CardStats cs=u.consumeLastSummarized();
							final double[] ldlcVals={cs.ldlc(), cs.dlcSbs(), cs.hc(),
								u.fgraEstimate(), cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2()};
							for(int e=0; e<7; e++){
								final double v=ldlcVals[e];
								final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
								ldlcSumErr[ti][e]+=lerr; ldlcSumAbsErr[ti][e]+=Math.abs(lerr); ldlcSumSqErr[ti][e]+=lerr*lerr;
							}
						}else if(ddl.getClass()==ProtoLogLog16c.class){
							final ProtoLogLog16c p=(ProtoLogLog16c)ddl;
							final double[] ldlcR=p.ldlcEstimate();
							final int[] ldlcIdx={0, 1, 2, 4, 5, 6, 7};
							for(int e=0; e<7; e++){
								final double v=ldlcR[ldlcIdx[e]];
								final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
								ldlcSumErr[ti][e]+=lerr; ldlcSumAbsErr[ti][e]+=Math.abs(lerr); ldlcSumSqErr[ti][e]+=lerr*lerr;
							}
						}
						{
							final int[] extraIdx={AbstractCardStats.LC_NOMICRO_IDX, AbstractCardStats.SBS_NOMICRO_IDX};
							for(int e=0; e<extraIdx.length; e++){
								final double v=(extraIdx[e]<est.length ? est[extraIdx[e]] : 0);
								final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
								ldlcSumErr[ti][7+e]+=lerr; ldlcSumAbsErr[ti][7+e]+=Math.abs(lerr); ldlcSumSqErr[ti][7+e]+=lerr*lerr;
							}
						}
						ti++;
						if(ti>=thresholds.length){break;}
					}
				}
			}
			if(out4Pw!=null){out4Pw.close();}
		}

		final long seed; final int numDDLs; final long[] thresholds;
		final int buckets; final int k; final long maxTrue;
		final String loglogtype; final String out4; final int dupFactor;
		final int[] n; final double[] occSum;
		final double[][] sumErr, sumAbsErr, sumSqErr;
		final double[][] ldlcSumErr, ldlcSumAbsErr, ldlcSumSqErr;
		volatile boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	int numDDLs=128;
	int buckets=2048;
	int k=31;
	double maxMult=10;
	double reportFrac=0.01;
	long masterSeed=12345L;
	int threads=Shared.threads();
	int dupFactor=1;
	String out1="stdout.txt";
	String out3=null;
	String out4=null;
	String loglogtype="ddl";
	String notes="";
	String cffile=null;
	String pllmode=null;
	boolean hcWeightExplicit=false;
	long maxTrue;
	long[] thresholds;

}
