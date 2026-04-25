package cardinality;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
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
	static final int NUM_LDLC=14;
	static final int HLDLC_IDX=11, MEAN16_IDX=12, DUALLC_IDX=13;
	static final int NUM_LDLC_BASE=7;
	static final String[] LDLC_NAMES={"LDLC", "DLC_L", "HC", "FGRA", "HLL+H", "Mean+H", "Hybrid+2", "LC_noMicro", "SBS_noMicro", "WordEst", "WordEstCV",
		"HLDLC", "Mean16", "DualLC"};

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
			Parser.printSetThreads=false;
			parse(args);
			loglogtype=Parser.loglogType.toLowerCase();
			threads=Shared.threads();
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
			else if(a.equals("card") || a.equals("maxcard") || a.equals("cardinality")){maxCard=Parse.parseKMG(b);}
			else if(a.equals("dupfactor") || a.equals("dupefactor") || a.equals("dupes")){dupFactor=Integer.parseInt(b);}
			else if(a.equals("reportfrac")){reportFrac=Double.parseDouble(b);}
			else if(a.equals("seed")){masterSeed=Long.parseLong(b);}
			else if(a.equals("valseed")){/*ignored*/}
			else if(a.equals("out") || a.equals("out1")){out1=b;}
			else if(a.equals("out3")){out3=b;}
			else if(a.equals("out4")){out4=b;}
			else if(a.equals("hcweight") || a.equals("ldlcweight")){CardinalityParser.parse(arg, a, b);}
			else if(CardinalityParser.parse(arg, a, b)){}
			else if(a.startsWith("hsbtable")){System.err.println("Note: hsbtable= is deprecated; HSB values are now hardcoded in StateTable.");}
			else if(a.equals("hsbttll4")){
				String[] parts=b.split(",");
				double[] tbl=new double[parts.length];
				for(int i=0; i<parts.length; i++){tbl[i]=Double.parseDouble(parts[i]);}
				StateTable.CF_TTLL_4_OVERRIDE=tbl;
				StateTable.USE_TTLL_HSB=true;
			}
			else if(a.equals("termcf")){StateTable.terminalCFOverride=Double.parseDouble(b);}
			else if(a.equals("sbsfile")){
				sbsFileOverride=b;
			}
			else if(a.equals("cffile")){cffile=b; CorrectionFactor.USE_CORRECTION=true;}
			else if(a.equals("out2")){System.err.println("Note: out2= is not supported by DDLCalibrationDriver2; ignoring.");
			}else if(a.equals("notes")){notes=b.replace('_',' ');
			}else if(a.equals("pllmode") || a.equals("pmode")){pllmode=b;
			}else if(a.equals("clamp") || a.equals("clamptoadded")){CLAMP_TO_ADDED=Parse.parseBoolean(b);
			}else if(a.equals("lcmode") || a.equals("lc")){lcMode=Parse.parseBoolean(b);
			}else if(!Parser.parseStatic(arg, a, b)){throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}
	}

	/** Configure global state: class init, CF loading, formula coefficients.
	 *  Must be called inside synchronized(DDLCalibrationDriver2.class). */
	private void initGlobalState(){
		maxTrue=(maxCard>0) ? maxCard : Math.max(1, (long)(buckets*maxMult));
		CardinalityParser.initializeAll(loglogtype, buckets, k, cffile, pllmode);
		if(sbsFileOverride!=null){
			CorrectionFactor.sbsFile=sbsFileOverride;
			CorrectionFactor.SBS_CF_TABLE=null;
			CorrectionFactor.loadSbsTable();
		}
		thresholds=DDLCalibrationDriver.computeThresholds(maxTrue, reportFrac);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void process(final long t0){
		{
			final CardinalityTracker tmp=DDLCalibrationDriver.makeInstance(loglogtype, buckets, k, 0, 0);
			System.err.println("terminalMeanCF="+tmp.terminalMeanCF()+" terminalMeanPlusCF="+tmp.terminalMeanPlusCF());
		}
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
		final AtomicLong nextDDL=new AtomicLong(0);
		final CalThread[] calThreads=new CalThread[numThreads];
		for(int t=0; t<numThreads; t++){
			calThreads[t]=new CalThread(masterSeed, numDDLs, nextDDL, thresholds,
				buckets, k, maxTrue, loglogtype, out4, dupFactor, lcMode);
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
		System.err.println("=== Calibration Summary ===");
		final int actualBk=DDLCalibrationDriver.makeInstance(loglogtype, buckets, k, 0, 0).actualBuckets();
		System.err.println("Type: "+loglogtype+"  Buckets: "+actualBk+"  DDLs: "+numDDLs
			+"  MaxCard: "+maxTrue+"  Rows: "+mergedRows.size()
			+"  Elapsed: "+String.format("%.1f", elapsed)+"s"
			+(dupFactor>1 ? "  DupFactor: "+dupFactor : ""));
		System.err.println("--- Log-Weighted and Width-Weighted Avg Absolute Error, Peak, Signed (lower = better) ---");
		final int rows=mergedRows.size();
		// Bucket widths: trailing difference of trueCard (row 0 uses trueCard[0]-0).
		final double[] widths=new double[rows];
		widths[0]=Math.max(1, mergedRows.get(0).trueCard);
		for(int i=1; i<rows; i++){
			widths[i]=Math.max(1, mergedRows.get(i).trueCard-mergedRows.get(i-1).trueCard);
		}

		final double[] totalAbsErr=new double[NUM_EST], widthWtAbsErr=new double[NUM_EST], countWtAbsErr=new double[NUM_EST];
		final double[] peakAbsErr=new double[NUM_EST], totalSignErr=new double[NUM_EST], totalCV=new double[NUM_EST];
		double widthWeightSum=0; double countWeightSum=0; int cvRows=0;
		int rowIdx=0;
		for(DDLCalibrationDriver.ReportRow row : mergedRows){
			final double w=widths[rowIdx]; widthWeightSum+=w; countWeightSum+=row.n;
			for(int e=0; e<NUM_EST; e++){
				final double meanErr=row.sumErr[e]/row.n;
				final double absAtRow=row.sumAbsErr[e]/row.n;
				totalAbsErr[e]+=absAtRow; widthWtAbsErr[e]+=absAtRow*w; countWtAbsErr[e]+=absAtRow*row.n; totalSignErr[e]+=meanErr;
				if(absAtRow>peakAbsErr[e]){peakAbsErr[e]=absAtRow;}
				final double var=row.sumSqErr[e]/row.n-meanErr*meanErr;
				final double sd=Math.sqrt(Math.max(0, var));
				final double denom=Math.abs(1.0+meanErr);
				if(denom>0){totalCV[e]+=sd/denom;}
			}
			cvRows++;
			rowIdx++;
		}
		System.err.println(String.format("%-12s %-13s %-13s %-13s %-13s %-13s %s",
			"", "LogWtAbsErr", "WidthWtAbsErr", "CountWtAbsErr", "PeakAbsErr", "AvgSignErr", "AvgCV"));
		for(int e=0; e<Math.min(17, NUM_EST); e++){
			System.err.println(String.format("%-12s %.8f   %.8f   %.8f   %.8f   %+.8f   %.8f",
				ESTIMATOR_NAMES[e], totalAbsErr[e]/rows,
				widthWeightSum>0 ? widthWtAbsErr[e]/widthWeightSum : 0,
				countWeightSum>0 ? countWtAbsErr[e]/countWeightSum : 0,
				peakAbsErr[e],
				rows>0 ? totalSignErr[e]/rows : 0, cvRows>0 ? totalCV[e]/cvRows : 0));
		}
		// LDLC summary — trailing-difference widths, gap-aware across skipped empty thresholds
		final double[] tWidths=new double[numThresholds];
		{
			long prev=0;
			for(int ti=0; ti<numThresholds; ti++){
				if(m.n[ti]<1){tWidths[ti]=0; continue;}
				tWidths[ti]=Math.max(1, thresholds[ti]-prev);
				prev=thresholds[ti];
			}
		}

		final double[] ldlcTotAbs=new double[NUM_LDLC], ldlcWidthWt=new double[NUM_LDLC], ldlcCountWt=new double[NUM_LDLC];
		final double[] ldlcPeak=new double[NUM_LDLC], ldlcTotSign=new double[NUM_LDLC];
		final double[] ldlcTotCV=new double[NUM_LDLC];
		double ldlcWidthSum=0; double ldlcCountSum=0; int ldlcCvRows=0;
		for(int ti=0; ti<numThresholds; ti++){
			if(m.n[ti]<1){continue;}
			final double w=tWidths[ti]; ldlcWidthSum+=w; ldlcCountSum+=m.n[ti];
			for(int e=0; e<NUM_LDLC; e++){
				final double meanErr=m.ldlcErr[ti][e]/m.n[ti];
				final double absAtRow=m.ldlcAbsErr[ti][e]/m.n[ti];
				ldlcTotAbs[e]+=absAtRow; ldlcWidthWt[e]+=absAtRow*w; ldlcCountWt[e]+=absAtRow*m.n[ti];
				ldlcTotSign[e]+=meanErr;
				if(absAtRow>ldlcPeak[e]){ldlcPeak[e]=absAtRow;}
				final double var=m.ldlcSqErr[ti][e]/m.n[ti]-meanErr*meanErr;
				final double sd=Math.sqrt(Math.max(0, var));
				final double denom=Math.abs(1.0+meanErr);
				if(denom>0){ldlcTotCV[e]+=sd/denom;}
			}
			ldlcCvRows++;
		}
		for(int e=0; e<NUM_LDLC; e++){
			System.err.println(String.format("%-12s %.8f   %.8f   %.8f   %.8f   %+.8f   %.8f",
				LDLC_NAMES[e], ldlcTotAbs[e]/rows,
				ldlcWidthSum>0 ? ldlcWidthWt[e]/ldlcWidthSum : 0,
				ldlcCountSum>0 ? ldlcCountWt[e]/ldlcCountSum : 0,
				ldlcPeak[e],
				rows>0 ? ldlcTotSign[e]/rows : 0,
				ldlcCvRows>0 ? ldlcTotCV[e]/ldlcCvRows : 0));
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

		CalThread(long masterSeed, int totalDDLs, AtomicLong nextDDL, long[] thresholds,
				int buckets, int k, long maxTrue, String loglogtype, String out4, int dupFactor, boolean lcMode){
			this.masterSeed=masterSeed; this.totalDDLs=totalDDLs; this.nextDDL=nextDDL;
			this.thresholds=thresholds;
			this.buckets=buckets; this.k=k; this.maxTrue=maxTrue;
			this.loglogtype=loglogtype; this.out4=out4; this.dupFactor=dupFactor; this.lcMode=lcMode;
			this.clampToAdded=CLAMP_TO_ADDED;
		}

		@Override
		public void run(){
			synchronized(DDLCalibrationDriver2.class){
				// Ensures visibility of all global state set before thread creation on NUMA
			}
			synchronized(this){
				final int nt=thresholds.length;
				n=new int[nt]; occSum=new double[nt];
				sumErr=new double[nt][NUM_EST]; sumAbsErr=new double[nt][NUM_EST]; sumSqErr=new double[nt][NUM_EST];
				ldlcSumErr=new double[nt][NUM_LDLC]; ldlcSumAbsErr=new double[nt][NUM_LDLC]; ldlcSumSqErr=new double[nt][NUM_LDLC];
				try{runInner(); success=true;}
				catch(Exception e){e.printStackTrace();}
			}
		}

		void runInner(){
			final FastRandomXoshiro rng=new FastRandomXoshiro(0);
			java.io.PrintWriter out4Pw=null;
			if(out4!=null){try{out4Pw=new java.io.PrintWriter(new java.io.FileOutputStream(out4));}catch(Exception e){e.printStackTrace();}}
			long di;
			while((di=nextDDL.getAndIncrement())<totalDDLs){
				rng.setSeed(masterSeed+di);
				final long ddlSeed=rng.nextLong()&Long.MAX_VALUE;
				final CardinalityTracker ddl=DDLCalibrationDriver.makeInstance(loglogtype, buckets, k, ddlSeed, 0);
				int ti=0;
				if(lcMode){
					runLC(ddl, ddlSeed, ti, (int)di, out4Pw);
					continue;
				}
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
							final double v=clampToAdded ? Math.min(est[e], trueCard) : est[e];
							final double err=(v-trueCard)/(double)trueCard;
							sumErr[ti][e]+=err; sumAbsErr[ti][e]+=Math.abs(err); sumSqErr[ti][e]+=err*err;
						}
						if(ddl.getClass()==UltraDynamicLogLog6.class){
							final UltraDynamicLogLog6 u=(UltraDynamicLogLog6)ddl;
							accumulateCardStatsLdlc(u.consumeLastSummarized(), u.fgraEstimate(), trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==UltraDynamicLogLog6m.class){
							final UltraDynamicLogLog6m u=(UltraDynamicLogLog6m)ddl;
							accumulateCardStatsLdlc(u.consumeLastSummarized(), u.fgraEstimate(), trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==CompressedDynamicLogLog4.class){
							accumulateCardStatsLdlc(((CompressedDynamicLogLog4)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==BankedDynamicLogLog4.class){
							accumulateCardStatsLdlc(((BankedDynamicLogLog4)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==BankedDynamicLogLog5.class){
							accumulateCardStatsLdlc(((BankedDynamicLogLog5)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==CompressedDynamicLogLog5.class){
							accumulateCardStatsLdlc(((CompressedDynamicLogLog5)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==BankedCompressedDynamicLogLog5.class){
							accumulateCardStatsLdlc(((BankedCompressedDynamicLogLog5)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==ArithmeticCompressedDynamicLogLog5.class){
							accumulateCardStatsLdlc(((ArithmeticCompressedDynamicLogLog5)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==VariableCompressedDynamicLogLog4.class){
							accumulateCardStatsLdlc(((VariableCompressedDynamicLogLog4)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==ArithmeticCompressedDynamicLogLog4.class){
							accumulateCardStatsLdlc(((ArithmeticCompressedDynamicLogLog4)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==HalfCappedDynamicLogLog4.class){
							accumulateCardStatsLdlc(((HalfCappedDynamicLogLog4)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==ExpandedDynamicLogLog8.class){
							accumulateCardStatsLdlc(((ExpandedDynamicLogLog8)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==ArithmeticUltraDynamicLogLog32.class){
							accumulateCardStatsLdlc(((ArithmeticUltraDynamicLogLog32)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==ArithmeticUltraDynamicLogLog33.class){
							accumulateCardStatsLdlc(((ArithmeticUltraDynamicLogLog33)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==ArithmeticVariableDynamicLogLog32.class){
							accumulateCardStatsLdlc(((ArithmeticVariableDynamicLogLog32)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==ExpandedDynamicLogLog9.class){
							accumulateCardStatsLdlc(((ExpandedDynamicLogLog9)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==CompressedDynamicLogLog3.class){
							accumulateCardStatsLdlc(((CompressedDynamicLogLog3)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==BankedCompressedDynamicLogLog3.class){
							accumulateCardStatsLdlc(((BankedCompressedDynamicLogLog3)ddl).consumeLastSummarized(), 0, trueCard, ti, ddl.hldlcWeight());
						}else if(ddl.getClass()==ProtoLogLog16c.class){
							final ProtoLogLog16c p=(ProtoLogLog16c)ddl;
							accumulateLdlcEstimate(p.ldlcEstimate(), trueCard, ti);
						}else if(ddl.getClass()==TwinTailLogLog.class){
							final TwinTailLogLog t=(TwinTailLogLog)ddl;
							final double[] ldlcR=t.ldlcEstimate();
							accumulateLdlcEstimate(ldlcR, trueCard, ti);
							{
								final float hw=ddl.hldlcWeight(); final double hldlc=hw*ldlcR[0]+(1-hw)*ldlcR[7];
								final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
								ldlcSumErr[ti][HLDLC_IDX]+=lerr; ldlcSumAbsErr[ti][HLDLC_IDX]+=Math.abs(lerr); ldlcSumSqErr[ti][HLDLC_IDX]+=lerr*lerr;
							}
							{
								final double mean16=ldlcR[8];
								final double lerr=(mean16>0 ? (mean16-trueCard)/(double)trueCard : -1.0);
								ldlcSumErr[ti][MEAN16_IDX]+=lerr; ldlcSumAbsErr[ti][MEAN16_IDX]+=Math.abs(lerr); ldlcSumSqErr[ti][MEAN16_IDX]+=lerr*lerr;
							}
							if(ldlcR.length>9){
								final double dualLC=ldlcR[9];
								final double lerr=(dualLC>0 ? (dualLC-trueCard)/(double)trueCard : -1.0);
								ldlcSumErr[ti][DUALLC_IDX]+=lerr; ldlcSumAbsErr[ti][DUALLC_IDX]+=Math.abs(lerr); ldlcSumSqErr[ti][DUALLC_IDX]+=lerr*lerr;
							}
						}else if(ddl.getClass()==TwinTailLogLog3.class){
							final TwinTailLogLog3 t=(TwinTailLogLog3)ddl;
							final double[] ldlcR=t.ldlcEstimate();
							accumulateLdlcEstimate(ldlcR, trueCard, ti);
							{
								final float hw=ddl.hldlcWeight(); final double hldlc=hw*ldlcR[0]+(1-hw)*ldlcR[7];
								final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
								ldlcSumErr[ti][HLDLC_IDX]+=lerr; ldlcSumAbsErr[ti][HLDLC_IDX]+=Math.abs(lerr); ldlcSumSqErr[ti][HLDLC_IDX]+=lerr*lerr;
							}
							{
								final double mean16=ldlcR[8];
								final double lerr=(mean16>0 ? (mean16-trueCard)/(double)trueCard : -1.0);
								ldlcSumErr[ti][MEAN16_IDX]+=lerr; ldlcSumAbsErr[ti][MEAN16_IDX]+=Math.abs(lerr); ldlcSumSqErr[ti][MEAN16_IDX]+=lerr*lerr;
							}
							if(ldlcR.length>9){
								final double dualLC=ldlcR[9];
								final double lerr=(dualLC>0 ? (dualLC-trueCard)/(double)trueCard : -1.0);
								ldlcSumErr[ti][DUALLC_IDX]+=lerr; ldlcSumAbsErr[ti][DUALLC_IDX]+=Math.abs(lerr); ldlcSumSqErr[ti][DUALLC_IDX]+=lerr*lerr;
							}
						}else if(ddl.getClass()==ErtlULL.class){
							final ErtlULL u=(ErtlULL)ddl;
							// Only FGRA applies to ULL; other LDLC slots stay 0 → sentinel -1.0
							final double fgra=u.fgraEstimatePublic();
							final double lerr=(fgra>0 ? (fgra-trueCard)/(double)trueCard : -1.0);
							ldlcSumErr[ti][3]+=lerr; ldlcSumAbsErr[ti][3]+=Math.abs(lerr); ldlcSumSqErr[ti][3]+=lerr*lerr;
						}
						{
							final int[] extraIdx={AbstractCardStats.LC_NOMICRO_IDX, AbstractCardStats.SBS_NOMICRO_IDX};
							for(int e=0; e<extraIdx.length; e++){
								final double v=(extraIdx[e]<est.length ? est[extraIdx[e]] : 0);
								final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
								ldlcSumErr[ti][7+e]+=lerr; ldlcSumAbsErr[ti][7+e]+=Math.abs(lerr); ldlcSumSqErr[ti][7+e]+=lerr*lerr;
							}
						}
						if(ddl.getClass()==DynamicLogLog4.class && DynamicLogLog4.wordTable!=null){
							// Pull CF-corrected WordEst from rawEstimates() array
							final int wIdx=DDLCalibrationDriver.WORDEST_RAW_IDX;
							final double wEst=(wIdx<est.length ? est[wIdx] : 0);
							final double wEstCV=(wIdx+1<est.length ? est[wIdx+1] : 0);
							final double[] wVals={wEst, wEstCV};
							for(int e=0; e<2; e++){
								final double v=wVals[e];
								final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
								ldlcSumErr[ti][9+e]+=lerr; ldlcSumAbsErr[ti][9+e]+=Math.abs(lerr); ldlcSumSqErr[ti][9+e]+=lerr*lerr;
							}
						}
						ti++;
						if(ti>=thresholds.length){break;}
					}
				}
			}
			if(out4Pw!=null){out4Pw.close();}

		}

		// ldlcEstimate() index mapping: skip slot 3 (FGRA is at a different position in legacy arrays)
		private static final int[] LDLC_ESTIMATE_IDX={0, 1, 2, 4, 5, 6, 7};

		/** Accumulate LDLC errors from a ldlcEstimate()-style array (PLL16c, TTLL). */
		private void accumulateLdlcEstimate(double[] ldlcR, long trueCard, int ti){
			for(int e=0; e<LDLC_ESTIMATE_IDX.length; e++){
				final double v=ldlcR[LDLC_ESTIMATE_IDX[e]];
				final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
				ldlcSumErr[ti][e]+=lerr; ldlcSumAbsErr[ti][e]+=Math.abs(lerr); ldlcSumSqErr[ti][e]+=lerr*lerr;
			}
		}

		/** Accumulate LDLC errors for a CardStats-based estimator (UDLL6, CDLL4, BDLL4, BDLL5, CDLL5). */
		private void accumulateCardStatsLdlc(CardStats cs, double fgra, long trueCard, int ti, float hw){
			final double[] ldlcVals={cs.ldlc(), cs.dlcSbs(), cs.hc(),
				fgra, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2()};
			for(int e=0; e<NUM_LDLC_BASE; e++){
				final double v=ldlcVals[e];
				final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
				ldlcSumErr[ti][e]+=lerr; ldlcSumAbsErr[ti][e]+=Math.abs(lerr); ldlcSumSqErr[ti][e]+=lerr*lerr;
			}
			final double hldlc=hw*ldlcVals[0]+(1-hw)*ldlcVals[6];
			final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
			ldlcSumErr[ti][HLDLC_IDX]+=lerr; ldlcSumAbsErr[ti][HLDLC_IDX]+=Math.abs(lerr); ldlcSumSqErr[ti][HLDLC_IDX]+=lerr*lerr;
		}

		/** LC mode: replay from same seed with increasing windows. */
		void runLC(CardinalityTracker ddl, long ddlSeed, int ti, int di, java.io.PrintWriter out4Pw){
			final long valSeed=ddlSeed*3+1; // Must differ from ddlSeed to avoid hash correlation
			final FastRandomXoshiro valRng=new FastRandomXoshiro(valSeed);
			long trueCard=0;
			long currentCount=0;
			long currentLimit=1;
			while(trueCard<maxTrue){
				valRng.setSeed(valSeed);
				currentCount=0;
				while(currentCount<currentLimit){
					final long val=valRng.nextLong();
					ddl.add(val);
					currentCount++;
					if(currentCount>trueCard){
						trueCard=currentCount;
						// Check threshold
						if(ti<thresholds.length && trueCard>=thresholds[ti]){
							final double occ=ddl.occupancy();
							final double[] est=ddl.rawEstimates();
							occSum[ti]+=occ; n[ti]++;
							for(int e=0; e<Math.min(NUM_EST, est.length); e++){
								final double v=clampToAdded ? Math.min(est[e], trueCard) : est[e];
								final double err=(v-trueCard)/(double)trueCard;
								sumErr[ti][e]+=err; sumAbsErr[ti][e]+=Math.abs(err); sumSqErr[ti][e]+=err*err;
							}
							if(ddl.getClass()==DynamicLogLog4.class && DynamicLogLog4.wordTable!=null){
								// Pull CF-corrected WordEst from rawEstimates() array
								final int wIdx=DDLCalibrationDriver.WORDEST_RAW_IDX;
								final double wEst=(wIdx<est.length ? est[wIdx] : 0);
								final double wEstCV=(wIdx+1<est.length ? est[wIdx+1] : 0);
								final double[] wVals={wEst, wEstCV};
								for(int e=0; e<2; e++){
									final double v=wVals[e];
									final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
									ldlcSumErr[ti][9+e]+=lerr; ldlcSumAbsErr[ti][9+e]+=Math.abs(lerr); ldlcSumSqErr[ti][9+e]+=lerr*lerr;
								}
							}
							ti++;
							if(ti>=thresholds.length){return;}
						}
					}
				}
				// Grow window for next round
				long growth=1+valRng.nextLong(Math.max(1, currentLimit*2));
				currentLimit=Math.min(maxTrue, currentLimit+growth);
			}
		}

		final long masterSeed; final int totalDDLs; final AtomicLong nextDDL;
		final long[] thresholds;
		final int buckets; final int k; final long maxTrue;
		final String loglogtype; final String out4; final int dupFactor; final boolean lcMode;
		final boolean clampToAdded;
		int[] n; double[] occSum;
		double[][] sumErr, sumAbsErr, sumSqErr;
		double[][] ldlcSumErr, ldlcSumAbsErr, ldlcSumSqErr;
		volatile boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	int numDDLs=128;
	int buckets=2048;
	int k=31;
	double maxMult=10;
	long maxCard=-1;
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
	String sbsFileOverride=null;
	String pllmode=null;
	boolean lcMode=false;
	long maxTrue;
	long[] thresholds;

}
