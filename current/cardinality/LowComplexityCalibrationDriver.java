package cardinality;

import java.util.ArrayList;
import java.util.BitSet;
import rand.FastRandomXoshiro;
import parse.Parse;
import parse.Parser;
import shared.Shared;

/**
 * Low-complexity cardinality calibration driver.
 * Tests estimator accuracy on datasets with bounded cardinality and repeated values.
 * Each estimator draws with replacement from a virtual value array, biased
 * toward lower indices via min(rand(), rand()) to simulate skewed frequency distributions.
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
		float iterations=4;
		int buckets=2048;
		int k=31;
		int threads=4;
		long masterSeed=1;
		double reportFrac=0.01;
		String loglogtype="dll4";
		boolean minRand=true;
		String cffile=null;
		CardinalityTracker.clampToAdded=false;

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
			else if(a.equals("t") || a.equals("threads")){threads=Integer.parseInt(b);}
			else if(a.equals("seed")){masterSeed=Long.parseLong(b);}
			else if(a.equals("reportfrac")){reportFrac=Double.parseDouble(b);}
			else if(a.equals("cffile")){cffile=b; CorrectionFactor.USE_CORRECTION=true;}
			else if(a.equals("minrand") || a.equals("skew")){minRand=Parse.parseBoolean(b);}
			else if(CardinalityParser.parse(arg, a, b)){}
			else if(a.equals("loglogtype") || a.equals("type")){Parser.loglogType=b.toLowerCase();}
			else if(a.startsWith("hsbtable")){}
			else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}
		loglogtype=Parser.loglogType;
		if(threads<1){threads=Shared.threads();}

		final int NUM_EST=DDLCalibrationDriver.NUM_EST;
		final int NUM_OUT1=DDLCalibrationDriver.NUM_OUT1;
		final int NUM_EXT=DDLCalibrationDriver2.NUM_EXTENDED_TYPES;
		final String[] EXT_NAMES=DDLCalibrationDriver2.EXTENDED_NAMES;

		CardinalityParser.initializeAll(loglogtype, buckets, k, cffile, null);

		final long totalAdds=(iterations>0 ? (long)(cardinality*iterations) : Long.MAX_VALUE);
		final boolean stopOnSaturation=(iterations<=0);

		final long[] thresholds=DDLCalibrationDriver.computeThresholds(cardinality, reportFrac);
		final int numThresholds=thresholds.length;
		System.err.println("Thresholds: "+numThresholds+" points from 1 to "+cardinality);

		final FastRandomXoshiro masterRng=new FastRandomXoshiro(masterSeed);
		int lowerLen=Integer.highestOneBit((int)Math.ceil(Math.sqrt(cardinality)));
		if((long)lowerLen*lowerLen<cardinality){lowerLen<<=1;}
		final int upperLen=(cardinality+lowerLen-1)/lowerLen;
		final int lowerMask=lowerLen-1;
		final int lowerShift=Integer.numberOfTrailingZeros(lowerLen);
		final int[] lowerArray=new int[lowerLen];
		final int[] upperArray=new int[upperLen];
		for(int i=0; i<lowerLen; i++){lowerArray[i]=masterRng.nextInt();}
		for(int i=0; i<upperLen; i++){upperArray[i]=masterRng.nextInt();}

		final long[] seeds=new long[numDDLs];
		for(int i=0; i<numDDLs; i++){seeds[i]=masterRng.nextLong();}

		/* Spawn worker threads with static range partitioning */
		final long t0=System.nanoTime();
		final int perThread=(numDDLs+threads-1)/threads;
		final LCWorker[] workers=new LCWorker[threads];
		for(int tid=0; tid<threads; tid++){
			final int estStart=tid*perThread;
			final int estEnd=Math.min(estStart+perThread, numDDLs);
			workers[tid]=new LCWorker(
				estStart, estEnd, cardinality, numThresholds,
				buckets, k, totalAdds, stopOnSaturation, minRand,
				loglogtype, seeds, thresholds,
				lowerArray, upperArray, lowerMask, lowerShift,
				NUM_OUT1, NUM_EST, NUM_EXT
			);
			workers[tid].start();
		}
		for(LCWorker w : workers){
			try{w.join();}catch(InterruptedException e){throw new RuntimeException(e);}
		}

		/* Merge accumulators */
		final double[][] sumErr=new double[numThresholds][NUM_OUT1];
		final double[][] sumAbsErr=new double[numThresholds][NUM_OUT1];
		final double[][] sumSqErr=new double[numThresholds][NUM_OUT1];
		final double[][] ldlcSumErr=new double[numThresholds][NUM_EXT];
		final double[][] ldlcSumAbsErr=new double[numThresholds][NUM_EXT];
		final double[][] ldlcSumSqErr=new double[numThresholds][NUM_EXT];
		final double[] occSum=new double[numThresholds];
		final int[] nArr=new int[numThresholds];
		for(LCWorker w : workers){
			for(int ti=0; ti<numThresholds; ti++){
				if(w.nArr[ti]==0){continue;}
				nArr[ti]+=w.nArr[ti];
				occSum[ti]+=w.occSum[ti];
				for(int e=0; e<NUM_OUT1; e++){
					sumErr[ti][e]+=w.sumErr[ti][e];
					sumAbsErr[ti][e]+=w.sumAbsErr[ti][e];
					sumSqErr[ti][e]+=w.sumSqErr[ti][e];
				}
				for(int e=0; e<NUM_EXT; e++){
					ldlcSumErr[ti][e]+=w.ldlcErr[ti][e];
					ldlcSumAbsErr[ti][e]+=w.ldlcAbsErr[ti][e];
					ldlcSumSqErr[ti][e]+=w.ldlcSqErr[ti][e];
				}
			}
		}

		/* Build report rows */
		final ArrayList<DDLCalibrationDriver.ReportRow> mergedRows=new ArrayList<>(numThresholds);
		for(int ti=0; ti<numThresholds; ti++){
			if(nArr[ti]<1){continue;}
			mergedRows.add(new DDLCalibrationDriver.ReportRow(
				thresholds[ti], nArr[ti], occSum[ti], sumErr[ti], sumAbsErr[ti], sumSqErr[ti]));
		}

		/* Stderr summary */
		printSummary(loglogtype, buckets, numDDLs, cardinality, iterations,
			mergedRows, nArr, thresholds, numThresholds,
			ldlcSumErr, ldlcSumAbsErr, ldlcSumSqErr,
			NUM_EST, NUM_EXT, EXT_NAMES, t0);

		/* Stdout: data rows with LDLC columns */
		{
			final StringBuilder hdr=new StringBuilder(DDLCalibrationDriver.header1());
			for(int e=0; e<NUM_EXT; e++){
				hdr.append('\t').append(EXT_NAMES[e]).append("_err");
				hdr.append('\t').append(EXT_NAMES[e]).append("_abs");
			}
			System.out.println(hdr);
		}
		int finalTi=0;
		for(DDLCalibrationDriver.ReportRow row : mergedRows){
			final StringBuilder sb=new StringBuilder(DDLCalibrationDriver.formatRow1(row));
			final int ti=finalTi;
			final int ni=row.n;
			for(int e=0; e<NUM_EXT; e++){
				sb.append('\t').append(String.format("%.6f", ldlcSumErr[ti][e]/ni));
				sb.append('\t').append(String.format("%.6f", ldlcSumAbsErr[ti][e]/ni));
			}
			System.out.println(sb);
			finalTi++;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Worker Thread         ----------------*/
	/*--------------------------------------------------------------*/

	static class LCWorker extends Thread {

		LCWorker(int estStart, int estEnd, int cardinality, int numThresholds,
				int buckets, int k, long totalAdds, boolean stopOnSaturation,
				boolean minRand, String loglogType, long[] seeds, long[] thresholds,
				int[] lowerArray, int[] upperArray, int lowerMask, int lowerShift,
				int numOut1, int numEst, int numExt){
			this.estStart=estStart;
			this.estEnd=estEnd;
			this.cardinality=cardinality;
			this.numThresholds=numThresholds;
			this.buckets=buckets;
			this.k=k;
			this.totalAdds=totalAdds;
			this.stopOnSaturation=stopOnSaturation;
			this.minRand=minRand;
			this.loglogType=loglogType;
			this.seeds=seeds;
			this.thresholds=thresholds;
			this.lowerArray=lowerArray;
			this.upperArray=upperArray;
			this.lowerMask=lowerMask;
			this.lowerShift=lowerShift;
			this.numEst=numEst;

			sumErr=new double[numThresholds][numOut1];
			sumAbsErr=new double[numThresholds][numOut1];
			sumSqErr=new double[numThresholds][numOut1];
			ldlcErr=new double[numThresholds][numExt];
			ldlcAbsErr=new double[numThresholds][numExt];
			ldlcSqErr=new double[numThresholds][numExt];
			occSum=new double[numThresholds];
			nArr=new int[numThresholds];
		}

		@Override
		public void run(){
			final FastRandomXoshiro rng=new FastRandomXoshiro(0);
			final BitSet seen=new BitSet(cardinality);

			for(int estIdx=estStart; estIdx<estEnd; estIdx++){
				rng.setSeed(seeds[estIdx]);
				seen.clear();
				int trueCard=0;
				int ti=0;
				final CardinalityTracker est=DDLCalibrationDriver.makeInstance(
					loglogType, buckets, k, seeds[estIdx], 0);

				for(long add=0; add<totalAdds; add++){
					final int pos=minRand ?
						Math.min(rng.nextInt(cardinality), rng.nextInt(cardinality)) :
						rng.nextInt(cardinality);
					final long val=(lowerArray[pos&lowerMask]&0xFFFFFFFFL)|
						((long)(upperArray[pos>>>lowerShift])<<32);
					est.add(val);
					if(!seen.get(pos)){
						seen.set(pos);
						trueCard++;
						while(ti<numThresholds && trueCard>thresholds[ti]){ti++;}
					}
					if(ti<numThresholds && trueCard==thresholds[ti]){
						recordAtThreshold(est, trueCard, ti);
					}
					if(stopOnSaturation && trueCard==cardinality){break;}
				}
			}
		}

		private void recordAtThreshold(CardinalityTracker est, int trueCard, int ti){
			final double occ=est.occupancy();
			final double[] estimates=est.rawEstimates();
			occSum[ti]+=occ;
			nArr[ti]++;
			for(int e=0; e<Math.min(numEst, estimates.length); e++){
				final double v=estimates[e];
				final double err=(v-trueCard)/(double)trueCard;
				sumErr[ti][e]+=err;
				sumAbsErr[ti][e]+=Math.abs(err);
				sumSqErr[ti][e]+=err*err;
			}
			recordLdlc(est, trueCard, ti, estimates);
		}

		private void recordLdlc(CardinalityTracker est, int trueCard, int ti, double[] estimates){
			final int NUM_BASE=DDLCalibrationDriver2.NUM_EXTENDED_BASE;
			final Class<?> cls=est.getClass();

			if(cls==UltraDynamicLogLog6.class){
				final UltraDynamicLogLog6 u=(UltraDynamicLogLog6)est;
				final CardStats cs=u.consumeLastSummarized();
				accLdlcBase(cs.ldlc(), cs.dlcSbs(), cs.hc(), u.fgraEstimate(), cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti, NUM_BASE);
				accHldlc(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==UltraDynamicLogLog6m.class){
				final UltraDynamicLogLog6m u=(UltraDynamicLogLog6m)est;
				final CardStats cs=u.consumeLastSummarized();
				accLdlcBase(cs.ldlc(), cs.dlcSbs(), cs.hc(), u.fgraEstimate(), cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti, NUM_BASE);
				accHldlc(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==UltraDynamicLogLog36.class){
				final UltraDynamicLogLog36 u=(UltraDynamicLogLog36)est;
				final CardStats cs=u.consumeLastSummarized();
				accLdlcBase(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti, NUM_BASE);
				accHldlc(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
				accVwMean(cs, est, cs.ldlc(), trueCard, ti);
			}else if(cls==BankedDynamicLogLog5.class){
				final CardStats cs=((BankedDynamicLogLog5)est).consumeLastSummarized();
				accLdlcBase(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti, NUM_BASE);
				accHldlc(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==CompressedDynamicLogLog4.class){
				final CardStats cs=((CompressedDynamicLogLog4)est).consumeLastSummarized();
				accLdlcBase(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti, NUM_BASE);
				accHldlc(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==CompressedDynamicLogLog5.class){
				final CardStats cs=((CompressedDynamicLogLog5)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==BankedCompressedDynamicLogLog5.class){
				final CardStats cs=((BankedCompressedDynamicLogLog5)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==ArithmeticCompressedDynamicLogLog5.class){
				final CardStats cs=((ArithmeticCompressedDynamicLogLog5)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==VariableCompressedDynamicLogLog4.class){
				final CardStats cs=((VariableCompressedDynamicLogLog4)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==ArithmeticCompressedDynamicLogLog4.class){
				final CardStats cs=((ArithmeticCompressedDynamicLogLog4)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==HalfCappedDynamicLogLog4.class){
				final CardStats cs=((HalfCappedDynamicLogLog4)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==ExpandedDynamicLogLog8.class){
				final CardStats cs=((ExpandedDynamicLogLog8)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==ExpandedDynamicLogLog9.class){
				final CardStats cs=((ExpandedDynamicLogLog9)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==CompressedDynamicLogLog3.class){
				final CardStats cs=((CompressedDynamicLogLog3)est).consumeLastSummarized();
				accLdlcBase(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti, NUM_BASE);
				accHldlc(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==BankedCompressedDynamicLogLog3.class){
				final CardStats cs=((BankedCompressedDynamicLogLog3)est).consumeLastSummarized();
				accLdlcBase(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti, NUM_BASE);
				accHldlc(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==ArithmeticUltraDynamicLogLog32.class){
				final CardStats cs=((ArithmeticUltraDynamicLogLog32)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==ArithmeticUltraDynamicLogLog33.class){
				final CardStats cs=((ArithmeticUltraDynamicLogLog33)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==ArithmeticVariableDynamicLogLog32.class){
				final CardStats cs=((ArithmeticVariableDynamicLogLog32)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==ArithmeticVariableDynamicLogLog36.class){
				final CardStats cs=((ArithmeticVariableDynamicLogLog36)est).consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
			}else if(cls==ArithmeticVariableDynamicLogLog34.class){
				final ArithmeticVariableDynamicLogLog34 c=(ArithmeticVariableDynamicLogLog34)est;
				final CardStats cs=c.consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
				accVwMean(cs, est, cs.ldlc(), trueCard, ti);
			}else if(cls==ArithmeticVariableDynamicLogLog64.class){
				final ArithmeticVariableDynamicLogLog64 c=(ArithmeticVariableDynamicLogLog64)est;
				final CardStats cs=c.consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
				accVwMean(cs, est, cs.ldlc(), trueCard, ti);
			}else if(cls==ArithmeticVariableLogLog.class){
				final ArithmeticVariableLogLog c=(ArithmeticVariableLogLog)est;
				final CardStats cs=c.consumeLastSummarized();
				accLdlcBase7(cs.ldlc(), cs.dlcSbs(), cs.hc(), 0, cs.hllRaw(), cs.meanHistCF(), cs.hybridPlus2(), trueCard, ti);
				accHldlc7(est, cs.ldlc(), cs.hybridPlus2(), trueCard, ti);
				accVwMean(cs, est, cs.ldlc(), trueCard, ti);
			}else if(cls==ErtlULL.class){
				final double fgra=((ErtlULL)est).fgraEstimatePublic();
				final double lerr=(fgra>0 ? (fgra-trueCard)/(double)trueCard : -1.0);
				ldlcErr[ti][3]+=lerr; ldlcAbsErr[ti][3]+=Math.abs(lerr); ldlcSqErr[ti][3]+=lerr*lerr;
			}else if(cls==ProtoLogLog16c.class){
				final double[] ldlcR=((ProtoLogLog16c)est).ldlcEstimate();
				if(ldlcR!=null){accLdlcFromArray(ldlcR, trueCard, ti);}
			}else if(cls==TwinTailLogLog.class){
				final double[] ldlcR=((TwinTailLogLog)est).ldlcEstimate();
				accLdlcFromArray(ldlcR, trueCard, ti);
				accHldlcFromArray(est, ldlcR, trueCard, ti);
				if(ldlcR.length>8){accSingle(ldlcR[8], trueCard, ti, DDLCalibrationDriver2.MEAN16_IDX);}
				if(ldlcR.length>9){accSingle(ldlcR[9], trueCard, ti, DDLCalibrationDriver2.DUALLC_IDX);}
			}else if(cls==TwinTailLogLog3.class){
				final double[] ldlcR=((TwinTailLogLog3)est).ldlcEstimate();
				accLdlcFromArray(ldlcR, trueCard, ti);
				accHldlcFromArray(est, ldlcR, trueCard, ti);
				if(ldlcR.length>8){accSingle(ldlcR[8], trueCard, ti, DDLCalibrationDriver2.MEAN16_IDX);}
				if(ldlcR.length>9){accSingle(ldlcR[9], trueCard, ti, DDLCalibrationDriver2.DUALLC_IDX);}
			}

			/* LC_noMicro and SBS_noMicro from rawEstimates */
			{
				final int lcIdx=AbstractCardStats.LC_NOMICRO_IDX;
				final int sbsIdx=AbstractCardStats.SBS_NOMICRO_IDX;
				if(lcIdx<estimates.length){accSingle(estimates[lcIdx], trueCard, ti, 7);}
				if(sbsIdx<estimates.length){accSingle(estimates[sbsIdx], trueCard, ti, 8);}
			}
		}

		/*--------------------------------------------------------------*/
		/*----------------      Accumulator Helpers     ----------------*/
		/*--------------------------------------------------------------*/

		private void accLdlcBase(double ldlc, double dlcSbs, double hc, double fgra,
				double hllRaw, double meanHistCF, double hybridPlus2,
				int trueCard, int ti, int count){
			final double[] vals={ldlc, dlcSbs, hc, fgra, hllRaw, meanHistCF, hybridPlus2};
			for(int e=0; e<count; e++){
				final double v=vals[e];
				final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
				ldlcErr[ti][e]+=lerr; ldlcAbsErr[ti][e]+=Math.abs(lerr); ldlcSqErr[ti][e]+=lerr*lerr;
			}
		}

		private void accLdlcBase7(double ldlc, double dlcSbs, double hc, double fgra,
				double hllRaw, double meanHistCF, double hybridPlus2,
				int trueCard, int ti){
			accLdlcBase(ldlc, dlcSbs, hc, fgra, hllRaw, meanHistCF, hybridPlus2, trueCard, ti, 7);
		}

		private void accHldlc(CardinalityTracker est, double ldlc, double hybridPlus2, int trueCard, int ti){
			final float hw=est.hldlcWeight();
			final double hldlc=hw*ldlc+(1-hw)*hybridPlus2;
			final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
			ldlcErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr;
			ldlcAbsErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=Math.abs(lerr);
			ldlcSqErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr*lerr;
		}

		private void accHldlc7(CardinalityTracker est, double ldlc, double hybridPlus2, int trueCard, int ti){
			final float hw=est.hldlcWeight();
			final double hldlc=hw*ldlc+(1-hw)*hybridPlus2;
			final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
			ldlcErr[ti][11]+=lerr; ldlcAbsErr[ti][11]+=Math.abs(lerr); ldlcSqErr[ti][11]+=lerr*lerr;
		}

		private void accHldlcFromArray(CardinalityTracker est, double[] ldlcR, int trueCard, int ti){
			final float hw=est.hldlcWeight();
			final double hldlc=hw*ldlcR[0]+(1-hw)*ldlcR[7];
			final double lerr=(hldlc>0 ? (hldlc-trueCard)/(double)trueCard : -1.0);
			ldlcErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr;
			ldlcAbsErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=Math.abs(lerr);
			ldlcSqErr[ti][DDLCalibrationDriver2.HLDLC_IDX]+=lerr*lerr;
		}

		private static final int[] LDLC_IDX={0, 1, 2, 4, 5, 6, 7};

		private void accLdlcFromArray(double[] ldlcR, int trueCard, int ti){
			for(int e=0; e<DDLCalibrationDriver2.NUM_EXTENDED_BASE; e++){
				final double v=ldlcR[LDLC_IDX[e]];
				final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
				ldlcErr[ti][e]+=lerr; ldlcAbsErr[ti][e]+=Math.abs(lerr); ldlcSqErr[ti][e]+=lerr*lerr;
			}
		}

		private void accSingle(double v, int trueCard, int ti, int idx){
			final double lerr=(v>0 ? (v-trueCard)/(double)trueCard : -1.0);
			ldlcErr[ti][idx]+=lerr; ldlcAbsErr[ti][idx]+=Math.abs(lerr); ldlcSqErr[ti][idx]+=lerr*lerr;
		}

		private void accVwMean(CardStats cs, CardinalityTracker est, double ldlc, int trueCard, int ti){
			if(CardStats.VW_STATE_TABLE==null){return;}
			final double vwPure=cs.vwMean();
			accSingle(vwPure, trueCard, ti, DDLCalibrationDriver2.MEAN16_IDX);
			double vw=vwPure;
			if(vwPure>0){
				final double sbs=cs.sbs(); final double zoneEst=cs.dlcRaw();
				final double hb0=DDLCalibrationDriver2.VW_BLEND_LO*buckets;
				final double hb1=DDLCalibrationDriver2.VW_BLEND_HI*buckets;
				if(zoneEst<=hb0){vw=sbs;}
				else if(zoneEst<hb1){final double t=Math.log(zoneEst/hb0)/Math.log(hb1/hb0); vw=(1-t)*sbs+t*vwPure;}
			}
			accSingle(vw, trueCard, ti, DDLCalibrationDriver2.VWMEAN_IDX);
			if(vw>0 && ldlc>0){
				final double vldlc=DDLCalibrationDriver2.VLDLC_LDLC*ldlc+DDLCalibrationDriver2.VLDLC_VW*vw;
				final double vlerr=(vldlc-trueCard)/(double)trueCard;
				ldlcErr[ti][DDLCalibrationDriver2.VLDLC_IDX]+=vlerr;
				ldlcAbsErr[ti][DDLCalibrationDriver2.VLDLC_IDX]+=Math.abs(vlerr);
				ldlcSqErr[ti][DDLCalibrationDriver2.VLDLC_IDX]+=vlerr*vlerr;
			}
		}

		/*--------------------------------------------------------------*/
		/*----------------           Fields            ----------------*/
		/*--------------------------------------------------------------*/

		final int estStart, estEnd;
		final int cardinality, numThresholds, buckets, k, numEst;
		final long totalAdds;
		final boolean stopOnSaturation, minRand;
		final String loglogType;
		final long[] seeds, thresholds;
		final int[] lowerArray, upperArray;
		final int lowerMask, lowerShift;

		final double[][] sumErr, sumAbsErr, sumSqErr;
		final double[][] ldlcErr, ldlcAbsErr, ldlcSqErr;
		final double[] occSum;
		final int[] nArr;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Summary Output        ----------------*/
	/*--------------------------------------------------------------*/

	private static void printSummary(String loglogtype, int buckets, int numDDLs,
			int cardinality, float iterations,
			ArrayList<DDLCalibrationDriver.ReportRow> mergedRows,
			int[] nArr, long[] thresholds, int numThresholds,
			double[][] ldlcSumErr, double[][] ldlcSumAbsErr, double[][] ldlcSumSqErr,
			int NUM_EST, int NUM_EXT, String[] EXT_NAMES, long t0){
		final double elapsed=(System.nanoTime()-t0)*1e-9;
		final int rows=mergedRows.size();
		System.err.println();
		System.err.println("=== LC Calibration Summary ===");
		final int actualBk=DDLCalibrationDriver.makeInstance(loglogtype, buckets, 31, 0, 0).actualBuckets();
		System.err.println("Type: "+loglogtype+"  Buckets: "+actualBk+"  DDLs: "+numDDLs
			+"  MaxCard: "+cardinality+"  Rows: "+rows
			+"  Elapsed: "+String.format("%.1f", elapsed)+"s"
			+"  LowComplexity: iter="+iterations);
		System.err.println("--- Log-Weighted and Width-Weighted Avg Absolute Error, Peak, Signed (lower = better) ---");

		final double[] widths=new double[rows];
		widths[0]=Math.max(1, mergedRows.get(0).trueCard);
		for(int i=1; i<rows; i++){widths[i]=Math.max(1, mergedRows.get(i).trueCard-mergedRows.get(i-1).trueCard);}

		final double[] totalAbsErr=new double[NUM_EST], widthWtAbsErr=new double[NUM_EST];
		final double[] countWtAbsErr=new double[NUM_EST], peakAbsErr=new double[NUM_EST];
		final double[] totalSignedErr=new double[NUM_EST], totalCV=new double[NUM_EST];
		double widthWeightSum=0, countWeightSum=0;
		int rowIdx=0;
		for(DDLCalibrationDriver.ReportRow row : mergedRows){
			final double w=widths[rowIdx];
			widthWeightSum+=w; countWeightSum+=row.n;
			for(int e=0; e<NUM_EST; e++){
				final double meanErr=row.sumErr[e]/row.n;
				final double meanAbs=row.sumAbsErr[e]/row.n;
				totalAbsErr[e]+=meanAbs; widthWtAbsErr[e]+=meanAbs*w; countWtAbsErr[e]+=meanAbs*row.n;
				totalSignedErr[e]+=meanErr;
				if(meanAbs>peakAbsErr[e]){peakAbsErr[e]=meanAbs;}
				final double var=row.sumSqErr[e]/row.n-meanErr*meanErr;
				final double denom=Math.abs(1.0+meanErr);
				if(denom>0){totalCV[e]+=Math.sqrt(Math.max(0, var))/denom;}
			}
			rowIdx++;
		}
		System.err.println(String.format("%-12s %-13s %-13s %-13s %-13s %-13s %s",
			"", "LogWtAbsErr", "WidthWtAbsErr", "CountWtAbsErr", "PeakAbsErr", "AvgSignErr", "AvgCV"));
		final String[] EST_NAMES=DDLCalibrationDriver.ESTIMATOR_NAMES;
		for(int e=0; e<Math.min(17, NUM_EST); e++){
			System.err.println(String.format("%-12s %.8f   %.8f   %.8f   %.8f   %+.8f   %.8f",
				EST_NAMES[e], totalAbsErr[e]/rows,
				widthWeightSum>0 ? widthWtAbsErr[e]/widthWeightSum : 0,
				countWeightSum>0 ? countWtAbsErr[e]/countWeightSum : 0,
				peakAbsErr[e], rows>0 ? totalSignedErr[e]/rows : 0,
				rows>0 ? totalCV[e]/rows : 0));
		}

		/* LDLC summary */
		final double[] tWidths=new double[numThresholds];
		long prev=0;
		for(int ti2=0; ti2<numThresholds; ti2++){
			if(nArr[ti2]<1){tWidths[ti2]=0; continue;}
			tWidths[ti2]=Math.max(1, thresholds[ti2]-prev); prev=thresholds[ti2];
		}
		final double[] ldlcTotalAbs=new double[NUM_EXT], ldlcWidthWt=new double[NUM_EXT];
		final double[] ldlcCountWt=new double[NUM_EXT], ldlcPeak=new double[NUM_EXT];
		final double[] ldlcTotalSigned=new double[NUM_EXT];
		double ldlcWwSum=0, ldlcCwSum=0;
		for(int ti2=0; ti2<numThresholds; ti2++){
			if(nArr[ti2]<1){continue;}
			final double w=tWidths[ti2]; ldlcWwSum+=w; ldlcCwSum+=nArr[ti2];
			for(int e=0; e<NUM_EXT; e++){
				final double mAbs=ldlcSumAbsErr[ti2][e]/nArr[ti2];
				ldlcTotalAbs[e]+=mAbs; ldlcWidthWt[e]+=mAbs*w; ldlcCountWt[e]+=mAbs*nArr[ti2];
				ldlcTotalSigned[e]+=ldlcSumErr[ti2][e]/nArr[ti2];
				if(mAbs>ldlcPeak[e]){ldlcPeak[e]=mAbs;}
			}
		}
		for(int e=0; e<NUM_EXT; e++){
			System.err.println(String.format("%-12s %.8f   %.8f   %.8f   %.8f   %+.8f   %.8f",
				EXT_NAMES[e], ldlcTotalAbs[e]/rows,
				ldlcWwSum>0 ? ldlcWidthWt[e]/ldlcWwSum : 0,
				ldlcCwSum>0 ? ldlcCountWt[e]/ldlcCwSum : 0,
				ldlcPeak[e], rows>0 ? ldlcTotalSigned[e]/rows : 0, 0.0));
		}
		System.err.println();
	}

}
