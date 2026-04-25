package cardinality;

import parse.Parse;
import parse.PreParser;
import shared.Shared;

import java.io.*;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

/**
 * SBS table trainer for TwinTailLogLog3 (3-bit tails).
 * Simulates TTLL3's symmetric encoding (two 3-bit tails selected by histBit)
 * to build per-state expected distinct counts at each occupancy level.
 *
 * State = (nlzBin, combined_h) where combined_h = (h1&lt;&lt;3)|h0.
 * 3 NLZ bins (0, 1, 2+) x 64 history patterns = 192 grid cells.
 * Bin 0: 3 reachable (MSB-only). Bin 1: 12 (MSB+mid). Bin 2+: 48. Total: 63.
 *
 * @author Chloe
 * @date April 25, 2026
 */
public class LCHistTrainerTTLL3 {

	static final int HBITS=6;
	static final int TAIL_BITS=3;
	static final int TAIL_MASK=0x7;
	static final int NUM_BINS=3;
	static final int NUM_CH=64;
	static final int MSB_MASK=0x24;     // h0 bit2 + h1 bit2
	static final int LSB_MASK=0x09;     // h0 bit0 + h1 bit0
	static final int MID_LSB_MASK=0x1B; // h0 bits0,1 + h1 bits0,1

	static final int SAMPLE_ENTRY=0;
	static final int SAMPLE_ALL=1;
	static final int SAMPLE_BOTH=2;
	static final int SAMPLE_PCTILE=3;
	static final String[] SAMPLE_NAMES={"entry","all","both","percentile"};

	static final int AVG_LIN=0;
	static final int AVG_GEO=1;
	static final int AVG_HARM=2;
	static final String[] AVG_NAMES={"lin","geo","harm"};

	static long hash64shift(long key){
		key=(~key)+(key<<21);
		key=key^(key>>>24);
		key=(key+(key<<3))+(key<<8);
		key=key^(key>>>14);
		key=(key+(key<<2))+(key<<4);
		key=key^(key>>>28);
		key=key+(key<<31);
		return key;
	}

	static int stateIndex(int rawNlz, int combinedH){
		return Math.min(rawNlz, NUM_BINS-1)*NUM_CH+(combinedH&0x3F);
	}

	static boolean isReachable(int bin, int ch){
		if((ch&MSB_MASK)==0){return false;}
		if(bin==0 && (ch&MID_LSB_MASK)!=0){return false;}
		if(bin==1 && (ch&LSB_MASK)!=0){return false;}
		return true;
	}

	static final class TrainerThread extends Thread {

		TrainerThread(long seed, int numTrials, int buckets, int sampleMode, int avgMode){
			this.seed=seed;
			this.numTrials=numTrials;
			this.buckets=buckets;
			this.sampleMode=sampleMode;
			this.avgMode=avgMode;
			this.nStates=NUM_BINS*NUM_CH;
			accumSums=new double[buckets+1][nStates];
			observations=new long[buckets+1][nStates];
		}

		@Override
		public void run(){
			try{ runInner(); success=true; }
			catch(Exception e){ e.printStackTrace(); }
		}

		void runInner(){
			final Random rng=new Random(seed);
			final int B=buckets;
			final int sm=sampleMode;

			for(int t=0; t<numTrials; t++){
				int[] rawNlzArr=new int[B];
				int[] ch=new int[B];
				int[] distinct=new int[B];
				boolean[] filled=new boolean[B];
				int filledCount=0;
				long totalAdds=0;
				long nextPctile=1;

				while(filledCount<B){
					final long number=rng.nextLong();
					final long key=hash64shift(number);
					final int nlz=Long.numberOfLeadingZeros(key);
					final int bucket=(int)(Long.remainderUnsigned(key, B));
					final int histBit=(int)(Long.divideUnsigned(key, B)&1);

					distinct[bucket]++;
					totalAdds++;
					boolean wasFill=false;

					if(!filled[bucket]){
						rawNlzArr[bucket]=nlz;
						final int bitPos=(histBit==0) ? (TAIL_BITS-1) : (TAIL_BITS-1+TAIL_BITS);
						ch[bucket]=(1<<bitPos);
						filled[bucket]=true;
						filledCount++;
						wasFill=true;
					}else{
						final int oldNlz=rawNlzArr[bucket];
						final int delta=nlz-oldNlz;

						if(delta<-(TAIL_BITS-1)){
							// ignore
						}else if(delta<=0){
							final int bitPos=(TAIL_BITS-1)+delta;
							final int regBitPos=(histBit==0) ? bitPos : (bitPos+TAIL_BITS);
							ch[bucket]|=(1<<regBitPos);
						}else{
							final int shiftAmt=Math.min(delta, TAIL_BITS);
							int h0=(ch[bucket]&TAIL_MASK)>>>shiftAmt;
							int h1=((ch[bucket]>>>TAIL_BITS)&TAIL_MASK)>>>shiftAmt;
							if(histBit==0){h0|=(1<<(TAIL_BITS-1));}
							else          {h1|=(1<<(TAIL_BITS-1));}
							ch[bucket]=(h1<<TAIL_BITS)|h0;
							rawNlzArr[bucket]=nlz;
						}
					}

					if(sm==SAMPLE_ENTRY){
						if(wasFill){snapshot(filledCount, rawNlzArr, ch, distinct, filled);}
					}else if(sm==SAMPLE_ALL){
						snapshot(filledCount, rawNlzArr, ch, distinct, filled);
					}else if(sm==SAMPLE_BOTH){
						if(wasFill){snapshot(filledCount, rawNlzArr, ch, distinct, filled);}
						snapshot(filledCount, rawNlzArr, ch, distinct, filled);
					}else{
						if(totalAdds>=nextPctile){
							snapshot(filledCount, rawNlzArr, ch, distinct, filled);
							nextPctile=Math.max(nextPctile+1, (long)Math.ceil(nextPctile*1.01));
						}
					}
				}
			}
		}

		void snapshot(int filledCount, int[] rawNlz, int[] combinedH,
				int[] distinct, boolean[] filled){
			final int am=avgMode;
			for(int b=0; b<buckets; b++){
				if(filled[b]){
					final int si=stateIndex(rawNlz[b], combinedH[b]);
					final double d=distinct[b];
					if(am==AVG_LIN){
						accumSums[filledCount][si]+=d;
					}else if(am==AVG_GEO){
						accumSums[filledCount][si]+=Math.log(d);
					}else{
						accumSums[filledCount][si]+=1.0/d;
					}
					observations[filledCount][si]++;
				}
			}
		}

		final long seed;
		final int numTrials;
		final int buckets;
		final int sampleMode;
		final int avgMode;
		final int nStates;
		final double[][] accumSums;
		final long[][] observations;
		boolean success=false;
	}

	public static void main(String[] args){
		{ PreParser pp=new PreParser(args, null, false); args=pp.args; }
		final long t0=System.nanoTime();
		int buckets=1536;
		long iters=100000;
		int threads=Shared.threads();
		long masterSeed=42;
		String out="stdout.txt";
		int sampleMode=SAMPLE_ENTRY;
		int avgMode=AVG_LIN;

		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<0){continue;}
			final String a=arg.substring(0, eq).toLowerCase();
			final String b=arg.substring(eq+1).toLowerCase();
			if(a.equals("buckets") || a.equals("b")){
				buckets=Parse.parseIntKMG(b);
			}else if(a.equals("iters") || a.equals("trials") || a.equals("i")){
				iters=Parse.parseIntKMG(b);
			}else if(a.equals("threads") || a.equals("t")){
				threads=Integer.parseInt(b);
			}else if(a.equals("seed")){
				masterSeed=Long.parseLong(b);
			}else if(a.equals("out")){
				out=arg.substring(eq+1);
			}else if(a.equals("sample") || a.equals("s")){
				if(b.equals("entry")){sampleMode=SAMPLE_ENTRY;}
				else if(b.equals("all")){sampleMode=SAMPLE_ALL;}
				else if(b.equals("both")){sampleMode=SAMPLE_BOTH;}
				else if(b.startsWith("pct") || b.startsWith("percent")){sampleMode=SAMPLE_PCTILE;}
			}else if(a.equals("avg") || a.equals("averaging")){
				if(b.equals("lin") || b.equals("linear")){avgMode=AVG_LIN;}
				else if(b.equals("geo") || b.equals("geometric")){avgMode=AVG_GEO;}
				else if(b.equals("harm") || b.equals("harmonic")){avgMode=AVG_HARM;}
			}
		}

		final int numThreads=Math.min(threads, (int)iters);
		System.err.println("LCHistTrainerTTLL3: buckets="+buckets+" bins="+NUM_BINS
			+" iters="+iters+" threads="+numThreads
			+" sample="+SAMPLE_NAMES[sampleMode]+" avg="+AVG_NAMES[avgMode]);

		final TrainerThread[] workers=new TrainerThread[numThreads];
		for(int t=0; t<numThreads; t++){
			final int threadTrials=(int)(iters/numThreads+(t<iters%numThreads ? 1 : 0));
			workers[t]=new TrainerThread(masterSeed+t, threadTrials, buckets, sampleMode, avgMode);
		}
		for(TrainerThread w : workers){w.start();}
		for(TrainerThread w : workers){
			try{w.join();}catch(InterruptedException e){Thread.currentThread().interrupt();}
			if(!w.success){System.err.println("Warning: trainer thread failed.");}
		}

		final int nStates=NUM_BINS*NUM_CH;
		final double[][] totalSums=new double[buckets+1][nStates];
		final long[][] totalObs=new long[buckets+1][nStates];
		for(TrainerThread w : workers){
			for(int f=1; f<=buckets; f++){
				for(int s=0; s<nStates; s++){
					totalSums[f][s]+=w.accumSums[f][s];
					totalObs[f][s]+=w.observations[f][s];
				}
			}
		}

		final double elapsed=(System.nanoTime()-t0)/1e9;
		System.err.println("Training done: "+String.format("%.1f", elapsed)+"s");

		int numValid=0;
		for(int bin=0; bin<NUM_BINS; bin++){
			for(int ch=0; ch<NUM_CH; ch++){
				if(isReachable(bin, ch)){numValid++;}
			}
		}
		final int[][] validStates=new int[numValid][2];
		final String[] stateNames=new String[numValid];
		int vi=0;
		for(int bin=0; bin<NUM_BINS; bin++){
			final String binLabel=(bin>=NUM_BINS-1) ? (bin+"+") : String.valueOf(bin);
			for(int ch=0; ch<NUM_CH; ch++){
				if(!isReachable(bin, ch)){continue;}
				validStates[vi][0]=bin*NUM_CH+ch;
				validStates[vi][1]=1;
				final int h1=(ch>>>3)&0x7;
				final int h0=ch&0x7;
				stateNames[vi]=binLabel+"."+
					((h1>>>2)&1)+((h1>>>1)&1)+(h1&1)+
					((h0>>>2)&1)+((h0>>>1)&1)+(h0&1);
				vi++;
			}
		}

		int zeros=0;
		for(int s=0; s<nStates; s++){
			if(totalObs[buckets][s]==0){zeros++;}
		}
		System.err.println("Valid states: "+numValid+" (zeros in final row: "+zeros
			+", expected "+(nStates-numValid)+")");

		double[][] avg=new double[buckets+1][nStates];
		for(int f=1; f<=buckets; f++){
			for(int s=0; s<nStates; s++){
				if(totalObs[f][s]>0){
					final double sum=totalSums[f][s];
					final long cnt=totalObs[f][s];
					if(avgMode==AVG_LIN){
						avg[f][s]=sum/cnt;
					}else if(avgMode==AVG_GEO){
						avg[f][s]=Math.exp(sum/cnt);
					}else{
						avg[f][s]=cnt/sum;
					}
				}
			}
		}

		try{
			OutputStream os;
			if(out.equals("stdout.txt")){
				os=System.out;
			}else if(out.endsWith(".gz")){
				os=new GZIPOutputStream(new FileOutputStream(out));
			}else{
				os=new FileOutputStream(out);
			}
			PrintWriter pw=new PrintWriter(new OutputStreamWriter(os));

			pw.println("#LC_HISTORY_TABLE");
			pw.println("#Buckets\t"+buckets);
			pw.println("#HistBits\t"+HBITS);
			pw.println("#Trials\t"+iters);
			pw.println("#Type\tTTLL3");
			pw.println("#Bins\t"+NUM_BINS);
			pw.println("#Sample\t"+SAMPLE_NAMES[sampleMode]);
			pw.println("#Averaging\t"+AVG_NAMES[avgMode]);
			pw.println("#ValidStates\t"+numValid);

			StringBuilder hdr=new StringBuilder("#Filled");
			for(String name : stateNames){hdr.append('\t').append(name);}
			pw.println(hdr);

			for(int f=1; f<=buckets; f++){
				StringBuilder row=new StringBuilder();
				row.append(f);
				for(int v=0; v<numValid; v++){
					int si=validStates[v][0];
					int minVal=validStates[v][1];
					double val=Math.max(avg[f][si], minVal);
					row.append('\t').append(String.format("%.8f", val));
				}
				pw.println(row);
			}

			pw.flush();
			if(os!=System.out){pw.close();}
		}catch(IOException e){
			e.printStackTrace();
		}

		System.err.println("Output: "+out+" ("+numValid+" valid states)");
	}
}
