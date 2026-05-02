package cardinality;

import parse.Parse;
import parse.PreParser;
import shared.Shared;

import java.io.*;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

/**
 * SBS table trainer for TwinTailLogLog.
 * Simulates TTLL's symmetric encoding (two 2-bit tails selected by histBit)
 * to build per-state expected distinct counts at each occupancy level.
 *
 * State = (nlzBin, combined_h) where combined_h = (h1&lt;&lt;2)|h0.
 * 3 NLZ bins (0, 1, 2+) x 16 history patterns = 48 grid cells.
 * Bin 0: 3 reachable (MSB-only). Bin 1, 2+: 12 each. Total: 27 valid states.
 *
 * Supports 4 sampling modes x 3 averaging modes = 12 combinations:
 *   sample=entry|all|both|percentile
 *   avg=lin|geo|harm
 *
 * Usage: java -ea -cp current/ cardinality.LCHistTrainerTTLL
 *        buckets=2048 t=4 iters=1m sample=all avg=lin out=sbsTTLL.tsv.gz
 *
 * @author Chloe
 * @date April 23, 2026
 */
public class LCHistTrainerTTLL {

	static int NUM_TAILS=2;
	static int HIST_LEN=2;
	static int HBITS=NUM_TAILS*HIST_LEN;
	static final int NUM_BINS=3;

	static int computeMsbMask(){
		int mask=0;
		for(int ti=0; ti<NUM_TAILS; ti++){mask|=(1<<(ti*HIST_LEN+HIST_LEN-1));}
		return mask;
	}
	static int computeLsbMask(){
		int mask=0;
		for(int ti=0; ti<NUM_TAILS; ti++){mask|=(1<<(ti*HIST_LEN));}
		return mask;
	}
	static int CH_MSB_MASK=computeMsbMask();
	static int CH_LSB_MASK=computeLsbMask();

	static boolean COMPRESSED=false;
	static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	// Sampling modes
	static final int SAMPLE_ENTRY=0;
	static final int SAMPLE_ALL=1;
	static final int SAMPLE_BOTH=2;
	static final int SAMPLE_PCTILE=3;
	static final String[] SAMPLE_NAMES={"entry","all","both","percentile"};

	// Averaging modes
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
		return Math.min(rawNlz, NUM_BINS-1)*(1<<HBITS)+(combinedH&((1<<HBITS)-1));
	}

	static boolean isReachable(int bin, int ch){
		if((ch&CH_MSB_MASK)==0){return false;}
		if(bin==0 && (ch&CH_LSB_MASK)!=0){return false;}
		return true;
	}

	static final class TrainerThread extends Thread {

		TrainerThread(long seed, int numTrials, int buckets, int sampleMode, int avgMode){
			this.seed=seed;
			this.numTrials=numTrials;
			this.buckets=buckets;
			this.bucketMask=buckets-1;
			this.sampleMode=sampleMode;
			this.avgMode=avgMode;
			this.nStates=NUM_BINS*(1<<HBITS);
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
					final int nlz;
					if(COMPRESSED){
						final int rawNlz=Long.numberOfLeadingZeros(key);
						final int mBits;
						if(rawNlz>=43){mBits=(int)((key<<(rawNlz-42))&0xFFFFF);}
						else{mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);}
						nlz=(2*rawNlz+((mBits>=MANTISSA_THRESHOLD)?1:0))/3;
					}else{
						nlz=Long.numberOfLeadingZeros(key);
					}
					final int bucket=(int)(key&bucketMask);
					final int tailIdx=(int)((key>>>Integer.numberOfTrailingZeros(B))&(NUM_TAILS-1));
					final int hl=HIST_LEN;
					final int tailMask=(1<<hl)-1;

					distinct[bucket]++;
					totalAdds++;
					boolean wasFill=false;

					if(!filled[bucket]){
						rawNlzArr[bucket]=nlz;
						ch[bucket]=(1<<(tailIdx*hl+hl-1));
						filled[bucket]=true;
						filledCount++;
						wasFill=true;
					}else{
						final int oldNlz=rawNlzArr[bucket];
						final int delta=nlz-oldNlz;

						if(delta<-(hl-1)){
							// ignore
						}else if(delta<=0){
							ch[bucket]|=(1<<(tailIdx*hl+(hl-1)+delta));
						}else{
							final int shiftAmt=Math.min(delta, hl);
							int newH=0;
							for(int ti=0; ti<NUM_TAILS; ti++){
								int tail=((ch[bucket]>>>(ti*hl))&tailMask)>>>shiftAmt;
								if(ti==tailIdx){tail|=(1<<(hl-1));}
								newH|=(tail<<(ti*hl));
							}
							ch[bucket]=newH;
							rawNlzArr[bucket]=nlz;
						}
					}

					// Sampling dispatch
					if(sm==SAMPLE_ENTRY){
						if(wasFill){snapshot(filledCount, rawNlzArr, ch, distinct, filled);}
					}else if(sm==SAMPLE_ALL){
						snapshot(filledCount, rawNlzArr, ch, distinct, filled);
					}else if(sm==SAMPLE_BOTH){
						if(wasFill){snapshot(filledCount, rawNlzArr, ch, distinct, filled);}
						snapshot(filledCount, rawNlzArr, ch, distinct, filled);
					}else{// SAMPLE_PCTILE
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
					}else{// AVG_HARM
						accumSums[filledCount][si]+=1.0/d;
					}
					observations[filledCount][si]++;
				}
			}
		}

		final long seed;
		final int numTrials;
		final int buckets;
		final int bucketMask;
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
		int buckets=2048;
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
			}else if(a.equals("compressed")){
				COMPRESSED=Parse.parseBoolean(b);
			}else if(a.equals("numtails") || a.equals("tails")){
				NUM_TAILS=Integer.parseInt(b);
				HBITS=NUM_TAILS*HIST_LEN;
				CH_MSB_MASK=computeMsbMask();
				CH_LSB_MASK=computeLsbMask();
			}else if(a.equals("histlen")){
				HIST_LEN=Integer.parseInt(b);
				HBITS=NUM_TAILS*HIST_LEN;
				CH_MSB_MASK=computeMsbMask();
				CH_LSB_MASK=computeLsbMask();
			}else if(a.equals("out")){
				out=arg.substring(eq+1); // preserve case for filename
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

		buckets=Integer.highestOneBit(buckets);

		final int numThreads=Math.min(threads, (int)iters);
		System.err.println("LCHistTrainerTTLL: buckets="+buckets+" bins="+NUM_BINS
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

		// Merge
		final int numCH=1<<HBITS;
		final int nStates=NUM_BINS*numCH;
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

		// Build valid state list
		int numValid=0;
		for(int bin=0; bin<NUM_BINS; bin++){
			for(int ch=0; ch<numCH; ch++){
				if(isReachable(bin, ch)){numValid++;}
			}
		}
		final int[][] validStates=new int[numValid][2];
		final String[] stateNames=new String[numValid];
		int vi=0;
		for(int bin=0; bin<NUM_BINS; bin++){
			final String binLabel=(bin>=NUM_BINS-1) ? (bin+"+") : String.valueOf(bin);
			for(int ch=0; ch<numCH; ch++){
				if(!isReachable(bin, ch)){continue;}
				validStates[vi][0]=bin*numCH+ch;
				validStates[vi][1]=1;
				StringBuilder bits=new StringBuilder();
				for(int b=HBITS-1; b>=0; b--){bits.append((ch>>>b)&1);}
				stateNames[vi]=binLabel+"."+bits;
				vi++;
			}
		}

		int zeros=0;
		for(int s=0; s<nStates; s++){
			if(totalObs[buckets][s]==0){zeros++;}
		}
		System.err.println("Valid states: "+numValid+" (zeros in final row: "+zeros
			+", expected "+(nStates-numValid)+")");

		// Compute per-state averages
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
					}else{// AVG_HARM
						avg[f][s]=cnt/sum;
					}
				}
			}
		}

		// Write output
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
			pw.println("#Type\tTTLL");
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
