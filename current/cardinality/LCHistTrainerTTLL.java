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
 * TTLL update rule:
 *   histBit = hash & 1 (selects h0 or h1 tail)
 *   delta = hashNLZ - storedNLZ
 *   delta &lt; -1 : ignore
 *   delta == -1: set LSB of history[histBit]
 *   delta == 0 : set MSB of history[histBit]
 *   delta &gt; 0 : advance NLZ; shift BOTH tails right by min(delta,2);
 *               set MSB of history[histBit]
 *
 * State = combined_h = (h1&lt;&lt;2)|h0, 4 bits, 16 states per NLZ bin.
 * State indexing uses the same stateIndex as LCHistTrainer (flat grid).
 *
 * Usage: java -ea -cp current/ cardinality.LCHistTrainerTTLL
 *        buckets=2048 t=4 iters=1m out=sbsTTLL.tsv.gz
 *
 * @author Chloe
 * @date April 23, 2026
 */
public class LCHistTrainerTTLL {

	static final int HBITS=4;

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

	/** Maps (rawNlz, combined_h) to flat grid index.
	 *  Grid: (HBITS+2) NLZ bins x 16 hist slots. */
	static int stateIndex(int rawNlz, int combinedH){
		return Math.min(rawNlz, HBITS+1)*16+(combinedH&0xF);
	}

	static final class TrainerThread extends Thread {

		TrainerThread(long seed, int numTrials, int buckets){
			this.seed=seed;
			this.numTrials=numTrials;
			this.buckets=buckets;
			this.bucketMask=buckets-1;
			this.nStates=(HBITS+2)*16;
			collisionSums=new long[buckets+1][nStates];
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

			for(int t=0; t<numTrials; t++){
				int[] rawNlzArr=new int[B]; // stored NLZ per bucket
				int[] ch=new int[B];         // combined_h per bucket: (h1<<2)|h0
				int[] distinct=new int[B];
				boolean[] filled=new boolean[B];
				int filledCount=0;

				while(filledCount<B){
					final long number=rng.nextLong();
					final long key=hash64shift(number);
					final int nlz=Long.numberOfLeadingZeros(key);
					final int bucket=(int)(key&bucketMask);
					final int histBit=(int)((key>>>Integer.numberOfTrailingZeros(B))&1);

					distinct[bucket]++;

					if(!filled[bucket]){
						// First element: set NLZ, set MSB of selected tail
						rawNlzArr[bucket]=nlz;
						final int bitPos=(histBit==0) ? 1 : 3;
						ch[bucket]=(1<<bitPos);
						filled[bucket]=true;
						filledCount++;
						snapshot(filledCount, rawNlzArr, ch, distinct, filled);
					}else{
						final int oldNlz=rawNlzArr[bucket];
						final int delta=nlz-oldNlz;

						if(delta<-1){
							// Below range: ignore
						}else if(delta==-1){
							// Set LSB of selected tail
							final int bitPos=(histBit==0) ? 0 : 2;
							ch[bucket]|=(1<<bitPos);
						}else if(delta==0){
							// Set MSB of selected tail
							final int bitPos=(histBit==0) ? 1 : 3;
							ch[bucket]|=(1<<bitPos);
						}else{
							// Advance: shift both tails right, set MSB of selected tail
							final int shiftAmt=Math.min(delta, 2);
							final int oldH=ch[bucket]&0xF;
							final int h0=(oldH&0x3)>>>shiftAmt;
							final int h1=((oldH>>>2)&0x3)>>>shiftAmt;
							ch[bucket]=(h1<<2)|h0;
							final int bitPos=(histBit==0) ? 1 : 3;
							ch[bucket]|=(1<<bitPos);
							rawNlzArr[bucket]=nlz;
						}
					}
				}
			}
		}

		void snapshot(int filledCount, int[] rawNlz, int[] combinedH,
				int[] distinct, boolean[] filled){
			for(int b=0; b<buckets; b++){
				if(filled[b]){
					int si=stateIndex(rawNlz[b], combinedH[b]);
					collisionSums[filledCount][si]+=distinct[b];
					observations[filledCount][si]++;
				}
			}
		}

		final long seed;
		final int numTrials;
		final int buckets;
		final int bucketMask;
		final int nStates;
		final long[][] collisionSums;
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

		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<0){continue;}
			final String a=arg.substring(0, eq).toLowerCase();
			final String b=arg.substring(eq+1);
			if(a.equals("buckets") || a.equals("b")){
				buckets=Parse.parseIntKMG(b);
			}else if(a.equals("iters") || a.equals("trials") || a.equals("i")){
				iters=Parse.parseIntKMG(b);
			}else if(a.equals("threads") || a.equals("t")){
				threads=Integer.parseInt(b);
			}else if(a.equals("seed")){
				masterSeed=Long.parseLong(b);
			}else if(a.equals("out")){
				out=b;
			}
		}

		buckets=Integer.highestOneBit(buckets);

		final int numThreads=Math.min(threads, (int)iters);
		System.err.println("LCHistTrainerTTLL: buckets="+buckets+
			" hbits="+HBITS+" iters="+iters+" threads="+numThreads);

		final TrainerThread[] workers=new TrainerThread[numThreads];
		for(int t=0; t<numThreads; t++){
			final int threadTrials=(int)(iters/numThreads+(t<iters%numThreads ? 1 : 0));
			workers[t]=new TrainerThread(masterSeed+t, threadTrials, buckets);
		}
		for(TrainerThread w : workers){w.start();}
		for(TrainerThread w : workers){
			try{w.join();}catch(InterruptedException e){Thread.currentThread().interrupt();}
			if(!w.success){System.err.println("Warning: trainer thread failed.");}
		}

		// Merge
		final int nStates=(HBITS+2)*16;
		final long[][] totalSums=new long[buckets+1][nStates];
		final long[][] totalObs=new long[buckets+1][nStates];
		for(TrainerThread w : workers){
			for(int f=1; f<=buckets; f++){
				for(int s=0; s<nStates; s++){
					totalSums[f][s]+=w.collisionSums[f][s];
					totalObs[f][s]+=w.observations[f][s];
				}
			}
		}

		final double elapsed=(System.nanoTime()-t0)/1e9;
		System.err.println("Training done: "+String.format("%.1f", elapsed)+"s");

		// Flat grid: all (HBITS+2) bins x 16 states = 96 entries.
		// No invalid-state filtering — TTLL's bit layout doesn't follow
		// standard tier-depth ordering, so all states are potentially reachable.
		final int numBins=HBITS+2; // 6 bins (0..4, 5+)
		final int numValid=numBins*16; // 96 states
		final int[][] validStates=new int[numValid][2];
		final String[] stateNames=new String[numValid];
		int vi=0;
		for(int bin=0; bin<numBins; bin++){
			final String binLabel=(bin>HBITS) ? (bin+"+") : String.valueOf(bin);
			for(int ch2=0; ch2<16; ch2++){
				final int gridIdx=bin*16+ch2;
				validStates[vi][0]=gridIdx;
				validStates[vi][1]=1; // min distinct = 1 for any filled bucket
				stateNames[vi]=binLabel+"."+
					((ch2>>>3)&1)+((ch2>>>2)&1)+((ch2>>>1)&1)+(ch2&1);
				vi++;
			}
		}

		// Compute per-state averages
		double[][] avg=new double[buckets+1][nStates];
		for(int f=1; f<=buckets; f++){
			for(int s=0; s<nStates; s++){
				if(totalObs[f][s]>0){
					avg[f][s]=(double)totalSums[f][s]/totalObs[f][s];
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
