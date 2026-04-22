package cardinality;

import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;

import java.io.*;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

/**
 * Multithreaded LC history CF table trainer.
 * Simulates HyperLogLog with 2-bit history, snapshots at filled-count transitions,
 * and outputs per-state expected distinct counts for each occupancy level.
 *
 * Each thread runs independent trials with its own RNG and accumulator arrays.
 * Results are merged after all threads finish.
 *
 * Usage: trainLCHist.sh buckets=2048 t=4 hbits=2 iters=1m out=table.tsv.gz
 *
 * @author Eru
 * @date March 29, 2026
 */
public class LCHistTrainer {

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact copy of Tools.hash64shift for self-contained simulation. */
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

	/** Maps (rawNlz, hist) to a flat index in the full NLZ x hist grid.
	 *  Grid has (hbits+2) NLZ bins x (1&lt;&lt;hbits) hist slots. */
	static int stateIndex(int rawNlz, int hist, int hbits){
		final int numHist=1<<hbits;
		return Math.min(rawNlz, hbits+1)*numHist+(hist&(numHist-1));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker Thread        ----------------*/
	/*--------------------------------------------------------------*/

	static final class TrainerThread extends Thread{

		TrainerThread(long seed, int numTrials, int buckets, int hbits){
			this.seed=seed;
			this.numTrials=numTrials;
			this.buckets=buckets;
			this.hbits=hbits;
			this.bucketMask=buckets-1;
			this.nStates=(hbits+2)*(1<<hbits);
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
			final int hb=hbits;
			final int hmask=(1<<hb)-1;

			for(int t=0; t<numTrials; t++){
				int[] rawNlzArr=new int[B];
				int[] hist=new int[B];
				int[] distinct=new int[B];
				boolean[] filled=new boolean[B];
				int filledCount=0;

				while(filledCount<B){
					final long number=rng.nextLong();
					final long key=hash64shift(number);
					final int nlz=Long.numberOfLeadingZeros(key);
					final int bucket=(int)(key&bucketMask);
					final int rawNlz=nlz;

					distinct[bucket]++;

					if(!filled[bucket]){
						rawNlzArr[bucket]=rawNlz;
						hist[bucket]=0;
						filled[bucket]=true;
						filledCount++;
						snapshot(filledCount, rawNlzArr, hist, distinct, filled);
					}else{
						int oldNlz=rawNlzArr[bucket];
						if(rawNlz>oldNlz){
							// Tier advance: shift old history and set carry bit.
							int oldH=hist[bucket]&hmask;
							int diff=rawNlz-oldNlz;
							hist[bucket]=((oldH|(1<<hb))>>diff)&hmask;
							rawNlzArr[bucket]=rawNlz;
						}else if(rawNlz<oldNlz){
							int diff=oldNlz-rawNlz;
							if(diff>=1 && diff<=hb){
								hist[bucket]|=1<<(hb-diff);
							}
						}
					}
				}
			}
		}

		void snapshot(int filledCount, int[] rawNlz, int[] hist,
				int[] distinct, boolean[] filled){
			for(int b=0; b<buckets; b++){
				if(filled[b]){
					int si=stateIndex(rawNlz[b], hist[b], hbits);
					collisionSums[filledCount][si]+=distinct[b];
					observations[filledCount][si]++;
				}
			}
		}

		/*--------------------------------------------------------------*/
		/*----------------           Fields            ----------------*/
		/*--------------------------------------------------------------*/

		final long seed;
		final int numTrials;
		final int buckets;
		final int hbits;
		final int bucketMask;
		final int nStates;
		final long[][] collisionSums;
		final long[][] observations;
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Main              ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		{ PreParser pp=new PreParser(args, null, false); args=pp.args; }
		final long t0=System.nanoTime();
		int buckets=2048;
		int hbits=2;
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
			}else if(a.equals("hbits") || a.equals("historybits")){
				hbits=Integer.parseInt(b);
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

		// Ensure buckets is power of 2
		buckets=Integer.highestOneBit(buckets);

		final int numThreads=Math.min(threads, (int)iters);
		System.err.println("LCHistTrainer: buckets="+buckets+" hbits="+hbits+
			" iters="+iters+" threads="+numThreads+" seed="+masterSeed);

		// Create and start worker threads
		final TrainerThread[] workers=new TrainerThread[numThreads];
		for(int t=0; t<numThreads; t++){
			final int threadTrials=(int)(iters/numThreads+(t<iters%numThreads ? 1 : 0));
			workers[t]=new TrainerThread(masterSeed+t, threadTrials, buckets, hbits);
		}
		for(TrainerThread w : workers){w.start();}

		// Wait for all threads
		for(TrainerThread w : workers){
			try{w.join();}catch(InterruptedException e){Thread.currentThread().interrupt();}
			if(!w.success){System.err.println("Warning: a trainer thread did not complete successfully.");}
		}

		// Merge accumulators
		final int nStates=(hbits+2)*(1<<hbits);
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

		// Build valid state list dynamically from hbits.
		// At NLZ bin k (0 through hbits+1), top min(k, hbits) history bits are valid.
		// Total valid states = 3 * (1 << hbits) - 1.
		final int numValid=3*(1<<hbits)-1;
		final int[][] validStates=new int[numValid][2]; // [i] = {gridIndex, minDistinct}
		final String[] stateNames=new String[numValid];
		int vi=0;
		for(int bin=0; bin<=hbits+1; bin++){
			final int validSlots=Math.min(bin, hbits);
			final int invalidBits=hbits-validSlots;
			final int numHist=1<<validSlots;
			final String binLabel=(bin>hbits) ? (bin+"+") : String.valueOf(bin);
			for(int j=0; j<numHist; j++){
				final int hist=j<<invalidBits; // valid history value
				final int gridIdx=bin*(1<<hbits)+hist;
				// minDistinct = 1 + number of set bits in valid portion
				final int minDistinct=1+Integer.bitCount(j);
				validStates[vi][0]=gridIdx;
				validStates[vi][1]=minDistinct;
				// Name: "bin.histbits" (e.g. "2.10", "3+.011")
				StringBuilder name=new StringBuilder(binLabel).append('.');
				for(int b=validSlots-1; b>=0; b--){
					name.append((j>>>b)&1);
				}
				if(validSlots==0){name.append('0');}
				stateNames[vi]=name.toString();
				vi++;
			}
		}

		// Verify: count zeros in full grid at final row (should equal nStates - numValid)
		int zeros=0;
		for(int s=0; s<nStates; s++){
			if(totalObs[buckets][s]==0){zeros++;}
		}
		final int expectedZeros=nStates-numValid;
		System.err.println("Zeros in final row: "+zeros+" (expected "+expectedZeros+")");

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
			pw.println("#HistBits\t"+hbits);
			pw.println("#Trials\t"+iters);

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

		System.err.println("Output: "+out);
	}
}
