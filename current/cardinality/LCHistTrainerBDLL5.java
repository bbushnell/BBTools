package cardinality;

import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;

import java.io.*;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

/**
 * SBS table trainer for BankedDynamicLogLog5.
 * Simulates the full BDLL5 mechanic: 6 buckets per word, 3-bit stored,
 * 2-bit history, 2-bit bank exponent, bank promotion, global promotion.
 * Uses Long.remainderUnsigned for non-power-of-2 bucket selection.
 *
 * @author Brian Bushnell, Chloe
 * @date April 2026
 */
public class LCHistTrainerBDLL5 {

	/*--------------------------------------------------------------*/
	/*----------------            Main              ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		{PreParser pp=new PreParser(args, null, false); args=pp.args;}
		final long t0=System.nanoTime();
		int buckets=2048;
		long iters=100000;
		int threads=Shared.threads();
		long masterSeed=42;
		String out="stdout.txt";
		float promoteFrac=0.004f;

		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<0){continue;}
			final String a=arg.substring(0, eq).toLowerCase();
			final String b=arg.substring(eq+1);
			if(a.equals("buckets")||a.equals("b")){buckets=Parse.parseIntKMG(b);}
			else if(a.equals("iters")||a.equals("trials")||a.equals("i")){iters=Parse.parseIntKMG(b);}
			else if(a.equals("threads")||a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("seed")){masterSeed=Long.parseLong(b);}
			else if(a.equals("out")){out=b;}
			else if(a.equals("promotefrac")||a.equals("pf")){promoteFrac=Float.parseFloat(b);}
		}

		final int words=Math.max(1, (buckets+3)/6);
		final int modBuckets=words*6;
		PROMOTE_THRESHOLD=(promoteFrac>0 ? (int)(modBuckets*promoteFrac) : 0);

		final int numThreads=Math.min(threads, (int)iters);
		System.err.println("LCHistTrainerBDLL5: modBuckets="+modBuckets+" words="+words
			+" hbits="+HBITS+" iters="+iters+" threads="+numThreads
			+" promoteFrac="+promoteFrac+" promoteThreshold="+PROMOTE_THRESHOLD);

		final TrainerThread[] workers=new TrainerThread[numThreads];
		for(int t=0; t<numThreads; t++){
			final int threadTrials=(int)(iters/numThreads+(t<iters%numThreads ? 1 : 0));
			workers[t]=new TrainerThread(masterSeed+t, threadTrials, modBuckets, words);
		}
		for(TrainerThread w : workers){w.start();}
		for(TrainerThread w : workers){
			try{w.join();}catch(InterruptedException e){Thread.currentThread().interrupt();}
			if(!w.success){System.err.println("Warning: thread failed.");}
		}

		final long[][] totalSums=new long[modBuckets+1][NSTATES];
		final long[][] totalObs=new long[modBuckets+1][NSTATES];
		for(TrainerThread w : workers){
			for(int f=1; f<=modBuckets; f++){
				for(int s=0; s<NSTATES; s++){
					totalSums[f][s]+=w.collisionSums[f][s];
					totalObs[f][s]+=w.observations[f][s];
				}
			}
		}

		final double elapsed=(System.nanoTime()-t0)/1e9;
		System.err.println("Training done: "+String.format("%.1f", elapsed)+"s");

		final int numValid=3*(1<<HBITS)-1;
		final int[][] validStates=new int[numValid][2];
		final String[] stateNames=new String[numValid];
		int vi=0;
		for(int bin=0; bin<=HBITS+1; bin++){
			final int validSlots=Math.min(bin, HBITS);
			final int invalidBits=HBITS-validSlots;
			final int numHist=1<<validSlots;
			final String binLabel=(bin>HBITS) ? (bin+"+") : String.valueOf(bin);
			for(int j=0; j<numHist; j++){
				final int h=j<<invalidBits;
				final int gridIdx=bin*(1<<HBITS)+h;
				final int minDistinct=1+Integer.bitCount(j);
				validStates[vi][0]=gridIdx;
				validStates[vi][1]=minDistinct;
				StringBuilder name=new StringBuilder(binLabel).append('.');
				for(int b2=validSlots-1; b2>=0; b2--){name.append((j>>>b2)&1);}
				if(validSlots==0){name.append('0');}
				stateNames[vi]=name.toString();
				vi++;
			}
		}

		int zeros=0;
		for(int s=0; s<NSTATES; s++){if(totalObs[modBuckets][s]==0){zeros++;}}
		System.err.println("Zeros in final row: "+zeros+" (expected "+(NSTATES-numValid)+")");

		double[][] avg=new double[modBuckets+1][NSTATES];
		for(int f=1; f<=modBuckets; f++){
			for(int s=0; s<NSTATES; s++){
				if(totalObs[f][s]>0){avg[f][s]=(double)totalSums[f][s]/totalObs[f][s];}
			}
		}

		try{
			OutputStream os;
			if(out.equals("stdout.txt")){os=System.out;}
			else if(out.endsWith(".gz")){os=new GZIPOutputStream(new FileOutputStream(out));}
			else{os=new FileOutputStream(out);}
			PrintWriter pw=new PrintWriter(new OutputStreamWriter(os));

			pw.println("#LC_HISTORY_TABLE");
			pw.println("#Buckets\t"+modBuckets);
			pw.println("#HistBits\t"+HBITS);
			pw.println("#Trials\t"+iters);

			StringBuilder hdr=new StringBuilder("#Filled");
			for(String name : stateNames){hdr.append('\t').append(name);}
			pw.println(hdr);

			for(int f=1; f<=modBuckets; f++){
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
		}catch(IOException e){e.printStackTrace();}

		System.err.println("Output: "+out);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker Thread        ----------------*/
	/*--------------------------------------------------------------*/

	static final class TrainerThread extends Thread {

		/** Construct a trainer thread for a single trial batch. */
		TrainerThread(final long seed, final int numTrials, final int modBuckets, final int words){
			this.seed=seed;
			this.numTrials=numTrials;
			this.modBuckets=modBuckets;
			this.words=words;
			nStates=NSTATES;
			collisionSums=new long[modBuckets+1][nStates];
			observations=new long[modBuckets+1][nStates];
		}

		@Override public void run(){
			try{runInner(); success=true;}
			catch(Exception e){e.printStackTrace();}
		}

		/** Run all trials, accumulating collision counts into collisionSums/observations. */
		void runInner(){
			final Random rng=new Random(seed);
			final int B=modBuckets;
			final int W=words;

			for(int t=0; t<numTrials; t++){
				int[] stored=new int[B];
				int[] hist=new int[B];
				int[] bank=new int[W];
				int[] distinct=new int[B];
				boolean[] filled=new boolean[B];
				int filledCount=0;
				int globalNLZ=-1;
				int minZeroCount=B;

				while(filledCount<B){
					final long number=rng.nextLong();
					final long key=hash64shift(number);
					final int nlz=Long.numberOfLeadingZeros(key);
					final int bucket=(int)(Long.remainderUnsigned(key, B));
					final int wordIdx=bucket/6;
					int bk=bank[wordIdx];
					int localRelNlz=nlz-globalNLZ-1-bk;

					distinct[bucket]++;

					// Bank promotion
					if(localRelNlz>=7 && bk<3){
						if(canPromote(stored, wordIdx)){
							promoteBank(stored, bank, wordIdx);
							bk=bank[wordIdx];
							localRelNlz=nlz-globalNLZ-1-bk;
						}
					}

					final int newStored=Math.min(localRelNlz+1, 7);
					final int oldStored=stored[bucket];
					final int oldHist=hist[bucket];

					if(newStored<=0){
						// Sub-floor hist update
						if(oldStored>0){
							final int delta=oldStored-newStored;
							if(delta==1||delta==2){hist[bucket]=oldHist|(1<<(2-delta));}
						}
						continue;
					}

					if(newStored==oldStored){continue;}
					if(newStored<oldStored){
						final int delta=oldStored-newStored;
						if(delta==1||delta==2){hist[bucket]=oldHist|(1<<(2-delta));}
						continue;
					}

					// Tier advance: carry-shift history
					final int newHist;
					if(oldStored==0){
						newHist=0;
					}else{
						final int delta=newStored-oldStored;
						newHist=(delta>=3) ? 0 : (((1<<2)|oldHist)>>>delta)&0x3;
					}
					stored[bucket]=newStored;
					hist[bucket]=newHist;

					final boolean wasFilled=filled[bucket];
					if(!wasFilled){
						filled[bucket]=true;
						filledCount++;
						snapshot(filledCount, stored, hist, bank, distinct, filled, globalNLZ);
					}

					// Global promotion
					if(oldStored==0 && bk==0){
						minZeroCount--;
						while(minZeroCount<=PROMOTE_THRESHOLD && globalNLZ<63){
							globalNLZ++;
							minZeroCount=countAndDecrement(stored, hist, bank, B, W);
						}
					}
				}
			}
		}

		/** Returns true if all 6 registers in the word are non-empty (stored > 0). */
		boolean canPromote(final int[] stored, final int wordIdx){
			final int base=wordIdx*6;
			for(int b=0; b<6; b++){
				if(stored[base+b]==0){return false;}
			}
			return true;
		}

		/** Decrement all 6 registers in the word by 1 and increment its bank exponent. */
		void promoteBank(final int[] stored, final int[] bank, final int wordIdx){
			final int base=wordIdx*6;
			for(int b=0; b<6; b++){stored[base+b]--;}
			bank[wordIdx]++;
		}

		/**
		 * Decrement all registers after a global floor advance.
		 * Words with bank > 0 absorb by decrementing bank only.
		 * Words at bank=0 decrement each stored value.
		 * @return New minZeroCount (buckets at stored=0 eligible for next advance)
		 */
		int countAndDecrement(final int[] stored, final int[] hist, final int[] bank, final int B, final int W){
			int count=0;
			for(int w=0; w<W; w++){
				if(bank[w]>0){
					bank[w]--;
					if(bank[w]==0){
						final int base=w*6;
						for(int b=0; b<6; b++){
							if(stored[base+b]==0){count++;}
						}
					}
				}else{
					final int base=w*6;
					for(int b=0; b<6; b++){
						if(stored[base+b]>0){
							stored[base+b]--;
							if(stored[base+b]==0){count++;}
						}else{
							count++;
						}
					}
				}
			}
			return count;
		}

		/** Record collision and observation counts for all filled buckets at this fill level. */
		void snapshot(final int filledCount, final int[] stored, final int[] hist, final int[] bank,
				final int[] distinct, final boolean[] filled, final int globalNLZ){
			for(int b=0; b<modBuckets; b++){
				if(filled[b]){
					final int bk=bank[b/6];
					final int s=stored[b];
					final int absNlz=s+globalNLZ+bk;
					if(absNlz<0){continue;}
					final int h=(absNlz==0) ? 0 : hist[b];
					final int si=stateIndex(absNlz, h, HBITS);
					collisionSums[filledCount][si]+=distinct[b];
					observations[filledCount][si]++;
				}
			}
		}

		final long seed;
		final int numTrials, modBuckets, words, nStates;
		final long[][] collisionSums, observations;
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Compute compressed-tier state index from absNlz and history bits. */
	static int stateIndex(final int rawNlz, final int hist, final int hbits){
		final int numHist=1<<hbits;
		return Math.min(rawNlz, hbits+1)*numHist+(hist&(numHist-1));
	}

	/** BBTools hash64shift: avalanche a 64-bit key. */
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

	/** Number of history bits per bucket. */
	static final int HBITS=2;
	/** Total number of (tier, hist) state slots. */
	static final int NSTATES=(HBITS+2)*(1<<HBITS);
	/** Minimum minZeroCount before triggering a global floor advance. */
	static int PROMOTE_THRESHOLD=0;

}
