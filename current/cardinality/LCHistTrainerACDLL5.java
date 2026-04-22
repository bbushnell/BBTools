package cardinality;

import parse.Parse;
import parse.PreParser;
import shared.Shared;

import java.io.*;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

/**
 * LC history CF table trainer for ACDLL5 (compressed tiers, 2-bit history).
 * Tier = (2*rawNlz + mantissa) / 3, with mantissa = key_fraction >= (2-sqrt(2)).
 * History update: standard 2-bit carry-shift ((oldH | carry) >> delta) & hmask.
 * No overflow modeling — with 10 exponent states, overflow is ~1/11585 (negligible).
 * No floor advancement — same reasoning.
 * Modulo bucket selection matching ACDLL5/CDLL5 production behavior.
 *
 * Forked from LCHistTrainerCDLL4 (Chloe, April 2026).
 *
 * Usage: java cardinality.LCHistTrainerACDLL5 buckets=2048 hbits=2 iters=1m out=table.tsv.gz
 *
 * @author UMP45
 * @date April 2026
 */
public class LCHistTrainerACDLL5{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods       ----------------*/
	/*--------------------------------------------------------------*/

	static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

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

	/** Compute compressed absolute tier from a hashed key. */
	static int tierOf(long key){
		final int rawNlz=Long.numberOfLeadingZeros(key);
		final int mBits;
		if(rawNlz>=43){
			mBits=(int)((key<<(rawNlz-42))&0xFFFFF);
		}else{
			mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
		}
		final int mant=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
		return (2*rawNlz+mant)/3;
	}

	/** Maps (tier, hist) to a flat grid index. Grid: (hbits+2) tier bins * (1<<hbits) hist slots. */
	static int stateIndex(int tier, int hist, int hbits){
		final int numHist=1<<hbits;
		return Math.min(tier, hbits+1)*numHist+(hist&(numHist-1));
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
			final int hmask=(1<<hbits)-1;

			for(int t=0; t<numTrials; t++){
				int[] tierArr=new int[B];
				int[] hist=new int[B];
				int[] distinct=new int[B];
				boolean[] filled=new boolean[B];
				int filledCount=0;

				while(filledCount<B){
					final long number=rng.nextLong();
					final long key=hash64shift(number);
					final int tier=tierOf(key);
					final int bucket=(int)(Long.remainderUnsigned(key, B));

					distinct[bucket]++;

					if(!filled[bucket]){
						tierArr[bucket]=tier;
						hist[bucket]=0;
						filled[bucket]=true;
						filledCount++;
						snapshot(filledCount, tierArr, hist, distinct, filled);
					}else{
						final int oldTier=tierArr[bucket];
						final int delta=tier-oldTier;
						if(delta>0){
							final int oldH=hist[bucket]&hmask;
							hist[bucket]=((oldH|(1<<hbits))>>delta)&hmask;
							tierArr[bucket]=tier;
						}else if(delta<0){
							final int diff=-delta;
							if(diff>=1 && diff<=hbits){
								hist[bucket]|=1<<(hbits-diff);
							}
						}
					}
				}
			}
		}

		void snapshot(int filledCount, int[] tier, int[] hist, int[] distinct, boolean[] filled){
			for(int b=0; b<buckets; b++){
				if(filled[b]){
					int si=stateIndex(tier[b], hist[b], hbits);
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
			if(eq<0) continue;
			final String a=arg.substring(0, eq).toLowerCase();
			final String b=arg.substring(eq+1);
			if(a.equals("buckets")||a.equals("b")){buckets=Parse.parseIntKMG(b);}
			else if(a.equals("hbits")||a.equals("historybits")){hbits=Integer.parseInt(b);}
			else if(a.equals("iters")||a.equals("trials")||a.equals("i")){iters=Parse.parseIntKMG(b);}
			else if(a.equals("threads")||a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("seed")){masterSeed=Long.parseLong(b);}
			else if(a.equals("out")){out=b;}
		}

		final int numThreads=Math.min(threads, (int)iters);
		System.err.println("LCHistTrainerACDLL5: buckets="+buckets+" hbits="+hbits+
			" iters="+iters+" threads="+numThreads+" seed="+masterSeed);

		final TrainerThread[] workers=new TrainerThread[numThreads];
		for(int t=0; t<numThreads; t++){
			final int threadTrials=(int)(iters/numThreads+(t<iters%numThreads?1:0));
			workers[t]=new TrainerThread(masterSeed+t, threadTrials, buckets, hbits);
		}
		for(TrainerThread w : workers){w.start();}
		for(TrainerThread w : workers){
			try{ w.join(); }catch(InterruptedException e){ Thread.currentThread().interrupt(); }
			if(!w.success){System.err.println("Warning: a trainer thread did not complete successfully.");}
		}

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

		final int numValid=3*(1<<hbits)-1;
		final int[][] validStates=new int[numValid][2];
		final String[] stateNames=new String[numValid];
		int vi=0;
		for(int bin=0; bin<=hbits+1; bin++){
			final int validSlots=Math.min(bin, hbits);
			final int invalidBits=hbits-validSlots;
			final int numHist=1<<validSlots;
			final String binLabel=(bin>hbits) ? (bin+"+") : String.valueOf(bin);
			for(int j=0; j<numHist; j++){
				final int hist=j<<invalidBits;
				final int gridIdx=bin*(1<<hbits)+hist;
				final int minDistinct=1+Integer.bitCount(j);
				validStates[vi][0]=gridIdx;
				validStates[vi][1]=minDistinct;
				StringBuilder name=new StringBuilder(binLabel).append('.');
				for(int bb=validSlots-1; bb>=0; bb--){ name.append((j>>>bb)&1); }
				if(validSlots==0) name.append('0');
				stateNames[vi]=name.toString();
				vi++;
			}
		}

		int zeros=0;
		for(int s=0; s<nStates; s++){ if(totalObs[buckets][s]==0) zeros++; }
		final int expectedZeros=nStates-numValid;
		System.err.println("Zeros in final row: "+zeros+" (expected "+expectedZeros+")");

		double[][] avg=new double[buckets+1][nStates];
		for(int f=1; f<=buckets; f++){
			for(int s=0; s<nStates; s++){
				if(totalObs[f][s]>0){ avg[f][s]=(double)totalSums[f][s]/totalObs[f][s]; }
			}
		}

		try{
			OutputStream os;
			if(out.equals("stdout.txt")){ os=System.out; }
			else if(out.endsWith(".gz")){ os=new GZIPOutputStream(new FileOutputStream(out)); }
			else{ os=new FileOutputStream(out); }
			PrintWriter pw=new PrintWriter(new OutputStreamWriter(os));

			pw.println("#LC_HISTORY_TABLE_ACDLL5");
			pw.println("#Buckets\t"+buckets);
			pw.println("#HistBits\t"+hbits);
			pw.println("#TierGeom\tACDLL5_half_nlz");
			pw.println("#Trials\t"+iters);

			StringBuilder hdr=new StringBuilder("#Filled");
			for(String name : stateNames) hdr.append('\t').append(name);
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
			if(os!=System.out){ pw.close(); }
		}catch(IOException e){ e.printStackTrace(); }

		System.err.println("Output: "+out);
	}
}
