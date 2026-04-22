package cardinality;

import parse.Parse;
import parse.PreParser;
import shared.Shared;

import java.io.*;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

/**
 * LC history CF table trainer for AUDLL32 (standard tiers, 3-state collapsed history).
 * Tier = rawNlz directly (TIER_SCALE=1.0, no mantissa compression).
 * History: 3 states collapsed from 2-bit model (10+11 → state 2):
 *   State 0 (00): no sub-tier observations
 *   State 1 (01): observed at -2 tiers only
 *   State 2 (1x): observed at -1 tier (possibly also -2)
 * Modulo bucket selection matching AUDLL32 production behavior.
 * No overflow modeling — 13 exponent states, overflow negligible.
 * No floor advancement — same reasoning.
 *
 * Forked from LCHistTrainerACDLL5 (UMP45, April 2026).
 *
 * Usage: java cardinality.LCHistTrainerAUDLL32 buckets=1536 hbits=2 iters=1m t=4
 *
 * @author Neptune
 * @date April 2026
 */
public class LCHistTrainerAUDLL32 {

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

	/** Standard tier: rawNlz directly (no compression). */
	static int tierOf(long key){
		return Long.numberOfLeadingZeros(key);
	}

	/*--------------------------------------------------------------*/
	/*----------------   3-State History Transitions ----------------*/
	/*--------------------------------------------------------------*/

	/** History state after tier advance by delta.
	 *  delta=1: any→3 (11, collapsed 1x), delta=2: any→1 (01), delta>=3: any→0 (00).
	 *  Note: state 2 (binary 10) is unused; state 3 (binary 11) is the collapsed state. */
	static int advanceHistory(int oldHist, int delta){
		if(delta==1){return 3;}
		if(delta==2){return 1;}
		return 0;
	}

	/** Set history for sub-tier observation at given distance below current tier.
	 *  diff=1 (-1 tier): always becomes state 3 (binary 11 = collapsed 1x).
	 *  diff=2 (-2 tiers): becomes state 1 unless already state 3 (preserve -1 info).
	 *  Note: internal state 2 (binary 10) is unused; state 3 (binary 11) represents
	 *  the collapsed "1x" state for HC bit-counting compatibility. */
	static int observeHistory(int oldHist, int diff){
		if(diff==1){return 3;}
		return Math.max(oldHist, 1);
	}

	/** Maps (tier, hist) to flat grid index.
	 *  Uses hbits=2 grid (4 hist slots per tier bin) for CardStats compatibility.
	 *  State 3 (11) is never produced, so those grid cells stay at zero. */
	static int stateIndex(int tier, int hist, int hbits){
		final int numHist=1<<hbits;
		return Math.min(tier, hbits+1)*numHist+(hist&(numHist-1));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker Thread        ----------------*/
	/*--------------------------------------------------------------*/

	static final class TrainerThread extends Thread {

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
							// Tier advance: use 3-state collapsed transition
							hist[bucket]=advanceHistory(hist[bucket], delta);
							tierArr[bucket]=tier;
						}else if(delta<0){
							// Sub-tier observation: set history
							final int diff=-delta;
							if(diff>=1 && diff<=2){
								hist[bucket]=observeHistory(hist[bucket], diff);
							}
						}
					}
				}
			}
		}

		void snapshot(int filledCount, int[] tier, int[] hist, int[] distinct, boolean[] filled){
			for(int b=0; b<buckets; b++){
				if(filled[b]){
					// Apply same emit logic as AUDLL32.emitHistMask:
					// tier 0: hist=0, tier 1: state 3→2 (only high bit), tier 2+: as-is
					final int emitHist;
					if(tier[b]==0){emitHist=0;}
					else if(tier[b]==1){emitHist=(hist[b]==3) ? 2 : 0;}
					else if(hist[b]==3){emitHist=3;} // tier 2+: collapsed as 11
					else{emitHist=hist[b];}
					int si=stateIndex(tier[b], emitHist, hbits);
					collisionSums[filledCount][si]+=distinct[b];
					observations[filledCount][si]++;
				}
			}
		}

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
		int buckets=1536;
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
		System.err.println("LCHistTrainerAUDLL32: buckets="+buckets+" hbits="+hbits+
			" iters="+iters+" threads="+numThreads+" seed="+masterSeed);

		final TrainerThread[] workers=new TrainerThread[numThreads];
		for(int t=0; t<numThreads; t++){
			final int threadTrials=(int)(iters/numThreads+(t<iters%numThreads?1:0));
			workers[t]=new TrainerThread(masterSeed+t, threadTrials, buckets, hbits);
		}
		for(TrainerThread w : workers) w.start();
		for(TrainerThread w : workers){
			try{ w.join(); }catch(InterruptedException e){ Thread.currentThread().interrupt(); }
			if(!w.success) System.err.println("Warning: a trainer thread did not complete successfully.");
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

			pw.println("#LC_HISTORY_TABLE_AUDLL32");
			pw.println("#Buckets\t"+buckets);
			pw.println("#HistBits\t"+hbits);
			pw.println("#TierGeom\tstandard_nlz");
			pw.println("#Trials\t"+iters);

			StringBuilder hdr=new StringBuilder("#Filled");
			for(String name : stateNames) hdr.append('\t').append(name);
			pw.println(hdr);

			for(int f=1; f<=buckets; f++){
				StringBuilder row=new StringBuilder();
				row.append(f);
				for(int v=0; v<numValid; v++){
					int si=validStates[v][0];
					double val=avg[f][si];
					row.append('\t').append(String.format("%.2f", val));
				}
				pw.println(row);
			}
			pw.flush();
			if(os!=System.out){ os.close(); }
		}catch(IOException e){ throw new RuntimeException(e); }

		System.err.println("SBS table written to "+(out.equals("stdout.txt")?"stdout":out));
	}
}
