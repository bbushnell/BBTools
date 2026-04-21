package cardinality;

import parse.Parse;
import parse.PreParser;
import shared.Shared;

import java.io.*;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

/**
 * LC history CF table trainer for BCDLL5 (banked compressed 2-bit history).
 * Models full BCDLL5 behavior: compressed tiers (halfNlz/3), 2-bit carry-shift
 * history, per-word 2-bit bank exponents, bank promotion, global floor
 * advancement with bank absorption, and HISTORY_MARGIN=2 eemask relaxation.
 *
 * Usage: java cardinality.LCHistTrainerBCDLL5 buckets=192 iters=1m out=table.tsv.gz
 *
 * @author Eru
 * @date April 2026
 */
public class LCHistTrainerBCDLL5 {

	/*--------------------------------------------------------------*/
	/*----------------            Main              ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		{PreParser pp=new PreParser(args, null, false); args=pp.args;}
		final long t0=System.nanoTime();
		int buckets=192;
		long iters=100000;
		int threads=Shared.threads();
		long masterSeed=42;
		String out="stdout.txt";

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
		}

		// Power of 2 for bucketMask (matching real BCDLL5 behavior)
		buckets=Integer.highestOneBit(buckets);
		if(buckets<2){buckets=2;}

		final int numThreads=Math.min(threads, (int)iters);
		System.err.println("LCHistTrainerBCDLL5: buckets="+buckets+" hbits="+HBITS+
			" iters="+iters+" threads="+numThreads+" seed="+masterSeed);

		final TrainerThread[] workers=new TrainerThread[numThreads];
		for(int i=0; i<numThreads; i++){
			final int threadTrials=(int)(iters/numThreads+(i<iters%numThreads?1:0));
			workers[i]=new TrainerThread(masterSeed+i, threadTrials, buckets);
		}
		for(TrainerThread w : workers){w.start();}
		for(TrainerThread w : workers){
			try{w.join();}catch(InterruptedException e){Thread.currentThread().interrupt();}
			if(!w.success){System.err.println("Warning: a trainer thread did not complete successfully.");}
		}

		final int nStates=(HBITS+2)*(1<<HBITS);
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

		// Build valid state list (same as LCHistTrainerCDLL4)
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
				for(int bb=validSlots-1; bb>=0; bb--){name.append((j>>>bb)&1);}
				if(validSlots==0){name.append('0');}
				stateNames[vi]=name.toString();
				vi++;
			}
		}

		int zeros=0;
		for(int s=0; s<nStates; s++){if(totalObs[buckets][s]==0){zeros++;}}
		final int expectedZeros=nStates-numValid;
		System.err.println("Zeros in final row: "+zeros+" (expected "+expectedZeros+")");

		double[][] avg=new double[buckets+1][nStates];
		for(int f=1; f<=buckets; f++){
			for(int s=0; s<nStates; s++){
				if(totalObs[f][s]>0){avg[f][s]=(double)totalSums[f][s]/totalObs[f][s];}
			}
		}

		try{
			OutputStream os;
			if(out.equals("stdout.txt")){os=System.out;}
			else if(out.endsWith(".gz")){os=new GZIPOutputStream(new FileOutputStream(out));}
			else{os=new FileOutputStream(out);}
			PrintWriter pw=new PrintWriter(new OutputStreamWriter(os));

			pw.println("#LC_HISTORY_TABLE_BCDLL5");
			pw.println("#Buckets\t"+buckets);
			pw.println("#HistBits\t"+HBITS);
			pw.println("#TierGeom\tBCDLL5_banked_half_nlz");
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
			System.err.println("Output: "+out);
		}catch(IOException e){
			e.printStackTrace();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker Thread        ----------------*/
	/*--------------------------------------------------------------*/

	static final class TrainerThread extends Thread {

		/** Construct a trainer thread for one batch of trials. */
		TrainerThread(final long seed, final int numTrials, final int buckets){
			this.seed=seed;
			this.numTrials=numTrials;
			this.buckets=buckets;
			this.words=(buckets+REGS_PER_WORD-1)/REGS_PER_WORD;
			this.bucketMask=buckets-1;
			this.nStates=(HBITS+2)*(1<<HBITS);
			collisionSums=new long[buckets+1][nStates];
			observations=new long[buckets+1][nStates];
		}

		@Override public void run(){
			try{runInner(); success=true;}
			catch(Exception e){e.printStackTrace();}
		}

		/** Run all trials, accumulating collision counts into collisionSums/observations. */
		void runInner(){
			final Random rng=new Random(seed);
			final int B=buckets;

			for(int t=0; t<numTrials; t++){
				// Per-bucket state
				int[] tierPart=new int[B];  // 0=empty/floor-level, 1-7=relTier+1
				int[] hist=new int[B];
				int[] distinct=new int[B];

				// Per-word bank exponents
				int[] bank=new int[words];

				// Global state
				int globalNLZ=-1;
				int filledCount=0;  // buckets with tierPart>0 at bank=0 level
				int minZeroCount=B;
				long eeMask=-1L;

				while(filledCount<B){
					final long number=rng.nextLong();
					final long key=hash64shift(number);

					if(Long.compareUnsigned(key, eeMask)>0){continue;}

					final int absTier=tierOf(key);
					final int bucket=(int)(key&bucketMask);
					final int wordIdx=bucket/REGS_PER_WORD;
					final int bk=bank[wordIdx];
					final int relTier=absTier-globalNLZ-1-bk;

					if(relTier<-HISTORY_MARGIN){continue;}

					distinct[bucket]++;

					final int oldTp=tierPart[bucket];
					final int oldHist=hist[bucket];

					if(relTier>=0){
						int newTp=Math.min(relTier+1, 7);

						// Bank promotion on overflow
						if(newTp>7 && bk<3){
							boolean canPromote=true;
							final int base=wordIdx*REGS_PER_WORD;
							for(int r=0; r<REGS_PER_WORD; r++){
								if(tierPart[base+r]==0){canPromote=false; break;}
							}
							if(canPromote){
								bank[wordIdx]++;
								for(int r=0; r<REGS_PER_WORD; r++){
									tierPart[base+r]--;
								}
								final int newBk=bank[wordIdx];
								final int newRelTier=absTier-globalNLZ-1-newBk;
								newTp=Math.min(newRelTier+1, 7);
								// Re-read bucket state after promotion
								final int promTp=tierPart[bucket];
								final int promHist=hist[bucket];
								if(newTp<=promTp){continue;}
								final int delta=(promTp>0) ? (newTp-promTp) : (newRelTier+1);
								final int carry=(promTp>0 || globalNLZ+newBk>=0) ? HIST_CARRY : 0;
								hist[bucket]=((promHist|carry)>>delta)&HMASK;
								tierPart[bucket]=newTp;
								continue;
							}
						}

						if(newTp>oldTp){
							final int delta=(oldTp>0) ? (newTp-oldTp) : (relTier+1);
							final int carry=(oldTp>0 || globalNLZ+bk>=0) ? HIST_CARRY : 0;
							hist[bucket]=((oldHist|carry)>>delta)&HMASK;
							tierPart[bucket]=newTp;

							if(oldTp==0 && bk==0){
								filledCount++;
								snapshot(filledCount, tierPart, hist, bank, distinct,
									globalNLZ, B);
								if(--minZeroCount<1){
									while(minZeroCount==0 && globalNLZ<64){
										globalNLZ++;
										final int relaxedTier=Math.max(0, globalNLZ-HISTORY_MARGIN);
										final int minNlz=(3*relaxedTier)/2;
										eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
										minZeroCount=countAndDecrement(
											tierPart, hist, bank, words, B, globalNLZ);
									}
									// Recount filled after floor advance
									filledCount=0;
									for(int i=0; i<B; i++){
										if(tierPart[i]>0 || bank[i/REGS_PER_WORD]>0 || globalNLZ>=0){
											filledCount++;
										}
									}
								}
							}
						}else if(newTp<oldTp){
							final int diff=oldTp-newTp;
							if(diff>=1 && diff<=HBITS){hist[bucket]|=1<<(HBITS-diff);}
						}
					}else if(relTier==-1){
						if(oldTp>0){
							final int diff=oldTp;
							if(diff>=1 && diff<=HBITS){hist[bucket]|=1<<(HBITS-diff);}
						}
					}else{
						// relTier == -2
						final int diff=(oldTp>0) ? (oldTp+1) : 1;
						if(diff>=1 && diff<=HBITS){hist[bucket]|=1<<(HBITS-diff);}
					}
				}
			}
		}

		/**
		 * Decrement all registers after a global floor advance.
		 * Words with bank > 0 absorb by decrementing bank only.
		 * Words at bank=0 decrement each tierPart value.
		 * @return New minZeroCount (buckets at tierPart=0 eligible for next advance)
		 */
		static int countAndDecrement(final int[] tierPart, final int[] hist, final int[] bank,
				final int words, final int buckets, final int globalNLZ){
			int newMinZeroCount=0;
			for(int w=0; w<words; w++){
				final int bk=bank[w];
				if(bk>0){
					bank[w]=bk-1;
					if(bk==1){
						final int base=w*REGS_PER_WORD;
						for(int r=0; r<REGS_PER_WORD; r++){
							if(base+r<buckets && tierPart[base+r]==0){newMinZeroCount++;}
						}
					}
				}else{
					final int base=w*REGS_PER_WORD;
					for(int r=0; r<REGS_PER_WORD; r++){
						final int idx=base+r;
						if(idx>=buckets){break;}
						if(tierPart[idx]>0){
							tierPart[idx]--;
							if(tierPart[idx]==0){newMinZeroCount++;}
						}else{
							newMinZeroCount++;
						}
					}
				}
			}
			return newMinZeroCount;
		}

		/** Record collision and observation counts for all filled buckets at this fill level. */
		void snapshot(final int filledCount, final int[] tierPart, final int[] hist, final int[] bank,
				final int[] distinct, final int globalNLZ, final int B){
			for(int b=0; b<B; b++){
				final int tp=tierPart[b];
				final int bk=bank[b/REGS_PER_WORD];
				final int h=hist[b];
				final int absTier=tp+globalNLZ+bk;
				if(absTier<0){continue;} // truly empty
				final int si=stateIndex(absTier, h, HBITS);
				collisionSums[filledCount][si]+=distinct[b];
				observations[filledCount][si]++;
			}
		}

		final long seed;
		final int numTrials;
		final int buckets;
		final int words;
		final int bucketMask;
		final int nStates;
		final long[][] collisionSums;
		final long[][] observations;
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Compute compressed-tier state index from absolute tier and history bits. */
	static int stateIndex(final int tier, final int hist, final int hbits){
		final int numHist=1<<hbits;
		return Math.min(tier, hbits+1)*numHist+(hist&(numHist-1));
	}

	/** Compute compressed absolute tier from a 64-bit hash key. */
	static int tierOf(final long key){
		final int rawNlz=Long.numberOfLeadingZeros(key);
		final int mBits;
		if(rawNlz>=43){
			mBits=(int)((key<<(rawNlz-42))&0xFFFFF); // zero-extend remaining bits
		}else{
			mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
		}
		final int mant=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
		return (2*rawNlz+mant)/3;
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

	/** Mantissa threshold for compressed tiers: (2-sqrt(2)) * 1048576. */
	static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);
	/** Number of history bits per bucket. */
	static final int HBITS=2;
	/** Bitmask for history: (1<<HBITS)-1 = 0x3. */
	static final int HMASK=(1<<HBITS)-1;
	/** Carry bit injected on tier advance: 1<<HBITS = 4. */
	static final int HIST_CARRY=1<<HBITS;
	/** eeMask relaxation: accept hashes this many tiers below floor for history. */
	static final int HISTORY_MARGIN=2;
	/** Registers (buckets) packed per word. */
	static final int REGS_PER_WORD=6;

}
