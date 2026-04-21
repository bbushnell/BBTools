package cardinality;

import parse.Parse;
import parse.PreParser;
import shared.Shared;

import java.io.*;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

/**
 * HSB (per-tier state-bias) table generator for BCDLL5.
 * Multi-bucket simulator that models full banked compressed tier behavior:
 * compressed tiers (halfNlz/3), 2-bit carry-shift history, per-word 2-bit
 * bank exponents, bank promotion, global floor advancement with bank
 * absorption, and HISTORY_MARGIN=2.
 *
 * For each (tier, hist-state), records the average true cardinality of buckets
 * occupying that state, then computes log2-space correction factors:
 *   CF(state) = log2(mean_card_in_state / mean_card_in_tier)
 *
 * Output format matches mantissacompare.sh mode=ctll for StateTable.loadHsbTable.
 *
 * Usage: java cardinality.MantissaCompareBCDLL5 buckets=256 outer=131072 inner=32768 maxtier=11 out=hsb.tsv.gz
 *
 * @author Eru
 * @date April 2026
 */
public class MantissaCompareBCDLL5 {

	static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);
	static final int HBITS=2;
	static final int HMASK=(1<<HBITS)-1;
	static final int HIST_CARRY=1<<HBITS;
	static final int HISTORY_MARGIN=2;
	static final int REGS_PER_WORD=6;
	static final int NUM_HIST=1<<HBITS; // 4

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

	static int tierOf(long key){
		final int rawNlz=Long.numberOfLeadingZeros(key);
		final int mant;
		final int mBits;
		if(rawNlz>=43){
			mBits=(int)((key<<(rawNlz-42))&0xFFFFF); // zero-extend remaining bits
		}else{
			mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
		}
		mant=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
		return (2*rawNlz+mant)/3;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker Thread        ----------------*/
	/*--------------------------------------------------------------*/

	static final class SimThread extends Thread {

		SimThread(long seed, int numTrials, int buckets, int maxTier, int innerSteps){
			this.seed=seed;
			this.numTrials=numTrials;
			this.buckets=buckets;
			this.words=(buckets+REGS_PER_WORD-1)/REGS_PER_WORD;
			this.bucketMask=buckets-1;
			this.maxTier=maxTier;
			this.innerSteps=innerSteps;
			// Accumulators: [tier][hist] -> sum of cardinalities and count
			cardSums=new double[maxTier+1][NUM_HIST];
			counts=new long[maxTier+1][NUM_HIST];
		}

		@Override
		public void run(){
			try{ runInner(); success=true; }
			catch(Exception e){ e.printStackTrace(); }
		}

		void runInner(){
			final Random rng=new Random(seed);
			final int B=buckets;

			for(int trial=0; trial<numTrials; trial++){
				int[] tierPart=new int[B];
				int[] hist=new int[B];
				int[] bank=new int[words];
				int globalNLZ=-1;
				int minZeroCount=B;
				long eeMask=-1L;
				long cardinality=0;

				for(int step=0; step<innerSteps; step++){
					final long number=rng.nextLong();
					final long key=hash64shift(number);

					if(Long.compareUnsigned(key, eeMask)>0){continue;}

					final int absTier=tierOf(key);
					final int bucket=(int)(key&bucketMask);
					final int wordIdx=bucket/REGS_PER_WORD;
					final int bk=bank[wordIdx];
					final int relTier=absTier-globalNLZ-1-bk;

					if(relTier<-HISTORY_MARGIN){continue;}

					cardinality++;

					final int oldTp=tierPart[bucket];
					final int oldHist=hist[bucket];

					if(relTier>=0){
						int newTp=Math.min(relTier+1, 7);

						// Bank promotion on overflow
						if(newTp>7 && bk<3){
							boolean canPromote=true;
							final int base=wordIdx*REGS_PER_WORD;
							for(int r=0; r<REGS_PER_WORD && base+r<B; r++){
								if(tierPart[base+r]==0){canPromote=false; break;}
							}
							if(canPromote){
								bank[wordIdx]++;
								for(int r=0; r<REGS_PER_WORD && base+r<B; r++){
									tierPart[base+r]--;
								}
								final int newBk=bank[wordIdx];
								final int newRelTier=absTier-globalNLZ-1-newBk;
								newTp=Math.min(newRelTier+1, 7);
								final int promTp=tierPart[bucket];
								final int promHist=hist[bucket];
								if(newTp<=promTp){
									sampleAll(tierPart, hist, bank, globalNLZ, B, cardinality);
									continue;
								}
								final int delta=(promTp>0) ? (newTp-promTp) : (newRelTier+1);
								final int carry=(promTp>0 || globalNLZ+newBk>=0) ? HIST_CARRY : 0;
								hist[bucket]=((promHist|carry)>>delta)&HMASK;
								tierPart[bucket]=newTp;
								sampleAll(tierPart, hist, bank, globalNLZ, B, cardinality);
								continue;
							}
						}

						if(newTp>oldTp){
							final int delta=(oldTp>0) ? (newTp-oldTp) : (relTier+1);
							final int carry=(oldTp>0 || globalNLZ+bk>=0) ? HIST_CARRY : 0;
							hist[bucket]=((oldHist|carry)>>delta)&HMASK;
							tierPart[bucket]=newTp;

							if(oldTp==0 && bk==0){
								if(--minZeroCount<1){
									while(minZeroCount==0 && globalNLZ<63){
										globalNLZ++;
										final int relaxedTier=Math.max(0, globalNLZ+1-HISTORY_MARGIN);
										final int minNlz=(3*relaxedTier)/2;
										eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
										minZeroCount=countAndDecrement(
											tierPart, hist, bank, words, B);
									}
								}
							}
						}else if(newTp<oldTp){
							final int diff=oldTp-newTp;
							if(diff>=1 && diff<=HBITS){
								hist[bucket]|=1<<(HBITS-diff);
							}
						}
					}else if(relTier==-1){
						if(oldTp>0){
							final int diff=oldTp;
							if(diff>=1 && diff<=HBITS){
								hist[bucket]|=1<<(HBITS-diff);
							}
						}
					}else{
						final int diff=(oldTp>0) ? (oldTp+1) : 1;
						if(diff>=1 && diff<=HBITS){
							hist[bucket]|=1<<(HBITS-diff);
						}
					}

					sampleAll(tierPart, hist, bank, globalNLZ, B, cardinality);
				}
			}
		}

		void sampleAll(int[] tierPart, int[] hist, int[] bank, int globalNLZ,
				int B, long cardinality){
			for(int i=0; i<B; i++){
				final int tp=tierPart[i];
				final int bk=bank[i/REGS_PER_WORD];
				final int h=hist[i];
				final int absTier=tp+globalNLZ+bk;
				if(absTier<0){continue;}
				if(absTier>=0 && absTier<=maxTier){
					cardSums[absTier][h]+=cardinality;
					counts[absTier][h]++;
				}
			}
		}

		static int countAndDecrement(int[] tierPart, int[] hist, int[] bank,
				int words, int buckets){
			int newMinZeroCount=0;
			for(int w=0; w<words; w++){
				final int bk=bank[w];
				if(bk>0){
					bank[w]=bk-1;
					if(bk==1){
						final int base=w*REGS_PER_WORD;
						for(int r=0; r<REGS_PER_WORD; r++){
							if(base+r<buckets && tierPart[base+r]==0){
								newMinZeroCount++;
							}
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

		final long seed;
		final int numTrials;
		final int buckets;
		final int words;
		final int bucketMask;
		final int maxTier;
		final int innerSteps;
		final double[][] cardSums;
		final long[][] counts;
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Main              ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		{ PreParser pp=new PreParser(args, null, false); args=pp.args; }
		final long t0=System.nanoTime();
		int buckets=256;
		int outerTrials=131072;
		int innerSteps=32768;
		int maxTier=11;
		int threads=Shared.threads();
		long masterSeed=42;
		String out="stdout.txt";
		boolean linCF=true;

		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<0) continue;
			final String a=arg.substring(0, eq).toLowerCase();
			final String b=arg.substring(eq+1);
			if(a.equals("buckets")||a.equals("b")){buckets=Parse.parseIntKMG(b);}
			else if(a.equals("outer")||a.equals("trials")){outerTrials=Parse.parseIntKMG(b);}
			else if(a.equals("inner")||a.equals("steps")){innerSteps=Parse.parseIntKMG(b);}
			else if(a.equals("maxtier")){maxTier=Integer.parseInt(b);}
			else if(a.equals("threads")||a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("seed")){masterSeed=Long.parseLong(b);}
			else if(a.equals("out")){out=b;}
			else if(a.equals("avg")){linCF=b.equalsIgnoreCase("lin");}
		}

		buckets=Integer.highestOneBit(buckets);
		if(buckets<2){buckets=2;}
		final int numThreads=Math.min(threads, outerTrials);
		System.err.println("MantissaCompareBCDLL5: buckets="+buckets+" outer="+outerTrials+
			" inner="+innerSteps+" maxTier="+maxTier+" threads="+numThreads+
			" avg="+(linCF?"lin":"geo")+" seed="+masterSeed);

		final SimThread[] workers=new SimThread[numThreads];
		for(int i=0; i<numThreads; i++){
			final int threadTrials=outerTrials/numThreads+(i<outerTrials%numThreads?1:0);
			workers[i]=new SimThread(masterSeed+i, threadTrials, buckets, maxTier, innerSteps);
		}
		for(SimThread w : workers) w.start();
		for(SimThread w : workers){
			try{ w.join(); }catch(InterruptedException e){ Thread.currentThread().interrupt(); }
			if(!w.success) System.err.println("Warning: a sim thread did not complete successfully.");
		}

		// Merge accumulators
		final double[][] totalCardSums=new double[maxTier+1][NUM_HIST];
		final long[][] totalCounts=new long[maxTier+1][NUM_HIST];
		for(SimThread w : workers){
			for(int t=0; t<=maxTier; t++){
				for(int h=0; h<NUM_HIST; h++){
					totalCardSums[t][h]+=w.cardSums[t][h];
					totalCounts[t][h]+=w.counts[t][h];
				}
			}
		}

		final double elapsed=(System.nanoTime()-t0)/1e9;
		System.err.println("Simulation done: "+String.format("%.1f", elapsed)+"s");

		// Compute per-state mean cardinalities
		final double[][] meanCard=new double[maxTier+1][NUM_HIST];
		for(int t=0; t<=maxTier; t++){
			for(int h=0; h<NUM_HIST; h++){
				if(totalCounts[t][h]>0){
					meanCard[t][h]=totalCardSums[t][h]/totalCounts[t][h];
				}
			}
		}

		// Compute per-tier average and CFs
		final double[] tierMean=new double[maxTier+1];
		final double[][] cf=new double[maxTier+1][NUM_HIST];
		for(int t=0; t<=maxTier; t++){
			if(linCF){
				// Linear average: weighted by observation count
				double sumCard=0; long sumCount=0;
				for(int h=0; h<NUM_HIST; h++){
					sumCard+=totalCardSums[t][h];
					sumCount+=totalCounts[t][h];
				}
				tierMean[t]=(sumCount>0) ? sumCard/sumCount : 0;
			}else{
				// Geometric average of state means
				double logSum=0; int n=0;
				for(int h=0; h<NUM_HIST; h++){
					if(meanCard[t][h]>0){logSum+=Math.log(meanCard[t][h]); n++;}
				}
				tierMean[t]=(n>0) ? Math.exp(logSum/n) : 0;
			}
			for(int h=0; h<NUM_HIST; h++){
				if(meanCard[t][h]>0 && tierMean[t]>0){
					cf[t][h]=Math.log(meanCard[t][h]/tierMean[t])/Math.log(2.0);
				}
			}
		}

		// Compute steady-state CF (average of high tiers)
		final double[] ssCF=new double[NUM_HIST];
		int ssTiers=0;
		for(int t=Math.max(8, maxTier-3); t<=maxTier; t++){
			boolean allNonzero=true;
			for(int h=0; h<NUM_HIST; h++){if(totalCounts[t][h]==0){allNonzero=false; break;}}
			if(allNonzero){
				for(int h=0; h<NUM_HIST; h++){ssCF[h]+=cf[t][h];}
				ssTiers++;
			}
		}
		if(ssTiers>0){for(int h=0; h<NUM_HIST; h++){ssCF[h]/=ssTiers;}}

		// Output in HSB format (matches StateTable.loadHsbTable)
		try{
			OutputStream os;
			if(out.equals("stdout.txt")){ os=System.out; }
			else if(out.endsWith(".gz")){ os=new GZIPOutputStream(new FileOutputStream(out)); }
			else{ os=new FileOutputStream(out); }
			PrintWriter pw=new PrintWriter(new OutputStreamWriter(os));

			pw.println("# BCDLL5 state-bias table, 2-bit history, banked CTLL compressed-tier geometry");
			pw.println("# Generated via MantissaCompareBCDLL5 buckets="+buckets+
				" outer="+outerTrials+" inner="+innerSteps+" maxtier="+maxTier+
				" avg="+(linCF?"lin":"geo"));
			pw.println("# LinCF values (avg=lin)");
			pw.println("# Columns: tier  CF(hist=00)  CF(hist=01)  CF(hist=10)  CF(hist=11)");

			for(int t=0; t<=maxTier; t++){
				StringBuilder row=new StringBuilder();
				row.append(t);
				for(int h=0; h<NUM_HIST; h++){
					row.append('\t').append(String.format("%+.8f", cf[t][h]));
				}
				pw.println(row);
			}
			// Steady-state row
			StringBuilder ssRow=new StringBuilder("ss");
			for(int h=0; h<NUM_HIST; h++){
				ssRow.append('\t').append(String.format("%+.8f", ssCF[h]));
			}
			pw.println(ssRow);

			pw.flush();
			if(os!=System.out){pw.close();}
			System.err.println("Output: "+out);
		}catch(IOException e){
			e.printStackTrace();
		}
	}
}
