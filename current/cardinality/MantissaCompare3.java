package cardinality;

import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;

import java.io.*;
import java.util.zip.GZIPOutputStream;

/**
 * Multi-bucket HSB table generator with per-bucket cardinality tracking.
 * Replaces MantissaCompareBCDLL5 (which used global cardinality — broken).
 *
 * Three modes:
 *   history — uncompressed 2x tiers, validates against MC2 mode=history
 *   ctll    — compressed tiers (halfNlz/3), validates against MC2 mode=ctll
 *   bcdll5  — compressed tiers + 2-bit banking per 6-register word
 *
 * Per-bucket cardinality: each bucket tracks trueCard[i]++. When sampling,
 * sum[tier][state] += trueCard[i], treating each bucket as an independent
 * single-register experiment.
 *
 * @author Eru
 * @date April 2026
 */
public class MantissaCompare3 {

	static final int MODE_HISTORY=0, MODE_CTLL=1, MODE_BCDLL5=2;
	static final String[] MODE_NAMES={"history", "ctll", "bcdll5"};
	static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);
	static final int HBITS=2;
	static final int HMASK=(1<<HBITS)-1;
	static final int HIST_CARRY=1<<HBITS;
	static final int HISTORY_MARGIN=2;
	static final int REGS_PER_WORD=6;
	static final int NUM_HIST=1<<HBITS;

	static int tierOf(long key, int mode){
		final int rawNlz=Long.numberOfLeadingZeros(key);
		if(mode==MODE_HISTORY){return rawNlz;}
		final int mBits;
		if(rawNlz>=43){
			mBits=(int)((key<<(rawNlz-42))&0xFFFFF);
		}else{
			mBits=(int)((key>>>(42-rawNlz))&0xFFFFF);
		}
		final int mant=(mBits>=MANTISSA_THRESHOLD) ? 1 : 0;
		return (2*rawNlz+mant)/3;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker Thread        ----------------*/
	/*--------------------------------------------------------------*/

	static final class SimThread extends Thread {

		SimThread(long seed, int numTrials, int buckets, int maxTier,
				int innerSteps, int mode, boolean geoAvg){
			this.seed=seed;
			this.numTrials=numTrials;
			this.buckets=buckets;
			this.maxTier=maxTier;
			this.innerSteps=innerSteps;
			this.mode=mode;
			this.geoAvg=geoAvg;
			if(mode==MODE_BCDLL5){
				words=(buckets+REGS_PER_WORD-1)/REGS_PER_WORD;
			}else{
				words=0;
			}
			cardSums=new double[maxTier+1][NUM_HIST];
			geoSums=new double[maxTier+1][NUM_HIST];
			counts=new long[maxTier+1][NUM_HIST];
		}

		@Override
		public void run(){
			try{ runInner(); success=true; }
			catch(Exception e){ e.printStackTrace(); }
		}

		void runInner(){
			final rand.FastRandomXoshiro rng=new rand.FastRandomXoshiro(seed);
			final int B=buckets;
			final boolean isBanked=(mode==MODE_BCDLL5);
			final boolean isHistory=(mode==MODE_HISTORY);
			long bankPromotions=0;

			for(int trial=0; trial<numTrials; trial++){
				int[] tierPart=new int[B];
				int[] hist=new int[B];
				long[] trueCard=new long[B];
				int[] bank=isBanked ? new int[words] : null;
				int minZeros=0;
				int minZeroCount=B;
				long eeMask=-1L;

				for(int step=0; step<innerSteps; step++){
					final long number=rng.nextLong();
					final long key=Tools.hash64shift(number);
					if(Long.compareUnsigned(key, eeMask)>0){continue;}

					final int absTier=tierOf(key, mode);
					final int bucket=(int)(Long.remainderUnsigned(key, B));
					trueCard[bucket]++;

					final int bk=isBanked ? bank[bucket/REGS_PER_WORD] : 0;
					final int relTier=absTier-minZeros-bk;
					if(relTier<-HISTORY_MARGIN){
						sampleBucket(bucket, tierPart, hist, bank, minZeros, trueCard);
						continue;
					}

					final int oldTp=tierPart[bucket];
					final int oldHist=hist[bucket];

					if(relTier>=0){
						// Bank promotion (bcdll5 only)
						if(isBanked && relTier+1>7 && bk<3){
							final int wordIdx=bucket/REGS_PER_WORD;
							boolean canPromote=true;
							final int base=wordIdx*REGS_PER_WORD;
							for(int r=0; r<REGS_PER_WORD && base+r<B; r++){
								if(tierPart[base+r]==0){canPromote=false; break;}
							}
							if(canPromote){
								bank[wordIdx]++;
								bankPromotions++;
								for(int r=0; r<REGS_PER_WORD && base+r<B; r++){
									tierPart[base+r]--;
								}
								final int newBk=bank[wordIdx];
								final int newRelTier=absTier-minZeros-newBk;
								final int newTp=Math.min(newRelTier+1, 7);
								final int promTp=tierPart[bucket];
								final int promHist=hist[bucket];
								if(newTp<=promTp){
									sampleBucket(bucket, tierPart, hist, bank, minZeros, trueCard);
									continue;
								}
								final int delta=(promTp>0) ? (newTp-promTp) : (newRelTier+1);
								final int carry=(promTp>0 || minZeros+newBk>0) ? HIST_CARRY : 0;
								hist[bucket]=((promHist|carry)>>delta)&HMASK;
								tierPart[bucket]=newTp;
								sampleBucket(bucket, tierPart, hist, bank, minZeros, trueCard);
								continue;
							}
						}

						final int newTp=Math.min(relTier+1, 7);
						if(newTp>oldTp){
							final int delta=(oldTp>0) ? (newTp-oldTp) : (relTier+1);
							final int carry=(oldTp>0 || minZeros+bk>0) ? HIST_CARRY : 0;
							hist[bucket]=((oldHist|carry)>>delta)&HMASK;
							tierPart[bucket]=newTp;

							final boolean shouldDecrement=isBanked
								? (oldTp==0 && bk==0)
								: (oldTp==0);
							if(shouldDecrement && --minZeroCount<1){
								while(minZeroCount==0 && minZeros<64){
									if(isHistory){
										minZeros++;
										eeMask>>>=1;
									}else{
										minZeros++;
										final int relaxedTier=Math.max(0, minZeros-HISTORY_MARGIN);
										final int minNlz=(3*relaxedTier)/2;
										eeMask=(minNlz>=64) ? 0 : ~0L>>>minNlz;
									}
									minZeroCount=countAndDecrement(
										tierPart, hist, bank, words, B, isBanked);
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

					sampleBucket(bucket, tierPart, hist, bank, minZeros, trueCard);
				}
			}

			this.bankPromotions=bankPromotions;
		}

		void sampleBucket(int i, int[] tierPart, int[] hist, int[] bank,
				int minZeros, long[] trueCard){
			final int tp=tierPart[i];
			final int bk=(bank!=null) ? bank[i/REGS_PER_WORD] : 0;
			final int h=hist[i];
			final int absTier;
			if(tp>0){
				absTier=(tp-1)+minZeros+bk;
			}else if(minZeros+bk>0){
				absTier=minZeros+bk-1;
			}else{
				return;
			}
			if(absTier>=0 && absTier<=maxTier && trueCard[i]>0){
				cardSums[absTier][h]+=trueCard[i];
				if(geoAvg){geoSums[absTier][h]+=Math.log(trueCard[i]);}
				counts[absTier][h]++;
			}
		}

		static int countAndDecrement(int[] tierPart, int[] hist, int[] bank,
				int words, int buckets, boolean isBanked){
			int newMinZeroCount=0;
			if(isBanked){
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
			}else{
				for(int i=0; i<buckets; i++){
					if(tierPart[i]>0){
						tierPart[i]--;
						if(tierPart[i]==0){newMinZeroCount++;}
					}else{
						newMinZeroCount++;
					}
				}
			}
			return newMinZeroCount;
		}

		final long seed;
		final int numTrials, buckets, words, maxTier, innerSteps, mode;
		final boolean geoAvg;
		final double[][] cardSums, geoSums;
		final long[][] counts;
		long bankPromotions;
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
		int mode=MODE_CTLL;
		boolean geoAvg=false;

		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<0) continue;
			final String a=arg.substring(0, eq).toLowerCase();
			final String b=arg.substring(eq+1);
			if(a.equals("buckets")||a.equals("b")){buckets=Parse.parseIntKMG(b);}
			else if(a.equals("outer")||a.equals("trials")){outerTrials=Parse.parseIntKMG(b);}
			else if(a.equals("inner")||a.equals("steps")){innerSteps=Parse.parseIntKMG(b);}
			else if(a.equals("maxtier")||a.equals("mt")){maxTier=Integer.parseInt(b);}
			else if(a.equals("threads")||a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("seed")){masterSeed=Long.parseLong(b);}
			else if(a.equals("avg")){geoAvg=b.equalsIgnoreCase("geo");}
			else if(a.equals("mode")){
				if(b.equals("history")){mode=MODE_HISTORY;}
				else if(b.equals("ctll")){mode=MODE_CTLL;}
				else if(b.equals("bcdll5")){mode=MODE_BCDLL5;}
				else{throw new RuntimeException("Unknown mode '"+b+"'");}
			}
		}

		if(mode==MODE_BCDLL5){
			buckets=((buckets+REGS_PER_WORD-1)/REGS_PER_WORD)*REGS_PER_WORD;
			if(buckets<REGS_PER_WORD*32){buckets=REGS_PER_WORD*32;}
		}
		if(mode!=MODE_HISTORY && maxTier>11){maxTier=11;}

		final int numThreads=Math.min(threads, outerTrials);
		System.err.println("MantissaCompare3: mode="+MODE_NAMES[mode]+" buckets="+buckets+
			" outer="+outerTrials+" inner="+innerSteps+" maxTier="+maxTier+
			" threads="+numThreads+" avg="+(geoAvg?"geo":"lin")+" seed="+masterSeed);

		final SimThread[] workers=new SimThread[numThreads];
		for(int i=0; i<numThreads; i++){
			final int threadTrials=outerTrials/numThreads+(i<outerTrials%numThreads?1:0);
			workers[i]=new SimThread(masterSeed+i, threadTrials, buckets, maxTier,
				innerSteps, mode, geoAvg);
		}
		for(SimThread w : workers) w.start();
		long totalBankPromotions=0;
		for(SimThread w : workers){
			try{ w.join(); }catch(InterruptedException e){ Thread.currentThread().interrupt(); }
			if(!w.success) System.err.println("Warning: a sim thread did not complete.");
			totalBankPromotions+=w.bankPromotions;
		}

		// Merge
		final double[][] totalCardSums=new double[maxTier+1][NUM_HIST];
		final double[][] totalGeoSums=new double[maxTier+1][NUM_HIST];
		final long[][] totalCounts=new long[maxTier+1][NUM_HIST];
		for(SimThread w : workers){
			for(int t=0; t<=maxTier; t++){
				for(int h=0; h<NUM_HIST; h++){
					totalCardSums[t][h]+=w.cardSums[t][h];
					totalGeoSums[t][h]+=w.geoSums[t][h];
					totalCounts[t][h]+=w.counts[t][h];
				}
			}
		}

		final double elapsed=(System.nanoTime()-t0)/1e9;
		System.err.println("Simulation done: "+String.format("%.1f", elapsed)+"s");
		if(mode==MODE_BCDLL5){
			System.err.println("Bank promotions: "+totalBankPromotions+
				" (avg "+String.format("%.1f", (double)totalBankPromotions/outerTrials)+" per trial)");
			if(totalBankPromotions==0){
				System.err.println("WARNING: Zero bank promotions! inner="+innerSteps+
					" may be too small. Try inner=262144 or larger.");
			}
		}

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
			if(geoAvg){
				double logSum=0; int n=0;
				for(int h=0; h<NUM_HIST; h++){
					if(totalCounts[t][h]>0){
						final double stateGeoMean=Math.exp(totalGeoSums[t][h]/totalCounts[t][h]);
						logSum+=Math.log(stateGeoMean); n++;
					}
				}
				tierMean[t]=(n>0) ? Math.exp(logSum/n) : 0;
				for(int h=0; h<NUM_HIST; h++){
					if(totalCounts[t][h]>0 && tierMean[t]>0){
						final double stateGeoMean=Math.exp(totalGeoSums[t][h]/totalCounts[t][h]);
						cf[t][h]=Math.log(stateGeoMean/tierMean[t])/Math.log(2.0);
					}
				}
			}else{
				double sumCard=0; long sumCount=0;
				for(int h=0; h<NUM_HIST; h++){
					sumCard+=totalCardSums[t][h];
					sumCount+=totalCounts[t][h];
				}
				tierMean[t]=(sumCount>0) ? sumCard/sumCount : 0;
				for(int h=0; h<NUM_HIST; h++){
					if(totalCounts[t][h]>10 && tierMean[t]>0){
						cf[t][h]=Math.log(meanCard[t][h]/tierMean[t])/Math.log(2.0);
					}
				}
			}
		}

		// Steady-state CF (average of high tiers)
		final int ssStart=(mode==MODE_HISTORY) ? Math.max(8, maxTier-3) : 3;
		final double[] ssCF=new double[NUM_HIST];
		int ssTiers=0;
		for(int t=ssStart; t<=maxTier; t++){
			boolean allNonzero=true;
			for(int h=0; h<NUM_HIST; h++){if(totalCounts[t][h]<100){allNonzero=false; break;}}
			if(allNonzero){
				for(int h=0; h<NUM_HIST; h++){ssCF[h]+=cf[t][h];}
				ssTiers++;
			}
		}
		if(ssTiers>0){for(int h=0; h<NUM_HIST; h++){ssCF[h]/=ssTiers;}}

		// Per-tier detailed report (to stderr)
		System.err.println("\nPer-tier results:");
		StringBuilder hdr=new StringBuilder("Tier\tTotal\tAvgCard");
		for(int h=0; h<NUM_HIST; h++){hdr.append("\tP("+h+")\tCF("+h+")");}
		System.err.println(hdr);
		for(int t=0; t<=maxTier; t++){
			long tierTotal=0;
			for(int h=0; h<NUM_HIST; h++){tierTotal+=totalCounts[t][h];}
			if(tierTotal<100) continue;
			StringBuilder row=new StringBuilder();
			row.append(t).append('\t').append(tierTotal);
			row.append('\t').append(String.format("%.2f", tierMean[t]));
			for(int h=0; h<NUM_HIST; h++){
				row.append('\t').append(String.format("%.6f", totalCounts[t][h]/(double)tierTotal));
				row.append('\t').append(String.format("%+.8f", cf[t][h]));
			}
			System.err.println(row);
		}
		System.err.println("\nSteady-state CFs (tiers "+ssStart+"-"+maxTier+", "+ssTiers+" tiers):");
		StringBuilder ssSb=new StringBuilder("ss");
		for(int h=0; h<NUM_HIST; h++){ssSb.append('\t').append(String.format("%+.8f", ssCF[h]));}
		System.err.println(ssSb);

		// HSB table output (to stdout, for StateTable.loadHsbTable)
		System.out.println("# MantissaCompare3 HSB table");
		System.out.println("# mode="+MODE_NAMES[mode]+" buckets="+buckets+
			" outer="+outerTrials+" inner="+innerSteps+" maxTier="+maxTier+
			" avg="+(geoAvg?"geo":"lin"));
		System.out.println("# Columns: tier  CF(hist=00)  CF(hist=01)  CF(hist=10)  CF(hist=11)");
		for(int t=0; t<=maxTier; t++){
			long tierTotal=0;
			for(int h=0; h<NUM_HIST; h++){tierTotal+=totalCounts[t][h];}
			if(tierTotal<100) continue;
			StringBuilder row=new StringBuilder();
			row.append(t);
			for(int h=0; h<NUM_HIST; h++){
				row.append('\t').append(String.format("%+.8f", cf[t][h]));
			}
			System.out.println(row);
		}
		StringBuilder ssRow=new StringBuilder("ss");
		for(int h=0; h<NUM_HIST; h++){
			ssRow.append('\t').append(String.format("%+.8f", ssCF[h]));
		}
		System.out.println(ssRow);

		System.err.println("\nDone. "+(geoAvg?"Geometric":"Linear")+" averaging.");
	}
}
