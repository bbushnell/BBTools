package cardinality;

import rand.FastRandomXoshiro;
import parse.Parse;
import parse.PreParser;

/**
 * BCLL (BankedCeilingLogLog) Simulator.
 *
 * Simulates a single 16-bit BCLL word to build per-tier correction factor tables.
 * Same physical layout as FLL2 (6 x 2-bit registers + 4-bit localExp per word),
 * but the localExp represents a CEILING (max NLZ) instead of a floor.
 *
 * Update rule:
 *   delta = hashNLZ - ceiling
 *   delta > 0  : raise ceiling. delta==1: MSBs to LSBs. delta>=2: wipe.
 *                Set MSB for arriving bucket.
 *   delta == 0 : set MSB (ceiling hit)
 *   delta == -1: set LSB (one below ceiling)
 *   delta <= -2: ignore
 *   Overflow: if newLocalExp > 15, ignore hash entirely.
 *
 * No promote() method. Ceiling advancement is inline in the simulation loop.
 *
 * Run: java -ea cardinality.BCLLSimulator [iters=N] [threads=N] [maxTier=N]
 *
 * @author Brian Bushnell, Chloe
 * @date April 9, 2026
 */
public class BCLLSimulator {

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	static final int BUCKETS_PER_WORD=6;
	static final int NUM_TIERS=16;
	static final int LSB_MASK=0x555;
	static final int MSB_MASK=0xAAA;
	static final int NUM_EQUIV=84;
	static final int[] OFFSET_A={0, 28, 49, 64, 74, 80, 83};

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public BCLLSimulator(int iters, int threads, int maxTier){
		this.iters=iters;
		this.threads=threads;
		this.maxTier=maxTier;
		this.counts=new long[NUM_TIERS][NUM_EQUIV];
		this.sums=new double[NUM_TIERS][NUM_EQUIV];
		this.sumSq=new double[NUM_TIERS][NUM_EQUIV];
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	static long hash64shift(long key){
		key=(~key)+(key<<21);
		key^=(key>>>24);
		key+=(key<<3)+(key<<8);
		key^=(key>>>14);
		key+=(key<<2)+(key<<4);
		key^=(key>>>28);
		key+=(key<<31);
		return key;
	}

	static int idx84(int history){
		final int h1=history>>>1;
		final int d=Integer.bitCount(history&h1&LSB_MASK);
		final int a=6-Integer.bitCount((history|h1)&LSB_MASK);
		final int c=Integer.bitCount(~history&h1&LSB_MASK);
		final int b=6-a-c-d;
		final int r=6-a;
		return OFFSET_A[a]+b*(r+1)-b*(b-1)/2+c;
	}

	static int[] idx84ToHistogram(int idx){
		for(int a=0; a<=6; a++){
			final int base=OFFSET_A[a];
			final int nextBase=(a<6) ? OFFSET_A[a+1] : NUM_EQUIV;
			if(idx<nextBase){
				final int r=6-a;
				final int within=idx-base;
				for(int b=0; b<=r; b++){
					final int bOffset=b*(r+1)-b*(b-1)/2;
					final int maxC=r-b;
					if(within>=bOffset && within<bOffset+maxC+1){
						final int c=within-bOffset;
						final int d=r-b-c;
						return new int[]{a, b, c, d};
					}
				}
			}
		}
		throw new IllegalArgumentException("Invalid idx84: "+idx);
	}

	static String idx84ToString(int idx){
		final int[] h=idx84ToHistogram(idx);
		return "{"+h[0]+","+h[1]+","+h[2]+","+h[3]+"}";
	}

	static int histogramToIdx84(int[] hist){
		final int a=hist[0], b=hist[1], c=hist[2];
		final int r=6-a;
		return OFFSET_A[a]+b*(r+1)-b*(b-1)/2+c;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	void simulate() throws InterruptedException{
		final long[][][]   threadCounts=new long[threads][NUM_TIERS][NUM_EQUIV];
		final double[][][] threadSums  =new double[threads][NUM_TIERS][NUM_EQUIV];
		final double[][][] threadSumSq =new double[threads][NUM_TIERS][NUM_EQUIV];
		final long[] threadOverflow=new long[threads];
		final long[] threadBelow   =new long[threads];
		final long[] threadTotal   =new long[threads];
		final long[] threadRaises  =new long[threads];
		final long[] threadWipes   =new long[threads];

		final int base=iters/threads;
		final int rem =iters%threads;

		final Thread[] threadArr=new Thread[threads];
		for(int t=0; t<threads; t++){
			final int tid    =t;
			final int myIters=base+(tid<rem ? 1 : 0);

			threadArr[t]=new Thread(()->{
				final FastRandomXoshiro rng=new FastRandomXoshiro(tid+1);
				final long[][]   lCounts=threadCounts[tid];
				final double[][] lSums  =threadSums[tid];
				final double[][] lSumSq =threadSumSq[tid];
				long lOver=0, lBelow=0, lTotal=0, lRaises=0, lWipes=0;

				for(int iter=0; iter<myIters; iter++){
					int word=0;// ceiling starts at 0
					// In single-word sim, localExp IS the absCeiling (no global)
					// eeMask rejects NLZ < ceiling - 1
					long eeMask=-1L;// ceiling=0 -> accept everything

					for(long trueCard=1; trueCard<=10_000_000L; trueCard++){
						final long hash=hash64shift(rng.nextLong());
						lTotal++;

						if(Long.compareUnsigned(hash, eeMask)>0){
							lBelow++;
						}else{
							final int hashNLZ=Long.numberOfLeadingZeros(hash);
							final int ceiling=(word>>>12)&0xF;
							final int delta  =hashNLZ-ceiling;

							if(delta<=-2){
								lBelow++;
							}else if(delta>0){
								// Raise ceiling
								final int newCeiling=ceiling+delta;
								if(newCeiling>15){
									lOver++;
								}else{
									final int bucket=(int)((hash&0x7FFF_FFFFL)%BUCKETS_PER_WORD);
									if(delta==1){
										// Shift: MSBs become LSBs
										word=((word&MSB_MASK)>>>1)|(newCeiling<<12);
										lRaises++;
									}else{
										// Wipe: delta >= 2
										word=(newCeiling<<12);
										lWipes++;
									}
									// Set MSB for arriving bucket
									word|=(1<<(1+bucket*2));
									// Update eeMask
									eeMask=(newCeiling<=1) ? -1L : (-1L)>>>(newCeiling-1);
								}
							}else{
								// delta == 0: set MSB; delta == -1: set LSB
								final int bucket=(int)((hash&0x7FFF_FFFFL)%BUCKETS_PER_WORD);
								final int bitToSet=(delta==0)
									? (1<<(1+bucket*2))
									: (1<<(bucket*2));
								word|=bitToSet;
							}
						}

						// Record into 84-slot table
						final int curTier=(word>>>12)&0xF;
						final int eqIdx  =idx84(word&0xFFF);
						lCounts[curTier][eqIdx]++;
						lSums[curTier][eqIdx]+=trueCard;
						lSumSq[curTier][eqIdx]+=(double)trueCard*trueCard;

						if(curTier>maxTier){break;}
					}
				}

				threadOverflow[tid]=lOver;
				threadBelow[tid]   =lBelow;
				threadTotal[tid]   =lTotal;
				threadRaises[tid]  =lRaises;
				threadWipes[tid]   =lWipes;
			});
			threadArr[t].start();
		}

		for(final Thread t : threadArr){t.join();}

		for(int t=0; t<threads; t++){
			for(int tier=0; tier<NUM_TIERS; tier++){
				for(int s=0; s<NUM_EQUIV; s++){
					counts[tier][s]+=threadCounts[t][tier][s];
					sums[tier][s]  +=threadSums[t][tier][s];
					sumSq[tier][s] +=threadSumSq[t][tier][s];
				}
			}
			overflowCount+=threadOverflow[t];
			belowCount   +=threadBelow[t];
			totalAdds    +=threadTotal[t];
			ceilingRaises+=threadRaises[t];
			wipeCount    +=threadWipes[t];
		}
	}

	/** Fills sparse states by borrowing from L1-distance neighbors in histogram space. */
	void smoothSparseStates(long minObs){
		final long[][]   origCounts=new long[NUM_TIERS][NUM_EQUIV];
		final double[][] origSums  =new double[NUM_TIERS][NUM_EQUIV];
		final double[][] origSumSq =new double[NUM_TIERS][NUM_EQUIV];
		for(int tier=0; tier<NUM_TIERS; tier++){
			System.arraycopy(counts[tier], 0, origCounts[tier], 0, NUM_EQUIV);
			System.arraycopy(sums[tier], 0, origSums[tier], 0, NUM_EQUIV);
			System.arraycopy(sumSq[tier], 0, origSumSq[tier], 0, NUM_EQUIV);
		}

		for(int tier=0; tier<NUM_TIERS; tier++){
			for(int s=0; s<NUM_EQUIV; s++){
				if(origCounts[tier][s]>=minObs){continue;}
				final int[] hist=idx84ToHistogram(s);

				for(int from=0; from<4; from++){
					if(hist[from]==0){continue;}
					hist[from]--;
					for(int to=0; to<4; to++){
						if(to==from){continue;}
						hist[to]++;
						final int nIdx=histogramToIdx84(hist);
						counts[tier][s]+=origCounts[tier][nIdx];
						sums[tier][s]  +=origSums[tier][nIdx];
						sumSq[tier][s] +=origSumSq[tier][nIdx];
						hist[to]--;
					}
					hist[from]++;
				}
			}
		}
	}

	double[] buildTierAverages(){
		final double[] tierAvg=new double[NUM_TIERS];
		for(int tier=0; tier<NUM_TIERS; tier++){
			long totalCount=0;
			double totalSum=0;
			for(int s=0; s<NUM_EQUIV; s++){
				totalCount+=counts[tier][s];
				totalSum  +=sums[tier][s];
			}
			tierAvg[tier]=totalCount>0 ? totalSum/totalCount : 0;
		}
		return tierAvg;
	}

	double[][] buildMultipliers(double[] tierAvg){
		final double[][] mult=new double[NUM_TIERS][NUM_EQUIV];
		for(int tier=0; tier<NUM_TIERS; tier++){
			for(int s=0; s<NUM_EQUIV; s++){mult[tier][s]=1.0;}
			if(tierAvg[tier]<=0){continue;}
			for(int s=0; s<NUM_EQUIV; s++){
				if(counts[tier][s]>0){
					mult[tier][s]=(sums[tier][s]/counts[tier][s])/tierAvg[tier];
				}
			}
		}
		return mult;
	}

	/** Build CV (coefficient of variation) table: CV = stddev / mean per state. */
	double[][] buildCV(){
		final double[][] cv=new double[NUM_TIERS][NUM_EQUIV];
		for(int tier=0; tier<NUM_TIERS; tier++){
			for(int s=0; s<NUM_EQUIV; s++){
				cv[tier][s]=1.0;// default: no information
				if(counts[tier][s]>1){
					final double mean=sums[tier][s]/counts[tier][s];
					if(mean>0){
						final double variance=sumSq[tier][s]/counts[tier][s]-mean*mean;
						if(variance>0){
							cv[tier][s]=Math.sqrt(variance)/mean;
						}else{
							cv[tier][s]=0.001;// near-zero variance
						}
					}
				}
			}
		}
		return cv;
	}

	void printResults(double[] tierAvg, double[][] mult, double[][] cv, long[][] origCounts){
		System.err.println("=== BCLL Simulator Results ===");
		System.err.println("iters="+iters+"  threads="+threads+"  maxTier="+maxTier);
		System.err.println("totalAdds="+totalAdds);
		System.err.printf("overflowRate =%.6f  (localExp would exceed 15)%n",
			totalAdds>0 ? (double)overflowCount/totalAdds : 0);
		System.err.printf("belowRate    =%.6f  (delta<=-2 or below eeMask)%n",
			totalAdds>0 ? (double)belowCount/totalAdds : 0);
		System.err.printf("ceilingRaises=%d  (delta==1 shifts)%n", ceilingRaises);
		System.err.printf("wipeCount    =%d  (delta>=2 wipes)%n", wipeCount);
		System.err.println();

		System.err.printf("%-6s %-16s %-10s %-12s%n", "Tier", "AvgCard", "Growth", "Observations");
		System.out.println("# BCLL Tier averages: Tier\tAvgCard\tGrowth\tObservations");
		for(int tier=0; tier<=maxTier; tier++){
			long obs=0;
			for(int s=0; s<NUM_EQUIV; s++){obs+=origCounts[tier][s];}
			final double growth=(tier>0 && tierAvg[tier-1]>0)
				? tierAvg[tier]/tierAvg[tier-1] : Double.NaN;
			System.err.printf("%-6d %-16.2f %-10.4f %-12d%n",
				tier, tierAvg[tier], growth, obs);
			System.out.printf("tierAvg\t%d\t%.8f\t%.8f\t%d%n",
				tier, tierAvg[tier], growth, obs);
		}
		System.err.println();

		final StringBuilder header=new StringBuilder("Tier\tType");
		for(int s=0; s<NUM_EQUIV; s++){
			header.append('\t').append(idx84ToString(s));
		}
		System.out.println(header);

		for(int tier=0; tier<=maxTier; tier++){
			final StringBuilder rowCounts=new StringBuilder();
			final StringBuilder rowMults =new StringBuilder();
			final StringBuilder rowCV    =new StringBuilder();
			rowCounts.append(tier).append("\tcounts");
			rowMults.append(tier).append("\tmult");
			rowCV.append(tier).append("\tcv");
			for(int s=0; s<NUM_EQUIV; s++){
				rowCounts.append('\t').append(origCounts[tier][s]);
				rowMults.append('\t').append(String.format("%.8f", mult[tier][s]));
				rowCV.append('\t').append(String.format("%.8f", cv[tier][s]));
			}
			System.out.println(rowCounts);
			System.out.println(rowMults);
			System.out.println(rowCV);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws InterruptedException{
		{PreParser pp=new PreParser(args, null, false); args=pp.args;}
		int  iters  =10000;
		int  threads=8;
		int  maxTier=12;
		long minObs =100;

		for(final String arg : args){
			final String[] kv=arg.split("=", 2);
			if(kv.length!=2){throw new RuntimeException("Unknown parameter '"+arg+"'");}
			switch(kv[0]){
				case "iters":   iters  =Parse.parseIntKMG(kv[1]); break;
				case "threads": threads=Integer.parseInt(kv[1]);   break;
				case "maxTier": maxTier=Integer.parseInt(kv[1]);   break;
				case "minObs":  minObs =Parse.parseKMG(kv[1]);     break;
				default: throw new RuntimeException("Unknown parameter '"+arg+"'");
			}
		}

		System.err.println("BCLLSimulator: iters="+iters+" threads="+threads
			+" maxTier="+maxTier+" minObs="+minObs);

		final long t0=System.currentTimeMillis();
		final BCLLSimulator sim=new BCLLSimulator(iters, threads, maxTier);
		sim.simulate();
		final long t1=System.currentTimeMillis();
		System.err.println("Simulation time: "+(t1-t0)+" ms");

		final long[][] origCounts=new long[NUM_TIERS][NUM_EQUIV];
		for(int t=0; t<NUM_TIERS; t++){
			System.arraycopy(sim.counts[t], 0, origCounts[t], 0, NUM_EQUIV);
		}

		final double[] tierAvg=sim.buildTierAverages();
		sim.smoothSparseStates(minObs);
		final double[][] mult=sim.buildMultipliers(tierAvg);
		final double[][] cv=sim.buildCV();
		sim.printResults(tierAvg, mult, cv, origCounts);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final int iters;
	final int threads;
	final int maxTier;

	long[][]   counts;
	double[][] sums;
	double[][] sumSq;

	long overflowCount;
	long belowCount;
	long totalAdds;
	long ceilingRaises;
	long wipeCount;

}
