package cardinality;

import rand.FastRandomXoshiro;
import parse.Parse;
import parse.PreParser;

/**
 * FLL2 (FutureLogLog 2-bit) Simulator.
 *
 * Simulates a single 16-bit FLL2 word to build per-tier correction factor tables.
 * The word encodes one bank of 6 buckets sharing a 4-bit local exponent.
 *
 * Word layout: [15:12]=localExp [11:0]=6x2-bit future bitmap
 * Per-bucket state (2 bits):
 *   Bit 0 (LSB) = "seen hash with NLZ == absNLZ"     (floor hit)
 *   Bit 1 (MSB) = "seen hash with NLZ == absNLZ+1"   (future hit)
 *   All four states {00,01,10,11} are reachable.
 *
 * Promotion fires when all 6 LSBs are set: (word &amp; 0x555) == 0x555.
 * Promotion: MSBs shift right to LSB positions, MSBs clear, localExp++.
 * Max cascade = 2 (all-{11} -> all-{01} -> all-{00}).
 *
 * States are stored in 84-slot equivalence class representation.
 * The 84 classes correspond to bucket histograms (a,b,c,d) where
 * a=count(00), b=count(01), c=count(10), d=count(11), a+b+c+d=6.
 * This eliminates permutation-inflation bugs and is more cache-friendly
 * (84 doubles per tier = 672 bytes vs 4096 doubles = 32KB).
 *
 * Run with: java -ea cardinality.FLL2Simulator [iters=N] [threads=N] [maxTier=N]
 *
 * @author Brian Bushnell, Chloe
 * @date April 8, 2026
 */
public class FLL2Simulator {

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	static final int BUCKETS_PER_WORD=6;
	static final int NUM_TIERS=16;
	/** Bits 0,2,4,6,8,10 -- floor-level bits (LSB of each bucket). */
	static final int LSB_MASK=0x555;
	/** Bits 1,3,5,7,9,11 -- future-level bits (MSB of each bucket). */
	static final int MSB_MASK=0xAAA;

	/** Number of equivalence classes for 6 buckets with 4 possible states = C(9,3). */
	static final int NUM_EQUIV=84;
	/** When true, overflow hashes (delta>1) clamp to delta=1 instead of being discarded. */
	static boolean CLAMP_OVERFLOW=false;

	/**
	 * Prefix sums for idx84 computation.
	 * OFFSET_A[a] = number of compositions of 6 into 4 parts where count(00) &lt; a.
	 */
	static final int[] OFFSET_A={0, 28, 49, 64, 74, 80, 83};

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public FLL2Simulator(int iters, int threads, int maxTier, int startState){
		assert iters>0   : "iters must be positive: "+iters;
		assert threads>0 : "threads must be positive: "+threads;
		assert maxTier>=0 && maxTier<NUM_TIERS : "maxTier out of range: "+maxTier;
		this.iters=iters;
		this.threads=threads;
		this.maxTier=maxTier;
		this.startState=startState;
		this.counts=new long[NUM_TIERS][NUM_EQUIV];
		this.sums=new double[NUM_TIERS][NUM_EQUIV];
		this.promotionHist=new long[4];
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Thomas Wang's 64-bit hash finalizer. */
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

	/**
	 * Converts a 12-bit history state to its 84-slot equivalence class index.
	 * Uses bitwise popcount to extract bucket histogram (a,b,c,d) without
	 * array allocation, then maps to a lex-ordered index via OFFSET_A.
	 *
	 * @param history 12-bit history from word &amp; 0xFFF
	 * @return index in [0, 83]
	 */
	static int idx84(int history){
		final int h1=history>>>1;
		final int d=Integer.bitCount(history&h1&LSB_MASK);       // state 11: both bits set
		final int a=6-Integer.bitCount((history|h1)&LSB_MASK);   // state 00: neither set
		final int c=Integer.bitCount(~history&h1&LSB_MASK);      // state 10: MSB only
		final int b=6-a-c-d;                                     // state 01: LSB only
		final int r=6-a;
		return OFFSET_A[a]+b*(r+1)-b*(b-1)/2+c;
	}

	/**
	 * Returns the bucket histogram (a,b,c,d) for a given idx84.
	 * Inverse of idx84().
	 */
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

	/** Human-readable string for an idx84, e.g. "{2,3,0,1}". */
	static String idx84ToString(int idx){
		final int[] h=idx84ToHistogram(idx);
		return "{"+h[0]+","+h[1]+","+h[2]+","+h[3]+"}";
	}

	/** Convert histogram (a,b,c,d) to idx84. */
	static int histogramToIdx84(int[] hist){
		final int a=hist[0], b=hist[1], c=hist[2];
		final int r=6-a;
		return OFFSET_A[a]+b*(r+1)-b*(b-1)/2+c;
	}

	/** Shifts MSBs to LSBs and increments localExp. Max cascade = 2. */
	static int promote(int word, long[] promHist){
		int promotionCount=0;
		while((word&LSB_MASK)==LSB_MASK){
			final int localExp=(word>>>12)&0xF;
			if(localExp>=15){break;}
			word=((word&MSB_MASK)>>>1)|((localExp+1)<<12);
			promotionCount++;
			assert promotionCount<=2 : "Cascade > 2: "+promotionCount;
		}
		if(promHist!=null){
			promHist[Math.min(promotionCount, promHist.length-1)]++;
		}
		return word;
	}

	/** Validates word invariants. */
	static boolean repOK(int word){
		final int localExp=(word>>>12)&0xF;
		assert localExp<=15 : "localExp="+localExp;
		assert (word&0xFFFF)==word : "Word > 16 bits";
		assert localExp>=15 || (word&LSB_MASK)!=LSB_MASK
			: "Promotable state: 0x"+Integer.toHexString(word);
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	void simulate() throws InterruptedException{
		final long[][][]   threadCounts   =new long[threads][NUM_TIERS][NUM_EQUIV];
		final double[][][] threadSums     =new double[threads][NUM_TIERS][NUM_EQUIV];
		final long[][]     threadPromHist =new long[threads][4];
		final long[]       threadOverflow =new long[threads];
		final long[]       threadUnderflow=new long[threads];
		final long[]       threadTotal    =new long[threads];

		final int base=iters/threads;
		final int rem =iters%threads;

		final Thread[] threadArr=new Thread[threads];
		for(int t=0; t<threads; t++){
			final int tid    =t;
			final int myIters=base+(tid<rem ? 1 : 0);
			final int myStart=startState&0xFFFF;

			threadArr[t]=new Thread(()->{
				final FastRandomXoshiro rng=new FastRandomXoshiro(tid+1);
				final long[][]   lCounts  =threadCounts[tid];
				final double[][] lSums    =threadSums[tid];
				final long[]     lPromHist=threadPromHist[tid];
				long lOver=0, lUnder=0, lTotal=0;

				for(int iter=0; iter<myIters; iter++){
					int word=myStart;
					long eeMask=-1L>>>((word>>>12)&0xF);

					for(long trueCard=1; trueCard<=10_000_000L; trueCard++){
						final long hash=hash64shift(rng.nextLong());
						lTotal++;

						if(Long.compareUnsigned(hash, eeMask)>0){
							lUnder++;
						}else{
							final int hashNLZ=Long.numberOfLeadingZeros(hash);
							final int tier   =(word>>>12)&0xF;
							final int delta  =hashNLZ-tier;

							if(CLAMP_OVERFLOW ? (delta<0) : (delta<0 || delta>1)){
								lOver++;
							}else{
								final int clampedDelta=Math.min(delta, 1);
								final int bucket  =(int)((hash&0x7FFF_FFFFL)%BUCKETS_PER_WORD);
								final int bitToSet=1<<(clampedDelta+bucket*2);
								if((word&bitToSet)==0){
									word|=bitToSet;
									if((word&LSB_MASK)==LSB_MASK){
										word=promote(word, lPromHist);
										eeMask=-1L>>>((word>>>12)&0xF);
									}else{
										lPromHist[0]++;
									}
								}else{
									lPromHist[0]++;
								}
							}
							assert repOK(word);
						}

						// Record into 84-slot table
						final int curTier=(word>>>12)&0xF;
						final int eqIdx  =idx84(word&0xFFF);
						lCounts[curTier][eqIdx]++;
						lSums[curTier][eqIdx]+=trueCard;

						if(curTier>maxTier){break;}
					}
				}

				threadOverflow[tid] =lOver;
				threadUnderflow[tid]=lUnder;
				threadTotal[tid]    =lTotal;
			});
			threadArr[t].start();
		}

		for(final Thread t : threadArr){t.join();}

		// Combine per-thread results
		for(int t=0; t<threads; t++){
			for(int tier=0; tier<NUM_TIERS; tier++){
				for(int s=0; s<NUM_EQUIV; s++){
					counts[tier][s]+=threadCounts[t][tier][s];
					sums[tier][s]  +=threadSums[t][tier][s];
				}
			}
			for(int i=0; i<4; i++){promotionHist[i]+=threadPromHist[t][i];}
			overflowCount +=threadOverflow[t];
			underflowCount+=threadUnderflow[t];
			totalAdds     +=threadTotal[t];
		}
	}

	/**
	 * For states with fewer than minObs observations, supplement with
	 * L1=2 neighbors in histogram space (one bucket-state change).
	 */
	void smoothSparseStates(long minObs){
		final long[][]   origCounts=new long[NUM_TIERS][NUM_EQUIV];
		final double[][] origSums  =new double[NUM_TIERS][NUM_EQUIV];
		for(int tier=0; tier<NUM_TIERS; tier++){
			System.arraycopy(counts[tier], 0, origCounts[tier], 0, NUM_EQUIV);
			System.arraycopy(sums[tier], 0, origSums[tier], 0, NUM_EQUIV);
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

	void printResults(double[] tierAvg, double[][] mult, long[][] origCounts){
		System.err.println("=== FLL2 Simulator Results ===");
		System.err.println("iters="+iters+"  threads="+threads+"  maxTier="+maxTier);
		System.err.println("totalAdds="+totalAdds);
		System.err.printf("overflowRate =%.4f  (delta>1, IOT-ignored)%n",
			totalAdds>0 ? (double)overflowCount/totalAdds : 0);
		System.err.printf("underflowRate=%.4f  (delta<0, ignored)%n",
			totalAdds>0 ? (double)underflowCount/totalAdds : 0);
		long totalProm=0;
		for(long v : promotionHist){totalProm+=v;}
		System.err.printf("Promotions per add  0=%.4f  1=%.4f  2=%.4f  3+=%.4f%n",
			totalProm>0 ? (double)promotionHist[0]/totalProm : 0,
			totalProm>0 ? (double)promotionHist[1]/totalProm : 0,
			totalProm>0 ? (double)promotionHist[2]/totalProm : 0,
			totalProm>0 ? (double)promotionHist[3]/totalProm : 0);
		System.err.println();

		System.err.printf("%-6s %-16s %-10s %-12s%n", "Tier", "AvgCard", "Growth", "Observations");
		System.out.println("# Tier averages: Tier\tAvgCard\tGrowth\tObservations");
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
		System.err.println("Equivalence classes: "+NUM_EQUIV);

		// TSV table: header with idx84 labels
		final StringBuilder header=new StringBuilder("Tier\tType");
		for(int s=0; s<NUM_EQUIV; s++){
			header.append('\t').append(idx84ToString(s));
		}
		System.out.println(header);

		// Two rows per tier: original counts and smoothed multipliers
		for(int tier=0; tier<=maxTier; tier++){
			final StringBuilder rowCounts=new StringBuilder();
			final StringBuilder rowMults =new StringBuilder();
			rowCounts.append(tier).append("\tcounts");
			rowMults.append(tier).append("\tmult");
			for(int s=0; s<NUM_EQUIV; s++){
				rowCounts.append('\t').append(origCounts[tier][s]);
				rowMults.append('\t').append(String.format("%.8f", mult[tier][s]));
			}
			System.out.println(rowCounts);
			System.out.println(rowMults);
		}
	}

	/**
	 * Runs fresh single-threaded simulations, computing estimate at each step
	 * using the provided CF table. Reports signed and absolute error.
	 * Validates that the CF table produces correct estimates in the
	 * exact same simulation environment that generated it.
	 */
	static void validateTable(double[] tierAvg, double[][] mult, int valIters, int maxTier){
		System.err.println();
		System.err.println("=== CF Table Validation (single-word estimate) ===");
		System.err.println("valIters="+valIters+"  maxTier="+maxTier);
		final int cfTiers=mult.length;

		final FastRandomXoshiro rng=new FastRandomXoshiro(999);
		double sumSignedErr=0, sumAbsErr=0;
		long totalObs=0;

		// Per-tier stats
		final double[] tierSignedSum=new double[NUM_TIERS];
		final double[] tierAbsSum=new double[NUM_TIERS];
		final long[] tierObs=new long[NUM_TIERS];

		for(int iter=0; iter<valIters; iter++){
			int word=0;
			long eeMask=-1L;

			for(long trueCard=1; trueCard<=10_000_000L; trueCard++){
				final long hash=hash64shift(rng.nextLong());

				if(Long.compareUnsigned(hash, eeMask)>0){
					// below floor, no state change
				}else{
					final int hashNLZ=Long.numberOfLeadingZeros(hash);
					final int tier   =(word>>>12)&0xF;
					final int delta  =hashNLZ-tier;

					if(delta<=1 && delta>=0){
						final int bucket  =(int)((hash&0x7FFF_FFFFL)%BUCKETS_PER_WORD);
						final int bitToSet=1<<(delta+bucket*2);
						if((word&bitToSet)==0){
							word|=bitToSet;
							if((word&LSB_MASK)==LSB_MASK){
								word=promote(word, null);
								eeMask=-1L>>>((word>>>12)&0xF);
							}
						}
					}
				}

				// Compute estimate using CF table
				final int curTier=(word>>>12)&0xF;
				final int history=word&0xFFF;
				final int eqIdx  =idx84(history);
				final int cfTier =Math.min(curTier, cfTiers-1);

				final double estimate;
				if(word==0){
					estimate=0;// match estimator's zero-word skip
				}else if(curTier<tierAvg.length){
					estimate=tierAvg[curTier]*mult[cfTier][eqIdx];
				}else{
					// extrapolate tierAvg
					final double lastGrowth=tierAvg[tierAvg.length-1]/tierAvg[tierAvg.length-2];
					final double base2=tierAvg[tierAvg.length-1];
					final int extra=curTier-(tierAvg.length-1);
					estimate=base2*Math.pow(lastGrowth, extra)*mult[cfTier][eqIdx];
				}

				final double err=(estimate-trueCard)/(double)trueCard;
				sumSignedErr+=err;
				sumAbsErr+=Math.abs(err);
				totalObs++;
				tierSignedSum[curTier]+=err;
				tierAbsSum[curTier]+=Math.abs(err);
				tierObs[curTier]++;

				if(curTier>maxTier){break;}
			}
		}

		System.err.printf("Overall: signedErr=%.6f  absErr=%.6f  observations=%d%n",
			sumSignedErr/totalObs, sumAbsErr/totalObs, totalObs);
		System.err.printf("%-6s %-14s %-14s %-12s%n", "Tier", "SignedErr", "AbsErr", "Observations");
		for(int t=0; t<NUM_TIERS; t++){
			if(tierObs[t]==0){continue;}
			System.err.printf("%-6d %-14.6f %-14.6f %-12d%n",
				t, tierSignedSum[t]/tierObs[t], tierAbsSum[t]/tierObs[t], tierObs[t]);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws InterruptedException{
		{PreParser pp=new PreParser(args, null, false); args=pp.args;}
		int  iters     =10000;
		int  threads   =8;
		int  maxTier   =15;
		long minObs    =100;
		int  startState=0;
		boolean validate=false;
		int  valIters  =1000;

		for(final String arg : args){
			final String[] kv=arg.split("=", 2);
			if(kv.length!=2){throw new RuntimeException("Unknown parameter '"+arg+"'");}
			switch(kv[0]){
				case "iters":      iters     =Parse.parseIntKMG(kv[1]);   break;
				case "threads":    threads   =Integer.parseInt(kv[1]);    break;
				case "maxTier":    maxTier   =Integer.parseInt(kv[1]);    break;
				case "minObs":     minObs    =Parse.parseKMG(kv[1]);      break;
				case "startState": startState=Integer.decode(kv[1]);      break;
				case "validate":   validate  =Boolean.parseBoolean(kv[1]); break;
				case "valIters":   valIters  =Parse.parseIntKMG(kv[1]);   break;
				case "clampoverflow": CLAMP_OVERFLOW=Boolean.parseBoolean(kv[1]); break;
				default: throw new RuntimeException("Unknown parameter '"+arg+"'");
			}
		}

		System.err.println("FLL2Simulator: iters="+iters+" threads="+threads
			+" maxTier="+maxTier+" minObs="+minObs
			+" startState=0x"+Integer.toHexString(startState));

		final long t0=System.currentTimeMillis();
		final FLL2Simulator sim=new FLL2Simulator(iters, threads, maxTier, startState);
		sim.simulate();
		final long t1=System.currentTimeMillis();
		System.err.println("Simulation time: "+(t1-t0)+" ms");

		// Snapshot raw counts before smoothing
		final long[][] origCounts=new long[NUM_TIERS][NUM_EQUIV];
		for(int t=0; t<NUM_TIERS; t++){
			System.arraycopy(sim.counts[t], 0, origCounts[t], 0, NUM_EQUIV);
		}

		final double[] tierAvg=sim.buildTierAverages();
		sim.smoothSparseStates(minObs);
		final double[][] mult=sim.buildMultipliers(tierAvg);
		sim.printResults(tierAvg, mult, origCounts);

		if(validate){
			validateTable(tierAvg, mult, valIters, maxTier);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final int iters;
	final int threads;
	final int maxTier;
	final int startState;

	long[][]   counts;
	double[][] sums;

	long[] promotionHist;
	long overflowCount;
	long underflowCount;
	long totalAdds;

}
