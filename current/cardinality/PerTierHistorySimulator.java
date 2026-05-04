package cardinality;

import java.io.PrintStream;
import java.util.concurrent.atomic.AtomicLong;
import jdk.incubator.vector.*;
import parse.Parse;
import parse.PreParser;
import rand.FastRandomXoshiro;

/**
 * Per-tier history state simulator for TTLL-family estimators.
 *
 * Simulates single buckets receiving random hashes until the bucket reaches
 * stopTier. Records (tier, combined_h, trueCard) under 4 sampling modes and
 * 3 averaging modes (12 combos). Self-tests by looking up the trained table
 * at 1%-spaced checkpoints and measuring signed/abs error per tier.
 *
 * Each simulation run is seeded by its run number via AtomicLong counter,
 * so results are deterministic regardless of thread count.
 *
 * @author Chloe, Brian
 */
public class PerTierHistorySimulator {

	static final int S_ENTRY=0, S_ALL=1, S_ENTRYEXIT=2, S_PCTILE=3;
	static final int A_LIN=0, A_GEO=1, A_HARM=2;
	static final int NS=4, NA=3;
	static final String[] SNAME={"entry","all","entryexit","pctile"};
	static final String[] ANAME={"lin","geo","harm"};

	static boolean COMPRESSED=false;
	static boolean EXPANDED=false;
	static boolean USE_SIMD=true;
	static int TIER_NUMER=1, TIER_DENOM=1;
	static int NUM_TAILS=2;
	static int PS_TARGET=0; // 0=harmonic(default), 1=arithmetic, 2=geometric
	static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	final int tailBits, numTails, tailWidth, tailMask, numStates, numCondensed, margin, stopTier, threads;

	long[][][] counts;
	double[][][] sums, sumsLog, sumsInv, sumsSquared;

	PerTierHistorySimulator(int tailBits_, int numTails_, int stopTier_, int threads_){
		tailBits=tailBits_; numTails=numTails_;
		tailWidth=tailBits/numTails; tailMask=(1<<tailWidth)-1;
		numStates=1<<tailBits; margin=tailWidth-1;
		int nc=1; for(int d=0; d<tailWidth; d++){nc*=(numTails+1);} numCondensed=nc;
		stopTier=stopTier_; threads=threads_;
		counts=new long[NS][stopTier][numStates];
		sums=new double[NS][stopTier][numStates];
		sumsLog=new double[NS][stopTier][numStates];
		sumsInv=new double[NS][stopTier][numStates];
		sumsSquared=new double[NS][stopTier][numStates];
	}

	/*--------------------------------------------------------------*/
	/*----------------     Single-Bucket Simulation  ----------------*/
	/*--------------------------------------------------------------*/

	/** Simulates one bucket from trueCard=1 until exp reaches stopTier.
	 *  Records into the provided accumulators under all 4 sampling modes. */
	void simulateOneBucket(long seed,
			long[][][] lc, double[][][] ls, double[][][] lsl, double[][][] lsi, double[][][] lss){
		final FastRandomXoshiro rng=new FastRandomXoshiro(seed);
		int exp=0, tail=0;
		long eeMask=-1L;
		int prevTier=-1, prevState=-1;
		double nextPctile=1.0;

		for(long tc=1; ; tc++){
			final long hash=rng.nextLong();

			if(Long.compareUnsigned(hash, eeMask)>0){
				record(lc[S_ALL], ls[S_ALL], lsl[S_ALL], lsi[S_ALL], lss[S_ALL], exp, tail, tc);
				if(tc>=nextPctile){
					record(lc[S_PCTILE], ls[S_PCTILE], lsl[S_PCTILE], lsi[S_PCTILE], lss[S_PCTILE], exp, tail, tc);
					nextPctile=Math.max(tc+1, nextPctile*1.01);
				}
				continue;
			}

			final int rawNLZ=Long.numberOfLeadingZeros(hash);
			final int hashNLZ;
			if(TIER_NUMER!=1){
				final int mBits;
				if(rawNLZ>=43){mBits=(int)((hash<<(rawNLZ-42))&0xFFFFF);}
				else{mBits=(int)((hash>>>(42-rawNLZ))&0xFFFFF);}
				final int mantissa=(mBits>=MANTISSA_THRESHOLD)?1:0;
				hashNLZ=(TIER_NUMER*rawNLZ+mantissa)/TIER_DENOM;
			}else{
				hashNLZ=rawNLZ;
			}
			final int delta=hashNLZ-exp;
			final int tailIdx=(int)((hash>>>1)&(numTails-1));

			if(delta<-margin){
				record(lc[S_ALL], ls[S_ALL], lsl[S_ALL], lsi[S_ALL], lss[S_ALL], exp, tail, tc);
				if(tc>=nextPctile){
					record(lc[S_PCTILE], ls[S_PCTILE], lsl[S_PCTILE], lsi[S_PCTILE], lss[S_PCTILE], exp, tail, tc);
					nextPctile=Math.max(tc+1, nextPctile*1.01);
				}
				continue;
			}

			if(delta<=0){
				tail|=(1<<(tailIdx*tailWidth+(tailWidth-1)+delta));
			}else{
				if(exp+delta>=stopTier){break;}
				final int shiftAmt=Math.min(delta, tailWidth);
				int newTail=0;
				for(int ti=0; ti<numTails; ti++){
					int t2=((tail>>>(ti*tailWidth))&tailMask)>>>shiftAmt;
					if(ti==tailIdx){t2|=(1<<(tailWidth-1));}
					newTail|=(t2<<(ti*tailWidth));
				}
				tail=newTail;
				exp+=delta;
				{
					final int target=exp-margin;
					if(target<=0){eeMask=-1L;}
					else{
						final int minRaw=(TIER_DENOM*target)/TIER_NUMER;
						eeMask=(minRaw>=64)?0:(-1L)>>>minRaw;
					}
				}
			}

			final boolean changed=(exp!=prevTier || tail!=prevState);

			record(lc[S_ALL], ls[S_ALL], lsl[S_ALL], lsi[S_ALL], lss[S_ALL], exp, tail, tc);

			if(changed){
				record(lc[S_ENTRY], ls[S_ENTRY], lsl[S_ENTRY], lsi[S_ENTRY], lss[S_ENTRY], exp, tail, tc);
			}

			if(changed && prevTier>=0){
				record(lc[S_ENTRYEXIT], ls[S_ENTRYEXIT], lsl[S_ENTRYEXIT], lsi[S_ENTRYEXIT], lss[S_ENTRYEXIT],
					prevTier, prevState, tc);
			}
			if(changed){
				record(lc[S_ENTRYEXIT], ls[S_ENTRYEXIT], lsl[S_ENTRYEXIT], lsi[S_ENTRYEXIT], lss[S_ENTRYEXIT],
					exp, tail, tc);
			}

			if(tc>=nextPctile){
				record(lc[S_PCTILE], ls[S_PCTILE], lsl[S_PCTILE], lsi[S_PCTILE], lss[S_PCTILE], exp, tail, tc);
				nextPctile=Math.max(tc+1, nextPctile*1.01);
			}

			prevTier=exp; prevState=tail;
		}
	}

	static void record(long[][] c, double[][] s, double[][] sl, double[][] si, double[][] ss,
			int tier, int state, long tc){
		c[tier][state]++;
		s[tier][state]+=tc;
		sl[tier][state]+=Math.log(tc);
		si[tier][state]+=1.0/tc;
		ss[tier][state]+=(double)tc*tc;
	}

	/** Entry-only simulation: only records on state changes. ~10x faster at high tiers. */
	void simulateOneBucketEntry(long seed,
			long[][] lc, double[][] ls, double[][] lsl, double[][] lsi, double[][] lss){
		final FastRandomXoshiro rng=new FastRandomXoshiro(seed);
		int exp=0, tail=0;
		long eeMask=-1L;

		for(long tc=1; ; tc++){
			final long hash=rng.nextLong();

			if(Long.compareUnsigned(hash, eeMask)>0){continue;}

			final int rawNLZ=Long.numberOfLeadingZeros(hash);
			final int hashNLZ;
			if(TIER_NUMER!=1){
				final int mBits;
				if(rawNLZ>=43){mBits=(int)((hash<<(rawNLZ-42))&0xFFFFF);}
				else{mBits=(int)((hash>>>(42-rawNLZ))&0xFFFFF);}
				hashNLZ=(TIER_NUMER*rawNLZ+((mBits>=MANTISSA_THRESHOLD)?1:0))/TIER_DENOM;
			}else{hashNLZ=rawNLZ;}
			final int delta=hashNLZ-exp;
			final int tailIdx=(int)((hash>>>1)&(numTails-1));

			if(delta<-margin){continue;}

			final int oldTail=tail;
			if(delta<=0){
				tail|=(1<<(tailIdx*tailWidth+(tailWidth-1)+delta));
			}else{
				if(exp+delta>=stopTier){break;}
				final int shiftAmt=Math.min(delta, tailWidth);
				int newTail=0;
				for(int ti=0; ti<numTails; ti++){
					int t2=((tail>>>(ti*tailWidth))&tailMask)>>>shiftAmt;
					if(ti==tailIdx){t2|=(1<<(tailWidth-1));}
					newTail|=(t2<<(ti*tailWidth));
				}
				tail=newTail;
				exp+=delta;
				final int target=exp-margin;
				if(target<=0){eeMask=-1L;}
				else{
					final int minRaw=(TIER_DENOM*target)/TIER_NUMER;
					eeMask=(minRaw>=64)?0:(-1L)>>>minRaw;
				}
			}

			if(tail!=oldTail || delta>0){
				lc[exp][tail]++;
				ls[exp][tail]+=tc;
				lsl[exp][tail]+=Math.log(tc);
				lsi[exp][tail]+=1.0/tc;
				lss[exp][tail]+=(double)tc*tc;
			}
		}
	}

	private static final VectorSpecies<Long> SPECIES=LongVector.SPECIES_256;
	private static final LongVector FLIP_BIT=LongVector.broadcast(SPECIES, Long.MIN_VALUE);

	/** Single-bucket SIMD: 4 RNG streams feed one bucket. 4× hash throughput, zero tail waste. */
	void simulateOneBucketEntrySIMD(long seed,
			long[][] lc, double[][] ls, double[][] lsl, double[][] lsi, double[][] lss){

		// 4 independent xoshiro256+ streams for one bucket
		long[] ss0=new long[4], ss1=new long[4], ss2=new long[4], ss3=new long[4];
		for(int i=0; i<4; i++){
			long x=seed*4+i;
			x=splitmix(x); ss0[i]=x;
			x=splitmix(x); ss1[i]=x;
			x=splitmix(x); ss2[i]=x;
			x=splitmix(x); ss3[i]=x;
		}
		LongVector vs0=LongVector.fromArray(SPECIES, ss0, 0);
		LongVector vs1=LongVector.fromArray(SPECIES, ss1, 0);
		LongVector vs2=LongVector.fromArray(SPECIES, ss2, 0);
		LongVector vs3=LongVector.fromArray(SPECIES, ss3, 0);

		for(int w=0; w<4; w++){
			LongVector t=vs1.lanewise(VectorOperators.LSHL, 17);
			vs2=vs2.lanewise(VectorOperators.XOR, vs0);
			vs3=vs3.lanewise(VectorOperators.XOR, vs1);
			vs1=vs1.lanewise(VectorOperators.XOR, vs2);
			vs0=vs0.lanewise(VectorOperators.XOR, vs3);
			vs2=vs2.lanewise(VectorOperators.XOR, t);
			vs3=vs3.lanewise(VectorOperators.ROL, 45);
		}

		int exp=0, tail=0;
		long eeMask=-1L;
		long tc=0;
		LongVector mFlip=LongVector.broadcast(SPECIES, eeMask^Long.MIN_VALUE);

		for(;;){
			// SIMD INNER LOOP: generate 4 hashes, test all against single eeMask
			LongVector result;
			long inner=0;
			do{
				result=vs0.add(vs3);
				LongVector t=vs1.lanewise(VectorOperators.LSHL, 17);
				vs2=vs2.lanewise(VectorOperators.XOR, vs0);
				vs3=vs3.lanewise(VectorOperators.XOR, vs1);
				vs1=vs1.lanewise(VectorOperators.XOR, vs2);
				vs0=vs0.lanewise(VectorOperators.XOR, vs3);
				vs2=vs2.lanewise(VectorOperators.XOR, t);
				vs3=vs3.lanewise(VectorOperators.ROL, 45);
				inner++;
			}while(result.lanewise(VectorOperators.XOR, FLIP_BIT)
				.compare(VectorOperators.GT, mFlip).allTrue());

			// Account for (inner-1) full-skip iterations × 4 hashes each
			tc+=(inner-1)*4;

			// Process 4 hashes from the last iteration sequentially into one bucket
			final long h0=result.lane(0), h1=result.lane(1), h2=result.lane(2), h3=result.lane(3);
			final long[] hh={h0, h1, h2, h3};
			for(int i=0; i<4; i++){
				tc++;
				final long hash=hh[i];
				if(Long.compareUnsigned(hash, eeMask)>0){continue;}

				final int rawNLZ=Long.numberOfLeadingZeros(hash);
				final int hashNLZ=(TIER_NUMER!=1)?(TIER_NUMER*rawNLZ+((((rawNLZ>=43)?(int)((hash<<(rawNLZ-42))&0xFFFFF):(int)((hash>>>(42-rawNLZ))&0xFFFFF))>=MANTISSA_THRESHOLD)?1:0))/TIER_DENOM:rawNLZ;
				final int delta=hashNLZ-exp;
				if(delta<-margin){continue;}

				final int tailIdx=(int)((hash>>>1)&(numTails-1));
				final int oldTail=tail;
				if(delta<=0){
					tail|=(1<<(tailIdx*tailWidth+(tailWidth-1)+delta));
				}else{
					if(exp+delta>=stopTier){return;}
					final int shiftAmt=Math.min(delta, tailWidth);
					int newTail=0;
					for(int ti=0; ti<numTails; ti++){
						int t2=((tail>>>(ti*tailWidth))&tailMask)>>>shiftAmt;
						if(ti==tailIdx){t2|=(1<<(tailWidth-1));}
						newTail|=(t2<<(ti*tailWidth));
					}
					tail=newTail; exp+=delta;
					final int target=exp-margin;
					if(target<=0){eeMask=-1L;}
					else{
						final int minRaw=(TIER_DENOM*target)/TIER_NUMER;
						eeMask=(minRaw>=64)?0:(-1L)>>>minRaw;
					}
					mFlip=LongVector.broadcast(SPECIES, eeMask^Long.MIN_VALUE);
				}
				if(tail!=oldTail || delta>0){
					lc[exp][tail]++; ls[exp][tail]+=tc;
					lsl[exp][tail]+=Math.log(tc); lsi[exp][tail]+=1.0/tc;
					lss[exp][tail]+=(double)tc*tc;
				}
			}
		}
	}

	private static long splitmix(long x){
		x+=0x9E3779B97F4A7C15L;
		x=(x^(x>>>30))*0xBF58476D1CE4E5B9L;
		x=(x^(x>>>27))*0x94D049BB133111EBL;
		return x^(x>>>31);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Training            ----------------*/
	/*--------------------------------------------------------------*/

	void train(int totalRuns, long seedOffset) throws InterruptedException {
		train(totalRuns, seedOffset, -1);
	}

	void train(int totalRuns, long seedOffset, int smFilter) throws InterruptedException {
		if(smFilter==S_ENTRY){
			trainEntry(totalRuns, seedOffset);
			return;
		}
		final AtomicLong counter=new AtomicLong(0);
		final long[][][][] tC=new long[threads][NS][stopTier][numStates];
		final double[][][][] tS=new double[threads][NS][stopTier][numStates];
		final double[][][][] tSL=new double[threads][NS][stopTier][numStates];
		final double[][][][] tSI=new double[threads][NS][stopTier][numStates];
		final double[][][][] tSS=new double[threads][NS][stopTier][numStates];

		final Thread[] arr=new Thread[threads];
		for(int t=0; t<threads; t++){
			final int tid=t;
			arr[t]=new Thread(()->{
				long runId;
				while((runId=counter.getAndIncrement())<totalRuns){
					simulateOneBucket(seedOffset+runId, tC[tid], tS[tid], tSL[tid], tSI[tid], tSS[tid]);
				}
			});
			arr[t].start();
		}
		for(final Thread th : arr){th.join();}

		for(int t=0; t<threads; t++){
			for(int sm=0; sm<NS; sm++){
				for(int tier=0; tier<stopTier; tier++){
					for(int s=0; s<numStates; s++){
						counts[sm][tier][s]+=tC[t][sm][tier][s];
						sums[sm][tier][s]+=tS[t][sm][tier][s];
						sumsLog[sm][tier][s]+=tSL[t][sm][tier][s];
						sumsInv[sm][tier][s]+=tSI[t][sm][tier][s];
						sumsSquared[sm][tier][s]+=tSS[t][sm][tier][s];
					}
				}
			}
		}
	}

	void trainEntry(int totalRuns, long seedOffset) throws InterruptedException {
		final AtomicLong counter=new AtomicLong(0);
		final long[][][] tC=new long[threads][stopTier][numStates];
		final double[][][] tS=new double[threads][stopTier][numStates];
		final double[][][] tSL=new double[threads][stopTier][numStates];
		final double[][][] tSI=new double[threads][stopTier][numStates];
		final double[][][] tSS=new double[threads][stopTier][numStates];

		final boolean useSIMD=USE_SIMD && totalRuns>=4;
		final Thread[] arr=new Thread[threads];
		for(int t=0; t<threads; t++){
			final int tid=t;
			arr[t]=new Thread(()->{
				long runId;
				if(useSIMD){
					while((runId=counter.getAndIncrement())<totalRuns){
						simulateOneBucketEntrySIMD(seedOffset+runId, tC[tid], tS[tid], tSL[tid], tSI[tid], tSS[tid]);
					}
				}else{
					while((runId=counter.getAndIncrement())<totalRuns){
						simulateOneBucketEntry(seedOffset+runId, tC[tid], tS[tid], tSL[tid], tSI[tid], tSS[tid]);
					}
				}
			});
			arr[t].start();
		}
		for(final Thread th : arr){th.join();}

		for(int t=0; t<threads; t++){
			for(int tier=0; tier<stopTier; tier++){
				for(int s=0; s<numStates; s++){
					counts[S_ENTRY][tier][s]+=tC[t][tier][s];
					sums[S_ENTRY][tier][s]+=tS[t][tier][s];
					sumsLog[S_ENTRY][tier][s]+=tSL[t][tier][s];
					sumsInv[S_ENTRY][tier][s]+=tSI[t][tier][s];
					sumsSquared[S_ENTRY][tier][s]+=tSS[t][tier][s];
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Build Tables           ----------------*/
	/*--------------------------------------------------------------*/

	double[][] buildTable(int sm, int am){
		final double[][] t=new double[stopTier][numStates];
		for(int tier=0; tier<stopTier; tier++){
			for(int s=0; s<numStates; s++){
				final long n=counts[sm][tier][s];
				if(n==0) continue;
				switch(am){
					case A_LIN:  t[tier][s]=sums[sm][tier][s]/n; break;
					case A_GEO:  t[tier][s]=Math.exp(sumsLog[sm][tier][s]/n); break;
					case A_HARM: t[tier][s]=(sumsInv[sm][tier][s]>0)?n/sumsInv[sm][tier][s]:0; break;
				}
			}
		}
		return t;
	}

	int condensedClass(int rawState){
		int cls=0, base=1;
		for(int d=0; d<tailWidth; d++){
			int count=0;
			for(int ti=0; ti<numTails; ti++){
				if(((rawState>>>(ti*tailWidth+d))&1)!=0) count++;
			}
			cls+=count*base;
			base*=(numTails+1);
		}
		return cls;
	}

	double[][] buildVarianceTable(int sm){
		final long[][] pc=new long[stopTier][numCondensed];
		final double[][] ps=new double[stopTier][numCondensed];
		final double[][] pss=new double[stopTier][numCondensed];
		for(int tier=0; tier<stopTier; tier++){
			for(int s=0; s<numStates; s++){
				final int cls=condensedClass(s);
				pc[tier][cls]+=counts[sm][tier][s];
				ps[tier][cls]+=sums[sm][tier][s];
				pss[tier][cls]+=sumsSquared[sm][tier][s];
			}
		}
		final double[][] v=new double[stopTier][numCondensed];
		for(int tier=0; tier<stopTier; tier++){
			for(int cls=0; cls<numCondensed; cls++){
				final long n=pc[tier][cls];
				if(n<2) continue;
				double mean=ps[tier][cls]/n;
				double meanSq=pss[tier][cls]/n;
				v[tier][cls]=Math.max(0, meanSq-mean*mean);
			}
		}
		return v;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Self-Test             ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns [tier][0=signedSum, 1=absSum, 2=count] */
	double[][] selfTest(double[][] table, int totalRuns, long seedOffset){
		final double[][] result=new double[stopTier][3];
		final AtomicLong counter=new AtomicLong(0);
		final double[][][] tResult=new double[threads][stopTier][3];

		final Thread[] arr=new Thread[threads];
		for(int t=0; t<threads; t++){
			final int tid=t;
			arr[t]=new Thread(()->{
				final double[][] lr=tResult[tid];
				long runId;
				while((runId=counter.getAndIncrement())<totalRuns){
					testOneBucket(seedOffset+runId, table, lr);
				}
			});
			arr[t].start();
		}
		try{for(final Thread th:arr){th.join();}}catch(InterruptedException e){throw new RuntimeException(e);}

		for(int t=0; t<threads; t++){
			for(int tier=0; tier<stopTier; tier++){
				result[tier][0]+=tResult[t][tier][0];
				result[tier][1]+=tResult[t][tier][1];
				result[tier][2]+=tResult[t][tier][2];
			}
		}
		return result;
	}

	void testOneBucket(long seed, double[][] table, double[][] result){
		final FastRandomXoshiro rng=new FastRandomXoshiro(seed);
		int exp=0, tail=0;
		long eeMask=-1L;
		double nextCheck=1.0;

		for(long tc=1; ; tc++){
			final long hash=rng.nextLong();
			if(Long.compareUnsigned(hash, eeMask)<=0){
				final int rawNLZ2=Long.numberOfLeadingZeros(hash);
				final int hashNLZ2;
				if(TIER_NUMER!=1){
					final int mB;
					if(rawNLZ2>=43){mB=(int)((hash<<(rawNLZ2-42))&0xFFFFF);}
					else{mB=(int)((hash>>>(42-rawNLZ2))&0xFFFFF);}
					hashNLZ2=(TIER_NUMER*rawNLZ2+((mB>=MANTISSA_THRESHOLD)?1:0))/TIER_DENOM;
				}else{hashNLZ2=rawNLZ2;}
				final int delta=hashNLZ2-exp;
				final int tailIdx2=(int)((hash>>>1)&(numTails-1));
				if(delta>=-margin){
					if(delta<=0){
						tail|=(1<<(tailIdx2*tailWidth+(tailWidth-1)+delta));
					}else{
						if(exp+delta>=stopTier){break;}
						final int shiftAmt=Math.min(delta, tailWidth);
						int newTail=0;
						for(int ti=0; ti<numTails; ti++){
							int t2=((tail>>>(ti*tailWidth))&tailMask)>>>shiftAmt;
							if(ti==tailIdx2){t2|=(1<<(tailWidth-1));}
							newTail|=(t2<<(ti*tailWidth));
						}
						tail=newTail;
						exp+=delta;
						{
							final int target=exp-margin;
							if(target<=0){eeMask=-1L;}
							else{
								final int minRaw=(TIER_DENOM*target)/TIER_NUMER;
								eeMask=(minRaw>=64)?0:(-1L)>>>minRaw;
							}
						}
					}
				}
			}

			if(tc>=nextCheck){
				if(exp<stopTier){
					final double est=table[exp][tail];
					if(est>0){
						final double err=(est-tc)/(double)tc;
						result[exp][0]+=err;
						result[exp][1]+=Math.abs(err);
						result[exp][2]++;
					}
				}
				nextCheck=Math.max(tc+1, nextCheck*1.01);
			}
		}
	}

	/** Tests multiple tables in a single simulation pass. Returns [tableIdx][tier][0=signed,1=abs,2=count] */
	double[][][] selfTestMulti(double[][][] tables, int totalRuns, long seedOffset){
		final int nt=tables.length;
		final double[][][] result=new double[nt][stopTier][3];
		final AtomicLong counter=new AtomicLong(0);
		final double[][][][] tResult=new double[threads][nt][stopTier][3];

		final boolean useSIMD=USE_SIMD;
		final Thread[] arr=new Thread[threads];
		for(int t=0; t<threads; t++){
			final int tid=t;
			arr[t]=new Thread(()->{
				final double[][][] lr=tResult[tid];
				long runId;
				while((runId=counter.getAndIncrement())<totalRuns){
					if(useSIMD){testOneBucketMultiSIMD(seedOffset+runId, tables, lr);}
					else{testOneBucketMulti(seedOffset+runId, tables, lr);}
				}
			});
			arr[t].start();
		}
		try{for(final Thread th:arr){th.join();}}catch(InterruptedException e){throw new RuntimeException(e);}

		for(int t=0; t<threads; t++)
			for(int ti=0; ti<nt; ti++)
				for(int tier=0; tier<stopTier; tier++){
					result[ti][tier][0]+=tResult[t][ti][tier][0];
					result[ti][tier][1]+=tResult[t][ti][tier][1];
					result[ti][tier][2]+=tResult[t][ti][tier][2];
				}
		return result;
	}

	void testOneBucketMulti(long seed, double[][][] tables, double[][][] result){
		final FastRandomXoshiro rng=new FastRandomXoshiro(seed);
		int exp=0, tail=0;
		long eeMask=-1L;
		double nextCheck=1.0;

		for(long tc=1; ; tc++){
			final long hash=rng.nextLong();
			if(Long.compareUnsigned(hash, eeMask)<=0){
				final int rawNLZ2=Long.numberOfLeadingZeros(hash);
				final int hashNLZ2;
				if(TIER_NUMER!=1){
					final int mB;
					if(rawNLZ2>=43){mB=(int)((hash<<(rawNLZ2-42))&0xFFFFF);}
					else{mB=(int)((hash>>>(42-rawNLZ2))&0xFFFFF);}
					hashNLZ2=(TIER_NUMER*rawNLZ2+((mB>=MANTISSA_THRESHOLD)?1:0))/TIER_DENOM;
				}else{hashNLZ2=rawNLZ2;}
				final int delta=hashNLZ2-exp;
				final int tailIdx2=(int)((hash>>>1)&(numTails-1));
				if(delta>=-margin){
					if(delta<=0){
						tail|=(1<<(tailIdx2*tailWidth+(tailWidth-1)+delta));
					}else{
						if(exp+delta>=stopTier){break;}
						final int shiftAmt=Math.min(delta, tailWidth);
						int newTail=0;
						for(int ti=0; ti<numTails; ti++){
							int t2=((tail>>>(ti*tailWidth))&tailMask)>>>shiftAmt;
							if(ti==tailIdx2){t2|=(1<<(tailWidth-1));}
							newTail|=(t2<<(ti*tailWidth));
						}
						tail=newTail;
						exp+=delta;
						final int target=exp-margin;
						if(target<=0){eeMask=-1L;}
						else{
							final int minRaw=(TIER_DENOM*target)/TIER_NUMER;
							eeMask=(minRaw>=64)?0:(-1L)>>>minRaw;
						}
					}
				}
			}

			if(tc>=nextCheck){
				if(exp<stopTier){
					final double invTc=1.0/tc;
					for(int ti=0; ti<tables.length; ti++){
						final double est=tables[ti][exp][tail];
						if(est>0){
							final double err=(est-tc)*invTc;
							result[ti][exp][0]+=err;
							result[ti][exp][1]+=Math.abs(err);
							result[ti][exp][2]++;
						}
					}
				}
				nextCheck=Math.max(tc+1, nextCheck*1.01);
			}
		}
	}

	/** SIMD version of testOneBucketMulti: 4 RNG streams, one bucket, checkpoint-aware. */
	void testOneBucketMultiSIMD(long seed, double[][][] tables, double[][][] result){
		long[] ss0=new long[4], ss1=new long[4], ss2=new long[4], ss3=new long[4];
		for(int i=0; i<4; i++){
			long x=seed*4+i;
			x=splitmix(x); ss0[i]=x;
			x=splitmix(x); ss1[i]=x;
			x=splitmix(x); ss2[i]=x;
			x=splitmix(x); ss3[i]=x;
		}
		LongVector vs0=LongVector.fromArray(SPECIES, ss0, 0);
		LongVector vs1=LongVector.fromArray(SPECIES, ss1, 0);
		LongVector vs2=LongVector.fromArray(SPECIES, ss2, 0);
		LongVector vs3=LongVector.fromArray(SPECIES, ss3, 0);
		for(int w=0; w<4; w++){
			LongVector t=vs1.lanewise(VectorOperators.LSHL, 17);
			vs2=vs2.lanewise(VectorOperators.XOR, vs0);
			vs3=vs3.lanewise(VectorOperators.XOR, vs1);
			vs1=vs1.lanewise(VectorOperators.XOR, vs2);
			vs0=vs0.lanewise(VectorOperators.XOR, vs3);
			vs2=vs2.lanewise(VectorOperators.XOR, t);
			vs3=vs3.lanewise(VectorOperators.ROL, 45);
		}

		int exp=0, tail=0;
		long eeMask=-1L;
		long tc=0;
		double nextCheck=1.0;
		LongVector mFlip=LongVector.broadcast(SPECIES, eeMask^Long.MIN_VALUE);

		for(;;){
			// How many SIMD iterations until next checkpoint?
			final long maxIter=Math.max(1, ((long)nextCheck-tc+3)/4);

			// SIMD INNER LOOP: generate 4 hashes, check eeMask, cap at checkpoint
			LongVector result_vec;
			long inner=0;
			do{
				result_vec=vs0.add(vs3);
				LongVector t=vs1.lanewise(VectorOperators.LSHL, 17);
				vs2=vs2.lanewise(VectorOperators.XOR, vs0);
				vs3=vs3.lanewise(VectorOperators.XOR, vs1);
				vs1=vs1.lanewise(VectorOperators.XOR, vs2);
				vs0=vs0.lanewise(VectorOperators.XOR, vs3);
				vs2=vs2.lanewise(VectorOperators.XOR, t);
				vs3=vs3.lanewise(VectorOperators.ROL, 45);
				inner++;
			}while(result_vec.lanewise(VectorOperators.XOR, FLIP_BIT)
				.compare(VectorOperators.GT, mFlip).allTrue() && inner<maxIter);

			// Account for (inner-1) full-skip iterations
			tc+=(inner-1)*4;

			// Process 4 hashes from last iteration
			final long h0=result_vec.lane(0), h1=result_vec.lane(1), h2=result_vec.lane(2), h3=result_vec.lane(3);
			final long[] hh={h0, h1, h2, h3};
			for(int i=0; i<4; i++){
				tc++;
				final long hash=hh[i];
				if(Long.compareUnsigned(hash, eeMask)<=0){
					final int rawNLZ=Long.numberOfLeadingZeros(hash);
					final int hashNLZ=(TIER_NUMER!=1)?(TIER_NUMER*rawNLZ+((((rawNLZ>=43)?(int)((hash<<(rawNLZ-42))&0xFFFFF):(int)((hash>>>(42-rawNLZ))&0xFFFFF))>=MANTISSA_THRESHOLD)?1:0))/TIER_DENOM:rawNLZ;
					final int delta=hashNLZ-exp;
					if(delta>=-margin){
						final int tailIdx=(int)((hash>>>1)&(numTails-1));
						if(delta<=0){
							tail|=(1<<(tailIdx*tailWidth+(tailWidth-1)+delta));
						}else{
							if(exp+delta>=stopTier){return;}
							final int shiftAmt=Math.min(delta, tailWidth);
							int newTail=0;
							for(int ti=0; ti<numTails; ti++){
								int t2=((tail>>>(ti*tailWidth))&tailMask)>>>shiftAmt;
								if(ti==tailIdx){t2|=(1<<(tailWidth-1));}
								newTail|=(t2<<(ti*tailWidth));
							}
							tail=newTail; exp+=delta;
							final int target=exp-margin;
							if(target<=0){eeMask=-1L;}
							else{
								final int minRaw=(TIER_DENOM*target)/TIER_NUMER;
								eeMask=(minRaw>=64)?0:(-1L)>>>minRaw;
							}
							mFlip=LongVector.broadcast(SPECIES, eeMask^Long.MIN_VALUE);
						}
					}
				}
				if(tc>=nextCheck){
					if(exp<stopTier){
						final double invTc=1.0/tc;
						for(int ti=0; ti<tables.length; ti++){
							final double est=tables[ti][exp][tail];
							if(est>0){
								final double err=(est-tc)*invTc;
								result[ti][exp][0]+=err;
								result[ti][exp][1]+=Math.abs(err);
								result[ti][exp][2]++;
							}
						}
					}
					nextCheck=Math.max(tc+1, nextCheck*1.01);
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Correction            ----------------*/
	/*--------------------------------------------------------------*/

	double[][] applyCorrection(double[][] table, double[] cf){
		final double[][] out=new double[table.length][numStates];
		for(int tier=0; tier<table.length; tier++)
			for(int s=0; s<numStates; s++)
				out[tier][s]=table[tier][s]*cf[tier];
		return out;
	}

	double[][][] selfTestPerState(double[][] table, int totalRuns, long seedOffset){
		final double[][][] result=new double[stopTier][numStates][3];
		final AtomicLong counter=new AtomicLong(0);
		final double[][][][] tResult=new double[threads][stopTier][numStates][3];

		final Thread[] arr=new Thread[threads];
		for(int t=0; t<threads; t++){
			final int tid=t;
			arr[t]=new Thread(()->{
				final double[][][] lr=tResult[tid];
				long runId;
				while((runId=counter.getAndIncrement())<totalRuns){
					testOneBucketPerState(seedOffset+runId, table, lr);
				}
			});
			arr[t].start();
		}
		try{for(final Thread th:arr){th.join();}}catch(InterruptedException e){throw new RuntimeException(e);}

		for(int t=0; t<threads; t++)
			for(int tier=0; tier<stopTier; tier++)
				for(int s=0; s<numStates; s++)
					for(int i=0; i<3; i++)
						result[tier][s][i]+=tResult[t][tier][s][i];
		return result;
	}

	void testOneBucketPerState(long seed, double[][] table, double[][][] result){
		final FastRandomXoshiro rng=new FastRandomXoshiro(seed);
		int exp=0, tail=0;
		long eeMask=-1L;
		double nextCheck=1.0;

		for(long tc=1; ; tc++){
			final long hash=rng.nextLong();
			if(Long.compareUnsigned(hash, eeMask)<=0){
				final int rawNLZ2=Long.numberOfLeadingZeros(hash);
				final int hashNLZ2;
				if(TIER_NUMER!=1){
					final int mB;
					if(rawNLZ2>=43){mB=(int)((hash<<(rawNLZ2-42))&0xFFFFF);}
					else{mB=(int)((hash>>>(42-rawNLZ2))&0xFFFFF);}
					hashNLZ2=(TIER_NUMER*rawNLZ2+((mB>=MANTISSA_THRESHOLD)?1:0))/TIER_DENOM;
				}else{hashNLZ2=rawNLZ2;}
				final int delta=hashNLZ2-exp;
				final int tailIdx2=(int)((hash>>>1)&(numTails-1));
				if(delta>=-margin){
					if(delta<=0){
						tail|=(1<<(tailIdx2*tailWidth+(tailWidth-1)+delta));
					}else{
						if(exp+delta>=stopTier){break;}
						final int shiftAmt=Math.min(delta, tailWidth);
						int newTail=0;
						for(int ti=0; ti<numTails; ti++){
							int t2=((tail>>>(ti*tailWidth))&tailMask)>>>shiftAmt;
							if(ti==tailIdx2){t2|=(1<<(tailWidth-1));}
							newTail|=(t2<<(ti*tailWidth));
						}
						tail=newTail;
						exp+=delta;
						{
							final int target=exp-margin;
							if(target<=0){eeMask=-1L;}
							else{
								final int minRaw=(TIER_DENOM*target)/TIER_NUMER;
								eeMask=(minRaw>=64)?0:(-1L)>>>minRaw;
							}
						}
					}
				}
			}

			if(tc>=nextCheck){
				if(exp<stopTier){
					final double est=table[exp][tail];
					if(est>0){
						if(PS_TARGET==0){
							final double err=(est-tc)/(double)tc;
							result[exp][tail][0]+=err;
						}else if(PS_TARGET==1){
							result[exp][tail][0]+=tc;
						}else{
							result[exp][tail][0]+=Math.log(tc);
						}
						result[exp][tail][1]+=Math.abs((est-tc)/(double)tc);
						result[exp][tail][2]++;
					}
				}
				nextCheck=Math.max(tc+1, nextCheck*1.01);
			}
		}
	}

	double[][] applyPerStateCorrection(double[][] table, double[][] cf){
		final double[][] out=new double[table.length][numStates];
		for(int tier=0; tier<table.length; tier++)
			for(int s=0; s<numStates; s++)
				out[tier][s]=table[tier][s]*cf[tier][s];
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------            I/O                ----------------*/
	/*--------------------------------------------------------------*/

	void writeTable(double[][] table, String name, PrintStream out){
		writeTable(table, name, out, numStates);
	}

	void writeTable(double[][] table, String name, PrintStream out, int cols){
		out.printf("#%s\t%d%n", name, table.length);
		StringBuilder hdr=new StringBuilder("#Tier");
		for(int s=0; s<cols; s++) hdr.append('\t').append(s);
		out.println(hdr);
		for(int tier=0; tier<table.length; tier++){
			StringBuilder row=new StringBuilder(Integer.toString(tier));
			for(int s=0; s<cols; s++) row.append('\t').append(String.format("%.6f", table[tier][s]));
			out.println(row);
		}
		out.println("#");
	}

	static long hash64shift(long key){
		key=(~key)+(key<<21); key^=(key>>>24);
		key+=(key<<3)+(key<<8); key^=(key>>>14);
		key+=(key<<2)+(key<<4); key^=(key>>>28);
		key+=(key<<31); return key;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Main               ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws Exception {
		{PreParser pp=new PreParser(args, null, false); args=pp.args;}

		int tailBits=4, stopTier=16, runs=1000, threads=4;
		long masterSeed=0; // 0=deterministic (seeds 0..runs-1), -1=random, >0=offset
		boolean variance=false, rawOutput=false, perStateCorrection=false, allowNullOut=false;
		int smFilter=-1; // -1=all, 0=entry, 1=all, 2=entryexit, 3=pctile
		String outFile=null;

		for(String arg : args){
			String[] kv=arg.split("=",2);
			if(kv.length!=2) throw new RuntimeException("Unknown: "+arg);
			switch(kv[0]){
				case "tailBits": case "tail": tailBits=Integer.parseInt(kv[1]); break;
				case "numTails": case "tails": NUM_TAILS=Integer.parseInt(kv[1]); break;
				case "stopTier": case "stoptier": case "tiers": stopTier=Integer.parseInt(kv[1]); break;
				case "runs": runs=Parse.parseIntKMG(kv[1]); break;
				case "threads": case "t": threads=Integer.parseInt(kv[1]); break;
				case "seed": masterSeed=Long.parseLong(kv[1]); break;
				case "out": outFile=kv[1]; break;
				case "compressed": COMPRESSED=Parse.parseBoolean(kv[1]); break;
				case "expanded": EXPANDED=Parse.parseBoolean(kv[1]); break;
				case "simd": USE_SIMD=Parse.parseBoolean(kv[1]); break;
				case "variance": case "var": variance=Parse.parseBoolean(kv[1]); break;
				case "raw": rawOutput=Parse.parseBoolean(kv[1]); break;
				case "perstate": perStateCorrection=Parse.parseBoolean(kv[1]); break;
				case "allownullout": allowNullOut=Parse.parseBoolean(kv[1]); break;
				case "sm": case "sample":
					switch(kv[1].toLowerCase()){
						case "0": case "entry": smFilter=0; break;
						case "1": case "all": smFilter=1; break;
						case "2": case "entryexit": smFilter=2; break;
						case "3": case "pctile": smFilter=3; break;
						default: smFilter=Integer.parseInt(kv[1]); break;
					} break;
				case "target":
					if(kv[1].startsWith("harm") || kv[1].equals("0")){PS_TARGET=0;}
					else if(kv[1].startsWith("arith") || kv[1].equals("1")){PS_TARGET=1;}
					else if(kv[1].startsWith("geo") || kv[1].equals("2")){PS_TARGET=2;}
					else{throw new RuntimeException("Unknown target: "+kv[1]);}
					break;
				default: throw new RuntimeException("Unknown: "+arg);
			}
		}
		if(COMPRESSED){TIER_NUMER=2; TIER_DENOM=3;}
		else if(EXPANDED){TIER_NUMER=2; TIER_DENOM=1;}
		else{TIER_NUMER=1; TIER_DENOM=1;}
		if(tailBits%NUM_TAILS!=0) throw new RuntimeException("tailBits must be divisible by numTails");
		if(outFile==null && !allowNullOut) throw new RuntimeException("out= parameter required (use allownullout=t for benchmark-only mode)");

		final long seedOffset;
		if(masterSeed==-1){seedOffset=System.nanoTime();}
		else{seedOffset=masterSeed;}

		System.err.printf("PerTierHistorySimulator: tailBits=%d numTails=%d stopTier=%d runs=%d threads=%d tierNumer=%d tierDenom=%d seed=%d%n",
			tailBits, NUM_TAILS, stopTier, runs, threads, TIER_NUMER, TIER_DENOM, seedOffset);
		System.err.printf("  %d runs for training (seeds %d..%d), %d runs for testing (seeds %d..%d)%n",
			runs, seedOffset, seedOffset+runs-1, runs, seedOffset+runs, seedOffset+2*runs-1);

		PerTierHistorySimulator sim=new PerTierHistorySimulator(tailBits, NUM_TAILS, stopTier, threads);

		// Phase 1: Train
		long t0=System.currentTimeMillis();
		sim.train(runs, seedOffset, smFilter);
		System.err.printf("Training: %d ms%n", System.currentTimeMillis()-t0);
		System.err.println();

		final int summSm=(smFilter>=0 ? smFilter : S_ALL);
		System.err.printf("%-6s %-14s %-12s%n", "Tier", "LinAvg("+SNAME[summSm]+")", "Obs");
		for(int tier=0; tier<stopTier; tier++){
			long obs=0; double s=0;
			for(int st=0; st<sim.numStates; st++){obs+=sim.counts[summSm][tier][st]; s+=sim.sums[summSm][tier][st];}
			System.err.printf("%-6d %-14.2f %-12d%n", tier, (obs>0?s/obs:0), obs);
		}

		// Per-tier bit occupancy
		System.err.println();
		{
			final int tw=tailBits/NUM_TAILS;
			StringBuilder hdr=new StringBuilder(String.format("%-6s %10s", "Tier", "Obs"));
			for(int b=0; b<tailBits; b++){
				int tailNum=b/tw; int depth=tw-1-(b%tw);
				hdr.append(String.format(" %8s", "t"+tailNum+"d"+depth));
			}
			System.err.println(hdr);
			for(int tier=0; tier<stopTier; tier++){
				long total=0;
				long[] bitSet=new long[tailBits];
				for(int s=0; s<sim.numStates; s++){
					long n=sim.counts[S_ALL][tier][s];
					total+=n;
					for(int b=0; b<tailBits; b++){
						if((s&(1<<b))!=0){bitSet[b]+=n;}
					}
				}
				if(total==0) continue;
				StringBuilder row=new StringBuilder(String.format("%-6d %10d", tier, total));
				for(int b=0; b<tailBits; b++){
					row.append(String.format(" %8.4f", (double)bitSet[b]/total));
				}
				System.err.println(row);
			}
		}
		System.err.println();

		// Build all 12 tables
		double[][][][] rawTables=new double[NS][NA][][];
		for(int sm=0; sm<NS; sm++)
			for(int am=0; am<NA; am++)
				rawTables[sm][am]=sim.buildTable(sm, am);

		// Phase 2: Test raw, measure per-tier bias, apply half+full correction, report
		System.err.printf("%nTesting (runs=%d, 1%%-spaced checkpoints, raw → half-corrected → full-corrected):%n", runs);
		System.err.printf("%-12s %-6s │%10s %10s│%10s %10s│%10s %10s%n",
			"Sample", "Avg", "RawSgn%", "RawAbs%", "HalfSgn%", "HalfAbs%", "FullSgn%", "FullAbs%");
		System.err.println("─────────────────────┼─────────────────────┼─────────────────────┼─────────────────────");

		t0=System.currentTimeMillis();
		final long testSeedOffset=seedOffset+runs;
		final int testSmStart=(smFilter>=0 ? smFilter : 0);
		final int testSmEnd=(smFilter>=0 ? smFilter+1 : NS);
		for(int sm=testSmStart; sm<testSmEnd; sm++){
			// Batch all averaging modes: test raw tables in one pass
			double[][][] rawBatch=new double[NA][][];
			for(int am=0; am<NA; am++) rawBatch[am]=rawTables[sm][am];
			double[][][] rawResults=sim.selfTestMulti(rawBatch, runs, testSeedOffset);

			// Compute corrections, then batch-test corrected tables
			double[][][] corrBatch=new double[NA*2][][];
			for(int am=0; am<NA; am++){
				double[] cfFull=new double[stopTier];
				double[] cfHalf=new double[stopTier];
				for(int tier=0; tier<stopTier; tier++){
					if(rawResults[am][tier][2]>0){
						double bias=rawResults[am][tier][0]/rawResults[am][tier][2];
						cfFull[tier]=(1.0+bias)!=0 ? 1.0/(1.0+bias) : 1.0;
						cfHalf[tier]=(1.0+0.5*bias)!=0 ? 1.0/(1.0+0.5*bias) : 1.0;
					}else{ cfFull[tier]=1.0; cfHalf[tier]=1.0; }
				}
				corrBatch[am*2]=sim.applyCorrection(rawBatch[am], cfHalf);
				corrBatch[am*2+1]=sim.applyCorrection(rawBatch[am], cfFull);
			}
			double[][][] corrResults=sim.selfTestMulti(corrBatch, runs, testSeedOffset);

			for(int am=0; am<NA; am++){
				double sR=0,aR=0,sH=0,aH=0,sF=0,aF=0; long nR=0,nH=0,nF=0;
				for(int tier=0; tier<stopTier; tier++){
					sR+=rawResults[am][tier][0]; aR+=rawResults[am][tier][1]; nR+=(long)rawResults[am][tier][2];
					sH+=corrResults[am*2][tier][0]; aH+=corrResults[am*2][tier][1]; nH+=(long)corrResults[am*2][tier][2];
					sF+=corrResults[am*2+1][tier][0]; aF+=corrResults[am*2+1][tier][1]; nF+=(long)corrResults[am*2+1][tier][2];
				}
				System.err.printf("%-12s %-6s │%+10.2f %10.2f│%+10.2f %10.2f│%+10.2f %10.2f%n",
					SNAME[sm], ANAME[am],
					nR>0?100*sR/nR:0, nR>0?100*aR/nR:0,
					nH>0?100*sH/nH:0, nH>0?100*aH/nH:0,
					nF>0?100*sF/nF:0, nF>0?100*aF/nF:0);
			}
		}
		System.err.printf("Testing: %d ms%n", System.currentTimeMillis()-t0);

		if(outFile!=null){
			PrintStream out=new PrintStream(outFile);
			if(rawOutput){
				for(int sm=0; sm<NS; sm++)
					for(int am=0; am<NA; am++)
						sim.writeTable(rawTables[sm][am], SNAME[sm]+"_"+ANAME[am], out);
				if(variance){
					for(int sm=0; sm<NS; sm++)
						sim.writeTable(sim.buildVarianceTable(sm), "var_"+SNAME[sm], out, sim.numCondensed);
				}
				out.close();
				System.err.printf("Raw tables written to %s%s%n", outFile, variance?" (with variance)":"");
			}else if(perStateCorrection){
				final int smStart=(smFilter>=0 ? smFilter : 0);
				final int smEnd=(smFilter>=0 ? smFilter+1 : NS);
				for(int sm=smStart; sm<smEnd; sm++){
					for(int am=0; am<NA; am++){
						double[][] raw=rawTables[sm][am];
						double[][][] rPS=sim.selfTestPerState(raw, runs, testSeedOffset);
						double[][] cfPS=new double[stopTier][sim.numStates];
						for(int tier=0; tier<stopTier; tier++){
							for(int s=0; s<sim.numStates; s++){
								if(rPS[tier][s][2]>0 && raw[tier][s]>0){
									if(PS_TARGET==0){
										double bias=rPS[tier][s][0]/rPS[tier][s][2];
										cfPS[tier][s]=(1.0+bias)!=0 ? 1.0/(1.0+bias) : 1.0;
									}else if(PS_TARGET==1){
										double meanTrue=rPS[tier][s][0]/rPS[tier][s][2];
										cfPS[tier][s]=meanTrue/raw[tier][s];
									}else{
										double meanLogTrue=rPS[tier][s][0]/rPS[tier][s][2];
										cfPS[tier][s]=Math.exp(meanLogTrue)/raw[tier][s];
									}
								}else{cfPS[tier][s]=1.0;}
							}
						}
						double[][] corrected=sim.applyPerStateCorrection(raw, cfPS);
						sim.writeTable(corrected, SNAME[sm]+"_"+ANAME[am], out);
						if(variance && am==0){
							long t1=System.currentTimeMillis();
							double[][][] rPS2=sim.selfTestPerState(corrected, runs, testSeedOffset+runs);
							double[][] varPS=new double[stopTier][sim.numStates];
							for(int tier=0; tier<stopTier; tier++){
								for(int s=0; s<sim.numStates; s++){
									if(rPS2[tier][s][2]>0){
										varPS[tier][s]=rPS2[tier][s][1]/rPS2[tier][s][2];
									}
								}
							}
							sim.writeTable(varPS, "var_"+SNAME[sm], out, sim.numStates);
							System.err.printf("  Phase C variance for %s_%s: %d ms%n", SNAME[sm], ANAME[am], System.currentTimeMillis()-t1);
						}
					}
				}
				out.close();
				System.err.printf("Per-state corrected tables written to %s%s%n", outFile, variance?" (with per-state variance)":"");
			}else{
				final int outSmStart=(smFilter>=0 ? smFilter : 0);
				final int outSmEnd=(smFilter>=0 ? smFilter+1 : NS);
				for(int sm=outSmStart; sm<outSmEnd; sm++){
					double[][][] rawBatch=new double[NA][][];
					for(int am=0; am<NA; am++) rawBatch[am]=rawTables[sm][am];
					double[][][] results=sim.selfTestMulti(rawBatch, runs, testSeedOffset);
					for(int am=0; am<NA; am++){
						double[] cf=new double[stopTier];
						for(int tier=0; tier<stopTier; tier++){
							if(results[am][tier][2]>0){
								double bias=results[am][tier][0]/results[am][tier][2];
								cf[tier]=(1.0+bias)!=0 ? 1.0/(1.0+bias) : 1.0;
							}else{cf[tier]=1.0;}
						}
						sim.writeTable(sim.applyCorrection(rawBatch[am], cf), SNAME[sm]+"_"+ANAME[am], out);
					}
				}
				if(variance){
					for(int sm=outSmStart; sm<outSmEnd; sm++)
						sim.writeTable(sim.buildVarianceTable(sm), "var_"+SNAME[sm], out, sim.numCondensed);
				}
				out.close();
				System.err.printf("Corrected tables written to %s%s%n", outFile, variance?" (with variance)":"");
			}
		}
	}
}
