package cardinality;

import java.io.PrintStream;
import java.util.concurrent.atomic.AtomicLong;
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
	static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*1048576);

	final int tailBits, halfWidth, halfMask, numStates, margin, stopTier, threads;

	long[][][] counts;
	double[][][] sums, sumsLog, sumsInv;

	PerTierHistorySimulator(int tailBits_, int stopTier_, int threads_){
		tailBits=tailBits_; halfWidth=tailBits/2; halfMask=(1<<halfWidth)-1;
		numStates=1<<tailBits; margin=halfWidth-1;
		stopTier=stopTier_; threads=threads_;
		counts=new long[NS][stopTier][numStates];
		sums=new double[NS][stopTier][numStates];
		sumsLog=new double[NS][stopTier][numStates];
		sumsInv=new double[NS][stopTier][numStates];
	}

	/*--------------------------------------------------------------*/
	/*----------------     Single-Bucket Simulation  ----------------*/
	/*--------------------------------------------------------------*/

	/** Simulates one bucket from trueCard=1 until exp reaches stopTier.
	 *  Records into the provided accumulators under all 4 sampling modes. */
	void simulateOneBucket(long seed,
			long[][][] lc, double[][][] ls, double[][][] lsl, double[][][] lsi){
		final FastRandomXoshiro rng=new FastRandomXoshiro(seed);
		int exp=0, tail=0;
		long eeMask=-1L;
		int prevTier=-1, prevState=-1;
		double nextPctile=1.0;

		for(long tc=1; ; tc++){
			final long hash=hash64shift(rng.nextLong());

			if(Long.compareUnsigned(hash, eeMask)>0){
				record(lc[S_ALL], ls[S_ALL], lsl[S_ALL], lsi[S_ALL], exp, tail, tc);
				if(tc>=nextPctile){
					record(lc[S_PCTILE], ls[S_PCTILE], lsl[S_PCTILE], lsi[S_PCTILE], exp, tail, tc);
					nextPctile=Math.max(tc+1, nextPctile*1.01);
				}
				continue;
			}

			final int rawNLZ=Long.numberOfLeadingZeros(hash);
			final int hashNLZ;
			if(COMPRESSED){
				final int mBits;
				if(rawNLZ>=43){mBits=(int)((hash<<(rawNLZ-42))&0xFFFFF);}
				else{mBits=(int)((hash>>>(42-rawNLZ))&0xFFFFF);}
				final int mantissa=(mBits>=MANTISSA_THRESHOLD)?1:0;
				hashNLZ=(2*rawNLZ+mantissa)/3;
			}else{
				hashNLZ=rawNLZ;
			}
			final int delta=hashNLZ-exp;
			final int histBit=(int)((hash>>>1)&1);

			if(delta<-margin){
				record(lc[S_ALL], ls[S_ALL], lsl[S_ALL], lsi[S_ALL], exp, tail, tc);
				if(tc>=nextPctile){
					record(lc[S_PCTILE], ls[S_PCTILE], lsl[S_PCTILE], lsi[S_PCTILE], exp, tail, tc);
					nextPctile=Math.max(tc+1, nextPctile*1.01);
				}
				continue;
			}

			if(delta<=0){
				tail|=(1<<((halfWidth-1)+delta+(histBit*halfWidth)));
			}else{
				if(exp+delta>=stopTier){break;}
				final int shiftAmt=Math.min(delta, halfWidth);
				final int h0=(tail&halfMask)>>>shiftAmt;
				final int h1=((tail>>>halfWidth)&halfMask)>>>shiftAmt;
				tail=(h1<<halfWidth)|h0;
				tail|=(1<<((histBit==0)?(halfWidth-1):(tailBits-1)));
				exp+=delta;
				if(COMPRESSED){
					final int target=exp-margin;
					if(target<=0){eeMask=-1L;}
					else{final int minRaw=(3*target)/2; eeMask=(minRaw>=64)?0:(-1L)>>>minRaw;}
				}else{
					eeMask=(exp<=margin)?-1L:(-1L)>>>(exp-margin);
				}
			}

			final boolean changed=(exp!=prevTier || tail!=prevState);

			record(lc[S_ALL], ls[S_ALL], lsl[S_ALL], lsi[S_ALL], exp, tail, tc);

			if(changed){
				record(lc[S_ENTRY], ls[S_ENTRY], lsl[S_ENTRY], lsi[S_ENTRY], exp, tail, tc);
			}

			if(changed && prevTier>=0){
				record(lc[S_ENTRYEXIT], ls[S_ENTRYEXIT], lsl[S_ENTRYEXIT], lsi[S_ENTRYEXIT],
					prevTier, prevState, tc);
			}
			if(changed){
				record(lc[S_ENTRYEXIT], ls[S_ENTRYEXIT], lsl[S_ENTRYEXIT], lsi[S_ENTRYEXIT],
					exp, tail, tc);
			}

			if(tc>=nextPctile){
				record(lc[S_PCTILE], ls[S_PCTILE], lsl[S_PCTILE], lsi[S_PCTILE], exp, tail, tc);
				nextPctile=Math.max(tc+1, nextPctile*1.01);
			}

			prevTier=exp; prevState=tail;
		}
	}

	static void record(long[][] c, double[][] s, double[][] sl, double[][] si,
			int tier, int state, long tc){
		c[tier][state]++;
		s[tier][state]+=tc;
		sl[tier][state]+=Math.log(tc);
		si[tier][state]+=1.0/tc;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Training            ----------------*/
	/*--------------------------------------------------------------*/

	void train(int totalRuns, long seedOffset) throws InterruptedException {
		final AtomicLong counter=new AtomicLong(0);
		final long[][][][] tC=new long[threads][NS][stopTier][numStates];
		final double[][][][] tS=new double[threads][NS][stopTier][numStates];
		final double[][][][] tSL=new double[threads][NS][stopTier][numStates];
		final double[][][][] tSI=new double[threads][NS][stopTier][numStates];

		final Thread[] arr=new Thread[threads];
		for(int t=0; t<threads; t++){
			final int tid=t;
			arr[t]=new Thread(()->{
				long runId;
				while((runId=counter.getAndIncrement())<totalRuns){
					simulateOneBucket(seedOffset+runId, tC[tid], tS[tid], tSL[tid], tSI[tid]);
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
					}
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
			final long hash=hash64shift(rng.nextLong());
			if(Long.compareUnsigned(hash, eeMask)<=0){
				final int rawNLZ2=Long.numberOfLeadingZeros(hash);
				final int hashNLZ2;
				if(COMPRESSED){
					final int mB;
					if(rawNLZ2>=43){mB=(int)((hash<<(rawNLZ2-42))&0xFFFFF);}
					else{mB=(int)((hash>>>(42-rawNLZ2))&0xFFFFF);}
					hashNLZ2=(2*rawNLZ2+((mB>=MANTISSA_THRESHOLD)?1:0))/3;
				}else{hashNLZ2=rawNLZ2;}
				final int delta=hashNLZ2-exp;
				final int histBit=(int)((hash>>>1)&1);
				if(delta>=-margin){
					if(delta<=0){
						tail|=(1<<((halfWidth-1)+delta+(histBit*halfWidth)));
					}else{
						if(exp+delta>=stopTier){break;}
						final int shiftAmt=Math.min(delta, halfWidth);
						final int h0=(tail&halfMask)>>>shiftAmt;
						final int h1=((tail>>>halfWidth)&halfMask)>>>shiftAmt;
						tail=(h1<<halfWidth)|h0;
						tail|=(1<<((histBit==0)?(halfWidth-1):(tailBits-1)));
						exp+=delta;
						if(COMPRESSED){
							final int target=exp-margin;
							if(target<=0){eeMask=-1L;}
							else{final int minRaw=(3*target)/2; eeMask=(minRaw>=64)?0:(-1L)>>>minRaw;}
						}else{
							eeMask=(exp<=margin)?-1L:(-1L)>>>(exp-margin);
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

	/*--------------------------------------------------------------*/
	/*----------------            I/O                ----------------*/
	/*--------------------------------------------------------------*/

	void writeTable(double[][] table, String name, PrintStream out){
		out.printf("#%s\t%d%n", name, table.length);
		StringBuilder hdr=new StringBuilder("#Tier");
		for(int s=0; s<numStates; s++) hdr.append('\t').append(s);
		out.println(hdr);
		for(int tier=0; tier<table.length; tier++){
			StringBuilder row=new StringBuilder(Integer.toString(tier));
			for(int s=0; s<numStates; s++) row.append('\t').append(String.format("%.6f", table[tier][s]));
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
		String outFile=null;

		for(String arg : args){
			String[] kv=arg.split("=",2);
			if(kv.length!=2) throw new RuntimeException("Unknown: "+arg);
			switch(kv[0]){
				case "tailBits": case "tail": tailBits=Integer.parseInt(kv[1]); break;
				case "stopTier": case "stoptier": case "tiers": stopTier=Integer.parseInt(kv[1]); break;
				case "runs": runs=Parse.parseIntKMG(kv[1]); break;
				case "threads": case "t": threads=Integer.parseInt(kv[1]); break;
				case "out": outFile=kv[1]; break;
				case "compressed": COMPRESSED=Parse.parseBoolean(kv[1]); break;
				default: throw new RuntimeException("Unknown: "+arg);
			}
		}
		if(tailBits%2!=0) throw new RuntimeException("tailBits must be even");

		System.err.printf("PerTierHistorySimulator: tailBits=%d stopTier=%d runs=%d threads=%d compressed=%b%n",
			tailBits, stopTier, runs, threads, COMPRESSED);
		System.err.printf("  %d runs for training (seeds 0..%d), %d runs for testing (seeds %d..%d)%n",
			runs, runs-1, runs, runs, 2*runs-1);

		PerTierHistorySimulator sim=new PerTierHistorySimulator(tailBits, stopTier, threads);

		// Phase 1: Train
		long t0=System.currentTimeMillis();
		sim.train(runs, 0);
		System.err.printf("Training: %d ms%n%n", System.currentTimeMillis()-t0);

		System.err.printf("%-6s %-14s %-12s%n", "Tier", "LinAvg(all)", "Obs(all)");
		for(int tier=0; tier<stopTier; tier++){
			long obs=0; double s=0;
			for(int st=0; st<sim.numStates; st++){obs+=sim.counts[S_ALL][tier][st]; s+=sim.sums[S_ALL][tier][st];}
			System.err.printf("%-6d %-14.2f %-12d%n", tier, (obs>0?s/obs:0), obs);
		}

		// Per-tier bit occupancy
		System.err.println();
		{
			final int hw=tailBits/2;
			StringBuilder hdr=new StringBuilder(String.format("%-6s %10s", "Tier", "Obs"));
			for(int b=0; b<tailBits; b++){
				int tailNum=b/hw; int depth=hw-1-(b%hw);
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
		final long testSeedOffset=runs;
		for(int sm=0; sm<NS; sm++){
			for(int am=0; am<NA; am++){
				double[][] raw=rawTables[sm][am];

				// Test raw
				double[][] rRaw=sim.selfTest(raw, runs, testSeedOffset);

				// Compute per-tier correction from raw results
				double[] cfFull=new double[stopTier];
				double[] cfHalf=new double[stopTier];
				for(int tier=0; tier<stopTier; tier++){
					if(rRaw[tier][2]>0){
						double bias=rRaw[tier][0]/rRaw[tier][2];
						cfFull[tier]=(1.0+bias)!=0 ? 1.0/(1.0+bias) : 1.0;
						cfHalf[tier]=(1.0+0.5*bias)!=0 ? 1.0/(1.0+0.5*bias) : 1.0;
					}else{
						cfFull[tier]=1.0; cfHalf[tier]=1.0;
					}
				}

				// Test corrected (same seed — measuring correction effect on same data)
				double[][] halfCorr=sim.applyCorrection(raw, cfHalf);
				double[][] rHalf=sim.selfTest(halfCorr, runs, testSeedOffset);
				double[][] fullCorr=sim.applyCorrection(raw, cfFull);
				double[][] rFull=sim.selfTest(fullCorr, runs, testSeedOffset);

				double sR=0,aR=0,sH=0,aH=0,sF=0,aF=0; long nR=0,nH=0,nF=0;
				for(int tier=0; tier<stopTier; tier++){
					sR+=rRaw[tier][0]; aR+=rRaw[tier][1]; nR+=(long)rRaw[tier][2];
					sH+=rHalf[tier][0]; aH+=rHalf[tier][1]; nH+=(long)rHalf[tier][2];
					sF+=rFull[tier][0]; aF+=rFull[tier][1]; nF+=(long)rFull[tier][2];
				}

				System.err.printf("%-12s %-6s │%+10.2f %10.2f│%+10.2f %10.2f│%+10.2f %10.2f%n",
					SNAME[sm], ANAME[am],
					nR>0?100*sR/nR:0, nR>0?100*aR/nR:0,
					nH>0?100*sH/nH:0, nH>0?100*aH/nH:0,
					nF>0?100*sF/nF:0, nF>0?100*aF/nF:0);
			}
		}
		System.err.printf("Testing: %d ms%n", System.currentTimeMillis()-t0);

		// Also output corrected tables
		if(outFile!=null){
			PrintStream out=new PrintStream(outFile);
			for(int sm=0; sm<NS; sm++){
				for(int am=0; am<NA; am++){
					double[][] raw=rawTables[sm][am];
					double[][] rRaw=sim.selfTest(raw, runs, testSeedOffset);
					double[] cf=new double[stopTier];
					for(int tier=0; tier<stopTier; tier++){
						if(rRaw[tier][2]>0){
							double bias=rRaw[tier][0]/rRaw[tier][2];
							cf[tier]=(1.0+bias)!=0 ? 1.0/(1.0+bias) : 1.0;
						}else{cf[tier]=1.0;}
					}
					sim.writeTable(sim.applyCorrection(raw, cf), SNAME[sm]+"_"+ANAME[am], out);
				}
			}
			out.close();
			System.err.printf("Corrected tables written to %s%n", outFile);
		}
		if(false){// disabled raw output
			PrintStream out=new PrintStream(outFile);
			for(int sm=0; sm<NS; sm++)
				for(int am=0; am<NA; am++)
					sim.writeTable(rawTables[sm][am], SNAME[sm]+"_"+ANAME[am], out);
			out.close();
			System.err.printf("Tables written to %s%n", outFile);
		}
	}
}
