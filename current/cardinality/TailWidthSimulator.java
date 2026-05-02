package cardinality;

import rand.FastRandomXoshiro;

/**
 * Simulates a single TTLL register with variable tail width to measure
 * per-bit set probabilities at each tier.
 *
 * With N-bit tails, each tail tracks:
 *   bit N-1 = "saw current tier"    (E-0)
 *   bit N-2 = "saw tier-1"          (E-1)
 *   ...
 *   bit 0   = "saw tier-(N-1)"      (E-(N-1))
 *
 * Update rule:
 *   delta < -(N-1): ignore
 *   delta in [-(N-1), 0]: set bit (N-1+delta) of history[histBit]
 *   delta > 0: advance exp, shift both tails right by min(delta,N),
 *              set MSB of history[histBit]
 *
 * Reports P(bit=1) for each bit position at each tier (one tail only,
 * since both tails are symmetric by construction).
 *
 * Usage: java cardinality.TailWidthSimulator [tailbits=2] [iters=100000]
 *        [maxcard=10000000] [threads=4]
 *
 * @author Chloe
 * @date April 22, 2026
 */
public class TailWidthSimulator {

	public static void main(String[] args) throws InterruptedException {
		int tailBits=2;
		int iters=100000;
		long maxCard=10_000_000L;
		int threads=4;

		for(String arg : args){
			String[] kv=arg.split("=",2);
			if(kv.length!=2){continue;}
			switch(kv[0]){
				case "tailbits": tailBits=Integer.parseInt(kv[1]); break;
				case "iters": iters=Integer.parseInt(kv[1]); break;
				case "maxcard": maxCard=Long.parseLong(kv[1]); break;
				case "threads": threads=Integer.parseInt(kv[1]); break;
			}
		}

		final int tb=tailBits;
		final int nIters=iters;
		final long mc=maxCard;
		final int nThreads=threads;
		final int NUM_TIERS=16;
		final int tailMask=(1<<tb)-1;

		// Per-thread accumulators
		final long[][][] tBitSet=new long[nThreads][NUM_TIERS][tb];
		final long[][] tTierObs=new long[nThreads][NUM_TIERS];

		final int base=nIters/nThreads;
		final int rem=nIters%nThreads;

		Thread[] threadArr=new Thread[nThreads];
		for(int t=0; t<nThreads; t++){
			final int tid=t;
			final int myIters=base+(tid<rem ? 1 : 0);
			threadArr[t]=new Thread(()->{
				final FastRandomXoshiro rng=new FastRandomXoshiro(tid+1);
				final long[][] lBitSet=tBitSet[tid];
				final long[] lTierObs=tTierObs[tid];

				for(int iter=0; iter<myIters; iter++){
					int exp=0;
					int h0=0, h1=0;

					for(long card=1; card<=mc; card++){
						final long hash=TTLLSimulator.hash64shift(rng.nextLong());
						final int nlz=Long.numberOfLeadingZeros(hash);
						final int hb=(int)(hash&1);
						final int delta=nlz-exp;

						if(delta<-(tb-1)){
							// ignore
						}else if(delta>0){
							final int newExp=Math.min(exp+delta, 15);
							final int shift=Math.min(delta, tb);
							h0=(h0>>>shift)&tailMask;
							h1=(h1>>>shift)&tailMask;
							if(hb==0){h0|=(1<<(tb-1));}
							else     {h1|=(1<<(tb-1));}
							exp=newExp;
						}else{
							// delta in [-(tb-1), 0]
							final int bitPos=(tb-1)+delta;
							if(hb==0){h0|=(1<<bitPos);}
							else     {h1|=(1<<bitPos);}
						}

						// Record h0 state (h1 is symmetric)
						if(exp<NUM_TIERS){
							lTierObs[exp]++;
							for(int b=0; b<tb; b++){
								if(((h0>>>b)&1)!=0){lBitSet[exp][b]++;}
							}
						}
						if(exp>14){break;}
					}
				}
			});
			threadArr[t].start();
		}
		for(Thread th : threadArr){th.join();}

		// Merge
		long[][] bitSet=new long[NUM_TIERS][tb];
		long[] tierObs=new long[NUM_TIERS];
		for(int t=0; t<nThreads; t++){
			for(int tier=0; tier<NUM_TIERS; tier++){
				tierObs[tier]+=tTierObs[t][tier];
				for(int b=0; b<tb; b++){
					bitSet[tier][b]+=tBitSet[t][tier][b];
				}
			}
		}

		// Print
		System.out.printf("TailWidthSimulator: tailbits=%d iters=%d maxcard=%d threads=%d%n",
			tb, nIters, mc, nThreads);
		System.out.println();

		StringBuilder hdr=new StringBuilder(String.format("%4s", "Tier"));
		for(int b=tb-1; b>=0; b--){
			hdr.append(String.format("  %8s", "E-"+(tb-1-b)));
		}
		System.out.println(hdr);

		for(int tier=0; tier<NUM_TIERS; tier++){
			if(tierObs[tier]==0){continue;}
			StringBuilder row=new StringBuilder(String.format("%4d", tier));
			for(int b=tb-1; b>=0; b--){
				double p=(double)bitSet[tier][b]/tierObs[tier];
				row.append(String.format("  %8.4f", p));
			}
			System.out.println(row);
		}

		// Steady-state
		System.out.println();
		System.out.print("  ss");
		for(int b=tb-1; b>=0; b--){
			long setSum=0, obsSum=0;
			for(int tier=3; tier<=8; tier++){
				setSum+=bitSet[tier][b];
				obsSum+=tierObs[tier];
			}
			System.out.printf("  %8.4f", obsSum>0 ? (double)setSum/obsSum : 0);
		}
		System.out.println();
	}
}
