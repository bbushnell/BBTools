package cardinality;

/**
 * WidthWt52: THE HONEST REFEREE TABLE.  1024 bytes each, every bit paid for:
 * AVLL (1408 regs) vs Fll52 (512 words = 2560 tails + in-word counters).
 * Estimators: AVLL+HLDLC; Fll52+MLE (likelihood, dup-blind); Fll52 BLENDED —
 * cHat from the modular window net (k=2 windows averaged), then
 * w = clamp((cHat-0.85)/0.10, 0, 1) and nHat = w*MLE + (1-w)*adds*cHat
 * (smooth handoff replaces the twitchy 0.95 threshold that cost FLL4 its HC).
 * Brian's metrics: LogWt/WidthWt/CountWt/Peak/AvgSigned over the reportFrac=0.01
 * ladder, HC + LC(minrand iter=4), clampToAdded=false.
 *
 * Usage: java cardinality.WidthWt52 [instances] [threads] [table.tsv] [net.txt] [maxMult]
 *
 * @author Amber
 * @date July 2026
 */
public class WidthWt52 {

	static final long GOLD=0x9E3779B97F4A7C15L;
	static final double LOG2=Math.log(2.0);
	// Sizes set from args: honest equal bytes (AVLL regs = words*11/8 exactly
	// when words is a multiple of 8: bytes = words*2 = regs*8/11).
	static int FLL_WORDS=512;
	static int FLL_BUCKETS=FLL_WORDS*5;
	static int AVLL_REGS=1408;
	static final int NUM_EST=3;
	static final String[] EST_NAMES={"AVLL_HLDLC", "FLL52_MLE", "FLL52_blend"};

	public static void main(String[] args) throws Exception{
		final int instances=(args.length>0) ? Integer.parseInt(args[0]) : 64;
		final int threads=(args.length>1) ? Integer.parseInt(args[1])
			: Runtime.getRuntime().availableProcessors();
		final String tablePath=(args.length>2) ? args[2] : "fll52mle_table.tsv";
		final String netPath=(args.length>3) ? args[3] : "netwin52.txt";
		final long maxMult=(args.length>4) ? Long.parseLong(args[4]) : 2048;
		if(args.length>5){
			FLL_WORDS=Integer.parseInt(args[5]);
			FLL_BUCKETS=FLL_WORDS*5;
			AVLL_REGS=(int)(FLL_WORDS*2L*11/8);   // equal bytes: bytes*11 regs per 8 bytes
		}
		System.err.println("HONEST bytes="+(FLL_WORDS*2)+"  fllWords="+FLL_WORDS
			+"  avllRegs="+AVLL_REGS);
		CardinalityTracker.clampToAdded=false;

		final Fll52MLE mle=new Fll52MLE(tablePath);
		final ComplexityFarm.BigNet net=ComplexityFarm.BigNet.load(netPath);
		final long maxTrue=maxMult*FLL_BUCKETS;

		final java.util.ArrayList<Long> tl=new java.util.ArrayList<>();
		long next=1;
		while(next<maxTrue){
			tl.add(next);
			next=Math.max(next+1, (long)(next*1.01));
		}
		tl.add(maxTrue);
		final long[] thresholds=new long[tl.size()];
		for(int i=0; i<thresholds.length; i++){thresholds[i]=tl.get(i);}
		final int NT=thresholds.length;

		for(int mode=0; mode<2; mode++){
			final boolean lc=(mode==1);
			final double[][] sumAbs=new double[NUM_EST][NT];
			final double[][] sumSigned=new double[NUM_EST][NT];
			final long[] cnt=new long[NT];
			final Object lock=new Object();

			FllSkewTest.runPool(threads, instances, inst->{
				final ArithmeticVariableLogLog avll=new ArithmeticVariableLogLog(AVLL_REGS, 31, -1, 0);
				final Fll52 fll=new Fll52(FLL_BUCKETS, 31, -1, 0);
				final java.util.Random seedr=new java.util.Random((lc ? 555000L : 111000L)+inst*104729L);
				final long seedBase=seedr.nextLong();
				final java.util.SplittableRandom rng=new java.util.SplittableRandom(seedr.nextLong());
				final double[][] lAbs=new double[NUM_EST][NT];
				final double[][] lSig=new double[NUM_EST][NT];
				final long[] lCnt=new long[NT];
				final boolean[] seen=lc ? new boolean[(int)maxTrue] : null;
				long distinct=0, draws=0;
				final long maxDraws=lc ? 4L*maxTrue : maxTrue;
				int ti=0;
				while(draws<maxDraws && ti<NT){
					final long key;
					if(lc){
						final int pos=Math.min(rng.nextInt((int)maxTrue), rng.nextInt((int)maxTrue));
						if(!seen[pos]){seen[pos]=true; distinct++;}
						key=seedBase+pos*GOLD;
					}else{
						distinct++;
						key=seedBase+distinct*GOLD;
					}
					avll.add(key);
					fll.add(key);
					draws++;
					if(ti<NT && distinct==thresholds[ti]){
						final double t=distinct;
						final double[] ests=new double[NUM_EST];
						ests[0]=avll.cardinality();
						ests[1]=mle.estimate(fll);
						// Blended: modular net cHat (k windows), smooth handoff
						double ySum=0;
						int k=0;
						for(int w0=0; w0+Fll52FarmW.WIN<=FLL_WORDS; w0+=Fll52FarmW.WIN){
							ySum+=net.predict(Fll52FarmW.features(fll, draws, w0));
							k++;
						}
						final double cHat=Math.pow(2.0, ySum/k);
						if(ests[1]<32.0*FLL_BUCKETS){
							ests[2]=ests[1];   // cold sketch: trust likelihood
						}else{
							final double w=Math.max(0, Math.min(1, (cHat-0.85)/0.10));
							ests[2]=w*ests[1]+(1-w)*draws*cHat;
						}
						for(int e=0; e<NUM_EST; e++){
							final double err=(ests[e]-t)/t;
							lAbs[e][ti]+=Math.abs(err);
							lSig[e][ti]+=err;
						}
						lCnt[ti]++;
						ti++;
					}
				}
				synchronized(lock){
					for(int e=0; e<NUM_EST; e++){
						for(int i=0; i<NT; i++){
							sumAbs[e][i]+=lAbs[e][i]; sumSigned[e][i]+=lSig[e][i];
						}
					}
					for(int i=0; i<NT; i++){cnt[i]+=lCnt[i];}
				}
			});

			System.out.println("\n=== "+(lc ? "LOW" : "HIGH")+" complexity, HONEST 1KB each, "
				+instances+" instances, maxTrue="+maxTrue+(lc ? " (pool), iter=4" : "")+" ===");
			System.out.println(String.format("%-12s %-13s %-13s %-13s %-13s %s",
				"", "LogWtAbsErr", "WidthWtAbsErr", "CountWtAbsErr", "PeakAbsErr", "AvgSignErr"));
			for(int e=0; e<NUM_EST; e++){
				double logWt=0, widthWt=0, cntWt=0, peak=0, signed=0;
				double widthSum=0, cntSum=0;
				int rows=0;
				long prev=0;
				for(int i=0; i<NT; i++){
					final double w=Math.max(1, thresholds[i]-prev);
					prev=thresholds[i];
					if(cnt[i]==0){continue;}
					final double meanAbs=sumAbs[e][i]/cnt[i];
					logWt+=meanAbs; widthWt+=meanAbs*w; cntWt+=meanAbs*cnt[i];
					signed+=sumSigned[e][i]/cnt[i];
					widthSum+=w; cntSum+=cnt[i];
					if(meanAbs>peak){peak=meanAbs;}
					rows++;
				}
				System.out.println(String.format("%-12s %.8f   %.8f   %.8f   %.8f   %+.8f",
					EST_NAMES[e], logWt/rows, widthWt/widthSum, cntWt/cntSum, peak, signed/rows));
			}
		}
	}
}
