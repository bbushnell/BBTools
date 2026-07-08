package cardinality;

/**
 * WidthWt53C: chronicle-blend referee. THE HONEST REFEREE TABLE.  1024 bytes each, every bit paid for:
 * AVLL (1408 regs) vs Fll53 (512 words = 2560 tails + in-word counters).
 * Estimators: AVLL+HLDLC; Fll53+MLE (likelihood, dup-blind); Fll53 BLENDED —
 * cHat from the modular window net (k=2 windows averaged), then
 * w = clamp((cHat-0.85)/0.10, 0, 1) and nHat = w*MLE + (1-w)*adds*cHat
 * (smooth handoff replaces the twitchy 0.95 threshold that cost FLL4 its HC).
 * Brian's metrics: LogWt/WidthWt/CountWt/Peak/AvgSigned over the reportFrac=0.01
 * ladder, HC + LC(minrand iter=4), clampToAdded=false.
 *
 * Usage: java cardinality.WidthWt53 [instances] [threads] [table.tsv] [net.txt] [maxMult]
 *        [orgs] [chronNet] [lcPatterns] [trajC 0/1] [totalBytes 0/1]
 * trajC=1: chron lane uses the 37-byte fp16 Fll53TrajC instead of the ~2KB
 *   Traj (ablation-validated, drift <<referee noise).  totalBytes=1: AVLL is
 *   sized to the chron lane's TOTAL memory (sketch + STATE_BYTES) — the
 *   equal-TOTAL-bytes crown (Brian 2026-07-05; AVLL's own ~17B of vestigial/
 *   optional state waived as tiny, per his sizing convention).  Both default
 *   0 = old behavior exactly.
 *
 * @author Amber
 * @date July 2026
 */
public class WidthWt53C {

	static final long GOLD=0x9E3779B97F4A7C15L;
	static final double LOG2=Math.log(2.0);
	// Sizes set from args: honest equal bytes (AVLL regs = words*11/8 exactly
	// when words is a multiple of 8: bytes = words*2 = regs*8/11).
	// 21/20-bit organisms pack 3 per long: bytes = orgs*8/3.
	static int FLL_ORGS=384;
	static int FLL_BUCKETS=FLL_ORGS*5;
	static int AVLL_REGS=1408;
	static final int NUM_EST=4;
	static final String[] EST_NAMES={"AVLL_HLDLC", "FLL53_MLE", "FLL53_blend",
		"FLL53_chron"};

	public static void main(String[] args) throws Exception{
		final int instances=(args.length>0) ? Integer.parseInt(args[0]) : 64;
		final int threads=(args.length>1) ? Integer.parseInt(args[1])
			: Runtime.getRuntime().availableProcessors();
		final String tablePath=(args.length>2) ? args[2] : "fll53mle_384o.tsv";
		final String netPath=(args.length>3) ? args[3] : "netwin53.txt";
		final long maxMult=(args.length>4) ? Long.parseLong(args[4]) : 2048;
		final String chronPath=(args.length>6) ? args[6] : "netchron53b.txt";
		// LC patterns, comma-separated.  minrand = the historical crown pattern
		// (identical seeds/streams to WidthWt53); the others exercise
		// time-VARYING complexity, which a single constant-complexity pool
		// never tests.  The crown is the average across patterns.
		final String[] lcPats=((args.length>7) ? args[7]
			: "minrand,unifpool,zipfish,onehot").split(",");
		if(args.length>5){
			FLL_ORGS=Integer.parseInt(args[5]);
			FLL_BUCKETS=FLL_ORGS*5;
			AVLL_REGS=(int)(FLL_ORGS*8L/3*11/8);   // equal bytes
		}
		final boolean trajC=(args.length>8) && Integer.parseInt(args[8])!=0;
		final boolean totalBytes=(args.length>9) && Integer.parseInt(args[9])!=0;
		if(totalBytes){   // equal TOTAL bytes: AVLL matches sketch+chronicle state
			AVLL_REGS=(int)((FLL_ORGS*8L/3+Fll53TrajC.STATE_BYTES)*11L/8);
		}
		System.err.println("HONEST bytes="+(FLL_ORGS*8/3)+"  orgs="+FLL_ORGS
			+"  avllRegs="+AVLL_REGS+"  trajC="+trajC+"  totalBytes="+totalBytes
			+(totalBytes ? " (chron total="+(FLL_ORGS*8/3+Fll53TrajC.STATE_BYTES)+"B)" : ""));
		CardinalityTracker.clampToAdded=false;

		final Fll53MLE mle=new Fll53MLE(tablePath);
		final ComplexityFarm.BigNet net=ComplexityFarm.BigNet.load(netPath);
		final ComplexityFarm.BigNet chron=ComplexityFarm.BigNet.load(chronPath);
		System.err.println("chron net "+chronPath+" (D="+chron.D+", act="+chron.act+")");
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

		final double[][] lcWW=new double[NUM_EST][lcPats.length];
		int lcDone=0;
		for(int mode=0; mode<1+lcPats.length; mode++){
			final boolean lc=(mode>=1);
			final String pat=lc ? lcPats[mode-1] : "HC";
			// minrand keeps the historical 555000 base (bit-comparable with all
			// prior crown tables); other patterns get distinct seed spaces.
			final int patSeed=!lc ? 111000
				: "minrand".equals(pat) ? 555000 : 555000+7777*mode;
			final double[][] sumAbs=new double[NUM_EST][NT];
			final double[][] sumSigned=new double[NUM_EST][NT];
			final long[] cnt=new long[NT];
			final Object lock=new Object();

			FllSkewTest.runPool(threads, instances, inst->{
				final ArithmeticVariableLogLog avll=new ArithmeticVariableLogLog(AVLL_REGS, 31, -1, 0);
				final Fll53 fll=new Fll53(FLL_BUCKETS, 31, -1, 0);
				final Fll53Chron.Traj tr=trajC ? null
					: new Fll53Chron.Traj(net, mle, 4L*FLL_BUCKETS);
				final Fll53TrajC trC=trajC
					? new Fll53TrajC(net, mle, 4L*FLL_BUCKETS) : null;
				final StreamSource src=StreamSource.make(pat, maxTrue,
					patSeed+inst*104729L);
				final double[][] lAbs=new double[NUM_EST][NT];
				final double[][] lSig=new double[NUM_EST][NT];
				final long[] lCnt=new long[NT];
				long draws=0;
				final long maxDraws=src.maxDraws(maxTrue);
				int ti=0;
				while(draws<maxDraws && ti<NT){
					final long key=src.next();
					avll.add(key);
					fll.add(key);
					draws++;
					if(trajC){trC.tick(fll, draws);}else{tr.tick(fll, draws);}
					if(ti<NT && src.distinct()==thresholds[ti]){
						final double t=src.distinct();
						final double[] ests=new double[NUM_EST];
						ests[0]=avll.cardinality();
						ests[1]=mle.estimate(fll);
						// Blended: modular net cHat (k windows), smooth handoff
						double ySum=0;
						int k=0;
						for(int w0=0; w0+Fll53FarmW.WIN<=FLL_ORGS; w0+=Fll53FarmW.WIN){
							ySum+=net.predict(Fll53FarmW.features(fll, draws, w0));
							k++;
						}
						if(k==0 || ests[1]<32.0*FLL_BUCKETS){
							ests[2]=ests[1];   // cold or sub-window sketch
							ests[3]=ests[1];
						}else{
							final double cHat=Math.pow(2.0, ySum/k);
							final double w=Math.max(0, Math.min(1, (cHat-0.85)/0.10));
							ests[2]=w*ests[1]+(1-w)*draws*cHat;
							// Chronicle lane: trajectory-fused complexity, same ramp
							final double[] ff=trajC
								? Fll53TrajC.fuseC(trC, fll, draws, FLL_BUCKETS)
								: Fll53Chron.fuse(tr, fll, draws, FLL_BUCKETS);
							final double cC=Math.pow(2.0, chron.predict(ff));
							final double wc=Math.max(0, Math.min(1, (cC-0.85)/0.10));
							ests[3]=wc*ests[1]+(1-wc)*draws*cC;
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

			System.out.println("\n=== "+(lc ? ("LOW complexity ["+pat+"]") : "HIGH complexity")
				+", HONEST equal bytes, "+instances+" instances, maxTrue="+maxTrue+" ===");
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
				if(lc){lcWW[e][mode-1]=widthWt/widthSum;}
			}
			if(lc){lcDone++;}
		}
		// The crown score: LC WidthWt averaged across patterns.
		System.out.println("\n=== LC AVERAGE (WidthWtAbsErr over "+lcDone+" patterns) ===");
		for(int e=0; e<NUM_EST; e++){
			double sum=0;
			for(int p=0; p<lcDone; p++){sum+=lcWW[e][p];}
			System.out.println(String.format("%-12s %.8f", EST_NAMES[e], sum/lcDone));
		}
	}
}
