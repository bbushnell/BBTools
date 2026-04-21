package cardinality;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.concurrent.atomic.AtomicInteger;

import rand.FastRandomXoshiro;

/**
 * HCDLCTierAccuracy - measure per-tier accuracy of DLC and HC estimates
 * across sampling modes and averaging modes for UDLL6.
 * <p>
 * Samples the internal state of a {@link UltraDynamicLogLog6} tracker at
 * selected points during a DDL run, computes per-tier DLC and HC
 * estimates, and bins absolute error by integer-valued axis:
 * <ul>
 *   <li>DLC: binned by V = B - filled_at_or_above_t</li>
 *   <li>HC (Beff axis): binned by hcBeff_t</li>
 *   <li>HC (Unseen axis): binned by hcUnseen_t</li>
 * </ul>
 * For each bin the simulator accumulates linear / geometric / harmonic
 * sums of |err|, and mean signed error.  All three flavours are written
 * to the output in separate sections so a single run produces all three
 * averaging modes.
 * <p>
 * Sampling modes:
 * <ul>
 *   <li><b>logcard</b> - record at log-spaced cardinality checkpoints
 *     (the current DLCTierAccuracy mode, fully supported in v1)</li>
 *   <li><b>entry</b> - record just after each state change
 *     (detected via {@code lastCardinality==-1})</li>
 *   <li><b>entryexit</b> - record twice per state change: the new state
 *     (entry) and the previous state about to be left (exit)</li>
 *   <li><b>add</b> - record at every insertion, weighted by dwell
 *     (implemented as entry-mode with dwell-count weighting)</li>
 * </ul>
 * Supports hbits in {0, 1, 2}; hbits=3 is accepted as a flag but not
 * executed in v1 (no UDLL7 yet).
 *
 * @author Chloe
 * @date April 14, 2026
 */
public class HCDLCTierAccuracy {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum tier index examined per sample. */
	static final int MAX_TIER=50;
	/** UDLL6 history margin (see UltraDynamicLogLog6.HISTORY_MARGIN). */
	static final int HISTORY_MARGIN=2;

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws Exception {
		int buckets=2048;
		int ddls=100_000;
		int maxmult=8192;
		long seed=-1;
		int hbits=2;
		int threads=Math.max(1, Runtime.getRuntime().availableProcessors());
		int points=300;
		String sample="logcard";
		String outPrefix="hcdlc_";
		double timeSec=0;// 0 = use ddls; >0 = run until wall time exceeded
		String tableinDlc=null;// path to V-indexed DLC accuracy table
		String tableinHc=null;// path to hcBeff-indexed HC accuracy table
		String tableinHcUnseen=null;// path to hcUnseen-indexed HC accuracy table
		double powerDlc=4.5;
		double powerHc=2.0;
		String avgCol="geo";// which avgAbsErr column to read: lin|geo|harm
		boolean dlcFormula=false;// true = Mode 2 formula, ignores tablein
		boolean hcFormula=false;// true = Mode 2 formula, ignores tablein_hc
		String hcWeight="auto";// HC weighting: auto|beff|unseen|geo|max|v|formula
		String trackerName="udll6";// tracker: udll6|cdll4

		for(String arg : args){
			final String[] kv=arg.split("=", 2);
			if(kv.length<2){throw new RuntimeException("Unknown parameter '"+arg+"'");}
			final String a=kv[0].toLowerCase(), b=kv[1];
			switch(a){
				case "buckets": case "b": buckets=Integer.parseInt(b); break;
				case "ddls":    ddls=Integer.parseInt(b); break;
				case "maxmult": maxmult=Integer.parseInt(b); break;
				case "seed":    seed=Long.parseLong(b); break;
				case "hbits":   hbits=Integer.parseInt(b); break;
				case "threads": case "t": threads=Integer.parseInt(b); break;
				case "points":  points=Integer.parseInt(b); break;
				case "sample":  sample=b.toLowerCase(); break;
				case "time":    timeSec=Double.parseDouble(b); break;
				case "tablein_dlc": case "tablein": case "tabledlc": tableinDlc=b; break;
				case "tablein_hc": case "tablein_hc_beff": case "tablehc": tableinHc=b; break;
				case "tablein_hc_unseen": case "tablehcunseen": tableinHcUnseen=b; break;
				case "power_dlc": case "pdlc": powerDlc=Double.parseDouble(b); break;
				case "power_hc": case "phc": powerHc=Double.parseDouble(b); break;
				case "power": powerDlc=powerHc=Double.parseDouble(b); break;
				case "avg": avgCol=b.toLowerCase(); break;
				case "dlcformula": case "formula_dlc":
					dlcFormula=parseBool(b); break;
				case "hcformula": case "formula_hc":
					hcFormula=parseBool(b); break;
				case "formula": dlcFormula=hcFormula=parseBool(b); break;
				case "hc_weight": case "hcweight": case "hc_combine":
					hcWeight=b.toLowerCase(); break;
				case "tracker": case "type": case "loglogtype":
					trackerName=b.toLowerCase(); break;
				case "out_prefix": case "out": case "outprefix": outPrefix=b; break;
				default: throw new RuntimeException("Unknown parameter '"+arg+"'");
			}
		}

		if(hbits<0 || hbits>3){
			throw new RuntimeException("hbits must be in [0,3]; got "+hbits);
		}
		if(hbits==3){
			throw new RuntimeException("hbits=3 is not supported in v1 "+
				"(no UDLL7 class yet).  Simulator parameterizes for hbits=3 "+
				"but execution is deferred.");
		}
		if(!sample.equals("logcard") && !sample.equals("entry") &&
		   !sample.equals("entryexit") && !sample.equals("add")){
			throw new RuntimeException("Unknown sample mode '"+sample+
				"'.  Expected logcard|entry|entryexit|add.");
		}

		// Tracker-specific parameters.  UDLL6: tier width 2^1, up to 2 history
		// bits per bucket.  CDLL4: tier width 2^1.5 (half-NLZ with (2-sqrt(2))
		// mantissa threshold), 1 history bit per bucket.
		final int trackerType;
		final double tierScale;
		final int maxHbitsSupported;
		switch(trackerName){
			case "udll6":
				trackerType=0; tierScale=1.0; maxHbitsSupported=2; break;
			case "cdll4":
				trackerType=1; tierScale=1.5; maxHbitsSupported=1; break;
			default:
				throw new RuntimeException("Unknown tracker='"+trackerName+
					"' (udll6|cdll4)");
		}
		if(hbits>maxHbitsSupported){
			throw new RuntimeException("tracker="+trackerName+" supports at most "+
				maxHbitsSupported+" history bits; got hbits="+hbits);
		}
		if(trackerType==1 && dlcFormula){
			throw new RuntimeException("dlcformula=true not supported for cdll4 "+
				"(Mode 2 formula is derived for UDLL6 tier width).  Use tablein_dlc= "+
				"instead.");
		}
		if(trackerType==1 && hcFormula){
			throw new RuntimeException("hcformula=true not supported for cdll4. "+
				"Use tablein_hc=/tablein_hc_unseen= instead.");
		}

		final long maxCard=(long)buckets*maxmult;
		final int nBins=buckets+1;
		final long masterSeed=(seed>=0 ? seed : System.nanoTime());

		// Precompute log-spaced cardinality checkpoints
		final long[] thresh=logThresh(maxCard, points);
		final int nThresh=thresh.length;

		// Time-limit mode: set ddls to effectively infinite; workers will exit
		// when wall time exceeds timeSec.
		final long timeLimitMs=(timeSec>0 ? (long)(timeSec*1000) : Long.MAX_VALUE);
		if(timeSec>0){ddls=Integer.MAX_VALUE;}

		// Load error tables for blend computation.
		// errByV[v] = DLC absolute error when V=v (from a previous sim).
		// errByBeff[b] = HC absolute error when hcBeff=b.
		// In formula mode the table is not consulted — Mode 2 theoretical formula is used.
		final double[] errByV=(tableinDlc!=null && !dlcFormula ? loadErrTable(tableinDlc, avgCol, nBins) : null);
		final double[] errByBeff=(tableinHc!=null && !hcFormula ? loadErrTable(tableinHc, avgCol, nBins) : null);
		final double[] errByUnseen=(tableinHcUnseen!=null && !hcFormula ? loadErrTable(tableinHcUnseen, avgCol, nBins) : null);
		final boolean doBlend=(errByV!=null || dlcFormula);
		final boolean doHcBlend=(errByBeff!=null || errByUnseen!=null || hcFormula);

		System.err.println("HCDLCTierAccuracy: type="+trackerName+" B="+buckets+
			" hbits="+hbits+" ddls="+(timeSec>0 ? "time="+timeSec+"s" : String.valueOf(ddls))+
			" maxmult="+maxmult+" sample="+sample+" threads="+threads+
			" points="+nThresh+" maxCard="+maxCard+" seed="+masterSeed+
			" tierScale="+tierScale+
			(doBlend ? " dlc="+(dlcFormula?"formula":"table")+
				" power_dlc="+powerDlc+
				(doHcBlend ? " hc="+(hcFormula?"formula":"table")+
					" power_hc="+powerHc : "")+
				" avg="+avgCol : " blend=OFF"));

		// Global accumulators [5 stats] × [nBins]:
		//   0: count, 1: sum(absErr), 2: sum(ln(absErr)), 3: sum(1/absErr),
		//   4: sum(signedErr)
		final double[][] dlcAcc=new double[5][nBins];
		final double[][] hcBeffAcc=new double[5][nBins];
		final double[][] hcUnseenAcc=new double[5][nBins];

		// Blend-error accumulators: log2(trueCard) bins (0..64).
		// Separate DLC-blend and HC-blend end-to-end error as f(cardinality).
		final int nCardBins=64;
		final double[][] dlcBlendAcc=new double[5][nCardBins];
		final double[][] hcBlendAcc=new double[5][nCardBins];

		final AtomicInteger ddlCounter=new AtomicInteger(0);
		final long t0=System.currentTimeMillis();

		final int finalBuckets=buckets;
		final int finalHbits=hbits;
		final int finalDdls=ddls;
		final long finalMaxCard=maxCard;
		final String finalSample=sample;
		final long finalSeed=masterSeed;
		final int finalTrackerType=trackerType;
		final double finalTierScale=tierScale;

		final double[] finalErrByV=errByV;
		final double[] finalErrByBeff=errByBeff;
		final double[] finalErrByUnseen=errByUnseen;
		final double finalPowerDlc=powerDlc;
		final double finalPowerHc=powerHc;
		final boolean finalDlcFormula=dlcFormula;
		// Resolve hcWeight enum: 0=auto 1=beff 2=unseen 3=geo 4=max 5=V 6=formula
		final int finalHcWeight;
		if(hcFormula){finalHcWeight=6;}
		else{
			switch(hcWeight){
				case "auto":    finalHcWeight=0; break;
				case "beff":    finalHcWeight=1; break;
				case "unseen":  finalHcWeight=2; break;
				case "geo":     finalHcWeight=3; break;
				case "max":     finalHcWeight=4; break;
				case "v":       finalHcWeight=5; break;
				case "formula": finalHcWeight=6; break;
				default: throw new RuntimeException("Unknown hc_weight='"+hcWeight+
					"' (auto|beff|unseen|geo|max|v|formula)");
			}
		}

		final Thread[] workers=new Thread[threads];
		for(int ti=0; ti<threads; ti++){
			final int tid=ti;
			workers[ti]=new Thread(){
				@Override public void run(){
					final double[][] lDlc=new double[5][nBins];
					final double[][] lBeff=new double[5][nBins];
					final double[][] lUns=new double[5][nBins];
					final double[][] lDlcBlend=new double[5][nCardBins];
					final double[][] lHcBlend=new double[5][nCardBins];
					int d;
					while((d=ddlCounter.getAndIncrement())<finalDdls){
						// Time-limit exit: check before starting a new DDL to avoid
						// partial (and thus biased) contribution.
						if(System.currentTimeMillis()-t0>=timeLimitMs){break;}
						final long instSeed=finalSeed+d*12345678901L+tid*97L;
						final CardinalityTracker ct=(finalTrackerType==1)
							? new CompressedDynamicLogLog4(finalBuckets, 31, instSeed, 0)
							: new UltraDynamicLogLog6(finalBuckets, 31, instSeed, 0);
						final FastRandomXoshiro rng=new FastRandomXoshiro(instSeed^0xDEADBEEFL);
						runOneDDL(ct, rng, finalBuckets, finalHbits, finalMaxCard,
							finalSample, thresh, lDlc, lBeff, lUns,
							finalErrByV, finalErrByBeff, finalErrByUnseen,
							finalPowerDlc, finalPowerHc,
							finalDlcFormula, finalHcWeight,
							finalTrackerType, finalTierScale,
							lDlcBlend, lHcBlend);

						if((d+1)%Math.max(1, finalDdls/20)==0){
							final long el=System.currentTimeMillis()-t0;
							System.err.println("  DDL "+(d+1)+"/"+finalDdls+
								" ("+String.format("%.1f", el/1000.0)+"s)");
						}
					}
					synchronized(dlcAcc){
						mergeInto(dlcAcc, lDlc, nBins);
						mergeInto(hcBeffAcc, lBeff, nBins);
						mergeInto(hcUnseenAcc, lUns, nBins);
						mergeInto(dlcBlendAcc, lDlcBlend, nCardBins);
						mergeInto(hcBlendAcc, lHcBlend, nCardBins);
					}
				}
			};
			workers[ti].start();
		}
		for(Thread w : workers){w.join();}

		final long elapsed=System.currentTimeMillis()-t0;
		final int completedDdls=ddlCounter.get();
		System.err.println("Simulation complete in "+
			String.format("%.1f", elapsed/1000.0)+"s ("+
			Math.min(completedDdls, ddls)+" DDLs completed)");

		// Emit 3 (or 5) output files
		final String hdr=String.format(
			"B=%d hbits=%d ddls=%d maxmult=%d sample=%s sim=UDLL6 "+
			"seed=%d simMs=%d%s",
			buckets, hbits, ddls, maxmult, sample, masterSeed, elapsed,
			(doBlend ? " tablein="+tableinDlc+
				(tableinHc!=null ? " tablein_hc="+tableinHc : "")+
				" power_dlc="+powerDlc+
				(errByBeff!=null ? " power_hc="+powerHc : "")+
				" avg="+avgCol : ""));
		writeTable(outPrefix+"dlc_v.tsv", hdr, "V", dlcAcc, nBins);
		if(hbits>=1){
			writeTable(outPrefix+"hc_beff.tsv", hdr, "hcBeff", hcBeffAcc, nBins);
			writeTable(outPrefix+"hc_unseen.tsv", hdr, "hcUnseen", hcUnseenAcc, nBins);
		}
		if(doBlend){
			writeTable(outPrefix+"dlc_blend.tsv", hdr, "log2Card", dlcBlendAcc, nCardBins);
			if(doHcBlend){
				writeTable(outPrefix+"hc_blend.tsv", hdr, "log2Card", hcBlendAcc, nCardBins);
			}
		}

		// Quick stderr summary at bin count midpoint
		summarize(System.err, "DLC  (V)       ", dlcAcc, nBins);
		if(hbits>=1){
			summarize(System.err, "HC   (Beff)    ", hcBeffAcc, nBins);
			summarize(System.err, "HC   (Unseen)  ", hcUnseenAcc, nBins);
		}
		if(doBlend){
			summarizeCard(System.err, "DLCblend(card) ", dlcBlendAcc, nCardBins);
			if(doHcBlend){
				summarizeCard(System.err, "HC blend(card) ", hcBlendAcc, nCardBins);
			}
			// Single-scalar terminal-weighted average absolute error:
			// sum of (weight * |err|) over all samples, divided by sum of weight.
			// Under dwell weighting this is the mean error across cardinality,
			// emphasizing high-card/steady-state where the dwell is longer.
			System.err.println();
			System.err.printf("TERMINAL  dlcBlend=%.4f%%   hcBlend=%.4f%%%n",
				100*scalarWeightedAvg(dlcBlendAcc),
				(doHcBlend ? 100*scalarWeightedAvg(hcBlendAcc) : 0));
		}
		// Also print per-tier scalar for reference (regardless of blend)
		System.err.printf("           dlcPerTier=%.4f%%   hcBeffPerTier=%.4f%%%n",
			100*scalarWeightedAvg(dlcAcc),
			(hbits>=1 ? 100*scalarWeightedAvg(hcBeffAcc) : 0));
	}

	/*--------------------------------------------------------------*/
	/*----------------           Per-DDL            ----------------*/
	/*--------------------------------------------------------------*/

	/** Run one DDL: feed up to maxCard random longs, sampling state per mode. */
	static void runOneDDL(CardinalityTracker ct, FastRandomXoshiro rng,
			int B, int hbits, long maxCard, String sample,
			long[] thresh, double[][] dlcAcc, double[][] beffAcc, double[][] unsAcc,
			double[] errByV, double[] errByBeff, double[] errByUnseen,
			double powerDlc, double powerHc,
			boolean dlcFormula, int hcWeightMode,
			int trackerType, double tierScale,
			double[][] dlcBlendAcc, double[][] hcBlendAcc){
		final boolean doBlend=(errByV!=null || dlcFormula);

		// Per-DDL extracted vectors (valid at last sample point)
		final int[] nlzBucketCount=new int[64];
		final int[][] nlzHbitSet=(hbits>0 ? new int[hbits][64] : null);

		if(sample.equals("logcard")){
			long added=0;
			long prevTarget=0;
			for(int ti=0; ti<thresh.length; ti++){
				final long target=Math.min(thresh[ti], maxCard);
				while(added<target){
					ct.hashAndStore(rng.nextLong());
					added++;
				}
				// Dwell weight: cardinality units since the previous checkpoint.
				// This makes every unit of true cardinality contribute equally,
				// which biases the aggregate toward steady-state (terminal) perf.
				final double w=Math.max(1, target-prevTarget);
				prevTarget=target;
				extractState(ct, trackerType, B, hbits, nlzBucketCount, nlzHbitSet);
				recordSampleWeighted(added, B, hbits, tierScale, nlzBucketCount, nlzHbitSet, w,
					dlcAcc, beffAcc, unsAcc);
				if(doBlend){
					recordBlend(added, B, hbits, tierScale, nlzBucketCount, nlzHbitSet,
						errByV, errByBeff, errByUnseen, powerDlc, powerHc, dlcFormula, hcWeightMode, w,
						dlcBlendAcc, hcBlendAcc);
				}
			}
			return;
		}

		// entry / entryexit / add modes all follow state-change detection
		// via lastCardinality==-1.  To use this, set lastCardinality>=0
		// before each add.  If after the add it is -1, the state changed.
		final boolean entry=sample.equals("entry");
		final boolean entryexit=sample.equals("entryexit");
		final boolean addMode=sample.equals("add");

		// For entryexit we need to remember the PREVIOUS extracted state
		// so we can emit it as an "exit" sample just before the new entry.
		final int[] prevNlzBucketCount=new int[64];
		final int[][] prevNlzHbitSet=(hbits>0 ? new int[hbits][64] : null);
		long prevAdded=0;// cardinality at which prev state was extracted
		boolean havePrev=false;

		long added=0;
		while(added<maxCard){
			ct.lastCardinality=0;// sentinel: flips to -1 on state change
			ct.hashAndStore(rng.nextLong());
			added++;
			if(ct.lastCardinality!=-1){continue;}// no state change

			// State changed.  Extract current state.
			extractState(ct, trackerType, B, hbits, nlzBucketCount, nlzHbitSet);

			if(entryexit && havePrev){
				// Emit EXIT sample: prev state, valid from prevAdded to added-1.
				recordSample(added-1, B, hbits, tierScale, prevNlzBucketCount, prevNlzHbitSet,
					dlcAcc, beffAcc, unsAcc);
				if(doBlend){
					recordBlend(added-1, B, hbits, tierScale, prevNlzBucketCount, prevNlzHbitSet,
						errByV, errByBeff, errByUnseen, powerDlc, powerHc, dlcFormula, hcWeightMode, 1.0,
						dlcBlendAcc, hcBlendAcc);
				}
			}
			if(addMode && havePrev){
				// Dwell-weighted sample of previous state.
				final long weight=added-1-prevAdded;
				if(weight>0){
					recordSampleWeighted(prevAdded+1, B, hbits, tierScale,
						prevNlzBucketCount, prevNlzHbitSet, weight,
						dlcAcc, beffAcc, unsAcc);
					if(doBlend){
						recordBlend(prevAdded+1, B, hbits, tierScale,
							prevNlzBucketCount, prevNlzHbitSet,
							errByV, errByBeff, errByUnseen, powerDlc, powerHc, dlcFormula, hcWeightMode, (double)weight,
							dlcBlendAcc, hcBlendAcc);
					}
				}
			}
			if(entry || entryexit){
				recordSample(added, B, hbits, tierScale, nlzBucketCount, nlzHbitSet,
					dlcAcc, beffAcc, unsAcc);
				if(doBlend){
					recordBlend(added, B, hbits, tierScale, nlzBucketCount, nlzHbitSet,
						errByV, errByBeff, errByUnseen, powerDlc, powerHc, dlcFormula, hcWeightMode, 1.0,
						dlcBlendAcc, hcBlendAcc);
				}
			}

			// Save as prev for next iteration
			System.arraycopy(nlzBucketCount, 0, prevNlzBucketCount, 0, 64);
			if(nlzHbitSet!=null){
				for(int d=0; d<hbits; d++){
					System.arraycopy(nlzHbitSet[d], 0, prevNlzHbitSet[d], 0, 64);
				}
			}
			prevAdded=added;
			havePrev=true;
		}

		// Flush tail for add mode: last state valid until maxCard
		if(addMode && havePrev){
			final long weight=maxCard-prevAdded;
			if(weight>0){
				recordSampleWeighted(prevAdded+1, B, hbits, tierScale,
					prevNlzBucketCount, prevNlzHbitSet, weight,
					dlcAcc, beffAcc, unsAcc);
				if(doBlend){
					recordBlend(prevAdded+1, B, hbits, tierScale,
						prevNlzBucketCount, prevNlzHbitSet,
						errByV, errByBeff, errByUnseen, powerDlc, powerHc, dlcFormula, hcWeightMode, (double)weight,
						dlcBlendAcc, hcBlendAcc);
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       State Extraction       ----------------*/
	/*--------------------------------------------------------------*/

	/** Dispatcher: extract state from either tracker type. */
	static void extractState(CardinalityTracker ct, int trackerType, int B, int hbits,
			int[] nlzBucketCount, int[][] nlzHbitSet){
		if(trackerType==1){
			extractStateCDLL4((CompressedDynamicLogLog4)ct, B, hbits,
				nlzBucketCount, nlzHbitSet);
		}else{
			extractStateUDLL6((UltraDynamicLogLog6)ct, B, hbits,
				nlzBucketCount, nlzHbitSet);
		}
	}

	/** Decode UDLL6 buckets into per-absNlz counts and per-hbit set counts. */
	static void extractStateUDLL6(UltraDynamicLogLog6 ct, int B, int hbits,
			int[] nlzBucketCount, int[][] nlzHbitSet){
		java.util.Arrays.fill(nlzBucketCount, 0);
		if(nlzHbitSet!=null){
			for(int d=0; d<hbits; d++){java.util.Arrays.fill(nlzHbitSet[d], 0);}
		}
		final int globalNLZ=ct.getMinZeros()-1;
		for(int i=0; i<B; i++){
			final int reg=ct.getReg(i);
			int absNlz, histPattern;
			if(reg==0){
				absNlz=globalNLZ; // floor level; -1 when nothing seen
				histPattern=0;
			}else{
				final int nlzPart=reg>>>2;
				final int fullHist=reg&3;
				if(nlzPart>=HISTORY_MARGIN){
					absNlz=nlzPart-HISTORY_MARGIN+globalNLZ+1;
					histPattern=fullHist;
				}else{
					// Sub-margin bucket: at floor level, history preserved
					absNlz=globalNLZ;
					histPattern=fullHist;
				}
			}
			if(absNlz<0 || absNlz>=64){continue;}
			nlzBucketCount[absNlz]++;
			if(hbits>0){
				// UDLL6 stores 2 history bits; if hbits<2, keep only the
				// top (hbits) bits by right-shifting the low-order bits away.
				final int masked=(hbits>=2 ? histPattern : (histPattern>>>(2-hbits)));
				for(int d=0; d<hbits; d++){
					// Per CardStats.java:286: d=0 is top bit.
					if((masked&(1<<(hbits-1-d)))!=0){
						nlzHbitSet[d][absNlz]++;
					}
				}
			}
		}
	}

	/**
	 * Decode CDLL4 nibbles into per-absTier counts and per-hbit set counts.
	 * Nibble layout: bits [3:1]=tierPart (0=empty/floor-level, 1-7=relTier 0-6),
	 * bit [0]=1-bit history.  absTier = tierPart + globalNLZ.
	 */
	static void extractStateCDLL4(CompressedDynamicLogLog4 ct, int B, int hbits,
			int[] nlzBucketCount, int[][] nlzHbitSet){
		java.util.Arrays.fill(nlzBucketCount, 0);
		if(nlzHbitSet!=null){
			for(int d=0; d<hbits; d++){java.util.Arrays.fill(nlzHbitSet[d], 0);}
		}
		final int globalNLZ=ct.getMinZeros()-1;
		for(int i=0; i<B; i++){
			final int nib=ct.readNibble(i);
			final int tp=nib>>>1;
			final int hist=nib&1;
			int absTier=tp+globalNLZ;
			if(absTier<0){continue;} // truly empty
			if(absTier>=64){continue;}
			nlzBucketCount[absTier]++;
			if(hbits>0 && hist!=0){
				// CDLL4 has 1-bit history; place into d=0 slot.
				nlzHbitSet[0][absTier]++;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Sample Recording       ----------------*/
	/*--------------------------------------------------------------*/

	/** Compute per-tier DLC/HC estimates and accumulate into bins (weight=1). */
	static void recordSample(long trueCard, int B, int hbits, double tierScale,
			int[] nlzBucketCount, int[][] nlzHbitSet,
			double[][] dlcAcc, double[][] beffAcc, double[][] unsAcc){
		recordSampleWeighted(trueCard, B, hbits, tierScale, nlzBucketCount, nlzHbitSet,
			1.0, dlcAcc, beffAcc, unsAcc);
	}

	/** Weighted sample recording. */
	static void recordSampleWeighted(long trueCard, int B, int hbits, double tierScale,
			int[] nlzBucketCount, int[][] nlzHbitSet, double weight,
			double[][] dlcAcc, double[][] beffAcc, double[][] unsAcc){
		if(trueCard<=0){return;}

		// Find max active tier (highest absNlz with any buckets)
		int maxNlz=0;
		for(int t=63; t>=0; t--){
			if(nlzBucketCount[t]>0){maxNlz=t; break;}
		}
		final int maxTierHC=Math.min(maxNlz+2, Math.min(63, MAX_TIER));

		// V_t = B - cum(filled at absNlz >= t)
		// Iterate t from top down, keep running cumulative
		int cumFilled=0;
		for(int t=63; t>=0; t--){cumFilled+=nlzBucketCount[t];}
		// cumFilled now total filled.  For DLC we need V_t = B - (filled at >= t)
		// but the definition used in CardStats.dlcPure / paper.md:248 is
		// V_t = V + sum(n_i for i=0..t-1) = count of buckets at absNlz < t.
		// That is B minus (filled at absNlz >= t).
		int filledGE=cumFilled;// filled at absNlz >= 0 = total filled (includes floor-level; negatives skipped)
		// Actually nlzBucketCount[k] is count at absNlz==k, so filledGE at tier 0 = sum_{k>=0}.
		// We need filled at absNlz >= t for each t; iterate top-down.
		// Reset and rebuild incrementally.
		int aboveOrEq=0;
		int[] filledAtOrAbove=new int[65];
		for(int t=63; t>=0; t--){
			aboveOrEq+=nlzBucketCount[t];
			filledAtOrAbove[t]=aboveOrEq;
		}
		filledAtOrAbove[64]=0;

		for(int t=0; t<=maxTierHC; t++){
			// DLC(t): V_t = B - filledAtOrAbove[t]
			final int V=B-filledAtOrAbove[t];
			if(V>=1 && V<B){
				final double scale=Math.pow(2.0, t*tierScale);// 2^(t*S)
				final double est=scale*(double)B*Math.log((double)B/V);
				final double signed=(est-trueCard)/(double)trueCard;
				final double abs=Math.abs(signed);
				final int bin=Math.min(V, B);
				if(bin>=0 && bin<=B && abs>0 && !Double.isInfinite(abs)){
					addToBin(dlcAcc, bin, abs, signed, weight);
				}
			}

			// HC(t): only when hbits >= 1
			if(hbits>=1){
				int hcBeff=0, hcUnseen=0;
				for(int d=0; d<hbits; d++){
					final int src=t+d+1;
					if(src>=64){continue;}
					hcBeff+=nlzBucketCount[src];
					hcUnseen+=(nlzBucketCount[src]-nlzHbitSet[d][src]);
				}
				if(hcBeff>=8 && hcUnseen>=1 && hcUnseen<hcBeff){
					// HC scale: 2^((t+1)*S) / (2^S - 1).  The (2^S-1) divisor
					// accounts for the tier-probability ratio; for UDLL6 S=1
					// it equals 1 (vanishes), for CDLL4 S=1.5 it equals ~1.828.
					final double tierDenom=Math.pow(2.0, tierScale)-1.0;
					final double est=Math.pow(2.0, (t+1)*tierScale)*(double)B*
						Math.log((double)hcBeff/hcUnseen)/tierDenom;
					if(est>0 && !Double.isInfinite(est) && !Double.isNaN(est)){
						final double signed=(est-trueCard)/(double)trueCard;
						final double abs=Math.abs(signed);
						if(abs>0 && !Double.isInfinite(abs)){
							final int beffBin=Math.min(hcBeff, B);
							final int unsBin=Math.min(hcUnseen, B);
							addToBin(beffAcc, beffBin, abs, signed, weight);
							addToBin(unsAcc, unsBin, abs, signed, weight);
						}
					}
				}
			}
		}
	}

	private static final double SQRT_2_OVER_PI=Math.sqrt(2.0/Math.PI);

	/**
	 * Compute table- or formula-weighted DLC (and optionally HC) tier-blend
	 * estimates and record the blend's |err| in a log2(trueCard)-indexed
	 * accumulator.  Blend formula: geometric mean of per-tier estimates
	 * weighted by (1/err_t)^power, where err_t is either a table lookup
	 * (default) or the Mode 2 theoretical formula (when *Formula=true).
	 */
	static void recordBlend(long trueCard, int B, int hbits, double tierScale,
			int[] nlzBucketCount, int[][] nlzHbitSet,
			double[] errByV, double[] errByBeff, double[] errByUnseen,
			double powerDlc, double powerHc,
			boolean dlcFormula, int hcWeightMode, double weight,
			double[][] dlcBlendAcc, double[][] hcBlendAcc){
		if(trueCard<=0){return;}

		int[] filledAtOrAbove=new int[65];
		int acc=0;
		for(int t=63; t>=0; t--){acc+=nlzBucketCount[t]; filledAtOrAbove[t]=acc;}
		filledAtOrAbove[64]=0;

		int maxNlz=0;
		for(int t=63; t>=0; t--){
			if(nlzBucketCount[t]>0){maxNlz=t; break;}
		}
		final int maxTier=Math.min(maxNlz+2, Math.min(63, MAX_TIER));

		// DLC blend (table OR Mode 2 formula)
		double sumWDlc=0, sumWLogEstDlc=0;
		for(int t=0; t<=maxTier; t++){
			final int V=B-filledAtOrAbove[t];
			if(V<1 || V>=B){continue;}
			final double scale=Math.pow(2.0, t*tierScale);
			final double est=scale*(double)B*Math.log((double)B/V);
			if(est<=0 || Double.isInfinite(est) || Double.isNaN(est)){continue;}
			final double err;
			if(dlcFormula){
				// Mode 2: err = √(2/π) · √(occ/(B·V)) / ln(B/V)
				final int occ=B-V;
				err=SQRT_2_OVER_PI*Math.sqrt((double)occ/((double)B*V))/Math.log((double)B/V);
			}else{
				err=errByV[Math.min(V, B)];
			}
			if(!(err>0) || Double.isInfinite(err) || Double.isNaN(err)){continue;}
			final double w=Math.pow(1.0/err, powerDlc);
			sumWDlc+=w;
			sumWLogEstDlc+=w*Math.log(est);
		}
		if(sumWDlc>0){
			final double blend=Math.exp(sumWLogEstDlc/sumWDlc);
			final double signed=(blend-trueCard)/(double)trueCard;
			final double abs=Math.abs(signed);
			final int cbin=log2CardBin(trueCard);
			if(abs>0 && !Double.isInfinite(abs) && cbin>=0 && cbin<dlcBlendAcc[0].length){
				addToBin(dlcBlendAcc, cbin, abs, signed, weight);
			}
		}

		// HC blend (table OR Mode 2 formula; requires history)
		// hcWeightMode: 0=auto, 1=beff, 2=unseen, 3=geo(beff,unseen), 4=max(beff,unseen),
		//               5=V (DLC table), 6=formula
		final boolean doHc=(hbits>=1) && (hcWeightMode==6 ||
			(hcWeightMode==5 && errByV!=null) ||
			(hcWeightMode==1 && errByBeff!=null) ||
			(hcWeightMode==2 && errByUnseen!=null) ||
			((hcWeightMode==3 || hcWeightMode==4) && errByBeff!=null && errByUnseen!=null) ||
			(hcWeightMode==0 && (errByBeff!=null || errByUnseen!=null)));
		if(doHc){
			double sumWHc=0, sumWLogEstHc=0;
			for(int t=0; t<=maxTier; t++){
				int hcBeff=0, hcUnseen=0;
				for(int d=0; d<hbits; d++){
					final int src=t+d+1;
					if(src>=64){continue;}
					hcBeff+=nlzBucketCount[src];
					hcUnseen+=(nlzBucketCount[src]-nlzHbitSet[d][src]);
				}
				if(hcBeff<8 || hcUnseen<1 || hcUnseen>=hcBeff){continue;}
				final double tierDenom=Math.pow(2.0, tierScale)-1.0;
				final double est=Math.pow(2.0, (t+1)*tierScale)*(double)B*
					Math.log((double)hcBeff/hcUnseen)/tierDenom;
				if(est<=0 || Double.isInfinite(est) || Double.isNaN(est)){continue;}
				final double err;
				switch(hcWeightMode){
					case 6: {// Mode 2 formula
						final int hcOcc=hcBeff-hcUnseen;
						err=SQRT_2_OVER_PI*Math.sqrt((double)hcOcc/((double)hcBeff*hcUnseen))
							/Math.log((double)hcBeff/hcUnseen);
						break;
					}
					case 5: {// V-axis lookup: uses DLC's V_t for the same tier
						final int V=B-filledAtOrAbove[t];
						err=(V>=1 && V<=B ? errByV[V] : 0);
						break;
					}
					case 1: err=errByBeff[Math.min(hcBeff, B)]; break;
					case 2: err=errByUnseen[Math.min(hcUnseen, B)]; break;
					case 3: {// geo mean of marginals
						final double eb=errByBeff[Math.min(hcBeff, B)];
						final double eu=errByUnseen[Math.min(hcUnseen, B)];
						err=(eb>0 && eu>0 ? Math.sqrt(eb*eu) : Math.max(eb, eu));
						break;
					}
					case 4: {// max of marginals (worst-case blend)
						final double eb=errByBeff[Math.min(hcBeff, B)];
						final double eu=errByUnseen[Math.min(hcUnseen, B)];
						err=Math.max(eb, eu);
						break;
					}
					default: {// 0=auto
						if(errByBeff!=null && errByUnseen!=null){
							final double eb=errByBeff[Math.min(hcBeff, B)];
							final double eu=errByUnseen[Math.min(hcUnseen, B)];
							err=(eb>0 && eu>0 ? Math.sqrt(eb*eu) : Math.max(eb, eu));
						}else if(errByBeff!=null){err=errByBeff[Math.min(hcBeff, B)];}
						else{err=errByUnseen[Math.min(hcUnseen, B)];}
					}
				}
				if(!(err>0) || Double.isInfinite(err) || Double.isNaN(err)){continue;}
				final double w=Math.pow(1.0/err, powerHc);
				sumWHc+=w;
				sumWLogEstHc+=w*Math.log(est);
			}
			if(sumWHc>0){
				final double blend=Math.exp(sumWLogEstHc/sumWHc);
				final double signed=(blend-trueCard)/(double)trueCard;
				final double abs=Math.abs(signed);
				final int cbin=log2CardBin(trueCard);
				if(abs>0 && !Double.isInfinite(abs) && cbin>=0 && cbin<hcBlendAcc[0].length){
					addToBin(hcBlendAcc, cbin, abs, signed, weight);
				}
			}
		}
	}

	static boolean parseBool(String b){
		if(b==null){return false;}
		final String s=b.toLowerCase();
		return s.equals("t") || s.equals("true") || s.equals("1") || s.equals("yes");
	}

	static int log2CardBin(long card){
		if(card<=0){return 0;}
		return 63-Long.numberOfLeadingZeros(card);
	}

	/** Add one sample to accumulator bin. */
	static void addToBin(double[][] acc, int bin, double absErr, double signedErr,
			double weight){
		acc[0][bin]+=weight;// count
		acc[1][bin]+=weight*absErr;
		acc[2][bin]+=weight*Math.log(absErr);
		acc[3][bin]+=weight/absErr;
		acc[4][bin]+=weight*signedErr;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Utilities          ----------------*/
	/*--------------------------------------------------------------*/

	static long[] logThresh(long maxCard, int points){
		final long[] tmp=new long[points];
		for(int p=0; p<points; p++){
			final double frac=(double)p/(points-1);
			tmp[p]=Math.max(1, (long)Math.pow((double)maxCard, frac));
		}
		int unique=1;
		for(int i=1; i<points; i++){if(tmp[i]>tmp[i-1]){unique++;}}
		final long[] out=new long[unique];
		out[0]=tmp[0];
		int j=1;
		for(int i=1; i<points; i++){
			if(tmp[i]>tmp[i-1]){out[j++]=tmp[i];}
		}
		return out;
	}

	static void mergeInto(double[][] dst, double[][] src, int nBins){
		for(int s=0; s<5; s++){
			for(int i=0; i<nBins; i++){dst[s][i]+=src[s][i];}
		}
	}

	/**
	 * Load a V-axis (or hcBeff-axis) error table from a previous sim.
	 * File format: `#` comments, then data rows
	 * `axis\tcount\tlinAvgAbsErr\tgeoAvgAbsErr\tharmAvgAbsErr\tavgSignedErr`.
	 * Returns double[nBins] indexed by axis value; unpopulated bins are 0
	 * (caller guards against 0).  avgCol picks which column: lin|geo|harm.
	 */
	static double[] loadErrTable(String path, String avgCol, int nBins)
			throws IOException {
		final int colIdx;
		switch(avgCol.toLowerCase()){
			case "lin": case "linear":   colIdx=2; break;
			case "geo": case "geometric": colIdx=3; break;
			case "harm": case "harmonic": colIdx=4; break;
			default: throw new RuntimeException("Unknown avg='"+avgCol+"'");
		}
		final double[] table=new double[nBins];
		int loaded=0;
		try(BufferedReader br=new BufferedReader(new FileReader(path))){
			String line;
			while((line=br.readLine())!=null){
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				final String[] parts=line.split("\t");
				if(parts.length<=colIdx){continue;}
				try{
					final int bin=Integer.parseInt(parts[0]);
					if(bin<0 || bin>=nBins){continue;}
					final double v=Double.parseDouble(parts[colIdx]);
					if(v>0 && !Double.isNaN(v) && !Double.isInfinite(v)){
						table[bin]=v;
						loaded++;
					}
				}catch(NumberFormatException e){/* skip malformed */}
			}
		}
		// Fill gaps by left-to-right interpolation so every bin has a value
		double last=0;
		for(int i=0; i<nBins; i++){
			if(table[i]>0){last=table[i];}else if(last>0){table[i]=last;}
		}
		// Back-fill left edge
		double first=0;
		for(int i=0; i<nBins; i++){
			if(table[i]>0){first=table[i]; break;}
		}
		if(first>0){
			for(int i=0; i<nBins; i++){if(table[i]==0){table[i]=first;}}
		}
		System.err.println("Loaded "+loaded+" rows from "+path+" (col="+avgCol+")");
		return table;
	}

	static void writeTable(String path, String hdr, String axisName,
			double[][] acc, int nBins) throws FileNotFoundException {
		try(PrintStream out=new PrintStream(path)){
			out.println("# "+hdr+" axis="+axisName);
			out.println("#"+axisName+"\tcount\tlinAvgAbsErr\tgeoAvgAbsErr\tharmAvgAbsErr\tavgSignedErr");
			for(int i=0; i<nBins; i++){
				final double c=acc[0][i];
				if(c<=0){continue;}
				final double lin=acc[1][i]/c;
				final double geo=Math.exp(acc[2][i]/c);
				final double harm=c/acc[3][i];
				final double signed=acc[4][i]/c;
				out.printf("%d\t%.0f\t%.8f\t%.8f\t%.8f\t%.8f%n",
					i, c, lin, geo, harm, signed);
			}
		}
		System.err.println("Wrote "+path);
	}

	/**
	 * Single-scalar terminal-weighted mean |err|.  Sums (weight * |err|) across
	 * all bins and divides by total weight.  Under dwell-based weighting (our
	 * default for logcard + add) this is the cardinality-averaged mean error,
	 * biased toward high-cardinality / steady-state by construction.
	 */
	static double scalarWeightedAvg(double[][] acc){
		double sumAbs=0, sumW=0;
		for(int i=0; i<acc[0].length; i++){
			sumW+=acc[0][i];// count (= weight if weighted)
			sumAbs+=acc[1][i];// sum of weighted |err|
		}
		return (sumW>0 ? sumAbs/sumW : 0);
	}

	/** Summary for log2(card)-binned blend accumulators. */
	static void summarizeCard(PrintStream out, String label, double[][] acc, int nBins){
		final StringBuilder sb=new StringBuilder(label+" | ");
		// Print at selected log2 bins: 8, 12, 16, 20, 24
		int[] targets={8, 12, 16, 20, 24};
		for(int tgt : targets){
			if(tgt>=nBins){continue;}
			final double c=acc[0][tgt];
			if(c<=0){sb.append("  L2="+tgt+":---"); continue;}
			final double lin=acc[1][tgt]/c;
			final double geo=Math.exp(acc[2][tgt]/c);
			sb.append(String.format("  L2=%2d lin=%.2f%% geo=%.2f%%",
				tgt, 100*lin, 100*geo));
		}
		out.println(sb);
	}

	static void summarize(PrintStream out, String label, double[][] acc, int nBins){
		// Report err at bin B/4, B/2, 3B/4 if present; else nearest
		final int B=nBins-1;
		int[] targets={B/8, B/4, B/2, 3*B/4, 7*B/8};
		final StringBuilder sb=new StringBuilder(label+" | ");
		for(int tgt : targets){
			int best=-1;
			double bestDist=Double.MAX_VALUE;
			for(int i=0; i<nBins; i++){
				if(acc[0][i]<=0){continue;}
				final double d=Math.abs(i-tgt);
				if(d<bestDist){bestDist=d; best=i;}
			}
			if(best<0){sb.append("  V="+tgt+":---"); continue;}
			final double c=acc[0][best];
			final double lin=acc[1][best]/c;
			final double geo=Math.exp(acc[2][best]/c);
			sb.append(String.format("  bin=%4d lin=%.2f%% geo=%.2f%%",
				best, 100*lin, 100*geo));
		}
		out.println(sb);
	}
}
