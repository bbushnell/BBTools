package cardinality;

import java.io.PrintStream;

import fileIO.ByteFile;
import parse.Parse;
import parse.PreParser;
import rand.FastRandomXoshiro;

/**
 * TTLL4Simulator — TwinTailLogLog4 Simulator (4-bit tails).
 *
 * Simulates a single 12-bit TTLL4 register to build per-(tier, combined_h)
 * correction factor tables.
 *
 * Register layout:
 *   [9:6] = 4-bit exp
 *   [5:3] = h1 (3-bit tail for histBit==1)
 *   [2:0] = h0 (3-bit tail for histBit==0)
 *
 * Update rule (symmetric only):
 *   delta = hashNLZ - exp
 *   delta &lt; -2 : ignore
 *   delta in [-2, 0]: set bit (2+delta) of history[histBit]
 *   delta &gt; 0 : advance exp; shift both tails right by min(delta,3);
 *                set MSB of history[histBit]
 *   Overflow (newExp &gt; 15): cap at 15.
 *
 * State: tier=exp (0-15), combined_h=(h1&lt;&lt;3)|h0 (0-63).
 *
 * @author Brian, Chloe
 * @date April 26, 2026
 */
public class TTLL4Simulator {

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int NUM_TIERS=16;
	static final int TAIL_BITS=4;
	static final int TAIL_MASK=0xF;
	static final int EXP_SHIFT=2*TAIL_BITS; // 8
	static final int TAILS_MASK=(1<<EXP_SHIFT)-1; // 0xFF
	static final int NUM_COMBINED=1<<EXP_SHIFT; // 256

	static final int AVG_LIN=0;
	static final int AVG_GEO=1;
	static final int AVG_HARM=2;
	static final int AVG_BLEND=3;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public TTLL4Simulator(int iters_, int threads_, int maxTier_){
		iters=iters_;
		threads=threads_;
		maxTier=maxTier_;
		counts  =new long[NUM_TIERS][NUM_COMBINED];
		sums    =new double[NUM_TIERS][NUM_COMBINED];
		sumsLog =new double[NUM_TIERS][NUM_COMBINED];
		sumsInv =new double[NUM_TIERS][NUM_COMBINED];
		sumSq   =new double[NUM_TIERS][NUM_COMBINED];
		tierLogErr    =new double[NUM_TIERS];
		tierLogSigned =new double[NUM_TIERS];
		tierLinNum    =new double[NUM_TIERS];
		tierLinSigned =new double[NUM_TIERS];
		tierLinDen    =new double[NUM_TIERS];
		tierTransitions=new long[NUM_TIERS];
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	void simulate() throws InterruptedException {
		final boolean trackErrors=(loadedTable!=null);
		final long[][][]   threadCounts  =new long[threads][NUM_TIERS][NUM_COMBINED];
		final double[][][] threadSums    =new double[threads][NUM_TIERS][NUM_COMBINED];
		final double[][][] threadSumsLog =new double[threads][NUM_TIERS][NUM_COMBINED];
		final double[][][] threadSumsInv =new double[threads][NUM_TIERS][NUM_COMBINED];
		final double[][][] threadSumSq   =new double[threads][NUM_TIERS][NUM_COMBINED];
		final long[]   threadOverflow=new long[threads];
		final long[]   threadBelow  =new long[threads];
		final long[]   threadTotal  =new long[threads];
		final long[]   threadAdv    =new long[threads];
		final double[] tLogErr    =new double[threads];
		final double[] tLogSigned =new double[threads];
		final double[] tLinNum    =new double[threads];
		final double[] tLinSigned =new double[threads];
		final double[] tLinDen    =new double[threads];
		final long[]   tTrans     =new long[threads];
		final double[][] tTierLogErr    =new double[threads][NUM_TIERS];
		final double[][] tTierLogSigned =new double[threads][NUM_TIERS];
		final double[][] tTierLinNum    =new double[threads][NUM_TIERS];
		final double[][] tTierLinSigned =new double[threads][NUM_TIERS];
		final double[][] tTierLinDen    =new double[threads][NUM_TIERS];
		final long[][]   tTierTrans     =new long[threads][NUM_TIERS];

		final int base=iters/threads;
		final int rem =iters%threads;

		final Thread[] threadArr=new Thread[threads];
		for(int t=0; t<threads; t++){
			final int tid    =t;
			final int myIters=base+(tid<rem ? 1 : 0);

			threadArr[t]=new Thread(()->{
				final FastRandomXoshiro rng=new FastRandomXoshiro(tid+1);
				final long[][]   lCounts  =threadCounts[tid];
				final double[][] lSums    =threadSums[tid];
				final double[][] lSumsLog =threadSumsLog[tid];
				final double[][] lSumsInv =threadSumsInv[tid];
				final double[][] lSumSq   =threadSumSq[tid];
				long lOver=0, lBelow=0, lTotal=0, lAdv=0;
				double lLogErr=0, lLogSigned=0, lLinNum=0, lLinSigned=0, lLinDen=0;
				long   lTrans=0;
				final double[] lTierLogErr    =tTierLogErr[tid];
				final double[] lTierLogSigned =tTierLogSigned[tid];
				final double[] lTierLinNum    =tTierLinNum[tid];
				final double[] lTierLinSigned =tTierLinSigned[tid];
				final double[] lTierLinDen    =tTierLinDen[tid];
				final long[]   lTierTrans     =tTierTrans[tid];

				for(int iter=0; iter<myIters; iter++){
					int word=0;
					long eeMask=-1L;
					int prevTier=-1, prevCh=-1;

					for(long trueCard=1; trueCard<=10_000_000L; trueCard++){
						final long hash=TTLLSimulator.hash64shift(rng.nextLong());
						lTotal++;

						if(Long.compareUnsigned(hash, eeMask)>0){
							lBelow++;
						}else{
							final int hashNLZ=Long.numberOfLeadingZeros(hash);
							final int exp    =(word>>>EXP_SHIFT)&0xF;
							final int delta  =hashNLZ-exp;
							final int histBit=(int)(hash&1);

							if(delta<-(TAIL_BITS-1)){
								lBelow++;
							}else if(delta<=0){
								final int bitPos=(TAIL_BITS-1)+delta;
								final int regBitPos=(histBit==0) ? bitPos : (bitPos+TAIL_BITS);
								word|=(1<<regBitPos);
							}else{
								final int newExp=Math.min(exp+delta, 15);
								final int shiftAmt=Math.min(delta, TAIL_BITS);
								int h0=(word&TAIL_MASK)>>>shiftAmt;
								int h1=((word>>>TAIL_BITS)&TAIL_MASK)>>>shiftAmt;
								if(histBit==0){h0|=(1<<(TAIL_BITS-1));}
								else          {h1|=(1<<(TAIL_BITS-1));}
								word=(newExp<<EXP_SHIFT)|(h1<<TAIL_BITS)|h0;
								eeMask=(newExp<=(TAIL_BITS-1)) ? -1L : (-1L)>>>(newExp-(TAIL_BITS-1));
								lAdv++;
							}
						}

						final int tier=(word>>>EXP_SHIFT)&0xF;
						final int ch  =word&TAILS_MASK;
						lCounts[tier][ch]++;
						lSums[tier][ch]   +=trueCard;
						lSumsLog[tier][ch]+=Math.log(trueCard);
						lSumsInv[tier][ch]+=1.0/trueCard;
						lSumSq[tier][ch]  +=(double)trueCard*trueCard;

						if(trackErrors && (tier!=prevTier || ch!=prevCh)){
							prevTier=tier; prevCh=ch;
							if(trueCard>10){
								final int t2=Math.min(tier, loadedTable.length-1);
								final double est=loadedTable[t2][ch];
								if(est>0){
									final double diff=est-trueCard;
									final double relErr=Math.abs(diff)/trueCard;
									lLogErr    +=relErr;
									lLogSigned +=diff/trueCard;
									lLinNum    +=Math.abs(diff);
									lLinSigned +=diff;
									lLinDen    +=trueCard;
									lTrans++;
									if(tier<NUM_TIERS){
										lTierLogErr[tier]    +=relErr;
										lTierLogSigned[tier] +=diff/trueCard;
										lTierLinNum[tier]    +=Math.abs(diff);
										lTierLinSigned[tier] +=diff;
										lTierLinDen[tier]    +=trueCard;
										lTierTrans[tier]++;
									}
								}
							}
						}

						if(tier>maxTier){break;}
					}
				}

				threadOverflow[tid]=lOver;
				threadBelow[tid]   =lBelow;
				threadTotal[tid]   =lTotal;
				threadAdv[tid]     =lAdv;
				tLogErr[tid]    =lLogErr;
				tLogSigned[tid] =lLogSigned;
				tLinNum[tid]    =lLinNum;
				tLinSigned[tid] =lLinSigned;
				tLinDen[tid]    =lLinDen;
				tTrans[tid]     =lTrans;
			});
			threadArr[t].start();
		}

		for(final Thread th : threadArr){th.join();}

		for(int t=0; t<threads; t++){
			for(int tier=0; tier<NUM_TIERS; tier++){
				for(int s=0; s<NUM_COMBINED; s++){
					counts[tier][s]  +=threadCounts[t][tier][s];
					sums[tier][s]    +=threadSums[t][tier][s];
					sumsLog[tier][s] +=threadSumsLog[t][tier][s];
					sumsInv[tier][s] +=threadSumsInv[t][tier][s];
					sumSq[tier][s]   +=threadSumSq[t][tier][s];
				}
			}
			overflowCount+=threadOverflow[t];
			belowCount   +=threadBelow[t];
			totalAdds    +=threadTotal[t];
			advanceCount +=threadAdv[t];
			logErrSum    +=tLogErr[t];
			logSignedSum +=tLogSigned[t];
			linErrNum    +=tLinNum[t];
			linSignedNum +=tLinSigned[t];
			linErrDen    +=tLinDen[t];
			transitionCount+=tTrans[t];
			for(int tier=0; tier<NUM_TIERS; tier++){
				tierLogErr[tier]    +=tTierLogErr[t][tier];
				tierLogSigned[tier] +=tTierLogSigned[t][tier];
				tierLinNum[tier]    +=tTierLinNum[t][tier];
				tierLinSigned[tier] +=tTierLinSigned[t][tier];
				tierLinDen[tier]    +=tTierLinDen[t][tier];
				tierTransitions[tier]+=tTierTrans[t][tier];
			}
		}
	}

	void loadTable(String fname, String sectionName){
		ByteFile bf=ByteFile.makeByteFile(fname, false);
		final int tableTiers=maxTier+2;
		loadedTable=new double[tableTiers][NUM_COMBINED];
		boolean inSection=false;
		int loaded=0;
		for(byte[] raw=bf.nextLine(); raw!=null; raw=bf.nextLine()){
			if(raw.length==0){continue;}
			final String line=new String(raw).trim();
			if(line.startsWith(sectionName)){
				inSection=true;
				continue;
			}
			if(inSection){
				if(line.startsWith("#")){
					if(line.startsWith("#Tier")){continue;}
					break;
				}
				final String[] parts=line.split("\t");
				if(parts.length<2){continue;}
				final int tier;
				try{tier=Integer.parseInt(parts[0]);}
				catch(NumberFormatException e){continue;}
				if(tier>=tableTiers){continue;}
				for(int s=1; s<parts.length && (s-1)<NUM_COMBINED; s++){
					try{loadedTable[tier][s-1]=Double.parseDouble(parts[s]);}
					catch(NumberFormatException e){/* leave 0 */}
				}
				loaded++;
			}
		}
		bf.close();
		System.err.println("Loaded "+sectionName+" table: "+loaded
			+" tiers from "+fname);
	}

	void smoothSparseStates(long minObs){
		final long[][]   origCounts =new long[NUM_TIERS][NUM_COMBINED];
		final double[][] origSums   =new double[NUM_TIERS][NUM_COMBINED];
		final double[][] origSumsLog=new double[NUM_TIERS][NUM_COMBINED];
		final double[][] origSumsInv=new double[NUM_TIERS][NUM_COMBINED];
		final double[][] origSumSq  =new double[NUM_TIERS][NUM_COMBINED];
		for(int tier=0; tier<NUM_TIERS; tier++){
			System.arraycopy(counts[tier],  0, origCounts[tier],  0, NUM_COMBINED);
			System.arraycopy(sums[tier],    0, origSums[tier],    0, NUM_COMBINED);
			System.arraycopy(sumsLog[tier], 0, origSumsLog[tier], 0, NUM_COMBINED);
			System.arraycopy(sumsInv[tier], 0, origSumsInv[tier], 0, NUM_COMBINED);
			System.arraycopy(sumSq[tier],   0, origSumSq[tier],   0, NUM_COMBINED);
		}
		for(int tier=0; tier<NUM_TIERS; tier++){
			for(int ch=0; ch<NUM_COMBINED; ch++){
				if(origCounts[tier][ch]>=minObs){continue;}
				for(int bit=0; bit<EXP_SHIFT; bit++){
					final int nb=ch^(1<<bit);
					counts[tier][ch] +=origCounts[tier][nb];
					sums[tier][ch]   +=origSums[tier][nb];
					sumsLog[tier][ch]+=origSumsLog[tier][nb];
					sumsInv[tier][ch]+=origSumsInv[tier][nb];
					sumSq[tier][ch]  +=origSumSq[tier][nb];
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Build Methods          ----------------*/
	/*--------------------------------------------------------------*/

	double[] buildLinTierAvg(){
		final double[] a=new double[NUM_TIERS];
		for(int tier=0; tier<NUM_TIERS; tier++){
			long n=0; double s=0;
			for(int ch=0; ch<NUM_COMBINED; ch++){n+=counts[tier][ch]; s+=sums[tier][ch];}
			a[tier]=(n>0) ? s/n : 0;
		}
		return a;
	}

	double[] buildGeoTierAvg(){
		final double[] a=new double[NUM_TIERS];
		for(int tier=0; tier<NUM_TIERS; tier++){
			long n=0; double s=0;
			for(int ch=0; ch<NUM_COMBINED; ch++){n+=counts[tier][ch]; s+=sumsLog[tier][ch];}
			a[tier]=(n>0) ? Math.exp(s/n) : 0;
		}
		return a;
	}

	double[] buildHarmTierAvg(){
		final double[] a=new double[NUM_TIERS];
		for(int tier=0; tier<NUM_TIERS; tier++){
			long n=0; double s=0;
			for(int ch=0; ch<NUM_COMBINED; ch++){n+=counts[tier][ch]; s+=sumsInv[tier][ch];}
			a[tier]=(s>0) ? n/s : 0;
		}
		return a;
	}

	double[] buildBlendTierAvg(double[] geoAvg, double[] harmAvg){
		final double[] a=new double[NUM_TIERS];
		for(int tier=0; tier<NUM_TIERS; tier++){
			a[tier]=(2.0*harmAvg[tier]+geoAvg[tier])/3.0;
		}
		return a;
	}

	double[][] buildLinAvgTable(){
		final double[][] t=new double[NUM_TIERS][NUM_COMBINED];
		for(int tier=0; tier<NUM_TIERS; tier++){
			for(int ch=0; ch<NUM_COMBINED; ch++){
				t[tier][ch]=(counts[tier][ch]>0) ? sums[tier][ch]/counts[tier][ch] : 0;
			}
		}
		return t;
	}

	double[][] buildGeoAvgTable(){
		final double[][] t=new double[NUM_TIERS][NUM_COMBINED];
		for(int tier=0; tier<NUM_TIERS; tier++){
			for(int ch=0; ch<NUM_COMBINED; ch++){
				t[tier][ch]=(counts[tier][ch]>0)
					? Math.exp(sumsLog[tier][ch]/counts[tier][ch]) : 0;
			}
		}
		return t;
	}

	double[][] buildHarmAvgTable(){
		final double[][] t=new double[NUM_TIERS][NUM_COMBINED];
		for(int tier=0; tier<NUM_TIERS; tier++){
			for(int ch=0; ch<NUM_COMBINED; ch++){
				t[tier][ch]=(sumsInv[tier][ch]>0)
					? counts[tier][ch]/sumsInv[tier][ch] : 0;
			}
		}
		return t;
	}

	double[][] buildBlendAvgTable(double[][] geo, double[][] harm){
		final double[][] t=new double[NUM_TIERS][NUM_COMBINED];
		for(int tier=0; tier<NUM_TIERS; tier++){
			for(int ch=0; ch<NUM_COMBINED; ch++){
				t[tier][ch]=(2.0*harm[tier][ch]+geo[tier][ch])/3.0;
			}
		}
		return t;
	}

	double[][] buildMult(double[][] avgTable, double[] tierAvg){
		final double[][] m=new double[NUM_TIERS][NUM_COMBINED];
		for(int tier=0; tier<NUM_TIERS; tier++){
			for(int ch=0; ch<NUM_COMBINED; ch++){
				m[tier][ch]=(tierAvg[tier]>0 && avgTable[tier][ch]>0)
					? avgTable[tier][ch]/tierAvg[tier] : 1.0;
			}
		}
		return m;
	}

	double[][] buildProb(long[][] origCounts){
		final double[][] p=new double[NUM_TIERS][NUM_COMBINED];
		for(int tier=0; tier<NUM_TIERS; tier++){
			long total=0;
			for(int ch=0; ch<NUM_COMBINED; ch++){total+=origCounts[tier][ch];}
			if(total<=0){continue;}
			for(int ch=0; ch<NUM_COMBINED; ch++){
				p[tier][ch]=(double)origCounts[tier][ch]/total;
			}
		}
		return p;
	}

	double[][] buildCV(){
		final double[][] cv=new double[NUM_TIERS][NUM_COMBINED];
		for(int tier=0; tier<NUM_TIERS; tier++){
			for(int ch=0; ch<NUM_COMBINED; ch++){
				cv[tier][ch]=1.0;
				if(counts[tier][ch]>1){
					final double mean=sums[tier][ch]/counts[tier][ch];
					if(mean>0){
						final double var=sumSq[tier][ch]/counts[tier][ch]-mean*mean;
						cv[tier][ch]=(var>0) ? Math.sqrt(var)/mean : 0.001;
					}
				}
			}
		}
		return cv;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Output             ----------------*/
	/*--------------------------------------------------------------*/

	private static void printTable(PrintStream out, String name,
			double[][] table, int maxTier){
		out.printf("#%s\t%d%n", name, maxTier+1);
		final StringBuilder colHdr=new StringBuilder("#Tier");
		for(int s=0; s<NUM_COMBINED; s++){colHdr.append('\t').append(s);}
		out.println(colHdr);
		for(int tier=0; tier<=maxTier; tier++){
			final StringBuilder row=new StringBuilder(Integer.toString(tier));
			for(int s=0; s<NUM_COMBINED; s++){
				row.append('\t').append(String.format("%.8f", table[tier][s]));
			}
			out.println(row);
		}
		out.println("#");
	}

	private static void printCountTable(PrintStream out, long[][] origCounts,
			int maxTier){
		out.printf("#StateCounts\t%d%n", maxTier+1);
		final StringBuilder colHdr=new StringBuilder("#Tier");
		for(int s=0; s<NUM_COMBINED; s++){colHdr.append('\t').append(s);}
		out.println(colHdr);
		for(int tier=0; tier<=maxTier; tier++){
			final StringBuilder row=new StringBuilder(Integer.toString(tier));
			for(int s=0; s<NUM_COMBINED; s++){row.append('\t').append(origCounts[tier][s]);}
			out.println(row);
		}
		out.println("#");
	}

	void printResults(double[] linAvg, double[] geoAvg, double[] harmAvg,
			double[][] linAvgT, double[][] geoAvgT, double[][] harmAvgT,
			double[][] blendAvgT, double[][] linMult, double[][] geoMult,
			double[][] harmMult, double[][] blendMult, double[][] prob,
			double[][] cv, long[][] origCounts, PrintStream out,
			String cmdLine, long simMs){
		final int rows=maxTier+1;

		out.println("#TTLL4Simulator");
		out.println("#"+cmdLine);
		out.printf("#iters=%d\tthreads=%d\tmaxTier=%d%n", iters, threads, maxTier);
		out.printf("#states=%d\ttiers=%d\ttailBits=%d%n", NUM_COMBINED, rows, TAIL_BITS);
		out.printf("#simulationTime=%dms%n", simMs);
		out.printf("#totalAdds=%d%n", totalAdds);
		out.printf("#overflowRate=%.6f\tbelowRate=%.6f\tadvanceCount=%d%n",
			totalAdds>0 ? (double)overflowCount/totalAdds : 0,
			totalAdds>0 ? (double)belowCount/totalAdds : 0,
			advanceCount);
		out.println("#");

		out.printf("#TierAvg\t%d%n", rows);
		out.println("#Tier\tlin\tgeo\tharm\tobs");
		for(int tier=0; tier<=maxTier; tier++){
			long obs=0;
			for(int s=0; s<NUM_COMBINED; s++){obs+=origCounts[tier][s];}
			out.printf("%d\t%.8f\t%.8f\t%.8f\t%d%n",
				tier, linAvg[tier], geoAvg[tier], harmAvg[tier], obs);
		}
		out.println("#");

		printTable(out, "StateLinAvg",   linAvgT,   maxTier);
		printTable(out, "StateGeoAvg",   geoAvgT,   maxTier);
		printTable(out, "StateHarmAvg",  harmAvgT,  maxTier);
		printTable(out, "StateBlendAvg", blendAvgT, maxTier);
		printTable(out, "StateLinMult",   linMult,   maxTier);
		printTable(out, "StateGeoMult",   geoMult,   maxTier);
		printTable(out, "StateHarmMult",  harmMult,  maxTier);
		printTable(out, "StateBlendMult", blendMult, maxTier);
		printTable(out, "StateProb", prob, maxTier);
		printTable(out, "StateCV",   cv,   maxTier);
		printCountTable(out, origCounts, maxTier);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws Exception {
		{PreParser pp=new PreParser(args, null, false); args=pp.args;}
		int    iters    =10000;
		int    threads  =8;
		int    maxTier  =14;
		long   minObs   =100;
		String outFile  =null;
		String tableFile=null;
		int    avgMode  =AVG_LIN;

		for(final String arg : args){
			final String[] kv=arg.split("=", 2);
			if(kv.length!=2){throw new RuntimeException("Unknown parameter '"+arg+"'");}
			switch(kv[0]){
				case "iters":   iters  =Parse.parseIntKMG(kv[1]); break;
				case "threads": threads=Integer.parseInt(kv[1]); break;
				case "maxTier": maxTier=Integer.parseInt(kv[1]); break;
				case "minObs":  minObs =Parse.parseKMG(kv[1]);    break;
				case "out":     outFile  =kv[1];                 break;
				case "table":   tableFile=kv[1];                 break;
				case "avg":
					switch(kv[1].toLowerCase()){
						case "lin":   case "linear":   avgMode=AVG_LIN;   break;
						case "geo":   case "geometric":avgMode=AVG_GEO;   break;
						case "harm":  case "harmonic": avgMode=AVG_HARM;  break;
						case "blend":                  avgMode=AVG_BLEND; break;
						default: throw new RuntimeException("Unknown avg mode '"+kv[1]+"'");
					}
					break;
				default: throw new RuntimeException("Unknown parameter '"+arg+"'");
			}
		}

		final String cmdLine="iters="+iters+" threads="+threads
			+" maxTier="+maxTier+" minObs="+minObs
			+(outFile!=null ? " out="+outFile : " out=stdout")
			+(tableFile!=null ? " table="+tableFile+" avg="+TTLLSimulator.avgName(avgMode) : "");
		System.err.println("TTLL4Simulator: "+cmdLine);

		final long t0=System.currentTimeMillis();
		final TTLL4Simulator sim=new TTLL4Simulator(iters, threads, maxTier);
		if(tableFile!=null){sim.loadTable(tableFile, "#"+TTLLSimulator.avgName(avgMode));}
		sim.simulate();
		final long simMs=System.currentTimeMillis()-t0;
		System.err.println("Simulation time: "+simMs+" ms");

		final long[][] origCounts=new long[NUM_TIERS][NUM_COMBINED];
		for(int t=0; t<NUM_TIERS; t++){
			System.arraycopy(sim.counts[t], 0, origCounts[t], 0, NUM_COMBINED);
		}

		final double[] linAvg =sim.buildLinTierAvg();
		final double[] geoAvg =sim.buildGeoTierAvg();
		final double[] harmAvg=sim.buildHarmTierAvg();

		System.err.println();
		System.err.printf("%-6s %-14s %-14s %-14s %-12s%n",
			"Tier", "LinAvg", "GeoAvg", "HarmAvg", "Observations");
		for(int tier=0; tier<=maxTier; tier++){
			long obs=0;
			for(int s=0; s<NUM_COMBINED; s++){obs+=origCounts[tier][s];}
			final double growth=(tier>0 && linAvg[tier-1]>0)
				? linAvg[tier]/linAvg[tier-1] : Double.NaN;
			System.err.printf("%-6d %-14.2f %-14.2f %-14.2f %-12d (x%.4f)%n",
				tier, linAvg[tier], geoAvg[tier], harmAvg[tier], obs, growth);
		}

		if(tableFile!=null && sim.transitionCount>0){
			final long tc=sim.transitionCount;
			System.err.println();
			System.err.println("=== Error Tracking ("+TTLLSimulator.avgName(avgMode)+") ===");
			System.err.printf("transitions   =%d%n", tc);
			System.err.printf("logAvgErr     =%.4f%%  (mean |est-true|/true)%n",
				100.0*sim.logErrSum/tc);
			System.err.printf("logSignedErr  =%+.4f%%  (mean (est-true)/true)%n",
				100.0*sim.logSignedSum/tc);
			System.err.printf("linAvgErr     =%.4f%%  (sum|est-true|/sumTrue)%n",
				100.0*sim.linErrNum/sim.linErrDen);
			System.err.printf("linSignedErr  =%+.4f%%  (sum(est-true)/sumTrue)%n",
				100.0*sim.linSignedNum/sim.linErrDen);
			System.err.println();
			System.err.printf("%-6s %-12s %-12s %-12s %-10s%n",
				"Tier", "LogErr%", "LogSigned%", "LinErr%", "Transitions");
			for(int tier=0; tier<=maxTier; tier++){
				final long tt=sim.tierTransitions[tier];
				if(tt==0){continue;}
				System.err.printf("%-6d %-12.4f %-12.4f %-12.4f %-10d%n",
					tier,
					100.0*sim.tierLogErr[tier]/tt,
					100.0*sim.tierLogSigned[tier]/tt,
					100.0*sim.tierLinNum[tier]/sim.tierLinDen[tier],
					tt);
			}
		}

		sim.smoothSparseStates(minObs);
		final double[][] linAvgT  =sim.buildLinAvgTable();
		final double[][] geoAvgT  =sim.buildGeoAvgTable();
		final double[][] harmAvgT =sim.buildHarmAvgTable();
		final double[][] blendAvgT=sim.buildBlendAvgTable(geoAvgT, harmAvgT);
		final double[][] linMult  =sim.buildMult(linAvgT,   linAvg);
		final double[][] geoMult  =sim.buildMult(geoAvgT,   geoAvg);
		final double[][] harmMult =sim.buildMult(harmAvgT,  harmAvg);
		final double[][] blendMult=sim.buildMult(blendAvgT,
			sim.buildBlendTierAvg(geoAvg, harmAvg));
		final double[][] prob     =sim.buildProb(origCounts);
		final double[][] cv       =sim.buildCV();

		final PrintStream out=(outFile!=null)
			? new PrintStream(outFile) : System.out;
		sim.printResults(linAvg, geoAvg, harmAvg,
			linAvgT, geoAvgT, harmAvgT, blendAvgT,
			linMult, geoMult, harmMult, blendMult,
			prob, cv, origCounts, out, cmdLine, simMs);
		if(out!=System.out){out.close();}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	final int iters;
	final int threads;
	final int maxTier;

	long[][]   counts;
	double[][] sums;
	double[][] sumsLog;
	double[][] sumsInv;
	double[][] sumSq;

	long overflowCount, belowCount, totalAdds, advanceCount;

	double[][] loadedTable;
	double logErrSum, logSignedSum, linErrNum, linSignedNum, linErrDen;
	long   transitionCount;
	double[] tierLogErr, tierLogSigned, tierLinNum, tierLinSigned, tierLinDen;
	long[]   tierTransitions;

}
