package cardinality;

import rand.FastRandomXoshiro;
import parse.Parse;
import parse.PreParser;
import shared.Tools;

/**
 * Multi-bucket per-state bias simulator for BDLL5 (banked 3-bit exponent
 * + 2-bit history, 6 buckets per 32-bit word).
 * <p>
 * Unlike {@link MantissaCompare2} which simulates a single 16-bit register
 * with no banking, this simulator runs a full BDLL5 tracker (default 1536
 * buckets = 256 words = 1 KB, the standard BDLL5 optimization size). Bank
 * promotion is a word-level operation that correlates 6 buckets, and
 * changes per-state collision dynamics vs plain UDLL6.
 * <p>
 * For each insertion, records the TARGET bucket's post-update (absNlz, hist)
 * and the current cardinality. Outputs:
 * <ul>
 *   <li>stderr: per-tier verbose breakdown (MantissaCompare2-style)
 *   <li>stdout: HSB table TSV ready for StateTable.loadHsbTable()
 * </ul>
 *
 * @author Chloe
 * @date April 2026
 */
public class BankedBiasSimulator {

	/*--------------------------------------------------------------*/
	/*----------------            Main              ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws InterruptedException {
		{PreParser pp=new PreParser(args, null, false); args=pp.args;}
		int buckets=1536, inner=4000000, outer=1024;
		int maxTier=11, threads=1, ssTier=11;
		int avgMode=AVG_GEO;
		float promoteFrac=0.004f;

		for(String arg : args){
			final String[] ab=arg.split("=");
			final String a=ab[0].toLowerCase(), b=ab.length>1 ? ab[1] : "";
			if(a.equals("buckets") || a.equals("b")){buckets=Parse.parseIntKMG(b);}
			else if(a.equals("inner") || a.equals("iters")){inner=Parse.parseIntKMG(b);}
			else if(a.equals("outer") || a.equals("trials")){outer=Parse.parseIntKMG(b);}
			else if(a.equals("maxtier") || a.equals("mt")){maxTier=Integer.parseInt(b);}
			else if(a.equals("sstier") || a.equals("ss")){ssTier=Integer.parseInt(b);}
			else if(a.equals("threads") || a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("pf") || a.equals("promotefrac")){promoteFrac=Float.parseFloat(b);}
			else if(a.equals("avg")){
				if(b.equals("lin")){avgMode=AVG_LIN;}
				else if(b.equals("geo")){avgMode=AVG_GEO;}
				else if(b.equals("harm")){avgMode=AVG_HARM;}
				else{throw new RuntimeException("Unknown avg mode '"+b+"'");}
			}
			else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}

		final int numStates=4; // 2-bit history only (4 states)
		System.err.println("BankedBiasSimulator: buckets="+buckets+" promoteFrac="+promoteFrac
				+" inner="+inner+" outer="+outer+" maxtier="+maxTier
				+" ssTier="+ssTier+" threads="+threads+" avg="+AVG_NAMES[avgMode]);
		final long t0=System.nanoTime();

		// Thread-local accumulators
		final long[][][] tCount=new long[threads][64][numStates];
		final double[][][] tSum=new double[threads][64][numStates];
		final double[][][] tGeo=new double[threads][64][numStates];
		final double[][][] tHarm=new double[threads][64][numStates];

		final int fInner=inner, fMaxTier=maxTier;
		final int fBuckets=buckets;
		final float fPF=promoteFrac;
		final int fOuter=outer;

		final Thread[] workers=new Thread[threads];

		for(int ti=0; ti<threads; ti++){
			final int threadIdx=ti;
			final int trialStart=ti*(fOuter/threads);
			final int trialEnd=(ti==threads-1) ? fOuter : (ti+1)*(fOuter/threads);
			workers[ti]=new Thread(()->{
				final long[][] lc=tCount[threadIdx];
				final double[][] ls=tSum[threadIdx];
				final double[][] lg=tGeo[threadIdx];
				final double[][] lh=tHarm[threadIdx];

				final int words=Math.max(1, (fBuckets+3)/6);
				final int modBuckets=words*6;
				final int[] maxArray=new int[words];
				final int promoteThreshold=(fPF>0 ? (int)(modBuckets*fPF) : 0);

				for(int trial=trialStart; trial<trialEnd; trial++){
					final FastRandomXoshiro rng=new FastRandomXoshiro(trial*0x9E3779B97F4A7C15L+0x1234567L);
					java.util.Arrays.fill(maxArray, 0);
					int globalNLZ=-1;
					int minZeroCount=modBuckets;
					long eeMask=-1L;

					for(int card=1; card<=fInner; card++){
						final long number=rng.nextLong();
						final long key=Tools.hash64shift(number);
						if(Long.compareUnsigned(key, eeMask)>0){continue;}
						final int nlz=Long.numberOfLeadingZeros(key);
						final int bucket=(int)(Long.remainderUnsigned(key, modBuckets));
						final int wordIdx=bucket/6;
						final int bucketShift=(bucket%6)*5;

						int word=maxArray[wordIdx];
						int bank=(word>>>30)&0x3;
						int localRelNlz=nlz-globalNLZ-1-bank;

						// Single-attempt bank promotion (matches BDLL3 and fixed
						// tracker). Cascading breaks minZeroCount tracking.
						if(localRelNlz>=7 && bank<3){
							boolean canPromote=true;
							for(int bb=0; bb<6; bb++){
								if(((word>>>(bb*5))&0x7)==0){canPromote=false; break;}
							}
							if(canPromote){
								int result=0;
								for(int bb=0; bb<6; bb++){
									final int shift=bb*5;
									final int nibble=(word>>>shift)&0x1F;
									int stored=nibble&0x7;
									final int hist=nibble&0x18;
									stored--;
									result|=((stored|hist)<<shift);
								}
								maxArray[wordIdx]=result|((bank+1)<<30);
								word=maxArray[wordIdx];
								bank++;
								localRelNlz=nlz-globalNLZ-1-bank;
							}
						}

						final int newStored=Math.min(localRelNlz+1, 7);
						final int oldNibble=(word>>>bucketShift)&0x1F;
						final int oldStored=oldNibble&0x7;
						final int oldHist=(oldNibble>>>3)&0x3;

						if(newStored<=0){
							if(oldStored>0){
								final int delta=oldStored-newStored;
								if(delta==1 || delta==2){
									final int bit=1<<(2-delta);
									final int nh=oldHist|bit;
									if(nh!=oldHist){
										final int nibble=(oldStored&0x7)|((nh&0x3)<<3);
										maxArray[wordIdx]=(maxArray[wordIdx]&~(0x1F<<bucketShift))|(nibble<<bucketShift);
									}
								}
							}else if(globalNLZ+bank>=0){
								// Floor-level hist update (mirrors tracker fix)
								final int delta=-1-localRelNlz;
								if(delta==1 || delta==2){
									final int bit=1<<(2-delta);
									final int nh=oldHist|bit;
									if(nh!=oldHist){
										final int nibble=(0&0x7)|((nh&0x3)<<3);
										maxArray[wordIdx]=(maxArray[wordIdx]&~(0x1F<<bucketShift))|(nibble<<bucketShift);
									}
								}
							}
						}else if(newStored==oldStored){
							// no-op
						}else if(newStored<oldStored){
							final int delta=oldStored-newStored;
							if(delta==1 || delta==2){
								final int bit=1<<(2-delta);
								final int nh=oldHist|bit;
								if(nh!=oldHist){
									final int nibble=(oldStored&0x7)|((nh&0x3)<<3);
									maxArray[wordIdx]=(maxArray[wordIdx]&~(0x1F<<bucketShift))|(nibble<<bucketShift);
								}
							}
						}else{
							// Advance: newStored > oldStored
							final int newHist;
							if(oldStored==0){
								if(globalNLZ+bank>=0){
									// Floor-level→real: carry-shift from floor-level hist (delta=newStored)
									newHist=(newStored>=3) ? 0 : (((1<<2)|oldHist)>>>newStored)&0x3;
								}else{
									newHist=0;
								}
							}else{
								final int delta=newStored-oldStored;
								newHist=(delta>=3) ? 0 : (((1<<2)|oldHist)>>>delta)&0x3;
							}
							final int nibble=(newStored&0x7)|((newHist&0x3)<<3);
							maxArray[wordIdx]=(maxArray[wordIdx]&~(0x1F<<bucketShift))|(nibble<<bucketShift);

							// Global promotion trigger (EARLY_PROMOTE)
							if(oldStored==0 && bank==0){
								if(--minZeroCount<=promoteThreshold){
									while(minZeroCount<=promoteThreshold && globalNLZ<64){
										globalNLZ++;
										final int exitThreshold=Math.max(0, globalNLZ-1);
										eeMask=(exitThreshold==0) ? -1L : -1L>>>exitThreshold;
										int newMinZeroCount=0;
										for(int w=0; w<words; w++){
											int bb=(maxArray[w]>>>30)&0x3;
											if(bb>0){
												final int newBank=bb-1;
												maxArray[w]=(maxArray[w]&0x3FFFFFFF)|((newBank&0x3)<<30);
												if(newBank==0){
													for(int bi=0; bi<6; bi++){
														int s=(maxArray[w]>>>(bi*5))&0x7;
														if(s==0){newMinZeroCount++;}
													}
												}
											}else{
												int w2=maxArray[w];
												if((w2&0x3FFFFFFF)==0){
													newMinZeroCount+=6;
													continue;
												}
												int result=0;
												for(int bi=0; bi<6; bi++){
													final int shift=bi*5;
													final int nibble2=(w2>>>shift)&0x1F;
													int s=nibble2&0x7;
													final int hist=nibble2&0x18;
													if(s>0){
														s--;
														if(s==0){newMinZeroCount++;}
													}
													result|=((s|hist)<<shift);
												}
												maxArray[w]=result;
											}
										}
										minZeroCount=newMinZeroCount;
									}
								}
							}
						}

						// Record a UNIFORMLY RANDOM bucket's current state (not the
						// hash target). Target-only sampling biases toward "just
						// updated" buckets (low-hist) and undersamples long-settled
						// buckets (high-hist). Uniform sampling removes that bias.
						final int sampleBucket=(int)((rng.nextLong()>>>1)%modBuckets);
						final int sampleShift=(sampleBucket%6)*5;
						final int sampleWordIdx=sampleBucket/6;
						final int curWord=maxArray[sampleWordIdx];
						final int curBank=(curWord>>>30)&0x3;
						final int curNibble=(curWord>>>sampleShift)&0x1F;
						final int curStored=curNibble&0x7;
						final int curHist=(curNibble>>>3)&0x3;
						final int absNlz=curStored+globalNLZ+curBank;
						if(absNlz<0){continue;} // truly empty
						final int emitHist=(absNlz==0) ? 0 : curHist;
						if(absNlz<=fMaxTier){
							lc[absNlz][emitHist]++;
							ls[absNlz][emitHist]+=card;
							lg[absNlz][emitHist]+=Math.log(card);
							lh[absNlz][emitHist]+=1.0/card;
						}
					}
				}
			});
			workers[ti].start();
		}
		for(Thread w : workers){w.join();}

		// Merge thread arrays
		final long[][] count=new long[64][numStates];
		final double[][] sum=new double[64][numStates];
		final double[][] geoSum=new double[64][numStates];
		final double[][] harmSum=new double[64][numStates];
		for(int ti=0; ti<threads; ti++){
			for(int t=0; t<64; t++){
				for(int s=0; s<numStates; s++){
					count[t][s]+=tCount[ti][t][s];
					sum[t][s]+=tSum[ti][t][s];
					geoSum[t][s]+=tGeo[ti][t][s];
					harmSum[t][s]+=tHarm[ti][t][s];
				}
			}
		}

		final double elapsed=(System.nanoTime()-t0)*1e-9;
		System.err.println(String.format("Elapsed: %.1fs", elapsed));

		// Verbose per-tier report to stderr
		StringBuilder header=new StringBuilder("Tier\tTotal\tAvgCard\tRatio");
		for(int s=0; s<numStates; s++){header.append("\tP("+s+")\tLinCF("+s+")\tGeoCF("+s+")\tHarmCF("+s+")");}
		System.err.println(header);

		double prevAvg=0;
		// Cache per-tier CFs for later extraction to HSB table
		final double[][] perTierLinCF=new double[maxTier+1][numStates];
		final double[][] perTierGeoCF=new double[maxTier+1][numStates];
		final double[][] perTierHarmCF=new double[maxTier+1][numStates];

		for(int t=0; t<=maxTier; t++){
			long tierTotal=0; double tierSum=0, tierGeo=0, tierHarm=0;
			for(int s=0; s<numStates; s++){
				tierTotal+=count[t][s]; tierSum+=sum[t][s];
				tierGeo+=geoSum[t][s]; tierHarm+=harmSum[t][s];
			}
			if(tierTotal<100){continue;}
			final double tierLinAvg=tierSum/tierTotal;
			final double tierGeoAvg=Math.exp(tierGeo/tierTotal);
			final double tierHarmAvg=tierTotal/tierHarm;
			StringBuilder sb=new StringBuilder();
			sb.append(t).append('\t').append(tierTotal);
			sb.append('\t').append(String.format("%.2f", tierLinAvg));
			sb.append('\t').append(prevAvg>0 ? String.format("%.4f", tierLinAvg/prevAvg) : "-");
			for(int s=0; s<numStates; s++){
				sb.append('\t').append(String.format("%.6f", count[t][s]/(double)tierTotal));
				if(count[t][s]>10){
					final double stateLinAvg=sum[t][s]/count[t][s];
					final double stateGeoAvg=Math.exp(geoSum[t][s]/count[t][s]);
					final double stateHarmAvg=count[t][s]/harmSum[t][s];
					final double linCF=Math.log(stateLinAvg/tierLinAvg)/Math.log(2);
					final double geoCF=Math.log(stateGeoAvg/tierGeoAvg)/Math.log(2);
					final double harmCF=Math.log(stateHarmAvg/tierHarmAvg)/Math.log(2);
					sb.append('\t').append(String.format("%+.8f", linCF));
					sb.append('\t').append(String.format("%+.8f", geoCF));
					sb.append('\t').append(String.format("%+.8f", harmCF));
					perTierLinCF[t][s]=linCF;
					perTierGeoCF[t][s]=geoCF;
					perTierHarmCF[t][s]=harmCF;
				}else{
					sb.append("\tN/A\tN/A\tN/A");
				}
			}
			System.err.println(sb);
			prevAvg=tierLinAvg;
		}

		// HSB table TSV to stdout (ready for StateTable.loadHsbTable)
		System.out.println("# BDLL5 2-bit history biases from BankedBiasSimulator");
		System.out.println("# buckets="+buckets+" promoteFrac="+promoteFrac+" inner="+inner+" outer="+outer
				+" avg="+AVG_NAMES[avgMode]+" ssTier="+ssTier);
		System.out.println("# Source: tiers 0-2 per-tier rows; ss = tier "+ssTier);

		final double[][] chosenCF=(avgMode==AVG_LIN) ? perTierLinCF
				: (avgMode==AVG_GEO) ? perTierGeoCF : perTierHarmCF;

		for(int t=0; t<3; t++){
			StringBuilder sb=new StringBuilder();
			sb.append(t);
			for(int s=0; s<numStates; s++){
				sb.append('\t').append(String.format("%+.8f", chosenCF[t][s]));
			}
			System.out.println(sb);
		}

		// ss row = ssTier values
		StringBuilder ssSb=new StringBuilder("ss");
		for(int s=0; s<numStates; s++){
			ssSb.append('\t').append(String.format("%+.8f", chosenCF[ssTier][s]));
		}
		System.out.println(ssSb);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Linear averaging mode: mean cardinality per state. */
	static final int AVG_LIN=0;
	/** Geometric averaging mode: exp(mean(log(card))) per state. */
	static final int AVG_GEO=1;
	/** Harmonic averaging mode: n/sum(1/card) per state. */
	static final int AVG_HARM=2;
	/** Display names for the three averaging modes. */
	static final String[] AVG_NAMES={"lin", "geo", "harm"};

}
