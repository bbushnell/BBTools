package cardinality;

import shared.Tools;

/**
 * Single-bucket 16-bit simulation comparing sub-NLZ extraction methods.
 * Each run tests ONE mode at a specified bit depth. Run multiple times
 * to compare modes.
 *
 * 16-bit register: upper 6 bits = nlzStored (absNlz+1), lower 10 bits available.
 *
 * Modes:
 *   mantissa:  Inverted top N bits after leading 1-bit. (~(key>>>shift))&mask
 *   andtissa:  Top N bits of (hash & (hash<<2)) after leading 1-bit. NOT inverted.
 *   nlz2:      min(cap, NLZ(key << (NLZ(key)+1))). Double NLZ.
 *   history:   Shift-register tracking sawMinus1..sawMinusN.
 *   luck:      min(cap, maxNlz - secondMaxNlz - 1). Gap between top 2 NLZ values.
 *
 * Usage: java cardinality.MantissaCompare [mode=mantissa] [bits=2]
 *        [inner=32768] [outer=131072] [maxtier=11]
 *
 * @author Brian Bushnell, Chloe
 * @date March 2026
 */
public class MantissaCompare {

	static final int MODE_MANTISSA=0, MODE_ANDTISSA=1, MODE_NLZ2=2, MODE_HISTORY=3, MODE_LUCK=4;
	static final String[] MODE_NAMES={"Mantissa", "Andtissa", "NLZ2", "History", "Luck"};

	public static void main(String[] args){
		int inner=32768;
		int outer=131072;
		int maxTier=11;
		int mode=MODE_MANTISSA;
		int bits=2;

		for(String arg : args){
			final String[] ab=arg.split("=");
			final String a=ab[0].toLowerCase(), b=ab.length>1 ? ab[1] : "";
			if(a.equals("inner")){inner=Integer.parseInt(b);}
			else if(a.equals("outer")){outer=Integer.parseInt(b);}
			else if(a.equals("maxtier") || a.equals("mt")){maxTier=Integer.parseInt(b);}
			else if(a.equals("bits")){bits=Integer.parseInt(b);}
			else if(a.equals("mode")){
				if(b.equals("mantissa")){mode=MODE_MANTISSA;}
				else if(b.equals("andtissa")){mode=MODE_ANDTISSA;}
				else if(b.equals("nlz2")){mode=MODE_NLZ2;}
				else if(b.equals("history")){mode=MODE_HISTORY;}
				else if(b.equals("luck")){mode=MODE_LUCK;}
				else{throw new RuntimeException("Unknown mode: "+b);}
			}
		}

		final int numStates=1<<bits;
		final int cap=numStates-1; // max value for capped modes
		System.err.println("MantissaCompare: mode="+MODE_NAMES[mode]+" bits="+bits
			+" states="+numStates+" inner="+inner+" outer="+outer+" maxtier="+maxTier);
		final long t0=System.nanoTime();

		// [tier][state]
		final long[][] count=new long[64][numStates];
		final double[][] sum=new double[64][numStates];

		for(int trial=0; trial<outer; trial++){
			long seed=Tools.hash64shift((long)trial*1234567891L+54321);

			// Single bucket state (16-bit register: upper 6 = nlzStored, lower 10 = extra)
			char stored=0; // unsigned 16-bit
			int maxNlz=-1;
			int extra=0;

			// For luck: track second-best NLZ
			int luckSecond=-1;

			// For history: the extra field IS the history, carried via shift register
			// We track it inside 'extra' and rebuild 'stored' on changes

			for(int card=1; card<=inner; card++){
				final long key=Tools.hash64shift(seed+card);
				final int nlz=Long.numberOfLeadingZeros(key);

				if(mode==MODE_MANTISSA || mode==MODE_ANDTISSA || mode==MODE_NLZ2){
					// Hash-based modes: compute value from this hash
					int val=0;
					if(mode==MODE_MANTISSA){
						final int shift=63-nlz-bits;
						val=(shift<0) ? 0 : (int)((~(key>>>shift))&cap);
					}else if(mode==MODE_ANDTISSA){
						final long anded=key&(key<<2);
						final int shift=63-nlz-bits;
						val=(shift<0) ? 0 : (int)((anded>>>shift)&cap);
					}else{ // NLZ2
						final int consumed=nlz+1;
						if(consumed>=64){val=cap;}
						else{val=Math.min(cap, Long.numberOfLeadingZeros(key<<consumed));}
					}

					// Update: higher NLZ always wins; same NLZ: higher extra wins
					if(nlz>maxNlz){
						maxNlz=nlz;
						extra=val;
					}else if(nlz==maxNlz && val>extra){
						extra=val;
					}

				}else if(mode==MODE_HISTORY){
					final int newNlzStored=Math.min(nlz+1, 63);
					final int oldNlzStored=(stored>0) ? (stored>>>10) : 0;

					if(newNlzStored>oldNlzStored){
						// Promotion: shift register carry
						final int k=newNlzStored-oldNlzStored;
						final int oldHist=(stored>0) ? (stored&cap) : 0;
						final int carry=(stored>0) ? (1<<bits) : 0;
						extra=((oldHist|carry)>>k)&cap;
						maxNlz=newNlzStored-1;
					}else if(newNlzStored<oldNlzStored){
						// Active: set history bit for NLZ difference
						final int diff=oldNlzStored-newNlzStored;
						if(diff>=1 && diff<=bits){
							// Bit (bits-1) = sawMinus1 (MSB = nearer)
							// Bit 0 = sawMinusN (LSB = farthest)
							extra|=(1<<(bits-diff));
						}
					}
					// Same NLZ: no change

				}else if(mode==MODE_LUCK){
					// Luck = raw gap between top two all-time NLZ values.
					// State 0 = small gap = established. State cap = big gap = outlier.
					// Step 1: Update leaderboard.
					if(nlz>maxNlz){
						luckSecond=maxNlz;
						maxNlz=nlz;
					}else if(nlz>luckSecond){
						luckSecond=nlz;
					}
					// Step 2: Compute raw gap.
					if(maxNlz>=0 && luckSecond>=0){
						extra=Math.min(cap, maxNlz-luckSecond);
					}else{
						extra=cap; // no second yet = max gap = outlier
					}
				}

				// Update stored register
				if(maxNlz>=0){
					stored=(char)(((maxNlz+1)<<10)|extra);
				}

				// Record observation
				if(maxNlz>=0 && maxNlz<=maxTier && extra<numStates){
					count[maxNlz][extra]++;
					sum[maxNlz][extra]+=card;
				}
			}
		}

		final double elapsed=(System.nanoTime()-t0)*1e-9;
		System.err.println(String.format("Elapsed: %.1fs", elapsed));

		// Report per-tier
		StringBuilder header=new StringBuilder("Tier\tTotal");
		for(int s=0; s<numStates; s++){header.append("\tP("+s+")\tCF("+s+")");}
		System.out.println(header);

		for(int t=0; t<=maxTier; t++){
			long tierTotal=0;
			double tierSum=0;
			for(int s=0; s<numStates; s++){
				tierTotal+=count[t][s];
				tierSum+=sum[t][s];
			}
			if(tierTotal<100){continue;}
			double tierAvg=tierSum/tierTotal;

			StringBuilder sb=new StringBuilder();
			sb.append(t).append('\t').append(tierTotal);
			for(int s=0; s<numStates; s++){
				sb.append('\t').append(String.format("%.4f", count[t][s]/(double)tierTotal));
				if(count[t][s]>10){
					double stateAvg=sum[t][s]/count[t][s];
					double cf=Math.log(stateAvg/tierAvg)/Math.log(2);
					sb.append('\t').append(String.format("%+.4f", cf));
				}else{sb.append("\tN/A");}
			}
			System.out.println(sb);
		}

		// Observation-weighted summary for tiers 3-maxTier
		long grandTotal=0;
		double[] grandCF=new double[numStates];
		for(int t=3; t<=maxTier; t++){
			long tierTotal=0;
			double tierSum=0;
			for(int s=0; s<numStates; s++){
				tierTotal+=count[t][s];
				tierSum+=sum[t][s];
			}
			if(tierTotal<100) continue;
			double tierAvg=tierSum/tierTotal;
			for(int s=0; s<numStates; s++){
				if(count[t][s]>10){
					double stateAvg=sum[t][s]/count[t][s];
					double cf=Math.log(stateAvg/tierAvg)/Math.log(2);
					grandCF[s]+=cf*count[t][s];
				}
				grandTotal+=count[t][s];
			}
		}

		System.out.println("\n# "+MODE_NAMES[mode]+" "+bits+"-bit weighted CFs (tiers 3-"+maxTier+"):");
		StringBuilder cfLine=new StringBuilder("# CF=[");
		double minCF=Double.MAX_VALUE, maxCF=-Double.MAX_VALUE;
		for(int s=0; s<numStates; s++){
			double cf=grandCF[s]/(grandTotal>0?grandTotal:1);
			if(s>0) cfLine.append(", ");
			cfLine.append(String.format("%+.6f", cf));
			if(cf<minCF) minCF=cf;
			if(cf>maxCF) maxCF=cf;
		}
		cfLine.append(String.format("]  range=%.4f", maxCF-minCF));
		System.out.println(cfLine);
	}
}
