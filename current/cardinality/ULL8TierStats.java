package cardinality;

import shared.Tools;

/**
 * Single-bucket simulation to measure empirical correction factors.
 * Outputs per-tier, per-state average cardinality ratios and CFs.
 * <p>
 * Designed to be re-run with increasing iterations for more precision.
 * Ignore tiers above maxtier (where observation counts stop doubling).
 * <p>
 * Usage: {@code java cardinality.ULL8TierStats [inner=32768] [outer=65536] [maxtier=12]}
 *
 * @author Brian Bushnell
 */
public class ULL8TierStats {

	/*--------------------------------------------------------------*/
	/*----------------           Main              ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		int inner=32768;
		int outer=65536;
		int maxTier=12;

		for(String arg : args){
			final String[] ab=arg.split("=");
			final String a=ab[0].toLowerCase(), b=ab.length>1 ? ab[1] : "";
			if(a.equals("inner")){inner=Integer.parseInt(b);}
			else if(a.equals("outer")){outer=Integer.parseInt(b);}
			else if(a.equals("maxtier") || a.equals("mt")){maxTier=Integer.parseInt(b);}
		}

		System.err.println("ULL8TierStats: inner="+inner+" outer="+outer+" maxtier="+maxTier);
		final long t0=System.nanoTime();

		final long[][] count=new long[64][4];
		final double[][] sum=new double[64][4];

		// Disable ULL corrections for measurement
		UltraLogLog8.STATE_POWER=0.0;

		for(int trial=0; trial<outer; trial++){
			long seed=Tools.hash64shift((long)trial*1234567891L+12345);
			int stored=0;

			for(int card=1; card<=inner; card++){
				final long key=Tools.hash64shift(seed+card);
				final int nlz=Long.numberOfLeadingZeros(key);
				final int newNlzStored=Math.min(nlz+1, 63);
				final int oldNlzStored=(stored>0) ? (stored>>2) : 0;

				if(newNlzStored>oldNlzStored){
					final int k=newNlzStored-oldNlzStored;
					final int oldHist=(stored>0) ? (stored&0x3) : 0;
					final int carry=(stored>0) ? 0x4 : 0;
					final int newHist=((oldHist|carry)>>k)&0x3;
					stored=(newNlzStored<<2)|newHist;
				}else if(newNlzStored<oldNlzStored){
					final int diff=oldNlzStored-newNlzStored;
					if(diff==1){stored|=0x02;}
					else if(diff==2){stored|=0x01;}
				}

				if(stored>0){
					final int absNlz=(stored>>2)-1;
					final int state=stored&0x3;
					if(absNlz>=0 && absNlz<=maxTier){
						count[absNlz][state]++;
						sum[absNlz][state]+=card;
					}
				}
			}
		}

		final double elapsed=(System.nanoTime()-t0)*1e-9;
		System.err.println(String.format("Elapsed: %.1fs", elapsed));

		// Header
		System.out.println("NLZ\tTotal\tAvgCard\t"+
			"P(00)\tP(01)\tP(10)\tP(11)\t"+
			"AvgC(00)\tAvgC(01)\tAvgC(10)\tAvgC(11)\t"+
			"Ratio(00)\tRatio(01)\tRatio(10)\tRatio(11)\t"+
			"CF(00)\tCF(01)\tCF(10)\tCF(11)");

		for(int nlz=0; nlz<=maxTier; nlz++){
			long total=0;
			double totalSum=0;
			for(int s=0; s<4; s++){
				total+=count[nlz][s];
				totalSum+=sum[nlz][s];
			}
			if(total<100){continue;}

			final double avgCard=totalSum/total;
			final StringBuilder sb=new StringBuilder();
			sb.append(nlz).append('\t').append(total).append('\t');
			sb.append(String.format("%.2f", avgCard));

			for(int s=0; s<4; s++){
				sb.append('\t').append(String.format("%.6f", count[nlz][s]/(double)total));
			}
			for(int s=0; s<4; s++){
				if(count[nlz][s]>0){
					sb.append('\t').append(String.format("%.2f", sum[nlz][s]/count[nlz][s]));
				}else{sb.append("\tN/A");}
			}
			for(int s=0; s<4; s++){
				if(count[nlz][s]>10){
					sb.append('\t').append(String.format("%.6f", (sum[nlz][s]/count[nlz][s])/avgCard));
				}else{sb.append("\tN/A");}
			}
			for(int s=0; s<4; s++){
				if(count[nlz][s]>10){
					final double ratio=(sum[nlz][s]/count[nlz][s])/avgCard;
					sb.append('\t').append(String.format("%+.6f", Math.log(ratio)/Math.log(2)));
				}else{sb.append("\tN/A");}
			}
			System.out.println(sb);
		}
	}
}
