package cardinality;

import shared.Tools;

/**
 * Single-bucket simulation for DDL8 (mantissa mode).
 * Measures average cardinality per (tier, mantissa_state) to derive
 * empirical correction factors, analogous to ULL8TierStats.
 * <p>
 * DDL8 encoding: stored = ((relNlz+1) &lt;&lt; 2) | invMantissa
 * where invMantissa = (~(key &gt;&gt;&gt; shift)) &amp; 3, giving 4 states per tier.
 * <p>
 * Unlike ULL8 (which has history carry on promotion), DDL8's mantissa
 * is purely local: it depends only on the current max hash, not on
 * previous tiers. So the mantissa resets on every promotion.
 * <p>
 * Usage: {@code java cardinality.DDL8TierStats [inner=32768] [outer=131072] [maxtier=12]}
 *
 * @author Brian Bushnell
 */
public class DDL8TierStats {

	/*--------------------------------------------------------------*/
	/*----------------        Constants            ----------------*/
	/*--------------------------------------------------------------*/

	static final int mantissaBits=2;
	static final int mantissaMask=(1<<mantissaBits)-1; // 3
	static final int wordlen=64;

	/*--------------------------------------------------------------*/
	/*----------------           Main              ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		int inner=32768;
		int outer=131072;
		int maxTier=12;

		for(String arg : args){
			final String[] ab=arg.split("=");
			final String a=ab[0].toLowerCase(), b=ab.length>1 ? ab[1] : "";
			if(a.equals("inner")){inner=Integer.parseInt(b);}
			else if(a.equals("outer")){outer=Integer.parseInt(b);}
			else if(a.equals("maxtier") || a.equals("mt")){maxTier=Integer.parseInt(b);}
		}

		System.err.println("DDL8TierStats: inner="+inner+" outer="+outer+" maxtier="+maxTier);
		final long t0=System.nanoTime();

		// count[absNlz][mantissa] and sum[absNlz][mantissa]
		// mantissa index = invMantissa (0-3)
		final long[][] count=new long[64][4];
		final double[][] sum=new double[64][4];

		for(int trial=0; trial<outer; trial++){
			final long seed=Tools.hash64shift((long)trial*1234567891L+54321);

			// Single bucket: track absNlz and invMantissa
			// No shared exponent (single bucket, no promotion needed)
			int maxNlz=-1;      // absolute NLZ of current max
			int maxInvMant=0;   // inverted mantissa of current max

			for(int card=1; card<=inner; card++){
				final long key=Tools.hash64shift(seed+card);
				final int nlz=Long.numberOfLeadingZeros(key);

				// Compute inverted mantissa: top 2 bits after NLZ, inverted
				final int shift=wordlen-nlz-mantissaBits-1;
				final int invMant=(shift<0) ? 0 : (int)((~(key>>>shift))&mantissaMask);

				// DDL8 comparison: higher NLZ wins; on same NLZ, higher invMantissa wins
				// (because inverted mantissa means higher = smaller hash = rarer)
				if(nlz>maxNlz || (nlz==maxNlz && invMant>maxInvMant)){
					maxNlz=nlz;
					maxInvMant=invMant;
				}

				// Record observation
				if(maxNlz>=0 && maxNlz<=maxTier){
					count[maxNlz][maxInvMant]++;
					sum[maxNlz][maxInvMant]+=card;
				}
			}
		}

		final double elapsed=(System.nanoTime()-t0)*1e-9;
		System.err.println(String.format("Elapsed: %.1fs", elapsed));

		// Header
		System.out.println("NLZ\tTotal\tAvgCard\t"+
			"P(m0)\tP(m1)\tP(m2)\tP(m3)\t"+
			"AvgC(m0)\tAvgC(m1)\tAvgC(m2)\tAvgC(m3)\t"+
			"Ratio(m0)\tRatio(m1)\tRatio(m2)\tRatio(m3)\t"+
			"CF(m0)\tCF(m1)\tCF(m2)\tCF(m3)");

		for(int nlz=0; nlz<=maxTier; nlz++){
			long total=0;
			double totalSum=0;
			for(int m=0; m<4; m++){
				total+=count[nlz][m];
				totalSum+=sum[nlz][m];
			}
			if(total<100){continue;}

			final double avgCard=totalSum/total;
			final StringBuilder sb=new StringBuilder();
			sb.append(nlz).append('\t').append(total).append('\t');
			sb.append(String.format("%.2f", avgCard));

			// Frequencies
			for(int m=0; m<4; m++){
				sb.append('\t').append(String.format("%.6f", count[nlz][m]/(double)total));
			}
			// Average cardinality per mantissa state
			for(int m=0; m<4; m++){
				if(count[nlz][m]>0){
					sb.append('\t').append(String.format("%.2f", sum[nlz][m]/count[nlz][m]));
				}else{sb.append("\tN/A");}
			}
			// Ratios
			for(int m=0; m<4; m++){
				if(count[nlz][m]>10){
					final double stateAvg=sum[nlz][m]/count[nlz][m];
					sb.append('\t').append(String.format("%.6f", stateAvg/avgCard));
				}else{sb.append("\tN/A");}
			}
			// CFs = log2(ratio)
			for(int m=0; m<4; m++){
				if(count[nlz][m]>10){
					final double stateAvg=sum[nlz][m]/count[nlz][m];
					final double ratio=stateAvg/avgCard;
					sb.append('\t').append(String.format("%+.6f", Math.log(ratio)/Math.log(2)));
				}else{sb.append("\tN/A");}
			}

			System.out.println(sb);
		}
	}
}
