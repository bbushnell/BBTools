package cardinality;

import shared.Tools;

/**
 * Single-bucket simulation for ProtoLogLog16c.
 * Measures avg cardinality per (tier, extra_state) for any mode/combination.
 * Supports up to 10 extra bits = 1024 states.
 *
 * Usage: java cardinality.ProtoTierStats [inner=32768] [outer=131072] [maxtier=11]
 *        [mode=mantissa|nlz2|history|luck|andtissa]
 *        [mbits=2] [nbits=2] [hbits=2] [lbits=2] [abits=2]
 *
 * Can combine modes: mode=mantissa+history
 */
public class ProtoTierStats {

	public static void main(String[] args){
		int inner=32768;
		int outer=131072;
		int maxTier=11;

		for(String arg : args){
			final String[] ab=arg.split("=");
			final String a=ab[0].toLowerCase(), b=ab.length>1 ? ab[1] : "";
			if(a.equals("inner")){inner=Integer.parseInt(b);}
			else if(a.equals("outer")){outer=Integer.parseInt(b);}
			else if(a.equals("maxtier") || a.equals("mt")){maxTier=Integer.parseInt(b);}
			else if(a.equals("mode")){
				ProtoLogLog16c.MODE=0;
				for(String m : b.split("\\+")){
					if(m.equals("mantissa")){ProtoLogLog16c.MODE|=ProtoLogLog16c.MODE_MANTISSA;}
					else if(m.equals("nlz2")){ProtoLogLog16c.MODE|=ProtoLogLog16c.MODE_NLZ2;}
					else if(m.equals("history")){ProtoLogLog16c.MODE|=ProtoLogLog16c.MODE_HISTORY;}
					else if(m.equals("luck")){ProtoLogLog16c.MODE|=ProtoLogLog16c.MODE_LUCK;}
					else if(m.equals("andtissa")){ProtoLogLog16c.MODE|=ProtoLogLog16c.MODE_ANDTISSA;}
					else if(m.equals("none")){ProtoLogLog16c.MODE=0;}
				}
			}
			else if(a.equals("mbits")){ProtoLogLog16c.MANTISSA_BITS=Integer.parseInt(b);}
			else if(a.equals("nbits")){ProtoLogLog16c.NLZ2_BITS=Integer.parseInt(b);}
			else if(a.equals("hbits")){ProtoLogLog16c.HISTORY_BITS=Integer.parseInt(b);}
			else if(a.equals("lbits")){ProtoLogLog16c.LUCK_BITS=Integer.parseInt(b);}
			else if(a.equals("abits")){ProtoLogLog16c.ANDTISSA_BITS=Integer.parseInt(b);}
			else if(a.equals("nlzbits")){ProtoLogLog16c.NLZ_BITS=Integer.parseInt(b);}
		}

		// Calculate total extra bits used
		int totalExtra=0;
		if((ProtoLogLog16c.MODE & ProtoLogLog16c.MODE_MANTISSA)!=0) totalExtra+=ProtoLogLog16c.MANTISSA_BITS;
		if((ProtoLogLog16c.MODE & ProtoLogLog16c.MODE_ANDTISSA)!=0) totalExtra+=ProtoLogLog16c.ANDTISSA_BITS;
		if((ProtoLogLog16c.MODE & ProtoLogLog16c.MODE_NLZ2)!=0) totalExtra+=ProtoLogLog16c.NLZ2_BITS;
		if((ProtoLogLog16c.MODE & ProtoLogLog16c.MODE_HISTORY)!=0) totalExtra+=ProtoLogLog16c.HISTORY_BITS;
		if((ProtoLogLog16c.MODE & ProtoLogLog16c.MODE_LUCK)!=0) totalExtra+=ProtoLogLog16c.LUCK_BITS;
		final int numStates=1<<totalExtra;

		System.err.println("ProtoTierStats: inner="+inner+" outer="+outer+" maxtier="+maxTier
			+" mode="+ProtoLogLog16c.MODE+" extraBits="+totalExtra+" states="+numStates);
		final long t0=System.nanoTime();

		// Single-bucket simulation using ProtoLogLog16c's hashAndStore
		final long[][] count=new long[64][numStates];
		final double[][] sum=new double[64][numStates];

		final int nlzShift=16-ProtoLogLog16c.NLZ_BITS;
		final int extraMask=(1<<nlzShift)-1;

		for(int trial=0; trial<outer; trial++){
			long seed=Tools.hash64shift((long)trial*1234567891L+54321);
			ProtoLogLog16c proto=new ProtoLogLog16c(4, 2, -1, 0); // 4 buckets (min valid)

			for(int card=1; card<=inner; card++){
				long val=Tools.hash64shift(seed+card);
				proto.add(val);

				int stored=proto.maxArray[0]&0xFFFF;
				if(stored>0){
					int absNlz=(stored>>>nlzShift)-1;
					int extra=stored&extraMask;
					if(absNlz>=0 && absNlz<=maxTier && extra<numStates){
						count[absNlz][extra]++;
						sum[absNlz][extra]+=card;
					}
				}
			}
		}

		final double elapsed=(System.nanoTime()-t0)*1e-9;
		System.err.println(String.format("Elapsed: %.1fs", elapsed));

		// Report: per-tier CFs for each state
		StringBuilder header=new StringBuilder("Tier\tTotal");
		for(int s=0; s<numStates; s++){
			header.append(String.format("\tP(%d)\tCF(%d)", s, s));
		}
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
				double freq=count[t][s]/(double)tierTotal;
				sb.append('\t').append(String.format("%.4f", freq));
				if(count[t][s]>10){
					double stateAvg=sum[t][s]/count[t][s];
					double cf=Math.log(stateAvg/tierAvg)/Math.log(2);
					sb.append('\t').append(String.format("%+.4f", cf));
				}else{sb.append("\tN/A");}
			}
			System.out.println(sb);
		}

		// Observation-weighted summary for tiers 3-maxTier
		System.out.println("\n# Weighted CFs (tiers 3-"+maxTier+"):");
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
