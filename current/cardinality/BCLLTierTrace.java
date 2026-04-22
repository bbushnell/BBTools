package cardinality;

/**
 * Diagnostic tool that traces per-tier DLC (distinct-level counting)
 * statistics for a BankedCeilingLogLog instance.
 * Adds a configurable number of elements, then prints ceiling distribution
 * and per-tier fill rates with entropy information.
 *
 * @author Brian Bushnell
 */
public class BCLLTierTrace {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		int numBuckets=96;
		if(args.length>0){numBuckets=Integer.parseInt(args[0]);}

		final BankedCeilingLogLog bcll=new BankedCeilingLogLog(numBuckets, 31, 42, 0);
		final int trueCard=numBuckets*10; // 10x buckets
		for(int i=1; i<=trueCard; i++){bcll.add(i);}

		// Access private tier data via rawEstimates to trigger collectTierCounts
		bcll.rawEstimates();

		System.err.println("buckets="+bcll.getModBuckets()+" words="+bcll.getNumWords()
			+" trueCard="+trueCard+" globalCeiling="+bcll.getGlobalCeiling());

		// Print ceiling distribution
		final int[] ceilHist=new int[64];
		for(int w=0; w<bcll.getNumWords(); w++){
			final int word=bcll.getWord(w);
			final int le=(word>>>12)&0xF;
			final int c=bcll.getGlobalCeiling()+le;
			if(c<64){ceilHist[c]++;}
		}
		System.err.print("Ceiling dist: ");
		for(int t=0; t<64; t++){
			if(ceilHist[t]>0){System.err.printf("C%d=%d ", t, ceilHist[t]);}
		}
		System.err.println();

		// Now manually compute per-tier DLC
		final int gc=bcll.getGlobalCeiling();
		final int maxC=gc+15;
		final int maxTier=maxC+2;
		final int[] tFilled=new int[maxTier];
		final int[] tTotal=new int[maxTier];
		final int modB=bcll.getModBuckets();

		for(int w=0; w<bcll.getNumWords(); w++){
			final int word=bcll.getWord(w);
			final int le=(word>>>12)&0xF;
			final int C=gc+le;
			final int history=word&0xFFF;

			int nonZero=0, msbSet=0;
			for(int r=0; r<6; r++){
				final int bits=(history>>>(r*2))&0x3;
				if(bits!=0){nonZero++;}
				if((bits&2)!=0){msbSet++;}
			}

			if(C<maxTier){tTotal[C]+=6; tFilled[C]+=msbSet;}
			if(C>=1 && C-1<maxTier){tTotal[C-1]+=6; tFilled[C-1]+=nonZero;}
			for(int t=C+1; t<maxTier; t++){tTotal[t]+=6;}
			for(int t=Math.max(0, C-2); t>=0; t--){tTotal[t]+=nonZero; tFilled[t]+=nonZero;}
		}

		System.err.printf("%-6s %8s %8s %8s %12s %12s%n", "Tier", "B", "Filled", "V/B", "DLC_T", "Info");
		for(int t=0; t<maxTier; t++){
			final int B=tTotal[t];
			if(B==0){continue;}
			final int filled=tFilled[t];
			final int V=B-filled;
			final double vb=(double)V/B;
			double dlc=-1;
			double info=0;
			if(V>0 && V<B){
				dlc=(1L<<t)*(-B*Math.log(vb));
				final double frac=vb;
				info=-B*(frac*Math.log(frac)+(1-frac)*Math.log(1-frac));
			}
			System.err.printf("%-6d %8d %8d %8.4f %12.1f %12.1f%n", t, B, filled, vb, dlc, info);
		}
	}
}
