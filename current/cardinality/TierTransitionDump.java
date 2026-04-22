package cardinality;

import rand.FastRandomXoshiro;

/**
 * Dumps average tier occupancy at each DLL3 tier transition,
 * with DLL4 comparison. Measures raw and corrected NLZ distributions
 * and stored overflow at each minZeros transition point.
 *
 * @author Brian Bushnell
 */
public class TierTransitionDump {

	/*--------------------------------------------------------------*/
	/*----------------           Main              ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final int B=512, N=512;
		final long maxTrue=(long)B*1000;
		DynamicLogLog3.CORRECT_OVERFLOW=true;
		DynamicLogLog3.USE_STORED_OVERFLOW=true;
		DynamicLogLog3.OVERFLOW_SCALE=1.0;

		for(int targetMZ=1; targetMZ<=6; targetMZ++){
			final double[] r3s=new double[20], c3s=new double[20], r4s=new double[20];
			final double[] ovfSum=new double[20]; // storedOverflow accumulator
			int cnt=0;
			for(int di=0; di<N; di++){
				final FastRandomXoshiro rng=new FastRandomXoshiro(42+di);
				final long s3=rng.nextLong()&Long.MAX_VALUE;
				final FastRandomXoshiro rng4=new FastRandomXoshiro(42+di);
				final long s4=rng4.nextLong()&Long.MAX_VALUE;
				final DynamicLogLog3 d3=new DynamicLogLog3(B, 31, s3, 0);
				final DynamicLogLog4 d4=new DynamicLogLog4(B, 31, s4, 0);
				int prev=0;
				boolean done=false;
				for(long c=1; c<=maxTrue && !done; c++){
					d3.add(rng.nextLong());
					d4.add(rng4.nextLong());
					final int mz=d3.getMinZeros();
					if(mz>=targetMZ && prev<targetMZ){
						d3.rawEstimates(); d4.rawEstimates();
						if(d3.lastRawNlz!=null){
							for(int t=0; t<20; t++){
								r3s[t]+=d3.lastRawNlz[t];
								c3s[t]+=d3.lastCorrNlz[t];
								r4s[t]+=d4.lastRawNlz[t];
							}
							final int[] so=d3.getStoredOverflow();
							for(int t=0; t<20; t++){ovfSum[t]+=so[t];}
							cnt++;
						}
						done=true;
					}
					prev=mz;
				}
			}
			// Compute cumulative raw sums (top down)
			final double[] cumRaw3=new double[20], cumRaw4=new double[20];
			cumRaw3[19]=r3s[19]/cnt; cumRaw4[19]=r4s[19]/cnt;
			for(int t=18; t>=0; t--){
				cumRaw3[t]=cumRaw3[t+1]+r3s[t]/cnt;
				cumRaw4[t]=cumRaw4[t+1]+r4s[t]/cnt;
			}

			System.out.println("=== Transition to minZeros="+targetMZ+" (tier "+(6+targetMZ)+" active, "+cnt+" DDLs) ===");
			System.out.println("Tier\tDLL3_raw\tDLL3_corr\tDLL4_raw\tcumRaw3\t\tcumRaw4\t\tstored_ovf");
			for(int t=0; t<20; t++){
				final double r=r3s[t]/cnt, c=c3s[t]/cnt, d=r4s[t]/cnt, o=ovfSum[t]/cnt;
				if(r>0.5 || c>0.5 || d>0.5){
					System.out.printf("%d\t%.1f\t\t%.1f\t\t%.1f\t\t%.1f\t\t%.1f\t\t%.1f%n",
						t, r, c, d, cumRaw3[t], cumRaw4[t], o);
				}
			}
			System.out.println();
		}
	}
}
