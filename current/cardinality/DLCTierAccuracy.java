package cardinality;

import rand.FastRandomXoshiro;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Measures DLC per-tier absolute error as a function of tier occupancy.
 * Uses LogLog6 (simplest class: no minZeros, no microIndex, no promotion).
 * Tier occupancy = number of buckets with NLZ >= tier (integer 0..B).
 *
 * Each thread creates one LL6, feeds it maxCard elements, and at each
 * log-spaced cardinality threshold samples the tier occupancy and error.
 * Then repeats with a new LL6 until ddls are exhausted.
 *
 * Output: TSV with columns: Occupancy, Count, AvgAbsErr, AvgSignedErr
 *
 * Usage: java cardinality.DLCTierAccuracy [buckets=2048] [ddls=2000000]
 *        [maxmult=512] [seed=1] [tier=3] [threads=128] [points=500]
 *
 * @author Eru
 */
public class DLCTierAccuracy {

	public static void main(String[] args){
		int buckets=2048;
		int ddls=100000;
		int maxmult=512;
		long seed=1;
		int targetTier=3;
		int threads=1;
		int numPoints=500;

		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : "";
			if(a.equals("buckets")){buckets=Integer.parseInt(b);}
			else if(a.equals("ddls")){ddls=Integer.parseInt(b);}
			else if(a.equals("maxmult")){maxmult=Integer.parseInt(b);}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("tier")){targetTier=Integer.parseInt(b);}
			else if(a.equals("threads") || a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("points")){numPoints=Integer.parseInt(b);}
		}

		final long maxCard=(long)buckets*maxmult;
		final int numBins=buckets+1;

		System.err.println("DLCTierAccuracy (LL6): buckets="+buckets+" ddls="+ddls+
			" maxmult="+maxmult+" tier="+targetTier+" threads="+threads+
			" points="+numPoints+" maxCard="+maxCard);

		// Precompute cardinality thresholds (log-spaced)
		final long[] thresholds=new long[numPoints];
		for(int p=0; p<numPoints; p++){
			final double frac=(double)p/(numPoints-1);
			thresholds[p]=Math.max(1, (long)Math.pow(maxCard, frac));
		}
		// Deduplicate (low end may have repeats)
		final long[] thresh;
		{
			int unique=1;
			for(int i=1; i<numPoints; i++){
				if(thresholds[i]>thresholds[i-1]){unique++;}
			}
			thresh=new long[unique];
			thresh[0]=thresholds[0];
			int j=1;
			for(int i=1; i<numPoints; i++){
				if(thresholds[i]>thresholds[i-1]){thresh[j++]=thresholds[i];}
			}
		}
		final int nThresh=thresh.length;

		// Global accumulators
		final long[] globalCounts=new long[numBins];
		final double[] globalAbsErr=new double[numBins];
		final double[] globalSignedErr=new double[numBins];

		// Shared DDL counter
		final AtomicInteger ddlCounter=new AtomicInteger(0);

		final int B=buckets;
		final int tier=targetTier;
		final long sd=seed;
		final int totalDdls=ddls;
		final long t0=System.currentTimeMillis();

		final Thread[] workers=new Thread[threads];
		for(int t=0; t<threads; t++){
			workers[t]=new Thread(){
				@Override
				public void run(){
					final long[] localCounts=new long[numBins];
					final double[] localAbsErr=new double[numBins];
					final double[] localSignedErr=new double[numBins];

					int d;
					while((d=ddlCounter.getAndIncrement())<totalDdls){
						final long instanceSeed=sd+d*12345678901L;
						final LogLog6 ll=new LogLog6(B, 31, instanceSeed, 0);
						final FastRandomXoshiro rng=new FastRandomXoshiro(instanceSeed);

						long added=0;
						for(int ti=0; ti<nThresh; ti++){
							final long target=thresh[ti];
							// Feed elements up to this threshold
							while(added<target){
								ll.add(rng.nextLong());
								added++;
							}

							// Sample tier occupancy and error
							ll.rawEstimates();
							final int[] nlzCounts=ll.lastRawNlz;
							if(nlzCounts==null){continue;}

							int filled=0;
							for(int k=0; k<nlzCounts.length; k++){filled+=nlzCounts[k];}
							int vk=B-filled;
							for(int k=0; k<tier && k<nlzCounts.length; k++){
								vk+=nlzCounts[k];
							}

							final int tierOcc=B-vk;
							if(tierOcc<0 || tierOcc>B){continue;}

							if(vk>=1 && vk<B){
								final double tierEst=(1L<<tier)*(double)B*Math.log((double)B/vk);
								final double absErr=Math.abs(tierEst-target)/(double)target;
								final double signedErr=(tierEst-target)/(double)target;

								localCounts[tierOcc]++;
								localAbsErr[tierOcc]+=absErr;
								localSignedErr[tierOcc]+=signedErr;
							}
						}

						// Progress every 1000 DDLs
						if((d+1)%1000==0){
							final long elapsed=System.currentTimeMillis()-t0;
							System.err.println("  DDL "+(d+1)+"/"+totalDdls+
								" ("+String.format("%.1f", elapsed/1000.0)+"s)");
						}
					}

					// Merge into global
					synchronized(globalCounts){
						for(int i=0; i<numBins; i++){
							globalCounts[i]+=localCounts[i];
							globalAbsErr[i]+=localAbsErr[i];
							globalSignedErr[i]+=localSignedErr[i];
						}
					}
				}
			};
			workers[t].start();
		}

		for(Thread w : workers){
			try{w.join();}catch(InterruptedException e){e.printStackTrace();}
		}

		// Output
		System.out.println("Occupancy\tCount\tAvgAbsErr\tAvgSignedErr");
		for(int occ=0; occ<numBins; occ++){
			if(globalCounts[occ]>0){
				System.out.println(occ+"\t"+globalCounts[occ]+"\t"+
					String.format("%.8f", globalAbsErr[occ]/globalCounts[occ])+"\t"+
					String.format("%.8f", globalSignedErr[occ]/globalCounts[occ]));
			}
		}

		long total=0;
		for(long c : globalCounts){total+=c;}
		final long elapsed=System.currentTimeMillis()-t0;
		System.err.println("Done. Total samples: "+total+
			" Elapsed: "+String.format("%.1f", elapsed/1000.0)+"s");
	}
}
