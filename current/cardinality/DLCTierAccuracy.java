package cardinality;

import rand.FastRandomXoshiro;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Measures DLC per-tier absolute error as a function of tier occupancy.
 * Supports any CardinalityTracker type via loglogtype flag.
 * Tier occupancy = number of buckets with NLZ >= tier (integer 0..B).
 *
 * Each thread creates one tracker, feeds it maxCard elements, and at each
 * log-spaced cardinality threshold samples the tier occupancy and error.
 * Then repeats with a new tracker until ddls are exhausted.
 *
 * Output: TSV with columns: Occupancy, Count, AvgAbsErr, AvgSignedErr
 *
 * Usage: java cardinality.DLCTierAccuracy [buckets=2048] [ddls=2000000]
 *        [maxmult=512] [seed=1] [tier=3] [threads=128] [points=500]
 *        [type=ll6]
 *
 * @author Eru, Chloe
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
		String loglogtype="ll6";

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
			else if(a.equals("loglogtype") || a.equals("type")){loglogtype=b.toLowerCase();}
			else if(a.equals("correctoverflow") || a.equals("co")){
				final boolean v=b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");
				DynamicLogLog3.CORRECT_OVERFLOW=v;
				BankedDynamicLogLog3.CORRECT_OVERFLOW=v;
				if(v){DynamicLogLog3.IGNORE_OVERFLOW=false; DynamicLogLog2.IGNORE_OVERFLOW=false; BankedDynamicLogLog3.IGNORE_OVERFLOW=false;}
			}else if(a.equals("ignoreoverflow") || a.equals("io")){
				final boolean v=b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");
				DynamicLogLog3.IGNORE_OVERFLOW=v;
				DynamicLogLog2.IGNORE_OVERFLOW=v;
				BankedDynamicLogLog3.IGNORE_OVERFLOW=v;
				if(v){DynamicLogLog3.CORRECT_OVERFLOW=false; BankedDynamicLogLog3.CORRECT_OVERFLOW=false;}
			}else if(a.equals("overflowscale") || a.equals("os")){
				DynamicLogLog3.OVERFLOW_SCALE=Double.parseDouble(b);
			}else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}

		final long maxCard=(long)buckets*maxmult;
		final int numBins=buckets+1;

		System.err.println("DLCTierAccuracy: type="+loglogtype+" buckets="+buckets+" ddls="+ddls+
			" maxmult="+maxmult+" tier="+targetTier+" threads="+threads+
			" points="+numPoints+" maxCard="+maxCard+
			" io="+DynamicLogLog3.IGNORE_OVERFLOW+
			" os="+DynamicLogLog3.OVERFLOW_SCALE);

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
		final String type=loglogtype;
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
						final CardinalityTracker ct=DDLCalibrationDriver.makeInstance(
							type, B, 31, instanceSeed, 0);
						final FastRandomXoshiro rng=new FastRandomXoshiro(instanceSeed);

						long added=0;
						for(int ti=0; ti<nThresh; ti++){
							final long target=thresh[ti];
							while(added<target){
								ct.add(rng.nextLong());
								added++;
							}

							// Trigger summarize to populate nlzCounts
							ct.rawEstimates();
							final int[] nlz=ct.nlzCounts;
							if(nlz==null){continue;}

							// Compute V_k for this tier: empties + buckets below tier
							// nlz format: [0]=empties, [k+1]=buckets at absNlz==k
							int vk=nlz[0]; // empties
							for(int k=0; k<tier && k+1<nlz.length; k++){
								vk+=nlz[k+1]; // buckets at absNlz==k
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
