package cardinality;

import rand.FastRandomXoshiro;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Maps V_x (buckets below tier x) to true cardinality for each DLC tier.
 *
 * For each DLL2 instance, adds elements incrementally and snapshots the state.
 * At each snapshot, for each tier 0-6, computes V_x = empties + buckets below tier x,
 * and accumulates trueCard into a [V_x][tier] table.
 *
 * Output: TSV with columns: Vx, Count, T0_avgCard, T1_avgCard, ..., T6_avgCard
 * where each T column is the average true cardinality observed at that V_x for that tier.
 *
 * Usage: java cardinality.DLCTierVxMap [buckets=2048] [ddls=100000] [maxmult=5000]
 *        [threads=4] [loglogtype=dll2]
 *
 * @author Eru, Chloe
 */
public class DLCTierVxMap{

	/*--------------------------------------------------------------*/
	/*----------------        Main Method          ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		int buckets=2048;
		int ddls=100000;
		int maxmult=5000;
		long seed=1;
		int threads=4;
		String loglogtype="dll2";

		for(String arg : args){
			final String[] split=arg.split("=");
			final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : "";
			if(a.equals("buckets")){buckets=Integer.parseInt(b);}
			else if(a.equals("ddls")){ddls=Integer.parseInt(b);}
			else if(a.equals("maxmult")){maxmult=Integer.parseInt(b);}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("threads") || a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("loglogtype") || a.equals("type")){loglogtype=b.toLowerCase();}
			else if(a.equals("ignoreoverflow") || a.equals("io")){
				final boolean v=b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");
				DynamicLogLog2.IGNORE_OVERFLOW=v;
				DynamicLogLog3.IGNORE_OVERFLOW=v;
			}else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}

		final int B=buckets;
		final long maxCard=(long)B*maxmult;
		final int numBins=B+1; // V_x ranges from 0 to B

		System.err.println("DLCTierVxMap: type="+loglogtype+" B="+B+" ddls="+ddls+
			" maxmult="+maxmult+" threads="+threads+" maxCard="+maxCard+
			" io="+DynamicLogLog2.IGNORE_OVERFLOW);

		// Precompute snapshot points: every 1% cardinality increment, log-spaced
		// Use ~2000 points to get good coverage
		final int numPoints=2000;
		final long[] thresholds;
		{
			final java.util.TreeSet<Long> set=new java.util.TreeSet<>();
			for(int p=0; p<=numPoints; p++){
				final double frac=(double)p/numPoints;
				final long t=Math.max(1, (long)Math.pow(maxCard, frac));
				set.add(t);
			}
			thresholds=new long[set.size()];
			int i=0;
			for(long v : set){thresholds[i++]=v;}
		}
		final int nThresh=thresholds.length;

		// Global accumulators: [V_x][tier] -> sum of trueCard, count
		// V_x ranges 0..B, tier ranges 0..NUM_TIERS-1
		final double[][] globalSumCard=new double[numBins][NUM_TIERS];
		final long[][] globalCounts=new long[numBins][NUM_TIERS];

		final AtomicInteger ddlCounter=new AtomicInteger(0);
		final int totalDdls=ddls;
		final long sd=seed;
		final String type=loglogtype;
		final long t0=System.currentTimeMillis();

		final Thread[] workers=new Thread[threads];
		for(int t=0; t<threads; t++){
			workers[t]=new Thread(){
				@Override
				public void run(){
					final double[][] localSum=new double[numBins][NUM_TIERS];
					final long[][] localCounts=new long[numBins][NUM_TIERS];

					int d;
					while((d=ddlCounter.getAndIncrement())<totalDdls){
						final long instanceSeed=sd+d*12345678901L;
						final CardinalityTracker ct=DDLCalibrationDriver.makeInstance(
							type, B, 31, instanceSeed, 0);
						final FastRandomXoshiro rng=new FastRandomXoshiro(instanceSeed);

						long added=0;
						for(int ti=0; ti<nThresh; ti++){
							final long target=thresholds[ti];
							while(added<target){
								ct.add(rng.nextLong());
								added++;
							}

							// Get NLZ histogram
							ct.rawEstimates();
							final int[] nlz=ct.nlzCounts;
							if(nlz==null){continue;}

							// Compute V_x for each tier and record
							int vx=nlz[0]; // start with empties
							for(int tier=0; tier<NUM_TIERS; tier++){
								// V_x for this tier = empties + buckets at tiers below
								if(vx>=0 && vx<=B){
									localSum[vx][tier]+=target;
									localCounts[vx][tier]++;
								}
								// Add this tier's count to get V_x for next tier
								if(tier+1<nlz.length){
									vx+=nlz[tier+1]; // nlz[k+1] = buckets at absNlz=k
								}
							}
						}

						if((d+1)%1000==0){
							final long elapsed=System.currentTimeMillis()-t0;
							System.err.println("  DDL "+(d+1)+"/"+totalDdls+
								" ("+String.format("%.1f", elapsed/1000.0)+"s)");
						}
					}

					// Merge into global
					synchronized(globalCounts){
						for(int v=0; v<numBins; v++){
							for(int tier=0; tier<NUM_TIERS; tier++){
								globalSumCard[v][tier]+=localSum[v][tier];
								globalCounts[v][tier]+=localCounts[v][tier];
							}
						}
					}
				}
			};
			workers[t].start();
		}

		for(Thread w : workers){
			try{w.join();}catch(InterruptedException e){e.printStackTrace();}
		}

		// Output header
		System.out.print("Vx\tCount");
		for(int tier=0; tier<NUM_TIERS; tier++){
			System.out.print("\tT"+tier+"_avgCard\tT"+tier+"_count");
		}
		System.out.println();

		// Output rows
		for(int vx=0; vx<=B; vx++){
			long totalCount=0;
			for(int tier=0; tier<NUM_TIERS; tier++){totalCount+=globalCounts[vx][tier];}
			if(totalCount==0){continue;}

			System.out.print(vx+"\t"+totalCount);
			for(int tier=0; tier<NUM_TIERS; tier++){
				if(globalCounts[vx][tier]>0){
					System.out.printf("\t%.2f\t%d",
						globalSumCard[vx][tier]/globalCounts[vx][tier],
						globalCounts[vx][tier]);
				}else{
					System.out.print("\tNaN\t0");
				}
			}
			System.out.println();
		}

		long total=0;
		for(int v=0; v<=B; v++){
			for(int tier=0; tier<NUM_TIERS; tier++){total+=globalCounts[v][tier];}
		}
		final long elapsed=System.currentTimeMillis()-t0;
		System.err.println("Done. Total samples: "+total+
			" Elapsed: "+String.format("%.1f", elapsed/1000.0)+"s");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Constants            ----------------*/
	/*--------------------------------------------------------------*/

	static final int NUM_TIERS=11; // tiers 0-10
}
