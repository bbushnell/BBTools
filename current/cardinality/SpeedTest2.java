package cardinality;

import rand.FastRandomXoshiro;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Cache-pressure throughput benchmark for cardinality estimators.
 * Each thread maintains N simultaneous estimators, cycling through them
 * every add. Simulates real-world deployment with many active sketches.
 *
 * @author Brian Bushnell, Chloe
 */
public class SpeedTest2 {

	public static void main(String[] args){
		CardinalityTracker.clampToAdded=false;

		int buckets=2048;
		int sim=128;
		long maxCard=40_000_000L;
		int threads=Runtime.getRuntime().availableProcessors();
		String typeFilter=null;

		for(int i=0; i<args.length; i++){
			final String arg=args[i].toLowerCase();
			final int eq=arg.indexOf('=');
			if(eq<0){continue;}
			final String key=arg.substring(0, eq);
			final String val=arg.substring(eq+1);
			if(key.equals("buckets") || key.equals("b")){buckets=Integer.parseInt(val);}
			else if(key.equals("sim") || key.equals("simultaneous")){sim=Integer.parseInt(val);}
			else if(key.equals("maxcard") || key.equals("card")){maxCard=Long.parseLong(val);}
			else if(key.equals("threads") || key.equals("t")){threads=Integer.parseInt(val);}
			else if(key.equals("type") || key.equals("loglogtype")){typeFilter=val;}
		}

		final int totalSim=threads*sim;
		final String[] allTypes={"avll", "exa", "udll6", "ull", "dll4", "ll6", "hll4", "hlll", "htc4", "htb"};
		final String[] allLabels={"AVLL", "EXA", "UDLL6", "ULL", "DLL4", "LL6", "HLL4", "HLLL", "HTC4", "HTB"};

		System.err.println("SpeedTest2: buckets="+buckets+" sim/thread="+sim+
			" totalSim="+totalSim+" maxCard="+maxCard+" threads="+threads+
			(typeFilter!=null ? " type="+typeFilter : " type=ALL"));
		System.out.println("Type\tM_adds/s\tAvgCard\tExpectedCard\tSimultaneous");

		for(int i=0; i<allTypes.length; i++){
			if(typeFilter!=null && !allTypes[i].equals(typeFilter)){continue;}
			runTest(allLabels[i], allTypes[i], buckets, sim, maxCard, threads);
		}
	}

	static void runTest(String label, String type, int buckets,
			int sim, long maxCard, int threads){
		final int k=31;
		final AtomicLong cardSum=new AtomicLong(0);
		final int totalSim=threads*sim;

		final Thread[] workers=new Thread[threads];
		final long t0=System.nanoTime();

		for(int ti=0; ti<threads; ti++){
			final int threadId=ti;
			workers[ti]=new Thread(()->{
				final FastRandomXoshiro rng=new FastRandomXoshiro(threadId*982451653L+123456789L);
				final CardinalityTracker[] cts=new CardinalityTracker[sim];
				final long[] masks=new long[sim];
				for(int i=0; i<sim; i++){
					cts[i]=DDLCalibrationDriver.makeInstance(type, buckets, k, 0, 0);
					masks[i]=rng.nextLong();
				}

				for(long c=0; c<maxCard; c++){
					final long base=rng.nextLong();
					for(int i=0; i<sim; i++){
						cts[i].add(base^masks[i]);
					}
				}

				long localCardSum=0;
				for(int i=0; i<sim; i++){
					cts[i].lastCardinality=-1;
					localCardSum+=cts[i].cardinality();
				}
				cardSum.addAndGet(localCardSum);
			});
			workers[ti].start();
		}

		for(Thread w : workers){
			try{w.join();}catch(InterruptedException e){Thread.currentThread().interrupt();}
		}

		final long t1=System.nanoTime();
		final double elapsedSec=(t1-t0)/1e9;
		final double totalAdds=(double)totalSim*maxCard;
		final double mAddsPerSec=totalAdds/elapsedSec/1e6;
		final double avgCard=(double)cardSum.get()/totalSim;

		System.out.printf("%s\t%.1f\t%.0f\t%d\t%d%n", label, mAddsPerSec, avgCard, maxCard, totalSim);
		System.err.printf("  %s: %.1f M adds/s (%.1f s, avg card=%.0f, %d simultaneous)%n",
			label, mAddsPerSec, elapsedSec, avgCard, totalSim);
	}

}
