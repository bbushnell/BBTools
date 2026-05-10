package cardinality;

import rand.FastRandomXoshiro;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Throughput benchmark for cardinality estimator types.
 * Measures million adds/second across many estimators in parallel.
 *
 * @author Brian Bushnell, Chloe
 */
public class SpeedTest {

	public static void main(String[] args){
		CardinalityTracker.clampToAdded=false;

		int buckets=2048;
		int numEstimators=16384;
		long maxCard=40_000_000L;
		int threads=Runtime.getRuntime().availableProcessors();

		for(int i=0; i<args.length; i++){
			final String arg=args[i].toLowerCase();
			final int eq=arg.indexOf('=');
			if(eq<0){continue;}
			final String key=arg.substring(0, eq);
			final String val=arg.substring(eq+1);
			if(key.equals("buckets") || key.equals("b")){buckets=Integer.parseInt(val);}
			else if(key.equals("estimators") || key.equals("e")){numEstimators=Integer.parseInt(val);}
			else if(key.equals("maxcard") || key.equals("card")){maxCard=Long.parseLong(val);}
			else if(key.equals("threads") || key.equals("t")){threads=Integer.parseInt(val);}
		}

		final String[] types={"udll6", "ull", "dll4", "ll6", "hll4", "hlll", "htc4", "htb"};
		final String[] labels={"UDLL6", "ULL", "DLL4", "LL6", "HLL4", "HLLL", "HTC4", "HTB"};

		System.err.println("SpeedTest: buckets="+buckets+" estimators="+numEstimators+
			" maxCard="+maxCard+" threads="+threads);
		System.out.println("Type\tM_adds/s\tAvgCard\tExpectedCard");

		for(int i=0; i<types.length; i++){
			runTest(labels[i], types[i], buckets, numEstimators, maxCard, threads);
		}
	}

	static void runTest(String label, String type, int buckets,
			int numEstimators, long maxCard, int threads){
		final int k=31;
		final AtomicInteger nextEst=new AtomicInteger(0);
		final AtomicLong cardSum=new AtomicLong(0);

		final Thread[] workers=new Thread[threads];
		final long t0=System.nanoTime();

		for(int ti=0; ti<threads; ti++){
			final int threadId=ti;
			workers[ti]=new Thread(()->{
				final FastRandomXoshiro rng=new FastRandomXoshiro(threadId*982451653L+123456789L);
				long localCardSum=0;
				int est;
				while((est=nextEst.getAndIncrement())<numEstimators){
					final CardinalityTracker ct=DDLCalibrationDriver.makeInstance(type, buckets, k, 0, 0);
					for(long c=0; c<maxCard; c++){
						ct.add(rng.nextLong());
					}
					ct.lastCardinality=-1;
					localCardSum+=ct.cardinality();
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
		final double totalAdds=(double)numEstimators*maxCard;
		final double mAddsPerSec=totalAdds/elapsedSec/1e6;
		final double avgCard=(double)cardSum.get()/numEstimators;

		System.out.printf("%s\t%.1f\t%.0f\t%d%n", label, mAddsPerSec, avgCard, maxCard);
		System.err.printf("  %s: %.1f M adds/s (%.1f s, avg card=%.0f)%n",
			label, mAddsPerSec, elapsedSec, avgCard);
	}

}
