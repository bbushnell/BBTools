package cardinality;

import java.util.concurrent.atomic.AtomicInteger;

/**
 * Measures merge accuracy degradation for DLL4 as a function of the number
 * of merged instances. Generates data for Fig 19 of the DLL paper.
 * Parallelizes across trials (each trial is thread-local, no contention).
 *
 * Usage: java cardinality.MergeAccuracyTest [card=1720000000] [buckets=2048] [trials=100]
 *
 * @author Brian Bushnell, Chloe
 */
public class MergeAccuracyTest {

	public static void main(String[] args){
		long totalCard=1720000000L;
		int buckets=2048;
		int trials=100;
		int threads=Runtime.getRuntime().availableProcessors();

		for(String arg : args){
			final String[] split=arg.split("=");
			final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : null;
			if(a.equals("card")){totalCard=Long.parseLong(b);}
			else if(a.equals("buckets") || a.equals("b")){buckets=Integer.parseInt(b);}
			else if(a.equals("trials") || a.equals("t")){trials=Integer.parseInt(b);}
			else if(a.equals("threads")){threads=Integer.parseInt(b);}
		}

		final int[] splitCounts={1, 2, 4, 8, 16, 32, 64};

		System.err.println("MergeAccuracyTest: card="+totalCard+" buckets="+buckets+
			" trials="+trials+" threads="+threads);
		CardinalityTracker.clampToAdded=false;

		System.out.println("Splits\tMean_abs\tMax_abs\tMean_signed\tStddev");

		for(int sc : splitCounts){
			final long card=totalCard;
			final int bk=buckets;
			final int nTrials=trials;
			final double[] errors=new double[nTrials];
			final AtomicInteger nextTrial=new AtomicInteger(0);

			Thread[] workers=new Thread[threads];
			for(int wi=0; wi<threads; wi++){
				workers[wi]=new Thread(()->{
					int t;
					while((t=nextTrial.getAndIncrement())<nTrials){
						errors[t]=runOneTrial(card, bk, sc, t*7919L+sc*104729L);
					}
				});
				workers[wi].start();
			}
			for(Thread w : workers){try{w.join();}catch(InterruptedException e){}}

			double sumAbs=0, sumSigned=0, sumSq=0, maxAbs=0;
			for(int t=0; t<nTrials; t++){
				double e=errors[t];
				sumAbs+=Math.abs(e);
				sumSigned+=e;
				sumSq+=e*e;
				maxAbs=Math.max(maxAbs, Math.abs(e));
			}
			double meanAbs=sumAbs/nTrials;
			double meanSigned=sumSigned/nTrials;
			double variance=sumSq/nTrials - (sumSigned/nTrials)*(sumSigned/nTrials);
			double stddev=Math.sqrt(Math.max(0, variance));
			System.out.printf("%d\t%.6f\t%.6f\t%.6f\t%.6f%n", sc, meanAbs, maxAbs, meanSigned, stddev);
			System.err.printf("splits=%d  mean_abs=%.4f%%  max_abs=%.4f%%  signed=%.4f%%%n",
				sc, meanAbs*100, maxAbs*100, meanSigned*100);
		}
	}

	/** Run a single trial: create 'splits' instances, fill round-robin, merge, return relative error. */
	static double runOneTrial(long totalCard, int buckets, int splits, long seed){
		final DynamicLogLog4[] instances=new DynamicLogLog4[splits];
		for(int s=0; s<splits; s++){
			instances[s]=new DynamicLogLog4(buckets, 31, -1, 0);
		}
		for(long i=0; i<totalCard; i++){
			instances[(int)(i%splits)].add(wangHash(seed+i));
		}
		for(int s=1; s<splits; s++){
			instances[0].add(instances[s]);
		}
		long estimate=instances[0].cardinality();
		return (double)(estimate-totalCard)/(double)totalCard;
	}

	static long wangHash(long key){
		key=(~key)+(key<<21);
		key=key^(key>>>24);
		key=(key+(key<<3))+(key<<8);
		key=key^(key>>>14);
		key=(key+(key<<2))+(key<<4);
		key=key^(key>>>28);
		key=key+(key<<31);
		return key;
	}
}
