package cardinality;

/**
 * FllSkewTest: does FLL2's duplicate signature (lsbFrac, msbFrac) stay
 * identifiable under SKEWED real-world frequency distributions?
 * <p>
 * Training is identical to FllDupCal (structured interleaved-sweep duplication,
 * 6-coefficient quadratic in (lsbFrac, msbFrac)).  Evaluation adds Zipf(s)
 * streams over a key pool: hot keys repeat constantly, cold keys are nearly
 * unique — a time structure the training never saw.  If the corrector trained
 * on interleaved sweeps still lands, the manifold is distribution-robust;
 * if not, the residual pattern tells us what FLL3's correction must key on.
 * True distinct counts tracked exactly via boolean[].
 *
 * Usage: java cardinality.FllSkewTest [buckets] [trainInst] [evalInst] [threads]
 *
 * @author Amber
 * @date July 2026
 */
public class FllSkewTest {

	static final long GOLD=0x9E3779B97F4A7C15L;
	static final double LOG2=Math.log(2.0);

	public static void main(String[] args) throws Exception{
		final int B=(args.length>0) ? Integer.parseInt(args[0]) : 2048;
		final int trainInst=(args.length>1) ? Integer.parseInt(args[1]) : 12;
		final int evalInst=(args.length>2) ? Integer.parseInt(args[2]) : 16;
		final int threads=(args.length>3) ? Integer.parseInt(args[3])
			: Runtime.getRuntime().availableProcessors();
		CardinalityTracker.clampToAdded=false;

		// ---- Training: same protocol as FllDupCal, sized to n/B in [24, 390] ----
		final long[] ns={24L*B, 49L*B, 98L*B, 195L*B, 390L*B};
		final double[] dups={1, 1.5, 2, 2.5, 3, 4, 8};
		final java.util.List<double[]> X=java.util.Collections.synchronizedList(new java.util.ArrayList<>());
		final java.util.List<Double> Y=java.util.Collections.synchronizedList(new java.util.ArrayList<>());
		final java.util.List<double[]> X1=java.util.Collections.synchronizedList(new java.util.ArrayList<>());
		final java.util.List<Double> Y1=java.util.Collections.synchronizedList(new java.util.ArrayList<>());

		final java.util.ArrayList<double[]> jobs=new java.util.ArrayList<>();
		for(long n : ns){for(double dup : dups){for(int i=0; i<trainInst; i++){
			jobs.add(new double[]{n, dup, i});
		}}}
		runPool(threads, jobs.size(), j->{
			final double[] job=jobs.get(j);
			final long n=(long)job[0];
			final double dup=job[1];
			final int inst=(int)job[2];
			final FutureLogLog2 f=new FutureLogLog2(B, 31, -1, 0);
			final long seedBase=new java.util.Random(4242L+inst*7919L).nextLong();
			final long total=(long)(n*dup);
			long added=0;
			outer:
			while(true){
				for(long i=0; i<n; i++){
					f.add(seedBase+i*GOLD);
					if(++added>=total){break outer;}
				}
			}
			final double[] s=FllDupCal.summarize(f);
			final double y=Math.log(n)/LOG2-s[0];
			final double[] x=FllDupCal.features(s[1], s[2]);
			synchronized(X){
				X.add(x); Y.add(y);
				if(dup==1){X1.add(x); Y1.add(y);}
			}
		});

		final double[] beta=FllDupCal.lsq(X, Y);
		// State-blind baseline: constant correction fit on unique streams only.
		// (A quadratic fit on dup=1 points is unstable off-manifold: the unique
		// manifold is nearly a point in (lsb,msb), so extrapolation is arbitrary.)
		double meanY1=0;
		for(double y : Y1){meanY1+=y;}
		meanY1/=Y1.size();
		final double[] beta1={meanY1, 0, 0, 0, 0, 0};
		System.out.println("# dup-aware beta: "+java.util.Arrays.toString(beta));
		System.out.println("# naive beta:     "+java.util.Arrays.toString(beta1));

		// ---- Evaluation ----
		System.out.println("#stream\tpool\tdraws\ttrueN\tdupRatio\tmeanExp-log2N\tlsbFrac\tmsbFrac\tnaive_err%\tdupaware_err%\tgated_err%");
		// Reference points: unique + Brian's LC protocol (sized to ~146*B and ~214*B)
		evalConfig("HC_unique", 146*B, 146L*B, 0, B, evalInst, threads, beta, beta1);
		evalConfig("LC_minrand", 244*B, 976L*B, -1, B, evalInst, threads, beta, beta1);
		// Zipf sweeps: skew s x pool x pressure (pools sized to B for regime match)
		for(double s : new double[]{0.5, 1.0, 1.5, 2.0}){
			// Zipf(2) distinct grows ~sqrt(draws): only a small-B run can reach the
			// calibrated n/B band at feasible cost, and it needs heavier multipliers.
			if(s>=2.0 && B>512){continue;}
			for(int P : new int[]{98*B, 488*B}){
				final int[] mults=(s<2.0) ? new int[]{4, 16} : new int[]{64, 256};
				for(int mult : mults){
					evalConfig(String.format("Zipf_%.1f", s), P, (long)mult*P, s, B,
						evalInst, threads, beta, beta1);
				}
			}
		}
	}

	/** Gate: unique-manifold pin is msbFrac~=0.33; beyond the gate, duplication is proven. */
	static final double MSB_GATE=0.37;

	/** skew>0: Zipf(skew). skew==0: unique stream of 'draws' keys. skew==-1: LC minrand. */
	static void evalConfig(String name, int pool, long draws, double skew, int B,
			int instances, int threads, double[] beta, double[] beta1) throws Exception{
		// Zipf cumulative weights (shared, read-only)
		final double[] cum;
		if(skew>0){
			cum=new double[pool];
			double t=0;
			for(int i=0; i<pool; i++){t+=Math.pow(i+1, -skew); cum[i]=t;}
			for(int i=0; i<pool; i++){cum[i]/=t;}
		}else{cum=null;}

		final double[] acc=new double[8]; // e_naive, e_aware, trueN, dupRatio, dExp, lsb, msb, e_gated
		final Object lock=new Object();
		runPool(threads, instances, inst->{
			final FutureLogLog2 f=new FutureLogLog2(B, 31, -1, 0);
			final java.util.Random seedr=new java.util.Random(999000L+inst*104729L);
			final long seedBase=seedr.nextLong();
			final java.util.SplittableRandom rng=new java.util.SplittableRandom(seedr.nextLong());
			long trueN;
			if(skew==0){
				for(long i=0; i<draws; i++){f.add(seedBase+i*GOLD);}
				trueN=draws;
			}else{
				final boolean[] seen=new boolean[pool];
				long distinct=0;
				for(long d=0; d<draws; d++){
					final int pos;
					if(skew<0){pos=Math.min(rng.nextInt(pool), rng.nextInt(pool));}
					else{
						final double u=rng.nextDouble();
						int lo=0, hi=pool-1;
						while(lo<hi){
							final int mid=(lo+hi)>>>1;
							if(cum[mid]<u){lo=mid+1;}else{hi=mid;}
						}
						pos=lo;
					}
					if(!seen[pos]){seen[pos]=true; distinct++;}
					f.add(seedBase+pos*GOLD);
				}
				trueN=distinct;
			}
			final double[] s=FllDupCal.summarize(f);
			final double[] x=FllDupCal.features(s[1], s[2]);
			final double estA=Math.pow(2.0, s[0]+FllDupCal.dot(x, beta));
			final double estN=Math.pow(2.0, s[0]+FllDupCal.dot(x, beta1));
			final double estG=(s[2]<MSB_GATE) ? estN : estA;
			synchronized(lock){
				acc[0]+=Math.abs(estN-trueN)/trueN;
				acc[1]+=Math.abs(estA-trueN)/trueN;
				acc[2]+=trueN;
				acc[3]+=(double)draws/trueN;
				acc[4]+=s[0]-Math.log(trueN)/LOG2;
				acc[5]+=s[1];
				acc[6]+=s[2];
				acc[7]+=Math.abs(estG-trueN)/trueN;
			}
		});
		System.out.println(String.format("%s\t%d\t%d\t%.0f\t%.2f\t%.3f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f",
			name, pool, draws, acc[2]/instances, acc[3]/instances, acc[4]/instances,
			acc[5]/instances, acc[6]/instances,
			100.0*acc[0]/instances, 100.0*acc[1]/instances, 100.0*acc[7]/instances));
	}

	interface Job{void run(int idx);}

	static void runPool(int threads, int count, Job job) throws Exception{
		final java.util.concurrent.atomic.AtomicInteger next=new java.util.concurrent.atomic.AtomicInteger(0);
		final Thread[] pool=new Thread[Math.min(threads, count)];
		for(int t=0; t<pool.length; t++){
			pool[t]=new Thread(()->{
				int i;
				while((i=next.getAndIncrement())<count){job.run(i);}
			});
			pool[t].start();
		}
		for(Thread th : pool){th.join();}
	}
}
