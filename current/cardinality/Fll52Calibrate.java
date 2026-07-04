package cardinality;

/**
 * Fll52Calibrate: class-distribution table for likelihood extraction from
 * Fll52 (5-antenna organisms).  Per-word class = min(localExp,7)*56 +
 * idx56(antenna bits): 5 buckets x 4 states -> C(8,3)=56 equivalence classes.
 * Table conditioned on (globalExp, phase m) exactly as Fll2Calibrate
 * (M_OFFSET for negative m, per-instance ladder jitter, G-conditioning).
 *
 * Usage: java cardinality.Fll52Calibrate [buckets] [instances] [maxMult] [threads] [out.tsv]
 *
 * @author Amber
 * @date July 2026
 */
public class Fll52Calibrate {

	static final int NCLASS=56;
	static final int NC=8*NCLASS;
	static final int BINS_PER_UNIT=20;
	static final int M_OFFSET=12;
	static final int MAX_BINS=440;
	static final int MAX_G=32;
	/** Prefix sums over a=count(00 buckets) for 5 buckets: r=5-a slots. */
	static final int[] OFFSET_A={0, 21, 36, 46, 52, 55};

	/** 10-bit antenna field (word bits 2..11 shifted down) -> 56-class index. */
	static int idx56(int hist){
		final int h1=hist>>>1;
		final int LSB5=0x155;   // bits 0,2,4,6,8 of the shifted field
		final int d=Integer.bitCount(hist&h1&LSB5);
		final int a=5-Integer.bitCount((hist|h1)&LSB5);
		final int c=Integer.bitCount(~hist&h1&LSB5);
		final int b=5-a-c-d;
		final int r=5-a;
		return OFFSET_A[a]+b*(r+1)-b*(b-1)/2+c;
	}

	static int classOf(int word){
		final int le=Math.min((word>>>12)&0xF, 7);
		return le*NCLASS+idx56((word>>>2)&0x3FF);
	}

	interface Job {void run(int idx);}

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

	public static void main(String[] args) throws Exception{
		final int B=(args.length>0) ? Integer.parseInt(args[0]) : 2560;
		final int instances=(args.length>1) ? Integer.parseInt(args[1]) : 512;
		final long maxMult=(args.length>2) ? Long.parseLong(args[2]) : 2048;
		final int threads=(args.length>3) ? Integer.parseInt(args[3])
			: Runtime.getRuntime().availableProcessors();
		final String outPath=(args.length>4) ? args[4] : "fll52mle_table.tsv";
		CardinalityTracker.clampToAdded=false;
		final double LOG2=Math.log(2.0);

		final java.util.ArrayList<Long> cpl=new java.util.ArrayList<>();
		long nextCp=1;
		while(nextCp<(long)maxMult*B){
			cpl.add(nextCp);
			nextCp=Math.max(nextCp+1, (long)(nextCp*1.03));
		}
		cpl.add((long)maxMult*B);
		final long[] cps=new long[cpl.size()];
		for(int i=0; i<cps.length; i++){cps[i]=cpl.get(i);}
		System.err.println("Fll52Calibrate B="+B+" inst="+instances+" checkpoints="+cps.length);

		final long[][][] hist=new long[MAX_G][MAX_BINS][NC];
		final Object lock=new Object();
		runPool(threads, instances, inst->{
			final long[][][] lh=new long[MAX_G][MAX_BINS][NC];
			final Fll52 f=new Fll52(B, 31, -1, 0);
			final java.util.SplittableRandom randy=new java.util.SplittableRandom(
				new java.util.Random(818000L+inst).nextLong());
			final int nw=f.getNumWords(), nb=nw*5;
			final double jit=Math.pow(2.0, randy.nextDouble());
			final long[] icps=new long[cps.length];
			for(int i=0; i<cps.length; i++){icps[i]=Math.max(1, (long)(cps[i]*jit));}
			int cpIdx=0;
			for(long n=1; cpIdx<icps.length; n++){
				f.add(randy.nextLong());
				if(n==icps[cpIdx]){
					cpIdx++;
					final int G=f.getGlobalExp();
					final double m=Math.log((double)n/nb)/LOG2-G;
					final int bin=(int)((m+M_OFFSET)*BINS_PER_UNIT);
					if(bin<0 || bin>=MAX_BINS || G<0 || G>=MAX_G){continue;}
					for(int w=0; w<nw; w++){
						lh[G][bin][classOf(f.getWord(w))]++;
					}
				}
			}
			synchronized(lock){
				for(int g=0; g<MAX_G; g++){
					for(int b=0; b<MAX_BINS; b++){
						for(int k=0; k<NC; k++){hist[g][b][k]+=lh[g][b][k];}
					}
				}
			}
		});

		final StringBuilder sb=new StringBuilder();
		sb.append("#G\tbin\tmCenter\ttotal\tcounts["+NC+"]\n");
		for(int g=0; g<MAX_G; g++){
			for(int b=0; b<MAX_BINS; b++){
				long tot=0;
				for(long v : hist[g][b]){tot+=v;}
				if(tot==0){continue;}
				sb.append(g).append('\t').append(b).append('\t')
					.append(String.format("%.4f", (b+0.5)/BINS_PER_UNIT-M_OFFSET)).append('\t')
					.append(tot);
				for(long v : hist[g][b]){sb.append('\t').append(v);}
				sb.append('\n');
			}
		}
		java.nio.file.Files.write(java.nio.file.Paths.get(outPath), sb.toString().getBytes());
		System.err.println("Wrote "+outPath);
	}
}
