package ddl;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import rand.FastRandomXoshiro;
import shared.Tools;

/**
 * Benchmark driver for DDL paper experiments. Generates random data internally,
 * builds sketches, and measures collision rates, similarity accuracy, and speed.
 * No file I/O — everything is synthetic.
 *
 * Modes:
 *   collision — N random sketches, all-pairs comparison, report false match rates
 *   ani      — generate pairs at known mutation rates, compare estimated vs true ANI
 *   speed    — time sketch construction over many reps
 *   scaling  — one reference sketch vs sketches of varying sizes
 *
 * Usage: DDLBenchmark mode=collision type=ddl n=50 size=5000000 buckets=2048 k=31
 *        DDLBenchmark mode=ani type=ddl size=5000000 buckets=2048 k=31
 *        DDLBenchmark mode=speed type=ddl size=5000000 buckets=2048 k=31 reps=20
 *        DDLBenchmark mode=scaling type=ddl refsize=5000000 buckets=2048 k=31
 *
 * @author Noire
 * @date May 12, 2026
 */
public class DDLBenchmark {

	public static void main(String[] args){
		String mode="collision";
		String type="ddl";
		int n=50, size=5000000, buckets=2048, k=31, reps=20, refsize=5000000;
		long seed=12345L;
		double minrate=-1, maxrate=-1, ratestep=-1;
		int trials=20;

		for(String arg : args){
			String[] ab=arg.split("=");
			String a=ab[0].toLowerCase(), b=ab.length>1 ? ab[1] : "";
			if(a.equals("mode")){mode=b;}
			else if(a.equals("type") || a.equals("loglogtype")){type=b;}
			else if(a.equals("n")){n=Integer.parseInt(b);}
			else if(a.equals("size")){size=Integer.parseInt(b);}
			else if(a.equals("refsize")){refsize=Integer.parseInt(b);}
			else if(a.equals("buckets")){buckets=Integer.parseInt(b);}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("reps")){reps=Integer.parseInt(b);}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("minrate")){minrate=Double.parseDouble(b);}
			else if(a.equals("maxrate")){maxrate=Double.parseDouble(b);}
			else if(a.equals("ratestep")){ratestep=Double.parseDouble(b);}
			else if(a.equals("trials")){trials=Integer.parseInt(b);}
			else if(a.equals("exponent") || a.equals("ebits")){DynamicDemiLog.setExponent(Integer.parseInt(b));}
		}

		if(mode.equals("collision")){collisionTest(type, n, size, buckets, k, seed);}
		else if(mode.equals("ani")){aniTest(type, size, buckets, k, seed, minrate, maxrate, ratestep, trials);}
		else if(mode.equals("speed")){speedTest(type, size, buckets, k, reps, seed);}
		else if(mode.equals("scaling")){scalingTest(type, refsize, buckets, k, seed, trials);}
		else{System.err.println("Unknown mode: "+mode);}
	}

	/** Builds a sketch from random 64-bit values (simulating unique k-mers). */
	private static CardinalityTracker buildSketch(String type, int buckets, int k,
			long sketchSeed, int numElements){
		CardinalityTracker ct=CardinalityTracker.makeTracker(type, buckets, k, sketchSeed, 0f);
		final FastRandomXoshiro rng=new FastRandomXoshiro(sketchSeed);
		for(int i=0; i<numElements; i++){
			ct.hashAndStore(rng.nextLong());
		}
		return ct;
	}

	/** Builds a sketch from a base sketch's elements plus mutations. */
	private static CardinalityTracker buildMutatedSketch(String type, int buckets, int k,
			long baseSeed, int numElements, double mutRate, long mutSeed){
		CardinalityTracker ct=CardinalityTracker.makeTracker(type, buckets, k, baseSeed, 0f);
		final FastRandomXoshiro baseRng=new FastRandomXoshiro(baseSeed);
		final FastRandomXoshiro mutRng=new FastRandomXoshiro(mutSeed);
		for(int i=0; i<numElements; i++){
			final long baseVal=baseRng.nextLong();
			if(mutRng.nextFloat()<mutRate){
				ct.hashAndStore(mutRng.nextLong());
			}else{
				ct.hashAndStore(baseVal);
			}
		}
		return ct;
	}

	private static int[] compareBuckets(int[] a, int[] b){
		int lower=0, equal=0, higher=0, bothEmpty=0;
		for(int i=0; i<a.length; i++){
			if(a[i]==0 && b[i]==0){bothEmpty++; continue;}
			if(a[i]<b[i]){lower++;}
			else if(a[i]>b[i]){higher++;}
			else{equal++;}
		}
		return new int[]{lower, equal, higher, bothEmpty};
	}

	private static float wkid(int lower, int equal, int higher){
		int div=Math.min(equal+lower, equal+higher);
		if(div<6){div=6;}
		if(equal<=0){return 0;}
		return Math.min(1f, (float)equal/div);
	}

	private static float ani(int lower, int equal, int higher, int k){
		float w=wkid(lower, equal, higher);
		return w>0 ? (float)Math.exp(Math.log(w)/k) : 0;
	}

	/** Mode: collision — N random unrelated sketches, all-pairs, report false match stats. */
	private static void collisionTest(String type, int n, int size, int buckets, int k, long seed){
		System.err.println("Collision test: type="+type+" n="+n+" size="+size+
			" buckets="+buckets+" k="+k);

		CardinalityTracker[] sketches=new CardinalityTracker[n];
		int[][] bucketVals=new int[n][];
		for(int i=0; i<n; i++){
			sketches[i]=buildSketch(type, buckets, k, seed+i*99991L, size);
			bucketVals[i]=sketches[i].bucketValues();
			if(bucketVals[i]==null){
				System.err.println("Type '"+type+"' does not support bucketValues().");
				return;
			}
		}
		System.err.println("Built "+n+" sketches.");

		long totalPairs=0, totalEqual=0;
		int maxEqual=0;
		double sumWkid=0, sumAni=0;

		System.out.println("IdxA\tIdxB\tEqual\tLower\tHigher\tBothEmpty\tWKID\tANI");
		for(int a=0; a<n; a++){
			for(int b=a+1; b<n; b++){
				int[] cmp=compareBuckets(bucketVals[a], bucketVals[b]);
				int lower=cmp[0], equal=cmp[1], higher=cmp[2], bothEmpty=cmp[3];
				float w=wkid(lower, equal, higher);
				float an=ani(lower, equal, higher, k);
				totalPairs++;
				totalEqual+=equal;
				if(equal>maxEqual){maxEqual=equal;}
				sumWkid+=w;
				sumAni+=an;
				if(equal>0){
					System.out.println(a+"\t"+b+"\t"+equal+"\t"+lower+"\t"+higher+"\t"+
						bothEmpty+"\t"+String.format("%.6f",w)+"\t"+String.format("%.6f",an));
				}
			}
		}

		System.err.println("Pairs: "+totalPairs);
		System.err.println("Total false matches: "+totalEqual);
		System.err.println("Avg false matches/pair: "+String.format("%.4f",(double)totalEqual/totalPairs));
		System.err.println("Max false matches in any pair: "+maxEqual);
		System.err.println("Avg collision rate/bucket: "+String.format("%.6f",
			(double)totalEqual/(totalPairs*buckets)));
		System.err.println("Avg WKID: "+String.format("%.6f",sumWkid/totalPairs));
		System.err.println("Avg ANI: "+String.format("%.6f",sumAni/totalPairs));
	}

	/** Mode: ani — generate pairs at known mutation rates, measure accuracy. */
	private static void aniTest(String type, int size, int buckets, int k, long seed,
			double minrate, double maxrate, double ratestep, int trials){
		final double[] rates;
		if(minrate>=0 && maxrate>0 && ratestep>0){
			int nRates=(int)Math.round((maxrate-minrate)/ratestep)+1;
			rates=new double[nRates];
			for(int i=0; i<nRates; i++){
				rates[i]=minrate+i*ratestep;
			}
		}else{
			rates=new double[]{0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07,
				0.10, 0.15, 0.20, 0.30, 0.50};
		}
		int trialsPerRate=trials;

		System.err.println("ANI accuracy test: type="+type+" size="+size+
			" buckets="+buckets+" k="+k+" trials="+trialsPerRate+
			" rates="+rates.length+" ("+String.format("%.4f",rates[0])+
			" to "+String.format("%.4f",rates[rates.length-1])+")");

		System.out.println("MutRate\tTrueANI\tTrial\tEqual\tLower\tHigher\tEstWKID\tEstANI\tANIerr");
		for(double rate : rates){
			double trueANI=1.0-rate;
			for(int trial=0; trial<trialsPerRate; trial++){
				long baseSeed=seed+trial*777L;
				CardinalityTracker base=buildSketch(type, buckets, k, baseSeed, size);
				CardinalityTracker mutant=buildMutatedSketch(type, buckets, k,
					baseSeed, size, rate, seed+trial*333L+1000000L);

				int[] bv1=base.bucketValues();
				int[] bv2=mutant.bucketValues();
				if(bv1==null || bv2==null){
					System.err.println("Type '"+type+"' does not support bucketValues().");
					return;
				}
				int[] cmp=compareBuckets(bv1, bv2);
				float w=wkid(cmp[0], cmp[1], cmp[2]);
				float an=ani(cmp[0], cmp[1], cmp[2], k);
				double err=an-trueANI;
				System.out.println(String.format("%.3f\t%.3f\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%+.6f",
					rate, trueANI, trial, cmp[1], cmp[0], cmp[2], w, an, err));
			}
		}
	}

	/** Mode: speed — time sketch construction. */
	private static void speedTest(String type, int size, int buckets, int k, int reps, long seed){
		System.err.println("Speed test: type="+type+" size="+size+
			" buckets="+buckets+" k="+k+" reps="+reps);

		// Warmup
		for(int i=0; i<3; i++){buildSketch(type, buckets, k, seed+i, size);}

		long t0=System.nanoTime();
		for(int i=0; i<reps; i++){
			buildSketch(type, buckets, k, seed+i*13L, size);
		}
		long elapsed=System.nanoTime()-t0;
		double seconds=elapsed*1e-9;
		double perRep=seconds/reps;
		double addsPerSec=(double)size/perRep;

		System.out.println("Type\tSize\tBuckets\tReps\tTotalSec\tPerRepSec\tAddsPerSec");
		System.out.println(String.format("%s\t%d\t%d\t%d\t%.3f\t%.6f\t%.0f",
			type, size, buckets, reps, seconds, perRep, addsPerSec));
	}

	/** Mode: scaling — ref sketch vs sketches of varying sizes, multi-trial averaged. */
	private static void scalingTest(String type, int refsize, int buckets, int k, long seed, int samples){
		int[] sizes={2000, 4000, 8000, 16000, 32000, 64000, 128000, 256000,
			512000, 1000000, 2000000, 4000000, 8000000, 16000000, 32000000, 64000000};

		System.err.println("Scaling test: type="+type+" refsize="+refsize+
			" buckets="+buckets+" k="+k+" samples="+samples);

		final int S=sizes.length;
		final double[][] rndStats=new double[S][7];
		final double[][] subStats=new double[S][7];
		final int[] rndGeoCount=new int[S];
		final int[] subGeoCount=new int[S];

		for(int t=0; t<samples; t++){
			final long trialSeed=seed+t*982451653L;
			final CardinalityTracker ref=buildSketch(type, buckets, k, trialSeed, refsize);
			final int[] refVals=ref.bucketValues();
			if(refVals==null){
				System.err.println("Type '"+type+"' does not support bucketValues().");
				return;
			}

			for(int si=0; si<S; si++){
				final int size=sizes[si];

				{// Random (unrelated)
					final CardinalityTracker rnd=buildSketch(type, buckets, k, trialSeed+size*7L+99999L, size);
					final int[] cmp=compareBuckets(refVals, rnd.bucketValues());
					final float w=wkid(cmp[0], cmp[1], cmp[2]);
					final float an=ani(cmp[0], cmp[1], cmp[2], k);
					rndStats[si][0]+=cmp[1]; rndStats[si][1]+=cmp[0];
					rndStats[si][2]+=cmp[2]; rndStats[si][3]+=cmp[3];
					rndStats[si][4]+=w; rndStats[si][5]+=an;
					if(w>0){rndStats[si][6]+=Math.log(w); rndGeoCount[si]++;}
				}

				{// Subset (shared elements — same seed as ref, fewer elements)
					final CardinalityTracker sub=buildSketch(type, buckets, k, trialSeed, size);
					final int[] cmp=compareBuckets(refVals, sub.bucketValues());
					final float w=wkid(cmp[0], cmp[1], cmp[2]);
					final float an=ani(cmp[0], cmp[1], cmp[2], k);
					subStats[si][0]+=cmp[1]; subStats[si][1]+=cmp[0];
					subStats[si][2]+=cmp[2]; subStats[si][3]+=cmp[3];
					subStats[si][4]+=w; subStats[si][5]+=an;
					if(w>0){subStats[si][6]+=Math.log(w); subGeoCount[si]++;}
				}
			}
			if((t+1)%10==0){System.err.println("  trial "+(t+1)+"/"+samples);}
		}

		final double n=samples;
		System.out.println("Size\tMode\tEqual\tLower\tHigher\tBothEmpty\tWKID\tANI\tGeoWKID");
		for(int si=0; si<S; si++){
			final double rndGeo=rndGeoCount[si]>0 ? Math.exp(rndStats[si][6]/rndGeoCount[si]) : 0;
			System.out.println(String.format("%d\trandom\t%.1f\t%.1f\t%.1f\t%.1f\t%.6f\t%.6f\t%.6f",
				sizes[si], rndStats[si][0]/n, rndStats[si][1]/n, rndStats[si][2]/n, rndStats[si][3]/n,
				rndStats[si][4]/n, rndStats[si][5]/n, rndGeo));
			final double subGeo=subGeoCount[si]>0 ? Math.exp(subStats[si][6]/subGeoCount[si]) : 0;
			System.out.println(String.format("%d\tsubset\t%.1f\t%.1f\t%.1f\t%.1f\t%.6f\t%.6f\t%.6f",
				sizes[si], subStats[si][0]/n, subStats[si][1]/n, subStats[si][2]/n, subStats[si][3]/n,
				subStats[si][4]/n, subStats[si][5]/n, subGeo));
		}
	}
}
