package cardinality;

/**
 * FllComplexity: Brian's 2026-07-04 idea — a predictor optimized for STREAM
 * COMPLEXITY (distinct/adds), graded on its own accuracy, decoupled from
 * cardinality estimation.  The estimator always knows adds; the question is
 * whether the FLL2 state distinguishes streams with the same adds and same
 * distinct but different duplication STRUCTURE (uniform 2x vs one hot key
 * added millions of times with a unique tail).
 * <p>
 * Baseline model: complexity from r = log2(adds) - meanExp only (the "rough
 * sense of duplicates" from counting adds against the state's scale).
 * State model: adds (lsbFrac, msbFrac) interactions.
 * Training families: uniform-dup sweeps, Zipf 0.5/1.5, one-hot(k=10, f=0.9).
 * HELD-OUT eval families: LC minrand, Zipf 1.0, one-hot(k=1, f=0.99) and
 * (k=100, f=0.5) — generalization, not memorization.
 * Reported: mean |relative error| of the complexity estimate per family,
 * and the cardinality estimate it implies (nHat = adds * cHat).
 *
 * Usage: java cardinality.FllComplexity [buckets] [trainInst] [evalInst] [threads]
 *
 * @author Amber (idea: Brian)
 * @date July 2026
 */
public class FllComplexity {

	static final long GOLD=0x9E3779B97F4A7C15L;
	static final double LOG2=Math.log(2.0);

	/** Stream spec: type 0=uniformDup, 1=zipf, 2=onehot, 3=lcMinrand. */
	static final class Spec {
		final int type;
		final long param1;    // uniformDup: n distinct; zipf/lc: pool; onehot: k hot keys
		final double param2;  // uniformDup: dup; zipf: s; onehot: hot fraction f
		final long draws;
		Spec(int t, long p1, double p2, long d){type=t; param1=p1; param2=p2; draws=d;}
		String name(){
			switch(type){
				case 0: return String.format("unif_dup%.1f", param2);
				case 1: return String.format("Zipf_%.1f", param2);
				case 2: return String.format("onehot_k%d_f%.2f", param1, param2);
				default: return "LC_minrand";
			}
		}
	}

	/** Runs a stream; returns {distinct, adds, meanExp, lsb, msb}. */
	static double[] runStream(Spec sp, int inst, int B){
		final FutureLogLog2 f=new FutureLogLog2(B, 31, -1, 0);
		final java.util.Random seedr=new java.util.Random(808000L+inst*104729L+sp.name().hashCode());
		final long seedBase=seedr.nextLong();
		final java.util.SplittableRandom rng=new java.util.SplittableRandom(seedr.nextLong());
		long distinct=0;
		if(sp.type==0){
			final long n=sp.param1;
			long added=0;
			outer:
			while(true){
				for(long i=0; i<n; i++){
					f.add(seedBase+i*GOLD);
					if(++added>=sp.draws){break outer;}
				}
			}
			distinct=Math.min(n, sp.draws);
		}else if(sp.type==1){
			final int pool=(int)sp.param1;
			final double[] cum=new double[pool];
			double t=0;
			for(int i=0; i<pool; i++){t+=Math.pow(i+1, -sp.param2); cum[i]=t;}
			for(int i=0; i<pool; i++){cum[i]/=t;}
			final boolean[] seen=new boolean[pool];
			for(long d=0; d<sp.draws; d++){
				final double u=rng.nextDouble();
				int lo=0, hi=pool-1;
				while(lo<hi){
					final int mid=(lo+hi)>>>1;
					if(cum[mid]<u){lo=mid+1;}else{hi=mid;}
				}
				if(!seen[lo]){seen[lo]=true; distinct++;}
				f.add(seedBase+lo*GOLD);
			}
		}else if(sp.type==2){
			// one-hot: with prob f draw one of k hot keys; else a FRESH unique key
			final int k=(int)sp.param1;
			final double hf=sp.param2;
			final boolean[] hotSeen=new boolean[k];
			long fresh=0;
			for(long d=0; d<sp.draws; d++){
				if(rng.nextDouble()<hf){
					final int h=rng.nextInt(k);
					if(!hotSeen[h]){hotSeen[h]=true; distinct++;}
					f.add(seedBase+h*GOLD);
				}else{
					fresh++;
					distinct++;
					f.add(seedBase+(k+fresh)*GOLD);
				}
			}
		}else{
			final int pool=(int)sp.param1;
			final boolean[] seen=new boolean[pool];
			for(long d=0; d<sp.draws; d++){
				final int pos=Math.min(rng.nextInt(pool), rng.nextInt(pool));
				if(!seen[pos]){seen[pos]=true; distinct++;}
				f.add(seedBase+pos*GOLD);
			}
		}
		final double[] s=FllDupCal.summarize(f);
		// Brian's counter observables, normalized per add
		final double d=sp.draws;
		// Per-word 2-bit saturating counter stats: state fractions + the mean
		// localExp deviation of saturated words (are the hot words high or low?)
		final int nw=f.getNumWords();
		final int[] cc=new int[4];
		double expAll=0, expSat=0;
		for(int w=0; w<nw; w++){
			final int c=(f.ctrX==null) ? 0 : f.ctrX[0][w];
			cc[c]++;
			final int le=(f.getWord(w)>>>12)&0xF;
			expAll+=le;
			if(c==3){expSat+=le;}
		}
		expAll/=nw;
		final double dExpSat=(cc[3]>0) ? expSat/cc[3]-expAll : 0;
		return new double[]{distinct, sp.draws, s[0], s[1], s[2],
			f.ctrRehit/d, f.ctrBelow/d, f.ctrHigh/d,
			(double)cc[1]/nw, (double)cc[2]/nw, (double)cc[3]/nw, dExpSat};
	}

	static double[] baseFeatures(double r){return new double[]{1, r, r*r};}

	static double[] linFeatures(double r, double lsb, double msb){
		return new double[]{1, r, lsb, msb};
	}

	static double[] stateFeatures(double r, double lsb, double msb){
		return new double[]{1, r, r*r, lsb, msb, lsb*lsb, msb*msb, lsb*msb, r*lsb, r*msb};
	}

	static double[] cubicFeatures(double r, double lsb, double msb){
		return new double[]{1, r, lsb, msb, r*r, lsb*lsb, msb*msb, r*lsb, r*msb, lsb*msb,
			r*r*r, lsb*lsb*lsb, msb*msb*msb, r*r*lsb, r*r*msb, lsb*lsb*r, lsb*lsb*msb,
			msb*msb*r, msb*msb*lsb, r*lsb*msb};
	}

	/** Brian's counters (log-domain rates) + state quadratic. */
	static double[] ctrFeatures(double r, double lsb, double msb,
			double rehit, double below, double high){
		final double lre=Math.log1p(1000*rehit), lbe=Math.log1p(1000*below),
			lhi=Math.log1p(1000*high);
		return new double[]{1, r, lsb, msb, r*r, lsb*msb,
			lre, lbe, lhi, lre*lre, lre*r, lbe*r, lre*lbe};
	}

	/** Brian's 2-bit saturating counter histogram features. */
	static double[] ctr2Features(double r, double lsb, double msb,
			double f1, double f2, double f3, double dExpSat){
		return new double[]{1, r, lsb, msb, f1, f2, f3, dExpSat,
			r*r, lsb*msb, f1*f1, f3*f3, f1*f3, r*f3, msb*f3, dExpSat*f3};
	}

	/** Tiny MLP: D -> H tanh -> 1 linear, full-batch Adam, weight decay. */
	static final class TinyNet {
		final int H;
		final int D;
		final double[][] w1;
		final double[] b1, w2;
		double b2;
		final double[] mean, sd;

		TinyNet(java.util.List<double[]> rawX, java.util.List<Double> rawY){
			this(rawX, rawY, 12);
		}

		TinyNet(java.util.List<double[]> rawX, java.util.List<Double> rawY, int hidden){
			H=hidden;
			D=rawX.get(0).length;
			mean=new double[D]; sd=new double[D];
			final int n=rawX.size();
			for(int d=0; d<D; d++){
				double m=0, s=0;
				for(double[] x : rawX){m+=x[d];}
				m/=n;
				for(double[] x : rawX){s+=(x[d]-m)*(x[d]-m);}
				mean[d]=m; sd[d]=Math.max(1e-9, Math.sqrt(s/n));
			}
			final java.util.Random rnd=new java.util.Random(42);
			w1=new double[H][D]; b1=new double[H]; w2=new double[H];
			for(int h=0; h<H; h++){
				for(int d=0; d<D; d++){w1[h][d]=rnd.nextGaussian()*0.5;}
				w2[h]=rnd.nextGaussian()*0.5;
			}
			// Adam, full batch
			final double lr=0.01, be1=0.9, be2=0.999, eps=1e-8, wd=1e-4;
			final int P=H*D+H+H+1;
			final double[] mAd=new double[P], vAd=new double[P];
			final double[][] xn=new double[n][D];
			for(int i=0; i<n; i++){
				for(int d=0; d<D; d++){xn[i][d]=(rawX.get(i)[d]-mean[d])/sd[d];}
			}
			final double[] y=new double[n];
			for(int i=0; i<n; i++){y[i]=rawY.get(i);}
			for(int ep=1; ep<=6000; ep++){
				final double[] g=new double[P];
				for(int i=0; i<n; i++){
					final double[] a=new double[H];
					double out=b2;
					for(int h=0; h<H; h++){
						double z=b1[h];
						for(int d=0; d<D; d++){z+=w1[h][d]*xn[i][d];}
						a[h]=Math.tanh(z);
						out+=w2[h]*a[h];
					}
					final double e=2.0*(out-y[i])/n;
					int p=0;
					for(int h=0; h<H; h++){
						final double dh=e*w2[h]*(1-a[h]*a[h]);
						for(int d=0; d<D; d++){g[p++]+=dh*xn[i][d];}
					}
					for(int h=0; h<H; h++){g[p++]+=e*w2[h]*(1-a[h]*a[h]);}   // b1
					for(int h=0; h<H; h++){g[p++]+=e*a[h];}                  // w2
					g[p]+=e;                                                 // b2
				}
				int p=0;
				for(int h=0; h<H; h++){
					for(int d=0; d<D; d++){
						w1[h][d]=adam(w1[h][d], g[p]+wd*w1[h][d], mAd, vAd, p, ep, lr, be1, be2, eps); p++;
					}
				}
				for(int h=0; h<H; h++){b1[h]=adam(b1[h], g[p], mAd, vAd, p, ep, lr, be1, be2, eps); p++;}
				for(int h=0; h<H; h++){w2[h]=adam(w2[h], g[p]+wd*w2[h], mAd, vAd, p, ep, lr, be1, be2, eps); p++;}
				b2=adam(b2, g[p], mAd, vAd, p, ep, lr, be1, be2, eps);
			}
		}

		static double adam(double w, double g, double[] m, double[] v, int p, int t,
				double lr, double b1, double b2, double eps){
			m[p]=b1*m[p]+(1-b1)*g;
			v[p]=b2*v[p]+(1-b2)*g*g;
			final double mh=m[p]/(1-Math.pow(b1, t));
			final double vh=v[p]/(1-Math.pow(b2, t));
			return w-lr*mh/(Math.sqrt(vh)+eps);
		}

		double predict(double... raw){
			final double[] x=new double[D];
			for(int d=0; d<D; d++){x[d]=(raw[d]-mean[d])/sd[d];}
			double out=b2;
			for(int h=0; h<H; h++){
				double z=b1[h];
				for(int d=0; d<D; d++){z+=w1[h][d]*x[d];}
				out+=w2[h]*Math.tanh(z);
			}
			return out;
		}
	}

	public static void main(String[] args) throws Exception{
		final int B=(args.length>0) ? Integer.parseInt(args[0]) : 2048;
		final int trainInst=(args.length>1) ? Integer.parseInt(args[1]) : 8;
		final int evalInst=(args.length>2) ? Integer.parseInt(args[2]) : 16;
		final int threads=(args.length>3) ? Integer.parseInt(args[3])
			: Runtime.getRuntime().availableProcessors();
		CardinalityTracker.clampToAdded=false;

		// ---- Training families ----
		final java.util.ArrayList<Spec> train=new java.util.ArrayList<>();
		for(long n : new long[]{49L*B, 98L*B, 195L*B, 390L*B}){
			for(double dup : new double[]{1, 1.5, 2, 3, 4, 8, 16}){
				if(n*dup<=30_000_000L){train.add(new Spec(0, n, dup, (long)(n*dup)));}
			}
		}
		for(double s : new double[]{0.5, 1.5}){
			for(long pool : new long[]{98L*B, 488L*B}){
				for(int mult : new int[]{4, 16}){
					train.add(new Spec(1, pool, s, mult*pool));
				}
			}
		}
		for(long pool : new long[]{98L*B, 488L*B}){
			train.add(new Spec(2, 10, 0.9, 10*pool));
			train.add(new Spec(2, 10, 0.5, 2*pool));
		}

		// ---- Held-out eval families ----
		final java.util.ArrayList<Spec> eval=new java.util.ArrayList<>();
		eval.add(new Spec(3, 244L*B, 0, 976L*B));                 // LC minrand
		eval.add(new Spec(1, 98L*B, 1.0, 4L*98*B));               // Zipf 1.0
		eval.add(new Spec(1, 488L*B, 1.0, 16L*488*B));
		eval.add(new Spec(2, 1, 0.99, 100L*98*B));                // ONE key, 99% of adds
		eval.add(new Spec(2, 100, 0.5, 2L*195*B));                // 100 hot keys, mild
		eval.add(new Spec(0, 146L*B, 1.0, 146L*B));               // unique control
		eval.add(new Spec(0, 98L*B, 2.0, 2L*98*B));               // uniform 2x control

		// ---- Collect training data (raw triples; bases derived) ----
		final java.util.List<double[]> RAW=java.util.Collections.synchronizedList(new java.util.ArrayList<>());
		final java.util.List<Double> Y=java.util.Collections.synchronizedList(new java.util.ArrayList<>());
		final java.util.ArrayList<int[]> jobs=new java.util.ArrayList<>();
		for(int si=0; si<train.size(); si++){
			for(int i=0; i<trainInst; i++){jobs.add(new int[]{si, i});}
		}
		FllSkewTest.runPool(threads, jobs.size(), j->{
			final int[] job=jobs.get(j);
			final double[] o=runStream(train.get(job[0]), job[1], B);
			final double r=Math.log(o[1])/LOG2-o[2];
			synchronized(RAW){
				RAW.add(new double[]{r, o[3], o[4], o[5], o[6], o[7], o[8], o[9], o[10], o[11]});
				Y.add(Math.log(o[0]/o[1])/LOG2);   // y = log2(complexity), <= 0
			}
		});
		System.err.println("training samples: "+RAW.size());
		final java.util.List<double[]> XB=new java.util.ArrayList<>();
		final java.util.List<double[]> XL=new java.util.ArrayList<>();
		final java.util.List<double[]> XS=new java.util.ArrayList<>();
		final java.util.List<double[]> XC=new java.util.ArrayList<>();
		for(double[] o : RAW){
			XB.add(baseFeatures(o[0]));
			XL.add(linFeatures(o[0], o[1], o[2]));
			XS.add(stateFeatures(o[0], o[1], o[2]));
			XC.add(cubicFeatures(o[0], o[1], o[2]));
		}
		final double[] betaB=FllDupCal.lsq(XB, Y);
		final double[] betaL=FllDupCal.lsq(XL, Y);
		final double[] betaS=FllDupCal.lsq(XS, Y);
		final double[] betaC=FllDupCal.lsq(XC, Y);
		final java.util.List<double[]> XCTR=new java.util.ArrayList<>();
		final java.util.List<double[]> RAW3=new java.util.ArrayList<>();
		final java.util.List<double[]> RAW6=new java.util.ArrayList<>();
		for(double[] o : RAW){
			XCTR.add(ctrFeatures(o[0], o[1], o[2], o[3], o[4], o[5]));
			RAW3.add(new double[]{o[0], o[1], o[2]});
			RAW6.add(new double[]{o[0], o[1], o[2],
				Math.log1p(1000*o[3]), Math.log1p(1000*o[4]), Math.log1p(1000*o[5])});
		}
		final double[] betaCtr=FllDupCal.lsq(XCTR, Y);
		final java.util.List<double[]> XC2=new java.util.ArrayList<>();
		final java.util.List<double[]> RAW10=new java.util.ArrayList<>();
		for(double[] o : RAW){
			XC2.add(ctr2Features(o[0], o[1], o[2], o[6], o[7], o[8], o[9]));
			RAW10.add(new double[]{o[0], o[1], o[2],
				Math.log1p(1000*o[3]), Math.log1p(1000*o[4]), Math.log1p(1000*o[5]),
				o[6], o[7], o[8], o[9]});
		}
		final double[] betaC2=FllDupCal.lsq(XC2, Y);
		final TinyNet net=new TinyNet(RAW3, Y);
		final TinyNet bigNet=new TinyNet(RAW10, Y, 24);
		System.err.println("models fitted (LSQ + ctr2 quad + NN3 + bigNN 10-24-1)");

		// ---- Evaluate on held-out families ----
		System.out.println("#family\tdraws\ttrueComplexity\tbase_r_only%\tquad%\tctr2Quad%\ttinyNN3%\tbigNN10%");
		final double[][] tot=new double[5][2];
		for(final Spec sp : eval){
			final double[] acc=new double[7];
			final Object lock=new Object();
			FllSkewTest.runPool(threads, evalInst, inst->{
				final double[] o=runStream(sp, 7777+inst, B);
				final double c=o[0]/o[1];
				final double r=Math.log(o[1])/LOG2-o[2];
				final double[] est={
					Math.pow(2.0, FllDupCal.dot(baseFeatures(r), betaB)),
					Math.pow(2.0, FllDupCal.dot(stateFeatures(r, o[3], o[4]), betaS)),
					Math.pow(2.0, FllDupCal.dot(ctr2Features(r, o[3], o[4], o[8], o[9], o[10], o[11]), betaC2)),
					Math.pow(2.0, net.predict(r, o[3], o[4])),
					Math.pow(2.0, bigNet.predict(r, o[3], o[4],
						Math.log1p(1000*o[5]), Math.log1p(1000*o[6]), Math.log1p(1000*o[7]),
						o[8], o[9], o[10], o[11]))};
				synchronized(lock){
					for(int e=0; e<5; e++){acc[e]+=Math.abs(est[e]-c)/c;}
					acc[5]+=c;
				}
			});
			final StringBuilder sb=new StringBuilder(String.format("%s\t%d\t%.4f",
				sp.name(), sp.draws, acc[5]/evalInst));
			for(int e=0; e<5; e++){
				final double err=100*acc[e]/evalInst;
				sb.append(String.format("\t%.2f", err));
				tot[e][0]+=err; tot[e][1]++;
			}
			System.out.println(sb);
		}
		final StringBuilder sb=new StringBuilder("#MEAN_OVER_FAMILIES\t-\t-");
		for(int e=0; e<5; e++){sb.append(String.format("\t%.2f", tot[e][0]/tot[e][1]));}
		System.out.println(sb);
	}
}
