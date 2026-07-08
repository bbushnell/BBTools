package cardinality;

/**
 * Fll53Chron: chronicle complexity estimation (Brian's idea, 2026-07-04).
 * The snapshot window net estimates complexity from the sketch's final state
 * only, discarding the history of how occupancy evolved — but saturation is a
 * one-way door, and the early era (when the sketch could still see individual
 * novelty) carries information the final state has destroyed.
 *
 * Mechanism: while streaming, at every doubling of adds (from 4*B on), the
 * BASE window net is evaluated on the sketch and (log2 adds, yHat) is
 * recorded — a trajectory of at most ~64 doubles per SKETCH (honest memory:
 * bytes per sketch, not per organism).  At estimation time a FUSION net reads
 * 16 trajectory-shape features (current estimate, lagged estimates, slopes,
 * curvature, extremes) and outputs a corrected log2-complexity.  A stream
 * whose complexity is falling (duplicates arriving late) and one that is flat
 * are different animals a snapshot cannot distinguish.
 *
 * The fusion net is trained with BigNet.train (parallel float trainer,
 * rslog output) on randomized continuous shape priors and random sizes,
 * and evaluated on held-out families against the base net side by side.
 *
 * Usage: java cardinality.Fll53Chron [streams] [evalInst] [threads] [hidden]
 *        [baseNet] [netSave] [onlyOrgs]
 *
 * @author Amber (idea: Brian)
 * @date July 2026
 */
public class Fll53Chron {

	static final long GOLD=0x9E3779B97F4A7C15L;
	static final double LOG2=Math.log(2.0);
	static final int WIN=192;
	static final int FNF=36;   // 16 adds-clock + 8 exp-clock + 4 MLE + 8 raw snapshot

	/** Per-stream trajectory recorder on TWO clocks (Brian's design): the
	 * external clock (every doubling of total adds) and the internal clock
	 * (every global-NLZ increment) — totally different patterns.  Promotion
	 * cadence relative to adds is itself a complexity signature.  Read-only
	 * use of the shared base net (thread-safe). */
	static final class Traj {
		final ComplexityFarm.BigNet base;
		final Fll53MLE mle;   // independent complexity witness (may be null)
		final double[] la=new double[64], y=new double[64];       // adds clock
		final double[] yM=new double[64];                          // MLE log2-c
		final double[] laE=new double[40], yE=new double[40];     // exp clock
		int n=0, nE=0, lastG=Integer.MIN_VALUE;
		long next;

		Traj(ComplexityFarm.BigNet base, long firstCk){this(base, null, firstCk);}

		Traj(ComplexityFarm.BigNet base, Fll53MLE mle, long firstCk){
			assert(firstCk>0) : "firstCk="+firstCk;
			this.base=base; this.mle=mle; next=firstCk;
		}

		/** Called after every add; fires either clock when due. */
		void tick(Fll53 f, long adds){
			if(adds>=next){
				if(n<64){
					la[n]=Math.log(adds)/LOG2;
					y[n]=yNow(base, f, adds);
					yM[n]=yMle(mle, f, adds);
					n++;
				}
				next*=2;
				if(lastG==Integer.MIN_VALUE){lastG=f.getGlobalExp();}
			}
			final int g=f.getGlobalExp();
			if(lastG!=Integer.MIN_VALUE && g!=lastG){
				lastG=g;
				if(nE<40){
					laE[nE]=Math.log(adds)/LOG2;
					yE[nE]=yNow(base, f, adds);
					nE++;
				}
			}
		}
	}

	/** MLE-implied log2 complexity, clamped to a sane range; 0-cost if mle==null. */
	static double yMle(Fll53MLE mle, Fll53 f, long adds){
		if(mle==null){return 0;}
		final double c=mle.estimate(f)/Math.max(1.0, adds);
		return Math.max(-8, Math.min(0.5, Math.log(Math.max(1e-3, c))/LOG2));
	}

	/** Mean base-net prediction over all windows (log2-complexity estimate). */
	static double yNow(ComplexityFarm.BigNet base, Fll53 f, long adds){
		double sum=0;
		int k=0;
		for(int w0=0; w0+WIN<=f.getNumOrgs(); w0+=WIN){
			sum+=base.predict(Fll53FarmW.features(f, adds, w0));
			k++;
		}
		assert(k>0) : "Sketch smaller than one window: "+f.getNumOrgs();
		return sum/k;
	}

	/** 24 fusion features: snapshot estimate + adds-clock trajectory shape
	 * (features 0-15) + exp-clock promotion cadence (features 16-23). */
	static double[] fuse(Traj tr, Fll53 f, long adds, int B){
		final double now=yNow(tr.base, f, adds);
		final double la=Math.log(adds)/LOG2;
		final int n=tr.n;
		// Lagged estimates, clamped to the oldest recorded checkpoint.
		final double y1=(n>0) ? tr.y[n-1] : now;
		final double y2=(n>1) ? tr.y[n-2] : y1;
		final double y3=(n>2) ? tr.y[n-3] : y2;
		final double y4=(n>3) ? tr.y[n-4] : y3;
		final double yF=(n>0) ? tr.y[0] : now;
		final double laF=(n>0) ? tr.la[0] : la;
		final double la1=(n>0) ? tr.la[n-1] : la;
		double mn=now, mx=now, mean=now;
		if(n>0){
			mn=mx=mean=0;
			mn=Double.MAX_VALUE; mx=-Double.MAX_VALUE;
			for(int i=0; i<n; i++){
				mn=Math.min(mn, tr.y[i]); mx=Math.max(mx, tr.y[i]); mean+=tr.y[i];
			}
			mean/=n;
		}
		// Exp-clock (promotion) features.
		final int nE=tr.nE;
		final double yE1=(nE>0) ? tr.yE[nE-1] : now;
		final double yE2=(nE>1) ? tr.yE[nE-2] : yE1;
		final double laE1=(nE>0) ? tr.laE[nE-1] : la;
		final double laE2=(nE>1) ? tr.laE[nE-2] : laE1;
		double meanE=now;
		if(nE>0){
			meanE=0;
			for(int i=0; i<nE; i++){meanE+=tr.yE[i];}
			meanE/=nE;
		}
		// Independent witnesses: MLE complexity now + trajectory + marginal rate.
		final double yMnow=yMle(tr.mle, f, adds);
		final double yM1=(n>0) ? tr.yM[n-1] : yMnow;
		// marginal novelty between last checkpoint and now, in log2 space:
		// (nHat_now - nHat_1)/(adds - adds_1), clamped like yMle.
		double mSlope=yMnow;
		if(tr.mle!=null && n>0){
			final double a1=Math.pow(2, tr.la[n-1]);
			if(adds>a1+1){
				final double dn=tr.mle.estimate(f)-Math.pow(2, yM1)*a1;
				mSlope=Math.max(-8, Math.min(0.5,
					Math.log(Math.max(1e-3, dn/(adds-a1)))/LOG2));
			}
		}
		// Raw snapshot summary: mean 16-vector over windows.
		final double[] raw=new double[16];
		{
			int k=0;
			for(int w0=0; w0+WIN<=f.getNumOrgs(); w0+=WIN){
				final double[] wf=Fll53FarmW.features(f, adds, w0);
				for(int i=0; i<16; i++){raw[i]+=wf[i];}
				k++;
			}
			for(int i=0; i<16; i++){raw[i]/=Math.max(1, k);}
		}
		return new double[]{
			la-Math.log(B)/LOG2,   // stream length in sketch units
			now,                   // snapshot estimate (the baseline)
			y1, y2, y3, y4,        // lagged estimates (newest..oldest)
			now-y1,                // slope, last doubling
			y1-y3,                 // slope, older era
			now-2*y1+y2,           // curvature
			la-la1,                // fractional position past last checkpoint
			mn, mx, mean,          // trajectory envelope
			n/32.0,                // trajectory length
			yF,                    // earliest estimate (least saturated era)
			laF-Math.log(B)/LOG2,  // when the record began
			yE1, yE2,              // estimates at the last two promotions
			la-laE1,               // doublings since last promotion
			laE1-laE2,             // last promotion interval (cadence)
			nE/16.0,               // promotion count
			meanE,                 // mean estimate across promotions
			now-yE1,               // drift since last promotion
			f.getGlobalExp()/32.0, // current global exponent
			yMnow,                 // MLE-implied complexity NOW (indep. witness)
			yM1,                   // ...at the last checkpoint
			yMnow-yM1,             // MLE drift
			mSlope,                // marginal MLE novelty rate (local complexity)
			raw[1], raw[2], raw[3],        // antenna tier fills (mean over windows)
			raw[5], raw[6], raw[7],        // counter occupancy c1,c2,c3
			raw[4], raw[15]};              // dExpSat, zeroFrac
	}

	/*--------------- Training-envelope guard (v5a, 2026-07-05) ---------------*/
	/* The fusion net is only trustworthy inside the feature region it trained
	 * on: netchron53d_b33 exploded +60% on a real 1.3B-add stream when
	 * unbounded features (log2 adds/B, trajectory length) left the farm's
	 * range.  The envelope (per-feature min/max over the training data) is
	 * saved beside the net; at inference, out-of-envelope means fall back to
	 * the snapshot blend. */

	/** Per-feature training envelope: [0]=min, [1]=max. */
	static double[][] envelope(java.util.List<double[]> data, int nf){
		final double[] mn=new double[nf], mx=new double[nf];
		java.util.Arrays.fill(mn, Double.MAX_VALUE);
		java.util.Arrays.fill(mx, -Double.MAX_VALUE);
		for(double[] row : data){
			for(int i=0; i<nf; i++){
				mn[i]=Math.min(mn[i], row[i]); mx[i]=Math.max(mx[i], row[i]);
			}
		}
		return new double[][]{mn, mx};
	}

	static void saveEnv(String path, double[][] env) throws Exception{
		try(java.io.PrintWriter pw=new java.io.PrintWriter(path)){
			pw.println("#feature\tmin\tmax");
			for(int i=0; i<env[0].length; i++){
				pw.println(i+"\t"+env[0][i]+"\t"+env[1][i]);
			}
		}
	}

	/** Loads an envelope file; returns null (guard disabled) if absent. */
	static double[][] loadEnv(String path){
		try(java.io.BufferedReader br=new java.io.BufferedReader(
				new java.io.FileReader(path))){
			final java.util.ArrayList<double[]> rows=new java.util.ArrayList<>();
			String line;
			while((line=br.readLine())!=null){
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				final String[] s=line.split("\t");
				rows.add(new double[]{Double.parseDouble(s[1]), Double.parseDouble(s[2])});
			}
			final double[] mn=new double[rows.size()], mx=new double[rows.size()];
			for(int i=0; i<rows.size(); i++){mn[i]=rows.get(i)[0]; mx[i]=rows.get(i)[1];}
			return new double[][]{mn, mx};
		}catch(Exception e){return null;}
	}

	/** True iff every feature lies within [min-m*range, max+m*range]. */
	static boolean inEnvelope(double[][] env, double[] f, double margin){
		if(env==null){return true;}
		assert(env[0].length==f.length) : "envelope "+env[0].length+" vs features "+f.length;
		for(int i=0; i<f.length; i++){
			final double r=Math.max(1e-6, env[1][i]-env[0][i]);
			if(f[i]<env[0][i]-margin*r || f[i]>env[1][i]+margin*r){return false;}
		}
		return true;
	}

	static boolean finite(double[] v){
		for(double x : v){if(!Double.isFinite(x)){return false;}}
		return true;
	}

	interface Sink {void emit(Fll53 f, long distinct, long adds);}

	/** Same stream generators as Fll53FarmW.runStream, plus the chronicle
	 * hook: tr.record fires at every doubling of adds. */
	static void runStream(ComplexityFarm.Spec sp, long seed, int B, Traj tr, Sink sink){
		final Fll53 f=new Fll53(B, 31, -1, 0);
		final java.util.Random seedr=new java.util.Random(seed);
		final long sb=seedr.nextLong();
		final java.util.SplittableRandom rng=new java.util.SplittableRandom(seedr.nextLong());
		long distinct=0, adds=0;
		long nextCp=24L*B;
		if(sp.type==0){
			final long n=sp.p1;
			final long total=(long)(n*sp.p2);
			outer:
			while(true){
				for(long i=0; i<n; i++){
					f.add(sb+i*GOLD);
					adds++;
					tr.tick(f, adds);
					if(adds<=n){distinct=adds;}
					if(distinct>=nextCp){sink.emit(f, Math.min(distinct, n), adds); nextCp=(long)(nextCp*1.5);}
					if(adds>=total){break outer;}
				}
			}
			sink.emit(f, Math.min(n, adds), adds);
		}else if(sp.type==1){
			final int pool=(int)sp.p1;
			final double[] cum=new double[pool];
			double t=0;
			for(int i=0; i<pool; i++){t+=Math.pow(i+1, -sp.p2); cum[i]=t;}
			for(int i=0; i<pool; i++){cum[i]/=t;}
			final boolean[] seen=new boolean[pool];
			final long total=(long)(pool*sp.p3);
			for(long d=0; d<total; d++){
				final double u=rng.nextDouble();
				int lo=0, hi=pool-1;
				while(lo<hi){
					final int mid=(lo+hi)>>>1;
					if(cum[mid]<u){lo=mid+1;}else{hi=mid;}
				}
				if(!seen[lo]){seen[lo]=true; distinct++;}
				f.add(sb+lo*GOLD);
				adds++;
				tr.tick(f, adds);
				if(distinct>=nextCp){sink.emit(f, distinct, adds); nextCp=(long)(nextCp*1.5);}
			}
			sink.emit(f, distinct, adds);
		}else if(sp.type==2){
			final int k=(int)sp.p1;
			final double hf=sp.p2;
			final long targetFresh=(long)sp.p3;
			final boolean[] hotSeen=new boolean[k];
			long fresh=0;
			while(fresh<targetFresh){
				if(rng.nextDouble()<hf){
					final int h=rng.nextInt(k);
					if(!hotSeen[h]){hotSeen[h]=true; distinct++;}
					f.add(sb+h*GOLD);
				}else{
					fresh++; distinct++;
					f.add(sb+(k+fresh)*GOLD);
				}
				adds++;
				tr.tick(f, adds);
				if(distinct>=nextCp){sink.emit(f, distinct, adds); nextCp=(long)(nextCp*1.5);}
			}
			sink.emit(f, distinct, adds);
		}else{
			final int pool=(int)sp.p1;
			final boolean[] seen=new boolean[pool];
			final long total=(long)(pool*sp.p2);
			for(long d=0; d<total; d++){
				final int pos=Math.min(rng.nextInt(pool), rng.nextInt(pool));
				if(!seen[pos]){seen[pos]=true; distinct++;}
				f.add(sb+pos*GOLD);
				adds++;
				tr.tick(f, adds);
				if(distinct>=nextCp){sink.emit(f, distinct, adds); nextCp=(long)(nextCp*1.5);}
			}
			sink.emit(f, distinct, adds);
		}
	}

	public static void main(String[] args) throws Exception{
		final int streams=(args.length>0) ? Integer.parseInt(args[0]) : 400000;
		final int evalInst=(args.length>1) ? Integer.parseInt(args[1]) : 64;
		final int threads=(args.length>2) ? Integer.parseInt(args[2])
			: Runtime.getRuntime().availableProcessors();
		final int hidden=(args.length>3) ? Integer.parseInt(args[3]) : 32;
		final String basePath=(args.length>4) ? args[4] : "netwin53.txt";
		final String netSave=(args.length>5) ? args[5] : null;
		final int onlyOrgs=(args.length>6) ? Integer.parseInt(args[6]) : 0;
		// lcBoost: probability of forcing an lcPool-type (minrand-family) spec,
		// on Brian's call that constant complexity is NOT a corner case.
		final double lcBoost=(args.length>7) ? Double.parseDouble(args[7]) : 0;
		CardinalityTracker.clampToAdded=false;

		final java.util.Map<Integer,Fll53MLE> MLES=new java.util.HashMap<>();
		for(int o : new int[]{192, 384, 768, 1536}){
			try{MLES.put(o, new Fll53MLE("fll53mle_"+o+"o.tsv"));}
			catch(Exception e){System.err.println("no MLE table for "+o+"o");}
		}

		final ComplexityFarm.BigNet base=ComplexityFarm.BigNet.load(basePath);
		System.err.println("base net "+basePath+" loaded (H="+base.H+", act="+base.act+")");

		final java.util.List<double[]> DATA=java.util.Collections.synchronizedList(
			new java.util.ArrayList<>());
		final int[] SIZES=(onlyOrgs>0) ? new int[]{onlyOrgs}
			: new int[]{192, 384, 768, 1536};
		FllSkewTest.runPool(threads, streams, j->{
			final java.util.SplittableRandom prng=new java.util.SplittableRandom(
				new java.util.Random(7474L+j).nextLong());
			final int orgs=SIZES[prng.nextInt(SIZES.length)];
			final int B=orgs*5;
			ComplexityFarm.Spec sp=ComplexityFarmW.randomSpec(prng, B);
			if(lcBoost>0 && prng.nextDouble()<lcBoost){
				final long pool=(long)(30L*B*Math.exp(prng.nextDouble()*Math.log(83)));
				sp=new ComplexityFarm.Spec(3, pool, 1.2+prng.nextDouble()*8.8, 0);
			}
			final Traj tr=new Traj(base, MLES.get(orgs), 4L*B);
			runStream(sp, 808000L+j*104729L, B, tr, (f, distinct, adds)->{
				if(distinct<24L*B){return;}
				final double[] row=new double[FNF+1];
				System.arraycopy(fuse(tr, f, adds, B), 0, row, 0, FNF);
				row[FNF]=Math.log((double)distinct/adds)/LOG2;
				if(!finite(row)){   // crash loud even under -da: NaN here poisons
					throw new AssertionError("non-finite fusion row (WIN/farm mismatch?)");
				}
				DATA.add(row);
			});
		});
		System.err.println("fusion samples="+DATA.size());
		final double[][] ENV=envelope(DATA, FNF);
		if(netSave!=null){saveEnv(netSave+".env", ENV); System.err.println("saved "+netSave+".env");}

		// LINEAR output: the target log2-complexity spans [-6.6, 0], where an
		// rslog final is gradient-starved (needs |z|~e^|y| pre-activation) —
		// measured 2026-07-04: rslog fusion lost to base everywhere, 398% on
		// onehot_k1.  rslog finals suit O(1) target ranges only.
		final ComplexityFarm.BigNet fusion=
			ComplexityFarm.BigNet.train(DATA, hidden, 40, threads, 0, 42);
		System.err.println("fusion net trained (linear output, parallel float)");
		if(netSave!=null){fusion.save(netSave); System.err.println("saved "+netSave);}

		System.out.println("#family\torgs\tk\ttrueComplexity\tbase_cErr%\tchron_cErr%");
		for(final int orgs : SIZES){
			final int B=orgs*5;
			final ComplexityFarm.Spec[] evals={
				new ComplexityFarm.Spec(3, 244L*B, 4, 0),
				new ComplexityFarm.Spec(1, 98L*B, 1.0, 4),
				new ComplexityFarm.Spec(2, 1, 0.99, 196L*B),
				new ComplexityFarm.Spec(2, 100, 0.5, 390L*B),
				new ComplexityFarm.Spec(0, 146L*B, 1.0, 0),
				new ComplexityFarm.Spec(0, 98L*B, 2.0, 0),
			};
			final String[] names={"LC_minrand", "Zipf1.0", "onehot_k1",
				"onehot_k100", "unif_dup1", "unif_dup2"};
			for(int e=0; e<evals.length; e++){
				final ComplexityFarm.Spec sp=evals[e];
				final double[] acc=new double[3];
				final Object lock=new Object();
				FllSkewTest.runPool(threads, evalInst, inst->{
					final long[] fin=new long[2];
					final Fll53[] sk=new Fll53[1];
					final Traj tr=new Traj(base, MLES.get(orgs), 4L*B);
					runStream(sp, 909000L+inst*104729L, B, tr,
						(f, distinct, adds)->{sk[0]=f; fin[0]=distinct; fin[1]=adds;});
					final double c=(double)fin[0]/fin[1];
					final double cBase=Math.pow(2.0, yNow(base, sk[0], fin[1]));
					final double cChron=Math.pow(2.0,
						fusion.predict(fuse(tr, sk[0], fin[1], B)));
					synchronized(lock){
						acc[0]+=Math.abs(cBase-c)/c;
						acc[1]+=Math.abs(cChron-c)/c;
						acc[2]+=c;
					}
				});
				System.out.println(String.format("%s\t%d\t%d\t%.4f\t%.2f\t%.2f",
					names[e], orgs, orgs/WIN, acc[2]/evalInst,
					100*acc[0]/evalInst, 100*acc[1]/evalInst));
			}
		}
	}
}
