package cardinality;

/**
 * ComplexityFarm: Brian's prescription — millions of training vectors, weights
 * in the hundreds.  Streams from a wide family grid (uniform-dup, LC pools,
 * Zipf, one-hot) emit a sample at every geometric distinct-checkpoint, so each
 * stream yields ~10 samples.  Features (16): r=log2(adds)-meanExp, lsbFrac,
 * msbFrac, dExpSat, and the occupancy vectors (f1,f2,f3) of all four 2-bit
 * saturating counter semantics (see FutureLogLog2.ctrX).  Target:
 * y=log2(distinct/adds).  Sample rows carry y as the LAST column — one atomic
 * append, no pairing races.
 * Models: quadratic LSQ on (r,lsb,msb) [the old champion, same data] and an
 * H-hidden tanh net on all 16 inputs, mini-batch Adam.
 * Eval: held-out families (never trained): LC_minrand iter4 at an untrained
 * pool, Zipf 1.0, onehot(k=1,f=0.99), onehot(k=100,f=0.5), unif dup 1 and 2.
 *
 * Usage: java cardinality.ComplexityFarm [buckets] [streamsPerSpec] [evalInst] [threads] [hidden]
 *
 * @author Amber (design objective: Brian)
 * @date July 2026
 */
public class ComplexityFarm {

	static final long GOLD=0x9E3779B97F4A7C15L;
	static final double LOG2=Math.log(2.0);
	static final int NF=16; // feature count; row = NF features + y

	/*------------------------ stream generation ------------------------*/

	/** type: 0=unifDup(n=p1, dup=p2) 1=zipf(pool=p1, s=p2, mult=p3)
	 *  2=onehot(k=p1, f=p2, targetFresh=p3) 3=lcPool(pool=p1, iter=p2) */
	static final class Spec {
		final int type; final long p1; final double p2, p3;
		Spec(int t, long a, double b, double c){type=t; p1=a; p2=b; p3=c;}
	}

	interface SampleSink {void emit(FutureLogLog2 f, long distinct, long adds);}

	static void runStream(Spec sp, long seed, int B, SampleSink sink){
		final FutureLogLog2 f=new FutureLogLog2(B, 31, -1, 0);
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
					if(i<n && adds<=n){distinct=adds;}
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
				if(distinct>=nextCp){sink.emit(f, distinct, adds); nextCp=(long)(nextCp*1.5);}
			}
			sink.emit(f, distinct, adds);
		}
	}

	/** 16 features from sketch state at snapshot time. */
	static double[] features(FutureLogLog2 f, long adds){
		final double[] s=FllDupCal.summarize(f);
		final double r=Math.log((double)adds)/LOG2-s[0];
		final int nw=f.getNumWords();
		final double[] out=new double[NF];
		out[0]=r; out[1]=s[1]; out[2]=s[2];
		double expAll=0, expSat=0;
		int satCnt=0;
		for(int w=0; w<nw; w++){
			final int le=(f.getWord(w)>>>12)&0xF;
			expAll+=le;
			if(f.ctrX!=null && f.ctrX[1][w]==3){expSat+=le; satCnt++;}
		}
		out[3]=(satCnt>0) ? expSat/satCnt-expAll/nw : 0;
		for(int sem=0; sem<4; sem++){
			final int[] cc=new int[4];
			for(int w=0; w<nw; w++){
				cc[(f.ctrX==null) ? 0 : f.ctrX[sem][w]]++;
			}
			for(int k=0; k<3; k++){out[4+sem*3+k]=(double)cc[k+1]/nw;}
		}
		return out;
	}

	/*------------------------------ main ------------------------------*/

	public static void main(String[] args) throws Exception{
		final int B=(args.length>0) ? Integer.parseInt(args[0]) : 2048;
		final int perSpec=(args.length>1) ? Integer.parseInt(args[1]) : 64;
		final int evalInst=(args.length>2) ? Integer.parseInt(args[2]) : 64;
		final int threads=(args.length>3) ? Integer.parseInt(args[3])
			: Runtime.getRuntime().availableProcessors();
		final int hidden=(args.length>4) ? Integer.parseInt(args[4]) : 32;
		CardinalityTracker.clampToAdded=false;

		// ---- Randomized CONTINUOUS training distribution ----
		// Discrete family grids make held-out shapes EXTRAPOLATIONS: the grid farm
		// measured bigNN 46% on onehot_k1 despite 756k samples while the rigid
		// quadratic held 4-5%.  Drawing every generator parameter from continuous
		// priors makes realistic shapes INTERIOR points of the training measure.
		final int totalStreams=perSpec*118;  // arg keeps old scale: perSpec~streams/118
		System.err.println("randomized streams="+totalStreams);

		// ---- Farm the samples (row = 16 features + y, single atomic append) ----
		final java.util.List<double[]> DATA=java.util.Collections.synchronizedList(
			new java.util.ArrayList<>());
		FllSkewTest.runPool(threads, totalStreams, j->{
			final java.util.SplittableRandom prng=new java.util.SplittableRandom(
				new java.util.Random(9090L+j).nextLong());
			final Spec sp;
			switch(prng.nextInt(4)){
				case 0: {  // uniform dup: n logU[24B,1200B], dup logU[1,32]
					final long n=(long)(24L*B*Math.exp(prng.nextDouble()*Math.log(50)));
					double dup=Math.exp(prng.nextDouble()*Math.log(32));
					if(n*dup>40_000_000L){dup=40_000_000.0/n;}
					sp=new Spec(0, n, dup, 0); break;}
				case 1: {  // zipf: s U[0.3,2.2], pool logU[40B,1000B], mult logU[2,20]
					final double zs=0.3+prng.nextDouble()*1.9;
					final long pool=(long)(40L*B*Math.exp(prng.nextDouble()*Math.log(25)));
					double mult=Math.exp(Math.log(2)+prng.nextDouble()*Math.log(10));
					if(pool*mult>40_000_000L){mult=40_000_000.0/pool;}
					sp=new Spec(1, pool, zs, mult); break;}
				case 2: {  // onehot: k logU[1,1000], f U[0.3,0.995], fresh logU[30B,600B]
					final long k=Math.max(1, (long)Math.exp(prng.nextDouble()*Math.log(1000)));
					final double hf=0.3+prng.nextDouble()*0.695;
					final long fresh=(long)(30L*B*Math.exp(prng.nextDouble()*Math.log(20)));
					sp=new Spec(2, k, hf, fresh); break;}
				default: {  // lcPool: pool logU[30B,2500B], iter U[1.2,10]
					final long pool=(long)(30L*B*Math.exp(prng.nextDouble()*Math.log(83)));
					double iter=1.2+prng.nextDouble()*8.8;
					if(pool*iter>40_000_000L){iter=40_000_000.0/pool;}
					sp=new Spec(3, pool, iter, 0); break;}
			}
			final long seed=101000L+j*104729L;
			runStream(sp, seed, B, (f, distinct, adds)->{
				if(distinct<24L*B){return;}
				final double[] row=new double[NF+1];
				System.arraycopy(features(f, adds), 0, row, 0, NF);
				row[NF]=Math.log((double)distinct/adds)/LOG2;
				DATA.add(row);
			});
		});
		System.err.println("samples="+DATA.size());

		// ---- Optional dump for Brian's trainer (train.sh / ml.Trainer) ----
		// Format: '#dims\t16\t1' then tab-delimited floats, target normalized
		// to [0,1] via y' = (y+13)/13  (y = log2 complexity in [-13, 0]).
		final String dumpPath=(args.length>5) ? args[5] : null;
		if(dumpPath!=null){
			final StringBuilder sb=new StringBuilder("#dims\t"+NF+"\t1\n");
			for(double[] row : DATA){
				for(int i=0; i<NF; i++){sb.append(String.format("%.6f", row[i])).append('\t');}
				final double yn=Math.max(0, Math.min(1, (row[NF]+13)/13));
				sb.append(String.format("%.6f", yn)).append('\n');
			}
			java.nio.file.Files.write(java.nio.file.Paths.get(dumpPath),
				sb.toString().getBytes());
			System.err.println("dumped "+DATA.size()+" vectors to "+dumpPath);
		}

		// ---- Baseline: quadratic on (r,lsb,msb), same data ----
		final java.util.List<double[]> XQ=new java.util.ArrayList<>();
		final java.util.List<Double> Y=new java.util.ArrayList<>();
		for(double[] row : DATA){
			XQ.add(FllComplexity.stateFeatures(row[0], row[1], row[2]));
			Y.add(row[NF]);
		}
		final double[] betaQ=FllDupCal.lsq(XQ, Y);

		// ---- Big net: 16 -> hidden -> 1, mini-batch Adam ----
		final BigNet net=new BigNet(DATA, hidden, 40);
		System.err.println("net trained: "+NF+"-"+hidden+"-1, "
			+((NF+2)*hidden+1)+" params, "+DATA.size()+" samples");
		if(args.length>6){
			net.save(args[6]);
			System.err.println("net saved to "+args[6]);
		}

		// ---- Held-out eval ----
		final Spec[] evals={
			new Spec(3, 244L*B, 4, 0),          // LC_minrand iter4 (untrained pool+iter)
			new Spec(1, 98L*B, 1.0, 4),         // Zipf 1.0 (untrained s)
			new Spec(1, 488L*B, 1.0, 16),
			new Spec(2, 1, 0.99, 196L*B),       // one key, 99% of adds
			new Spec(2, 100, 0.5, 390L*B),      // 100 hot keys, mild
			new Spec(0, 146L*B, 1.0, 0),        // unique control
			new Spec(0, 98L*B, 2.0, 0),         // uniform 2x control
		};
		final String[] evalNames={"LC_minrand", "Zipf1.0_small", "Zipf1.0_big",
			"onehot_k1_f0.99", "onehot_k100_f0.5", "unif_dup1", "unif_dup2"};
		System.out.println("#family\ttrueComplexity\tquad%\tbigNN%");
		double tq=0, tn=0;
		for(int e=0; e<evals.length; e++){
			final Spec sp=evals[e];
			final double[] acc=new double[3];
			final Object lock=new Object();
			final java.util.List<String> evalRows=(dumpPath!=null)
				? java.util.Collections.synchronizedList(new java.util.ArrayList<>()) : null;
			FllSkewTest.runPool(threads, evalInst, inst->{
				final long seed=717000L+inst*104729L;
				final double[][] last=new double[1][];
				final long[] fin=new long[2];
				runStream(sp, seed, B, (f, distinct, adds)->{
					last[0]=features(f, adds);
					fin[0]=distinct; fin[1]=adds;
				});
				final double c=(double)fin[0]/fin[1];
				if(evalRows!=null){
					final StringBuilder rb=new StringBuilder();
					for(int i=0; i<NF; i++){rb.append(String.format("%.6f", last[0][i])).append('\t');}
					rb.append(String.format("%.6f",
						Math.max(0, Math.min(1, (Math.log(c)/LOG2+13)/13))));
					evalRows.add(rb.toString());
				}
				final double cQ=Math.pow(2.0, FllDupCal.dot(
					FllComplexity.stateFeatures(last[0][0], last[0][1], last[0][2]), betaQ));
				final double cN=Math.pow(2.0, net.predict(last[0]));
				synchronized(lock){
					acc[0]+=Math.abs(cQ-c)/c;
					acc[1]+=Math.abs(cN-c)/c;
					acc[2]+=c;
				}
			});
			final double eq=100*acc[0]/evalInst, en=100*acc[1]/evalInst;
			tq+=eq; tn+=en;
			System.out.println(String.format("%s\t%.4f\t%.2f\t%.2f",
				evalNames[e], acc[2]/evalInst, eq, en));
			if(evalRows!=null){
				final StringBuilder eb=new StringBuilder("#dims\t"+NF+"\t1\n");
				for(String row : evalRows){eb.append(row).append('\n');}
				java.nio.file.Files.write(java.nio.file.Paths.get(
					dumpPath+".eval_"+evalNames[e]+".tsv"), eb.toString().getBytes());
			}
		}
		System.out.println(String.format("#MEAN\t-\t%.2f\t%.2f",
			tq/evals.length, tn/evals.length));
	}

	/*------------------------- mini-batch net -------------------------*/

	static final class BigNet {
		final int D, H;
		final double[][] w1;
		final double[] b1, w2, mean, sd;
		double b2;
		/** Output activation: 0=linear (legacy), 1=rslog (sign(z)*log(1+|z|)).
		 * rslog matches BBTools ml.Function.RSLOG so rslog nets export exactly
		 * to BBNet format (see BigNetToBBNet). */
		final int act;

		/** Raw-field constructor for load(). */
		BigNet(int hidden){this(hidden, 0, NF);}

		BigNet(int hidden, int outAct){this(hidden, outAct, NF);}

		BigNet(int hidden, int outAct, int inputs){
			assert(hidden>0) : "hidden="+hidden;
			assert(outAct==0 || outAct==1) : "act="+outAct+" not in {0=linear,1=rslog}";
			assert(inputs>0) : "inputs="+inputs;
			H=hidden; act=outAct; D=inputs;
			w1=new double[H][D]; b1=new double[H]; w2=new double[H];
			mean=new double[D]; sd=new double[D];
		}

		void save(String path) throws Exception{
			final StringBuilder sb=new StringBuilder();
			sb.append(H).append('\t').append(act).append('\t').append(D).append('\n');
			for(double[] arr : new double[][]{mean, sd, b1, w2}){
				for(double v : arr){sb.append(v).append('\t');}
				sb.append('\n');
			}
			for(double[] row : w1){
				for(double v : row){sb.append(v).append('\t');}
				sb.append('\n');
			}
			sb.append(b2).append('\n');
			java.nio.file.Files.write(java.nio.file.Paths.get(path), sb.toString().getBytes());
		}

		static BigNet load(String path) throws Exception{
			final java.util.List<String> L=java.nio.file.Files.readAllLines(
				java.nio.file.Paths.get(path));
			// header: "H" (legacy, linear) or "H\tact" or "H\tact\tD"
			final String[] head=L.get(0).trim().split("\t");
			final BigNet n=new BigNet(Integer.parseInt(head[0].trim()),
				(head.length>1) ? Integer.parseInt(head[1].trim()) : 0,
				(head.length>2) ? Integer.parseInt(head[2].trim()) : NF);
			final double[][] heads={n.mean, n.sd, n.b1, n.w2};
			for(int i=0; i<4; i++){
				final String[] f=L.get(1+i).trim().split("\t");
				for(int j=0; j<heads[i].length; j++){heads[i][j]=Double.parseDouble(f[j]);}
			}
			for(int h=0; h<n.H; h++){
				final String[] f=L.get(5+h).trim().split("\t");
				for(int j=0; j<n.D; j++){n.w1[h][j]=Double.parseDouble(f[j]);}
			}
			n.b2=Double.parseDouble(L.get(5+n.H).trim());
			return n;
		}

		/** Legacy serial double-math trainer (linear output).  DO NOT alter its
		 * numerics — in-flight cluster experiments depend on byte-identical
		 * behavior for comparability.  New work: use train() below. */
		BigNet(java.util.List<double[]> data, int hidden, int epochs){
			H=hidden; act=0; D=data.get(0).length-1;
			mean=new double[D]; sd=new double[D];
			final int n=data.size();
			for(int d=0; d<D; d++){
				double m=0, s=0;
				for(double[] row : data){m+=row[d];}
				m/=n;
				for(double[] row : data){s+=(row[d]-m)*(row[d]-m);}
				mean[d]=m; sd[d]=Math.max(1e-9, Math.sqrt(s/n));
			}
			final java.util.Random rnd=new java.util.Random(42);
			w1=new double[H][D]; b1=new double[H]; w2=new double[H];
			for(int h=0; h<H; h++){
				for(int d=0; d<D; d++){w1[h][d]=rnd.nextGaussian()*0.3;}
				w2[h]=rnd.nextGaussian()*0.3;
			}
			final int P=H*D+H+H+1;
			final double[] mAd=new double[P], vAd=new double[P];
			final int BATCH=8192;
			final double be1=0.9, be2=0.999, eps=1e-8, wd=1e-5;
			final int[] idx=new int[n];
			for(int i=0; i<n; i++){idx[i]=i;}
			final java.util.Random shuf=new java.util.Random(7);
			int t=0;
			for(int ep=0; ep<epochs; ep++){
				final double lr=0.005*Math.pow(0.93, ep);
				for(int i=n-1; i>0; i--){
					final int j=shuf.nextInt(i+1);
					final int tmp=idx[i]; idx[i]=idx[j]; idx[j]=tmp;
				}
				for(int start=0; start+BATCH<=n; start+=BATCH){
					t++;
					final double[] g=new double[P];
					for(int bi=start; bi<start+BATCH; bi++){
						final double[] row=data.get(idx[bi]);
						final double[] x=new double[D];
						for(int d=0; d<D; d++){x[d]=(row[d]-mean[d])/sd[d];}
						final double[] a=new double[H];
						double out=b2;
						for(int h=0; h<H; h++){
							double z=b1[h];
							for(int d=0; d<D; d++){z+=w1[h][d]*x[d];}
							a[h]=Math.tanh(z);
							out+=w2[h]*a[h];
						}
						final double e=2.0*(out-row[D])/BATCH;
						int p=0;
						for(int h=0; h<H; h++){
							final double dh=e*w2[h]*(1-a[h]*a[h]);
							for(int d=0; d<D; d++){g[p++]+=dh*x[d];}
						}
						for(int h=0; h<H; h++){g[p++]+=e*w2[h]*(1-a[h]*a[h]);}
						for(int h=0; h<H; h++){g[p++]+=e*a[h];}
						g[p]+=e;
					}
					int p=0;
					for(int h=0; h<H; h++){
						for(int d=0; d<D; d++){
							w1[h][d]=FllComplexity.TinyNet.adam(w1[h][d], g[p]+wd*w1[h][d],
								mAd, vAd, p, t, lr, be1, be2, eps); p++;
						}
					}
					for(int h=0; h<H; h++){
						b1[h]=FllComplexity.TinyNet.adam(b1[h], g[p], mAd, vAd, p, t, lr, be1, be2, eps); p++;
					}
					for(int h=0; h<H; h++){
						w2[h]=FllComplexity.TinyNet.adam(w2[h], g[p]+wd*w2[h], mAd, vAd, p, t, lr, be1, be2, eps); p++;
					}
					b2=FllComplexity.TinyNet.adam(b2, g[p], mAd, vAd, p, t, lr, be1, be2, eps);
				}
			}
		}

		double predict(double[] raw){
			double out=b2;
			for(int h=0; h<H; h++){
				double z=b1[h];
				for(int d=0; d<D; d++){z+=w1[h][d]*(raw[d]-mean[d])/sd[d];}
				out+=w2[h]*Math.tanh(z);
			}
			return (act==1) ? rslog(out) : out;
		}

		static double rslog(double z){
			return (z<0) ? -Math.log(1.0-z) : Math.log(1.0+z);
		}

		/*------------------------------------------------------------*/
		/*  Parallel float trainer (Brian's conventions: strictly      */
		/*  floats for training and evaluation; <=8 threads per net;   */
		/*  beyond that, parallel random restarts, keep the best).     */
		/*------------------------------------------------------------*/

		/** Train with up to totalThreads.  R=totalThreads/8 nets train in
		 * parallel from different random starts (8 data-parallel threads
		 * each); the net with the best held-out MSE wins.  R==1 uses all
		 * threads on one net with no holdout (matches legacy data usage). */
		static BigNet train(final java.util.List<double[]> data, final int hidden,
				final int epochs, final int totalThreads, final int outAct,
				final long seed) throws Exception{
			return train(data, hidden, epochs, totalThreads, outAct, seed, 0.005f);
		}

		static BigNet train(final java.util.List<double[]> data, final int hidden,
				final int epochs, final int totalThreads, final int outAct,
				final long seed, final float lr0) throws Exception{
			final int n=data.size();
			assert(n>0) : "Empty training data";
			assert(totalThreads>0) : "totalThreads="+totalThreads;
			final int dim=data.get(0).length-1;
			// One-time conversion to packed float rows [dim features + y].
			final float[][] F=new float[n][];
			for(int i=0; i<n; i++){
				final double[] r=data.get(i);
				assert(r.length==dim+1) : "row "+i+" length "+r.length+" != "+(dim+1);
				final float[] f=new float[dim+1];
				for(int d=0; d<=dim; d++){f[d]=(float)r[d];}
				F[i]=f;
			}
			final int R=Math.max(1, totalThreads/8);
			if(R==1){
				return new BigNet(F, hidden, epochs, Math.min(8, totalThreads),
					outAct, seed, lr0);
			}
			// Strided holdout (2%) for restart selection; rest trains.
			final int nVal=n/50, nTrain=n-nVal;
			final float[][] TR=new float[nTrain][], VA=new float[nVal][];
			for(int i=0, tr=0, va=0; i<n; i++){
				if(i%50==0 && va<nVal){VA[va++]=F[i];}else{TR[tr++]=F[i];}
			}
			final int tpn=Math.max(1, Math.min(8, totalThreads/R));
			final BigNet[] cand=new BigNet[R];
			final Thread[] th=new Thread[R];
			for(int i=0; i<R; i++){
				final int ci=i;
				th[i]=new Thread(()->{
					try{
						cand[ci]=new BigNet(TR, hidden, epochs, tpn, outAct,
							seed+104729L*ci, lr0);
					}catch(Exception e){throw new RuntimeException(e);}
				});
				th[i].start();
			}
			for(Thread t : th){t.join();}
			BigNet best=null; double bestMse=Double.MAX_VALUE;
			for(int i=0; i<R; i++){
				final double mse=cand[i].mseF(VA);
				System.err.println(String.format(
					"restart %d (seed %d): heldout MSE=%.6f", i,
					seed+104729L*i, mse));
				if(mse<bestMse){bestMse=mse; best=cand[i];}
			}
			if(best==null){   // all-NaN MSE (poisoned data); crash loud even under -da —
				// returning null here cost a 1.7h farm run before the NPE (2026-07-05)
				throw new AssertionError("No candidate survived (all heldout MSE NaN?); R="+R);
			}
			return best;
		}

		/** Float-path MSE over packed rows (evaluation in float, per house
		 * convention — matches what an exported BBNet would compute). */
		double mseF(float[][] rows){
			double sum=0;
			for(float[] row : rows){
				float out=(float)b2;
				for(int h=0; h<H; h++){
					float z=(float)b1[h];
					for(int d=0; d<D; d++){
						z+=(float)w1[h][d]*((row[d]-(float)mean[d])/(float)sd[d]);
					}
					out+=(float)w2[h]*(float)Math.tanh(z);
				}
				if(act==1){out=(out<0) ? (float)-Math.log(1f-out) : (float)Math.log(1f+out);}
				final float e=out-row[D];
				sum+=e*(double)e;
			}
			return sum/Math.max(1, rows.length);
		}

		/** Data-parallel float-math trainer: one net, tpn threads splitting
		 * each minibatch.  Per-thread gradient partials are summed in thread
		 * order, so results are deterministic for a fixed tpn. */
		private BigNet(final float[][] rows, final int hidden, final int epochs,
				final int tpn, final int outAct, final long seed, final float lr0)
				throws Exception{
			H=hidden; act=outAct; D=rows[0].length-1;
			mean=new double[D]; sd=new double[D];
			final int n=rows.length;
			assert(n>0) : "Empty training rows";
			for(int d=0; d<D; d++){
				double m=0, s=0;
				for(float[] row : rows){m+=row[d];}
				m/=n;
				for(float[] row : rows){s+=(row[d]-m)*(row[d]-m);}
				mean[d]=m; sd[d]=Math.max(1e-9, Math.sqrt(s/n));
			}
			final float[] fm=new float[D], fs=new float[D];
			for(int d=0; d<D; d++){fm[d]=(float)mean[d]; fs[d]=(float)sd[d];}
			final java.util.Random rnd=new java.util.Random(seed);
			final float[] wf1=new float[H*D], bf1=new float[H], wf2=new float[H];
			for(int h=0; h<H; h++){
				for(int d=0; d<D; d++){wf1[h*D+d]=(float)(rnd.nextGaussian()*0.3);}
				wf2[h]=(float)(rnd.nextGaussian()*0.3);
			}
			float bf2=0;
			final int P=H*D+H+H+1;
			final float[] mAd=new float[P], vAd=new float[P];
			final int BATCH=8192;
			final float be1=0.9f, be2=0.999f, eps=1e-8f, wd=1e-5f;
			final int[] idx=new int[n];
			for(int i=0; i<n; i++){idx[i]=i;}
			final java.util.Random shuf=new java.util.Random(seed*31+7);
			final java.util.concurrent.ExecutorService ex=
				java.util.concurrent.Executors.newFixedThreadPool(tpn);
			try{
				final float[][] gPart=new float[tpn][P];
				final java.util.List<java.util.concurrent.Callable<Void>> tasks=
					new java.util.ArrayList<>(tpn);
				int t=0;
				for(int ep=0; ep<epochs; ep++){
					final float lr=lr0*(float)Math.pow(0.93, ep);
					for(int i=n-1; i>0; i--){
						final int j=shuf.nextInt(i+1);
						final int tmp=idx[i]; idx[i]=idx[j]; idx[j]=tmp;
					}
					for(int start=0; start+BATCH<=n; start+=BATCH){
						t++;
						final int base=start;
						final float bf2c=bf2;
						tasks.clear();
						final int per=BATCH/tpn;
						for(int k=0; k<tpn; k++){
							final int lo=base+k*per;
							final int hi=(k==tpn-1) ? base+BATCH : lo+per;
							final float[] g=gPart[k];
							tasks.add(()->{
								java.util.Arrays.fill(g, 0f);
								final float[] x=new float[D], a=new float[H];
								for(int bi=lo; bi<hi; bi++){
									final float[] row=rows[idx[bi]];
									for(int d=0; d<D; d++){x[d]=(row[d]-fm[d])/fs[d];}
									float zo=bf2c;
									for(int h=0; h<H; h++){
										float z=bf1[h];
										final int off=h*D;
										for(int d=0; d<D; d++){z+=wf1[off+d]*x[d];}
										a[h]=(float)Math.tanh(z);
										zo+=wf2[h]*a[h];
									}
									float out=zo;
									if(outAct==1){
										out=(zo<0) ? (float)-Math.log(1f-zo)
											: (float)Math.log(1f+zo);
									}
									float e=2f*(out-row[D])/BATCH;
									if(outAct==1){e/=(1f+Math.abs(zo));}
									int p=0;
									for(int h=0; h<H; h++){
										final float dh=e*wf2[h]*(1f-a[h]*a[h]);
										for(int d=0; d<D; d++){g[p++]+=dh*x[d];}
									}
									for(int h=0; h<H; h++){
										g[p++]+=e*wf2[h]*(1f-a[h]*a[h]);
									}
									for(int h=0; h<H; h++){g[p++]+=e*a[h];}
									g[p]+=e;
								}
								return null;
							});
						}
						for(java.util.concurrent.Future<Void> f : ex.invokeAll(tasks)){
							f.get();   // rethrow worker exceptions
						}
						final float[] g=gPart[0];
						for(int k=1; k<tpn; k++){
							final float[] gk=gPart[k];
							for(int p=0; p<P; p++){g[p]+=gk[p];}
						}
						int p=0;
						for(int h=0; h<H; h++){
							for(int d=0; d<D; d++){
								wf1[h*D+d]=adamF(wf1[h*D+d], g[p]+wd*wf1[h*D+d],
									mAd, vAd, p, t, lr, be1, be2, eps); p++;
							}
						}
						for(int h=0; h<H; h++){
							bf1[h]=adamF(bf1[h], g[p], mAd, vAd, p, t, lr, be1, be2, eps); p++;
						}
						for(int h=0; h<H; h++){
							wf2[h]=adamF(wf2[h], g[p]+wd*wf2[h], mAd, vAd, p, t, lr, be1, be2, eps); p++;
						}
						bf2=adamF(bf2, g[p], mAd, vAd, p, t, lr, be1, be2, eps);
					}
				}
			}finally{ex.shutdown();}
			w1=new double[H][D]; b1=new double[H]; w2=new double[H];
			for(int h=0; h<H; h++){
				for(int d=0; d<D; d++){w1[h][d]=wf1[h*D+d];}
				b1[h]=bf1[h]; w2[h]=wf2[h];
			}
			b2=bf2;
			assert(Double.isFinite(b2)) : "Training diverged: b2="+b2;
		}

		static float adamF(float x, float g, float[] m, float[] v, int p, int t,
				float lr, float be1, float be2, float eps){
			m[p]=be1*m[p]+(1f-be1)*g;
			v[p]=be2*v[p]+(1f-be2)*g*g;
			final float mh=m[p]/(1f-(float)Math.pow(be1, t));
			final float vh=v[p]/(1f-(float)Math.pow(be2, t));
			return x-lr*mh/((float)Math.sqrt(vh)+eps);
		}
	}
}
