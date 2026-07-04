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
		final int D=NF, H;
		final double[][] w1;
		final double[] b1, w2, mean, sd;
		double b2;

		/** Raw-field constructor for load(). */
		BigNet(int hidden){
			H=hidden;
			w1=new double[H][D]; b1=new double[H]; w2=new double[H];
			mean=new double[D]; sd=new double[D];
		}

		void save(String path) throws Exception{
			final StringBuilder sb=new StringBuilder();
			sb.append(H).append('\n');
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
			final BigNet n=new BigNet(Integer.parseInt(L.get(0).trim()));
			final double[][] heads={n.mean, n.sd, n.b1, n.w2};
			for(int i=0; i<4; i++){
				final String[] f=L.get(1+i).trim().split("\t");
				for(int j=0; j<heads[i].length; j++){heads[i][j]=Double.parseDouble(f[j]);}
			}
			for(int h=0; h<n.H; h++){
				final String[] f=L.get(5+h).trim().split("\t");
				for(int j=0; j<NF; j++){n.w1[h][j]=Double.parseDouble(f[j]);}
			}
			n.b2=Double.parseDouble(L.get(5+n.H).trim());
			return n;
		}

		BigNet(java.util.List<double[]> data, int hidden, int epochs){
			H=hidden;
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
			return out;
		}
	}
}
