package cardinality;

/**
 * Fll53bFarmW: modular window farm for the 21-bit species.  Window = 192
 * organisms = 64 longs = 512 bytes = 960 tails, honest.  Same randomized
 * continuous shape priors and random sketch sizes as Fll52FarmW.
 *
 * Usage: java cardinality.Fll53bFarmW [streams] [evalInst] [threads] [hidden] [netSave]
 *
 * @author Amber (design: Brian)
 * @date July 2026
 */
public class Fll53bFarmW {

	static final long GOLD=0x9E3779B97F4A7C15L;
	static final double LOG2=Math.log(2.0);
	static final int WIN=192;              // organisms per window: 512 bytes (3/long)
	static final int NF=16;

	interface Sink {void emit(Fll53b f, long distinct, long adds);}

	static void runStream(ComplexityFarm.Spec sp, long seed, int B, Sink sink){
		final Fll53b f=new Fll53b(B, 31, -1, 0);
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

	/** 16 window features over organisms [w0, w0+WIN). */
	static double[] features(Fll53b f, long adds, int w0){
		final int nw=f.getNumOrgs();
		final double winAdds=Math.max(1.0, adds*(double)WIN/nw);
		long floorB=0, fut1=0, fut2=0;
		double expSum=0, expSat=0;
		int satCnt=0, zeroCnt=0;
		final int[] ca=new int[4], cb=new int[4];
		for(int w=w0; w<w0+WIN; w++){
			final int org=f.getOrg(w);
			final int le=(org>>>Fll53b.EXP_SHIFT)&0xF;
			expSum+=le;
			final int field=(org>>>Fll53b.FIELD_SHIFT)&Fll53b.FIELD_MASK;
			floorB+=Integer.bitCount(field&0x249);
			fut1+=Integer.bitCount(field&0x492);
			fut2+=Integer.bitCount(field&0x924);
			final int a=f.ctrA(w), b=f.ctrB(w);
			ca[a]++; cb[b]++;
			if(a==3){expSat+=le; satCnt++;}
			if(field==0){zeroCnt++;}
		}
		final double meanExp=f.getGlobalExp()+expSum/WIN;
		final double ff=(double)floorB/(WIN*4), f1f=(double)fut1/(WIN*4), f2f=(double)fut2/(WIN*4);
		final double a1=(double)ca[1]/WIN, a2=(double)ca[2]/WIN, a3=(double)ca[3]/WIN;
		final double b1=(double)cb[1]/WIN, b2=(double)cb[2]/WIN, b3=(double)cb[3]/WIN;
		final double r=Math.log(winAdds)/LOG2-meanExp;
		final double dExpSat=(satCnt>0) ? expSat/satCnt-expSum/WIN : 0;
		return new double[]{r, ff, f1f, f2f, dExpSat, a1, a2, a3, b1, b2, b3,
			(a1+2*a2+3*a3)/3, (b1+2*b2+3*b3)/3, r*a3, f1f*a3,
			(double)zeroCnt/WIN};
	}

	public static void main(String[] args) throws Exception{
		final int streams=(args.length>0) ? Integer.parseInt(args[0]) : 400000;
		final int evalInst=(args.length>1) ? Integer.parseInt(args[1]) : 64;
		final int threads=(args.length>2) ? Integer.parseInt(args[2])
			: Runtime.getRuntime().availableProcessors();
		final int hidden=(args.length>3) ? Integer.parseInt(args[3]) : 32;
		final String netSave=(args.length>4) ? args[4] : null;
		CardinalityTracker.clampToAdded=false;

		for(final boolean flow : new boolean[]{true, false}){
			Fll53b.FLOW=flow;
			System.out.println("### FLOW="+flow+" ###");
			runMode(streams, evalInst, threads, hidden,
				(netSave!=null) ? netSave+(flow ? ".flow" : ".reset") : null);
		}
	}

	static void runMode(int streams, int evalInst, int threads, int hidden,
			String netSave) throws Exception{
		final java.util.List<double[]> DATA=java.util.Collections.synchronizedList(
			new java.util.ArrayList<>());
		final int[] SIZES={192, 384, 768, 1536};   // organisms: 512B..4KB
		FllSkewTest.runPool(threads, streams, j->{
			final java.util.SplittableRandom prng=new java.util.SplittableRandom(
				new java.util.Random(6363L+j).nextLong());
			final int orgs=SIZES[prng.nextInt(SIZES.length)];
			final int B=orgs*5;
			final ComplexityFarm.Spec sp=ComplexityFarmW.randomSpec(prng, B);
			runStream(sp, 707000L+j*104729L, B, (f, distinct, adds)->{
				if(distinct<24L*B){return;}
				final double y=Math.log((double)distinct/adds)/LOG2;
				final int nw=f.getNumOrgs();
				for(int w0=0; w0+WIN<=nw; w0+=WIN){
					final double[] row=new double[NF+1];
					System.arraycopy(features(f, adds, w0), 0, row, 0, NF);
					row[NF]=y;
					DATA.add(row);
				}
			});
		});
		System.err.println("window samples="+DATA.size());

		final ComplexityFarm.BigNet net=new ComplexityFarm.BigNet(DATA, hidden, 40);
		System.err.println("Fll53b window net trained");
		if(netSave!=null){net.save(netSave); System.err.println("saved "+netSave);}

		System.out.println("#family\torgs\tk\ttrueComplexity\tmodular_cErr%");
		for(final int orgs : new int[]{192, 384, 768, 1536}){
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
				final double[] acc=new double[2];
				final Object lock=new Object();
				FllSkewTest.runPool(threads, evalInst, inst->{
					final long[] fin=new long[2];
					final Fll53b[] sk=new Fll53b[1];
					runStream(sp, 909000L+inst*104729L, B,
						(f, distinct, adds)->{sk[0]=f; fin[0]=distinct; fin[1]=adds;});
					final double c=(double)fin[0]/fin[1];
					double ySum=0;
					int k=0;
					for(int w0=0; w0+WIN<=sk[0].getNumOrgs(); w0+=WIN){
						ySum+=net.predict(features(sk[0], fin[1], w0));
						k++;
					}
					final double cHat=Math.pow(2.0, ySum/k);
					synchronized(lock){
						acc[0]+=Math.abs(cHat-c)/c;
						acc[1]+=c;
					}
				});
				System.out.println(String.format("%s\t%d\t%d\t%.4f\t%.2f",
					names[e], orgs, orgs/WIN, acc[1]/evalInst,
					100*acc[0]/evalInst));
			}
		}
	}
}
