package cardinality;

/**
 * Fll52FarmW: modular window farm for the Fll52 species (honest bytes:
 * 256 words = 512B = 1280 tails, counter in-word).  One net predicts
 * log2(complexity) on a 256-word window; k*512B estimators average k
 * predictions in log space.  Training: random sketch sizes x randomized
 * continuous shape priors; every window of every snapshot = one row.
 * Eval: held-out families at 512B/1KB/2KB/4KB, one net.
 *
 * Usage: java cardinality.Fll52FarmW [streams] [evalInst] [threads] [hidden] [netSave]
 *
 * @author Amber (design: Brian)
 * @date July 2026
 */
public class Fll52FarmW {

	static final long GOLD=0x9E3779B97F4A7C15L;
	static final double LOG2=Math.log(2.0);
	static final int WIN=256;              // words per window: exactly 512 bytes
	static final int NF=16;                // matches ComplexityFarm.BigNet

	interface Sink {void emit(Fll52 f, long distinct, long adds);}

	/** Same stream generators as ComplexityFarm.runStream, feeding Fll52. */
	static void runStream(ComplexityFarm.Spec sp, long seed, int B, Sink sink){
		final Fll52 f=new Fll52(B, 31, -1, 0);
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

	/** 16 window features from Fll52 state. */
	static double[] features(Fll52 f, long adds, int w0){
		final int nw=f.getNumWords();
		final double winAdds=Math.max(1.0, adds*(double)WIN/nw);
		long lsb=0, msb=0;
		double expSum=0, expSat=0;
		int satCnt=0, zeroCnt=0;
		final int[] cc=new int[4];
		for(int w=w0; w<w0+WIN; w++){
			final int word=f.getWord(w);
			final int le=(word>>>12)&0xF;
			expSum+=le;
			lsb+=Integer.bitCount(word&Fll52.LSB_MASK);
			msb+=Integer.bitCount(word&Fll52.MSB_MASK);
			final int c=word&Fll52.CTR_MASK;
			cc[c]++;
			if(c==3){expSat+=le; satCnt++;}
			if((word&0xFFC)==0){zeroCnt++;}
		}
		final double meanExp=f.getGlobalExp()+expSum/WIN;
		final double lf=(double)lsb/(WIN*5), mf=(double)msb/(WIN*5);
		final double f1=(double)cc[1]/WIN, f2=(double)cc[2]/WIN, f3=(double)cc[3]/WIN;
		final double r=Math.log(winAdds)/LOG2-meanExp;
		final double dExpSat=(satCnt>0) ? expSat/satCnt-expSum/WIN : 0;
		return new double[]{r, lf, mf, dExpSat, f1, f2, f3,
			(f1+2*f2+3*f3)/3, f2+f3, lf*mf, r*f3, mf*f3, lf*lf, mf*mf, f3*f3,
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

		final java.util.List<double[]> DATA=java.util.Collections.synchronizedList(
			new java.util.ArrayList<>());
		final int[] SIZES={256, 512, 1024, 2048};
		FllSkewTest.runPool(threads, streams, j->{
			final java.util.SplittableRandom prng=new java.util.SplittableRandom(
				new java.util.Random(5252L+j).nextLong());
			final int words=SIZES[prng.nextInt(SIZES.length)];
			final int B=words*5;
			final ComplexityFarm.Spec sp=ComplexityFarmW.randomSpec(prng, B);
			runStream(sp, 606000L+j*104729L, B, (f, distinct, adds)->{
				if(distinct<24L*B){return;}
				final double y=Math.log((double)distinct/adds)/LOG2;
				final int nw=f.getNumWords();
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
		System.err.println("Fll52 window net trained");
		if(netSave!=null){net.save(netSave); System.err.println("saved "+netSave);}

		System.out.println("#family\twords\tk\ttrueComplexity\tmodular_cErr%\tmodular_nErr%");
		for(final int words : new int[]{256, 512, 1024, 2048}){
			final int B=words*5;
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
					final Fll52[] sk=new Fll52[1];
					runStream(sp, 909000L+inst*104729L, B,
						(f, distinct, adds)->{sk[0]=f; fin[0]=distinct; fin[1]=adds;});
					final double c=(double)fin[0]/fin[1];
					double ySum=0;
					int k=0;
					for(int w0=0; w0+WIN<=sk[0].getNumWords(); w0+=WIN){
						ySum+=net.predict(features(sk[0], fin[1], w0));
						k++;
					}
					final double cHat=Math.pow(2.0, ySum/k);
					synchronized(lock){
						acc[0]+=Math.abs(cHat-c)/c;
						acc[1]+=c;
						acc[2]+=Math.abs(cHat*fin[1]-fin[0])/fin[0];
					}
				});
				System.out.println(String.format("%s\t%d\t%d\t%.4f\t%.2f\t%.2f",
					names[e], words, words/WIN, acc[1]/evalInst,
					100*acc[0]/evalInst, 100*acc[2]/evalInst));
			}
		}
	}
}
