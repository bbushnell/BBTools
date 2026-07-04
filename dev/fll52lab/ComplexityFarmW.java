package cardinality;

/**
 * ComplexityFarmW: Brian's MODULAR window design — one complexity network that
 * predicts on a fixed 256-word (512-byte) window of FLL2 state, applied k times
 * for a k*512-byte estimator and averaged in log space.
 * <p>
 * Training streams use RANDOM sketch sizes (256..2048 words) and randomized
 * continuous shape priors (the grid-vs-continuous lesson); every 256-word
 * window of every snapshot emits one training row.  The window's share of
 * adds is estimated as adds*256/numWords (unbiased under uniform hashing).
 * Eval: held-out families at 512B/1KB/2KB/4KB estimator sizes, same single
 * net, k=1/2/4/8 window predictions averaged — the modularity demonstration.
 *
 * Usage: java cardinality.ComplexityFarmW [streams] [evalInst] [threads] [hidden] [netSave]
 *
 * @author Amber (design: Brian)
 * @date July 2026
 */
public class ComplexityFarmW {

	static final long GOLD=0x9E3779B97F4A7C15L;
	static final double LOG2=Math.log(2.0);
	static final int WIN=256;              // words per window (512 bytes)
	static final int NF=16;

	/** Window features over words [w0, w0+WIN): same 16 features as
	 * ComplexityFarm but window-local, with adds scaled to the window share. */
	static double[] features(FutureLogLog2 f, long adds, int w0){
		final int nw=f.getNumWords();
		final double winAdds=Math.max(1.0, adds*(double)WIN/nw);
		long lsb=0, msb=0;
		double expSum=0, expSat=0;
		int satCnt=0;
		for(int w=w0; w<w0+WIN; w++){
			final int word=f.getWord(w);
			expSum+=(word>>>12)&0xF;
			lsb+=Integer.bitCount(word&0x555);
			msb+=Integer.bitCount(word&0xAAA);
			if(f.ctrX!=null && f.ctrX[1][w]==3){expSat+=(word>>>12)&0xF; satCnt++;}
		}
		final double meanExp=f.getGlobalExp()+expSum/WIN;
		final double[] out=new double[NF];
		out[0]=Math.log(winAdds)/LOG2-meanExp;
		out[1]=(double)lsb/(WIN*6);
		out[2]=(double)msb/(WIN*6);
		out[3]=(satCnt>0) ? expSat/satCnt-expSum/WIN : 0;
		for(int sem=0; sem<4; sem++){
			final int[] cc=new int[4];
			for(int w=w0; w<w0+WIN; w++){
				cc[(f.ctrX==null) ? 0 : f.ctrX[sem][w]]++;
			}
			for(int k=0; k<3; k++){out[4+sem*3+k]=(double)cc[k+1]/WIN;}
		}
		return out;
	}

	public static void main(String[] args) throws Exception{
		final int streams=(args.length>0) ? Integer.parseInt(args[0]) : 400000;
		final int evalInst=(args.length>1) ? Integer.parseInt(args[1]) : 64;
		final int threads=(args.length>2) ? Integer.parseInt(args[2])
			: Runtime.getRuntime().availableProcessors();
		final int hidden=(args.length>3) ? Integer.parseInt(args[3]) : 32;
		final String netSave=(args.length>4) ? args[4] : null;
		CardinalityTracker.clampToAdded=false;

		// ---- Farm: random sketch size x random continuous shape ----
		final java.util.List<double[]> DATA=java.util.Collections.synchronizedList(
			new java.util.ArrayList<>());
		final int[] SIZES={256, 512, 1024, 2048};   // words
		FllSkewTest.runPool(threads, streams, j->{
			final java.util.SplittableRandom prng=new java.util.SplittableRandom(
				new java.util.Random(4141L+j).nextLong());
			final int words=SIZES[prng.nextInt(SIZES.length)];
			final int B=words*6;
			final ComplexityFarm.Spec sp=randomSpec(prng, B);
			ComplexityFarm.runStream(sp, 505000L+j*104729L, B, (f, distinct, adds)->{
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

		// ---- Train the single window net ----
		final ComplexityFarm.BigNet net=new ComplexityFarm.BigNet(DATA, hidden, 40);
		System.err.println("window net trained: "+NF+"-"+hidden+"-1");
		if(netSave!=null){net.save(netSave); System.err.println("saved "+netSave);}

		// ---- Modularity eval: same net, k windows averaged, 4 sizes ----
		System.out.println("#family\twords\tk\ttrueComplexity\tmodular_cErr%");
		for(final int words : new int[]{256, 512, 1024, 2048}){
			final int B=words*6;
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
					final double[][] lastF=new double[1][];
					final long[] fin=new long[2];
					final FutureLogLog2[] sk=new FutureLogLog2[1];
					ComplexityFarm.runStream(sp, 909000L+inst*104729L, B,
						(f, distinct, adds)->{sk[0]=f; fin[0]=distinct; fin[1]=adds;});
					final double c=(double)fin[0]/fin[1];
					// k window predictions, averaged in log space
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
					}
				});
				System.out.println(String.format("%s\t%d\t%d\t%.4f\t%.2f",
					names[e], words, words/WIN, acc[1]/evalInst,
					100*acc[0]/evalInst));
			}
		}
	}

	static ComplexityFarm.Spec randomSpec(java.util.SplittableRandom prng, int B){
		switch(prng.nextInt(4)){
			case 0: {
				final long n=(long)(24L*B*Math.exp(prng.nextDouble()*Math.log(50)));
				double dup=Math.exp(prng.nextDouble()*Math.log(32));
				if(n*dup>40_000_000L){dup=40_000_000.0/n;}
				return new ComplexityFarm.Spec(0, n, dup, 0);}
			case 1: {
				final double zs=0.3+prng.nextDouble()*1.9;
				final long pool=(long)(40L*B*Math.exp(prng.nextDouble()*Math.log(25)));
				double mult=Math.exp(Math.log(2)+prng.nextDouble()*Math.log(10));
				if(pool*mult>40_000_000L){mult=40_000_000.0/pool;}
				return new ComplexityFarm.Spec(1, pool, zs, mult);}
			case 2: {
				final long k=Math.max(1, (long)Math.exp(prng.nextDouble()*Math.log(1000)));
				final double hf=0.3+prng.nextDouble()*0.695;
				final long fresh=(long)(30L*B*Math.exp(prng.nextDouble()*Math.log(20)));
				return new ComplexityFarm.Spec(2, k, hf, fresh);}
			default: {
				final long pool=(long)(30L*B*Math.exp(prng.nextDouble()*Math.log(83)));
				double iter=1.2+prng.nextDouble()*8.8;
				if(pool*iter>40_000_000L){iter=40_000_000.0/pool;}
				return new ComplexityFarm.Spec(3, pool, iter, 0);}
		}
	}
}
