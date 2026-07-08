package cardinality;

/**
 * FllDupCal: duplicate-aware estimation for FutureLogLog2, trained and tested.
 * <p>
 * Model: log2(n) = meanExp + f(lsbFrac, msbFrac), f = quadratic (6 coefficients)
 * fit by least squares on structured-duplication training streams.  The servo
 * pins (lsbFrac, msbFrac) on unique streams; duplication moves the pair off the
 * unique manifold in a direction n cannot mimic, so f absorbs the duplicate
 * distortion without knowing the rate.  NAIVE variant: f fit on dup=1 only
 * (a constant-offset-style corrector, the state-blind baseline).
 * Evaluation: held-out unique streams (HC) and Brian's low-complexity protocol
 * (pool draws via min(rand,rand), iterations=4), true distinct tracked exactly.
 *
 * @author Amber
 * @date July 2026
 */
public class FllDupCal {

	static double[] features(double lsb, double msb){
		return new double[]{1, lsb, msb, lsb*lsb, msb*msb, lsb*msb};
	}

	/** Returns {meanExp, lsbFrac, msbFrac} for a sketch. */
	static double[] summarize(FutureLogLog2 f){
		final int nw=f.getNumWords();
		long lsb=0, msb=0, expSum=0;
		for(int w=0; w<nw; w++){
			final int word=f.getWord(w);
			expSum+=(word>>>12)&0xF;
			lsb+=Integer.bitCount(word&0x555);
			msb+=Integer.bitCount(word&0xAAA);
		}
		return new double[]{f.getGlobalExp()+(double)expSum/nw,
			(double)lsb/(nw*6), (double)msb/(nw*6)};
	}

	public static void main(String[] args) throws Exception{
		final int B=(args.length>0) ? Integer.parseInt(args[0]) : 2048;
		final int trainInst=(args.length>1) ? Integer.parseInt(args[1]) : 12;
		CardinalityTracker.clampToAdded=false;
		final double LOG2=Math.log(2.0);

		// ---- Training: structured sweeps at various duplication levels ----
		final long[] ns={50000, 100000, 200000, 400000, 800000};
		final double[] dups={1, 1.5, 2, 2.5, 3, 4, 8};
		final java.util.ArrayList<double[]> X=new java.util.ArrayList<>();
		final java.util.ArrayList<Double> Y=new java.util.ArrayList<>();
		final java.util.ArrayList<double[]> X1=new java.util.ArrayList<>(); // dup=1 only
		final java.util.ArrayList<Double> Y1=new java.util.ArrayList<>();

		for(long n : ns){
			for(double dup : dups){
				for(int inst=0; inst<trainInst; inst++){
					final FutureLogLog2 f=new FutureLogLog2(B, 31, -1, 0);
					final long seedBase=new java.util.Random(4242L+inst*7919L).nextLong();
					final long total=(long)(n*dup);
					long added=0;
					outer:
					while(true){
						for(long i=0; i<n; i++){
							f.add(seedBase+i*0x9E3779B97F4A7C15L);
							if(++added>=total){break outer;}
						}
					}
					final double[] s=summarize(f);
					final double y=Math.log(n)/LOG2-s[0];
					final double[] x=features(s[1], s[2]);
					X.add(x); Y.add(y);
					if(dup==1){X1.add(x); Y1.add(y);}
				}
			}
		}

		final double[] beta=lsq(X, Y);
		final double[] beta1=lsq(X1, Y1);
		System.out.println("# dup-aware beta: "+java.util.Arrays.toString(beta));
		System.out.println("# naive beta:     "+java.util.Arrays.toString(beta1));

		// ---- Evaluation ----
		System.out.println("#stream\ttrueN\tnaive_err%\tdupaware_err%");
		// (a) held-out unique streams
		for(long n : new long[]{70000, 300000, 600000}){
			evalPoint("HC_unique", null, n, 0, B, beta, beta1);
		}
		// (b) Brian's LC protocol: pool P, draws=4P, pos=min(rand,rand)
		for(int P : new int[]{200000, 500000, 1000000}){
			evalPoint("LC_minrand", null, 0, P, B, beta, beta1);
		}
	}

	static void evalPoint(String name, Object unused, long nUnique, int pool,
			int B, double[] beta, double[] beta1){
		final double LOG2=Math.log(2.0);
		final int INST=16;
		double e0=0, e1=0; long tn=0;
		for(int inst=0; inst<INST; inst++){
			final FutureLogLog2 f=new FutureLogLog2(B, 31, -1, 0);
			final java.util.Random seedr=new java.util.Random(999000L+inst*104729L);
			long trueN;
			if(pool==0){
				final long seedBase=seedr.nextLong();
				for(long i=0; i<nUnique; i++){f.add(seedBase+i*0x9E3779B97F4A7C15L);}
				trueN=nUnique;
			}else{
				final long seedBase=seedr.nextLong();
				final java.util.SplittableRandom rng=new java.util.SplittableRandom(seedr.nextLong());
				final boolean[] seen=new boolean[pool];
				long distinct=0;
				final long draws=4L*pool;
				for(long d=0; d<draws; d++){
					final int pos=Math.min(rng.nextInt(pool), rng.nextInt(pool));
					if(!seen[pos]){seen[pos]=true; distinct++;}
					f.add(seedBase+pos*0x9E3779B97F4A7C15L);
				}
				trueN=distinct;
			}
			final double[] s=summarize(f);
			final double[] x=features(s[1], s[2]);
			final double estA=Math.pow(2.0, s[0]+dot(x, beta));
			final double estN=Math.pow(2.0, s[0]+dot(x, beta1));
			e1+=Math.abs(estA-trueN)/trueN;
			e0+=Math.abs(estN-trueN)/trueN;
			tn=trueN;
		}
		System.out.println(String.format("%s\t%d\t%.2f\t%.2f",
			name, tn, 100.0*e0/INST, 100.0*e1/INST));
	}

	static double dot(double[] a, double[] b){
		double s=0;
		for(int i=0; i<a.length; i++){s+=a[i]*b[i];}
		return s;
	}

	/** Least squares via normal equations + Gaussian elimination. */
	static double[] lsq(java.util.List<double[]> X, java.util.List<Double> Y){
		final int m=X.get(0).length;
		final double[][] A=new double[m][m+1];
		for(int r=0; r<X.size(); r++){
			final double[] x=X.get(r);
			final double y=Y.get(r);
			for(int i=0; i<m; i++){
				for(int j=0; j<m; j++){A[i][j]+=x[i]*x[j];}
				A[i][m]+=x[i]*y;
			}
		}
		for(int c=0; c<m; c++){
			int piv=c;
			for(int r=c+1; r<m; r++){if(Math.abs(A[r][c])>Math.abs(A[piv][c])){piv=r;}}
			final double[] t=A[c]; A[c]=A[piv]; A[piv]=t;
			for(int r=0; r<m; r++){
				if(r==c || A[c][c]==0){continue;}
				final double fct=A[r][c]/A[c][c];
				for(int j=c; j<=m; j++){A[r][j]-=fct*A[c][j];}
			}
		}
		final double[] beta=new double[m];
		for(int i=0; i<m; i++){beta[i]=(A[i][i]!=0) ? A[i][m]/A[i][i] : 0;}
		return beta;
	}
}
