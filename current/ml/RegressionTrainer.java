package ml;

import java.util.ArrayList;
import java.util.Random;

/**
 * RegressionTrainer: trains a feed-forward network for CONTINUOUS outputs and
 * saves it in BBNet format (via CellNet's own serializer, so the file is
 * byte-compatible with everything that consumes BBNets).
 * <p>
 * Differences from ml.Trainer: pure MSE regression (no balancing, cutoffs,
 * FPR/FNR machinery), Adam optimizer with cosine learning-rate decay,
 * mini-batches, input standardization learned from the data and FOLDED INTO
 * the first layer's weights at export — the saved net takes raw inputs.
 * Hidden layers are tanh, output is sigmoid (targets must be in [0,1];
 * normalize upstream like '#dims N 1' data always is).
 * After writing, the net is reloaded through CellNetParser and checked
 * against the in-memory model on training samples (round-trip verification).
 *
 * Usage: java ml.RegressionTrainer in=<data.tsv> out=<net.bbnet> dims=16,32,1
 *          [epochs=60] [batch=8192] [lr=0.003] [wd=1e-4] [seed=1] [vfraction=0.1]
 *
 * Input format: same as train.sh — '#dims <in> <out>' header, then
 * tab-delimited floats, inputs first, target(s) last.  One output supported.
 *
 * @author Amber (for Brian)
 * @date July 2026
 */
public class RegressionTrainer {

	public static void main(String[] args) throws Exception{
		String in=null, out=null;
		int[] dims=null;
		int epochs=60, batch=8192;
		double lr=0.003, wd=1e-4;
		long seed=1;
		double vfraction=0.1;
		boolean rslogFinal=true;
		for(String arg : args){
			final String[] kv=arg.split("=", 2);
			final String k=kv[0].toLowerCase();
			final String v=(kv.length>1) ? kv[1] : "";
			if(k.equals("in")){in=v;}
			else if(k.equals("out")){out=v;}
			else if(k.equals("dims")){
				final String[] p=v.split("[,x]");
				dims=new int[p.length];
				for(int i=0; i<p.length; i++){dims[i]=Integer.parseInt(p[i]);}
			}
			else if(k.equals("epochs")){epochs=Integer.parseInt(v);}
			else if(k.equals("batch")){batch=Integer.parseInt(v);}
			else if(k.equals("lr") || k.equals("alpha")){lr=Double.parseDouble(v);}
			else if(k.equals("wd")){wd=Double.parseDouble(v);}
			else if(k.equals("seed")){seed=Long.parseLong(v);}
			else if(k.equals("vfraction")){vfraction=Double.parseDouble(v);}
			else if(k.equals("final")){rslogFinal=v.equalsIgnoreCase("rslog");}
			else{throw new IllegalArgumentException("unknown arg "+arg);}
		}
		if(in==null || out==null || dims==null){
			throw new IllegalArgumentException("required: in=, out=, dims=");
		}
		if(dims[dims.length-1]!=1){
			throw new IllegalArgumentException("one output supported; last dim must be 1");
		}

		// ---- Load ----
		final ArrayList<float[]> X=new ArrayList<>();
		final ArrayList<Float> Y=new ArrayList<>();
		final int NI=dims[0];
		for(String line : java.nio.file.Files.readAllLines(java.nio.file.Paths.get(in))){
			if(line.isEmpty() || line.charAt(0)=='#'){continue;}
			final String[] f=line.split("\t");
			if(f.length<NI+1){continue;}
			final float[] x=new float[NI];
			for(int i=0; i<NI; i++){x[i]=Float.parseFloat(f[i]);}
			X.add(x);
			Y.add(Float.parseFloat(f[NI]));
		}
		final int n=X.size();
		System.err.println("loaded "+n+" samples, dims="+java.util.Arrays.toString(dims));

		// ---- Standardize inputs (folded into layer 1 at export) ----
		final double[] mean=new double[NI], sd=new double[NI];
		for(float[] x : X){for(int i=0; i<NI; i++){mean[i]+=x[i];}}
		for(int i=0; i<NI; i++){mean[i]/=n;}
		for(float[] x : X){for(int i=0; i<NI; i++){sd[i]+=(x[i]-mean[i])*(x[i]-mean[i]);}}
		for(int i=0; i<NI; i++){sd[i]=Math.max(1e-9, Math.sqrt(sd[i]/n));}

		// ---- Split train/validation ----
		final Random rnd=new Random(seed);
		final int[] perm=new int[n];
		for(int i=0; i<n; i++){perm[i]=i;}
		for(int i=n-1; i>0; i--){
			final int j=rnd.nextInt(i+1);
			final int t=perm[i]; perm[i]=perm[j]; perm[j]=t;
		}
		final int nVal=(int)(n*vfraction);
		final int nTrain=n-nVal;

		// ---- Model arrays: W[l][i][j], B[l][i]; layer l maps dims[l]->dims[l+1] ----
		final int L=dims.length-1;
		final double[][][] W=new double[L][][];
		final double[][] Bs=new double[L][];
		for(int l=0; l<L; l++){
			W[l]=new double[dims[l+1]][dims[l]];
			Bs[l]=new double[dims[l+1]];
			final double scale=1.0/Math.sqrt(dims[l]);
			for(double[] row : W[l]){
				for(int j=0; j<row.length; j++){row[j]=rnd.nextGaussian()*scale;}
			}
		}
		// Adam state
		final double[][][] mW=new double[L][][], vW=new double[L][][];
		final double[][] mB=new double[L][], vB=new double[L][];
		for(int l=0; l<L; l++){
			mW[l]=new double[dims[l+1]][dims[l]]; vW[l]=new double[dims[l+1]][dims[l]];
			mB[l]=new double[dims[l+1]]; vB[l]=new double[dims[l+1]];
		}

		// ---- Train ----
		final double be1=0.9, be2=0.999, eps=1e-8;
		long step=0;
		double bestVal=Double.MAX_VALUE;
		double[][][] bestW=null;
		double[][] bestB=null;
		final double[][] act=new double[dims.length][];
		final double[][] delta=new double[dims.length][];
		for(int i=0; i<dims.length; i++){act[i]=new double[dims[i]]; delta[i]=new double[dims[i]];}
		final int[] order=new int[nTrain];
		for(int i=0; i<nTrain; i++){order[i]=perm[i];}

		for(int ep=1; ep<=epochs; ep++){
			// shuffle
			for(int i=nTrain-1; i>0; i--){
				final int j=rnd.nextInt(i+1);
				final int t=order[i]; order[i]=order[j]; order[j]=t;
			}
			final double lrNow=lr*0.5*(1+Math.cos(Math.PI*(ep-1)/epochs));
			double trainMse=0;
			for(int start=0; start<nTrain; start+=batch){
				final int end=Math.min(nTrain, start+batch);
				final int bs=end-start;
				// zero grads (reuse Adam arrays for grads via local accumulation)
				final double[][][] gW=new double[L][][];
				final double[][] gB=new double[L][];
				for(int l=0; l<L; l++){
					gW[l]=new double[dims[l+1]][dims[l]];
					gB[l]=new double[dims[l+1]];
				}
				for(int s=start; s<end; s++){
					final float[] x=X.get(order[s]);
					final double y=Y.get(order[s]);
					for(int i=0; i<NI; i++){act[0][i]=(x[i]-mean[i])/sd[i];}
					for(int l=0; l<L; l++){
						for(int i=0; i<dims[l+1]; i++){
							double z=Bs[l][i];
							final double[] wr=W[l][i];
							for(int j=0; j<dims[l]; j++){z+=wr[j]*act[l][j];}
							act[l+1][i]=(l==L-1) ? fin(z, rslogFinal) : Math.tanh(z);
						}
					}
					final double outv=act[L][0];
					final double err=outv-y;
					trainMse+=err*err;
					// MSE x final-activation derivative (rslog: d/dz = exp(-|fx|))
					delta[L][0]=2*err*(rslogFinal ? Math.exp(-Math.abs(outv)) : outv*(1-outv));
					for(int l=L-1; l>=0; l--){
						for(int i=0; i<dims[l+1]; i++){
							final double d=delta[l+1][i];
							gB[l][i]+=d;
							final double[] wr=W[l][i], gr=gW[l][i];
							for(int j=0; j<dims[l]; j++){gr[j]+=d*act[l][j];}
						}
						if(l>0){
							for(int j=0; j<dims[l]; j++){
								double s2=0;
								for(int i=0; i<dims[l+1]; i++){s2+=delta[l+1][i]*W[l][i][j];}
								delta[l][j]=s2*(1-act[l][j]*act[l][j]);   // tanh'
							}
						}
					}
				}
				// Adam update
				step++;
				final double c1=1-Math.pow(be1, step), c2=1-Math.pow(be2, step);
				for(int l=0; l<L; l++){
					for(int i=0; i<dims[l+1]; i++){
						for(int j=0; j<dims[l]; j++){
							final double g=gW[l][i][j]/bs+wd*W[l][i][j];
							mW[l][i][j]=be1*mW[l][i][j]+(1-be1)*g;
							vW[l][i][j]=be2*vW[l][i][j]+(1-be2)*g*g;
							W[l][i][j]-=lrNow*(mW[l][i][j]/c1)/(Math.sqrt(vW[l][i][j]/c2)+eps);
						}
						final double g=gB[l][i]/bs;
						mB[l][i]=be1*mB[l][i]+(1-be1)*g;
						vB[l][i]=be2*vB[l][i]+(1-be2)*g*g;
						Bs[l][i]-=lrNow*(mB[l][i]/c1)/(Math.sqrt(vB[l][i]/c2)+eps);
					}
				}
			}
			// validation
			double valMse=0;
			for(int s=nTrain; s<n; s++){
				final double p=predict(X.get(perm[s]), mean, sd, W, Bs, dims, rslogFinal);
				final double e=p-Y.get(perm[s]);
				valMse+=e*e;
			}
			valMse/=Math.max(1, nVal);
			if(valMse<bestVal){
				bestVal=valMse;
				bestW=deepCopy(W);
				bestB=new double[L][];
				for(int l=0; l<L; l++){bestB[l]=Bs[l].clone();}
			}
			if(ep%5==0 || ep==1 || ep==epochs){
				System.err.println(String.format(
					"epoch %d lr=%.5f trainMSE=%.6f valMSE=%.6f best=%.6f",
					ep, lrNow, trainMse/nTrain, valMse, bestVal));
			}
		}
		final double[][][] Wf=(bestW!=null) ? bestW : W;
		final double[][] Bf=(bestB!=null) ? bestB : Bs;

		// ---- Export: build a real CellNet, fold standardization into layer 1 ----
		final CellNet net=new CellNet(dims, Math.max(0, seed), 1.0f, 0f, 1, new ArrayList<>());
		net.randomize();   // allocates dense edge topology + weight arrays
		for(int l=0; l<L; l++){
			final Cell[] layer=net.net[l+1];
			for(int i=0; i<dims[l+1]; i++){
				final Cell c=layer[i];
				c.function=Function.getFunction((l==L-1)
					? (rslogFinal ? Function.RSLOG : Function.SIG) : Function.TANH);
				double b=Bf[l][i];
				for(int j=0; j<dims[l]; j++){
					double w=Wf[l][i][j];
					if(l==0){       // fold (x-mean)/sd into weights and bias
						w/=sd[j];
						b-=Wf[l][i][j]*mean[j]/sd[j];
					}
					c.weights[j]=(float)w;
				}
				c.setBias((float)b, true);
			}
		}
		final structures.ByteBuilder bb=net.toBytes();
		java.nio.file.Files.write(java.nio.file.Paths.get(out), bb.toBytes());
		System.err.println("wrote "+out+" (valMSE="+String.format("%.6f", bestVal)+")");

		// ---- Round-trip verification ----
		final CellNet check=CellNetParser.load(out);
		double maxDiff=0;
		for(int s=0; s<Math.min(1000, n); s++){
			final float[] x=X.get(s);
			check.applyInput(x);
			final double a=check.feedForward();
			final double b=predict(x, mean, sd, Wf, Bf, dims, rslogFinal);
			maxDiff=Math.max(maxDiff, Math.abs(a-b));
		}
		System.err.println(String.format(
			"round-trip check vs in-memory model: maxDiff=%.2e over %d samples %s",
			maxDiff, Math.min(1000, n), (maxDiff<1e-3 ? "(OK)" : "(SUSPICIOUS!)")));
	}

	static double fin(double z, boolean rslog){
		if(!rslog){return 1.0/(1.0+Math.exp(-z));}
		return (z<0) ? -Math.log(1-z) : Math.log(1+z);
	}

	static double predict(float[] x, double[] mean, double[] sd,
			double[][][] W, double[][] Bs, int[] dims, boolean rslogFinal){
		final int L=dims.length-1;
		double[] a=new double[dims[0]];
		for(int i=0; i<dims[0]; i++){a[i]=(x[i]-mean[i])/sd[i];}
		for(int l=0; l<L; l++){
			final double[] nx=new double[dims[l+1]];
			for(int i=0; i<dims[l+1]; i++){
				double z=Bs[l][i];
				for(int j=0; j<dims[l]; j++){z+=W[l][i][j]*a[j];}
				nx[i]=(l==L-1) ? fin(z, rslogFinal) : Math.tanh(z);
			}
			a=nx;
		}
		return a[0];
	}

	static double[][][] deepCopy(double[][][] w){
		final double[][][] o=new double[w.length][][];
		for(int l=0; l<w.length; l++){
			o[l]=new double[w[l].length][];
			for(int i=0; i<w[l].length; i++){o[l][i]=w[l][i].clone();}
		}
		return o;
	}
}
