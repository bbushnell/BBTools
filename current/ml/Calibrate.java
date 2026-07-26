package ml;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.PreParser;
import shared.Timer;
import shared.Tools;
import structures.FloatList;

/**
 * Fits probability-calibration constants for the monotone map
 *   p = K * sigmoid(a*logit(x) + b)^c
 * so that calibrate(rawScore) approximates P(label==1 | rawScore). This is the
 * conversion from a network's raw continuous output (or any raw score) into a
 * calibrated probability. It is fit per dataset, independently, so callers with
 * many levels/nets simply run it once each (see calibrate.sh).
 *
 * Generic input (either form):
 *   net= in=vectors.tsv : load the net, run it on each vector (a "#dims IN OUT"
 *                         header set), pairing the net output with the vector's
 *                         first target column -> (rawScore, label).
 *   in=pairs.tsv        : read (score,label) straight from two columns (default
 *                         col 0 and col 1); no net. Works for any floats and
 *                         booleans (label parsed as >=0.5 -> 1).
 *
 * Scales to any row count: raw scores are streamed into a fixed histogram
 * (fitbins buckets), so the fit cost is O(iterations*bins), independent of N.
 *
 * Output: a "#cal K a b c" line (paste-ready into a .bbnets file for
 * clade.SerialNNLoader / CladeConfidence.calibrate), plus log-loss and ECE
 * (expected calibration error) before and after, so the improvement is visible.
 *
 * Also emits a RELIABILITY TABLE (per probability bin: n, mean predicted,
 * empirical rate, gap). A scalar ECE can look excellent while the tail is badly
 * miscalibrated; the table is what makes that visible.
 *
 * And it fits a nonparametric alternative -- an isotonic (pool-adjacent-
 * violators) monotone lookup table -- then compares it against the 4-parameter
 * form on HELD-OUT rows, so the comparison is not biased toward the model with
 * more free parameters. The 4-parameter form is what CladeConfidence already
 * implements, so it stays the default; the LUT is there to show what is being
 * given up, and can be dumped with lutout= if a level needs it.
 *
 * Crashes loudly on malformed input rather than emitting a silent wrong fit.
 *
 * @author Noire
 */
public class Calibrate {

	public static void main(String[] args){
		Timer t=new Timer();
		Calibrate x=new Calibrate(args);
		x.process(t);
	}

	public Calibrate(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=", 2);
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("in") || a.equals("in1")){in=b;}
			else if(a.equals("net") || a.equals("netin")){netFile=b;}
			else if(a.equals("out")){out=b;}
			else if(a.equals("reliability") || a.equals("rel")){relOut=b;}
			else if(a.equals("lutout")){lutOut=b;}
			else if(a.equals("scorecol")){scoreCol=Integer.parseInt(b);}
			else if(a.equals("labelcol")){labelCol=Integer.parseInt(b);}
			else if(a.equals("bins")){bins=Integer.parseInt(b);}
			else if(a.equals("fitbins")){fitBins=Integer.parseInt(b);}
			else if(a.equals("holdout")){holdout=Integer.parseInt(b);}
			else if(a.equals("iters") || a.equals("maxiters")){maxIters=Integer.parseInt(b);}
			else if(a.equals("restarts")){restarts=Integer.parseInt(b);}
			else if(a.equals("label") || a.equals("name")){label=b;}
			else if(b==null && !arg.startsWith("-") && in==null){in=arg;}
			else{throw new RuntimeException("Unknown parameter '"+arg+"'");}
		}
		if(in==null){throw new RuntimeException("Required: in=<vectors-or-score/label-pairs.tsv>");}
		if(bins<1 || fitBins<16){throw new RuntimeException("bins>=1 and fitbins>=16 required");}
		if(holdout<0 || holdout==1){throw new RuntimeException("holdout must be 0 (disabled) or >=2 (use every Nth row for evaluation)");}
	}

	public void process(Timer t){
		final Accum[] acc=(netFile!=null) ? loadFromNet() : loadFromPairs();
		final Accum accFit=acc[0], accEval=acc[1];
		final Hist hAll=Accum.merge(accFit, accEval).finish();
		if(hAll.totalN<4){throw new RuntimeException("Only "+hAll.totalN+" usable rows in "+in+"; need at least 4 to fit.");}

		// The SHIPPED constants are fit on every row; the holdout split below is only for the comparison.
		final float[] cal=fit(hAll);

		final double[] pRaw=predictRaw(hAll), pCal=predictCal(hAll, cal);
		if(label!=null){outstream.println("Dataset:   "+label);}
		outstream.println("Rows="+hAll.totalN+"  mean_label="+String.format("%.5f", hAll.totalPos/(double)hAll.totalN)
				+"  occupied_bins="+hAll.meanRaw.length+"/"+fitBins);
		outstream.println(String.format("log-loss:  raw=%.6f  ->  calibrated=%.6f", logLoss(hAll, pRaw), logLoss(hAll, pCal)));
		outstream.println(String.format("ECE(%d):     raw=%.6f  ->  calibrated=%.6f", bins, ece(hAll, pRaw), ece(hAll, pCal)));
		outstream.println(String.format("params:    K=%.6f  a=%.6f  b=%.6f  c=%.6f", cal[0], cal[1], cal[2], cal[3]));
		final String calLine="#cal "+cal[0]+" "+cal[1]+" "+cal[2]+" "+cal[3];
		outstream.println(calLine);

		compare(accFit, accEval);

		final String table=reliabilityTable(hAll, pCal, pRaw);
		outstream.println();
		outstream.println(table);

		if(out!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(out, true, false, false);
			bsw.start();
			bsw.print((calLine+"\n").getBytes());
			bsw.poisonAndWait();
			outstream.println("Wrote "+out);
		}
		if(relOut!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(relOut, true, false, false);
			bsw.start();
			if(label!=null){bsw.print(("#dataset\t"+label+"\n").getBytes());}
			bsw.print((calLine+"\n").getBytes());
			bsw.print((table+"\n").getBytes());
			bsw.poisonAndWait();
			outstream.println("Wrote "+relOut);
		}
		if(lutOut!=null){
			final Lut lut=buildLut(hAll);
			ByteStreamWriter bsw=new ByteStreamWriter(lutOut, true, false, false);
			bsw.start();
			bsw.print(("#lut\t"+lut.x.length+"\n#rawScore\tprobability\n").getBytes());
			for(int i=0; i<lut.x.length; i++){
				bsw.print((lut.x[i]+"\t"+lut.y[i]+"\n").getBytes());
			}
			bsw.poisonAndWait();
			outstream.println("Wrote "+lutOut+" ("+lut.x.length+" knots)");
		}
		t.stop();
		outstream.println("Time: "+t);
	}

	/**
	 * Fits BOTH models on the fit split and scores BOTH on the held-out split.
	 * In-sample the isotonic table always wins (it has one free value per occupied
	 * bin against four parameters); only a holdout says whether that is real.
	 */
	private void compare(Accum accFit, Accum accEval){
		if(holdout<2){return;}
		final Hist hFit=accFit.finish(), hEval=accEval.finish();
		if(hFit.totalN<4 || hEval.totalN<4){
			outstream.println("Holdout comparison SKIPPED: fit="+hFit.totalN+" eval="+hEval.totalN+" rows.");
			return;
		}
		final float[] calF=fit(hFit);
		final Lut lut=buildLut(hFit);

		final double[] pRaw=predictRaw(hEval), pPar=predictCal(hEval, calF), pLut=predictLut(hEval, lut);
		outstream.println();
		outstream.println("Held-out comparison (fit on "+hFit.totalN+" rows, scored on "+hEval.totalN+" rows;"
				+" isotonic knots="+lut.x.length+"):");
		outstream.println(String.format("           %-12s %-12s %-12s", "raw", "4-param", "isotonic"));
		outstream.println(String.format("log-loss:  %-12.6f %-12.6f %-12.6f",
				logLoss(hEval, pRaw), logLoss(hEval, pPar), logLoss(hEval, pLut)));
		outstream.println(String.format("ECE(%d):     %-12.6f %-12.6f %-12.6f", bins,
				ece(hEval, pRaw), ece(hEval, pPar), ece(hEval, pLut)));
	}

	/*---------- input loaders (stream directly into the histograms) ----------*/

	/** Load a net, run it on each vector; label = first target column (>=0.5 -> 1). */
	private Accum[] loadFromNet(){
		final CellNet net=CellNetParser.load(netFile);
		if(net==null){throw new RuntimeException("Failed to load net "+netFile);}
		final int nIn=net.numInputs(), nOut=net.numOutputs();
		if(nOut!=1){throw new RuntimeException("Calibrate expects a single-output net; "+netFile+" has "+nOut+".");}
		final Accum fitAcc=new Accum(fitBins), evalAcc=new Accum(fitBins);
		final FloatList fl=new FloatList(nIn);
		final ByteFile bf=ByteFile.makeByteFile(in, true);
		byte[] line;
		long row=0;
		while((line=bf.nextLine())!=null){
			if(line.length==0 || line[0]=='#'){continue;}
			final String[] cols=new String(line).split("\t", -1);
			if(cols.length<nIn+nOut){
				throw new RuntimeException("Vector row has "+cols.length+" cols; net needs "+(nIn+nOut)+" (in+out).");
			}
			fl.size=0;
			for(int j=0; j<nIn; j++){fl.add(Float.parseFloat(cols[j]));}
			net.applyInput(fl);
			final float raw=net.feedForward();
			final boolean lab=Float.parseFloat(cols[nIn])>=0.5f;
			if(holdout>=2 && row%holdout==0){evalAcc.add(raw, lab);}else{fitAcc.add(raw, lab);}
			row++;
		}
		bf.close();
		return new Accum[]{fitAcc, evalAcc};
	}

	/** Read (score,label) directly from two columns. */
	private Accum[] loadFromPairs(){
		final int need=Math.max(scoreCol, labelCol)+1;
		final Accum fitAcc=new Accum(fitBins), evalAcc=new Accum(fitBins);
		final ByteFile bf=ByteFile.makeByteFile(in, true);
		byte[] line;
		long row=0;
		while((line=bf.nextLine())!=null){
			if(line.length==0 || line[0]=='#'){continue;}
			final String[] cols=new String(line).split("\t", -1);
			if(cols.length<need){
				throw new RuntimeException("Row has "+cols.length+" cols; need >="+need
						+" (scorecol="+scoreCol+", labelcol="+labelCol+"). For net mode, pass net=.");
			}
			final float raw=Float.parseFloat(cols[scoreCol]);
			final boolean lab=Float.parseFloat(cols[labelCol])>=0.5f;
			if(holdout>=2 && row%holdout==0){evalAcc.add(raw, lab);}else{fitAcc.add(raw, lab);}
			row++;
		}
		bf.close();
		return new Accum[]{fitAcc, evalAcc};
	}

	/*---------- histogram of (rawScore -> count,positives) ----------*/

	/** Fills fixed [0,1) bins while streaming; finish() compacts to non-empty bins. */
	private static final class Accum {
		Accum(int b){B=b; count=new long[b]; pos=new long[b]; sumRaw=new double[b];}
		void add(float raw, boolean label){
			double r=Math.min(0.9999999, Math.max(0, raw));
			int bi=(int)(r*B); if(bi>=B){bi=B-1;}
			count[bi]++; sumRaw[bi]+=raw; totalN++;
			if(label){pos[bi]++; totalPos++;}
		}
		Hist finish(){return new Hist(B, count, pos, sumRaw, totalN, totalPos);}
		/** Element-wise union of two same-width accumulators; neither input is modified. */
		static Accum merge(Accum x, Accum y){
			assert(x.B==y.B) : "Accumulator width mismatch: "+x.B+" vs "+y.B;
			final Accum z=new Accum(x.B);
			for(int i=0; i<x.B; i++){
				z.count[i]=x.count[i]+y.count[i];
				z.pos[i]=x.pos[i]+y.pos[i];
				z.sumRaw[i]=x.sumRaw[i]+y.sumRaw[i];
			}
			z.totalN=x.totalN+y.totalN; z.totalPos=x.totalPos+y.totalPos;
			return z;
		}
		final int B; final long[] count, pos; final double[] sumRaw; long totalN, totalPos;
	}

	/** Compacted histogram: one entry per non-empty bin (meanRaw, count, positives), ascending in meanRaw. */
	private static final class Hist {
		Hist(int B, long[] c, long[] p, double[] s, long tN, long tP){
			int m=0; for(int i=0; i<B; i++){if(c[i]>0){m++;}}
			meanRaw=new float[m]; count=new long[m]; pos=new long[m];
			int j=0;
			for(int i=0; i<B; i++){
				if(c[i]>0){meanRaw[j]=(float)(s[i]/c[i]); count[j]=c[i]; pos[j]=p[i]; j++;}
			}
			totalN=tN; totalPos=tP;
		}
		final float[] meanRaw; final long[] count, pos; final long totalN, totalPos;
	}

	/*---------- calibration map ----------*/

	/** p = K * sigmoid(a*logit(x)+b)^c ; identical to clade.CladeConfidence.calibrate. */
	static float calibrate(float x, float[] p){
		x=Tools.mid(0.0001f, x, 0.9999f);
		final double lx=Math.log(x/(1.0-x));
		final double s=1.0/(1.0+Math.exp(-(p[1]*lx+p[2])));
		return (float)Math.min(1.0, p[0]*Math.pow(s, p[3]));
	}

	/*---------- nonparametric alternative: isotonic (PAV) lookup table ----------*/

	/** Monotone lookup table over raw score, linearly interpolated between knots. */
	static final class Lut {
		Lut(float[] x_, float[] y_){x=x_; y=y_;}
		float apply(float raw){
			final int last=x.length-1;
			if(raw<=x[0]){return y[0];}
			if(raw>=x[last]){return y[last];}
			int lo=0, hi=last;
			while(lo+1<hi){
				final int mid=(lo+hi)>>>1;
				if(x[mid]<=raw){lo=mid;}else{hi=mid;}
			}
			final float dx=x[hi]-x[lo];
			final float f=(dx<=0 ? 0 : (raw-x[lo])/dx);
			return y[lo]+f*(y[hi]-y[lo]);
		}
		final float[] x, y;
	}

	/**
	 * Pool-adjacent-violators isotonic regression of empirical rate on raw score.
	 * Hist bins are already ascending in meanRaw, so a single left-to-right pass
	 * with backward pooling gives the least-squares monotone fit.
	 */
	private static Lut buildLut(Hist h){
		final int m=h.meanRaw.length;
		//Track pooled positives and counts rather than rates, so the final value can be
		//smoothed from the block's actual sample size.
		final double[] sPos=new double[m], sCnt=new double[m];
		final int[] end=new int[m];
		int nb=0;
		for(int k=0; k<m; k++){
			sPos[nb]=h.pos[k]; sCnt[nb]=h.count[k]; end[nb]=k; nb++;
			while(nb>1 && sPos[nb-2]/sCnt[nb-2] > sPos[nb-1]/sCnt[nb-1]){
				sPos[nb-2]+=sPos[nb-1]; sCnt[nb-2]+=sCnt[nb-1];
				end[nb-2]=end[nb-1]; nb--;
			}
		}
		final float[] xs=new float[m], ys=new float[m];
		int start=0;
		for(int b=0; b<nb; b++){
			//Jeffreys smoothing: a block that happened to observe zero positives means
			//"very unlikely", never "impossible". An exact 0 or 1 asserts certainty that
			//finite data cannot support, and downstream it prints as a 0.0% confidence
			//that no amount of contrary evidence can move. Negligible for large blocks.
			final float v=(float)((sPos[b]+0.5)/(sCnt[b]+1.0));
			for(int k=start; k<=end[b]; k++){xs[k]=h.meanRaw[k]; ys[k]=v;}
			start=end[b]+1;
		}
		return new Lut(xs, ys);
	}

	/*---------- fit (Nelder-Mead on log-loss, in unconstrained space) ----------*/

	/** Fits [K,a,b,c]. Optimizes unconstrained z=[logit(K),ln(a),b,ln(c)] so K in (0,1), a>0, c>0. */
	private float[] fit(Hist h){
		float[] best=null;
		double bestLoss=Double.MAX_VALUE;
		final double[][] starts={
			{0, 0, 0, 0}, {4, 0, 0, 0}, {2, 0.7, 0, 0}, {4, 0, -1, 0.4}, {1, -0.5, 1, -0.4},
		};
		final int R=Math.min(restarts, starts.length);
		for(int r=0; r<R; r++){
			final float[] p=toParams(nelderMead(h, starts[r]));
			final double loss=logLossCal(h, p);
			if(loss<bestLoss){bestLoss=loss; best=p;}
		}
		return best;
	}

	private static float[] toParams(double[] z){
		final double K=1.0/(1.0+Math.exp(-z[0]));
		return new float[]{(float)K, (float)Math.exp(z[1]), (float)z[2], (float)Math.exp(z[3])};
	}

	/** Standard Nelder-Mead simplex minimizing log-loss over the 4 unconstrained params. */
	private double[] nelderMead(Hist h, double[] start){
		final int N=4;
		final double[][] simplex=new double[N+1][];
		final double[] fval=new double[N+1];
		simplex[0]=start.clone();
		for(int i=0; i<N; i++){
			final double[] v=start.clone();
			v[i]+=(v[i]==0 ? 0.5 : 0.05*Math.abs(v[i])+0.5);
			simplex[i+1]=v;
		}
		for(int i=0; i<=N; i++){fval[i]=logLossCal(h, toParams(simplex[i]));}

		final double alpha=1.0, gamma=2.0, rho=0.5, sigma=0.5;
		for(int iter=0; iter<maxIters; iter++){
			for(int i=1; i<=N; i++){
				final double fk=fval[i]; final double[] sk=simplex[i]; int j=i-1;
				while(j>=0 && fval[j]>fk){fval[j+1]=fval[j]; simplex[j+1]=simplex[j]; j--;}
				fval[j+1]=fk; simplex[j+1]=sk;
			}
			if(Math.abs(fval[N]-fval[0])<=1e-9*(Math.abs(fval[0])+1e-12)){break;}

			final double[] cen=new double[N];
			for(int i=0; i<N; i++){for(int d=0; d<N; d++){cen[d]+=simplex[i][d];}}
			for(int d=0; d<N; d++){cen[d]/=N;}

			final double[] xr=new double[N];
			for(int d=0; d<N; d++){xr[d]=cen[d]+alpha*(cen[d]-simplex[N][d]);}
			final double fr=logLossCal(h, toParams(xr));

			if(fr<fval[0]){
				final double[] xe=new double[N];
				for(int d=0; d<N; d++){xe[d]=cen[d]+gamma*(xr[d]-cen[d]);}
				final double fe=logLossCal(h, toParams(xe));
				if(fe<fr){simplex[N]=xe; fval[N]=fe;}else{simplex[N]=xr; fval[N]=fr;}
			}else if(fr<fval[N-1]){
				simplex[N]=xr; fval[N]=fr;
			}else{
				final double[] xc=new double[N];
				for(int d=0; d<N; d++){xc[d]=cen[d]+rho*(simplex[N][d]-cen[d]);}
				final double fc=logLossCal(h, toParams(xc));
				if(fc<fval[N]){simplex[N]=xc; fval[N]=fc;}
				else{
					for(int i=1; i<=N; i++){
						for(int d=0; d<N; d++){simplex[i][d]=simplex[0][d]+sigma*(simplex[i][d]-simplex[0][d]);}
						fval[i]=logLossCal(h, toParams(simplex[i]));
					}
				}
			}
		}
		int bi=0;
		for(int i=1; i<=N; i++){if(fval[i]<fval[bi]){bi=i;}}
		return simplex[bi];
	}

	/*---------- predictions (one probability per occupied histogram bin) ----------*/

	private static double[] predictRaw(Hist h){
		final double[] p=new double[h.meanRaw.length];
		for(int k=0; k<p.length; k++){p[k]=h.meanRaw[k];}
		return p;
	}

	private static double[] predictCal(Hist h, float[] cal){
		final double[] p=new double[h.meanRaw.length];
		for(int k=0; k<p.length; k++){p[k]=calibrate(h.meanRaw[k], cal);}
		return p;
	}

	private static double[] predictLut(Hist h, Lut lut){
		final double[] p=new double[h.meanRaw.length];
		for(int k=0; k<p.length; k++){p[k]=lut.apply(h.meanRaw[k]);}
		return p;
	}

	/*---------- metrics (over the compacted histogram) ----------*/

	/**
	 * Allocation-free log-loss for the Nelder-Mead inner loop; the reporting path
	 * uses the double[] form instead so it can score any model, not just the
	 * 4-parameter one.
	 */
	private static double logLossCal(Hist h, float[] cal){
		double sum=0;
		for(int k=0; k<h.meanRaw.length; k++){
			double p=calibrate(h.meanRaw[k], cal);
			p=Math.min(1.0-1e-7, Math.max(1e-7, p));
			sum+=-h.pos[k]*Math.log(p)-(h.count[k]-h.pos[k])*Math.log(1.0-p);
		}
		return sum/h.totalN;
	}

	/** Count-weighted mean negative log-likelihood of an arbitrary per-bin prediction. */
	private static double logLoss(Hist h, double[] pred){
		double sum=0;
		for(int k=0; k<h.meanRaw.length; k++){
			final double p=Math.min(1.0-1e-7, Math.max(1e-7, pred[k]));
			sum+=-h.pos[k]*Math.log(p)-(h.count[k]-h.pos[k])*Math.log(1.0-p);
		}
		return sum/h.totalN;
	}

	/** Expected calibration error over equal-width probability bins, count-weighted. */
	private double ece(Hist h, double[] pred){
		final double[] sp=new double[bins], sl=new double[bins];
		final long[] cnt=new long[bins];
		for(int k=0; k<h.meanRaw.length; k++){
			final double p=Math.min(0.9999999, Math.max(0, pred[k]));
			int bi=(int)(p*bins); if(bi>=bins){bi=bins-1;}
			sp[bi]+=p*h.count[k]; sl[bi]+=h.pos[k]; cnt[bi]+=h.count[k];
		}
		double err=0;
		for(int i=0; i<bins; i++){
			if(cnt[i]==0){continue;}
			err+=(cnt[i]/(double)h.totalN)*Math.abs(sp[i]/cnt[i]-sl[i]/cnt[i]);
		}
		return err;
	}

	/**
	 * Per-bin reliability diagram, binned by CALIBRATED probability. A scalar ECE
	 * averages the gap over bins weighted by population, so a large error in a
	 * sparse high-confidence bin -- exactly where a confidence number gets trusted
	 * -- barely moves it. This table shows that bin explicitly.
	 */
	private String reliabilityTable(Hist h, double[] pCal, double[] pRaw){
		final double[] sp=new double[bins], sr=new double[bins], sl=new double[bins];
		final long[] cnt=new long[bins];
		for(int k=0; k<h.meanRaw.length; k++){
			final double p=Math.min(0.9999999, Math.max(0, pCal[k]));
			int bi=(int)(p*bins); if(bi>=bins){bi=bins-1;}
			sp[bi]+=p*h.count[k]; sr[bi]+=pRaw[k]*h.count[k]; sl[bi]+=h.pos[k]; cnt[bi]+=h.count[k];
		}
		final StringBuilder sb=new StringBuilder();
		sb.append("#reliability (binned by calibrated probability)\n");
		sb.append("#bin\tlo\thi\tn\tfrac\tmeanRaw\tmeanPred\tempirical\tgap\n");
		boolean first=true;
		for(int i=0; i<bins; i++){
			if(cnt[i]==0){continue;}
			final double mp=sp[i]/cnt[i], mr=sr[i]/cnt[i], emp=sl[i]/(double)cnt[i];
			if(!first){sb.append('\n');}
			first=false;
			sb.append(i).append('\t');
			sb.append(String.format("%.4f\t%.4f\t", i/(double)bins, (i+1)/(double)bins));
			sb.append(cnt[i]).append('\t');
			sb.append(String.format("%.6f\t%.6f\t%.6f\t%.6f\t%+.6f",
					cnt[i]/(double)h.totalN, mr, mp, emp, mp-emp));
		}
		return sb.toString();
	}

	/*---------- fields ----------*/

	private String in=null;
	private String netFile=null;
	private String out=null;
	private String relOut=null;
	private String lutOut=null;
	private String label=null;
	private int scoreCol=0;
	private int labelCol=1;
	private int bins=20;
	private int fitBins=4096;
	private int holdout=5;
	private int maxIters=2000;
	private int restarts=5;
	private java.io.PrintStream outstream=System.err;
}
