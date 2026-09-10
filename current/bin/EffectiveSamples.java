package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.LongList;

/**
 * Measurement scaffold: estimates the number of LOGICAL (effective) samples in a
 * QuickBin cov file, for auto-calibrating CovMaker's condense target.
 * Computes the sample-sample correlation matrix of log-depths over large contigs,
 * then several effective-count estimators:
 *   PR    - eigenvalue participation ratio (sum L)^2/(sum L^2)
 *   erank - spectral-entropy effective rank exp(-sum p ln p), p=L/sum  [Roy & Vetterli 2007]
 *   v95   - eigenvalues needed to reach 95% of trace
 *   clusters(t) - single-linkage cluster count at correlation threshold t
 *   nGap  - single-linkage merge-distance elbow (the CovMaker-integrable stopping rule)
 * Ground-truth validation data: RandomReadsMG depthseed/jitter arms + CovCorrSim synthetics
 * (see plans/sample_reduction_design_2026-09-09.md).
 * @author Amber
 * @date September 9, 2026
 */
public class EffectiveSamples {

	public static void main(String[] args){
		Timer t=new Timer();
		EffectiveSamples x=new EffectiveSamples(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public EffectiveSamples(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("in")){in=b;}
			else if(a.equals("minsize")){minSize=Parse.parseIntKMG(b);}
			else if(a.equals("weighted")){weighted=Parse.parseBoolean(b);}
			else if(a.equals("verbose")){verbose=Parse.parseBoolean(b);}
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		assert(in!=null) : "in=<cov file> is required";
	}

	void process(Timer t){
		loadCov();
		final int n=numDepths;
		assert(n>=2) : "Need >=2 samples, found "+n;
		double[][] corr=correlationMatrix();
		printMatrix(corr, "Correlation");
		double[] eig=jacobiEigenvalues(corr);
		Arrays.sort(eig);
		for(int i=0, j=eig.length-1; i<j; i++, j--){double tmp=eig[i]; eig[i]=eig[j]; eig[j]=tmp;}
		outstream.println("Eigenvalues:\t"+fmt(eig));
		//Pairwise-deleted correlation matrices can be slightly non-PSD (small negative
		//eigenvalues); assert on the RAW sum (trace is exact), then clamp for the stats.
		double rawSum=0;
		for(double L : eig){rawSum+=L;}
		assert(Math.abs(rawSum-n)<0.01*n) : "Eigenvalue sum "+rawSum+" != trace "+n;
		double sum=0, sum2=0;
		for(double L : eig){L=Math.max(L, 0); sum+=L; sum2+=L*L;}
		double pr=sum*sum/sum2;
		double ent=0;
		for(double L : eig){
			if(L>1e-12){double p=L/sum; ent-=p*Math.log(p);}
		}
		double erank=Math.exp(ent);
		int v95=0;
		double acc=0;
		for(double L : eig){acc+=L; v95++; if(acc>=0.95*sum){break;}}
		outstream.println("PR:   \t"+String.format("%.3f", pr));
		outstream.println("erank:\t"+String.format("%.3f", erank));
		outstream.println("v95:  \t"+v95);
		for(double tau : new double[] {0.70, 0.80, 0.90, 0.95}){
			outstream.println("clusters(r>"+tau+"):\t"+clusterCount(corr, tau));
		}
		outstream.println("nGap: \t"+gapEstimate(corr));
		t.stop("Total time:");
	}

	private void loadCov(){
		ByteFile bf=ByteFile.makeByteFile(in, true);
		LineParser1 lp=new LineParser1('\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length<1){continue;}
			if(line[0]=='#'){
				lp.set(line);
				if(lp.termEquals("#Depths", 0)){numDepths=lp.parseInt(1);}
				continue;
			}
			lp.set(line);
			//Format: ShortName ID Size Cov_0..Cov_{d-1} [Edge Weight ...]
			long size=lp.parseLong(2);
			if(size<minSize){continue;}
			float[] row=new float[numDepths];
			for(int s=0; s<numDepths; s++){row[s]=lp.parseFloat(3+s);}
			rows.add(row);
			weights.add(weighted ? size : 1L);
		}
		bf.close();
		outstream.println("Loaded "+rows.size()+" contigs >= "+minSize+" bp, "+numDepths+" samples.");
		assert(rows.size()>=100) : "Too few contigs ("+rows.size()+") for a stable estimate; lower minsize.";
	}

	/**
	 * Weighted PAIRWISE Pearson correlation of log(depth+0.5) across contigs.
	 * Rows where BOTH samples are absent (depth==0) are EXCLUDED from that pair:
	 * shared absence is not evidence of correlation (NEON-like data is mostly
	 * zeros in every library — naive global Pearson reads co-absence as high r
	 * and would wrongly condense independent sparse samples). Present-in-one,
	 * absent-in-other rows ARE included: informative discordance.
	 */
	private double[][] correlationMatrix(){
		final int n=numDepths, m=rows.size();
		double[][] x=new double[n][m];
		boolean[][] present=new boolean[n][m];
		double[] w=new double[m];
		for(int c=0; c<m; c++){
			w[c]=weights.get(c);
			float[] row=rows.get(c);
			for(int s=0; s<n; s++){
				x[s][c]=Math.log(row[s]+0.5);
				present[s][c]=(row[s]>0);
			}
		}
		double[][] corr=new double[n][n];
		for(int i=0; i<n; i++){
			corr[i][i]=1;
			for(int j=i+1; j<n; j++){
				double sw=0, sx=0, sy=0, sxx=0, syy=0, sxy=0;
				for(int c=0; c<m; c++){
					if(!present[i][c] && !present[j][c]){continue;}
					final double wc=w[c], xi=x[i][c], yj=x[j][c];
					sw+=wc; sx+=wc*xi; sy+=wc*yj;
					sxx+=wc*xi*xi; syy+=wc*yj*yj; sxy+=wc*xi*yj;
				}
				double r=0;
				if(sw>0){
					double cov=sxy/sw-(sx/sw)*(sy/sw);
					double vx=sxx/sw-(sx/sw)*(sx/sw), vy=syy/sw-(sy/sw)*(sy/sw);
					double d=Math.sqrt(vx*vy);
					r=(d>0 ? cov/d : 0);
				}
				corr[i][j]=corr[j][i]=r;
			}
		}
		return corr;
	}

	/** Cyclic Jacobi eigenvalue iteration; fine for small symmetric matrices (n<=~50). */
	static double[] jacobiEigenvalues(double[][] a0){
		final int n=a0.length;
		double[][] a=new double[n][n];
		for(int i=0; i<n; i++){a[i]=a0[i].clone();}
		for(int sweep=0; sweep<100; sweep++){
			double off=0;
			for(int i=0; i<n; i++){
				for(int j=i+1; j<n; j++){off+=a[i][j]*a[i][j];}
			}
			if(off<1e-18){break;}
			for(int p=0; p<n; p++){
				for(int q=p+1; q<n; q++){
					if(Math.abs(a[p][q])<1e-15){continue;}
					double theta=(a[q][q]-a[p][p])/(2*a[p][q]);
					double t=Math.signum(theta)/(Math.abs(theta)+Math.sqrt(theta*theta+1));
					if(theta==0){t=1;}
					double c=1/Math.sqrt(t*t+1), s=t*c;
					for(int k=0; k<n; k++){
						double akp=a[k][p], akq=a[k][q];
						a[k][p]=c*akp-s*akq;
						a[k][q]=s*akp+c*akq;
					}
					for(int k=0; k<n; k++){
						double apk=a[p][k], aqk=a[q][k];
						a[p][k]=c*apk-s*aqk;
						a[q][k]=s*apk+c*aqk;
					}
				}
			}
		}
		double[] eig=new double[n];
		for(int i=0; i<n; i++){eig[i]=a[i][i];}
		return eig;
	}

	/** Single-linkage: union samples whose correlation exceeds tau; count components. */
	private int clusterCount(double[][] corr, double tau){
		final int n=corr.length;
		int[] parent=new int[n];
		for(int i=0; i<n; i++){parent[i]=i;}
		for(int i=0; i<n; i++){
			for(int j=i+1; j<n; j++){
				if(corr[i][j]>tau){union(parent, i, j);}
			}
		}
		int count=0;
		for(int i=0; i<n; i++){if(find(parent, i)==i){count++;}}
		return count;
	}

	private static int find(int[] p, int x){while(p[x]!=x){p[x]=p[p[x]]; x=p[x];} return x;}
	private static void union(int[] p, int a, int b){p[find(p, a)]=find(p, b);}

	/**
	 * Merge-cost elbow (the CovMaker-integrable stopping rule): single-linkage
	 * agglomeration on d=1-corr records n-1 merge distances in ascending order;
	 * the largest relative jump marks the boundary between within-group and
	 * between-group merges. N = n - (merges before the jump). If no jump exceeds
	 * GAP_MIN (structureless case, all distances comparable), returns n.
	 */
	static int gapEstimate(double[][] corr){
		final int n=corr.length;
		//Kruskal-style single linkage: sort pair distances, count component-joining merges.
		int pairs=n*(n-1)/2;
		double[][] edges=new double[pairs][3];
		int e=0;
		for(int i=0; i<n; i++){
			for(int j=i+1; j<n; j++){edges[e++]=new double[] {1-corr[i][j], i, j};}
		}
		Arrays.sort(edges, (x, y) -> Double.compare(x[0], y[0]));
		int[] parent=new int[n];
		for(int i=0; i<n; i++){parent[i]=i;}
		double[] merges=new double[n-1];
		int m=0;
		for(double[] ed : edges){
			int a=find(parent, (int)ed[1]), b=find(parent, (int)ed[2]);
			if(a!=b){parent[a]=b; merges[m++]=ed[0];}
		}
		assert(m==n-1) : m;
		final double GAP_MIN=3.0, EPS=1e-4;
		int cut=-1;
		double bestRatio=GAP_MIN;
		for(int i=0; i<m-1; i++){
			double ratio=(merges[i+1]+EPS)/(merges[i]+EPS);
			if(ratio>bestRatio){bestRatio=ratio; cut=i;}
		}
		return (cut<0 ? n : n-(cut+1));
	}

	private void printMatrix(double[][] m, String name){
		outstream.println(name+" matrix:");
		for(double[] row : m){
			StringBuilder sb=new StringBuilder();
			for(double v : row){sb.append(String.format("%7.3f", v));}
			outstream.println(sb.toString());
		}
	}

	private static String fmt(double[] a){
		StringBuilder sb=new StringBuilder();
		for(double v : a){sb.append(String.format("%.3f ", v));}
		return sb.toString();
	}

	private String in=null;
	private int minSize=5000;
	private boolean weighted=true;
	private int numDepths=-1;
	private final ArrayList<float[]> rows=new ArrayList<float[]>();
	private final LongList weights=new LongList();

	public static boolean verbose=false;
	private PrintStream outstream=System.err;
}
