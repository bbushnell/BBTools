package cardinality;

/**
 * Fll53MLE: composite-likelihood cardinality extraction for Fll52, using the
 * (G, m)-conditioned class table from Fll53Calibrate.  Identical machinery to
 * Fll2MLE with the 56-class 5-antenna geometry.
 *
 * @author Amber
 * @date July 2026
 */
public class Fll53MLE {

	static final int NC=Fll53Calibrate.NC;
	static final int MAX_G=32;
	static final double ALPHA=0.5;

	/*---------------- Construction (path may be .gz) ----------------*/

	public Fll53MLE(String tablePath) throws Exception{
		final java.util.List<String> lines=readLines(tablePath);
		final java.util.TreeMap<Integer, java.util.ArrayList<double[]>> rowsByG=new java.util.TreeMap<>();
		final java.util.TreeMap<Integer, java.util.ArrayList<Double>> centersByG=new java.util.TreeMap<>();
		for(String line : lines){
			if(line.isEmpty() || line.charAt(0)=='#'){continue;}
			final String[] f=line.split("\t");
			final int g=Integer.parseInt(f[0]);
			final double total=Double.parseDouble(f[3]);
			final double[] lp=new double[NC];
			final double denom=total+ALPHA*NC;
			for(int k=0; k<NC; k++){
				lp[k]=Math.log((Double.parseDouble(f[4+k])+ALPHA)/denom);
			}
			rowsByG.computeIfAbsent(g, x->new java.util.ArrayList<>()).add(lp);
			centersByG.computeIfAbsent(g, x->new java.util.ArrayList<>())
				.add(Double.parseDouble(f[2]));
		}
		for(int g : rowsByG.keySet()){
			final java.util.ArrayList<double[]> rows=rowsByG.get(g);
			final java.util.ArrayList<Double> centers=centersByG.get(g);
			mCenterByG[g]=new double[centers.size()];
			logPByG[g]=new double[rows.size()][];
			for(int i=0; i<rows.size(); i++){
				mCenterByG[g][i]=centers.get(i);
				logPByG[g][i]=rows.get(i);
			}
		}
	}

	private final double[][] mCenterByG=new double[MAX_G][];
	private final double[][][] logPByG=new double[MAX_G][][];

	/** Reads a text table, transparently handling .gz (via fileIO.ByteFile). */
	private static java.util.List<String> readLines(String path) throws Exception{
		final fileIO.FileFormat ff=fileIO.FileFormat.testInput(path, fileIO.FileFormat.TEXT, null, false, true);
		final fileIO.ByteFile bf=fileIO.ByteFile.makeByteFile(ff, 1);
		final java.util.ArrayList<String> lines=new java.util.ArrayList<>();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			lines.add(new String(line));
		}
		bf.close();
		return lines;
	}

	public double estimate(Fll53 f){
		final int nw=f.getNumOrgs();
		final int[] counts=new int[NC];
		for(int w=0; w<nw; w++){
			counts[Fll53Calibrate.classOf(f.getOrg(w))]++;
		}
		final int G=f.getGlobalExp();
		int g=Math.min(G, MAX_G-1);
		int lo=g, hi=g;
		while(logPByG[g]==null){
			lo--; hi++;
			if(lo>=0 && logPByG[lo]!=null){g=lo; break;}
			if(hi<MAX_G && logPByG[hi]!=null){g=hi; break;}
			if(lo<0 && hi>=MAX_G){throw new RuntimeException("empty table");}
		}
		final double[] mCenter=mCenterByG[g];
		final double[][] logP=logPByG[g];

		int best=0;
		double bestL=Double.NEGATIVE_INFINITY;
		final double[] L=new double[mCenter.length];
		for(int b=0; b<mCenter.length; b++){
			double s=0;
			final double[] lp=logP[b];
			for(int k=0; k<NC; k++){
				if(counts[k]!=0){s+=counts[k]*lp[k];}
			}
			L[b]=s;
			if(s>bestL){bestL=s; best=b;}
		}
		double m=mCenter[best];
		if(best>0 && best<mCenter.length-1){
			final double lm=L[best-1], l0=L[best], lpn=L[best+1];
			final double d=(lm-2*l0+lpn);
			if(d<0){
				final double off=0.5*(lm-lpn)/d;
				if(Math.abs(off)<=1.0){
					m=mCenter[best]+off*(mCenter[best+1]-mCenter[best]);
				}
			}
		}
		return nw*5.0*Math.pow(2.0, G+m);
	}
}
