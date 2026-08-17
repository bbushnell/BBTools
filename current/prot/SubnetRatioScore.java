package prot;

import java.util.ArrayList;

import fileIO.ByteFile;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;

/**
 * Evaluates a "denominator" subnet by the quality of the DERIVED subset-relative
 * completeness, not just its raw target. Such a subnet predicts a gene-set's native
 * COMPLEMENT (e.g. an organism's total ncRNA count); the useful quantity is
 * completeness = observed / predicted-native. This tool reports:
 *
 * <ul>
 * <li>the raw native-count R2 (how well the denominator itself is predicted), and</li>
 * <li>the derived-completeness R2: R2 of (observed / predicted-native, clamped to
 * [0,1]) against the TRUE completeness (observed / true-native), on the same rows.</li>
 * </ul>
 *
 * The first {@code numobs} input columns are the observed counts (their sum is the
 * numerator); the single output is the predicted native complement; the first target
 * column is the TRUE native complement.
 *
 * Usage: java prot.SubnetRatioScore net=net_ncrna.bbnet in=ncrna_val.tsv numobs=5
 *
 * @author Eru
 */
public class SubnetRatioScore {

	public static void main(String[] args) throws Exception{
		String netFile=null, in=null;
		int numObs=5;
		for(String a : args){
			int e=a.indexOf('=');
			if(e<0){continue;}
			String k=a.substring(0, e), v=a.substring(e+1);
			if(k.equals("net")){netFile=v;}
			else if(k.equals("in")){in=v;}
			else if(k.equals("numobs")){numObs=Integer.parseInt(v);}
		}
		if(netFile==null || in==null){throw new RuntimeException("Required: net= in= [numobs=5]");}

		CellNet net=CellNetParser.load(netFile);
		final ByteFile bf=ByteFile.makeByteFile(in, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		lp.set(bf.nextLine());
		final int numIn=lp.parseInt(1);
		final int numOut=lp.parseInt(2);

		ArrayList<double[]> rows=new ArrayList<double[]>();//{trueNat, predNat, trueComp, predComp}
		float[] input=new float[numIn];
		double sseN=0, sseC=0, maeC=0;
		long skipped=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				lp.set(line);
				if(lp.terms()<numIn+numOut){
					//Malformed row: logged, not silently dropped (same fix as MagQCScore).
					skipped++;
					continue;
				}
				double obs=0;
				for(int i=0; i<numIn; i++){input[i]=lp.parseFloat(i); if(i<numObs){obs+=input[i];}}
				net.applyInput(input);
				net.feedForward();
				double pnat=net.getOutput(0);
				double tnat=lp.parseDouble(numIn);//first output column = native target
				double tcomp=(tnat>0 ? obs/tnat : 0);
				double pcomp=Math.min(1.0, obs/Math.max(0.5, pnat));
				rows.add(new double[]{tnat, pnat, tcomp, pcomp});
				sseN+=(pnat-tnat)*(pnat-tnat);
				sseC+=(pcomp-tcomp)*(pcomp-tcomp);
				maeC+=Math.abs(pcomp-tcomp);
			}
		}finally{
			bf.close();
		}
		if(skipped>0){System.err.println("WARNING: skipped "+skipped+" malformed row(s) in "+in);}
		if(rows.isEmpty()){throw new RuntimeException("No valid rows scored from "+in+" ("+skipped+" skipped).");}

		final int n=rows.size();
		double mN=0, mC=0;
		for(double[] r : rows){mN+=r[0]; mC+=r[2];}
		mN/=Math.max(1, n); mC/=Math.max(1, n);
		double bseN=0, bseC=0;
		for(double[] r : rows){bseN+=sq(r[0]-mN); bseC+=sq(r[2]-mC);}
		double r2N=1.0-sseN/Math.max(1e-12, bseN);
		double r2C=1.0-sseC/Math.max(1e-12, bseC);

		System.out.printf("Scored %d rows (%s)%n", n, netFile);
		System.out.printf("native-count:                                 R2=%.3f  RMSE=%.3f%n",
			r2N, Math.sqrt(sseN/Math.max(1, n)));
		System.out.printf("derived completeness (obs/pred, clamp[0,1]):  R2=%.3f  MAE=%.4f  RMSE=%.4f%n",
			r2C, maeC/Math.max(1, n), Math.sqrt(sseC/Math.max(1, n)));
		System.out.printf("(true completeness mean=%.3f; native mean=%.2f)%n", mC, mN);
	}

	static double sq(double x){return x*x;}
}
