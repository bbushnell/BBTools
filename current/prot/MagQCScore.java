package prot;

import java.util.ArrayList;

import fileIO.ByteFile;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;

/**
 * Scores a trained MAG-QC regression net against a held-out validation TSV and
 * reports per-output MAE/RMSE alongside a predict-the-mean baseline, so the net's
 * accuracy is judged honestly (does it beat guessing the average?) rather than by
 * a single summed MSE.
 *
 * Usage: java prot.MagQCScore net=magqc_bac.bbnet in=val_bac.tsv
 * The TSV is the same format the trainer consumes: '#dims <in> <out>' header, then
 * rows of <in floats> <out floats>. Assumes 2 outputs (completeness, contamination).
 */
public class MagQCScore {

	public static void main(String[] args) throws Exception{
		String netFile=null, in=null;
		for(String a : args){
			int e=a.indexOf('=');
			if(e<0){continue;}
			String k=a.substring(0, e), v=a.substring(e+1);
			if(k.equals("net")){netFile=v;}
			else if(k.equals("in")){in=v;}
		}
		if(netFile==null || in==null){throw new RuntimeException("Required: net= in=");}

		CellNet net=CellNetParser.load(netFile);
		final ByteFile bf=ByteFile.makeByteFile(in, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		final byte[] header=bf.nextLine();
		lp.set(header);
		final int numIn=lp.parseInt(1);
		final int numOut=lp.parseInt(2);

		double[] sae=new double[numOut], sse=new double[numOut];
		ArrayList<double[]> trues=new ArrayList<double[]>();
		float[] input=new float[numIn];
		long n=0, skipped=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				lp.set(line);
				if(lp.terms()<numIn+numOut){
					//Malformed row (not just a blank line): logged, not silently dropped -- a
					//large skip count would otherwise hide most of the validation set vanishing.
					skipped++;
					continue;
				}
				for(int i=0; i<numIn; i++){input[i]=lp.parseFloat(i);}
				net.applyInput(input);
				net.feedForward();
				double[] t=new double[numOut];
				for(int k=0; k<numOut; k++){
					double pred=net.getOutput(k);
					double tru=lp.parseDouble(numIn+k);
					t[k]=tru;
					double err=pred-tru;
					sae[k]+=Math.abs(err);
					sse[k]+=err*err;
				}
				trues.add(t);
				n++;
			}
		}finally{
			bf.close();
		}
		if(skipped>0){System.err.println("WARNING: skipped "+skipped+" malformed row(s) in "+in);}
		if(n==0){throw new RuntimeException("No valid rows scored from "+in+" ("+skipped+" skipped).");}

		// predict-the-mean baseline (best constant predictor = the mean of the val targets)
		double[] mean=new double[numOut];
		for(double[] t : trues){for(int k=0; k<numOut; k++){mean[k]+=t[k];}}
		for(int k=0; k<numOut; k++){mean[k]/=Math.max(1, n);}
		double[] bae=new double[numOut], bse=new double[numOut];
		for(double[] t : trues){
			for(int k=0; k<numOut; k++){
				double e=mean[k]-t[k];
				bae[k]+=Math.abs(e);
				bse[k]+=e*e;
			}
		}

		String[] names={"completeness", "contamination"};
		System.out.println("Scored "+n+" held-out validation samples ("+netFile+")");
		System.out.printf("%-14s  %8s %8s   %8s %8s   %6s%n", "output", "MAE", "RMSE", "baseMAE", "baseRMSE", "R2");
		for(int k=0; k<numOut; k++){
			double mae=sae[k]/n, rmse=Math.sqrt(sse[k]/n);
			double bmae=bae[k]/n, brmse=Math.sqrt(bse[k]/n);
			double r2=1.0-(sse[k]/Math.max(1e-12, bse[k]));//1 - net_SSE/baseline_SSE
			String nm=(k<names.length ? names[k] : "out"+k);
			System.out.printf("%-14s  %8.4f %8.4f   %8.4f %8.4f   %6.3f%n", nm, mae, rmse, bmae, brmse, r2);
		}
	}
}
