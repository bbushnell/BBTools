package prot;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

import ml.CellNet;
import ml.CellNetParser;

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
		BufferedReader br=new BufferedReader(new FileReader(in));
		String header=br.readLine();
		String[] h=header.split("\t");
		final int numIn=Integer.parseInt(h[1]);
		final int numOut=Integer.parseInt(h[2]);

		double[] sae=new double[numOut], sse=new double[numOut];
		ArrayList<double[]> trues=new ArrayList<double[]>();
		float[] input=new float[numIn];
		String line;
		long n=0;
		while((line=br.readLine())!=null){
			if(line.length()==0){continue;}
			String[] f=line.split("\t");
			if(f.length<numIn+numOut){continue;}
			for(int i=0; i<numIn; i++){input[i]=Float.parseFloat(f[i]);}
			net.applyInput(input);
			net.feedForward();
			double[] t=new double[numOut];
			for(int k=0; k<numOut; k++){
				double pred=net.getOutput(k);
				double tru=Double.parseDouble(f[numIn+k]);
				t[k]=tru;
				double err=pred-tru;
				sae[k]+=Math.abs(err);
				sse[k]+=err*err;
			}
			trues.add(t);
			n++;
		}
		br.close();

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
