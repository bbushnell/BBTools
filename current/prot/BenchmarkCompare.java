package prot;

import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ByteFile;
import parse.LineParser1;

/**
 * Compares a genome-QC tool's per-bin completeness/contamination predictions against the benchmark
 * ground truth, reporting MAE / RMSE / R2 (vs predict-the-mean baseline) / bias for each target.
 * Truth and predictions are joined on bin id, so CheckM1, CheckM2, and our net are scored the SAME
 * way on the SAME bins. The CheckM outputs are first normalized (by a small awk) to a 3-column
 * (binID, completeness, contamination) TSV; this tool does the metric math (Java, per project RULES).
 *
 * truth=     bench_lean_truth.tsv  (binID, tid, completeness, contamination, ...; 0-1 fractions)
 * pred=      normalized TSV        (binID, completeness, contamination)
 * predscale= multiply pred values by this (0.01 for CheckM's 0-100 percentages; default 1.0)
 * name=      label for the output rows (e.g. checkm2)
 * @author UMP45
 */
public class BenchmarkCompare {

	public static void main(String[] args){
		String truthFile=null, predFile=null, name="tool";
		double predScale=1.0;
		for(final String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<0){continue;}
			final String a=arg.substring(0, eq).toLowerCase(), b=arg.substring(eq+1);
			if(a.equals("truth")){truthFile=b;}
			else if(a.equals("pred")){predFile=b;}
			else if(a.equals("predscale")){predScale=Double.parseDouble(b);}
			else if(a.equals("name")){name=b;}
		}
		if(truthFile==null || predFile==null){throw new RuntimeException("Required: truth= pred=");}

		// truth: binID -> {completeness, contamination}
		final HashMap<String,double[]> truth=new HashMap<String,double[]>();
		{
			final ByteFile bf=ByteFile.makeByteFile(truthFile, true);
			final LineParser1 lp=new LineParser1((byte)'\t');
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(lp.terms()<4){continue;}
				truth.put(lp.parseString(0), new double[]{lp.parseDouble(2), lp.parseDouble(3)});
			}
			bf.close();
		}
		assert(!truth.isEmpty()) : "no truth rows in "+truthFile;

		// pred: join on binID, accumulate errors
		double sae0=0, sae1=0, sse0=0, sse1=0, sbias0=0, sbias1=0;
		int n=0, unmatched=0;
		final ArrayList<double[]> tvals=new ArrayList<double[]>();
		{
			final ByteFile bf=ByteFile.makeByteFile(predFile, true);
			final LineParser1 lp=new LineParser1((byte)'\t');
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(lp.terms()<3){continue;}
				final double[] t=truth.get(lp.parseString(0));
				if(t==null){unmatched++; continue;}
				final double e0=lp.parseDouble(1)*predScale-t[0], e1=lp.parseDouble(2)*predScale-t[1];
				sae0+=Math.abs(e0); sae1+=Math.abs(e1);
				sse0+=e0*e0; sse1+=e1*e1;
				sbias0+=e0; sbias1+=e1;
				tvals.add(t);
				n++;
			}
			bf.close();
		}
		if(n==0){throw new RuntimeException("No bins joined between "+truthFile+" and "+predFile
			+" (check binID match; "+unmatched+" pred bins unmatched)");}

		// predict-the-mean baseline SSE for R2
		double m0=0, m1=0;
		for(final double[] t : tvals){m0+=t[0]; m1+=t[1];}
		m0/=n; m1/=n;
		double b0=0, b1=0;
		for(final double[] t : tvals){b0+=(m0-t[0])*(m0-t[0]); b1+=(m1-t[1])*(m1-t[1]);}
		final double r2c=1.0-sse0/Math.max(1e-12, b0), r2k=1.0-sse1/Math.max(1e-12, b1);

		System.out.printf("%-10s %-14s %8s %8s %8s %8s%n", "tool", "target", "MAE", "RMSE", "R2", "bias");
		System.out.printf("%-10s %-14s %8.4f %8.4f %8.3f %+8.4f%n", name, "completeness", sae0/n, Math.sqrt(sse0/n), r2c, sbias0/n);
		System.out.printf("%-10s %-14s %8.4f %8.4f %8.3f %+8.4f%n", name, "contamination", sae1/n, Math.sqrt(sse1/n), r2k, sbias1/n);
		System.err.println(name+": joined "+n+" bins ("+unmatched+" pred bins unmatched to truth)");
	}
}
