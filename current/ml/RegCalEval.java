package ml;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import structures.FloatList;

/**
 * Ranks continuous RegressionTrainer networks by calibration rather than MSE.
 * Loads accepted validation rows once, then reuses them for each net_*.bbnet in a directory.
 * Each tab-delimited row supplies raw features followed by one target, normally in [0,1].
 * Comment headers, including #dims, are ignored; the first data row establishes feature width.
 * Rows with a different field count or unparseable numeric fields are skipped.
 * A positive bin count and at least one accepted validation row are required.
 * Raw inputs suit networks whose input standardization was folded into the exported first layer.
 * Reports MSE and depth-ECE: the population-weighted absolute difference between mean prediction
 * and mean target in equal-width prediction buckets. Only bucket indices are clamped; scores
 * themselves are not clipped. Nonfinite targets or returned predictions stop evaluation;
 * finite targets outside [0,1] remain accepted.
 * MSE and calibration are different properties (Amber): a net can fit well on average yet be
 * miscalibrated, so the fleet is ranked on ECE and the MSE-rank/ECE-rank disagreement is reported.
 *
 * The optional manifest maps network IDs to display-only dimension strings. Results and rank
 * inversions go to stdout; loading diagnostics go to stderr. main owns its readers and
 * evaluates networks serially, with normal-path reader cleanup.
 * Reported reader-completion errors stop evaluation before scoring begins.
 *
 * Usage: regcaleval.sh nets=&lt;dir&gt; val=&lt;val.tsv&gt; [manifest=&lt;tsv&gt;] [bins=20] [maxrows=500000]
 *
 * @author Noire
 */
public class RegCalEval{

	/** Loads at most maxrows accepted rows, scores matching network files, and prints rankings.
	 * A model absent at load time is skipped; parsing or inference exceptions otherwise propagate.
	 * @param args Required nets=/val= paths and optional manifest=/bins=/maxrows= settings */
	public static void main(String[] args){
		String dir=null, val=null, manifest=null;
		int bins=20, maxrows=500000;
		for(String a : args){
			final String[] s=a.split("=", 2);
			final String k=s[0].toLowerCase(), v=(s.length>1 ? s[1] : null);
			if(k.equals("nets") || k.equals("dir")){dir=v;}
			else if(k.equals("val") || k.equals("in")){val=v;}
			else if(k.equals("manifest")){manifest=v;}
			else if(k.equals("bins")){bins=Integer.parseInt(v);}
			else if(k.equals("maxrows")){maxrows=Integer.parseInt(v);}
			else{throw new RuntimeException("Unknown parameter: "+a);}
		}
		if(dir==null || val==null){throw new RuntimeException("nets=<dir> and val=<file> are required");}
		if(bins<1){throw new IllegalArgumentException("bins must be positive: "+bins);}

		//--- load val once: raw features + target ---
		final ArrayList<float[]> X=new ArrayList<float[]>();
		final FloatList Y=new FloatList();
		int nin=-1;
		final ByteFile bf=ByteFile.makeByteFile(val, true);
		byte[] line;
		long lineNumber=0;
		while((line=bf.nextLine())!=null && X.size()<maxrows){
			lineNumber++;
			if(line.length==0 || line[0]=='#'){continue;}
			final String[] f=new String(line).split("\t");
			if(nin<0){nin=f.length-1;}
			if(f.length!=nin+1){continue;}
			final float[] x=new float[nin];
			boolean ok=true;
			for(int i=0; i<nin; i++){try{x[i]=Float.parseFloat(f[i]);}catch(Exception e){ok=false; break;}}
			if(!ok){continue;}
			final float target;
			try{target=Float.parseFloat(f[nin]);}catch(Exception e){continue;}
			if(!Float.isFinite(target)){
				final IllegalArgumentException error=new IllegalArgumentException(
					"Non-finite target in "+val+" at line "+lineNumber+": "+target);
				if(bf.close()){error.addSuppressed(new RuntimeException("Error reading or closing "+val));}
				throw error;
			}
			Y.add(target);
			X.add(x);
		}
		if(bf.close()){throw new RuntimeException("Error reading or closing "+val);}
		final int n=X.size();
		if(n<1){throw new IllegalArgumentException("No usable validation rows in "+val);}
		System.err.println("val rows="+n+"  inputs="+nin);

		//--- manifest: net_id -> dims (for display) ---
		final java.util.HashMap<String,String> dimsById=new java.util.HashMap<String,String>();
		if(manifest!=null){
			final ByteFile mb=ByteFile.makeByteFile(manifest, true);
			while((line=mb.nextLine())!=null){
				if(line.length==0 || line[0]=='#'){continue;}
				final String[] f=new String(line).split("\t");
				if(f.length>=2){dimsById.put(f[0], f[1]);}
			}
			if(mb.close()){throw new RuntimeException("Error reading or closing "+manifest);}
		}

		//--- evaluate each net ---
		final File[] nets=new File(dir).listFiles((d, name)->name.startsWith("net_") && name.endsWith(".bbnet"));
		if(nets==null || nets.length==0){throw new RuntimeException("no net_*.bbnet in "+dir);}
		Arrays.sort(nets);
		final ArrayList<Result> results=new ArrayList<Result>();
		final FloatList fl=new FloatList(nin);
		for(File nf : nets){
			final CellNet net=CellNetParser.load(nf.getPath());
			if(net==null){System.err.println("skip (load failed): "+nf.getName()); continue;}
			final double[] binPred=new double[bins], binTrue=new double[bins];
			final long[] binCount=new long[bins];
			double sse=0;
			for(int r=0; r<n; r++){
				fl.size=0;
				final float[] x=X.get(r);
				for(int i=0; i<nin; i++){fl.add(x[i]);}
				net.applyInput(fl);
				final float p=net.feedForward();
				if(!Float.isFinite(p)){
					throw new IllegalArgumentException("Non-finite prediction for "+nf.getName()+
						" at validation row "+(r+1)+": "+p);
				}
				final float y=Y.get(r);
				final double error=(double)p-y;
				sse+=error*error;
				int b=(int)(p*bins); if(b<0){b=0;} if(b>=bins){b=bins-1;}
				binPred[b]+=p; binTrue[b]+=y; binCount[b]++;
			}
			double ece=0;
			for(int b=0; b<bins; b++){
				if(binCount[b]==0){continue;}
				final double mp=binPred[b]/binCount[b], mt=binTrue[b]/binCount[b];
				ece+=(binCount[b]/(double)n)*Math.abs(mp-mt);
			}
			final String id=nf.getName().replace("net_", "").replace(".bbnet", "");
			results.add(new Result(id, dimsById.getOrDefault(id, "?"), sse/n, ece));
		}

		//--- rank by ECE; report MSE-rank vs ECE-rank divergence ---
		results.sort((a, b)->Double.compare(a.ece, b.ece));
		final PrintStream out=System.out;
		out.println("REGRESSION-CALIBRATION RANKING (lower ECE = better calibrated). MSE shown for contrast.\n");
		out.println(String.format("%-6s %-22s %10s %10s", "net", "dims", "depthECE", "MSE"));
		for(Result r : results){
			out.println(String.format("net_%-2s %-22s %10.5f %10.5f", r.id, r.dims, r.ece, r.mse));
		}
		//rank correlation: how often does MSE order disagree with ECE order?
		final ArrayList<Result> byMse=new ArrayList<Result>(results);
		byMse.sort((a, b)->Double.compare(a.mse, b.mse));
		int inversions=0;
		for(int i=0; i<results.size(); i++){
			for(int j=i+1; j<results.size(); j++){
				if(byMse.indexOf(results.get(i))>byMse.indexOf(results.get(j))){inversions++;}
			}
		}
		final int pairs=results.size()*(results.size()-1)/2;
		out.println(String.format("%nMSE-rank vs ECE-rank inversions: %d/%d (%.1f%%) -- high = MSE is a poor proxy for calibration",
				inversions, pairs, 100.0*inversions/Math.max(pairs, 1)));
		//TODO: Probable bug (ML046) - if every listed model disappears before loading, no best result remains to index.
		out.println("BEST by ECE: net_"+results.get(0).id+" ("+results.get(0).dims+")  ECE="+String.format("%.5f", results.get(0).ece));
	}

	/** One model's displayed identity and metrics; identity equality supports rank lookup. */
	static class Result{
		/** Network filename-derived ID and optional manifest dimension description. */
		final String id, dims;
		/** Mean squared error and population-weighted depth calibration error. */
		final double mse, ece;
		/** Stores one evaluated model; metric values are retained without normalization. */
		Result(String id, String dims, double mse, double ece){this.id=id; this.dims=dims; this.mse=mse; this.ece=ece;}
	}
}
