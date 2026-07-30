package ml;

import java.util.ArrayList;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;

/**
 * Subsamples a training vector file by prediction error: keep 100% of the
 * highest-error rows and taper down to a floor fraction for the lowest-error rows.
 * Produces a training set enriched for hard cases without duplicating any row.
 *
 * Usage: java ml.ErrorTaper vectors=<file> net=<bbnet> out=<file>
 *          [floor=0.1] [seed=42]
 *
 * The taper is linear: row i (sorted by descending error) is kept with probability
 * 1.0 - (1.0 - floor) * i / n. The hardest row is always kept; the easiest is
 * kept with probability floor.
 *
 * @author Neptune
 * @date July 2026
 */
public class ErrorTaper{

	public static void main(String[] args) throws Exception{
		String vecFile=null, netFile=null, outFile=null;
		float floor=0.1f;
		long seed=42;
		for(String arg : args){
			final String[] kv=arg.split("=", 2);
			final String k=kv[0].toLowerCase();
			final String v=(kv.length>1) ? kv[1] : "";
			if(k.equals("vectors") || k.equals("in")){vecFile=v;}
			else if(k.equals("net")){netFile=v;}
			else if(k.equals("out")){outFile=v;}
			else if(k.equals("floor")){floor=Float.parseFloat(v);}
			else if(k.equals("seed")){seed=Long.parseLong(v);}
		}
		if(vecFile==null || netFile==null || outFile==null){
			throw new IllegalArgumentException("required: vectors= net= out=");
		}

		final CellNet net=CellNetParser.load(netFile);
		int NI=33;

		// Load all rows and score
		final ArrayList<String> lines=new ArrayList<>();
		final ArrayList<Float> errors=new ArrayList<>();
		String header=null;

		final ByteFile bf=ByteFile.makeByteFile(vecFile, true);
		byte[] line;
		while((line=bf.nextLine())!=null){
			final String s=new String(line);
			if(s.isEmpty()){continue;}
			if(s.charAt(0)=='#'){
				header=s;
				if(s.startsWith("#dims")){
					final String[] hf=s.split("\t");
					if(hf.length>=2){NI=Integer.parseInt(hf[1].trim());}
				}
				continue;
			}
			final String[] f=s.split("\t");
			if(f.length<NI+1){continue;}
			final float[] x=new float[NI];
			for(int i=0; i<NI; i++){x[i]=Float.parseFloat(f[i]);}
			final float label=Float.parseFloat(f[NI]);

			net.applyInput(x);
			final float pred=(float)net.feedForward();
			final float err=Math.abs(pred-label);

			lines.add(s);
			errors.add(err);
		}
		bf.close();

		final int n=lines.size();
		System.err.println("Loaded "+n+" rows, NI="+NI);

		// Sort indices by error descending
		final Integer[] idx=new Integer[n];
		for(int i=0; i<n; i++){idx[i]=i;}
		java.util.Arrays.sort(idx, (a,b)->Float.compare(errors.get(b), errors.get(a)));

		// Taper: row at rank i kept with prob 1.0 - (1.0-floor)*i/n
		final Random rnd=new Random(seed);
		final ByteStreamWriter bsw=new ByteStreamWriter(outFile, true, false, true);
		bsw.start();
		if(header!=null){bsw.println(header);}

		int kept=0;
		for(int rank=0; rank<n; rank++){
			final float prob=1.0f-(1.0f-floor)*rank/(float)(n-1);
			if(rnd.nextFloat()<prob){
				bsw.println(lines.get(idx[rank]));
				kept++;
			}
		}
		bsw.poisonAndWait();

		final float topErr=errors.get(idx[0]), botErr=errors.get(idx[n-1]);
		final float medErr=errors.get(idx[n/2]);
		System.err.println(String.format("Kept %d of %d (%.1f%%). floor=%.2f. Error range: %.4f - %.4f (median %.4f)",
			kept, n, 100.0*kept/n, floor, botErr, topErr, medErr));
	}
}
