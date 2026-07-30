package ml;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/**
 * Enriches a ranking vector file by duplicating queries where the net's
 * top-1 pick is not the best available hit. Outputs a new vector file
 * with the hard cases duplicated N times.
 *
 * Usage: java ml.EnrichFailures vectors=<file> net=<bbnet> out=<file>
 *          [copies=3] [threshold=0.1]
 *
 * A query is a "failure" when the best available LCA depth minus the
 * net's top-1 pick's LCA depth exceeds threshold. The query's rows
 * are written (copies+1) times: once normally plus N duplicates.
 *
 * @author Neptune
 * @date July 2026
 */
public class EnrichFailures{

	public static void main(String[] args) throws Exception{
		String vecFile=null, netFile=null, outFile=null;
		int copies=3;
		float threshold=0.1f;
		for(String arg : args){
			final String[] kv=arg.split("=", 2);
			final String k=kv[0].toLowerCase();
			final String v=(kv.length>1) ? kv[1] : "";
			if(k.equals("vectors") || k.equals("in")){vecFile=v;}
			else if(k.equals("net")){netFile=v;}
			else if(k.equals("out")){outFile=v;}
			else if(k.equals("copies")){copies=Integer.parseInt(v);}
			else if(k.equals("threshold")){threshold=Float.parseFloat(v);}
		}
		if(vecFile==null || netFile==null || outFile==null){
			throw new IllegalArgumentException("required: vectors= net= out=");
		}

		final CellNet net=CellNetParser.load(netFile);

		// Load all rows, detect query boundaries, score with net
		final ArrayList<String> lines=new ArrayList<>();
		final ArrayList<Float> labels=new ArrayList<>();
		final ArrayList<Double> scores=new ArrayList<>();
		final ArrayList<Integer> qstarts=new ArrayList<>();
		String header=null;

		int NI=33;
		final ByteFile bf=ByteFile.makeByteFile(vecFile, true);
		byte[] line;
		int rowInQuery=0;
		float prevDim29=-1;
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
			final float dim29=x[28];

			if(lines.isEmpty() || rowInQuery>=10 || (dim29!=prevDim29 && !lines.isEmpty())){
				qstarts.add(lines.size());
				rowInQuery=0;
			}

			net.applyInput(x);
			final double score=net.feedForward();

			lines.add(s);
			labels.add(label);
			scores.add(score);
			prevDim29=dim29;
			rowInQuery++;
		}
		bf.close();

		final int n=lines.size();
		final int nq=qstarts.size();

		// Find failure queries
		int failures=0;
		final boolean[] isFailure=new boolean[nq];
		for(int qi=0; qi<nq; qi++){
			final int s=qstarts.get(qi);
			final int e=(qi+1<nq ? qstarts.get(qi+1) : n);

			// Find net's top pick and best available
			int netPick=s;
			float bestLabel=labels.get(s);
			for(int i=s; i<e; i++){
				if(scores.get(i)>scores.get(netPick)){netPick=i;}
				if(labels.get(i)>bestLabel){bestLabel=labels.get(i);}
			}
			final float gap=bestLabel-labels.get(netPick);
			if(gap>threshold){
				isFailure[qi]=true;
				failures++;
			}
		}

		// Write output: all rows, with failure queries duplicated
		final ByteStreamWriter bsw=new ByteStreamWriter(outFile, true, false, true);
		bsw.start();
		if(header!=null){bsw.println(header);}

		int normalRows=0, enrichedRows=0;
		for(int qi=0; qi<nq; qi++){
			final int s=qstarts.get(qi);
			final int e=(qi+1<nq ? qstarts.get(qi+1) : n);
			final int reps=(isFailure[qi] ? copies+1 : 1);
			for(int r=0; r<reps; r++){
				for(int i=s; i<e; i++){
					bsw.println(lines.get(i));
					if(r==0){normalRows++;}
					else{enrichedRows++;}
				}
			}
		}
		bsw.poisonAndWait();

		System.err.println("Queries: "+nq+"  failures: "+failures
			+" ("+String.format("%.1f", 100.0*failures/nq)+"%)"
			+"  threshold: "+threshold);
		System.err.println("Output rows: "+(normalRows+enrichedRows)
			+" ("+normalRows+" normal + "+enrichedRows+" enriched, "
			+copies+"x copies of failures)");
	}
}
