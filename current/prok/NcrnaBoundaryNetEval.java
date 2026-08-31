package prok;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import ml.CellNet;

/**
 * Evaluates an ncRNA boundary network on genome-split vector rows.  Vectors are
 * emitted as one START candidate block followed by one STOP candidate block per
 * locus; the label-1 row is the exact endpoint.  This reports the metric used by
 * production selection (argmax within each boundary block), not row-wise error.
 */
public class NcrnaBoundaryNetEval {

	public static void main(String[] args) throws IOException{
		String netPath=null, inPath=null;
		for(String arg : args){
			String[] split=arg.split("=", 2);
			if(split.length<2){continue;}
			if(split[0].equalsIgnoreCase("net")){netPath=split[1];}
			else if(split[0].equalsIgnoreCase("in")){inPath=split[1];}
		}
		if(netPath==null || inPath==null){
			throw new IllegalArgumentException("Usage: net=<model.bbnet> in=<vectors.tsv>");
		}
		CellNet net=NcrnaBoundaryScorer.load(netPath);
		Result result=evaluate(net, inPath);
		System.out.println("loci\tstartExact\tstopExact\tbothExact\tstartPct\tstopPct\tbothPct");
		System.out.printf("%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f%n",
			result.loci, result.startExact, result.stopExact, result.bothExact,
			pct(result.startExact, result.loci), pct(result.stopExact, result.loci),
			pct(result.bothExact, result.loci));
	}

	static Result evaluate(CellNet net, String path) throws IOException{
		Result result=new Result();
		BoundaryBlock start=new BoundaryBlock(false);
		BoundaryBlock stop=new BoundaryBlock(true);
		boolean sawStop=false;
		int dims=-1;
		try(BufferedReader reader=new BufferedReader(new FileReader(path))){
			String line;
			while((line=reader.readLine())!=null){
				if(line.isEmpty()){continue;}
				if(line.charAt(0)=='#'){
					String[] fields=line.split("\t");
					if(fields.length>=3 && fields[0].equals("#dims")){
						dims=Integer.parseInt(fields[1]);
						if(dims!=NcrnaBoundaryScorer.NUM_FEATURES || Integer.parseInt(fields[2])!=1){
							throw new IllegalArgumentException("Expected #dims\t"+NcrnaBoundaryScorer.NUM_FEATURES+"\t1: "+line);
						}
					}
					continue;
				}
				if(dims<0){throw new IllegalArgumentException("Missing #dims header in "+path);}
				String[] fields=line.split("\t");
				if(fields.length!=dims+1){throw new IllegalArgumentException("Expected "+(dims+1)+" columns: "+line);}
				float[] input=new float[dims];
				for(int i=0; i<dims; i++){input[i]=Float.parseFloat(fields[i]);}
				int label=Integer.parseInt(fields[dims]);
				if(label!=0 && label!=1){throw new IllegalArgumentException("Boundary label must be 0/1: "+label);}
				boolean isStop=input[4]>=0.5f;
				if(!isStop && sawStop){finishLocus(start, stop, result); sawStop=false;}
				if(isStop){sawStop=true;}
				BoundaryBlock block=(isStop ? stop : start);
				net.applyInput(input);
				block.accept(net.feedForward(), label);
			}
		}
		if(start.rows>0 || stop.rows>0){finishLocus(start, stop, result);}
		return result;
	}

	private static void finishLocus(BoundaryBlock start, BoundaryBlock stop, Result result){
		start.validate();
		stop.validate();
		boolean startOk=start.bestLabel==1;
		boolean stopOk=stop.bestLabel==1;
		result.loci++;
		if(startOk){result.startExact++;}
		if(stopOk){result.stopExact++;}
		if(startOk && stopOk){result.bothExact++;}
		start.clear();
		stop.clear();
	}

	private static double pct(long numerator, long denominator){
		return denominator==0 ? 0 : 100.0*numerator/denominator;
	}

	static class BoundaryBlock {
		BoundaryBlock(boolean stop_){stop=stop_;}
		void accept(float score, int label){
			if(label==1){positiveRows++;}
			// Production's candidate sweep also keeps the first candidate on an exact tie.
			if(rows==0 || score>bestScore){bestScore=score; bestLabel=label;}
			rows++;
		}
		void validate(){
			if(rows<1 || positiveRows!=1){
				throw new IllegalArgumentException((stop ? "STOP" : "START")+" block has rows="+rows
					+" label-1 rows="+positiveRows);
			}
		}
		void clear(){rows=positiveRows=bestLabel=0; bestScore=-Float.MAX_VALUE;}
		final boolean stop;
		int rows=0, positiveRows=0, bestLabel=0;
		float bestScore=-Float.MAX_VALUE;
	}

	static class Result {
		long loci, startExact, stopExact, bothExact;
	}
}
