package prok;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import consensus.BaseGraph;
import fileIO.FileFormat;
import parse.Parse;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;

/** Stop-boundary vectors centered on actual accepted caller spans and winning models. */
public class NcrnaCallerStateVectorGen {
	public static void main(String[] args) throws IOException{
		String fasta=null,outPath=null,stopTablePath=null;
		String r58Kmers=null,r58Consensus=null,r58Models=null,lsuKmers=null,lsuConsensus=null,lsuModels=null;
		boolean leaveOneOut=false;
		for(String a:args){String[] kv=a.split("=",2); if(kv.length<2){continue;}
			if(kv[0].equalsIgnoreCase("fasta")){fasta=kv[1];}
			else if(kv[0].equalsIgnoreCase("out")){outPath=kv[1];}
			else if(kv[0].equalsIgnoreCase("tablestop")){stopTablePath=kv[1];}
			else if(kv[0].equalsIgnoreCase("leaveoneout")||kv[0].equalsIgnoreCase("loo")){leaveOneOut=Parse.parseBoolean(kv[1]);}
			else if(kv[0].equalsIgnoreCase("r58kmers")){r58Kmers=kv[1];}
			else if(kv[0].equalsIgnoreCase("r58consensus")){r58Consensus=kv[1];}
			else if(kv[0].equalsIgnoreCase("r58models")){r58Models=kv[1];}
			else if(kv[0].equalsIgnoreCase("lsukmers")){lsuKmers=kv[1];}
			else if(kv[0].equalsIgnoreCase("lsuconsensus")){lsuConsensus=kv[1];}
			else if(kv[0].equalsIgnoreCase("lsumodels")){lsuModels=kv[1];}
		}
		if(fasta==null||outPath==null||stopTablePath==null){throw new IllegalArgumentException("Required: fasta= out= tablestop=");}
		CallGenes.setR58LsuEnabled(true); CallGenes.NCRNA_BOUNDARY_NN_ENABLED=false;
		CallGenes.R58_KMERS_OVERRIDE=r58Kmers; CallGenes.R58_CONSENSUS_OVERRIDE=r58Consensus; CallGenes.R58_MODELS_OVERRIDE=r58Models;
		CallGenes.LSU_KMERS_OVERRIDE=lsuKmers; CallGenes.LSU_CONSENSUS_OVERRIDE=lsuConsensus; CallGenes.LSU_MODELS_OVERRIDE=lsuModels;
		CallGenes.loadNcrnaResources(); NcrnaFamily fam=null;
		for(NcrnaFamily f:GeneCaller.ncrnaFamilies){if(f.name.equals("lsu")){fam=f;break;}}
		if(fam==null||fam.models==null){throw new IllegalArgumentException("LSU resources unavailable");}
		final TrnaNinemerTableBuilder.LoadedTable lt=TrnaNinemerTableBuilder.loadTable(stopTablePath);
		if(lt.type!=TrnaBoundaryFeatures.BoundaryType.STOP){throw new IllegalArgumentException("Not a STOP table");}
		final float meanLen=medianLength(fam.library);
		Shared.TRIM_READ_DESCRIPTION=false; Shared.TRIM_RNAME=true; Read.TO_UPPER_CASE=true;
		final ArrayList<Read> reads=ReadInputStream.toReads(fasta,FileFormat.FA,-1);
		long written=0;
		try(PrintStream out=new PrintStream(new FileOutputStream(outPath))){
			out.println("#dims\t10\t1");
			for(Read r:reads){
				final int s=parseInt(r.id,"call_start="), e0=parseInt(r.id,"call_stop="), truth=parseInt(r.id,"truth_stop=");
				final String modelName=parseString(r.id,"model="); final float gc=parseFloat(r.id,"contig_gc=");
				int mi=-1; for(int i=0;i<fam.modelNames.length;i++){if(firstToken(fam.modelNames[i]).equals(firstToken(modelName))){mi=i;break;}}
				if(mi<0){throw new IllegalArgumentException("Unknown model in "+r.id);}
				int positives=0,rows=0;
				for(int off:fam.boundaryStopOffsets){
					final int e=e0+off; if(s<0||e>=r.length()||e-s<15){throw new IllegalArgumentException("Candidate out of bounds: "+r.id);}
					final byte[] cand=java.util.Arrays.copyOfRange(r.bases,s,e+1);
					final float ani=TrnaBoundaryFeatures.aniFeature(cand,fam.library[mi]);
					final float[] prof=(leaveOneOut
						? TrnaBoundaryFeatures.enrichmentProfile(r.bases,e,truth,TrnaBoundaryFeatures.BoundaryType.STOP,lt.insideCount,lt.outsideCount,lt.table)
						: TrnaBoundaryFeatures.enrichmentProfile(r.bases,e,TrnaBoundaryFeatures.BoundaryType.STOP,lt.insideCount,lt.outsideCount,lt.table));
					final BaseGraph model=fam.models[mi]; final float[] fuzz=TrnaBoundaryFeatures.tipFuzzinessFeature(cand,model,false);
					final float lengthRatio=(e-s+1)/meanLen; final int label=(e==truth?1:0); positives+=label; rows++;
					out.printf("%.6f\t%.6f\t%.6f\t%.6f\t1.000000\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d%n",
						ani,prof[0],prof[1],prof[2],fuzz[0],fuzz[1],fuzz[2],lengthRatio,gc,label);
				}
				if(rows!=9||positives!=1){throw new IllegalArgumentException("Bad capture block rows="+rows+" positives="+positives+": "+r.id);}
				written+=rows;
			}
		}
		System.err.println("NcrnaCallerStateVectorGen wrote="+written+" loci="+reads.size()+" leaveOneOut="+leaveOneOut);
	}

	private static int parseInt(String s,String k){return Integer.parseInt(parseString(s,k));}
	private static float parseFloat(String s,String k){return Float.parseFloat(parseString(s,k));}
	private static String parseString(String s,String k){int a=s.indexOf(k);if(a<0){throw new IllegalArgumentException("Missing "+k+" in "+s);}a+=k.length();int b=s.indexOf(';',a);return b<0?s.substring(a):s.substring(a,b);}
	private static String firstToken(String s){int x=s.indexOf(' ');return x<0?s:s.substring(0,x);}
	private static float medianLength(byte[][] lib){int[] x=new int[lib.length];for(int i=0;i<x.length;i++){x[i]=lib[i].length;}java.util.Arrays.sort(x);return x[x.length/2];}
}
