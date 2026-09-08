package prok;

import java.io.PrintWriter;
import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.LongHashSet;
import parse.Parse;
import shared.Shared;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;

/** Development-only direct runner for a seed-and-align rRNA experiment.
 * Deliberately bypasses CallGenes, so it is suitable for isolated ncRNA
 * calibration but must never be construed as a production calling path. */
public class RrnaScavengerHarness {

	public static void main(String[] args) throws Exception {
		String in=null, out=null, consensus=null, kmers=null, hbmPath=null, workloadOut=null, family="5S";
		int k=15, minLen=90, windowPad=100, indexTopN=12, quantumThresh=120, fixedMinHits=12;
		boolean adaptive=true, normalizeCase=false;
		float adaptFloor=11f, adaptTopFrac=.48f, adaptQFrac=.072f;
		float idPass=.60f, idBorderline=.60f, hbmPass=Float.NaN, collapseFrac=.85f;
		for(String arg : args) {
			final int eq=arg.indexOf('=');
			if(eq<1) {throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0, eq).toLowerCase(), value=arg.substring(eq+1);
			if(key.equals("in")) {in=value;} else if(key.equals("out")) {out=value;}
			else if(key.equals("consensus")) {consensus=value;} else if(key.equals("kmers")) {kmers=value;}
			else if(key.equals("hbm")) {hbmPath=value;} else if(key.equals("hbmpass")) {hbmPass=Float.parseFloat(value);}
			else if(key.equals("workloadout")) {workloadOut=value;} else if(key.equals("family")) {family=value;}
			else if(key.equals("k")) {k=Integer.parseInt(value);} else if(key.equals("minlen")) {minLen=Integer.parseInt(value);}
			else if(key.equals("windowpad")) {windowPad=Integer.parseInt(value);} else if(key.equals("indextopn")) {indexTopN=Integer.parseInt(value);}
			else if(key.equals("quantumthresh")) {quantumThresh=Integer.parseInt(value);} else if(key.equals("adaptive")) {adaptive=Parse.parseBoolean(value);}
			else if(key.equals("adaptfloor")) {adaptFloor=Float.parseFloat(value);} else if(key.equals("adapttopfrac")) {adaptTopFrac=Float.parseFloat(value);}
			else if(key.equals("adaptqfrac")) {adaptQFrac=Float.parseFloat(value);} else if(key.equals("fixedminhits")) {fixedMinHits=Integer.parseInt(value);}
			else if(key.equals("normalizecase")) {normalizeCase=Parse.parseBoolean(value);} else if(key.equals("idpass")) {idPass=Float.parseFloat(value);}
			else if(key.equals("idborderline")) {idBorderline=Float.parseFloat(value);} else if(key.equals("collapsefrac")) {collapseFrac=Float.parseFloat(value);}
			else {throw new IllegalArgumentException("Unknown parameter: "+key);}
		}
		if(in==null || out==null || consensus==null || kmers==null) {throw new IllegalArgumentException("Required: in= out= consensus= kmers=");}
		if(k<1 || k>31 || minLen<1 || windowPad<0 || indexTopN<1 || quantumThresh<1 || collapseFrac<0 || collapseFrac>1 || idPass<0 || idPass>1 || idBorderline<0 || idBorderline>idPass) {throw new IllegalArgumentException("Invalid scavenger parameters");}
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(consensus);
		final String[] names=TrnaConsensusBuilder.lastLoadedNames;
		if(library==null || library.length==0) {throw new IllegalArgumentException("No consensus records: "+consensus);}
		final consensus.BaseGraph[] models=(hbmPath==null ? null : TrnaConsensusBuilder.loadModels(hbmPath));
		if(models!=null) {CallGenes.requireNcrnaLibraryModelAlignment(family, library, models, names, consensus, hbmPath); if(Float.isNaN(hbmPass)){hbmPass=.75f;}}
		else if(!Float.isNaN(hbmPass)) {throw new IllegalArgumentException("hbmpass requires hbm=");}
		final LongHashSet seedSet=ProkObject.loadLongKmers(kmers, k);
		if(seedSet==null || seedSet.size()<1) {throw new IllegalArgumentException("No length-"+k+" seeds: "+kmers);}
		final NcrnaScavenger scavenger=new NcrnaScavenger(library, models, names, seedSet, k, minLen, windowPad, 7, indexTopN, adaptive, adaptFloor, adaptTopFrac, adaptQFrac, fixedMinHits, 0f, 1f, idPass, idBorderline);
		scavenger.hbmPass=(models==null ? 1.01f : hbmPass); scavenger.collapseFrac=collapseFrac; scavenger.quantumThresh=quantumThresh;
		long reads=0, calls=0;
		final PrintWriter workloadPw=(workloadOut==null ? null : new PrintWriter(workloadOut));
		if(workloadPw!=null){
			workloadPw.println("event\tcontig\tstrand\tpass\tstart\tstop\tseed_hits");
			scavenger.setWorkloadSink(new NcrnaWorkloadInstrumentSink() {
				@Override public void seedHits(String contig, int strand, int[] hitPositions) {workloadPw.println("seed\t"+contig+"\t"+strand+"\t.\t.\t.\t"+hitPositions.length);}
				@Override public void scheduledWindow(String contig, int strand, int pass, int start, int stop) {workloadPw.println("window\t"+contig+"\t"+strand+"\t"+pass+"\t"+start+"\t"+stop+"\t.");}
			});
		}
		final FileFormat ffin=FileFormat.testInput(in, FileFormat.FA, null, true, true);
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ffin, null);
		try(PrintWriter pw=new PrintWriter(out)) {
			pw.println("##gff-version 3"); cris.start();
			for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()) {
				for(Read read : ln) {
					if(read.bases==null || read.length()<minLen) {continue;} reads++;
					final byte[] bases=(normalizeCase ? read.bases.clone() : read.bases); if(normalizeCase){Tools.toUpperCase(bases);}
					calls+=writeCalls(pw, scavenger.scavenge(read.id, bases, Shared.PLUS, new ArrayList<int[]>()), family, false);
					calls+=writeCalls(pw, scavenger.scavenge(read.id, AminoAcid.reverseComplementBases(bases), Shared.MINUS, new ArrayList<int[]>()), family, true);
				}
				cris.returnList(ln);
			}
		} finally {ReadWrite.closeStream(cris); if(workloadPw!=null){workloadPw.close();}}
		System.err.println("RrnaScavengerHarness: reads="+reads+" calls="+calls+" consensus="+library.length+" seeds="+seedSet.size());
	}

	private static int writeCalls(PrintWriter pw, ArrayList<Orf> calls, String family, boolean flip) {
		if(calls==null) {return 0;}
		for(Orf orf : calls) {if(flip) {orf.flip();} orf.trnaModel=family; pw.println(orf.toGff());}
		return calls.size();
	}
}
