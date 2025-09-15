package jgi;

import fileIO.TextFile;
import shared.KillSwitch;
import shared.Tools;

public class RQCFilterStats {
	
	long readsIn;
	long basesIn;
	
	long readsOut;
	long basesOut;
	
	long readsDuplicate;
	long basesDuplicate;
	
	long readsLowQuality;
	long basesLowQuality;

	/** These are already counted under low quality */
	long readsPolyG;
	long basesPolyG;

	/** These are already counted under low quality */
	long readsN;
	long basesN;

	long readsArtifact;
	long basesArtifact;

	long readsFTrimmed;
	long basesFTrimmed;

	long readsAdapter;
	long basesAdapter;

	long readsSpikin;
	long basesSpikin;

	long readsRiboMap;
	long basesRiboMap;

	long readsChloroMap;
	long basesChloroMap;

	long readsMitoMap;
	long basesMitoMap;

	long readsSipMap;
	long basesSipMap;

	long readsRiboKmer;
	long basesRiboKmer;

	long readsMouse;
	long basesMouse;

	long readsCat;
	long basesCat;

	long readsDog;
	long basesDog;

	long readsHuman;
	long basesHuman;

	long readsMicrobe;
	long basesMicrobe;

	long readsOther;
	long basesOther;
	
	double gcPolymerRatio;
	
	long totalReadsRemoved(){
		return readsLowQuality/*+readsN*/+readsArtifact/*+readsFTrimmed*/
				+readsAdapter+readsSpikin+readsDuplicate
				+readsRiboMap+readsChloroMap+readsMitoMap+readsSipMap+readsRiboKmer
				+readsMouse+readsCat+readsDog+readsHuman+readsMicrobe+readsOther;
	}
	long totalBasesRemoved(){
		return basesLowQuality/*+basesN*/+basesArtifact+basesFTrimmed
				+basesAdapter+basesSpikin+basesDuplicate
				+basesRiboMap+basesChloroMap+basesMitoMap+basesSipMap+basesRiboKmer
				+basesMouse+basesCat+basesDog+basesHuman+basesMicrobe+basesOther;
	}
	
	@Override
	public String toString(){
		return toString(false);
	}
	
	public String toString(boolean skipAssertion){
		StringBuilder sb=new StringBuilder(1000);
		sb.append("#Class\tReads\tBases\tReadPct\tBasePct\tNotes\n");
		sb.append(format("Input", readsIn, basesIn, readsIn, basesIn));
		sb.append(format("Output", readsOut, basesOut, readsIn, basesIn));
		sb.append(format("Duplicate", readsDuplicate, basesDuplicate, readsIn, basesIn));
		sb.append(format("LowQuality", readsLowQuality, basesLowQuality, readsIn, basesIn));
		sb.append(format("PolyG", readsPolyG, basesPolyG, readsIn, basesIn, "\tSubsetOfLowQuality"));
		sb.append(format("N", readsN, basesN, readsIn, basesIn, "\tSubsetOfLowQuality"));
		sb.append(format("Artifact", readsArtifact, basesArtifact, readsIn, basesIn));
		sb.append(format("Spike-in", readsSpikin, basesSpikin, readsIn, basesIn));
		sb.append(format("ForceTrim", /*readsFTrimmed*/0, basesFTrimmed, readsIn, basesIn));
		sb.append(format("Adapter", readsAdapter, basesAdapter, readsIn, basesIn));
		sb.append(format("SipMap", readsSipMap, basesSipMap, readsIn, basesIn));
		sb.append(format("ChloroMap", readsChloroMap, basesChloroMap, readsIn, basesIn));
		sb.append(format("MitoMap", readsMitoMap, basesMitoMap, readsIn, basesIn));
		sb.append(format("RiboMap", readsRiboMap, basesRiboMap, readsIn, basesIn));
		sb.append(format("RiboKmer", readsRiboKmer, basesRiboKmer, readsIn, basesIn));
		sb.append(format("Human", readsHuman, basesHuman, readsIn, basesIn));
		sb.append(format("Mouse", readsMouse, basesMouse, readsIn, basesIn));
		sb.append(format("Cat", readsCat, basesCat, readsIn, basesIn));
		sb.append(format("Dog", readsDog, basesDog, readsIn, basesIn));
		sb.append(format("Microbe", readsMicrobe, basesMicrobe, readsIn, basesIn));
		sb.append(format("Other", readsOther, basesOther, readsIn, basesIn));
		
		assert(skipAssertion || readsIn>=readsOut) : toString(true)+"\n\nsb:\n"+sb+"\n";
		assert(skipAssertion || basesIn>=basesOut) : toString(true)+"\n\n"+sb+"\n";
		assert(skipAssertion || readsIn-totalReadsRemoved()==readsOut) : toString(true)+"\n\ntrr="+totalReadsRemoved()+"\nri="+readsIn+"\nro="+readsOut+"\n\nsb:\n"+sb+"\n";
		assert(skipAssertion || basesIn-totalBasesRemoved()==basesOut) : toString(true)+"\n"+basesIn+"-"+totalBasesRemoved()+"!="+basesOut+"\n\nsb:\n"+sb+"\n";
		
		return sb.toString();
	}
	
	StringBuilder format(String name, long reads, long bases, long rtot, long btot){
		return format(name, reads, bases, rtot, btot, null);
	}
	
	StringBuilder format(String name, long reads, long bases, long rtot, long btot, String suffix){
		assert(bases>=reads) : name+", "+reads+", "+bases+", "+rtot+", "+btot+", "+suffix;
		assert(btot>=rtot) : name+", "+reads+", "+bases+", "+rtot+", "+btot+", "+suffix;
		assert(reads<=rtot) : name+", "+reads+", "+bases+", "+rtot+", "+btot+", "+suffix;
		assert(bases<=btot) : name+", "+reads+", "+bases+", "+rtot+", "+btot+", "+suffix;
		StringBuilder sb=new StringBuilder();
		sb.append(name).append('\t');
		sb.append(reads).append('\t');
		sb.append(bases).append('\t');
		sb.append(toPercent(reads, rtot, 3)).append('\t');
		sb.append(toPercent(bases, btot, 3));
		if(suffix!=null){sb.append(suffix);}
		sb.append('\n');
		return sb;
	}
	
	private static String toPercent(long numerator, long denominator, int decimals){
		if(denominator<1){denominator=1;}
		return Tools.format("%."+decimals+"f",numerator*100.0/denominator);
	}
	
	void parseHuman(String fname){
		if(fname==null){return;}
		TextFile tf=new TextFile(fname);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(!line.startsWith("#")){
				long[] ret=parseStatsLine(line);
				if(line.startsWith("human")){
					readsHuman+=ret[0];
					basesHuman+=ret[1];
				}else if(line.startsWith("cat")){
					readsCat+=ret[0];
					basesCat+=ret[1];
				}else if(line.startsWith("dog")){
					readsDog+=ret[0];
					basesDog+=ret[1];
				}else if(line.startsWith("mouse")){
					readsMouse+=ret[0];
					basesMouse+=ret[1];
				}else{
					assert(false) : line;
				}
			}
		}
		tf.close();
	}
	
//	long[] parseStatsLine(String line){
//		String[] split=line.split("\t");
//		long[] ret=KillSwitch.allocLong1D(2);
//		ret[0]=Long.parseLong(split[5])+Long.parseLong(split[6])/2;
//		ret[1]=(long)(1000000L*(Double.parseDouble(split[2])+Long.parseLong(split[4])/2));
//		return ret;
//	}
	
	long[] parseStatsLine(String line){
		String[] split=line.split("\t");
		long[] ret=KillSwitch.allocLong1D(2);
		ret[0]=Long.parseLong(split[7]);
		ret[1]=Long.parseLong(split[8]);
		return ret;
	}
	
	void parseChloro(String fname){
		if(fname==null){return;}
		TextFile tf=new TextFile(fname);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(!line.startsWith("#")){
				long[] ret=parseStatsLine(line);
				if(line.contains("lcl|SSU_") || line.contains("lcl|LSU_")){
					readsRiboMap+=ret[0];
					basesRiboMap+=ret[1];
				}else if(line.contains("mitochondrion")){
					readsMitoMap+=ret[0];
					basesMitoMap+=ret[1];
				}else if(line.startsWith("plastid") || line.startsWith("chloroplast")){
					readsChloroMap+=ret[0];
					basesChloroMap+=ret[1];
				}else{
					readsOther+=ret[0];
					basesOther+=ret[1];
				}
			}
		}
		tf.close();
	}
	
	void parseSip(String fname){
		if(fname==null){return;}
		TextFile tf=new TextFile(fname);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(!line.startsWith("#")){
				long[] ret=parseStatsLine(line);
				if(true){
					readsSipMap+=ret[0];
					basesSipMap+=ret[1];
				}else{
					readsOther+=ret[0];
					basesOther+=ret[1];
				}
			}
		}
		tf.close();
	}
	
//	//Ribo:
//	tid|39138|lcl|SSU_GQ398331.1.1377 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter;bacterium A1(2009)
//	//Chloro:
//	tid|48534|ref|NC_016471.1| Neottia nidus-avis plastid, complete genome
//	tid|13331|ref|NC_033913.1| Akebia quinata chloroplast, complete genome
//	//Mito:
//	tid|935657|ref|NC_026218.1| Colletes gigas mitochondrion, complete genome
	
	/*
	Input Reads: 53768366
	Input Bases: 8119023266

	Reads Removed for:
	Low Quality: xxx
	N's: xxx
	Too Short after trimming: (minlen param)
	Artifact/synthetic contamination: xxx
	Adapter:
	Spike-in
	Ribosomal RNA (also Mito, Chloro)
	Microbial contamination:
	Cat DNA
	Dog DNA
	Mouse DNA
	Human DNA
	Total reads removed:

	Bases removed from reads by trimming:
	low quality on end: xxx
	adapters: xxx
	Total bases removed: xxx

	Remaining Reads: xxx (yy%)
	Remaining Bases: xxxx (yy%)
	*/
	
}
