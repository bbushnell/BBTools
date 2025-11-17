package jgi;

import fileIO.TextFile;
import shared.KillSwitch;
import shared.Tools;

/**
 * Tracks and formats statistics for RQC (Read Quality Control) filtering operations.
 * Accumulates counts of reads and bases removed by various quality control steps
 * including quality trimming, contamination removal, and deduplication.
 * Provides formatted output for reporting filtering results.
 *
 * @author Brian Bushnell
 */
public class RQCFilterStats {
	
	/** Total number of input reads */
	long readsIn;
	/** Total number of input bases */
	long basesIn;
	
	/** Number of reads remaining after filtering */
	long readsOut;
	/** Number of bases remaining after filtering */
	long basesOut;
	
	/** Number of duplicate reads removed */
	long readsDuplicate;
	/** Number of bases from duplicate reads removed */
	long basesDuplicate;
	
	/** Number of low quality reads removed */
	long readsLowQuality;
	/** Number of bases from low quality reads removed */
	long basesLowQuality;

	/** These are already counted under low quality */
	long readsPolyG;
	/** Number of bases from PolyG reads removed (subset of low quality) */
	long basesPolyG;

	/** These are already counted under low quality */
	long readsN;
	/** Number of bases from N-containing reads removed (subset of low quality) */
	long basesN;

	/** Number of reads removed as artifacts or synthetic contamination */
	long readsArtifact;
	/** Number of bases from artifact reads removed */
	long basesArtifact;

	/** Number of reads that underwent force trimming */
	long readsFTrimmed;
	/** Number of bases removed by force trimming */
	long basesFTrimmed;

	/** Number of reads removed for adapter contamination */
	long readsAdapter;
	/** Number of bases from adapter-contaminated reads removed */
	long basesAdapter;

	/** Number of reads removed as spike-in sequences */
	long readsSpikin;
	/** Number of bases from spike-in reads removed */
	long basesSpikin;

	/** Number of reads removed by ribosomal RNA mapping */
	long readsRiboMap;
	/** Number of bases from ribosomal RNA reads removed */
	long basesRiboMap;

	/** Number of reads removed by chloroplast mapping */
	long readsChloroMap;
	/** Number of bases from chloroplast reads removed */
	long basesChloroMap;

	/** Number of reads removed by mitochondrial mapping */
	long readsMitoMap;
	/** Number of bases from mitochondrial reads removed */
	long basesMitoMap;

	/** Number of reads removed by SIP (Stable Isotope Probing) mapping */
	long readsSipMap;
	/** Number of bases from SIP reads removed */
	long basesSipMap;

	/** Number of reads removed by ribosomal k-mer filtering */
	long readsRiboKmer;
	/** Number of bases from ribosomal k-mer filtered reads removed */
	long basesRiboKmer;

	/** Number of reads removed as mouse DNA contamination */
	long readsMouse;
	/** Number of bases from mouse DNA reads removed */
	long basesMouse;

	/** Number of reads removed as cat DNA contamination */
	long readsCat;
	/** Number of bases from cat DNA reads removed */
	long basesCat;

	/** Number of reads removed as dog DNA contamination */
	long readsDog;
	/** Number of bases from dog DNA reads removed */
	long basesDog;

	/** Number of reads removed as human DNA contamination */
	long readsHuman;
	/** Number of bases from human DNA reads removed */
	long basesHuman;

	/** Number of reads removed as microbial contamination */
	long readsMicrobe;
	/** Number of bases from microbial contamination reads removed */
	long basesMicrobe;

	/** Number of reads removed as other/unclassified contamination */
	long readsOther;
	/** Number of bases from other contamination reads removed */
	long basesOther;
	
	/** Ratio related to GC content and polymer detection */
	double gcPolymerRatio;
	
	/**
	 * Calculates the total number of reads removed across all filtering categories.
	 * Excludes readsFTrimmed and readsN from the total as noted in comments.
	 * @return Total count of reads removed by all filtering steps
	 */
	long totalReadsRemoved(){
		return readsLowQuality/*+readsN*/+readsArtifact/*+readsFTrimmed*/
				+readsAdapter+readsSpikin+readsDuplicate
				+readsRiboMap+readsChloroMap+readsMitoMap+readsSipMap+readsRiboKmer
				+readsMouse+readsCat+readsDog+readsHuman+readsMicrobe+readsOther;
	}
	/**
	 * Calculates the total number of bases removed across all filtering categories.
	 * Includes basesFTrimmed but excludes basesN from the total as noted in comments.
	 * @return Total count of bases removed by all filtering steps
	 */
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
	
	/**
	 * Generates a comprehensive formatted report of RQC filtering statistics.
	 * Creates tab-separated table with read counts, base counts, and percentages
	 * for each filtering category. Includes validation assertions unless skipped.
	 *
	 * @param skipAssertion If true, skips validation assertions for debugging
	 * @return Formatted multi-line report with headers and statistics
	 */
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
	
	/**
	 * Formats a single statistics line for the report.
	 *
	 * @param name Category name
	 * @param reads Read count for this category
	 * @param bases Base count for this category
	 * @param rtot Total reads for percentage calculation
	 * @param btot Total bases for percentage calculation
	 * @return Formatted StringBuilder with tab-separated values
	 */
	StringBuilder format(String name, long reads, long bases, long rtot, long btot){
		return format(name, reads, bases, rtot, btot, null);
	}
	
	/**
	 * Formats a single statistics line with optional suffix annotation.
	 * Creates tab-separated line with counts and percentages, validates inputs.
	 *
	 * @param name Category name
	 * @param reads Read count for this category
	 * @param bases Base count for this category
	 * @param rtot Total reads for percentage calculation
	 * @param btot Total bases for percentage calculation
	 * @param suffix Optional annotation text (e.g., "SubsetOfLowQuality")
	 * @return Formatted StringBuilder with tab-separated values and newline
	 */
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
	
	/**
	 * Converts fraction to formatted percentage string.
	 *
	 * @param numerator Numerator value
	 * @param denominator Denominator value (minimum 1 to avoid division by zero)
	 * @param decimals Number of decimal places in output
	 * @return Formatted percentage string
	 */
	private static String toPercent(long numerator, long denominator, int decimals){
		if(denominator<1){denominator=1;}
		return Tools.format("%."+decimals+"f",numerator*100.0/denominator);
	}
	
	/**
	 * Parses contamination statistics from human/animal screening results file.
	 * Reads tab-separated file and accumulates counts for human, cat, dog, and mouse
	 * DNA contamination based on line prefixes.
	 * @param fname Path to contamination screening results file
	 */
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
	
	/**
	 * Parses read and base counts from a tab-separated statistics line.
	 * Extracts values from columns 7 and 8 (0-based indexing).
	 * @param line Tab-separated input line
	 * @return Array containing [read_count, base_count]
	 */
	long[] parseStatsLine(String line){
		String[] split=line.split("\t");
		long[] ret=KillSwitch.allocLong1D(2);
		ret[0]=Long.parseLong(split[7]);
		ret[1]=Long.parseLong(split[8]);
		return ret;
	}
	
	/**
	 * Parses organellar contamination statistics from mapping results file.
	 * Categorizes reads as ribosomal (SSU/LSU), mitochondrial, chloroplast,
	 * or other based on reference sequence identifiers in the mapping results.
	 * @param fname Path to organellar mapping results file
	 */
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
	
	/**
	 * Parses SIP (Stable Isotope Probing) contamination statistics from file.
	 * Currently adds all reads to SIP category unconditionally due to
	 * condition always being true.
	 * @param fname Path to SIP mapping results file
	 */
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
