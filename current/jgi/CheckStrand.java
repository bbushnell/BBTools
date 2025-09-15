package jgi;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import dna.Data;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import gff.GffLine;
import prok.CallGenes;
import prok.GeneCaller;
import prok.GeneModel;
import prok.GeneModelParser;
import prok.Orf;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.Sketch;
import sketch.SketchMakerMini;
import sketch.SketchObject;
import sketch.SketchTool;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.Read;
import stream.ReadInputStream;
import structures.DoubleList;
import structures.Feature;
import structures.ListNum;

/**
 * Checks the strandedness of RNA-seq reads.
 *  
 * @author Brian Bushnell
 * @date Aug 4, 2023
 *
 */
public class CheckStrand {

	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		CheckStrand x=new CheckStrand(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CheckStrand(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		parser.out1=out1;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("size") || a.equals("len") || a.equals("length") || a.equals("sketchsize")){
				sketchSize=Parse.parseIntKMG(b);
			}else if(a.equals("samplerate")){
				samplerate=Float.parseFloat(b);
			}else if(a.equals("sampleseed") || a.equals("seed")){
				sampleseed=Long.parseLong(b);
			}else if(a.equals("normalize")){
				normalize=Parse.parseBoolean(b);
			}else if(a.equals("ref") || a.equals("fna")){
				fna=b;
			}else if(a.equals("gff")){
				gff=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				//				throw new RuntimeException("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
	}
	
	void process(Timer t){
		if(verbose){System.err.println("Setting sketch params.");}
		setSketchStatics();
		
		SketchTool canonTool=new SketchTool(sketchSize, 0, true, true, true);
		SketchTool forwardTool=new SketchTool(sketchSize, 0, true, true, false);
		
		FASTQ.PAIR_READS=false;
		final int threads=Tools.min(Shared.threads(), 16); //Can't seem to scale beyond around 8 threads...
		
		System.err.println("Making canonical sketch.");
		Sketch canonSketch=canonTool.processReadsMT(in1, threads, maxReads, SketchObject.ONE_SKETCH,
				samplerate, 0, 0, (byte)0, false);
		
		System.err.println("Making forward sketch.");
		Sketch fwdSketch=forwardTool.processReadsMT(in1, threads, maxReads, SketchObject.ONE_SKETCH,
				samplerate, 0, 0, (byte)0, false);
		
		final double[] refResults=calcPMRatioWithRef(canonTool, forwardTool, canonSketch, fwdSketch);
		
		FASTQ.PAIR_READS=true;
		if(verbose){outstream.println("Finished reading data.");}
		
		outstream.println();
		double[] results=calcStrandedness(canonSketch, fwdSketch);
		if(normalize) {results=calcStrandednessNormalized(canonSketch, fwdSketch);}
		
		outputResults(results, refResults);
		
		outstream.println();
		t.stop("Time:\t");
		assert(!errorState) : "An error was encountered.";
	}
	
	private void setSketchStatics() {
		SketchObject.AUTOSIZE=false;
		SketchObject.k=32;
		SketchObject.k2=-1;//Otherwise the estimates are too high
		SketchObject.setK=true;
		SketchObject.sampleseed=sampleseed;
//		SketchObject.defaultParams.minKeyOccuranceCount=2;
		
		SketchObject.AUTOSIZE=false;
		SketchObject.AUTOSIZE_LINEAR=false;
		SketchObject.targetSketchSize=sketchSize;
		SketchObject.SET_TARGET_SIZE=true;
		SketchObject.processSSU=false;
		
		SketchObject.defaultParams.parse("trackcounts", "trackcounts", null);
		SketchObject.defaultParams.samplerate=samplerate;
		SketchObject.postParse();
	}
	
	//Old version
//	@Deprecated
//	private double[] calcPMRatio_old(SketchTool canonTool, SketchTool fwdTool, 
//			Sketch canonSketch, Sketch fwdSketch) {
//		
//		FASTQ.PAIR_READS=false;
//		SketchObject.rcomp=false;
//		
//		if(fna==null) {return null;}
//
//		ArrayList<Read> genes=null;
//		if(gff!=null) {
//			ArrayList<GffLine> lines=getGffLines(gff, types);
//			HashMap<String, Read> seqMap=getSequenceMap(fna);
//			genes=CheckStrand.grabGenes(lines, seqMap);
//		}else{
//			System.err.println("Calling genes.");
//			genes=callGenes(fna);
//		}
//		if(genes==null || genes.isEmpty()) {return null;}
//		System.err.println("Processing "+genes.size()+" genes.");
//		
//		Sketch[] geneSketches=sketchGenes(genes, canonTool, fwdTool);
//		Sketch canonGeneSketch=geneSketches[0], plusSketch=geneSketches[1], minusSketch=geneSketches[2];
//		
//		double[] p=CheckStrand.countSharedSum(fwdSketch, plusSketch);
//		double[] m=CheckStrand.countSharedSum(fwdSketch, minusSketch);
//		
//		double ratio=p[3]/(p[3]+m[3]);
//		double[] abFractions=CheckStrand.calcCoverage(canonSketch, canonGeneSketch);
//		
//		double[] ret={ratio, abFractions[0], abFractions[1]};
//		return ret;
//	}
	
	/**
	 * Calculate the Plus/Minus mapping ratio using transcriptome kmers.
	 * @param canonTool SketchTool for canonical kmers.
	 * @param fwdTool SketchTool for forward kmers.
	 * @param canonSketch Read Sketch using canonical kmers.
	 * @param fwdSketch Read Sketch using forward kmers.
	 * @return Array of results (see calcPMRatio).
	 */
	private double[] calcPMRatioWithRef(SketchTool canonTool, SketchTool fwdTool, 
			Sketch canonSketch, Sketch fwdSketch) {
		ArrayList<Read> genes=grabGenes();
		if(genes==null) {return null;}
		outstream.println("Processing "+genes.size()+" genes.");
		Sketch[] geneSketches=sketchGenes(genes, canonTool, fwdTool);
		Sketch canonGeneSketch=geneSketches[0], plusSketch=geneSketches[1], minusSketch=geneSketches[2];
		return calcPMRatioWithGeneSketches(canonSketch, fwdSketch, canonGeneSketch, plusSketch, minusSketch);
	}
	
	/** 
	 * Generate a list of gene sequences from the reference.
	 * @return The list.
	 */
	private ArrayList<Read> grabGenes(){
		if(fna==null) {return null;}

		ArrayList<Read> genes=null;
		if(gff!=null) {
			ArrayList<GffLine> lines=getGffLines(gff, types);
			HashMap<String, Read> map=getSequenceMap(fna);
			genes=CheckStrand.grabGenes(lines, map);
		}else{
			System.err.println("Calling genes.");
			genes=callGenes(fna);
		}
		return (genes==null || genes.isEmpty()) ? null : genes;
	}
	
	/**
	 * Calculate the Plus/Minus mapping ratio using transcriptome kmers.
	 * @param canonSketch Read Sketch using canonical kmers.
	 * @param fwdSketch Read Sketch using forward kmers.
	 * @param canonGeneSketch Transcriptome Sketch using canonical kmers.
	 * @param plusSketch Transcriptome Sketch using plus-strand forward kmers.
	 * @param minusSketch Transcriptome Sketch using minus-strand forward kmers.
	 * @return Results vector: {ratio, abFractions[0], abFractions[1]}
	 */
	static double[] calcPMRatioWithGeneSketches(Sketch canonSketch, Sketch fwdSketch, 
			Sketch canonGeneSketch, Sketch plusSketch, Sketch minusSketch) {
		double[] p=CheckStrand.countSharedSum(fwdSketch, plusSketch);
		double[] m=CheckStrand.countSharedSum(fwdSketch, minusSketch);
		
		double ratio=p[3]/(p[3]+m[3]);
		double[] abFractions=CheckStrand.calcCoverage(canonSketch, canonGeneSketch);
		
		double[] ret={ratio, abFractions[0], abFractions[1]};
		return ret;
	}
	
	/**
	 * Create sketches of the transcriptome.
	 * @param genes List of gene sequences (sense on plus strand).
	 * @param canonTool SketchTool for canonical kmers.
	 * @param fwdTool SketchTool for forward kmers.
	 * @return Vector of: {canonGeneSketch, plusSketch, minusSketch}
	 */
	static Sketch[] sketchGenes(ArrayList<Read> genes, SketchTool canonTool, SketchTool fwdTool) {
		SketchMakerMini smmPlus=new SketchMakerMini(fwdTool, SketchObject.ONE_SKETCH, SketchObject.defaultParams);
		SketchMakerMini smmMinus=new SketchMakerMini(fwdTool, SketchObject.ONE_SKETCH, SketchObject.defaultParams);
		SketchMakerMini smmCanon=new SketchMakerMini(canonTool, SketchObject.ONE_SKETCH, SketchObject.defaultParams);
		for(Read r : genes) {
			smmCanon.processRead(r);
			smmPlus.processRead(r);
			r.reverseComplement();
			r.setStrand(0);
			smmMinus.processRead(r);
		}
		Sketch canonGeneSketch=smmCanon.toSketch(0);
		Sketch plusSketch=smmPlus.toSketch(0);
		Sketch minusSketch=smmMinus.toSketch(0);
		return new Sketch[] {canonGeneSketch, plusSketch, minusSketch};
	}
	
	/**
	 * Print the final program results.
	 * @param results Results from read kmer depth analysis
	 * @param refResults Results based on transcriptome kmer comparison
	 */
	private void outputResults(double[] results, double[] refResults){
		if(ffout1==null) {return;}
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();
		
		double strandedness=results[4];
		double depth=results[5];
		
		//Write stuff to the bsw
		bsw.println(String.format("Strandedness:\t%.2f%%", strandedness*100));
		bsw.println(String.format("AvgKmerDepth:\t%.2f", depth));
		if(refResults!=null) {
			double pmRatio=refResults[0];
			double aFraction=refResults[1];
			double bFraction=refResults[2];
			bsw.println(String.format("P/(P+M)_Ratio:\t%.6f", pmRatio));
			bsw.println("MajorStrand:\t"+((pmRatio>=0.5) ? "Plus" : "Minus"));
			bsw.println(String.format("GeneCoverage:\t%.4f", aFraction));
			bsw.println(String.format("GenePrecision:\t%.4f", bFraction));
		}
		
		errorState=bsw.poisonAndWait() | errorState;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Determine the strandedness of a set of reads by comparing a Sketch of
	 * canonical kmers to forward kmers.
	 * @param saCanon Sketch of canonical kmers.
	 * @param sbFwd Sketch of forward kmers.
	 * @return Results vector: {totalSum, minSum, expectedMinSum, maxPossibleMinSum, 
	 * strandedness, depth, matches, nonUniqueFraction}
	 */
	static double[] calcStrandedness(Sketch saCanon, Sketch sbFwd) {
		final long[] a=saCanon.keys, b=sbFwd.keys;
		final int[] aCounts=saCanon.keyCounts, bCounts=sbFwd.keyCounts;
		
		int matches=0;
		int totalCount=0;
		long sharedSum=0;
		long nonUniqueCount=0;
		double totalSum=0;
		double minSum=0;
		double expectedMinSum=0;
		double maxPossibleMinSum=0;
		
		//Here we walk down the Sketches and find where they share kmers.
		//In those cases, the counts are compared to determine balance.
		for(int i=0, j=0; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){//Match
				matches++;
				final int ca=aCounts[i];
				final int cb=bCounts[j];
				final int cmin=Tools.min(cb, ca-cb);
				assert(cmin>=0) : ca+", "+cb;
				nonUniqueCount+=(ca>1 ? 1 : 0);
				totalSum+=ca;
				totalCount++;
				sharedSum+=ca;
				minSum+=cmin;
				maxPossibleMinSum+=(ca/2);
				expectedMinSum+=expectedMinorAlleleCount(ca);
				i++;
				j++;
			}else if(ka<kb){//kb was missing; thus it had a minor count of 0
				final int ca=aCounts[i];
				final int cb=0;
				final int cmin=Tools.min(cb, ca-cb);//Should always be cb
				assert(cmin>=0);//Should always be 0
				nonUniqueCount+=(ca>1 ? 1 : 0);
				totalSum+=ca;
				totalCount++;
				minSum+=cmin;
				maxPossibleMinSum+=(ca/2);
				expectedMinSum+=expectedMinorAlleleCount(ca);
				i++;
			}else{//ka was missing; thus this is a noncanonical key
				j++;
			}
		}
		
		//Strandedness will be (0.5-1.0) for normal (fully unstranded-fully stranded) libraries.
		//Synthetic or binned libraries with a perfectly flat distribution will get 0.0,
		//but anything between 0.0 and 0.5 would be unusual.
		//Basically 1.0 is preference for a strand, 0.5 is no preference for a strand,
		//and 0.0 is preference for perfect balance between strands - as you would get
		//if you treated paired reads as single-ended.
		double strandedness;
		if(minSum<=expectedMinSum) {//Normal case; strandedness between 0 and 50%
//			strandedness=1-(minSum/(expectedMinSum+minSum));//Not sure about the proper formula; needs thought
			strandedness=0.5+(1-(minSum/expectedMinSum))*0.5;
		}else{//Odd case; distribution is more even than expected by chance
			assert(minSum<=maxPossibleMinSum) : minSum+", "+maxPossibleMinSum;
			double range=(maxPossibleMinSum-expectedMinSum);
			double delta=minSum-expectedMinSum;
			assert(delta>=0 && delta<=range);
			double x=0.5*(1-(delta/range));
			strandedness=x;//Not really sure about this either
		}
		double depth=totalSum/totalCount;
		double nonUniqueFraction=nonUniqueCount/(1.0*totalCount);
		return new double[] {totalSum, minSum, expectedMinSum, maxPossibleMinSum, 
				strandedness, depth, matches, nonUniqueFraction};
	}
	
	/**
	 * Determine the strandedness of a set of reads by comparing a Sketch of
	 * canonical kmers to forward kmers.
	 * In this mode, each kmer contributes 0-1.0 strandedness regardless of its depth.
	 * @param saCanon Sketch of canonical kmers.
	 * @param sbFwd Sketch of forward kmers.
	 * @return Results vector: {totalSum, minSum, expectedMinSum, maxPossibleMinSum, 
	 * strandedness, depth, matches, nonUniqueFraction}
	 */
	static double[] calcStrandednessNormalized(Sketch saCanon, Sketch sbFwd) {
		final long[] a=saCanon.keys, b=sbFwd.keys;
		final int[] aCounts=saCanon.keyCounts, bCounts=sbFwd.keyCounts;
		
		int matches=0;
		int totalCount=0;
		long sharedSum=0;
		long nonUniqueCount=0;
		double totalSum=0;
		double minSum=0;
		double expectedMinSum=0;
		double maxPossibleMinSum=0;
		
		//Here we walk down the Sketches and find where they share kmers.
		//In those cases, the counts are compared to determine balance.
		for(int i=0, j=0; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){//Match
				final int ca=aCounts[i];
				final int cb=bCounts[j];
				if(ca>1) {
					matches++;
					final int cmin=Tools.min(cb, ca-cb);
					assert(cmin>=0) : ca+", "+cb;
					nonUniqueCount+=(ca>1 ? 1 : 0);
					double expectedMinor=expectedMinorAlleleCount(ca);
					assert(expectedMinor>0);
					double ratio=cmin/(double)ca;
					totalSum+=ca;
					totalCount++;
					sharedSum+=ca;
					minSum+=ratio;
					maxPossibleMinSum+=((ca/2)/(double)ca);
					expectedMinSum+=(expectedMinor/ca);
				}
				i++;
				j++;
			}else if(ka<kb){//kb was missing; thus it had a minor count of 0
				final int ca=aCounts[i];//TODO
				final int cb=0;
				if(ca>1) {
					final int cmin=Tools.min(cb, ca-cb);//Should always be cb
					assert(cmin>=0);//Should always be 0
					nonUniqueCount+=(ca>1 ? 1 : 0);
					double expectedMinor=expectedMinorAlleleCount(ca);
					double ratio=cmin/(double)ca;
					totalSum+=ca;
					totalCount++;
					sharedSum+=ca;
					minSum+=ratio;
					maxPossibleMinSum+=((ca/2)/(double)ca);
					expectedMinSum+=(expectedMinor/ca);
				}
				i++;
			}else{//ka was missing; thus this is a noncanonical key
				j++;
			}
		}
		
		//Strandedness will be (0.5-1.0) for normal (fully unstranded-fully stranded) libraries.
		//Synthetic or binned libraries with a perfectly flat distribution will get 0.0,
		//but anything between 0.0 and 0.5 would be unusual.
		//Basically 1.0 is preference for a strand, 0.5 is no preference for a strand,
		//and 0.0 is preference for perfect balance between strands - as you would get
		//if you treated paired reads as single-ended.
		double strandedness;
		if(minSum<=expectedMinSum) {//Normal case; strandedness between 0 and 50%
//			strandedness=1-(minSum/(expectedMinSum+minSum));//Not sure about the proper formula; needs thought
			strandedness=0.5+(1-(minSum/expectedMinSum))*0.5;
		}else{//Odd case; distribution is more even than expected by chance
			assert(minSum<=maxPossibleMinSum) : minSum+", "+maxPossibleMinSum;
			double range=(maxPossibleMinSum-expectedMinSum);
			double delta=minSum-expectedMinSum;
			assert(delta>=0 && delta<=range);
			double x=0.5*(1-(delta/range));
			strandedness=x;//Not really sure about this either
		}
		double depth=totalSum/totalCount;
		double nonUniqueFraction=nonUniqueCount/(1.0*totalCount);
		return new double[] {totalSum, minSum, expectedMinSum, maxPossibleMinSum, 
				strandedness, depth, matches, nonUniqueFraction};
	}
	
	static float strandedness(long plus, long minus) {
		long sum=plus+minus;
		long maxPossibleMinor=sum/2;
		float expectedMinor=(float)expectedMinorAlleleCount(sum);
		assert(expectedMinor>0);
		
		return strandedness(plus, minus, maxPossibleMinor, expectedMinor);
	}
	
	static float strandedness(long plus, long minus, long maxPossibleMinor, float expectedMinor) {
		long sum=plus+minus;
		long min=Math.min(plus, minus);
		assert(expectedMinor>0);
		
		if(min<=expectedMinor) {//Expected case
			float strandedness=1-0.5f*(min/expectedMinor);//0.5-1
			assert(strandedness>=0 && strandedness<=1.01) : strandedness+
				", "+plus+", "+minus+", "+sum+", "+min+", "+expectedMinor;
			return strandedness;
		}else{//Rare case; overly unstranded
			float range=maxPossibleMinor-expectedMinor;
			float dif=min-expectedMinor;
			float strandedness=0.5f*(1-(dif/range));
			assert(strandedness>=0 && strandedness<=1.01) : strandedness+", "+plus+", "+minus+
				", "+sum+", "+min+", "+expectedMinor+", "+range+", "+dif;
			return strandedness;
		}
	}
	
	/**
	 * Counts the fraction of total kmers shared between sketches, which includes their counts.
	 * @param saFwd Forward sketch of reads.
	 * @param sbTranscriptStrand Forward (or reverse) sketch of transcriptome.
	 * @return Results vector: {totalSum, depth, matches, sharedSum}
	 */
	static double[] countSharedSum(Sketch saFwd, Sketch sbTranscriptStrand) {
		final long[] a=saFwd.keys, b=sbTranscriptStrand.keys;
		final int[] aCounts=saFwd.keyCounts, bCounts=sbTranscriptStrand.keyCounts;
		
		int matches=0;
		int totalCount=0;
		long sharedSum=0;
		double totalSum=0;
		
		for(int i=0, j=0; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){//Match
				matches++;
				final int ca=aCounts[i];
				totalSum+=ca;
				totalCount++;
				sharedSum+=ca;
				i++;
				j++;
			}else if(ka<kb){//kb was missing
				final int ca=aCounts[i];
				totalSum+=ca;
				totalCount++;
				i++;
			}else{//ka was missing
				j++;
			}
		}
		
		double depth=sharedSum/(double)Tools.max(matches, 1);
		return new double[] {totalSum, depth, matches, sharedSum};
	}
	
	/**
	 * Calculate each sketch's fractional coverage of the other sketch, ignoring counts.
	 * @return Results vector: {aFraction, bFraction},
	 * where aFraction is sketch a's coverage of sketch b.
	 */
	static double[] calcCoverage(Sketch sa, Sketch sb) {
		final long[] a=sa.keys, b=sb.keys;
		final int[] aCounts=sa.keyCounts, bCounts=sb.keyCounts;
		
		int matches=0;
		
		int i=0, j=0;
		for(; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){//Match
				matches++;
				i++;
				j++;
			}else if(ka<kb){//kb was missing
				i++;
			}else{//ka was missing
				j++;
			}
		}

		double aFraction=matches/(double)j;
		double bFraction=matches/(double)i;
		
		return new double[] {aFraction, bFraction};
	}
	
	/*--------------------------------------------------------------*/
	
	/** 
	 * This works for a diploid het allele, or a fair coin, or other things with 50/50 outcomes.
	 * Returns the expected minor allele frequency for a given depth.
	 * @param depth i.e., number of observations.
	 * @return Expected number of observations of the minor allele.
	 */
	public static final double expectedMinorAlleleFrequency(final long depth) {
		if(depth<expectedMinorAlleleFreq.length) {return depth>1 ? expectedMinorAlleleFreq[(int)depth] : 0;}
		long d=depth;
		double mult=1.0;
		
		//If the depth is greater than array length, downscale it.
		//This could alternatively use power-of-two multipliers
		//like 2 and 4, or 8 and 64, where one is the square of the other.
		//Or probably some formula using logs.
		//TODO: Consider switching to a log-based formula instead of a loop
		while(d>=expectedMinorAlleleFreq.length) {
			d=d/100;
			mult*=0.1;
		}
		double maf=expectedMinorAlleleFreq[(int)d];
		double dif=0.5-maf;
		return 0.5-(dif*mult);
	}
	
	public static final double expectedMinorAlleleCount(final long depth) {
		if(depth<expectedMinorAlleleCount.length) {return depth>1 ? expectedMinorAlleleCount[(int)depth] : 0;}
		return depth*expectedMinorAlleleFrequency(depth);
	}
	
	/** This is for generating the stats file for loading later; only needed once. */
	private static void printMinorAlleleCount() {
		System.out.println("#Expected minor allele count for N coin flips, starting at 0, 10m simulations.");
		for(int i=0; i<expectedMinorAlleleCount.length; i++) {
			double c=expectedMinorAlleleCount[i];
			int decimals=Tools.max(1, 7-Integer.toString((int)c).length());
			System.out.println(String.format("%."+decimals+"f", expectedMinorAlleleCount[i]));
		}
	}
	
	/** 
	 * Runs a simulation.  This will happen automatically if minorAlleleCount.txt is not found,
	 * but it will be less precise due to fewer trials.  The total number of coin flips
	 * is maxDepth*trials.
	 * @param maxDepth Max total allele count (coin flips in a series).
	 * @param trials Number of simulated series.
	 * @return Array of average minor allele counts.
	 */
	public static final double[] makeExpectedMinorAlleleArray(final int maxDepth, final int trials) {
		{
			double[] d=loadExpectedMinorAlleleArray();
			if(d!=null) {return d;}
		}
		final long[] minorSum=new long[maxDepth+1];
		final Random randy=Shared.threadLocalRandom();
		final int[] headsTails=new int[2];
		for(int trial=0; trial<trials; trial++) {
			headsTails[0]=headsTails[1]=0;
			for(int i=1; i<maxDepth+1; i++) {
				int bit=randy.nextInt()&1;
				headsTails[bit]++;
				minorSum[i]+=Tools.min(headsTails[0], headsTails[1]);
			}
		}
		final double[] expected=new double[maxDepth+1];
		final double mult=1.0/trials;
		for(int i=0; i<maxDepth+1; i++) {
			expected[i]=minorSum[i]*mult;
		}
		return expected;
	}
	
	/** 
	 * Load minor allele counts from a file.
	 * It should be in bbmap/resources/minorAlleleCount.txt
	 * */
	public static double[] loadExpectedMinorAlleleArray() {
		String path=Data.findPath("?minorAlleleCount.txt", true);
		if(path==null) {return null;}
		File f=new File(path);
		if(!f.exists() || !f.canRead()) {return null;}
		
		DoubleList dl=new DoubleList();
		String[] lines=TextFile.toStringLines(path);
		for(String s : lines) {
			if(!Tools.startsWith(s, '#')) {
				double d=Parse.parseDouble(s, 0, s.length());
				dl.add(d);
			}
		}
		return dl.toArray();
	}
	
	/** Make the frequency array from the count array */
	private static double[] makeExpectedMinorAlleleFreq(double[] counts) {
		double[] freq=new double[counts.length];
		for(int i=1; i<counts.length; i++) {
			freq[i]=counts[i]/i;
		}
		return freq;
	}
	
	/*--------------------------------------------------------------*/

	/**
	 * Load a gff file.
	 * @param gff The file path.
	 * @param types Types of features to load, such as "CDS,rRNA".
	 * @return
	 */
	static ArrayList<GffLine> getGffLines(String gff, String types){
		return GffLine.loadGffFile(gff, types, false);
	}

	/**
	 * Load a fasta file as a HashMap of names to sequences.
	 * @param fna Fasta file.
	 * @return Map of names to sequences.
	 */
	static HashMap<String, Read> getSequenceMap(String fna){
		ArrayList<Read> list=ReadInputStream.toReads(fna, FileFormat.FA, -1);
		return getSequenceMap(list);
	}
	
	/**
	 * Generate a map of names to sequences from a list of sequences.
	 * Also maps name prefix up to the first whitespace.
	 * @param list List of sequences.
	 * @return The map.
	 */
	static HashMap<String, Read> getSequenceMap(ArrayList<Read> list){
		HashMap<String, Read> map=new HashMap<String, Read>(1+list.size()*3);
		
		for(Read r : list){
			map.put(r.id, r);
			//Faster to use a lineparser but they don't support whitespaceplus
			String id2=Tools.whitespacePlus.split(r.id)[0];
			if(!id2.equals(r.id)) {map.put(id2, r);}
		}
		return map;
	}
	
	/**
	 * Generates a list of sequences by cutting out the specified regions.
	 * Intended for generating gene sequences given a list of GffLines.
	 * @param <K> A Feature such as a GffLine.
	 * @param lines List of features.
	 * @param map Map of name to sequence (for the reference genome).
	 * @return A list of sequences of the input features, named by the features.
	 */
	static <K extends Feature> ArrayList<Read> grabGenes(ArrayList<? extends K> lines, HashMap<String, Read> map){
		ArrayList<Read> list=null; 
		
//		HashSet<String> set=new HashSet<String>();
//		for(String s : types) {set.add(s);}
		
		for(K gline : lines){
//			if(set.contains(gline.type)){
				Read scaf=map.get(gline.seqid());
				assert(scaf!=null) : "Can't find "+gline.seqid()+" in "+map.keySet();
				
				final int start=gline.start();
				final int stop=gline.stop();
				
				
				if(start>=0 && stop<scaf.length()){
					String id=gline.name();
					Read r=new Read(Arrays.copyOfRange(scaf.bases, start, stop+1), null, id, 1);

//					assert(!r.containsLowercase()) : r.toFasta()+"\n"
//					+ "validated="+r.validated()+", scaf.validated="+scaf.validated()+", tuc="+Read.TO_UPPER_CASE+", vic="+Read.VALIDATE_IN_CONSTRUCTOR;
					if(r!=null){
						if(gline.strand()==1){r.reverseComplement();}
						if(list==null){list=new ArrayList<Read>(8);}
						list.add(r);
					}
//				}
					
			}
		}
		return list;
	}
	
	/**
	 * Call genes from a reference file and return the gene sequences.
	 * @param fna Fasta reference.
	 * @return Gene sequences.
	 */
	private ArrayList<Read> callGenes(String fna){
		final ConcurrentReadInputStream cris=makeFastaCris(fna);
		
		if(pgmFile==null){
			pgmFile=Data.findPath("?model.pgm");
		}
		GeneModel pgm=GeneModelParser.loadModel(pgmFile);
		GeneCaller gCaller=CallGenes.makeGeneCaller(pgm);
		ArrayList<Read> genes=callGenes(cris, gCaller);
		
		//Close the input stream
		errorState|=ReadWrite.closeStream(cris);
		return genes;
	}
	
	/**
	 * Makes a read input stream for a file (assumed to be fasta).
	 * @param fname File path.
	 * @return The read input stream.
	 */
	private ConcurrentReadInputStream makeFastaCris(String fname){
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ffin, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/**
	 * Call genes from a read stream and return the gene sequences.
	 * @param cris Read stream.
	 * @param gCaller The gene caller.
	 * @return Gene sequences.
	 */
	static ArrayList<Read> callGenes(ConcurrentReadInputStream cris, GeneCaller gCaller){
		ArrayList<Read> genes=new ArrayList<Read>();

		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();

		CallGenes.callCDS=true;
		CallGenes.calltRNA=CallGenes.call16S=CallGenes.call23S
				=CallGenes.call5S=CallGenes.call18S=true;
		
		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){

			for(Read r : ln) {
				ArrayList<Orf> orfs=gCaller.callGenes(r);
				if(orfs!=null) {
					for(Orf orf : orfs) {
						Read gene=CallGenes.fetch(orf, r);
						genes.add(gene);
					}
				}
//				System.err.println(r.length()+", "+orfs.size()+", "+orfs);
			}
			//Fetch a new list
			ln=cris.nextList();
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		return genes;
	}
	
	/*--------------------------------------------------------------*/
	
	private int sketchSize=20000;
	
	private String in1=null;
	private String out1="stdout.txt";
	String fna=null;
	String gff=null;
	String pgmFile=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/** Features to pull from gff files */
	private String types="CDS,rRNA,tRNA,ncRNA,exon,5S,16S,23S";
	
	/*--------------------------------------------------------------*/

	private boolean normalize=false;
	private long maxReads=-1;
	private float samplerate=1;
	private long sampleseed=17;
	private boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	
	static final double[] expectedMinorAlleleCount=
			makeExpectedMinorAlleleArray(10000, 100000);
	static final double[] expectedMinorAlleleFreq=
			makeExpectedMinorAlleleFreq(expectedMinorAlleleCount);
	
	/*--------------------------------------------------------------*/
	
	/** Output screen messages here */
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
