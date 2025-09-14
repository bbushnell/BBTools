package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import gff.GffLine;
import prok.CallGenes;
import prok.GeneCaller;
import prok.GeneModel;
import prok.GeneModelParser;
import prok.Orf;
import prok.ProkObject;
import shared.LineParser4;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.Sketch;
import sketch.SketchHeap;
import sketch.SketchMakerMini;
import sketch.SketchObject;
import sketch.SketchTool;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.ReadInputStream;
import stream.SamLine;
import structures.ListNum;
import structures.LongPair;
import structures.Range;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Checks the strandedness of RNA-seq reads.
 * 
 * @author Brian Bushnell
 * @date August 7, 2023
 *
 */
public class CheckStrand2 implements Accumulator<CheckStrand2.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		CheckStrand2 x=new CheckStrand2(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CheckStrand2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;
			extin=parser.extin;

			out1=parser.out1;
			extout=parser.extout;
		}

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, extout, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		//Allows unchanged sam output for sam input.
		SamLine.SET_FROM_OK=true;
		
		//If there is a gff file that will not be used, nullify it.
		if(gff!=null) {
			if(transcriptome) {
				System.err.println("Ignoring gff file due to transcriptome mode.");
				gff=null;
			}else if(!ffin1.samOrBam() && fna==null) {
				System.err.println("Ignoring gff file due to using unaligned reads with no reference.");
				gff=null;
			}
		}
		
		//Decide whether gene-calling is needed, and pick a gene model.
		if((fna!=null && (gff==null || passes>1) && !transcriptome && !rnaContigs) || scoreReadGenes) {
			if(pgmFile==null){pgmFile=Data.findPath("?model.pgm");}
			final GeneModel pgm0=GeneModelParser.loadModel(pgmFile);
			if(passes<2 || fna==null || transcriptome) {pgm=pgm0;}
			else {
				Timer t=new Timer();
				System.err.print("Refining gene model: ");
				pgm=CallGenes.makeMultipassModel(pgm0, genomeSequence(), gff, passes);
				if(streamGenome) {genomeSequenceCache=null;}//To save memory
				t.stopAndPrint();
			}
			gCaller=CallGenes.makeGeneCaller(pgm);
		}else {
			pgm=null;
			gCaller=null;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
		parser.out1=out1;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(a.equals("size") || a.equals("len") || a.equals("length") || 
					a.equals("sketchsize") || a.equals("sketchlen")){
				sketchSize=Parse.parseIntKMG(b);
			}else if(a.equals("samplerate")){
				samplerate=Float.parseFloat(b);
			}else if(a.equals("sampleseed") || a.equals("seed")){
				sampleseed=Long.parseLong(b);
			}else if(a.equals("ref") || a.equals("fna")){
				fna=b;
			}else if(a.equals("gff")){
				gff=b;
			}else if(a.equals("merge")){
				mergePairs=Parse.parseBoolean(b);
			}else if(a.equals("bin") || a.equals("binlen")){
				binlen=Parse.parseIntKMG(b);
			}else if(a.equals("scoregenes") || a.equals("scorereads") || a.equals("callgenes")
					|| a.equals("callreads") || a.equals("kfa")){
				scoreReadGenes=Parse.parseBoolean(b);
			}else if(a.equals("pgm") || a.equals("model")){
				pgmFile=b;
			}else if(a.equals("streamgenome")){
				streamGenome=Parse.parseBoolean(b);
			}else if(a.equals("stopanalysis") || a.equals("orf")){
				doStopAnalysis=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("polya")){
				testPolyA=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("transcriptome")){
				transcriptome=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("rnacontigs")){
				rnaContigs=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("normalize")){
				normalize=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("minreads")){
				minReads=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("genome")){
				transcriptome=!Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("outp") || a.equalsIgnoreCase("outplus")){
				outPlus=b;
			}else if(a.equalsIgnoreCase("outm") || a.equalsIgnoreCase("outminus")){
				outMinus=b;
			}else if(a.equalsIgnoreCase("scaffoldReport") || a.equalsIgnoreCase("scafReport")){
				scaffoldReport=b;
			}else if(a.equalsIgnoreCase("firstorf")){
				useFirstORF=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("minpolya") || a.equalsIgnoreCase("polyalen")
					|| a.equalsIgnoreCase("polyalength")){
				minPolyA=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("passes")){
				passes=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("2pass") || a.equalsIgnoreCase("twopass")){
				if(Parse.parseBoolean(b)) {passes=2;}
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				boolean parsed=false;
				if(b==null && arg.indexOf('=')<0) {
					File f=new File(arg);
					String ext=ReadWrite.rawExtension(a);
					if(f.isFile() && f.canRead()) {
						if(parser.in1==null && (FileFormat.isSamOrBamExt(ext) || FileFormat.isFastqExt(ext))) {
							outstream.println("Set input to "+f.getName());
							parser.in1=arg;
							parsed=true;
						}else if(gff==null && FileFormat.isGffExt(ext)) {
							outstream.println("Set gff to "+f.getName());
							gff=arg;
							parsed=true;
						}else if(fna==null && FileFormat.isFastaExt(ext)) {
							outstream.println("Set ref to "+f.getName());
							fna=arg;
							parsed=true;
						}
					}
				}
				if(!parsed) {
					outstream.println("Unknown parameter "+args[i]);
					assert(false) : "Unknown parameter "+args[i];
				}
			}
		}
		
		return parser;
	}
	
	/** Replace # with 1 and 2 in headers */
	private void doPoundReplacement(){
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		qfin1=Tools.fixExtension(qfin1);
		qfin2=Tools.fixExtension(qfin2);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, outPlus, outMinus)){
			outstream.println((out1==null)+", "+(outPlus==null)+", "+(outMinus==null)+
					", "+out1+", "+outPlus+", "+outMinus);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; "
					+ "Can't write to output files "+out1+", "+outPlus+", "+outMinus+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, outPlus, outMinus)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Make sure interleaving agrees with number of input and output files */
	private void adjustInterleaving(){
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
//				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		makeTools();
		makeGeneSketches();
		
		if(ffin1.samOrBam()) {//Aligned reads
			if(gffLines==null) {//Features not yet loaded
				if(gff!=null) {
					gffLines=CheckStrand.getGffLines(gff, types);
				}else if(fna!=null && !transcriptome && !rnaContigs) {
					ArrayList<Read> reads=genomeSequence();
					System.err.print("Calling genes from genome: ");
					ArrayList<Orf> orfs=gCaller.callGenes(reads);
					gffLines=Orf.toGffLines(orfs);
				}
			}
			if(gffLines!=null) {
				SamLine.RNAME_AS_BYTES=false;//Must come before cris starts
				rangeMap=GffLine.makeRangeMap(gffLines);
			}
			if(binlen>0 || scaffoldReport!=null) {
				if(binlen<1) {binlen=1000;}
				binMap=Seq.loadSamHeader(ffin1);
			}
		}
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris();
		final ConcurrentReadOutputStream crosP=makeCros(outPlus);
		final ConcurrentReadOutputStream crosM=makeCros(outMinus);
		
		//Reset counters
		readsIn=readsProcessed=0;
		basesIn=basesProcessed=0;
		readsMerged=0;
		
		//Process the reads in separate threads
		spawnThreads(cris, crosP, crosM);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, crosP, crosM);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;

		outstream.println();
		analyze();
		
		//Report timing and results
		t.stop();
		
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
//		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Do math on the raw results to determine strandedness and major strand */
	void analyze() {
		double[] results=CheckStrand.calcStrandedness(canonSketch, fwdSketch);
		double[] resultsN=CheckStrand.calcStrandednessNormalized(canonSketch, fwdSketch);
		
		double[] refResults=(canonGeneSketch==null ? null : CheckStrand.calcPMRatioWithGeneSketches(
				canonSketch, fwdSketch, canonGeneSketch, plusSketch, minusSketch));
		outputResults(results, resultsN, refResults);
	}
	
	/**
	 * Start a read input stream from the read files
	 * @return The stream.
	 */
	private ConcurrentReadInputStream makeCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(
				maxReads, true, ffin1, ffin2, qfin1, qfin2);
		cris.setSampleRate(samplerate, sampleseed);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){//Announce whether the data is being procesed as paired
			outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));
		}
		return cris;
	}

	/** 
	 * Start a read output stream for some destination
	 * @param fname Destination file for stream
	 * @return The stream.
	 */
	private ConcurrentReadOutputStream makeCros(String fname){
		if(fname==null) {return null;}
		
		//Select output buffer size based on whether it needs to be ordered
		final boolean ordered=true;
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		final FileFormat ffout=FileFormat.testOutput(
				fname, ffin1.format(), null, true, overwrite, false, ordered);
		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(
				ffout, null, buff, null, true);
		ros.start(); //Start the stream
		return ros;
	}
	
	/**
	 * Set Sketch static parameters and generate SketchTools
	 * for forward and canonical kmers.
	 */
	private void makeTools() {

		if(verbose){System.err.println("Setting sketch params.");}
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
		
		canonTool=new SketchTool(sketchSize, 0, true, true, true);
		fwdTool=new SketchTool(sketchSize, 0, true, true, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, 
			ConcurrentReadOutputStream crosP, ConcurrentReadOutputStream crosM){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int maxThreads=mergePairs ? 64 : 16; //This is to limit memory use by sketches
		final int threads=Tools.min(maxThreads, Shared.threads());
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, crosP, crosM, i));
		}
		
		//Start the threads and wait for them to finish
		ThreadWaiter.startThreads(alpt);
//		
//		//While that's going, process the ref
//		makeGeneSketches();
		
		boolean success=ThreadWaiter.waitForThreadsToFinish(alpt, this);
		outstream.println("Finished processing "+readsIn+" reads and "+basesIn+" bases.");
		errorState&=!success;
		
		//Do anything necessary after processing
		//This transforms SketchHeaps (used for building Sketches) into finished Sketches.
		canonSketch=new Sketch(canonHeap, false, true, null);
		fwdSketch=new Sketch(fwdHeap, false, true, null);
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt) {
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			readsMerged+=pt.readsMergedT;

			fCountSum+=pt.fCountSumT;
			rCountSum+=pt.rCountSumT;
			fPosSum+=pt.fPosSumT;
			rPosSum+=pt.rPosSumT;

			Tools.add(fBestFrame, pt.fBestFrameT);
			Tools.add(rBestFrame, pt.rBestFrameT);
			Tools.add(ACGTCount, pt.ACGTCountT);

			polyACount+=pt.polyACountT;
			polyTCount+=pt.polyTCountT;

			gCallerPlusCount+=pt.gCallerPlusCountT;
			gCallerMinusCount+=pt.gCallerMinusCountT;
			gCallerPlusCalled+=pt.gCallerPlusCalledT;
			gCallerMinusCalled+=pt.gCallerMinusCalledT;
			gCallerPlusScore+=pt.gCallerPlusScoreT;
			gCallerMinusScore+=pt.gCallerMinusScoreT;

			plusAlignedReads+=pt.plusAlignedReadsT;
			minusAlignedReads+=pt.minusAlignedReadsT;
			allAlignedReads+=pt.allAlignedReadsT;

			samLinesProcessed+=pt.samLinesProcessedT;
			samLinesAlignedToFeatures+=pt.samLinesAlignedToFeaturesT;
			
			errorState|=(!pt.success);

			if(!pt.canonSmm.isEmpty()) {
				if(canonHeap==null){canonHeap=pt.canonSmm.heap();}
				else{canonHeap.add(pt.canonSmm.heap());}
			}
			if(!pt.fwdSmm.isEmpty()) {
				if(fwdHeap==null){fwdHeap=pt.fwdSmm.heap();}
				else{fwdHeap.add(pt.fwdSmm.heap());}
			}
		}
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Generate sketches of the reference transcriptome.
	 * @return True if successful.
	 */
	private boolean makeGeneSketches() {
		genes=grabGenes();
		if(genes==null) {return false;}
		outstream.println("Processing "+genes.size()+" genes.");
		Sketch[] geneSketches=CheckStrand.sketchGenes(genes, canonTool, fwdTool);
		if(geneSketches==null) {return false;}
		canonGeneSketch=geneSketches[0];
		plusSketch=geneSketches[1];
		minusSketch=geneSketches[2];
		return true;
	}
	
	/**
	 * Load the reference transcriptome; requires a reference fasta file.
	 * In transcriptome mode, the fasta file is simply loaded and returned.
	 * In genome mode, uses a gff if provided, otherwise the genes are called
	 * using CallGenes.  Then the genes are cut from the reference like CutGff.
	 * Does not currently fuse exons together, and is untested on Euk annotations,
	 * but they are expected to work fine (as long as a gff is provided).
	 * @return Gene sequences, sense strand as plus.
	 */
	private ArrayList<Read> grabGenes(){
		if(fna==null) {return null;}
		
		final ArrayList<Read> genes;
		Timer t=new Timer();
		
		if(transcriptome || rnaContigs) {
			System.err.print("Loading genes from transcriptome: ");
			genes=genomeSequence(); //Skips loading if already cached
		}else if(gff!=null) {
			System.err.print("Loading genes from genome and gff: ");
			HashMap<String, Read> map=(genomeSequenceCache==null ? 
					CheckStrand.getSequenceMap(fna) : CheckStrand.getSequenceMap(genomeSequence()));
			genes=CheckStrand.grabGenes(gffLines, map);
		}else{
			ArrayList<Read> reads=genomeSequence();
			if(gffLines==null) {
				System.err.print("Calling genes from genome: ");
				ArrayList<Orf> orfs=gCaller.callGenes(reads);
				gffLines=Orf.toGffLines(orfs);
			}
			HashMap<String, Read> map=CheckStrand.getSequenceMap(reads);
			genes=CheckStrand.grabGenes(gffLines, map);
		}
		t.stopAndPrint();
		
		return (genes==null || genes.isEmpty()) ? null : genes;
	}
	
	/**
	 * Calls genes while streaming the fasta to save memory.
	 * Works fine but for some reason is extremely slow.
	 * @param fna Genome fasta.
	 * @return Gene sequences.
	 */
	@Deprecated
	private ArrayList<Read> callGenes(String fna){
		final ConcurrentReadInputStream cris=makeFastaCris(fna);
		
		ArrayList<Read> genes=CheckStrand.callGenes(cris, gCaller);
		
		//Close the input stream
		errorState|=ReadWrite.closeStream(cris);
		return genes;
	}
	
	/**
	 * Creates a read input stream for the fasta reference.
	 * @param fname Fasta path.
	 * @return The stream.
	 */
	private ConcurrentReadInputStream makeFastaCris(String fname){
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ffin, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/**
	 * Print the final program results.
	 * @param results Results from read kmer depth analysis
	 * @param refResults Results based on transcriptome kmer comparison
	 */
	private void outputResults(double[] results, double[] resultsN, double[] refResults){
		if(ffout1==null) {return;}
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();

		final double invReads=1.0/readsProcessed;
		final double invBases=1.0/basesProcessed;
		
		//Results based on read kmer depth
		if(doDepthAnalysis){
			double strandedness=results[4];
			double strandednessN=resultsN[4];
			double depth=results[5];
			double nonUniqueFraction=results[7];
			
			bsw.println("Depth_Analysis:");
			bsw.println(String.format("Strandedness:\t%.2f%%", strandedness*100));
			bsw.println(String.format("StrandednessN:\t%.2f%%", strandednessN*100));
			bsw.println(String.format("AvgKmerDepth:\t%.3f", depth));
			bsw.println(String.format("Kmers>Depth1:\t%.4f", nonUniqueFraction));
		}

		//Results based on stop codons
		if(doStopAnalysis){
			bsw.println();
			boolean plus=(fPosSum>=rPosSum);
			double gc=(ACGTCount[1]+ACGTCount[2])/(1.0*shared.Vector.sum(ACGTCount));
			bsw.println("Stop_Codon_Analysis:");
			bsw.println(String.format("MajorStrandORF:\t"+(plus ? "Plus" : "Minus")));
			bsw.println(String.format("AvgReadLen:\t%.2f", basesProcessed*invReads));
			bsw.println(String.format("AvgORFLen+:\t%.2f", fPosSum*invReads));
			bsw.println(String.format("AvgORFLen-:\t%.2f", rPosSum*invReads));
			bsw.println(String.format("AvgStopCount+:\t%.4f", fCountSum*invReads));
			bsw.println(String.format("AvgStopCount-:\t%.4f", rCountSum*invReads));
			bsw.println(String.format("GC_Content:\t%.4f", gc));
			
			//Frame analysis seemed useless so I suppressed it,
			//but this prints which frame has the longest ORF
			long[] bestFrame=(plus ? fBestFrame : rBestFrame);
			double invSum=1.0/shared.Vector.sum(bestFrame);
//			bsw.println(String.format("FrameStats"/*+(plus ? "+" : "-")*/+":\t%.4f\t%.4f\t%.4f",
//					bestFrame[0]*invSum, bestFrame[1]*invSum, bestFrame[2]*invSum));
		}

		//Results based on terminal poly-As versus poly-Ts
		//TODO: Read 2 should also be used for this analysis... otherwise the insert size interferes...
		if(testPolyA) {
			bsw.println();
			boolean plus=(polyACount>=polyTCount);
			bsw.println("PolyA_Analysis:");
			bsw.println(String.format("MajorStrandPA:\t"+(plus ? "Plus" : "Minus")));
			bsw.println(String.format("PolyA/(PA+PT):\t%.6f", (polyACount/(1.0*(polyACount+polyTCount)))));
			bsw.println(String.format("PolyAFraction:\t%.6f", polyACount*invReads));
		}
		
		//Results based on comparing transcriptome kmers to read kmers
		if(refResults!=null) {
			bsw.println();
			double pmRatio=refResults[0];
			double aFraction=refResults[1];
			double bFraction=refResults[2];
			boolean plus=(pmRatio>=0.5);
			bsw.println("Ref_Analysis:");
			bsw.println("MajorStrandREF:\t"+(plus ? "Plus" : "Minus"));
			bsw.println(String.format("P/(P+M)_Ratio:\t%.6f", pmRatio));
			bsw.println(String.format("GeneCoverage:\t%.4f", aFraction));
			bsw.println(String.format("GenePrecision:\t%.4f", bFraction));
		}
		
//		if(scoreReadGenes) {
//			bsw.println();
//			double readsCalled=(gCallerPlusCount+gCallerMinusCount);
//			double pmRatio=gCallerPlusCount/(1.0*readsCalled);
//			double pScore=gCallerPlusScore/(1.0*readsCalled);
//			double mScore=gCallerMinusScore/(1.0*readsCalled);
//			boolean plus=(pmRatio>=0.5);
//			bsw.println("Kmer_Frequency_Analysis:");
//			bsw.println("MajorStrandKFA:\t"+(plus ? "Plus" : "Minus"));
//			bsw.println(String.format("P/(P+M)_Ratio:\t%.6f", pmRatio));
//			bsw.println(String.format("AvgScorePlus:\t%.6f", pScore));
//			bsw.println(String.format("AvgScoreMinus:\t%.6f", mScore));
//			bsw.println(String.format("UsedFraction:\t%.6f", readsCalled*invReads));
//		}
		
		//Results based on calling and scoring genes on each side of the read
		if(scoreReadGenes) {
			bsw.println();
			double readsCalled=(gCallerPlusCount+gCallerMinusCount);
			double pmRatio=gCallerPlusCount/(1.0*readsCalled);
			double pScore=gCallerPlusScore/(1.0*gCallerPlusCalled);
			double mScore=gCallerMinusScore/(1.0*gCallerMinusCalled);
			boolean plus=(pmRatio>=0.5);
			bsw.println("Read_Gene_Calling_Analysis:");
			bsw.println("MajorStrandRGC:\t"+(plus ? "Plus" : "Minus"));
			bsw.println(String.format("P/(P+M)_Ratio:\t%.6f", pmRatio));
			bsw.println(String.format("AvgScorePlus:\t%.6f", pScore));
			bsw.println(String.format("AvgScoreMinus:\t%.6f", mScore));
			bsw.println(String.format("UsedFraction:\t%.6f", readsCalled*invReads));
		}
		
		//Alignment-based results
		//TODO: Add a flag; should probably be disabled for euks with no gff
		if(ffin1.samOrBam() && (transcriptome || rnaContigs || gff!=null || fna!=null)) {
			long genesAligned=0;
			long genesGT1=0;
			double strandednessSumNorm=0;
			double pmRatioSum=0;
			long plusGeneCount=0;
			double minorSum=0, majorSum=0, maxMinorSum=0, expectedSum=0;
			{//Normalized Statistics
				for(LongPair p : geneMap.values()) {
					if(p.a+p.b>0) {
						float pmRatio=p.a/(float)(p.a+p.b);
						genesAligned++;
						pmRatioSum+=pmRatio;
						plusGeneCount+=(pmRatio>0.5f ? 1 : 0);
						if(p.a+p.b>=minReads) {
							float strandedness=CheckStrand.strandedness(p.a, p.b);
							assert(strandedness>=0 && strandedness<=1.01) : strandedness+", "+p.a+", "+p.b;
							strandednessSumNorm+=strandedness;
							genesGT1++;
							assert(strandednessSumNorm<=genesGT1) : strandednessSumNorm+", "+genesGT1;
							minorSum+=p.min();
							majorSum+=p.max();
							maxMinorSum+=p.sum()/2;
							expectedSum+=CheckStrand.expectedMinorAlleleCount(p.sum());
						}
					}
				}
			}
			final double[] binStrandedness=binStrandedness();
			binStrandednessPerSequence(scaffoldReport);
			
			final boolean plus=(plusAlignedReads>=minusAlignedReads);
			final long readsAligned0=plusAlignedReads+minusAlignedReads;
			final long readsAligned=Tools.max(1, readsAligned0);
			final double alignmentRate=allAlignedReads/(1.0*samLinesProcessed);
			final double alignmentRate2=samLinesAlignedToFeatures/(1.0*samLinesProcessed);
			final double pmRatio=plusAlignedReads/(1.0*readsAligned);
			double strandednessALN=Tools.max(plusAlignedReads, minusAlignedReads)/(1.0*readsAligned);
			if(rnaContigs) {
				//strandednessALN=CheckStrand.strandedness(minorSum, majorSum);
				if(minorSum<=expectedSum) {//Expected case
					strandednessALN=1-0.5f*(minorSum/expectedSum);//0.5-1
				}else{//Rare case; overly unstranded
					double range=maxMinorSum-expectedSum;
					double dif=minorSum-expectedSum;
					strandednessALN=0.5f*(1-(dif/range));
				}
			}
			
			bsw.println();
			bsw.println("Alignment_Results:");
			bsw.println(String.format("StrandednessAL:\t%.2f%%", strandednessALN*100));
			assert(strandednessSumNorm<=genesGT1) : strandednessSumNorm+", "+genesGT1;
			bsw.println(String.format("StrandednessAN:\t%.2f%%", 100*strandednessSumNorm/genesGT1));
//			assert(false) : strandednessSumNorm+", "+genesGT1;
			if(binStrandedness!=null) {
				bsw.println(String.format("StrandednessB:\t%.2f%%", binStrandedness[0]*100));
				bsw.println(String.format("StrandednessBN:\t%.2f%%", binStrandedness[1]*100));
				bsw.println(String.format("StrandednessBS:\t%.2f%%", binStrandedness[2]*100));
				bsw.println(String.format("StrandednessBNS:\t%.2f%%", binStrandedness[3]*100));
			}
			if(!rnaContigs) {
				bsw.println("MajorStrandAL:\t"+(plus ? "Plus" : "Minus"));
				bsw.println(String.format("P/(P+M)_Ratio:\t%.6f", pmRatio));
				bsw.println(String.format("P/(P+M)_RatioN:\t%.6f", pmRatioSum/genesAligned));
			}
			if(!transcriptome) {
				bsw.println(String.format("PlusFeatures:\t%.6f", plusGeneCount/(float)genesAligned));
			}
			bsw.println(String.format("AlignmentRate:\t%.6f", alignmentRate));
			if(!transcriptome && !rnaContigs) {bsw.println(String.format("Feature-Mapped:\t%.6f", alignmentRate2));}
			
//			bsw.println(String.format("Plus-mapped:\t%d", plusAlignedReads));
//			bsw.println(String.format("Minus-mapped:\t%d", minusAlignedReads));
		}
		
		//Displays fraction of reads merged
		if(mergePairs && readsProcessed<readsIn) {
			bsw.println();
			bsw.println(String.format("PairMergeRate:\t%.4f", readsMerged/(1.0*readsProcessed)));
		}
		
		errorState=bsw.poisonAndWait() | errorState;
	}
	
	double[] binStrandedness() {
		if(binMap==null || binMap.isEmpty()) {return null;}
		long binsAligned=0;
		long binsGT1=0;
		double strandednessSumNorm=0;
		double pmRatioSum=0;
		double simpleSumNorm=0;
		long plusBinCount=0;
		double minorSum=0, majorSum=0, maxMinorSum=0, expectedSum=0;
//		int i=0;
		for(Seq seq : binMap.values()) {
			for(LongPair p : seq.bins) {
				if(p!=null && p.a+p.b>0) {
					float pmRatio=p.a/(float)(p.a+p.b);
					binsAligned++;
					pmRatioSum+=pmRatio;
					plusBinCount+=(pmRatio>0.5f ? 1 : 0);
					if(p.a+p.b>=minReads) {
						float strandedness=CheckStrand.strandedness(p.a, p.b);
						assert(strandedness>=0 && strandedness<=1.01) : strandedness+", "+p.a+", "+p.b;
						strandednessSumNorm+=strandedness;
						binsGT1++;
						assert(strandednessSumNorm<=binsGT1) : strandednessSumNorm+", "+binsGT1;
						minorSum+=p.min();
						majorSum+=p.max();
						maxMinorSum+=p.sum()/2;
						expectedSum+=CheckStrand.expectedMinorAlleleCount(p.sum());
						simpleSumNorm+=p.max()/(double)p.sum();
					}
				}
			}
//			if(i++<10) {System.err.println("i="+i+", name="+seq.name+", majorSum="+majorSum);}
		}
		{System.err.println("majorSum="+majorSum);}
		double strandedness=CheckStrand.strandedness((long)minorSum, (long)majorSum, (long)maxMinorSum, (float)expectedSum);
		double strandednessN=strandednessSumNorm/binsGT1;
		double strandednessS=majorSum/(majorSum+minorSum);
		double strandednessNS=simpleSumNorm/binsGT1;
		return new double[] {strandedness, strandednessN, strandednessS, strandednessNS};
	}
	
	void binStrandednessPerSequence(String fname) {
		if(binMap==null || binMap.isEmpty() || fname==null) {return;}
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, append, true);
		bsw.println("#Name\tLength\tReads\tStrandednessB\tStrandednessBN\tStrandednessBS\tStrandednessBNS");
//		int i=0;
//		long majorSum2=0;
		for(Seq seq : binMap.values()) {
			long binsGT1=0;
			double strandednessSumNorm=0;
			double simpleSumNorm=0;
			double minorSum=0, majorSum=0, maxMinorSum=0, expectedSum=0;
			for(LongPair p : seq.bins) {
				if(p!=null && p.a+p.b>0) {
					if(p.a+p.b>=minReads) {
						float strandedness=CheckStrand.strandedness(p.a, p.b);
						assert(strandedness>=0 && strandedness<=1.01) : strandedness+", "+p.a+", "+p.b;
						strandednessSumNorm+=strandedness;
						binsGT1++;
						assert(strandednessSumNorm<=binsGT1) : strandednessSumNorm+", "+binsGT1;
						minorSum+=p.min();
						majorSum+=p.max();
//						majorSum2+=p.max();
						maxMinorSum+=p.sum()/2;
						expectedSum+=CheckStrand.expectedMinorAlleleCount(p.sum());
						simpleSumNorm+=p.max()/(double)p.sum();
					}
				}
			}
//			if(i++<10) {System.err.println("i="+i+", name="+seq.name+", majorSum="+majorSum);}
			double strandedness=0, strandednessN=0, strandednessS=0, strandednessNS=0;
			if(majorSum>0) {
				strandedness=CheckStrand.strandedness((long)minorSum, (long)majorSum, (long)maxMinorSum, (float)expectedSum);
				strandednessN=strandednessSumNorm/binsGT1;
				strandednessS=majorSum/(majorSum+minorSum);
				strandednessNS=simpleSumNorm/binsGT1;
			}
			bsw.print(seq.name).tab().print(seq.len).tab().print((long)(majorSum+minorSum));
			bsw.tab().print(strandedness, 6).tab().print(strandednessN, 6).tab().print(strandednessS, 6).tab().print(strandednessNS, 6).nl();
		}
//		{System.err.println("majorSum="+majorSum2);}
		bsw.poison();
	}
	
	/**
	 * Manages a cached copy of the ref fasta to prevent reading it multiple times.
	 * @return The reference, as Read objects.
	 */
	private ArrayList<Read> genomeSequence(){
		assert(fna!=null);
		if(genomeSequenceCache==null) {
			genomeSequenceCache=ReadInputStream.toReads(fna, FileFormat.FASTA, -1);
		}
		return genomeSequenceCache;
	}
	
	/*--------------------------------------------------------------*/
	
	void incrementGeneMap(String gene, ConcurrentHashMap<String, LongPair> map, int senseStrand, int amt){
		if(gene==null) {return;}
		LongPair pair=map.get(gene);
		if(pair==null) {
			map.putIfAbsent(gene, new LongPair());
			pair=map.get(gene);
		}
		synchronized(pair) {
			if(senseStrand==0) {pair.a+=amt;}
			else {pair.b+=amt;}
		}
	}
	
	void incrementBinMap(String gene, ConcurrentHashMap<String, Seq> map, int pos, int strand, int amt){
		if(gene==null || map==null) {return;}
		Seq seq=map.get(gene);
		if(seq==null) {
			map.putIfAbsent(gene, new Seq(gene, 0));
			seq=map.get(gene);
		}
		int idx=Tools.max(0, pos)/binlen;
		increment(seq.bins, idx, strand, amt);
	}
	
	void increment(ArrayList<LongPair> list, int idx, int strand, int amt) {
		LongPair bin;
		synchronized(list) {
			while(list.size()<=idx) {list.add(null);}
			bin=list.get(idx);
			if(bin==null) {
				bin=new LongPair();
				list.set(idx, bin);
			}
		}
		synchronized(bin) {
			if(strand==0) {bin.a+=amt;}
			else {bin.b+=amt;}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Handles all operations on input reads, 
	 * such as merging, sketching, and stop-codon finding. */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, 
				ConcurrentReadOutputStream crosP_, ConcurrentReadOutputStream crosM_, final int tid_){
			cris=cris_;
			crosP=crosP_;
			crosM=crosM_;
			tid=tid_;//Thread ID
			canonSmm=new SketchMakerMini(canonTool, SketchObject.ONE_SKETCH, minEntropy, minProb, minQual);
			fwdSmm=new SketchMakerMini(fwdTool, SketchObject.ONE_SKETCH, minEntropy, minProb, minQual);
			gCallerT=(scoreReadGenes ? CallGenes.makeGeneCaller(pgm) : null);
//			if(gCallerT!=null) {gCallerT.keepAtLeastOneOrf=true;}//Increases usage rate, but decreases precision
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Fetch lists of reads from the input stream. */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
				if(verbose){outstream.println("Fetched "+ln.size()+" reads.");}
				
				processList(ln);
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
				
				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		/** Iterate through the reads */
		void processList(ListNum<Read> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;
			pReads=(crosP==null ? null : new ArrayList<Read>());
			mReads=(crosM==null ? null : new ArrayList<Read>());
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				final Read r2=r1.mate;
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}
				if(r2!=null && !r2.validated()){r2.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r1.length();
				final int initialLength2=r1.mateLength();

				//Increment counters
				readsInT+=r1.pairCount();
				basesInT+=initialLength1+initialLength2;
				
				processReadPair(r1, r2);
			}
			
			//Output reads to the output streams
			if(crosP!=null){crosP.add(pReads, ln.id);}
			if(crosM!=null){crosM.add(mReads, ln.id);}
			pReads=mReads=null;
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 */
		void processReadPair(Read r1, Read r2){
			if(r1.bases==null) {return;} //Especially for sam lines
			final SamLine sl=r1.samline;
			
			if(sl!=null) {//This means the input was a sam or bam file
				processSamLine(r1, sl);
				if(sl.pairnum()>0) {return;}//Don't process read 2 from sam input
			}
			
			//Merge pairs here if needed.
			if(mergePairs && r2!=null){
				final int insert=BBMerge.findOverlapStrict(r1, r2, false);
				if(insert>0){
					r2.reverseComplement();
					r1=r1.joinRead(insert);
					r2=null;
					readsMergedT++;
				}
			}
			
			readsProcessedT++;
			basesProcessedT+=r1.length();
			r1.mate=null;
			if(doStopAnalysis) {countStopCodons(r1);}
			canonSmm.processReadNucleotide(r1);
			fwdSmm.processReadNucleotide(r1);
			if(testPolyA) {analyzePolyA(r1);}
			if(scoreReadGenes) {scoreGenes2(r1);}
		}
		
		/**
		 * Counts stop codons in each frame, and the longest ORF per frame.
		 */
		void countStopCodons(Read r) {
			final byte[] bases=r.bases;
			if(bases.length<5) {return;}
			
			final int k=3;
			final int shift=2*k;
			final int shift2=shift-2;
			final int mask=0b111111;
			
			int kmer=0;
			int rkmer=0;
			int len=0;

			Arrays.fill(fCount, 0);
			Arrays.fill(rCount, 0);
			Arrays.fill(fPos, bases.length);
			Arrays.fill(rPos, bases.length);
			
			Arrays.fill(fMaxOrf, -1);
			Arrays.fill(rMaxOrf, -1);
			Arrays.fill(fLastStop, 0);
			Arrays.fill(rLastStop, 0);
			
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				int x=AminoAcid.baseToNumber[b];
				int x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(x<0){
					len=0;
					rkmer=0;
				}else{
					len++;
					ACGTCountT[x]++;
				}
				if(len>=k){
					final int j=i%3;
					final byte aa=AminoAcid.codeToByte[kmer];
					final byte raa=AminoAcid.codeToByte[rkmer];
					if(aa=='*'){
//						System.err.println("aa="+Character.toString(aa));
						fCount[j]++;
						fPos[j]=Tools.min(fPos[j], i);//TODO: This only tracks the first ORF, not the longest
						
						final int orfLen=i-fLastStop[j];
						fLastStop[j]=i;
						fMaxOrf[j]=Tools.max(fMaxOrf[j], orfLen);
					}
					if(raa=='*') {
						final int i2=bases.length-i-2;
						rCount[j]++;
						rPos[j]=Tools.min(rPos[j], i2);
						
						final int orfLen=rLastStop[j]-i2;
						rLastStop[j]=i2;
						rMaxOrf[j]=Tools.max(rMaxOrf[j], orfLen);
//						System.err.println(i2+", "+rLastStop[j]+", "+rMaxOrf[j]);
					}
				}
			}
			fCountSumT+=Tools.min(fCount);
			rCountSumT+=Tools.min(rCount);
			
			for(int i=0; i<3; i++) {
				if(fMaxOrf[i]<0) {fMaxOrf[i]=bases.length;}
				if(rMaxOrf[i]<0) {rMaxOrf[i]=bases.length;}
			}
			
			final int[] fORF=(useFirstORF ? fPos : fMaxOrf);
			final int[] rORF=(useFirstORF ? rPos : rMaxOrf);
			final int mfp=Tools.max(fORF);
			final int mrp=Tools.max(rORF);
			fPosSumT+=mfp;
			rPosSumT+=mrp;

//			fBestFrameT[Tools.maxIndex(fORF)]++;
//			rBestFrameT[Tools.maxIndex(rORF)]++;
			for(int i=0; i<3; i++) {//This is better since it increments all equal ones
				if(fORF[i]==mfp){fBestFrameT[i]++;}
				if(rORF[i]==mrp){rBestFrameT[i]++;}
			}
//			assert(false) : "\n"+Arrays.toString(fCount)+"\n"+Arrays.toString(rCount)+"\n"+Arrays.toString(fPos)+"\n"+Arrays.toString(rPos);
		}
		
		/** 
		 * Counts tip poly-As and poly-Ts.
		 * TODO: Should take into account pairnum so read 2 can be used.
		 */
		void analyzePolyA(Read r) {
			final int trailingPolyA=r.countTrailing('A'), leadingPolyT=r.countLeading('T');
//			final int trailingPolyT, leadingPolyA;//These two should be spurious, but could be used to infer whether the poly-A is due to real poly-A tails.
			
			if(trailingPolyA>leadingPolyT && trailingPolyA>minPolyA) {polyACountT++;}
			else if(trailingPolyA<leadingPolyT && leadingPolyT>minPolyA) {polyTCountT++;}
			
		}
		
		/**
		 * Calls genes on both strands, and records which strand
		 * had the highest-scoring gene.
		 */
		void scoreGenes2(Read r) {
			ArrayList<Orf> list;
			list=gCallerT.callGenes(r, true);
//			System.err.print((list==null ? 0 : list.size())+", ");
//			Orf bestPlus=null;
//			Orf bestMinus=null;
			final float min=-99999;
			float plusScore=min, minusScore=min;
			if(list!=null && !list.isEmpty()) {
				for(Orf orf : list) {
//					if(orf.strand==0) {
//						if(bestPlus==null || bestPlus.score()<orf.score()) {bestPlus=orf;}
//					}else {
//						if(bestMinus==null || bestMinus.score()<orf.score()) {bestMinus=orf;}
//					}
					if(orf.strand==0) {
						plusScore=Tools.max(plusScore, orf.score());
					}else {
						minusScore=Tools.max(minusScore, orf.score());
					}
				}
			}

			//This block is to ensure the uncalled strand gets a lower score,
			//in case the called strand has a negative score
			if(plusScore==min && minusScore==min) {
				plusScore=minusScore=0;
			}else if(plusScore==min) {
				plusScore=Tools.min(0, minusScore*1.5f-1);
			}else if(minusScore==min) {
				minusScore=Tools.min(0, plusScore*1.5f-1);
			}
			
//			if(plusScore<=0 && minusScore<=0) {return;}
			
			gCallerPlusScoreT+=plusScore;
			gCallerMinusScoreT+=minusScore;
			if(plusScore>minusScore) {gCallerPlusCountT++;}
			else if(minusScore>plusScore){gCallerMinusCountT++;}
			if(plusScore>0){gCallerPlusCalledT++;}
			if(minusScore>0){gCallerMinusCalledT++;}
			
		}
		
		@Deprecated
		/** 
		 * Old version; just looked at enriched interior kmers instead of
		 * doing normal gene-calling.  Didn't work very will.
		 */
		void scoreGenes(Read r) {
			byte[] bases=r.bases;
			double plusScoreCDS=gCaller.scoreFeature(bases, ProkObject.CDS);//These scores are suspiciously low; I wonder if frame tracking is working correctly?
			double plusScore16S=gCaller.scoreFeature(bases, ProkObject.r16S);
			double plusScore5S=gCaller.scoreFeature(bases, ProkObject.r5S);
			AminoAcid.reverseComplementBasesInPlace(bases);
			double minusScoreCDS=gCaller.scoreFeature(bases, ProkObject.CDS);
			double minusScore16S=gCaller.scoreFeature(bases, ProkObject.r16S);
			double minusScore5S=gCaller.scoreFeature(bases, ProkObject.r5S);
			AminoAcid.reverseComplementBasesInPlace(bases);

			double maxCDS=Tools.max(plusScoreCDS, minusScoreCDS)*1.6;
			double max16S=Tools.max(plusScore16S, minusScore16S)*1.2;
			double max5S=Tools.max(plusScore5S, minusScore5S)*0.9;
			boolean cds=maxCDS>=max16S;
			boolean r5s=max5S>=max16S && max5S>=maxCDS;
			double plusScore=r5s ? plusScore5S : cds ? plusScoreCDS : plusScore16S;
			double minusScore=r5s ? minusScore5S : cds ? minusScoreCDS : minusScore16S;
			
			//Doesn't fit the model well...
			if(plusScore<0.97 && minusScore<0.97) {return;}
			
			gCallerPlusScoreT+=plusScore;
			gCallerMinusScoreT+=minusScore;
			if(plusScore>=minusScore) {gCallerPlusCountT++;}
			else {gCallerMinusCountT++;}
		}
		
		/**
		 * Associate an aligned read with the plus or minus strand.
		 * @param r1 The read.
		 * @param sl The read's SamLine.
		 */
		void processSamLine(final Read r1, final SamLine sl) {
			if(!sl.primary() || sl.supplementary()) {return;}
			samLinesProcessedT++;
			if(!sl.mapped()) {return;}
			allAlignedReadsT++;
			
			//Flip read 2's associated strand.
			final int mappedStrand=(sl.strand()^sl.pairnum());
			
			//The strand on which the gene is located.
			int senseStrand=-1;
			String gene=null;
			if(transcriptome) {
				senseStrand=mappedStrand;
				gene=sl.rnameS();
			}else if(rnaContigs) {
				senseStrand=mappedStrand;
				gene=sl.rnameS();
			}else if(gffLines!=null){//Try to figure out the sense strand
				final int start=sl.start(false, false);
//				final int stop=sl.stop(start, false, false);
//				final int p=(start+stop)/2;
				final int p=start;//A point used as proxy for the read location
				
				//Get the ranges for the contig
				String rname=sl.rnameS();
				Range[] ranges=rangeMap.get(rname);
				
				//Retry the name, trimmed to the first whitespace
				if(ranges==null) {ranges=rangeMap.get(sl.rnamePrefix());}
//				assert(ranges!=null) : "Can't find rname "+sl.rnameS()+":\n"+sl+"\n";// in "+rangeMap.keySet()+"\n"+sl;
				
				if(ranges!=null) {
					//Find the range containing the proxy point
					int idx=Range.findIndex(p, ranges);
					if(idx<0) {idx=-idx-1;}
					if(idx<0 || idx>=ranges.length){return;}
					final Range r=ranges[idx];

					//List of features overlapping the point
					final ArrayList<GffLine> lines=(ArrayList<GffLine>) r.obj1;
					for(int i=lines.size()-1; i>=0; i--) {//Now to select a line...
						final GffLine line=lines.get(i);
						//TODO: With multiple overlapping features, this may or may not be the best one
						if(r.intersects(line.start, line.stop)) {
							senseStrand=mappedStrand^line.strand;
							gene=sl.rnameS();
							break;
						}
					}
				}
			}
			
			if(gene!=null) {
				incrementGeneMap(gene, geneMap, senseStrand, 1);
				incrementBinMap(gene, binMap, sl.pos, mappedStrand, 1);

				//Track results according to the determined strand
				if(senseStrand==0) {
					samLinesAlignedToFeaturesT++;
					plusAlignedReadsT++;
					if(crosP!=null) {pReads.add(r1);}
				}else if(senseStrand==1){
					samLinesAlignedToFeaturesT++;
					minusAlignedReadsT++;
					if(crosM!=null) {mReads.add(r1);}
				}else {
					assert(senseStrand==-1);
				}
			}
			
		}

		/** Number of reads in to this thread */
		protected long readsInT=0;
		/** Number of bases in to this thread */
		protected long basesInT=0;

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;

		/** Number of reads merged by this thread */
		protected long readsMergedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		
		private final ConcurrentReadOutputStream crosP;
		private final ConcurrentReadOutputStream crosM;
		
		private ArrayList<Read> pReads, mReads;

		final SketchMakerMini canonSmm;
		final SketchMakerMini fwdSmm;
		
		//These are just buffers
		/** Counts of stop codons on forward strand */
		private final int[] fCount=new int[3];
		/** Counts of stop codons on reverse strand */
		private final int[] rCount=new int[3];
		/** Longest initial ORF on forward strand */
		private final int[] fPos=new int[3];
		/** Longest initial ORF on reverse strand */
		private final int[] rPos=new int[3];
		
		/** Longest ORF on forward strand */
		private final int[] fMaxOrf=new int[3];
		/** Longest ORF on reverse strand */
		private final int[] rMaxOrf=new int[3];
		/** Last stop seen on forward strand */
		private final int[] fLastStop=new int[3];
		/** Last stop seen on reverse strand */
		private final int[] rLastStop=new int[3];

		private long fCountSumT;
		private long rCountSumT;
		private long fPosSumT;
		private long rPosSumT;

		private final long[] fBestFrameT=new long[3];
		private final long[] rBestFrameT=new long[3];
		private final long[] ACGTCountT=new long[4];
		
		private long polyACountT=0;
		private long polyTCountT=0;

		private long gCallerPlusCountT=0;
		private long gCallerMinusCountT=0;
		private long gCallerPlusCalledT=0;
		private long gCallerMinusCalledT=0;
		private double gCallerPlusScoreT=0;
		private double gCallerMinusScoreT=0;
		
		private long plusAlignedReadsT=0;
		private long minusAlignedReadsT=0;
		private long allAlignedReadsT=0;

		private long samLinesProcessedT=0;
		private long samLinesAlignedToFeaturesT=0;
		
		private GeneCaller gCallerT;
		
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	
	static class Seq {
		
		public Seq(String name_, int len_) {
			name=name_;
			len=len_;
		}
		
		String name;
		int len;
		ArrayList<LongPair> bins=new ArrayList<LongPair>(1);
		
		public static ConcurrentHashMap<String, Seq> loadSamHeader(FileFormat ff){
			ConcurrentHashMap<String, Seq> map=new ConcurrentHashMap<String, Seq>();
			ByteFile bf=ByteFile.makeByteFile(ff);
			LineParser4 lp=new LineParser4("\t:\t:");
			byte[] line=bf.nextLine();
			while(line!=null && line.length>0){
				lp.set(line);
				if(lp.startsWith("@SQ\t")){
					Seq seq=new Seq(lp.parseString(2), lp.parseInt(4));
					map.put(seq.name, seq);
				}else if(line[0]!='@'){
					break;
				}
				line=bf.nextLine();
			}
			bf.close();

			return map;
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	String fna=null;
	String gff=null;
	String pgmFile=null;
	private String types="CDS,rRNA,tRNA,ncRNA,exon,5S,16S,23S";
	private ArrayList<GffLine> gffLines=null;
	
	private ConcurrentHashMap<String, Seq> binMap;
	
	//TODO:
	//scaffold_name  scaffold_length  scaffold_mapped_reads  scaffold_strandedness
	
	/** 
	 * Map associating with contig names with ordered arrays of nonoverlapping ranges
	 * Each range has attached a list of features that fully contain it. */
	private HashMap<String, Range[]> rangeMap;

	/** Primary output file path */
	private String out1="stdout.txt";
	
	/** Output file path for plus-aligned reads */
	private String outPlus=null;
	
	/** Output file path for minus-aligned reads */
	private String outMinus=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	private String scaffoldReport=null;
	
	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;
	
	private int binlen=0;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsIn=0;
	/** Number of bases processed */
	protected long basesIn=0;
	
	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	
	protected long readsMerged=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;

	private float samplerate=1;
	private long sampleseed=17;
	private int sketchSize=80000;
	
	private float minEntropy=0;
	private float minProb=0;
	private byte minQual=0;
	private boolean mergePairs=false;
	private final boolean doDepthAnalysis=true;
	private boolean doStopAnalysis=true;
	private boolean testPolyA=true;
	
	private ArrayList<Read> genes;
	
	private SketchTool canonTool;
	private SketchTool fwdTool;
	
	private SketchHeap canonHeap=null;
	private SketchHeap fwdHeap=null;
	
	private Sketch canonSketch=null;
	private Sketch fwdSketch=null;

	private Sketch canonGeneSketch=null;
	private Sketch plusSketch=null;
	private Sketch minusSketch=null;

	private long fCountSum=0;
	private long rCountSum=0;
	private long fPosSum=0;
	private long rPosSum=0;

	private boolean useFirstORF=false;
	private final long[] fBestFrame=new long[3];
	private final long[] rBestFrame=new long[3];
	private final long[] ACGTCount=new long[4];
	
	private long polyACount=0;
	private long polyTCount=0;
	private int minPolyA=6;

	private long gCallerPlusCount=0;
	private long gCallerMinusCount=0;
	private long gCallerPlusCalled=0;
	private long gCallerMinusCalled=0;
	private double gCallerPlusScore=0;
	private double gCallerMinusScore=0;
	
	private long plusAlignedReads=0;
	private long minusAlignedReads=0;
	private long allAlignedReads=0;
	private long samLinesProcessed=0;
	private long samLinesAlignedToFeatures=0;
	
	//Possibly should change this to scoring the best orf per strand
	//This gave the wrong answer for a metatranscriptome
	private boolean scoreReadGenes=true;
	private final GeneModel pgm;
	private final GeneCaller gCaller;
	private int passes=2;
	private boolean streamGenome=false;
	private ArrayList<Read> genomeSequenceCache=null;

	ConcurrentHashMap<String, LongPair> geneMap=new ConcurrentHashMap<String, LongPair>();
	
	/** Whether the input reference is a transcriptome or genome */
	private boolean transcriptome=false;
	/** True for unknown sense contigs from RNA-seq data */
	private boolean rnaContigs=false;
	
	private boolean normalize=false;
	private int minReads=2;//Don't calculate strandedness from fewer reads
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output file */
	private final FileFormat ffout1;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
