package bbduk;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import aligner.SideChannel3;
import cardinality.CardinalityTracker;
import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.CalcTrueQuality;
import json.JsonObject;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import stream.FASTQ;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Separates, trims, or masks sequences based on matching kmers in a reference.
 * Supports Hamming and and edit distance.
 * Supports K 1-31 and emulated K>31.
 * @author Brian Bushnell
 * @date Aug 30, 2013
 *
 */
public class BBDukS {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Create a new BBDuk instance
		BBDukS x=new BBDukS(args);
		
		//And run it
		x.process();
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
		assert(WAYS==7 && !verbose) : "Undo WAYS and verbose";
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public BBDukS(String[] args){
		
		/* Parse arguments */
		
		BBDukParser p=new BBDukParser(args, getClass());
		parser=p;
		
		silent=p.silent;
		json=p.json;
		khistIn=p.khistIn;
		khistOut=p.khistOut;
		ref=p.ref;
		literal=p.literal;
		in1=p.in1;
		in2=p.in2;
		qfin1=p.qfin1;
		qfin2=p.qfin2;
		qfout1=p.qfout1;
		qfout2=p.qfout2;
		out1=p.out1;
		out2=p.out2;
		outb1=p.outb1;
		outb2=p.outb2;
		outsingle=p.outsingle;
		outstats=p.outstats;
		outrqc=p.outrqc;
		outrpkm=p.outrpkm;
		outrefstats=p.outrefstats;
		polymerStatsFile=p.polymerStatsFile;
		tossJunk=p.tossJunk;
		maxReads=p.maxReads;
		samplerate=p.samplerate;
		sampleseed=p.sampleseed;
		ordered=p.ordered;
		samFile=p.samFile;
		jsonStats=p.jsonStats;
		
		recalibrateQuality=p.recalibrateQuality;
		skipreads=p.skipreads;	
		threadsIn=p.threadsIn;
		threadsOut=p.threadsOut;
		
		outstream=BBDukParser.outstream;
		overwrite=BBDukParser.overwrite;
		append=BBDukParser.append;
		showSpeed=BBDukParser.showSpeed;
		DISPLAY_PROGRESS=BBDukParser.DISPLAY_PROGRESS;
		THREADS=BBDukParser.THREADS;
		STATS_COLUMNS=BBDukParser.STATS_COLUMNS;
		REPLICATE_AMBIGUOUS=BBDukParser.REPLICATE_AMBIGUOUS;

//		WAYS=p.WAYS;
//		initialSizeDefault=p.initialSizeDefault;
		
		loglogIn=(p.loglog ? CardinalityTracker.makeTracker(p.parser) : null);
		loglogOut=(p.loglogOut ? CardinalityTracker.makeTracker(p.parser) : null);
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, qfout1, qfout2, outb1, outb2, outsingle, outstats, 
				outrpkm, outrqc, outrefstats, polymerStatsFile, khistIn, khistOut, p.alignOut)){
			throw new RuntimeException("\nCan't write to some output files; overwrite="+overwrite+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2, qfin1, qfin2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(!Tools.testInputFiles(true, true, ref)){
			throw new RuntimeException("\nCan't read to some reference files.\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, qfin1, qfin2, qfout1, qfout2,
				out1, out2, outb1, outb2, outsingle, outstats, outrpkm, outrqc, outrefstats, polymerStatsFile, khistIn, khistOut, p.alignOut)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		assert(THREADS>0) : "THREADS must be greater than 0.";

		assert(in1==null || in1.toLowerCase().startsWith("stdin") || in1.toLowerCase().startsWith("standardin") || new File(in1).exists()) : "Can't find "+in1;
		assert(in2==null || in2.toLowerCase().startsWith("stdin") || in2.toLowerCase().startsWith("standardin") || new File(in2).exists()) : "Can't find "+in2;
		
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
		
		final int defaultFormat=(ffin1==null ? FileFormat.FASTQ : ffin1.format());
		ffout1=FileFormat.testOutput(out1, defaultFormat, null, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, defaultFormat, null, true, overwrite, append, ordered);
		ffoutb1=FileFormat.testOutput(outb1, defaultFormat, null, true, overwrite, append, ordered);
		ffoutb2=FileFormat.testOutput(outb2, defaultFormat, null, true, overwrite, append, ordered);
		ffouts=FileFormat.testOutput(outsingle, defaultFormat, null, true, overwrite, append, ordered);

		p.parser.validateStdio(ffin1, ffout1, ffoutb1, ffouts);
		
		index=new BBDukIndex(p);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Main processing method that coordinates the entire workflow.
	 * Loads variants if specified, initializes quality recalibration,
	 * calls process2 for core processing, and calculates final statistics.
	 */
	public void process(){
		
		{//TODO: State sync can be mostly eliminated once Index is used for complete encapsulation 
			//1. Load Index
			index.loadRef(in1);

			//2. Sync State - References and Arrays
			sidechannel=index.sidechannel;
			ref=index.ref;
			scaffoldNames=index.scaffoldNames;
			
			errorState|=index.errorState;  // Merge error states
		}
	    
		if(recalibrateQuality){
			if(samFile!=null){
				CalcTrueQuality.main2(new String[] {"in="+samFile, "showstats=f"});
			}
			CalcTrueQuality.initializeMatrices();
		}
		
		/* Check for output file collisions */
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outb1, outb2, outstats, outrpkm, outrqc, outrefstats)){
			throw new RuntimeException("One or more output files were duplicate or could not be written to.  Check the names or set the 'overwrite=true' flag.");
		}
		
		process2();
		
		if(outstream!=System.err && outstream!=System.out){outstream.close();}
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	/**
	 * Core processing method that loads reference kmers and processes input reads.
	 * Fills kmer tables from reference sequences, then spawns threads to match
	 * reads against reference kmers and perform filtering/trimming operations.
	 */
	public void process2(){
		
		/* Start phase timer */
		Timer t0=new Timer();
		Timer t1=new Timer();
		final long start=t0.time1;

		if(DISPLAY_PROGRESS && !json){
			outstream.println("Initial:");
			Shared.printMemory(outstream);
			outstream.println();
		}
		
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=THREADS<4;
		
		/* Do kmer matching of input reads */
		BBDukProcessorS proc=spawnProcessThreads(parser, t1);
		
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		proc.printOutput(start);
		
		/* Stop timer and calculate speed statistics */
		t0.stop();
		t1.stop();
		lastReadsOut=proc.readsOut;
		
		if(showSpeed && !json){
			outstream.println();
			outstream.println(Tools.timeReadsBasesProcessed(t0, proc.readsIn, proc.basesIn, 8));
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	

	/**
	 * Match reads against reference kmers, using multiple ProcessThread.
	 * @param t
	 */
	private BBDukProcessorS spawnProcessThreads(BBDukParser p, Timer t){
		t.start();
		
		/* Create read input stream */
		final Streamer cris;
		final boolean paired;
		{
			cris=StreamerFactory.getReadInputStream(maxReads, ffin1.samOrBam(), 
				ffin1, ffin2, qfin1, qfin2, threadsIn);
			cris.setSampleRate(samplerate, sampleseed);
			cris.start();
			paired=cris.paired();
			assert(ffin2!=null || !paired || (FASTQ.FORCE_INTERLEAVED || FASTQ.TEST_INTERLEAVED)) : 
				paired+", "+FASTQ.FORCE_INTERLEAVED+", "+FASTQ.TEST_INTERLEAVED+", "+ffin1.interleaved();
			if(!ffin1.samOrBam() && !silent){
				if(json){
					jsonStats.add("paired", paired);
				}else{
					outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));
				}
			}
		}
		
		/* Create read output streams */
		final Writer ros, rosb, ross;
		{
			final int buff=(!ordered ? 12 : Tools.max(32, 2*Shared.threads()));
			if(out1!=null){
				ros=WriterFactory.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, true, threadsOut);
				ros.start();
			}else{ros=null;}
			if(outb1!=null){
				rosb=WriterFactory.getStream(ffoutb1, ffoutb2, null, null, buff, null, true, threadsOut);
				rosb.start();
			}else{rosb=null;}
			if(outsingle!=null){
				ross=WriterFactory.getStream(ffouts, null, null, null, buff, null, true, threadsOut);
				ross.start();
			}else{ross=null;}
			if(sidechannel!=null) {sidechannel.start();}
			if(ros!=null || rosb!=null || ross!=null || sidechannel!=null){
				t.stop();
				if(!silent && !json){outstream.println("Started output streams:\t"+t);}
				t.start();
			}
		}
		
		/* Optionally skip the first reads, since initial reads may have lower quality */
		if(skipreads>0){
			long skipped=0;

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(skipped<skipreads && ln!=null && !ln.poison()){
				skipped+=reads.size();
				
				if(rosb!=null){rosb.add(new ArrayList<Read>(1), ln.id);}
				if(ros!=null){ros.add(new ArrayList<Read>(1), ln.id);}
				if(ross!=null){ross.add(new ArrayList<Read>(1), ln.id);}
				
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln==null || ln.poison() || ln.last()){
				ReadWrite.closeStreams(cris, ros, rosb, ross);
				if(sidechannel!=null) {sidechannel.shutdown();}
				outstream.println("Skipped all of the reads.");
				System.exit(0);
			}
		}
		
		BBDukProcessorS processor=new BBDukProcessorS(p, index, cris, ros, rosb, ross);
		
		/* Create ProcessThreads */
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new ProcessThread(processor));}
		for(ProcessThread pt : alpt){pt.start();}
		
		/* Wait for threads to die, and gather statistics */
		for(ProcessThread pt : alpt){
			
			/* Wait for a thread to die */
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			
			/* Accumulate data from per-thread counters */
			processor.add(pt.processor);
			errorState|=(!pt.processor.finishedSuccessfully);
		}
		
		/* Shut down I/O streams; capture error status */
		{
			//Prevent a spurious error message in the event of a race condition when maxReads is set.
			boolean b=ReadWrite.closeStream(cris);
			if(maxReads<1 || maxReads==Long.MAX_VALUE || (maxReads!=processor.readsIn && maxReads*2!=processor.readsIn && samplerate<1)){errorState|=b;}
		}
		errorState|=ReadWrite.closeOutputStreams(ros, rosb, ross);
		if(sidechannel!=null) {errorState|=sidechannel.shutdown();}
		errorState|=ReadStats.writeAll();
		
		t.stop();
		if(showSpeed && !json){
			outstream.println("Processing time:   \t\t"+t);
		}
		return processor;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Matches read kmers against reference kmers, performs binning and/or trimming, and writes output.
	 */
	private class ProcessThread extends Thread{
		
		/** Constructor */
		public ProcessThread(BBDukProcessorS processor_){
			processor=processor_.clone();
		}
		
		@Override
		public void run(){
			processor.process();
		}
		
		final BBDukProcessorS processor;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses literal sequence argument into array of sequences.
	 * @param arg Comma-separated list of literal sequences
	 * @return Array of processed literal sequences
	 */
	public static final String[] processLiteralArg(String arg) {
		if(arg==null) {return null;}
		String[] split=arg.split(",");
		ArrayList<String> list=new ArrayList<String>(split.length);
		for(String b : split) {
			String c=processLiteralTerm(b);
			if(c!=null) {list.add(c);}
		}
		String[] ret=list.isEmpty() ? null : list.toArray(new String[0]);
		return ret;
	}
	
	/**
	 * Processes individual literal sequence term, expanding polymer shortcuts.
	 * @param b Single literal sequence or polymer specification
	 * @return Processed literal sequence
	 */
	public static final String processLiteralTerm(String b) {
		assert(b.length()>0) : "Invalid literal sequence: '"+b+"'";
		if(AminoAcid.isACGTN(b)) {
			return b;
		}
		b=b.toUpperCase().replaceAll("-", "");
		if(b.startsWith("POLY")) {
			b=b.replace("POLY", "");
			assert(AminoAcid.isACGTN(b)) : "Invalid literal sequence: '"+b+"'";
			StringBuilder sb=new StringBuilder(40);
			final int minlen=Tools.max(31+b.length(), b.length()*3);
			while(sb.length()<minlen) {sb.append(b);}
			System.err.println("Adding literal polymer "+sb);
			return sb.toString();
		}else {
			assert(false) : "Invalid literal sequence: '"+b+"'";
		}
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Has this class encountered errors while processing? */
	public boolean errorState;
	
	/** Stores JSON output */
	private final JsonObject jsonStats;
	
	/*--------------------------------------------------------------*/
	
	private final BBDukParser parser;
	private final BBDukIndex index;
	
	/** A scaffold's name is stored at scaffoldNames.get(id).
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	private ArrayList<String> scaffoldNames;
	/** Array of reference files from which to load kmers */
	private String[] ref;
	/** Array of literal strings from which to load kmers */
	private final String[] literal;
	
	/*--------------------------------------------------------------*/
	/*----------------          Immutable           ----------------*/
	/*--------------------------------------------------------------*/

	//TODO: Document
	private final boolean silent;
	private final boolean json;
	
	/** For calculating kmer cardinality in input */
	private final CardinalityTracker loglogIn;
	/** For calculating kmer cardinality in output */
	private final CardinalityTracker loglogOut;
	/** Requires (and sets) cardinality tracking.  This is for input kmers. */
	private final String khistIn;
	/** Requires (and sets) cardinality tracking.  This is for output kmers. */
	private final String khistOut;

	/** Input reads */
	private final String in1, in2;
	/** Input FileFormats */
	private final FileFormat ffin1, ffin2;
	/** Input qual files */
	private final String qfin1, qfin2;
	/** Output qual files */
	private final String qfout1, qfout2;
	/** Output reads (unmatched and at least minlen) */
	private final String out1, out2;
	/** Output reads (matched or shorter than minlen) */
	private final String outb1, outb2;
	/** Output FileFormats */
	private final FileFormat ffout1, ffout2, ffoutb1, ffoutb2, ffouts;
	/** Output reads whose mate was discarded */
	private final String outsingle;
	/** Statistics output files */
	private final String outstats, outrqc, outrpkm, outrefstats, polymerStatsFile;	
	
	private int threadsIn=-1;
	private int threadsOut=-1;
	
	final boolean tossJunk;
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	private final long maxReads;
	/** Process this fraction of input reads. */
	private final float samplerate;
	/** Set samplerate seed to this value. */
	private final long sampleseed;
	
	/** Output reads in input order.  May reduce speed. */
	private final boolean ordered;
	
	/*--------------------------------------------------------------*/
	/*----------------       Variant-Related        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Optional file for quality score recalibration */
	private final String samFile;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Skip this many initial input reads */
	private final long skipreads;
	/** Recalibrate quality scores using matrices */
	private final boolean recalibrateQuality;
	
	/*--------------------------------------------------------------*/
	/*----------------         Side Channel         ----------------*/
	/*--------------------------------------------------------------*/
	
	private SideChannel3 sidechannel;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of tables (and threads, during loading) */
	private static final int WAYS=7;
	/** Verbose messages */
	public static final boolean verbose=false;
	
	/** Number of reads output in the last run */
	public static long lastReadsOut;
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Print speed statistics upon completion */
	public static boolean showSpeed=true;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.threads();
	/** Indicates end of input stream */
	static final ArrayList<Read> POISON=new ArrayList<Read>(0);
	/** Number of columns for statistics output, 3 or 5 */
	public static int STATS_COLUMNS=3;
	/** Release memory used by kmer storage after processing reads */
	public static boolean RELEASE_TABLES=true;
	/** Make unambiguous copies of ref sequences with ambiguous bases */
	public static boolean REPLICATE_AMBIGUOUS=false;
	
	/** Stores some data for statistics when running RQCFilter; not used otherwise. */
	public static HashMap<String, Long> RQC_MAP=null;
	
}
