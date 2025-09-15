package bin;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import aligner.Factory;
import aligner.IDAligner;
import aligner.SingleStateAlignerFlat2;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import prok.GeneCaller;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Command-line tool for searching genomic sequences against a reference database.
 * Identifies the closest taxonomic match for query genomes using k-mer profiling.
 * Supports multithreaded processing for efficient searching of large databases.
 * 
 * @author Brian Bushnell
 * @date April 12, 2025
 */
public class CladeSearcher extends BinObject implements Accumulator<CladeSearcher.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();

		//Set command-line entrance defaults
		GeneCaller.useIDAligner=true;
		
		//Create an instance of this class
		CladeSearcher x=new CladeSearcher(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Default constructor.
	 * Ensures the entropy adjustment model is loaded.
	 */
	public CladeSearcher(){
		if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150) {
			AdjustEntropy.load(4, 150);
		}
	}
	
	/**
	 * Constructor with command line arguments.
	 * Parses arguments, sets up parameters, and validates input/output paths.
	 * 
	 * @param args Command line arguments
	 */
	public CladeSearcher(String[] args){
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
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			extin=parser.extin;

			if(parser.setOut) {out=parser.out1;}
			extout=parser.extout;
		}

		if(ref.isEmpty()) {
			String s=defaultRef();
			if(s!=null) {ref.add(s);}
		}
		CladeIndex.heapSize=Math.max(CladeIndex.heapSize, maxHitsToPrint);
		validateParams();
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, extout, true, overwrite, append, ordered);
	}
	
	/**
	 * Returns the default reference database path based on the execution environment.
	 * 
	 * @return Path to the default reference database, or null if unavailable
	 */
	public static String defaultRef() {
		if(Shared.DORI) {
			return ("/clusterfs/jgi/groups/gentech/homes/bbushnell/clade/refseqA48_with_ribo.spectra.gz");
		}else if(Shared.PERLMUTTER) {
			return ("/global/cfs/cdirs/bbtools/clade/refseqA48_with_ribo.spectra.gz");
		}else {
			return Data.findPath("?refseqA48_with_ribo.spectra.gz");
		}
	}
	
	/**
	 * Sets up the searcher by loading necessary components.
	 * Loads entropy model and taxonomy tree if needed.
	 */
	void setup() {
		if(calcCladeEntropy && (AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150)) {
			AdjustEntropy.load(4, 150);
		}
		
		if(useTree) {BinObject.loadTree();}
		if(Clade.callSSU) {
			GeneTools.loadPGM();
			GeneCaller.call23S=GeneCaller.call5S=GeneCaller.calltRNA=GeneCaller.callCDS=false;
		}
	}
	
	/**
	 * Loads and indexes reference clades from the specified files.
	 */
	void loadIndex() {
		Timer t=new Timer(outstream, false);
		index=CladeIndex.loadIndex(ref);
		t.stop("Indexed "+index.size()+" spectra in ");
	}
	
	/**
	 * Loads query clades from input files.
	 * Can process individual contigs separately if perContig is true.
	 */
	void loadQueries() {
		Timer t=new Timer(outstream, false);
		CladeLoaderMF loaderMF=new CladeLoaderMF();
		queries=loaderMF.loadFiles(in, perContig, minContig);
		t.stopAndStart("Loaded "+queries.size()+" queries in ");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Parses command line arguments.
	 * Handles program-specific parameters and delegates standard flags to Parser.
	 * 
	 * @param args Command line arguments
	 * @return Parser with parsed standard flags
	 */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("multithreaded") || a.equalsIgnoreCase("maxCompareThreads") || 
					a.equals("comparethreads")){
				if(Tools.isNumeric(b)) {
					maxCompareThreads=Integer.parseInt(b);
					multithreaded=maxCompareThreads>1;
				}else {
					multithreaded=Parse.parseBoolean(b);
				}
			}else if(a.equals("hits") || a.equals("maxhits") || a.equals("records")){
				maxHitsToPrint=Integer.parseInt(b);
			}else if(a.equals("percontig") || a.equals("persequence")){
				perContig=Parse.parseBoolean(b);
			}else if(a.equals("format")){
				assert(b!=null) : "'format' flag requires and option, like 'format=oneline'";
				if(Tools.isNumeric(b)) {
					format=Integer.parseInt(b);
				}else if(b.equals("machine") || b.equals("oneline")){
					format=MACHINE;
				}else if(b.equals("human")){
					format=HUMAN;
				}else {
					assert(false) : "Unknown format "+b;
				}
			}else if(a.equals("machine") || a.equals("machineout") || a.equals("oneline")){
				if(Parse.parseBoolean(b)) {format=MACHINE;}
			}else if(a.equals("qtid") || a.equals("printqtid")){
				printQTID=Parse.parseBoolean(b);
			}else if(a.equals("mincontig") || a.equals("minlen") || a.equals("minlength")){
				minContig=Integer.parseInt(b);
			}else if(a.equals("printmetrics")){
				printMetrics=Parse.parseBoolean(b);
			}else if(a.equals("usetree")){
				useTree=Parse.parseBoolean(b);
			}else if(a.equals("ref")){
				Tools.getFileOrFiles(b, ref, true, false, false, false);
			}else if(a.equals("in")){
				Tools.getFileOrFiles(b, in, true, false, false, false);
			}else if(CladeIndex.parse(arg, a, b)){
				//Do nothing
			}else if(b==null && new File(arg).isFile()){
				in.add(arg);
			}else if(b==null && new File(arg).isDirectory()){
				Tools.getFileOrFiles(arg, in, true, false, false, false);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/** 
	 * Ensures input files can be read and output files can be written.
	 * Throws an exception if file access issues are detected.
	 */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in, ref)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in.get(0), ref.get(0), out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** 
	 * Adjusts file-related static fields as needed for this program.
	 * Optimizes file reading mode based on thread count.
	 */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/** 
	 * Ensures parameter ranges are within bounds and required parameters are set.
	 * 
	 * @return true if parameters are valid
	 */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
//		assert(false) : "TODO";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** 
	 * Main processing method.
	 * Creates read streams, processes all data, and writes results.
	 * 
	 * @param t0 Timer for overall execution time tracking
	 */
	void process(Timer t0){
		Timer t=new Timer(outstream, false);
		
		//Reset counters
		readsProcessed=basesProcessed=0;
		setup();
		loadIndex();
		loadQueries();
		ArrayList<Object> results;
		t.start();
		if(multithreaded) {
			results=spawnThreads(queries, index);
		}else {
			results=searchST(queries, index);
		}
		
		if(printMetrics) {
			evaluate(results);
		}
		outstream.println("Made "+index.comparisons+" fast and "+index.slowComparisons+" slow comparisons.");
		t.stop("Searched "+in.size()+" queries in ");
		
		if(ffout!=null) {
//			outstream.println();
			write(ffout, results);
		}
		
		//Report timing and results
		t0.stop();
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t0, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.things("Queries", queries.size(), 8));
		outstream.println(Tools.things("Bytes Out", bytesOut, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+
					" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Performs single-threaded search of queries against the index.
	 * 
	 * @param queries List of query Clades to search
	 * @param index The CladeIndex to search against
	 * @return List of Comparison results in the same order as queries
	 */
	public ArrayList<Object> searchST(ArrayList<Clade> queries, CladeIndex index) {
		ArrayList<Object> results=new ArrayList<Object>();
		for(Clade query : queries) {
			Object c=index.findBest(query);
			results.add(c);
		}
		return results;
	}
	
	/**
	 * Evaluates search results and prints performance metrics.
	 * Calculates correctness statistics and taxonomic level accuracy.
	 * 
	 * @param results List of Comparison results to evaluate
	 */
	public void evaluate(ArrayList<Object> results) {
		int correct=0, wrong=0, nohit=0, other=0;
		int[] levels=new int[TaxTree.numTaxLevelNames];
		
		double[][] metrics=new double[7][TaxTree.numTaxLevelNames];
		String[] mNames= {"r34", "r45", "strand", "gc", "3mer", "4mer", "5mer"};
		double[] ratio34=metrics[0];
		double[] ratio45=metrics[1];
		double[] strandDif=metrics[2];
		double[] gcDif=metrics[3];
		double[] triDif=metrics[4];
		double[] tetDif=metrics[5];
		double[] pentDif=metrics[6];
		
		for(int i=0; i<queries.size(); i++) {
			Object o=results.get(i);
			@SuppressWarnings("unchecked")
			Comparison c=(o.getClass()==Comparison.class ? (Comparison)o : 
				((ArrayList<Comparison>)o).get(0));
			if(c==null || c.ref==null) {nohit++;}
			else if(c.correct()) {correct++;}
			else if(c.incorrect()) {wrong++;}
			else {other++;}
			int level=levels.length-1;
			
			if(useTree) {level=c.correctLevel();}
			levels[level]++;
			
			if(c!=null && c.ref!=null) {
				ratio34[level]+=c.k3dif/Math.max(c.k4dif, 0.000001f);
				ratio45[level]+=c.k4dif/Math.max(c.k5dif, 0.000001f);
				strandDif[level]+=c.strdif;
				gcDif[level]+=c.gcdif;
				triDif[level]+=c.k3dif;
				tetDif[level]+=c.k4dif;
				pentDif[level]+=c.k5dif;
			}
		}
		outstream.println("Results:\n");
		outstream.println("Correct:  \t"+correct);
		outstream.println("Incorrect:\t"+wrong);
		outstream.println("No-Hit:   \t"+nohit);
		outstream.println("Other:    \t"+other);
		
		if(useTree) {
			outstream.println("\nBest Hit Accuracy:");
			long score=0, sum=0;
			for(int level=0; level<levels.length; level++) {
				int x=levels[level];
				if(x>0) {
					sum+=x;
					score+=x*(levels.length-level-1);
					outstream.println(Tools.padRight(TaxTree.levelToString(level)+":", 9)+"\t"+x);
				}
			}
			outstream.println("\n"+Tools.padRight("Score:", 9)+String.format("\t%.3f", score/(float)sum));
		}

		if(printMetrics) {
			outstream.println("\nMetrics:");
			outstream.print("Level");
			for(int level=0; level<levels.length; level++) {outstream.print("\t"+TaxTree.levelToStringShort(level));}
			outstream.print("\nCount");
			for(int level=0; level<levels.length; level++) {outstream.print("\t"+levels[level]);}
			outstream.println();
			for(int metric=0; metric<metrics.length; metric++) {
				outstream.print(mNames[metric]);
				for(int level=0; level<levels.length; level++) {
					int count=levels[level];
					if(count==0) {outstream.print('\t');}
					else {outstream.print(String.format("\t%.4f", metrics[metric][level]/count));}
				}
				outstream.println();
			}
		}
	}
	
	/**
	 * Writes comparison results to the specified output format.
	 * 
	 * @param ff Output file format
	 * @param comps Collection of comparison results to write
	 */
	void write(FileFormat ff, Collection<Object> comps) {
		if(ff==null) {return;}
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ff);
		ByteBuilder bb=new ByteBuilder(1024);
		if(format==MACHINE) {bsw.println(Comparison.machineHeader(printQTID));}
		for(Object c : comps) {
			bb.clear();
			if(maxHitsToPrint>0) {
				appendResult(c, bb);
				bsw.print(bb);
			}
		}
		bsw.poison();
	}
	
	ByteBuilder appendResult(Object o, ByteBuilder bb) {
		if(o.getClass()==Comparison.class) {
			return appendResult((Comparison)o, bb, 0);
		}
		@SuppressWarnings("unchecked")
		Collection<Comparison> coll=(Collection<Comparison>)o;
		int i=0;
		for(Comparison c : coll) {
			if(i>=maxHitsToPrint) {break;}
			appendResult(c, bb, i);
			i++;
		}
		return bb;
	}
	
	/**
	 * Appends a comparison result to a ByteBuilder in the configured format.
	 * 
	 * @param c Comparison result to append
	 * @param bb ByteBuilder to append to
	 * @return The ByteBuilder with appended result
	 */
	ByteBuilder appendResult(Comparison c, ByteBuilder bb, int hitNum) {
		if(format==MACHINE) {return appendResultMachine(c, bb);}
		else {return appendResultHuman(c, bb, hitNum);}
	}
	
	/**
	 * Appends a comparison result in human-readable format.
	 * 
	 * @param c Comparison result to append
	 * @param bb ByteBuilder to append to
	 * @return The ByteBuilder with appended result
	 */
	ByteBuilder appendResultHuman(Comparison c, ByteBuilder bb, int hitNum) {
		if(bb==null) {bb=new ByteBuilder();}
		c.appendResultHuman(bb, hitNum);
		bb.nl().nl();
		bytesOut+=bb.length;
		linesOut++;
		return bb;
	}
	
	/**
	 * Appends a comparison result in machine-readable format.
	 * 
	 * @param c Comparison result to append
	 * @param bb ByteBuilder to append to
	 * @return The ByteBuilder with appended result
	 */
	ByteBuilder appendResultMachine(Comparison c, ByteBuilder bb) {
		if(bb==null) {bb=new ByteBuilder();}
		c.appendResultMachine(printQTID, bb);
		bb.nl();
		bytesOut+=bb.length;
		linesOut++;
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Spawns multiple threads to process queries in parallel.
	 * 
	 * @param queries List of query Clades to process
	 * @param index The CladeIndex to search against
	 * @return List of Comparison results in the same order as queries
	 */
	private ArrayList<Object> spawnThreads(final ArrayList<Clade> queries, CladeIndex index){
		
		//Do anything necessary prior to processing
		ArrayList<Object> results=new ArrayList<Object>(queries.size());
		
		//Determine how many threads may be used
		final int threads=Tools.mid(1, Shared.threads(), Tools.min(maxCompareThreads, queries.size()/16));
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(queries, index, i, threads));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		while(results.size()<queries.size()) {results.add(null);}
		for(int tnum=0; tnum<threads; tnum++) {
			ProcessThread pt=alpt.get(tnum);
			synchronized(pt) {
				for(int i=0, j=tnum; i<pt.results.size(); i++, j+=threads) {
					Object o=pt.results.get(i);
					results.set(j, o);
				}
			}
		}
		return results;
	}
	
	/**
	 * Accumulates statistics from a finished ProcessThread.
	 * Implementation of Accumulator interface.
	 * 
	 * @param pt ProcessThread to accumulate statistics from
	 */
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt) {
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			index.comparisons+=pt.index.comparisons;
			index.slowComparisons+=pt.index.slowComparisons;
			errorState|=(!pt.success);
		}
	}
	
	/**
	 * Returns whether all processing was successful.
	 * Implementation of Accumulator interface.
	 * 
	 * @return true if no errors occurred during processing
	 */
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Worker thread for parallel clade search processing.
	 * Each thread processes a subset of queries assigned by thread ID.
	 */
	static class ProcessThread extends Thread {
		
		/**
		 * Constructs a new ProcessThread.
		 * 
		 * @param queries_ List of query Clades to process
		 * @param index_ The CladeIndex to search against
		 * @param tid_ Thread ID
		 * @param threads_ Total number of threads
		 */
		ProcessThread(ArrayList<Clade> queries_, CladeIndex index_, final int tid_, final int threads_){
			queries=queries_;
			index=index_.clone(); // Clone to avoid synchronization issues
			tid=tid_;
			threads=threads_;
		}
		
		/**
		 * Main thread execution method.
		 * Processes queries assigned to this thread and stores results.
		 */
		@Override
		public synchronized void run(){
			//Do anything necessary prior to processing
			IDAligner ssa=(Clade.callSSU ? aligner.Factory.makeIDAligner() : null);
			//Run queries
			for(int i=tid; i<queries.size(); i+=threads) {
				Clade clade=queries.get(i);
				synchronized(clade) {
					readsProcessedT+=clade.contigs;
					basesProcessedT+=clade.bases;
					ArrayList<Comparison> list=index.findBest(clade);
					results.add(list);
					if(list!=null && Clade.callSSU) {
						for(Comparison comp : list) {comp.align(ssa);}
					}
					Collections.sort(list);
				}
			}
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		final CladeIndex index;
		
		/** Shared data source */
		private final ArrayList<Clade> queries;
		/** Clade storage */
		private final ArrayList<Object> results=new ArrayList<Object>();
		/** Thread ID */
		final int tid;
		final int threads;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private ArrayList<String> in=new ArrayList<String>();

	/** Primary output file path */
	private String out="stdout.txt";
	
	/** Reference database paths */
	private ArrayList<String> ref=new ArrayList<String>();
	/** Loaded query clades */
	private ArrayList<Clade> queries;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** The clade index for reference searching */
	public CladeIndex index;
	
	/** Whether to print detailed comparison metrics */
	boolean printMetrics=false;
	/** Whether to use multiple threads for searching */
	boolean multithreaded=true;
	/** Whether to process each contig separately */
	boolean perContig=false;
	/** Minimum contig length to process */
	int minContig=0;
	/** Maximum number of threads to use for comparison */
	int maxCompareThreads=9999;
	
	/** Whether to print query taxon IDs in output */
	boolean printQTID=false;
	/** Output format (HUMAN or MACHINE) */
	int format=1;
	/** Format constants */
	public static final int HUMAN=1, MACHINE=2;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of lines out */
	protected long linesOut=0;
	/** Number of bytes out */
	protected long bytesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	private int maxHitsToPrint=1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary output file */
	private FileFormat ffout;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Whether to use taxonomy tree for evaluation */
	static boolean useTree=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public static boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	/** Reads are output in input order */
	private boolean ordered=false;
	
}