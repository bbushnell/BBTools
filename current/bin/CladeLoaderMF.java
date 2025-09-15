package bin;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import prok.CallGenes;
import prok.GeneCaller;
import prok.GeneModelParser;
import prok.Orf;
import shared.LineParser1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;
import tracker.ReadStats;

/**
 * Designed to load one clade per file in metagenomic bins.
 * Unlike CladeLoader, this class is optimized for multiple
 * files with one taxonomic group per file
 * (rather than a single file with multiple taxonomic IDs).
 * Uses one thread per file for parallel processing.
 * 
 * @author Brian Bushnell
 * @date April 12, 2025
 */
public class CladeLoaderMF extends BinObject implements Accumulator<CladeLoaderMF.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		
		//Create an instance of this class
		CladeLoaderMF x=new CladeLoaderMF(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Default constructor.
	 * Loads entropy model if needed.
	 */
	public CladeLoaderMF(){
		if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150) {
			AdjustEntropy.load(4, 150);
		}
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CladeLoaderMF(String[] args){
		
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
			
			extin=parser.extin;

			out=parser.out1;
			extout=parser.extout;
		}

		validateParams();
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout=FileFormat.testOutput(out, FileFormat.CLADE, extout, true, overwrite, append, ordered);
		
		if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150) {
			AdjustEntropy.load(4, 150);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Parse arguments from the command line 
	 * @param args Arguments to parse
	 * @return Parser with parsed values
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
			
			if(a.equals("tree")){
				if(b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true")) {
					treePath="auto";
				}else if(b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false")) {
					treePath=null;
				}else {
					treePath=b;
				}
			}else if(a.equals("in")){
				Tools.getFileOrFiles(b, in, true, false, false, false);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("percontig")){
				perContig=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("a48")){
				Clade.outputCoding=Parse.parseBoolean(b) ? Clade.A48 : Clade.DECIMAL;
			}else if(a.equals("maxk") || a.equals("kmax")){
				Comparison.maxK=Clade.MAXK=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("callSSU")){
				Clade.callSSU=Parse.parseBoolean(b);
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
	 * Ensure files can be read and written 
	 * @throws RuntimeException If files cannot be read or written
	 */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in.get(0), out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** 
	 * Adjust file-related static fields as needed for this program 
	 */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/** 
	 * Ensure parameter ranges are within bounds and required parameters are set 
	 * @return True if parameters are valid
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
	 * Create read streams and process all data 
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		//Reset counters
		readsProcessed=basesProcessed=0;
		
		ArrayList<Clade> list=loadFiles(in, perContig, minContig);
		
		if(ffout!=null) {
			write(ffout, list);
		}
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.things("Clades", list.size(), 8));
		outstream.println(Tools.things("Bytes Out", bytesOut, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+
					" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Write a collection of clades to a file.
	 * @param ff FileFormat for output
	 * @param coll Collection of clades to write
	 */
	void write(FileFormat ff, Collection<Clade> coll) {
		if(ff==null) {return;}
		ArrayList<Clade> list=new ArrayList<Clade>(coll);
		Collections.sort(list);
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ff);
		ByteBuilder bb=new ByteBuilder(1024);
		for(Clade c : list) {
			bb.clear();
			c.toBytes(bb);
			bsw.print(bb);
			bytesOut+=bb.length;
			linesOut++;//Actually, clades out
		}
		bsw.poison();
	}
	
	/**
	 * Load a single clade from a file.
	 * @param fname File name to load
	 * @param et Entropy tracker
	 * @return Clade loaded from the file
	 */
	public static Clade loadOneClade(String fname, EntropyTracker et) {
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, false);
		if(ff.clade()) {
			return loadOneCladeFromClade(ff);
		}else {
			return loadOneCladeFromSequence(ff, et);
		}
	}
	
	/**
	 * Load a single clade from a clade-format file.
	 * @param ff FileFormat for input
	 * @return Clade loaded from the file
	 */
	public static Clade loadOneCladeFromClade(FileFormat ff) {
		ArrayList<byte[]> lines=ByteFile.toLines(ff);
		Clade c=Clade.parseClade(lines, new LineParser1('\t'));
		return c;
	}
	
	/**
	 * Load a single clade from a sequence file.
	 * @param ff FileFormat for input
	 * @param et Entropy tracker
	 * @return Clade loaded from the file
	 */
	private static Clade loadOneCladeFromSequence(FileFormat ff, EntropyTracker et) {
		ArrayList<Read> list=ConcurrentReadInputStream.getReads(-1, false, ff, null, null, null);
		if(list==null || list.isEmpty()) {return null;}
		
		final GeneCaller caller=Clade.callSSU ? GeneTools.makeGeneCaller() : null;
		final int tid=resolveTaxID(list.get(0).id);
		Clade c=new Clade(tid, -1, ff.simpleName());
		if(tree!=null) {c.level=tree.toLevel(tid);}
		synchronized(c) {
			for(Read r : list) {
				c.add(r, et, caller);
			}
			c.finish();
		}
		return c;
	}
	
	/**
	 * Load clades from a file.
	 * @param fname File name to load
	 * @param et Entropy tracker
	 * @param perContig Whether to treat each contig as a separate clade
	 * @param minContig Minimum contig length to process
	 * @return List of clades loaded from the file
	 */
	public static ArrayList<Clade> loadClades(String fname, EntropyTracker et, 
			boolean perContig, int minContig) {
//		assert(false) : perContig+", "+fname;
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, false);
		if(ff.clade()) {
			return loadCladesFromClade(ff);
		}else if(!perContig) {
			ArrayList<Clade> list=new ArrayList<Clade>(1);
			Clade c=loadOneCladeFromSequence(ff, et);
			list.add(c);
			return list;
		}else {
			return loadCladesFromSequence(ff, et, minContig);
		}
	}
	
	/**
	 * Load clades from a clade-format file.
	 * @param ff FileFormat for input
	 * @return List of clades loaded from the file
	 */
	public static ArrayList<Clade> loadCladesFromClade(FileFormat ff) {
		ByteFile bf=ByteFile.makeByteFile(ff);
		LineParser1 lp=new LineParser1('\t');
		ArrayList<Clade> out=new ArrayList<Clade>();
		final ArrayList<byte[]> set=new ArrayList<byte[]>(20); //Should be 13 lines
		for(ListNum<byte[]> list=bf.nextList(); list!=null; list=bf.nextList()) {
			for(byte[] line : list) {
				if(Tools.startsWith(line, '#') && set.size()>5) {//New record
					addClade(set, out, lp);
					set.clear();
				}
				set.add(line);
			}
		}
		addClade(set, out, lp);
		return out;
	}
	
	private static boolean addClade(final ArrayList<byte[]> set, 
			ArrayList<Clade> clades, LineParser1 lp) {
		if(set.size()<=5) {return false;}
		Clade c=Clade.parseClade(set, lp);
		synchronized(c) {
			c.finish();
			clades.add(c);
		}
		return true;
	}
	
	/**
	 * Create a ConcurrentReadInputStream for a FileFormat.
	 * @param ff FileFormat to create stream from
	 * @return ConcurrentReadInputStream for the FileFormat
	 */
	private static ConcurrentReadInputStream makeCris(FileFormat ff){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null);
		cris.start(); //Start the stream
		if(verbose){System.err.println("Started cris");}
		boolean paired=cris.paired();
		if(ff.fastq()){System.err.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		else {assert(!paired);}
		return cris;
	}
	
	/**
	 * Load clades from a sequence file.
	 * @param ff FileFormat for input
	 * @param et Entropy tracker
	 * @param minContig Minimum contig length to process
	 * @return List of clades loaded from the file
	 */
	private static ArrayList<Clade> loadCladesFromSequence(FileFormat ff, 
			EntropyTracker et, int minContig) {
		ConcurrentReadInputStream cris=makeCris(ff);
		ArrayList<Clade> list=processReads(cris, et, minContig);
		ReadWrite.closeStream(cris);
		return list;
	}
	
	/**
	 * Process reads from a stream into clades.
	 * @param cris Input read stream
	 * @param et Entropy tracker
	 * @param minContig Minimum contig length to process
	 * @return List of clades processed from the reads
	 */
	private static ArrayList<Clade> processReads(ConcurrentReadInputStream cris, 
			EntropyTracker et, int minContig) {
		ArrayList<Clade> out=new ArrayList<Clade>();

		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();

		//Check to ensure pairing is as expected
		if(ln!=null && !ln.isEmpty()){
			Read r=ln.get(0);
		}

		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){
			for(Read r : ln) {
				if(r.length()<minContig) {continue;}
				int tid=resolveTaxID(r.id);
				Clade c=new Clade(tid, -1, r.id);
				synchronized(c) {
					if(tree!=null) {c.level=tree.toLevel(tid);}
					c.add(r.bases, et);
					c.finish();
				}
				out.add(c);
			}

			//Notify the input stream that the list was used
			cris.returnList(ln);
			//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

			//Fetch a new list
			ln=cris.nextList();
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		return out;
	}
	
	/**
	 * Load clades from multiple files.
	 * @param in List of file paths to load
	 * @param perContig Whether to treat each contig as a separate clade
	 * @param minContig Minimum contig length to process
	 * @return List of all clades loaded from the files
	 */
	public ArrayList<Clade> loadFiles(ArrayList<String> in, boolean perContig, int minContig){
		ArrayList<Clade> list=spawnThreads(in, perContig, minContig);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Spawn process threads to handle multiple files in parallel.
	 * @param files List of file paths to process
	 * @param perContig Whether to treat each contig as a separate clade
	 * @param minContig Minimum contig length to process
	 * @return List of all clades loaded from the files
	 */
	private ArrayList<Clade> spawnThreads(final ArrayList<String> files, boolean perContig, int minContig){
		
		//Do anything necessary prior to processing
		ConcurrentHashMap<Integer, ArrayList<Clade>> cladeMap=new ConcurrentHashMap<Integer, ArrayList<Clade>>(files.size());
		
		//Determine how many threads may be used
		final int threads=Tools.min(files.size(), Shared.threads());
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(files, cladeMap, i, threads, perContig, minContig));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		ArrayList<Clade> out=new ArrayList<Clade>(files.size());
		for(int i=0; i<files.size(); i++) {
			ArrayList<Clade> list=cladeMap.get(i);
			for(Clade c : list) {
				synchronized(c) {
					assert(c.finished());
					out.add(c);
				}
			}
		}
		return out;
	}
	
	/**
	 * Accumulate statistics from a ProcessThread.
	 * @param pt ProcessThread to accumulate from
	 */
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt) {
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			errorState|=(!pt.success);
		}
	}
	
	/**
	 * Check if processing was successful.
	 * @return True if no errors occurred
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
	 * Thread for processing files in parallel.
	 * This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier.
	 */
	static class ProcessThread extends Thread {
		
		/**
		 * Constructor.
		 * @param files_ List of files to process
		 * @param cladeMap_ Map to store processed clades
		 * @param tid_ Thread ID
		 * @param threads_ Total number of threads
		 * @param perContig_ Whether to treat each contig as a separate clade
		 * @param minContig_ Minimum contig length to process
		 */
		ProcessThread(ArrayList<String> files_, ConcurrentHashMap<Integer, ArrayList<Clade>> cladeMap_, 
				final int tid_, final int threads_, boolean perContig_, int minContig_){
			files=files_;
			cladeMap=cladeMap_;
			tid=tid_;
			threads=threads_;
			perContig=perContig_;
			minContig=minContig_;
		}
		
		/**
		 * Process assigned files.
		 * Called by start()
		 */
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			for(int i=tid; i<files.size(); i+=threads) {
				ArrayList<Clade> list=loadClades(files.get(i), et, perContig, minContig);
				cladeMap.put(i, list);
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
		
		private EntropyTracker et=new EntropyTracker(entropyK, entropyWindow, false);
		private final boolean perContig;
		private final int minContig;
		
		/** Shared data source */
		private final ArrayList<String> files;
		/** Clade storage */
		private final ConcurrentHashMap<Integer, ArrayList<Clade>> cladeMap;
		/** Thread ID */
		final int tid;
		/** Total number of threads */
		final int threads;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private ArrayList<String> in=new ArrayList<String>();

	/** Primary output file path */
	private String out=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** Whether to treat each contig as a separate clade */
	boolean perContig=false;
	/** Minimum contig length to process */
	int minContig=0;
	
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
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary output file */
	private FileFormat ffout;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
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