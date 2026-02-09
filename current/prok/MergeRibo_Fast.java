package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.IDAligner;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tax.GiToTaxid;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Picks one ribosomal 16S sequence per taxID from input files using a global
 * consensus sequence. This legacy implementation does not build per-taxID
 * consensus sequences; it simply retains the read with the best identity to
 * the global consensus for each taxonomic identifier, using multi-threaded
 * processing for speed.
 * @author Brian Bushnell
 * @date November 19, 2015
 */
public class MergeRibo_Fast implements Accumulator<MergeRibo_Fast.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Command-line entry point.
	 * @param args Command line arguments */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		MergeRibo_Fast x=new MergeRibo_Fast(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Parses command line arguments, initializes IO formats, validates parameters,
	 * and configures processing options for a MergeRibo_Fast run.
	 * @param args Command line arguments containing file paths and options
	 */
	public MergeRibo_Fast(String[] args){
		
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

			out1=parser.out1;
			extout=parser.extout;
		}

		validateParams();
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin=new ArrayList<FileFormat>(in.size());
		ffalt=FileFormat.testInput(alt, FileFormat.FASTA, extin, true, true);
		for(String s : in){
			FileFormat ff=FileFormat.testInput(s, FileFormat.FASTA, extin, true, true);
			ffin.add(ff);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses command line arguments into a Parser configuration.
	 * Handles input files, output paths, and basic processing flags.
	 * @param args Array of command line arguments
	 * @return Configured Parser with parsed parameters
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
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("in")){
				Tools.addFiles(b, in);
			}else if(a.equals("alt")){
				alt=b;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		assert(!in.isEmpty()) : "No input file.";
		return parser;
	}
	
	/** Validates that input files exist and are readable and that output files can
	 * be written, throwing a runtime exception if any check fails. */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
//		//Ensure that no file was specified multiple times
//		if(!Tools.testForDuplicateFiles(true, out1, in.toArray(new String[0]))){
//			throw new RuntimeException("\nSome file names were specified multiple times.\n");
//		}
	}
	
	/** Configures interleaving settings for FASTQ input and disables forced
	 * interleaving since this tool processes single-ended ribosomal reads. */
	private void adjustInterleaving(){
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
	}
	
	/** Adjusts static file-handling settings for this program, enabling BF2 mode
	 * for multi-threaded file reading when appropriate. */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/**
	 * Validates configuration parameters and serves as an extension point for
	 * future checks; currently always returns true.
	 * @return true if parameters are considered valid
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
	 * Main processing pipeline that loads the global 16S (and optional 18S)
	 * consensus sequences, processes input files in parallel threads, and writes
	 * the best ribosomal sequence per taxID.
	 * @param t Timer for tracking execution performance
	 */
	void process(Timer t){

		if(process16S){
			Read[] data=ProkObject.loadConsensusSequenceType("16S", true, true);
			consensus16S=data[0].bases;
		}
		if(process18S){
			Read[] data=ProkObject.loadConsensusSequenceType("18S", true, true);
			consensus18S=data[0].bases;
		}
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		for(FileFormat ff : ffin) {
			//Create a read input stream
			final ConcurrentReadInputStream cris=makeCris(ff);

			//Process the reads in separate threads
			spawnThreads(cris);
			errorState|=ReadWrite.closeStream(cris);
		}
		
		//Do anything necessary after processing
		if(ffout1!=null){
			//Optionally create a read output stream
			final ConcurrentReadOutputStream ros=makeCros(false);
			long num=0;
			for(Entry<Integer, Ribo> e : bestMap.entrySet()){
				Read r=e.getValue().r;
				readsOut++;
				basesOut+=r.length();
				ArrayList<Read> list=new ArrayList<Read>(1);
				list.add(r);
				ros.add(list, num);
				num++;
			}
			//Close the read streams
			errorState|=ReadWrite.closeStream(ros);
		}
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private ConcurrentReadInputStream makeCris(FileFormat ff){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		assert(!paired) : "This should not be paired input.";
		return cris;
	}
	
	private ConcurrentReadOutputStream makeCros(boolean pairedInput){
		if(ffout1==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates and runs worker threads to process reads from the input stream.
	 * Spawns one ProcessThread per available worker thread, starts them, waits
	 * for completion, and records any error state.
	 * @param cris Input stream to be processed by worker threads
	 */
	private void spawnThreads(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, i));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
	}
	
	/**
	 * Accumulates statistics from a completed worker thread by updating global
	 * read and base counters and aggregating error status.
	 * @param pt Completed ProcessThread containing statistics
	 */
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		errorState|=(!pt.success);
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Worker thread that processes ribosomal sequences from the shared input
	 * stream, aligns each read to consensus sequences, extracts taxonomic IDs,
	 * and tracks the best-scoring sequence per taxID.
	 */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final int tid_){
			cris=cris_;
			tid=tid_;
		}
		
		//Called by start()
		/** Main thread body that repeatedly fetches lists of reads from the input
		 * stream, delegates to processInner(), and marks completion status. */
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/**
		 * Core processing loop that fetches read lists from the input stream, calls
		 * processList(...) on each list, and returns lists to the stream until no
		 * more reads are available.
		 */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(ln);
				
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
		}
		
		void processList(ListNum<Read> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r1.length();

				//Increment counters
				readsProcessedT++;
				basesProcessedT+=initialLength1;
				
				boolean keep=processRead(r1);
			}
		}
		
		/**
		 * Processes a single ribosomal read by extracting its taxonomic ID, aligning
		 * it to the consensus sequence(s), and updating the best sequence map.
		 * @param r Ribosomal sequence read to process
		 * @return true if this read becomes the current best for its taxID
		 */
		boolean processRead(final Read r){
			Integer key=GiToTaxid.parseTaxidNumber(r.id, '|');
			if(key==null || key==-1){return false;}
			float id=align(r);
			float product=id*r.length();
			Ribo ribo=new Ribo(r, key, id);
			
			
			synchronized(bestMap){
				Ribo old=bestMap.get(key);
				if(old==null || old.product<product){
					bestMap.put(key, ribo);
					return true;
				}
			}
			return false;
		}
		
		float align(Read r){
			float a=(process16S ? ssa.align(r.bases, consensus16S) : 0);
			float b=(process18S ? ssa.align(r.bases, consensus18S) : 0);
			return Tools.max(a, b);
		}
		
		IDAligner ssa=idaligner.Factory.makeIDAligner();

		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		
		boolean success=false;
		
		private final ConcurrentReadInputStream cris;
		final int tid;
	}
	
	private static class Ribo implements Comparable<Ribo>{
		Ribo(Read r_, int tid_, float identity_){
			r=r_;
			tid=tid_;
			identity=identity_;
			product=r.length()*identity;
		}
		
		/**
		 * Compares Ribo objects for sorting by quality metrics, first by product
		 * (length × identity) and then by sequence length.
		 * @param o Other Ribo object to compare against
		 * @return Negative if this is better, positive if other is better, 0 if equal
		 */
		@Override
		public int compareTo(Ribo o) {
			if(o.product>product){return -1;}
			else if(o.product<product){return 1;}
			else if(o.r.length()>r.length()){return -1;}
			else if(o.r.length()<r.length()){return 1;}
			return 0;
		}
		
		Read r;
		int tid;
		float identity;
		float product;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file paths. */
	private ArrayList<String> in=new ArrayList<String>();
	
	private String alt=null;
	
	private String out1=null;
	
	private String extin=null;
	private String extout=null;

	HashMap<Integer, Ribo> bestMap=new HashMap<Integer, Ribo>(10000000);
	HashMap<Integer, ArrayList<Ribo>> listMap=new HashMap<Integer, ArrayList<Ribo>>(10000000);
	
	static byte[] consensus16S;
	byte[] consensus18S;
	
	/*--------------------------------------------------------------*/

	protected long readsProcessed=0;
	protected long basesProcessed=0;

	protected long readsOut=0;
	protected long basesOut=0;

	private long maxReads=-1;

	private boolean process16S=true;
	private boolean process18S=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file formats. */
	private final ArrayList<FileFormat> ffin;
	private final FileFormat ffalt;
	
	/** Primary output file format. */
	private final FileFormat ffout1;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean ordered=false;
	
}
