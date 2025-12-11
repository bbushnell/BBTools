package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import aligner.IDAligner;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Splits a mix of ribosomal sequences (such as Silva) into different files per
 * type (16S, 18S, etc).
 * Processes FASTA/FASTQ input and sorts sequences by ribosomal RNA type using alignment
 * to consensus sequences. Outputs separate files for each type with configurable thresholds.
 *
 * @author Brian Bushnell
 * @date November 19, 2015
 */
public class SplitRibo implements Accumulator<SplitRibo.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Program entry point for ribosomal sequence splitting.
	 * @param args Command line arguments */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		SplitRibo x=new SplitRibo(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command-line arguments and initializes processing parameters.
	 * Sets up input/output file formats, loads consensus sequences, and validates configuration.
	 * @param args Command-line arguments containing file paths and processing options
	 */
	public SplitRibo(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		Shared.capBufferLen(50);
		ReadWrite.ZIPLEVEL=9;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			qfin1=parser.qfin1;
			extin=parser.extin;

			outPattern=parser.out1;
			extout=parser.extout;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		
		numTypes=sequenceTypes.length;
		readsOut=new long[numTypes];
		basesOut=new long[numTypes];
		consensusSequences=loadConsensusSequenceFromFile();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses command-line arguments using the BBTools Parser framework.
	 * Handles parameters for input/output files, sequence types, and alignment thresholds.
	 * @param args Command-line argument array
	 * @return Configured Parser object with processed arguments
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
			}else if(a.equalsIgnoreCase("minid")){
				minID=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minid2") || a.equalsIgnoreCase("refineid")){
				refineID=Float.parseFloat(b);
			}else if(a.equals("out") || a.equals("pattern") || a.equals("outpattern")){
				parser.out1=b;
			}else if(a.equals("type") || a.equals("types")){
				parseTypes(b);
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
	
	private void parseTypes(String b){
		sequenceTypes=null;
		if(b==null){
			assert(false) : "'types' flag requires a list of types, such as 'types=16S,18S'";
			sequenceTypes=new String[] {"Other"};
		}else{
			String[] split=b.split(",");
			sequenceTypes=new String[split.length+1];
			sequenceTypes[0]="Other";
			for(int i=0; i<split.length; i++){
				String s=split[i].replace('s', 'S');
				if(s.startsWith("its")){s=s.replaceFirst("its", "ITS");}
				sequenceTypes[i+1]=s;
			}
		}
	}
	
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		qfin1=Tools.fixExtension(qfin1);
	}
	
	private void checkFileExistence(){
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		if(outPattern==null){return;}
		
		if(!outPattern.contains("#")){
			throw new RuntimeException("OutPattern must contain '#' symbol: "+outPattern);
		}
		
		for(String type : sequenceTypes) {
			String out=outPattern.replaceFirst("#", type);
			
			//Ensure output files can be written
			if(!Tools.testOutputFiles(overwrite, append, false, out)){
				outstream.println((outPattern==null)+", "+(out==null)+", "+outPattern+", "+out);
				throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
			}

			//Ensure that no file was specified multiple times
			if(!Tools.testForDuplicateFiles(true, in1, out)){
				throw new RuntimeException("\nSome file names were specified multiple times.\n");
			}
		}
	}
	
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/**
	 * Validates that required parameters are set and within acceptable ranges.
	 * @return true if all parameters are valid
	 * @throws RuntimeException if required input file is not specified
	 */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		return true;
	}
	
	private final Read[][] loadConsensusSequenceFromFile(){
		Read[][] seqs=new Read[numTypes][];
		m16S_index=Tools.find("m16S", sequenceTypes);
		m18S_index=Tools.find("m18S", sequenceTypes);
		p16S_index=Tools.find("p16S", sequenceTypes);
		boolean stripM16S=(m16S_index>=0);
		boolean stripM18S=(m18S_index>=0);
		boolean stripP16S=(p16S_index>=0);
		for(int st=1; st<numTypes; st++){
			String name=sequenceTypes[st];
			boolean is16S=name.equalsIgnoreCase("16S");
			boolean is18S=name.equalsIgnoreCase("18S");
			seqs[st]=ProkObject.loadConsensusSequenceType(name, ((is16S && stripM16S) || (is18S && stripM18S)), (is16S && stripP16S));
		}
		return seqs;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main processing method that orchestrates the sequence splitting pipeline.
	 * Creates input/output streams, spawns processing threads, and reports results.
	 * @param t Timer for tracking execution time and performance metrics
	 */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris();
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream[] rosa=makeCrosArray();
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		Arrays.fill(readsOut, 0);
		Arrays.fill(basesOut, 0);
		
		//Process the reads in separate threads
		spawnThreads(cris, rosa);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//assert(!errorState);
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, rosa);
		//assert(!errorState);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;

		long readsOut2=shared.Vector.sum(readsOut)-readsOut[0];
		long basesOut2=shared.Vector.sum(basesOut)-basesOut[0];
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut2, basesOut2, 8, true));

		outstream.println();
		outstream.println(Tools.string("Type", "Count", 8));
		for(int type=0; type<numTypes; type++){
			outstream.println(Tools.number(sequenceTypes[type], readsOut[type], 8));
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private ConcurrentReadInputStream makeCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, qfin1, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	private ConcurrentReadOutputStream[] makeCrosArray(){
		ConcurrentReadOutputStream[] rosa=new ConcurrentReadOutputStream[numTypes];
		for(int i=0; i<numTypes; i++){
			String type=sequenceTypes[i];
			final ConcurrentReadOutputStream ros=makeCros(type);
			rosa[i]=ros;
		}
		return rosa;
	}
	
	private ConcurrentReadOutputStream makeCros(String type){
		if(outPattern==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(2, 16, (Shared.threads()*2)/3) : 4);
		final String fname=outPattern.replaceFirst("#", type);
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.FASTA, extout, true, overwrite, append, ordered);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates and manages worker threads for processing sequences.
	 * Spawns multiple ProcessThread instances and waits for completion.
	 * @param cris Input stream for reading sequences
	 * @param rosa Array of output streams for each sequence type
	 */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream[] rosa){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, rosa, i));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		//assert(!errorState);
		
		//Do anything necessary after processing
		
	}
	
	/**
	 * Accumulates processing statistics from completed worker threads.
	 * Combines read/base counts and error states from all threads.
	 * @param pt ProcessThread containing results to accumulate
	 */
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		Tools.add(readsOut, pt.readsOutT);
		Tools.add(basesOut, pt.basesOutT);
		errorState|=(!pt.success);
		//assert(!errorState);
	}
	
	/** Reports overall processing success status.
	 * @return true if no errors occurred during processing */
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Worker thread for processing sequences and classifying them by ribosomal RNA type.
	 * Each thread reads batches of sequences, aligns them to consensus references,
	 * and routes them to appropriate output streams based on best alignment identity.
	 */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream[] rosa_, final int tid_){
			cris=cris_;
			rosa=rosa_;
			tid=tid_;
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
		
		/** Main processing loop that reads sequence batches and processes them.
		 * Continues until input stream is exhausted, processing each batch through processList. */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
				assert(r.mate==null);
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
			
			@SuppressWarnings("unchecked")
			final ArrayList<Read>[] out=new ArrayList[numTypes];
			for(int i=0; i<numTypes; i++){
				ArrayList<Read> list=new ArrayList<Read>(50);
				out[i]=list;
			}
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r1.length();
				final int initialLength2=r1.mateLength();

				//Increment counters
				readsProcessedT+=r1.pairCount();
				basesProcessedT+=initialLength1+initialLength2;
				
				{
					//Reads are processed in this block.
					final int type=processRead(r1);
					readsOutT[type]+=r1.pairCount();
					basesOutT[type]+=r1.pairLength();
					out[type].add(r1);
				}
			}

			//Output reads to the output stream
			if(rosa!=null){
				for(int type=0; type<numTypes; type++){
					rosa[type].add(out[type], ln.id);
				}
			}
		}
		
		/**
		 * Classifies a single read by aligning it to consensus sequences for each type.
		 * Uses two-stage alignment: first to overall consensus, then to clade-specific
		 * consensuses if needed. Returns the best matching type index.
		 *
		 * @param r The read to classify
		 * @return Index of best matching sequence type, or 0 for unclassified
		 */
		private int processRead(final Read r){
			int bestType=0;
			float bestID=-1;
			for(int type=1; type<numTypes; type++){//Align to only the overall consensus
				Read[] refs=consensusSequences[type];
				float id=align(r, refs, 0, 1);
				if(id>bestID && id>=minID){
					bestType=type;
					bestID=id;
				}
			}
			if(bestType<1 || bestID<refineID || bestType==p16S_index){//If nothing met minID, or if it matched chloro, align to clade-specific consensuses
				for(int type=1; type<numTypes; type++){
					Read[] refs=consensusSequences[type];
					float id=align(r, refs, 1, refs.length);
					if(id>bestID && id>=minID){
						bestType=type;
						bestID=id;
					}
				}
			}
			r.obj=bestID;//If desired...  in actuality, more info might be useful, like alignment length
			return bestID<minID ? 0 : bestType;
		}
		
		private float align(Read r, Read[] refs, int minRef, int maxRef){
			float bestID=-1;
			if(refs!=null){
				for(int i=minRef; i<maxRef; i++){
					Read ref=refs[i];
					float id=ssa.align(r.bases, ref.bases);
					bestID=Tools.max(id,  bestID);
				}
			}
			return bestID;
		}
		
		IDAligner ssa=aligner.Factory.makeIDAligner();

		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		
		/** Number of reads retained by this thread for each sequence type */
		protected long[] readsOutT=new long[numTypes];
		/** Number of bases retained by this thread for each sequence type */
		protected long[] basesOutT=new long[numTypes];
		
		boolean success=false;
		
		private final ConcurrentReadInputStream cris;
		private final ConcurrentReadOutputStream[] rosa;
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;
	
	private String qfin1=null;

	/** Output file pattern with '#' placeholder for sequence type names */
	private String outPattern=null;
	
	private String extin=null;
	private String extout=null;

	float minID=0.59f; //This could be a per-type value
	float refineID=0.70f; //Refine alignment if best is less than this
	
	private int m16S_index=-2;
	private int m18S_index=-2;
	private int p16S_index=-2;
	
	/*--------------------------------------------------------------*/
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	
	private long maxReads=-1;
	
	private String[] sequenceTypes=new String[] {"Other", "16S", "18S", "23S", "5S", "m16S", "m18S", "p16S"};
	private final int numTypes;//=sequenceTypes.length;
	final Read[][] consensusSequences;
	
	/** Count of reads output to each sequence type file */
	final long[] readsOut;
	/** Count of bases output to each sequence type file */
	final long[] basesOut;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** FileFormat object for primary input file */
	private final FileFormat ffin1;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean ordered=true;
	
}
