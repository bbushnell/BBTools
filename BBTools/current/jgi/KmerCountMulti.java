package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import cardinality.MultiLogLog;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Counts unique kmers in a file.
 * Tracks multiple kmer lengths independently.
 * 
 * @author Brian Bushnell
 * @date December 30, 2016
 *
 */
public class KmerCountMulti {
	
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
		KmerCountMulti x=new KmerCountMulti(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerCountMulti(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_UNPIGZ=true;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.out1="stdout.txt";
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] split2=b.split(",");
				for(String k : split2){
					parser.loglogKlist.add(Integer.parseInt(k));
				}
			}else if(a.equals("ways") || a.equals("hashes")){
				ways=Integer.parseInt(b);
			}else if(a.equals("showstdev") || a.equals("showstddev") || a.equals("stdev") || a.equals("stddev")){
				showStdev=Parse.parseBoolean(b);
			}else if(a.equals("wavg")){
				useWavg=Parse.parseBoolean(b);
			}else if(a.equals("seed")){
				parser.parse(null, "loglogseed", b);
			}else if(a.equals("buckets")){
				parser.parse(null, "loglogbuckets", b);
			}else if(a.equals("minprob")){
				parser.parse(null, "loglogminprob", b);
			}else if(a.equals("sweep")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] split2=b.split(",");
				int mink=Integer.parseInt(split2[0]);
				int maxk=Integer.parseInt(split2[1]);
				int incr=Integer.parseInt(split2[2]);
				for(int k=mink; k<=maxk; k+=incr){
					parser.loglogKlist.add(k);
				}
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
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out=parser.out1;
			
			extin=parser.extin;
			
			overwrite=parser.overwrite;
			append=parser.append;
		}
		
		mlogArray=new MultiLogLog[ways];
		for(int i=0; i<ways; i++){
			mlogArray[i]=new MultiLogLog(parser);
			if(parser.loglogseed>=0){
				parser.loglogseed++;
			}
		}
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		for(ProcessThread pt : alpt){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		if(ffout!=null){
			writeOutput();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void writeOutput0(){
		TextStreamWriter tsw=new TextStreamWriter(ffout);
		tsw.start();
		if(ffout.stdio()){tsw.println();}
		tsw.print("#K\tCount\n");
		int numK=mlogArray[0].counters.length;
		for(int knum=0; knum<numK; knum++){
			long sum=0;
			for(MultiLogLog mlog : mlogArray){
				sum+=mlog.counters[knum].cardinality();
			}
			tsw.print(mlogArray[0].counters[knum].k+"\t"+((sum+ways-1)/ways)+"\n");
		}
		if(ffout.stdio()){tsw.println();}
		errorState|=tsw.poisonAndWait();
	}
	
	private void writeOutput(){
		TextStreamWriter tsw=new TextStreamWriter(ffout);
		tsw.start();
		if(ffout.stdio()){tsw.println();}
		tsw.print("#K\tCount      "+(ways>1 && showStdev ? "\tStdDev" : "")+"\n");
		
		int numK=mlogArray[0].counters.length;
		long[][] counts=new long[numK][ways];
				
		for(int way=0; way<ways; way++){
			MultiLogLog mlog=mlogArray[way];
			for(int knum=0; knum<numK; knum++){
				counts[knum][way]=mlog.counters[knum].cardinality();
			}
		}
		
		for(int knum=0; knum<numK; knum++){
			long[] array=counts[knum];
			Arrays.sort(array);
			double avg=shared.Vector.sum(array)/(double)ways;
			double wavg=Tools.weightedAverage(array);
			double stdev=Tools.standardDeviation(array)*100/avg;
			int k=mlogArray[0].counters[knum].k;
			String avgs=""+(long)Math.round(useWavg ? wavg : avg);
			while(avgs.length()<11){avgs=avgs+" ";}
			if(ways>1 && showStdev){
				tsw.print(Tools.format("%d\t%s\t%.2f%%\n",k,avgs,stdev));
			}else{
				tsw.print(Tools.format("%d\t%s\n",k,avgs));
			}
				
		}
		if(ffout.stdio()){tsw.println();}
		errorState|=tsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final int tid_){
			cris=cris_;
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
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

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
					readsProcessedT+=r1.pairCount();
					basesProcessedT+=initialLength1+initialLength2;
					
					{
						//Reads are processed in this block.
						processReadPair(r1, r2);
					}
				}

				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 */
		void processReadPair(final Read r1, final Read r2){
			for(MultiLogLog mlog : mlogArray){
				mlog.hash(r1);
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Thread ID */
		final int tid;
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

	/** Primary output file path */
	private String out=null;
	
	/** Override input file extension */
	private String extin=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/** Number of hash functions to use */
	int ways=1;
	
	boolean showStdev=false;
	boolean useWavg=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output file */
	private final FileFormat ffout;
	
	private final MultiLogLog[] mlogArray;
	
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
