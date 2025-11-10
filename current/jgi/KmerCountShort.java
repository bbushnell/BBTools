package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
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
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * This class does nothing.
 * It is designed to be easily modified into a program
 * that processes reads in multiple threads, by
 * filling in the processReadPair method.
 * 
 * @author Brian Bushnell
 * @date November 19, 2015
 *
 */
public class KmerCountShort implements Accumulator<KmerCountShort.ProcessThread> {
	
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
		KmerCountShort x=new KmerCountShort(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerCountShort(String[] args){
		
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
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;
			extin=parser.extin;

			if(outKmers==null) {outKmers=parser.out1;}
		}

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		counts=new long[1<<(2*k)];
		
		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		ffout=FileFormat.testOutput(outKmers, FileFormat.TXT, null, true, overwrite, append, false);
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
			}else if(a.equals("out") || a.equals("out1") || a.equals("outkmers") || a.equals("outk") || a.equals("dump")){
				outKmers=b;
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("headermark") || a.equals("comment")){
				headerMark=(b==null ? "" : b);
			}else if(a.equals("mincount") || a.equals("min")){
				minCount=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("skip") || a.equals("rskip") || a.equals("advance")){
				skip=Integer.parseInt(b);
			}else if(a.equals("reverse")){
				reverse=Parse.parseBoolean(b);
			}else if(a.equals("out") || a.equals("out1") || a.equals("outkmers") || a.equals("outk") || a.equals("dump")){
				outKmers=b;
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
		if(!Tools.testOutputFiles(overwrite, append, false, outKmers)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+outKmers+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, outKmers)){
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
		assert(k>0 && k<=15) : "k="+k;
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
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris();
		
		//Reset counters
		readsProcessed=basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads(cris);
		
		if(ffout!=null) {
			if(ffout.fasta()) {
				printKmersFasta(ffout, counts, k, minCount, rcomp, reverse);
			}else {
				printKmersTSV(ffout, headerMark+"Kmer\tCount", counts, k, minCount, rcomp, reverse);
			}
		}
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
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
	
	private ConcurrentReadInputStream makeCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, i, k));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt) {
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			Tools.add(counts, pt.countsT);
			errorState|=(!pt.success);
		}
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void addKmers(final byte[] bases, final long[] counts, int k){
		if(bases==null || bases.length<k){return;}
		assert(k<=15);
		final int bits=2*k;
		final int mask=~((-1)<<bits);
		
		int kmer=0;
		int len=0;
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				counts[kmer]+=(len>=k ? 1 : 0);
			}else{len=kmer=0;}
			
			//Branchless version; slower
//			len=(x>=0 ? len+1 : 0);
//			counts[kmer]+=(len>=k ? 1 : 0);
		}
//		assert(false) : bits+", "+mask+", "+Arrays.toString(counts);
	}
	
	public static void addKmers(final byte[] bases, final long[] counts, final int k, final int skip){
		if(bases==null || bases.length<k){return;}
		assert(k<=15);
		final int bits=2*k;
		final int mask=~((-1)<<bits);
		
		int kmer=0;
		int len=0;
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			
			if(x>=0){
				len++;
				counts[kmer]+=(len>=k && len%skip==0 ? 1 : 0);
			}else{len=kmer=0;}
		}
	}
	
	public static void printKmersTSV(FileFormat ff, String header, long[] counts, int k, int minCount, boolean rcomp, boolean reverse) {
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ff);
		if(header!=null) {bsw.println(header);}
		ByteBuilder bb=new ByteBuilder();
		for(int kmer=0; kmer<counts.length; kmer++) {
			final int rkmer=AminoAcid.reverseComplementBinaryFast(kmer, k);
			final long a=counts[kmer], b=counts[rkmer];
			final long count=(rcomp && kmer!=rkmer) ? a+b : a;
			if(count<minCount || (rcomp && kmer>rkmer)) {continue;}
			bb.appendKmer(reverse ? rkmer : kmer, k).tab();
			bb.append(count).nl();
			bsw.print(bb);
			bb.clear();
		}
		bsw.poison();
	}
	
	public static void printKmersFasta(FileFormat ff, long[] counts, int k, int minCount, boolean rcomp, boolean reverse) {
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ff);
		ByteBuilder bb=new ByteBuilder();
		for(int kmer=0; kmer<counts.length; kmer++) {
			final int rkmer=AminoAcid.reverseComplementBinaryFast(kmer, k);
			final long a=counts[kmer], b=counts[rkmer];
			final long count=(rcomp && kmer!=rkmer) ? a+b : a;
			if(count<minCount || (rcomp && kmer>rkmer)) {continue;}
			bb.append('>').append(count).nl();
			bb.appendKmer(reverse ? rkmer : kmer, k).nl();
			bsw.print(bb);
			bb.clear();
		}
		bsw.poison();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	static class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final int tid_, final int k_){
			cris=cris_;
			tid=tid_;
			k=k_;
		}
		
		//Called by start()
		@Override
		public void run(){
			synchronized(this) {
				//Do anything necessary prior to processing

				//Process the reads
				countsT=new long[1<<(2*k)];
				processInner();

				//Do anything necessary after processing

				//Indicate successful exit status
				success=true;
			}
		}
		
		/** Iterate through the reads */
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
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 */
		void processReadPair(final Read r1, final Read r2){
			if(skip<2) {
				addKmers(r1.bases, countsT, k);
				if(r2!=null) {addKmers(r2.bases, countsT, k);}
			}else {
				addKmers(r1.bases, countsT, k, skip);
				if(r2!=null) {addKmers(r2.bases, countsT, k, skip);}
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		private final int k;
		private long[] countsT;
		
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
	
	/** Kmer count output file */
	private String outKmers="stdout.tsv";
	
	/** Override input file extension */
	private String extin=null;
	
	/*--------------------------------------------------------------*/

	private static int skip=1;
	private int k=4;
	private int minCount=0;
	private boolean rcomp=true;
	private boolean reverse=false;
	private final long[] counts;
	private String headerMark="";
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	private final FileFormat ffout;
	
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
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
