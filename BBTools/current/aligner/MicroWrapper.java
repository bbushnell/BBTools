package aligner;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBDuk;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.ReadStreamByteWriter;
import stream.SamLine;
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
public class MicroWrapper implements Accumulator<MicroWrapper.ProcessThread> {
	
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
		MicroWrapper x=new MicroWrapper(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public MicroWrapper(String[] args){
		
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
			extin=parser.extin;

			out1=parser.out1;
			out2=parser.out2;
			extout=parser.extout;
		}
		
		makeReadStats=ReadStats.collectingStats();

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.SAM, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.SAM, extout, true, overwrite, append, ordered);
		ffoutu1=FileFormat.testOutput(outu1, FileFormat.SAM, extout, true, overwrite, append, ordered);
		ffoutu2=FileFormat.testOutput(outu2, FileFormat.SAM, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		samOutput=(ffout1!=null && ffout1.samOrBam()) || (ffoutu1!=null && ffoutu1.samOrBam());
		if(samOutput) {ReadStreamByteWriter.USE_ATTACHED_SAMLINE=true;}
		fixK();
		ref=SideChannel3.fixRefPath(ref);
		minIdentity1=SideChannel3.fixID(minIdentity1);
		minIdentity2=SideChannel3.fixID(minIdentity2);
		assert(k1>0);
		
		final Read r=MicroIndex3.loadRef(ref, samOutput);
		index1=new MicroIndex3(k1, mm1, r);
		index2=(k2<1 ? null : new MicroIndex3(k2, mm1, r));
	}
	
	void fixK() {
		if(k1<k2) {
			int temp=k1;
			k1=k2;
			k2=temp;
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
			}else if(a.equals("mappedonly")){
				mappedOnly=Parse.parseBoolean(b);
			}
			
			else if(a.equals("outu")){
				outu1=b;
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("k") || a.equals("k1") || a.equals("kbig")){
				k1=Integer.parseInt(b);
			}else if(a.equals("k2") || a.equals("ksmall")){
				k2=Integer.parseInt(b);
			}else if(a.equals("minid") || a.equals("minid1")){
				minIdentity1=Float.parseFloat(b);
			}else if(a.equals("minid2")){
				minIdentity2=Float.parseFloat(b);
			}else if(a.equals("mm")){
				mm1=mm2=Integer.parseInt(b);
			}else if(a.equals("mm1")){
				mm1=Integer.parseInt(b);
			}else if(a.equals("mm2")){
				mm2=Integer.parseInt(b);
			}
			
			else if(a.equals("parse_flag_goes_here")){
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

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}

		//Do output file # replacement
		if(outu1!=null && outu2==null && outu1.indexOf('#')>-1){
			outu2=outu1.replace("#", "2");
			outu1=outu1.replace("#", "1");
		}
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
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
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
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
//		assert(false) : "TODO";
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
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros(ffout1, ffout2, cris.paired());
		final ConcurrentReadOutputStream rosu=makeCros(ffoutu1, ffoutu2, cris.paired());
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads(cris, ros, rosu);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros, rosu);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		outstream.println(stats(readsProcessed, basesProcessed));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private ConcurrentReadInputStream makeCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}
	
	private ConcurrentReadOutputStream makeCros(FileFormat ff1, FileFormat ff2, boolean pairedInput){
		if(ff1==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

//		//Notify user of output mode
//		if(pairedInput && out2==null && (in1!=null && !ffin1.samOrBam() && !ff1.samOrBam())){
//			outstream.println("Writing interleaved.");
//		}

		final ConcurrentReadOutputStream ros=
				ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, 
			final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosu){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, rosu, i));
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
			
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			
			readsMapped+=pt.readsMappedT;
			basesMapped+=pt.basesMappedT;
			
			identitySum+=pt.identitySumT;
			errorState|=(!pt.success);
		}
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public String stats(long readsIn, long basesIn) {
		long ro=readsOut, ri=readsIn; 
		long bo=basesOut, bi=basesIn; 
		long rm=readsMapped, bm=basesMapped;
		double idsum=identitySum;
		ByteBuilder bb=new ByteBuilder();
		bb.append("Aligned Reads:     \t"+rm+"    \t"+BBDuk.toPercent(rm, ri)).nl();
		bb.append("Aligned Bases:     \t"+bm+"  \t"+BBDuk.toPercent(bm, bi)).nl();
		bb.append("Avg Identity:      \t"+Tools.format("%.4f%%", idsum*100.0/(rm)));
		return bb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_, 
				final ConcurrentReadOutputStream rosu_, final int tid_){
			cris=cris_;
			ros=ros_;
			rosu=rosu_;
			tid=tid_;
			readstats=makeReadStats ? new ReadStats() : null;
			mapper1=new MicroAligner3(index1, minIdentity1, false);
			mapper2=(k2<1 ? null : new MicroAligner3(index2, minIdentity2, false));
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
			final ArrayList<Read> ulist=new ArrayList<Read>(ln.size());
			
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
					boolean mapped=processReadPair(r1, r2);
					
					
					if(mapped || !mappedOnly){
						readsOutT+=r1.pairCount();
						basesOutT+=r1.pairLength();
					}
					if(!mapped) {
						ulist.add(r1);
						if(mappedOnly) {reads.set(idx, null);}
					}

					boolean output=(!mapped && rosu!=null) || (ros!=null && (mapped || !mappedOnly));
					if(samOutput && output) {
						r1.samline=(r1.samline!=null ? r1.samline : new SamLine(r1, 0));
						if(r2!=null) {r2.samline=(r2.samline!=null ? r2.samline : new SamLine(r2, 1));}
					}
				}
			}

			//Output reads to the output stream
			if(ros!=null){ros.add(reads, ln.id);}
			if(rosu!=null){rosu.add(ulist, ln.id);}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		boolean processReadPair(final Read r1, final Read r2){
			
			boolean mapped=map(r1, r2);

			if(readstats!=null){
				readstats.addToHistograms(r1);
			}
			return mapped;
		}
		/*--------------------------------------------------------------*/
		
		public boolean map(Read r1, Read r2) {
			float id1=mapper1.map(r1);
			float id2=mapper1.map(r2);
			if(id1+id2<=0) {return false;}//Common case
			
			if(r2!=null) {
				if(mapper2!=null) {
					if(r1.mapped() && !r2.mapped()) {id2=mapper2.map(r2);}
					else if(r2.mapped() && !r1.mapped()) {id1=mapper2.map(r1);}
				}
				boolean properPair=(r1.mapped() && r2.mapped() && r1.chrom==r2.chrom && 
						r1.strand()!=r2.strand() && Tools.absdif(r1.start, r2.start)<=1000);
				r1.setPaired(properPair);
				r2.setPaired(properPair);
			}
			
			if(r1.mapped()) {
				readsMappedT++;
				basesMappedT+=r1.length();
				identitySumT+=id1;
			}else {id1=0;}
			
			if(r2!=null && r2.mapped()) {
				readsMappedT++;
				basesMappedT+=r2.length();
				identitySumT+=id2;
			}else {id2=0;}
			
			return id1>0 || id2>0;
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;

		/** Number of reads mapped by this thread */
		protected long readsMappedT=0;
		/** Number of bases mapped by this thread */
		protected long basesMappedT=0;
		
		protected double identitySumT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		final ReadStats readstats;
		final MicroAligner3 mapper1;
		final MicroAligner3 mapper2;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Shared unmapped output stream */
		private final ConcurrentReadOutputStream rosu;
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

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;

	/** Primary output file path */
	private String outu1=null;
	/** Secondary output file path */
	private String outu2=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	/** Number of reads retained */
	protected long readsMapped=0;
	/** Number of bases retained */
	protected long basesMapped=0;
	
	protected double identitySum=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	public int k1=17;
	public int k2=13;
	public float minIdentity1=0.66f;
	public float minIdentity2=0.56f;
	public int mm1=1;
	public int mm2=1;
	public final MicroIndex3 index1;
	public final MicroIndex3 index2;
	
	public String ref;
	boolean mappedOnly=false;
	final boolean makeReadStats;
	final boolean samOutput;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean TRACK_STATS=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	
	/** Primary output file */
	private final FileFormat ffoutu1;
	/** Secondary output file */
	private final FileFormat ffoutu2;
	
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
	/** Reads are output in input order */
	private boolean ordered=false;
	
}
