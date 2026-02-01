package var2;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.ReadStreamWriter;
import stream.SamLine;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Realigns SAM records to a reference using Realigner.
 * Takes SAM/BAM input, optionally trims reads, and attempts local realignment
 * against a loaded reference to improve alignment positions and CIGAR strings.
 * @author Brian Bushnell
 * @date April 26, 2017
 */
public class Realign {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		Realign x=new Realign(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Parses command-line arguments, configures SAM filtering, quality trimming,
	 * input/output formats, and loads the reference into ScafMap.
	 * @param args Command line arguments
	 */
	public Realign(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set static variables
		SamLine.PARSE_OPTIONAL_MD_ONLY=true; //I only need the MD tag..
		SamLine.RNAME_AS_BYTES=false;
		SamLine.SET_FROM_OK=true;
		ReadStreamWriter.USE_ATTACHED_SAMLINE=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		Shared.TRIM_READ_COMMENTS=Shared.TRIM_RNAME=true;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.qtrimLeft=qtrimLeft;
		parser.qtrimRight=qtrimRight;
		parser.trimq=trimq;

		samFilter.includeUnmapped=false;
		samFilter.includeSupplementary=false;
		samFilter.includeDuplicate=true;
		samFilter.includeNonPrimary=true;
		samFilter.includeQfail=false;
		samFilter.minMapq=4;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}
			
			
			else if(a.equals("unclip")){
				unclip=Parse.parseBoolean(b);
			}else if(a.equals("realignrows") || a.equals("rerows")){
				Realigner.defaultMaxrows=Integer.parseInt(b);
			}else if(a.equals("realigncols") || a.equals("recols")){
				Realigner.defaultColumns=Integer.parseInt(b);
			}else if(a.equals("realignpadding") || a.equals("repadding") || a.equals("padding")){
				Realigner.defaultPadding=Integer.parseInt(b);
			}else if(a.equals("msa")){
				Realigner.defaultMsaType=b;
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("border")){
				border=Integer.parseInt(b);
			}
			
			else if(samFilter.parse(arg, a, b)){
				//do nothing
			}
			
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		samFilter.setSamtoolsFilter();
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
			
			trimq=parser.trimq;
			trimE=parser.trimE();
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null || ref==null){throw new RuntimeException("Error - one input file and one reference are required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, ref)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, ref, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.SAM, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.SAM, extin, true, true);
		
		loadReference();
	}

	/** Loads the reference genome into ScafMap for realignment and sets the
	 * static Realigner map, skipping work if the reference has already been loaded. */
	private void loadReference(){
		if(loadedRef){return;}
		assert(ref!=null);
		ScafMap.loadReference(ref, scafMap, samFilter, true);
		Realigner.map=scafMap;
		loadedRef=true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates read streams, spawns worker threads, processes all reads, and
	 * prints realignment statistics and timing information.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, null, null);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros;
		if(ffout1!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			//Notify user of output mode
			if(cris.paired() && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
				outstream.println("Writing interleaved.");
			}
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, null, null, buff, null, true);
			ros.start(); //Start the stream
		}else{ros=null;}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads(cris, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println("Realignments:  \t"+realignmentsAttempted);
		outstream.println("Successes:     \t"+realignmentsSucceeded);
		outstream.println("Improvements:  \t"+realignmentsImproved);
		outstream.println("Retained:      \t"+realignmentsRetained);
		outstream.println("Bases trimmed: \t"+basesTrimmed);
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Spawns worker threads to consume reads from the input stream, optionally
	 * write to the output stream, and accumulates per-thread statistics and
	 * realignment counts.
	 * @param cris Concurrent input stream for reading
	 * @param ros Concurrent output stream for writing (may be null)
	 */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, i));
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
			
			trimmedBasesProcessed+=pt.trimmedBasesProcessedT;
			basesTrimmed+=pt.basesTrimmedT;
			readsDiscarded+=pt.readsDiscardedT;
			pairedInSequencingReadsProcessed+=pt.pairedInSequencingReadsProcessedT;
			properlyPairedReadsProcessed+=pt.properlyPairedReadsProcessedT;

			realignmentsAttempted+=pt.realigner.realignmentsAttempted;
			realignmentsImproved+=pt.realigner.realignmentsImproved;
			realignmentsSucceeded+=pt.realigner.realignmentsSucceeded;
			realignmentsRetained+=pt.realigner.realignmentsRetained;
			
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Worker thread that validates reads, invokes processRead(...) for realignment,
	 * and tracks per-thread statistics and trimming counters. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_, final int tid_){
			cris=cris_;
			ros=ros_;
			tid=tid_;
			realigner=new Realigner();
		}
		
		//Called by start()
		/** Main thread body that repeatedly fetches read lists, calls processInner(),
		 * and marks successful completion when processing finishes. */
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
		 * Iterates through lists of reads from the input stream, processes each read
		 * via processRead(...), discards failed reads, and forwards retained reads to
		 * the output stream.
		 */
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
						boolean keep=processRead(r1);
						if(!keep){
							reads.set(idx, null);
							readsDiscardedT++;
						}
					}
				}

				//Output reads to the output stream
				if(ros!=null){ros.add(reads, ln.id);}

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
		 * Processes a read or read pair by applying SAM filters, attempting
		 * realignment with Realigner, trimming with TrimRead if configured, and
		 * updating per-thread counters.
		 * @param r Primary read (mate may be attached)
		 * @return true if the reads should be kept, false if they should be discarded
		 */
		boolean processRead(final Read r){
			if(r.bases==null || r.length()<=1){return false;}
			final SamLine sl=r.samline;
			final int len0=r.length();
			
			if(samFilter!=null && !samFilter.passesFilter(sl)){return false;}
			
			if(sl.properPair()){properlyPairedReadsProcessedT++;}
			if(sl.hasMate()){pairedInSequencingReadsProcessedT++;}
			final Scaffold scaf=scafMap.getScaffold(sl);
			final int scafnum=scaf.number;
			
//			r.toLongMatchString(false); //Not necessary if scoring can be done on short match string
			assert(sl.cigar!=null) : sl;
			boolean realigned=realigner.realign(r, sl, scaf, unclip);
			if(!realigned && r.match!=null){sl.cigar=SamLine.toCigar14(r.match, r.start, r.stop, scaf.length, sl.seq);}
			assert(sl.cigar!=null) : sl;
			
			int leftTrimAmount=border, rightTrimAmount=border;
			if(qtrimLeft || qtrimRight){
				long packed=TrimRead.testOptimal(r.bases, r.quality, trimE);
				if(qtrimLeft){leftTrimAmount=Tools.max(leftTrimAmount, (int)((packed>>32)&0xFFFFFFFFL));}
				if(qtrimRight){rightTrimAmount=Tools.max(rightTrimAmount, (int)((packed)&0xFFFFFFFFL));}
			}
			if(len0-leftTrimAmount-rightTrimAmount<2){return false;}
			
			int trimmed=(leftTrimAmount<1 && rightTrimAmount<1 ? 0 : TrimRead.trimReadWithMatch(r, sl, leftTrimAmount, rightTrimAmount, 0, scaf.length, false));
			if(trimmed<0 || r.length()<2 || sl.cigar==null){return false;}
			basesTrimmedT+=trimmed;
			if(trimmed>0){sl.optional=null;}
			
			assert(sl.cigar!=null) : leftTrimAmount+", "+rightTrimAmount+", "+len0+"\n"+sl;
			return true;
		}
		
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		/** Number of trimmed, mapped bases processed by this thread. */
		protected long trimmedBasesProcessedT=0;
		/** Number of bases trimmed by this thread. */
		protected long basesTrimmedT=0;
		protected long readsDiscardedT=0;
		protected long pairedInSequencingReadsProcessedT=0;
		protected long properlyPairedReadsProcessedT=0;
		
		boolean success=false;
		
		private final ConcurrentReadInputStream cris;
		private final ConcurrentReadOutputStream ros;
		final Realigner realigner;
		
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;

	private String ref=null;
	
	private String out1=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	protected long readsProcessed=0;
	protected long basesProcessed=0;

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	
	private boolean loadedRef=false;

	private boolean qtrimLeft=false;
	private boolean qtrimRight=true;
	private float trimq=10;
	private final float trimE;
	
	/*--------------------------------------------------------------*/
	
	protected long trimmedBasesProcessed=0;
	protected long readsDiscarded=0;
	protected long pairedInSequencingReadsProcessed=0;
	protected long properlyPairedReadsProcessed=0;
	protected long basesTrimmed=0;
	
	protected long realignmentsAttempted;
	protected long realignmentsImproved;
	protected long realignmentsSucceeded;
	protected long realignmentsRetained;
	
	public final ScafMap scafMap=new ScafMap();
	
	public int border=0;
	public boolean unclip=false;
	
	public final SamFilter samFilter=new SamFilter();
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private final FileFormat ffin1;
	
	private final FileFormat ffout1;
	
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
