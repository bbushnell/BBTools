package consensus;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;
import var2.SamFilter;

/**
 * Resizes scaffold gaps to represent the best estimate 
 * based on the insert size distribution of paired reads.
 * 
 * @author Brian Bushnell
 * @date September 11, 2019
 *
 */
public class FixScaffoldGaps implements Accumulator<FixScaffoldGaps.ProcessThread> {
	
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
		FixScaffoldGaps x=new FixScaffoldGaps(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public FixScaffoldGaps(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		samFilter.includeUnmapped=false;
		samFilter.includeSupplementary=false;
//		samFilter.includeDuplicate=false;
		samFilter.includeNonPrimary=false;
		samFilter.includeQfail=false;
//		samFilter.minMapq=4;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in=parser.in1;
			extin=parser.extin;

			out=parser.out1;
			extout=parser.extout;
		}
		
		{
//			if("auto".equalsIgnoreCase(atomic)){Scaffold.setCA3A(Shared.threads()>8);}
//			else{Scaffold.setCA3A(Parse.parseBoolean(atomic));}
			samFilter.setSamtoolsFilter();
			
			streamerThreads=Tools.max(1, Tools.min(streamerThreads, Shared.threads()));
			assert(streamerThreads>0) : streamerThreads;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout=FileFormat.testOutput(out, FileFormat.FASTA, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin=FileFormat.testInput(in, FileFormat.SAM, extin, true, true);
		ffref=FileFormat.testInput(ref, FileFormat.FASTA, null, true, true);
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
			}else if(a.equals("ref") || a.equals("scaffolds")){
				ref=b;
			}else if(a.equals("insertlist")){
				insertList=b;
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("ns") || a.equalsIgnoreCase("n") || a.equalsIgnoreCase("scaffoldbreak") || a.equalsIgnoreCase("gap")){
				scaffoldBreakNs=Integer.parseInt(b);
				assert(scaffoldBreakNs>0);
			}else if(a.equalsIgnoreCase("mindepth")){
				minDepth=Integer.parseInt(b);
				assert(minDepth>0);
			}else if(a.equalsIgnoreCase("trim") || a.equalsIgnoreCase("trimFraction") || a.equalsIgnoreCase("border")){
				trimFraction=Float.parseFloat(b);
				assert(trimFraction>=0 && trimFraction<=1) : "trimFraction should be between 0 and 1";
			}else if(a.equals("clearfilters") || a.equals("clearfilter")){
				if(Parse.parseBoolean(b)){
					samFilter.clear();
				}
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(samFilter.parse(arg, a, b)){
				//do nothing
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in=Tools.fixExtension(in);
		ref=Tools.fixExtension(ref);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){

		//Ensure there is an input file
		if(in==null){throw new RuntimeException("Error - an input file is required.");}

		//Ensure there is an input file
		if(ref==null){throw new RuntimeException("Error - a reference file is required.");}
		
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
		if(!Tools.testForDuplicateFiles(true, in, ref, out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
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
		
		//Create a read input stream
		final Streamer ss=makeStreamer(ffin);
		
		//Load reference
		loadReferenceCustom();
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads(ss);

		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros();
		
		if(verbose){outstream.println("Fixing reference.");}
		
		fixScaffolds(ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStream(ros);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, scaffoldsOut, scaffoldLengthOut, 8, false));
		
		outstream.println();
		outstream.println(Tools.number("Average Insert", totalAverageInsert, 2, 8));
		outstream.println(Tools.number("Gaps Unchanged", gapsUnchanged, 8));
		outstream.println(Tools.number("Gaps Widened  ", gapsWidened, 8));
		outstream.println(Tools.number("Gaps Narrowed ", gapsNarrowed, 8));
		outstream.println(Tools.number("Ns Total      ", nsTotal, 8));
		outstream.println(Tools.number("Ns Added      ", nsAdded, 8));
		outstream.println(Tools.number("Ns Removed    ", nsRemoved, 8));
		
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/** Loads reference sequences into scaffold maps.
	 * Creates Scaffold objects and populates both full and trimmed name maps. */
	private synchronized void loadReferenceCustom(){
		assert(!loadedRef);
		ConcurrentReadInputStream cris=makeRefCris();
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()) {
			for(Read r : ln){
				String name=r.id;
				String name2=Tools.trimToWhitespace(r.id);
				Scaffold scaf=new Scaffold(name, r.bases, r.numericID);
				refMap.put(name, scaf);
				refMap2.put(name2, scaf);
			}
		}
		loadedRef=true;
	}
	
	/** Creates input stream for reading reference sequences.
	 * @return Configured ConcurrentReadInputStream for reference file */
	private ConcurrentReadInputStream makeRefCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffref, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		assert(!paired) : "References should not be paired.";
		return cris;
	}
	
	private Streamer makeStreamer(FileFormat ff){
		if(ff==null){return null;}
		Streamer ss=StreamerFactory.makeSamOrBamStreamer(ff, streamerThreads, true, false, maxReads, true);
		ss.start(); //Start the stream
		if(verbose){outstream.println("Started Streamer");}
		return ss;
	}
	
	/** Creates output stream for writing fixed scaffolds.
	 * @return Configured ConcurrentReadOutputStream or null if no output */
	private ConcurrentReadOutputStream makeCros(){
		if(ffout==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final Streamer ss){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(ss, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for threads to finish
		boolean success=ThreadWaiter.waitForThreadsToFinish(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		totalAverageInsert=totalInsertSum/(double)totalInsertCount;
		insertByPercentile=Tools.makeHistogram(insertCounts, buckets);
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		readsOut+=pt.readsOutT;
		basesOut+=pt.basesOutT;
		
		totalInsertSum+=pt.totalInsertSumT;
		totalInsertCount+=pt.totalInsertCountT;
		
		errorState|=(!pt.success);
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Fixes all scaffolds based on collected insert size data.
	 * Iterates through reference scaffolds and applies gap corrections.
	 * @param ros Output stream for writing corrected scaffolds
	 */
	private void fixScaffolds(ConcurrentReadOutputStream ros){
		ByteBuilder bb=new ByteBuilder(1000000);

		ArrayList<Read> list=new ArrayList<Read>(200);
		long num=0;
		long lengthSum=0;
		for(Entry<String, Scaffold> e : refMap.entrySet()){
			Scaffold scaf=e.getValue();
			Read r=scaf.fixScaffold(bb);
			lengthSum+=r.length();
			list.add(r);
			scaffoldsOut++;
			scaffoldLengthOut+=r.length();
			
			if(list.size()>=200 || lengthSum>=100000){
				if(ros!=null){ros.add(list, num);}
				list=new ArrayList<Read>(200);
				num++;
				lengthSum=0;
			}
		}
		if(list.size()>0){
			if(ros!=null){ros.add(list, num);}
		}
	}
	
	/**
	 * Calculates insert size from SAM alignment.
	 * @param sl SAM line representing a paired read alignment
	 * @return Absolute insert size
	 */
	private static int calcInsertSize(SamLine sl) {
		assert(sl.mapped() && sl.pairedOnSameChrom());
		assert(sl.primary());
		assert(!sl.supplementary());
		assert(sl.leftmost());
		
		assert(sl.tlen>0) : sl.tlen+"\n\n"+sl;
		return sl.tlen>0 ? sl.tlen : -sl.tlen;
		
//		final int insertSize;
//		String insertTag=null;
//		if(sl.optional!=null){
//			for(String s : sl.optional){
//				if(s.startsWith("X8:Z:")){
//					insertTag=s;
//					break;
//				}
//			}
//		}
//		if(insertTag!=null){
//			insertSize=Integer.parseInt(insertTag.substring(5));
//		}else{
//			insertSize=sl.tlen;//This is unsafe due to indels.
//			assert(false) : "Reads need insert size tags.";
//		}
//		assert(insertSize>0) : sl;
//		
//		return insertSize;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final Streamer ss_, final int tid_){
			ss=ss_;
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
			
			//Grab and process all lists
			for(ListNum<Read> ln=ss.nextList(); ln!=null; ln=ss.nextList()){
//				if(verbose){outstream.println("Got list of size "+list.size());} //Disabled due to non-static access
				
				processList(ln);
			}
			
		}
		
		/** Processes a list of reads for insert size analysis.
		 * @param ln ListNum containing reads to process */
		void processList(ListNum<Read> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r=reads.get(idx);
				
				//Validate reads in worker threads
				if(!r.validated()){r.validate(true);}

				//Track the initial length for statistics
				final int initialLength=r.length();

				//Increment counters
				readsProcessedT+=r.pairCount();
				basesProcessedT+=initialLength;
				
				processRead(r);
			}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		void processRead(final Read r){
			final SamLine sl=r.samline;
			assert(sl!=null) : sl;
			if(samFilter!=null && !samFilter.passesFilter(sl)){return;}
			
			//sl.nextMapped();
			if(sl.mapped() && sl.pairedOnSameChrom() && sl.properPair() && sl.primary() && !sl.supplementary() && sl.leftmost()){
				final String rname=sl.rnameS();
				Scaffold scaf=refMap.get(rname);
				if(scaf==null){scaf=refMap2.get(Tools.trimToWhitespace(rname));}
				assert(scaf!=null) : "Can't find graph for "+rname;
				
				if(scaf!=null){
					final int insertSize=calcInsertSize(sl);
					insertCounts.incrementAndGet(Tools.mid(0, insertSize, insertCounts.length()));
					scaf.add(sl, insertSize);

					readsOutT++;
					basesOutT+=r.length();

					totalInsertSumT+=insertSize;
					totalInsertCountT++;
				}
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		/** Sum of all insert sizes observed by this thread */
		protected long totalInsertSumT=0;
		/** Count of insert sizes observed by this thread */
		protected long totalInsertCountT=0;
		
		/** Local sum of insert sizes for this thread */
		long insertSum=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final Streamer ss;
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Represents a scaffold sequence with coverage and insert size tracking.
	 * Maintains atomic arrays for thread-safe accumulation of coverage depth and insert sizes. */
	private class Scaffold {
		
		/**
		 * Creates a scaffold with coverage tracking arrays.
		 * @param name_ Scaffold name/identifier
		 * @param bases_ Sequence bases
		 * @param numericID_ Numeric identifier
		 */
		Scaffold(String name_, byte[] bases_, long numericID_){
			name=name_;
			bases=bases_;
			numericID=(int)numericID_;
			depthArray=new AtomicIntegerArray(bases.length);
			insertArray=new AtomicLongArray(bases.length);
		}
		
		/**
		 * Adds coverage and insert size data from an aligned read pair.
		 * Updates depth and insert arrays for positions covered by the read pair.
		 * @param sl SAM line representing the alignment
		 * @param insertSize Insert size of the read pair
		 */
		void add(SamLine sl, int insertSize){
			assert(sl.mapped() && sl.pairedOnSameChrom());
			assert(sl.primary());
			assert(!sl.supplementary());
			assert(sl.leftmost());

//			final int insertSize=calcInsertSize(sl);
			
			int start=sl.pos-1;
			int stop=start+sl.tlen;

			int trim=(int)(sl.length()*trimFraction);
			start+=trim;
			stop-=trim;
			
			for(int i=start; i<stop; i++){
				if(i>=0 && i<bases.length){
					depthArray.incrementAndGet(i);
					insertArray.addAndGet(i, insertSize);
				}
			}
			
		}
		
		/**
		 * Generates corrected scaffold sequence with adjusted gap sizes.
		 * Uses insert size statistics to resize N-gaps based on expected insert sizes.
		 * @param bb ByteBuilder for sequence construction
		 * @return Read containing the corrected scaffold sequence
		 */
		Read fixScaffold(ByteBuilder bb){
			int streak=0;
			bb.clear();
			
			if(insertList!=null){
				for(int i=0; i<bases.length; i++) {
					bb.append(i).tab().append(depthArray.get(i)).tab().append(insertArray.get(i)/(Tools.max(1, depthArray.get(i)))).nl();
				}
				ReadWrite.writeString(bb, insertList, false);
				bb.clear();
			}
			
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				if(b=='N'){
					streak++;
				}else{
					if(streak>=scaffoldBreakNs && i-streak>300 && i<bases.length-300){
						int pivot=i-streak/2-1;
						long depthSum=depthArray.get(pivot);
						long insertSum=insertArray.get(pivot);
						double avgInsert=(insertSum/(double)depthSum);
						
						int avgDepth=((depthArray.get(i-200-streak)+depthArray.get(i+200))/2);
						int percentile=(int)(buckets*Tools.max(0.5, 1.0-depthSum/(double)(avgDepth+depthSum)));
						int insertProxy=insertByPercentile[percentile];
						
//						assert(false) : Arrays.toString(insertByPercentile);
						
						int dif=(int)Math.round(insertProxy-avgInsert);
						int toAdd=Tools.max(scaffoldBreakNs, streak+dif);

//						System.err.println("totalAverageInsert="+(int)totalAverageInsert+", avg="+(int)avgInsert+", dif="+dif);
//						System.err.println("proxy="+insertProxy+", percentile="+percentile+", avgDepth="+(int)avgDepth+", depth="+depthSum);
//						System.err.println("pivot="+pivot+", depthSum="+depthSum+", avg="+(int)avgInsert+", streak="+streak+", dif="+dif+", toAdd="+toAdd);
						
						if(dif>0){
							//TODO:  This is tricky because long-insert reads will be self-selected.
							//Should be compared to average coverage, or nearby coverage, and then consult a size distribution histogram.
							gapsWidened++;
							nsAdded+=dif;
						}else if(dif<0){
							gapsNarrowed++;
							nsRemoved-=dif;
						}else{
							gapsUnchanged++;
						}
						nsTotal+=toAdd;
						for(int j=0; j<toAdd; j++){
							bb.append('N');
						}
					}
					streak=0;
					bb.append(b);
				}
			}
			
			return new Read(bb.toBytes(), null, name, numericID);
		}
		
		/** Numeric identifier for this scaffold */
		final int numericID;
		/** Scaffold name/identifier */
		final String name;
		/** Original scaffold sequence bases */
		final byte[] bases;
		/** Thread-safe array tracking coverage depth at each position */
		final AtomicIntegerArray depthArray;
		/** Thread-safe array tracking cumulative insert sizes at each position */
		final AtomicLongArray insertArray;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in=null;
	/** Secondary input file path */
	private String ref=null;

	/** Primary output file path */
	private String out=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** Path to write insert size list for debugging */
	private String insertList=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	/** Number of scaffolds output */
	protected long scaffoldsOut=0;
	/** Total length of scaffolds output */
	protected long scaffoldLengthOut=0;
	
	/** Count of gaps that were not modified */
	protected long gapsUnchanged=0;
	/** Count of gaps that were made larger */
	protected long gapsWidened=0;
	/** Count of gaps that were made smaller */
	protected long gapsNarrowed=0;
	/** Total number of N bases added to gaps */
	protected long nsAdded=0;
	/** Total number of N bases removed from gaps */
	protected long nsRemoved=0;
	/** Total number of N bases in final output */
	protected long nsTotal=0;
	
	/** Sum of all observed insert sizes */
	protected long totalInsertSum=0;
	/** Count of all observed insert sizes */
	protected long totalInsertCount=0;
	/** Average insert size across all observations */
	protected double totalAverageInsert;
	
	/** Histogram of insert size frequencies for distribution analysis */
	protected AtomicLongArray insertCounts=new AtomicLongArray(20000);
	/** Insert sizes by percentile for gap size estimation */
	protected int[] insertByPercentile;
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	/** Threads dedicated to reading the sam file */
	private int streamerThreads=-1;
	
	/** Flag indicating whether reference has been loaded */
	private boolean loadedRef=false;
	
	/** Minimum number of consecutive Ns to consider a scaffold break */
	private int scaffoldBreakNs=10;
	
	/** Number of buckets for insert size percentile calculation */
	int buckets=1000;
	
	/** Minimum coverage depth required for gap analysis */
	private int minDepth=10;
	
	/** Fraction of read ends to trim when calculating coverage */
	private float trimFraction=0.4f;
	
	/** Filter for SAM records to include in analysis */
	public final SamFilter samFilter=new SamFilter();
	
	/** Uses full ref names */
	public LinkedHashMap<String, Scaffold> refMap=new LinkedHashMap<String, Scaffold>();
	/** Uses truncated ref names */
	public LinkedHashMap<String, Scaffold> refMap2=new LinkedHashMap<String, Scaffold>();;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin;
	/** Secondary input file */
	private final FileFormat ffref;
	
	/** Primary output file */
	private final FileFormat ffout;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	/** Read-write lock for thread synchronization */
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
