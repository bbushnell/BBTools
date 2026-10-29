package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashSet;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Filters sequences according to their taxonomy as determined by sequence names.
 * Sequences should be labeled with gi numbers or NCBI taxIDs in their headers.
 * Supports both include and exclude filtering modes based on taxonomic hierarchy.
 *
 * @author Brian Bushnell
 * @date November 23, 2015
 */
public class FilterByTaxa{

	/**
	 * Program entry point from command line.
	 * Creates FilterByTaxa instance and processes input with timing.
	 * @param args Command line arguments for configuration
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		FilterByTaxa x=new FilterByTaxa(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructor that parses command-line arguments and initializes the filter.
	 * Sets up input/output files, creates TaxFilter, and configures processing options.
	 * Handles file validation, format detection, and interleaving settings.
	 * @param args Command-line arguments containing file paths and filter parameters
	 */
	public FilterByTaxa(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		boolean setInterleaved=false; //Whether interleaved was explicitly set.

		//Set shared static variables
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		//Create a parser object
		Parser parser=new Parser();

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("besteffort")){
				bestEffort=Parse.parseBoolean(b);
			}else if(a.equals("results") || a.equals("result")){
				resultsFile=b;
			}else if(TaxFilter.validArgument(a)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		if(resultsFile!=null){nodes=new LinkedHashSet<TaxNode>();}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;

			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;

			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;

			extin=parser.extin;
			extout=parser.extout;
		}

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

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}

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

		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, resultsFile)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2+", "+resultsFile);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, resultsFile)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTA, extout, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTA, extin, true, true);

		//Make the actual filter
		filter=TaxFilter.makeFilter(args);

		//Widen filter to ensure matches
		if(bestEffort){filter.reviseByBestEffort(in1);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates read streams and processes all sequence data through the taxonomic filter.
	 * Sets up concurrent input/output streams, resets input counters only, and reports
	 * accepted reads and bases even when no output stream is requested.
	 * @param t Timer for tracking execution time and performance metrics
	 */
	public void process(Timer t){

		//Create a read input stream
		final Streamer cris;
		{
			cris=StreamerFactory.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2, 1);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}

		//Optionally create a read output stream
		final Writer ros;
		if(ffout1!=null){
			final int buff=4;

			if(cris.paired() && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
				outstream.println("Writing interleaved.");
			}

			ros=WriterFactory.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false, 1);
			ros.start(); //Start the stream
		}else{ros=null;}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Process the read stream
		Throwable processFailure=null;
		try{
			processInner(cris, ros);

			if(verbose){outstream.println("Finished; closing streams.");}

			//Write anything that was accumulated by ReadStats
			errorState|=ReadStats.writeAll();
		}catch(RuntimeException|Error e){
			processFailure=e;
			throw e;
		}finally{
			errorState|=closeProcessStreams(cris, ros, processFailure);
		}

		//Report timing and results
		{
			t.stop();

			String ri="Reads In:               \t"+readsProcessed+" reads";
			String ro="Reads Out:              \t"+readsOut+" reads";
			while(ro.length()<ri.length()){ro=ro+" ";}

			outstream.println(ri+"\t"+basesProcessed+" bases");
			outstream.println(ro+"\t"+basesOut+" bases");
			outstream.println();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}

		//Throw an exception if there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/** Closes input and output independently after processing.
	 * Attempts output cleanup even when input cleanup throws, preserves a primary
	 * processing failure, and suppresses distinct cleanup failures onto it.
	 * @param cris Input read stream to close, may be null
	 * @param ros Output writer to close, may be null
	 * @param primary Primary processing failure, or null when processing finished
	 * @return true if normal cleanup observed an error state
	 */
	private static boolean closeProcessStreams(final Streamer cris, final Writer ros,
			final Throwable primary){
		boolean errorState=false;
		Throwable cleanupFailure=null;
		try{
			errorState|=ReadWrite.closeStream(cris);
		}catch(RuntimeException|Error e){
			cleanupFailure=suppressDistinct(cleanupFailure, e);
		}
		try{
			errorState|=ReadWrite.closeOutputStreams(ros);
		}catch(RuntimeException|Error e){
			cleanupFailure=suppressDistinct(cleanupFailure, e);
		}
		if(cleanupFailure!=null){
			if(primary==null){throwUnchecked(cleanupFailure);}
			suppressDistinct(primary, cleanupFailure);
		}
		return errorState;
	}

	/** Suppresses one throwable onto another unless they are the same object.
	 * @param target Throwable that should receive {@code suppressed}, or null
	 * @param suppressed Throwable to attach, or null
	 * @return {@code target}, or {@code suppressed} when target was null
	 */
	private static Throwable suppressDistinct(final Throwable target, final Throwable suppressed){
		if(target==null){return suppressed;}
		if(suppressed!=null && suppressed!=target){target.addSuppressed(suppressed);}
		return target;
	}

	/** Rethrows a cleanup failure known to be unchecked.
	 * @param failure RuntimeException or Error to rethrow */
	private static void throwUnchecked(final Throwable failure){
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		throw (Error)failure;
	}

	/** Streams list batches, filters reads in place, writes accepted reads when requested,
	 * and records accepted taxonomic nodes when results output is enabled.
	 * @param cris Input read stream
	 * @param ros Optional output read stream */
	void processInner(final Streamer cris, final Writer ros){

		//Do anything necessary prior to processing

		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;

					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());

					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;

					boolean keep=processReadPair(r1, r2);
					if(!keep){reads.set(idx, null);}
					else{
						readsOut+=r1.pairCount();
						basesOut+=initialLength1+initialLength2;
					}
				}

				//Output reads to the output stream
				if(ros!=null){ros.add(reads, ln.id);}

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
		}

		//Do anything necessary after processing

		if(resultsFile!=null){
			TextStreamWriter tsw=new TextStreamWriter(resultsFile, overwrite, append, false);
			tsw.start();
			for(TaxNode tn : nodes){
				tsw.println(tn.id+"\t"+tn.levelStringExtended(false)+"\t"+tn.name);
			}
			errorState|=tsw.poisonAndWait();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Processes one root read through the taxonomic filter and applies that decision
	 * to the whole pair. The mate is accepted or rejected with {@code r1}; its header
	 * is not inspected.
	 * @param r1 Primary read supplying the taxonomy-bearing header
	 * @param r2 Mate read, may be null for unpaired data
	 * @return true if the read or pair passes the taxonomic filter and should be kept
	 */
	boolean processReadPair(final Read r1, final Read r2){
		boolean b=filter.passesFilter(r1.id);
		if(b && nodes!=null){
			TaxNode tn=filter.tree().parseNodeFromHeader(r1.id, true);
			if(tn!=null){nodes.add(tn);}
		}
		return b;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path. */
	private String in1=null;
	/** Secondary input file path for paired reads. */
	private String in2=null;

	/** Primary input quality file path. */
	private String qfin1=null;
	/** Secondary input quality file path. */
	private String qfin2=null;

	/** Primary output file path. */
	private String out1=null;
	/** Secondary output file path. */
	private String out2=null;

	/** Primary output quality file path. */
	private String qfout1=null;
	/** Secondary output quality file path. */
	private String qfout2=null;

	/** Input extension override. */
	private String extin=null;
	/** Output extension override. */
	private String extout=null;

	/** Header-based taxonomy filter. */
	private final TaxFilter filter;

	/** Recur at higher taxonomic levels until some sequence matches; intended for include mode. */
	public boolean bestEffort=false;

	/** Optional path for emitting matched taxonomic nodes. */
	public String resultsFile=null;

	/** Matched taxonomic nodes retained in encounter order for results output. */
	public LinkedHashSet<TaxNode> nodes=null;

	/*--------------------------------------------------------------*/

	/** Number of input reads processed. */
	protected long readsProcessed=0;
	/** Number of input bases processed. */
	protected long basesProcessed=0;

	/** Number of reads accepted by the filter, regardless of whether output is written. */
	public long readsOut=0;
	/** Number of bases accepted by the filter, regardless of whether output is written. */
	public long basesOut=0;

	/** Maximum number of root fragments requested from the Streamer, or -1 for no limit. */
	private long maxReads=-1;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file format handler. */
	private final FileFormat ffin1;
	/** Secondary input file format handler. */
	private final FileFormat ffin2;

	/** Primary output file format handler. */
	private final FileFormat ffout1;
	/** Secondary output file format handler. */
	private final FileFormat ffout2;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print stream for status messages and logging output. */
	private PrintStream outstream=System.err;
	/** Enable verbose status messages during processing. */
	public static boolean verbose=false;
	/** Flag indicating if an error was encountered during processing. */
	public boolean errorState=false;
	/** Allow overwriting existing output files. */
	private boolean overwrite=true;
	/** Append to existing output files instead of overwriting. */
	private boolean append=false;
	/** Unused retained ordering flag. */
	private final boolean ordered=false;

}
