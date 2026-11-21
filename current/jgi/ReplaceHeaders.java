package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Replaces read headers with headers from a separate file.
 *
 * <p>Replacement headers are read from one or two header streams in step with the
 * sequence reads.  Each consumed read keeps its sequence and quality; only its ID
 * changes, either by direct replacement or by prefixing the replacement ID.</p>
 *
 * @author Brian Bushnell
 * @date May 23, 2016
 */
public class ReplaceHeaders{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Program entry point that constructs ReplaceHeaders and runs processing.
	 * @param args Command line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		ReplaceHeaders x=new ReplaceHeaders(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Parses command line arguments, prepares read and header descriptors, validates
	 * path combinations, and records pairing/interleaving choices.
	 *
	 * <p>Shared IO settings are configured here, including the intentionally short
	 * shared read-list length.  The read/header streams themselves are created in
	 * {@link #process(Timer)}.</p>
	 *
	 * @param args Command line arguments specifying input/output files and options
	 */
	public ReplaceHeaders(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Configure shared IO settings; the shared read-list length is intentionally one.
		Shared.setBufferLen(1);//Tools.min(200, Shared.bufferLen());
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Only the first '=' separates the option; filenames may contain more.
			final int equals=arg.indexOf('=');
			String a=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
			String b=(equals<0 || equals==arg.length()-1 ? null : arg.substring(equals+1));
			
			if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(a.equals("hin") || a.equals("hin1") || a.equals("inh") || a.equals("inh1") || a.equals("headerin") || a.equals("headerin1")){
				hin1=b;
			}else if(a.equals("hin2") || a.equals("inh2") || a.equals("headerin2")){
				hin2=b;
			}else if(a.equals("prefix") || a.equals("addprefix")){
				prefix=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Dormant branch reserved for a future local flag.
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
		
		//Do header input file # replacement
		if(hin1!=null && hin2==null && hin1.indexOf('#')>-1 && !new File(hin1).exists()){
			hin2=hin1.replace("#", "2");
			hin1=hin1.replace("#", "1");
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
		
		//Ensure there is a read input file.
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Ensure there is a header input file.
		if(hin1==null){throw new RuntimeException("Error - at least one header input file is required.");}
		
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
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2, hin1, hin2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, hin1, hin2, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		//Create read input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

		//Create header input FileFormat objects
		ffhin1=FileFormat.testInput(hin1, FileFormat.HEADER, null, true, true);
		ffhin2=FileFormat.testInput(hin2, FileFormat.HEADER, null, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates read, header, and optional output streams, delegates read/header ID
	 * replacement, then closes streams on normal completion.
	 *
	 * <p>The borrowed streams passed to {@link #processInner(ConcurrentReadInputStream,
	 * ConcurrentReadInputStream, ConcurrentReadOutputStream)} remain owned here:
	 * after processing, accumulated {@link ReadStats} are written and the read/header
	 * streams and nullable output stream are closed.</p>
	 *
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		//Create the read input stream.
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		final boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Create the replacement-header input stream and check pairedness.
		final ConcurrentReadInputStream hcris;
		{
			hcris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffhin1, ffhin2, null, null);
			hcris.start(); //Start the stream
			boolean paired2=hcris.paired();
			if(verbose){outstream.println("Started hcris");}
			if(paired!=paired2){
				KillSwitch.kill("Reads and heades must agree on pairedness: reads="+paired+", headers="+paired2);
			}
		}
		
		//Optionally create and start the read output stream.
		final ConcurrentReadOutputStream ros;
		if(ffout1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
				outstream.println("Writing interleaved.");
			}
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
			ros.start(); //Start the stream
		}else{ros=null;}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process paired read/header batches.
		processInner(cris, hcris, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write accumulated statistics and close the streams we own.
		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(hcris);
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		//Throw an exception if there was an error in a thread.
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Iterates through read and header streams synchronously for all consumed read
	 * batches.
	 *
	 * <p>The first read/header batches must be equal in size.  Later consumed read
	 * batches require at least as many replacement headers as reads, which permits an
	 * unused tail in a fetched header batch and does not drain additional surplus
	 * header batches once the read stream ends.</p>
	 *
	 * @param cris Input stream for sequence reads
	 * @param hcris Input stream for replacement headers
	 * @param ros Output stream for modified reads (may be null)
	 */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadInputStream hcris, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		{
			//Grab the first ListNum of reads.
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum.
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//Grab the first ListNum of headers.
			ListNum<Read> hln=hcris.nextList();
			//Grab the actual header list from the ListNum.
			ArrayList<Read> hreads=(hln!=null ? hln.list : null);
			
			if(hreads==null || hreads.size()!=reads.size()){
				KillSwitch.kill("Wrong number of headers for the reads: "+reads.size()+"!="+hreads.size());
			}
			
			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				Read h=hreads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
				assert((r.mate!=null) == (h.mate!=null));
			}
			
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//Later read batches need at least the headers they consume.
				if(hreads==null || hreads.size()<reads.size()){
					KillSwitch.kill("Missing replacement headers in batch "+ln.id+": expected "+reads.size()+", found "+(hreads==null ? 0 : hreads.size())+".");
				}
				if(verbose){
					outstream.println("Fetched "+reads.size()+" reads: "+reads.get(0).id+"\n"+reads.get(0).length());
					outstream.println("Fetched "+hreads.size()+" headers.");
				}
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final Read h1=hreads.get(idx);
					final Read h2=h1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;
					
					boolean keep=processReadPair(r1, r2, h1, h2);
					if(!keep){reads.set(idx, null);}
				}
				
				//Output reads to the output stream.
				if(ros!=null){ros.add(reads, ln.id);}
				
				//Return the consumed read list to its stream.
				cris.returnList(ln);
				if(verbose){outstream.println("Returned a read list.");}
				
				//Return the consumed header list to its stream.
				hcris.returnList(hln.id, hln.list.isEmpty());
				if(verbose){outstream.println("Returned a header list.");}
				
				//Fetch the next read list.
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
				
				//Fetch the next header list.
				hln=hcris.nextList();
				hreads=(hln!=null ? hln.list : null);
			}
			
			//Return the final fetched read list, which was not consumed above.
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
			
			//Return the final fetched header list, which may contain unused headers.
			if(hln!=null){
				hcris.returnList(hln.id, hln.list==null || hln.list.isEmpty());
			}
		}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Replaces only the ID fields for one read or read pair.
	 *
	 * <p>Prefix mode joins the replacement ID, a single space, and the original ID.
	 * Non-prefix mode assigns the replacement ID directly.  Sequences, qualities,
	 * pairing, and all other read state are left unchanged.</p>
	 *
	 * @param r1 First read in the pair
	 * @param r2 Second read in the pair (may be null for unpaired)
	 * @param h1 Replacement header for first read
	 * @param h2 Replacement header for second read (may be null for unpaired)
	 * @return true; reads are always kept
	 */
	boolean processReadPair(final Read r1, final Read r2, final Read h1, final Read h2){
		if(prefix){
			r1.id=h1.id+" "+r1.id;
			if(r2!=null){r2.id=h2.id+" "+r2.id;}
		}else{
			r1.id=h1.id;
			if(r2!=null){r2.id=h2.id;}
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** First and optional second read-input paths. */
	private String in1=null;
	private String in2=null;

	/** First and optional second replacement-header paths. */
	private String hin1=null;
	private String hin2=null;
	
	/** Optional quality-sidecar input paths for the read inputs. */
	private String qfin1=null;
	private String qfin2=null;

	/** First and optional second read-output paths. */
	private String out1=null;
	private String out2=null;

	/** Optional quality-sidecar output paths for the read outputs. */
	private String qfout1=null;
	private String qfout2=null;
	
	/** Read input and output extension overrides. */
	private String extin=null;
	private String extout=null;
	
	/** True to prepend replacement IDs to existing IDs rather than replacing them. */
	private boolean prefix=false;
	
	/*--------------------------------------------------------------*/

	/** Read and base counters for consumed reads. */
	protected long readsProcessed=0;
	protected long basesProcessed=0;

	/** Maximum reads or pairs to request from both read and header streams. */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Read input descriptors. */
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	/** Replacement-header input descriptors. */
	private final FileFormat ffhin1;
	private final FileFormat ffhin2;
	
	/** Read output descriptors; the first may be null to suppress output. */
	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Diagnostic output stream, usually stderr unless redirected by PreParser. */
	private PrintStream outstream=System.err;
	/** Enables verbose progress diagnostics. */
	public static boolean verbose=false;
	/** Set when statistics or stream closure reports an error. */
	public boolean errorState=false;
	/** Output overwrite flag. */
	private boolean overwrite=true;
	/** Output append flag. */
	private boolean append=false;
	
}
