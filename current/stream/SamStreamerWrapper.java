package stream;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.bam.BamIndexWriter;
import structures.ListNum;
import var2.SamFilter;
import var2.ScafMap;

/**
 * Streams SAM/BAM files with optional filtering and format conversion.
 * 
 * Supports:
 * - SAM/BAM to FASTQ/FASTA conversion
 * - SAM/BAM to SAM/BAM conversion with filtering
 * - CIGAR normalization for different SAM versions
 * - SamFilter for quality/mapping filters
 * - Reference-based coordinate lookups
 * 
 * Usage examples:
 *   SamStreamerWrapper in=mapped.bam out=reads.fq.gz
 *   SamStreamerWrapper in=input.sam out=filtered.bam mapped=t
 *   SamStreamerWrapper in=aligned.bam out=output.sam samversion=1.4
 * 
 * @author Brian Bushnell, Isla
 * @date November 6, 2025
 */
public class SamStreamerWrapper{

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
		SamStreamerWrapper x=new SamStreamerWrapper(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	SamStreamerWrapper(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null/*getClass()*/, false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();

			in1=parser.in1;
			out1=parser.out1;
		}

		//Do input/output setup
		fixExtensions();
		checkFileExistence();

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.SAM, null, true, true);
		ffout1=FileFormat.testOutput(out1, FileFormat.SAM, null, true, true, false, true);

		//Determine if we need to parse SAM fields or can skip for performance
		if(!forceParse && !fixCigar && (ffout1==null || !ffout1.samOrBam())){
			SamLine.PARSE_2=false;
			SamLine.PARSE_5=false;
			SamLine.PARSE_6=false;
			SamLine.PARSE_7=false;
			SamLine.PARSE_8=false;
			SamLine.PARSE_OPTIONAL=false;
		}

		//Enable optimizations for SAM->SAM conversion
		ReadStreamByteWriter.USE_ATTACHED_SAMLINE=true;
	}

	/** 
	 * Parse arguments from the command line.
	 * @param args Command line arguments
	 * @return Parser object with standard flags processed
	 */
	private Parser parse(String[] args){

		//Create filter
		filter=new SamFilter();
		filter.includeNonPrimary=true;
		filter.includeLengthZero=true;
		boolean doFilter=true;

		//Create a parser object
		Parser parser=new Parser();

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
			}else if(a.equals("forceparse")){
				forceParse=Parse.parseBoolean(b);
			}else if(a.equals("ref") || a.equals("reference")){
				ref=b;
			}else if(a.equals("rnameasbytes")){
				SamLine.RNAME_AS_BYTES=Parse.parseBoolean(b);
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("samversion") || a.equals("samv") || a.equals("sam")){
				Parser.parseSam(arg, a, b);
				fixCigar=true;
			}else if(a.equals("threadsin") || a.equals("tin")){
				threadsIn=Integer.parseInt(b);
			}else if(a.equals("threadsout") || a.equals("tout")){
				threadsOut=Integer.parseInt(b);

			//Filter parameters
			}else if(a.equals("filter")){
				doFilter=Parse.parseBoolean(b);
			}else if(filter.parse(arg, a, b)){
				//do nothing

			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(i==0 && !arg.contains("=") && parser.in1==null &&
					FileFormat.isSamOrBamFile(arg) && new File(arg).isFile()){
				parser.in1=arg;
			}else if(i==1 && !arg.contains("=") && parser.out1==null && parser.in1!=null &&
					(FileFormat.isSequence(arg) || FileFormat.isBaiFile(arg))){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		//Disable filter if requested
		if(!doFilter){filter=null;}

		return parser;
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Add or remove .gz or .bz2 extensions as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
	}

	/** Ensure input files can be read and output files can be written */
	private void checkFileExistence(){
		//Ensure input file exists
		if(in1==null){
			throw new RuntimeException("Error - at least one input file is required.");
		}

		//Ensure output files can be written
		if(out1!=null && !Tools.testOutputFiles(true, false, false, out1)){
			throw new RuntimeException("\nCan't write to output file "+out1+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read input file "+in1+"\n");
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Primary Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Create Streamer and Writer, then process all data.
	 * @param t Timer for tracking elapsed time
	 */
	void process(Timer t){

		//Determine processing mode
		final boolean inputSam=(ffin1!=null && ffin1.samOrBam());
		final boolean outputSam=(ffout1!=null && ffout1.samOrBam());
		final boolean outputReads=(ffout1!=null && !ffout1.samOrBam());
		final boolean outputBai=(ffout1!=null && ffout1.bai());
		final boolean useSharedHeader=inputSam && outputSam;
		final boolean makeReads=(outputReads || !inputSam);
		if(!inputSam) {
			System.err.println("Input is "+ffin1.formatString()+"; sam filter disabled.");
			filter=null;
			ref=null;
		}
		
		if(outputBai){
			assert(ffin1.bam()) : "bai output requires bam input.";
			try{
				BamIndexWriter.writeIndex(in1, out1);
			}catch(Throwable e){
				KillSwitch.exceptionKill(e);
			}
			t.stop();
			outstream.println("Time:                         \t"+t);
			return;
		}

		//Load reference if specified
		if(ref!=null){
			ScafMap.loadReference(ref, true);
			SamLine.RNAME_AS_BYTES=false;
		}

		//Create streamer and writer
		Streamer st=StreamerFactory.makeStreamer(ffin1, null, ordered, maxReads, useSharedHeader, makeReads, threadsIn);
		Writer fw=(ffout1==null ? null : WriterFactory.makeWriter(ffout1, null, threadsOut, null, useSharedHeader));

		//Process data
		st.start();
		if(fw!=null){fw.start();}

		if(outputReads || !inputSam){
			processAsReads(st, fw);
		}else{
			processAsSam(st, fw);
		}

		//Wait for writer to finish
		if(fw!=null){
			errorState|=fw.poisonAndWait();
			readsOut=fw.readsWritten();
			basesOut=fw.basesWritten();
		}

		//Check for errors
		errorState|=st.errorState();

		//Print statistics
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+String.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		outstream.println("Bases Processed:    "+basesProcessed+" \t"+String.format("%.2f Mbp/sec", (basesProcessed/(double)(t.elapsed))*1000));
		if(ffout1!=null){
			outstream.println("Reads Out:          "+readsOut);
			outstream.println("Bases Out:          "+basesOut);
		}

		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getSimpleName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Process data as Read objects (for SAM/BAM -> FASTQ/FASTA conversion).
	 * @param st Input Streamer
	 * @param fw Output Writer (may be null)
	 */
	private void processAsReads(Streamer st, Writer fw){
		for(ListNum<Read> ln=st.nextList(); ln!=null && ln.size()>0; ln=st.nextList()){
			ArrayList<Read> list=ln.list;
			if(verbose){outstream.println("Got list of size "+ln.size());}

			ArrayList<Read> out=(filter==null ? list : new ArrayList<Read>(list.size()));

			for(Read r : list){
				final int len=r.length();
				readsProcessed++;
				basesProcessed+=len;

				//Apply filter if present
				boolean keep=(filter==null || filter.passesFilter(r.samline));
				if(keep && filter!=null){
					out.add(r);
				}
			}

			if(fw!=null){fw.addReads(new ListNum<Read>(out, ln.id));}
		}
		if(verbose){outstream.println("Finished.");}
	}

	/**
	 * Process data as SamLine objects (for SAM/BAM -> SAM/BAM conversion).
	 * @param st Input Streamer
	 * @param fw Output Writer (may be null)
	 */
	private void processAsSam(Streamer st, Writer fw){
		for(ListNum<SamLine> ln=st.nextLines(); ln!=null && ln.size()>0; ln=st.nextLines()){
			ArrayList<SamLine> list=ln.list;
			if(verbose){outstream.println("Got list of size "+ln.size());}

			ArrayList<SamLine> out=(filter==null && !fixCigar ? list : new ArrayList<SamLine>(list.size()));

			for(SamLine sl : list){
				final int len=sl.lengthOrZero();
				readsProcessed++;
				basesProcessed+=len;

				//Apply filter if present
				boolean keep=(filter==null || filter.passesFilter(sl)) && 
					(len>0 || ffout1==null || ffout1.samOrBam());
				if(!keep){continue;}

				//Fix CIGAR if needed
				if(fixCigar && sl.cigar!=null){
					if(SamLine.VERSION==1.3f){
						sl.cigar=SamLine.toCigar13(sl.cigar);
					}else{
						byte[] shortMatch=sl.toShortMatch(true);
						byte[] longMatch=Read.toLongMatchString(shortMatch);
						int start=sl.pos-1;
						int stop=start+Read.calcMatchLength(longMatch)-1;
						sl.cigar=SamLine.toCigar14(longMatch, start, stop, Integer.MAX_VALUE, sl.seq);
					}
				}

				if(filter!=null || fixCigar){
					out.add(sl);
				}
			}

			if(fw!=null){fw.addLines(new ListNum<SamLine>(out, ln.id));}
		}
		if(verbose){outstream.println("Finished.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** SAM/BAM filter for quality/mapping criteria */
	private SamFilter filter;

	/** Primary input file path */
	private String in1=null;
	/** Primary output file path */
	private String out1=null;
	/** Reference genome path for coordinate lookups */
	private String ref=null;

	/** Primary input file format */
	private FileFormat ffin1;
	/** Primary output file format */
	private FileFormat ffout1;

	/** Number of threads for input streaming (-1 = auto) */
	private int threadsIn=-1;
	/** Number of threads for output writing (-1 = auto) */
	private int threadsOut=-1;

	/** Number of reads processed */
	private long readsProcessed=0;
	/** Number of reads written */
	private long readsOut=0;
	/** Number of bases processed */
	private long basesProcessed=0;
	/** Number of bases written */
	private long basesOut=0;

	/*--------------------------------------------------------------*/

	/** True if an error was encountered */
	public boolean errorState=false;
	/** Maintain input order in output */
	public boolean ordered=true;
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	/** Force parsing of all SAM fields even if not needed */
	private boolean forceParse;
	/** Normalize CIGAR strings to specified SAM version */
	private boolean fixCigar;

	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;

}