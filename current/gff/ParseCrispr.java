package gff;

import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Command-line utility for parsing and extracting CRISPR annotation sequences
 * from GFF files. Processes input GFF files to extract direct repeat CRISPR
 * sequences, converting them into FASTA format with sequential identifiers.
 * Supports filtering valid CRISPR annotations and handling file compression.
 *
 * @author Brian Bushnell
 */
public class ParseCrispr {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Program entry point that initializes a ParseCrispr instance and processes
	 * the specified input files.
	 * @param args Command line arguments including input/output file paths
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		ParseCrispr x=new ParseCrispr(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command line arguments, validates parameters,
	 * and initializes file format objects for input and output processing.
	 * Sets up compression settings and prepares file streams.
	 * @param args Command line arguments for configuration
	 */
	public ParseCrispr(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, /*getClass()*/null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses command line arguments into configuration parameters.
	 * Supports parameters for input files, output files, maximum line limits,
	 * verbose mode, and invalid entry handling.
	 *
	 * @param args Array of command line arguments in key=value format
	 * @return Configured Parser object with processed arguments
	 */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
//		parser.out1="stdout";
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	/** Adjusts file extensions for compression formats (.gz, .bz2) and validates
	 * that required input files are specified. */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Validates that input files can be read and output files can be written.
	 * Checks for duplicate file specifications to prevent data corruption. */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/**
	 * Configures static file processing settings based on available thread count.
	 * Forces ByteFile2 mode for multi-threaded processing when more than 2 threads
	 * are available.
	 */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
//		if(!ByteFile.FORCE_MODE_BF2){
//			ByteFile.FORCE_MODE_BF2=false;
//			ByteFile.FORCE_MODE_BF1=true;
//		}
	}
	
	/**
	 * Validates parameter ranges and ensures required parameters are set.
	 * Currently contains a TODO marker for future parameter validation logic.
	 * @return true if parameters are valid
	 */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
		assert(false) : "TODO";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main processing method that creates input/output streams and processes
	 * all CRISPR data. Handles file I/O setup, calls inner processing method,
	 * and reports statistics including processing time, lines processed,
	 * and valid/invalid line counts.
	 *
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		ByteStreamWriter bsw=makeBSW(ffout1);
		ByteStreamWriter bswInvalid=makeBSW(ffoutInvalid);
		
//		assert(false) : "Header goes here.";
		if(bsw!=null){
//			assert(false) : "Header goes here.";
		}
		
		processInner(bf, bsw, bswInvalid);
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Valid Lines:       \t"+linesOut);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesOut));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	//tid|1000566|NZ_JACGWZ010000001.1	RefSeq	direct_repeat	1135444	1135673	.	+	.	
	//ID=id-NZ_JACGWZ010000001.1:1135444..1135673;gbkey=repeat_region;
	//inference=COORDINATES: alignment:pilercr:v1.02;rpt_family=CRISPR;
	//rpt_type=direct;rpt_unit_range=1135453..1135473;rpt_unit_seq=gcccctacgagggatcgcaac
	
	private void processInner(ByteFile bf, ByteStreamWriter bsw, ByteStreamWriter bswInvalid){
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();
		final int taglen="rpt_unit_seq=".length();
		
		long number=1;
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				
				
				boolean valid=(line[0]!='#' && Tools.indexOf(line, "rpt_family=CRISPR;rpt_type=direct;", 0)>=0);
				
				if(valid){
					int idx=Tools.indexOf(line, "rpt_unit_seq=", 0);
					assert(idx>=0) : new String(line);
					
					bb.clear().append('>').append(number).nl();
					number++;
					for(int i=idx+taglen; i<line.length; i++) {
						byte b=line[i];
						if(!Tools.isLetter(b)) {break;}
						bb.append(Tools.toUpperCase(b));
					}
					bb.nl();
					
					linesOut+=2;
					bytesOut+=(bb.length);
					bsw.print(bb);
					bb.clear();
				}else{
					if(bswInvalid!=null){
						bswInvalid.println(line);
					}
				}
			}
			line=bf.nextLine();
		}
	}
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;

	private String out1=null;

	/** Optional output file path for invalid GFF entries */
	private String outInvalid=null;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input file format configuration for primary input file */
	private final FileFormat ffin1;
	/** Output file format configuration for primary output file */
	private final FileFormat ffout1;
	private final FileFormat ffoutInvalid;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	/** Controls whether to append to existing output files */
	private boolean append=false;
	
}
