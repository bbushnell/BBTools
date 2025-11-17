package bloom;

import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Reads a text file.
 * Prints it to another text file.
 * Filters out invalid lines and prints them to an optional third file.
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class ParseBloomFilter {
	
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
		ParseBloomFilter x=new ParseBloomFilter(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public ParseBloomFilter(String[] args){
		
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

			out1=(parser.out1==null ? "stdout.txt" : parser.out1);
		}
		
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
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
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
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
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
	
	/** Adjust file-related static fields as needed for this program */
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
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that coordinates file reading, parsing, and output.
	 * Opens input/output streams, processes the file content, and reports
	 * execution statistics including valid/invalid line counts and timing.
	 * @param t Timer for tracking execution time and performance metrics
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
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Core parsing logic that processes each line of the input file.
	 * Extracts thread counts from "Executing bloom.BloomFilterWrapper" lines,
	 * parses numeric values from "Keys Counted:" and "Increments:" lines,
	 * and extracts timing from "Filter creation:" lines. Invalid lines
	 * are written to the invalid output stream if specified.
	 *
	 * @param bf Input file reader
	 * @param bsw Output stream for valid parsed data
	 * @param bswInvalid Output stream for invalid lines (may be null)
	 */
	private void processInner(ByteFile bf, ByteStreamWriter bsw, ByteStreamWriter bswInvalid){
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();
		
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				final boolean header=(line[0]=='#');
				String s=new String(line);
				
				if(header){
					bb.nl();
					bb.append(line);
				}else if(s.startsWith("Executing bloom.BloomFilterWrapper")){
					//parse threads
					bb.nl();
					String[] split=s.replace(' ', '\t').replace(',', '\t').split("\t");
					for(String sub : split){
						if(sub.startsWith("t=")){
							bb.append(sub.substring(2)).tab();
							break;
						}
					}
				}else if(s.startsWith("Keys Counted:") || s.startsWith("Increments:")){
					String[] split=Tools.whitespacePlus.split(s);
					bb.append(split[split.length-1]).tab();
				}else if(s.startsWith("Filter creation:")){
					String[] split=Tools.whitespacePlus.split(s);
					bb.append(split[split.length-2]).tab();
				}else{
					if(bswInvalid!=null){
						bswInvalid.println(line);
					}
				}
			}
			line=bf.nextLine();
		}

		bb.nl();
		bytesOut+=bb.length;
		bsw.print(bb.toBytes());
		bb.clear();
	}
	
	/**
	 * Creates and starts a ByteStreamWriter for the given file format.
	 * Returns null if the file format is null.
	 * @param ff File format specification for the output stream
	 * @return Started ByteStreamWriter instance, or null if ff is null
	 */
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input file path */
	private String in1=null;
	/** Output file path for valid parsed lines */
	private String out1=null;
	/** Output file path for invalid lines that couldn't be parsed */
	private String outInvalid=null;
	
	/*--------------------------------------------------------------*/
	
	/** Total number of input lines processed */
	private long linesProcessed=0;
	/** Number of valid output lines written */
	private long linesOut=0;
	/** Total number of input bytes processed */
	private long bytesProcessed=0;
	/** Total number of output bytes written */
	private long bytesOut=0;
	
	/** Maximum number of lines to process before stopping */
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** File format specification for the input file */
	private final FileFormat ffin1;
	/** File format specification for the main output file */
	private final FileFormat ffout1;
	/** File format specification for the invalid lines output file */
	private final FileFormat ffoutInvalid;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and statistics */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging and detailed progress reporting */
	public static boolean verbose=false;
	/** Flag indicating if an error occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
