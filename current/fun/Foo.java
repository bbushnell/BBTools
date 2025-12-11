package fun;

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
import structures.LongList;

/**
 * Text file processor that reads delimited data and computes file size statistics.
 * Parses pipe-delimited files containing filesystem metadata including inode numbers,
 * generation numbers, snapshot IDs, file sizes, and file types. Filters files by
 * type and computes statistical distributions of file sizes.
 *
 * Expected input format per line:
 * INODE|GENERATION|SNAPSHOT_ID|FILE_SIZE|MEP_ID|TIMESTAMP|FILE_TYPE|...
 *
 * Only processes files marked with type 'F' and outputs comprehensive size statistics
 * including median, mean, mode, and various percentiles.
 *
 * @author Brian Bushnell
 */
public class Foo {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Program entry point that processes delimited text files and computes statistics.
	 * Creates Foo instance, runs processing with timing, and ensures proper cleanup.
	 * @param args Command line arguments for file paths and processing options
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		Foo x=new Foo(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs Foo processor by parsing arguments and initializing file formats.
	 * Sets up input/output streams, validates parameters, and prepares file formats
	 * for reading input and writing valid/invalid output streams.
	 * @param args Command line arguments including file paths and options
	 */
	public Foo(String[] args){
		
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
	 * Recognizes 'invalid' (junk output path), 'lines' (max lines to process),
	 * 'verbose' (detailed logging), and standard parser arguments.
	 * @param args Array of command line arguments in key=value format
	 * @return Configured Parser instance with file paths and options set
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
	
	/** Adds or removes compression extensions (.gz, .bz2) as needed for input files.
	 * Ensures at least one input file is specified and properly formatted. */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/**
	 * Validates that all input files are readable and output files are writable.
	 * Checks for duplicate file specifications and enforces overwrite/append policies.
	 * Throws RuntimeException if any file access issues are detected.
	 */
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
	 * Adjusts file-related static configuration for optimal performance.
	 * Forces ByteFile2 mode when using more than 2 threads for improved
	 * multi-threaded file reading performance.
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
	 * Currently returns true as placeholder for future parameter validation.
	 * @return Always true in current implementation
	 */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
//		assert(false) : "TODO";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main processing method that creates streams and processes all input data.
	 * Sets up ByteFile reader and ByteStreamWriter outputs, processes data through
	 * inner method, closes streams properly, and reports timing and statistics.
	 * @param t Timer instance for performance measurement
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
	
	/**
	 * Core processing logic that parses pipe-delimited lines and extracts file sizes.
	 * Processes lines with format: INODE|GENERATION|SNAPSHOT_ID|FILE_SIZE|MEP_ID|...
	 * Filters for files marked with type 'F', collects sizes into LongList,
	 * sorts the data, and computes comprehensive statistics including total,
	 * mean, median, mode, and various percentiles (P80, P90, P95, P99).
	 * @param bf ByteFile reader for input data
	 * @param bsw ByteStreamWriter for valid output (currently unused)
	 * @param bswInvalid ByteStreamWriter for invalid data (currently unused)
	 */
	private void processInner(ByteFile bf, ByteStreamWriter bsw, ByteStreamWriter bswInvalid){
		byte[] line=bf.nextLine();
		
		byte delimiter=(byte)'|';
//		SuperLongList sizes=new SuperLongList(10000000);
		LongList sizes=new LongList(100000000);
		
//		INODE           - Specifies the file's inode number
//		GENERATION Number - Specifies a number that is incremented whenever an INODE number is reused.
//		SMNAPSHOT ID    - Specifies the snapshot ID
//		FILE_SIZE       - Specifies the current size or length of the file, in bytes.

		
		while(line!=null){
			linesProcessed++;
			bytesProcessed+=(line.length+1);
			if(line.length>0){
				int a=0, b=0;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing term : '"+new String(line)+"'";
//				long w=Parse.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing term : '"+new String(line)+"'";
//				long y=Parse.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing term : '"+new String(line)+"'";
//				long z=Parse.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing term : '"+new String(line)+"'";
				long size=Parse.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing term : '"+new String(line)+"'";
//				long z=Parse.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing term : '"+new String(line)+"'";
//				long z=Parse.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing term : '"+new String(line)+"'";
				if(line[a]=='F') {
					linesOut++;
					sizes.add(size);
				}
				b++;
				a=b;
			}
			if(linesProcessed>=maxLines) {break;}
			line=bf.nextLine();
		}
		sizes.sort();
		
		long total=sizes.sumLong();
		long mode=sizes.mode();
		long mean=sizes.sumLong()/sizes.size;
		long median=sizes.get((int)(sizes.size*0.5));
		System.out.println("total lines:\t"+Tools.padKMB(linesProcessed, 0)+"\t("+linesProcessed+")");
		System.out.println("total size: \t"+Tools.padKMB(total, 0)+"\t("+total+")"+"\t"+"("+((total/(1024L*1024L*1024L*1024L)))+" tebibytes)");
		System.out.println("median size:\t"+median);
		System.out.println("mean size:  \t"+mean);
		System.out.println("mode size:  \t"+mode);
		System.out.println("P80 size:   \t"+sizes.get((int)(sizes.size*0.8)));
		System.out.println("P90 size:   \t"+sizes.get((int)(sizes.size*0.9)));
		System.out.println("P95 size:   \t"+sizes.get((int)(sizes.size*0.95)));
		System.out.println("P99 size:   \t"+sizes.get((int)(sizes.size*0.99)));
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
	
	/** Input file format specification */
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	private final FileFormat ffoutInvalid;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	/** Allow overwriting existing output files */
	private boolean overwrite=true;
	/** Append to existing output files instead of overwriting */
	private boolean append=false;
	
}
