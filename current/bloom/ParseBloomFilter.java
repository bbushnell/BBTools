package bloom;

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
 * Extracts a compact table from the historical BloomFilterWrapper log format.
 * Hash-prefixed labels and execution lines start new output lines; thread count,
 * counted keys, increments and filter-creation time are appended in encounter
 * order, separated by tabs. Unrecognized nonempty lines optionally go to a
 * separate file. This is a format-specific summarizer, not a general log parser.
 * The complete summary is buffered until input ends; empty input lines are ignored.
 *
 * @author Brian Bushnell
 * @date May 9, 2016
 */
public class ParseBloomFilter{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Program entry point.
	 * Creates an instance of ParseBloomFilter, runs the processing pipeline,
	 * and handles output stream cleanup.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		
		ParseBloomFilter x=new ParseBloomFilter(args);
		
		x.process(t);
		
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command-line arguments and initializes file formats.
	 * Sets up input/output files, validates file existence, and configures
	 * compression settings and file access modes.
	 * @param args Command-line arguments containing input/output file paths and options
	 */
	public ParseBloomFilter(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
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
	
	/**
	 * Parses command-line arguments into configuration settings.
	 * Handles options for invalid output file, maximum lines to process,
	 * and verbose mode settings.
	 *
	 * @param args Command-line argument array
	 * @return Parser object containing parsed configuration
	 */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=", 2);//Values may contain '='; split only the flag delimiter.
			String a=split[0].toLowerCase();
			String b=(split.length>1 && !split[1].isEmpty() ? split[1] : null);
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
			}
		}
		
		return parser;
	}
	
	/**
	 * Adds or removes compression extensions (.gz, .bz2) as needed for input files.
	 * Ensures at least one input file is specified.
	 * @throws RuntimeException if no input file is provided
	 */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/**
	 * Checks input and both output paths before any processing stream is opened.
	 * Invalid-line output obeys the same overwrite and alias rules as the summary.
	 * @throws RuntimeException if files cannot be accessed or are duplicated
	 */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, outInvalid)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+outInvalid+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1, outInvalid)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/**
	 * Selects ByteFile2 when more than two threads are available and neither backend
	 * is already forced. This process-wide override is not restored after processing.
	 */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Opens owned streams, processes input, drains outputs and reports counters.
	 * Normal completion folds close/write errors into errorState and throws if set.
	 * Exceptions during processing currently bypass the normal close sequence.
	 * @param t Started timer to stop after stream closure */
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		ByteStreamWriter bsw=makeBSW(ffout1);
		ByteStreamWriter bswInvalid=makeBSW(ffoutInvalid);
		
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
	 * Builds one summary in memory and copies unrecognized nonempty lines separately.
	 * Preserves the historical prefix matching, field positions and trailing tabs.
	 * Counts input bytes as line length plus one, independent of the original newline.
	 * @param bf Borrowed input positioned at the first unconsumed line
	 * @param bsw Started summary writer; must be non-null
	 * @param bswInvalid Optional started writer for unrecognized lines
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
				boolean valid=true;
				
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
					valid=false;
					if(bswInvalid!=null){
						bswInvalid.println(line);
					}
				}
				if(valid){linesOut++;}
			}
			line=bf.nextLine();
		}

		bb.nl();
		bytesOut+=bb.length;
		bsw.print(bb.toBytes());
		bb.clear();
	}
	
	/** Starts a writer for the supplied descriptor, or returns null for no output. */
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Required input log path after extension resolution. */
	private String in1=null;
	/** Summary destination; defaults to stdout.txt. */
	private String out1=null;
	/** Optional destination for unrecognized nonempty lines. */
	private String outInvalid=null;
	
	/*--------------------------------------------------------------*/
	
	/** Nonempty input lines processed, excluding a line beyond the limit. */
	private long linesProcessed=0;
	/** Recognized nonempty input lines contributing to the summary. */
	private long linesOut=0;
	/** Input line bytes plus one assumed newline per processed nonempty line. */
	private long bytesProcessed=0;
	/** Summary bytes submitted to the writer, excluding invalid-line output. */
	private long bytesOut=0;
	
	/** Nonempty-line limit; current loop treats zero as unlimited too. */
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Resolved input descriptor. */
	private final FileFormat ffin1;
	/** Resolved summary descriptor. */
	private final FileFormat ffout1;
	/** Optional invalid-line descriptor. */
	private final FileFormat ffoutInvalid;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Diagnostic stream selected by PreParser. */
	private PrintStream outstream=System.err;
	/** Shared diagnostic setting, also forwarded to the byte I/O classes. */
	public static boolean verbose=false;
	/** Sticky input/output error reported on the normal close path. */
	public boolean errorState=false;
	/** Allow replacement of existing outputs. */
	private boolean overwrite=true;
	/** Append to outputs when requested by Parser. */
	private boolean append=false;
	
}
