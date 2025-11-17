package driver;

import java.io.PrintStream;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 17, 2014
 *
 */
public class RenameNcbiToTid {
	
	/**
	 * Program entry point for NCBI to TID header conversion.
	 * Creates a RenameNcbiToTid instance and processes the input file.
	 * @param args Command-line arguments containing input/output file paths
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		RenameNcbiToTid x=new RenameNcbiToTid(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command-line arguments and initializes file formats.
	 * Sets up input and output streams, validates file paths, and configures
	 * processing parameters including overwrite settings and verbosity.
	 *
	 * @param args Command-line arguments containing file paths and options
	 * @throws RuntimeException if input file is not specified or output cannot be written
	 */
	public RenameNcbiToTid(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, true);
	}
	
	/**
	 * Main processing method that reads input file and writes converted output.
	 * Processes each line through the header conversion function and tracks
	 * statistics including lines processed and execution time.
	 *
	 * @param t Timer for tracking execution time
	 * @throws RuntimeException if processing encounters errors or corruption
	 */
	void process(Timer t){
		
		final TextFile tf;
		{
			tf=new TextFile(ffin1);
			if(verbose){outstream.println("Started tf");}
		}
		
		final TextStreamWriter tsw;
		{
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
			if(verbose){outstream.println("Started tsw");}
		}
		
		long linesProcessed=0;
		long charsProcessed=0;
		
		{
			String line;
			while((maxReads<0 || linesProcessed<maxReads) && (line=tf.nextLine())!=null){
				linesProcessed++;
				charsProcessed+=line.length();
				String result=processLine(line);
				if(tsw!=null && result!=null){tsw.println(result);}
			}
		}
		
		errorState|=tsw.poisonAndWait();
		errorState|=tf.close();
		
		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, charsProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	/**
	 * Converts NCBI header format to TID format for a single line.
	 * Changes ">ncbi" prefix to ">tid" and inserts pipe character between
	 * identifier and description. Non-header lines pass through unchanged.
	 *
	 * @param line Input line to process
	 * @return Converted line with TID format, or original line if not an NCBI header
	 */
	private static String processLine(String line){
		if(line.startsWith(">ncbi")){
			line=line.replaceFirst(">ncbi", ">tid");
			int firstSpace=line.indexOf(' ');
			line=line.substring(0, firstSpace)+"|"+line.substring(firstSpace+1);
		}
		return line;
	}
	
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** Input file path */
	private String in1=null;
	/** Output file path */
	private String out1=null;
	
	/*--------------------------------------------------------------*/

	/** Maximum number of lines to process (-1 for unlimited) */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	/** FileFormat for input file handling */
	private final FileFormat ffin1;
	/** FileFormat for output file handling */
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	/** Controls verbosity of status output during processing */
	public static boolean verbose=false;
	/** Tracks whether processing encountered errors */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
