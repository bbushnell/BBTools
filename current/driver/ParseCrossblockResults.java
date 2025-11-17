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
 * @date June 9, 2016
 *
 */
public class ParseCrossblockResults {
	
	/**
	 * Program entry point for parsing crossblock results.
	 * Creates ParseCrossblockResults instance and processes input file with timing.
	 * @param args Command-line arguments including input file and options
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		ParseCrossblockResults x=new ParseCrossblockResults(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs parser from command-line arguments.
	 * Initializes file formats, validates input/output paths, and configures processing options.
	 * @param args Command-line arguments containing file paths and settings
	 * @throws RuntimeException if required input file is missing or output validation fails
	 */
	public ParseCrossblockResults(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		
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

		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, false, false);
	}
	
	/**
	 * Main processing loop that reads input file and parses each result line.
	 * Skips comment lines starting with '#' and processes data lines through processLine().
	 * Reports timing statistics and throws exception if errors occurred during processing.
	 *
	 * @param t Timer for tracking execution time and reporting performance statistics
	 * @throws RuntimeException if errorState is true after processing
	 */
	void process(Timer t){
		
		final TextFile tf;
		{
			tf=new TextFile(ffin1);
			if(verbose){outstream.println("Started tf");}
		}
		
		long linesProcessed=0;
		long charsProcessed=0;
		
		{
			String line;
			while((maxReads<0 || linesProcessed<maxReads) && (line=tf.nextLine())!=null){
				linesProcessed++;
				charsProcessed+=line.length();
				if(!line.startsWith("#")){
					processLine(line);
				}
			}
		}
		errorState|=tf.close();

		if(ffout1!=null){
			final TextStreamWriter tsw;
			{
				tsw=new TextStreamWriter(ffout1);
				tsw.start();
				if(verbose){outstream.println("Started tsw");}
				errorState|=tsw.poisonAndWait();
			}
		}
		
		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, charsProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	/**
	 * Processes a single tab-separated result line from crossblock output.
	 * Parses line into ResultsLine object and updates statistics based on removal status.
	 * Increments appropriate counters for bases/contigs kept or discarded.
	 * @param line Tab-separated line containing crossblock result data
	 */
	private void processLine(String line){
		ResultsLine rl=new ResultsLine(line);
		if(rl.removed){
			basesDiscarded+=rl.length;
			contigsDiscarded++;
		}else{
			basesKept+=rl.length;
			contigsKept++;
		}
	}
	
	
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Represents a parsed line from crossblock results file.
	 * Contains length and removal status extracted from tab-separated values.
	 * Expected format: field0 \t field1 \t removed_flag \t length_value
	 */
	private static class ResultsLine{
		
		/**
		 * Parses tab-separated result line into structured data.
		 * Extracts length from field 3 and removal status from field 2 (1=removed, 0=kept).
		 * @param s Tab-separated line from crossblock results file
		 */
		public ResultsLine(String s){
			String[] split=s.split("\t");
			length=Integer.parseInt(split[3]);
			removed=Integer.parseInt(split[2])==1;
		}
		
		/** Length of the sequence from field 3 of tab-separated line */
		final int length;
		/** Whether sequence was removed (true if field 2 equals 1) */
		final boolean removed;
		
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** Input file path for crossblock results file */
	private String in1=null;
	/** Output file path (currently unused in processing logic) */
	private String out1=null;
	
	/*--------------------------------------------------------------*/

	/** Maximum number of lines to process (-1 for unlimited) */
	private long maxReads=-1;

	/** Total number of bases in sequences that were kept */
	private long basesKept=0;
	/** Total number of bases in sequences that were removed */
	private long basesDiscarded=0;
	/** Count of sequences that were kept after filtering */
	private long contigsKept=0;
	/** Count of sequences that were removed during filtering */
	private long contigsDiscarded=0;

	/** Gets total bases in kept sequences */
	public long basesKept(){return basesKept;}
	/** Gets total bases in discarded sequences */
	public long basesDiscarded(){return basesDiscarded;}
	/** Gets count of kept sequences */
	public long contigsKept(){return contigsKept;}
	/** Gets count of discarded sequences */
	public long contigsDiscarded(){return contigsDiscarded;}
	
	/** Gets total sequence count (kept + discarded) */
	public long contigs(){return contigsKept+contigsDiscarded;}
	/** Gets total base count (kept + discarded) */
	public long bases(){return basesKept+basesDiscarded;}
	
	/*--------------------------------------------------------------*/
	
	/** Input file format wrapper for reading crossblock results */
	private final FileFormat ffin1;
	/** Output file format wrapper (currently unused) */
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and results */
	private PrintStream outstream=System.err;
	/** Controls verbose output for debugging and detailed logging */
	public static boolean verbose=false;
	/** Tracks whether errors occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files rather than overwriting */
	private boolean append=false;
	
}
