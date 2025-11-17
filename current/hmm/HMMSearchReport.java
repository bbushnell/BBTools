package hmm;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

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

/**
 * Loads output of HMMSearch.
 * @author Brian Bushnell
 * @date April 9, 2020
 *
 */
public class HMMSearchReport {
	
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
		HMMSearchReport x=new HMMSearchReport(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public HMMSearchReport(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, /*getClass()*/null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		
		{//Parse the arguments
			final Parser parser=parse(args);
			parser.out1="stdout.txt";
			overwrite=parser.overwrite;
			append=parser.append;
			
			in=parser.in1;

//			out=parser.out1;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

//		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, append, false);
		ffin=FileFormat.testInput(in, FileFormat.TXT, null, true, true);
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

			if(a.equals("verbose")){
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
		in=Tools.fixExtension(in);
		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
//		if(!Tools.testOutputFiles(overwrite, append, false, out)){
//			outstream.println((out==null)+", "+out);
//			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
//		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
//		//Ensure that no file was specified multiple times
//		if(!Tools.testForDuplicateFiles(true, in, out)){
//			throw new RuntimeException("\nSome file names were specified multiple times.\n");
//		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that executes HMM search result analysis.
	 * Creates input file reader, processes all HMM search lines,
	 * handles cleanup, and reports timing statistics.
	 * @param t Timer for tracking execution performance
	 */
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin);
		ByteStreamWriter bsw=null;//makeBSW(ffout);
		
//		assert(false) : "Header goes here.";
		if(bsw!=null){
//			assert(false) : "Header goes here.";
		}
		
		processInner(bf, bsw);
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
//		outstream.println();
//		outstream.println("Valid Lines:       \t"+linesOut);
//		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesOut));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Internal processing logic for HMM search results.
	 * Loads all HMM search lines and adds them to the protein summary map.
	 * Outputs each processed line to stderr for monitoring.
	 *
	 * @param bf ByteFile reader for input HMM search results
	 * @param bsw ByteStreamWriter for output (currently unused)
	 */
	private void processInner(ByteFile bf, ByteStreamWriter bsw){
		ArrayList<HMMSearchLine> lines=load(bf);
		for(HMMSearchLine line : lines){
			addToMap(line);
			System.err.println(line);
		}
	}
	
	/**
	 * Adds an HMM search line to the protein summary mapping.
	 * Creates new ProteinSummary if none exists for the protein name,
	 * otherwise adds the hit to the existing summary.
	 * @param line HMMSearchLine containing search hit information
	 */
	private void addToMap(HMMSearchLine line){
		ProteinSummary ps=map.get(line.name);
		if(ps==null){
			ps=new ProteinSummary(line.name);
			map.put(line.name, ps);
		}
		ps.add(line);
	}
	
	/**
	 * Loads and parses all HMM search lines from input file.
	 * Skips comment lines starting with '#' and creates HMMSearchLine
	 * objects for each valid result line. Updates processing statistics.
	 *
	 * @param bf ByteFile reader for the input file
	 * @return ArrayList of parsed HMMSearchLine objects
	 */
	private ArrayList<HMMSearchLine> load(ByteFile bf){
		byte[] line=bf.nextLine();
		
		ArrayList<HMMSearchLine> lines=new ArrayList<HMMSearchLine>();
		while(line!=null){
			if(line.length>0){
				linesProcessed++;
				bytesProcessed+=(line.length+1);

				if(line[0]!='#'){
					HMMSearchLine hline=new HMMSearchLine(line);
					lines.add(hline);
				}
			}
			line=bf.nextLine();
		}
		return lines;
	}
	
	/**
	 * Creates and starts a ByteStreamWriter for output formatting.
	 * Returns null if no output format is specified.
	 * @param ff FileFormat specification for output
	 * @return Started ByteStreamWriter or null if no output needed
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
	
	/** Input file path for HMM search results */
	private String in=null;
//	private String out=null;
	
	/** Protein summary mapping from protein names to their best hit summaries */
	public HashMap<String, ProteinSummary> map=new HashMap<String, ProteinSummary>();
	
	/*--------------------------------------------------------------*/
	
	/** Total number of lines processed from input file */
	private long linesProcessed=0;
	/** Total number of bytes processed from input file */
	private long bytesProcessed=0;
//	private long linesOut=0;
//	private long bytesOut=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input file format specification */
	private final FileFormat ffin;
//	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and results */
	private PrintStream outstream=System.err;
	/** Controls verbose output for debugging and detailed logging */
	public static boolean verbose=false;
	/** Tracks whether any errors occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files rather than overwrite */
	private boolean append=false;
	
}
