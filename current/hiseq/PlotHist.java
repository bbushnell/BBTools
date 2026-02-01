package hiseq;

import java.io.PrintStream;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Generates histograms from data matrix files, specifically designed for plotting
 * TileDump data. Processes input text files containing numerical data, creating
 * histogram distributions by binning values across multiple columns.
 *
 * The tool dynamically calculates bin ranges based on the maximum values in each
 * column and generates output TSV files representing the histogram counts. Each
 * column from the input file produces a separate histogram output file.
 *
 * @author Brian Bushnell
 */
public class PlotHist {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		PlotHist x=new PlotHist(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command line arguments and initializes the PlotHist
	 * instance. Sets up file formats, validates parameters, and prepares for
	 * data processing.
	 * @param args Command line arguments
	 */
	public PlotHist(String[] args){
		
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
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses arguments from the command line using a Parser object.
	 * Supports configuration of bins count and verbose output settings.
	 * @param args Command line arguments to parse
	 * @return Configured Parser instance with parsed settings
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

			if(a.equals("bins")){
				bins=Integer.parseInt(b);
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
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	/**
	 * Adds or removes .gz or .bz2 extensions as needed for input files.
	 * Validates that at least one input file is specified.
	 * @throws RuntimeException if no input file is provided
	 */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensures input files can be read by testing file accessibility.
	 * @throws RuntimeException if input files cannot be read */
	private void checkFileExistence(){
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
	}
	
	/**
	 * Adjusts file-related static fields as needed for this program.
	 * Configures ByteFile mode based on available thread count for optimal
	 * I/O performance.
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
	 * Ensures parameter ranges are within bounds and required parameters are set.
	 * Validates that input file and bins count are properly configured.
	 * @return true if all parameters are valid
	 */
	private boolean validateParams(){
		assert(in1!=null);
		assert(bins>0);
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates streams and processes all data through the histogram generation
	 * pipeline. Coordinates file reading, data gathering, histogram calculation,
	 * and output writing.
	 * @param t Timer for tracking processing duration and performance metrics
	 */
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		
		processInner(bf);
		writeData();
		
		errorState|=bf.close();
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Lines Out:         \t"+linesOut);
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void gatherData(ByteFile bf) {
		byte[] line=bf.nextLine();
		
		LineParser1 lp=new LineParser1('\t');
		while(line!=null){
			if(line.length>0){
				final boolean valid=(line[0]!='#');
				if(valid) {
					lp.set(line);
					int lineTerms=Tools.min(terms, lp.terms());
					for(int term=0; term<lineTerms; term++) {
						double d=lp.parseDouble(term);
						maxArray[term]=Tools.max(maxArray[term], d);
					}
				}else {
					header=new String(line).substring(1).split("\t");
					terms=header.length;
					maxArray=new double[terms];
				}
			}
			line=bf.nextLine();
		}
		bf.reset();
		countMatrix=new long[terms][bins+1];
	}
	
	private void processInner(ByteFile bf){
		gatherData(bf);
		
		byte[] line=bf.nextLine();

		LineParser1 lp=new LineParser1('\t');
		while(line!=null){
			if(line.length>0){
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				final boolean valid=(line[0]!='#');
				
				if(valid){
					lp.set(line);
					int lineTerms=Tools.min(terms, lp.terms());
					for(int term=0; term<lineTerms; term++) {
						double d=lp.parseDouble(term);
						int bin=Tools.max(0, (int)((d/maxArray[term])*bins));
						assert(bin<=bins && term<countMatrix.length) : 
							"\n"+new String(line)+"\nterm="+term+", bin="+bin+"/"+bins+", d="+d+
							", max="+maxArray[term]+"\n"+Arrays.toString(maxArray);
						countMatrix[term][bin]++;
					}
				}
			}
			line=bf.nextLine();
		}
	}
	
	private void writeData() {
		ByteBuilder bb=new ByteBuilder();
		for(int term=0; term<terms; term++) {
			String name=header[term];
			long[] counts=countMatrix[term];
			ByteStreamWriter bsw=ByteStreamWriter.makeBSW(name+".tsv", overwrite, false, false);
			bb.clear().append(name+"\tcount\n");
			linesOut++;
			double incr=maxArray[term]/bins;
			double binStart=0;
			for(int bin=0; bin<counts.length; bin++) {
				bb.append(binStart, 6, true).tab().append(counts[bin]).nl();
				linesOut++;
				binStart+=incr;
			}
			bytesOut+=bb.length();
			bsw.print(bb);
			bsw.poison();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	
	String[] header;
	long[][] countMatrix;
	double[] maxArray;
	int bins=1000;
	int terms=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
