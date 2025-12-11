package barcode;

import java.io.PrintStream;
import java.util.Collections;

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
 * Analyzes paired barcode counts from CountBarcodes2 output.
 * Generates statistics comparing expected vs actual barcode pair distributions
 * and identifies problematic barcode pairs that may indicate technical issues.
 *
 * @author Brian Bushnell
 * @date July 19, 2023
 */
public class AnalyzeBarcodes {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Program entry point.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		AnalyzeBarcodes x=new AnalyzeBarcodes(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command-line arguments and initializes the analyzer.
	 * Sets up input/output files, validates parameters, and configures processing options.
	 * @param args Command-line arguments containing file paths and options
	 */
	public AnalyzeBarcodes(String[] args){
		
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
	 * Parses command-line arguments into program configuration.
	 * Processes parameters for input files, expected barcode counts, output files,
	 * and various analysis options.
	 *
	 * @param args Array of command-line arguments in key=value format
	 * @return Configured Parser object with parsed settings
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

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("expected") || a.equals("expectedbarcodes")){
				if(Tools.isNumeric(b)) {
					expectedBarcodeCount=Integer.parseInt(b);
				}else {
					expectedBarcodeFile=b;
				}
			}else if(a.equals("outleft") || a.equals("left") || a.equals("leftstats")){
				outLeft=b;
			}else if(a.equals("outright") || a.equals("right") || a.equals("rightstats")){
				outRight=b;
			}else if(a.equals("outpair") || a.equals("pair") || a.equals("pairstats")){
				outPair=b;
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
	
	/** Adds or removes file extensions (.gz, .bz2) as needed for input files.
	 * Ensures input file path has correct compression extension. */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");} //Possible bug: misleading error message since expectedBarcodeFile is also required
	}
	
	/**
	 * Validates that input files exist and output files can be created.
	 * Checks for file accessibility and prevents duplicate file specifications.
	 * @throws RuntimeException if files cannot be read/written or duplicates exist
	 */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, outLeft, outRight, outPair)){
//			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to some output files.\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, expectedBarcodeFile)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, expectedBarcodeFile, outLeft, outRight, outPair)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjusts static file handling settings based on thread count.
	 * Enables multi-threaded file reading when sufficient threads are available. */
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
	 * Validates that required parameters are set and within acceptable ranges.
	 * Ensures input file and expected barcode file are specified.
	 * @return true if all parameters are valid
	 */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
		assert(in1!=null && expectedBarcodeFile!=null) : "Barcode stats and expected barcodes must be specified.";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main processing method that loads barcode statistics and generates analysis output.
	 * Loads barcode data, calculates statistics for left/right/pair distributions,
	 * and writes analysis results to specified output files.
	 * @param t Timer for tracking execution time and performance metrics
	 */
	void process(Timer t){
		
		bs=BarcodeStats.loadStatic(in1, expectedBarcodeFile, expectedBarcodeCount);
		linesProcessed+=bs.linesProcessed;
		bytesProcessed+=bs.bytesProcessed;
		bs.leftStats=bs.makeLeft();
		bs.rightStats=bs.makeRight();
		bs.calcStats();
		bs.leftStats.calcStats();
		bs.rightStats.calcStats();

//		float errorFraction1=bs.leftStats.errorFraction();
		float badPairFraction=bs.badPairFraction;
		
//		HashMap<String, Barcode> leftBadPairs
		
		if(outLeft!=null){
			ByteStreamWriter bsw=makeBSW(outLeft);
			printMismatchVsTotal(bs.leftStats, bsw);
			if(bsw!=null){errorState|=bsw.poisonAndWait();}
		}
		
		if(outRight!=null){
			ByteStreamWriter bsw=makeBSW(outRight);
			printMismatchVsTotal(bs.rightStats, bsw);
			if(bsw!=null){errorState|=bsw.poisonAndWait();}
		}
		
		if(outPair!=null){
			ByteStreamWriter bsw=makeBSW(outPair);
			printActualVsExpected(bs, bsw);
			if(bsw!=null){errorState|=bsw.poisonAndWait();}
		}
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Lines Out:         \t"+linesOut);
		outstream.println("Bytes Out:         \t"+bytesOut);
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Prints barcode statistics comparing total counts vs bad pair counts.
	 * Outputs tab-separated data showing each expected barcode's total occurrences
	 * and how many were part of problematic pairs.
	 *
	 * @param bs BarcodeStats containing the barcode count data
	 * @param bsw ByteStreamWriter for output file writing
	 */
	void printMismatchVsTotal(BarcodeStats bs, ByteStreamWriter bsw) {
		ByteBuilder bb=new ByteBuilder();
		bb.append("#Name\tTotal\tBadPair\n");
		linesOut++;
//		ArrayList<Barcode> list=new ArrayList<Barcode>(bs.expectedCodeList.size());
		for(Barcode e : bs.expectedCodeList) {
			Barcode bTotal=bs.codeMap.get(e.name);
			if(bTotal!=null) {e.setCount(bTotal.count());}
		}
		Collections.sort(bs.expectedCodeList);
		for(Barcode e : bs.expectedCodeList) {
			Barcode bTotal=bs.codeMap.get(e.name);
			Barcode bBad=bs.badPairMap.get(e.name);
			bb.append(e.name).tab().append((bTotal==null ? 0 : bTotal.count())).tab().append(bBad==null ? 0 : bBad.count()).nl();
			bsw.print(bb);
			bytesOut+=bb.length();
			bb.clear();
		}
	}
	
	/**
	 * Prints analysis of unexpected barcode pairs found in the data.
	 * Calculates expected frequencies for barcode pairs not in the expected set
	 * and compares them to actual observed counts to identify anomalies.
	 *
	 * @param bs BarcodeStats containing barcode pair distribution data
	 * @param bsw ByteStreamWriter for output file writing
	 */
	void printActualVsExpected(BarcodeStats bs, ByteStreamWriter bsw) {
		ByteBuilder bb=new ByteBuilder();
		bb.append("#Name\tExpectedCount\tActualCount\n");
		linesOut++;
		String dString=(bs.delimiter<1 ? "" : Character.toString((char)bs.delimiter));
		final double badRate=bs.badPairFraction;
		final double validLeft=bs.leftStats.expectedCodes;
		final double validRight=bs.rightStats.expectedCodes;
		
		//For two barcodes, a.count*b.count*badMult will give the expect number of bad pairs
		final double badMult=bs.validArray[0]/(validLeft*validRight);
		
		for(Barcode ebLeft : bs.leftStats.expectedCodeList) {
			final double leftCount;
			{
				Barcode bLeft=bs.leftStats.codeMap.get(ebLeft.name);
				leftCount=(bLeft==null ? 0 : bLeft.count());
			}
			for(Barcode ebRight : bs.rightStats.expectedCodeList) {
				final double rightCount;
				{
					Barcode bRight=bs.rightStats.codeMap.get(ebRight.name);
					rightCount=(bRight==null ? 0 : bRight.count());
				}
				String pair=ebLeft.name+dString+ebRight.name;
				if(!bs.expectedCodeMap.containsKey(pair)) {
					Barcode bPair=bs.codeMap.get(pair);
					final double pairCount=(bPair==null ? 0 : bPair.count());
					final double expectedCount=leftCount*rightCount*badMult;

					bb.append(pair).tab().append(expectedCount,2).tab().append(pairCount,0).nl();
					bsw.print(bb);
					bytesOut+=bb.length();
					linesOut++;
					bb.clear();
				}
			}
		}
	}
	
	private ByteStreamWriter makeBSW(String fname){
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TXT, null, true, overwrite, append, false);
		return makeBSW(ff);
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

	/** Primary input file path containing barcode statistics data */
	private String in1=null;
	
	private String expectedBarcodeFile=null;
	private int expectedBarcodeCount=-1;

	/** Output file path for left barcode statistics */
	private String outLeft="leftStats.txt";
	private String outRight="rightStats.txt";
	private String outPair="pairStats.txt";
	
	private BarcodeStats bs;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long bytesProcessed=0;
	private long linesOut=0;
	private long bytesOut=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** FileFormat object for the primary input file */
	private final FileFormat ffin1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	/** Enable verbose logging output */
	public static boolean verbose=false;
	/** Flag indicating whether an error occurred during processing */
	public boolean errorState=false;
	/** Allow overwriting existing output files */
	private boolean overwrite=true;
	/** Append output to existing files instead of overwriting */
	private boolean append=false;
	
}
