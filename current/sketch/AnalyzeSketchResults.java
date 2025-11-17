package sketch;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicInteger;

import fileIO.ByteFile;
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
import structures.FloatList;
import tax.TaxTree;
import template.ThreadWaiter;

/**
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class AnalyzeSketchResults {
	
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
		AnalyzeSketchResults x=new AnalyzeSketchResults(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public AnalyzeSketchResults(String[] args){
		
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
			in2=parser.in2;

			out1=parser.out1;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutMap=FileFormat.testOutput(outMap, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutAccuracy=FileFormat.testOutput(outAccuracy, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutBad=FileFormat.testOutput(outBad, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.TXT, null, true, true);
		
		recordSets=(ffoutAccuracy==null && ffoutBad==null ? null : new ArrayList<RecordSet>());
		tree=(treeFile==null) ? null : TaxTree.loadTaxTree(treeFile, outstream, true, false);
		SSUMap.load(outstream);
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

			if(a.equals("map") || a.equals("outmap")){
				outMap=b;
			}else if(a.equals("accuracy") || a.equals("outacc") || a.equals("outaccuracy")){
				outAccuracy=b;
			}else if(a.equals("outbad")){
				outBad=b;
			}else if(a.equals("tree")){
				treeFile=b;
			}else if(a.equals("shrinkonly")){
				shrinkOnly=Parse.parseBoolean(b);
			}else if(a.equals("ssu") || a.equals("ssufile")){
				assert(false) : "ssu and ssufile are deprecated; please specify 16S or 18S independently";
			}else if(a.equalsIgnoreCase("16S") || a.equalsIgnoreCase("16Sfile")){
				SSUMap.r16SFile=b;
			}else if(a.equalsIgnoreCase("18S") || a.equalsIgnoreCase("18Sfile")){
				SSUMap.r18SFile=b;
			}else if(a.equals("mash")){
				mode=MASH_MODE;
			}else if(a.equals("sourmash")){
				mode=SOURMASH_MODE;
			}else if(a.equals("blast")){
				mode=BLAST_MODE;
			}else if(a.equals("bbsketch")){
				mode=BBSKETCH_MODE;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("minsamples")){
				minSamples=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				RecordSet.verbose=Record.verbose=verbose=Parse.parseBoolean(b);
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
		in2=Tools.fixExtension(in2);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, outMap)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, outMap)){
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
	 * Main processing pipeline for sketch result analysis.
	 * Handles shrink-only mode, processes input files through ResultLineParser,
	 * generates accuracy statistics with SSU validation, creates ANI/AAI mappings,
	 * and outputs formatted results with timing information.
	 *
	 * @param t Timer for tracking execution performance
	 */
	void process(Timer t){
		if(verbose){System.err.println("process()");}

		ByteFile bf1=ByteFile.makeByteFile(ffin1);
		ByteFile bf2=(ffin2==null ? null : ByteFile.makeByteFile(ffin2));

		if(shrinkOnly){
			runShrinkOnly(bf1);
			return;
		}
		
		bswBad=makeBSW(ffoutBad);
		ResultLineParser parser=new ResultLineParser(mode, tree, bswBad, recordSets, false);
		
		processInner(bf1, parser, ffin2==null ? null : aniMap);
		ByteStreamWriter bsw=makeBSW(ffout1);
		printResults(parser, bsw);
		if(bsw!=null){errorState|=bsw.poisonAndWait();}

		ByteStreamWriter bswAcc=makeBSW(ffoutAccuracy);
		if(recordSets!=null){
			if(SSUMap.hasMap()){
				processSetsThreaded();
			}
			printAccuracy(bswAcc);
		}
		if(bswAcc!=null){errorState|=bswAcc.poisonAndWait();}
		if(bswBad!=null){errorState|=bswBad.poisonAndWait();}

		ByteStreamWriter bswMap=makeBSW(ffoutMap);
		if(bf2!=null){
			processInner(bf2, parser, aaiMap);
			printMap(bswMap);
		}
		if(bswMap!=null){errorState|=bswMap.poisonAndWait();}
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Lines Out:         \t"+linesOut);
		outstream.println("Bytes Out:         \t"+bytesOut);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Processes input in shrink-only mode for data reduction.
	 * Parses sketch results, sorts and filters records, then outputs
	 * only the essential data without full analysis or statistics.
	 * @param bf Input ByteFile containing sketch comparison results
	 */
	private void runShrinkOnly(ByteFile bf){
		if(verbose){System.err.println("runShrinkOnly("+bf.name()+")");}
		ResultLineParser parser=new ResultLineParser(mode, tree, null, null, true);
		ByteStreamWriter bsw=makeBSW(ffout1);

		byte[] line=bf.nextLine();

		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				parser.parse(line);
				RecordSet rs=parser.processData(null, true);
				if(rs!=null){
					rs.sortAndSweep();
					for(Record r : rs.records){
						if(bsw!=null){
							bsw.println(r.text);
						}
					}
				}
			}
			line=bf.nextLine();
		}
		errorState|=bf.close();
		
		RecordSet rs=parser.currentSet;
		if(rs!=null){
			rs.sortAndSweep();
			for(Record r : rs.records){
				if(bsw!=null){
					bsw.println(r.text);
				}
			}
		}
		
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
	}
	
	/**
	 * Core processing loop for reading and parsing sketch result lines.
	 * Reads input file line by line, parses each result through ResultLineParser,
	 * and populates the provided map with processed data.
	 *
	 * @param bf Input ByteFile to process
	 * @param parser ResultLineParser for processing sketch result lines
	 * @param map HashMap to store processed results (ANI or AAI mappings)
	 */
	private void processInner(ByteFile bf, ResultLineParser parser, HashMap<Long, Float> map){
		if(verbose){System.err.println("processInner("+bf.name()+")");}
		byte[] line=bf.nextLine();
		
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				parser.parse(line);
				parser.processData(map, recordSets!=null);
			}
			line=bf.nextLine();
		}
		errorState|=bf.close();
	}
	
	/**
	 * Multi-threaded processing of RecordSet collection for SSU analysis.
	 * Creates worker threads to process SSU rRNA data in parallel,
	 * with each thread handling record sorting and SSU validation.
	 */
	private void processSetsThreaded(){
		if(verbose){System.err.println("processSetsThreaded("+Shared.threads()+")");}
		final int threads=Shared.threads();
//		final ThreadWaiter waiter=new ThreadWaiter();
		final ArrayList<SSUThread> list=new ArrayList<SSUThread>(threads);
		final AtomicInteger atom=new AtomicInteger(0);
		for(int i=0; i<threads; i++){
			list.add(new SSUThread(atom));
		}
		boolean success=ThreadWaiter.startAndWait(list);
		if(!success){errorState=true;}
	}
	
	/**
	 * Generates accuracy statistics across taxonomic levels.
	 * Evaluates correctness of taxonomic assignments by comparing sketch results
	 * with SSU rRNA validation, reporting correct, incorrect, and missing classifications.
	 * @param bswAcc Output stream for accuracy results, or null for stdout
	 */
	private void printAccuracy(ByteStreamWriter bswAcc){
		if(verbose){System.err.println("printAccuracy("+bswAcc+")");}
		int[][] results=new int[taxLevels][16];
		for(RecordSet rs : recordSets){
			if(mode==MASH_MODE){
				rs.sortAndSweep();
				rs.processSSU();
			}
			int[] statusArray=rs.test(bswBad);
			for(int level=0; level<statusArray.length; level++){
				int status=statusArray[level];
				results[level][status]++;
				assert(status==NOHIT || status==CORRECT || 
						status==INCORRECT_TAX_CORRECT_SSU || status==INCORRECT_TAX_INCORRECT_SSU || 
						status==INCORRECT_TAX_MISSING_SSU) : status;
			}
		}
		ByteBuilder bb=new ByteBuilder();
		bb.append("#Level          \tCorrect\tbadTaxGoodSSU\tBadTaxNoSSU\tbadTaxBadSSU\tNoHit\n");
		if(bswAcc!=null){bswAcc.print(bb);}
		for(String levelName : printLevels){
			int level=TaxTree.stringToLevelExtended(levelName);
			bb.clear();
			bb.append(TaxTree.levelToStringExtended(level));
			while(bb.length<"species subgroup".length()){bb.append(' ');}
			bb.tab();
			bb.append(results[level][CORRECT]).tab();
			bb.append(results[level][INCORRECT_TAX_CORRECT_SSU]).tab();
			bb.append(results[level][INCORRECT_TAX_MISSING_SSU]).tab();
			bb.append(results[level][INCORRECT_TAX_INCORRECT_SSU]).tab();
			bb.append(results[level][NOHIT]).nl();
			if(bswAcc!=null){bswAcc.print(bb);}
		}
	}
	
	/**
	 * Outputs ANI/AAI mapping results in tabular format.
	 * Correlates query and reference IDs with their calculated ANI and AAI values
	 * for comparative genomics analysis.
	 * @param bsw Output stream for mapping results, or null for stdout
	 */
	private void printMap(ByteStreamWriter bsw){
		if(verbose){System.err.println("printMap("+bsw+")");}
		ByteBuilder bb=new ByteBuilder();
		bb.append("#qID\trID\tANI\tAAI\n");
		if(bsw!=null){bsw.print(bb);}
		for(Entry<Long, Float> e : aniMap.entrySet()){
			long key=e.getKey();
			int qID=((int)(key>>>32));
			int rID=(int)(key&Integer.MAX_VALUE);
			float ani=e.getValue();
			Float aai=aaiMap.get(key);
			if(aai!=null){
				bb.clear();
				bb.append(qID).tab().append(rID).tab().append(ani, 3).tab().append(aai, 3).nl();
				if(bsw!=null){bsw.print(bb);}
				linesOut++;
				bytesOut+=bb.length();
			}
		}
	}
	
	/**
	 * Outputs aggregated statistics across taxonomic levels.
	 * Reports average ANI and SSU values, standard deviations, and sample counts
	 * for each taxonomic rank that meets minimum sample requirements.
	 *
	 * @param parser ResultLineParser containing accumulated statistics
	 * @param bsw Output stream for results, or null for stdout
	 */
	private void printResults(ResultLineParser parser, ByteStreamWriter bsw){
		if(verbose){System.err.println("printResults("+bsw+")");}
		ByteBuilder bb=new ByteBuilder();
		bb.append("#Level    \tRank\tANI_AVG\tSSU_AVG\tANI_STD\tSSU_STD\tSamples");
		if(bsw!=null){bsw.println(bb);}
		for(int level=0; level<taxLevels; level++){
			long aniCount=parser.levelCounts[level];
			long ssuCount=parser.levelCountsSSU[level];
			if(aniCount>=minSamples && ((1L<<level)&printLevelsMask)!=0){
				bb.clear();
				String name=TaxTree.levelToStringExtended(level);
				double aniSum=parser.levelAniSums[level];
				double ssuSum=parser.levelSSUSums[level];
				FloatList aniList=parser.aniLists[level];
				FloatList ssuList=parser.ssuLists[level];
				aniList.sort();
				ssuList.sort();

				double aniAvg=aniSum/aniCount;
				double ssuAvg=ssuSum/ssuCount;
				double aniStd=aniList.stdev();
				double ssuStd=ssuList.stdev();

				bb.append(name);
				while(bb.length<9){bb.append(' ');}
				bb.tab().append(level);
				bb.tab().append(aniAvg, 3);
				bb.tab().append(ssuAvg, 3);
				bb.tab().append(aniStd, 3);
				bb.tab().append(ssuStd, 3);
				bb.tab().append(aniCount);
				if(bsw!=null){bsw.println(bb);}
				linesOut++;
				bytesOut+=bb.length();
			}
		}
	}
	
	/**
	 * Creates and starts a ByteStreamWriter for the specified FileFormat.
	 * Returns null if FileFormat is null, otherwise initializes and starts
	 * the writer for output operations.
	 *
	 * @param ff FileFormat specification for output, may be null
	 * @return Started ByteStreamWriter or null if no output specified
	 */
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(verbose){System.err.println("makeBSW("+ff+")");}
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/**
	 * Creates a bitmask for taxonomic levels to include in output.
	 * Converts taxonomic level names to bit positions and combines
	 * them into a single long value for efficient level filtering.
	 *
	 * @param printLevelsArray Array of taxonomic level names to include
	 * @return Bitmask representing selected taxonomic levels
	 */
	private static long makePrintLevelsMask(String[] printLevelsArray){
		long mask=0;
		for(String s : printLevelsArray){
			int level=TaxTree.stringToLevelExtended(s);
			long bit=1L<<level;
			assert(Long.bitCount(bit)==1);
			mask|=bit;
		}
		return mask;
	}
	
	/*--------------------------------------------------------------*/
	
	private class SSUThread extends Thread{
		
		/** Constructor for SSU processing worker thread.
		 * @param atom_ Atomic counter for thread-safe work distribution */
		SSUThread(AtomicInteger atom_){
			atom=atom_;
		}
		
		@Override
		public void run(){
			for(int idx=atom.getAndIncrement(); idx<recordSets.size(); idx=atom.getAndIncrement()){
				RecordSet rs=recordSets.get(idx);
				rs.sortAndSweep();
				rs.processSSU();
			}
		}
		
		/** Thread-safe counter for distributing work across SSU processing threads */
		private final AtomicInteger atom;
	}
	
	/*--------------------------------------------------------------*/
	
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path for sketch comparison results */
	private String in1=null;
	/** Secondary input file path for additional sketch results (optional) */
	private String in2=null;
	/** Primary output file path for aggregated results */
	private String out1="stdout.txt";
	/** Output file path for ANI/AAI mapping results */
	private String outMap=null;
	/** Output file path for accuracy statistics and validation results */
	private String outAccuracy=null;
	/** Output file path for records failing validation or accuracy checks */
	private String outBad=null;
	/** Path to taxonomic tree file for classification validation */
	private String treeFile=null;
	
	/** Output stream writer for records failing validation checks */
	private ByteStreamWriter bswBad=null;
	
	/*--------------------------------------------------------------*/
	
	/** Count of input lines processed from sketch result files */
	private long linesProcessed=0;
	/** Count of lines written to output files */
	private long linesOut=0;
	/** Total bytes read from input files during processing */
	private long bytesProcessed=0;
	/** Total bytes written to output files */
	private long bytesOut=0;
	
	/** Flag to enable shrink-only mode for data reduction without full analysis */
	private boolean shrinkOnly=false;
	
	/**
	 * Minimum number of samples required to report statistics for a taxonomic level
	 */
	int minSamples=1;
	/**
	 * Maximum number of input lines to process, used for limiting analysis scope
	 */
	private long maxLines=Long.MAX_VALUE;
	
	/** Number of taxonomic levels supported by the classification system */
	final static int taxLevels=TaxTree.numTaxLevelNamesExtended;

	/** Mapping of query-reference ID pairs to Average Nucleotide Identity values */
	final HashMap<Long, Float> aniMap=new HashMap<Long, Float>();
	/** Mapping of query-reference ID pairs to Average Amino acid Identity values */
	final HashMap<Long, Float> aaiMap=new HashMap<Long, Float>();
	
	/**
	 * Taxonomic tree structure for classification validation and hierarchy traversal
	 */
	final TaxTree tree;
	/** Collection of RecordSet objects for accuracy validation and SSU analysis */
	final ArrayList<RecordSet> recordSets;
	
//	static HashMap<Integer, byte[]> ssuMap;
	
	/** Array of taxonomic level names to include in output reports */
	static final String[] printLevels=new String[] {"strain", "species", "genus", "family", "order", "class", "phylum", "superkingdom", "life"};
	/** Bitmask for efficient filtering of taxonomic levels in output generation */
	static final long printLevelsMask=makePrintLevelsMask(printLevels);
	
	/*--------------------------------------------------------------*/

	/** Status code indicating no classification hit was found */
	static int NOHIT=0;
	/** Status code indicating correct taxonomic classification */
	static int CORRECT=1;
	/** Status flag for incorrect taxonomic assignment */
	static int INCORRECT_TAX=2;
	/** Status flag for incorrect SSU rRNA classification */
	static int INCORRECT_SSU=4;
	/** Status flag indicating missing SSU rRNA data for validation */
	static int MISSING_SSU=8;
	/** Composite status for incorrect taxonomy but correct SSU classification */
	private static int INCORRECT_TAX_CORRECT_SSU=INCORRECT_TAX;
	/** Composite status for incorrect taxonomy and incorrect SSU classification */
	private static int INCORRECT_TAX_INCORRECT_SSU=INCORRECT_TAX|INCORRECT_SSU;
	/** Composite status for incorrect taxonomy with missing SSU data */
	private static int INCORRECT_TAX_MISSING_SSU=INCORRECT_TAX|MISSING_SSU;
	
	/** Analysis mode constant for BBSketch output format */
	static final int BBSKETCH_MODE=0;
	/** Analysis mode constant for Mash output format */
	static final int MASH_MODE=1;
	/** Analysis mode constant for Sourmash output format */
	static final int SOURMASH_MODE=2;
	/** Analysis mode constant for BLAST output format */
	static final int BLAST_MODE=3;
	/** Current analysis mode determining input format parsing strategy */
	static int mode=BBSKETCH_MODE;
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** FileFormat specification for primary input file */
	private final FileFormat ffin1;
	/** FileFormat specification for secondary input file */
	private final FileFormat ffin2;
	/** FileFormat specification for primary output file */
	private final FileFormat ffout1;
	/** FileFormat specification for ANI/AAI mapping output file */
	private final FileFormat ffoutMap;
	/** FileFormat specification for accuracy statistics output file */
	private final FileFormat ffoutAccuracy;
	/** FileFormat specification for failed validation records output file */
	private final FileFormat ffoutBad;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and logging information */
	private PrintStream outstream=System.err;
	/** Flag enabling detailed diagnostic output during processing */
	public static boolean verbose=false;
	/** Flag indicating whether any errors occurred during processing */
	public boolean errorState=false;
	/** Flag allowing existing output files to be overwritten */
	private boolean overwrite=true;
	/** Flag enabling append mode for output files instead of overwriting */
	private boolean append=false;
	
}
