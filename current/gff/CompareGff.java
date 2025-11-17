package gff;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import prok.ProkObject;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.StringNum;

/**
 * Compares gff files for the purpose of grading gene-calling.
 * @author Brian Bushnell
 * @date October 3, 2018
 *
 */
public class CompareGff {
	
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
		CompareGff x=new CompareGff(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CompareGff(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
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
			
			in=parser.in1;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffin=FileFormat.testInput(in, FileFormat.GFF, null, true, true);
		ffref=FileFormat.testInput(ref, FileFormat.GFF, null, true, true);
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

			if(a.equals("ref")){
				ref=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
//				ByteFile1.verbose=verbose;
//				ByteFile2.verbose=verbose;
//				ReadWrite.verbose=verbose;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(i==0 && arg.indexOf('=')<0){
				parser.in1=arg;
			}else if(i==1 && arg.indexOf('=')<0 && ref==null){
				ref=arg;
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
		ref=Tools.fixExtension(ref);
		if(in==null || ref==null){throw new RuntimeException("Error - at least two input files are required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(true, true, in, ref)){
			throw new RuntimeException("\nCan't read some input files.\n");  
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
	 * Main processing method for GFF comparison analysis.
	 * Loads reference GFF, processes query lines, calculates statistics,
	 * and outputs comprehensive comparison metrics including SNR.
	 * @param t Timer for performance tracking
	 */
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin);
		
		processInner(bf);
		
		errorState|=bf.close();
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Ref count:           \t"+refCount);
		outstream.println("Query count:         \t"+queryCount);

		outstream.println();
		outstream.println("Ref-relative counts:");
		outstream.println("True Positive Start: \t"+truePositiveStart+"\t"+(Tools.format("%.3f%%", truePositiveStart*100.0/refCount)));
		outstream.println("True Positive Stop:  \t"+truePositiveStop+"\t"+(Tools.format("%.3f%%", truePositiveStop*100.0/refCount)));
//		outstream.println("False Positive Start:\t"+falsePositiveStart+"\t"+(Tools.format("%.3f%%", falsePositiveStart*100.0/refCount)));
//		outstream.println("False Positive Stop: \t"+falsePositiveStop+"\t"+(Tools.format("%.3f%%", falsePositiveStop*100.0/refCount)));
		outstream.println("False Negative Start:\t"+falseNegativeStart+"\t"+(Tools.format("%.3f%%", falseNegativeStart*100.0/refCount)));
		outstream.println("False Negative Stop: \t"+falseNegativeStop+"\t"+(Tools.format("%.3f%%", falseNegativeStop*100.0/refCount)));

		outstream.println();
		outstream.println("Query-relative counts:");
		outstream.println("True Positive Start: \t"+truePositiveStart2+"\t"+(Tools.format("%.3f%%", truePositiveStart2*100.0/queryCount)));
		outstream.println("True Positive Stop:  \t"+truePositiveStop2+"\t"+(Tools.format("%.3f%%", truePositiveStop2*100.0/queryCount)));
		outstream.println("False Positive Start:\t"+falsePositiveStart2+"\t"+(Tools.format("%.3f%%", falsePositiveStart2*100.0/queryCount)));
		outstream.println("False Positive Stop: \t"+falsePositiveStop2+"\t"+(Tools.format("%.3f%%", falsePositiveStop2*100.0/queryCount)));
		
		outstream.println();
		outstream.println("SNR: \t"+Tools.format("%.4f", 10*Math.log10((truePositiveStart2+truePositiveStop2+0.1)/(falsePositiveStart2+falsePositiveStop2+0.1))));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Core comparison logic between reference and query GFF files.
	 * Builds hash maps from reference file for efficient lookup, then processes
	 * each query line to determine true/false positives and negatives.
	 * @param bf ByteFile reader for the query GFF file
	 */
	@SuppressWarnings("unchecked")
	private void processInner(ByteFile bf){
		byte[] line=bf.nextLine();
		
		{
			ArrayList<GffLine> refLines=GffLine.loadGffFile(ffref, "CDS,rRNA,tRNA", true);

			refCount=refLines.size();
			lineMap=new HashMap<StringNum, GffLine>();
			startCountMap=new HashMap<StringNum, Integer>();
			stopCountMap=new HashMap<StringNum, Integer>();
			
			for(GffLine gline : refLines){
				final int stop=gline.trueStop();
				StringNum sn=new StringNum(gline.seqid, stop);
				lineMap.put(sn, gline);
				startCountMap.put(sn, 0);
				stopCountMap.put(sn, 0);
				assert(lineMap.get(sn)==gline);
//				assert(false) : "\n\nsn='"+sn+"'\n"+lineMap.containsKey(sn)+"\n"+lineMap.keySet();
			}
			if(verbose){
				System.err.println(lineMap);
				System.err.println(startCountMap);
				System.err.println(stopCountMap);
			}
		}

		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				final boolean valid=(line[0]!='#');
				if(valid){
					queryCount++;
					GffLine gline=new GffLine(line);
					processLine(gline);
				}
			}
			line=bf.nextLine();
		}
		
		for(Entry<StringNum, Integer> e : startCountMap.entrySet()){
			if(e.getValue()<1){
				falseNegativeStart++;
			}
		}
		for(Entry<StringNum, Integer> e : stopCountMap.entrySet()){
			if(e.getValue()<1){
				falseNegativeStop++;
			}
		}
	}
	
	/**
	 * Processes individual query GFF line against reference data.
	 * Compares feature type, strand, start, and stop positions to classify
	 * as true positive, false positive, or false negative matches.
	 * @param gline Query GFF line to evaluate against reference
	 */
	private void processLine(GffLine gline){
//		boolean cds=gline.type.equals("CDS");
//		boolean trna=gline.type.equals("tRNA");
//		boolean rrna=gline.type.equals("rRNA");
//		if(!cds && !trna && !rrna){return;}
//		if(cds && !ProkObject.callCDS){return;}
//		if(trna && !ProkObject.calltRNA){return;}
//		if(rrna){
//			int type=gline.prokType();
//			if(ProkObject.processType(type)){return;}
//		}
		int type=gline.prokType();
		if(!ProkObject.processType(type)){return;}
		
		final int stop=gline.trueStop();
		final int start=gline.trueStart();
		
//		System.err.println("Considering "+start+", "+stop);

		StringNum sn=new StringNum(gline.seqid, stop);
		GffLine refline=lineMap.get(sn);
		
		boolean fail=(refline==null || refline.strand!=gline.strand || !refline.type.equals(gline.type));
		if(fail){
			if(verbose){
				System.err.println("Can't find "+sn+"\n"+gline+"\n"+refline);
				assert(false) : "\n\nsn='"+sn+"'\n"+lineMap.containsKey(sn)+"\n"+lineMap.keySet();
			}
			falsePositiveStart++;
			falsePositiveStop++;
			falsePositiveStart2++;
			falsePositiveStop2++;
		}else{
			assert(stop==refline.trueStop());
			truePositiveStop++;
			truePositiveStop2++;
			stopCountMap.put(sn, stopCountMap.get(sn)+1);
			if(start==refline.trueStart()){
				truePositiveStart++;
				truePositiveStart2++;
				startCountMap.put(sn, startCountMap.get(sn)+1);
			}else{
				falsePositiveStart++;
				falsePositiveStart2++;
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input query GFF file path */
	private String in=null;
	/** Reference GFF file path for comparison */
	private String ref=null;
	
	
	/*--------------------------------------------------------------*/

	/** Maps sequence ID and stop position to reference GFF lines for lookup */
	private HashMap<StringNum, GffLine> lineMap;
	/** Tracks count of start position matches for each reference feature */
	private HashMap<StringNum, Integer> startCountMap;
	/** Tracks count of stop position matches for each reference feature */
	private HashMap<StringNum, Integer> stopCountMap;
	
//	private HashMap<Integer, ArrayList<GffLine>> map;
//	private HashSet<Integer> stopSet;
//	private HashSet<Integer> startSet;
//	private HashSet<Integer> stopSetM;
//	private HashSet<Integer> startSetM;
	
	/** Total number of input lines processed */
	private long linesProcessed=0;
	/** Number of lines written to output */
	private long linesOut=0;
	/** Total bytes read from input files */
	private long bytesProcessed=0;
	/** Total bytes written to output files */
	private long bytesOut=0;
	
	/** Maximum number of lines to process (configurable limit) */
	private long maxLines=Long.MAX_VALUE;

	/** Count of incorrect start positions in query (reference-relative) */
	private long falsePositiveStart=0;
	/** Count of incorrect stop positions in query (reference-relative) */
	private long falsePositiveStop=0;
	/** Count of correct start positions in query (reference-relative) */
	private long truePositiveStart=0;
	/** Count of correct stop positions in query (reference-relative) */
	private long truePositiveStop=0;
	/** Count of missing start positions in query (reference features not found) */
	private long falseNegativeStart=0;
	/** Count of missing stop positions in query (reference features not found) */
	private long falseNegativeStop=0;
	
	/** Count of incorrect start positions in query (query-relative) */
	private long falsePositiveStart2=0;
	/** Count of incorrect stop positions in query (query-relative) */
	private long falsePositiveStop2=0;
	/** Count of correct start positions in query (query-relative) */
	private long truePositiveStart2=0;
	/** Count of correct stop positions in query (query-relative) */
	private long truePositiveStop2=0;
	
	/** Total number of features in reference GFF file */
	private long refCount=0;
	/** Total number of features in query GFF file */
	private long queryCount=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** File format handler for input query GFF file */
	private final FileFormat ffin;
	/** File format handler for reference GFF file */
	private final FileFormat ffref;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for results and error messages */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging */
	public static boolean verbose=false;
	/** Tracks whether processing encountered errors */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;
	
}
