package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map.Entry;

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
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date January 14, 2017
 *
 */
public class CompareVCF {
	
	/** Program entry point for VCF comparison operations.
	 * @param args Command-line arguments specifying input files and operation mode */
	public static void main(String[] args){
		Timer t=new Timer();
		CompareVCF x=new CompareVCF(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs a CompareVCF instance from command-line arguments.
	 * Parses operation mode (difference/union/intersection), input files,
	 * quality filters, and variant splitting options.
	 * @param args Command-line arguments containing file paths and options
	 */
	public CompareVCF(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		int mode_=DIFFERENCE;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("difference") || a.equals("minus") || a.equals("dif") || a.equals("diff") || a.equals("subtraction") || a.equals("subtract")){
				mode_=DIFFERENCE;
			}else if(a.equals("union") || a.equals("plus")){
				mode_=UNION;
			}else if(a.equals("intersection") || a.equals("shared")){
				mode_=INTERSECTION;
			}else if(a.equals("addsamples")){
				addSamples=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("splitalleles")) {
				splitAlleles=Parse.parseBoolean(b);
			}else if(a.equals("splitsubs") || a.equals("splitsnps")) {
				splitSubs=Parse.parseBoolean(b);
			}else if(a.equals("splitcomplex")) {
				splitComplex=Parse.parseBoolean(b);
			}else if(a.equals("sass") || a.equals("split")) {
				splitAlleles=splitSubs=Parse.parseBoolean(b);
			}else if(a.equals("splitall") || a.equals("sascsss")) {
				splitAlleles=splitComplex=splitSubs=Parse.parseBoolean(b);
			}else if(a.equals("minscore") || a.equals("minqual") || a.equals("minq")) {
				minScore=Double.parseDouble(b);
			}else if(a.equals("trimtocanonical") || a.equals("canonicalize") || a.equals("canonicize") || a.equals("canonize")){
				VCFLine.TRIM_TO_CANONICAL=Parse.parseBoolean(b);
			}else if(a.equals("ref")){
				ref=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		mode=mode_;
		{//Process parser fields
			in1=(parser.in1==null ? null : parser.in1.split(","));
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null || in1.length<2){throw new RuntimeException("Error - at least two input files are required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}
		ffin1=new FileFormat[in1.length];
		for(int i=0; i<in1.length; i++){
			ffin1[i]=FileFormat.testInput(in1[i], FileFormat.VCF, null, true, true);
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.VCF, null, true, overwrite, append, false);
		if(ffout1!=null && ffout1.type()==FileFormat.VAR){outputVar=true;}
		if(ref!=null){ScafMap.loadReference(ref, null, null, true);}
	}
	
	/**
	 * Loads VCF variants from a file into a HashSet.
	 * Optionally splits complex variants and filters by minimum quality score.
	 * Updates processing statistics and header information.
	 *
	 * @param ff FileFormat for the input VCF file
	 * @param set Existing HashSet to add variants to, or null to create new set
	 * @return HashSet containing variants from the file
	 */
	public HashSet<VCFLine> getSet(FileFormat ff, HashSet<VCFLine> set){
		if(set==null){set=new HashSet<VCFLine>();}
		VCFFile vfile=new VCFFile(ff);
		samples.addAll(vfile.sampleNames);
		if(header==null){
			header=vfile.header;
			if(ScafMap.defaultScafMap()==null){
				ScafMap.setDefaultScafMap(vfile.toScafMap(null), ff.name());
			}
		}
		for(Entry<VCFLine, VCFLine> e : vfile.map.entrySet()){
			VCFLine v=e.getValue();
			ArrayList<VCFLine> list=null;
			if(splitAlleles || splitComplex || splitSubs){list=v.split(splitAlleles, splitComplex, splitSubs);}
			if(list==null || list.isEmpty()){
				if(!set.contains(v) && v.qual>=minScore){set.add(v);}
			}else{
				for(VCFLine line : list){
					if(!set.contains(line) && v.qual>=minScore){set.add(line);}
				}
			}
		}
		
		linesProcessed+=vfile.linesProcessed();
		headerLinesProcessed+=vfile.header.size();
		variantLinesProcessed+=vfile.map.size();
		bytesProcessed+=vfile.bytesProcessed();
		
		errorState|=vfile.errorState;
		return set;
	}
	
	/**
	 * Creates the union of all variants from input VCF files.
	 * Combines all unique variants from all input files into a single set.
	 * @return HashSet containing all unique variants from all input files
	 */
	public HashSet<VCFLine> union(){
		final HashSet<VCFLine> set=new HashSet<VCFLine>();
		for(FileFormat ff : ffin1){
			getSet(ff, set);
		}
		return set;
	}
	
	/**
	 * Creates the intersection of variants across all input VCF files.
	 * Returns only variants that are present in every input file.
	 * @return HashSet containing variants common to all input files
	 */
	public HashSet<VCFLine> intersection(){
		HashSet<VCFLine> set0=null;
		for(FileFormat ff : ffin1){
			HashSet<VCFLine> set=getSet(ff, null);
			if(set0==null){set0=set;}
			else{set0.retainAll(set);}
		}
		return set0;
	}
	
	/**
	 * Creates the difference between the first VCF file and all subsequent files.
	 * Returns variants present in the first file but absent from all others.
	 * @return HashSet containing variants unique to the first input file
	 */
	public HashSet<VCFLine> difference(){
		HashSet<VCFLine> set0=null;
		for(FileFormat ff : ffin1){
			HashSet<VCFLine> set=getSet(ff, null);
			if(set0==null){set0=set;}
			else{set0.removeAll(set);}
		}
		return set0;
	}
	
	/**
	 * Performs the specified set operation and returns results as a sorted list.
	 * Executes difference, union, or intersection based on configured mode.
	 * @return Sorted ArrayList of variants after applying the set operation
	 */
	ArrayList<VCFLine> toList(){
		final HashSet<VCFLine> set;
		if(mode==DIFFERENCE){
			set=difference();
		}else if(mode==UNION){
			set=union();
		}else if(mode==INTERSECTION){
			set=intersection();
		}else{
			throw new RuntimeException("Unknown mode "+mode);
		}
		ArrayList<VCFLine> list=new ArrayList<VCFLine>(set.size());
		list.addAll(set);
		Shared.sort(list);
		return list;
	}
	
	/**
	 * Main processing method that executes the VCF comparison and writes output.
	 * Performs the configured set operation and writes results to output file.
	 * Prints processing statistics and timing information.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		ArrayList<VCFLine> list=toList();
		
		ByteStreamWriter bsw=null;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}

		{//Processing block
			for(byte[] line : header){
				headerLinesOut++;
				if(bsw!=null){bsw.println(line);}
			}
			ByteBuilder bb=new ByteBuilder(33000);
			for(VCFLine line : list){
				variantLinesOut++;
				if(outputVar){
					assert(false) : "TODO";
//					Var v=line.toVar();
//					v.toText(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, ScafMap.defaultScafMap());
				}else{
					line.toText(bb);
				}
				bb.nl();
				if(bb.length>=32000){
					if(bsw!=null){bsw.print(bb);}
					bb.clear();
				}
			}
			if(bb.length>0){
				if(bsw!=null){bsw.print(bb);}
				bb.clear();
			}
		}
		
		if(bsw!=null){errorState|=bsw.poisonAndWait();}

		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Header Lines In:   \t"+headerLinesProcessed);
		outstream.println("Variant Lines In:  \t"+variantLinesProcessed);
		outstream.println("Header Lines Out:  \t"+headerLinesOut);
		outstream.println("Variant Lines Out: \t"+variantLinesOut);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/

	/** Total number of lines processed across all input files */
	private long linesProcessed=0;
	/** Number of header lines processed from input files */
	private long headerLinesProcessed=0;
	/** Number of variant lines processed from input files */
	private long variantLinesProcessed=0;
	/** Number of header lines written to output */
	private long headerLinesOut=0;
	/** Number of variant lines written to output */
	private long variantLinesOut=0;
	/** Total bytes processed from all input files */
	private long bytesProcessed=0;
	
	/** Maximum number of lines to process (configurable limit) */
	private long maxLines=Long.MAX_VALUE;

	/** VCF header lines from input files */
	public ArrayList<byte[]> header=null;
	/** Sample names collected from VCF headers */
	public ArrayList<String> samples=new ArrayList<String>();
	
	/*--------------------------------------------------------------*/
	
	/** Array of input VCF file paths */
	private String in1[]=null;
	/** Output file path */
	private String out1=null;
	/** Reference genome file path */
	private String ref=null;

	/** FileFormat objects for input VCF files */
	private final FileFormat ffin1[];
	/** FileFormat object for output file */
	private final FileFormat ffout1;
	
	/** Operation mode: DIFFERENCE, UNION, or INTERSECTION */
	public final int mode;
	
	/** Whether to add sample information to output */
	public boolean addSamples=true;
	/** Whether to output in VAR format instead of VCF */
	private boolean outputVar=false;

	/** Whether to split multi-allelic variants into separate lines */
	boolean splitAlleles=false;
	/** Whether to split substitution variants */
	boolean splitSubs=false;
	/** Whether to split complex variants into components */
	boolean splitComplex=false;
	/** Minimum quality score for variants to include in comparison */
	double minScore=-99999;
	
	/*--------------------------------------------------------------*/
	
	/** Constant for intersection operation mode (value 2) */
	/** Constant for union operation mode (value 1) */
	/** Constant for difference operation mode (value 0) */
	public static int DIFFERENCE=0, UNION=1, INTERSECTION=2;
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for logging and status messages */
	private PrintStream outstream=System.err;
	/** Flag for verbose output and debugging information */
	public static boolean verbose=false;
	/** Flag indicating whether an error occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;
	
}
