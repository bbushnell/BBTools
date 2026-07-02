package driver;

import java.io.PrintStream;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tax.TaxFilter;

/**
 * Filters assembly summary files based on taxonomic criteria.
 * Processes tab-delimited assembly summary files and retains only lines
 * matching configured taxonomic filters. Designed for filtering NCBI-style
 * assembly summaries by taxonomic identifiers.
 *
 * @author Brian Bushnell
 * @date May 1, 2016
 */
public class FilterAssemblySummary {
	
	public static void main(String[] args){
		Timer t=new Timer();
		FilterAssemblySummary x=new FilterAssemblySummary(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public FilterAssemblySummary(String[] args){
		
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
			}else if(TaxFilter.validArgument(a)){
				//do nothing
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
		
		//Make the actual filter
		filter=TaxFilter.makeFilter(args);
	}
	
	void process(Timer t){
		
		final TextFile tf;
		{
			tf=new TextFile(ffin1);
			if(verbose){outstream.println("Started tf");}
		}
		
		final TextStreamWriter tsw;
		{
			//TODO: Possible bug [driver/FilterAssemblySummary#001] out= omitted => ffout1==null (FileFormat.testOutput
			//returns null for a null name), and new TextStreamWriter(null) NPEs immediately (TextStreamWriter ctor
			//dereferences ff.fastq() on line 1). in1 gets a clean "input required" throw (L85) but out1 gets none, so
			//running with no out= yields a raw NPE instead of a clean error or stdout. The dead `if(tsw!=null)` guard at
			//L129 shows optional-output was INTENDED, but only half-wired: construction here and poisonAndWait below are
			//unguarded. INTENT UNCLEAR (require-out-with-clean-error vs make-out-optional) => escalate, not fixing in place.
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
			if(verbose){outstream.println("Started tsw");}
		}

		long linesProcessed=0;
		long linesRetained=0;
		long charsProcessed=0;
		
		{
			String line;
			while((maxReads<0 || linesProcessed<maxReads) && (line=tf.nextLine())!=null){
				linesProcessed++;
				charsProcessed+=line.length();
				String result=processLine(line);
				if(result!=null){
					linesRetained++;
					if(tsw!=null){tsw.println(result);} //tsw is never null here (built unconditionally L112) => this guard is DEAD/vestigial; it is the fingerprint of the intended-optional-output that #001 left half-wired.
				}
			}
		}
		
		errorState|=tsw.poisonAndWait(); //GOOD: writer drain return captured into errorState (contrast SummarizeCoverage#002 which dropped it). [If #001 makes out optional, this too needs an `if(tsw!=null)` guard.]
		errorState|=tf.close(); //GOOD: reader close return also captured.
		
		t.stop();
		
		String kpstring=(linesRetained<100000 ? ""+linesRetained : linesRetained<100000000 ? (linesRetained/1000)+"k" : (linesRetained/1000000)+"m");
		while(kpstring.length()<8){kpstring=" "+kpstring;}
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, charsProcessed, 8));
		outstream.println("Lines Retained:     "+kpstring);
		
		
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	private String processLine(String line){
//		System.out.println("Processing line "+line);
		if(line.startsWith("#")){return null;} //header lines (NCBI assembly_summary starts data-header rows with '#') skipped.
		String[] split=line.split("\t");
		//TODO: Possible bug [driver/FilterAssemblySummary#002] LOW/format-contract: guard is assert-ONLY. With -da (the
		//shell default disables nothing, but production runs often use -da), a short/blank line (<7 cols) makes split[6]
		//throw AIOOBE, and a non-numeric col-7 makes Integer.parseInt throw NumberFormatException — both uncaught => the
		//whole run dies on ONE malformed line. Crash-loud, no corruption; same shape as EstherFilter#001 / SummarizeSealStats#002.
		assert(split.length>6) : split.length+"\n"+"'"+line+"'";
		String id=split[6]; //split[6] = column 7 = species_taxid in NCBI assembly_summary.txt (col 6/index 5 is 'taxid'). Using species_taxid is plausibly intentional (species-level filtering); NOT flagged as a bug — comprehension note only.
		int number=Integer.parseInt(id);
//		System.out.println("Found number "+number+" from "+id+"; node: "+filter.tree().getNode(number));
		boolean b=filter.passesFilter((int)number); //(int) cast redundant: number is already int. Cosmetic.
//		System.out.println("passesFilter? "+b);
		return b ? line : null;
	}
	
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	private final TaxFilter filter;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
