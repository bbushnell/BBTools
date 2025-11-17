package driver;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import jgi.DecontaminateByNormalization;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date June 1, 2016
 *
 */
public class SummarizeCrossblock {
	
	/** Program entry point that creates SummarizeCrossblock instance and processes files.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		SummarizeCrossblock x=new SummarizeCrossblock(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs SummarizeCrossblock instance with argument parsing and file setup.
	 * Disables compression utilities, parses command-line arguments, and configures
	 * input/output files with validation.
	 * @param args Command-line arguments for configuration
	 */
	public SummarizeCrossblock(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_UNBGZIP=false;
		
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
		
		String in=null;
		{//Process parser fields
			Parser.processQuality();
			
			overwrite=parser.overwrite;
			append=parser.append;
			
			in=parser.in1;

			out1=parser.out1;
		}
		
		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(in.contains(",")){
			for(String s : in.split(",")){
				inList.add(s);
			}
		}else{
			inList.add(in);
			DecontaminateByNormalization.parseStringsFromFiles(inList);
		}
		
		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, false);
	}
	
	/**
	 * Main processing method that summarizes crossblock results for all input files.
	 * Creates output stream, processes each input file with ParseCrossblockResults,
	 * and writes tab-delimited summary including filename, copies, contigs, and bases
	 * statistics. Handles processing errors gracefully.
	 *
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		final TextStreamWriter tsw;
		tsw=ffout1!=null ? new TextStreamWriter(ffout1) : null;
		if(tsw!=null){tsw.start();}
		if(tsw!=null){tsw.print("#fname\tcopies\tcontigs\tcontigsDiscarded\tbases\tbasesDiscarded\n");}
		
		int i=1;
		for(String fname : inList){
			ParseCrossblockResults pcr=null;
			try{
				pcr=new ParseCrossblockResults(new String[] {"in="+fname});
				Timer t2=new Timer();
				pcr.process(t2);
				if(tsw!=null){tsw.print(fname+"\t"+i+"\t"+pcr.contigs()+"\t"+pcr.contigsDiscarded()+"\t"+pcr.bases()+"\t"+pcr.basesDiscarded()+"\n");}
			}catch(Throwable e){
				System.err.println(e);
				if(tsw!=null){tsw.print(fname+"\tERROR\n");}
			}
			i++;
		}
		if(tsw!=null){errorState|=tsw.poisonAndWait();}
	}

	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** List of input file paths to process */
	private ArrayList<String> inList=new ArrayList<String>();
	/** Output file path for summary results */
	private String out1=null;
	
	/*--------------------------------------------------------------*/

	/** Total count of bases retained after processing */
	private long basesKept=0;
	/** Total count of bases discarded during processing */
	private long basesDiscarded=0;
	/** Total count of contigs retained after processing */
	private long contigsKept=0;
	/** Total count of contigs discarded during processing */
	private long contigsDiscarded=0;

	/** Returns the total number of bases kept after processing */
	public long basesKept(){return basesKept;}
	/** Returns the total number of bases discarded during processing */
	public long basesDiscarded(){return basesDiscarded;}
	/** Returns the total number of contigs kept after processing */
	public long contigsKept(){return contigsKept;}
	/** Returns the total number of contigs discarded during processing */
	public long contigsDiscarded(){return contigsDiscarded;}
	
	/*--------------------------------------------------------------*/
	
	/** FileFormat for output file configuration */
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for messages, defaults to System.err */
	private PrintStream outstream=System.err;
	/** Global verbosity flag for detailed output */
	public static boolean verbose=false;
	/** Flag indicating whether an error occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
