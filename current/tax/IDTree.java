package tax;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;

/**
 * Builds a greedy maximum-identity clustering tree from a tab-delimited matrix.
 * Each non-comment row starts with a label, followed by identities to preceding
 * rows. Values beyond that lower triangle are ignored. IDNode merges clusters
 * using maximum similarities and emits a Newick tree; no taxonomy database is used.
 *
 * @author Brian Bushnell
 * @date July 1, 2016
 */
public class IDTree{
	
	/** Program entry point.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		IDTree x=new IDTree(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor that parses command-line arguments and initializes file formats.
	 * Configures input/output files, validates file accessibility, and sets up
	 * processing parameters including overwrite and append modes.
	 * @param args Command-line arguments for configuration
	 */
	public IDTree(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
			
			maxLines=parser.maxReads;
		}
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Reads lower-triangular identity rows, merges clusters and optionally writes Newick.
	 * Lines beginning with # are skipped. maxLines limits accepted data rows, not
	 * physical lines. Input closure is attempted before clustering, including on
	 * parse failure. A close failure is thrown, or suppressed onto an earlier failure.
	 * @param t Started timer stopped after output
	 */
	void process(Timer t){
		
		ArrayList<IDNode> list=new ArrayList<IDNode>();
		TextFile tf=new TextFile(in1);
		Throwable failure=null;
		try{
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				if(maxLines>=0 && list.size()>=maxLines){break;}
				if(!line.startsWith("#")){
					String[] split=line.split("\t");
					double[] array=new double[list.size()];
					for(int i=0; i<array.length; i++){
						array[i]=Double.parseDouble(split[i+1]);
					}
					IDNode n=new IDNode(array, list.size(), split[0]);
					list.add(n);
				}
			}
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			try{
				if(tf.close()){throw new RuntimeException("Failed to close input "+in1);}
			}catch(RuntimeException | Error e){
				if(failure==null){throw e;}
				if(e!=failure){failure.addSuppressed(e);}
			}
		}
		
		IDNode[] nodes=list.toArray(new IDNode[0]);
		IDNode head=IDNode.makeTree(nodes);
		
		StringBuilder sb=head.toNewick();
		sb.append(';');
		
		if(out1!=null){
			ReadWrite.writeString(sb, out1);
			outstream.println("Wrote tree to "+out1);
		}
		
		t.stop();
		outstream.println("Time: \t"+t);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Unused template placeholder; always throws without processing either read.
	 *
	 * @param r1 Read 1
	 * @param r2 Read 2 (may be null)
	 * @return Never returns normally
	 * @throws RuntimeException Always, because this placeholder is unimplemented
	 */
	boolean processReadPair(final Read r1, final Read r2){
		throw new RuntimeException("TODO");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input matrix path; required by the constructor. */
	private String in1=null;

	/** Optional Newick output path; null still builds the tree without writing it. */
	private String out1=null;
	
	/*--------------------------------------------------------------*/

	/** Reserved template counter; the current processing loop does not update it. */
	protected long linesProcessed=0;
	
	/** Maximum accepted matrix rows, or a negative value for no limit. */
	protected long maxLines=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Validated input descriptor; processing currently opens in1 directly. */
	private final FileFormat ffin1;
	
	/** Validated output descriptor; writing currently uses out1 directly. */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Preparser-selected diagnostics stream, closed by main after normal completion. */
	private PrintStream outstream=System.err;
	/** Parsed verbosity flag; currently not consulted by processing. */
	public static boolean verbose=false;
	/** Reserved template failure flag; currently not updated or checked. */
	public boolean errorState=false;
	/** Whether an existing output is permitted by constructor validation. */
	private boolean overwrite=true;
	/** Append flag used for output validation; not forwarded to writeString. */
	private boolean append=false;
	
}
