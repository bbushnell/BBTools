package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.LinkedHashMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;

/**
 * Constructs a directory and file tree of sequences corresponding to a taxonomic tree.
 * Reads FASTA sequences with taxonomic headers and organizes them into directory
 * structures matching the taxonomic hierarchy, with one file per taxonomic node.
 *
 * @author Brian Bushnell
 * @date December 12, 2017
 */
public class ExplodeTree{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Program entry point for command-line execution.
	 * Creates ExplodeTree instance, processes arguments, and executes the tree explosion.
	 * @param args Command line arguments for input/output paths and options
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		ExplodeTree x=new ExplodeTree(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor that parses command-line arguments and initializes the ExplodeTree.
	 * Sets up input/output paths, loads taxonomic tree, and validates file accessibility.
	 * @param args Command-line arguments containing input files, output paths, and options
	 */
	public ExplodeTree(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		//Create a parser object
		Parser parser=new Parser();

		//Parse each argument
		boolean useTree=true;
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("tree") || a.equals("usetree") || a.equals("taxtree")){
				useTree=TaxTree.parseTreeFlag(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("out") || a.equals("path") || a.equals("outpath")){
				outPath=b;
			}else if(a.equals("prefix")){
				prefix=b;
			}else if(a.equals("results") || a.equals("result")){
				resultsFile=b;
			}else if(a.equals("makedirectories") || a.equals("mkdirs") || a.equals("mkdir")){
				makeDirectories=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		if(prefix==null){prefix="";}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;

			overwrite=parser.overwrite;

			in1=parser.in1;

			extin=parser.extin;
		}

		if(outPath==null || outPath.trim().length()==0){outPath="";}
		else{
			outPath=outPath.trim().replace('\\', '/').replaceAll("/+", "/");
			if(!outPath.endsWith("/")){outPath=outPath+"/";}
		}

		assert(FastaReadInputStream.settingsOK());

		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, false, false, resultsFile)){
			outstream.println(resultsFile);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+resultsFile+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, resultsFile)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);

		if(!useTree){
			throw new RuntimeException("ExplodeTree requires a taxonomy tree; tree=f is not supported.");
		}
		tree=TaxTree.loadTaxTree(outstream, true, false);
		if(tree==null){
			throw new RuntimeException("ExplodeTree requires a taxonomy tree; none could be loaded.");
		}
	}

	/** Constructor for finite tests that supply their own tree and input reader.
	 * @param ffin1_ Primary input descriptor
	 * @param tree_ Taxonomy tree used for header parsing and output placement
	 * @param outPath_ Root output directory path
	 * @param prefix_ Output filename prefix
	 * @param overwrite_ Whether output files may be overwritten */
	protected ExplodeTree(final FileFormat ffin1_, final TaxTree tree_, final String outPath_,
			final String prefix_, final boolean overwrite_){
		ffin1=ffin1_;
		tree=tree_;
		outPath=outPath_;
		prefix=prefix_;
		overwrite=overwrite_;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates one directory per tax node and optional name marker files.
	 * @param root Root output directory
	 * @param writeNames Whether to write one name marker per tax node */
	public void makeDirectoryTree(String root, boolean writeNames){
		for(TaxNode node : tree.nodes){
			if(node!=null){
				String dir=tree.toDir(node, root);
				File df=new File(dir);
				if(!df.isDirectory()){errorState|=!df.mkdirs() && !df.isDirectory();}
				if(writeNames){
					try{
						String fname=node.simpleName()+".name";
						//FIXED [tax/ExplodeTree#001]: use dir+fname for the exists()
						//check as well as the write target. The old CWD check missed
						//target .name files and could skip on coincidental CWD files;
						//the pre-fix redundant writes were idempotent because content is node.name.
						File nf=new File(dir+fname);
						if(!nf.exists()){
							ReadWrite.writeString(node.name, dir+fname);
						}
					}catch(Exception e){
						e.printStackTrace();
						errorState=true;
					}
				}
			}
		}
	}

	/** Creates optional taxon directories, distributes reads, and reports counts.
	 * @param t Timer used to report elapsed processing time */
	public void process(Timer t){

		Timer t2=new Timer();
		if(makeDirectories){
			makeDirectoryTree(outPath, true);
			t2.stop("Finished making directories. ");
			t2.start();
		}
		processInner();
		t2.stop();
		t2.stop("Finished writing data. ");

		//Do anything necessary after processing

		if(resultsFile!=null){
			TextStreamWriter tsw=new TextStreamWriter(resultsFile, overwrite, false, false);
			tsw.start();
			for(TaxNode tn : nodes.keySet()){
				Long data=nodes.get(tn);
				if(data==null){data=0L;}
				tsw.println(tn.id+"\t"+data+"\t"+tn.levelStringExtended(false)+"\t"+tn.name);
			}
			errorState|=tsw.poisonAndWait();
		}

		//Report timing and results
		{
			t.stop();

			//Calculate units per nanosecond
			double rpnano=readsProcessed/(double)(t.elapsed);
			double lpnano=linesProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);

			//Add "k" and "m" for large numbers
			String rpstring=Tools.padKMB(readsProcessed, 8);
			String lpstring=Tools.padKMB(linesProcessed, 8);
			String bpstring=Tools.padKMB(basesProcessed, 8);

			String li="Lines In:               \t"+linesProcessed+" lines";
			String lo="Lines Out:              \t"+linesOut+" lines";
			while(lo.length()<li.length()){lo=lo+" ";}

			String ri="Reads In:               \t"+readsProcessed+" reads";
			String ro="Reads Out:              \t"+readsOut+" reads";
			while(ro.length()<ri.length()){ro=ro+" ";}

			outstream.println(ri+"\t"+basesProcessed+" bases");
			outstream.println(ro+"\t"+basesOut+" bases");
			outstream.println(li);
			outstream.println(lo);
			outstream.println();

			outstream.println("Time:                         \t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+Tools.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Lines Processed:    "+lpstring+" \t"+Tools.format("%.2fk reads/sec", lpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
		}

		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Core sequence processing loop that reads FASTA input and distributes sequences.
	 * Parses taxonomic information from headers, creates output files for each
	 * taxonomic node, and writes sequences to appropriate directory locations.
	 */
	void processInner(){
		ByteFile bf=makeByteFile(ffin1);
		TaxNode currentNode=null;
		long currentSize=0;
		ByteStreamWriter bsw=null;
		Throwable failure=null;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				linesProcessed++;
				if(line.length>0){
					final boolean header=(line[0]=='>');
					if(header){
						if(maxReads>0 && readsProcessed>=maxReads){break;}
						readsProcessed++;
						if(currentNode!=null){nodes.put(currentNode, nodes.get(currentNode)+currentSize);}

						final TaxNode tn=tree.parseNodeFromHeader(new String(line, 1, line.length-1), false);

						if(tn==null || tn!=currentNode){
							if(bsw!=null){
								final ByteStreamWriter old=bsw;
								bsw=null;
								errorState=old.poisonAndWait()|errorState;
							}
						}
						if(tn!=null && tn!=currentNode){
							String dir=tree.toDir(tn, outPath);
							final boolean found=nodes.containsKey(tn);
							if(!found){nodes.put(tn, 0L);}
							FileFormat ff=FileFormat.testOutput(dir+prefix+tn.id+".fa.gz",
									FileFormat.FASTA, null, true, overwrite && !found, found, false);
							bsw=makeWriter(ff);
							bsw.start();
						}

						currentNode=tn;
						currentSize=0;
						if(bsw!=null){readsOut++;}
					}else{
						basesProcessed+=line.length;
						currentSize+=line.length;
					}
					if(bsw!=null){
						linesOut++;
						if(!header){basesOut+=line.length;}
						bsw.println(line);
					}
				}
			}
			if(bsw!=null){
				final ByteStreamWriter old=bsw;
				bsw=null;
				errorState=old.poisonAndWait()|errorState;
				if(currentNode!=null){nodes.put(currentNode, nodes.get(currentNode)+currentSize);}
			}
		}catch(RuntimeException|Error e){
			failure=e;
			throw e;
		}finally{
			errorState|=closeProcessStreams(bsw, bf, failure);
		}
	}

	/** Closes the active writer and reader after normal or exceptional processing.
	 * Writers that failed before worker startup are closed directly; writers with
	 * a started worker are poisoned and joined.
	 * @param bsw Writer to clean up, or null
	 * @param bf Reader to close
	 * @param primary Primary processing failure, or null
	 * @return True if a stream reported a close error
	 */
	private boolean closeProcessStreams(final ByteStreamWriter bsw, final ByteFile bf, final Throwable primary){
		boolean errorState=false;
		Throwable failure=primary;
		try{
			if(bsw!=null){errorState|=closeWriter(bsw);}
		}catch(RuntimeException|Error e){
			failure=suppressDistinct(failure, e);
		}
		try{
			errorState|=bf.close();
		}catch(RuntimeException|Error e){
			failure=suppressDistinct(failure, e);
		}
		if(primary==null && failure!=null){throwUnchecked(failure);}
		return errorState;
	}

	/** Closes an output writer whether or not the worker started.
	 * @param bsw Writer whose constructor already opened the output stream
	 * @return True if writer cleanup reported an error */
	private static boolean closeWriter(final ByteStreamWriter bsw){
		return bsw.getState()==Thread.State.NEW ? bsw.closeBeforeStart() :
				bsw.poisonAndWait();
	}

	/** Suppresses one throwable onto another unless they are the same object.
	 * @param target Throwable that should receive {@code suppressed}, or null
	 * @param suppressed Throwable to attach, or null
	 * @return {@code target}, or {@code suppressed} when target was null
	 */
	private static Throwable suppressDistinct(final Throwable target, final Throwable suppressed){
		if(target==null){return suppressed;}
		if(suppressed!=null && suppressed!=target){target.addSuppressed(suppressed);}
		return target;
	}

	/** Throws a cleanup failure without wrapping.
	 * @param failure Cleanup failure to throw
	 */
	private static void throwUnchecked(final Throwable failure){
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		throw (Error)failure;
	}

	/** Opens the input reader; finite tests can inject a close-failing reader.
	 * @param ff Input file descriptor
	 * @return Byte reader for {@code ff} */
	protected ByteFile makeByteFile(final FileFormat ff){
		return ByteFile.makeByteFile(ff);
	}

	/** Creates an output writer; finite tests can inject start failures.
	 * @param ff Output file descriptor
	 * @return Writer for {@code ff} */
	protected ByteStreamWriter makeWriter(final FileFormat ff){return new ByteStreamWriter(ff);}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path containing sequences to be distributed */
	private String in1=null;

	/** Root output directory path where taxonomic structure will be created */
	private String outPath=null;

	/** Prefix prepended to per-taxon output filenames */
	private String prefix;

	/** Override input file extension for format detection */
	private String extin=null;

	/** Optional summary output path for per-taxon base counts */
	public String resultsFile=null;

	/** Whether to create the taxonomy directory tree before writing reads */
	public boolean makeDirectories=true;

	/** Accumulated output bases by taxon node */
	public LinkedHashMap<TaxNode, Long> nodes=new LinkedHashMap<TaxNode, Long>();

	/*--------------------------------------------------------------*/

	/** Number of sequence reads processed from input */
	protected long readsProcessed=0;
	/** Number of input lines processed including headers and sequences */
	protected long linesProcessed=0;
	/** Number of sequence bases processed from input */
	protected long basesProcessed=0;

	/** Number of sequence reads submitted to output files */
	public long readsOut=0;
	/** Number of input lines submitted to output files */
	public long linesOut=0;
	/** Total number of sequence bases submitted to output files */
	public long basesOut=0;

	/** Maximum number of input reads to process; -1 means no limit */
	private long maxReads=-1;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** FileFormat object for the primary input file */
	private final FileFormat ffin1;

	/** Taxonomy tree used for header parsing and output placement */
	private final TaxTree tree;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Output stream for status messages and progress reports */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging and detailed progress information */
	public static boolean verbose=false;
	/** Flag indicating whether an error was encountered during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;

}
