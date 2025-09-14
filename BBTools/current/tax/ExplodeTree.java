package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.LinkedHashMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;

/**
 * Constructs a directory and file tree of sequences
 * corresponding to a taxonomic tree.
 * 
 * @author Brian Bushnell
 * @date December 12, 2017
 *
 */
public class ExplodeTree {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		ExplodeTree x=new ExplodeTree(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
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
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("out") || a.equals("path") || a.equals("outpath")){
				outPath=b;
			}else if(a.equals("prefix")){
				prefix=b;
			}else if(a.equals("results") || a.equals("result")){
				resultsFile=b;
			}else if(a.equals("makedirectories") || a.equals("mkdirs") || a.equals("mkdir")){
				makeDirectories=Parse.parseBoolean(b);
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		if(prefix==null){prefix="";}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		
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
		
		tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void makeDirectoryTree(String root, boolean writeNames){
		for(TaxNode node : tree.nodes){
			if(node!=null){
				String dir=tree.toDir(node, root);
				File df=new File(dir);
				if(!df.exists()){df.mkdirs();}
				if(writeNames){
					try {
						String fname=node.simpleName()+".name";
						File nf=new File(fname);
						if(!nf.exists()){
							ReadWrite.writeString(node.name, dir+fname);
						}
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
	}
	
	/** Create read streams and process all data */
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
	
	/** Iterate through the reads */
	void processInner(){
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		TaxNode currentNode=null;
		long currentSize=0;
		ByteStreamWriter bsw=null;
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
						if(bsw!=null){errorState=bsw.poisonAndWait()|errorState; bsw=null;}
					}
					if(tn!=null && tn!=currentNode){
						String dir=tree.toDir(tn, outPath);
						final boolean found=nodes.containsKey(tn);
						if(!found){nodes.put(tn, 0L);}
						FileFormat ff=FileFormat.testOutput(dir+prefix+tn.id+".fa.gz", FileFormat.FASTA, null, true, overwrite && !found, found, false);
						bsw=new ByteStreamWriter(ff);
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
			errorState=bsw.poisonAndWait()|errorState; bsw=null;
			if(currentNode!=null){nodes.put(currentNode, nodes.get(currentNode)+currentSize);}
		}
		bf.close();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String outPath=null;
	
	private String prefix;
	
	/** Override input file extension */
	private String extin=null;
	
	/** For listing what is present in the output */
	public String resultsFile=null;
	
	public String taxTreeFile=null;
	
	public boolean makeDirectories=true;
	
	public LinkedHashMap<TaxNode, Long> nodes=new LinkedHashMap<TaxNode, Long>();
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of lines processed */
	protected long linesProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads out */
	public long readsOut=0;
	/** Number of lines out */
	public long linesOut=0;
	/** Number of bases out */
	public long basesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	
	private final TaxTree tree;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	
}
