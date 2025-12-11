package tax;

import java.io.PrintStream;
import java.util.Arrays;

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
import stream.FastaReadInputStream;
import structures.IntLongHashMap;

/**
 * Calculates the sizes of sequences corresponding to TaxIDs.
 * Processes FASTA files containing taxonomy-labeled sequences and computes
 * size statistics including base counts, sequence counts, and node counts
 * with cumulative values propagated up the taxonomic hierarchy.
 *
 * @author Brian Bushnell
 * @date December 13, 2017
 */
public class TaxSize {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		Timer t=new Timer();
		TaxSize x=new TaxSize(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * Parses command-line arguments, initializes I/O settings, loads taxonomy tree,
	 * and validates input/output files.
	 * @param args Command line arguments
	 */
	public TaxSize(String[] args){
		
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
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			
			in1=parser.in1;
			out=parser.out1;
			
			extin=parser.extin;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, false, false, out)){
			outstream.println( out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+ out+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1,  out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
		
		tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Create read streams and process all data.
	 * Executes the main processing pipeline, generates output statistics,
	 * and reports timing information.
	 * @param t Timer for tracking execution time
	 */
	public void process(Timer t){
		
		processInner();
		
		//Do anything necessary after processing
		
		if(out!=null){
			percolateUp();
			ByteStreamWriter bsw=new ByteStreamWriter(out, overwrite, false, false);
			bsw.start();
			int[] keys=(printEmptyNodes ? cNodeMap.toArray() : sizeMap.toArray());
			Arrays.sort(keys);
			bsw.print("#taxID\tbases\tbasesC\tseqs\tseqsC\tnodesC\n".getBytes());
			for(int key : keys){
				final long size=Tools.max(0, sizeMap.get(key));
				if(size>0 || printEmptyNodes){
					final long csize=Tools.max(0, cSizeMap.get(key));
					final long seqs=Tools.max(0, seqMap.get(key));
					final long cseqs=Tools.max(0, cSeqMap.get(key));
					final long cnodes=Tools.max(0, cNodeMap.get(key));
					bsw.print(key).print('\t');
					bsw.print(size).print('\t');
					bsw.print(csize).print('\t');
					bsw.print(seqs).print('\t');
					bsw.print(cseqs).print('\t');
					bsw.print(cnodes).print('\n');
				}
			}
			errorState|=bsw.poisonAndWait();
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
			String lo="Lines Out:              \t"+linesAssigned+" lines";
			while(lo.length()<li.length()){lo=lo+" ";}

			String ri="Reads In:               \t"+readsProcessed+" reads";
			String ro="Reads Out:              \t"+readsAssigned+" reads";
			while(ro.length()<ri.length()){ro=ro+" ";}

			outstream.println(ri+"\t"+basesProcessed+" bases");
			outstream.println(ro+"\t"+basesAssigned+" bases");
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
	
	void percolateUp(){
		for(TaxNode tn : tree.nodes){
			if(tn!=null){
				final int tid0=tn.id;
				final long size=sizeMap.get(tid0);
				final long seqs=seqMap.get(tid0);
				int tid=tid0;
				while(true){
					cNodeMap.increment(tid, 1);
					if(size>0){cSizeMap.increment(tid, size);}
					if(seqs>0){cSeqMap.increment(tid, seqs);}
					int pid=tree.getParentID(tid);
					if(pid==tid){break;}
					tid=pid;
				}
			}
		}
	}
	
	void processInner(){
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		TaxNode currentNode=null;
		long currentSize=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0){
				linesProcessed++;
				final boolean header=(line[0]=='>');
				if(header){
					if(maxReads>0 && readsProcessed>=maxReads){break;}
					readsProcessed++;
					
					if(currentNode!=null){
						sizeMap.increment(currentNode.id, currentSize);
						seqMap.increment(currentNode.id, 1);
					}
					
					currentNode=tree.parseNodeFromHeader(new String(line, 1, line.length-1), false);
					currentSize=0;
					if(currentNode!=null){readsAssigned++;}
				}else{
					basesProcessed+=line.length;
					currentSize+=line.length;
				}
				if(currentNode!=null){
					linesAssigned++;
					if(!header){basesAssigned+=line.length;}
				}
			}
		}
		if(currentNode!=null){
			sizeMap.increment(currentNode.id, currentSize);
			seqMap.increment(currentNode.id, 1);
		}
		bf.close();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;

	private String out=null;
	
	private String extin=null;
	
	public String taxTreeFile=null;

	public IntLongHashMap sizeMap=new IntLongHashMap();
	public IntLongHashMap cSizeMap=new IntLongHashMap();
	public IntLongHashMap seqMap=new IntLongHashMap();
	public IntLongHashMap cSeqMap=new IntLongHashMap();
	public IntLongHashMap cNodeMap=new IntLongHashMap();
	
	/*--------------------------------------------------------------*/

	protected long readsProcessed=0;
	protected long linesProcessed=0;
	protected long basesProcessed=0;

	/** Number of reads successfully assigned to taxonomy nodes */
	public long readsAssigned=0;
	public long linesAssigned=0;
	public long basesAssigned=0;

	private long maxReads=-1;
	
	private boolean printEmptyNodes=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private final FileFormat ffin1;
	
	private final TaxTree tree;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	
}
