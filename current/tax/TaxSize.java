package tax;

import java.io.PrintStream;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.IntLongHashMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;

/**
 * Generates a taxsize TSV from FASTA records whose headers encode taxonomy IDs.
 * The output reports direct base/read counts per taxID, cumulative base/read
 * counts propagated through the loaded taxonomy tree, and cumulative node
 * counts including each emitted taxon itself and all descendants.
 *
 * @author Brian Bushnell
 * @date December 13, 2017
 */
public class TaxSize{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Runs the command-line tool.
	 * @param args command-line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		TaxSize x=new TaxSize(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Parses command-line arguments, configures FASTA input and optional TSV output,
	 * and loads the taxonomy tree used to resolve record headers.
	 * @param args command-line arguments
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
				// TODO: Throw RuntimeException so -da cannot ignore unknown arguments.
				assert(false) : "Unknown parameter "+args[i];
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
			outstream.println(out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out)){
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
	 * Scans the FASTA input, writes the optional taxsize TSV, and reports timing
	 * and assignment counts to the configured status stream.
	 * @param t timer for tracking execution time
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

	/**
	 * Fills cumulative maps by adding every loaded taxonomy node's direct sequence
	 * counts to that node and each ancestor up to the root.
	 */
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

	/**
	 * Reads the input FASTA one line at a time, counting nonheader line bytes as
	 * sequence length without validating residue contents, and assigns each
	 * sequence's bases to the taxID parsed from that sequence's primary header.
	 * Closes the acquired reader on normal and exceptional exits; a processing
	 * failure remains primary if cleanup also throws.
	 */
	void processInner(){
		ByteFile bf=makeByteFile(ffin1);
		Throwable thrown=null;
		try{
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
		}catch(RuntimeException|Error e){
			thrown=e;
			throw e;
		}finally{
			try{
				errorState|=bf.close();
			}catch(RuntimeException|Error closeFailure){
				if(thrown==null){throw closeFailure;}
				if(closeFailure!=thrown){thrown.addSuppressed(closeFailure);}
			}
		}
	}

	/** Opens the configured FASTA input for line-oriented reading.
	 * @param ff Validated FASTA input descriptor
	 * @return ByteFile reader for the input */
	protected ByteFile makeByteFile(final FileFormat ff){
		return ByteFile.makeByteFile(ff);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input FASTA file whose headers contain taxonomic identifiers. */
	private String in1=null;

	/** Optional taxsize TSV output path. */
	private String out=null;

	/** Optional input extension override. */
	private String extin=null;

	/** Taxonomy tree path, or auto for the default tree. */
	public String taxTreeFile=null;

	/** Direct bases assigned to each parsed taxID. */
	public IntLongHashMap sizeMap=new IntLongHashMap();
	/** Cumulative bases assigned to each taxID and all descendants. */
	public IntLongHashMap cSizeMap=new IntLongHashMap();
	/** Direct sequence count assigned to each parsed taxID. */
	public IntLongHashMap seqMap=new IntLongHashMap();
	/** Cumulative sequence count assigned to each taxID and all descendants. */
	public IntLongHashMap cSeqMap=new IntLongHashMap();
	/** Cumulative count of each taxID itself plus descendant nodes. */
	public IntLongHashMap cNodeMap=new IntLongHashMap();

	/*--------------------------------------------------------------*/

	/** FASTA records observed before the read limit. */
	protected long readsProcessed=0;
	/** Nonempty input lines observed, including the stopping header when {@link #maxReads} is reached. */
	protected long linesProcessed=0;
	/** Sequence bases observed before the read limit. */
	protected long basesProcessed=0;

	/** Number of reads successfully assigned to taxonomy nodes. */
	public long readsAssigned=0;
	/** Number of lines belonging to reads assigned to taxonomy nodes. */
	public long linesAssigned=0;
	/** Number of bases belonging to reads assigned to taxonomy nodes. */
	public long basesAssigned=0;

	/** Maximum number of FASTA records to count, or a nonpositive value for no limit. */
	private long maxReads=-1;

	/** Whether to emit taxa with no direct sequence bases. */
	private boolean printEmptyNodes=true;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Validated input descriptor. */
	private final FileFormat ffin1;

	/** Taxonomy tree used to resolve headers and propagate cumulative counts. */
	private final TaxTree tree;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Status stream from the preparser; closed by {@link #main(String[])}. */
	private PrintStream outstream=System.err;
	/** Verbose logging flag. */
	public static boolean verbose=false;
	/** True if input or output closure reported an error. */
	public boolean errorState=false;
	/** Whether an existing output path may be overwritten. */
	private boolean overwrite=true;

}
