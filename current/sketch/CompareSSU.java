package sketch;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

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
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.FloatList;
import tax.TaxNode;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Compares SSU (Small Subunit) ribosomal RNA sequences in all-to-all or fractional matrix format.
 * Performs sequence identity comparisons between SSU sequences using taxonomic
 * hierarchy information.
 * Results are grouped by taxonomic levels and identity statistics are computed.
 *
 * @author Brian Bushnell
 * @date December 2, 2019
 */
public class CompareSSU implements Accumulator<CompareSSU.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Program entry point that creates a CompareSSU instance and processes sequences.
	 * @param args Command line arguments */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		CompareSSU x=new CompareSSU(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses arguments and initializes data structures.
	 * Loads SSU sequences from the input file, creates taxonomic tree,
	 * and prepares comparison matrices.
	 * @param args Command line arguments containing input files and parameters
	 */
	public CompareSSU(String[] args){
		
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
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}

		validateParams();
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, ordered);
		
		tree=(treeFile==null) ? null : TaxTree.loadTaxTree(treeFile, outstream, true, false);
		
		SSUMap.r16SFile=in1;
		if(SSUMap.r16SFile!=null){
			SSUMap.load(outstream);
			HashMap<Integer, byte[]> ssuMap=SSUMap.r16SMap;
			ssuList=new ArrayList<Read>(ssuMap.size());
			for(Entry<Integer, byte[]> e : ssuMap.entrySet()){
				int id=e.getKey();
				byte[] value=e.getValue();
				if(value.length>=minlen && value.length<=maxlen){
					Read r=new Read(value, null, id);//Sets numeric ID to TaxID.
					if(maxns<0 || r.countNocalls()<=maxns){
						ssuList.add(r);
					}
				}
			}
			Collections.shuffle(ssuList);
		}
		for(int i=0; i<idLists.length; i++){
			idLists[i]=new FloatList();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses command-line arguments into program parameters.
	 * Handles parameters like verbose, tree file, length filters, and comparison modes.
	 * @param args Command line argument array
	 * @return Configured Parser object with standard settings
	 */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("tree")){
				treeFile=b;
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("ata") || a.equals("alltoall")){
				allToAll=Parse.parseBoolean(b);
			}else if(a.equals("store") || a.equals("storeresults")){
				storeResults=Parse.parseBoolean(b);
			}else if(a.equals("minlen") || a.equals("maxlength")){
				minlen=Parse.parseIntKMG(b);
			}else if(a.equals("maxlen") || a.equals("maxlength")){
				maxlen=Parse.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("maxns")){
				maxns=Parse.parseIntKMG(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/**
	 * Validates that input files can be read and output files can be written.
	 * Checks for duplicate file specifications and file accessibility.
	 * @throws RuntimeException if files cannot be accessed or are duplicated
	 */
	private void checkFileExistence(){
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
	}
	
	/** Adjusts static file-related settings for optimal performance.
	 * Configures ByteFile modes based on thread count and validates stream settings. */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/** Validates parameter ranges and required settings.
	 * @return true if all parameters are valid */
	private boolean validateParams(){
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main processing method that executes SSU sequence comparisons.
	 * Spawns worker threads, collects results, and outputs statistics by taxonomic level.
	 * @param t Timer for tracking execution performance
	 */
	void process(Timer t){
		
		ByteStreamWriter bsw=makeBSW(ffout1);
		if(bsw!=null){
			bsw.forcePrint("#Level\tIdentity\tQueryID\tRefID\n");
		}
		
		//Reset counters
		queriesProcessed=0;
		comparisons=0;
		
		//Process the reads in separate threads
		spawnThreads(bsw);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		{
			ByteBuilder bb=new ByteBuilder();
			bb.append("\nLevel       \tCount\tMean"+(storeResults ? "\tMedian\t90%ile\t10%ile\tSTDev" : "")+"\n");
			outstream.print(bb);
			final int minlen="superkingdom".length();
			for(int level=0; level<taxLevels; level++){
				if(counts[level]>0){
					bb.clear();
					bb.append(TaxTree.levelToStringExtended(level));
					while(bb.length()<minlen){bb.space();}
					bb.tab().append(counts[level]).tab();
					bb.append(sums[level]/counts[level]*100, 3);
					if(storeResults){
						FloatList list=idLists[level];
						list.sort();
						double stdev=list.stdev();
						double median=list.median();
//						double mode=list.mode();
						double percent90=list.percentile(0.9f);
						double percent10=list.percentile(0.1f);
						bb.tab().append(median*100, 3).tab().append(percent90*100, 3).tab().append(percent10*100, 3).tab().append(stdev*100, 3);
					}
					bb.nl();
					outstream.print(bb);
				}
			}
		}
		
		//Report timing and results
		t.stop();
		outstream.println();
		outstream.println(Tools.timeQueriesComparisonsProcessed(t, queriesProcessed, comparisons, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates and manages ProcessThread instances for parallel sequence comparison.
	 * @param bsw Output stream writer for results */
	private void spawnThreads(ByteStreamWriter bsw){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(bsw, i, threads));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	/**
	 * Accumulates results from a completed ProcessThread.
	 * Merges query counts, comparison counts, identity lists, and error states.
	 * @param pt Completed ProcessThread containing results to accumulate
	 */
	@Override
	public final void accumulate(ProcessThread pt){
		queriesProcessed+=pt.querysProcessedT;
		comparisons+=pt.comparisonsT;
		errorState|=(!pt.success);

		for(int i=0; i<taxLevels; i++){
			idLists[i].addAll(pt.idListsT[i]);
			counts[i]+=pt.countsT[i];
			sums[i]+=pt.sumsT[i];
		}
	}
	
	/** Reports whether all processing completed successfully.
	 * @return true if no error state was encountered, false otherwise */
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Worker thread that performs pairwise sequence comparisons.
	 * Processes a subset of query sequences against reference sequences,
	 * computing identity scores and grouping by taxonomic levels.
	 */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(ByteStreamWriter bsw_, final int tid_, final int threads_){
			bsw=bsw_;
			threadID=tid_;
			threads=threads_;
			listCopy=new ArrayList<Read>(ssuList.size());
			listCopy.addAll(ssuList);
			for(int i=0; i<idListsT.length; i++){
				idListsT[i]=new FloatList();
			}
		}
		
		//Called by start()
		/** Thread execution method that processes assigned sequence comparisons.
		 * Called automatically when thread starts. */
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Main processing loop that iterates through assigned sequences.
		 * Uses atomic counter to distribute work among threads. */
		void processInner(){
			final long limit=Tools.min(ssuList.size(), (maxReads>0 ? maxReads : Integer.MAX_VALUE));
			
			for(int num=next.getAndIncrement(); num<limit; num=next.getAndIncrement()){
				Read r=ssuList.get(num);
				processRead(r);
			}
		}
		
		void processRead(Read query){
			if(query.numericID<1){return;}//invalid TID
			final int qid=(int)query.numericID;
			if(tree.getNode(qid)==null) {return;}
			if(querysProcessedT%5==0) {
				Collections.shuffle(listCopy);
			}
			querysProcessedT++;
			
			ByteBuilder bb=new ByteBuilder();
			
			long seen=0;
			for(Read ref :listCopy){
				int rid=(int)ref.numericID;
				if(rid!=qid && rid>0 && tree.getNode(rid)!=null){
					int aid=tree.commonAncestor(qid, rid);
					if(aid>0){
						TaxNode tn=tree.getNode(aid);
						if(tn.isRanked()){
							int level=tn.levelExtended;
							long mask=1L<<level;
							if(allToAll || ((mask&printLevels)!=0 && (mask&seen)==0)){
								seen|=mask;
								float identity=compare(query, ref, level);
								bb.append(TaxTree.levelToStringExtended(level)).tab().append(identity, 6);
								bb.tab().append(qid).tab().append(rid).nl();
							}
						}
					}
				}
			}
			if(bsw!=null){
				bsw.addJob(bb);
			}
		}
		
		float compare(Read query, Read ref, int level){
			comparisonsT++;
			float identity=SketchObject.align(query.bases, ref.bases);
			if(storeResults){idListsT[level].add(identity);}
			countsT[level]++;
			sumsT[level]+=identity;
			return identity;
		}

		/** Number of queries processed by this thread */
		protected long querysProcessedT=0;
		/** Number of comparisons performed by this thread */
		protected long comparisonsT=0;
		
		/** Whether this thread completed successfully */
		boolean success=false;
		
		/** Output stream writer for this thread's results */
		private final ByteStreamWriter bsw;
		/** Unique identifier for this processing thread */
		final int threadID;
		
		final int threads;
		
		ArrayList<Read> listCopy;
		
		final FloatList[] idListsT=new FloatList[taxLevels];
		long[] countsT=new long[taxLevels];
		double[] sumsT=new double[taxLevels];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path for SSU sequences */
	private String in1=null;
	
	private String treeFile="auto";

	/** Primary output file path for comparison results */
	private String out1=null;
	
	public static ArrayList<Read> ssuList=null;

	final static int taxLevels=TaxTree.numTaxLevelNamesExtended;
	static final String[] printLevelsArray=new String[] {"strain", "species", "genus", "family", "order", "class", "phylum", "superkingdom", "life"};
	static final long printLevels=makePrintLevels(printLevelsArray);
	
	private final TaxTree tree;
	
	private static final long makePrintLevels(String[] names){
		long mask=0;
		for(String s : names){
			int level=TaxTree.stringToLevelExtended(s);
			mask|=(1L<<level);
		}
		return mask;
	}
	
	private FloatList[] idLists=new FloatList[taxLevels];
	private long[] counts=new long[taxLevels];
	private double[] sums=new double[taxLevels];
	
	private int minlen=0;
	private int maxlen=Integer.MAX_VALUE;
	private int maxns=-1;
	
	/*--------------------------------------------------------------*/

	/** Total number of query sequences processed */
	protected long queriesProcessed=0;
	/** Total number of pairwise sequence comparisons performed */
	protected long comparisons=0;

	/** Maximum number of reads to process, -1 for no limit */
	private long maxReads=-1;
	
	private AtomicInteger next=new AtomicInteger(0);

	private boolean allToAll=false;
	private boolean storeResults=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output file format specification for comparison results */
	private final FileFormat ffout1;
	
	/** Returns the read-write lock for thread synchronization.
	 * @return ReadWriteLock instance used by this object */
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	/** Whether to print verbose status messages */
	public static boolean verbose=false;
	/** Whether an error was encountered during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;
	/** Whether results should be output in input order */
	private boolean ordered=false;
	
}
