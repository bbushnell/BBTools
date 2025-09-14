package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import aligner.SingleStateAlignerFlat2;
import consensus.BaseGraph;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.IntHashSet;
import structures.ListNum;
import tax.CanonicalLineage;
import tax.GiToTaxid;
import tax.PrintTaxonomy;
import tax.TaxNode;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Picks one ribosomal (16S) sequence per taxID.
 * 
 * @author Brian Bushnell
 * @date November 19, 2015
 *
 */
public class MergeRibo implements Accumulator<MergeRibo.ProcessThread> {
	
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
		MergeRibo x=new MergeRibo(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public MergeRibo(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
//		Shared.capBufferLen(40);//This does not help; the slowness comes from unevenness in list length during pickBest.
		//To fix it, long lists should be sorted to be first.
		
		BaseGraph.MAF_sub=0.251f;
		BaseGraph.MAF_del=0.0f;
		BaseGraph.MAF_ins=0.0f;
		BaseGraph.MAF_noref=0.0f;
		BaseGraph.trimDepthFraction=0.3f;
		BaseGraph.trimNs=true;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			extin=parser.extin;

			out1=parser.out1;
			extout=parser.extout;
		}

		validateParams();
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin=new ArrayList<FileFormat>(in.size());
		ffalt=FileFormat.testInput(alt, FileFormat.FASTA, extin, true, true);
		for(String s : in){
			FileFormat ff=FileFormat.testInput(s, FileFormat.FASTA, extin, true, true);
			ffin.add(ff);
		}
		
		//Determine how many threads may be used
		threads=Shared.threads();
		useTree|=(dada2 || taxLevelE>0);
		if(useTree) {loadTree();}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
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
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("consensus")){
				useConsensus=Parse.parseBoolean(b);
			}else if(a.equals("best")){
				useConsensus=!Parse.parseBoolean(b);
			}else if(a.equals("fast")){
				fast=Parse.parseBoolean(b);
			}else if(a.equals("minid")){
				minID=Float.parseFloat(b);
			}else if(a.equals("maxns")){
				maxns=Integer.parseInt(b);
			}else if(a.equals("minlen")){
				minlen=Integer.parseInt(b);
			}else if(a.equals("maxlen")){
				maxlen=Integer.parseInt(b);
			}else if(a.equals("in")){
				Tools.addFiles(b, in);
			}else if(a.equals("alt")){
				alt=b;
			}else if(a.equals("tree")){
				treePath=b;
				useTree=(treePath!=null);
			}else if(a.equals("usetree")){
				useTree=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("dada2")){
				dada2=Parse.parseBoolean(b);
			}else if(a.equals("level") || a.equals("taxlevel")){
				taxLevelE=TaxTree.parseLevelExtended(b);
				if(taxLevelE>0) {useTree=true;}
			}else if(a.equalsIgnoreCase("process16S") || a.equalsIgnoreCase("16S")){
				process16S=Parse.parseBoolean(b);
				process18S=!process16S;
			}else if(a.equalsIgnoreCase("process18S") || a.equalsIgnoreCase("18S")){
				process18S=Parse.parseBoolean(b);
				process16S=!process18S;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		assert(!in.isEmpty()) : "No input file.";
		return parser;
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
//		//Ensure that no file was specified multiple times
//		if(!Tools.testForDuplicateFiles(true, out1, in.toArray(new String[0]))){
//			throw new RuntimeException("\nSome file names were specified multiple times.\n");
//		}
	}
	
	/** Make sure interleaving agrees with number of input and output files */
	private void adjustInterleaving(){
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
//		assert(false) : "TODO";
		assert(process16S || process18S) : "16S or 18S must be selected.";
		assert(!process16S || !process18S) : "16S or 18S are both selected; only one may be active.";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){

		if(process16S){
			Read[] data=ProkObject.loadConsensusSequenceType("16S", true, true);
			consensus16S=data[0].bases;
			if(verbose){System.err.println("process16S: Loaded 16S consensus, length "+consensus16S.length+": "+new String(consensus16S));}
		}
		if(process18S){
			Read[] data=ProkObject.loadConsensusSequenceType("18S", true, true);
			consensus18S=data[0].bases;
			if(verbose){System.err.println("process18S: Loaded 18S consensus, length "+consensus18S.length+": "+new String(consensus18S));}
		}
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Align everything to global consensus
		for(FileFormat ff : ffin) {
			//Create a read input stream
			final ConcurrentReadInputStream cris=makeCris(ff);

			//Process the reads in separate threads
			spawnThreads(cris, false);
			errorState|=ReadWrite.closeStream(cris);
		}
		
		if(ffalt!=null){
			//Create a read input stream
			final ConcurrentReadInputStream cris=makeCris(ffalt);

			//Process the reads in separate threads
			spawnThreads(cris, true);
			errorState|=ReadWrite.closeStream(cris);
		}
		
//		queue=new ConcurrentLinkedQueue<ArrayList<Ribo>>();
//		for(Entry<Integer, ArrayList<Ribo>> e : listMap.entrySet()){
//			queue.add(e.getValue());
//		}
//		listMap=null;
		queue=makeQueue();
		
		//Run a second pass to pick the best SSU per taxID
		spawnThreads(null, false);
		
		//Do anything necessary after processing
		if(ffout1!=null){
			//Optionally create a read output stream
			final ConcurrentReadOutputStream ros=makeCros();
			long num=0;
			for(Ribo ribo : bestList){
				Read r=ribo.r;
				readsOut++;
				basesOut+=r.length();
				ArrayList<Read> list=new ArrayList<Read>(1);
				list.add(r);
				ros.add(list, num);
				num++;
			}
			//Close the read streams
			errorState|=ReadWrite.closeStream(ros);
		}
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private ConcurrentLinkedQueue<ArrayList<Ribo>> makeQueue(){
		ArrayList<ArrayList<Ribo>> listList=new ArrayList<ArrayList<Ribo>>(listMap.size());
		for(Entry<Integer, ArrayList<Ribo>> e : listMap.entrySet()){
			listList.add(e.getValue());
		}
		listMap=null;
		Collections.sort(listList, new ListComparator());
		assert(listList.isEmpty() || listList.get(0).size()>=listList.get(listList.size()-1).size());
		ConcurrentLinkedQueue<ArrayList<Ribo>> q=new ConcurrentLinkedQueue<ArrayList<Ribo>>();
		for(ArrayList<Ribo> x : listList){
			q.add(x);
		}
		return q;
	}
	
	private ConcurrentReadInputStream makeCris(FileFormat ff){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		assert(!paired) : "This should not be paired input.";
		return cris;
	}
	
	private ConcurrentReadOutputStream makeCros(){
		if(ffout1==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, boolean altData){
		
		//Do anything necessary prior to processing
		
		//Fill a list with ProcessThreads
		if(verbose){System.err.println("Spawning "+threads+" threads.");}
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, i, altData));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		if(verbose){System.err.println("Threads finished with success="+success+".");}
		errorState&=!success;
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		errorState|=(!pt.success);
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final int tid_, boolean alt_){
			cris=cris_;
			threadID=tid_;
			processInput=(cris!=null);
			altData=alt_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			if(processInput){
				//Process the reads
				processInner();
			}else{
				pickBest();
			}
//			assert(false) : processInput;
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			if(verbose && threadID==0){System.err.println("processInner() for tid="+threadID);}
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processInput(ln);
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
				
				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		void processInput(ListNum<Read> ln){
			if(verbose && threadID==0){System.err.println("processInput() for tid="+threadID);}

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r1.length();

				//Increment counters
				readsProcessedT++;
				basesProcessedT+=initialLength1;
				
				processRead(r1);
			}
		}
		
		void pickBest(){
			if(verbose && threadID==0){System.err.println("pickBest() for tid="+threadID);}
			for(ArrayList<Ribo> list=queue.poll(); list!=null; list=queue.poll()){
				Ribo best=pickBest(list);
				list.clear();
				if(best!=null) {//For example, fails length test or no taxonomy
					synchronized(bestList){
						bestList.add(best);
					}
				}
			}
		}
		
		/** Get the best sequence and possibly rename it */
		Ribo pickBest(ArrayList<Ribo> list) {
			//TODO: Customize all names, not just for dada2 output 
			Ribo best=pickBestInner(list);
			if(dada2) {
//				TaxNode tn=tree.getNode(best.tid);
////				System.err.println("tid="+best.tid+", node="+tn);
//				Object id=best.r.id;
//				if(tn==null) {
//					id=(best.tid+" No_Node");
//				}else if(tn.levelExtended==TaxTree.SPECIES_E){
////					bb.append(best.tid).space().append(tn.name);//Just the "id genus species".
//					id=PrintTaxonomy.makeTaxLine(tree, tn, TaxTree.SPECIES_E, TaxTree.KINGDOM_E, true, true);	
//				}else {//Some other level; use the full taxonomy from kingdom to species
//					id=PrintTaxonomy.makeTaxLine(tree, tn, TaxTree.SPECIES_E, TaxTree.KINGDOM_E, true, true);		
//				}
//				best.r.id=id.toString();
				CanonicalLineage cn=new CanonicalLineage(best.tid, tree);
				if(cn.levelsDefined<1) {return null;}
				best.r.id=cn.toString();
			}
			return best;
		}
		
		Ribo pickBestInner(ArrayList<Ribo> list){
			if(verbose && threadID==0){System.err.println("pickBest(list[="+list.size()+"]) for tid="+threadID);}
			assert(list!=null && list.size()>0);
			if(list.size()==1){return list.get(0);}
			Collections.sort(list);
			Collections.reverse(list);
			assert(list.get(0).product>=list.get(1).product);
			if(list.size()<3 || fast){return list.get(0);}
			
			Ribo base=list.get(0);
			int pad=Tools.max(10, (1600-base.r.length()));
			BaseGraph bg=new BaseGraph(base.r.name(), base.r.bases, base.r.quality, base.r.numericID, pad);
			for(Ribo r : list){
				bg.alignAndGenerateMatch(r.r, ssa);
			}
			Read consensus=bg.traverse();
			Ribo best;
			if(useConsensus){
				best=new Ribo(consensus, base.tid, 1);
			}else{
				for(Ribo r : list){
					float id=ssa.align(r.r.bases, consensus.bases);
					r.identity=id;
					r.product=score(r.length(), r.identity);
				}
				Collections.sort(list);
				Collections.reverse(list);
				assert(list.get(0).product>=list.get(1).product);
				best=list.get(0);
			}
			return best;
		}
		
		/**
		 * Process a read or a read pair.
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		void processRead(final Read r){
			if(verbose && threadID==0){System.err.println("processRead()");}
			if(r.length()<minlen || r.length()>maxlen){return;}
			if(maxns>=0 && r.countNocalls()>maxns){return;}
			Integer key=GiToTaxid.parseTaxidNumber(r.id, '|');
			if(verbose && threadID==0){System.err.println("key="+key);}
			
			if(key!=null && key>0 && taxLevelE>0) {
				int levelID=tree.getIdAtLevelExtended(key, taxLevelE);
				if(levelID>0 && levelID!=key) {key=levelID;}
			}
			
			if(key==null || key==-1 || (altData && seenTaxID.contains(key))){return;}
			float id=align(r);
			if(id<minID){return;}
			Ribo ribo=new Ribo(r, key, id);
			
			synchronized(listMap){
				ArrayList<Ribo> list=listMap.get(key);
				if(list==null){
					list=new ArrayList<Ribo>(8);
					listMap.put(key, list);
				}
				list.add(ribo);
				if(!altData){seenTaxID.add(key);}
			}
		}
		
		float align(Read r){
			float a=(process16S ? ssa.align(r.bases, consensus16S) : 0);
			float b=(process18S ? ssa.align(r.bases, consensus18S) : 0);
			if(verbose && threadID==0){System.err.println("Aligned; a="+a+", b="+b);}
			return Tools.max(a, b);
		}
		
		SingleStateAlignerFlat2 ssa=new SingleStateAlignerFlat2();

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Thread ID */
		final int threadID;
		
		//Run mode
		final boolean processInput;
		final boolean altData;
	}
	
	private class Ribo implements Comparable<Ribo>{
		Ribo(Read r_, int tid_, float identity_){
			r=r_;
			tid=tid_;
			identity=identity_;
			product=score(r.length(), identity);
		}
		
		@Override
		public int compareTo(Ribo b) {
			if(b.product>product){return -1;}
			else if(b.product<product){return 1;}
			else if(b.r.length()>r.length()){return -1;}
			else if(b.r.length()<r.length()){return 1;}
			return 0;
		}
		
		int length(){return r.length();}
		
		public String toString() {return "tid="+tid+", header="+r.id;}
		
		Read r;
		int tid;
		float identity;
		float product;
	}
	
	private class ListComparator implements Comparator<ArrayList<Ribo>> {

		@Override
		public int compare(ArrayList<Ribo> a, ArrayList<Ribo> b) {
			return a.size()>b.size() ? -1 : a.size()<b.size() ? 1 : 0;
		}
		
	}
	
	private float lengthMult(int len){
		int idealLength=idealLength();
		int max=Tools.max(len, idealLength, 1);
		int min=Tools.min(len, idealLength);
		return min/(float)max;
	}
	
	private float score(int len, float identity){
		return lengthMult(len)*identity;
	}
	
	private TaxTree loadTree() {
		if("auto".equals(treePath)){treePath=TaxTree.defaultTreeFile();}
		if(treePath!=null) {
			tree=TaxTree.loadTaxTree(treePath, System.err, false, false);
		}
		return tree;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file path */
	private ArrayList<String> in=new ArrayList<String>();
	
	/** Alternate input file path */
	private String alt=null;
	
	/** Primary output file path */
	private String out1=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;

	ArrayList<Ribo> bestList=new ArrayList<Ribo>();
	HashMap<Integer, ArrayList<Ribo>> listMap=new HashMap<Integer, ArrayList<Ribo>>(100000);
	ConcurrentLinkedQueue<ArrayList<Ribo>> queue;
	
	private TaxTree tree=null;
	private String treePath="auto";
	private boolean useTree=false;
	private int taxLevelE=-1;
	private boolean dada2=false;
	
	IntHashSet seenTaxID=new IntHashSet(1000000);
	
	byte[] consensus16S;
	byte[] consensus18S;
	
	int idealLength(){
		if(process16S){return consensus16S.length;}
		return consensus18S.length;
	}
	
	boolean useConsensus=false;
	boolean fast=false;
	int maxns=-1;
	int minlen=1;
	int maxlen=4000;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;

	private float minID=0.62f;
	
	private boolean process16S=true;
	private boolean process18S=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final ArrayList<FileFormat> ffin;
	private final FileFormat ffalt;
	
	/** Primary output file */
	private final FileFormat ffout1;
	
	final int threads;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
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
	/** Append to existing output files */
	private boolean append=false;
	/** Reads are output in input order */
	private boolean ordered=false;
	
}
