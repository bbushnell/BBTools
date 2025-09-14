package bin;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import prok.GeneCaller;
import shared.LineParser1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;
import tracker.ReadStats;

/**
 * Loads fasta files with TID-labeled contigs,
 * to produce Clade record output with kmer frequencies, etc.
 * 
 * @author Brian Bushnell
 * @date April 12, 2025
 *
 */
public class CladeLoader extends BinObject implements Accumulator<CladeLoader.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		
		//Create an instance of this class
		CladeLoader x=new CladeLoader(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Default constructor.
	 */
	public CladeLoader(){
		if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150) {
			AdjustEntropy.load(4, 150);
		}
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CladeLoader(String[] args){
		
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
			
			extin=parser.extin;

			out=parser.out1;
			extout=parser.extout;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout=FileFormat.testOutput(out, FileFormat.CLADE, extout, true, overwrite, append, ordered);
		
		if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150) {
			AdjustEntropy.load(4, 150);
		}
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
			
			if(a.equals("tree")){
				if(b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true")) {
					treePath="auto";
				}else if(b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false")) {
					treePath=null;
				}else {
					treePath=b;
				}
			}else if(a.equals("in")){
				Tools.getFileOrFiles(b, in, true, false, false, false);
			}else if(a.equals("mergedupes")){
				mergeDuplicateTaxIDs=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("dummy") || a.equals("usedummy")){
				useDummy=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("a48")){
				Clade.outputCoding=Parse.parseBoolean(b) ? Clade.A48 : Clade.DECIMAL;
			}else if(a.equals("maxk") || a.equals("kmax")){
				Comparison.maxK=Clade.MAXK=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("callSSU")){
				Clade.callSSU=Parse.parseBoolean(b);
			}else if(a.equals("aligner") || a.equals("idaligner")){
				GeneCaller.useIDAligner=(b==null || !("f".equals(b) || "false".equals(b)));
				if(GeneCaller.useIDAligner) {aligner.Factory.setType(b);}
			}else if(a.equalsIgnoreCase("useTree")){
				useTree=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("16S")){
				r16sFile.clear();
				Tools.getFileOrFiles(b, r16sFile, true, false, false, false);
//				assert(false): r16sFile;
			}else if(a.equalsIgnoreCase("18S")){
				r18sFile.clear();
				Tools.getFileOrFiles(b, r18sFile, true, false, false, false);
//				assert(false): r18sFile;
			}else if(a.equalsIgnoreCase("replaceribo")){
				replaceRibo=Parse.parseBoolean(b);
			}
			
			else if(b==null && new File(arg).isFile()){
				in.add(arg);
			}else if(b==null && new File(arg).isDirectory()){
				Tools.getFileOrFiles(arg, in, true, false, false, false);
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
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in=Tools.fixExtension(in);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in.get(0), out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
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
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Reset counters
		readsProcessed=basesProcessed=0;
		
		ConcurrentHashMap<Integer, Clade> map=null;
		map=load(in, map);
		
		if(!r16sFile.isEmpty()) {
			for(String fname : r16sFile) {
				r16sAdded+=addRibo(map, fname, true);
			}
			System.err.println("Added "+r16sAdded+" 16S.");
		}
		if(!r18sFile.isEmpty()) {
			for(String fname : r18sFile) {
				r18sAdded+=addRibo(map, fname, true);
			}
			System.err.println("Added "+r18sAdded+" 18S.");
		}
		if(ffout!=null) {
			write(ffout, map.values());
		}
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.things("Clades", map.size(), 8));
		outstream.println(Tools.things("Bytes Out", bytesOut, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+
					" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Load clades from a collection of files.
	 * @param fnames Collection of file names to load
	 * @param map Map to store clades in (created if null)
	 * @return Map containing loaded clades
	 */
	public ConcurrentHashMap<Integer, Clade> load(Collection<String> fnames, 
			ConcurrentHashMap<Integer, Clade> map){
		assert(fnames!=null && !fnames.isEmpty());
		//Loop through input files
		for(String fname : fnames) {
			map=load(fname, map);
			assert(map!=null);
		}
		for(Clade c : map.values()) {
			Integer key=c.taxID;
			if(map.get(key)!=c) {
				System.err.println("key "+key+" mapped to: "+map.get(key)+"\ninstead of: "+c);
//				assert(false);//This happened ~3 times for RefSeq Bacteria output...
			}
			c.finish();
		}
		return map;
	}
	
	/**
	 * Load clades from a single file.
	 * @param fname File name to load
	 * @param map Map to store clades in (created if null)
	 * @return Map containing loaded clades
	 */
	public ConcurrentHashMap<Integer, Clade> load(String fname, ConcurrentHashMap<Integer, Clade> map){
		//Create input FileFormat object
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, extin, true, true);
		//Process the reads in separate threads
//		assert(ff.clade()) : ff.rawExtension()+", "+ff.type();
		if(ff.clade() || ff.extensionEquals("tsv")) {
			if(tree==null && treePath!=null && useTree) {tree=loadTree();}
			map=loadFromClade(ff, map);
		}else {
			if(tree==null && treePath!=null) {tree=loadTree();}
			map=loadFromSequence(ff, map);
		}
		return map;
	}

	/**
	 * Load clades from sequences in a file.
	 * @param ff FileFormat for input
	 * @param cladeMap Map to store clades in (created if null)
	 * @return Map containing loaded clades
	 */
	ConcurrentHashMap<Integer, Clade> loadFromSequence(FileFormat ff, 
			ConcurrentHashMap<Integer, Clade> cladeMap){
		outstream.println("Loading "+ff.name());
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris(ff);
		
		if(cladeMap==null) {cladeMap=new ConcurrentHashMap<Integer, Clade>(16000);}
		
		//Process the reads in separate threads
		spawnThreads(cris, cladeMap);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		return cladeMap;
	}

	/**
	 * Load clades from a clade file.
	 * @param ff FileFormat for input
	 * @param cladeMap Map to store clades in (created if null)
	 * @return Map containing loaded clades
	 */
	ConcurrentHashMap<Integer, Clade> loadFromClade(FileFormat ff, 
			ConcurrentHashMap<Integer, Clade> cladeMap){
		outstream.println("Loading "+ff.name());
		
		if(cladeMap==null) {cladeMap=new ConcurrentHashMap<Integer, Clade>(16000);}
		
		ByteFile bf=ByteFile.makeByteFile(ff);
		loadClades(bf, cladeMap);
		bf.close();
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		return cladeMap;
	}
	
	/**
	 * Load clades from a ByteFile.
	 * @param bf ByteFile containing clade data
	 * @param map Map to store clades in
	 */
	static void loadClades(ByteFile bf, ConcurrentHashMap<Integer, Clade> map) {
		LineParser1 lp=new LineParser1('\t');
		final ArrayList<byte[]> set=new ArrayList<byte[]>(20); //Should be 13 lines
		for(ListNum<byte[]> list=bf.nextList(); list!=null; list=bf.nextList()) {
			for(byte[] line : list) {
				if(Tools.startsWith(line, '#') && set.size()>5) {//New record
					addClade(set, map, lp);
					set.clear();
				}
				set.add(line);
			}
		}
		addClade(set, map, lp);//Possible last entry if not terminated with a '#'.
	}
	
	private static boolean addClade(final ArrayList<byte[]> set, 
			ConcurrentHashMap<Integer, Clade> map, LineParser1 lp) {
		if(set.size()<=5) {return false;}
		Clade c=Clade.parseClade(set, lp);
		Integer key=c.taxID;
		Clade old=map.get(key);
		if(old==null) {
			map.put(c.taxID, c);
			return true;
		}else if(c.bases>0){
//			System.err.println("Duplicate tid "+c.taxID);
			if(mergeDuplicateTaxIDs) {
				old.add(c);
			}else if(c.bases>old.bases) {
				map.put(c.taxID, c);
			}
		}
		return false;
	}
	
	/**
	 * Load clades from a clade-format file.
	 * @param ff FileFormat for input
	 * @return List of clades loaded from the file
	 */
	public static ArrayList<Clade> loadCladesFromClade(FileFormat ff) {
		ByteFile bf=ByteFile.makeByteFile(ff);
		LineParser1 lp=new LineParser1('\t');
		ArrayList<Clade> out=new ArrayList<Clade>();
		final ArrayList<byte[]> set=new ArrayList<byte[]>(20); //Should be 13 lines
		for(ListNum<byte[]> list=bf.nextList(); list!=null; list=bf.nextList()) {
			for(byte[] line : list) {
//				if(Tools.startsWith(line, '#')) {continue;}
//				assert(set.size()<=13);
				if(Tools.startsWith(line, '#') && !set.isEmpty()) {//New record
					Clade c=Clade.parseClade(set, lp);
					synchronized(c) {
						c.finish();
						out.add(c);
					}
					set.clear();
				}
				set.add(line);
			}
		}
		if(set.size()>1) {//New record
			Clade c=Clade.parseClade(set, lp);
			synchronized(c) {
				c.finish();
				out.add(c);
				set.clear();
			}
		}
		return out;
	}
	
	static int addRibo(ConcurrentHashMap<Integer, Clade> cladeMap, String ssuFile, boolean r16s) {
		FileFormat ff=FileFormat.testInput(ssuFile, FileFormat.FASTA, null, true, false);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null);
		cris.start(); //Start the stream
		int added=addRibo(cladeMap, cris, r16s);
		boolean errorState=false;
		errorState|=ReadWrite.closeStreams(cris);
		return added;
	}
	
	static int addRibo(ConcurrentHashMap<Integer, Clade> map, ConcurrentReadInputStream cris, boolean r16s) {
		ListNum<Read> ln=cris.nextList();
		int added=0;
		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){
			for(Read r : ln) {
				added+=addRibo(map, r, r16s);
			}
			cris.returnList(ln);
			ln=cris.nextList();
		}
		
		if(ln!=null){cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());}
		return added;
	}
	
	static int addRibo(ConcurrentHashMap<Integer, Clade> map, Read r, boolean r16s) {
		int tid=TaxTree.parseHeaderStatic(r.id);
		if(tid<0) {return 0;}
		Clade c=map.get(tid);
		if(c==null) {return 0;}
		if(r16s && (replaceRibo || c.r16S==null)) {c.r16S=r.bases;return 1;}
		if(!r16s && (replaceRibo || c.r18S==null)) {c.r18S=r.bases;return 1;}
		return 0;
	}
	
	/**
	 * Write a collection of clades to a file.
	 * @param ff FileFormat for output
	 * @param coll Collection of clades to write
	 */
	void write(FileFormat ff, Collection<Clade> coll) {
		if(ff==null) {return;}
		ArrayList<Clade> list=new ArrayList<Clade>(coll);
		Collections.sort(list);
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ff);
		ByteBuilder bb=new ByteBuilder(1024);
		for(Clade c : list) {
			bb.clear();
			c.toBytes(bb);
			bsw.print(bb);
			bytesOut+=bb.length;
			linesOut++;//Actually, clades out
		}
		bsw.poison();
	}
	
	/*--------------------------------------------------------------*/

	/**
	 * Create a ConcurrentReadInputStream for a FileFormat.
	 * @param ff FileFormat to create stream from
	 * @return ConcurrentReadInputStream for the FileFormat
	 */
	private ConcurrentReadInputStream makeCris(FileFormat ff){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(ff.fastq()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		else {assert(!paired);}
		return cris;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, ConcurrentHashMap<Integer, Clade> cladeMap){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, cladeMap, i));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt) {
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			errorState|=(!pt.success);
		}
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
	static class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, 
				ConcurrentHashMap<Integer, Clade> cladeMap_, final int tid_){
			cris=cris_;
			cladeMap=cladeMap_;
			tid=tid_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			synchronized(dummy) {
				processInner();
				if(useDummy) {emitDummy();}
			}
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			
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
				
				processList(ln);
				
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
		
		/**
		 * Process a list of reads.
		 * @param ln ListNum containing reads to process
		 */
		void processList(ListNum<Read> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				final Read r2=r1.mate;
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}
				if(r2!=null && !r2.validated()){r2.validate(true);}

				//Increment counters
				readsProcessedT+=r1.pairCount();
				basesProcessedT+=r1.pairLength();
				
				processContig(r1);
				if(r2!=null) {processContig(r2);}
			}
		}
		
		/**
		 * Process a single contig.
		 * @param r Read containing contig data
		 * @return True if successfully processed
		 */
		boolean processContig(final Read r){
			if(useDummy) {return processContigDummy(r);}
			int tid=resolveTaxID(r.id);
			if(tid<1) {return false;}
			Clade c=getOrMakeClade(tid);
			assert(c.taxID==tid);
			c.add(r.bases, et);
			return true;
		}
		
		/**
		 * Process a single contig using the dummy approach.
		 * @param r Read containing contig data
		 * @return True if successfully processed
		 */
		boolean processContigDummy(final Read r){
			int tid=resolveTaxID(r.id);
			if(tid<1) {return false;}
			
			if(tid!=dummy.taxID) {
				emitDummy();
				assert(dummy.taxID<0 && dummy.bases==0) : tid+", "+dummy.bases+", "+dummy;
				dummy.taxID=tid;
			}
			assert(dummy.taxID==tid) : tid+", "+dummy.bases+", "+dummy;
			dummy.add(r.bases, et);
			return true;
		}
		
		/**
		 * Emit the current dummy clade to the main clade map.
		 */
		private void emitDummy() {
			final int tid=dummy.taxID;
			if(tid<1) {return;}
			Clade c=getOrMakeClade(tid);
			assert(c.taxID==tid);
			c.add(dummy);
			dummy.clear();
		}
		
		/**
		 * Get a clade from the map or create a new one if it doesn't exist.
		 * @param tid Tax ID for the clade
		 * @return Clade object
		 */
		private Clade getOrMakeClade(final int tid) {
			Integer key=Integer.valueOf(tid);
			Clade c=cladeMap.get(key);
			if(c!=null) {
				assert(c.taxID==key.intValue()) : key+", "+c;
				return c;
			}
			synchronized(cladeMap) {
				c=cladeMap.get(key);
				if(c!=null) {
					assert(c.taxID==key.intValue()) : key+", "+c;
					return c;
				}
				c=Clade.makeClade(tid);
				assert(c.taxID==tid);
				Clade old=cladeMap.putIfAbsent(key, c);
				assert(old==null);
			}
			assert(cladeMap.get(key)==c) : key+", "+c;
			assert(cladeMap.get(key).taxID==key.intValue()) : key+", "+c;
			return c;
		}
		
		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		private EntropyTracker et=new EntropyTracker(entropyK, entropyWindow, false);
		
		private final Clade dummy=new Clade(-1, -1, null);
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Clade storage */
		private final ConcurrentHashMap<Integer, Clade> cladeMap;
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private ArrayList<String> in=new ArrayList<String>();

	private ArrayList<String> r16sFile=new ArrayList<String>();
	private ArrayList<String> r18sFile=new ArrayList<String>();
	int r16sAdded=0, r18sAdded=0;
	
	/** Primary output file path */
	private String out=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of lines out */
	protected long linesOut=0;
	/** Number of bytes out */
	protected long bytesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/** Whether to use a dummy clade for temporary storage */
	static boolean useDummy=true;
	
	/** Whether to merge duplicate tax IDs */
	static boolean mergeDuplicateTaxIDs=false;
	
	static boolean replaceRibo=false;
	static boolean useTree=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary output file */
	private FileFormat ffout;
	
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