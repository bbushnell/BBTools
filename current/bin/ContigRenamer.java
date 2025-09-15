package bin;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.LineParser1;
import shared.LineParser4;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.SamLineStreamer;
import stream.SamReadInputStream;
import stream.SamStreamer;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntList;
import structures.IntLongHashMap;
import structures.ListNum;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;
import var2.SamFilter;

/**
 * Renames contigs based on a sam file.
 * Appends coverage (cov_#).
 * Also appends taxid (taxid_#) if reads have taxids,
 * and the contig does not already have a taxid.
 * 
 * @author Brian Bushnell
 * @date September 6, 2019
 *
 */
public class ContigRenamer implements Accumulator<ContigRenamer.ProcessThread> {
	
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
		ContigRenamer x=new ContigRenamer(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public ContigRenamer(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		samFilter.includeUnmapped=false;
		samFilter.includeSupplimentary=false;
		samFilter.includeDuplicate=false;
		samFilter.includeNonPrimary=false;
		samFilter.includeQfail=false;
//		samFilter.minMapq=4;
		
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
		
		{
			samFilter.setSamtoolsFilter();
			
			streamerThreads=Tools.max(1, Tools.min(streamerThreads, Shared.threads()));
			assert(streamerThreads>0) : streamerThreads;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout=FileFormat.testOutput(out, FileFormat.FASTA, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffref=FileFormat.testInput(ref, FileFormat.FASTA, null, true, true);
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
			}else if(a.equals("ref") || a.equals("contigs") || a.equals("in")){
				ref=b;
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("delimiter")){
				delimiter=Parse.parseSymbol(b);
			}else if(a.equals("wipe")){
				wipe=Parse.parseBoolean(b);
			}else if(a.equals("depth") || a.equals("cov") || a.equals("adddepth")){
				addDepth=Parse.parseBoolean(b);
			}else if(a.equals("tid") || a.equalsIgnoreCase("taxid")){
				addTid=Parse.parseBoolean(b);
			}else if(a.equals("cami")){
				cami=b;
			}else if(a.equals("clearfilters")){
				if(Parse.parseBoolean(b)){
					samFilter.clear();
				}
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(a.equals("sam") || a.equals("insam") || a.equals("samin")){
				Tools.getFileOrFiles(b, in, false, false, true, false);
			}else if(samFilter.parse(arg, a, b)){
				//do nothing
			}else if(new File(arg).isFile()){
				in.add(arg);
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
		ref=Tools.fixExtension(ref);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){

		//Ensure there is an input file
		if(in==null){throw new RuntimeException("Error - an input file is required.");}

		//Ensure there is an input file
		if(ref==null){throw new RuntimeException("Error - a reference file is required.");}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		ArrayList<String> x=new ArrayList<String>();
		x.addAll(in);
		x.add(ref);
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, x)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		x.add(out);
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, x.toArray(new String[0]))){
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
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		if(in.size()>0) {
			processSam();
		}else if(cami!=null) {
			processCami();
		}else {
			assert(false) : "No input.";
		}
		processReference();
		
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
	
	private void processSam() {

		for(int i=0; i<in.size(); i++) {
			FileFormat ffin=FileFormat.testInput(in.get(i), FileFormat.SAM, extin, true, true);

			//Create a read input stream
			SamStreamer ss=makeStreamer(ffin);

			if(scafMap.isEmpty()) {scafMap=makeScafMap();}

			//Process the reads in separate threads
			spawnThreads(ss);
			
			for(Scaf scaf : scafMap.values()) {
				scaf.process();
				scaf.map.clear();
			}
		}
	}
	
	private void processCami() {
		LineParser1 lp=new LineParser1('\t');
		ArrayList<byte[]> lines=ByteFile.toLines(cami);
		scafMap=new HashMap<String, Scaf>();
		addDepth=false;
		addTid=true;
		for(byte[] line : lines) {
			if(!Tools.startsWith(line, '@')){
				lp.set(line);
				String name=lp.parseString(0);
				int taxID=lp.parseInt(2);
				Scaf scaf=new Scaf(name, 1);
				Scaf old=scafMap.putIfAbsent(name, scaf);
				assert(old==null) : "Name conflict for scaffold "+name;
				scaf.tidList.add(taxID);
				scaf.covList.add(1);
			}
		}
	}
	
	private ConcurrentReadInputStream makeRefCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffref, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		assert(!paired) : "References should not be paired.";
		return cris;
	}
	
	private SamStreamer makeStreamer(FileFormat ff){
		if(ff==null){return null;}
		SamStreamer ss=new SamLineStreamer(ff, streamerThreads, true, maxReads);
		ss.start(); //Start the stream
		if(verbose){outstream.println("Started Streamer");}
		return ss;
	}
	
	private static HashMap<String, Scaf> makeScafMap() {
		//@SQ	SN:foo	LN:999
		LineParser4 lp=new LineParser4("\t:\t:");
		ArrayList<byte[]> header=SamReadInputStream.getSharedHeader(true);
		HashMap<String, Scaf> map=new HashMap<String, Scaf>(header.size()+1);
		for(byte[] line : header) {
			if(Tools.startsWith(line, "@SQ\t")){
				lp.set(line);
				assert(lp.termEquals("SN", 1)) : "\n"+lp.parseString(1)+"\n"+lp+"\n"+new String(line);
				String name=lp.parseString(2);
				assert(lp.termEquals("LN", 3));
				int len=lp.parseInt(4);
				Scaf scaf=new Scaf(name, len);
				Scaf old=map.putIfAbsent(name, scaf);
				assert(old==null) : "Name conflict for scaffold "+name;
			}
		}
		return map;
	}
	
	private ConcurrentReadOutputStream makeCros(){
		if(ffout==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final SamStreamer ss){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(ss, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for threads to finish
		boolean success=ThreadWaiter.waitForThreadsToFinish(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
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

	private void processReference(){
		ConcurrentReadInputStream cris=makeRefCris();
		ConcurrentReadOutputStream ros=makeCros();
		ListNum<Read> ln=null;
		for(ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()) {
			for(Read r : ln) {
				r.id=(rename(r.id, r.length(), r.numericID));
				readsOut++;
				basesOut+=r.length();
				readsProcessed++;
				basesProcessed+=r.length();
			}
			if(ros!=null) {ros.add(ln.list, ln.id);}
			cris.returnList(ln);
		}

		//Notify the input stream that the list was used
		cris.returnList(ln);

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		errorState|=ReadWrite.closeStreams(cris, ros);
	}
	
	static String toShortName(String s) {
		int idx=Tools.indexOfWhitespace(s);
		return idx<0 ? s : s.substring(0, idx);
	}
	
	static String toShortName(byte[] s) {
		int idx=Tools.indexOfWhitespace(s);
		return idx<0 ? new String(s) : new String(s, 0, idx);
	}
	
	private String rename(final String old, final int scaflen, long numericID) {
		final String name=old;
		final String shortName=toShortName(name);
		final int oldTaxid=TaxTree.parseHeaderStatic2(old, null);
		Scaf scaf=scafMap.get(shortName);
		if(scaf==null) {scaf=scafMap.get(name);}
		assert(scaf!=null) : "Can't find "+shortName+" in map of size "+scafMap.size();
		if(scaf==null) {return old+delimiter+"cov_0";}
		
		IntList tids=scaf.tidList;
		FloatList covs=scaf.covList;
		float maxCov=0;
		int maxTid=-1;
		for(int i=0; i<covs.size(); i++) {
			int tid=tids.get(i);
			float cov=covs.get(i);
			if(maxTid<=0 || cov>maxCov) {
				maxCov=cov;
				maxTid=tid;
			}
		}
		
		ByteBuilder bb=new ByteBuilder(old.length()+(covs.size*9)+20);
		if(wipe) {
			bb.append("contig_").append(numericID);
		}else {
			bb.append(old);
		}
		
		if(oldTaxid<0 && maxTid>0 && addTid) {bb.append(delimiter).append("tid_").append(maxTid);}
		for(int i=0; i<covs.size && addDepth; i++) {
			bb.append(delimiter).append("cov_").append(covs.get(i), 2);
		}
//		assert(false) : oldTaxid+", "+maxTid+", "+bb;
		return bb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final SamStreamer ss_, final int tid_){
			ss=ss_;
			tid=tid_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab and process all lists
			for(ListNum<SamLine> ln=ss.nextLines(); ln!=null; ln=ss.nextLines()){
//				if(verbose){outstream.println("Got list of size "+list.size());} //Disabled due to non-static access
				
				processList(ln);
			}
			
		}
		
		void processList(ListNum<SamLine> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<SamLine> reads=ln.list;
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final SamLine r=reads.get(idx);

				//Track the initial length for statistics
				final int initialLength=r.length();

				//Increment counters
				readsProcessedT+=r.length();
				basesProcessedT+=initialLength;
				
				addLine(r);
			}
		}
		
		private void addLine(SamLine sl) {
			if(samFilter!=null && !samFilter.passesFilter(sl)){return;}
			String key=sl.rnameS();
			Scaf scaf=scafMap.get(key);
			assert(scaf!=null) : "Can't find scaffold "+key+" in map of size "+scafMap.size();
			scaf.add(sl);
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final SamStreamer ss;
		/** Thread ID */
		final int tid;
	}
	
	private static class Scaf {
		
		Scaf(String shortName_, int length_) {
			name=shortName_;
			length=length_;
		}
		
		void add(SamLine sl) {
			int taxid=TaxTree.parseHeaderStatic2(sl.qname, null);
			int length=sl.length();
			synchronized(this) {map.increment(taxid, length);}
		}
		
		void process() {
			
			final int[] keys=map.keys();
			final long[] values=map.values();
			final int invalid=map.invalid();
			int taxid=-1;
			long totalMappedBases=0;
			long maxMappedBases=0;
			for(int i=0; i<keys.length; i++) {
				int key=keys[i];
				long value=values[i];
				if(key!=invalid) {
					totalMappedBases+=value;
					if(value>maxMappedBases) {
						maxMappedBases=value;
						taxid=key;
					}
				}
			}
			assert(length>0);
			final float cov=totalMappedBases/(float)length;
			tidList.add(taxid);
			covList.add(cov);
		}
		
		final String name;
		final int length;
		
		IntLongHashMap map=new IntLongHashMap(7);
		IntList tidList=new IntList(2);
		FloatList covList=new FloatList(2);
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Sam input file path */
	private ArrayList<String> in=new ArrayList<String>();
	/** Contig input file path */
	private String ref=null;

	private String cami=null;
	
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

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;

	private String delimiter=" ";
	private boolean wipe=false;
	private boolean addDepth=true;
	private boolean addTid=true;
	
	/*--------------------------------------------------------------*/
	
	/** Threads dedicated to reading the sam file */
	private int streamerThreads=SamStreamer.DEFAULT_THREADS;
	
	public final SamFilter samFilter=new SamFilter();
	
	public HashMap<String, Scaf> scafMap=new HashMap<String, Scaf>();
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Secondary input file */
	private final FileFormat ffref;
	
	/** Primary output file */
	private final FileFormat ffout;
	
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
