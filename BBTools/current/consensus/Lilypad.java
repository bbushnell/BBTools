package consensus;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.AminoAcid;
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
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.SamReadStreamer;
import stream.SamStreamer;
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;
import var2.SamFilter;

/**
 * Scaffolds contigs based on paired read mapping.
 * 
 * @author Brian Bushnell
 * @date September 11, 2019
 *
 */
public class Lilypad implements Accumulator<Lilypad.ProcessThread> {
	
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
		Lilypad x=new Lilypad(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Lilypad(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		SamLine.RNAME_AS_BYTES=false;
		
		samFilter.includeUnmapped=false;
		samFilter.includeSupplimentary=false;
//		samFilter.includeDuplicate=false;
		samFilter.includeNonPrimary=false;
		samFilter.includeQfail=false;
		samFilter.minMapq=4;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in=parser.in1;
			extin=parser.extin;

			out=parser.out1;
			extout=parser.extout;
		}
		
		{
//			if("auto".equalsIgnoreCase(atomic)){Scaffold.setCA3A(Shared.threads()>8);}
//			else{Scaffold.setCA3A(Parse.parseBoolean(atomic));}
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
		ffin=FileFormat.testInput(in, FileFormat.SAM, extin, true, true);
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
			}else if(a.equals("ref") || a.equals("scaffolds")){
				ref=b;
			}else if(a.equals("insertlist")){
				insertList=b;
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("sameStrandPairs")){
				sameStrandPairs=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("ns") || a.equalsIgnoreCase("n") || a.equalsIgnoreCase("scaffoldbreak") || a.equalsIgnoreCase("gap") || a.equalsIgnoreCase("mingap")){
				scaffoldBreakNs=Integer.parseInt(b);
				assert(scaffoldBreakNs>0);
			}else if(a.equalsIgnoreCase("mindepth")){
				minDepth=Integer.parseInt(b);
				assert(minDepth>0);
			}else if(a.equalsIgnoreCase("maxinsert")){
				maxPairDist=Parse.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("minWeightRatio") || a.equalsIgnoreCase("minwr")){
				minWeightRatio=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minStrandRatio") || a.equalsIgnoreCase("minsr")){
				minStrandRatio=Float.parseFloat(b);
			}else if(a.equals("clearfilters") || a.equals("clearfilter")){
				if(Parse.parseBoolean(b)){
					samFilter.clear();
				}
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(samFilter.parse(arg, a, b)){
				//do nothing
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
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in, ref)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in, ref, out)){
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
		
		//Create a read input stream
		final SamStreamer ss=makeStreamer(ffin);
		
		//Load reference
		loadReferenceCustom();
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads(ss);

		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros();
		
		if(verbose){outstream.println("Fixing reference.");}
		
		makeScaffolds(ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStream(ros);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, scaffoldsOut, scaffoldLengthOut, 8, false));
		
		outstream.println();
		outstream.println(Tools.number("Average Insert", totalAverageInsert, 2, 8));
		outstream.println(Tools.number("Joins Made    ", gapsAdded, 8));
		outstream.println(Tools.number("Ns Added      ", nsAdded, 8));
		outstream.println(Tools.number("Contigs In    ", refMap.size(), 8));
		outstream.println(Tools.number("Scaffolds Out ", scaffoldsOut, 8));
		
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	private synchronized void loadReferenceCustom(){
		assert(!loadedRef);
		ConcurrentReadInputStream cris=makeRefCris();
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()) {
			for(Read r : ln){
				String name=r.id;
				String name2=Tools.trimToWhitespace(r.id);
				Contig cont=new Contig(name, r.bases, r.numericID);
				refMap.put(name, cont);
				refMap2.put(name2, cont);
			}
			cris.returnList(ln);
		}
		ReadWrite.closeStream(cris);
		loadedRef=true;
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
		SamStreamer ss=new SamReadStreamer(ff, streamerThreads, true, maxReads);
		ss.start(); //Start the stream
		if(verbose){outstream.println("Started Streamer");}
		return ss;
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
		totalAverageInsert=totalInsertSum/(double)totalInsertCount;
		insertByPercentile=Tools.makeHistogram(insertCounts, buckets);
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		readsOut+=pt.readsOutT;
		basesOut+=pt.basesOutT;
		
		totalInsertSum+=pt.totalInsertSumT;
		totalInsertCount+=pt.totalInsertCountT;
		
		errorState|=(!pt.success);
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void makeScaffolds(ConcurrentReadOutputStream ros){
		ByteBuilder bb=new ByteBuilder(1000000);

		ArrayList<Read> list=new ArrayList<Read>(200);
		long num=0;
		long lengthSum=0;
		for(Entry<String, Contig> e : refMap.entrySet()){
			Contig cont=e.getValue();
			if(!cont.processed()){
				Read r=cont.makeScaffold(bb);
				assert(r!=null);

				lengthSum+=r.length();
				list.add(r);
				scaffoldsOut++;
				scaffoldLengthOut+=r.length();

				if(list.size()>=200 || lengthSum>=100000){
					if(ros!=null){ros.add(list, num);}
					list=new ArrayList<Read>(200);
					num++;
					lengthSum=0;
				}
				assert(cont.processed());
			}
		}
		if(list.size()>0){
			if(ros!=null){ros.add(list, num);}
		}
	}
	
	private static int calcInsertSize(SamLine sl) {
		assert(sl.mapped() && sl.pairedOnSameChrom());
		assert(sl.primary());
		assert(!sl.supplementary());
		assert(sl.leftmost());
		
		assert(sl.tlen>0) : sl.tlen+"\n\n"+sl;
		return sl.tlen>0 ? sl.tlen : -sl.tlen;
		
//		final int insertSize;
//		String insertTag=null;
//		if(sl.optional!=null){
//			for(String s : sl.optional){
//				if(s.startsWith("X8:Z:")){
//					insertTag=s;
//					break;
//				}
//			}
//		}
//		if(insertTag!=null){
//			insertSize=Integer.parseInt(insertTag.substring(5));
//		}else{
//			insertSize=sl.tlen;//This is unsafe due to indels.
//			assert(false) : "Reads need insert size tags.";
//		}
//		assert(insertSize>0) : sl;
//		
//		return insertSize;
	}
	
	private Contig getScaffold(String rname){
		Contig scaf=refMap.get(rname);
		if(scaf==null){scaf=refMap2.get(Tools.trimToWhitespace(rname));}
		assert(scaf!=null) : "Can't find graph for "+rname;
		return scaf;
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
			for(ListNum<Read> ln=ss.nextReads(); ln!=null; ln=ss.nextReads()){
//				if(verbose){outstream.println("Got list of size "+list.size());} //Disabled due to non-static access
				
				processList(ln);
			}
			
		}
		
		void processList(ListNum<Read> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r=reads.get(idx);
				
				//Validate reads in worker threads
				if(!r.validated()){r.validate(true);}

				//Track the initial length for statistics
				final int initialLength=r.length();

				//Increment counters
				readsProcessedT+=r.pairCount();
				basesProcessedT+=initialLength;
				
				processRead(r);
			}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		void processRead(final Read r){
			final SamLine sl=r.samline;
			assert(sl!=null) : sl;
			if(samFilter!=null && !samFilter.passesFilter(sl)){return;}
			
			//sl.nextMapped();
			if(sl.mapped() && sl.primary() && !sl.supplementary()){
				final String rname=sl.rnameS();
				Contig scaf=getScaffold(rname);
				if(scaf!=null){
					if(sl.pairedOnSameChrom() && sl.properPair() && sl.leftmost()){
						final int insertSize=calcInsertSize(sl);
						insertCounts.incrementAndGet(Tools.mid(0, insertSize, insertCounts.length()-1));
						totalInsertSumT+=insertSize;
						totalInsertCountT++;
					}
					scaf.add(sl);

					readsOutT++;
					basesOutT+=r.length();
				}
			}
			if(sl.mapped() && sl.pairedOnSameChrom() && sl.properPair() && sl.primary() && !sl.supplementary() && sl.leftmost()){
				final String rname=sl.rnameS();
				
				Contig scaf=getScaffold(rname);
				if(scaf!=null){
					final int insertSize=calcInsertSize(sl);
					insertCounts.incrementAndGet(Tools.mid(0, insertSize, insertCounts.length()-1));
					scaf.add(sl);

					readsOutT++;
					basesOutT+=r.length();

					totalInsertSumT+=insertSize;
					totalInsertCountT++;
				}
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		protected long totalInsertSumT=0;
		protected long totalInsertCountT=0;
		
		long insertSum=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final SamStreamer ss;
		/** Thread ID */
		final int tid;
	}
	
	Contig findLeftmost(Contig source){
		if(verbose){System.err.println("findLeftmost("+source.name+")");}
		while(true) {
			assert(!source.processed());
			if(source.processed()){return null;}
			source.processedLeft=true;
			Edge se=source.bestEdge(true);
			if(verbose){System.err.println("Found source edge "+se);}
			if(se==null){return source;}
			Contig dest=se.dest;
			if(dest.processed()){
				if(verbose){System.err.println("Dest was processed; returning.");}
				return source;
			}
			if(se.sameStrand()){
				if(source.strand==dest.strand){

				}else{
					if(verbose){System.err.println("Flipping "+dest.name);}
					dest.flip();
				}
			}else{
				if(source.strand==dest.strand){
					if(verbose){System.err.println("Flipping "+dest.name);}
					dest.flip();
				}else{

				}
			}
			Edge de=dest.bestEdge(false);
			if(verbose){System.err.println("Found dest edge "+de);}
			if(de==null || de.dest!=source){
				if(dest.strand==1){dest.flip();}
				if(verbose){System.err.println("Dest edge did not match; returning.");}
				return source;
			}
			source=dest;
			if(verbose){System.err.println("Migrated to next node.");}
		}
	}
	
	Read expandRight(final Contig source0, ByteBuilder bb){
		if(verbose){System.err.println("expandRight("+source0.name+")");}
		bb.clear();
		Contig source=source0;
		while(true) {
			assert(!source.processedRight);
			if(source.processedRight){return null;}
			if(source.strand==1){
				Tools.reverseInPlace(source.depthArray);
			}
			source.processedRight=true;
			bb.append(source.bases);
			Edge se=source.bestEdge(false);
			if(verbose){System.err.println("Found source edge "+se);}
			if(se==null){break;}
			Contig dest=se.dest;
			if(dest.processedRight){
				if(verbose){System.err.println("Dest was processed; returning.");}
				break;
			}
			if(se.sameStrand()){
				if(source.strand==dest.strand){

				}else{
					if(verbose){System.err.println("Flipping "+dest.name);}
					dest.flip();
				}
			}else{
				if(source.strand==dest.strand){
					if(verbose){System.err.println("Flipping "+dest.name);}
					dest.flip();
				}else{

				}
			}
			Edge de=dest.bestEdge(true);
			if(verbose){System.err.println("Found dest edge "+de);}
			if(de==null || de.dest!=source){
				if(verbose){System.err.println("Dest edge did not match; returning.");}
				if(dest.strand==1){dest.flip();}
				break;
			}
			
			//Now append Ns
			int observedLength=(int)(se.distanceSum/se.count());
			long depth=se.count();
			int depthProxyIndex=(source.length()-Tools.min(source.length()/2, 300));
			long depthProxy=source.depthArray.get(depthProxyIndex);
			int percentile=(int)(buckets*depth/(float)(depth+depthProxy));
			int inferredLength=insertByPercentile[percentile];
			
			int Ns=(Tools.max(scaffoldBreakNs, inferredLength-observedLength));
			for(int i=0; i<Ns; i++){bb.append('N');}
			source=dest;
			
			gapsAdded++;
			nsAdded+=Ns;
		}
		Read r=new Read(bb.toBytes(), null, source0.name, source0.numericID);
		return r;
	}
	
	/*--------------------------------------------------------------*/
	
	private class Contig {
		
		Contig(String name_, byte[] bases_, long numericID_){
			name=name_;
			bases=bases_;
			numericID=(int)numericID_;
			depthArray=new AtomicIntegerArray(bases.length);
		}
		
		public Read makeScaffold(ByteBuilder bb) {
			assert(!processed());
			Contig leftmost=findLeftmost(this);
			return expandRight(leftmost, bb);
		}
		
		Edge bestEdge(boolean left) {
			final LinkedHashMap<String, Edge> map=(left ? leftEdgeMap : rightEdgeMap);
			if(map.isEmpty()){return null;}
			long weightSum=0;
			long countSum=0;
			Edge best=null;
			for(Entry<String, Edge> entry : map.entrySet()){
				Edge e=entry.getValue();
				weightSum+=e.weight;
				countSum+=e.count();
				if(best==null || e.weight>best.weight){best=e;}
			}
			if(best.count()<minDepth){return null;}
			if(weightSum*minWeightRatio>best.weight){return null;}
			if(best.strandRatio()<minStrandRatio){return null;}
			if(best.badCount>0.5*best.count()){return null;}
			return best;
		}

		void add(SamLine sl){
			assert(sl.mapped() && sl.primary() && !sl.supplementary());
			if(sl.nextMapped()){
				if(sl.pairedOnSameChrom()){
					if(!sl.properPair()){
						addCoverageSingleton(sl);
					}else if(sl.leftmost()){
						addCoveragePaired(sl);
					}
				}else{
					addCoverageSingleton(sl);
					handleMixedPair(sl);
				}
			}else{
				addCoverageSingleton(sl);
			}
		}
		
		private void addCoverageSingleton(SamLine sl){
			assert(sl.cigar!=null);
			int start=sl.pos-1;
			int stop=start+SamLine.calcCigarLength(sl.cigar, false, false);
			
			for(int i=start; i<stop; i++){
				if(i>=0 && i<bases.length){
					depthArray.incrementAndGet(i);
				}
			}
		}
		
		private void addCoveragePaired(SamLine sl){
			assert(sl.cigar!=null);
			assert(sl.leftmost() && sl.pairedOnSameChrom() && sl.nextMapped());
			int start=sl.pos-1;
			int stop=start+sl.tlen;
			
			for(int i=start; i<stop; i++){
				if(i>=0 && i<bases.length){
					depthArray.incrementAndGet(i);
				}
			}
		}
		
		/** Reads mapping to different contigs */
		private void handleMixedPair(SamLine sl){
			assert(sl.mapped() && sl.nextMapped() && !sl.pairedOnSameChrom());
			String rname=sl.rnameS();
			String rnext=sl.rnextS();
			if(rname.equals(rnext)){return;}
			
			final boolean left=(sl.strand()==1);
			LinkedHashMap<String, Edge> map=(left ? leftEdgeMap : rightEdgeMap);
			Edge e=map.get(rnext);
			Contig dest=null;
			if(e==null){
				dest=getScaffold(rnext);
				if(dest==null){return;}
				e=new Edge(this, dest, left);
				map.put(rnext, e);
			}
			e.add(sl);
		}
		
		void flip(){//Be careful with this
			AminoAcid.reverseComplementBasesInPlace(bases);
			strand^=1;
			LinkedHashMap<String, Edge> temp=leftEdgeMap;
			leftEdgeMap=rightEdgeMap;
			rightEdgeMap=temp;
		}
		
		int length(){return bases.length;}
		
		final int numericID;
		final String name;
		final byte[] bases;
		final AtomicIntegerArray depthArray;
		int strand=0;

		boolean processedLeft=false;
		boolean processedRight=false;
		boolean processed(){return processedLeft || processedRight;}

		LinkedHashMap<String, Edge> leftEdgeMap=new LinkedHashMap<String, Edge>();
		LinkedHashMap<String, Edge> rightEdgeMap=new LinkedHashMap<String, Edge>();
	}
	
	private class Edge{
		
		Edge(Contig source_, Contig dest_, boolean left_){
			source=source_;
			dest=dest_;
			leftEdge=left_;
		}

		void add(SamLine sl){
			final boolean sameStrandReads=(sl.strand()==sl.mateStrand());
			final boolean sameStrandContigs=(sameStrandPairs==sameStrandReads);
			final int spos, dpos;
			if(leftEdge){
				spos=sl.pos+sl.calcCigarLength(true, false)-1;
				dpos=(sameStrandContigs ? dest.length()-sl.pnext : sl.pnext+sl.length())-1;
			}else{
				spos=source.length()-sl.pos-1;
				dpos=(sameStrandContigs ? sl.pnext+sl.length() : dest.length()-sl.pnext)-1;
			}
			final int distance=spos+dpos;
			
			if(distance>maxPairDist){

//				assert(false) : "distance="+distance+", spos="+spos+", dpos="+dpos+", sameStrandContigs="+sameStrandContigs+
//					"\nsl.pos="+sl.pos+", sl.pnext="+sl.pnext+", strand="+sl.strand()+", nextStrand="+sl.mateStrand()+", left="+leftEdge
//					+"\n"+sl;
//				badCount++;
				return;
			}
			
			distanceSum+=distance;
			
			weight+=sl.mapq;
			if(sameStrandContigs){
				sameStrandCount++;
			}else{
				difStrandCount++;
			}
//			assert(false) : weight;
		}
		
		public float strandRatio() {
			return Tools.max(sameStrandCount, difStrandCount)/(float)(sameStrandCount+difStrandCount);
		}
		
		public boolean sameStrand(){
			return sameStrandCount>=difStrandCount;
		}
		
		@Override
		public String toString(){
			return "("+source.name+"->"+dest.name+", "+(leftEdge ? "left" : "right")+", weight="+weight+
					", same="+sameStrandCount+", dif="+difStrandCount+", bad="+badCount+")";
		}
		
		long count(){return sameStrandCount+difStrandCount;}
		
		final Contig source;
		final Contig dest;
		long sameStrandCount;
		long difStrandCount;
		long distanceSum;
		long weight;
		long badCount;
		final boolean leftEdge;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in=null;
	/** Secondary input file path */
	private String ref=null;

	/** Primary output file path */
	private String out=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	private String insertList=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	protected long scaffoldsOut=0;
	protected long scaffoldLengthOut=0;
	
	protected long totalInsertSum=0;
	protected long totalInsertCount=0;
	protected double totalAverageInsert;
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	boolean sameStrandPairs=false;
	
	int gapsAdded=0;
	long nsAdded=0;
	
	/*--------------------------------------------------------------*/
	
	/** Threads dedicated to reading the sam file */
	private int streamerThreads=SamStreamer.DEFAULT_THREADS;
	
	private boolean loadedRef=false;
	
	private int minDepth=4;

	private float minWeightRatio=0.8f;
	private float minStrandRatio=0.8f;
	
	private int scaffoldBreakNs=10;
	
	private int maxPairDist=3000;
	
	private int buckets=1000;
	protected AtomicLongArray insertCounts=new AtomicLongArray(20000);
	protected int[] insertByPercentile;
	
	public final SamFilter samFilter=new SamFilter();
	
	/** Uses full ref names */
	public LinkedHashMap<String, Contig> refMap=new LinkedHashMap<String, Contig>();
	/** Uses truncated ref names */
	public LinkedHashMap<String, Contig> refMap2=new LinkedHashMap<String, Contig>();
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin;
	/** Secondary input file */
	private final FileFormat ffref;
	
	/** Primary output file */
	private final FileFormat ffout;
	
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
