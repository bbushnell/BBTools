package consensus;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import prok.ProkObject;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.SketchObject;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
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
import var2.Realigner;
import var2.SamFilter;

/**
 * Alters a reference to represent the consensus of aligned reads.
 * 
 * @author Brian Bushnell
 * @date September 6, 1019
 *
 */
public class ConsensusMaker extends ConsensusObject implements Accumulator<ConsensusMaker.ProcessThread> {
	
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
		ConsensusMaker x=new ConsensusMaker(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public ConsensusMaker(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
		samFilter.includeUnmapped=false;
//		samFilter.includeSupplimentary=false;
//		samFilter.includeDuplicate=false;
		samFilter.includeNonPrimary=false;
		samFilter.includeQfail=false;
//		samFilter.minMapq=4;
		
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
			silent=Parser.silent;
		}
		
		{
//			if("auto".equalsIgnoreCase(atomic)){Scaffold.setCA3A(Shared.threads()>8);}
//			else{Scaffold.setCA3A(Parse.parseBoolean(atomic));}

			if(ploidy<1){System.err.println("WARNING: ploidy not set; assuming ploidy=1."); ploidy=1;}
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
		ffmodel=FileFormat.testOutput(outModel, FileFormat.ALM, null, true, overwrite, false, ordered);
		
		//Create input FileFormat objects
		ffin=FileFormat.testInput(in, FileFormat.SAM, extin, true, true);
		ffref=FileFormat.testInput(ref, FileFormat.FASTA, null, true, true);
		
		if(inModelFile!=null){
			ArrayList<BaseGraph> list=(ArrayList<BaseGraph>)ReadWrite.readObject(inModelFile, false);
			inModel=list.get(0);
			inModel.calcProbs();
			inModel.makeWeights();
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
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				BaseGraph.verbose=verbose;
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("outm") || a.equals("outmodel") || a.equals("model") || a.equals("alm")){
				outModel=b;
			}else if(a.equals("inm") || a.equals("inmodel")){
				inModelFile=b;
			}else if(a.equals("hist") || a.equals("histogram")){
				outHistFile=b;
			}else if(a.equals("realign")){
				realign=Parse.parseBoolean(b);
			}else if(a.equals("printscores")){
				printScores=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("useMapq")){
				useMapq=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("identityCeiling") || a.equalsIgnoreCase("idceiling") || a.equalsIgnoreCase("ceiling")){
				double d=Double.parseDouble(b);
				if(d<=2){d=d*100;}
				identityCeiling=(int)d;
				invertIdentity=true;
			}else if(a.equalsIgnoreCase("invertIdentity")){
				invertIdentity=Parse.parseBoolean(b);
			}else if(a.equals("name") || a.equals("rename") || a.equals("header")){
				name=b;
			}else if(a.equals("noindels")){
				noIndels=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("onlyConvertNs") || a.equalsIgnoreCase("nOnly") || a.equalsIgnoreCase("onlyN")){
				onlyConvertNs=Parse.parseBoolean(b);
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("mindepth")){
				minDepth=Integer.parseInt(b);
			}else if(a.equals("trimdepth") || a.equals("trimdepthfraction")){
				trimDepthFraction=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("trimNs")){
				trimNs=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("mafn") || a.equals("mafnoref")){
				MAF_noref=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("mafsub")){
				MAF_sub=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("mafdel")){
				MAF_del=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("mafins")){
				MAF_ins=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("maf")){
				MAF_ins=MAF_noref=MAF_del=MAF_sub=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("mafindel")){
				MAF_ins=MAF_del=Float.parseFloat(b);
			}else if(a.equals("clearfilters")){
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
		
		if(ref!=null && !Tools.existsInput(ref)){
			specialRef=ProkObject.isSpecialType(ref);
		}
		
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in=Tools.fixExtension(in);
		if(!specialRef){ref=Tools.fixExtension(ref);}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){

		//Ensure there is an input file
		if(in==null){throw new RuntimeException("Error - an input file is required.");}

		//Ensure there is an input file
		if(ref==null){throw new RuntimeException("Error - a reference file is required.");}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out, outModel, outHistFile)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+" or "+outModel+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(true, true, in, inModelFile, (specialRef ? null : ref))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in, ref, out, outModel, outHistFile, inModelFile)){
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
	
	private void writeHist(String fname){
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, false);
		bsw.start();
		bsw.print("#Value\tIdentity");
		if(inModel!=null){
			bsw.print("\tScore");
		}
		bsw.nl();
		for(int i=0; i<idHist.length; i++){
			bsw.print(0.01f*i, 2).tab().print(idHist[i]);
			if(inModel!=null){
				bsw.tab().print(scoreHist[i]);
			}
			bsw.nl();
		}
		errorState=bsw.poisonAndWait()|errorState;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Load reference;
		refMap=loadReferenceCustom();
		makeRefMap2();
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		alignedReads=0;

		if(ffin.samOrBam()){
			//Create a read input stream
			final SamStreamer ss=makeStreamer(ffin);
			//Process the reads in separate threads
			spawnThreads(ss);
		}else{
			Shared.capBufferLen(40);
			//Create a read input stream
			final ConcurrentReadInputStream cris=makeCris(ffin);
			//Process the reads in separate threads
			spawnThreads(cris);
		}
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros();
		
		outputConsensus(ros);
		
		if(outHistFile!=null){
			writeHist(outHistFile);
		}
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStream(ros);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		if(!silent){
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
			outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
			outstream.println();
		}
		outstream.println(Tools.number("Ref Count:", refCount, 8));
		outstream.println(Tools.number("Sub Count:", subCount, 8));
		outstream.println(Tools.number("Del Count:", delCount, 8));
		outstream.println(Tools.number("Ins Count:", insCount, 8));
		outstream.println(Tools.number("Avg Identity:", 100*identitySum/alignedReads, 3, 8));
		if(scoreSum>0){
			outstream.println(Tools.number("Avg Score:", 100*scoreSum/alignedReads, 3, 8));
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	private synchronized LinkedHashMap<String, BaseGraph> loadReferenceCustom(){
		assert(!loadedRef);
		if(specialRef){return loadReferenceSpecial();}
		ConcurrentReadInputStream cris=makeRefCris();
		LinkedHashMap<String, BaseGraph> map=new LinkedHashMap<String, BaseGraph>();
		ListNum<Read> ln=cris.nextList();
		for(; ln!=null && !ln.isEmpty(); ln=cris.nextList()){
			if(verbose){outstream.println("Fetched "+ln.size()+" reference sequences.");}
			for(Read r : ln){
				if(r.bases!=null && r.bases.length>0){
					BaseGraph bg=new BaseGraph(r.id, r.bases, r.quality, r.numericID, 0);
					map.put(r.name(), bg);
//					if(inModel!=null && inModel.name.equals(bg.name)){
////						assert(false) : Arrays.toString(bg.refWeights)+"\n"+Arrays.toString(inModel.refWeights);
//						bg.refWeights=inModel.refWeights;
//						bg.insWeights=inModel.insWeights;
//						bg.delWeights=inModel.delWeights;
//					}else{
//						bg.insWeights=new float[bg.ref.length];
//						bg.delWeights=new float[bg.ref.length];
//						Arrays.fill(bg.insWeights, 1);
//						Arrays.fill(bg.delWeights, 1);
////						assert(false) : Arrays.toString(bg.delWeights);
//					}
				}
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		if(verbose){outstream.println("Loaded "+map.size()+" reference sequences.");}
		loadedRef=true;
		return map;
	}
	
	//For ribo subunits in resources directory
	private synchronized LinkedHashMap<String, BaseGraph> loadReferenceSpecial(){
		assert(!loadedRef);
		Read[] array=ProkObject.loadConsensusSequenceType(ref, false, false);
		Read r=array[0];
		BaseGraph bg=new BaseGraph(r.id, r.bases, r.quality, r.numericID, 0);
		LinkedHashMap<String, BaseGraph> map=new LinkedHashMap<String, BaseGraph>();
		map.put(r.name(), bg);
		return map;
	}

	private synchronized void makeRefMap2(){
		assert(refMap!=null && refMap2==null);
		if(verbose){outstream.println("Making refMap2.");}
		refMap2=new LinkedHashMap<String, BaseGraph>((refMap.size()*3)/2);
		for(Entry<String, BaseGraph> e : refMap.entrySet()){
			String key=e.getKey();
			if(verbose){outstream.println("Considering "+key);}
			BaseGraph value=e.getValue();
			String key2=Tools.trimToWhitespace(key);
			if(verbose){outstream.println("key2="+key2);}
//			if(verbose){outstream.println("put "+key2+", "+value);}
			refMap2.put(key2, value);
//			if(verbose){outstream.println("putted "+key2+", "+value);}
			
			if(defaultRname==null){defaultRname=key;}
		}
//		assert(false) : refMap+"\n"+refMap2;
		if(verbose){outstream.println("Made refMap2.");}
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
	
	private ConcurrentReadInputStream makeCris(FileFormat ff){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		assert(!paired);
		return cris;
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
		if(verbose){outstream.println("Started worker threads.");}
		
		//Wait for threads to finish
		boolean success=ThreadWaiter.waitForThreadsToFinish(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started worker threads.");}
		
		//Wait for threads to finish
		boolean success=ThreadWaiter.waitForThreadsToFinish(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		ReadWrite.closeStreams(cris);
	}
	
	private void outputConsensus(ConcurrentReadOutputStream ros){
		if(verbose){outstream.println("Making consensus.");}
		ArrayList<Read> consensusList=new ArrayList<Read>(200);
		ArrayList<BaseGraph> graphList=new ArrayList<BaseGraph>(200);
		long num=0;
		for(Entry<String, BaseGraph> e : refMap.entrySet()){
			BaseGraph bg=e.getValue();
			graphList.add(bg);
			
			Read r=bg.traverse();
			if(name!=null){r.id=name;}
			consensusList.add(r);
			
			refCount+=bg.refCount;
			subCount+=bg.subCount;
			delCount+=bg.delCount;
			insCount+=bg.insCount;
			
			readsOut++;
			basesOut+=r.length();
			
			if(consensusList.size()>=200){
				if(ros!=null){ros.add(consensusList, num);}
				consensusList=new ArrayList<Read>(200);
				num++;
			}
		}
		if(consensusList.size()>0){
			if(ros!=null){ros.add(consensusList, num);}
			consensusList=new ArrayList<Read>(200);
			num++;
		}
		if(graphList.size()>0 && ffmodel!=null){
//			ReadWrite.writeObjectInThread(graphList, ffmodel.name(), true);
			ReadWrite.writeAsync(graphList, ffmodel.name(), true);
		}
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		readsProcessed+=pt.readsProcessedT;
		basesProcessed+=pt.basesProcessedT;
		alignedReads+=pt.alignedReadsT;
		identitySum+=pt.identitySumT;
		scoreSum+=pt.scoreSumT;
		for(int i=0; i<idHist.length; i++){
			idHist[i]+=pt.idHistT[i];
			scoreHist[i]+=pt.scoreHistT[i];
		}
		errorState|=(!pt.success);
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public ByteBuilder toText() {
		// TODO Auto-generated method stub
		return null;
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
			cris=null;
			tid=tid_;
			realigner=(realign ? new Realigner() : null);
		}
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final int tid_){
			ss=null;
			cris=cris_;
			tid=tid_;
			realigner=null;//(realign ? new Realigner() : null);
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
			
			if(ss!=null){
				for(ListNum<Read> ln=ss.nextReads(); ln!=null; ln=ss.nextReads()){
					processList(ln);
				}
			}else{
				for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
					processList(ln);
					cris.returnList(ln);
				}
			}
			
		}
		
		void processList(ListNum<Read> ln){
			
			//Loop through each read in the list
			for(Read r : ln){
				
				//Validate reads in worker threads
				if(!r.validated()){r.validate(true);}

				//Track the initial length for statistics
				final int initialLength=r.length();

				//Increment counters
				readsProcessedT++;
				basesProcessedT+=initialLength;
				
				{
					//Reads are processed in this block.
					processRead(r);
				}
			}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		void processRead(final Read r){
			if(r.bases==null || r.length()<=1){return;}
			
			final SamLine sl=r.samline;
			if(sl!=null && !sl.mapped()){return;}
			final String rname=(sl==null ? defaultRname : sl.rnameS());
			BaseGraph bg=refMap.get(rname);
			if(bg==null){bg=refMap2.get(Tools.trimToWhitespace(rname));}
			assert(bg!=null) : "Can't find graph for "+rname;
			float identity;
			float score=0;//unused

//			assert(false) : r;
			if(sl==null){
//				identity=SketchObject.alignAndMakeMatch(r, bg.original);
//				assert(false) : bg.refWeights[0];
//				identity=SketchObject.alignAndMakeMatch(r, bg.original, bg.refWeights, bg.insWeights, bg.delWeights);
				identity=SketchObject.alignAndMakeMatch(r, bg.original, bg.refWeights);
				if(identity<=0){return;}
				assert(r.match!=null) : identity+", "+r;
			}else{
				assert(sl!=null && sl.mapped() && sl.seq!=null) : sl;
				assert(r.match!=null);
				if(samFilter!=null && !samFilter.passesFilter(sl)){return;}
				identity=sl.calcIdentity();
			}
			assert(r.match!=null && identity>0) : identity+", "+r;
			
			if(inModel!=null){

				if(verbose){System.err.println(r.start+"\t"+r.stop+"\t"+new String(r.match));}
				score=inModel.score(r, false, true);
				if(printScores){
//					System.err.println(identity+"\t"+score+"\t"+inModel.score(r, false, true));
					System.err.println(String.format("%.5f\t%.5f", identity, score));
				}
//				assert(false);
//				if(score<0){
//					verbose=true;
//					inModel.score(r, false, false);
//				}
//				assert(!verbose);
			}
			
			alignedReadsT++;
			identitySumT+=identity;
			scoreSumT+=score;
			idHistT[Tools.mid(0, Math.round(100*identity), 100)]++;
			scoreHistT[Tools.mid(0, Math.round(100*score), 100)]++;
			
			
			
//			assert(sl.calcIdentity()<=0.78) : sl.calcIdentity()+", "+samFilter.maxId;
//			assert(false) : sl.cigar+", "+sl.calcIdentity();
			
//			System.err.println(rname);
			
			if(realigner!=null) {
				assert(false);
				realigner.realign(r, sl, bg.original, true);
			}
			
			bg.add(r);
		}

		long[] idHistT=new long[101];
		long[] scoreHistT=new long[101];

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		protected long alignedReadsT=0;

		double identitySumT=0;
		double scoreSumT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;

		/** Shared input stream */
		private final SamStreamer ss;
		/** Alternate input stream */
		private final ConcurrentReadInputStream cris;
		/** Thread ID */
		final int tid;
		/** For realigning reads */
		final Realigner realigner;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Read input file path */
	private String in=null;
	/** Reference input file path */
	private String ref=null;

	private String inModelFile;
	private BaseGraph inModel;
	private String outHistFile;
	
	/** Consensus output file path */
	private String out=null;
	/** Model output file path */
	private String outModel=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	
	protected long alignedReads=0;

	protected double identitySum=0;
	protected double scoreSum=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	public long subCount=0;
	public long refCount=0;
	public long delCount=0;
	public long insCount=0;

	long[] idHist=new long[101];
	long[] scoreHist=new long[101];
	boolean printScores=false;
	
	/*--------------------------------------------------------------*/
	
	/** Threads dedicated to reading the sam file */
	private int streamerThreads=SamStreamer.DEFAULT_THREADS;
	
	private String name=null;

	private boolean loadedRef=false;
	private boolean specialRef=false;
	
	private boolean realign=false;
	
	private int ploidy=1;
	
	private final boolean silent;
	
	public final SamFilter samFilter=new SamFilter();
	/** Uses full ref names */
	public LinkedHashMap<String, BaseGraph> refMap;
	/** Uses truncated ref names */
	public LinkedHashMap<String, BaseGraph> refMap2;
	
	private String defaultRname=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Read input file */
	private final FileFormat ffin;
	/** Reference input file */
	private final FileFormat ffref;

	/** Consensus output file */
	private final FileFormat ffout;
	/** Model output file */
	private final FileFormat ffmodel;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	/** Reads ARE output in input order, even though this is false */
	private final boolean ordered=false;
	
}
