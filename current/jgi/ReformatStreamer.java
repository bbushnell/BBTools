package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import processor.ReformatProcessor;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FASTQ;
import stream.Read;
import stream.ReadStreamByteWriter;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Reformat using Streamer/Writer interfaces with multithreading.
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 15, 2025
 */
public class ReformatStreamer implements Accumulator<ReformatStreamer.ProcessThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		ReformatStreamer x=new ReformatStreamer(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public ReformatStreamer(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		//Create the processor and parse arguments
		processor=new ReformatProcessor();
		final Parser parser=parse(args);

		//Force single-threading if necessary
		if(processor.uniqueNames || processor.sampleReadsExact || processor.sampleBasesExact){
			Shared.setThreads(1);
		}

		validateParams();
		doPoundReplacement();
		adjustInterleaving();
		fixExtensions();
		checkFileExistence();
		checkStatics();

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffoutsingle=FileFormat.testOutput(outsingle, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

		final boolean samIn=(ffin1!=null && ffin1.samOrBam());
		final boolean samOut=(ffout1!=null && ffout1.samOrBam());
		SamLine.SET_FROM_OK=samIn;
		ReadStreamByteWriter.USE_ATTACHED_SAMLINE=samIn && samOut;
		//Determine if we need to parse SAM fields or can skip for performance
		if(!forceParse && samIn && !samOut){
			SamLine.PARSE_2=false;
			SamLine.PARSE_5=false;
			SamLine.PARSE_6=false;
			SamLine.PARSE_7=false;
			SamLine.PARSE_8=false;
			SamLine.PARSE_OPTIONAL=false;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	private Parser parse(String[] args){

		Parser parser=new Parser();

		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("outsingle") || a.equals("outs")){
				outsingle=b;
			}else if(a.equals("qfoutsingle") || a.equals("qfouts")){
				qfoutsingle=b;
			}else if(a.equals("skipreads")){
				skipreads=Parse.parseKMG(b);
			}else if(a.equals("deleteinput")){
				deleteInput=Parse.parseBoolean(b);
			}else if(a.equals("workers") || a.equals("workerthreads") || a.equals("wt")){
				workers=Integer.parseInt(b);
			}else if(a.equals("threadsin") || a.equals("tin")){
				threadsIn=Integer.parseInt(b);
			}else if(a.equals("threadsout") || a.equals("tout")){
				threadsOut=Integer.parseInt(b);
			}else if(a.equals("forceparse")){
				forceParse=Parse.parseBoolean(b);
			}else if(processor.parse(arg, a, b)){
				//Argument consumed by processor
			}else if(parser.parse(arg, a, b)){
				//Argument consumed by parser
			}else if(i==0 && parser.in1==null && Tools.looksLikeInputSequenceStream(arg)){
				parser.in1=arg;
			}else if(i==1 && parser.in1!=null && parser.out1==null && Tools.looksLikeOutputSequenceStream(arg)){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		//Set processor fields from parser
		processor.setFromParser(parser);

		//Get fields from parser that the harness needs
		maxReads=parser.maxReads;
		breakLength=parser.breakLength;
		overwrite=ReadStats.overwrite=parser.overwrite;
		append=ReadStats.append=parser.append;
		testsize=parser.testsize;
		setInterleaved=parser.setInterleaved;
		if(workers<0){workers=Tools.mid(1, Shared.threads(), 8);}
		else if(workers==0){workers=1;}

		in1=parser.in1;
		in2=parser.in2;
		qfin1=parser.qfin1;
		qfin2=parser.qfin2;
		extin=parser.extin;

		out1=parser.out1;
		out2=parser.out2;
		qfout1=parser.qfout1;
		qfout2=parser.qfout2;
		extout=parser.extout;

		//Finalize processor settings
		processor.postParse();

		return parser;
	}

	private void doPoundReplacement(){
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}

	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		qfin1=Tools.fixExtension(qfin1);
		qfin2=Tools.fixExtension(qfin2);
	}

	private void checkFileExistence(){
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outsingle)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outsingle)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	private void adjustInterleaving(){
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null));
			if(in2!=null){
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}
	}

	private static void checkStatics(){
		//Empty for now
	}

	private boolean validateParams(){
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Primary Method        ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){

		final boolean inputReads=(ffin1!=null && !ffin1.samOrBam());
		final boolean inputSam=(ffin1!=null && ffin1.samOrBam());
		final boolean outputReads=(ffout1!=null && !ffout1.samOrBam());
		final boolean outputSam=(ffout1!=null && ffout1.samOrBam());
		final boolean saveHeader=inputSam && outputSam;

		//Create Streamer and Writers
		Streamer st=StreamerFactory.makeStreamer(ffin1, ffin2, ordered, maxReads,
			saveHeader, true, threadsIn);
		st.setSampleRate(processor.samplerate, processor.sampleseed);
		if(!ffin1.samOrBam()){
			outstream.println("Input is being processed as "+(st.paired() ? "paired" : "unpaired"));
		}

		Writer fw=WriterFactory.makeWriter(ffout1, ffout2, threadsOut, null, saveHeader);
		Writer fwb=WriterFactory.makeWriter(ffoutsingle, null, threadsOut, null, saveHeader);

		//Start streams
		st.start();
		if(fw!=null){fw.start();}
		if(fwb!=null){fwb.start();}

		//Process data
		if(workers>1){
			spawnThreads(st, fw, fwb, inputReads || outputReads);
		}else{
			processSingleThreaded(st, fw, fwb, inputReads || outputReads);
		}
		readsProcessed=processor.readsProcessedT;
		basesProcessed=processor.basesProcessedT;

		//Close writers
		if(fw!=null){
			fw.poisonAndWait();
			readsOut=fw.readsWritten();
			basesOut=fw.basesWritten();
		}
		if(fwb!=null){
			fwb.poisonAndWait();
			readsOut+=fwb.readsWritten();
			basesOut+=fwb.basesWritten();
		}

		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(st, fw, fwb);
		
		//Delete input files if requested
		if(deleteInput && !errorState && out1!=null && in1!=null){
			try {
				new File(in1).delete();
				if(in2!=null){new File(in2).delete();}
			} catch (Exception e) {
				outstream.println("WARNING: Failed to delete input files.");
			}
		}

		t.stop();
		printStats(t);

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Single-Threaded Mode      ----------------*/
	/*--------------------------------------------------------------*/

	private void processSingleThreaded(Streamer st, Writer fw, Writer fwb, boolean readMode){
		if(readMode){
			for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
				processReadList(ln, processor, fw, fwb);
			}
		}else{
			for(ListNum<SamLine> ln=st.nextLines(); ln!=null; ln=st.nextLines()){
				processLineList(ln, processor, fw, fwb);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Thread Management       ----------------*/
	/*--------------------------------------------------------------*/

	private void spawnThreads(Streamer st, Writer fw, Writer fwb, boolean readMode){

		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(workers);
		for(int i=0; i<workers; i++){
			alpt.add(new ProcessThread(st, fw, fwb, processor.clone(), readMode, i));
		}

		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
	}

	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt){
			errorState|=(!pt.success);
			processor.add(pt.processorT);
		}
	}

	@Override
	public final boolean success(){return !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	class ProcessThread extends Thread {

		ProcessThread(Streamer st_, Writer fw_, Writer fwb_,
			ReformatProcessor processorT_, boolean readMode_, int tid_){
			st=st_;
			fw=fw_;
			fwb=fwb_;
			processorT=processorT_;
			readMode=readMode_;
			tid=tid_;
		}

		@Override
		public void run(){
			processInner();
			success=true;
		}

		void processInner(){
			if(readMode){
				for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
					processReadList(ln, processorT, fw, fwb);
				}
			}else{
				for(ListNum<SamLine> ln=st.nextLines(); ln!=null; ln=st.nextLines()){
					processLineList(ln, processorT, fw, fwb);
				}
			}
		}
		
		boolean success=false;

		private final Streamer st;
		private final Writer fw;
		private final Writer fwb;
		final ReformatProcessor processorT;
		final boolean readMode;
		final int tid;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Processing Methods      ----------------*/
	/*--------------------------------------------------------------*/

	private void processReadList(ListNum<Read> ln, ReformatProcessor proc, Writer fw, Writer fwb){
		final ArrayList<Read> reads=ln.list;

		//Handle skipreads
		if(skipreads>0){
			int removed=0;
			for(int i=0; i<reads.size(); i++){
				Read r=reads.get(i);
				if(r.numericID<skipreads){
					reads.set(i, null);
					removed++;
				}else{
					skipreads=-1;
					break;
				}
			}
			if(removed>0){
				Tools.condenseStrict(reads);
			}
		}

		//Handle breakLength
		if(breakLength>0){
			Tools.breakReads(reads, breakLength, proc.minReadLength, (verbose ? outstream : null));
		}

		final ArrayList<Read> singles=(fwb==null ? null : new ArrayList<Read>());

		for(int idx=0; idx<reads.size(); idx++){
			final Read r1=reads.get(idx);
			final Read r2=r1.mate;

			if(!r1.validated()){r1.validate(true);}
			if(r2!=null && !r2.validated()){r2.validate(true);}

			final int keep=proc.processReadPair(r1, r2);

			if(keep==3){
				//Keep pair
			}else if(keep==1){
				//Keep r1 as singleton
				if(singles!=null){
					r1.mate=null;
					singles.add(r1);
				}
				reads.set(idx, null);
			}else if(keep==2){
				//Keep r2 as singleton
				if(singles!=null){
					r2.mate=null;
					singles.add(r2);
				}
				reads.set(idx, null);
			}else{
				//Discard both
				reads.set(idx, null);
			}
		}

		if(fw!=null){fw.addReads(ln);}
		if(fwb!=null && singles!=null && !singles.isEmpty()){
			ListNum<Read> lnb=new ListNum<Read>(singles, ln.id);
			fwb.addReads(lnb);
		}
	}

	private void processLineList(ListNum<SamLine> ln, ReformatProcessor proc, Writer fw, Writer fwb){
		final ArrayList<SamLine> lines=ln.list;

		//Handle skipreads
		if(skipreads>0){
			int removed=0;
			for(int i=0; i<lines.size(); i++){
				SamLine sl=lines.get(i);
				Read r=(Read)sl.obj;
				if(r.numericID<skipreads){
					lines.set(i, null);
					removed++;
				}else{
					skipreads=-1;
					break;
				}
			}
			if(removed>0){
				Tools.condenseStrict(lines);
			}
		}

		for(int i=0; i<lines.size(); i++){
			SamLine sl=lines.get(i);
			
			//Get the attached Read object if it exists
			Read r1=(sl.obj instanceof Read ? (Read)sl.obj : null);
			
			if(r1!=null){
				//Validate the read
				if(!r1.validated()){r1.validate(true);}
				
				//Process using the processor (no mate for SamLines)
				final int keep=proc.processReadPair(r1, null);
				
				//SamLines are never paired, so keep is either 1 (keep r1) or 0 (discard)
				if(keep==0){
					lines.set(i, null);
				}
			}else{
				//No Read object attached, just count it
				proc.readsProcessedT++;
				proc.basesProcessedT+=sl.lengthOrZero();
			}
		}

		if(fw!=null){fw.addLines(ln);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Statistics           ----------------*/
	/*--------------------------------------------------------------*/

	private void printStats(Timer t){

		long readsOut=processor.pairsOut*2+processor.singlesOut;
		long basesOut=processor.pairBasesOut+processor.singleBasesOut;

		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		processor.printStats(outstream);

		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null, in2=null;
	private String qfin1=null, qfin2=null;
	private String out1=null, out2=null, outsingle=null;
	private String qfout1=null, qfout2=null, qfoutsingle=null;//Unsupported currently
	private String extin=null, extout=null;
	private boolean setInterleaved=false;

	private int workers=-1;
	private int threadsIn=-1;
	private int threadsOut=-1;
	private boolean forceParse=false;

	private final ReformatProcessor processor;

	protected long readsProcessed=0;
	protected long basesProcessed=0;
	protected long readsOut=0;
	protected long basesOut=0;

	private long maxReads=-1;
	private long skipreads=-1;
	private int breakLength=0;

	private final FileFormat ffin1, ffin2;
	private final FileFormat ffout1, ffout2, ffoutsingle;

	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean ordered=false;
	private boolean testsize=false;
	private boolean deleteInput=false;

}