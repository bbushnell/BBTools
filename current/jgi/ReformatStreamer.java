package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random; // Added
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import processor.ReformatProcessor;
import shared.MetadataWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream; // Added
import stream.FASTQ;
import stream.Read;
import stream.ReadStreamByteWriter;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ListNum;
import structures.LongList; // Added
import structures.SuperLongList; // Added
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Reformat using Streamer/Writer interfaces with multithreading.
 * @author Brian Bushnell
 * @contributor Isla & Gemini
 * @date November 16, 2025
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

		// Grab fields from processor that are needed for setup
		sampleReadsExact=processor.sampleReadsExact;
		sampleBasesExact=processor.sampleBasesExact;
		prioritizeLength=processor.prioritizeLength;
		sampleReadsTarget=processor.sampleReadsTarget;
		sampleBasesTarget=processor.sampleBasesTarget;
		allowUpsample=processor.allowUpsample;

		//Force single-threading if necessary
		if(processor.uniqueNames || sampleReadsExact || sampleBasesExact){workers=1;}
		ordered=(workers>1);

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
			}else if(a.equals("workers") || a.equals("workerthreads") || a.equals("wt") || a.equals("w")){
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
		workers=(workers<1 ? processor.recommendedWorkers() : workers);

		//Copy sampling fields from processor
		sampleReadsExact=processor.sampleReadsExact;
		sampleBasesExact=processor.sampleBasesExact;
		prioritizeLength=processor.prioritizeLength;
		sampleReadsTarget=processor.sampleReadsTarget;
		sampleBasesTarget=processor.sampleBasesTarget;
		allowUpsample=processor.allowUpsample;

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

		//Handle sampling pre-pass
		if(sampleReadsExact || sampleBasesExact){
			handleSamplingPrepass();
		}

		//Create Streamer and Writers
		Streamer st=StreamerFactory.makeStreamer(ffin1, ffin2, ordered, maxReads,
			saveHeader, true, threadsIn);

		//Set samplerate only if NOT doing exact sampling
		if(!sampleReadsExact && !sampleBasesExact) {
			st.setSampleRate(processor.samplerate, processor.sampleseed);
		}

		if(!ffin1.samOrBam()){
			outstream.println("Input is being processed as "+(st.paired() ? "paired" : "unpaired"));
		}

		Writer fw=WriterFactory.makeWriter(ffout1, ffout2, threadsOut, null, saveHeader);
		Writer fwb=WriterFactory.makeWriter(ffoutsingle, null, threadsOut, null, saveHeader);
//		System.err.println("fw class: "+(fw==null ? "null" : fw.getClass()));
		
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
		if(fw!=null){fw.poisonAndWait();}
		if(fwb!=null){fwb.poisonAndWait();}
		if(verbose){System.err.println("Finished poisonAndWait().");}

		//Get output counts
		if(sampleReadsExact || sampleBasesExact) {
			//Trust the sampler's counts for pairs
			readsOut=sampledReadsOut;
			basesOut=sampledBasesOut;
			if (fwb!=null){//Add un-sampled singles
				readsOut += fwb.readsWritten();
				basesOut += fwb.basesWritten();
			}
		} else {
			//Trust the writers' counts
			if(fw!=null) {
				readsOut=fw.readsWritten();
				basesOut=fw.basesWritten();
			}
			if(fwb!=null) {
				readsOut += fwb.readsWritten();
				basesOut += fwb.basesWritten();
			}
		}

		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(st, fw, fwb);

		//Delete input files if requested
		if(deleteInput && !errorState && out1!=null && in1!=null){
			try{
				new File(in1).delete();
				if(in2!=null){new File(in2).delete();}
			}catch (Exception e){
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
		System.err.println("Spawning "+Tools.plural("worker", workers)+".");
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(workers);
		for(int i=0; i<workers; i++){
			alpt.add(new ProcessThread(st, fw, fwb, processor.clone(), readMode, i));
		}

		boolean success=ThreadWaiter.startAndWait(alpt, this);
		if(verbose){System.err.println("Finished waiting for threads.");}
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

	class ProcessThread extends Thread{

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
			if(verbose){System.err.println("Worker "+tid+" starting.");}
			if(readMode){
				for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
					assert(!ln.poison());
					assert(!ln.last());
					processReadList(ln, processorT, fw, fwb);
//					System.err.print(".");
				}
			}else{
				for(ListNum<SamLine> ln=st.nextLines(); ln!=null; ln=st.nextLines()){
					processLineList(ln, processorT, fw, fwb);
				}
			}
			if(verbose){System.err.println("Worker "+tid+" finished.");}
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

			if(keep==3 || (keep==1 && r2==null)){//Common case
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

		// Handle exact sampling
		if(sampleReadsExact || sampleBasesExact){
			final ArrayList<Read> listOut=new ArrayList<Read>();
			final ArrayList<Read> singlesOut=(fwb==null ? null : singles);
			if(singlesOut!=null){singlesOut.clear();}// Clear list from filtering step

			sampleExact(reads, listOut, singlesOut); // Fills listOut

			// Replace the original list with the new sampled list
			ln.list.clear();
			ln.list.addAll(listOut);

			// Send the modified ListNum (with the sampled list) to output
			if(fw!=null){ fw.addReads(ln);}
			if(fwb!=null && singles!=null && !singles.isEmpty()){
				// Note: Replicating original bug/feature of not sampling singles
				fwb.addReads(new ListNum<Read>(singles, ln.id)); 
			}
		} else {
			// No exact sampling; send the original filtered lists
			if(fw!=null){ fw.addReads(ln); }
			if(fwb!=null && singles!=null && !singles.isEmpty()){
				fwb.addReads(new ListNum<Read>(singles, ln.id));
			}
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
			Read r1=(Read)sl.obj;
			if(r1==null){r1=sl.toRead(false);}

			//Validate the read
			if(!r1.validated()){r1.validate(true);}

			//Process using the processor (no mate for SamLines)
			final boolean keep=proc.processSamLine(sl);

			if(!keep){lines.set(i, null);}
		}

		if(fw!=null){fw.addLines(ln);}
	}

	/*--------------------------------------------------------------*/
	/*----------------  Sampling Helper Methods     ----------------*/
	/*--------------------------------------------------------------*/

	/** Do a pre-pass to calculate exact sampling parameters */
	private void handleSamplingPrepass(){
		if(prioritizeLength){
			SuperLongList sll=makeLengthHist(maxReads);
			LongList list=sll.list();
			long[] array=sll.array();
			int minReadLength=processor.minReadLength;
			if(sampleReadsExact){
				long sum=0;
				for(int i=list.size()-1; i>=0 && sum<sampleReadsTarget; i--){
					long num=list.get(i);
					sum++;
					if(sum>=sampleReadsTarget){
						minReadLength=Tools.max(minReadLength, (int)num);
					}
				}
				for(int i=array.length-1; i>=0 && sum<sampleReadsTarget; i--){
					long count=array[i];
					sum+=count;
					if(sum>=sampleReadsTarget){
						minReadLength=Tools.max(minReadLength, i);
					}
				}
			}else{
				long sum=0;
				for(int i=list.size()-1; i>=0 && sum<sampleBasesTarget; i--){
					long num=list.get(i);
					sum+=num;
					if(sum>=sampleBasesTarget){
						minReadLength=Tools.max(minReadLength, (int)num);
					}
				}
				for(int i=array.length-1; i>=0 && sum<sampleBasesTarget; i--){
					long count=array[i];
					sum+=(count*i);
					if(sum>=sampleBasesTarget){
						minReadLength=Tools.max(minReadLength, i);
					}
				}
			}
			//Set the processor's minReadLength and turn off exact sampling
			processor.minReadLength=minReadLength;
			sampleReadsExact=processor.sampleReadsExact=false;
			sampleBasesExact=processor.sampleBasesExact=false;
		}else{
			long[] counts=countReads(maxReads);
			readsRemaining=counts[0]; //This is number of pairs/singletons
			basesRemaining=counts[2]; //This is total bases
			randy=Shared.threadLocalRandom(processor.sampleseed);
			double prob=(sampleReadsExact ? sampleReadsTarget/(double)(readsRemaining)
				: sampleBasesTarget/(double)(basesRemaining));
			if(!allowUpsample){prob=Tools.min(prob, 1.0);}
			outstream.println("Initial samplerate set to "+String.format("%.6f", prob));
		}
	}

	private void sampleExact(ArrayList<Read> reads, ArrayList<Read> listOut, ArrayList<Read> singlesOut){
		if(sampleReadsExact){
			// Sample pairs
			for(Read r : reads){
				if(r!=null){
					int bases=r.pairLength();
					assert(readsRemaining>0) : readsRemaining;
					double prob=sampleReadsTarget/(double)(readsRemaining);
					while(allowUpsample && prob>1){
						listOut.add(r);
						sampleReadsTarget--;
						sampledReadsOut++;
						sampledBasesOut+=bases;
						prob--;
					}
					if(randy.nextDouble()<prob){
						listOut.add(r);
						sampleReadsTarget--;
						sampledReadsOut++;
						sampledBasesOut+=bases;
					}
				}
				readsRemaining--; // Decrement for every read (pair) looked at
			}

			// Sample singles (replicates bug from original, only samples from main list)
			// This block is intentionally left empty to match original logic.

		}else{// sampleBasesExact
			// Sample pairs
			for(Read r : reads){
				int bases=0;
				if(r!=null){
					bases=r.pairLength();
					assert(basesRemaining>0) : basesRemaining;
					double prob=sampleBasesTarget/(double)(basesRemaining);
					while(allowUpsample && prob>1){
						listOut.add(r);
						sampleBasesTarget-=bases;
						sampledReadsOut++;
						sampledBasesOut+=bases;
						prob--;
					}
					if(randy.nextDouble()<prob){
						listOut.add(r);
						sampleBasesTarget-=bases;
						sampledReadsOut++;
						sampledBasesOut+=bases;
					}
				}
				basesRemaining-=bases; // Decrement by 0 if null
			}
			// As above, no sampling for singles to match original
		}
	}

	/** Copied from ReformatReads.java */
	private long[] countReads(long maxReads){
		if(ffin1.stdio()){
			throw new RuntimeException("Can't precount reads from standard in, only from a file.");
		}

		final ConcurrentReadInputStream cris;
		{
			//Use ConcurrentReadInputStream because it's available and thread-safe
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			if(verbose){outstream.println("Counting Reads");}
			cris.start();
		}

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		long count=0, count2=0, bases=0;

		while(ln!=null && reads!=null && reads.size()>0){
			count+=reads.size();
			for(Read r : reads){
				bases+=r.length();
				count2++;
				if(r.mate!=null){
					bases+=r.mateLength();
					count2++;
				}
			}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		errorState|=ReadWrite.closeStream(cris);
		return new long[]{count, count2, bases};
	}

	/** Copied from ReformatReads.java */
	private SuperLongList makeLengthHist(long maxReads){
		if(ffin1.stdio()){
			throw new RuntimeException("Can't precount reads from standard in, only from a file.");
		}

		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			if(verbose){outstream.println("Counting Reads");}
			cris.start();
		}

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		SuperLongList sll=new SuperLongList(200000);

		while(ln!=null && reads!=null && reads.size()>0){
			for(Read r : reads){
				sll.add(r.length());
				if(r.mate!=null){
					sll.add(r.mateLength());
				}
			}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		errorState|=ReadWrite.closeStream(cris);
		sll.sort();
		return sll;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Statistics           ----------------*/
	/*--------------------------------------------------------------*/

	private void printStats(Timer t){

		//Determine final output counts
		long readsOutPairs=processor.pairsOut*2;
		long basesOutPairs=processor.pairBasesOut;
		long readsOutSingles=processor.singlesOut;
		long basesOutSingles=processor.singleBasesOut;

		long readsOut=readsOutPairs+readsOutSingles;
		long basesOut=basesOutPairs+basesOutSingles;
		
		if(sampleReadsExact || sampleBasesExact) {
			readsOut=sampledReadsOut; // Sampler only tracks pairs
			basesOut=sampledBasesOut;
		}

		//Exact counts and fractions
		outstream.println(Tools.typeReadsBases("Input", readsProcessed, basesProcessed, 30, 25));
		outstream.println(Tools.typeReadsBases("Output", readsProcessed, readsOut, 
			basesProcessed, basesOut, 30, 25));
		outstream.println();

		processor.printStats(outstream);

		if(testsize){
			long bytesProcessed=(new File(in1).length()+(in2==null ? 0 : new File(in2).length())+
				(qfin1==null ? 0 : new File(qfin1).length())+
				(qfin2==null ? 0 : new File(qfin2).length()));//*passes
			double xpnano=bytesProcessed/(double)(t.elapsed);
			String xpstring=(bytesProcessed<100000 ? ""+bytesProcessed : 
				bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");
			while(xpstring.length()<8){xpstring=" "+xpstring;}
			outstream.println("Bytes Processed:    "+xpstring+" \t"+Tools.format("%.2fm bytes/sec", xpnano*1000));
		}
		MetadataWriter.write(null, readsProcessed, basesProcessed, readsOut, basesOut, false);

		//Final output stats
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));

		if(processor.verifypairing){outstream.println("Names appear to be correctly paired.");}
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
	protected long readsOut=0; //This is the final count, from writer or sampler
	protected long basesOut=0; //This is the final count, from writer or sampler

	private long maxReads=-1;
	private long skipreads=-1;
	private int breakLength=0;

	/*--------------------------------------------------------------*/

	//Sampling fields, mirrored from ReformatProcessor for harness logic
	private boolean sampleReadsExact=false;
	private boolean sampleBasesExact=false;
	private boolean allowUpsample=false;
	private boolean prioritizeLength=false;
	private long sampleReadsTarget=0;
	private long sampleBasesTarget=0;
	private long sampledReadsOut=0;
	private long sampledBasesOut=0;

	//Sampling state fields
	private long readsRemaining=0;
	private long basesRemaining=0;
	private Random randy=null;

	/*--------------------------------------------------------------*/

	private final FileFormat ffin1, ffin2;
	private final FileFormat ffout1, ffout2, ffoutsingle;

	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public static boolean verbose2=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean ordered=true;
	private boolean testsize=false;
	private boolean deleteInput=false;

}