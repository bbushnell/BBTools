package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
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
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.LongList;
import structures.SuperLongList;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * This class is designed to refactor ReformatReads by
 * using a ReformatProcessor to handle per-read logic
 * while this class manages I/O and multithreading.
 * @author Brian Bushnell
 * @contributor Gemini
 * @date November 16, 2025
 */
public class Reformat2 implements Accumulator<Reformat2.ProcessThread>{

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
		Reformat2 x=new Reformat2(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Reformat2(String[] args){

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
		final Parser parser=parse(args); // This calls processor.parse()

		// Grab fields from processor that are needed for setup
		sampleReadsExact=processor.sampleReadsExact;
		sampleBasesExact=processor.sampleBasesExact;
		prioritizeLength=processor.prioritizeLength;
		sampleReadsTarget=processor.sampleReadsTarget;
		sampleBasesTarget=processor.sampleBasesTarget;
		allowUpsample=processor.allowUpsample;

		//Force single-threading if necessary
		if(processor.uniqueNames || sampleReadsExact || sampleBasesExact){workers=1;}
		ordered=ordered && (workers>1);

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program	

		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffoutsingle=FileFormat.testOutput(outsingle, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Parse arguments from the command line */
	private Parser parse(String[] args){

		//Create a parser object
		Parser parser=new Parser();

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			//Pass arguments to the processor first
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
			}else if(processor.parse(arg, a, b)){
				//Argument was consumed by the processor
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(i==0 && parser.in1==null && Tools.looksLikeInputSequenceStream(arg)){
				parser.in1=arg;
			}else if(i==1 && parser.in1!=null && parser.out1==null && Tools.looksLikeOutputSequenceStream(arg)){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		//Set processor fields from the parser
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

		processor.postParse();
		workers=(workers<1 ? processor.recommendedWorkers() : workers);

		return parser;
	}

	/** Replace # with 1 and 2 in headers */
	private void doPoundReplacement(){
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}

		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}

	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		qfin1=Tools.fixExtension(qfin1);
		qfin2=Tools.fixExtension(qfin2);
	}

	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outsingle)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");	
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outsingle)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/** Make sure interleaving agrees with number of input and output files */
	private void adjustInterleaving(){
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){//If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{//There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
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
		// assert(false) : "TODO"; // Add any harness-level validation
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Primary Method        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){

		//Handle sampling pre-pass
		if(sampleReadsExact || sampleBasesExact){
			handleSamplingPrepass();
		}

		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris(processor.samplerate, processor.sampleseed);

		//Create read output streams
		final ConcurrentReadOutputStream ros=makeCros(cris.paired());
		final ConcurrentReadOutputStream rosb=makeRosb(cris.paired());

		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;

		//Process the reads in separate threads
		spawnThreads(cris, ros, rosb);

		if(verbose){outstream.println("Finished; closing streams.");}

		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros, rosb);

		//Delete input files if requested
		if(deleteInput && !errorState && out1!=null && in1!=null){
			try{
				new File(in1).delete();
				if(in2!=null){new File(in2).delete();}
			}catch (Exception e){
				outstream.println("WARNING: Failed to delete input files.");
			}
		}

		//Report timing and results
		t.stop();

		//Print stats
		printStats(t);

		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

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

	private ConcurrentReadInputStream makeCris(float samplerate, long sampleseed){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);

		//Set samplerate only if NOT doing exact sampling
		if(!sampleReadsExact && !sampleBasesExact){
			cris.setSampleRate(samplerate, sampleseed);
		}

		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}

	private ConcurrentReadOutputStream makeCros(boolean pairedInput){
		if(ffout1==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		//Notify user of output mode
		if(pairedInput && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
			outstream.println("Writing interleaved.");
		}

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}

	private ConcurrentReadOutputStream makeRosb(boolean pairedInput){
		if(ffoutsingle==null){return null;}
		final int buff=8;
		final ConcurrentReadOutputStream rosb=ConcurrentReadOutputStream.getStream(ffoutsingle, null, qfoutsingle, null, buff, null, false);
		rosb.start(); //Start the stream
		return rosb;
	}

	private void printStats(Timer t){

		long readsOut=processor.pairsOut*2+processor.singlesOut;
		long basesOut=processor.pairBasesOut+processor.singleBasesOut;
		if(sampleReadsExact || sampleBasesExact) {
			readsOut=sampledReadsOut;
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
		//outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false)); //This is redundant with the Tools.typeReadsBases call

		if(processor.verifypairing){outstream.println("Names appear to be correctly paired.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Thread Management       ----------------*/
	/*--------------------------------------------------------------*/

	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosb){
		System.err.println("Spawning "+Tools.plural("worker", workers)+".");
		//Do anything necessary prior to processing

		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(workers);
		for(int i=0; i<workers; i++){
			alpt.add(new ProcessThread(cris, ros, rosb, processor.clone(), i));
		}

		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;

		//Do anything necessary after processing

	}

	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt){
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			//readsOut+=pt.readsOutT; //No longer used, processor handles this
			//basesOut+=pt.basesOutT; //No longer used, processor handles this
			errorState|=(!pt.success);

			//Add processor stats
			processor.add(pt.processorT);
		}
	}

	@Override
	public final boolean success(){return !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	class ProcessThread extends Thread{

		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_, final ConcurrentReadOutputStream rosb_,
			final ReformatProcessor processorT_, final int tid_){
			cris=cris_;
			ros=ros_;
			rosb=rosb_;
			processorT=processorT_;
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

			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired());
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){

				processList(ln);

				//Notify the input stream that the list was used
				cris.returnList(ln);

				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		void processList(ListNum<Read> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;

			//Handle per-list operations from original ReformatReads
			if(skipreads>0){
				int removed=0;
				for(int i=0; i<reads.size(); i++){
					Read r=reads.get(i);
					if(r.numericID<skipreads){
						reads.set(i, null);
						removed++;
					}else{
						skipreads=-1; //Stop skipping
						break;
					}
				}
				if(removed>0){
					Tools.condenseStrict(reads);
				}
			}

			if(breakLength>0){
				Tools.breakReads(reads, breakLength, processorT.minReadLength, (verbose ? outstream : null));
			}

			//Create a list for reads that become singletons
			final ArrayList<Read> singles=(rosb==null ? null : new ArrayList<Read>());

			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				final Read r2=r1.mate;

				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}
				if(r2!=null && !r2.validated()){r2.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r1.length();
				final int initialLength2=r1.mateLength();

				//Increment harness counters (reads IN)
				readsProcessedT+=r1.pairCount();
				basesProcessedT+=initialLength1+initialLength2;

				{
					//Reads are processed by the processor.
					//The processor is responsible for ALL filtering and stats.
					final int keep=processorT.processReadPair(r1, r2);

					if(keep==3 || (keep==1 && r2==null)){//Common case
						//Keep pair; do nothing, r1 remains in list
					}else if(keep==1){
						assert(r2!=null);
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
						//Keep==0; discard both
						reads.set(idx, null);
					}
				}
			}

			// Handle exact sampling
			if(sampleReadsExact || sampleBasesExact){
				// Re-use 'singles' list for exact-sampled singletons
				final ArrayList<Read> listOut=new ArrayList<Read>();
				final ArrayList<Read> singlesOut=(rosb == null ? null : singles);
				if(singlesOut!=null){singlesOut.clear();}// Clear list, as it was populated above
				sampleExact(reads, listOut, singlesOut);

				// Send the *new, sampled* lists to output
				if(ros!=null){ros.add(listOut, ln.id);}// Send the sampled list
				// Send the UNsampled singles
				if(rosb!=null && singles!=null && !singles.isEmpty()){rosb.add(singles, ln.id);}

			}else{
				// No exact sampling; send the original lists
				if(ros!=null){ros.add(reads, ln.id);}
				if(rosb!=null && singles!=null && !singles.isEmpty()){rosb.add(singles, ln.id);}
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;

		/** True only if this thread has completed successfully */
		boolean success=false;

		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Shared output stream for singletons */
		private final ConcurrentReadOutputStream rosb;
		/** Per-thread processor instance */
		final ReformatProcessor processorT;
		/** Thread ID */
		final int tid;
	}

	/*--------------------------------------------------------------*/
	/*----------------  Sampling Helper Methods     ----------------*/
	/*--------------------------------------------------------------*/

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
			// If you want to fix the bug, you would iterate over 'singles' here
			// and apply the same variable-rate logic.

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
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			if(verbose){outstream.println("Counting Reads");}
			cris.start(); //4567
		}

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		long count=0, count2=0, bases=0;

		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
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
			cris.start(); //4567
		}

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		SuperLongList sll=new SuperLongList(200000);

		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
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
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;

	private String qfin1=null;
	private String qfin2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;
	/** Output file path for singletons */
	private String outsingle=null;

	private String qfout1=null;
	private String qfout2=null;
	private String qfoutsingle=null;

	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;

	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;

	private int workers=-1;

	/*--------------------------------------------------------------*/

	/** The main processor for accumulating stats */
	private final ReformatProcessor processor;

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0; // No longer used
	/** Number of bases retained */
	protected long basesOut=0; // No longer used

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	/** Skip this many input reads */
	private long skipreads=-1;
	/** Break reads into chunks of this length */
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
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;

	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	/** Singleton output file */
	private final FileFormat ffoutsingle;

	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/*--------------------------------------------------------------*/
	/*----------------         Common Fields        ----------------*/
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
	private boolean ordered=true;
	private boolean testsize=false;
	private boolean deleteInput=false;

}