package ml;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
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
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.LongList;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Filters sequences based on matching a neural network.
 * 
 * @author Brian Bushnell
 * @date November 1, 2023
 *
 */
public class NetFilter implements Accumulator<NetFilter.ProcessThread> {
	
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
		NetFilter x=new NetFilter(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public NetFilter(String[] args){
		
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
		}

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffoutu1=FileFormat.testOutput(outu1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffoutu2=FileFormat.testOutput(outu2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		ffnet=FileFormat.testOutput(netFile, FileFormat.BBNET, null, true, true, false, false);
		net0=CellNetParser.load(netFile);
		assert(net0!=null) : netFile;
		if(width<0) {width=(net0.numInputs()-4)/4;}
		else {assert(width==(net0.numInputs()-4)/4) : width+", "+net0.numInputs();}
		if(autoCutoff) {cutoff=net0.cutoff;}
		if(stepsize==Integer.MIN_VALUE) {stepsize=width-overlap;}
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
			}else if(a.equals("net")){
				netFile=b;
			}else if(a.equals("hist")){
				histFile=b;
			}else if(a.equals("width")){
				width=Integer.parseInt(b);
			}else if(a.equals("step") || a.equals("stepsize")){
				stepsize=Integer.parseInt(b);
				overlap=Integer.MIN_VALUE;
			}else if(a.equals("overlap")){
				overlap=Integer.parseInt(b);
				stepsize=Integer.MIN_VALUE;
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("parse")){
				parseHeader=Parse.parseBoolean(b);
			}else if(a.equals("cutoff")){
				if("auto".equalsIgnoreCase(b)) {
					autoCutoff=true;
				}else {
					autoCutoff=false;
					cutoff=Float.parseFloat(b);
				}
				filter=true;
			}else if(a.equals("highpass")){
				highpass=Parse.parseBoolean(b);
				filter=true;
			}else if(a.equals("lowpass")){
				highpass=!Parse.parseBoolean(b);
				filter=true;
			}else if(a.equals("mode") || a.equals("scoremode")){
				SCORE_MODE=parseMode(b);
			}else if(a.equals("pairmode")){
				PAIR_MODE=parseMode(b);
			}else if(a.equals("filter")){
				filter=Parse.parseBoolean(b);
			}else if(a.equals("outu") || a.equals("outu1")){
				outu1=b;
			}else if(a.equals("outu2")){
				outu2=b;
			}else if(a.equals("annotate") || a.equals("rename")){
				annotate=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
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

		//Do output file # replacement
		if(outu1!=null && outu2==null && outu1.indexOf('#')>-1){
			outu2=out1.replace("#", "2");
			outu1=out1.replace("#", "1");
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
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
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
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
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
		assert(stepsize>0) : "stepsize must be greater than 0: s="+stepsize+", o="+overlap+", w="+width;
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
		final ConcurrentReadInputStream cris=makeCris();
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros(cris.paired(), ffout1, ffout2, qfout1, qfout2);
		
		//Unmatched read output stream
		final ConcurrentReadOutputStream rosu=makeCros(cris.paired(), ffoutu1, ffoutu2, null, null);
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads(cris, ros, rosu);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros, rosu);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		float mult=100f/readsProcessed;
		outstream.println();
		outstream.println(String.format("Average Score:    \t%.4f", scoreSum/scoreCount));
		if(parseHeader) {
			outstream.println(String.format("Average Positive: \t%.4f\t\t%.4f%%", scoreSumPositive/scoreCountPositive, scoreCountPositive*mult));
			outstream.println(String.format("Average Negative: \t%.4f\t\t%.4f%%", scoreSumNegative/scoreCountNegative, scoreCountNegative*mult));
		}
		if(filter || true) {
			outstream.println(String.format("Average Pass:     \t%.4f\t\t%.4f%%", scoreSumPass/scoreCountPass, scoreCountPass*mult));
			outstream.println(String.format("Average Fail:     \t%.4f\t\t%.4f%%", scoreSumFail/scoreCountFail, scoreCountFail*mult));
		}
		if(parseHeader) {
			outstream.println(String.format("True Positive:    \t%d\t\t%.4f%%", tpCount, tpCount*mult));
			outstream.println(String.format("True Negative:    \t%d\t\t%.4f%%", tnCount, tnCount*mult));
			outstream.println(String.format("False Positive:   \t%d\t\t%.4f%%", fpCount, fpCount*mult));
			outstream.println(String.format("False Negative:   \t%d\t\t%.4f%%", fnCount, fnCount*mult));
		}
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private ConcurrentReadInputStream makeCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}
	
	private ConcurrentReadOutputStream makeCros(boolean pairedInput, FileFormat ff1, FileFormat ff2, String qf1, String qf2){
		if(ff1==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

//		//Notify user of output mode
//		if(pairedInput && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
//			outstream.println("Writing interleaved.");
//		}

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff1, ff2, qf1, qf2, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, 
			final ConcurrentReadOutputStream rosu){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, rosu, i));
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
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			
			scoreSum+=pt.scoreSumT;
			scoreSumPositive+=pt.scoreSumPositiveT;
			scoreSumNegative+=pt.scoreSumNegativeT;
			scoreSumPass+=pt.scoreSumPassT;
			scoreSumFail+=pt.scoreSumFailT;

			scoreCount+=pt.scoreCountT;
			scoreCountPositive+=pt.scoreCountPositiveT;
			scoreCountNegative+=pt.scoreCountNegativeT;
			scoreCountPass+=pt.scoreCountPassT;
			scoreCountFail+=pt.scoreCountFailT;
			
			fpCount+=pt.fpCountT;
			fnCount+=pt.fnCountT;
			tpCount+=pt.tpCountT;
			tnCount+=pt.tnCountT;
			
			phist.incrementBy(pt.phistT);
			mhist.incrementBy(pt.mhistT);
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
	
//	/** This class is static to prevent accidental writing to shared variables.
//	 * It is safe to remove the static modifier. */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_, 
				final ConcurrentReadOutputStream rosu_, final int tid_){
			cris=cris_;
			ros=ros_;
			rosu=rosu_;
			tid=tid_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			net=net0.copy(false);
			vec=new float[net.numInputs()];
			
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
		
		void processList(ListNum<Read> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;
			final ArrayList<Read> ureads=(rosu==null ? null : new ArrayList<Read>(reads.size()));
			
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

				//Increment counters
				readsProcessedT+=r1.pairCount();
				basesProcessedT+=initialLength1+initialLength2;
				
				{
					//Reads are processed in this block.
					boolean keep=processReadPair(r1, r2);
					
					if(keep){
						readsOutT+=r1.pairCount();
						basesOutT+=r1.pairLength();
					}else{
						reads.set(idx, null);
						if(ureads!=null) {ureads.add(r1);}
					}
				}
			}

			//Output reads to the output stream
			if(ros!=null){ros.add(reads, ln.id);}
			if(rosu!=null){rosu.add(ureads, ln.id);}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		boolean processReadPair(final Read r1, final Read r2){
			float score1=processRead(r1, vec);
			float score2=(r2==null ? score1 : processRead(r2, vec));
			float score=(r2==null || PAIR_MODE==SINGLE ? score1 : 
				PAIR_MODE==AVERAGE ? (score1+score2)*0.5f : 
				PAIR_MODE==MAX ? Tools.max(score1, score2) : 
				Tools.min(score1, score2));
			boolean pass=(!filter ? true : (score>=cutoff)==highpass);
			return pass;
		}
		
		private float processRead(Read r, float[] vec) {
			if(r==null) {return -1;}
			float result=(parseHeader ? Parse.parseFloat(r.id, "result=", '\t') : 0);
			float score=score(r.bases, vec, k, net, rcomp);

			boolean pass=(!filter ? true : (score>=cutoff)==highpass);
			boolean positive=(result>=0.5f);
			boolean correct=(positive==(score>=cutoff));
			scoreSumT+=score;
			scoreCountT++;
			
			if(positive) {
				scoreSumPositiveT+=score;
				scoreCountPositiveT++;
				tpCountT+=(correct ? 1 : 0);
				fnCountT+=(correct ? 0 : 1);
			}else {
				scoreSumNegativeT+=score;
				scoreCountNegativeT++;
				tnCountT+=(correct ? 1 : 0);
				fpCountT+=(correct ? 0 : 1);
			}
			
			if(pass) {
				scoreSumPassT+=score;
				scoreCountPassT++;
			}else {
				scoreSumFailT+=score;
				scoreCountFailT++;
			}
			
			int iscore=Tools.mid(0, Math.round(score*100), 100);
			if(result<0.5f) {mhistT.increment(iscore);}
			else {phistT.increment(iscore);}
			if(annotate) {r.id+=String.format("\tscore=%.4f", score);}
			return score;
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;

		double scoreSumT=0;
		double scoreSumPositiveT=0;
		double scoreSumNegativeT=0;
		double scoreSumPassT=0;
		double scoreSumFailT=0;
		long scoreCountT=0;
		long scoreCountPositiveT=0;
		long scoreCountNegativeT=0;
		long scoreCountPassT=0;
		long scoreCountFailT=0;
		
		long fpCountT=0;
		long fnCountT=0;
		long tpCountT=0;
		long tnCountT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;

		private LongList phistT=new LongList(101);
		private LongList mhistT=new LongList(101);
		
		private CellNet net;
		private float[] vec;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Shared unmatched output stream */
		private final ConcurrentReadOutputStream rosu;
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	
	float score(byte[] bases, float[] vec, int k, CellNet net, boolean rcomp) {
		final float f=score(bases, vec, k, net), r;
		if(rcomp) {
			AminoAcid.reverseComplementBasesInPlace(bases);
			r=score(bases, vec, k, net);
			AminoAcid.reverseComplementBasesInPlace(bases);
		}else {r=f;}
		return Tools.max(r, f);
	}
	
	float score(byte[] bases, float[] vec, int k, CellNet net) {
		if(SCORE_MODE==SINGLE) {
			return scoreSingle(bases, vec, k, net);
		}else {
			return scoreFrames(bases, vec, k, net, SCORE_MODE, stepsize, width);
		}
	}
	
	public static float scoreSingle(byte[] bases, float[] vec, int k, CellNet net) {
		assert(vec!=null);
		SequenceToVector.fillVector(bases, vec, k);
		net.applyInput(vec);
		final float f=net.feedForward();
		return f;
	}
	
	public static float scoreFrames(byte[] bases, float[] vec, int k, CellNet net, int SCORE_MODE, int stepsize, int width) {
		assert(vec!=null);
		assert(SCORE_MODE>=0 && SCORE_MODE!=SINGLE) : SCORE_MODE;
		float min=9999, max=-9999;
		double sum=0;
		long frames=1;
		int start=0, stop=width;
		do {
			SequenceToVector.fillVector(bases, vec, k, start, stop);
			net.applyInput(vec);
			final float f=net.feedForward();
			frames++;
			sum+=f;
			min=Tools.min(min, f);
			max=Tools.max(max, f);
			start+=stepsize;
			stop+=stepsize;
		}while(stop<bases.length-1);
		return SCORE_MODE==AVERAGE ? (float)(sum/frames) : SCORE_MODE==MIN ? min : max;
	}
	
	public static float score(byte[] bases, float[] vec, int k, CellNet net, int from, int to) {
		SequenceToVector.fillVector(bases, vec, k, from, to);
		net.applyInput(vec);
		final float f=net.feedForward();
		return f;
	}
	
	private static int parseMode(String b) {
		int x=Tools.find(b.toLowerCase(), MODES);
		if(x<0) {throw new RuntimeException("Unknown mode "+b+"; must be one of "+Arrays.toString(MODES));}
		return x;
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

	private String outu1=null;
	private String outu2=null;

	private String qfout1=null;
	private String qfout2=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;
	
	private String netFile=null;
	private String histFile=null;
	private final FileFormat ffnet;

	private LongList phist=new LongList(101);
	private LongList mhist=new LongList(101);
	
	/*--------------------------------------------------------------*/
	

	private boolean parseHeader=false;
	private final CellNet net0;
	
	/** For raw sequence classification */
	private boolean rcomp=false;
	private int width=-1;
	
	/** For kmer spectrum classification */
	private int k=0;
	private int dims; //asdf//TODO
	
	private boolean filter=true;
	private boolean highpass=true;
	private boolean annotate=false;
	private float cutoff=0.5f;
	private boolean autoCutoff=true;
	private int stepsize=1;
	private int overlap=Integer.MIN_VALUE;

	private int SCORE_MODE=SINGLE;
	private int PAIR_MODE=AVERAGE;
	private static final int SINGLE=0, AVERAGE=1, MAX=2, MIN=3;
	private static final String[] MODES= {"single", "average", "max", "min"};
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;
	
	double scoreSum=0;
	double scoreSumPositive=0;
	double scoreSumNegative=0;
	double scoreSumPass=0;
	double scoreSumFail=0;
	long scoreCount=0;
	long scoreCountPositive=0;
	long scoreCountNegative=0;
	long scoreCountPass=0;
	long scoreCountFail=0;
	long fpCount=0;
	long fnCount=0;
	long tpCount=0;
	long tnCount=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
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
	
	/** Unmatched output file */
	private final FileFormat ffoutu1;
	/** Secondary unmatched output file */
	private final FileFormat ffoutu2;
	
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
	private boolean ordered=true;
	
}
