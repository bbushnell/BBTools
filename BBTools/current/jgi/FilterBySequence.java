package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashSet;

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
import tracker.ReadStats;

/**
 * Filters by exact sequence matches.
 * Similar to Dedupe.
 * 
 * @author Brian Bushnell
 * @date December 18, 2015
 *
 */
public class FilterBySequence {
	
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
		FilterBySequence x=new FilterBySequence(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public FilterBySequence(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
//		FASTQ.FORCE_INTERLEAVED=false;
//		FASTQ.TEST_INTERLEAVED=false;
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("storebases") || a.equals("keepbases") || a.equals("sb")){
				storeBases=Parse.parseBoolean(b);
			}else if(a.equals("include")){
				include=Parse.parseBoolean(b);
			}else if(a.equals("exclude")){
				include=!Parse.parseBoolean(b);
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("casesensitive") || a.equals("case")){
				toUpperCase=!Parse.parseBoolean(b);
			}else if(a.equals("hdist") || a.equals("subs") || a.equals("s")){
				maxSubs=Integer.parseInt(b);
			}else if(a.equals("mismatchfraction") || a.equals("mf") || a.equals("f") || a.equals("fraction")
					 || a.equals("sf") || a.equals("subfraction")){
				mismatchFraction=Float.parseFloat(b);
			}else if(a.equals("lengthdif") || a.equals("maxlengthdif") || a.equals("lendif")){
				maxLengthDif=Integer.parseInt(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(a.equals("ref")){
				if(b==null){ref=null;}
				else{ref=b.split(",");}
			}else if(a.equals("literal")){
				if(b==null){literal=null;}
				else{literal=b.split(",");}
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
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
		
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
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
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure ref files can be read
		if(!Tools.testInputFiles(true, true, ref)){
			throw new RuntimeException("\nCan't read to some reference files.\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		assert(ref!=null || literal!=null) : "No reference sequences.";
		
		if(ref!=null){
			ffref=new FileFormat[ref.length];
			for(int i=0; i<ref.length; i++){
				ffref[i]=FileFormat.testInput(ref[i], FileFormat.FASTQ, null, true, true);
			}
		}else{
			ffref=null;
		}
		
		refSet=new HashSet<Code>();
		refList=new ArrayList<byte[]>();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		{
			if(ffref!=null){
				for(FileFormat ff : ffref){
					processReference(ff);
				}
			}
			if(literal!=null){
				for(String s : literal){
					addToRef(new Code(s.getBytes()));
				}
			}

			System.err.println("Loaded "+refSet.size()+" unique reference sequence"+(refSet.size()==1 ? "." : "s."));
		}
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros;
		if(ffout1!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			//Notify user of output mode
			if(cris.paired() && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
				outstream.println("Writing interleaved.");
			}
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
			ros.start(); //Start the stream
		}else{ros=null;}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnProcessThreads(cris, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 12));
		outstream.println();
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 12, false));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Load ref sequences */
	private void processReference(final FileFormat ff){
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null, null, null);
			cris.start(); //Start the stream
		}
		spawnLoadThreads(cris);
		ReadWrite.closeStream(cris);
	}
	
	
	/** Spawn process threads */
	private void spawnLoadThreads(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with LoadThreads
		ArrayList<LoadThread> alpt=new ArrayList<LoadThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new LoadThread(cris, i));
		}
		
		//Start the threads
		for(LoadThread pt : alpt){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		for(LoadThread pt : alpt){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			readsLoaded+=pt.readsProcessedT;
			basesLoaded+=pt.basesProcessedT;
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
	}
	
	/** Spawn process threads */
	private void spawnProcessThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		for(ProcessThread pt : alpt){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void addToRef(Collection<Code> codes) {
		synchronized(refSet) {
			for(Code c : codes) {addToRef(c);}
		}
	}
	
	private void addToRef(Code c) {
		boolean novel=refSet.add(c);
		if(novel) {refList.add(c.bases);}
//		System.err.println(new String(c.bases)+": "+novel);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_, final int tid_){
			cris=cris_;
			ros=ros_;
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
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

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
						if(!keep){reads.set(idx, null);}
						else{
							readsOutT+=r1.pairCount();
							basesOutT+=initialLength1+initialLength2;
						}
					}
				}

				//Output reads to the output stream
				if(ros!=null){ros.add(reads, ln.id);}

				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		boolean processReadPair(final Read r1, final Read r2){
			boolean contained=containsBoth(r1, r2);
			return contained==include;
		}
		
		boolean containsBoth(final Read r1, final Read r2){
			boolean match1=contains(r1);
			if(r2==null || !match1) {return match1;}
			boolean match2=contains(r2);
			return match1 && match2;
		}
		
		boolean contains(final Read r){
//			System.err.println(new String(r.bases)+": Query");
			if(r==null || r.bases==null || r.bases.length<1) {return true;}//Empty string is always contained
			final Code c=new Code(r.bases);
//			System.err.println(refSet.contains(c)+", "+include);
			if(refSet.contains(c)){return true;}
//			System.err.println("b");
			if(maxSubs<1 && maxLengthDif<1 && mismatchFraction<=0) {return false;}
//			System.err.println("c");
			return bruteForce(r.bases);
		}
		
		boolean bruteForce(byte[] bases) {
			if(bases==null) {return false;}
			boolean match=compare(bases);
//			System.err.println("g: "+match);
			if(match || !rcomp) {return match;}
			AminoAcid.reverseComplementBasesInPlace(bases);
			match=compare(bases);
//			System.err.println("h: "+match);
			AminoAcid.reverseComplementBasesInPlace(bases);
			return match;
		}
		
		boolean compare(byte[] bases) {
			int minlen=bases.length-maxLengthDif, maxlen=bases.length+maxLengthDif;
//			System.err.println("d: "+minlen+"-"+maxlen);
//			System.err.println(new String(bases)+": Query2");
			boolean match=false;
			for(byte[] ref : refList) {
				if(ref.length>=minlen && ref.length<=maxlen) {
//					System.err.println(new String(ref)+": Ref");
					final byte[] longer, shorter;
					if(bases.length<ref.length) {
						shorter=bases; longer=ref;
					}else {
						shorter=ref; longer=bases;
					}
					match=compare(longer, shorter);
					if(match) {break;}
				}
			}
//			System.err.println("compare returned "+match);
			return match;
		}
		
		boolean compare(byte[] a, byte[] b) {
			final int maxSubs0=Tools.max(maxSubs, 
					(int)Math.round(Tools.min(a.length, b.length)*mismatchFraction));
			final int aStartMax=a.length-b.length;
			for(int aStart=0; aStart<=aStartMax; aStart++) {
				int subs=0;
				for(int apos=aStart, bpos=0; bpos<b.length && subs<=maxSubs0; apos++, bpos++) {
					subs+=(a[apos]==b[bpos] ? 0 : 1);
				}
				if(subs<=maxSubs0) {return true;}
//				System.err.println("Subs="+subs+", maxSubs="+maxSubs);
			}
			return false;
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;

		/** Number of reads output by this thread */
		protected long readsOutT=0;
		/** Number of bases output by this thread */
		protected long basesOutT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Thread ID */
		final int tid;
	}
	
	private class LoadThread extends Thread {
		
		//Constructor
		LoadThread(final ConcurrentReadInputStream cris_, final int tid_){
			cris=cris_;
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
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}
			
			LinkedHashSet<Code> codes=new LinkedHashSet<Code>(4000);
			
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

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

					if(r1!=null){codes.add(new Code(r1.bases));}
					if(r2!=null){codes.add(new Code(r2.bases));}
				}
				
				if(codes.size()>2000){
					synchronized(refSet){
						addToRef(codes);
						codes.clear();
					}
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			if(codes.size()>0){
				synchronized(refSet){
					addToRef(codes);
					codes.clear();
				}
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
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
		/** Thread ID */
		final int tid;
	}
	
	private class Code {
		
		Code(byte[] bases_){
			long fwd=Dedupe.hash(bases_);
			long rev=Dedupe.hashReversed(bases_);
			a=(rcomp ? Tools.max(fwd, rev) : fwd);
			b=(rcomp ? Tools.min(fwd, rev) : rev);
			
			if(storeBases){
				if(a==fwd && !toUpperCase){
					bases=bases_;
				}else{
					bases=bases_.clone();
					if(a!=fwd){AminoAcid.reverseComplementBasesInPlace(bases);}
					for(int i=0; i<bases.length; i++){
						bases[i]=(byte) Tools.toUpperCase(bases[i]);
					}
				}
			}else{
				bases=null;
			}
		}
		
		@Override
		public boolean equals(Object o){
			return equals((Code)o);
		}
		
		public boolean equals(Code c){
			if(a!=c.a || b!=c.b){return false;}
			if(bases==null || c.bases==null){return true;}
			return Tools.equals(bases, c.bases);
		}
		
//		public boolean revEquals(Code c){
//			if(a!=c.b || b!=c.a){return false;}
//			if(bases==null || c.bases==null){return true;}
//			return AminoAcid.equalsReverseComp(bases, c.bases);
//		}
		
		@Override
		public int hashCode(){
			return (int)(a&0x7FFFFFFF);
		}
		
		final long a, b;
		final byte[] bases;
		
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

	private String qfout1=null;
	private String qfout2=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** Ref input file path */
	private String[] ref=null;
	
	/** Literals */
	private String[] literal=null;

	private HashSet<Code> refSet;
	private ArrayList<byte[]> refList;
	
	private boolean storeBases=true;
	
	private boolean include=true;
	
	private boolean rcomp=true;
	
	private boolean toUpperCase=true;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	protected long readsLoaded=0;
	protected long basesLoaded=0;

	/** Number of reads output*/
	protected long readsOut=0;
	/** Number of bases output*/
	protected long basesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;

	private int maxSubs=0;
	private float mismatchFraction=0f;
	private int maxLengthDif=0;
	
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
	
	/** Reference Files */
	private final FileFormat[] ffref;
	
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
