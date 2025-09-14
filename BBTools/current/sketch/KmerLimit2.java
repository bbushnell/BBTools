package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

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
import structures.IntMap;
import structures.ListNum;
import tracker.ReadStats;

/**
 * 
 * @author Brian Bushnell
 * @date July 30, 2018
 *
 */
public class KmerLimit2 extends SketchObject {
	
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
		KmerLimit2 x=new KmerLimit2(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerLimit2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		SketchObject.setKeyFraction(0.1);
		defaultParams.minEntropy=0;
		defaultParams.minProb=0.2f;
		
		boolean setHeapSize=false;
		int heapSize_=8091;
		long targetKmers_=0;
		int k_=32;
		int minCount_=1;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.overwrite=true;
		
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
			}else if(a.equals("size") || a.equals("heapsize")){
				heapSize_=Parse.parseIntKMG(b);
				setHeapSize=true;
			}else if(a.equals("kmers") || a.equals("target") || a.equals("limit")){
				targetKmers_=Parse.parseKMG(b);
			}else if(a.equals("mincount")){
				minCount_=Parse.parseIntKMG(b);
			}else if(a.equals("maxexpandedlength") || a.equals("maxlength") || a.equals("maxlen")){
				maxExpandedLength=Parse.parseIntKMG(b);
			}else if(a.equals("seed")){
				seed=Parse.parseKMG(b);
			}else if(a.equals("trials")){
				trials=Parse.parseIntKMG(b);
			}else if(parseSketchFlags(arg, a, b)){
				parser.parse(arg, a, b);
			}else if(defaultParams.parse(arg, a, b)){
				parser.parse(arg, a, b);
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
		
		if(!setHeapSize && minCount_>1){heapSize_=32000;}
		heapSize=heapSize_;
		targetKmers=targetKmers_;
		k=k_;
		minCount=minCount_;
		assert(targetKmers>0) : "Must set a kmer limit.";
		assert(heapSize>0) : "Heap size must be positive.";
		assert(k>0 && k<=32) : "0<k<33; k="+k;
		postParse();
		
//		if(minCount>1){
//			Shared.setBufferLen(800);
//		}
		
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

		minProb=defaultParams.minProb;
		minQual=defaultParams.minQual;
		
		shift=2*k;
		shift2=shift-2;
		mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		sharedHeap=new SketchHeap(heapSize, 0, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
//		//Optionally create a read output stream
//		final ConcurrentReadOutputStream ros;
//		if(ffout1!=null){
//			//Select output buffer size based on whether it needs to be ordered
//			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
//			
//			//Notify user of output mode
//			if(cris.paired() && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
//				outstream.println("Writing interleaved.");
//			}
//			
//			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
//			ros.start(); //Start the stream
//		}else{ros=null;}
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads0();
		
//		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;

		
		Sketch sketch=new Sketch(sharedHeap, true, true, null);
		sketch=capLengthAtCountSum(sketch, maxExpandedLength);
		final long reads=Tools.max(1, sketch.genomeSequences);
		final long targetReads=calcTargetReads(sketch, targetKmers, minCount, trials, seed);
		final double targetRate=Tools.min(1, targetReads/(double)reads);
		final String targetRateS=Tools.format("%.4f%%",targetRate*100);
		
		//Report timing and results
		t.stop();
		outstream.println("Finished counting kmers.");
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		String kstring0=Tools.padKMB(sketch.genomeSizeEstimate(minCount), 8);
		String rstring0=Tools.padKMB(targetReads, 8);
		outstream.println("Unique Kmers:       "+kstring0);
		outstream.println("Target Reads:       "+rstring0+"\t"+targetRateS);
		
//		outstream.println("Reads:        \t"+reads);
//		outstream.println("Unique Kmers: \t"+sketch.genomeSizeEstimate(minCount));
//		outstream.println("Target Reads: \t"+targetReads);
//		outstream.println("Sample Rate:  \t"+targetRateS);
//		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		t.start();
		outstream.println("\nSubsampling reads.");
		
//		String kstring=Tools.padKM(sharedHeap.genomeSizeEstimate(minCount), 8);
//		outstream.println("Unique Kmers Out:   "+kstring);
		
		
//		ArrayList<String> args=new ArrayList<String>();
//		args.add("in="+in1);
//		if(in2!=null){args.add("in2="+in2);}
//		args.add("out="+out1);
//		if(out2!=null){args.add("out2="+out2);}
//		args.add("ordered="+ordered);
//		args.add("ow="+(overwrite ? "t" : "f"));
//		if(targetRate<1){args.add("samplerate="+targetRateS);}
//		args.add("loglogout");
//		args.add("loglogk="+k);
//		args.add("loglogminprob="+minProb);
//		BBDukF.main(args.toArray(new String[0]));
		
//		Sketch sk=new Sketch(sharedHeap, true, true, null);
//		outstream.println(sk.genomeSizeEstimate());
		spawnThreads2(targetRate);
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		String kstring=Tools.padKMB(sharedHeap.genomeSizeEstimate(minCount), 8);
		outstream.println("Unique Kmers Out:   "+kstring);
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Spawn process threads */
	private void spawnThreads0(){
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Determine how many threads may be used
		final int threads=Tools.min(10, Shared.threads());
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, null, i, heapSize));
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
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
	}
	
	/** Spawn process threads */
	private void spawnThreads2(double rate){
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.setSampleRate((float)rate, seed);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
//		paired=cris.paired();
//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
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
		
		//Determine how many threads may be used
		final int threads=Tools.min(10, Shared.threads());
		
		sharedHeap.clear();			
//		readsProcessed=0;
//		basesProcessed=0;
		readsOut=0;
		basesOut=0;
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, i, heapSize));
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
//			readsProcessed+=pt.readsProcessedT;
//			basesProcessed+=pt.basesProcessedT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros);
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public static Sketch capLengthAtCountSum(Sketch sketch0, int max) {
		int len=0;
		long sum=0;
		for(; len<sketch0.keyCounts.length; len++){
			sum=sum+sketch0.keyCounts[len];
			if(sum>max){break;}
		}
		if(len>=sketch0.length()){return sketch0;}
		
		long[] keys=Arrays.copyOf(sketch0.keys, len);
		int[] counts=Arrays.copyOf(sketch0.keyCounts, len);
		
//		long[] array_, int[] counts_, int taxID_, long imgID_, long gSizeBases_, long gSizeKmers_, long gSequences_, double probCorrect_,
//		String taxName_, String name0_, String fname_, ArrayList<String> meta_
		
		Sketch sk=new Sketch(keys, counts, null, null, null, -1, -1, 
				sketch0.genomeSizeBases, sketch0.genomeSizeKmers, sketch0.genomeSequences, sketch0.probCorrect,
				null, null, null, null);
		
		return sk;
	}
	
	public static long calcTargetReads(Sketch sketch, long targetKmers, int minCount, int trials, long seed){
		final int[] counts0=sketch.keyCounts;
		final int[] counts=Arrays.copyOf(counts0, counts0.length);
		final long size=sketch.genomeSizeEstimate(minCount);
		final long reads=sketch.genomeSequences;
		final double targetKmerFraction=targetKmers/(double)size;
		if(targetKmerFraction>=1){return reads;}
		
		final int targetKeys=(int)(targetKmerFraction*counts.length);
		final long countSum=shared.Vector.sum(counts0);
		assert(countSum<Shared.MAX_ARRAY_LEN) : countSum;
//		System.err.println("countsum: "+countSum);
		
		final IntMap map=new IntMap(0, counts0.length);
		final int[] expanded=new int[(int)countSum];
		
		long roundSum=0;
		final Random randy=Shared.threadLocalRandom(seed);
		for(int i=0; i<trials; i++){
			Tools.fill(counts, counts0);
//			long rounds=reduceRounds(counts0, counts, minCount, targetKeys, randy);
			long rounds=reduceRoundsIM(counts0, expanded, minCount, targetKeys, randy, map);
			roundSum+=rounds;
		}
		double avgRounds=roundSum/(double)trials;
//		System.err.println("avgRounds: "+avgRounds);
		double targetCountFraction=1-(avgRounds/countSum);
//		System.err.println("targetFraction: "+targetCountFraction);
		return (long)(targetCountFraction*reads);
	}
	
//	public static int reduceRoundsOld(final int[] counts, final int minCount, final int targetKeys, final Random randy){
//		assert(minCount>=0) : minCount;
//		int rounds=0;
//		int valid=0;
//		for(int x : counts){
//			if(x>=minCount){valid++;}
//		}
//		
//		int len=counts.length;
//		System.err.println(targetKeys+", "+counts.length+", "+valid+", "+len+", "+rounds+", "+shared.Vector.sum(counts)+", "+Arrays.toString(counts));
//		for(; valid>targetKeys; rounds++){
//			int pos=randy.nextInt(len);
////			assert(counts[pos]>0) : pos+"/"+len+": "+targetKeys+", "+counts.length+", "+valid+", "+len+", "+rounds+", "+Arrays.toString(counts);
//			if(counts[pos]==minCount){valid--;}
//			counts[pos]--;
//			if(counts[pos]==0){
//				len--;//shrink the array
//				System.err.println("len="+len+", counts[len]="+counts[len]);
//				System.err.println("pos="+pos+", counts[pos]="+counts[pos]);
//				counts[pos]=counts[len];//move the last element to the empty slot 
//				counts[len]=0;
//				if(pos!=len && len>0){
//					assert(counts[pos]>0) : pos+"/"+len+": "+targetKeys+", "+counts.length+", "+valid+", "+len+", "+rounds+", "+Arrays.toString(counts);
//				}
//			}
//			System.err.println(len+", "+pos+": "+Arrays.toString(counts));
//		}
//		
//		System.err.println(targetKeys+", "+counts.length+", "+valid+", "+len+", "+rounds+", "+shared.Vector.sum(counts));
//		
//		return rounds;
//	}
	
	//This can be done faster with bins.
	//Each bin contains all kmers with count x.  When a bin is hit, one kmer moves to the next bin lower.
	//Alternately, expand the array into one physical kmer per count.  Store the current counts in an IntMap. Remove key each time.
	public static long reduceRounds(final int[] counts0, final int[] counts, final int minCount, final int targetKeys, final Random randy){
		assert(minCount>=0) : minCount;
		long rounds=0;
		int valid=0;
		for(int x : counts){
			if(x>=minCount){valid++;}
		}
		
		int len=counts.length;
		final long sum0=shared.Vector.sum(counts);
		long sum=sum0;
//		System.err.println(targetKeys+", "+counts.length+", "+valid+", "+len+", "+rounds+", "+shared.Vector.sum(counts)+", "+Arrays.toString(counts));
		for(; valid>targetKeys; rounds++){
			long posNum=(Long.MAX_VALUE&randy.nextLong())%sum;
			long sum2=0;
			int pos=0;
			
			for(int i=0; i<counts.length; i++){
				int x=counts[i];
				if(x>0){
					sum2+=x;
					if(sum2>=posNum){
						pos=i;
						break;
					}
				}
			}
			
//			for(int i=0; i<counts0.length; i++){
//				int x=counts0[i];
//				if(x>0){
//					sum2+=x;
//					if(sum2>=posNum){
//						pos=i;
//						break;
//					}
//				}
//			}
			
			sum--;
			
			assert(counts[pos]>0) : pos+"/"+len+": "+targetKeys+", "+counts.length+", "+valid+", "+len+", "+rounds+", "+Arrays.toString(counts);
			if(counts[pos]==minCount){valid--;}
			counts[pos]--;
			if(counts[pos]==0){
				len--;//shrink the array
			}
//			System.err.println(len+", "+pos+": "+Arrays.toString(counts));
		}
		
//		System.err.println(targetKeys+", "+counts.length+", "+valid+", "+len+", "+rounds+", "+shared.Vector.sum(counts));
		
		return rounds;
	}
	
	//This can be done faster with bins.
	//Each bin contains all kmers with count x.  When a bin is hit, one kmer moves to the next bin lower.
	//Alternately, expand the array into one physical kmer per count.  Store the current counts in an IntMap. Remove key each time.
	public static long reduceRoundsIM(final int[] counts0, final int[] expanded, final int minCount, final int targetKeys, final Random randy, final IntMap map){
		assert(minCount>=0) : minCount;
		long rounds=0;
		int valid=0;
		map.clear();
		for(int i=0, k=0; i<counts0.length; i++){
			int x=counts0[i];
//			counts[i]=counts0[i];
			if(x>=minCount){valid++;}
			map.put(i, x);
			for(int j=0; j<x; j++, k++){
				expanded[k]=i;
			}
		}
		assert(expanded.length==shared.Vector.sum(counts0));
		
		int len=expanded.length;
//		System.err.println(targetKeys+", "+counts.length+", "+valid+", "+len+", "+rounds+", "+shared.Vector.sum(counts)+", "+Arrays.toString(counts));
		for(; valid>targetKeys; rounds++){
			final int pos=randy.nextInt(len);
			final int key=expanded[pos];
			final int x=map.get(key);
			assert(x>0);
			
			
			if(x==minCount){valid--;}
			map.put(key, x-1);
			
			len--;//shrink the array
			//				System.err.println("len="+len+", counts[len]="+counts[len]);
			//				System.err.println("pos="+pos+", counts[pos]="+counts[pos]);
			expanded[pos]=expanded[len];//move the last element to the empty slot 
			expanded[len]=0;
			
//			System.err.println(len+", "+pos+": "+Arrays.toString(counts));
		}
		
//		System.err.println(targetKeys+", "+counts.length+", "+valid+", "+len+", "+rounds+", "+shared.Vector.sum(counts));
		
		return rounds;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_, final int tid_, final int size){
			cris=cris_;
			ros=ros_;
			tid=tid_;
			localHeap=new SketchHeap(size, 0, true);
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			dumpHeap();
			
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
					
					//Reads are processed in this block.
					processReadPair(r1, r2);
				}
				
				if(ros!=null){
					for(Read r1 : reads){
						readsOutT+=r1.pairCount();
						basesOutT+=r1.pairLength();
					}

					//Output reads to the output stream
					if(ros!=null){ros.add(reads, ln.id);}
				}

				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				if(ln.list!=null){ln.list.clear();}
				cris.returnList(ln.id, true);
			}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 */
		void processReadPair(final Read r1, final Read r2){
			processReadNucleotide(r1);
			if(r2!=null){processReadNucleotide(r2);}
		}
		
		void processReadNucleotide(final Read r){
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			long kmer=0;
			long rkmer=0;
			int len=0;
			assert(!r.aminoacid());
			
			final long min=minHashValue;
			localHeap.genomeSizeBases+=r.length();
			localHeap.genomeSequences++;
			
			if(quals==null || (minProb<=0 && minQual<2)){
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=AminoAcid.baseToNumber[b];
					long x2=AminoAcid.baseToComplementNumber[b];
					
					kmer=((kmer<<2)|x)&mask;
					rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
					
					if(x<0){len=0; rkmer=0;}else{len++;}
					if(len>=k){
						localHeap.genomeSizeKmers++;
						final long hashcode=hash(kmer, rkmer);
						if(hashcode>min){localHeap.add(hashcode);}
					}
				}
			}else{
				float prob=1;
				for(int i=0; i<bases.length; i++){
					final byte b=bases[i];
					final long x=AminoAcid.baseToNumber[b];
					final long x2=AminoAcid.baseToComplementNumber[b];
					
					{//Quality-related stuff
						final byte q=quals[i];
						assert(q>=0) : Arrays.toString(quals)+"\n"+minProb+", "+minQual;
						prob=prob*align2.QualityTools.PROB_CORRECT[q];
						if(len>k){
							byte oldq=quals[i-k];
							prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
						}
						if(x<0 || q<minQual){
							len=0;
							kmer=rkmer=0;
							prob=1;
						}else{
							len++;
						}
					}
					
					kmer=((kmer<<2)|x)&mask;
					rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
					
					if(len>=k && prob>=minProb){
						localHeap.genomeSizeKmers++;
						localHeap.probSum+=prob;
						final long hashcode=hash(kmer, rkmer);
						if(hashcode>min){localHeap.checkAndAdd(hashcode);}
					}
				}
			}
		}
		
		private void dumpHeap(){
			synchronized(sharedHeap){
				sharedHeap.add(localHeap);
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
		
		/** True only if this thread has completed successfully */
		boolean success=false;

		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Thread ID */
		final int tid;
		
		final SketchHeap localHeap;
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
	
	private boolean paired=false;
	private int trials=25;
	private long seed=-1;
	private int maxExpandedLength=50000000;
	
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
	
	private final SketchHeap sharedHeap;
	private final int heapSize;
	private final long targetKmers;
	private final int minCount;

	final int shift;
	final int shift2;
	final long mask;
	
	final float minProb;
	final byte minQual;
	
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
	/** Reads are output in input order (not enabled) */
	private boolean ordered=true;
	
}
