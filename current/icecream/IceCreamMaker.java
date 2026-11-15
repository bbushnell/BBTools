package icecream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import shared.Random;
import java.util.concurrent.atomic.AtomicLong;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Simulates PacBio subreads and truth-labelled adapter failures for IceCreamFinder.
 * Uses a synthetic genome or concatenated input sequences, samples molecule and
 * movie lengths, and alternates strand on successive passes. Missing adapters
 * join passes directly; hidden adapters join them with simulated adapter bases.
 * Emitted headers carry simulation metadata for IceCreamGrader, not measured
 * sequencing accuracy. Reads have no computed base qualities.
 * The ccs option retains one internal median-length simulated read per well;
 * it does not compute a consensus. Instances change shared I/O settings and are
 * intended for one CLI invocation, not concurrent independent runs in one JVM.
 * @author Brian Bushnell, Shinobu
 * @date June 8, 2019
 *
 */
public class IceCreamMaker{
	
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
		IceCreamMaker x=new IceCreamMaker(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Parses configuration, validates ranges and prepares file descriptions.
	 * Does not load the reference or start workers. Also changes shared I/O,
	 * fake-quality and read-format settings used by this CLI invocation.
	 * @param args BBTools key=value options, including a positive zmws count
	 */
	public IceCreamMaker(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		Shared.FAKE_QUAL=8;
//		Shared.FASTA_WRAP=511;
		SamLine.SET_FROM_OK=true;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			extin=parser.extin;

			out1=parser.out1;
			qfout1=parser.qfout1;
			extout=parser.extout;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		
		ffIdHist=FileFormat.testOutput(outIdHist, FileFormat.TXT, null, false, overwrite, append, false);
		
		insThresh=insFraction;
		delThresh=delFraction+insThresh;
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
			}else if(a.equals("idhist")){
				outIdHist=b;
			}else if(a.equals("minlength") || a.equals("minlen")){
				minMoleculeLength=Parse.parseIntKMG(b);
			}else if(a.equals("maxlength") || a.equals("maxlen")){
				maxMoleculeLength=Parse.parseIntKMG(b);
			}else if(a.equals("length") || a.equals("len")){
				minMoleculeLength=maxMoleculeLength=Parse.parseIntKMG(b);
			}else if(a.equals("minmovie") || a.equals("minmov")){
				minMovie=Parse.parseIntKMG(b);
			}else if(a.equals("maxmovie") || a.equals("maxmov")){
				maxMovie=Parse.parseIntKMG(b);
			}else if(a.equals("movie") || a.equals("mov")){
				minMovie=maxMovie=Parse.parseIntKMG(b);
			}else if(a.equals("missingrate") || a.equals("missing")){
				missingRate=Double.parseDouble(b);
				assert(missingRate<=1);
			}else if(a.equals("hiddenrate") || a.equals("hidden")){
				hiddenRate=Double.parseDouble(b);
				assert(hiddenRate<=1);
			}else if(a.equals("bothends")){
				allowBothEndsBad=Parse.parseBoolean(b);
				assert(false) : "TODO";
			}else if(a.equals("gc")){
				genomeGC=(float)Double.parseDouble(b);
				assert(genomeGC>=0 && genomeGC<=1);
			}else if(a.equals("genomesize")){
				genomeSize=Parse.parseKMG(b);
			}else if(a.equals("addns") || a.equals("ns")){
				addNs=Parse.parseBoolean(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("zmws") || a.equals("maxzmws") || a.equals("reads")){
				maxZMWs=Parse.parseIntKMG(b);
			}else if(a.equals("ccs")){
				makeCCS=Parse.parseBoolean(b);
			}else if(a.equals("invertedrepeatrate") || a.equals("invertrepeatrate") || a.equals("irrate")){
				invertedRepeatRate=Double.parseDouble(b);
				assert(invertedRepeatRate>=0);
			}else if(a.equals("invertedrepeatminlen") || a.equals("invertrepeatminlen") || a.equals("irminlen")){
				invertedRepeatMinLength=Parse.parseIntKMG(b);
			}else if(a.equals("invertedrepeatmaxlen") || a.equals("invertrepeatmaxlen") || a.equals("irmaxlen")){
				invertedRepeatMaxLength=Parse.parseIntKMG(b);
			}else if(a.equals("invertedrepeatlen") || a.equals("invertrepeatlen") || a.equals("irlen")){
				invertedRepeatMinLength=invertedRepeatMaxLength=Parse.parseIntKMG(b);
			}else if(a.equals("miner") || a.equals("minerrorrate")){
				minErrorRate=(float)Double.parseDouble(b);
				assert(minErrorRate>=0 && minErrorRate<=1);
			}else if(a.equals("maxer") || a.equals("maxerrorrate")){
				maxErrorRate=(float)Double.parseDouble(b);
				assert(maxErrorRate>=0 && maxErrorRate<=1);
			}else if(a.equals("er") || a.equals("errorrate")){
				minErrorRate=maxErrorRate=(float)Double.parseDouble(b);
				assert(minErrorRate>=0 && minErrorRate<=1);
			}else if(a.equals("minid") || a.equals("minidentity")){
				maxErrorRate=1-(float)Double.parseDouble(b);
				assert(maxErrorRate>=0 && maxErrorRate<=1);
			}else if(a.equals("maxid") || a.equals("maxidentity")){
				minErrorRate=1-(float)Double.parseDouble(b);
				assert(minErrorRate>=0 && minErrorRate<=1);
			}else if(a.equals("id") || a.equals("identity")){
				minErrorRate=maxErrorRate=1-(float)Double.parseDouble(b);
				assert(minErrorRate>=0 && minErrorRate<=1);
			}else if(a.equals("adderrors")){
				addErrors=Parse.parseBoolean(b);
			}else if(a.equals("ref")){
				parser.in1=b;
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
		in1=Tools.fixExtension(in1);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, outIdHist)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1, outIdHist)){
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
		assert(minMoleculeLength>0 && minMoleculeLength<=maxMoleculeLength) : minMoleculeLength+", "+maxMoleculeLength;
		assert(minMovie>0 && minMovie<=maxMovie) : minMovie+", "+maxMovie;

		assert(missingRate>=0 && missingRate<=1) : missingRate;
		assert(hiddenRate>=0 && hiddenRate<=1) : hiddenRate;
		assert(genomeGC>=0 && genomeGC<=1) : genomeGC;
		assert(in1!=null || genomeSize>maxMoleculeLength) : genomeSize;
		assert(in1!=null || invertedRepeatMaxLength*2<genomeSize) : genomeSize;
		assert(maxZMWs>0) : "zmsw="+maxZMWs+"; please set to a positive number.";

		assert(invertedRepeatRate>=0) : invertedRepeatRate;
		assert(invertedRepeatMinLength>0 && invertedRepeatMinLength<=invertedRepeatMaxLength) : invertedRepeatMinLength+", "+invertedRepeatMaxLength;

		assert(minErrorRate>=0 && minErrorRate<=maxErrorRate) : minErrorRate+", "+maxErrorRate;
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads or generates the reference, starts simulation workers and writes reports.
	 * Closes normally completed streams and checks accumulated worker/writer errors.
	 * @param t Running timer stopped before the final report
	 */
	void process(Timer t){
		
		//Enable input validation while loading the reference; restore it after normal completion.
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=true;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris();
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		Random randy=Shared.threadLocalRandom(seed);
		if(cris==null){
			ref=genSynthGenome(randy);
		}else{
			ref=loadData(cris, randy);
		}
		
		try{
			if(invertedRepeatRate>0){addInvertedRepeats(ref, randy);}
			validateSimulation();
		}catch(RuntimeException e){
			try{errorState|=ReadWrite.closeStream(cris);}
			finally{Read.VALIDATE_IN_CONSTRUCTOR=vic;}
			throw e;
		}
		//Do not start an output worker until reference and movie feasibility checks pass.
		final ConcurrentReadOutputStream ros=makeCros(false);
		
		//Process the reads in separate threads
		spawnThreads(cris, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		writeIdHist();
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsOut, basesOut, 8));
//		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Rejects impossible molecule/CCS requests without drawing from the simulator RNG.
	 * These necessary bounds do not guarantee a short wait for a very unlikely
	 * accepted candidate. The final transformed reference is checked here.
	 */
	private void validateSimulation(){
		final int minFragment=Tools.min(minMoleculeLength, ref.length);
		if(minFragment<1 || !hasUsableRun(ref, minFragment, true)){
			throw new IllegalArgumentException("Reference has no fully defined fragment of at least "+
					minFragment+" bases; simulation cannot generate a well.");
		}
		if(makeCCS){
			//A float draw cannot exceed nextDown(1f); strict inequality matters at that endpoint.
			final int passes=missingRate>Math.nextDown(1.0f) ? 4 : 3;
			final boolean deletions=addErrors && maxErrorRate>0 && delThresh>insThresh;
			final long minimumMovie=(passes-1L)*adapterLen+
					(deletions ? 1 : 2+(passes-2L)*minFragment);
			if(maxMovie<minimumMovie){
				throw new IllegalArgumentException("CCS requires at least "+passes+
						" raw passes; maxmovie="+maxMovie+" is below the minimum possible budget "+minimumMovie+".");
			}
		}
	}

	/** Stops at the first eligible run; normal references need only a short prefix scan.
	 * @param bases Reference to inspect without mutation
	 * @param length Positive required run length
	 * @param fullyDefined True for molecule sampling; false for repeat insertion's literal N exclusion
	 * @return Whether at least one eligible run exists
	 */
	private static boolean hasUsableRun(byte[] bases, int length, boolean fullyDefined){
		assert(length>0) : "A candidate window must contain at least one base: "+length;
		int run=0;
		for(byte b : bases){
			run=(fullyDefined ? AminoAcid.isFullyDefined(b) : b!='N') ? run+1 : 0;
			if(run>=length){return true;}
		}
		return false;
	}

	/**
	 * Writes one histogram observation per generated well, using its requested
	 * error rate rather than recounting realized substitutions or indels.
	 */
	private void writeIdHist(){
		if(ffIdHist==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(ffIdHist);
		bsw.start();
		bsw.print("#Identity\tCount\n".getBytes());
		for(int i=0; i<idHist.length; i++){
			bsw.print(i*100.0/(ID_BINS-1), 3).print('\t').println(idHist[i]);
		}
		errorState|=bsw.poisonAndWait();
	}
	
	/** Starts the optional reference input backend, or returns null for a synthetic reference. */
	private ConcurrentReadInputStream makeCris(){
		if(ffin1==null){return null;}
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/** Starts the optional unordered subread output backend, or returns null when output is disabled. */
	private ConcurrentReadOutputStream makeCros(boolean pairedInput){
		if(ffout1==null){return null;}

		//Set output buffer size
		final int buff=4;

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(
				ffout1, null, qfout1, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/** Starts simulation workers, joins them, and accumulates their counters and success flags. */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(ros, i, nextZmwID, seed));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for threads to finish
		waitForThreads(alpt);
		
		//Do anything necessary after processing
		
	}
	
	/** Wait until all worker threads are finished, then return */
	private void waitForThreads(ArrayList<ProcessThread> alpt){
		
		//Wait for completion of all threads
		boolean success=true;
		for(ProcessThread pt : alpt){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try{
					//Attempt a join operation
					pt.join();
				}catch(InterruptedException e){
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			Tools.add(idHist, pt.idHistT);
			
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Generate random base according to genome GC content.
	 * @param randy Random number generator
	 * @return Random nucleotide base (A, T, G, C)
	 */
	private byte randomBase(Random randy){
		float rGC=randy.nextFloat();
		if(rGC>=genomeGC){//AT
			return (byte)(randy.nextBoolean() ? 'A' : 'T');
		}else{//GC
			return (byte)(randy.nextBoolean() ? 'G' : 'C');
		}
	}
	
	/**
	 * Generate random length with bias toward shorter values.
	 * @param min Minimum length
	 * @param max Maximum length
	 * @param randy Random number generator
	 * @return Random length between min and max (inclusive)
	 */
	private static int randomLength(int min, int max, Random randy){
		if(min==max){return min;}
		int range=max-min+1;
		int x=min+Tools.min(randy.nextInt(range), randy.nextInt(range));
//		System.err.println(x+", "+min+", "+max);
//		new Exception().printStackTrace();
//		System.err.println(randy.getClass());
//		System.err.println(randy.nextLong());
		return x;
	}
	
	/**
	 * Generate random rate with bias toward lower values using multiple random draws.
	 * @param min Minimum rate
	 * @param max Maximum rate
	 * @param randy Random number generator
	 * @return Random rate between min and max
	 */
	private static float randomRate(float min, float max, Random randy){
		if(min==max){return min;}
		float range=max-min;
//		float a=(randy.nextFloat()+randy.nextFloat());
		float b=(randy.nextFloat()+randy.nextFloat());
		float c=(1.6f*randy.nextFloat()+0.4f*randy.nextFloat());
		float x=min+range*0.5f*Tools.min(b, c);
//		assert(false) : "x="+x+", a="+a+", b="+b+", c="+c+", min="+min+", max="+max;
//		System.err.println(x+", "+min+", "+max);
		return x;
	}
	
	/**
	 * Generate synthetic genome with random bases according to GC content.
	 * @param randy Random number generator
	 * @return Synthetic genome sequence
	 */
	private byte[] genSynthGenome(Random randy){
		assert(genomeSize<=MAX_GENOME_LENGTH) : genomeSize;
		byte[] array=new byte[(int)genomeSize];
		for(int i=0; i<genomeSize; i++){
			array[i]=randomBase(randy);
		}
		return array;
	}
	
	/**
	 * Concatenates primary input sequences up to MAX_GENOME_LENGTH. With addNs, inserts an
	 * N separator and preserves undefined bases; otherwise replaces undefined bases
	 * in the input arrays with sampled nucleotides before appending them.
	 * @param cris Input stream containing reference sequences
	 * @param randy Random number generator for base replacement
	 * @return Concatenated reference sequence
	 */
	private byte[] loadData(ConcurrentReadInputStream cris, Random randy){

		ByteBuilder bb=new ByteBuilder();
		
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();

		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){
//			if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
			
			for(Read r : ln){
				
//				System.err.println("Fetched read len="+r.length());//123
				
				//Increment counters
				readsProcessed+=r.pairCount();
				basesProcessed+=r.pairLength();
				//TODO: Probable bug [PIC-015] - paired reference input counts mates here but appends only primary bases below.
				//Normal reference FASTA is unpaired; paired-reference support needs a separate explicit contract and finite case.
				
				

//				System.err.println("Filter: addNs="+addNs+", (r.length()>maxMoleculeLength && (invertedRepeatRate<=0 || r.length()>2*invertedRepeatMaxLength);//123
				
				//Filter disabled; it causes short sequences to be ignored.
				//Optional filter criteria
//				if(!addNs || (r.length()>maxMoleculeLength && (invertedRepeatRate<=0 || r.length()>2*invertedRepeatMaxLength))){
					final byte[] bases=r.bases;
					
					if(addNs){
						if(bb.length()>0){bb.append('N');}
					}else{
						for(int i=0; i<bases.length; i++){
							if(!AminoAcid.isFullyDefined(bases[i])){
								bases[i]=randomBase(randy);
							}
						}
					}
					
					if(bb.length()+bases.length<=MAX_GENOME_LENGTH){
						bb.append(bases);
//						System.err.println("Appended "+r.length());//123
					}else{
//						System.err.println("Appending partial.");//123
						for(byte b : bases){
							bb.append(b);
							if(bb.length>=MAX_GENOME_LENGTH){
								cris.returnList(ln.id, true);
//								System.err.println("Returning partial "+bb.length);//123
								return bb.toBytes();
							}
						}
					}
//				}
			}

			//Notify the input stream that the list was used
			cris.returnList(ln);
//			if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
			
			//Fetch a new list
			ln=cris.nextList();
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}

//		System.err.println("Returning full "+bb.length);//123
		return bb.toBytes();
	}
	
	/**
	 * Add inverted repeats to reference sequence at random positions. Creates
	 * chimeric sequences by copying and reverse-complementing regions.
	 * @param bases Reference sequence to modify
	 * @param randy Random number generator for positioning
	 */
	private void addInvertedRepeats(byte[] bases, Random randy){
		long added=0;
		long toAdd=(long) (bases.length*invertedRepeatRate);
		if(toAdd<=0){return;}
		if(2L*invertedRepeatMaxLength>bases.length){
			throw new IllegalArgumentException("Reference length "+bases.length+
					" cannot fit every requested inverted-repeat pair; irmaxlen="+invertedRepeatMaxLength+".");
		}
		if(!hasUsableRun(bases, 2*invertedRepeatMinLength, false)){
			throw new IllegalArgumentException("Reference has no N-free window for an inverted-repeat pair of "+
					(2L*invertedRepeatMinLength)+" bases.");
		}
		while(added<toAdd){
			int len=randomLength(invertedRepeatMinLength, invertedRepeatMaxLength, randy);
			int start=randy.nextInt(bases.length-2*len+1);//PIC-011: include the last fitting window and exact fits.
			int stop=start+len;
			boolean OK=true;
			for(int i=0; i<len && OK; i++){
				OK&=(bases[start+i]!='N' && bases[stop+i]!='N');
			}
			if(OK){
				for(int i=0; i<len; i++){
					byte b=bases[stop-i-1];
					bases[stop+i]=AminoAcid.baseToComplementExtended[b];
				}
				added+=2*len;
//				System.err.println("added="+added+"/"+toAdd);
			}else{
//				
			}
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Owns one RNG and counters; claims well-ID batches through the shared atomic counter. */
	private class ProcessThread extends Thread{
		
		//Constructor
		/**
		 * Constructor for ProcessThread worker.
		 * @param ros_ Output stream for generated reads
		 * @param tid_ Thread ID
		 * @param nextZmwID_ Atomic counter for ZMW IDs
		 * @param seed Legacy unused argument; run uses the enclosing instance's configured seed
		 */
		ProcessThread(final ConcurrentReadOutputStream ros_, final int tid_, 
				final AtomicLong nextZmwID_, final long seed){
			ros=ros_;
			tid=tid_;
			atomicZmwID=nextZmwID_;
//			assert(false) : randy.getClass()+", "+randy.nextLong();
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			randy=Shared.threadLocalRandom(seed<0 ? seed : seed+(tid+1)*tid*1000L);
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Claims well IDs and emits unordered batches containing all their generated subreads. */
		void processInner(){

			//As long as there is a nonempty read list...
			for(long generated=atomicZmwID.getAndAdd(readsPerList); generated<maxZMWs; 
					generated=atomicZmwID.getAndAdd(readsPerList)){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

				long toGenerate=Tools.min(readsPerList, maxZMWs-generated);
				
				ArrayList<Read> reads=generateList((int)toGenerate, generated);

				//Output reads to the output stream
				if(ros!=null){ros.add(reads, 0);}
			}
		}
		
		/** Generates the requested number of successful wells, retrying rejected candidates.
		 * @param toGenerate Number of wells, not the number of emitted subreads
		 * @param nextID First well ID in this batch
		 * @return Concatenated subreads for the batch
		 */
		private ArrayList<Read> generateList(int toGenerate, long nextID){
			//TODO PIC-005: preflight excludes known impossible requests, but rare valid candidates can still take arbitrarily many retries.

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=new ArrayList<Read>(toGenerate);
			
			//Loop through each read in the list
			for(int i=0; i<toGenerate; i++, nextID++){
				ArrayList<Read> zmw=generateZMW(nextID);
				if(zmw==null){
					i--;
					nextID--;
				}else{reads.addAll(zmw);}
			}

			return reads;
		}
		
		/**
		 * Find median-length read from list, excluding first and last elements.
		 * @param list List of ReadBuilder objects
		 * @return First internal read at the lower median length, or null for fewer than three reads
		 */
		private ReadBuilder median(ArrayList<ReadBuilder> list){
			if(list.size()<3){return null;}
			IntList lengths=new IntList(list.size()-2);
			
			for(int i=1; i<list.size()-1; i++){
				lengths.add(list.get(i).length());
			}
			lengths.sort();
			int median=lengths.get((lengths.size-1)/2);
			
			for(int i=1; i<list.size()-1; i++){
				ReadBuilder rb=list.get(i);
				if(rb.length()==median){
					return rb;
				}
			}
			
			assert(false);
			return null;
		}
		
		/**
		 * Simulates a well, applies adapter joins, and optionally retains one internal
		 * median-length subread. Successful wells update output counts and the requested
		 * identity histogram. Rejected fragment/CCS candidates return null for retry.
		 * @param nextID Well identifier stored in every generated header
		 * @return Generated subreads, or null when this candidate is rejected
		 */
		private ArrayList<Read> generateZMW(final long nextID){
			
			final int movieLength=randomLength(minMovie, maxMovie, randy);
			final float errorRate=addErrors ? randomRate(minErrorRate, maxErrorRate, randy) : 0;
			
			ArrayList<ReadBuilder> baseCalls=baseCallAllPasses(movieLength, errorRate, nextID);
			if(baseCalls==null){
//				System.err.println(movieLength+", "+errorRate+", "+nextID);//123
				return null;
			}
			
			final boolean missing=randy.nextFloat()<missingRate;
			if(missing){
				final int missingMod=randy.nextInt2();
				ArrayList<ReadBuilder> temp=new ArrayList<ReadBuilder>();
				
				ReadBuilder current=null;
				for(int i=0; i<baseCalls.size(); i++){
					ReadBuilder rb=baseCalls.get(i);
					assert(rb.subreads==1);
					assert(rb.missing==0);
					assert(rb.adapters==0);
					if(current!=null && (i&1)==missingMod){
						current.add(rb);
						current.missing++;
						assert(current.subreads>1);
						assert(current.missing>0);
						temp.add(current);
						current=null;
					}else{
						if(current!=null){temp.add(current);}
						current=rb;
					}
				}
				if(current!=null){temp.add(current);}
				baseCalls=temp;
			}
			
			if(makeCCS){
				ReadBuilder median=median(baseCalls);
				if(median==null){return null;}
				baseCalls.clear();
				baseCalls.add(median);
			}else if(hiddenRate>0){
				ArrayList<ReadBuilder> temp=new ArrayList<ReadBuilder>();
				
				ReadBuilder current=null;
				for(int i=0; i<baseCalls.size(); i++){
					ReadBuilder rb=baseCalls.get(i);
					assert(rb.adapters==0);
					if(current!=null && randy.nextFloat()<hiddenRate){
						current.add(baseCallAdapter(errorRate));
						current.add(rb);
						assert(current.adapters>0);
						assert(current.subreads>1);
					}else{
						if(current!=null){temp.add(current);}
						current=rb;
						assert(current.adapters==0);
					}
				}
				if(current!=null){temp.add(current);}
				baseCalls=temp;
			}
			
			
			ArrayList<Read> reads=new ArrayList<Read>();
			for(ReadBuilder rb : baseCalls){
				Read r=rb.toRead();
				readsOutT+=r.pairCount();
				basesOutT+=r.length();
				reads.add(r);
			}
			
			idHistT[(int)((1-errorRate)*(ID_BINS-1)+0.5f)]++;
			return reads;
		}
		
		/**
		 * Simulate base calling for all passes through a DNA fragment, alternating
		 * strand direction for each pass.
		 * @param movieLength Emitted-base budget plus fixed adapter cost per pass
		 * @param errorRate Error rate for this ZMW
		 * @param zmw ZMW identifier
		 * @return Pass builders, or null if ten fragment attempts find no eligible region
		 */
		private ArrayList<ReadBuilder> baseCallAllPasses(final int movieLength, final float errorRate, long zmw){
			byte[] frag=null;
			
			for(int i=0; i<10 && frag==null; i++){//retry several times
				frag=fetchBases(ref);
			}
			
			if(frag==null){
//				System.err.println("Failed baseCallAllPasses("+movieLength+", "+errorRate+", "+zmw);//123
				return null;
			}
			
			ArrayList<ReadBuilder> list=new ArrayList<ReadBuilder>();
			int movieRemaining=movieLength;
			int moviePos=0;
			assert(frag.length>0) : frag.length;
			int start=randy.nextInt(frag.length);
			while(movieRemaining>0){
				ReadBuilder rb=baseCallOnePass(frag, errorRate, start, moviePos, movieRemaining, zmw);
				list.add(rb);
				start=0;
				int elapsed=rb.length()+adapterLen;
				moviePos+=elapsed;
				movieRemaining-=elapsed;
				Vector.reverseComplementInPlaceFast(frag);
			}
			return list;
		}
		
		/** Calls one oriented fragment suffix with a bounded emitted-base budget.
		 * An insertion consumes budget without advancing the fragment; a deletion
		 * advances the fragment without consuming budget. The truth pass fraction
		 * uses traversed fragment bases, and qualities are not synthesized here.
		 * @param frag Nonempty oriented fragment, read without mutation
		 * @param errorRate Per-base error-event probability in [0,1]
		 * @param start First fragment position
		 * @param movieStart Starting coordinate recorded in the synthetic header
		 * @param movieRemaining Remaining emitted-base budget
		 * @param zmw Well identifier
		 * @return Mutable sequence builder with simulation metadata
		 */
		private ReadBuilder baseCallOnePass(final byte[] frag, final float errorRate, final int start, final int movieStart, int movieRemaining, long zmw){
			ByteBuilder bb=new ByteBuilder();
			int fpos=start;
			for(; fpos<frag.length && movieRemaining>0; fpos++){
				float f=randy.nextFloat();
				byte b=frag[fpos];
				if(f>=errorRate){
					bb.append(b);
					movieRemaining--;
				}else{
					//PIC-003: normalize within the error event; the former remap exceeded one and suppressed all indels.
					f/=errorRate;//This branch is unreachable at errorRate=0; division also avoids reciprocal overflow.
					if(f<insThresh){//ins
						b=AminoAcid.numberToBase[randy.nextInt()&3];
						bb.append(b);
						fpos--;
						movieRemaining--;
					}else if(f<delThresh){//del
						
					}else{//sub
						int x=AminoAcid.baseToNumber[b]+randy.nextInt3()+1;
						b=AminoAcid.numberToBase[x&3];
						bb.append(b);
						movieRemaining--;
					}
				}
			}
			
			float passes=(fpos-start)*1.0f/frag.length;
			ReadBuilder rb=new ReadBuilder(bb, passes, movieStart, zmw);
			rb.errorRate=errorRate;
			return rb;
		}
		
		/** Call bases for an adapter sequence pass */
		private ReadBuilder baseCallAdapter(final float errorRate){
			ReadBuilder rb=baseCallOnePass(pacbioAdapter, errorRate, 0, 0, 999999, 0);
			rb.passes=0;
			rb.fullPasses=0;
			rb.subreads=0;
			rb.adapters=1;
			return rb;
		}
		
		/**
		 * Extract random fragment from reference sequence, checking for defined bases
		 * and optionally reverse-complementing.
		 * @param source Reference sequence
		 * @return Random fragment or null if no valid region found
		 */
		private byte[] fetchBases(byte[] source){
			
			final int len=randomLength(minMoleculeLength, maxMoleculeLength, randy);
			int start=len>=source.length ? 0 : randy.nextInt(source.length-len+1);//PIC-010: the final start is eligible too.
			int stop=Tools.min(start+len, source.length);
			
//			System.err.println("fetchBases(len="+len+", slen="+source.length+", start="+start+", stop="+stop);//123
			
			for(int i=start; i<stop; i++){
				if(!AminoAcid.isFullyDefined(source[i])){
					return null;
				}
			}
			if(stop-start<1){return null;}
			byte[] frag=Arrays.copyOfRange(source, start, stop);
			if(randy.nextBoolean()){Vector.reverseComplementInPlaceFast(frag);}
			return frag;
		}
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Atomic ZMW ID counter shared across threads */
		private final AtomicLong atomicZmwID;
		/** Wells claimed per batch; each successful well may emit several subreads. */
		private final int readsPerList=Shared.bufferLen();
		
		/** Random number source */
		private Random randy;
		
		/** Requested-identity observations, accumulated once per successful well. */
		private final long[] idHistT=new long[ID_BINS];
		
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String out1=null;

	/** Requested-identity histogram output path. */
	private String outIdHist=null;

	/** Quality file output path */
	private String qfout1=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	/** Minimum requested fragment length; a shorter reference can truncate it. */
	private int minMoleculeLength=500;
	/** Maximum requested fragment length, inclusive. */
	private int maxMoleculeLength=10000;
	/** Minimum movie budget in emitted bases plus fixed adapter costs. */
	private int minMovie=500;
	/** Maximum movie budget, inclusive. */
	private int maxMovie=40000;
	
	/** Probability that a well receives alternating missing-adapter joins. */
	private double missingRate=0.0;
	/** Probability of each hidden-adapter join when ccs mode is disabled. */
	private double hiddenRate=0.0;
	/** Unimplemented option; explicitly parsing bothends currently asserts. */
	private boolean allowBothEndsBad=false;
	
	/** GC probability for synthetic reference bases and undefined-base replacement. */
	private float genomeGC=0.6f;
	/** Synthetic reference length; unused when reference input is supplied. */
	private long genomeSize=10000000;
	/** Preserve undefined input bases and insert N separators rather than replacing them. */
	private boolean addNs=true;
	
	/** Target reference fraction counted as both halves of inserted inverted repeats. */
	private double invertedRepeatRate=0.0;
	/** Minimum length of one inverted-repeat half. */
	private int invertedRepeatMinLength=100;
	/** Maximum length of one inverted-repeat half, inclusive. */
	private int invertedRepeatMaxLength=10000;
	
	/** Lower endpoint of the per-well error-rate distribution. */
	private float minErrorRate=0.05f;
	/** Upper endpoint of the per-well error-rate distribution. */
	private float maxErrorRate=0.25f;
	/** When false, generates error-free calls regardless of the configured rate range. */
	private boolean addErrors=true;
	/** Retain one internal median simulated subread; does not compute a consensus. */
	private boolean makeCCS=false;
	
	/** Configured seed; negative values request nondeterministic RNG initialization. */
	private long seed=-1;
	
	/** Requested identity histogram, one observation per successful well. */
	private long[] idHist=new long[ID_BINS]; 
	
	//These should add to 1
	/** Fraction of errors that are insertions */
	private float insFraction=0.40f;
	/** Fraction of errors that are deletions */
	private float delFraction=0.35f;
	/** Fraction of errors that are substitutions */
	private float subFraction=0.25f;

	/** Threshold for insertion errors */
	private final float insThresh;
	/** Threshold for deletion errors */
	private final float delThresh;
	
	/*--------------------------------------------------------------*/

	/** Reference input records counted, including mates */
	protected long readsProcessed=0;
	/** Reference input bases counted, including mates */
	protected long basesProcessed=0;

	/** Number of simulated subreads emitted */
	protected long readsOut=0;
	/** Number of simulated bases emitted */
	protected long basesOut=0;

	/** Maximum reference input fragments requested from the reader; negative means unlimited */
	private long maxReads=-1;
	
	/** Quit after generating this many OUTPUT zmws */
	private long maxZMWs=-1;
	
	/** Reference genome, max 2Gbp */
	private byte[] ref;
	
	/** Thread-safe counter for ZMW identifiers */
	private AtomicLong nextZmwID=new AtomicLong(0);
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	
	/** Primary output file */
	private final FileFormat ffout1;
	
	/** Output descriptor for the requested-identity histogram */
	private final FileFormat ffIdHist;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of bins for identity histogram */
	private static final int ID_BINS=201;
	
	/** Maximum reference genome length (2Gbp) */
	private static final long MAX_GENOME_LENGTH=2000000000;

	/** Standard PacBio adapter sequence */
	public static final byte[] pacbioAdapter="ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT".getBytes();
	/** Length of PacBio adapter sequence */
	public static final int adapterLen=pacbioAdapter.length;
	
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
	
}
