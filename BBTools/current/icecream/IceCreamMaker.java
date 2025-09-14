package icecream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.atomic.AtomicLong;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
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
import stream.SamLine;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Generates chimeric PacBio reads containing inverted repeats
 * due to missing adapters.
 * 
 * @author Brian Bushnell
 * @date June 8, 2019
 *
 */
public class IceCreamMaker {
	
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
	 * Constructor.
	 * @param args Command line arguments
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
			}
			
			else if(a.equals("idhist")){
				outIdHist=b;
			}
			
			else if(a.equals("minlength") || a.equals("minlen")){
				minMoleculeLength=Parse.parseIntKMG(b);
			}else if(a.equals("maxlength") || a.equals("maxlen")){
				maxMoleculeLength=Parse.parseIntKMG(b);
			}else if(a.equals("length") || a.equals("len")){
				minMoleculeLength=maxMoleculeLength=Parse.parseIntKMG(b);
			}
			
			else if(a.equals("minmovie") || a.equals("minmov")){
				minMovie=Parse.parseIntKMG(b);
			}else if(a.equals("maxmovie") || a.equals("maxmov")){
				maxMovie=Parse.parseIntKMG(b);
			}else if(a.equals("movie") || a.equals("mov")){
				minMovie=maxMovie=Parse.parseIntKMG(b);
			}
			
			else if(a.equals("missingrate") || a.equals("missing")){
				missingRate=Double.parseDouble(b);
				assert(missingRate<=1);
			}else if(a.equals("hiddenrate") || a.equals("hidden")){
				hiddenRate=Double.parseDouble(b);
				assert(hiddenRate<=1);
			}else if(a.equals("bothends")){
				allowBothEndsBad=Parse.parseBoolean(b);
				assert(false) : "TODO";
			}
			
			else if(a.equals("gc")){
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
			}
			
			else if(a.equals("invertedrepeatrate") || a.equals("invertrepeatrate") || a.equals("irrate")){
				invertedRepeatRate=Double.parseDouble(b);
				assert(invertedRepeatRate>=0);
			}else if(a.equals("invertedrepeatminlen") || a.equals("invertrepeatminlen") || a.equals("irminlen")){
				invertedRepeatMinLength=Parse.parseIntKMG(b);
			}else if(a.equals("invertedrepeatmaxlen") || a.equals("invertrepeatmaxlen") || a.equals("irmaxlen")){
				invertedRepeatMaxLength=Parse.parseIntKMG(b);
			}else if(a.equals("invertedrepeatlen") || a.equals("invertrepeatlen") || a.equals("irlen")){
				invertedRepeatMinLength=invertedRepeatMaxLength=Parse.parseIntKMG(b);
			}
			
			else if(a.equals("miner") || a.equals("minerrorrate")){
				minErrorRate=(float)Double.parseDouble(b);
				assert(minErrorRate>=0 && minErrorRate<=1);
			}else if(a.equals("maxer") || a.equals("maxerrorrate")){
				maxErrorRate=(float)Double.parseDouble(b);
				assert(maxErrorRate>=0 && maxErrorRate<=1);
			}else if(a.equals("er") || a.equals("errorrate")){
				minErrorRate=maxErrorRate=(float)Double.parseDouble(b);
				assert(minErrorRate>=0 && minErrorRate<=1);
			}
			
			else if(a.equals("minid") || a.equals("minidentity")){
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
			}
			
			else if(a.equals("parse_flag_goes_here")){
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

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=true;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris();
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros(false);
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		Random randy=Shared.threadLocalRandom(seed);
		if(cris==null){
			ref=genSynthGenome(randy);
		}else{
			ref=loadData(cris, randy);
		}
		
		if(invertedRepeatRate>0){
			addInvertedRepeats(ref, randy);
		}
		
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
	
	/** Create a Read Input Stream */
	private ConcurrentReadInputStream makeCris(){
		if(ffin1==null){return null;}
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/** Create a Read Output Stream */
	private ConcurrentReadOutputStream makeCros(boolean pairedInput){
		if(ffout1==null){return null;}

		//Set output buffer size
		final int buff=4;

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(
				ffout1, null, qfout1, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/** Spawn process threads */
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
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
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
	
	private byte randomBase(Random randy) {
		float rGC=randy.nextFloat();
		if(rGC>=genomeGC){//AT
			return (byte)(randy.nextBoolean() ? 'A' : 'T');
		}else{//GC
			return (byte)(randy.nextBoolean() ? 'G' : 'C');
		}
	}
	
	private static int randomLength(int min, int max, Random randy) {
		if(min==max){return min;}
		int range=max-min+1;
		int x=min+Tools.min(randy.nextInt(range), randy.nextInt(range));
//		System.err.println(x+", "+min+", "+max);
//		new Exception().printStackTrace();
//		System.err.println(randy.getClass());
//		System.err.println(randy.nextLong());
		return x;
	}
	
	private static float randomRate(float min, float max, Random randy) {
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
	
	private byte[] genSynthGenome(Random randy){
		assert(genomeSize<=MAX_GENOME_LENGTH) : genomeSize;
		byte[] array=new byte[(int)genomeSize];
		for(int i=0; i<genomeSize; i++){
			array[i]=randomBase(randy);
		}
		return array;
	}
	
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
	
	private void addInvertedRepeats(byte[] bases, Random randy){
		
		long added=0;
		long toAdd=(long) (bases.length*invertedRepeatRate);
		while(added<toAdd){
			int len=randomLength(invertedRepeatMinLength, invertedRepeatMaxLength, randy);
			int start=randy.nextInt(bases.length-2*len);
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
	
	/** */
	private class ProcessThread extends Thread {
		
		//Constructor
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
		
		/** Iterate through the reads */
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
		
		/** Generate the next list of reads */
		private ArrayList<Read> generateList(int toGenerate, long nextID){

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
		 * Generate a single read.
		 */
		private ArrayList<Read> generateZMW(final long nextID){
			
			final int movieLength=randomLength(minMovie, maxMovie, randy);
			final float errorRate=randomRate(minErrorRate, maxErrorRate, randy);
			
			ArrayList<ReadBuilder> baseCalls=baseCallAllPasses(movieLength, errorRate, nextID);
			if(baseCalls==null){
//				System.err.println(movieLength+", "+errorRate+", "+nextID);//123
				return null;
			}
			
			final boolean missing=randy.nextFloat()<missingRate;
			if(missing){
				final int missingMod=randy.nextInt(2);
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
				AminoAcid.reverseComplementBasesInPlace(frag);
			}
			return list;
		}
		
		/** Call bases for one pass */
		private ReadBuilder baseCallOnePass(final byte[] frag, final float errorRate, final int start, final int movieStart, int movieRemaining, long zmw){
			final float mult=1/(1-errorRate);
			ByteBuilder bb=new ByteBuilder();
			int fpos=start;
			for(; fpos<frag.length && movieRemaining>0; fpos++){
				float f=randy.nextFloat();
				byte b=frag[fpos];
				if(f>=errorRate){
					bb.append(b);
					movieRemaining--;
				}else{
					f=mult*(1-f);
					if(f<insThresh){//ins
						b=AminoAcid.numberToBase[randy.nextInt(4)];
						bb.append(b);
						fpos--;
						movieRemaining--;
					}else if(f<delThresh){//del
						
					}else{//sub
						int x=AminoAcid.baseToNumber[b]+randy.nextInt(3)+1;
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
		
		private byte[] fetchBases(byte[] source){
			
			final int len=randomLength(minMoleculeLength, maxMoleculeLength, randy);
			int start=len>=source.length ? 0 : randy.nextInt(source.length-len);
			int stop=Tools.min(start+len, source.length);
			
//			System.err.println("fetchBases(len="+len+", slen="+source.length+", start="+start+", stop="+stop);//123
			
			for(int i=start; i<stop; i++){
				if(!AminoAcid.isFullyDefined(source[i])){
					return null;
				}
			}
			if(stop-start<1){return null;}
			byte[] frag=Arrays.copyOfRange(source, start, stop);
			if(randy.nextBoolean()){AminoAcid.reverseComplementBasesInPlace(frag);}
			return frag;
		}
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		private final AtomicLong atomicZmwID;
		private final int readsPerList=Shared.bufferLen();
		
		/** Random number source */
		private Random randy;
		
		/** Random number source */
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

	/** Primary output file path */
	private String outIdHist=null;

	private String qfout1=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	/**  */
	private int minMoleculeLength=500;
	/**  */
	private int maxMoleculeLength=10000;
	/**  */
	private int minMovie=500;
	/**  */
	private int maxMovie=40000;
	
	/**  */
	private double missingRate=0.0;
	/**  */
	private double hiddenRate=0.0;
	/**  */
	private boolean allowBothEndsBad=false;
	
	/**  */
	private float genomeGC=0.6f;
	/**  */
	private long genomeSize=10000000;
	/**  */
	private boolean addNs=true;
	
	/**  */
	private double invertedRepeatRate=0.0;
	/**  */
	private int invertedRepeatMinLength=100;
	/**  */
	private int invertedRepeatMaxLength=10000;
	
	/**  */
	private float minErrorRate=0.05f;
	/**  */
	private float maxErrorRate=0.25f;
	/**  */
	private boolean addErrors=true;
	/** One read per ZMW */
	private boolean makeCCS=false;
	
	/** */
	private long seed=-1;
	
	private long[] idHist=new long[ID_BINS]; 
	
	//These should add to 1
	private float insFraction=0.40f;
	private float delFraction=0.35f;
	private float subFraction=0.25f;

	private final float insThresh;
	private final float delThresh;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	/** Quit after processing this many INPUT reads */
	private long maxReads=-1;
	
	/** Quit after generating this many OUTPUT zmws */
	private long maxZMWs=-1;
	
	/** Reference genome, max 2Gbp */
	private byte[] ref;
	
	private AtomicLong nextZmwID=new AtomicLong(0);
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	
	/** Primary output file */
	private final FileFormat ffout1;
	
	/** Primary output file */
	private final FileFormat ffIdHist;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final int ID_BINS=201;
	
	private static final long MAX_GENOME_LENGTH=2000000000;

	public static final byte[] pacbioAdapter="ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT".getBytes();
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
