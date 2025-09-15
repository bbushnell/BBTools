package repeat;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.AbstractKmerTableSet;
import kmer.KmerTableSet;
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
import structures.ByteBuilder;
import structures.CRange;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;
import tracker.ReadStats;

/**
 * This class finds probable repeats using kmer-matching,
 * and prints them out ordered by depth, then length.
 * 
 * The repeats are only "probable" because no alignment is done;
 * a sequence is considered a probable repeat of depth D if
 * if all kmers within it have a depth of at least D.
 * 
 * @author Brian Bushnell
 * @date June 15, 2023
 *
 */
public class RepeatFinder implements Accumulator<RepeatFinder.ProcessThread> {
	
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
		RepeatFinder x=new RepeatFinder(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public RepeatFinder(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			extin=parser.extin;

			outr=parser.out1;
			extout=parser.extout;
			
			k=parser.k;
			assert(k>=1 && k<=32) : "1<=k<=31: k="+k;
			amino=Shared.AMINO_IN;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffoutr=FileFormat.testOutput(outr, FileFormat.TXT, null, true, overwrite, append, false);
		ffouts=FileFormat.testOutput(outs, FileFormat.FA, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
		
		tables=(processDepth ? new KmerTableSet(args, 12) : null);
		
		{//set some constants
			bitsPerBase=(amino ? 5 : 2);
			maxSymbol=(amino ? 20 : 3);
			symbols=maxSymbol+1;
			symbolArrayLen=(64+bitsPerBase-1)/bitsPerBase;
			symbolSpace=(1<<bitsPerBase);
			symbolMask=symbolSpace-1;
			
			symbolToNumber=AminoAcid.symbolToNumber(amino);
			symbolToNumber0=AminoAcid.symbolToNumber0(amino);
			symbolToComplementNumber0=AminoAcid.symbolToComplementNumber0(amino);
			
			clearMasks=new long[symbolArrayLen];
			leftMasks=new long[symbolArrayLen];
			rightMasks=new long[symbolArrayLen];
			setMasks=new long[symbols][symbolArrayLen];
			for(int i=0; i<symbolArrayLen; i++){
				clearMasks[i]=~(symbolMask<<(bitsPerBase*i));
				leftMasks[i]=((-1L)<<(bitsPerBase*i));
				rightMasks[i]=~((-1L)<<(bitsPerBase*i));
				for(long j=0; j<symbols; j++){
					setMasks[(int)j][i]=(j<<(bitsPerBase*i));
				}
			}
			
//			minlen=k-1;
			shift=bitsPerBase*k;
			shift2=shift-bitsPerBase;
			mask=(shift>63 ? -1L : ~((-1L)<<shift));
		}
		middleMask=(((k|1)==1) ? ~(3<<(k/2)) : ~(15<<(k/2-1)));//Works for even kmer lengths.
		
		baseMaskArray=new byte[127];
		for(int i=0; i<baseMaskArray.length; i++) {
			final byte b;
			if(!maskRepeats) {b=(byte)i;}
			else if(softMask) {b=Tools.toLowerCase((byte)i);
			}else{b=maskSymbol;}
			baseMaskArray[i]=b;
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
		parser.out1=outr;
		
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
			}else if(a.equals("displayprogress")){
				DISPLAY_PROGRESS=Parse.parseBoolean(b);
			}else if(a.equals("qhdist") || a.equals("hdist")){
				qHammingDist=Integer.parseInt(b);
			}else if(a.equals("maskmiddle") || a.equals("mm")){
				maskMiddle=AbstractKmerTableSet.MASK_MIDDLE=Parse.parseBoolean(b);
			}else if(a.equals("minrepeat") || a.equals("minlen") || a.equals("minlength")){
				minRepeat=Integer.parseInt(b);
			}else if(a.equals("maxgap") || a.equals("gap")){
				maxGap=Integer.parseInt(b);
			}else if(a.equals("maxdepth") || a.equals("maxcount")){
				maxDepth=Integer.parseInt(b);
				if(maxDepth<1) {maxDepth=Shared.MAX_ARRAY_LEN;}
			}else if(a.equals("mindepth") || a.equals("mincount")){
				minDepth=Tools.max(2, Integer.parseInt(b));
			}else if(a.equals("printlen") || a.equals("seqlen") || a.equals("println") || 
					a.equals("preview") || a.equals("printsequence") || a.equals("printseq")){
				if(b==null || Character.isLetter(b.charAt(0))){
					if(Parse.parseBoolean(b)) {
						Repeat.SEQ_AFFIX_LEN=12;
					}else{
						Repeat.SEQ_AFFIX_LEN=0;
					}
				}else{
					int x=Integer.parseInt(b);
					Repeat.SEQ_AFFIX_LEN=(x-3)/2;
				}
			}else if(a.equals("depth") || a.equals("processdepth")){
				processDepth=Parse.parseBoolean(b);
			}
			
			else if(a.equals("mask") || a.equals("maskrepeats")){
				if("lc".equalsIgnoreCase(b) || "lowercase".equalsIgnoreCase(b) ||
						"soft".equalsIgnoreCase(b)) {
					maskRepeats=true;
					softMask=true;
				}else if("hard".equalsIgnoreCase(b)) {
					maskRepeats=true;
					softMask=false;
				}else if("t".equalsIgnoreCase(b) || "f".equalsIgnoreCase(b) || b==null || b.length()>1) {
					maskRepeats=Parse.parseBoolean(b);
				}else{
					maskRepeats=true;
					softMask=false;
					maskSymbol=(byte)b.charAt(0);
				}
			}else if(a.equals("softmask") || a.equals("soft")){
				softMask=Parse.parseBoolean(b);
				if(softMask) {maskRepeats=true;}
			}else if(a.equals("hardmask") || a.equals("hard")){
				softMask=!Parse.parseBoolean(b);
				if(!softMask) {maskRepeats=true;}
			}else if(a.equals("symbol") || a.equals("masksymbol")){
				softMask=false;
				maskSymbol=(byte)b.charAt(0);
			}else if(a.equals("weak") || a.equals("weaksubsume") || a.equals("weaksubsumes") || a.equals("ignoregaps")){
				weakSubsumes=Parse.parseBoolean(b);
			}
			
			else if(a.equals("maskentropy") || a.equals("entropy") || a.equals("minentropy")){
				if(b==null){
					processEntropy=true;
				}else if(Tools.startsWithLetter(b)){
					processEntropy=Parse.parseBoolean(b);
				}else{
					entropyMaskRatio=Float.parseFloat(b);
					processEntropy=entropyMaskRatio>=0;
				}
			}else if(a.equals("maskentropywindow") || a.equals("entropymaskwindow") || a.equals("ew")){
				entropyWindow=Integer.parseInt(b);
			}else if(a.equals("entropyk") || a.equals("entropymaskk") || a.equals("maskentropyk") || a.equals("ek") || a.equals("ke")){
				entropyK=Integer.parseInt(b);
			}
			
			else if(a.equals("shorttandem") || a.equals("tandem") || a.equals("str")){
				processShortTandem=Parse.parseBoolean(b);
			}else if(a.equals("shorttandemmink") || a.equals("strmink") || a.equals("stmink")){
				shortTandemMinK=Integer.parseInt(b);
			}else if(a.equals("shorttandemmaxk") || a.equals("strmaxk") || a.equals("stmaxk")){
				shortTandemMaxK=Integer.parseInt(b);
			}else if(a.equals("shorttandemmincount") || a.equals("strmincount") || a.equals("stmincount")){
				shortTandemMincount=Integer.parseInt(b);
			}else if(a.equals("shorttandemminlen") || a.equals("strminlen") || a.equals("stminlen")){
				shortTandemMinlen=Integer.parseInt(b);
			}
			
			else if(a.equals("printrepeats") || a.equals("pr") || a.equals("print")){
				printRepeats=Parse.parseBoolean(b);
				maskRepeats=!printRepeats;
			}
			
			else if(a.equals("outs") || a.equals("outm") || a.equals("outmasked")){
				outs=b;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(parser.parseK(arg, a, b)){//Parse k
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
		if(!Tools.testOutputFiles(overwrite, append, false, outr, outs)){
			outstream.println((outr==null)+", "+outr+", "+outs);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+outr+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, outr, outs)){
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
//		assert(false) : "TODO";
		assert(k>=1 && k<=32) : k;
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
		
		AbstractKmerTableSet.DISPLAY_STATS=false;
		
		if(processDepth) {
			/* Fill tables with kmers */
			tables.process(t);

			if(DISPLAY_PROGRESS){
				outstream.println("After loading:");
				Shared.printMemory();
				outstream.println();
			}
		}
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris();
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros();
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads(cris, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Show the stuff the person involved cares about
		printRepeats(ffoutr, masterListOfRepeats);
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		outstream.println(Tools.number("Repeat Bases:       ", repeatBases, 8)+String.format(" \t%.2f%%", repeatBases*100.0/basesProcessed));
		if(processDepth){
			outstream.println(Tools.number("   High Depth:      ", highDepthBases, 8)+String.format(" \t%.2f%%", highDepthBases*100.0/basesProcessed));
		}
		if(processEntropy){
			outstream.println(Tools.number("   Low Entropy:     ", lowEntropyBases, 8)+String.format(" \t%.2f%%", lowEntropyBases*100.0/basesProcessed));
		}
		if(processShortTandem){
			outstream.println(Tools.number("   Short Tandem:    ", shortTandemBases, 8)+String.format(" \t%.2f%%", shortTandemBases*100.0/basesProcessed));
		}
		if(maskRepeats && ros!=null) {
			outstream.println(Tools.number("Bases Masked:       ", basesMasked, 8)+String.format(" \t%.2f%%", basesMasked*100.0/basesProcessed));
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	void printRepeats(FileFormat ff, ArrayList<Repeat> list) {
		Collections.sort(list);
		Collections.reverse(list);
		if(ff==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println(Repeat.tsvHeader());
		ByteBuilder bb=new ByteBuilder();
		for(Repeat r : list) {
			r.appendTo(bb);
			bsw.println(bb);
			bb.clear();
		}
		bsw.poisonAndWait();
	}
	
	private ConcurrentReadInputStream makeCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
		cris.start(); //Start the stream
		if(verbose){
			outstream.println("Started cris");
			boolean paired=cris.paired();
			if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		}
		return cris;
	}
	
	private ConcurrentReadOutputStream makeCros(){
		if(ffouts==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffouts, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, i));
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
			repeatBases+=pt.repeatBasesT;
			basesMasked+=pt.basesMaskedT;
			highDepthBases+=pt.highDepthBasesT;
			lowEntropyBases+=pt.lowEntropyBasesT;
			shortTandemBases+=pt.shortTandemBasesT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			errorState|=(!pt.success);
			masterListOfRepeats.addAll(pt.repeatSet.oldRepeats);
		}
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	final long rcomp(long kmer, int len){
		return amino ? kmer : AminoAcid.reverseComplementBinaryFast(kmer, len);
	}
	
//	final int maskRepeat(Repeat r) {
//		return maskRepeat(r.contig.bases, r.start, r.stop);
//	}
	
	final int maskRange(CRange r) {
		return maskRepeat(((Read)r.obj).bases, r.a, r.b);
	}
	
	final int maskRepeat(byte[] bases, int start, int stop) {
		int masked=0;
		for(int i=start; i<=stop; i++) {
			byte b=bases[i];
			byte b2=baseMaskArray[b];
			bases[i]=b2;
			masked+=(b2==b ? 0 : 1);
		}
		return masked;
	}
	
	/**
	 * Transforms a kmer into all canonical values for a given Hamming distance.
	 * Returns the maximal count stored in the tables.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param qHDist Hamming distance
	 * @return Value stored in table, or -1
	 */
	private final int getValue(final long kmer, final long rkmer, final int qHDist){
		if(verbose){outstream.println("getValue()");}
		final int count=tables.getCount(kmer, rkmer);
		int maxCount=count;
		if(/*count<2 &&*/ qHDist>0){
			final int qHDist2=qHDist-1;

			//Sub
			for(int j=0; j<symbols/* && count<2*/ && maxCount<maxDepth; j++){
				for(int i=0; i<k /*&& id<1*/; i++){
					final long temp=(kmer&clearMasks[i])|setMasks[j][i];
					//					outstream.println(i+", "+j+", "+setMasks[j][i]+", "+qHDist);
					if(temp!=kmer){
						long rtemp=rcomp(temp, k);
						int keyCount=getValue(temp, rtemp, qHDist2);
						maxCount=Tools.max(keyCount, maxCount);
					}
				}
			}
		}
		return maxCount;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_, final int tid_){
			cris=cris_;
			ros=ros_;
			tid=tid_;
			
			repeatSet=new RepeatSet(k, minDepth, maxDepth, minRepeat, maxGap, weakSubsumes, amino, entropyK, entropyWindow);
			et=(processEntropy ? new EntropyTracker(entropyK, entropyWindow, amino, entropyMaskRatio, true) : null);
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing

			synchronized(this) {
				//Process the reads
				processInner();

				//Do anything necessary after processing

				//Indicate successful exit status
				success=true;
			}
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			
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
			
//			int maxDepthSeen=0;
			assert(repeatSet.closedRepeats.isEmpty());
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r=reads.get(idx);
				
				//Validate reads in worker threads
				if(!r.validated()){r.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r.length();
				final int initialLength2=r.mateLength();

				//Increment counters
				readsProcessedT+=r.pairCount();
				basesProcessedT+=initialLength1+initialLength2;
				
				{//Reads are processed in this block.
					if(processDepth) {
//						int maxDepthThisRead=processReadDepth(r);
//						maxDepthSeen=Tools.max(maxDepthSeen, maxDepthThisRead);
						int x=processReadDepth(r);
						highDepthBasesT+=x;
					}
					if(processEntropy){
						int x=processReadEntropy(r);
						lowEntropyBasesT+=x;
					}
					if(processShortTandem){
						int x=processReadShortTandem(r, shortTandemMinK, shortTandemMaxK, shortTandemMincount, shortTandemMinlen);
						shortTandemBasesT+=x;
					}
					repeatSet.recent.clear();
				}
			}
			
//			repeatSet.collectResidual(maxDepthSeen);//Now performed per read
			
			//This is occasionally needed due to lazy collection
			repeatSet.subsumeClosed(weakSubsumes);
			
			ArrayList<CRange> ranges=repeatSet.closedToRanges(true);
			
			//Output reads to the output stream
			if(ros!=null){
				ArrayList<Read> seqList=reads;
				if(maskRepeats) {
					if(!softMask && Repeat.SEQ_AFFIX_LEN>0) {makePreview(repeatSet.closedRepeats);}
					basesMaskedT+=maskRepeats(ranges);
				}else if(printRepeats) {seqList=repeatSet.fetchRepeatSequence();}
				ros.add(seqList, ln.id);
				for(Read r : seqList) {
					readsOutT+=r.pairCount();
					basesOutT+=r.pairLength();
				}
			}
			repeatBasesT+=countRepeatBases(ranges);
			
			//Optionally clear sequence to save memory (but this interferes with printing repeat sequence)
//			for(Repeat r : closedRepeats) {r.contig=null;}
			repeatSet.retireClosed();
		}
		
		long maskRepeats(ArrayList<CRange> ranges) {
			long x=0;
			for(CRange r : ranges) {x+=maskRange(r);}
			return x;
		}
		
		long countRepeatBases(ArrayList<CRange> ranges) {
			long x=0;
			for(CRange r : ranges) {x+=r.length();}
			return x;
		}
		
		void makePreview(ArrayList<Repeat> list) {
			assert(maskRepeats && !softMask);
			if(list.isEmpty()){return;}
			ByteBuilder bb=new ByteBuilder(Repeat.SEQ_AFFIX_LEN*2+3);
			for(Repeat r : list) {r.setSeq(bb);}
		}
		
		//Deprecated and replaced by masking ranges
//		void maskRepeats() {
//			//Masking can be done on closedRepeats, but making a weakSubsumes list can save rework
//			ArrayList<Repeat> list=closedRepeats;
//			if(list.isEmpty()) {return;}
//			if(!softMask) {
//				//If hard masking, copies of the sequence need to be made first.
//				//Skipping this will make the preview come out in lower case,
//				//but it prevents memory crashes in the worst case.
//				//Although it uses more memory in the expected case.
//				ByteBuilder bb=new ByteBuilder(Repeat.SEQ_AFFIX_LEN*2+3);
//				for(Repeat r : list) {r.setSeq(bb);}
//			}
//			if(list.size()>1 && !weakSubsumes) {
//				list=(ArrayList<Repeat>) list.clone();
//				subsume(list, true);
//			}
//			for(Repeat r : list) {
//				maskRepeat(r);
//			}
//		}
		
		/**
		 * Process a read.
		 * @param r1 Read 1
		 * @return max depth observed.
		 */
		private int processReadDepth(final Read r){
			if(r==null || r.length()<k) {return 0;}
			final byte[] bases=r.bases;
			
			long kmer=0;
			long rkmer=0;
			int len=0;
			
			final int start=0;
			final int stop=bases.length;
			int max=0;
			
			assert(repeatSet.recent.isEmpty());
			repeatSet.recent.clear();
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			for(int i=start; i<stop; i++){
				byte b=bases[i];
				long x=symbolToNumber0[b];
				long x2=symbolToComplementNumber0[b];
				kmer=((kmer<<bitsPerBase)|x)&mask;
				rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
				if(forbidNs && !isFullyDefined(b)){len=0; rkmer=0;}else{len++;}
				//if(verbose){outstream.println("Scanning6b i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k){
					final int count=getValue(kmer, rkmer, qHammingDist);
					if(verbose){outstream.println("Testing kmer "+kmer+"; count="+count);}
					//assert(count>=1) : count;
					if(count>=minDepth) {
						repeatSet.increment(r, i, count);
						max=Tools.max(max, count);
					}
				}
			}
			
			repeatSet.collectResidual(max);
			ArrayList<CRange> recent=repeatSet.recentToRanges(true);
			repeatSet.recent.clear();
			
			int sum=0;
			for(CRange range : recent) {
				sum+=range.length();
			}
			return sum;
		}
		
		private int processReadEntropy(final Read rd){
			final int window=et.windowBases();
			if(rd==null || rd.length()<window){return 0;}
			final byte[] bases=rd.bases;
			
			et.clear();
			Repeat current=new Repeat(null, -1, minDepth, window, Tools.max(0, maxGap-k)+window, minRepeat, 'E'); //TODO: Use a buffer.
			int sum=0;
			for(int i=0, min=window-1; i<bases.length; i++){
				et.add(bases[i]);
				if(i>=min && et.ns()<1 && !et.passes()){
					int a=et.leftPos(), b=et.rightPos();
					Repeat old=current.increment(rd, b, 2);
					if(old!=null){
						repeatSet.addRepeat(old);
						sum+=old.length();
					}
				}
			}
			if(current.start>-1 && current.length()>=minRepeat){
				assert(current.contig==rd);
				repeatSet.addRepeat(current);
			}
			return sum;
		}
		
		private int processReadShortTandem(Read rd, int mink, int maxk, int mincount, int minlen){
			final byte[] bases=rd.bases;
			final ArrayList<CRange> ranges=new ArrayList<CRange>();//TODO: Make a buffer

			for(int k=mink; k<=maxk; k++){
				findShortTandem(bases, ranges, k, Tools.max(minlen, k*mincount), rd.numericID);
			}
			if(ranges.isEmpty()) {return 0;}
			CRange.mergeList(ranges, maxk>mink);
			int sum=0;
			for(CRange range : ranges){
				sum+=range.length();
				//Some info is lost by waiting until here to convert to a Repeat
				Repeat repeat=new Repeat(rd, range.a, minDepth, mink, maxGap, minRepeat, 'T');
				repeat.stop=range.b;
				repeat.maxDepth=minDepth;
				repeatSet.addRepeat(repeat);
			}
			return sum;
		}
		
//		assert(currentDepth>=depth);//TODO: This can be disabled, but then the currentDepth<depth case needs to be handled; currently the function is not called inside gaps
//		final int gap=pos-stop-1;
//		if(contigNum==currentContig.numericID && gap<=maxGap) {//advance
////			System.err.println("A:"+this);
//			stop=pos;
//			gapLen+=gap;
//			gapCount+=(gap>0 ? 1 : 0);
//			gapBP+=(gap>=k ? gap-k+1 : 0);
//			depthSum+=currentDepth;
//			maxDepth=Tools.max(currentDepth, maxDepth);
//			return null;
//		}
//		
//		Repeat r=null;
//		if(contigNum>=0 && length()>=minRepeat) {r=this.clone();}
////		System.err.println("B:"+r+", "+this);
//		clear();
//		contigNum=currentContig.numericID;
//		start=pos-k+1;
//		assert(start>=0) : start+", "+pos+", "+k;
//		stop=pos;
//		assert(stop<currentContig.bases.length) : stop+", "+currentContig.length();
//		minDepth=maxDepth=depth;
//		
//		//These are not *strictly* needed and can use a lot of memory.
//		//They could be cleared after calculating entropy and gc or removed entirely.
//		contig=currentContig;
//		contigName=currentContig.name();
//		


		private void findShortTandem(final byte[] bases, final ArrayList<CRange> ranges, final int k, final int minlen, final long rid){
			final int lim=bases.length-k;
			final int mask=(k>15 ? -1 : ~((-1)<<(2*k)));
			CRange last=null;
			for(int loc=0; loc<lim; loc++){
				int len=shortTandemLength(bases, k, mask, loc);
				if(len>=minlen){
					int a=loc-k, b=loc-k+len-1;
					if(last!=null && last.touches(a, b)) {
						last.b=b;
					}else{
						CRange r=new CRange(rid, a, b);
						ranges.add(r);
						loc=Tools.max(loc, b-minlen);
						last=r;
					}
				}else{
					//System.err.println("len="+len+" < minlen="+minlen);
				}
			}
		}

		private int shortTandemLength(final byte[] bases, final int k, final int mask, final int loc){

			final int lim=bases.length;
			final int key=getInitialKey(bases, loc, k);
			if(key<0){return 0;}
			int kmer=key;
			int gap=0, last=-1;
			for(int i=loc; i<lim && gap<k; i++){
				final byte b=bases[i];
				final int n=symbolToNumber[b];
				
				kmer=(((kmer<<2)&mask)|n);
				if(kmer==key){
					last=i;
					gap=0;
				}else if(kmer>=0){
					gap++;
				}else{
					break;//Undefined symbol encountered
				}
			}

			return (last<0 ? 0 : last-loc+k+1);
		}
		
		private int getInitialKey(byte[] bases, int loc, int k){
			assert(k<16);
			int start=loc-k;
			int key=0;
			if(start<0){return -1;}
			for(int i=start; i<loc && key>=0; i++){
				final byte b=bases[i];
				final int n=symbolToNumber[b];
				key=(key<<2)|n;
			}
			return key;
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		/** Number of repeat bases found by this thread */
		protected long repeatBasesT=0;
		
		/** Number of bases masked by this thread */
		protected long basesMaskedT=0;
		
		/** Number of high depth bases found by this thread */
		protected long highDepthBasesT=0;
		
		/** Number of low-entropy bases found by this thread */
		protected long lowEntropyBasesT=0;
		
		/** Number of short tandem repeat bases found by this thread */
		protected long shortTandemBasesT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Thread ID */
		final int tid;
		
		final RepeatSet repeatSet;
		final EntropyTracker et;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String outr="stdout.txt";
	
	/** Sequence output */
	private String outs=null;
	
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

	protected long repeatBases=0;
	protected long basesMasked=0;
	
	/** Number of high-depth bases found */
	protected long highDepthBases=0;
	
	/** Number of low-entropy bases found */
	protected long lowEntropyBases=0;
	
	/** Number of short tandem repeat bases */
	protected long shortTandemBases=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;

	boolean forbidNs=true;

	boolean printRepeats;
	boolean maskRepeats=true;
	boolean softMask=true;
	byte maskSymbol='N';
	boolean weakSubsumes=false;
	
	int minRepeat=0;
	int maxGap=0;

	boolean processDepth=true;
	int maxDepth=Shared.MAX_ARRAY_LEN;
	int minDepth=2;
	
	boolean processEntropy=false;
	int entropyWindow=80;
	int entropyK=5;
	float entropyMaskRatio=0.70f;

	boolean processShortTandem=false;
	int shortTandemMinK=2;
	int shortTandemMaxK=15;
	int shortTandemMincount=4;
	int shortTandemMinlen=32;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	
	/** Primary output file */
	private final FileFormat ffoutr;
	
	/** Sequence output file */
	private final FileFormat ffouts;
	
	private final int k;
	private final KmerTableSet tables;
	
	private final ArrayList<Repeat> masterListOfRepeats=new ArrayList<Repeat>();
	
	private final byte[] baseMaskArray;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*-----------        Symbol-Specific Constants        ----------*/
	/*--------------------------------------------------------------*/

	/** True for amino acid data, false for nucleotide data */
	final boolean amino;
	boolean maskMiddle=false;
	int qHammingDist=0;
	
//	final int maxSupportedK;
	final int bitsPerBase;
	final int maxSymbol;
	final int symbols;
	final int symbolArrayLen;
	final int symbolSpace;
	final long symbolMask;
	
//	final int minlen;
	final int shift;
	final int shift2;
	final long mask;
	final long middleMask;
	
	/** x&clearMasks[i] will clear base i */
	final long[] clearMasks;
	/** x|setMasks[i][j] will set base i to j */
	final long[][] setMasks;
	/** x&leftMasks[i] will clear all bases to the right of i (exclusive) */
	final long[] leftMasks;
	/** x&rightMasks[i] will clear all bases to the left of i (inclusive) */
	final long[] rightMasks;
	
	/** Symbol code; -1 for undefined */
	final byte[] symbolToNumber;
	/** Symbol code; 0 for undefined */
	final byte[] symbolToNumber0;
	/** Complementary symbol code; 0 for undefined */
	final byte[] symbolToComplementNumber0;
	
	/** For verbose / debugging output */
	final String kmerToString(long kmer, int k){
		return amino ? AminoAcid.kmerToStringAA(kmer, k) : AminoAcid.kmerToString(kmer, k);
	}
	
	/** Returns true if the symbol is not degenerate (e.g., 'N') for the alphabet in use. */
	final boolean isFullyDefined(byte symbol){
		return symbol>=0 && symbolToNumber[symbol]>=0;
	}
	
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
