package bloom;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import assemble.ErrorTracker;
import cardinality.CardinalityTracker;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import shared.KillSwitch;
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
import stream.SamLine;
import structures.IntList;
import structures.ListNum;
import structures.LongList;
import tracker.ReadStats;

/**
 * Wraps a BloomFilter to filter or correct reads.
 * 
 * @author Brian Bushnell
 * @date May 14, 2018
 *
 */
public class BloomFilterCorrectorWrapper {
	
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
		BloomFilterCorrectorWrapper x=new BloomFilterCorrectorWrapper(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public BloomFilterCorrectorWrapper(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Tools.max(Shared.threads()>1 ? 2 : 1, Shared.threads()>20 ? Shared.threads()/2 : Shared.threads()));
		BBMerge.strict=true;
		KCountArray.LOCKED_INCREMENT=true;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.loglog=parser.loglogOut=true;
		
		boolean setBits=false;
		int k_=31;
		int ksmall_=-1;
		int hashes_=3;
		int bits_=2;
		int minCount_=0;
		boolean rcomp_=true;
		boolean requireBothToPass_=true;
		boolean ecc_=true;
		boolean ecco_=false;
		boolean merge_=true;
		boolean testMerge_=true;
		boolean tossjunk_=false;
		boolean vstrict_=true;
		boolean ustrict_=false;
		float highCountFraction_=1.0f;
		corrector=new BloomFilterCorrector(null, k_, k_);
		
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
			}
			
			else if(a.equalsIgnoreCase("k") || a.equalsIgnoreCase("bloomK") || a.equalsIgnoreCase("bloomFilterK") || a.equalsIgnoreCase("kbig")){
				k_=Integer.parseInt(b);
				assert(k_>0 && k_<=31) : "K must be between 1 and 31, inclusive.";
			}else if(a.equalsIgnoreCase("ksmall") || a.equalsIgnoreCase("bloomKsmall") || a.equalsIgnoreCase("bloomFilterKsmall")){
				ksmall_=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("hashes") || a.equalsIgnoreCase("bloomHashes") || a.equalsIgnoreCase("bloomFilterHashes")){
				hashes_=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp_=Parse.parseBoolean(b);
			}else if(a.equals("bits")){
				setBits=true;
				bits_=Integer.parseInt(b);
				assert(bits_>0);
			}else if(a.equals("mincount")){
				minCount_=Integer.parseInt(b);
			}else if(a.equals("minprob")){
				ReadCounter.minProb=Float.parseFloat(b);
			}else if(a.equals("maxload")){
				maxLoad=Float.parseFloat(b);
			}else if(a.equals("requireboth")){
				requireBothToPass_=Parse.parseBoolean(b);
			}else if(a.equals("ecc")){
				ecc_=Parse.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco_=Parse.parseBoolean(b);
			}else if(a.equals("merge")){
				merge_=Parse.parseBoolean(b);
			}else if(a.equals("testmerge")){
				testMerge_=Parse.parseBoolean(b);
			}else if(a.equals("testmergewidth")){
				testMergeWidth=Integer.parseInt(b);
			}else if(a.equals("testmergethresh")){
				testMergeThresh=Integer.parseInt(b);
			}else if(a.equals("testmergemult")){
				testMergeMult=Parse.parseKMG(b);
			}else if(a.equals("vstrict")){
				vstrict_=Parse.parseBoolean(b);
			}else if(a.equals("ustrict")){
				ustrict_=Parse.parseBoolean(b);
			}else if(a.equals("tossjunk")){
				tossjunk_=Parse.parseBoolean(b);
			}else if(a.equals("memfraction") || a.equals("memmult") || a.equals("memratio")){
				memFraction=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("highCountFraction") || a.equalsIgnoreCase("mincountfraction") || a.equals("hcf")){
				highCountFraction_=Float.parseFloat(b);
			}
			
			else if(a.equals("pincer")){
				corrector.ECC_PINCER=Parse.parseBoolean(b);
			}else if(a.equals("tail")){
				corrector.ECC_TAIL=Parse.parseBoolean(b);
			}else if(a.equals("reassemble")){
				corrector.ECC_REASSEMBLE=Parse.parseBoolean(b);
			}else if(a.equals("smooth")){
				if(b!=null && Character.isDigit(b.charAt(0))){
					corrector.smoothWidth=Integer.parseInt(b);
					corrector.smooth=corrector.smoothWidth>0;
				}else{
					corrector.smooth=Parse.parseBoolean(b);
				}
//				assert(false) : corrector.smooth+", "+corrector.smoothWidth;
			}else if(a.equals("smoothwidth")){
				corrector.smoothWidth=Integer.parseInt(b);
			}else if(a.equals("cells")){
				BloomFilter.OVERRIDE_CELLS=Parse.parseKMG(b);
			}else if(a.equals("seed")){
				KCountArray7MTA.setSeed(Parse.parseKMG(b));
			}
			
			else if(a.equals("ref")){
				addFiles(b, ref);
			}else if(a.equals("extra")){
				addFiles(b, extra);
			}else if(a.equals("outm") || a.equals("outm1") || a.equals("out") || a.equals("out1")){
				out1=b;
			}else if(a.equals("outm2") || a.equals("out2")){
				out2=b;
			}else if(a.equals("outb") || a.equals("outb1") || a.equals("outbad") || a.equals("outbad1") || a.equals("outlow") || a.equals("outlow1")){
				outbad1=b;
			}else if(a.equals("outb2") || a.equals("outbad2") || a.equals("outlow2")){
				outbad2=b;
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
		
		while(minCount_>0 && (1L<<bits_)-1<minCount_){bits_*=2;}
		if(!setBits && ecc_ && bits_<4){bits_=4;}
		
		if(ksmall_<=0){ksmall_=k_;}
		assert(ksmall_<=k_) : k_+", "+ksmall_;
		
		k=k_;
		ksmall=Tools.min(k, ksmall_);
		corrector.k=k_;
		corrector.ksmall=ksmall_;
		bits=bits_;
		hashes=hashes_;
		minCount=minCount_;
		rcomp=rcomp_;
		requireBothToPass=requireBothToPass_;
		ecc=ecc_;
		ecco=ecco_;
		merge=merge_;
		testMerge=testMerge_;
		vstrict=vstrict_;
		ustrict=ustrict_;
		tossjunk=tossjunk_;
		highCountFraction=highCountFraction_;
		KmerCountAbstract.CANONICAL=rcomp; //rcomp, or true, perhaps...  hmmm.
		outstream.println("Using "+bits+" bits per cell.");
		junkWidth=Tools.max(corrector.smoothWidth, 1);
		
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
			
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
			
//			parser.loglogk=k;
			parser.loglogMinprob=ReadCounter.minProb;
			loglogOut=((parser.loglog&parser.loglogOut) ? CardinalityTracker.makeTracker(parser) : null);
		}
		
		SamLine.SET_FROM_OK=!ecc;
		
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
		
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		qfin1=Tools.fixExtension(qfin1);
		qfin2=Tools.fixExtension(qfin2);
		ref=Tools.fixExtension(ref);
		extra=Tools.fixExtension(extra);
		
		//Do output file # replacement
		if(outbad1!=null && outbad2==null && outbad1.indexOf('#')>-1){
			outbad2=outbad1.replace("#", "2");
			outbad1=outbad1.replace("#", "1");
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
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outbad1, outbad2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outbad1, outbad2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		
		//Create output FileFormat objects
		ffoutm1=FileFormat.testOutput(outbad1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffoutm2=FileFormat.testOutput(outbad2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

		{
			Timer t=new Timer(outstream, true);
			if(ref.isEmpty()){
				filter=new BloomFilter(in1, in2, extra, ksmall, k, bits, hashes, 1,
						rcomp, ecco, merge, memFraction);
			}else{
				ref.addAll(extra);
				filter=new BloomFilter(null, null, ref, ksmall, k, bits, hashes, 1,
						rcomp, ecco, merge, memFraction);
			}
			t.stop("Filter creation: \t\t");
			outstream.println(filter.filter.toShortString());
		}
		
		if(ecc){
			corrector.filter=filter;
		}
		
		{
			double a=filter.filter.estimateUniqueKmers(hashes);
			outstream.println("Estimated kmers of depth 1+: \t"+(long)a);
			if(bits>1){
				double usedFraction2=filter.filter.usedFraction(2);
				double b=filter.filter.estimateUniqueKmersFromUsedFraction(hashes, usedFraction2);
				outstream.println("Estimated kmers of depth 2+: \t"+(long)b);
				outstream.println("Used fraction for depth 2+:  \t"+Tools.format("%.3f%%", usedFraction2*100));
				if(usedFraction2>maxLoad){
					KillSwitch.kill("Max load exceeded, quitting: "+usedFraction2+" > "+maxLoad);
				}
			}else if(maxLoad<1.0f){
				double usedFraction=filter.filter.usedFraction(2);
				if(usedFraction>maxLoad){
					KillSwitch.kill("Max load exceeded, quitting: "+usedFraction+" > "+maxLoad);
				}
			}
		}
	}
	
	private static void addFiles(String b, ArrayList<String> list){
		if(b==null){list.clear();}
		else{
			if(new File(b).exists()){list.add(b);}
			else{
				for(String s : b.split(",")){list.add(s);}
			}
		}
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
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Optionally create read output streams
		final ConcurrentReadOutputStream ros, rosb;
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
		
		if(ffoutm1!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			rosb=ConcurrentReadOutputStream.getStream(ffoutm1, ffoutm2, qfoutbad1, qfoutbad2, buff, null, false);
			rosb.start(); //Start the stream
		}else{rosb=null;}
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		Timer t2=new Timer(outstream, true);
		
		//Process the reads in separate threads
		spawnThreads(cris, ros, rosb);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros, rosb);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		t2.stop("\nFiltering Time:  \t\t");
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, true));
		if(loglogOut!=null){
			outstream.println("Unique "+loglogOut.k+"-mers out:     \t"+loglogOut.cardinality());
		}
		
		if(ecc){
			final long corrected=(basesCorrectedTail+basesCorrectedPincer+basesCorrectedReassemble+basesCorrectedEcco);
			final long partial=(readsCorrected-readsFullyCorrected);
			outstream.println();
			outstream.println("Errors detected:            \t"+(basesDetected+basesCorrectedEcco));
			{
				StringBuilder sb=new StringBuilder();
				sb.append("Errors corrected:           \t"+Tools.padRight(corrected, 7)+" \t(");
				
				String comma="";
				if(corrector.ECC_PINCER){
					sb.append(comma).append(basesCorrectedPincer+" pincer");
					comma=", ";
				}
				if(corrector.ECC_TAIL || corrector.ECC_ALL){
					sb.append(comma).append(basesCorrectedTail+" tail");
					comma=", ";
				}
				if(corrector.ECC_REASSEMBLE){
					sb.append(comma).append(basesCorrectedReassemble+" reassemble");
					comma=", ";
				}
				if(ecco || merge){
					sb.append(comma).append(basesCorrectedEcco+" overlap");
					comma=", ";
				}
				
				sb.append(")");
				outstream.println(sb);
			}
			if(ecco || merge){outstream.println("Reads merged:               \t"+Tools.padRight(readsMerged, 7)+
					Tools.format(" \t(%.2f%%)", readsMerged*200.0/readsProcessed));}
			outstream.println("Reads with errors detected: \t"+Tools.padRight(readsDetected, 7)+
					Tools.format(" \t(%.2f%%)", readsDetected*100.0/readsProcessed));
			outstream.println("Reads fully corrected:      \t"+Tools.padRight(readsFullyCorrected, 7)+
					Tools.format(" \t(%.2f%% of detected)", readsFullyCorrected*100.0/readsDetected));
			outstream.println("Reads partly corrected:     \t"+Tools.padRight(partial, 7)+
					Tools.format(" \t(%.2f%% of detected)", partial*100.0/readsDetected));
			if(corrector.ECC_ROLLBACK || rollbacks>0){
				outstream.println("Rollbacks:                  \t"+Tools.padRight(rollbacks, 7)+
						Tools.format(" \t(%.2f%% of detected)", rollbacks*100.0/readsDetected));
			}
			
		}
		
		MetadataWriter.write(null, readsProcessed, basesProcessed, readsOut, basesOut, false);
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosb){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, rosb, i));
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
			
			readsExtended+=pt.readsExtendedT;
			basesExtended+=pt.basesExtendedT;
			readsCorrected+=pt.readsCorrectedT;
			basesCorrectedPincer+=pt.basesCorrectedPincerT;
			basesCorrectedTail+=pt.basesCorrectedTailT;
			basesCorrectedReassemble+=pt.basesCorrectedReassembleT;
			readsFullyCorrected+=pt.readsFullyCorrectedT;
			rollbacks+=pt.rollbacksT;
			readsDetected+=pt.readsDetectedT;
			basesDetected+=pt.basesDetectedT;
			readsMarked+=pt.readsMarkedT;
			basesMarked+=pt.basesMarkedT;

			readsMerged+=pt.readsMergedT;
			readsCorrectedEcco+=pt.readsCorrectedEccoT;
			basesCorrectedEcco+=pt.basesCorrectedEccoT;
			
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	final int findOverlap(Read r1, Read r2, boolean ecc){
		if(ustrict){
			return BBMerge.findOverlapUStrict(r1, r2, ecc);
		}else if(vstrict){
			return BBMerge.findOverlapVStrict(r1, r2, ecc);
		}else{
			return BBMerge.findOverlapStrict(r1, r2, ecc);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_, final ConcurrentReadOutputStream rosb_, final int tid_){
			cris=cris_;
			ros=ros_;
			rosb=rosb_;
			tid=tid_;
		}
		
		//Called by start()
		@Override
		public void run(){
			if(ecc){
				corrector.initializeThreadLocals();
				localTracker=corrector.localTracker.get();
				kmers=corrector.localLongList.get();
				counts=corrector.localIntList.get();
			}
			
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
				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				ArrayList<Read> keepList=new ArrayList<Read>(reads.size());
				ArrayList<Read> tossList=new ArrayList<Read>(reads.size());

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					Read r1=reads.get(idx);
					Read r2=r1.mate;
					final String r2id=r1.mateId();
					
					//Validate reads in worker threads
					if(!r1.validated()){r1.validate(true);}
					if(r2!=null && !r2.validated()){r2.validate(true);}

					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=r1.mateLength();

					//Increment counters
					readsProcessedT+=r1.pairCount();
					basesProcessedT+=initialLength1+initialLength2;

					final Read r1_0=r1, r2_0=r2;
					if(ecc){
						if(r2!=null && (merge || ecco)){
							final int insert=findOverlap(r1, r2, false);
							if(merge){
								if(insert>0){
									r2.reverseComplement();
									r1=r1.joinRead(insert);
									r2.reverseComplement();
									r2=null;
									if(testMerge && !corrector.mergeOK(r1, initialLength1, initialLength2, kmers, testMergeWidth, testMergeThresh, testMergeMult)){
										r1=r1_0;
										r2=r2_0;
									}else{
										r2_0.reverseComplement();
										int errors=BBMerge.countErrors(r1_0, r2_0, r1);
										r2_0.reverseComplement();
										basesCorrectedEccoT+=errors;
										readsCorrectedEccoT+=(errors>0 ? 1 : 0);
										readsMergedT++;
									}
								}
							}else if(ecco){
//								findOverlap(r1, r2, true);
								if(insert>0){
									r2.reverseComplement();
									Read merged=r1.joinRead(insert);
									if(!testMerge || corrector.mergeOK(merged, initialLength1, initialLength2, kmers, testMergeWidth, testMergeThresh, testMergeMult)){
										int errors=BBMerge.errorCorrectWithInsert(r1, r2, insert);
										basesCorrectedEccoT+=errors;
										readsCorrectedEccoT+=(errors>0 ?1 : 0);
										readsMergedT++;
									}
									r2.reverseComplement();
								}
							}
						}
						
						errorCorrect(r1);
						errorCorrect(r2);
						if(merge && r2==null && r2id!=null){
							if(!localTracker.rollback){
								//unmerge
								final int to=r1.length()-1;
								final int len=Tools.min(r1.length(), initialLength2);
								r2=r1.subRead(to-len+1, to);
								r2.setPairnum(1);
								r2.reverseComplement();
								r2.mate=r1;
								r1.mate=r2;
								r2.id=r2id;
								
								if(r1.length()>initialLength1){
									r1.bases=Arrays.copyOf(r1.bases, initialLength1);
									if(r1.quality!=null){r1.quality=Arrays.copyOf(r1.quality, initialLength1);}
								}
							}else{
								r1=r1_0;
								r2=r2_0;
							}
						}
					}

					boolean keep=true;
					if(minCount>0){
						if(r2!=null && !requireBothToPass){
							keep=filter.hasHighCountFraction(r1, minCount, highCountFraction) || filter.hasHighCountFraction(r2, minCount, highCountFraction);
						}else{
							boolean keep1=filter.hasHighCountFraction(r1, minCount, highCountFraction);
							boolean keep2=(r2==null || !keep1 ? keep1 : filter.hasHighCountFraction(r2, minCount, highCountFraction));
							keep=keep1 && keep2;
						}
					}
					if(tossjunk && keep){
						keep=!filter.isJunk(r1, r2, junkWidth);
					}
					
					assert(r2!=null || r2_0==null);
					assert(r1.mate==r2);
					
					if(keep){
						if(loglogOut!=null){loglogOut.hash(r1);}
						readsOutT+=r1.pairCount();
						basesOutT+=r1.pairLength();
						keepList.add(r1);
					}else{
						tossList.add(r1);
					}
				}

				//Output reads to the output streams
				if(ros!=null){ros.add(keepList, ln.id);}
				if(rosb!=null){rosb.add(tossList, ln.id);}

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
		
		void errorCorrect(Read r){
			if(r==null){return;}
			int corrected=corrector.errorCorrect(r);
			final int detected=localTracker.detected();
			final int correctedPincer=localTracker.correctedPincer;
			final int correctedTail=localTracker.correctedTail;
			final int correctedReassemble=localTracker.correctedReassemble();
			final int marked=localTracker.marked;
			assert(corrected==correctedPincer+correctedTail+correctedReassemble) : corrected+", "+localTracker;
			if(marked>0){
				readsMarkedT++;
				basesMarkedT+=marked;
			}
			if(localTracker.rollback){rollbacksT++;}
			if(detected>0){
				readsDetectedT++;
				basesDetectedT+=detected;
				if(corrected>0){
					readsCorrectedT++;
					basesCorrectedPincerT+=correctedPincer;
					basesCorrectedTailT+=correctedTail;
					basesCorrectedReassembleT+=correctedReassemble;
				}
				if(corrected==detected){
					readsFullyCorrectedT++;
				}
				
//				else if(discardUncorrectable){
//					r.setDiscarded(true);
//					if(r.mate!=null && !requireBothBad){
//						r.mate.setDiscarded(true);
//						return;
//					}else if(r.mate==null){return;}
//				}
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
		/** Matched output stream */
		private final ConcurrentReadOutputStream rosb;
		/** Thread ID */
		final int tid;
		
		long readsExtendedT=0;
		long basesExtendedT=0;
		long readsCorrectedT=0;
		long basesCorrectedPincerT=0;
		long basesCorrectedTailT=0;
		long basesCorrectedReassembleT=0;
		long readsFullyCorrectedT=0;
		long rollbacksT=0;
		long readsDetectedT=0;
		long basesDetectedT=0;
		long readsMarkedT=0;
		long basesMarkedT=0;

		long readsMergedT=0;
		long readsCorrectedEccoT=0;
		long basesCorrectedEccoT=0; //These numbers are not necessarily correct in the case of rollbacks

		ErrorTracker localTracker;
		LongList kmers;
		IntList counts;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<String> ref=new ArrayList<String>();
	private ArrayList<String> extra=new ArrayList<String>();
	
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

	/** Output file path for bad reads */
	private String outbad1=null;
	/** Secondary output file path for bad reads */
	private String outbad2=null;

	private String qfoutbad1=null;
	private String qfoutbad2=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** For calculating kmer cardinality in output */
	final CardinalityTracker loglogOut;
	
	/*--------------------------------------------------------------*/
	
	long readsExtended=0;
	long basesExtended=0;
	long readsCorrected=0;
	long basesCorrectedPincer=0;
	long basesCorrectedTail=0;
	long basesCorrectedReassemble=0;
	long readsFullyCorrected=0;
	long rollbacks=0;
	long readsDetected=0;
	long basesDetected=0;
	long readsMarked=0;
	long basesMarked=0;
	
	long readsMerged=0;
	long readsCorrectedEcco=0;
	long basesCorrectedEcco=0;
	
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
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	final FileFormat ffin1;
	/** Secondary input file */
	final FileFormat ffin2;
	
	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	
	/** Primary output file for matching reads */
	private final FileFormat ffoutm1;
	/** Secondary output file for matching reads */
	private final FileFormat ffoutm2;
	
	final BloomFilter filter;
	
	final BloomFilterCorrector corrector;

	final int k;
	final int ksmall;
	final int hashes;
	final int bits;
	final boolean rcomp;
	final boolean requireBothToPass;
	final boolean ecc;
	final boolean ecco;
	final boolean merge;
	final boolean testMerge;
	final boolean tossjunk;
	final int minCount;
	final float highCountFraction;
	final boolean vstrict;
	final boolean ustrict;
	final int junkWidth;
	float memFraction=1.0f;
	float maxLoad=1.0f; //Crash if this is exceeded
	
	int testMergeWidth=4;
	long testMergeMult=80L;
	int testMergeThresh=3;
	
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
