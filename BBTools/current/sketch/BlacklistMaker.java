package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicInteger;

import bloom.KCountArray;
import bloom.KmerCountAbstract;
import bloom.ReadCounter;
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
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.IntList;
import structures.IntListCompressor;
import structures.ListNum;
import structures.LongList;
import tax.AccessionToTaxid;
import tax.GiToTaxid;
import tax.TaxNode;
import tax.TaxTree;
import tracker.EntropyTracker;

/**
 * Makes blacklists in a taxa-aware manner.
 * 
 * @author Brian Bushnell
 * @date April 3, 2017
 *
 */
public class BlacklistMaker extends SketchObject {
	
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
		
		ReadCounter.maxShortKmerLength=32;
		
		//Create an instance of this class
		BlacklistMaker x=new BlacklistMaker(args);
		
		//Run the object
		x.process(t);
		
		ReadCounter.maxShortKmerLength=31;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	@SuppressWarnings("unchecked")
	public BlacklistMaker(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		if(Shared.threads()>=32){Shared.capThreads(Shared.threads()/2);}//Deals with hyperthreading, which causes excessive memory usage
		Shared.capBuffers(Tools.max(4, Shared.threads()+1));
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		KmerCountAbstract.SKETCH_MODE=true;
		KmerCountAbstract.STORE_HASHED=true;
		KmerCountAbstract.KEEP_DUPLICATE_KMERS=true;
		
		//Create a parser object
		Parser parser=new Parser();
		
		int mode_=PER_SEQUENCE;
		hashNames=true;
		
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
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}
			
			
			else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				if(b==null || "ignore".equalsIgnoreCase(b) || "skip".equalsIgnoreCase(b)){
					giTableFile=null;
					TaxTree.CRASH_IF_NO_GI_TABLE=false;
				}else{giTableFile=b;}
			}else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
			}else if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("imgfile") || a.equals("imgdump")){
				imgFile=b;
			}
			
			
			else if(a.equals("mincount") || a.equals("mintaxcount")){
				minTaxCount=Parse.parseIntKMG(b);
			}else if(a.equals("prefilter")){
				prefilter=Parse.parseBoolean(b);
			}else if(a.equals("prepasses") || a.equals("passes")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(b.equalsIgnoreCase("auto")){
					prepasses=2;
					autoPasses=true;
				}else{
					prepasses=Integer.parseInt(b);
					autoPasses=false;
				}
			}else if(a.equals("prehashes") || a.equals("hashes")){
				prehashes=Integer.parseInt(b);
			}else if(a.equals("prebits") || a.equals("bits")){
				prebits=Integer.parseInt(b);
			}else if(a.equals("name")){
				outName=b;
			}else if(a.equalsIgnoreCase("name0") || a.equalsIgnoreCase("nm0")){
				sketchName=b;
			}else if(a.equals("taxid") || a.equals("tid")){
				outTaxid=Integer.parseInt(b);
			}else if(parseMode(arg, a, b)>-1){
				mode_=parseMode(arg, a, b);
			}else if(a.equals("hist")){
				outHist=b;
			}
			
			else if(a.equals("taxlevel") || a.equals("tl") || a.equals("level") || a.equals("lv")){
				if(b==null){taxLevel=-1;}
				else if(Tools.isDigit(b.charAt(0))){
					taxLevel=Integer.parseInt(b);
				}else{
					taxLevel=TaxTree.parseLevel(b);
				}
			}
			
			else if(parseSketchFlags(arg, a, b)){
				//do nothing
			}else if(defaultParams.parse(arg, a, b)){
				//do nothing
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		if("auto".equalsIgnoreCase(imgFile)){imgFile=TaxTree.defaultImgFile();}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		
		mode=mode_;
		assert((mode!=PER_TAXA && mode!=PER_IMG) || taxTreeFile!=null);
		assert(mode!=PER_IMG || imgFile!=null);
		assert(mode==PER_TAXA || mode==PER_SEQUENCE || mode==PER_IMG);
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;
			in2=parser.in2;

			outSketch=parser.out1;
			
			extin=parser.extin;
		}
		
		postParse();
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
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
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, outSketch)){
			outstream.println((outSketch==null)+", "+outSketch);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+outSketch+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, outSketch, outHist)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		if(sketchName==null && outSketch!=null){
			sketchName=outSketch;
			sketchName=ReadWrite.stripToCore(sketchName);
		}
		
		//Create output FileFormat objects
		ffsketch=FileFormat.testOutput(outSketch, FileFormat.SKETCH, null, true, overwrite, append, false);
		ffhist=FileFormat.testOutput(outHist, FileFormat.TXT, null, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
//		tool=new SketchTool(targetSketchSize, defaultParams.minKeyOccuranceCount, false, false);
		
		if(taxTreeFile!=null){setTaxtree(taxTreeFile, outstream);}
		
		if(giTableFile!=null){
			loadGiToTaxid();
		}
		if(accessionFile!=null){
			AccessionToTaxid.tree=taxtree;
			outstream.println("Loading accession table.");
			AccessionToTaxid.load(accessionFile);
			System.gc();
		}
		if(imgFile!=null){
			TaxTree.loadIMG(imgFile, false, outstream);
		}
		
		maps=new HashMap[ways];
		for(int i=0; i<ways; i++){
			maps[i]=new HashMap<Long, IntListCompressor>();
		}
		
		calcMemory();
		
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		if(prefilter){
			makePrefilter();
		}
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads(cris);
		
		prefilterArray=null;
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println("Blacklist size: \t"+resultingSize+"\n");
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, i));
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
			success&=pt.success;
		}
		
		shrinkListsAndWriteHist();
		writeSketch(true);
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void writeSketch(boolean destroy){
		Sketch sk=toSketch(destroy);
		if(ffsketch!=null){errorState|=SketchTool.write(sk, ffsketch);}
	}
	
	private void shrinkListsAndWriteHist(){
		int max=1000000;
		long[] counts=new long[max+1];
		for(int i=0; i<ways; i++){
			for(Entry<Long, IntListCompressor> entry : maps[i].entrySet()){
				IntList value=entry.getValue().list;
				value.sort();
				value.shrinkToUnique();
				int index=Tools.min(max, value.size);
				counts[index]++;
			}
		}
		if(ffhist!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(ffhist);
			bsw.start();
			bsw.print("#count\tkmers\n".getBytes());
			for(int i=0; i<counts.length; i++){
				long count=counts[i];
				if(count>0){
					bsw.print(i);
					bsw.print('\t');
					bsw.print(count);
					bsw.print('\n');
				}
			}
			errorState|=bsw.poisonAndWait();
		}
	}

	private Sketch toSketch(boolean destroy){
		long[] array=toArray(destroy);
		hashArrayToSketchArray(array);
		ArrayList<String> meta=new ArrayList<String>();
		meta.add("minTaxCount:"+minTaxCount);
		meta.add("taxLevel:"+taxLevel);
		Sketch sk=new Sketch(array, null, null, null, null, outTaxid, -1, -1, -1, -1, -1, "blacklist", sketchName, ffin1.simpleName(), meta);
		return sk;
	}
	
	private long[] toArray(boolean destroy){
		LongList list=new LongList();
		for(int i=0; i<ways; i++){
			for(Entry<Long, IntListCompressor> entry : maps[i].entrySet()){
				Long key=entry.getKey();
				IntList value=entry.getValue().list;
				if(value.size()>=minTaxCount){
					list.add(key.longValue());
				}
			}
			if(destroy){maps[i]=null;}
		}
		resultingSize=list.size();
		return list.toArray();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Prefilter          ----------------*/
	/*--------------------------------------------------------------*/

	private void calcMemory(){
		long usableMemory=Shared.memAvailableAdvanced();
		
//		System.err.println(memory/1000000+", "+usableMemory/1000000+", "+usableMemory/1000000+", "+prepasses);
		
		if(prepasses==0 || !prefilter){
			filterMemory0=filterMemory1=0;
		}else{
			double low=Tools.min(prefilterFraction, 1-prefilterFraction);
			double high=(autoPasses ? low : 1-low);
			if((prepasses&1)==1){//odd passes
				filterMemory0=(long)(usableMemory*low);
				filterMemory1=(long)(usableMemory*high);
			}else{//even passes
				filterMemory0=(long)(usableMemory*high);
				filterMemory1=(long)(usableMemory*low);
			}
		}
	}
	
	private long filterMemory(int pass){return ((pass&1)==0) ? filterMemory0 : filterMemory1;}
	
	private void makePrefilter(){
		prefilterArray=makePrefilter_inner(new KCountArray[1], 0, minTaxCount);
		if(prefilterArray!=null){
			prefilterArray.purgeFilter();
//			filterMaxValue=Tools.min(filterMax, prefilterArray.maxValue-1);
		}
	}
	
	public final KCountArray makePrefilter_inner(final KCountArray[] filter, int currentPass, int overallFilterMax){
//		assert(false) : lastFilter+", "+prefilter+", "+filterMax()+", "+currentPass+", "+filterMemory(currentPass);
		if(!prefilter){return null;}
		
		if(filter[0]!=null){
			filter[0].purgeFilter();
			assert(filter[0].prefilter()==null);
		}
		
		KmerCountAbstract.CANONICAL=true;

		long precells=-1;
		int cbits=2;
		int filterMax=3;
		
		if(currentPass==0 && prebits>0){
			cbits=prebits;
			assert(Integer.bitCount(prebits)==1 && prebits<=32);
			filterMax=Tools.min((1<<cbits)-1, overallFilterMax);
		}else if(currentPass>0 || prepasses<2){
			filterMax=overallFilterMax;
			while(filterMax>=(1<<cbits)){cbits*=2;}
		}
		
		byte minq=0;
		if(precells<1){
			long prebits=(filterMemory(currentPass)-10)*8;
			
//			System.err.println("prebits="+prebits+", currentPass="+currentPass);
			
			precells=prebits/cbits;
			if(precells<100000){ //Not enough memory - no point.
				prefilter=false;
				return null;
			}
		}
		if(prehashes<1){prehashes=2;}

		{
			Timer ht=new Timer();

			ArrayList<String> extra=null;
			ReadCounter rc=new ReadCounter(k, rcomp, false, false, amino);
			filter[0]=rc.makeKca(in1, in2, extra, cbits, precells, prehashes, minq, maxReads, 1, 1, 1, filter[0], filterMax);
			assert(filterMax<=filter[0].maxValue || (currentPass==0)) : filterMax+", "+currentPass+", "+filter[0].maxValue;
			outstream.println("Made prefilter:   \t"+filter[0].toShortString(prehashes));
			double uf=filter[0].usedFraction();
//			System.err.println("cellsUsed: "+filter[0].cellsUsed(1)+" //123"); //123
			if(uf>0.5){
				outstream.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.95 ? "incredibly" : uf>0.9 ? "extremely" : uf>0.8 ? "very" :
					uf>0.7 ? "rather" : uf>0.6 ? "fairly" : "somewhat")+" full.  Ideal load is under 50% used." +
						"\nFor better accuracy, run on a node with more memory; quality-trim or error-correct reads; or increase prefiltersize.");
			}
			ht.stop();
			
			final double kmers=filter[0].estimateUniqueKmers(prehashes, Tools.min(overallFilterMax, filter[0].maxValue));
			outstream.println("Estimated valid kmers: \t\t"+(long)kmers);
			
//			outstream.println("Estimated valid kmers 1+: "+(long)filter[0].estimateUniqueKmers(prehashes, 1));
//			outstream.println("Estimated valid kmers 2+: "+(long)filter[0].estimateUniqueKmers(prehashes, 2));
//			outstream.println("Estimated valid kmers 3+: "+(long)filter[0].estimateUniqueKmers(prehashes, 3));
//			outstream.println("Estimated valid kmers 4+: "+(long)filter[0].estimateUniqueKmers(prehashes, 4));
			
			if(autoPasses && kmers<1000000){
				prepasses=currentPass;
			}
			
			if(currentPass+1<prepasses){
				return makePrefilter_inner(filter, currentPass+1, overallFilterMax);
			}
			
//			if(DISPLAY_PROGRESS){
				outstream.println("Prefilter time:\t"+ht);
				outstream.println("After prefilter:");
				Shared.printMemory();
				outstream.println();
//			}
		}
		
		return filter[0];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Tax Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private void loadGiToTaxid(){
		Timer t=new Timer();
		outstream.println("Loading gi to taxa translation table.");
		GiToTaxid.initialize(giTableFile);
		t.stop();
		if(true){
			outstream.println("Time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final int tid_){
			cris=cris_;
			tid=tid_;
			
			if(defaultParams.minEntropy>0){
				eTracker=new EntropyTracker(entropyK, k, amino, defaultParams.minEntropy, true);
			}else{
				eTracker=null;
			}
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
			assert(!translate); //Translate first
			
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
					
					if((initialLength1>=k || initialLength2>=k) && (!tossJunk || !r1.junk())){
						
						int taxID=-1;
						
						if(mode==PER_SEQUENCE){
							taxID=(int)(r1.numericID & ((long)Integer.MAX_VALUE));
						}else{
							assert(mode==PER_TAXA || mode==PER_IMG) : mode;
							TaxNode tn=null;
							tn=taxtree.parseNodeFromHeader(r1.id, bestEffort);
							while(tn!=null && tn.pid!=tn.id && tn.level<taxLevel){
								TaxNode temp=taxtree.getNode(tn.pid);
								if(temp==null || temp==tn || temp.level>=TaxTree.LIFE || temp.level>taxLevel){break;}
								tn=temp;
							}
							if(tn!=null){taxID=tn.id;}
							if(taxID<0){
								taxID=nextUnknown.getAndIncrement();
							}
						}
						assert(taxID>0) : r1.id+", "+mode;
//						System.err.println("TID="+taxID);

						//Process reads.
						if(amino){
							processReadAmino(r1, taxID);
							if(r2!=null){processReadAmino(r2, taxID);}
						}else{
							processReadNucleotide(r1, taxID);
							if(r2!=null){processReadNucleotide(r2, taxID);}
						}
					}
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
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		/**
		 * Process a read.
		 * @param r Read
		 */
		void processReadNucleotide(final Read r, final int value){
			assert(k<=32);

			final byte[] bases=r.bases;

			if(bases==null || bases.length<k){return;}

			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
			long kmer=0;
			long rkmer=0;
			int len=0;
			if(eTracker!=null){eTracker.clear();}

			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(eTracker!=null){eTracker.add(b);}

				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{
					len++;
					if(len>=k){
						if(eTracker==null || eTracker.passes()){
							final long hashcode=hash(kmer, rkmer);
//							System.err.println("kmer="+kmer+", rkmer="+rkmer+", mask="+mask+", hashcode="+hashcode+", count="+prefilterArray.read(hashcode));
							if(hashcode>=minHashValue){
								int precount=prefilter ? prefilterArray.read(hashcode) : minTaxCount;
								assert(!prefilter || precount>0) : "kmer="+kmer+", rkmer="+rkmer+", mask="+mask+", hashcode="+hashcode+", count="+precount;
								//							System.err.println("Sufficient hash for "+key+", "+value);
								if(precount>=minTaxCount){
									//								if(prefilter){System.err.println("Passed prefilter.");}
									addToMap(hashcode, value);
									keysAddedT++;
								}
							}
						}
					}
				}
			}
		}
		
		/**
		 * Process a read.
		 * @param r Read
		 */
		void processReadAmino(final Read r, final int value){
			final int aminoShift=AminoAcid.AMINO_SHIFT;
			assert(k*aminoShift<=64);
			
			final byte[] bases=r.bases;
			long kmer=0;
			int len=0;
			assert(r.aminoacid());
			
			final int shift=aminoShift*k;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			final long min=minHashValue;
			
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=AminoAcid.acidToNumberNoStops[b];
				kmer=((kmer<<aminoShift)|x)&mask;

//				System.err.println((char)b+" -> "+x);
				
				if(x<0){
					len=0;
					kmer=0;
//					System.err.println("Reset.");
				}else{
					len++;
//					System.err.println("Incremented length to "+len);
					if(len>=k){
						final long hashcode=hash(kmer, kmer);
						if(hashcode>=minHashValue){
//							System.err.println("Sufficient hash for "+key+", "+value);
							if(!prefilter || prefilterArray.read(hashcode)>=minTaxCount){
//								if(prefilter){System.err.println("Passed prefilter.");}
								addToMap(hashcode, value);
								keysAddedT++;
							}
						}
					}
				}
			}
//			System.err.println(keysAddedT);
//			assert(false);
		}
		
		
		void addToMap(long key0, int value){
			Long key=Long.valueOf(key0);
//			Long key=new Long(key0);
			HashMap<Long, IntListCompressor> map=maps[(int)(key%ways)];
			IntListCompressor lh=map.get(key);
			if(lh==null){
				synchronized(map){
					lh=map.get(key);
					if(lh==null){
						lh=new IntListCompressor();
						map.put(key, lh);
					}
				}
			}
			synchronized(lh){
				lh.add(value);
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		protected long keysAddedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Thread ID */
		final int tid;
		
		private final EntropyTracker eTracker;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final int mode;
	
	private String giTableFile=null;
	private String taxTreeFile=null;
	private String accessionFile=null;
	private String imgFile=null;

	private String outName=null;
	private String sketchName=null;
	private int outTaxid=-1;
	
	private int taxLevel=1;
	private boolean prefilter=true;
	private boolean tossJunk=true;
	private boolean bestEffort=true;
	private int minTaxCount=100;
	
	private int prepasses=2;
	private int prehashes=2;
	private int prebits=-1;
	private boolean autoPasses=false;
	
	double prefilterFraction=0.2;
	
	long filterMemory0;
	long filterMemory1;
	
	private HashMap<Long, IntListCompressor>[] maps;
	
	public KCountArray prefilterArray=null;
	
	final int ways=63;
	
	int resultingSize=-1;
	
	private final AtomicInteger nextUnknown=new AtomicInteger(SketchObject.minFakeID);
	
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;

	/** Primary output file path */
	private String outSketch=null;

	/** Histogram output file path */
	private String outHist=null;
	
	/** Override input file extension */
	private String extin=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

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
	private final FileFormat ffsketch;
	/** Histogram output file */
	private final FileFormat ffhist;
	
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
