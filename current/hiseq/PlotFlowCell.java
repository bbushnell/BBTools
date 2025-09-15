package hiseq;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import aligner.MicroAligner3;
import aligner.MicroIndex3;
import barcode.BarcodeStats;
import bloom.BloomFilter;
import bloom.BloomFilterCorrector;
import bloom.KmerCountAbstract;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.Dedupe;
import kmer.AbstractKmerTable;
import kmer.HashArray1D;
import kmer.ScheduleMaker;
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
import structures.ListNum;
import structures.LongList;
import template.Accumulator;
import template.ThreadWaiter;

//Example command:  plotflowcell.sh in=lane1.fq.gz size=999999 dump=dump.txt expected=expected.txt bloom multithreaded -Xmx40g

/**
 * Analyzes a flow cell for low-quality areas.
 * Removes reads in the low-quality areas.
 * 
 * @author Brian Bushnell
 * @date August 31, 2016
 *
 */
public class PlotFlowCell implements Accumulator<PlotFlowCell.WorkerThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		PlotFlowCell x=new PlotFlowCell(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public PlotFlowCell(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		//Create a parser object
		Parser parser=new Parser();
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("divisor") || a.equals("size")){
				Tile.xSize=Tile.ySize=Parse.parseIntKMG(b);
			}else if(a.equals("xdivisor") || a.equals("xsize")){
				Tile.xSize=Parse.parseIntKMG(b);
			}else if(a.equals("ydivisor") || a.equals("ysize")){
				Tile.ySize=Parse.parseIntKMG(b);
			}else if(a.equals("target")){
				targetAverageReads=Parse.parseIntKMG(b);
			}else if(a.equals("dump") || a.equals("out")){
				dump=b;
			}else if(a.equals("indump") || a.equals("ind") || a.equals("dumpin")){
				dumpIn=b;
			}else if(a.equals("pound")){
				pound=Parse.parseBoolean(b);
			}else if(a.equals("loadkmers") || a.equals("usekmers")){
				loadKmers=Parse.parseBoolean(b);
			}else if(a.equals("allkmers")){
				kmersPerRead=(Parse.parseBoolean(b) ? 0 : 1);
			}else if(a.equals("kmersperread")){
				kmersPerRead=Integer.parseInt(b);
			}else if(a.equals("multithreaded")){
				multithreadedLoad=multithreadedFill=Parse.parseBoolean(b);
			}else if(a.equals("multiload")){
				multithreadedLoad=Parse.parseBoolean(b);
			}else if(a.equals("multifill")){
				multithreadedFill=Parse.parseBoolean(b);
			}else if(a.equals("bloom")){
				bloomMode=Parse.parseBoolean(b);
			}else if(a.equals("bits") || a.equals("cbits")){
				cbits=Integer.parseInt(b);
			}else if(a.equals("hashes")){
				hashes=Integer.parseInt(b);
			}else if(a.equals("ref")){
				refPath=b;
			}else if(a.equals("minid")){
				minIdentity=Float.parseFloat(b);
				if(minIdentity>1) {minIdentity/=100;}
			}else if(a.equals("kalign") || a.equals("alignk")){
				alignK=Integer.parseInt(b);
			}else if(a.equals("expectedbarcodes") || a.equals("expected") || a.equals("barcodes")){
				expectedBarcodes=b;
			}else if(a.equals("shortheader")){
				MicroTile.shortHeader=Parse.parseBoolean(b);
			}else if(a.equals("longheader")){
				MicroTile.shortHeader=!Parse.parseBoolean(b);
			}
			
			else if(a.equals("minpolyg")){
				MicroTile.MIN_POLY_G=Integer.parseInt(b);
			}else if(a.equals("trackcycles")){
				MicroTile.TRACK_CYCLES=Parse.parseBoolean(b);
			}else if(a.equals("extra")){
				if(b!=null) {
					for(String s : b.split(",")) {
						extra.add(s);
					}
				}
			}
			
			else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(b==null && new File(arg).exists()){
				extra.add(b);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;

			if(dump==null && parser.out1!=null){dump=parser.out1;}
			
			extin=parser.extin;
		}
		
		if("phix".equalsIgnoreCase(refPath)) {
			refPath=Data.findPath("?phix2.fa.gz");
		}
		
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
		
		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null) : "\nin1="+in1+"\nin2="+in2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, dump)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+dump+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, dump)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	public void process(Timer t){
		
		final int oldThreads=Shared.threads();
		Shared.capThreads(kmersPerRead>0 ? 64 : 64);
		
		barcodeStats=loadBarcodes(expectedBarcodes);
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		MicroIndex3 index=(refPath==null ? null : new MicroIndex3(alignK, 0, refPath, false));
		mapper=(index==null ? null : new MicroAligner3(index, minIdentity, true));

		if(dumpIn==null){
			flowcell=new FlowCell(k);
			if(loadKmers){loadKmers();}
			fillTiles();
			keySets=null;
			bloomFilter=null;
			bloomCorrector=null;
		}else{
			flowcell=new FlowCell(dumpIn);
			
			if(flowcell.avgReads<targetAverageReads){
				flowcell=flowcell.widenToTargetReads(targetAverageReads);
			}
		}
		
		Shared.setThreads(oldThreads);
	}

	/** Create read streams and process all data */
	void loadKmers(){
		Timer t2=new Timer();
		outstream.print("Loading kmers:  \t");
		if(bloomMode) {
			BloomFilter.printMem=false;
			KmerCountAbstract.KMERS_PER_READ=kmersPerRead;
		}
		if(bloomMode && multithreadedLoad) {
			KmerCountAbstract.CANONICAL=true;
//			ReadCounter rc=new ReadCounter(k, true, false, true, false);
			bloomFilter=new BloomFilter(in1, in2, extra, k, k, cbits, hashes, 1, true, false, true, 0.7f);
		}else {
			//Create a read input stream
			final ConcurrentReadInputStream cris;
			{
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
				cris.start(); //Start the stream
				if(verbose){outstream.println("Started cris");}
			}
			boolean paired=cris.paired();
//			if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
			
			//Process the read stream
			loadKmersInner(cris);
			
			if(verbose){outstream.println("Finished; closing streams.");}
			
			//Close the read streams
			errorState|=ReadWrite.closeStreams(cris);
		}
		
		t2.stop();
		outstream.println(t2);
		
		if(bloomFilter!=null) {
			bloomCorrector=new BloomFilterCorrector(bloomFilter, k, k);
			double used=bloomFilter.filter.usedFraction();
			long unique=(long)bloomFilter.filter.estimateUniqueKmersFromUsedFraction(hashes, used);
			System.err.println(String.format("Bloom Occupancy:\t%.2f%%", 100*used));
			System.err.println("Unique Kmers:   \t"+unique);
//			System.err.println("Cells:          \t"+bloomFilter.filter.cells);
//			System.err.println("Cells Used:     \t"+bloomFilter.filter.cellsUsed());
			assert(keySets==null);
			BloomFilter.printMem=true;
		}
		KmerCountAbstract.KMERS_PER_READ=0;
	}
	
	BarcodeStats loadBarcodes(String expectedBarcodesFile) {
		if(delimiter<0) {
			delimiter=(byte)ffin1.barcodeDelimiter();
			barcodesPerRead=ffin1.barcodesPerRead();
		}
		BarcodeStats bs=new BarcodeStats(delimiter, barcodesPerRead, extin);
		if(bs.length1<1 && bs.length2<1) {
			bs.length1=ffin1.barcodeLength(1);
			bs.length2=ffin1.barcodeLength(2);
		}
		if(expectedBarcodesFile!=null) {
			bs.loadBarcodeList(expectedBarcodesFile, barcodesPerRead>1 ? delimiter : 0, false, false);
		}
		return bs;
	}

	/** Create read streams and process all data */
	void fillTiles(){
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
//		boolean paired=cris.paired();
//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}

		Timer t2=new Timer();
		
		//Process the read stream
		outstream.print("Filling tiles:  \t");
		if(multithreadedFill) {
			spawnThreads(cris);
		}else {
			fillTilesInner(cris);
		}
		
		t2.stop();
		outstream.println(t2);
		
		ArrayList<MicroTile> mtList=flowcell.calcStats();
		if(flowcell.avgReads<targetAverageReads){
			flowcell=flowcell.widenToTargetReads(targetAverageReads);
			mtList=flowcell.toList();
		}
		
		if(dump!=null){
			flowcell.dump(dump, overwrite);
		}
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Iterate through the reads */
	private void loadKmersInner(final ConcurrentReadInputStream cris){

		if(bloomMode) {
			bloomFilter=new BloomFilter(k, k, cbits, hashes, 1, true, false, true, 0.7f);
		}else {
			keySets=new AbstractKmerTable[WAYS];

			//Initialize tables
			ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, 12, false, 0.8);
			int[] schedule=scheduleMaker.makeSchedule();
			for(int j=0; j<WAYS; j++){
				keySets[j]=new HashArray1D(schedule, -1L);
			}
		}
		//Do anything necessary prior to processing
		loadKmersST(cris);
		if(bloomMode) {bloomFilter.filter.shutdown();}
		//Do anything necessary after processing
	}
	
	private void loadKmersST(final ConcurrentReadInputStream cris) {
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();
		//Grab the actual read list from the ListNum
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		//Check to ensure pairing is as expected
		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
		}
		
		LongList kmers=new LongList(300);

		//As long as there is a nonempty read list...
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				final Read r2=r1.mate;

				loadKmers(r1, kmers);
				loadKmers(r2, kmers);
			}

			//Notify the input stream that the list was used
			cris.returnList(ln);
			if(verbose){outstream.println("Returned a list.");}

			//Fetch a new list
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}
	
	private void loadKmers(Read r, LongList kmers) {
		if(r==null || r.length()<k) {return;}
//		if(!randy.nextBoolean()) {return;}//Speed optimization, I guess
		if(kmersPerRead<1) {
			kmers.clear();
			BloomFilter.toKmers(r, kmers, k, 0, 0, true);
			if(keySets!=null) {
				for(int i=0; i<kmers.size; i++) {
					long kmer=kmers.array[i];
					AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
					table.increment(kmer, 1);
				}
			}else {
				for(int i=0; i<kmers.size; i++) {
					long kmer=kmers.array[i];
					bloomFilter.filter.increment(kmer);
				}
			}
		}else {
//			final long kmer=toKmer(r.bases, randy.nextInt(r.length()-k2), k);
			final long kmer=toKmer(r.bases, (int)(r.numericID%(r.length()-k2)), k);
			if(kmer>=0){
				if(keySets!=null) {
					AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
					table.increment(kmer, 1);
				}else {
					bloomFilter.filter.increment(kmer);
				}
			}
		}
	}
	
	/** Iterate through the reads */
	private void fillTilesInner(final ConcurrentReadInputStream cris){
		
		LongList kmers=new LongList();
		IntList counts=new IntList();
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=r1.mateLength();
					
					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;
					
					final MicroTile mt=flowcell.getMicroTile(r1.id);
					
					if(loadKmers){
						processTileKmers(r1, mt, kmers, counts);
						processTileKmers(r2, mt, kmers, counts);
					}
					
					if(mapper!=null) {
						mapper.map(r1);
						mapper.map(r2);
					}
					
					if(barcodeStats!=null) {
						String code=r1.barcode(true);
						int hdist=barcodeStats.calcHdist(code);
						mt.barcodes+=barcodesPerRead;
						mt.barcodeHDistSum+=hdist;
					}

					mt.add(r1);
					mt.add(r2);
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
	}
	
	void processTileKmers(Read r, MicroTile mt, LongList kmers, IntList counts) {
		if(r==null || r.length()<k) {return;}
		
		if(kmersPerRead<1 && bloomMode) {
			processAllKmers(r, mt, kmers, counts);
		}else {
			processOneKmer(r, mt);
		}
	}
	
	void processOneKmer(Read r, MicroTile mt) {
//		final long kmer=toKmer(r.bases, randy.nextInt(r.length()-k2), k);
		final long kmer=toKmer(r.bases, (int)((r.numericID+1)%(r.length()-k2)), k);
		if(kmer<0) {
			mt.misses++;
			return;
		}
		final long key=(kmersPerRead<1 ? Tools.max(kmer, AminoAcid.reverseComplementBinaryFast(kmer, k)) : kmer);
		final int value=(keySets==null ? bloomFilter.getCount(key) : keySets[(int)(key%WAYS)].getValue(key));
		mt.depthSum+=value;
		final int cutoff=(kmersPerRead<1 ? 2 : 1);
		if(value>=cutoff) {mt.hits++;}
		else {mt.misses++;}
	}
	
	void processAllKmers(Read r, MicroTile mt, LongList kmers, IntList counts) {
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		
		int valid=bloomCorrector.fillKmers(bases, kmers);
		if(valid<2){return;}
//		if(!r.containsUndefined() && !bloomCorrector.hasErrorsFast(kmers)){return;}
		
		
		bloomCorrector.fillCounts(bases, kmers, counts);
		final int depth=(counts.size<1 ? 0 : counts.get((int)((r.numericID+1)%counts.size)));
		mt.depthSum+=depth;
		if(depth>1) {mt.hits++;}
		else {mt.misses++;}
		final int possibleErrors=bloomCorrector.countErrors(counts, quals);
		mt.kmerBaseErrorCount+=possibleErrors;
		mt.kmerReadErrorCount+=(possibleErrors>0 ? 1 : 0);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Threading           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<WorkerThread> alpt=new ArrayList<WorkerThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new WorkerThread(cris, i));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public final void accumulate(WorkerThread pt){
//		System.err.println("Accumulating "+pt.tid+": "+pt.readsProcessedT);
		synchronized(pt) {
			flowcell.add(pt.fc);
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			errorState|=(!pt.success);
		}
//		System.err.println("readsProcessed: "+readsProcessed+", "+flowcell.readsProcessed);
//		pt.fc.dump(pt.tid+".txt", overwrite);
	}
	
	@Override
	public final boolean success() {return !errorState;}
	
	class WorkerThread extends Thread {
		
		WorkerThread(ConcurrentReadInputStream cris_, int tid_){
			cris=cris_;
			tid=tid_;
			fc=new FlowCell(k);
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
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				final Read r2=r1.mate;
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}
				if(r2!=null && !r2.validated()){r2.validate(true);}
				
				{
					//Reads are processed in this block.
					processReadPair(r1, r2);
				}
			}

			//Output reads to the output stream
//			if(ros!=null){ros.add(reads, ln.id);}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		void processReadPair(final Read r1, final Read r2){
			//Increment counters
			readsProcessedT+=r1.pairCount();
			basesProcessedT+=r1.pairLength();
			
			final MicroTile mt=fc.getMicroTile(r1.id);
			
			if(loadKmers){
				processTileKmers(r1, mt, kmers, counts);
				processTileKmers(r2, mt, kmers, counts);
			}
			
			if(mapper!=null) {
				mapper.map(r1);
				mapper.map(r2);
			}
			
			if(barcodeStats!=null) {
				String code=r1.barcode(true);
				int hdist=barcodeStats.calcHdist(code);
				mt.barcodes+=barcodesPerRead;
				mt.barcodeHDistSum+=hdist;
			}
			
			mt.add(r1);
			mt.add(r2);
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
		
		final FlowCell fc;
		
		final LongList kmers=new LongList();
		final IntList counts=new IntList();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Generate a kmer from specified start location
	 * @param bases
	 * @param start
	 * @param klen kmer length
	 * @return kmer
	 */
	private static final long toKmer(final byte[] bases, final int start, final int klen){
		final int stop=start+klen;
		assert(stop<=bases.length) : klen+", "+bases.length;
		long kmer=0;
		
		for(int i=start; i<stop; i++){
			final byte b=bases[i];
			final long x=Dedupe.baseToNumber[b];
			kmer=((kmer<<2)|x);
		}
		return kmer;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private ArrayList<String> extra=new ArrayList<String>();
	
	/** Override input file extension */
	private String extin=null;
	
	private boolean pound=true;
	private String dump=null;
	private String dumpIn=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	public long readsProcessed=0;
	/** Number of bases processed */
	public long basesProcessed=0;
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	private AbstractKmerTable[] keySets;
	
	private BloomFilter bloomFilter;
	private BloomFilterCorrector bloomCorrector;
	
	private boolean bloomMode=false;
	private boolean multithreadedLoad=false;
	private boolean multithreadedFill=false;

	private boolean loadKmers=true;
	private int kmersPerRead=1;
	private int cbits=4;
	private int hashes=3;
	
	private int targetAverageReads=800;
	
	private static final int WAYS=31;
	private static final int k=31, k2=30;
	
//	private final Random randy=Shared.threadLocalRandom();
	private FlowCell flowcell;

	private String refPath="phix";
//	private byte[] ref;
//	private LongHashMap refIndex;
	private MicroAligner3 mapper;
	float minIdentity=0.65f;
	int alignK=19; //Phix is unique down to k=13
	
	
	private byte delimiter='+';
	private int barcodesPerRead=2;
	private String expectedBarcodes;
	private BarcodeStats barcodeStats;
//	private long minCountToUse=0;
//	
//	private boolean warned=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
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
	
}
