package hiseq;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import aligner.SideChannel3;
import barcode.Barcode;
import barcode.BarcodeStats;
import bloom.BloomFilter;
import bloom.KmerCountAbstract;
import bloom.PolyFilter;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import jgi.CalcTrueQuality;
import jgi.Dedupe;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.SamLineStreamer;
import structures.AtomicStringNum;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import structures.LongList;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * Analyzes a flow cell for low-quality areas.
 * Removes reads in the low-quality areas.
 * 
 * @author Brian Bushnell
 * @date August 31, 2016
 *
 */
public class AnalyzeFlowCell implements Accumulator<AnalyzeFlowCell.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		AnalyzeFlowCell x=new AnalyzeFlowCell(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public AnalyzeFlowCell(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=parse(args);
		
		if(gToN || discardG){MicroTile.TRACK_CYCLES=true;}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			out2=parser.out2;
			
			extin=parser.extin;
			extout=parser.extout;
			

			trimq=parser.trimq;
			trimE=parser.trimE();
			minlen=parser.minReadLength;
			trimLeft=parser.qtrimLeft;
			trimRight=parser.qtrimRight;
		}
		

		Read.VALIDATE_IN_CONSTRUCTOR=true;
		align=(align || alignOut!=null) && samInput==null && dumpIn==null;
		
		if(recalibrate) {
			CalcTrueQuality.initializeMatrices();
		}
//		if(mapper!=null) {
//			assert(mapper.getMap().size()>0) : refPath+", "+alignK+", "+minIdentity;
//		}
		checkFiles();
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutbad=FileFormat.testOutput(outbad, FileFormat.FASTQ, extout, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	private Parser parse(String[] args){
		//Create a parser object
		Parser parser=new Parser();
		parser.qtrimRight=trimRight;
		parser.trimq=trimq;
		parser.minReadLength=minlen;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("outb") || a.equals("outbad")){
				outbad=b;
			}else if(a.equals("sam") || a.equals("insam") || 
					a.equals("samin") || a.equals("saminput")){
				samInput=b;
			}else if(a.equals("sammt")){
//				processSamMT=Parse.parseBoolean(b);
			}else if(a.equals("divisor") || a.equals("size")){
				targetX=targetY=Tile.xSize=Tile.ySize=Parse.parseIntKMG(b);
			}else if(a.equals("xdivisor") || a.equals("xsize") || a.equals("x")){
				targetX=Tile.xSize=Parse.parseIntKMG(b);
			}else if(a.equals("ydivisor") || a.equals("ysize") || a.equals("y")){
				targetY=Tile.ySize=Parse.parseIntKMG(b);
			}else if(a.equals("target") || a.equals("targetreads")){
				targetAverageReads=Parse.parseIntKMG(b);
			}else if(a.equals("targetalignedreads") || a.equals("alignedreads")){
				targetAlignedReads=Parse.parseIntKMG(b);
			}else if(a.equals("dump") || a.equals("dumpout") || a.equals("outdump")){
				dumpOut=b;
			}else if(a.equals("indump") || a.equals("ind") || a.equals("dumpin")){
				dumpIn=b;
			}else if(a.equals("filterout") || a.equals("filterlist") || a.equals("coordinates")
					 || a.equals("coords") || a.equals("coordsout") || a.equals("coordinatesout")){
				coordsOut=b;
			}else if(a.equals("loadkmers") || a.equals("usekmers")){
				loadKmers=Parse.parseBoolean(b);
			}else if(a.equals("loadthreads")){
				loadThreads=Integer.parseInt(b);
			}else if(a.equals("fillthreads")){
				fillThreads=Integer.parseInt(b);
			}else if(a.equals("minprob")){
				minProb=Float.parseFloat(b);
			}else if(a.equals("bits") || a.equals("cbits")){
				cbits=Integer.parseInt(b);
			}else if(a.equals("hashes")){
				hashes=Integer.parseInt(b);
			}
			
			else if(a.equals("smooth") || a.equals("fixspikes")){
				if(!Tools.startsWithDigit(b)) {
					smoothDepths=Parse.parseBoolean(b) ? 3 : 0;
				}else {
					smoothDepths=Integer.parseInt(b);
				}
			}else if(a.equals("deblur") || a.equals("deconvolute")){
				deblurDepths=Parse.parseBoolean(b);
			}
			
			else if(a.equals("blur") || a.equals("blurtiles") || a.equals("smoothtiles")){
				blurTiles=Parse.parseBoolean(b);
			}
			
			else if(a.equals("recalibrate") || a.equals("recal")) {
				recalibrate=Parse.parseBoolean(b);
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
			}else if(a.equals("expectedbarcodes") || a.equals("expected") || a.equals("barcodes") || a.equals("barcodesin")){
				expectedBarcodes=b;
			}else if(a.equals("barcodesout") || a.equals("barcodecounts") || a.equals("counts")){
				barcodeCounts=b;
			}
			
			else if(a.equals("merge")){
				merge=Parse.parseBoolean(b);
			}else if(a.equals("strict")){
				strictmerge=Parse.parseBoolean(b);
			}else if(a.equals("loose")){
				strictmerge=!Parse.parseBoolean(b);
			}
			
			else if(a.equals("alignout") || a.equals("sideout") || a.equals("phixout") || a.equals("outphix") || 
					a.equals("outsam") || a.equals("samout")){
				alignOut=b;
			}else if(a.equals("align")){
				align=Parse.parseBoolean(b);
			}else if(a.equals("ref") || a.equals("alignref") || a.equals("sideref")){
				alignRef=b;
			}else if(a.equals("alignk") || a.equals("sidek") || a.equals("alignk1") || a.equals("sidek1")){
				alignK1=Integer.parseInt(b);
			}else if(a.equals("alignk2") || a.equals("sidek2")){
				alignK2=Integer.parseInt(b);
			}else if(a.equals("alignminid") || a.equals("alignminid1") || 
					a.equals("sideminid") || a.equals("sideminid1") || a.equals("minid") || a.equals("minid1")){
				alignMinid1=Float.parseFloat(b);
			}else if(a.equals("alignminid2") || a.equals("sideminid2") || a.equals("minid2")){
				alignMinid2=Float.parseFloat(b);
			}else if(a.equals("alignmm1") || a.equals("alignmidmask1") || a.equals("sidemm1") || 
					a.equals("sidemidmask1")){
				alignMM1=Integer.parseInt(b);
			}else if(a.equals("alignmm2") || a.equals("alignmidmask2") || a.equals("sidemm2") || 
					a.equals("sidemidmask2")){
				alignMM2=Integer.parseInt(b);
			}
			
			else if(a.equals("lqo") || a.equals("lowqualityonly")){
				discardOnlyLowQuality=Parse.parseBoolean(b);
			}else if(a.equals("dmult")){
				dmult=Float.parseFloat(b);
			}else if(a.equals("idmaskread")){
				idmask_read=Integer.parseInt(b);
			}else if(a.equals("idmaskwrite")){
				idmask_write=Integer.parseInt(b);
			}else if(a.equals("allkmers")){
//				kmersPerRead=(Parse.parseBoolean(b) ? 0 : 1);
				boolean x=Parse.parseBoolean(b);
				idmask_write=(x ? 0 : 15);
			}
//			else if(a.equals("kmersperread")){
//				kmersPerRead=Integer.parseInt(b);
//			}
			
			else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(TileDump.parseStatic(arg, a, b)){
				//do nothing
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
		return parser;
	}
	
	private void checkFiles(){
		doPoundReplacement();
		adjustInterleaving();
		checkFileExistence();
		checkStatics();
	}
	
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
	
	private void checkFileExistence(){
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outbad, dumpOut, barcodeCounts)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2+", "+outbad);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+", "+outbad+", "+dumpOut+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2, samInput, dumpIn, expectedBarcodes)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outbad, samInput, 
				dumpIn, dumpOut, expectedBarcodes, barcodeCounts)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
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
	
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	public void process(Timer t){
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		if(dumpIn==null){
			barcodeStats=loadBarcodes(expectedBarcodes);
			if(barcodeCounts!=null) {
				barcodeMap=new ConcurrentHashMap<String, AtomicStringNum>();
			}
			if(loadKmers){loadKmers();}
			flowcell=new FlowCell(k);
			fillTiles();
			
			bloomFilter=null;//Clearing before widen saves memory
			ArrayList<MicroTile> mtList;
//			if(!processSamMT) {loadSam_ST(samInput);}
			Timer t2=new Timer();
			boolean showTime=false;
			synchronized(flowcell) {
				mtList=flowcell.calcStats();
				if(showTime) {t2.stopAndStart("calcStats: ");}
				if(flowcell.avgReads<targetAverageReads){
					flowcell=flowcell.widenToTargetReads(targetAverageReads);
					if(showTime) {t2.stopAndStart("Widen: ");}
					mtList=flowcell.toList();
					if(showTime) {t2.stopAndStart("toList: ");}
				}
				if(blurTiles) {
					flowcell.blur();
					if(showTime) {t2.stopAndStart("Blur: ");}
				}
				//Temporarily widen to calculate a regression
				if(flowcell.avgAlignedReads<targetAlignedReads && flowcell.readsAligned>targetAlignedReads) {
					final int oldX=Tile.xSize, oldY=Tile.ySize;
					FlowCell temp=flowcell.widenToTargetAlignedReads(targetAlignedReads);
					temp.calcStats();
					flowcell.uniqueToReadErrorRateFormula=temp.uniqueToReadErrorRateFormula;
					flowcell.uniqueToBaseErrorRateFormula=temp.uniqueToBaseErrorRateFormula;
					Tile.xSize=oldX;
					Tile.ySize=oldY;
					if(showTime) {t2.stopAndStart("Widen2: ");}
				}
//				if(MicroTile.trackDepth) {flowcell.summarizeDepth();}
			}

			long readsToDiscard=TileDump.markTiles(flowcell, mtList, outstream);

			if(dumpOut!=null){
				flowcell.dump(dumpOut, overwrite);
				if(showTime) {t2.stopAndStart("Dump: ");}
			}
			if(barcodeCounts!=null) {
				dumpBarcodes(barcodeMap.values(), barcodeCounts, overwrite);
				if(showTime) {t2.stopAndStart("Dump Barcodes: ");}
			}
		}else{
			flowcell=new FlowCell(dumpIn);
//			if(loadKmers){loadKmers();}
			if(targetX>Tile.xSize || targetY>Tile.ySize) {
				flowcell=flowcell.widen(targetX, targetY, true);
			}
			if(flowcell.avgReads<targetAverageReads){
				flowcell.calcStats();//May be necessary for calculating average reads
				flowcell=flowcell.widenToTargetReads(targetAverageReads);
			}
			if(blurTiles) {flowcell.blur();}
			ArrayList<MicroTile> mtList=flowcell.calcStats();
			long readsToDiscard=TileDump.markTiles(flowcell, mtList, outstream);
		}
		System.err.println("Avg quality:     \t"+Tools.format("%.3f", flowcell.avgQuality));
		System.err.println("STDev quality:   \t"+Tools.format("%.4f", flowcell.stdQuality));
		System.err.println("Avg uniqueness:  \t"+Tools.format("%.4f", flowcell.avgUnique));
		System.err.println("STDev uniqueness:\t"+Tools.format("%.4f", flowcell.stdUnique));
		System.err.println("Avg depth:       \t"+Tools.format("%.4f", flowcell.avgDepth));
		System.err.println("STDev depth:     \t"+Tools.format("%.4f", flowcell.stdDepth));
		System.err.println("Avg poly-G:      \t"+Tools.format("%.4f", flowcell.avgPolyG));
		System.err.println("STDev poly-G:    \t"+Tools.format("%.4f", flowcell.stdPolyG));
		System.err.println("Alignment Rate:  \t"+Tools.format("%.8f", flowcell.alignmentRate()));
		System.err.println("Base Error Rate: \t"+Tools.format("%.8f", flowcell.baseErrorRate()));
		
		if(sidechannel!=null && sidechannel.samOut) {Data.unloadAll();}
		
		processReads(t);
	}

	/** Create read streams and process all data */
	void loadKmers(){
		Timer t2=new Timer();
		outstream.print("Loading kmers:  \t");
		
		loadThreads=Tools.min(loadThreads, Shared.threads());
		final int oldMCT=KmerCountAbstract.MAX_COUNT_THREADS;
		final float oldProb=KmerCountAbstract.minProb;
//		KmerCountAbstract.KMERS_PER_READ=kmersPerRead;
		KmerCountAbstract.IDMASK=idmask_write;
		KmerCountAbstract.minProb=minProb;
		KmerCountAbstract.CANONICAL=true;
		
		if(loadThreads>1) {
			if(idmask_write>7) {
				KmerCountAbstract.MAX_COUNT_THREADS=loadThreads;
			}
			bloomFilter=new BloomFilter(in1, in2, extra, k, k, cbits, hashes, 
					1, true, false, false, 0.65f);
			bloomFilter.filter.shutdown();
		}else {
			//Create a read input stream
			final ConcurrentReadInputStream cris;
			{
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
				cris.start(); //Start the stream
				if(verbose){outstream.println("Started cris");}
			}
			boolean paired=cris.paired();
			//Process the read stream
			loadKmersInner(cris);
			
//			if(verbose){outstream.println("Finished; closing streams.");}
			
			//Close the read streams
			errorState|=ReadWrite.closeStreams(cris);
			bloomFilter.filter.shutdown();
		}
		
		KmerCountAbstract.KMERS_PER_READ=0;
		KmerCountAbstract.IDMASK=0;
		KmerCountAbstract.MAX_COUNT_THREADS=oldMCT;
		KmerCountAbstract.minProb=oldProb;
		
		double used=bloomFilter.filter.usedFraction();
		long unique=(long)bloomFilter.filter.estimateUniqueKmersFromUsedFraction(hashes, used);
		System.err.println(String.format("Bloom Occupancy:\t%.2f%%", 100*used));
		System.err.println("Unique Kmers:   \t"+unique);
		
		t2.stop();
		outstream.println(t2);
	}

	/** Create read streams and process all data */
	void fillTiles(){
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		
		SamLineStreamer ss=null;
		if(samInput!=null && processSamMT) {
			outstream.println("Loading sam file.");
			FileFormat ff=FileFormat.testInput(samInput, FileFormat.SAM, null, true, false);
			final int streamerThreads=Tools.min(4, Shared.threads());
			ss=new SamLineStreamer(ff, streamerThreads, false, maxReads);
			ss.start();
		}
		

		if(align) {
			sidechannel=new SideChannel3(alignRef, alignOut, null, alignK1, alignK2, 
					alignMinid1, alignMinid2, alignMM1, alignMM2, overwrite, ordered);
			sidechannel.start();
		}else {
			sidechannel=null;
		}
		
		//Process the read stream
		fillTilesInner(cris, ss);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		if(sidechannel!=null) {errorState|=sidechannel.shutdown();}
	}
	
//	/** Singlethreaded version. */
//	private void loadSam_ST(String fname) {
//		if(fname==null || processSamMT) {return;}
//		Timer t=new Timer();
//		outstream.println("Loading sam file.");
//		final SamReadStreamer ss;
//		FileFormat ff=FileFormat.testInput(fname, FileFormat.SAM, null, true, false);
//		final int streamerThreads=Tools.min(4, Shared.threads());
//		
//		ss=new SamReadStreamer(ff, streamerThreads, false, maxReads);
//		ss.start();
//
//		ListNum<Read> ln=ss.nextList();
//		ArrayList<Read> reads=(ln==null ? null : ln.list);
//		final IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
//
//		while(ln!=null && reads!=null && reads.size()>0){
//
//			for(int idx=0; idx<reads.size(); idx++){
//				Read r=reads.get(idx);
//				assert(r.mate==null);
//				processSamLine(r, ihp);
//			}
//			ln=ss.nextList();
//			reads=(ln==null ? null : ln.list);
//		}
//		t.stopAndPrint();
//	}
//	
//	private void processSamLine(Read r, IlluminaHeaderParser2 ihp) {
//		if(r==null){return;}
//		final SamLine sl=r.samline;
////		if(!sl.mapped() && !sl.nextMapped()) {return;} //Probably not PhiX
//		if(!sl.mapped()) {return;}//TODO: Track unmapped info; requires modifying dump format
//		
//		
//		if(!sl.mapped() || r.bases==null || r.match==null) {return;}
//		final int pairnum=sl.pairnum();
//		assert(sl.strand()==r.strand());
//		
////		final boolean needsFixing=(varMap!=null && Read.containsVars(r.match));
//		
//		if(r.shortmatch()){r.toLongMatchString(false);}
//		final byte[] match=r.match;
//
//		int subs=0, inss=0, dels=0;
////		final AtomicLongArray matchCounts=lane.matchCounts[pairnum];
////		final AtomicLongArray subCounts=lane.subCounts[pairnum];
//		for(int mpos=0, qpos=0; mpos<match.length; mpos++){
//			byte m=match[mpos];
//			if(m=='m'){
////				matchCounts.incrementAndGet(qpos);
//				qpos++;
//			}else if(m=='S' || m=='N'){
////				subCounts.incrementAndGet(qpos);
//				subs++;
//				qpos++;
//			}else if(m=='I'){
////				subCounts.incrementAndGet(qpos);
//				inss++;
//				qpos++;
//			}else if(m=='X' || m=='Y' || m=='C'){
//				qpos++;
//			}else if(m=='D'){
//				dels++;
//			}else{
//				assert(false) : "Unhandled symbol "+m;
//			}
//		}
//		
//		final MicroTile mt;
//		final Lane lane;
//		ihp.parse(r.id);
//		final int lnum=ihp.lane(), tile=ihp.tile(), x=ihp.xPos(), y=ihp.yPos();
//		synchronized(flowcell) {
//			lane=flowcell.getLane(lnum);
//			mt=lane.getMicroTile(tile, x, y);
//		}
//
//		synchronized(mt) {
//			mt.alignedReadCount++;
//			mt.alignedBaseCount+=r.countAlignedBases();
//			mt.readErrorCount+=((subs+inss+dels>0) ? 1 : 0);
//			mt.baseErrorCount+=(subs+inss);
//			mt.readInsCount+=(inss>0 ? 1 : 0);
//			mt.readDelCount+=(dels>0 ? 1 : 0);
//		}
//	}
	
	//This version is not really any faster,
	//though it does use 15% less CPU-time.
	private void processSamLine(SamLine sl, IlluminaHeaderParser2 ihp) {
		if(sl==null){return;}
//		if(!sl.mapped() && !sl.nextMapped()) {return;} //Probably not PhiX
		if(!sl.mapped()) {return;}//TODO: Track unmapped info; requires modifying dump format
		
		if(!sl.mapped() || sl.seq==null || !sl.hasCigar()) {return;}
		final int pairnum=sl.pairnum();
		
//		final boolean needsFixing=(varMap!=null && Read.containsVars(r.match));
		
		final byte[] shortmatch=sl.toShortMatch(false);
		final byte[] match=Read.toLongMatchString(shortmatch);

		int subs=0, inss=0, dels=0;
		//Atomics are too slow here
//		final AtomicLongArray matchCounts=lane.matchCounts[pairnum];
//		final AtomicLongArray subCounts=lane.subCounts[pairnum];
		for(int mpos=0, qpos=0; mpos<match.length; mpos++){
			byte m=match[mpos];
			if(m=='m'){
//				matchCounts.incrementAndGet(qpos);
				qpos++;
			}else if(m=='S' || m=='N'){
//				subCounts.incrementAndGet(qpos);
				subs++;
				qpos++;
			}else if(m=='I'){
//				subCounts.incrementAndGet(qpos);
				inss++;
				qpos++;
			}else if(m=='X' || m=='Y' || m=='C'){
				qpos++;
			}else if(m=='D'){
				dels++;
			}else{
				assert(false) : "Unhandled symbol "+m;
			}
		}
		
		final MicroTile mt;
		final Lane lane;
		ihp.parse(sl.qname);
		final int lnum=ihp.lane(), tile=ihp.tile(), x=ihp.xPos(), y=ihp.yPos();
		synchronized(flowcell) {
			lane=flowcell.getLane(lnum);
			mt=lane.getMicroTile(tile, x, y);
		}

		synchronized(mt) {
			mt.alignedReadCount++;
			mt.alignedBaseCount+=Read.countAlignedBases(match);
			mt.readErrorCount+=((subs+inss+dels>0) ? 1 : 0);
			mt.baseErrorCount+=(subs+inss);
			mt.readInsCount+=(inss>0 ? 1 : 0);
			mt.readDelCount+=(dels>0 ? 1 : 0);
		}
	}

	/** Create read streams and process all data */
	void processReads(Timer t){
		
		if(ffout1!=null || ffoutbad!=null || coordsOut!=null){
			Timer t2=new Timer();
			outstream.print("Filtering reads:\t");

			//Create a read input stream
			final ConcurrentReadInputStream cris;
			{
				Read.VALIDATE_IN_CONSTRUCTOR=true;
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
				cris.start(); //Start the stream
				if(verbose){outstream.println("Started cris");}
			}
			boolean paired=cris.paired();
			//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}

			//Optionally read output streams
			final ConcurrentReadOutputStream ros, rosb;
			final int buff=4;
			if(ffout1!=null){
				ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, buff, null, false);
				ros.start(); //Start the stream
			}else{ros=null;}
			
			if(ffoutbad!=null){
				rosb=ConcurrentReadOutputStream.getStream(ffoutbad, null, null, null, buff, null, false);
				rosb.start(); //Start the stream
			}else{rosb=null;}
			
			final ByteStreamWriter coords=(coordsOut==null ? null : 
				ByteStreamWriter.makeBSW(coordsOut, overwrite, append, true));
			
			//Process the read stream
			processInner(cris, ros, rosb, coords);

			if(verbose){outstream.println("Finished; closing streams.");}

			//Close the read streams
			errorState|=ReadWrite.closeStreams(cris, ros, rosb);
			
			if(coords!=null) {
				errorState=coords.poisonAndWait()|errorState;
			}

			t2.stop();
			outstream.println(t2);
		}
		
		//Report timing and results
		{
			t.stop();
			lastReadsOut=readsProcessed-readsDiscarded;
			outstream.println();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
			
			if(ffout1!=null || ffoutbad!=null || coordsOut!=null){

				String rpstring=Tools.padKMB(readsDiscarded, 8);
				String bpstring=Tools.padKMB(basesDiscarded, 8);
				String gpstring=Tools.padKMB(gsTransformedToN, 8);
				outstream.println();
				outstream.println("Reads Discarded:    "+rpstring+" \t"+Tools.format("%.3f%%", readsDiscarded*100.0/readsProcessed));
				outstream.println("Bases Discarded:    "+bpstring+" \t"+Tools.format("%.3f%%", basesDiscarded*100.0/basesProcessed));
				if(gToN){outstream.println("Gs Masked By N:     "+gpstring+" \t"+Tools.format("%.3f%%", gsTransformedToN*100.0/basesProcessed));}
				outstream.println();
			}
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the reads */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, 
			final ConcurrentReadOutputStream rosb, final ByteStreamWriter coords){
		
		//Do anything necessary prior to processing
		readsProcessed=0;
		basesProcessed=0;
		final IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
		final ByteBuilder bb=new ByteBuilder();
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

				ArrayList<Read> keepList=new ArrayList<Read>(reads.size());
				ArrayList<Read> tossList=new ArrayList<Read>(4);
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;
					
					boolean keep=processReadPair(r1, r2);
					if(keep){
						keepList.add(r1);
					}else{
						tossList.add(r1);
						readsDiscarded+=r1.pairCount();
						basesDiscarded+=initialLength1+initialLength2;
					}
				}
				
				//Output reads to the output stream
				if(ros!=null){ros.add(keepList, ln.id);}
				if(rosb!=null){rosb.add(tossList, ln.id);}
				if(coords!=null) {
					bb.clear();
					for(Read r : tossList) {
						ihp.parse(r);
						ihp.appendCoordinates(bb);
						bb.nl();
					}
					if(!bb.isEmpty()) {coords.print(bb);}
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
		
		//Do anything necessary after processing
		
	}
	
	/** Iterate through the reads */
	public void loadKmersInner(final ConcurrentReadInputStream cris){
		
		bloomFilter=new BloomFilter(k, k, cbits, hashes, 1, true, false, true, 0.7f);
		
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
			
			LongList kmers=new LongList(300);
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

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
	}
	
	private void loadKmers(Read r, LongList kmers) {
		if(r==null || r.length()<k) {return;}
//		if(!randy.nextBoolean()) {return;}//Speed optimization, I guess
		
		kmers.clear();
		BloomFilter.toKmers(r, kmers, k, 0, 0, true);
		final int idmod=(int)((r.numericID+3)&idmask_write);
		for(int i=0; i<kmers.size; i++) {
			if((i&idmask_write)==idmod) {
				long kmer=kmers.array[i];
				bloomFilter.filter.increment(kmer);
			}
		}
	}
	
	/** Iterate through the reads */
	public void fillTilesInner(final ConcurrentReadInputStream cris, final SamLineStreamer ss){
		Timer t2=new Timer();
		
		if(merge && loadKmers) {fillThreads=Tools.max(fillThreads, fillThreadsM);}
		fillThreads=Tools.min(fillThreads, Shared.threads());
		outstream.print("Filling tiles with "+fillThreads+" threads:  \t");
		spawnThreads(cris, ss);
		
		t2.stop();
		outstream.println(t2);
	}
	
	private void processTileKmers(Read r, IntList klist, IntList blist) {
		fillDepthList(r, klist);
		if(smoothDepths>3) {smooth5(klist);}
		else if(smoothDepths>=3) {smooth3(klist);}
		if(deblurDepths) {deblur(klist, blist);}
	}
	
	//Samples all kmers
	private void fillDepthList(Read r, IntList depths) {
		if(r==null || r.length()<k) {return;}
		final byte[] bases=r.bases;
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(-1L)>>>(64-shift);

		long kmer=0, rkmer=0;
		
		depths.clear();
		
		final int idmod=(int)(r.numericID&idmask_read);
		
		//TODO: Debranch using 2 loops; one for 1st k-1 bases, and one for the rest.
		for(int i=0, len=0; i<bases.length; i++) {
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=x<0 ? 0 : ((rkmer>>>2)|(x2<<shift2))&mask;
			
			len=(x<0 ? 0 : len+1);
			rkmer=(x<0 ? 0 : rkmer);//is this necessary? TODO: Check.
//			if(x<0){len=0; rkmer=0;}else{len++;}//old branchy version
			if(i>=k && ((i&idmask_read)==idmod)) {
				final long key=toKey(kmer, rkmer);
				final int value=(len>=k ? bloomFilter.getCount(key) : 0);
				depths.add(value);
			}
		}
	}
	
	private void smooth3(IntList list) {
		if(list.size<3) {return;}
		final int max=list.size-1;
		{
			int a=list.get(0), b=list.get(1), c=list.get(2);
			a=Tools.min(a, b, c);
			list.set(0, a);
		}
		{
			int a=list.get(max), b=list.get(max-1), c=list.get(max-2);
			a=Tools.min(a, b, c);
			list.set(0, a);
		}
		for(int i=1; i<max; i++) {
			int a=list.get(i-1), b=list.get(i), c=list.get(i+1);
			list.set(i, Tools.min(b, Tools.max(a, c)));
		}
	}
	
	private void smooth5(IntList list) {
		if(list.size<5) {return;}
		final int max=list.size-1;
		{
			int a=list.get(0), b=list.get(1), c=list.get(2), d=list.get(3);
			a=Tools.min(a, b, c);
			b=Tools.min(b, Tools.max(a, Tools.min(c, d)));
			list.set(0, a);
			list.set(1, b);
		}
		{
			int a=list.get(max), b=list.get(max-1), c=list.get(max-2), d=list.get(max-3);
			a=Tools.min(a, b, c);
			b=Tools.min(b, Tools.max(a, Tools.min(c, d)));
			list.set(0, a);
			list.set(1, b);
		}
		for(int i=2; i<max-1; i++) {
			int a=list.get(i-2), b=list.get(i-1), c=list.get(i), d=list.get(i+1), e=list.get(i+2);
			int left=Tools.min(a, b);
			int right=Tools.min(d, e);
			list.set(i, Tools.min(c, Tools.max(left, right)));
		}
	}
	
	void deblur(IntList klist, IntList blist) {
		blist.clearFull();
		blist.set(klist.size+k-1, 0);
		for(int i=0; i<klist.size; i++) {
			int kdepth=klist.get(i);
			//Since max is only 3 this could be much more efficient, without a double loop
			if(kdepth>0) {
				for(int j=i, max=i+k; j<max; j++) {
					blist.set(j, Tools.max(blist.get(i), kdepth));
				}
			}
		}
	}
	
//	//Samples multiple kmers
//	private void processTileKmersSampled(Read r, MicroTile mt, int samples) {
//		if(r.length()<=3*k+1) {samples=1;}
//		samples=Tools.min(samples, r.length()-k2);
//		final int cutoff=(kmersPerRead<1 ? 2 : 1);
//		
//		long depthSum=0;
//		for(int i=0; i<samples; i++) {
//			int value=getKmerCount(r.bases, getPos(r.numericID+74+k3*i, r.length()));
//			depthSum+=value;
//			if(value>=cutoff) {mt.hits++;}
//			else {mt.misses++;}
//		}
//		mt.depthSum+=depthSum;
//	}
//	
//	private int getKmerCount(byte[] bases, int pos) {
//		final int lim=bases.length-k2;
//		pos=pos%(bases.length-k2);
//		assert(pos>=0 && pos<=lim);
//		final long kmer=toKmer(bases, pos, k);
//		if(kmer<0) {return 0;}
//		final long key=toKey(kmer);
//		final int value=bloomFilter.getCount(key);
////		if(maxReads==1) {System.err.println("Got "+value+" for kmer "+key);}
//		return value;
//	}
//	
//	private int getPos(long id, int len) {
//		return (deterministic ? (int)((id)%(len-k2)) : randy.nextInt(len-k2));
//	}
	
	private final long toKey(long kmer) {
		return Tools.max(kmer, AminoAcid.reverseComplementBinaryFast(kmer, k));
//		return (kmersPerRead==1 || kmer==-1 ? kmer : 
//			Tools.max(kmer, AminoAcid.reverseComplementBinaryFast(kmer, k)));
	}
	
	private final long toKey(long kmer, long rkmer) {
		return Tools.max(kmer, rkmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	boolean processReadPair(final Read r1, final Read r2){
		boolean passes=processReadPair_inner(r1, r2);
		if(passes){return true;}
		if(trimq>0){
			TrimRead.trimFast(r1, trimLeft, trimRight, trimq, trimE, 0);
			if(r2!=null){TrimRead.trimFast(r2, trimLeft, trimRight, trimq, trimE, 0);}
			return r1.length()>=minlen && (r2==null || r2.length()>=minlen);
		}else{
			return false;
		}
	}
	
	/**
	 * Process a single read pair.
	 * @param r1 Read 1
	 * @param r2 Read 2 (may be null)
	 * @return True if the reads should be kept, false if they should be discarded.
	 */
	boolean processReadPair_inner(final Read r1, final Read r2){
		
		MicroTile mt=flowcell.getMicroTile(r1.id);
		if(mt==null){
			if(!warned){
				outstream.println("\nWarning - a read was found with no corresponding MicroTile:\n"+r1.id);
				warned=true;
			}
			return true;
		}
		if(mt.discard<discardLevel){return true;}
		if(!discardOnlyLowQuality){return false;}
		
		if(shouldDiscard(r1, mt)) {return false;}
		if(gToN){gsTransformedToN+=doGToN(r1, mt);}
		
		if(shouldDiscard(r2, mt)) {return false;}
		if(gToN){gsTransformedToN+=doGToN(r2, mt);}
		
		return true;
	}
	
	private boolean shouldDiscard(Read r, MicroTile mt) {
		if(r==null || r.length()<1) {return false;}
		final int len=r.length();
		double qual=r.avgQualityByProbabilityDouble(true, len);
		double prob=100*r.probabilityErrorFree(true, len);
		if(qual<=flowcell.avgQuality-(dmult*TileDump.qDeviations*flowcell.stdQuality)){return true;}
		if(prob<=flowcell.avgErrorFree-(dmult*TileDump.eDeviations*flowcell.stdErrorFree)){return true;}
		if(PolyFilter.polymerLen(r.bases, (byte)'G', 0.16f)>15) {return true;}
		if(discardG && shouldDiscardG(r, mt)){return true;}
		return false;
	}
	
	private boolean shouldDiscardG(Read r, MicroTile mt){
		final byte[] bases=r.bases;
		final float[] gArray=mt.tracker.cycleAverages[2];
		
		final float thresh=(float)(flowcell.avgG+Tools.max(TileDump.gDeviations*flowcell.stdG, 
				flowcell.avgG*TileDump.gFraction, TileDump.gAbs));
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(b=='G' && gArray[i]>thresh){
				return true;
			}
		}
		return false;
	}
	
	private int doGToN(Read r, MicroTile mt){
		if(r==null || r.length()<1) {return 0;}
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		final float[] gArray=mt.tracker.cycleAverages[2];
		
		final float thresh=(float)(flowcell.avgG+Tools.max(TileDump.gDeviations*flowcell.stdG, 
				flowcell.avgG*TileDump.gFraction, TileDump.gAbs));
		int changes=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(b=='G' && gArray[i]>thresh){
				bases[i]='N';
				changes++;
				if(quals!=null){quals[i]=0;}
			}
		}
		return changes;
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
	/*----------------           Barcodes           ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	private static void dumpBarcodes(Collection<AtomicStringNum> counts, String fname, boolean overwrite) {
		System.err.println("Writing barcode counts.");
		if(fname==null || counts==null) {return;}
		ArrayList<AtomicStringNum> list=new ArrayList<AtomicStringNum>(counts);
		Collections.sort(list);
		Collections.reverse(list);
		long sum=0;
		for(AtomicStringNum asn : list) {sum+=asn.n.get();}
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, false, true);
		bsw.print("#Barcodes\t").print(sum).nl();
		bsw.print("#Unique\t").print(list.size()).nl();
		bsw.print("#Code\tCount\n");
		for(AtomicStringNum count : list) {
			bsw.print(count.s).tab().print(count.n.get()).nl();
		}
		bsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final SamLineStreamer ss){
		
		//Do anything necessary prior to processing
		
		//Increases concurrency by making flowcell copies
		//2 seems to be sufficient when merging is used. 
		FlowCell[] fca=new FlowCell[Tools.mid(1, (merge || loadKmers ? 2 : 3), (fillThreads+4)/8)];
		fca[0]=flowcell;
		for(int i=1; i<fca.length; i++) {fca[i]=new FlowCell(k);}
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(fillThreads);
		for(int i=0; i<fillThreads; i++){
			alpt.add(new ProcessThread(cris, ss, i, fca[(i%fca.length)]));
		}
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
		//Combine flowcell copies
		synchronized(flowcell) {
			for(int i=1; i<fca.length; i++) {
				synchronized(fca[i]) {
					flowcell.add(fca[i]);
				}
			}
		}
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt) {
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			errorState|=(!pt.success);
			synchronized(flowcell) {
				flowcell.xMin=(flowcell.xMin<0 ? pt.xmin : Tools.min(pt.xmin, flowcell.xMin));
				flowcell.xMax=(flowcell.xMax<0 ? pt.xmax : Tools.max(pt.xmax, flowcell.xMax));
				flowcell.yMin=(flowcell.yMin<0 ? pt.ymin : Tools.min(pt.ymin, flowcell.yMin));
				flowcell.yMax=(flowcell.yMax<0 ? pt.ymax : Tools.max(pt.ymax, flowcell.yMax));
				flowcell.tMin=(flowcell.tMin<0 ? pt.tmin : Tools.min(pt.tmin, flowcell.tMin));
				flowcell.tMax=(flowcell.tMax<0 ? pt.tmax : Tools.max(pt.tmax, flowcell.tMax));
				if(pt.flowcellName!=null) {
					assert(flowcell.name==null || flowcell.name.equals(pt.flowcellName));
					flowcell.name=pt.flowcellName;
				}
//				for(Lane lane : flowcell.lanes) {
//					if(lane!=null) {
//						for(int pairnum=0; pairnum<2; pairnum++) {
//							long[] counts=pt.laneDepthCounts[lane.lane][pairnum];
//							long[] sums=pt.laneDepthSums[lane.lane][pairnum];
//							for(int i=0; i<counts.length; i++) {
//								lane.depthCounts[pairnum].addAndGet(i, counts[i]);
//								lane.depthSums[pairnum].addAndGet(i, sums[i]);
//							}
//						}
//					}
//				}
			}
		}
	}
	
	@Override
	public final boolean success(){return !errorState;}	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final SamLineStreamer ss_, final int tid_, final FlowCell flowcell_){
			cris=cris_;
			ss=ss_;
			tid=tid_;
			flowcellT=flowcell_;
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
			
			if(ss!=null) {processSam();}
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

				//Track the initial length for statistics
				final int initialLength1=r1.length();
				final int initialLength2=r1.mateLength();

				//Increment counters
				readsProcessedT+=r1.pairCount();
				basesProcessedT+=initialLength1+initialLength2;
				
				processReadPair(r1, r2);
			}
			if(sidechannel!=null) {sidechannel.writeByStatus(reads, ln.id);}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 */
		void processReadPair(final Read r1, final Read r2){
			final int cutoff=(idmask_write<=3 && cbits>1) ? 2 : 1;
			if(recalibrate) {
				CalcTrueQuality.recalibrate(r1);
				CalcTrueQuality.recalibrate(r2);
			}
			
			ihp.parse(r1.id);
			final int lnum=ihp.lane(), tile=ihp.tile(), x=ihp.xPos(), y=ihp.yPos();
			xmin=Tools.min(x, xmin);
			xmax=Tools.max(x, xmax);
			ymin=Tools.min(y, ymin);
			ymax=Tools.max(y, ymax);
			tmin=Tools.min(tile, tmin);
			tmax=Tools.max(tile, tmax);
			if(flowcellName==null) {
				flowcellName=ihp.machine()+":"+ihp.run()+":"+ihp.flowcell();
			}
			
			long hits1=0, hits2=0, misses=0, depthSum=0;
			if(loadKmers){//All kmer processing is outside of the sync block
				processTileKmers(r1, kmerDepths0, baseDepths0);
				processTileKmers(r2, kmerDepths1, baseDepths1);
				
				IntList list0=(deblurDepths ? baseDepths0 : kmerDepths0);
//				long[] counts0=laneDepthCounts[lnum][0];
//				long[] sums0=laneDepthSums[lnum][0];
				for(int i=0; i<list0.size; i++) {
					int d=list0.get(i);
					int hit=(d>=cutoff ? 1 : 0);
					hits1+=hit;
					misses+=(hit^1);//This is clever.  No conditionals!
					depthSum+=d;
					//TODO: This atomic increment is super slow; replace it with locals.
//					lane.depthSums[0].addAndGet(i, d);
//					lane.depthCounts[0].incrementAndGet(i);
//					counts0[i]++;
//					sums0[i]+=d;
				}
//				assert(sums0[0]==0) : Arrays.toString(counts0)+"\n"+Arrays.toString(sums0);
				IntList list1=(deblurDepths ? baseDepths1 : kmerDepths1);
//				long[] counts1=laneDepthCounts[lnum][1];
//				long[] sums1=laneDepthSums[lnum][1];
				for(int i=0; i<list1.size; i++) {
					int d=list1.get(i);
					int hit=(d>=cutoff ? 1 : 0);
					hits2+=hit;
					misses+=(hit^1);
					depthSum+=d;
//					lane.depthSums[1].addAndGet(i, d);
//					lane.depthCounts[1].incrementAndGet(i);
//					counts1[i]++;
//					sums1[i]+=d;
				}
			}
			
			int merged=0;
			int insert=0;
			int overlap=0;
			int mergeErrors=0;
			if(merge && r2!=null) {
				if(strictmerge) {
					insert=BBMerge.findOverlapStrict(r1, r2, false);
				}else {
					insert=BBMerge.findOverlapLoose(r1, r2, false);
				}
				if(insert>0) {
					merged=2;
					overlap=Tools.min(insert, r1.length()+r2.length()-insert);
					mergeErrors=BBMerge.countErrors(r1.bases, r2.bases, insert);
				}else {insert=0;}
			}
			
			if(sidechannel!=null) {
				sidechannel.map(r1, r2);
			}
			
			int bchdist=0;
			int barcodePolymers=0;
			if(barcodeStats!=null) {
				String code=ihp.barcode();
				if(!barcodeStats.expectedCodeList.isEmpty()) {
					bchdist=barcodeStats.calcHdist(code);
				}
				barcodePolymers=Barcode.countPolymers(code);
//				assert(barcodePolymers==0 && !code.contains("GGGGGGGG")) : code+", "+barcodePolymers;
				if(barcodeMap!=null) {
					AtomicStringNum count=barcodeMap.get(code);
					if(count==null) {
						count=new AtomicStringNum(code, 0);
						barcodeMap.put(code, count);
					}
					count.increment();
				}
			}
			
			//Changes hits to misses if the read was a poly-G read.
			//Unwise for short-insert libraries.
			//Probably not necessary, either
			if(changePolyGHitsToMisses && r1!=null && r1.discarded()) {misses+=hits1; hits1=0;}
			if(changePolyGHitsToMisses && r2!=null && r2.discarded()) {misses+=hits2; hits2=0;}
			
			final MicroTile mt;
			synchronized(flowcellT) {
				mt=flowcellT.getMicroTile(lnum, tile, x, y, true);
			}
			
//			boolean addQuick=true; //This has a very minor impact on speed, <5% at t=64
//			double readQualityByProbSum=0, probErrorFreeSum=0, baseErrorProbSum=0;
//			int alignedReadCount=0, alignedBaseCount=0, readErrorCount=0, baseErrorCount=0, readInsCount=0, readDelCount=0;
//			if(addQuick){//Moving things out of the synchronized block to improve threading; changes add to addQuick
//				if(r1!=null) {
//					int len=r1.length();
//					readQualityByProbSum+=r1.avgQualityByProbabilityDouble(true, len);
//					probErrorFreeSum+=100*r1.probabilityErrorFree(true, len);
//					baseErrorProbSum+=r1.expectedErrors(true, len);
//					
//					if(r1.match!=null) {
//						int bc=r1.countAlignedBases();
//						if(bc>0) {
//							alignedReadCount++;
//							alignedBaseCount+=bc;
//							int errors=r1.countErrors();
//							readErrorCount+=(errors>0 ? 1 : 0);
//							baseErrorCount+=errors;
//							int[] mSCNID=Read.countMatchEvents(r1.match);
//							readInsCount+=(mSCNID[4]>0 ? 1 : 0);
//							readDelCount+=(mSCNID[5]>0 ? 1 : 0);
//						}
//					}else if(r1.samline!=null && r1.samline.mapped()) {
//						alignedReadCount++;
//					}
//					
//				}
//				if(r2!=null) {
//					int len=r2.length();
//					readQualityByProbSum+=r2.avgQualityByProbabilityDouble(true, len);
//					probErrorFreeSum+=100*r2.probabilityErrorFree(true, len);
//					baseErrorProbSum+=r2.expectedErrors(true, len);
//					
//					if(r2.match!=null) {
//						int bc=r2.countAlignedBases();
//						if(bc>0) {
//							alignedReadCount++;
//							alignedBaseCount+=bc;
//							int errors=r2.countErrors();
//							readErrorCount+=(errors>0 ? 1 : 0);
//							baseErrorCount+=errors;
//							int[] mSCNID=Read.countMatchEvents(r2.match);
//							readInsCount+=(mSCNID[4]>0 ? 1 : 0);
//							readDelCount+=(mSCNID[5]>0 ? 1 : 0);
//						}
//					}else if(r2.samline!=null && r2.samline.mapped()) {
//						alignedReadCount++;
//					}
//				}
//			}
			
			synchronized(mt) {
//				if(addQuick) {
//					mt.addQuick(r1);
//					mt.addQuick(r2);
//
//					mt.readQualityByProbSum+=readQualityByProbSum;
//					mt.probErrorFreeSum+=probErrorFreeSum;
//					mt.baseErrorProbSum+=baseErrorProbSum;
//
//					mt.alignedReadCount+=alignedReadCount;
//					mt.alignedBaseCount+=alignedBaseCount;
//					mt.readErrorCount+=readErrorCount;
//					mt.baseErrorCount+=baseErrorCount;
//					mt.readInsCount+=readInsCount;
//					mt.readDelCount+=readDelCount;
//				}else {
					mt.add(r1);
					mt.add(r2);
//				}
				
				mt.hits+=(hits1+hits2);
				mt.misses+=misses;
				mt.depthSum+=depthSum;
				
				mt.barcodes+=barcodesPerRead;
				mt.barcodeHDistSum+=bchdist;
				mt.barcodePolymers+=barcodePolymers;
				
				mt.mergedReads+=merged;
				mt.insertSum+=insert;
				mt.overlapSum+=overlap;
				mt.mergeErrorSum+=mergeErrors;
			}
		}
		
		private void processSam() {
			ListNum<SamLine> ln=ss.nextLines();
			ArrayList<SamLine> reads=(ln==null ? null : ln.list);
			final IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();

			while(ln!=null && reads!=null && reads.size()>0){

				for(int idx=0; idx<reads.size(); idx++){
					SamLine sl=reads.get(idx);
					processSamLine(sl, ihp);
				}
				ln=ss.nextLines();
				reads=(ln==null ? null : ln.list);
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;

		int xmin=Integer.MAX_VALUE, ymin=Integer.MAX_VALUE, tmin=Integer.MAX_VALUE;
		int xmax=-1, ymax=-1, tmax=-1;
		public String flowcellName=null;
		
		private IntList kmerDepths0=new IntList(151);
		private IntList kmerDepths1=new IntList(151);
		private IntList baseDepths0=new IntList(151);
		private IntList baseDepths1=new IntList(151);

//		long[][][] laneDepthSums=new long[9][2][500];
//		long[][][] laneDepthCounts=new long[9][2][500];
		
		private IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
		
		/** True only if this thread has completed successfully */
		boolean success=false;

		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Optional sam input stream */
		private final SamLineStreamer ss;
		/** Thread ID */
		final int tid;
		final FlowCell flowcellT;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private ArrayList<String> extra=new ArrayList<String>();

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;

	/** Discard output file path */
	private String outbad=null;

	/** Optional aligned reads (e.g. PhiX) */
	private String samInput=null;
	private final boolean processSamMT=true;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/*--------------------------------------------------------------*/	
	
//	private boolean pound=true;
	private String dumpOut=null;
	private String dumpIn=null;
	private String coordsOut=null;
	
	/*--------------------------------------------------------------*/	
	
	private byte delimiter='+';
	private int barcodesPerRead=2;
	private String expectedBarcodes;
	private String barcodeCounts;
	private BarcodeStats barcodeStats;
	
	private ConcurrentHashMap<String, AtomicStringNum> barcodeMap;
	
	/*--------------------------------------------------------------*/
	/*----------------         Side Channel         ----------------*/
	/*--------------------------------------------------------------*/
	
	boolean align=false;
	String alignOut=null;
	String alignRef="phix";
	float alignMinid1=0.66f;
	float alignMinid2=0.56f;
	int alignK1=17; //Phix is unique down to k=13 solid; unsure about gapped
	int alignK2=13;
	int alignMM1=1;
	int alignMM2=1;
	SideChannel3 sidechannel;
	
	/*--------------------------------------------------------------*/	

	/** Number of reads processed */
	public long readsProcessed=0;
	/** Number of bases processed */
	public long basesProcessed=0;

	/** Number of reads discarded */
	public long readsDiscarded=0;
	/** Number of bases discarded */
	public long basesDiscarded=0;
	
	protected long gsTransformedToN=0; 
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;

	//Test this; may need to be higher.  Was 16, now 24.
	//Seems unused though, since switching to a mask...
	//And Perlmutter seems to saturate at ~2000% CPU
	private int loadThreads=24;
	//32 was not quite enough on Perlmutter; maybe OK on Dori though
	//64 is not enough with merge enabled
	private int fillThreads=64;
	//128 threads only ran at 5600% utilization...  may be too many.
	private static final int fillThreadsM=96;//threads when merge is enabled
	private BloomFilter bloomFilter;

	int idmask_read=7;
	int idmask_write=15;
	private int smoothDepths=0;
	private boolean deblurDepths=false;
	
	private boolean blurTiles=false;
	
	private boolean loadKmers=true;
//	private int kmersPerRead=0;//for loading
	private float minProb=0;
	private boolean deterministic=true;
	private int cbits=2;
	private int hashes=2;
	
	private boolean recalibrate=false;
	private boolean merge=false;
	private boolean strictmerge=false;
	
	private int targetAverageReads=1600;
	private int targetAlignedReads=250;
	private int targetX=Tile.xSize;
	private int targetY=Tile.ySize;
	
	private static final int k=31;
	
	private FlowCell flowcell;
	
	private boolean changePolyGHitsToMisses=false;
	
	private float dmult=-0.2f;
	
	private boolean discardOnlyLowQuality=true;
	private int discardLevel=1;
	private boolean gToN=false;
	private boolean discardG=false;
	
	private int minlen=30;
	private float trimq=-1;
	private final float trimE;
	private boolean trimLeft=false;
	private boolean trimRight=true;
	
	private boolean warned=false;
	
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
	
	/** Output for discarded reads */
	private final FileFormat ffoutbad;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of reads output in the last run */
	public static long lastReadsOut;
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
	/** Output reads in input order (possibly todo) */
	private final boolean ordered=false;
	
}
