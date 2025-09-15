package barcode;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map.Entry;
import java.util.Set;
import java.util.regex.Pattern;

import aligner.MicroAligner2;
import barcode.stub.PCRMatrixProb;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import hiseq.IlluminaHeaderParser2;
import shared.KillSwitch;
import shared.LineParserS2;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.BufferedMultiCros;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;


/**
 * This class is designed to demultiplex barcoded reads using PCRMatrix.
 * 
 * @author Brian Bushnell
 * @date April 2, 2024
 *
 */
public class NovaDemux {

	/** Code entrance from the command line */
	public static void main(String[] args){
		
		//Capture values of static variables that might be modified in case this is called by another class.
		final int oldCap=Shared.numBuffers(), oldZipThreads=ReadWrite.MAX_ZIP_THREADS(), oldZl=ReadWrite.ZIPLEVEL;
		
		//External compressor selection
		final boolean oldPigz=ReadWrite.USE_PIGZ, oldUnpigz=ReadWrite.USE_UNPIGZ;
		final boolean oldBgzip=ReadWrite.USE_BGZIP, oldPreferBgzip=ReadWrite.PREFER_BGZIP;
		
		//Preserve quality scores even for reads with with incorrect quality scores
		final boolean oldCQ=Read.CHANGE_QUALITY;
		
		//Create and start a timer, for statistics
		Timer t=new Timer();
		
		//Create the demultiplexer object
		NovaDemux demultiplexer=new NovaDemux(args);
		
		//Process all data
		try {
			demultiplexer.process(t);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
		
		//Restore values of static variables.
		Shared.setBuffers(oldCap);
		ReadWrite.ZIPLEVEL=oldZl;
		ReadWrite.USE_PIGZ=oldPigz;
		ReadWrite.USE_BGZIP=oldBgzip;
		ReadWrite.PREFER_BGZIP=oldPreferBgzip;
		ReadWrite.USE_UNPIGZ=oldUnpigz;
		ReadWrite.setZipThreads(oldZipThreads);
		Read.CHANGE_QUALITY=oldCQ;
		
		//Close the print stream if it was redirected
		Shared.closeStream(demultiplexer.outstream);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public NovaDemux(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set some static variables
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=true;
		ReadWrite.USE_UNPIGZ=true;
		SamLine.SET_FROM_OK=true;
		
		//Reduce default compression level to increase speed
		//This automatically gets bumped up to 4 if pigz/bgzip is detected
		ReadWrite.ZIPLEVEL=2;
		Read.CHANGE_QUALITY=false;
		
		
		{//Argument-parsing block
			
			//Parse arguments
			final Parser parser=parse(args);
			
			//Process parser fields
			//This is in the constructor rather than a function to allow final fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			trackCardinality=parser.loglog;
			
			in1=parser.in1;
			in2=parser.in2;
//			qfin1=parser.qfin1;
//			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
//			qfout1=parser.qfout1;
//			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		PCRMatrix.postParseStatic();
		
		if(!ReadWrite.SET_ZIP_THREADS) {
			ReadWrite.setZipThreads(Shared.threads());
			if(ReadWrite.MAX_ZIP_THREADS()>=8) {
				float numerator=(ReadWrite.ZIPLEVEL>7 ? 32f : 8f);
				ReadWrite.setZipThreadMult(numerator/ReadWrite.MAX_ZIP_THREADS());
			}
		}
		
		//Ensure the input and output files are specified correctly
		validate();

		//Create FileFormat objects for input files
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		barcodeDelimiter=ffin1.barcodeDelimiter();
		barcodeLength1=ffin1.barcodeLength(1);
		barcodeLength2=ffin1.barcodeLength(2);
		
		//If the input and output are in sam/bam format, print the sam header
		if(ffin1!=null && out1!=null && ffin1.samOrBam()){
			useSharedHeader=FileFormat.isSamOrBamFile(out1);
		}
		
		//Set up data structures
		processHeaderDelimiter();
		
		if(expectedSet==null && sampleMapFile!=null) {
			expectedSet=new LinkedHashSet<String>();
			expectedSet.add(sampleMapFile);//This should work due to stripping of tabs.
		}
		assert(expectedSet!=null && expectedSet.size()>0) : "Expected barcodes are a required parameter.";
		BarcodeStats.loadBarcodeSet(expectedSet, (byte)barcodeDelimiter, rcIndex1, rcIndex2);
		BarcodeStats.loadBarcodeSet(outSubset, (byte)barcodeDelimiter, rcIndex1, rcIndex2);
		
//		//TODO: Activate this if stat dumps are needed
//		readstats=ReadStats.collectingStats() ? new ReadStats() : null;
		
		legacyWriter=(legacyPath==null ? null : new LegacyFileWriter());
	}
	
	/** 
	 * Parse command-line arguments
	 * @param args Command-line arguments
	 * @return A Parser object with fields filled from the parameters
	 */
	private Parser parse(final String[] args){
		
		//Handles parsing of common flags
		Parser parser=new Parser();
		parser.overwrite=overwrite;
		
		//Handles parsing of flags specific to this program
		for(int i=0; i<args.length; i++){//Parsing loop
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				BufferedMultiCros.verbose=verbose;
			}else if(a.equals("verboseclient")){
				verboseClient=Parse.parseBoolean(b);
			}else if(a.equals("verbosemcros")){
				BufferedMultiCros.verbose=Parse.parseBoolean(b);
			}else if(a.equals("names") || a.equals("expected") || a.equals("barcodes")){
				expectedSet=new LinkedHashSet<String>();
				for(String s : b.split(",")){
					expectedSet.add(s);
				}
			}else if(a.equals("subset") || a.equals("outset") || a.equals("outsubset")){
				outSubset=new LinkedHashSet<String>();
				for(String s : b.split(",")){
					outSubset.add(s);
				}
			}else if(a.equals("barcode") || a.equals("barcodemode") || a.equals("index")){
				if(Parse.parseBoolean(b)){mode=BARCODE_MODE;}
			}else if(a.equals("perheader") || a.equals("persequence") || a.equals("header") || a.equals("headermode")){
				if(Parse.parseBoolean(b)){mode=HEADER_MODE;}
			}else if(a.equals("suffixmode") || a.equals("suffix") || a.equals("sm")){
				if(Parse.parseBoolean(b)){mode=SUFFIX_MODE;}
			}else if(a.equals("prefixmode") || a.equals("prefix") || a.equals("pm")){
				if(Parse.parseBoolean(b)){mode=PREFIX_MODE;}
			}else if(a.equals("delimiter")){
				assert(false) : "Please specify headerdelimiter or barcodedelimiter instead of 'delimiter'";
			}else if(a.equals("hdelimiter") || a.equals("headerdelimiter")){
				headerDelimiter=Parse.parseSymbol(b);
			}else if(a.equals("bdelimiter") || a.equals("barcodedelimiter")){
				String s=Parse.parseSymbol(b);
				assert(s!=null && s.length()==1) : "Invalid barcode delimiter: "+b;
				barcodeDelimiter=s.charAt(0);
			}else if(a.equals("rc1") || a.equals("rcindex1")){
				rcIndex1=Parse.parseBoolean(b);
			}else if(a.equals("rc2") || a.equals("rcindex2")){
				rcIndex2=Parse.parseBoolean(b);
			}else if(a.equals("column")){
				column=Integer.parseInt(b);
				assert(column>0 || column==-1) : "Column is 1-based; must be 1+ or else -1 to disable.";
				column--;
			}else if(a.equalsIgnoreCase("length") || a.equalsIgnoreCase("len") || 
					a.equalsIgnoreCase("affixlength") || a.equalsIgnoreCase("affixlen")){
				fixedLength=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("stats") || a.equalsIgnoreCase("report") || a.equalsIgnoreCase("results")){
				stats=b;
			}else if(a.equalsIgnoreCase("minreadstodump") || a.equalsIgnoreCase("minreads")){
				minReadsToDump=Parse.parseKMG(b);
			}else if(a.equalsIgnoreCase("printRetireTime")){
				printRetireTime=Parse.parseBoolean(b);
			}else if(a.equals("outu") || a.equals("outu1")){
				outu1=b;
			}else if(a.equals("outu2")){
				outu2=b;
			}else if(a.equals("pattern")){
				parser.out1=b;
			}else if(a.equals("legacy") || a.equals("outlegacy") || a.equals("legacyout") || a.equals("legacypath")){
				legacyPath=b;
			}else if(a.equals("samplemapfile") || a.equals("samplemap")){
				sampleMapFile=b;
			}else if(a.equals("lane")){
				lane=Integer.parseInt(b);
			}else if(a.equals("statsonly")){
				statsOnly=Parse.parseBoolean(b);
			}else if(a.equals("remap") || a.equals("replace")){
				symbolRemap=Parse.parseRemap(b);
			}else if(a.equalsIgnoreCase("writeEmptyFiles") || a.equalsIgnoreCase("writeEmpty")) {
				writeEmptyFiles=Parse.parseBoolean(b);
			}else if(a.equals("countsin")){
				countsIn=b;
			}else if(a.equals("countsout")){
				countsOut=b;
			}else if(a.equals("barcodesout")){
				barcodesOut=b;
			}else if(a.equals("mincount0")){
				minCount0=Long.parseLong(b);
			}else if(a.equals("mincounta")){
				minCountA=Long.parseLong(b);
			}else if(a.equals("mincountr")){
				minCountR=Long.parseLong(b);
			}else if(a.equals("mapin")){
				mapIn=b;
			}else if(a.equals("mapout")){
				mapOut=b;
			}else if(a.equals("client") || a.equals("server") || a.equals("useserver")){
				if("auto".equalsIgnoreCase(b)) {
					useServer=false;
					setUseServer=false;
				}else{
					useServer=Parse.parseBoolean(b);
					setUseServer=true;
				}
			}
			
			else if(a.equals("rename")){
				rename=Parse.parseBoolean(b);
			}else if(a.equals("nosplit")){
				nosplit=Parse.parseBoolean(b);
			}
			

			else if(a.equals("refpath") || a.equals("spikepath") || a.equals("phixpath")){
				refPath=b;
			}else if(a.equals("spikelabel") || a.equals("phixlabel") 
					|| a.equals("spikeindex") || a.equals("phixindex")){
				spikeLabel=b;
			}else if(a.equals("kspike") || a.equals("kphix")){
				kSpike=Integer.parseInt(b);
			}else if(a.equals("skipmask")){
				skipMask=Integer.parseInt(b);
			}else if(a.equals("minid") || a.equals("minspikeid") || a.equals("phixid")){
				minSpikeIdentity=Float.parseFloat(b);
				if(minSpikeIdentity>1) {minSpikeIdentity/=100;}
			}else if(a.equalsIgnoreCase("mapUnexpected")){
				mapUnexpectedToSpike=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("mapUnknown")){
				mapUnexpectedToSpike=!Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("mapAll")){
				mapAllToSpike=Parse.parseBoolean(b);
			}
			
			else if(PCRMatrix.parseStatic(arg, a, b)){
				//Flag was captured by PCRMatrix; do nothing
			}else if(BufferedMultiCros.parseStatic(arg, a, b)){
				//Flag was captured by BufferedMultiCros; do nothing
			}else if(parser.parse(arg, a, b)){
				//Flag was captured by the parser; do nothing
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/** 
	 * Ensure input and output files are specified correctly,
	 * and that interleaved mode is consistent with the number of input files.
	 */
	private void validate(){
		
		//Validate patterns
		assert(out1==null || out1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(out2==null || out2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		
		//Perform # replacement for twin files
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(outu1!=null && outu2==null && outu1.indexOf('#')>-1){
			outu2=outu1.replace("#", "2");
			outu1=outu1.replace("#", "1");
		}

		//Disable interleaving if in2 is specified
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		{//Perform various file validation
			assert(FastaReadInputStream.settingsOK());
			if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
			if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
			
			//ByteFile2 is multithreaded; only use it if there are plenty of threads
			if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
				ByteFile.FORCE_MODE_BF2=true;
			}
			
			//Do some interleaving-detection logic if unspecified
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
			
			//If the user explicitly used the word "null", change it to null.
			if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
			if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
			
			//Ensure output files have no duplicates and can be written
			//Unfortunately, this does not prevent patterns expanding into a matching string
			if(!Tools.testOutputFiles(overwrite, append, false, outu1, outu2, stats)){
				outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2+", "+outu1+", "+outu2+", "+stats);
				throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+", "+outu1+", "+outu2+", "+stats+"\n");
			}

			//This does not prevent expansion into a matching string
			assert(out1==null || (!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1))) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
		}
	}
	
	private void processHeaderDelimiter() {
		//For length-1 delimiters, use a char instead of a String
		assert(column<0 || headerDelimiter!=null) : "Column may not be set if there is no delimiter.";
		if(headerDelimiter!=null){
			if(headerDelimiter.length()==1){
				delimiterChar=headerDelimiter.charAt(0);
			}else if(headerDelimiter.length()==2 && headerDelimiter.charAt(0)=='\\'){
				delimiterChar=headerDelimiter.charAt(1);
			}else{
				delimiterChar=0;
			}
			delimiterPattern=Pattern.compile(headerDelimiter);
			lp=new LineParserS2(delimiterChar);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private HashMap<String, String> makeAssignmentMap() {
		if(mapIn!=null) {return BarcodeCounter.loadAssignmentMap(mapIn);}
		
		Timer t=new Timer();
		final Collection<Barcode> counts;
		
		System.err.print("Counting Barcodes:\t");
		if(countsIn!=null) {counts=BarcodeCounter.loadCounts(countsIn, minCount0);}
		else {
			HashMap<String, Barcode> codeMap=BarcodeCounter.countBarcodes(ffin1, maxReads, PCRMatrix.byTile);
			counts=codeMap.values();
			if(countsOut!=null) {BarcodeCounter.writeCounts(counts, minCount0, countsOut, true, overwrite, append);}
		}
		if(legacyWriter!=null) {legacyWriter.counts=counts;}
		t.stop("");
		t.start();
		final HashMap<String, String> map;
		
		if(PCRMatrix.matrixType0==PCRMatrix.PROB_TYPE && PCRMatrixProb.clientside() && 
				!setUseServer && !useServer) {
			useServer=true;
		}
		if(useServer) {
			System.err.println("Using client-server mode for barcode analysis.");
			DemuxData dd=new DemuxData(barcodeLength1, barcodeLength2, barcodeDelimiter);
			dd.codeCounts=counts;
			dd.expectedList=expectedSet;
			DemuxClient client=new DemuxClient();
			map=client.getMap(dd, verboseClient);
			if(mapOut!=null) {
				PCRMatrix.printAssignmentMapStatic(map, mapOut, counts, overwrite, append);
			}
		}else {
			PCRMatrix pcrMatrix=PCRMatrix.create(barcodeLength1, barcodeLength2, barcodeDelimiter);
			assert(expectedSet!=null && !expectedSet.isEmpty());
			pcrMatrix.populateExpected(expectedSet);
			pcrMatrix.populateSplitCodes();
			pcrMatrix.initializeData();
			pcrMatrix.refine(counts, minCountR);
			map=pcrMatrix.makeAssignmentMap(counts, minCountA);
			if(mapOut!=null) {
				pcrMatrix.printAssignmentMap(map, mapOut, counts, overwrite, append);
			}
			pcrMatrix=null;
		}
		
		if(barcodesOut!=null) {
			ByteStreamWriter bsw=new ByteStreamWriter(barcodesOut, overwrite, append, true);
			bsw.start();
			ArrayList<Barcode> list=Barcode.summateAssignments(map, expectedSet, counts);
			for(Barcode b : list) {
				bsw.print(b.name).tab().print(b.frequency, 5).tab().print(b.count()).nl();
			}
			bsw.poisonAndWait();
		}
		
		t.stop("Assignment Time:\t\t");
//		if(mapOut!=null) {BarcodeCounter.writeAssignmentMap(map, mapOut, overwrite, append);}
		
//		map=filterAssignmentMap(map, outSubset);//Sends known things to unknown
		
		return map;
	}
	
	private HashMap<String, String> filterAssignmentMap(HashMap<String, String> map, Set<String> set) {
		if(map==null || set==null || map.isEmpty() || set.isEmpty()) {return map;}
		HashMap<String, String> filtered=new HashMap<String, String>();
		for(Entry<String, String> e : map.entrySet()) {
			if(set.contains(e.getValue())) {
				filtered.put(e.getKey(), e.getValue());
			}
		}
		return filtered;
	}
	
	/** 
	 * Primary method.
	 * Starts and stops I/O streams, processes the data, and prints results.
	 * @param t A timer that was already started.
	 */
	void process(Timer t){
		
		assignmentMap=makeAssignmentMap();
		
		Timer t2=new Timer();
		
		//Create stream for input reads
		final ConcurrentReadInputStream cris=makeInputStream();
		
		//Create streams for output reads other than unmatched
		final BufferedMultiCros mcros=(nosplit ? null : makeMatchedOutputStream(cris.paired()));
		
		//Create stream for unmatched output reads
		//TODO: Consider adding this to mcros
		final ConcurrentReadOutputStream rosu=makeUnmatchedOutputStream();
		
		spikeMapper=(refPath==null || spikeLabel==null ? null : 
			new MicroAligner2(kSpike, minSpikeIdentity, refPath));
		if(spikeMapper!=null) {spikeMapper.skipmask=skipMask;}//Only lookup some kmers for speed
		
		//Streams are set up, so process the reads
		processInner(cris, mcros, rosu);
		//At this point processing has finished.
		
		//Close streams
		cleanup(cris, mcros, rosu);
		
		t2.stop("Writing Time:\t\t\t");
		t2.start();
		
		//Report statistics to file
		if(stats!=null){printReport(mcros);}
		
		if(legacyPath!=null) {
			LinkedHashMap<String, String> sampleMap=legacyWriter.loadSampleMap(sampleMapFile, 
					(byte)barcodeDelimiter, rcIndex1, rcIndex2);
			if(legacyPath.equals(".")) {legacyPath="";}
			if(legacyPath.length()>0 && !legacyPath.endsWith("/")) {
				legacyPath+="/";
			}

			if(sampleMap==null || sampleMap.isEmpty()) {
				System.err.println("Warning: sample map is empty; "
						+ "sample names will be absent in legacy files.");
			}
			legacyWriter.writeTopUnknownBarcodes(expectedSet, assignmentMap, 
					(byte)barcodeDelimiter, 
					legacyPath+"Top_Unknown_Barcodes.csv", lane, 1000, overwrite);
			legacyWriter.writeQualityMetrics(expectedSet, sampleMap, 
					(byte)barcodeDelimiter, legacyPath+"Quality_Metrics.csv", lane, overwrite);
			legacyWriter.writeDemultiplexStats(legacyPath+"Demultiplex_Stats.csv", assignmentMap, 
					expectedSet, sampleMap, lane, overwrite);
			legacyWriter.writeIndexHoppingCounts(expectedSet, sampleMap, (byte)barcodeDelimiter,
					legacyPath+"Index_Hopping_Counts.csv", lane, overwrite);
			legacyWriter.writeDemultiplexTileStats(legacyPath+"Demultiplex_Tile_Stats.csv", 
					sampleMap, (sampleMap==null ? expectedSet : null), lane, (byte)barcodeDelimiter, overwrite);
		}
		if(stats!=null || legacyPath!=null) {
			t2.stop("Stats Writing Time:\t\t");
		}

		if(spikeMapper!=null) {
			System.err.println("Mapping Attempts:\t"+spikeMapper.mapCount);
			System.err.println("Quick Aligns:    \t"+spikeMapper.quickAligns);
			System.err.println("Slow Aligns:     \t"+spikeMapper.slowAligns);
			System.err.println("Align Success:   \t"+spikeMapper.metCutoff);
			double avgId=spikeMapper.idSum/(Tools.max(spikeMapper.metCutoff,1));
			System.err.println("Avg Identity:    \t"+String.format("%.6f", avgId));
		}
		
		if(printRetireTime && mcros!=null) {
			outstream.println("\n"+mcros.printRetireTime()+mcros.printCreateTime()+"\n");
		}
		
		//Stop the timer
		t.stop();
		
		//Print results
		printResultsToScreen(t);
		
		if(errorState){//Exit with a message and error code
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Create and start the input stream (cris). */
	private ConcurrentReadInputStream makeInputStream(){
		//Stream for input reads
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ffin1, ffin2);
		if(verbose){outstream.println("Started cris");}
		cris.start();
		
		//Report whether the input is being processed as paired.
		final boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}
	
	/** Create and start the normal output streams (mcros). */
	private BufferedMultiCros makeMatchedOutputStream(boolean pairedInput){
		//Streams for output reads other than unmatched
		final BufferedMultiCros mcros;
		if(out1!=null){
			mcros=BufferedMultiCros.make(out1, out2, overwrite, append,
					true, useSharedHeader, FileFormat.FASTQ);
			
			//Set mcros fields not available in constructor
			mcros.minReadsToDump=minReadsToDump;
			mcros.trackCardinality=trackCardinality;
			
			//Report output format
			if(pairedInput && out2==null && (in1==null || !ffin1.samOrBam())){
				outstream.println("Writing interleaved.");
			}
			
			//Start the mcros thread if necessary
			if(mcros.threaded){mcros.start();}
		}else{
			mcros=null;
		}
		return mcros;
	}
		
	/** Create and start the unmatched output stream (rosu). */
	private ConcurrentReadOutputStream makeUnmatchedOutputStream(){
		//Stream for unmatched output reads
		final ConcurrentReadOutputStream rosu;
		if(outu1!=null){
			//Number of output buffers; does not need to be high since access is single-threaded
			final int buff=4;
			
			FileFormat ffout1=FileFormat.testOutput(outu1, FileFormat.FASTQ, extout, true, overwrite, append, false);
			FileFormat ffout2=(outu2==null ? null : FileFormat.testOutput(outu2, FileFormat.FASTQ, extout, true, overwrite, append, false));
			rosu=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, true);
			rosu.start();
		}else{
			rosu=null;
		}
		
		return rosu;
	}
	
	/** 
	 * Process all reads.
	 * @param cris Input stream (required).
	 * @param mcros Matched read output stream (optional).
	 * @param rosu Unmatched read output stream (optional).
	 */
	private void processInner(final ConcurrentReadInputStream cris, final BufferedMultiCros mcros, final ConcurrentReadOutputStream rosu) {
		
		//Fetch the first list
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		//Ensure pairing seems correct
		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
		}
		
		//While there are more reads to process...
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			//List to accumulate unmatched reads
			ArrayList<Read> unmatched=new ArrayList<Read>();
			
			//For each read...
			for(Read r : reads){
				
				//Get the target file identifier
				String name=getValue(r);
				final int pairCount=r.pairCount(), pairLen=r.pairLength();
				if(legacyWriter!=null) {
					r.obj=name;
					legacyWriter.add(r, name);
					legacyWriter.add(r.mate, name);
				}
				if(name!=null){
					//Set the name so that the mcros will send it to the right place
					if(outSubset==null || outSubset.contains(name)) {
						r.obj=(symbolRemap==null ? name : Tools.remap(symbolRemap, name));
						readsOut+=pairCount;
						basesOut+=pairLen;
					}
				}else{
					//Send it to unmatched
					unmatched.add(r);
					readsUnmatched+=pairCount;
					basesUnmatched+=pairLen;
				}
				if(rename) {rename(r, name);}

				readsProcessed+=pairCount;
				basesProcessed+=pairLen;
			}
			if(rosu!=null){rosu.add(nosplit ? reads : unmatched, ln.id);}//Send unmatched reads to rosu
			if(mcros!=null){mcros.add(reads);}//Send matched reads to mcros
			
			//Notify the input stream that the list has been processed
			cris.returnList(ln);
			
			//Fetch a new list
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		
		//Notify the input stream that this thread is done processing
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}
	
	private void rename(Read r, String name) {
		if(name==null) {name="UNKNOWN";}
		r.id=r.id+"\t"+name;
		if(r.mate!=null) {r.mate.id=r.mate.id+"\t"+name;}
	}
	
	/** 
	 * Finish writing and close streams.
	 * @param cris Input stream (required).
	 * @param mcros Matched read output stream (optional).
	 * @param rosu Unmatched read output stream (optional).
	 */
	private void cleanup(final ConcurrentReadInputStream cris, final BufferedMultiCros mcros, 
			final ConcurrentReadOutputStream rosu){
		
		//Shut down mcros
		if(mcros!=null){
			mcros.close();
			errorState|=mcros.errorState();
		}

		if(minReadsToDump>0 && mcros!=null){
			//Dump the residual reads into the unmatched file
			mcros.dumpResidual(rosu);

			//Adjust statistics to reflect that residuals did not get demultiplexed
			readsOut-=mcros.residualReads;
			basesOut-=mcros.residualBases;
		}
		
		if(mcros!=null && writeEmptyFiles) {
			Set<String> observed=mcros.getKeys();
			for(String expected : expectedSet) {
				String key=(symbolRemap==null ? expected : Tools.remap(symbolRemap, expected));
				if(!observed.contains(key)) {

					String s1=out1.replaceFirst("%", key);
					String s2=out2==null ? null : out2.replaceFirst("%", key);
					ReadWrite.writeString("", s1, overwrite, false);
					if(s2!=null) {ReadWrite.writeString("", s2, overwrite, false);}
				}
			}
		}

		//errorState|=ReadStats.writeAll(); //Currently unused; allows tracking statistics.
		
		//Shut down normal streams
		errorState|=ReadWrite.closeStreams(cris, rosu);
	}
	
	/**
	 * @param t A Timer that has already been stopped
	 */
	private void printResultsToScreen(Timer t) {
		//Report statistics to screen
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);
		
		outstream.println("Time:               "+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+basesProcessed+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
		outstream.println("Reads Out:          "+readsOut);
		outstream.println("Bases Out:          "+basesOut);
		outstream.println("Yield:              "+String.format("%.5f", readsOut*1.0/readsProcessed));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns the value this read maps to; meaning, the variable part of the output filename */
	private String getValue(Read r){
		String key=getKey(r);
//		assert(key==null) : key;
//		if(verbose){System.err.println("Got key "+key);}
		String value=(key==null ? null : assignmentMap.get(key));
		//Note - generally...  unexpected are not in the assignment map anyway.
		if(spikeMapper!=null && (mapAllToSpike || value==null || 
				(mapUnexpectedToSpike && !expectedSet.contains(value)))) {
			float id=spikeMapper.map(r);
			if(id>=minSpikeIdentity) {value=spikeLabel;}
		}
		return value;
	}
	
	/** Generates a key from the read header */
	private String getKey(Read r){
		final String id=r.id;
		final int idlen=id.length();
		final String key;
		
		//Simplest case; every header is its own key
		if(mode==HEADER_MODE){return id;}
		
		//Parse the barcode
		if(mode==BARCODE_MODE){
			ihp.parse(id);
			String bc=ihp.barcode();
			return PCRMatrix.byTile ? bc+ihp.tile() : bc;
		}
		
		//Return a prefix or suffix
		if(mode==PREFIX_MODE || mode==SUFFIX_MODE){
			assert(fixedLength>0);
			key=(id.length()<=fixedLength ? id : mode==PREFIX_MODE ? 
					id.substring(0, fixedLength) : id.substring(idlen-fixedLength));
			return key;
		}
		
		//Return a delimited substring
		if(mode==DELIMITER_MODE){
			
			if(delimiterChar>0){//Use LineParser
				lp.set(id);
//				int col=Tools.min(column, lp.terms()-1); //TODO: Make LineParserS
				key=lp.parseString(column);
			}else{//Use a split operation, which is slow
				String[] split=delimiterPattern.split(id);//Faster; uses a precompiled pattern
				assert(split.length>1) : "Delimiter '"+headerDelimiter+"' was not found in name '"+id+"'";
				
				int col=Tools.min(column, split.length-1);
				key=split[col];
				if(col!=column && !warned){
					outstream.println("*** WARNING! ***\n"
							+ "Only "+(col+1)+" columns for record "+id+"\n"
							+ "Further warnings will be suppressed.\n");
					warned=true;
					assert(errorState=true); //Change error state to true if assertions are enabled.
				}
			}
			return key;
		}
		
		throw new RuntimeException("No mode is being used: "+mode); //Should be unreachable
	}
	
	/** 
	 * Print statistics about demultiplexing to a file.
	 * @param mcros Output stream, after processing is completely finished.
	 */
	void printReport(BufferedMultiCros mcros){
		if(stats==null){return;}
		
		//Make a writer for the stats file
		ByteStreamWriter bsw=new ByteStreamWriter(stats, overwrite, append, true);
		bsw.start();
		
		{
			ByteBuilder bb=new ByteBuilder();

			//Print the header
			bb.append("#ReadsIn\t").append(readsProcessed).nl();
			bb.append("#BasesIn\t").append(basesProcessed).nl();
			bb.append("#ReadsOut\t").append(readsOut).nl();
			bb.append("#BasesOut\t").append(basesOut).nl();
			bb.append("#Name\tReads\tBases"+(trackCardinality ? "\tCardinality" : "")+"\n");
			
			//Print results of unmatched reads, which is not captured by the mcros
			if(assignmentMap!=null && !assignmentMap.isEmpty()){
				bb.append("Unmatched\t").append(readsUnmatched).tab().append(basesUnmatched).nl();
			}
			bsw.print(bb);
		}
		
		if(mcros!=null){//Print results from the mcros
			ByteBuilder bb=mcros.report();
			bsw.print(bb);
		}
		
		//Finish writing
		bsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file */
	private String in1=null;
	/** Optional input file for read 2 of paired reads in twin files */
	private String in2=null;

	/** Output file pattern 1 */
	private String out1=null;
	/** Output file pattern 2 */
	private String out2=null;
	
	/** Unmatched read output file 1 */
	private String outu1=null;
	/** Unmatched read output file 2 */
	private String outu2=null;
	
	/** File extension override for input files */
	private String extin=null;
	/** File extension override for output files */
	private String extout=null;
	
	/** Pre-counted barcodes */
	private String countsIn=null;
	/** Output for barcode counts */
	private String countsOut=null;
	
	/** Output for barcode assignment counts */
	private String barcodesOut=null;
	
	private boolean rename=false;
	private boolean nosplit=false;
	
	/** 
	 * Ignore barcodes occuring fewer times than this.
	 * Saves memory and speeds processing. */
	private long minCount0=0;
	
	/** Use stricter cutoffs during refinement for barcodes with lower count. */
	private long minCountR=4;
	
	/** Use stricter cutoffs during assignment for barcodes with lower count. */
	private long minCountA=4;
	
	/** Premade assignment map */
	private String mapIn=null;
	/** Assignment map output */
	private String mapOut=null;

	private boolean useServer=false;
	private boolean setUseServer=false;
	
	//This is for a spike-in like PhiX
	private String spikeLabel=null;
	private String refPath="phix";
	private MicroAligner2 spikeMapper;
	private int kSpike=27;
	private int skipMask=3;
	float minSpikeIdentity=0.7f;
	boolean mapUnexpectedToSpike=false;
	boolean mapAllToSpike=false;
	
	/*--------------------------------------------------------------*/
	
	/** Primary input file */
	private final FileFormat ffin1;
	
	/** Read 2 input file */
	private final FileFormat ffin2;
	
	/*--------------------------------------------------------------*/
	
	/** Input reads */
	long readsProcessed=0;
	long basesProcessed=0;
	
	/** Demultiplexed output reads */
	long readsOut=0;
	long basesOut=0;
	
	/** Output reads that did not get demultiplexed */
	long readsUnmatched=0;
	long basesUnmatched=0;

	/** Stop after this many input reads */
	private long maxReads=-1;

	/** File to print number of reads sent to each output file */
	private String stats=null;

	/** For splitting headers on a symbol */
	private String headerDelimiter=null;
	
	/** For splitting headers faster if the delimiter is just one character */
	private char delimiterChar=0;
	
	private LineParserS2 lp=null;
	
	/** Precompiled pattern matching the delimiter */
	private Pattern delimiterPattern=null;
	
	/** For replacing symbols in the filename */
	private byte[] symbolRemap=null;
	
	/** If there is a delimiter, use this column after splitting.
	 * Column is 1-based for parsing parameters but 0-based internally. */
	private int column=-1;
	
	/** How to select a read's key; see static modes array */
	private int mode=BARCODE_MODE;

	/** Track per-output-file cardinality; may be slow */
	private boolean trackCardinality=false;
	
	/** Prevents issuing warnings multiple times */
	private boolean warned=false;
	
	private boolean printRetireTime=false;
	
	/** 
	 * Do not create files with under this many reads.
	 * Increases memory usage since reads must be retained until processing is finished.
	 * Upon running out of memory the program may go very slowly.
	 */
	private long minReadsToDump=0;
	
	/** Affix length if fixed */
	private int fixedLength=-1;
	
	/** 
	 * A set of recognized names and the output names they map to.
	 * The keys and values will be identical unless a Hamming distance is used.
	 * This is only filled if explicit names are provided.
	 */
	private HashMap<String, String> assignmentMap=null;//new HashMap<String, String>();

	private LinkedHashSet<String> expectedSet=null;
	private LinkedHashSet<String> outSubset=null;
	
	/** Whether input file interleaving was explicitly set */
	private boolean setInterleaved=false;
	
	private IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();

	//Legacy file stuff
//	final ReadStats readstats;
	private int lane=0;
	private String legacyPath;
	private String sampleMapFile;
	private final LegacyFileWriter legacyWriter;
	
	/*--------------------------------------------------------------*/
	
	private int barcodeDelimiter=0;
	private int barcodeLength1=0;
	private int barcodeLength2=0;
	private boolean rcIndex1=false;
	private boolean rcIndex2=false;
	
	/*--------------------------------------------------------------*/
	
	/** Print messages here */
	private PrintStream outstream=System.err;
	
	/** True if errors were encountered */
	public boolean errorState=false;
	
	/** Create empty files for barcodes that are expected but absent. */
	private boolean writeEmptyFiles=true;
	
	/** Permission to overwrite existing files. */
	private boolean overwrite=true;
	
	/** 
	 * Append to existing files rather than overwriting.
	 * It is not advisable to set this flag for this class. 
	 * */
	private boolean append=false;
	
	/** Retain header of sam/bam files. */
	private boolean useSharedHeader=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Operation mode for selection of key */
	public static final int BARCODE_MODE=1, PREFIX_MODE=2, SUFFIX_MODE=3, DELIMITER_MODE=4, HEADER_MODE=5;
	
	private static final String empty="";
	
	/** Verbose messages for debugging */
	public static boolean verbose=false;
	public static boolean verboseClient=true;
	
	public static boolean statsOnly=false;
	
}
