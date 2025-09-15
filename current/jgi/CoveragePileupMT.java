package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import align2.QualityTools;
import dna.Data;
import dna.Scaffold;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.KillSwitch;
import shared.LineParser1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.Read;
import stream.SamLine;
import stream.SamLineStreamer;
import stream.SamReadInputStream;
import structures.CoverageArray;
import structures.ListNum;
import structures.LongList;
import structures.StandardDeviator;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date April 17, 2024
 *
 */
public class CoveragePileupMT implements Accumulator<CoveragePileupMT.LoadThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		Timer t=new Timer();
		
		CoveragePileupMT x=new CoveragePileupMT(args);
		
		x.process();
		
		t.stop();
		if(!KEY_VALUE){
			x.outstream.println();
			x.outstream.println("Time: \t"+t);
		}
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public CoveragePileupMT(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, printCommand ? getClass() : null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		int vectorMode=-1;
		ReadWrite.USE_UNPIGZ=true;
//		SamLine.RNAME_AS_BYTES=false;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("ref") || a.equals("reference") || a.equals("fasta")){
				reference=b;
			}else if(a.equals("streams") || a.equals("maxstreams")){
				streams=Integer.parseInt(b);
			}else if(a.equals("atomic")){
				atomic=Parse.parseBoolean(b);
			}else if(a.equals("prealloc")){
				prealloc=Parse.parseBoolean(b);
			}else if(a.equals("addfromref")){
				ADD_FROM_REF=Parse.parseBoolean(b);
			}else if(a.equals("in") || a.equals("in1")){
				inputFiles.clear();
				if(b!=null) {
					for(String f : Tools.commaPattern.split(b)) {inputFiles.add(f);}
				}
			}else if(a.equals("out") || a.equals("coveragestats") || a.equals("covstats") || a.equals("stats")){
				covstats=b;
			}else if(a.equals("minscaf") || a.equals("covminscaf")){
				minscaf=Integer.parseInt(b);
			}else if(a.equals("mindepth") || a.equals("mincov")){
				minDepthToBeCovered=Integer.parseInt(b);
			}else if(a.equals("border")){
				border=Integer.parseInt(b);
			}else if(a.equals("qtrim")/* || a.equals("trim")*/){
				if(b==null || b.length()==0){qtrimRight=qtrimLeft=true;}
				else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){qtrimLeft=true;qtrimRight=false;}
				else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){qtrimLeft=false;qtrimRight=true;}
				else if(b.equalsIgnoreCase("both") || b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr")){qtrimLeft=qtrimRight=true;}
				else if(Tools.isDigit(b.charAt(0))){
					trimq=Float.parseFloat(b);
					qtrimRight=trimq>0;
				}else{qtrimRight=qtrimLeft=Parse.parseBoolean(b);}
			}else if(a.equals("trimq") || a.equals("trimquality")){
				trimq=Float.parseFloat(b);
			}else if(a.equals("minq") || a.equals("minmapq")){
				minMapq=Integer.parseInt(b);
			}else if(a.equals("rpkm") || a.equals("fpkm") || a.equals("outrpkm")){
				outrpkm=b;
			}else if(a.equals("outorf")){
				outorf=b;
			}else if(a.equals("orffasta") || a.equals("fastaorf")){
				orffasta=b;
			}else if(a.equals("basecov") || a.equals("outcov")){
				basecov=b;
			}else if(a.equals("bincov") || a.equals("outbinned")){
				bincov=b;
			}else if(a.equals("normcov") || a.equals("outnormalized")){
				normcov=b;
			}else if(a.equals("normcovo") || a.equals("outnormalizedoverall")){
				normcovOverall=b;
			}else if(a.equals("delta")){
				DELTA_ONLY=Parse.parseBoolean(b);
			}else if(a.equals("physical") || a.equals("physicalcoverage") || a.equals("physcov")){
				PHYSICAL_COVERAGE=Parse.parseBoolean(b);
			}else if(a.equals("tlen")){
				USE_TLEN=Parse.parseBoolean(b);
			}else if(a.equals("hist") || a.equals("histogram") || a.equals("covhist")){
				histogram=b;
			}else if(a.equals("histmax")){
				HISTMAX=Parse.parseIntKMG(b);
			}else if(a.equals("reads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("scafs") || a.equals("scaffolds")){
				initialScaffolds=Tools.mid(128, Integer.parseInt(b), 2000000000);
			}else if(a.equals("binsize")){
				binsize=Integer.parseInt(b);
			}else if(a.equals("32bit")){
				bits32=Parse.parseBoolean(b);
			}else if(a.equals("bitset") || a.equals("usebitset") || a.equals("bitsets") || a.equals("usebitsets")){
//				if(Parse.parseBoolean(b)){arrayMode=BITSET_MODE;}
				vectorMode=Parse.parseBoolean(b) ? BITSET_MODE : NOTHING_MODE;
			}else if(a.equals("array") || a.equals("arrays") || a.equals("usearrays")){
				vectorMode=Parse.parseBoolean(b) ? ARRAY_MODE : NOTHING_MODE;
			}else if(a.equals("median") || a.equals("calcmedian")){
				if(Parse.parseBoolean(b)){
					vectorMode=ARRAY_MODE;
				}
			}else if(a.startsWith("nonzero") || a.equals("nzo")){
				NONZERO_ONLY=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set NONZERO_ONLY to "+NONZERO_ONLY);}
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set overwrite to "+overwrite);}
			}else if(a.equalsIgnoreCase("twocolumn")){
				TWOCOLUMN=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set TWOCOLUMN to "+TWOCOLUMN);}
			}else if(a.equalsIgnoreCase("keyvalue") || a.equalsIgnoreCase("machineout")){
				KEY_VALUE=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("countgc")){
				COUNT_GC=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set COUNT_GC to "+COUNT_GC);}
			}else if(a.equals("secondary") || a.equals("usesecondary")){
				USE_SECONDARY=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set USE_SECONDARY_ALIGNMENTS to "+USE_SECONDARY);}
			}else if(a.equals("softclip") || a.equals("includesoftclip")){
				INCLUDE_SOFT_CLIP=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set INCLUDE_SOFT_CLIP to "+INCLUDE_SOFT_CLIP);}
			}else if(a.equals("keepshortbins") || a.equals("ksb")){
				KEEP_SHORT_BINS=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set KEEP_SHORT_BINS to "+KEEP_SHORT_BINS);}
			}else if(a.equals("strandedcoverage") || a.equals("strandedcov") || a.equals("covstranded") || a.equals("stranded")){
				STRANDED=Parse.parseBoolean(b);
			}else if(a.equals("startcov") || a.equals("covstart") || a.equals("startonly")){
				START_ONLY=Parse.parseBoolean(b);
			}else if(a.equals("stopcov") || a.equals("covstop") || a.equals("stoponly")){
				STOP_ONLY=Parse.parseBoolean(b);
			}else if(a.equals("concise")){
				CONCISE=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("verbosetime")){
				verboseTime=StandardDeviator.verbose=Parse.parseBoolean(b);
			}else if(a.equals("normc") || a.equals("normalizecoverage")){
				NORMALIZE_COVERAGE=Parse.parseBoolean(b);
			}else if(a.equals("header") || a.equals("hdr")){
				printHeader=Parse.parseBoolean(b);
			}else if(a.equals("headerpound") || a.equals("#")){
				headerPound=Parse.parseBoolean(b);
			}else if(a.equals("stdev")){
				calcCovStdev=Parse.parseBoolean(b);
			}else if(a.equals("delcov") || a.equals("dels") || a.equals("includedels") || a.equals("includedeletions") || a.equals("delcoverage")){
				INCLUDE_DELETIONS=Parse.parseBoolean(b);
			}else if(a.equals("dupecoverage") || a.equals("dupecov") || a.equals("dupes") || a.equals("duplicates") || a.equals("includeduplicates")){
				INCLUDE_DUPLICATES=Parse.parseBoolean(b);
			}else if(a.equals("ignoredupes") || a.equals("ignoreduplicates")){
				INCLUDE_DUPLICATES=!Parse.parseBoolean(b);
			}else if(a.equals("normb") || a.equals("normalizebins")){
				try {
					NORMALIZE_LENGTH_BINS=Integer.parseInt(b);
				} catch (NumberFormatException e) {
					boolean x=Parse.parseBoolean(b);
					NORMALIZE_LENGTH_BINS=x ? 100 : -1;
				}
			}else if(a.equals("covwindow")){
				if(b==null || b.length()<1 || Character.isLetter(b.charAt(0))){
					USE_WINDOW=Parse.parseBoolean(b);
				}else{
					LOW_COV_WINDOW=Integer.parseInt(b);
					USE_WINDOW=(LOW_COV_WINDOW>0);
				}
			}else if(a.equals("covwindowavg")){
				LOW_COV_DEPTH=Double.parseDouble(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(arg.indexOf('=')<0 && new File(arg).exists()){
				inputFiles.add(arg);
			}else{
				throw new RuntimeException("Unknown parameter: "+args[i]);
			}
			
		}
		trimE=(float)QualityTools.phredToProbError(trimq);
//		assert(false) : qtrimLeft+", "+qtrimRight+", "+trimq+", "+trimE;
		
		ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_UNMAPPED;
		if(!USE_SECONDARY){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_SECONDARY;}
		if(!INCLUDE_DUPLICATES){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_DUPLICATE;}
		
		SamLine.PARSE_0=false;
		SamLine.PARSE_6=false;
		SamLine.PARSE_7=false;
		SamLine.PARSE_8=false;
		if(k<1 && trimq<1){SamLine.PARSE_10=false;}
		SamLine.PARSE_OPTIONAL=false;
		
		prealloc=prealloc || atomic;
		caType=CoverageArray.getType(atomic, bits32);
		
		if(vectorMode>-1){
			USE_BITSETS=(vectorMode==BITSET_MODE);
			USE_COVERAGE_ARRAYS=(vectorMode==ARRAY_MODE);
		}else{
			if(histogram==null && basecov==null && bincov==null && normcov==null &&
					normcovOverall==null && outorf==null && !calcCovStdev){//No need for coverage array!
				USE_COVERAGE_ARRAYS=false;
				if(TWOCOLUMN){//No need for bitset, either!
					USE_BITSETS=false;
				}else{
					USE_BITSETS=true;
				}
			}
		}
		
		if(verbose){
			outstream.println("Set USE_COVERAGE_ARRAYS to "+USE_COVERAGE_ARRAYS);
			outstream.println("Set USE_BITSETS to "+USE_BITSETS);
		}
		
		if(maxReads<0){maxReads=Long.MAX_VALUE;}
		assert(inputFiles!=null && inputFiles.size()>0);
		
		if(STRANDED){
			assert(basecov==null || basecov.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
			assert(bincov==null || bincov.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
			assert(normcov==null || normcov.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
			assert(normcovOverall==null || normcovOverall.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
			assert(histogram==null || histogram.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
			assert(covstats==null || covstats.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
		}
		
		if(!Tools.testInputFiles(true, true, inputFiles)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(inputFiles.isEmpty()) {throw new RuntimeException("\nAt least one input file is required.\n");}
		
		if(!Tools.testInputFiles(false, true, reference)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, basecov, bincov, normcov, normcovOverall, histogram, covstats, outrpkm)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					basecov+", "+bincov+", "+normcov+", "+normcovOverall+", "+histogram+", "+covstats+", "+outrpkm+"\n");
		}
	}

	
	/*--------------------------------------------------------------*/
	/*----------------       Data Structures        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/** The goal of this is to garbage-collect unnecessary objects, not really for reusing the object */
	public void clear(){
		list=null;
		table=null;
		pairTable=null;
		
		program=null;
		version=null;
		
		inputFiles=null;
		covstats=null;
		outorf=null;
		outrpkm=null;
		reference=null;
		histogram=null;
		basecov=null;
		bincov=null;
		normcov=null;
		normcovOverall=null;
		orffasta=null;
		
		error=false;

		refBases=0;
		mappedBases=0;
		mappedNonClippedBases=0;
		mappedBasesWithDels=0;
		mappedReads=0;
		properPairs=0;
		readsProcessed=0;
		basesProcessed=0;
		kmersProcessed=0;
		mappedKmers=0;
		totalCoveredBases1=0;
		totalCoveredBases2=0;
		scaffoldsWithCoverage1=0;
		scaffoldsWithCoverage2=0;
		totalScaffolds=0;
	}
	
	public void createDataStructures(){
		refBases=0;
		mappedBases=0;
		mappedNonClippedBases=0;
		mappedBasesWithDels=0;
		mappedReads=0;
		properPairs=0;
		readsProcessed=0;
		basesProcessed=0;
		kmersProcessed=0;
		mappedKmers=0;
		totalCoveredBases1=0;
		totalCoveredBases2=0;
		scaffoldsWithCoverage1=0;
		scaffoldsWithCoverage2=0;
		totalScaffolds=0;
		error=false;
		list=new ArrayList<Scaffold>(initialScaffolds);
		table=new HashMap<String, Scaffold>(initialScaffolds);
		
		if(PHYSICAL_COVERAGE){
			pairTable=new HashMap<String, SamLine>();
			if(COUNT_GC){
				COUNT_GC=false;
				outstream.println("COUNT_GC disabled for physical coverage mode.");
			}
			if(USE_SECONDARY){
				USE_SECONDARY=false;
				outstream.println("USE_SECONDARY disabled for physical coverage mode.");
			}
			
			SamLine.PARSE_0=true;
			SamLine.PARSE_6=true;
			SamLine.PARSE_7=true;
			SamLine.PARSE_8=true;
			SamLine.PARSE_10=false;
			SamLine.PARSE_OPTIONAL=false;
		}
	}
	
	final CoverageArray makeCA(int len) {
		CoverageArray ca=CoverageArray.makeArray(1, len, caType);
		assert(ca.length()==len);
		return ca;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Read and process all input data. */
	public void process(){
		Timer t=new Timer();
		createDataStructures();
		
		if(inputFiles!=null && FileFormat.isBamFile(inputFiles.get(0))){
			if(Data.SAMTOOLS()){ReadWrite.USE_SAMBAMBA=false;} //Disable because it takes forever to read the header
		}
		
		ByteFile tf=ByteFile.makeByteFile(inputFiles.get(0), false);

		processHeader(tf);
		errorState=tf.close()|errorState;
		if(verboseTime) {t.stopAndStart("Process Header:");}
		
		ReadWrite.USE_SAMBAMBA=true;
		
		if(reference!=null) {
			processReference();
			if(verboseTime) {t.stopAndStart("Process Reference:");}
		}
		if(maxReads<0){maxReads=Long.MAX_VALUE;}
		
		spawnThreads();
		if(verboseTime) {t.stopAndStart("Process Reads:");}

		StandardDeviator sd=new StandardDeviator(STRANDED, 0);
		if(USE_COVERAGE_ARRAYS && list.size()>1) {
			boolean calcStd=true, calcMedian=true, calcHist=(histogram!=null);
			final int histmax=HISTMAX>0 ? HISTMAX : (bits32 ? 1000000 : Character.MAX_VALUE);
			sd.calculateStuff(Shared.threads(), list, calcStd, calcMedian, calcHist, 
					minDepthToBeCovered, histmax, USE_WINDOW ? LOW_COV_DEPTH : -1, LOW_COV_WINDOW);
			if(verboseTime) {t.stopAndStart("MT Calculation:");}
		}
		
		printOutput(sd.hist0, sd.hist1);
//		if(verboseTime) {t.stopAndStart("Print Output:");}
		
		if(orffasta!=null){
			processOrfsFasta(orffasta, outorf, table);
			if(verboseTime) {t.stopAndStart("Process ORFs:");}
		}
	}
	
	/** Spawn process threads */
	private void spawnThreads(){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		if(streams<1) {streams=(Shared.threads()+1)/2;}
		final int threads=Tools.mid(1, streams, inputFiles.size());
		
		//Fill a list with ProcessThreads
		ArrayList<LoadThread> alpt=new ArrayList<LoadThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new LoadThread(i));
		}
		
		//Start the threads and wait for them to finish
		ThreadWaiter.startThreads(alpt);
		
		for(String fname : inputFiles) {
			addToQueue(fname);
		}
		addToQueue(POISON);
		
		boolean success=ThreadWaiter.waitForThreadsToFinish(alpt, this);
		outstream.println("Finished processing "+readsProcessed+" reads and "+basesProcessed+" bases in "+Tools.plural("file", inputFiles.size())+".");
		errorState&=!success;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Setup            ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/** Process all sam header lines from the tf.
	 * Once a non-header line is encountered, return it.
	 * If non-null, print all lines to the tsw. */
	public void processHeader(ByteFile tf){
		byte[] line=null;
		for(line=tf.nextLine(); line!=null && (line.length==0 || line[0]=='@'); line=tf.nextLine()){
			processHeaderLine(line);
		}
		if(line!=null){tf.pushBack(line);}
//		return line;
	}
	
	void processHeaderLine(byte[] line) {
		if(line==null || line.length<3) {return;}
		final byte a=line[1], b=line[2];

		if(a=='S' && b=='Q'){
			lp.set(line);
			Scaffold scaf=new Scaffold(lp);
			if(prealloc){
				if(USE_COVERAGE_ARRAYS) {
					scaf.obj0=makeCA(scaf.length);
					if(STRANDED){scaf.obj1=makeCA(scaf.length);}
				}else if(USE_BITSETS) {
					scaf.obj0=new BitSet(scaf.length);
					if(STRANDED){scaf.obj1=new BitSet(scaf.length);}
				}
			}
			if(COUNT_GC){scaf.basecount=KillSwitch.allocLong1D(8);}
			assert(!table.containsKey(scaf.name)) : "\nDuplicate scaffold name!\n"+scaf+"\n\n"+table.get(scaf.name);
			table.put(scaf.name, scaf);
			list.add(scaf);
			refBases+=scaf.length;
			refKmers+=scaf.length-k+1;

			//				sc.obj=new CoverageArray2(table.size(), sc.length);
			//				outstream.println("Made scaffold "+sc.name+" of length "+sc.length);
		}else if(a=='P' && b=='G'){
			lp.set(line);
			for(int i=1, terms=lp.terms(); i<terms; i++) {
				if(lp.termStartsWith("PN:", i)){
					if(program==null){
						lp.incrementA(3);
						program=lp.parseStringFromCurrentField();
					}
				}else if(lp.termStartsWith("VN:", i)){
					if(version==null){
						lp.incrementA(3);
						version=lp.parseStringFromCurrentField();
					}
				}
			}
		}else if(a=='R' && b=='G'){
			//Do nothing
		}else if(a=='H' && b=='D'){
			//Do nothing
		}else if(a=='C' && b=='O'){
			//Do nothing
		}else{
			//				assert(false) : line;
		}
	}
	
	public void processReference(){
		if(reference==null){return;}

		ByteFile bf=ByteFile.makeByteFile(reference, true);
		Scaffold scaf=null;
		int len=0;
		final long[] acgtn=KillSwitch.allocLong1D(8);
		boolean addLen=false;
		for(byte[] s=bf.nextLine(); s!=null; s=bf.nextLine()){
			if(s.length>0 && s[0]=='>'){
				if(scaf!=null){//Then evict the old scaf and start a new one
					if(scaf.length>0 && scaf.length!=len){
						outstream.println("ERROR: Scaffold "+scaf.name+" has contradictory lengths of "+scaf.length+" and "+len+"\n"
								+ "This probably indicates a corrupt or incorrect reference.");
						errorState=true;
						if(Shared.EA()){KillSwitch.kill();}
					}
//					else{outstream.println(scaf.name+", "+scaf.length+", "+len);}
					scaf.length=len;
					if(addLen){
						refBases+=scaf.length;
						refKmers+=scaf.length-k+1;
						
						addLen=false;
					}
					scaf.gc=(float)((acgtn[1]+acgtn[2])*1d/Data.max(1, acgtn[0]+acgtn[1]+acgtn[2]+acgtn[3]));
					scaf=null;
					len=0;
					Arrays.fill(acgtn, 0);
				}
				
				String name=new String(s, 1, s.length-1);
				String shortName=Tools.trimToWhitespace(name);
				
				scaf=table.get(shortName);
				if(ADD_FROM_REF && scaf==null){
					scaf=new Scaffold(name, 0);
					if(!warned){
						outstream.println("Warning - SAM header did not include "+name+"\nAbsent scaffolds will be added; further warnings will be suppressed.");
						warned=true;
					}
					if(COUNT_GC){scaf.basecount=KillSwitch.allocLong1D(8);}
					table.put(shortName, scaf);
					list.add(scaf);
					addLen=true;
				}
			}else{
				len+=s.length;
				for(int i=0; i<s.length; i++){
					acgtn[charToNum[s[i]]]++;
				}
			}
		}
		if(scaf!=null){
			if(scaf.length>0 && scaf.length!=len){
				outstream.println("ERROR: Scaffold "+scaf.name+" has contradictory lengths of "+scaf.length+" and "+len+"\n"
						+ "This probably indicates a corrupt or incorrect reference.");
				errorState=true;
				if(Shared.EA()){KillSwitch.kill();}
			}
//			else{outstream.println(scaf.name+", "+scaf.length+", "+len);}
			scaf.length=len;
			if(addLen){
				refBases+=scaf.length;
				refKmers+=scaf.length-k+1;
				addLen=false;
			}
			scaf.gc=(float)((acgtn[1]+acgtn[2])*1d/Data.max(1, acgtn[0]+acgtn[1]+acgtn[2]+acgtn[3]));
			scaf=null;
			len=0;
			Arrays.fill(acgtn, 0);
		}
	}
	
	
	public void processOrfsFasta(String fname_in, String fname_out, HashMap<String, Scaffold> map){
		TextFile tf=new TextFile(fname_in, false);
		assert(!fname_in.equalsIgnoreCase(fname_out));
		TextStreamWriter tsw=new TextStreamWriter(fname_out, overwrite, false, true);
		tsw.start();
		
		if(printHeader){
			String pound=(headerPound ? "#" : "");
			tsw.print(pound+"mappedBases="+mappedBases+"\n");
			tsw.print(pound+"mappedNonClippedBases="+mappedNonClippedBases+"\n");
			tsw.print(pound+"mappedBasesWithDels="+mappedBasesWithDels+"\n");
			tsw.print(pound+"mappedReads="+mappedReads+"\n");
			tsw.print(pound+"name\tlength\tdepthSum\tavgDepth\tavgDepth/mappedBases\tminDepth\tmaxDepth\tmedianDepth\tstdDevDepth\tfractionCovered\n");
		}
		
		String line;
		final StringBuilder sb=new StringBuilder(500);
//		Formatter formatter=new Formatter(sb);
		
		while((line=tf.nextLine())!=null){
			if(line.length()>1 && line.charAt(0)=='>'){
				
				String[] split=line.split(" # "); //' # ' used as delimiters
				
				String orfname=split[0].substring(1).trim(); //In case there are spaces around the ' # ' delimiters
				String scafname=orfname;
				if(scafname.contains("_")){//PRODIGAL pads _1 to the name of the first orf of a scaffold, and etc
					int last=scafname.lastIndexOf('_');
					boolean numeric=false;
					for(int i=last+1; i<scafname.length(); i++){
						if(Tools.isDigit(scafname.charAt(i))){numeric=true;}
						else{numeric=false; break;}
					}
					if(numeric){scafname=scafname.substring(0, last);}
				}
				
				int start=Integer.parseInt(split[1].trim());
				int stop=Integer.parseInt(split[2].trim());
				int strand=Integer.parseInt(split[3].trim());
				if(strand==1){strand=Shared.PLUS;}else{strand=Shared.MINUS;}
				Orf orf=new Orf(orfname, start, stop, (byte)strand);
				
				Scaffold scaf=map.get(scafname);
				if(scaf==null){
					outstream.println("Can't find scaffold for ("+orf+")\nfrom line\n"+line+"\nscafname='"+scafname+"'\norfname='"+orfname+"'");
					if(ABORT_ON_ERROR){
						tsw.poison();
						throw new RuntimeException("Aborting.");
					}
				}else{
					if(orf.start<0 && orf.stop>=scaf.length){
						outstream.println("orf goes out of scaffold bounds.\n"+orf+"\n"+scaf);
						if(ABORT_ON_ERROR){
							tsw.poison();
							throw new RuntimeException("Aborting.");
						}
					}
					
					CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0); //TODO:  Strand logic here depends on stranding protocol.
					orf.readCoverageArray(ca);
				}
				
//				{
//					tsw.print(Tools.format("%s\t%d\t", args));
//				}
				
				sb.append(orf.name).append('\t');
				sb.append(orf.length()).append('\t');
				sb.append(orf.baseDepth).append('\t');
				sb.append(Tools.format("%.4f", orf.avgCoverage())).append('\t');
				sb.append(orf.avgCoverage()/mappedNonClippedBases);

				sb.append('\t');
				sb.append(orf.minDepth).append('\t');
				sb.append(orf.maxDepth).append('\t');
				sb.append(orf.medianDepth).append('\t');
				sb.append(Tools.format("%.4f",orf.stdevDepth)).append('\t');
				sb.append(Tools.format("%.4f",orf.fractionCovered()));

				sb.append('\n');
				tsw.print(sb.toString());
				sb.setLength(0);
			}
		}
		
		tsw.poisonAndWait();
	}
	
	private int trim(SamLine sl){
		assert(border>0 || (trimq>=0 && (qtrimLeft || qtrimRight)));
		Read r=null;
		Scaffold sc=null;
		int leftTrimAmount=border, rightTrimAmount=border;
		if(border>0){
			r=sl.toRead(false);
			sc=table.get(sl.rnameS());
			assert(sc!=null) : sl+"\n\n"+sl.rnameS()+"\n\n"+table;
			int skipTrimRange=Tools.max(10, border+5);
			if(r.start<skipTrimRange){
				if(r.strand()==Shared.PLUS){leftTrimAmount=0;}
				else{rightTrimAmount=0;}
			}
			if(r.stop>sc.length-skipTrimRange){
				if(r.strand()==Shared.PLUS){rightTrimAmount=0;}
				else{leftTrimAmount=0;}
			}
		}
		final int len0=sl.length();
		if(qtrimLeft || qtrimRight){
			long packed=TrimRead.testOptimal(sl.seq, sl.qual, trimE);
			if(qtrimLeft){leftTrimAmount=Tools.max(leftTrimAmount, (int)((packed>>32)&0xFFFFFFFFL));}
			if(qtrimRight){rightTrimAmount=Tools.max(rightTrimAmount, (int)((packed)&0xFFFFFFFFL));}
//			assert(false) : qtrimLeft+", "+qtrimRight+", "+trimq+", "+trimE+", "+rightTrimAmount+"\n"+sl+"\n";
		}
		final int trimmed;
		if(leftTrimAmount<1 && rightTrimAmount<1){trimmed=0;}
		else{
			if(r==null){r=sl.toRead(false);}
			if(sc==null){sc=table.get(sl.rnameS());}
			int scaflen=(sc==null ? 1999999999 : sc.length);
			trimmed=TrimRead.trimReadWithMatch(r, sl, leftTrimAmount, rightTrimAmount, 0, scaflen, false);
		}
//		assert(trimmed==len0-sl.length()) : trimmed+", "+len0+", "+sl.length();
//		assert(rightTrimAmount>0) : qtrimLeft+", "+qtrimRight+", "+trimq+", "+trimE+", "+rightTrimAmount+"\n"+sl+"\n";
		return trimmed;
	}

	
	/*--------------------------------------------------------------*/
	/*----------------        Output Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public void printOutput(LongList hist0, LongList hist1){
		Timer t=new Timer();
		
		totalScaffolds=list.size();
		
		String basecov1=(basecov==null ? null : (STRANDED ? basecov.replaceFirst("#", "1") : basecov));
		String bincov1=(bincov==null ? null : (STRANDED ? bincov.replaceFirst("#", "1") : bincov));
		String normcov1=(normcov==null ? null : (STRANDED ? normcov.replaceFirst("#", "1") : normcov));
		String normcovOverall1=(normcovOverall==null ? null : (STRANDED ? normcovOverall.replaceFirst("#", "1") : normcovOverall));
		String histogram1=(histogram==null ? null : (STRANDED ? histogram.replaceFirst("#", "1") : histogram));
		String stats1=(covstats==null ? null : (STRANDED ? covstats.replaceFirst("#", "1") : covstats));

		String basecov2=(basecov==null || !STRANDED ? null : basecov.replaceFirst("#", "2"));
		String bincov2=(bincov==null || !STRANDED ? null : bincov.replaceFirst("#", "2"));
		String normcov2=(normcov==null || !STRANDED ? null : normcov.replaceFirst("#", "2"));
		String normcovOverall2=(normcovOverall==null ? null : (STRANDED ? normcovOverall.replaceFirst("#", "2") : normcovOverall));
		String histogram2=(histogram==null || !STRANDED ? null : histogram.replaceFirst("#", "2"));
		String stats2=(covstats==null || !STRANDED ? null : covstats.replaceFirst("#", "2"));
		
		if(basecov!=null) {
			if(CONCISE){
				writeCoveragePerBaseConcise(basecov1, list, 0, minscaf);
				writeCoveragePerBaseConcise(basecov2, list, 1, minscaf);
			}else{
				writeCoveragePerBase(basecov1, list, DELTA_ONLY, 0, minscaf);
				writeCoveragePerBase(basecov2, list, DELTA_ONLY, 1, minscaf);
			}
			if(verboseTime) {t.stopAndStart("Write Coverage:");}
		}
		if(bincov!=null) {
			if(KEEP_SHORT_BINS){
				writeCoveragePerBaseBinned2(bincov1, list, binsize, 0, minscaf);
				writeCoveragePerBaseBinned2(bincov2, list, binsize, 1, minscaf);
			}else{
				writeCoveragePerBaseBinned(bincov1, list, binsize, 0, minscaf);
				writeCoveragePerBaseBinned(bincov2, list, binsize, 1, minscaf);
			}
			if(verboseTime) {t.stopAndStart("Write BinCov:");}
		}
		if(normcov!=null){
			writeCoveragePerBaseNormalized(normcov1, list, binsize, 0, minscaf);
			writeCoveragePerBaseNormalized(normcov2, list, binsize, 1, minscaf);
			if(verboseTime) {t.stopAndStart("Write NormCov:");}
		}
		if(normcovOverall!=null){
			writeCoveragePerBaseNormalizedOverall(normcovOverall1, list, binsize, 0, minscaf);
			writeCoveragePerBaseNormalizedOverall(normcovOverall2, list, binsize, 1, minscaf);
			if(verboseTime) {t.stopAndStart("Write NormCov2:");}
		}
		if(outrpkm!=null){
			writeRPKM(outrpkm, inputFiles, readsProcessed, NONZERO_ONLY,list);
			if(verboseTime) {t.stopAndStart("Write RPKM:");}
		}
		
		{
			writeStats(stats1, 0);
			if(STRANDED){writeStats(stats2, 1);}
			if(verboseTime) {t.stopAndStart("Write Stats:");}
		}
		
		if(histogram!=null) {
			if(hist0!=null){writeHist(histogram1, hist0.array);}
			if(STRANDED && hist1!=null){writeHist(histogram2, hist1.array);}
			if(verboseTime) {t.stopAndStart("Write Histogram:");}
		}
		
		final double mult=1.0/refBases;
		double depthCovered=mappedBases*mult;
		double depthCovered2=mappedNonClippedBases*mult;
		double depthCovered3=mappedBasesWithDels*mult;
		double pctScaffoldsWithCoverage=scaffoldsWithCoverage1*100.0/totalScaffolds;
		double pctCovered=totalCoveredBases1*100*mult;

		StandardDeviator sdTool=new StandardDeviator(STRANDED, 0);
		
		if(!KEY_VALUE){
			outstream.println("Reads:                               \t"+readsProcessed);
			outstream.println("Mapped reads:                        \t"+mappedReads);
			//		outstream.println("Mapped bases:                        \t"+mappedBases);
			//		outstream.println("Mapped non-clipped bases:            \t"+mappedNonClippedBases);
			outstream.println("Mapped bases:                        \t"+mappedNonClippedBases);
			outstream.println("Ref scaffolds:                       \t"+totalScaffolds);
			outstream.println("Ref bases:                           \t"+refBases);
			
			if(k>0){
				String kcovS=k+"-mer coverage:";
				String kcorrectS="Percent correct "+k+"-mers:";
				while(kcovS.length()<26){kcovS=kcovS+" ";}
				while(kcovS.length()<26){kcorrectS=kcorrectS+" ";}
				outstream.println(Tools.format("\n"+kcovS+"           \t%.3f", mappedKmers*1.0/refKmers));
				outstream.println(Tools.format(kcorrectS+"           \t%.3f", 100*correctKmers/kmersProcessed));
				//			outstream.println(kmersProcessed+", "+correctKmers);
			}
			
			outstream.println(Tools.format("\nPercent mapped:                      \t%.3f", mappedReads*100f/readsProcessed));
			outstream.println(Tools.format("Percent proper pairs:                \t%.3f", properPairs*100f/readsProcessed));
//			outstream.println(depthCovered);
			outstream.println(Tools.format("Average coverage:                    \t%.3f", depthCovered2));
			outstream.println(Tools.format("Average coverage with deletions:     \t%.3f", depthCovered3));
			if(USE_COVERAGE_ARRAYS && calcCovStdev){
				double[] stdev=sdTool.standardDeviation(list, minscaf);
				outstream.println(Tools.format("Standard deviation:                    \t%.3f", stdev[1]));
			}
			outstream.println(Tools.format("Percent scaffolds with any coverage: \t%.2f", pctScaffoldsWithCoverage));
			if(USE_COVERAGE_ARRAYS || USE_BITSETS){
				outstream.println(Tools.format("Percent of reference bases covered:  \t%.2f", pctCovered));
			}
		}else{
			outstream.println("reads="+readsProcessed);
			outstream.println("mappedReads="+mappedReads);
			outstream.println("mappedBases="+mappedNonClippedBases);
			outstream.println("mappedBasesWithDels="+mappedBasesWithDels);
			outstream.println("refScaffolds="+totalScaffolds);
			outstream.println("refBases="+refBases);
			outstream.println(Tools.format("percentMapped=%.3f", mappedReads*100f/readsProcessed));
			outstream.println(Tools.format("percentPaired=%.3f", properPairs*100f/readsProcessed));
			outstream.println(Tools.format("averageCoverage=%.3f", depthCovered2));
			if(USE_COVERAGE_ARRAYS && calcCovStdev){
				double[] stdev=sdTool.standardDeviation(list, minscaf);
				outstream.println(Tools.format("standardDeviation=%.3f", stdev[1]));
			}
			outstream.println(Tools.format("percentCoveredScaffolds=%.2f", pctScaffoldsWithCoverage));
			if(USE_COVERAGE_ARRAYS || USE_BITSETS){
				outstream.println(Tools.format("percentCoveredBases=%.2f", pctCovered));
			}
		}
		if(verboseTime) {t.stopAndStart("Calc StdDev:");}
	}
	
	public void writeStats(String fname, int strand){
//		outstream.println("Writing stats for "+fname+", "+strand);
		final ByteStreamWriter tsw=(fname==null ? null : new ByteStreamWriter(fname, overwrite, false, true));
		
		if(tsw!=null){
			tsw.start();
			if(printHeader){
				String pound=(headerPound ? "#" : "");
				if(TWOCOLUMN){
					tsw.println(pound+"ID\tAvg_fold");
				}else{
					tsw.println(pound+"ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads"+
							(COUNT_GC ? "\tRead_GC" : "")+
							(USE_COVERAGE_ARRAYS ? ("\tMedian_fold\tStd_Dev") : "")+
							(USE_WINDOW ? "\tUnder_"+Tools.format("%.0f",LOW_COV_DEPTH)+"/"+LOW_COV_WINDOW : ""));
				}
			}
			//Maximally:
			//"ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tRead_GC\tMedian_fold\tStd_Dev\tUnder_X"
		}
		
		long coveredScafTemp=0;
		long coveredBaseTemp=0;
		for(Scaffold scaf : list){
			synchronized(scaf) {
				final long sum=scaf.basehits;
				int covered=0;
				int median=0;
				int underWindowAverage=0;
				final double stdev;
				if(USE_COVERAGE_ARRAYS){
					CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
					if(ca!=null){
						covered=ca.covered(minDepthToBeCovered);
						stdev=ca.standardDeviation();
						underWindowAverage=(USE_WINDOW ? 
								ca.basesUnderAverageCoverage(LOW_COV_DEPTH, LOW_COV_WINDOW) : -1);
						median=ca.median();//This one is the most expensive, especially with atomics.
					}else{
						stdev=0;
					}
				}else if(USE_BITSETS){
					BitSet bs=(BitSet)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
					covered=(bs==null ? 0 : bs.cardinality());
					stdev=-1;
				}else{
					stdev=-1;
				}

				if(sum>0){
					coveredScafTemp++;
				}
				//			pw.print(scaf.name);
				if(tsw!=null && (sum>0 || !NONZERO_ONLY) && scaf.length>=minscaf){
					if(TWOCOLUMN){
//						tsw.print(Tools.format("%s\t%.4f\n", scaf.name, sum/(double)scaf.length));
						tsw.print(scaf.name).tab().print(sum/(double)scaf.length, 4).nl();
					}else {
						//Variable portion:
						//"\tRead_GC\tMedian_fold\tStd_Dev\tUnder_X\n"
						//\t%.4f\t%d\t%.2f\t%d\n
						
						tsw.print(scaf.name).tab().print(sum/(double)scaf.length, 4).tab().print(scaf.length).tab();
						tsw.print(scaf.gc, 4).tab().print(covered*100d/scaf.length, 3).tab().print(covered).tab();
						tsw.print(scaf.readhits-scaf.readhitsMinus).tab().print(scaf.readhitsMinus);
						if(COUNT_GC) {
							long[] bc=scaf.basecount;
							double gc=(bc[1]+bc[2])*1d/Data.max(1, bc[0]+bc[1]+bc[2]+bc[3]);
							tsw.tab().print(gc, 4);
						}
						if(USE_COVERAGE_ARRAYS) {tsw.tab().print(median).tab().print(stdev, 2);}
						if(USE_WINDOW) {tsw.tab().print(underWindowAverage);}
						tsw.nl();
					}
				}
				coveredBaseTemp+=covered;
			}
		}
		
		if(strand==0){
			scaffoldsWithCoverage1+=coveredScafTemp;
			totalCoveredBases1+=coveredBaseTemp;
		}else{
			scaffoldsWithCoverage2+=coveredScafTemp;
			totalCoveredBases2+=coveredBaseTemp;
		}
		
		if(tsw!=null){tsw.poisonAndWait();}
	}
	
	/**
	 * Write a histogram of number of bases covered to each depth
	 * @param fname Output filename
	 * @param counts counts[X] stores the number of bases with coverage X
	 */
	public static void writeHist(String fname, long[] counts){
		if(fname==null){return;}
		assert(counts!=null) : "Can't write a histogram with null counts.";
		ByteStreamWriter tsw=new ByteStreamWriter(fname, overwrite, false, false);
		tsw.start();
		if(printHeader){
			if(headerPound){tsw.print('#');}
			tsw.println("Coverage\tnumBases");
		}
		int max=0;
		for(max=counts.length-1; max>0 && counts[max]==0; max--){}
		for(int i=0; i<=max; i++){
			long x=counts[i];
			tsw.print(i);
			tsw.print('\t');
			tsw.println(x);
		}
		
		tsw.poisonAndWait();
	}
	
	/**
	 * Prints coverage in this format:
	 * scafname TAB position TAB coverage
	 * scafname TAB position TAB coverage
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param deltaOnly Only write lines when coverage changes
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public void writeCoveragePerBase(String fname, ArrayList<Scaffold> list, boolean deltaOnly, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		
		if(verbose){outstream.println("Starting tsw "+fname);}
		ByteStreamWriter tsw=new ByteStreamWriter(fname, overwrite, false, true);
		if(verbose){outstream.println("Created tsw "+fname);}
		tsw.start();
//		if(verbose){outstream.println("Started tsw "+fname);}
		if(printHeader){
			if(headerPound){tsw.print('#');}
			tsw.println("RefName\tPos\tCoverage");
		}
		
		for(Scaffold scaf : list){
			int last=-1;
			CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			if(scaf.length>=minscaf){
				for(int i=0, len=scaf.length; i<len; i++){
					int x=(ca==null ? 0 : ca.get(i));
					if(!deltaOnly || x!=last){
						if(x>0 || !NONZERO_ONLY){
							tsw.print(scaf.name).tab().print(i).tab().print(x).nl();
						}
						last=x;
					}
				}
			}
		}

		if(verbose){outstream.println("Closing tsw "+fname);}
		tsw.poisonAndWait();
		if(verbose){outstream.println("Closed tsw "+fname);}
	}
	
	//TODO:  Add a super-concise version where pos and depth are blank if
	//they can be logically inferred from the previous line.
	
	/**
	 * Prints coverage in this format, skipping zero-coverage positions:
	 * #scafname
	 * position TAB coverage
	 * position TAB coverage
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public void writeCoveragePerBaseConcise(String fname, ArrayList<Scaffold> list, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		
		if(verbose){outstream.println("Starting tsw "+fname);}
		ByteStreamWriter tsw=new ByteStreamWriter(fname, overwrite, false, true);
		tsw.start();
		if(verbose){outstream.println("Started tsw "+fname);}
//		tsw.print(pound+"RefName\tPos\tCoverage\n");
		
		for(Scaffold scaf : list){
			tsw.print('#').println(scaf.name);
			CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			if(scaf.length>=minscaf){
				for(int i=0; i<scaf.length; i++){
					int x=(ca==null ? 0 : ca.get(i));
					if(x>0){
						tsw.print(i).tab().print(x).nl();
					}
				}
			}
		}

		if(verbose){outstream.println("Closing tsw "+fname);}
		tsw.poisonAndWait();
		if(verbose){outstream.println("Closed tsw "+fname);}
//		assert(false);
	}
	
	/**
	 * Note.  As written, this will truncate all trailing bases of each scaffold's length modulo binsize.
	 * For example, with binsize 1000, the last 500 bases of a 1500 base scaffold will be ignored.
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param binsize Width of coverage bins in bp
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public static void writeCoveragePerBaseBinned(String fname, ArrayList<Scaffold> list, int binsize, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		ByteStreamWriter tsw=new ByteStreamWriter(fname, overwrite, false, false);
		tsw.start();
		if(printHeader){
			String pound=(headerPound ? "#" : "");
			if(calcCovStdev){
				StandardDeviator sdTool=new StandardDeviator(STRANDED, strand);
				double[] stdev=sdTool.standardDeviation(list, minscaf);
//				double[] stdev=sdTool.standardDeviationBinned(list, binsize, minscaf);//Singlethreaded
				if(stdev!=null){
					tsw.print(pound+"Mean\t"+Tools.format("%.3f", stdev[0])+"\n");
					tsw.print(pound+"STDev\t"+Tools.format("%.3f", stdev[1])+"\n");
				}
			}
			tsw.print(pound+"RefName\tCov\tPos\tRunningPos\n");
		}
		
		long running=0;
		final float invbin=1f/binsize;
		for(Scaffold scaf : list){
			if(scaf.length>=binsize && scaf.length>=minscaf){
				CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
				int lastPos=-1, nextPos=binsize-1;
				long sum=0;
				for(int i=0; i<scaf.length; i++){
					int x=(ca==null ? 0 : ca.get(i));
					sum+=x;
					if(i>=nextPos){
						if(sum>0 || !NONZERO_ONLY){
							tsw.print(scaf.name).tab().print(sum*invbin, 2);
							tsw.tab().print(i+1).tab().print(running).nl();
						}
						lastPos=i;
						running+=binsize;
						nextPos+=binsize;
						sum=0;
					}
				}
			}
		}
		
		tsw.poisonAndWait();
	}
	
	/**
	 * This version will NOT truncate all trailing bases of each scaffold's length modulo binsize.
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param binsize Width of coverage bins in bp
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public static void writeCoveragePerBaseBinned2(String fname, ArrayList<Scaffold> list, int binsize, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		ByteStreamWriter tsw=new ByteStreamWriter(fname, overwrite, false, false);
		tsw.start();
		if(printHeader){
			String pound=(headerPound ? "#" : "");
			if(calcCovStdev){
				StandardDeviator sdTool=new StandardDeviator(STRANDED, strand);
				double[] stdev=sdTool.standardDeviation(list, minscaf);
//				double[] stdev=sdTool.standardDeviationBinned(list, binsize, minscaf);//Singlethreaded
				if(stdev!=null){
					tsw.print(pound+"Mean\t"+Tools.format("%.3f", stdev[0])+"\n");
					tsw.print(pound+"STDev\t"+Tools.format("%.3f", stdev[1])+"\n");
				}
			}
			tsw.print(pound+"RefName\tCov\tPos\tRunningPos\n");
		}
		
		long running=0;
		for(Scaffold scaf : list){
			CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			int lastPos=-1, nextPos=binsize-1;
			long sum=0;
			final int lim=scaf.length-1;
			for(int i=0; i<scaf.length; i++){
				int x=(ca==null ? 0 : ca.get(i));
				sum+=x;
				if(i>=nextPos || i==lim){
					int bin=(i-lastPos);
					if(scaf.length>=minscaf){
						if(sum>0 || !NONZERO_ONLY){
							tsw.print(scaf.name).tab().print(sum/(float)bin, 2);
							tsw.tab().print(i+1).tab().print(running).nl();
						}
					}
					running+=bin;
					nextPos+=binsize;
					lastPos=i;
					sum=0;
				}
			}
		}
		
		tsw.poisonAndWait();
	}
	
	
	
	/**
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param binsize Width of coverage bins in bp
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public static void writeCoveragePerBaseNormalized(String fname, ArrayList<Scaffold> list, double binsize, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		if(printHeader){
			String pound=(headerPound ? "#" : "");
			tsw.print(pound+"RefName\tBin\tCov\tPos\tRunningPos\n");
		}
		
		double running=0;
		double invbin=1.0/binsize;
		final double invbincount=1.0/NORMALIZE_LENGTH_BINS;
		for(Scaffold scaf : list){
			if(NORMALIZE_LENGTH_BINS>0){
				binsize=scaf.length*invbincount;
				invbin=1.0/binsize;
			}
			
			if(scaf.length>=binsize && scaf.length>=minscaf){
				long max=-1;
				
				final CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
				double lastPos=-1, nextPos=binsize-1;
				long sum=0;

				if(NORMALIZE_COVERAGE){
					for(int i=0; i<scaf.length; i++){
						int x=(ca==null ? 0 : ca.get(i));
						sum+=x;
						if(i>=nextPos){
							max=Tools.max(sum, max);
							running+=binsize;
							nextPos+=binsize;
							sum=0;
						}
					}
					lastPos=-1;
					nextPos=binsize-1;
					sum=0;
					assert(max>-1) : max;
				}
				max=Tools.max(max, 1);
				final double binmult=(NORMALIZE_COVERAGE ? 1d/max : invbin);
				
//				assert(false) : NORMALIZE_COVERAGE+", "+binmult+", "+invbin+", "+max+", "+binsize;
				
				final String formatString=NORMALIZE_COVERAGE ? "%s\t%d\t%.5f\t%d\t%d\n" : "%s\t%d\t%.2f\t%d\t%d\n";
				int bin=1;
				for(int i=0; i<scaf.length; i++){
					int x=(ca==null ? 0 : ca.get(i));
					sum+=x;
					if(i>=nextPos){
//						outstream.println(x+", "+i+", "+nextPos+", "+sum+", "+(sum*binmult));
						tsw.print(String.format(formatString, scaf.name, bin, sum*binmult, (i+1), (long)running));
						bin++;
						lastPos=i;
						running+=binsize;
						nextPos+=binsize;
						sum=0;
					}
				}
			}
		}
		
		tsw.poisonAndWait();
	}
	

	
	/**
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param binsize Width of coverage bins in bp
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public static void writeCoveragePerBaseNormalizedOverall(String fname, ArrayList<Scaffold> list, double binsize, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		
		assert(NORMALIZE_LENGTH_BINS>0) : "Must set 'normalizebins' flag to a positive integer.";
		double running=0;
		double invbin=1.0/binsize;
		long usedScafs=0;
		final double invbincount=1.0/NORMALIZE_LENGTH_BINS;

		double[] normalized=new double[NORMALIZE_LENGTH_BINS+1];
		double[] absolute=new double[NORMALIZE_LENGTH_BINS+1];
		
		for(Scaffold scaf : list){
			if(NORMALIZE_LENGTH_BINS>0){
				binsize=scaf.length*invbincount;
				invbin=1.0/binsize;
			}
			
			if(scaf.length>=binsize && scaf.length>=minscaf){
				usedScafs++;

				if(scaf.readhits>0){
					long max=-1;
					final CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
					double lastPos=-1, nextPos=binsize-1;
					long sum=0;

					{
						for(int i=0; i<scaf.length; i++){
							int x=(ca==null ? 0 : ca.get(i));
							sum+=x;
							if(i>=nextPos){
								max=Tools.max(sum, max);
								running+=binsize;
								nextPos+=binsize;
								sum=0;
							}
						}
						lastPos=-1;
						nextPos=binsize-1;
						sum=0;
						assert(max>-1) : max;
					}
					max=Tools.max(max, 1);
					final double binmult=1d/max;
					
					//				assert(false) : NORMALIZE_COVERAGE+", "+binmult+", "+invbin+", "+max+", "+binsize;

					int bin=1;
					for(int i=0; i<scaf.length; i++){
						int x=(ca==null ? 0 : ca.get(i));
						sum+=x;
						if(i>=nextPos){
							normalized[bin]+=(sum*binmult);
							absolute[bin]+=(sum*invbin);
							bin++;
							lastPos=i;
							running+=binsize;
							nextPos+=binsize;
							sum=0;
						}
					}
				}
			}
		}
		
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		if(printHeader){
			String pound=(headerPound ? "#" : "");
			tsw.print(pound+"RefName\tBin\tAbsCov\tNormCov\n");
		}
		double invScafs=1d/Tools.max(1, usedScafs);
		
		final double maxNorm=Tools.max(normalized);
		final double normMult=1/maxNorm;
		
		for(int bin=1; bin<normalized.length; bin++){
//			assert((absolute[bin]*invScafs)!=Double.NaN && (normalized[bin]*invScafs)!=Double.NaN) : invScafs+", "+absolute[bin]+", "+normalized[bin];
//			assert(false) : invScafs+", "+absolute[bin]+", "+normalized[bin]+", "+(absolute[bin]*invScafs)+", "+(normalized[bin]*invScafs);
			tsw.print(Tools.format("%s\t%d\t%.5f\t%.5f\n", "all", bin, absolute[bin]*invScafs, normalized[bin]*normMult));
		}
		
		tsw.poisonAndWait();
	}
	

	
	/**
	 * Write RPKM statistics.
	 */
	public static void writeRPKM(String out, ArrayList<String> in1, long readsIn, boolean printNonZeroOnly, ArrayList<Scaffold> list){
		if(out==null){return;}
		final ByteStreamWriter tsw=new ByteStreamWriter(out, overwrite, false, false);
		tsw.start();

		/* Count mapped reads */
		long mappedReads=0;
		long mappedFrags=0;
		for(Scaffold scaf : list){
			mappedReads+=scaf.readhits;
			mappedFrags+=scaf.fraghits;
		}
		mappedFrags/=2;
		
		/* Print header */
		tsw.print("#File\t"+(in1==null ? "" : in1)+"\n");
		tsw.print(Tools.format("#Reads\t%d\n",readsIn));
		tsw.print(Tools.format("#Mapped\t%d\n",mappedReads));
		tsw.print(Tools.format("#RefSequences\t%d\n",list.size()));
		tsw.print("#Name\tLength\tBases\tCoverage\tReads\tRPKM\tFrags\tFPKM\n");
		
		final float readMult=1000000000f/Tools.max(1, mappedReads);
		final float fragMult=1000000000f/Tools.max(1, mappedFrags);
		
		/* Print data */
		for(final Scaffold scaf : list){
			final long reads=scaf.readhits;
			final long frags=scaf.fraghits/2;
			final long bases=scaf.basehits;
			final String s=scaf.name;
			final int len=scaf.length;
			final double invlen=1.0/Tools.max(1, len);
			final double readMult2=readMult*invlen;
			final double fragMult2=fragMult*invlen;
			if(reads>0 || !printNonZeroOnly){
//				tsw.print(Tools.format("%s\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\n",s,len,bases,bases*invlen,reads,reads*readMult2,frags,frags*fragMult2));
				tsw.print(s).tab().print(len).tab().print(bases).tab();
				tsw.print(bases*invlen, 4).tab().print(reads).tab().print(reads*readMult2, 4).tab();
				tsw.print(frags).tab().print(frags*fragMult2, 4).nl();
			}
		}
		tsw.poisonAndWait();
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------          Threading           ----------------*/
	/*--------------------------------------------------------------*/
	
	class LoadThread extends Thread {
		
		LoadThread(int tid_){
			tid=tid_;
		}
		
		@Override
		public void run() {
			synchronized(this) {
				for(String next=getNext(); next!=POISON; next=getNext()){
					processViaStreamer(next);
				}
				addToQueue(POISON);//Send it to the next one
			}
		}
		
		private void processViaStreamer(String fname){
			final boolean processHeader=false;
//			assert(false) : processHeader+", "+fname;//Should be 
			SamLineStreamer ss=new SamLineStreamer(fname, streamerThreads, processHeader, maxReads);
			ss.start();
			ListNum<SamLine> ln=ss.nextLines();
			
//			false, true, true, true, stdin.sam
//			200
//			if(stdin) {
//				//Something odd is happening; ss.header is null and shouldn't be.
//				System.err.println((ln==null)+", "+(ss.header==null)+", "+stdin+", "+processHeader+", "+fname);
//				System.err.println(ln.size());
//				System.err.println(ss.header.size());
//			}
			if(processHeader) {
				processHeader(SamReadInputStream.getSharedHeader(false));
				SamReadInputStream.setSharedHeader(null);
			}
			for(; ln!=null && ln.size()>0; ln=ss.nextLines()){
				ArrayList<SamLine> list=(ln==null ? null : ln.list);
				for(SamLine sl : list){
					processSamLine(sl);
				}
			}
		}
		
		private void processHeader(ArrayList<byte[]> lines) {
			for(byte[] line : lines) {processHeaderLine(line);}
		}
		
		public boolean processSamLine(SamLine sl){
			readsProcessed++;
			basesProcessed+=sl.length();
			if(sl.duplicate() && !INCLUDE_DUPLICATES){return false;}
			
			if(border>0 || (trimq>=0 && (qtrimLeft || qtrimRight))){
				if(sl.mapped() && sl.seq!=null){
					int trimmed=trim(sl);
					if(sl.length()<1 || trimmed<0){return false;}//trimmed<0 implies everything was trimmed
				}
			}
			
			final int kmers=Tools.countKmers(sl.seq, k);
			final double cKmers=sl.qual==null ? kmers : Tools.countCorrectKmers(sl.qual, k);
			kmersProcessed+=kmers;
			correctKmers+=cKmers;
			final boolean properPair=(sl.hasMate() && sl.mapped() && sl.primary() && sl.properPair());
			if(PHYSICAL_COVERAGE && properPair){
				SamLine mate=pairTable.remove(sl.qname);
				if(mate==null){pairTable.put(sl.qname, sl);}
				else{
					final int start1=sl.start(INCLUDE_SOFT_CLIP, false);
					final int stop1=sl.stop(start1, INCLUDE_SOFT_CLIP, false);
					final int start2=mate.start(INCLUDE_SOFT_CLIP, false);
					final int stop2=mate.stop(start2, INCLUDE_SOFT_CLIP, false);
					final int strand=(sl.pairnum()==0 ? sl.strand() : mate.strand());
					final int length=USE_TLEN ? sl.tlen : Tools.max(stop1, stop2)-Tools.min(start1, start2)+1;
					mappedKmers+=kmers;
					addCoverage(sl.rnameS(), null, null, Tools.min(start1, start2), Tools.max(stop1, stop2), length, sl.mappedNonClippedBases(), strand, 2, sl.properPair(), sl);
				}
			}else if(sl.mapped() && (USE_SECONDARY || sl.primary()) && sl.mapq>=minMapq){
				assert(sl.seq!=null || sl.cigar!=null) : "This program requires bases or a cigar string for every sam line.  Problem line:\n"+sl+"\n";
//				assert(sl.seq!=null) : sl.toString();
				final int length=sl.length();
				final int start=sl.start(INCLUDE_SOFT_CLIP, false);
				final int stop=sl.stop(start, INCLUDE_SOFT_CLIP, false);
//				assert(false && length==stop-start+1) : length+", "+start+", "+stop+", "+(stop-start+1);
//				assert(false) : "'"+new String(sl.rname())+"', '"+sl.rnameS()+"'";
//				assert(false) : "'"+sl.rnameS()+"'";
				final byte[] match=(INCLUDE_DELETIONS ? null : sl.toShortMatch(true));
				mappedKmers+=kmers;
				return addCoverage(sl.rnameS(), sl.seq, match, start, stop, length, sl.mappedNonClippedBases(), sl.strand(), sl.hasMate() ? 1 : 2, sl.properPair(), sl);
			}
			return false;
		}
		
		public boolean addCoverage(final String scafName, final byte[] seq, byte[] match, final int start0, final int stop0, final int readlen, 
				final int nonClippedBases, final int strand, int incrementFrags, boolean properPair, SamLine sl){//sl is optional
			Scaffold scaf=table.get(scafName);
			if(scaf==null){
				if(EA){
					KillSwitch.kill("ERROR: A read was mapped to unknown reference sequence "+scafName);
				}else if(!warned){
					outstream.println("Warning: Can't find "+scafName+"\nThis scaffold will not be included.  Further warnings will be suppressed.");
					warned=true;
					error=true;
				}
				return false;
			}
			return addCoverage(scaf, seq, match, start0, stop0, readlen, nonClippedBases, strand, incrementFrags, properPair, sl);
		}
		
		public boolean addCoverage(final Scaffold scaf, final byte[] seq, byte match[], final int start0, final int stop0, final int readlen, final int nonClippedBases, 
				final int strand, int incrementFrags, boolean properPair, SamLine sl){//sl is optional
			if(scaf==null){
//				if(EA){
//					KillSwitch.kill("ERROR: Adding coverage to a null Scaffold "+(sl==null ? "" : sl.rnameS()));
//				}
				assert(false) : "Adding coverage to a null Scaffold.";
				return false;
			}
			
			final int start=Tools.max(start0, 0);
			final int stop=Tools.min(stop0, scaf.length-1);
			
			assert(start>=0 && stop>=0) : "\nAn error was encountered when processing a read. Output will not be valid.\n"+
				"\nscafName="+scaf.name+"\nseq="+new String(seq)+"\nstart="+start+
				"\nstop="+stop+"\nreadlen="+readlen+"\nstrand="+strand+"\nscaf.length="+scaf.length+"\nscaf="+scaf;
			
			mappedBases+=readlen;
			mappedNonClippedBases+=nonClippedBases;
			mappedBasesWithDels+=(stop-start+1);
			mappedReads++;
			if(properPair){properPairs++;}

			synchronized(scaf) {
				scaf.readhits++;
				scaf.fraghits+=incrementFrags;
				if(strand==1){scaf.readhitsMinus++;}

				if(seq!=null && scaf.basecount!=null){
					final long[] counts=scaf.basecount;
					for(int i=0; i<seq.length; i++){
						counts[charToNum[seq[i]]]++;
					}
				}
			}
			
			if(!INCLUDE_DELETIONS && !START_ONLY && !STOP_ONLY){
				assert(match!=null) : "Coverage excluding deletions cannot be calculated without a match string.";
				return addCoverageIgnoringDeletions(scaf, seq, match, start, stop, readlen, strand, incrementFrags);
			}
			
			synchronized(scaf) {
				final int basehits=stop-start+1;
				scaf.basehits+=basehits;
			}
			
			if(USE_COVERAGE_ARRAYS){
				if(atomic) {
					assert(scaf.obj0!=null) : "Preallocate this.";
					CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
					if(START_ONLY){ca.increment(start);}
					else if(STOP_ONLY){ca.increment(stop);}
					else{ca.incrementRange(start, stop);}
				}else {
					synchronized(scaf) {
						if(scaf.obj0==null){
							scaf.obj0=makeCA(scaf.length);
							if(STRANDED){scaf.obj1=makeCA(scaf.length);}
						}
						CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
						if(START_ONLY){ca.increment(start);}
						else if(STOP_ONLY){ca.increment(stop);}
						else{ca.incrementRange(start, stop);}
					}
				}
			}else if(USE_BITSETS){
				synchronized(scaf) {
					if(scaf.obj0==null){
						scaf.obj0=new BitSet(scaf.length);
						if(STRANDED){scaf.obj1=new BitSet(scaf.length);}
					}
					BitSet bs=(BitSet)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
					if(START_ONLY){bs.set(start);}
					else if(STOP_ONLY){bs.set(stop);}
					else{bs.set(start, stop+1);}
				}
			}
			
			return true;
		}
		
		private boolean addCoverageIgnoringDeletions(final Scaffold scaf, final byte[] seq, byte match[], final int start, final int stop, final int readlen, final int strand, int incrementFrags){
			assert(!INCLUDE_DELETIONS && !START_ONLY && !STOP_ONLY);
			assert(match!=null) : "Coverage excluding deletions cannot be calculated without a match string.";
			
			if(Read.isShortMatchString(match)){
				match=Read.toLongMatchString(match);
			}
			
			int basehits=0;
			synchronized(scaf) {
				if(USE_COVERAGE_ARRAYS){
					if(atomic) {
						if(scaf.obj0==null){
							scaf.obj0=makeCA(scaf.length);
							if(STRANDED){scaf.obj1=makeCA(scaf.length);}
						}
						CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
						addCoverageIgnoringDeletions(start, stop, match, ca);
					}else {
						synchronized(scaf) {
							if(scaf.obj0==null){
								scaf.obj0=makeCA(scaf.length);
								if(STRANDED){scaf.obj1=makeCA(scaf.length);}
							}
							CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
							addCoverageIgnoringDeletions(start, stop, match, ca);
						}
					}
				}else if(USE_BITSETS){
					synchronized(scaf) {
						if(scaf.obj0==null){
							scaf.obj0=new BitSet(scaf.length);
							if(STRANDED){scaf.obj1=new BitSet(scaf.length);}
						}
						BitSet bs=(BitSet)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
						addCoverageIgnoringDeletions(start, stop, match, bs);
					}
				}
				scaf.basehits+=basehits;
			}
			
			return true;
		}
		
		int addCoverageIgnoringDeletions(int start, int stop, byte[] match, CoverageArray ca) {
			int basehits=0;
			for(int rpos=start, mpos=0; mpos<match.length && rpos<=stop; mpos++){
				byte m=match[mpos];
				if(m=='m' || m=='S' || m=='N'){
					ca.increment(rpos, 1);
					basehits++;
					rpos++;
				}else if(m=='X' || m=='Y' || m=='C' || m=='I'){
					//do nothing
				}else if(m=='D'){
					rpos++;
				}else{
					assert(false) : "Unhandled symbol "+m;
				}
			}
			return basehits;
		}
		
		int addCoverageIgnoringDeletions(int start, int stop, byte[] match, BitSet bs) {
			int basehits=0;
			for(int rpos=start, mpos=0; mpos<match.length && rpos<=stop; mpos++){
				byte m=match[mpos];
				if(m=='m' || m=='S' || m=='N'){
					bs.set(rpos);
					basehits++;
					rpos++;
				}else if(m=='X' || m=='Y' || m=='C' || m=='I'){
					//do nothing
				}else if(m=='D'){
					rpos++;
				}else{
					assert(false) : "Unhandled symbol "+m;
				}
			}
			return basehits;
		}
		
		final int tid;
		
		public long mappedBases=0;
		public long mappedNonClippedBases=0;
		public long mappedBasesWithDels=0;
		public long mappedReads=0;
		public long properPairs=0;
		public long readsProcessed=0;
		public long basesProcessed=0;
		public long kmersProcessed=0;
		public long mappedKmers=0;
		public double correctKmers=0;
		
	}
	
	private String getNext() {
		String next=null;
		while(next==null) {
			try {
				next=queue.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return next;
	}
	
	private void addToQueue(String s) {
		while(s!=null) {
			try {
				queue.put(s);
				s=null;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	@Override
	public void accumulate(LoadThread t) {
		synchronized(t) {
			mappedBases+=t.mappedBases;
			mappedNonClippedBases+=t.mappedNonClippedBases;
			mappedBasesWithDels+=t.mappedBasesWithDels;
			mappedReads+=t.mappedReads;
			properPairs+=t.properPairs;
			readsProcessed+=t.readsProcessed;
			basesProcessed+=t.basesProcessed;
			kmersProcessed+=t.kmersProcessed;
			mappedKmers+=t.mappedKmers;
			correctKmers+=t.correctKmers;
		}
	}

	@Override
	public boolean success() {
		// TODO Auto-generated method stub
		return false;
	}
	
	private final ArrayBlockingQueue<String> queue=new ArrayBlockingQueue(8);
	private final String POISON="POISON_NOT_A_FILE";
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private boolean atomic=false;
	private boolean prealloc=false;
	private Class<? extends CoverageArray> caType;
	
	/** The list of all scaffolds */
	private ArrayList<Scaffold> list;
	/** Maps names to scaffolds */
	private HashMap<String, Scaffold> table;
	
	/** Mapping program name */
	private String program=null;
	/** Mapping program version */
	private String version=null;

	//Inputs
	/** Primary input files (typically sam) */
	public ArrayList<String> inputFiles=new ArrayList<String>();
	public int streams=-1;
	
	/** Optional, for calculating GC */
	public String reference=null;
	public String orffasta=null;

	//Outputs
	/** Coverage statistics, one line per scaffold */
	public String covstats=null;
	public String outorf=null;
	/** Coverage histogram, one line per depth and one point per base */
	public String histogram=null;
	/** Coverage with one line per base */
	public String basecov=null;
	/** Coverage with one file per scaffold */
	public String basecov_ps=null;
	/** Coverage with one line per bin */
	public String bincov=null;
	/** Coverage with one line per bin, normalized by length and/or height */
	public String normcov=null;
	/** Coverage with one line per bin, normalized by length and/or height, for combined reference */
	public String normcovOverall=null;
	/** rpkm/fpkm output, similar to Seal */
	public String outrpkm=null;
	
	/** Typically indicates that a header line was encountered in an unexpected place, e.g. with concatenated sam files. */
	private boolean error=false;
	
	private boolean warned=false;
	private final boolean EA=Shared.EA();
	
	/** Total length of reference */
	public long refBases=0;
	public long refKmers=0;
	public long mappedBases=0;
	public long mappedNonClippedBases=0;
	public long mappedBasesWithDels=0;
	public long mappedReads=0;
	public long properPairs=0;
	public long readsProcessed=0;
	public long basesProcessed=0;
	public long kmersProcessed=0;
	public long mappedKmers=0;
	public double correctKmers=0;
	public long totalCoveredBases1=0;
	public long totalCoveredBases2=0;
	public long scaffoldsWithCoverage1=0;
	public long scaffoldsWithCoverage2=0;
	public long totalScaffolds=0;

	public int k=0;
	
	//Don't reset these variables when clearing.
	public long maxReads=-1;
	public int initialScaffolds=4096;
	public int binsize=1000;
	public boolean bits32=false;
	public int minMapq=0;
	public int streamerThreads=2;
	public int minDepthToBeCovered=1;
	
	private boolean qtrimLeft=false;
	private boolean qtrimRight=false;
	private float trimq=-1;
	private final float trimE;
	private int border=0;
	
	/** Don't print coverage info for scaffolds shorter than this */
	public int minscaf=0;
	
	public HashMap<String, SamLine> pairTable=new HashMap<String, SamLine>();
	
	public PrintStream outstream=System.err;
	
	private boolean errorState=false;
	
	private final LineParser1 lp=new LineParser1('\t');
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	

	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Print verbose log messages */
	public static boolean verbose=false;
	/** Print timing messages */
	public static boolean verboseTime=false;
	/** Print the arguments to main */
	public static boolean printCommand=true;
	
	/** Print headers in output files */
	public static boolean printHeader=true;
	/** Prepend '#' symbol to header lines */
	public static boolean headerPound=true;
	/** Calculate standard deviation of coverage */
	public static boolean calcCovStdev=true;
	
	/** Window size to use when calculating average coverage,
	 * for detecting contiguous low-coverage areas */
	public static int LOW_COV_WINDOW=500;
	/** Min average coverage to not be classified as low-depth */
	public static double LOW_COV_DEPTH=5;
	/** Print number of bases below a certain average coverage in a window */
	public static boolean USE_WINDOW=false;
	
	public static int HISTMAX=-1;
	
	/** Track base composition of reads covering each scaffold */
	public static boolean COUNT_GC=true;
	/** Output in 2-column format ("#ID\tAvg_fold\n") */
	public static boolean TWOCOLUMN=false;
	/** Track coverage for strands independently */
	public static boolean STRANDED=false;
	/** Add scaffold information from the reference (in addition to sam header) */
	public static boolean ADD_FROM_REF=true;
	/** Only print scaffolds with nonzero coverage */
	public static boolean NONZERO_ONLY=false;
	/** Store coverage info in numeric arrays */
	public static boolean USE_COVERAGE_ARRAYS=true;
	/** Store coverage info in bitsets */
	public static boolean USE_BITSETS=false;
	/** Only print lines when coverage changes (for compatibility with Jigsaw) */
	public static boolean DELTA_ONLY=false;
	/** Process secondary alignments */
	public static boolean USE_SECONDARY=true;
	/** Include coverage of unsequenced middle portion of pairs */
	public static boolean PHYSICAL_COVERAGE=false;
	/** Use 'tlen' field when calculating physical coverage */
	public static boolean USE_TLEN=true;
	/** Abort on error; otherwise, errors may be ignored */
	public static boolean ABORT_ON_ERROR=true;
	/** Print coverage for the last bin of a scaffold, even if it is shorter than binsize */
	public static boolean KEEP_SHORT_BINS=true;
	/** Only track coverage for start location */
	public static boolean START_ONLY=false;
	/** Only track coverage for stop location */
	public static boolean STOP_ONLY=false;
	/** This appears to be the same as nonzeroonly... */
	public static boolean CONCISE=false;
	/** Normalize coverage by expression contig coverage as a fraction of its max coverage */
	public static boolean NORMALIZE_COVERAGE=false;
	/** Normalize contig length by binning into this many bins per contig */
	public static int NORMALIZE_LENGTH_BINS=100;
	/** Include soft-clipped bases in coverage */
	public static boolean INCLUDE_SOFT_CLIP=false;
	/** Include deletions/introns in coverage */
	public static boolean INCLUDE_DELETIONS=true;
	/** Include reads flagged as duplicates in coverage */
	public static boolean INCLUDE_DUPLICATES=true;

	public static boolean KEY_VALUE;
	
	/** Translation array for tracking base counts */
	private static final byte[] charToNum=AssemblyStats2.makeCharToNum();
	
	private static final int NOTHING_MODE=0, BITSET_MODE=1, ARRAY_MODE=2;
	
}
