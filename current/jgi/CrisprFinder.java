package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import ml.CellNet;
import ml.CellNetParser;
import ml.ScoreSequence;
import repeat.Palindrome;
import repeat.PalindromeFinder;
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
import structures.Crispr;
import structures.ListNum;
import structures.LongList;
import structures.LongLongListHashMap;
import structures.Range;
import structures.SeqCount;
import structures.SeqCountM;
import structures.SeqMap;
import structures.SeqPosM;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.CrisprTracker;
import tracker.EntropyTracker;
import tracker.PalindromeTracker;
import tracker.ReadStats;

/**
 * Designed to find Crispr repeats within reads.
 * 
 * Helpful Tadpole command:
 * java -ea -Xmx8g assemble.Tadpole mode=extend k=93 in=merged.fq.gz 
 * int=f out=extended.fq.gz el=50 er=50 ow ibb=f mce=5 bm2=5 bm1=30
 * 
 * java -Xmx8g clump.Clumpify in=merged.fq.gz out=clumped.fq.gz ecc passes=12 int=f ow conservative
 * 
 * java -Xmx8g assemble.Tadpole in=clumped.fq.gz out=eccClumped.fq.gz ecc conservative k=72 int=f
 * 
 * @author Brian Bushnell
 * @date August 29, 2023
 *
 */
public class CrisprFinder implements Accumulator<CrisprFinder.ProcessThread> {
	
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
		
		Shared.capThreads(64);//to prevent using 256 threads on Perlmutter head nodes
		
		//Create an instance of this class
		CrisprFinder x=new CrisprFinder(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CrisprFinder(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			extin=parser.extin;

			out1=parser.out1;
			out2=parser.out2;
			extout=parser.extout;
		}

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		
		ffoutu1=FileFormat.testOutput(outu1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffoutu2=FileFormat.testOutput(outu2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		long mask=0;
//		if(maskMiddle>0){
//			//This is clever but not needed here because only forward kmers are used.
//			if((kRepeat&1)!=(maskMiddle&1)) {
//				maskMiddle++;
//				System.err.println("Set middle mask to "+maskMiddle+" base"+
//						(maskMiddle==1 ? "" : "s")+".");
//			}
//			assert((kRepeat&1)==(maskMiddle&1));
//			assert(kRepeat>maskMiddle);
//			int bits=maskMiddle*2;
//			mask=~((-1L)<<bits);
//			mask<<=(kRepeat-maskMiddle);
//			mask=~mask;
//		}
		if(maskMiddle>0){
			assert(kRepeat>maskMiddle+1);
			int bits=maskMiddle*2;
			int shift=(kRepeat-maskMiddle)&(~1);//Equivalent to (x/2)*2
			mask=~((-1L)<<bits);
			mask<<=(shift);
//			mask=~mask;
		}
		minTailRepeat=Tools.min(minTailRepeat, minRepeat);
		midmaskRepeat=mask;
		repeatMap=(outRepeat==null && outRefAndRepeats==null && !forceMaps ? null 
				: new HashMap<SeqCount,SeqCountM>());
		spacerMap=(outSpacer==null && !forceMaps ? null : new HashMap<SeqCount,SeqCountM>());
		refUseSet=(ref==null || (outRef==null && outRefAndRepeats==null) ? null 
				: new HashMap<SeqCountM, AtomicLong>());
		incrementRef=refUseSet!=null;
		minPeriod=minSpacer+minRepeat;
		minCrispr=minSpacer+2*minRepeat;
		minTailPeriod=minSpacer+minTailRepeat;
		minTailCrispr=minSpacer+minRepeat+minTailRepeat;
		bruteForce=bruteForceLeft||bruteForceRight||bruteForceMiddle;
		loadNet();
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
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("merge")){
				merge=Parse.parseBoolean(b);
			}
			
			else if(a.equals("outu") || a.equals("outu1") || a.equals("outnonmatch") ||
					a.equals("outnonmatch1") || a.equals("outunnmatch") || a.equals("outunmatch1")
					|| a.equals("outunnmatched") || a.equals("outunmatched1")){
				outu1=b;
			}else if(a.equals("outu2") || a.equals("outnonmatch2") || a.equals("outunmatch2") ||
					a.equals("outnonmatched2") || a.equals("outunmatched2")){
				outu2=b;
			}
			
			
			else if(a.equals("masked")){
				masked=Parse.parseBoolean(b);
			}else if(a.equals("consensus")){
				consensus=Parse.parseBoolean(b);
			}else if(a.equals("minrepeats") || a.equals("repeats")){
				minRepeats=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("outCrisprHist") || a.equals("crisprhist") || a.equals("chist")){
				outCrisprHist=b;
			}else if(a.equalsIgnoreCase("outPalindromeHist") || a.equals("palhist") || a.equals("phist")){
				outPalindromeHist=b;
			}else if(a.equals("krepeat") || a.equals("kr")){
				kRepeat=Integer.parseInt(b);
			}else if(a.equals("k")){
				kRepeat=kRef=Integer.parseInt(b);
			}else if(a.equals("maskmiddlerepeat") || a.equals("mmrepeat") || a.equals("mmr") || a.equals("mm")){
				maskMiddle=(b==null ? 1 : Tools.startsWithDigit(b) ? 
					Integer.parseInt(b) : Parse.parseBoolean(b) ? 1 : 0);
			}else if(a.equals("rqhdist") || a.equals("rhdist")){
				rqhdist=Integer.parseInt(b);
			}else if(a.equals("minspacer")){
				minSpacer=Integer.parseInt(b);
			}else if(a.equals("maxspacer")){
				maxSpacer=Integer.parseInt(b);
			}else if(a.equals("minrepeat")){
				minRepeat=Integer.parseInt(b);
			}else if(a.equals("minrepeat0")){
				minRepeat0=Integer.parseInt(b);
			}else if(a.equals("maxrepeat")){
				maxRepeat=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minTailRepeat") || a.equalsIgnoreCase("minRepeatTail") || 
					a.equals("mtr") || a.equals("minrepeattip") || a.equals("mintiprepeat")){
				minTailRepeat=Integer.parseInt(b);
			}else if(a.equals("minrgc")){
				minRGC=Float.parseFloat(b);
			}else if(a.equals("maxrgc")){
				maxRGC=Float.parseFloat(b);
			}else if(a.equals("pad")){
				crisprOutPad=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxRepeatMismatches") || a.equals("rmismatches") || 
					a.equals("maxrmismatches") || a.equals("rmm") || a.equals("mrmm")){
				maxRepeatMismatches=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxRepeatMismatchesTail") || a.equals("rmmt")){
				maxRepeatMismatchesTail=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("rmmpad") || a.equalsIgnoreCase("rmmtpad") || a.equals("mrmmpad")){
				mrmmPad=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minOverlapConsensus")){
				minOverlapConsensus=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxTrimConsensus")){
				maxTrimConsensus=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("grow")){
				grow=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("growlookahead") || a.equalsIgnoreCase("lookahead")){
				growLookahead=Integer.parseInt(b);
			}
			
			else if(a.equals("minentropy") || a.equals("entropy")){
				minEntropy=Float.parseFloat(b);
			}else if(a.equals("entropywindow") || a.equals("ewindow")){
				entropyWindow=Integer.parseInt(b);
			}else if(a.equals("entropyk") || a.equals("kentropy") || a.equals("ke")){
				entropyK=Integer.parseInt(b);
			}
			
			else if(a.equals("maxns")){
				maxNs=Integer.parseInt(b);
			}else if(a.equals("maxpoly")){
				maxPoly=Integer.parseInt(b);
			}
			
			else if(a.equals("net")){
				netFile=b;
			}else if(a.equals("cutoff")){
				netCutoff=Float.parseFloat(b);
				setCutoff=true;
			}else if(a.equals("nn")){
				useNet=Parse.parseBoolean(b);
			}else if(a.equals("maxns")){
				maxNs=Integer.parseInt(b);
			}
			
			else if(a.equals("forcemaps")){
				forceMaps=Parse.parseBoolean(b);
			}else if(a.equals("queues") || a.equals("usequeues")){
				useQueues=Parse.parseBoolean(b);
			}else if(a.equals("outcrispr") || a.equals("outc")){
				outCrispr=b;
			}else if(a.equals("outrepeat") || a.equals("outr")){
				outRepeat=b;
			}else if(a.equals("outspacer") || a.equals("outs")){
				outSpacer=b;
			}else if(a.equals("outref")) {
				outRef=b;
			}else if(a.equals("outrr")) {
				outRefAndRepeats=b;
			}
			
			else if(a.equalsIgnoreCase("bruteForce")){
				bruteForceLeft=bruteForceRight=bruteForceMiddle=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("bruteForceLeft")){
				bruteForceLeft=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("bruteForceRight")){
				bruteForceRight=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("bruteForceMiddle")){
				bruteForceMiddle=Parse.parseBoolean(b);
			}

			else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("kref")){
				kRef=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maskMiddleRef") || a.equals("mmref")){
				maskMiddleRef=(b==null ? 1 : Tools.startsWithDigit(b) ? 
					Integer.parseInt(b) : Parse.parseBoolean(b) ? 1 : 0);
			}else if(a.equalsIgnoreCase("minRefMatches") || a.equals("minrefm")){
				minRefMatches=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxRefMismatches") || a.equalsIgnoreCase("refMismatches")
					|| a.equals("maxrefmm") || a.equals("refmm")){
				maxRefMismatches=Integer.parseInt(b);
			}else if(a.equals("refmmf") || a.equalsIgnoreCase("refmismatchfraction")){
				maxRefMismatchFraction=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minRefCount") || a.equals("minrefc") || a.equals("minrefcopies")
					 || a.equals("mincount")){
				minRefCount=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minCountPrint") || a.equals("mincounttoprint") || a.equals("minctp")){
				minCountToPrint=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("revertAlignmentFailures") || a.equals("revertaf")){
				revertAlignmentFailures=Parse.parseBoolean(b);
				if(revertAlignmentFailures) {discardAlignmentFailures=shrinkAlignmentFailures=false;}
			}else if(a.equalsIgnoreCase("shrinkAlignmentFailures") || a.equals("shrinkaf") 
					|| a.equalsIgnoreCase("trimAlignmentFailures") || a.equals("trimaf")){
				shrinkAlignmentFailures=Parse.parseBoolean(b);
				if(shrinkAlignmentFailures) {revertAlignmentFailures=discardAlignmentFailures=false;}
			}else if(a.equalsIgnoreCase("discardAlignmentFailures") || a.equals("discardaf")){
				discardAlignmentFailures=Parse.parseBoolean(b);
				if(discardAlignmentFailures) {revertAlignmentFailures=shrinkAlignmentFailures=false;}
			}else if(a.equalsIgnoreCase("ignoreAlignmentFailures") || a.equals("ignoreaf")){
				boolean x=Parse.parseBoolean(b);
				if(x) {discardAlignmentFailures=revertAlignmentFailures=shrinkAlignmentFailures=false;}
			}else if(a.equalsIgnoreCase("discardUnaligned") || a.equals("discardua") || a.equals("dua")){
				discardUnaligned=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("doubleAlign") || a.equals("doublea") || a.equals("aa")){
				doubleAlign=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("doubleFetch") || a.equals("doublef") || a.equals("ff")){
				doubleFetch=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("double")){
				doubleAlign=doubleFetch=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("minRefOverlapFractionQ") || a.equals("minrefofq")){
				minRefOverlapFractionQ=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minRefOverlapFractionR") || a.equals("minrefofr")){
				minRefOverlapFractionR=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minRefOverlap") || a.equals("minrefo")){
				minRefOverlap=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxRefSkew") || a.equals("maxreflop")){
				maxRefSkew=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxRefTrim") || a.equals("maxreft")){
				maxRefTrim=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("sortRefCandidates") || a.equals("sortref")){
				if("auto".equalsIgnoreCase(b)) {
					autoSortRefCandidates=true;
				}else {
					sortRefCandidates=Parse.parseBoolean(b);
					autoSortRefCandidates=false;
				}
			}
			
			else if(a.equalsIgnoreCase("printpals")){
				printPals=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("printscore")){
				printScore=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("printstats")){
				printStats=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("printhist")){
				printHist=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("fullLengthPalStats") || a.equals("fullpals")){
				fullLengthPalStats=Parse.parseBoolean(b);
			}
			
			else if(a.equals("minpalindrome") || a.equals("palindrome") || a.equals("minpal") || a.equals("plen")){
				minPal=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minpalindromeloop") || a.equals("minploop") || a.equals("minloop")){
				minLoop=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxPalindromeLoop") || a.equals("maxploop") || a.equals("maxloop")){
				maxLoop=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("requirePalindrome") || a.equals("requirepal") || a.equals("reqpal")){
				requirePalindrome=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("forcesymmetry") || a.equals("symmetric")){
				forceSymmetry=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("minPalTail") || a.equals("minptail") || a.equals("mintail")){
				minTail=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxPalTail") || a.equals("maxptail") || a.equals("maxtail")){
				maxTail=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxPalTailDif") || a.equals("maxptaildif") || a.equals("maxtaildif")){
				maxTailDif=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minPalMatches") || a.equals("pmatches") || a.equals("minpmatches")){
				minPalMatches=Integer.parseInt(b);
//				assert(false) : minPalMatches;
			}else if(a.equalsIgnoreCase("maxPalMismatches") || a.equals("pmismatches") || a.equals("maxpmismatches")){
				maxPalMismatches=Integer.parseInt(b);
			}
			
			else if(a.equalsIgnoreCase("annotate")){
				annotate=Parse.parseBoolean(b);
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
	
	/** Replace # with 1 and 2 in headers */
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

		//Do output file # replacement
		if(outu1!=null && outu2==null && outu1.indexOf('#')>-1){
			outu2=outu1.replace("#", "2");
			outu1=outu1.replace("#", "1");
		}
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		
		if("crisprs".equalsIgnoreCase(ref) && !new File(ref).exists()) {
			ref=Data.findPath("?crisprs.fa.gz");
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outu1, outu2, outRef, 
				outRefAndRepeats, outCrispr, outRepeat, outSpacer)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to some output files.\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2, ref)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outu1, outu2, ref, outRef, outRefAndRepeats, 
				outCrispr, outRepeat, outSpacer)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Make sure interleaving agrees with number of input and output files */
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
		
		if(ref!=null) {
			ArrayList<Read> refSeqs=FastaReadInputStream.toReads(ref, FileFormat.FASTA, -1);
			//This takes about 0.05 seconds for a 1Mbp dataset.
			refMap=SeqMap.load(refSeqs, kRef, maskMiddleRef, minRefCount, true, net0);
			
			if(incrementRef) {
				for(Read r : refSeqs) {
					SeqCountM scm=new SeqCountM(r.bases);
					refUseSet.put(scm, new AtomicLong(0));
				}
			}
//			assert()
			maxRefCount=refMap.sort();
			if(autoSortRefCandidates) {sortRefCandidates=maxRefCount>1;}
//			refMapSum=validateRefMap();
//			assert(refMapSum>=refMap.size() || refMap.isEmpty()) : refMapSum+", "+refMap.size();
		}
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris();
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros(ffout1, ffout2, cris.paired(), ordered);
		final ConcurrentReadOutputStream rosu=makeCros(ffoutu1, ffoutu2, cris.paired(), ordered);
		final ConcurrentReadOutputStream rosCrispr=makeCrosCrispr();
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads(cris, ros, rosu, rosCrispr);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		if(outCrisprHist!=null) {
			ReadWrite.writeStringInThread(cTracker.toString(), outCrisprHist, false);
		}
		
		if(outPalindromeHist!=null) {
			PalindromeTracker pt=fullLengthPalStats ? pTrackerFull : pTracker;
			ReadWrite.writeStringInThread(pt.toString(), outPalindromeHist, false);
		}
		
		PalindromeFinder palFinder=new PalindromeFinder(minPal, minLoop, maxLoop,
				minPalMatches, maxPalMismatches, minTail, maxTail, maxTailDif);
		if(outRepeat!=null) {
			ArrayList<SeqCountM> list=toSortedList(repeatMap, net0, minCountToPrint);
			printList(list, outRepeat, palFinder, net0, overwrite, minCountToPrint);
		}
		
		if(outSpacer!=null) {
			ArrayList<SeqCountM> list=toSortedList(spacerMap, null, 0);
			printList(list, outSpacer, null, null, overwrite, 0);
		}
		
		if(outRef!=null && refUseSet!=null) {
			ArrayList<SeqCountM> list=toSortedListSA(refUseSet, net0, 0);
			printList(list, outRef, (printPals ? palFinder : null), (printScore ? net0 : null), overwrite, 0);
		}
		
		if(outRefAndRepeats!=null) {
			@SuppressWarnings("unchecked")
			HashMap<SeqCountM,AtomicLong> map=(HashMap<SeqCountM, AtomicLong>) refUseSet.clone();
			{
				ArrayList<SeqCountM> list0=toSortedList(repeatMap, net0, minCountToPrint);
				for(SeqCountM scm : list0) {
					assert(scm.count>=minCountToPrint);
					if(!map.containsKey(scm)) {
						map.put(scm, new AtomicLong(scm.count));
					}
				}
			}
			ArrayList<SeqCountM> list=toSortedListSA(map, net0, 0);
			printList(list, outRefAndRepeats, palFinder, net0, overwrite, 0);
		}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros, rosu, rosCrispr);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		

		if(printStats || printHist) {
			outstream.println();
			ByteBuilder bbc=crisprStats(readsProcessed-readsMerged, printStats, printHist);
			outstream.println(bbc);
			if(minPal>0) {
				ByteBuilder bbp=palindromeStats(cTracker.crisprsFound-cTracker.partialTipRepeats,
						printStats, printHist);
				outstream.println(bbp);
			}
			if(merge) {
				ByteBuilder bb=new ByteBuilder();
				float mergeRate=readsMerged*200f/readsProcessed;
				bb.append("Merge Rate:       \t").append(mergeRate,2).percent().nl();
				outstream.println(bb);
			}
		}
		
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private ByteBuilder crisprStats(long readsTotal, boolean stats, boolean hist) {
		ByteBuilder bb=new ByteBuilder();
		CrisprTracker c=cTracker;
		final long crisprs=c.crisprsFound;
		final double invCrisprs=1.0/crisprs;
		final double invReads=1.0/readsTotal;
		long readsWithCrisprs=c.readsWithCrisprs;
		long trimmedByConsensus=c.trimmedByConsensus;
//		System.err.println("readsTotal="+readsTotal+", readsWithCrisprs="+readsWithCrisprs);
		long copies=c.copyList.sumHist();
		long clusters=c.copyList.sumLong()-c.copyList.get(0);
		double avgCopies=copies/(1.0*clusters);
		double fraction=readsWithCrisprs*invReads;
		double avgRlen=c.rlenList.meanHist();
		double avgSlen=c.slenList.meanHist();
		double avgRGC=c.rgcList.meanHist();
		double avgSGC=c.sgcList.meanHist();
		double avgMM=c.mismatchList.meanHist();
		double avgRefMM=c.refMismatchList.meanHist();
		double avgRefMMV=c.refMismatchListValid.meanHist();
		
		if(stats){
			bb.append("Crisprs Found:    \t").append(crisprs).nl();
			bb.append("Reads With Crisprs:\t").append(readsWithCrisprs).tab().append(fraction*100, 2).append("%").nl();
			if(ref!=null) {
				final double alignmentsPerRepeat=c.alignments/(double)c.alignmentRequested;
				final double queriesPerRepeat=refMap.setQueries.get()/(double)refMap.queries.get();
//				bb.append("Total Alignments: \t").append(c.alignments).nl();
				bb.append("Queries Per Repeat:\t").append(queriesPerRepeat, 2).nl();
				bb.append("Alignments Per Repeat:\t").append(alignmentsPerRepeat, 2).nl();
				bb.append("Successful Alignments:\t").append(c.alignedToRef).nl();
				bb.append("Modified By Ref:  \t").append(c.modifiedByRef).nl();
				String x=(revertAlignmentFailures ? "\t(Reverted)" :
					shrinkAlignmentFailures ? "\t(Shrank)" :
					discardAlignmentFailures ? "\t(Discarded)" : "");
				bb.append("Bad Alignment:    \t").append(c.failedAlignment).append(x).nl();
				bb.append("Avg Ref Mismatches:\t").append(avgRefMMV,4).nl();
			}
			if(consensus) {bb.append("Trimmed By Consensus:\t").append(trimmedByConsensus).nl();}
			if(bruteForce || true) {bb.append("Partial Tip Repeats:\t").append(c.partialTipRepeats).nl();}
			bb.append("Avg Repeat Copies:\t").append(avgCopies,2).nl();
			bb.append("Avg Repeat Length:\t").append(avgRlen,2).nl();
			bb.append("Avg Spacer Length:\t").append(avgSlen,2).nl();
			bb.append("Avg Repeat GC:    \t").append(avgRGC,2).percent().nl();
			bb.append("Avg Spacer GC:    \t").append(avgSGC,2).percent().nl();
			bb.append("Avg Repeat Mismatches:\t").append(avgMM,4).nl();
			
			if(repeatMap!=null) {bb.append("Unique Repeats:    \t").append(repeatMap.size()).nl();}
			if(spacerMap!=null) {bb.append("Unique Spacers:    \t").append(spacerMap.size()).nl();}
		}
		if(hist) {
			bb.append("#Crispr Histogram").nl();
			cTracker.appendTo(bb);
		}
		return bb;
	}
	
	private ByteBuilder palindromeStats(long crisprsTotal, boolean stats, boolean hist) {
		ByteBuilder bb=new ByteBuilder();
		PalindromeTracker p=fullLengthPalStats ? pTrackerFull : pTracker;
		final long pals=p.found;
		final double invPals=1.0/pals;
		final double invCrisprs=1.0/crisprsTotal;
		double fraction=pals*invCrisprs;
		double avgPlen=p.plenList.meanHist();
		double avgLoop=p.loopList.meanHist();
		double avgTail=p.tailList.meanHist();
		double avgTDif=p.tailDifList.meanHist();
		double avgM=p.matchList.meanHist();
		double avgMM=p.mismatchList.meanHist();
		
		if(stats){
			bb.append("Palindromes Found:\t").append(pals).append("\t(In Full-Length Repeats)").nl();
			bb.append("Crisprs With Pals:\t").append(fraction*100, 2).percent().nl();
			bb.append("Avg Pal Length: \t").append(avgPlen,2).nl();
			bb.append("Avg Loop Length:\t").append(avgLoop,2).nl();
			bb.append("Avg Tail Length:\t").append(avgTail,2).nl();
			bb.append("Avg Tail Len Dif:\t").append(avgTDif,2).nl();
			bb.append("Avg Matches:    \t").append(avgM,2).nl();
			bb.append("Avg Mismatches: \t").append(avgMM,2).nl();
		}
		if(hist){
			bb.append("#Palindrome Histogram").nl();
			p.appendTo(bb);
		}
		
		return bb;
	}
	
	private ConcurrentReadInputStream makeCris(){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}
	
	private static ConcurrentReadOutputStream makeCros(
			FileFormat ff1, FileFormat ff2, boolean pairedInput, boolean ordered){
		if(ff1==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

//		//Notify user of output mode
//		if(pairedInput && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
//			outstream.println("Writing interleaved.");
//		}

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	private ConcurrentReadOutputStream makeCrosCrispr(){
		if(outCrispr==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);

		FileFormat ffout=FileFormat.testOutput(outCrispr, FileFormat.FASTA, null, true, overwrite, append, ordered);
		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, 
			final ConcurrentReadOutputStream rosu, final ConcurrentReadOutputStream rosCrispr){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, rosu, rosCrispr, i));
		}
		
		//Start the threads and wait for them to finish
		ThreadWaiter.startThreads(alpt);
		
		if(useQueues) {
			//Ooops...  main thread can't do this easily since there are 2 queues to manage.
			//I could set some kind of flag though.
		}
		
		boolean success=ThreadWaiter.waitForThreadsToFinish(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt) {
//			System.err.println("Accumulating tid "+pt.tid);
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			readsMerged+=pt.readsMergedT;
			errorState|=(!pt.success);
			cTracker.add(pt.cTrackerT);
			pTracker.add(pt.palFinder.tracker);
			pTrackerFull.add(pt.palFinder.trackerFull);
			
			if(repeatMap!=null && !useQueues) {
				synchronized(pt.repeatMapT) {
					for(Entry<SeqCount, SeqCountM> e : pt.repeatMapT.entrySet()) {
						SeqCountM sc=e.getValue();
						addToMap(sc, repeatMap, true);
					}
				}
			}
			
			if(spacerMap!=null && !useQueues) {
				synchronized(pt.spacerMapT) {
					for(Entry<SeqCount, SeqCountM> e : pt.spacerMapT.entrySet()) {
						SeqCountM sc=e.getValue();
						addToMap(sc, spacerMap, true);
					}
				}
			}
		}
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final ArrayList<Crispr> findRanges(final Read r, final int k, final int hdist, 
			final long midMask, final int minSpacer, final int maxSpacer, final int minRepeat,
			final int maxRepeat, final int maxMismatches, final boolean grow, final int lookahead,
			final float minRGC, final float maxRGC, 
			final LongLongListHashMap kmerPositionMap, final LongList allKmers, final int[] acgtn){
		final int maxPeriod=maxRepeat+maxSpacer;
		final int minPeriod=minRepeat+minSpacer;
		final int kMinus1=k-1;
		final byte[] bases=r.bases;
		final long[] kmers=allKmers.array;
		final long size=allKmers.size;
		
		ArrayList<Crispr> pairs=null;
		for(int i=0; i<size; i++) {
			final long kmer=kmers[i];
			final long mmKmer=kmer|midMask;
			LongList positions=(kmer==-1 ? null : kmerPositionMap.get(mmKmer));
			if(kmer==-1 || positions==null || (hdist<1 && positions.size()<2)){
				continue;//possibly my first use of this keyword
			}
			
			if(verbose2) {System.err.println("A: i="+i+"; found list: "+positions);}
			final int idx1, idx2, aStop0=(i+kMinus1);
			if(positions.size()>1) {
				idx1=positions.findIndex(aStop0);
				idx2=idx1+1;
				if(verbose2) {System.err.println("B: idx1="+idx1+", idx2="+idx2);}
			}else {
				assert(hdist>0);
				//Now we find a new positions array with hdist.
				positions=fetchListWithHdist(kmer, midMask, hdist, aStop0, kmerPositionMap);
				if(positions==null) {continue;}
				idx1=-1;//Not really needed
				idx2=positions.findIndexAfter(aStop0);
				if(verbose2) {System.err.println("C: idx1="+idx1+", idx2="+idx2);}
			}
			if(idx2>=positions.size) {continue;}
			
			int aStop=aStop0;
			int bStop=(int)positions.get(idx2);
			final int period=bStop-aStop;
			int gap=period-k;
			if(verbose2) {System.err.println("D: aStop="+aStop+", bStop="+bStop+", period="+period+", gap="+gap);}
			if(period>maxPeriod) {continue;}
			else if(period<minPeriod){//Some short repeat region, skip ahead
				i+=k/2;
				if(verbose2) {System.err.println("Skipping to i="+i);}
				continue;
			}
			
			//Inserted block
			{//Extend
				int aStart=i, bStart=bStop-kMinus1;
				assert(aStart==aStop-kMinus1);
				int matches=(hdist<1 && midMask==0) ? k :
					countMatches(bases, aStart, aStop, bStart, bStop);
//				assert(false) : hdist+", "+midMask;
				int mismatches=k-matches;
				if(verbose2) {System.err.println("E: aStart="+aStart+", aStop="+aStop+", bStart="+bStart+", bStop="+bStop+", matches="+matches+", mismatches="+mismatches);}
				if(mismatches>0) {
					if(hdist>0){//shrink through mismatches
						//Shrink left
						final int lookaheadShrink=3;
						while(aStart<=aStop && mismatches>0) {
							if(bases[aStart]!=bases[bStart] || !AminoAcid.isFullyDefined(bases[aStart])) {
								aStart++;
								bStart++;
								mismatches--;
							}else {
								int x=checkRight(bases, aStart+1, bStart+1, lookaheadShrink);
								if(x<lookaheadShrink) {
									aStart+=x+1;
									bStart+=x+1;
									matches-=x;
									mismatches--;
								}else {break;}
							}
						}
						//Shrink right
						while(aStart<=aStop && mismatches>0) {
							if(bases[aStop]!=bases[bStop] || !AminoAcid.isFullyDefined(bases[aStop])) {
								aStop--;
								bStop--;
								mismatches--;
							}else {
								int x=checkLeft(bases, aStop-1, bStop-1, lookaheadShrink);
								if(x<lookaheadShrink) {
									aStop-=(x+1);
									bStop-=(x+1);
									matches-=x;
									mismatches--;
								}else {break;}
							}
						}
						if(aStart>aStop) {continue;}
					}else{//This block does not allow shrinking through mismatches
						//Shrink left
						while(bases[aStart]!=bases[bStart] || !AminoAcid.isFullyDefined(bases[aStart])) {
							aStart++;
							bStart++;
							mismatches--;
						}
						//Shrink right
						while(bases[aStop]!=bases[bStop] || !AminoAcid.isFullyDefined(bases[aStop])) {
							aStop--;
							bStop--;
							mismatches--;
						}
					}
				}
				if(verbose2) {System.err.println("F: aStart="+aStart+", aStop="+aStop+", bStart="+bStart+", bStop="+bStop+", matches="+matches+", mismatches="+mismatches);}
				
				//Extend left
				for(aStart--, bStart--; aStart>=0 && bases[aStart]==bases[bStart] && 
						mismatches<=maxMismatches && AminoAcid.isFullyDefined(bases[aStart]);
						aStart--, bStart--) {
					matches++;
				}
				aStart++;
				bStart++;
				if(verbose2) {System.err.println("G: aStart="+aStart+", aStop="+aStop+", bStart="+bStart+", bStop="+bStop+", matches="+matches+", mismatches="+mismatches);}
				
				//Extend right
				for(aStop++, bStop++; bStop<bases.length && bases[aStop]==bases[bStop] && 
						mismatches<=maxMismatches && AminoAcid.isFullyDefined(bases[aStop]);
						aStop++, bStop++) {
					matches++;
				}
				aStop--;
				bStop--;
				if(verbose2) {System.err.println("H: aStart="+aStart+", aStop="+aStop+", bStart="+bStart+", bStop="+bStop+", matches="+matches+", mismatches="+mismatches);}
				if(aStop>=bStart) {
					if(aStop>aStop0) {
						i+=(aStop-aStop0);//Skip the extension kmers
					}else {
						i+=k/2;
					}
					if(verbose2) {System.err.println("Skipping to i="+i);}
					continue;
				}
				
				while(grow && mismatches<maxMismatches) {
					//Can assert here that the next bases are unequal
					//NOTE: These two assertions fired on long homopolymers.
					assert(aStart<=0 || bases[aStart-1]!=bases[bStart-1] || 
							symbolToNumber[bases[aStart-1]]<0) : 
						aStart+", "+Character.toString((char)(bases[aStart-1]))+", "+
						Character.toString((char)(bases[bStart-1]))+"\n"+
						new String(bases, (Tools.max(aStart, 0)), 
								Tools.min(bStop, bases.length-1)-Tools.max(0, aStart)+1)+"\n"+
						new String(bases, (Tools.max(aStart-50, 0)), 
								Tools.min(bStop+50, bases.length-1)-Tools.max(0, aStart-50)+1)+"\n"+
								aStart+"-"+aStop+", "+bStart+"-"+bStop+", len="+bases.length;
					assert(bStop>=bases.length-1 || bases[aStop+1]!=bases[bStop+1] || 
							symbolToNumber[bases[aStop+1]]<0) : 
						aStart+", "+Character.toString((char)(bases[aStop+1]))+", "+
						Character.toString((char)(bases[bStop+1]))+"\n"+
						new String(bases, (Tools.max(aStart, 0)), 
								Tools.min(bStop, bases.length-1)-Tools.max(0, aStart)+1)+"\n"+
						new String(bases, (Tools.max(aStart-50, 0)), 
								Tools.min(bStop+50, bases.length-1)-Tools.max(0, aStart-50)+1)+"\n"+
								aStart+"-"+aStop+", "+bStart+"-"+bStop+", len="+bases.length;
					//Technically I only need to do one of these scans per loop
					int left=checkLeft(bases, aStart-2, bStart-2, 999);
					int right=checkRight(bases, aStop+2, bStop+2, 999);
					if(left>=right && left>=lookahead) {
						aStart=aStart-left-1;
						bStart=bStart-left-1;
						mismatches++;
						matches+=left;
					}else if(right>=left && right>=lookahead) {
						aStop=aStop+right+1;
						bStop=bStop+right+1;
						mismatches++;
						matches+=right;
					}else {break;}
				}
				
//				//Extend left through mismatches
//				while(aStart>0 && mismatches<=maxMismatches) {
//					int x=symbolToNumber[bases[aStart-1]], y=symbolToNumber[bases[bStart-1]];
//					if(x==y && x>=0) {
//						aStart--;
//						bStart--;
//						matches++;
//					}else if(grow && mismatches<maxMismatches && 
//							checkLeft(bases, aStart-2, bStart-2, lookahead)==lookahead){
//						aStart--;
//						bStart--;
//						mismatches++;
//					}else {
//						break;
//					}
//				}
//				
//				//Extend right through mismatches
//				while(bStop<bases.length-1 && mismatches<=maxMismatches) {
////					int x=symbolToNumber[bases[aStop+1]], y=symbolToNumber[bases[bStop+1]];
//					final int x=bases[aStop+1], y=bases[bStop+1];
//					if(x==y && symbolToNumber[x]>=0) {
//						aStop++;
//						bStop++;
//						matches++;
//					}else if(grow && mismatches<maxMismatches && 
//							checkRight(bases, aStop+2, bStop+2, lookahead)==lookahead){
//						aStop++;
//						bStop++;
//						mismatches++;
//					}else {
//						break;
//					}
//				}
				
				gap=bStart-aStop-1;
				int repeat=aStop-aStart+1;
//				assert(mismatches<maxMismatches) : new RangePair(aStart, aStop, bStart, bStop).toString(bases);
				if(verbose2) {System.err.println("I: gap="+gap+", repeat="+repeat+", minSpacer="+minSpacer+", maxSpacer="+maxSpacer+", minRepeat="+minRepeat+", maxRepeat="+maxRepeat);}
				if(gap>=minSpacer && gap<=maxSpacer && repeat>=minRepeat && repeat<=maxRepeat) {
					final float rgc=Tools.calcGC(bases, aStart, aStop);
					final float sgc=Tools.calcGC(bases, aStop+1, bStart-1);
					if(rgc>=minRGC && rgc<=maxRGC && sgc>=0.05f && sgc<=0.95f 
							&& mismatches<=maxMismatches) {
						final Crispr rp=new Crispr(aStart, aStop, bStart, bStop);
//						System.err.println("Created "+rp.toString(bases));
						rp.matches=matches;
						rp.mismatches=mismatches;
						
						//Not generally needed unless hdist was enabled
//						if(rp.mismatches>0) {
//							int x=shrink(bases, rp);
//							gap=rp.gap();
//							repeat=rp.maxLength();
//							assert(bases[rp.a.a]==bases[rp.b.a]);
//							assert(bases[rp.a.b]==bases[rp.b.b]);
//							assert(rp.a.a==0 || bases[rp.a.a+1]==bases[rp.b.a+1]) : rp.toString(bases);
//							assert(rp.b.b>=bases.length-1 || bases[rp.a.b-1]==bases[rp.b.b-1]) : rp.toString(bases);
//							assert(rp.a.a==0 || bases[rp.a.a+2]==bases[rp.b.a+2]) : rp.toString(bases);
//							assert(rp.b.b>=bases.length-1 || bases[rp.a.b-2]==bases[rp.b.b-2]) : rp.toString(bases);
//						}
						
						//No longer needed since grow happens earlier
//						if(grow && maxMismatches>rp.mismatches) {
//							int x=grow(bases, rp, maxMismatches, lookahead);
//							gap=rp.gap();
//							repeat=rp.maxLength();
//						}
						
						if(gap>=minSpacer && gap<=maxSpacer && repeat>=minRepeat && repeat<=maxRepeat) {
							if(pairs==null) {pairs=new ArrayList<Crispr>();}
							pairs.add(rp);
						}
					}
				}
				if(aStop>aStop0) {
					i+=(aStop-aStop0);//Skip the extension kmers
					if(verbose2) {System.err.println("Skipping to i="+i);}
				}
			}
		}
		return pairs;
	}
	
	//TODO: Add capability to trim using an Alignment, particularly where thre is a 3-way mismatch
	private static int shrink(byte[] bases, Crispr rp) {
		return shrinkLeft(bases, rp)+shrinkRight(bases, rp);
	}
	
	//Trims mismatched bases off the ends
	private static int shrinkLeft(byte[] bases, Crispr rp) {
		Range left=rp.a, right=rp.b;
		if(left.a<=0) {return 0;}
		int shrank=0;
		for(boolean loop=true; loop; ) {
			loop=false;
			while(left.a<=left.b && right.a<bases.length) {
				int x=symbolToNumber[bases[left.a]], y=symbolToNumber[bases[right.a]];
				if(x==y || x<0 || y<0) {break;}
				left.a++;
				right.a++;
				shrank++;
				rp.mismatches--;
			}
			if(left.a+1<left.b && right.a<bases.length-2) {
				assert(right.a+1<bases.length) : rp+", "+bases.length;
				int x=symbolToNumber[bases[left.a+1]], y=symbolToNumber[bases[right.a+1]];
				if(x==y || x<0 || y<0) {break;}
				left.a+=2;
				right.a+=2;
				shrank+=2;
				loop=true;
				rp.matches--;
				rp.mismatches--;
			}else if(left.a+2<left.b && right.a<bases.length-3) {
				assert(right.a+1<bases.length) : rp+", "+bases.length;
				int x=symbolToNumber[bases[left.a+2]], y=symbolToNumber[bases[right.a+2]];
				if(x==y || x<0 || y<0) {break;}
				left.a+=3;
				right.a+=3;
				shrank+=3;
				loop=true;
				rp.matches-=2;
				rp.mismatches--;
			}
		}
		return shrank;
	}
	
	private static int shrinkRight(byte[] bases, Crispr rp) {
		Range left=rp.a, right=rp.b;
		if(right.b>=bases.length-1) {return 0;}
		int shrank=0;
		for(boolean loop=true; loop; ) {
			loop=false;
			while(right.b>=right.a && left.b>=0) {
				int x=symbolToNumber[bases[left.b]], y=symbolToNumber[bases[right.b]];
				if(x==y || x<0 || y<0) {break;}
				left.b--;
				right.b--;
				shrank++;
				rp.mismatches--;
			}
			if(right.b>right.a+1 && left.b>1) {
				int x=symbolToNumber[bases[left.b-1]], y=symbolToNumber[bases[right.b-1]];
				if(x==y || x<0 || y<0) {break;}
				left.b-=2;
				right.b-=2;
				shrank+=2;
				loop=true;
				rp.matches--;
				rp.mismatches--;
			}else if(right.b>right.a+2 && left.b>2) {
				int x=symbolToNumber[bases[left.b-2]], y=symbolToNumber[bases[right.b-2]];
				if(x==y || x<0 || y<0) {break;}
				left.b-=3;
				right.b-=3;
				shrank+=3;
				loop=true;
				rp.matches-=2;
				rp.mismatches--;
			} 
		}
		return shrank;
	}
	

	private static int grow(byte[] bases, Crispr rp, int maxMismatches, int lookahead) {
		return growLeft(bases, rp, maxMismatches, lookahead)+growRight(bases, rp, maxMismatches, lookahead);
	}
	
	//Extend through mismatches
	private static int growLeft(byte[] bases, Crispr rp, int maxMismatches, int lookahead) {
		Range left=rp.a, right=rp.b;
		final int a0=left.a, b0=left.b, mm0=rp.mismatches;
//		System.err.print("a");
		if(a0<=0) {return 0;}
		
		while(left.a>0 && rp.mismatches<=maxMismatches) {
			int x=symbolToNumber[bases[left.a-1]], y=symbolToNumber[bases[right.a-1]];
			if(x==y && x>=0) {
				left.a--;
				right.a--;
				rp.matches++;
//				System.err.print("b");
			}else if(rp.mismatches<maxMismatches && checkLeft(bases, left.a-2, right.a-2, lookahead)==lookahead){
				left.a--;
				right.a--;
				rp.mismatches++;
//				System.err.print("c");
			}else {
//				System.err.print("d");
				break;
			}
		}
//		if(left.a<a0 && rp.mismatches>mm0) {shrinkLeft(bases, rp);}
		assert(rp.mismatches<=maxMismatches);
		assert(rp.matches+rp.mismatches==rp.minLength());
		return a0-left.a;
	}
	
	//Extend through mismatches
	private static int growRight(byte[] bases, Crispr rp, int maxMismatches, int lookahead) {
		Range left=rp.a, right=rp.b;
		final int a0=left.a, b0=left.b, mm0=rp.mismatches;
		if(b0>=bases.length-1) {return 0;}
		
		while(right.b<bases.length-1 && rp.mismatches<=maxMismatches) {
			int x=symbolToNumber[bases[left.b+1]], y=symbolToNumber[bases[right.b+1]];
			if(x==y && x>=0) {
				left.b++;
				right.b++;
				rp.matches++;
			}else if(rp.mismatches<maxMismatches && checkRight(bases, left.b+2, right.b+2, lookahead)==lookahead){
				left.b++;
				right.b++;
				rp.mismatches++;
			}else {
				break;
			}
		}
//		if(left.b>b0 && rp.mismatches>mm0) {shrinkRight(bases, rp);}
		assert(rp.mismatches<=maxMismatches);
		assert(rp.matches+rp.mismatches==rp.minLength());
		return left.b-b0;
	}
	
	private static int checkLeft(final byte[] bases, int a, int b, final int limit) {
		int matches=0;
		for(; matches<limit && a>=0; matches++, a--, b--) {
			final byte x=bases[a], y=bases[b];
			if(x!=y || symbolToNumber[x]<0) {break;}
		}
		return matches;
	}
	
	private static int checkRight(final byte[] bases, int a, int b, final int limit) {
		int matches=0;
		for(; matches<limit && b<bases.length; matches++, a++, b++) {
			final byte x=bases[a], y=bases[b];
			if(x!=y || symbolToNumber[x]<0) {break;}
		}
		return matches;
	}
	
	public static final LongList fetchListWithHdist(long kmer, long midMask, int hdist, int pos, 
			LongLongListHashMap map) {
		LongList list=map.get(kmer&midMask);
		assert(list!=null);
		if(list.size>1) {return list;}//TODO: Ultimately, findIndexAfter should be used
		//Now generate mutants
		assert(false) : "TODO";
		return null;
	}
	
	public static final int countMatches(final byte[] s, int a1, final int b1, int a2, final int b2) {
//		return Vector.countMatches(s, s, a1, b1, a2, b2);
		int matches=0;
		for(; a1<=b1 && a2<=b2; a1++, a2++) {
			matches+=(s[a1]==s[a2] && AminoAcid.isFullyDefined(s[a1])) ? 1 : 0;
		}
		return matches;
	}
	
	public static final int fillMap(final Read r, int k, long midMask,
			LongLongListHashMap kmerPositionMap, LongList allKmers){
		if(r==null || r.length()<k){return 0;}
		assert(kmerPositionMap.isEmpty());
		assert(allKmers.isEmpty());
		
		final byte[] bases=r.bases;
		final int shift=2*k;
//		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		final int kMinus1=k-1;
		long kmer=0;
//		long rkmer=0;
		int repeats=0;
		int len=0;
		
		final int start=0;
		final int stop=bases.length;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
//		for(int i=start; i<stop; i++){
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=symbolToNumber[b];
//			long x2=symbolToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
//			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
//				rkmer=0;
			}else{len++;}
			if(len>=k){
				long mmKmer=kmer|midMask;
				boolean added=kmerPositionMap.put(mmKmer, i);
				allKmers.add(kmer);
				if(added) {
//					uniqueKmers.add(kmer);
				}else {
					repeats++;
//					repeatKmers.add(kmer);
				}
			}else if(i>=kMinus1) {
				allKmers.add(-1L);
			}
		}
		return repeats;
	}
	
	public static void findPalindromes(byte[] bases, ArrayList<Crispr> pairs, 
			PalindromeFinder palFinder, boolean forceSymmetry, boolean requirePalindrome, int minRepeat) {
		if(bases==null || pairs==null || pairs.isEmpty()) {return;}
		int removed=0;
		for(int i=0; i<pairs.size(); i++) {
			Crispr pair=pairs.get(i);
			Range ra=pair.a, rb=pair.b;
			Palindrome pa=(ra.length()>=rb.length() ? palFinder.longestPalindrome(bases, ra.a, ra.b) : null);
			Palindrome pb=(rb.length()>=ra.length() ? palFinder.longestPalindrome(bases, rb.a, rb.b) : null);
			if(pa!=null) {palFinder.tracker.add(pa, ra.a, ra.b);}
			else if(pb!=null) {palFinder.tracker.add(pb, rb.a, rb.b);}
			pair.pa=pa;
			pair.pb=pb;
			
			if(pair.internal(bases.length)) {//If not touching the edge
				if(pa!=null) {palFinder.trackerFull.add(pa, ra.a, ra.b);}
				else if(pb!=null) {palFinder.trackerFull.add(pb, rb.a, rb.b);}

				//Don't do this, or if you do, fix the palindromes.
				if(forceSymmetry && pa!=null && ra.length()==rb.length()) {
					int left=pa.a-ra.a;
					int right=ra.b-pa.b;
					if(left>right) {//Trim left
						int dif=left-right;
						pair.a.a+=dif;
						pair.b.a+=dif;
					}else if(right>left) {
						int dif=right-left;
						pair.a.b-=dif;
						pair.b.b-=dif;
					}
				}
			}
			
			if(pa==null && pb==null && requirePalindrome) {
//				System.err.println("Removed "+rng);
				pairs.set(i, null);
				removed++;
			}
		}
		if(removed>0) {Tools.condenseStrict(pairs);}
//		System.err.println("Returning "+pairs);
	}

	private final int scoreRepeats(byte[] bases, ArrayList<Crispr> pairs, float[] vec, CellNet net) {
		int removed=0;
		for(int i=0; i<pairs.size(); i++) {
			Crispr rp=pairs.get(i);
			if(rp.internal(bases.length) || rp.lengthDif()>0 || rp.spans(bases.length)) {
				//Don't bother scoring edge repeats unless it spans the whole thing
				float f=scoreRepeat(bases, rp, vec, net);
				if(f<netCutoff) {
					pairs.set(i, null);
					removed++;
				}
			}
		}
		if(removed>0) {Tools.condenseStrict(pairs);}
		return removed;
	}
	
	private final float scoreRepeat(byte[] bases, Crispr rp, float[] vec, CellNet net) {
		if(rp.sameLength() && rp.mismatches==0) {
			return rp.scoreA=rp.scoreB=ScoreSequence.score(bases, vec, 0, net, rp.a.a, rp.a.b);
		}
		if(rp.a.length()>=rp.b.length()) {
			rp.scoreA=ScoreSequence.score(bases, vec, 0, net, rp.a.a, rp.a.b);
		}
		if(rp.b.length()>=rp.a.length()) {
			rp.scoreB=ScoreSequence.score(bases, vec, 0, net, rp.b.a, rp.b.b);
		}
		return rp.maxScore();
	}
	
	private static int addToMap(SeqCountM s, HashMap<SeqCount,SeqCountM> map, boolean copy) {
//		synchronized(map) {
//			synchronized(s) {
				SeqCountM old=map.get(s);
				if(old==null) {
					final SeqCountM clone=(copy ? s.clone() : s);
//					synchronized(clone) {
						//This is crucial!
						//Key must NOT be the same as value since SeqCount is mutable
						old=map.put(new SeqCount(s.bases), clone);

						assert(old==null);
						return 1;
//					}
				}else {
//					synchronized(old) {
						old.add(s);
						assert(old.equals(s));
						return 0;
//					}
				}
//			}
//		}
	}
	
	private static ArrayList<SeqCountM> toSortedList(HashMap<SeqCount,SeqCountM> map, 
			CellNet net, long minCount) {
		ArrayList<SeqCountM> list=new ArrayList<SeqCountM>(map.size());
		for(Entry<SeqCount, SeqCountM> e : map.entrySet()) {
			SeqCountM value=e.getValue();
			if(value.count>=minCount) {list.add(value);}
		}
		if(net!=null) {
			float[] vec=new float[net.numInputs()];
			for(SeqCountM scm : list) {
				scm.score=ScoreSequence.score(scm.bases, vec, 0, net);
			}
		}
		Collections.sort(list);
		Collections.reverse(list);
		return list;
	}
	
	private static ArrayList<SeqCountM> toSortedListSA(HashMap<SeqCountM,AtomicLong> map, 
			CellNet net, long minCount) {
		ArrayList<SeqCountM> list=new ArrayList<SeqCountM>(map.size());
		for(Entry<SeqCountM, AtomicLong> e : map.entrySet()) {
			final SeqCountM scm=e.getKey();
			final long count=e.getValue().get();
			scm.count=(int)count;//TODO: Count should be a long
			if(count>=minCount) {list.add(scm);}
		}
		if(net!=null) {
			float[] vec=new float[net.numInputs()];
			for(SeqCountM scm : list) {
				scm.score=ScoreSequence.score(scm.bases, vec, 0, net);
			}
		}
		Collections.sort(list);
		Collections.reverse(list);
		return list;
	}
	
	private static <X extends SeqCount> void printList(ArrayList<X> list, String fname, 
			PalindromeFinder palFinder, CellNet net, boolean overwrite, int minCount) {
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.FASTA, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		ByteBuilder bb=new ByteBuilder(240);
		float[] vec=(net==null ? null : new float[net.numInputs()]);
		int i=1;
		for(SeqCount s : list) {
//			Read r=new Read(s.bases);
			if(s.count()>=minCount) {
				bb.append('>').append(i).tab().append("count=").append(s.count());
				if(palFinder!=null) {
					Palindrome pal=palFinder.longestPalindrome(s.bases, 0, s.bases.length-1);
					if(pal!=null) {pal.appendTo(bb.tab(), 0, s.bases.length-1);}
				}
				if(net!=null) {
					float score=ScoreSequence.score(s.bases, vec, 0, net, false);
					bb.tab().append("score=").append(score,4);
				}
				bb.nl().append(s.bases).nl();
				bsw.print(bb);
				bb.clear();
			}
			i++;
		}
		bsw.poisonAndWait();
	}
	
	private static int cullLowCountRepeats(ArrayList<Crispr> list, int minRepeats) {
		final int minPairs=minRepeats-1;
		int removed=0, full=0, partial=0;
		if(minPairs>list.size()) {
			removed=list.size();
			list.clear();
			return removed;
		}
		
		Crispr prev=null;
		for(int i=0; i<list.size(); i++) {
			Crispr rp=list.get(i);
//			if(rp.lengthDif()!=0) {continue;}
			
			if(prev==null) {
				full=partial=0;
			}else if(prev.b.equals(rp.a)) {
				//do nothing
			}else {
				if(full+partial/2<minPairs) {//Remove this repeat
					for(int j=i-full-partial; j<i; j++) {
						Crispr x=list.set(j, null);
						assert(x!=null);
						removed++;
					}
				}else {
					//do nothing
				}
				full=partial=0;
			}
			
			prev=rp;
			if(rp.lengthDif()==0) {full++;}
			else {partial++;}
		}
		
		if(full+partial/2<minPairs) {//Remove the last repeat
			for(int j=list.size()-full-partial; j<list.size(); j++) {
				Crispr x=list.set(j, null);
				assert(x!=null);
				removed++;
			}
		}
		if(removed>0) {Tools.condenseStrict(list);}
		return removed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_,  
				final ConcurrentReadOutputStream rosu_,
				final ConcurrentReadOutputStream rosCrispr_, final int tid_){
			cris=cris_;
			ros=ros_;
			rosu=rosu_;
			rosCrispr=rosCrispr_;
			tid=tid_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			synchronized(this) {
				//The initializations are here instead of the constructor to use thread-local memory.
				repeatMapT=(useQueues || (outRepeat==null && !forceMaps) ? null : 
					new HashMap<SeqCount,SeqCountM>());
				spacerMapT=(useQueues || (outSpacer==null && !forceMaps) ? null : 
					new HashMap<SeqCount,SeqCountM>());

				Object o1=(repeatMapT==null ? this : repeatMapT);
				Object o2=(spacerMapT==null ? this : spacerMapT);
				synchronized(o1) {synchronized(o2) {
				kmerPositionMap=new LongLongListHashMap(400);
				allKmers=new LongList(400);
				acgtn=new int[5];
				palFinder=new PalindromeFinder(minPal, minLoop, maxLoop,
						minPalMatches, maxPalMismatches, minTail, maxTail, maxTailDif);
				cTrackerT=new CrisprTracker();
				refMapT=refMap;//SeqMap.load(ref, kRef, maskMiddleRef, minRefCount, true);
				eTracker=(minEntropy<=0 ? null : new EntropyTracker(entropyK, entropyWindow, false, minEntropy, true));
				net=(net0==null ? null : net0.copy(false));
				vec=(net0==null ? null : new float[net0.numInputs()]);
				
				processInner();

				//Do anything necessary after 
				}}

				//Indicate successful exit status
				success=true;
			}
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

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
			final ArrayList<Read> rejects=(rosu==null ? null : new ArrayList<Read>());
			final ArrayList<Read> crisprs=(rosCrispr==null ? null : new ArrayList<Read>());
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				Read r1=reads.get(idx);
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
				
				{
					//Reads are processed in this block.
					boolean keep=processReadPair(r1, r2, crisprs);
					
					if(keep){
						readsOutT+=r1.pairCount();
						basesOutT+=r1.pairLength();
					}else{
						if(rejects!=null) {rejects.add(r1);}
						reads.set(idx, null);
					}
				}
			}

			//Output reads to the output stream
			if(ros!=null){ros.add(reads, ln.id);}
			if(rosu!=null){rosu.add(rejects, ln.id);}
			if(rosCrispr!=null) {rosCrispr.add(crisprs, ln.id);}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		boolean processReadPair(Read r1, Read r2, ArrayList<Read> crisprs){
			final byte[] b1=r1.bases, q1=r1.quality;
			final byte[] b2=(r2==null ? null : r2.bases), q2=(r2==null ? null : r2.quality);
			if(merge) {
				Read r=BBMerge.tryToMerge(r1, r2);
				if(r!=null) {
					r1.bases=r.bases;
					r1.quality=r.quality;
					r2.bases=null;
					r2.quality=null;
					readsMergedT++;
				}
			}
			boolean found=processRead(r1, crisprs);
			found=processRead(r2, crisprs)|found;
			if(merge && r2!=null) {//Unmerge before output
				r1.bases=b1; r2.bases=b2; r1.quality=q1; r2.quality=q2;
			}
			return found;
		}
		
		boolean processRead(Read r, ArrayList<Read> crisprs) {
			if(r==null || r.length()<kRepeat){return false;}
			int repeats=fillMap(r, kRepeat, midmaskRepeat, kmerPositionMap, allKmers);
			if(repeats<1) {
				clearMap();
				return false;
			}
			
			if(r.length()<800) {//Decide whether to discard short, low-quality reads entirely
				if(eTracker!=null && eTracker.averageEntropy(r.bases, true)<minEntropy) {
					clearMap();
					return false;
				}
//				if(maxNs>=0 && maxNs<r.length() && r.countUndefined()>maxNs) {
//					clearMap();
//					return false;
//				}
//				if(maxPoly>1 && maxPoly<r.length() && r.longestHomopolymer()>maxPoly) {
//					clearMap();
//					return false;
//				}
			}
			
//			assert(false) : kmerPositionMap.toStringSetView();
			ArrayList<Crispr> pairs=findRanges(r, kRepeat, rqhdist, midmaskRepeat, minSpacer, maxSpacer, 
					minRepeat0, maxRepeat, maxRepeatMismatches, grow, growLookahead, minRGC, maxRGC,
					kmerPositionMap, allKmers, acgtn);
			
//			System.err.println("a");
			if(pairs!=null && !pairs.isEmpty()) {
//				for(RangePair rp : pairs) {assert(rp.gap()>=minSpacer) : rp.toString(r.bases);}
				while(removeMistakes(pairs, r.length())>0) {}
//				System.err.println("b");
//				for(RangePair rp : pairs) {assert(rp.gap()>=minSpacer) : rp.toString(r.bases);}
				if(ref!=null) {alignToRef(r.bases, pairs);}
				

//				assert(false) : "TODO: Shrink here if over max read mismatches.";
				
//				System.err.println("c");
//				for(RangePair rp : pairs) {assert(rp.gap()>=minSpacer) : rp.toString(r.bases);}
				if(masked) {shrinkToMasked(r.bases, pairs);}
//				System.err.println("d");
				if(consensus) {consensus(r.bases, pairs);}
//				System.err.println("e");
//				for(RangePair rp : pairs) {assert(rp.gap()>=minSpacer) : rp.toString(r.bases);}
				if(bruteForce) {//This could be wrapped in a while loop
					int x=bruteForce(r.bases, pairs);
//					System.err.println("e2");
					if(x>0) {
						Collections.sort(pairs);
						while(removeMistakes(pairs, r.length())>0) {}
//						assert(false) : pairs;
						//Can potentially trim more with the new information
						if(consensus) {consensus(r.bases, pairs);}
//						System.err.println("e3");
					}
				}
				
//				filterPairs(r.bases, pairs);
				
//				for(RangePair rp : pairs) {assert(rp.gap()>=minSpacer) : rp.toString(r.bases);}
//				System.err.println(pairs.size());
//				System.err.println(pairs.size());
//				assert(false) : minPal+", "+requirePalindrome;
			}
//			System.err.println("f");

			if(pairs!=null) {
				filterPairs(r.bases, pairs);
//				for(Crispr rp : pairs) {
//					assert(rp.a.length()>=minRepeat || rp.b.length()>=minRepeat);
//				}
			}
			if(pairs!=null && minRepeats>2) {
				cullLowCountRepeats(pairs, minRepeats);
			}
//			System.err.println("g");

			if(pairs!=null && !pairs.isEmpty() && net!=null && netCutoff>-1){//Score with net
				scoreRepeats(r.bases, pairs, vec, net);
			}
			
			if(pairs!=null && !pairs.isEmpty() && minPal>0){//Process palindromes
				findPalindromes(r.bases, pairs, palFinder, forceSymmetry, requirePalindrome, minRepeat);
			}
//			System.err.println("h");

			if(pairs!=null && !pairs.isEmpty()){//Annotate and output
				trackCopies(pairs);
				if(repeatMapT!=null || spacerMapT!=null) {addToMaps(pairs, r.bases);}
				ByteBuilder bb=new ByteBuilder(), bbTotal=new ByteBuilder();
				bbTotal.append(r.id);
				for(Crispr rp : pairs) {
					cTrackerT.add(rp, r.bases);
					if(annotate) {rp.appendTo(bbTotal.tab(),r.length(),null);}
					if(crisprs!=null) {
						final int a=Tools.max(rp.a.a-crisprOutPad, 0);
						final int b=Tools.min(rp.b.b+crisprOutPad, r.length()-1);
						final int leftPad=rp.a.a-a, rightPad=b-rp.b.b;
						byte[] bases=Arrays.copyOfRange(r.bases, a, b+1);
						byte[] quals=(r.quality==null ? null : 
							Arrays.copyOfRange(r.quality, a, b+1));
						bb.clear().append(r.id);
						if(annotate || true) {rp.appendTo(bb.tab(),r.length(),null);}
						Read crispr=new Read(bases, quals, bb.toString(), r.numericID);
						final int alen=rp.a.length();
						final int blen=rp.b.length();
//						assert(len==rp.b.length());
						for(int i=0; i<alen; i++) {
							bases[i+leftPad]=Tools.toLowerCase(bases[i+leftPad]);
						}
						for(int i=0; i<blen; i++) {
							int x=bases.length-1-rightPad-i;
							bases[x]=Tools.toLowerCase(bases[x]);
						}
						crisprs.add(crispr);
					}
				}
				r.id=bbTotal.toString();
			}
//			else{//Not sure why this was here but it messed up the histogram
//				cTracker.copyList.increment(0);
//			}
			clearMap();
			return (pairs!=null && !pairs.isEmpty());
		}
		
		int filterPairs(byte[] bases, ArrayList<Crispr> pairs) {
			int removed=0;
			for(int i=0; i<pairs.size(); i++) {
				Crispr rp=pairs.get(i);
				if(maxNs>=0 && maxNs<bases.length && Read.countUndefined(bases, rp.a.a, rp.b.b)>maxNs) {
					pairs.set(i, null);
					removed++;
//				}else if(maxPoly>1 && maxPoly<bases.length && Read.longestHomopolymer(bases, rp.a.a, rp.b.b)>maxPoly) {
				}else if(maxPoly>1 && maxPoly<bases.length && (Read.longestHomopolymer(bases, rp.a.a, rp.a.b)>maxPoly ||
						Read.longestHomopolymer(bases, rp.b.a, rp.b.b)>maxPoly)) {
					pairs.set(i, null);
					removed++;
				}else if(rp.a.length()<minRepeat && rp.b.length()<minRepeat /*&& (ref==null || discardAlignmentFailures)*/) {//TODO: Make sure this was 'modified by reference'
					pairs.set(i, null);
					removed++;
				}else if((rp.gap()<minSpacer || rp.gap()>maxSpacer) /*&& (ref==null || discardAlignmentFailures)*/) {
					pairs.set(i, null);
					removed++;
				}
			}
			if(removed>0) {Tools.condenseStrict(pairs);}
			return removed;
		}
		
		int removeMistakes(ArrayList<Crispr> list, int length) {
			int removed=0, swapped=0;
			for(int i=1; i<list.size(); i++) {
				Crispr left=list.get(i-1), right=list.get(i);
				if(left!=null && right!=null) {
					if(left.containsInGap(right.a) || left.containsInGap(right.b)) {
						list.set(i-1, null);
						removed++;
					}else if(right.containsInGap(left.a) || right.containsInGap(left.b)) {
						list.set(i, null);
						removed++;
					}else if(left.a.includes(right.a) && left.b.includes(right.b)) {
						list.set(i, null);
						removed++;
					}else if(right.a.includes(left.a) && right.b.includes(left.b)) {
						list.set(i-1, null);
						removed++;
					}else if(left.a.overlap(right.a)>=left.a.length()/2 && 
							left.b.overlap(right.b)>=left.b.length()/2) {
						if(left.minLength()>=right.minLength()) {
							list.set(i-1, null);
							removed++;
						}else {
							list.set(i, null);
							removed++;
						}
					}else if(left.compareTo(right)>0) {
						swapped++;
						list.set(i, left);
						list.set(i-1, right);
					}
				}
			}
			if(removed>0) {Tools.condenseStrict(list);}
			if(swapped>0) {Collections.sort(list);}
			return removed+swapped;
		}
		
		void trackCopies(ArrayList<Crispr> pairs) {
			if(pairs==null || pairs.isEmpty()) {return;}
			Crispr prev=null;
			int copies=0;
			for(Crispr p : pairs) {
				if(prev==null) {
					copies=2;
				}else if(prev.b.overlaps(p.a)) {
					copies++;
				}else{
					assert(copies!=0) : p+"\n"+pairs;
					cTrackerT.copyList.increment(copies);
					copies=2;
				}
				prev=p;
			}
			assert(copies!=0) : pairs;
			cTrackerT.copyList.increment(copies);
			cTrackerT.readsWithCrisprs++;
		}
		
		private <X> void put(X x, ArrayBlockingQueue<X> queue) {
			for(success=false; !success; ) {
				try {
					queue.put(x);
					success=true;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		void addToMaps(ArrayList<Crispr> pairs, byte[] bases) {
			final int blen=bases.length;
			if(pairs==null || pairs.isEmpty()) {return;}
			Crispr prev=null;
			for(Crispr p : pairs) {
				if(p.internal(blen) && p.sameLength()) {
//					assert(p.a.length()>=minRepeat || p.b.length()>=minRepeat) : p;
					if(repeatMapT!=null) {
						if(p.a.length()>=p.b.length() && 
								(prev==null || !prev.b.equals(p.a) || !prev.internal(blen))){
							assert(p.a.length()>=minRepeat) : p.a.length()+", "+p.b.length()+", "+minRepeat;
							boolean passesNet=(net==null || p.scoreA>=netCutoff);
							boolean passesPal=(minPal<1 || (p.pa!=null && p.pa.matches>=minPalMatches));
							
							if(passesNet && passesPal) {
								SeqCountM scm=new SeqCountM(bases, p.a.a, p.a.b+1);

								if(useQueues) {put(scm, repeatQueue);}
								else {addToMap(scm, repeatMapT, false);}
							}
						}
						if(p.b.length()>=p.a.length()) {
							assert(p.b.length()>=minRepeat);
							boolean passesNet=(net==null || p.scoreB>=netCutoff);
							boolean passesPal=(minPal<1 || (p.pb!=null && p.pb.matches>=minPalMatches));
							
							if(passesNet && passesPal) {
								SeqCountM scm=new SeqCountM(bases, p.b.a, p.b.b+1);

								if(useQueues) {put(scm, repeatQueue);}
								else {addToMap(scm, repeatMapT, false);}
							}
						}
					}
					if(spacerMapT!=null) {
						SeqCountM s3=new SeqCountM(bases, p.a.b+1, p.b.a);
						if(useQueues) {put(s3, spacerQueue);}
						else {addToMap(s3, spacerMapT, false);}
					}
					prev=p;
				}else{
					if(spacerMapT!=null && ((!p.sameLength() && p.touchesEdge(blen)) || 
							(p.spans(blen) && p.sameLength()))) {//Should be OK since it was extended
						SeqCountM s3=new SeqCountM(bases, p.a.b+1, p.b.a);
						if(useQueues) {put(s3, spacerQueue);}
						else {addToMap(s3, spacerMapT, false);}
					}
					prev=null;
				}
			}
		}
		
		@Deprecated
		public void shrinkToMasked(byte[] bases, ArrayList<Crispr> list) {
			if(list==null) {return;}
			int removed=0;
			for(int i=0; i<list.size(); i++) {
				Crispr p=list.get(i);
				while(p.a.a<=p.a.b && Tools.isUpperCase(bases[p.a.a])) {
					if(bases[p.a.a]==bases[p.b.a]) {p.matches--;}
					else {p.mismatches--;}
					p.a.a++;
					p.b.a++;
				}
				while(p.a.a<=p.a.b && Tools.isUpperCase(bases[p.a.b])) {
					if(bases[p.a.b]==bases[p.b.b]) {p.matches--;}
					else {p.mismatches--;}
					p.a.b--;
					p.b.b--;
				}
				int repeat=p.a.length();
				int gap=p.gap();
				if(gap<minSpacer || gap>maxSpacer || repeat<minRepeat || repeat>maxRepeat) {
					removed++;
					list.set(i, null);
				}
			}
			if(removed>0) {
				Tools.condenseStrict(list);
			}
//			Tools.toUpperCase(bases);
		}
		
		public void alignToRef(byte[] bases, ArrayList<Crispr> list) {
			int changes=0;
			int removed=0;
			for(int i=0; i<list.size(); i++) {
				final Crispr rp=list.get(i);
				final int a1=rp.a.a, b1=rp.a.b, a2=rp.b.a, b2=rp.b.b; //Initials
				int leftDist=a1, rightDist=bases.length-1-b2;
				
				final int mismatches=rp.mismatches;
				//TODO: This assertion fires when the sequence has lowercase bases.
				assert(mismatches==countMismatches(rp, bases)) : mismatches+", "+countMismatches(rp, bases);
//				System.err.println("\nldist="+leftDist+", rdist="+rightDist+", "+rp.toString(bases));
				if(doubleAlign && leftDist>8 && rightDist>8 && rp.sameLength()) {
					Alignment aln=doubleAlignRefToRange(bases, rp.a, rp.b, mismatches);
					boolean aligned=(aln!=null);//doubleAlignRefToRange(bases, rp.a, rp.b, mismatches);
					if(aln!=null) {rp.a.a=aln.a2; rp.a.b=aln.b2;}
					int aDif=rp.a.a-a1, bDif=rp.a.b-b1;
					if(aDif!=0 || bDif!=0) {
						if(rp.a.a<=0) {
							//revert.
							rp.a.a=a1;
							rp.a.b=b1;
						}else {
							rp.b.a+=aDif;
							rp.b.b+=bDif;
							rp.fixBounds(bases.length);
							countMatches(rp, bases);
							boolean fail=(rp.minLength()<minTailRepeat || rp.maxLength()<minRepeat 
									|| rp.gap()<minSpacer || rp.gap()>maxSpacer
									|| rp.mismatches>maxRepeatMismatches 
									|| (rp.lengthDif()>0 && rp.mismatches>maxRepeatMismatchesTail));
							if(fail) {cTrackerT.failedAlignment++;}
							if(fail && revertAlignmentFailures){
								rp.set(a1, b1, a2, b2);
								countMatches(rp, bases);
							}else if(fail && shrinkAlignmentFailures){
								shrink(bases, rp);
								countMatches(rp, bases);
								changes++;
							}else if(fail && discardAlignmentFailures){
								list.set(i, null);
								removed++;
							}else {
								cTrackerT.modifiedByRef++;
								changes++;
							}
						}
					}else if(discardUnaligned && !aligned) {
						list.set(i, null);
						removed++;
					}
				}else if(leftDist>=rightDist && leftDist>0) {//use left side
					Alignment aln=alignRefToRange(bases, rp.a, rp.b, mismatches);
					boolean aligned=(aln!=null);
					if(aln!=null) {rp.a.a=aln.a2; rp.a.b=aln.b2;}
					int aDif=rp.a.a-a1, bDif=rp.a.b-b1;
					if(aDif!=0 || bDif!=0) {
						if(rp.a.a<=0) {
							//revert.  These are just always going to be problematic since they span the read.
							rp.a.a=a1;
							rp.a.b=b1;
						}else {
							rp.b.a+=aDif;
							rp.b.b+=bDif;
							rp.fixBounds(bases.length);
							countMatches(rp, bases);
							boolean fail=(rp.minLength()<minTailRepeat || rp.maxLength()<minRepeat 
									|| rp.gap()<minSpacer || rp.gap()>maxSpacer
									|| rp.mismatches>maxRepeatMismatches 
									|| (rp.lengthDif()>0 && rp.mismatches>maxRepeatMismatchesTail));
							if(fail) {cTrackerT.failedAlignment++;}
							if(fail && revertAlignmentFailures){
								rp.set(a1, b1, a2, b2);
								countMatches(rp, bases);
							}else if(fail && shrinkAlignmentFailures){
								shrink(bases, rp);
								countMatches(rp, bases);
								changes++;
							}else if(fail && discardAlignmentFailures){
								list.set(i, null);
								removed++;
							}else {
								cTrackerT.modifiedByRef++;
								changes++;
							}
						}
					}else if(discardUnaligned && !aligned) {
						list.set(i, null);
						removed++;
					}
				}else if(rightDist>0){//use right side
					Alignment aln=alignRefToRange(bases, rp.b, rp.a, mismatches);
					boolean aligned=(aln!=null);
					if(aln!=null) {rp.b.a=aln.a2; rp.b.b=aln.b2;}
					int aDif=rp.b.a-a2, bDif=rp.b.b-b2;
					if(aDif!=0 || bDif!=0) {
						if(rp.b.b>=bases.length-1) {
							//revert
							rp.b.a=a2;
							rp.b.b=b2;
						}else {
							rp.a.a+=aDif;
							rp.a.b+=bDif;
							rp.fixBounds(bases.length);
							countMatches(rp, bases);
							boolean fail=(rp.minLength()<minTailRepeat || rp.maxLength()<minRepeat 
									|| rp.gap()<minSpacer || rp.gap()>maxSpacer
									|| rp.mismatches>maxRepeatMismatches 
									|| (rp.lengthDif()>0 && rp.mismatches>maxRepeatMismatchesTail));
							if(fail) {cTrackerT.failedAlignment++;}
							if(fail && revertAlignmentFailures){
								rp.set(a1, b1, a2, b2);
								countMatches(rp, bases);
							}else if(fail && shrinkAlignmentFailures){
								shrink(bases, rp);
								countMatches(rp, bases);
								changes++;
							}else if(fail && discardAlignmentFailures){
								list.set(i, null);
								removed++;
							}else {
								cTrackerT.modifiedByRef++;
								changes++;
							}
						}
					}else if(discardUnaligned && !aligned) {
						list.set(i, null);
						removed++;
					}
				}
			}
			if(removed>0) {Tools.condenseStrict(list);}
			if(changes>0) {
				removeShort(bases, list);
				Collections.sort(list);
				while(removeMistakes(list, bases.length)>0) {}
//				for(RangePair rp : list) {assert(rp.gap()>=minSpacer) : rp.toString(bases);}
			}
//			for(RangePair rp : list) {assert(rp.gap()>=minSpacer) : rp.toString(bases);}
		}
		
		public Alignment alignRefToRange(byte[] bases, Range r, Range alt, int mm) {
			cTrackerT.alignmentRequested++;
			final ArrayList<SeqPosM> list;
			int fetchMM=Tools.max(maxRefMismatches, Math.round(bases.length*maxRefMismatchFraction));
			if(doubleFetch && mm>0 && r.length()==alt.length()) {
				list=refMapT.doubleFetch(bases, r.a, r.b, alt.a, alt.b, minRefOverlap, fetchMM, 
						maxRefTrim, maxRefSkew, minRefOverlapFractionQ, sortRefCandidates);
			}else {
				list=refMapT.fetch(bases, r.a, r.b, minRefOverlap, fetchMM, 
						maxRefTrim, maxRefSkew, minRefOverlapFractionQ, sortRefCandidates);
			}
			if(list==null) {return null;}
			
			Alignment[] alnv=localAlignments.get();
			final Alignment best=alnv[0], current=alnv[1];
			best.clear();
			current.clear();
			
			int alignments=0;
//			System.err.println("list.size()="+list.size());
			final int a0=r.a, b0=r.b;
			final int mmLimit=Tools.max(8, fetchMM+3);//+3 was tested as optimal at maxrefmm=5.
			best.mismatches=mmLimit;
			final int breakMismatches=0;//Tools.mid(2, maxRefMismatches-1, 0);
			for(SeqPosM sq : list) {
				alignments++;
				final byte[] seq;
				final int count, pos;
//				synchronized(sq) {synchronized(sq.seq()) {
					seq=sq.seq();
					count=sq.count; pos=sq.pos();
//				}}
				final int a2=pos, b2=pos+seq.length-1;
				int matches=0, mismatches=0;
				for(int i=Tools.max(-a2, 0), j=Tools.max(a2, 0), lim=Tools.min(seq.length, bases.length-j); 
						i<lim && mismatches<mmLimit; i++, j++) {
					final byte x=seq[i], y=bases[j];//Branchless version
					final int m=((x==y) & (symbolToNumber[x]>=0) ? 1 : 0);
					matches+=m;
					mismatches+=(m^1);
				}
				current.set(seq, pos, a0, b0, matches, mismatches, count);
//				final int score=matches-mismatches*3;
//				final int shift=Tools.absdif(a0, a2)+Tools.absdif(b0, b2);
//				boolean valid=(matches>=minRefMatches && mismatches<=maxRefMismatches);
				if(current.compareTo(best)>0) {
					assert(current.valid);
					best.setFrom(current);
					
					if(best.shift==0 && best.mismatches<=breakMismatches) {break;}
					if(best.mismatches==0 && best.a2>=0 && best.b2<bases.length && 
							(sortRefCandidates || 
									(best.a2<=a0 && best.b2>=b0 && (best.a2==a0 || best.b2==b0)))) {break;}
				}
			}
			cTrackerT.refMismatchList.increment(best.mismatches);
			cTrackerT.alignments+=alignments;
			
			if(!best.valid) {return null;}
			cTrackerT.refMismatchListValid.increment(best.mismatches);
//			cTrackerT.validAlignments++;
			
			if(incrementRef) {
				AtomicLong al=refUseSet.get(new SeqCount(best.seq.clone()));
				al.incrementAndGet();
			}
			
			assert(best.count>=minRefCount || best.count==1) : best.count+", "+minRefCount;
			cTrackerT.alignedToRef++;
			return new Alignment().setFrom(best);
		}

		//TODO: Add bounds constraints to prevent violation of minspacer/maxspacer
		public Alignment doubleAlignRefToRange(byte[] bases, Range left, Range right, int repeatMismatches) {
			cTrackerT.alignmentRequested++;
			final ArrayList<SeqPosM> list;
			int fetchMM=Tools.max(maxRefMismatches, Math.round(bases.length*maxRefMismatchFraction));
			if(doubleFetch && repeatMismatches>0 && left.length()==right.length()) {
				list=refMapT.doubleFetch(bases, left.a, left.b, right.a, right.b, minRefOverlap, fetchMM, 
						maxRefTrim, maxRefSkew, minRefOverlapFractionQ, sortRefCandidates);
			}else {
				list=refMapT.fetch(bases, left.a, left.b, minRefOverlap, fetchMM, 
						maxRefTrim, maxRefSkew, minRefOverlapFractionQ, sortRefCandidates);
			}
			if(list==null) {return null;}
			
			Alignment[] alnv=localAlignments.get();
			final Alignment best=alnv[0], current=alnv[1];
			best.clear();
			current.clear();
			
			int alignments=0;
//			System.err.println("list.size()="+list.size());
			final int a0=left.a, b0=left.b, a1=right.a, b1=right.b;
			final int period=a1-a0;
			final int mmLimit=2*Tools.max(7, fetchMM+2);
			final int breakMismatches=0;//Tools.mid(2, maxRefMismatches-1, 0);
			for(SeqPosM sq : list) {
				alignments++;
				final byte[] seq;
				final int count, pos;
				seq=sq.seq();
				count=sq.count; pos=sq.pos();
				final int a3=pos, b3=pos+seq.length-1;
				final int a4=a3+period, b4=b3+period;
				int matches1=0, mismatches1=0;
				int matches2=0, mismatches2=0;
				int mismatches3=0;//Case where all 3 differ; probably indicates spacer sequence
				
				if(a3>=0 && b4<bases.length) {//Middle version; common case
					for(int i=0, j=a3, k=a4; i<seq.length && (mismatches1+mismatches2)<mmLimit; i++, j++, k++) {
						final byte x=seq[i], y=bases[j], z=bases[k];
						final int m=((x==y) & (symbolToNumber[x]>=0) ? 1 : 0), mm=m^1;
						matches1+=m;
						mismatches1+=mm;
						final int n=((x==z) & (symbolToNumber[x]>=0) ? 1 : 0), nn=n^1;
						matches2+=n;
						mismatches2+=nn;
						mismatches3+=(mm&nn&(y==z ? 0 : 1));
					}
					
//					//SIMD mode.  Barely faster, slightly different answer.
//					int smax=seq.length-1;
//					matches1=Vector.countMatches(seq, bases, 0, smax, a3, b3);
//					matches2=Vector.countMatches(seq, bases, 0, smax, a4, b4);
//					mismatches1=seq.length-matches1;
//					mismatches2=seq.length-matches2;
////					mismatches3=seq.length-Vector.countMatches(bases, bases, a3, b3, a4, b4);//Could use just where the two differ for mismatches3
//					assert(matches1>=0 && matches1<=seq.length) : matches1+", "+seq.length+", "+pos+", "+a3+", "+b3;
//					assert(matches2>=0 && matches2<=seq.length);
//					assert(mismatches1>=0 && mismatches1<=seq.length);
//					assert(mismatches2>=0 && mismatches2<=seq.length);
					
				}else {//Tip version; rare
					for(int i=0, j=a3, k=a4; i<seq.length; i++, j++, k++) {
						final byte x=seq[i], y, z;
						final int m, mm, n, nn;
						if(j>=0 && j<bases.length) {
							y=bases[j];
							m=((x==y) & (symbolToNumber[x]>=0) ? 1 : 0);
							mm=m^1;
							matches1+=m;
							mismatches1+=mm;
						}else{y='?'; m=0; mm=0;}
						if(k>=0 && k<bases.length) {
							z=bases[k];
							n=((x==z) & (symbolToNumber[x]>=0) ? 1 : 0);
							nn=n^1;
							matches2+=n;
							mismatches2+=nn;
						}else{z='?'; n=0; nn=0;}
						mismatches3+=(mm&nn&(y==z ? 0 : 1));
					}
				}
				final int score1=matches1-mismatches1*3;
				final int score2=matches2-mismatches2*3;
				final int len1=matches1+mismatches1, len2=matches2+mismatches2;
				final int matches=(len1>len2 ? matches1 : len2>len1 ? matches2 : (score1>=score2) ? matches1 : matches2);
				final int mismatches=(len1>len2 ? mismatches1 : len2>len1 ? mismatches2 : (score1>=score2) ? mismatches1 : mismatches2);
				current.set(seq, pos, a0, b0, matches, mismatches, count);
				current.score=score1+score2-3*mismatches3;
				
				if(current.compareTo(best)>0) {
					assert(current.valid);
					best.setFrom(current);
					if(best.perfect) {break;}
					
					if(best.shift==0 && Tools.max(mismatches1, mismatches2)<=breakMismatches) {break;}
					if((best.mismatches==0) && best.a2>=0 && best.b2<bases.length && 
							(sortRefCandidates || (best.a2<=a0 && best.b2>=b0 && 
							(best.a2==a0 || best.b2==b0)))) {break;}
				}
			}
			cTrackerT.refMismatchList.increment(best.mismatches);
			cTrackerT.alignments+=alignments;
			
			if(!best.valid) {return null;}
			cTrackerT.refMismatchListValid.increment(best.mismatches);
			
			if(incrementRef) {
				AtomicLong al=refUseSet.get(new SeqCount(best.seq.clone()));
				al.incrementAndGet();
			}
			assert(best.count>=minRefCount || best.count==1) : best.count+", "+minRefCount;
			cTrackerT.alignedToRef++;
			return new Alignment().setFrom(best);
		}
		
		public void consensus(byte[] bases, ArrayList<Crispr> list) {
			if(list==null || list.size()<2) {return;}
//			System.err.println("Running consensus on "+list+"\n"+new String(bases));
			int delta=1;
			int changed=0;
//			System.err.println("ba");
//			System.err.println(list);
			for(int tries=0; delta>0 && tries<2; tries++) {
				//TODO: This went into an infinite loop once with:
				//AAAAGCGGCTCCATTGAAGCTCGAAAAGCTTACCCGTTATCTCCGTTCGA
				//CCCCACATTTTCCGCTTTGAAAAAAGCGGCTCCATTGAAGCAAACTTTTT
				//GCTTGACGGCGCGTATAGGATGGTGTAATTTTCCGCTTTGAAAAAAAGCG
				//GCTCCATTGAAGCACGAACGGGGCGACTGTTCCAGCGTCCGTGCCTTCCATTTTCCGCT
				//and command line:
				//java -ea -Xmx2g jgi.CrisprFinder in=merged.fq.gz consensus=t 
				//outc=foundCrisprs.fa outr=foundRepeats.fa outs=foundSpacers.fa 
				//ow int=f t=8 rmmpad=3 rmm=4 rmmt=1 k=15 mm=1 minrepeattail=9 grow=t 
				//lookahead=5 minrepeats=2 minOverlapConsensus=18 ref=Repeats.fna t=8 kref=17 
				//minrepeat=20 minrepeat0=19 reads=2 t=1 in=x.fq
				//...specifically due to setting minrepeat0=19 while using the ref.
				//It looks like the problem is that they span the whole sequence though, so each end
				//kept getting adjusted back to the way it had been.
//				System.err.println(list);
				delta=0;
				Crispr prev=null;
				int disordered=0;
				for(int i=0; i<list.size(); i++) {
					Crispr p=list.get(i);
//					System.err.println(p.toString(bases));
					if(prev==null) {prev=p;}
					else {
//						System.err.println(prev.toString(bases));
						int d=consensus(prev, p, bases);
						delta+=d;
						changed+=(d==0 ? 0 : 1);
//						System.err.println("delta="+delta);
						disordered+=(prev.compareTo(p)>0 ? 1 : 0);
					}
				}
//				System.err.println("bb");
				
				if(delta==0) {break;}
				if(disordered>0) {
					Collections.sort(list);
					disordered=0;
				}
				delta=0;
				prev=null;
				for(int i=list.size()-1; i>=0; i--) {
					Crispr p=list.get(i);
//					System.err.println(p.toString(bases));
					if(prev==null) {prev=p;}
					else {
//						System.err.println(prev.toString(bases));
						int d=consensus(p, prev, bases);
						delta+=d;
						changed+=(d==0 ? 0 : 1);
//						System.err.println("delta="+delta);
						disordered+=(prev.compareTo(p)>0 ? 1 : 0);
					}
				}
				if(disordered>0) {
					Collections.sort(list);
					disordered=0;
				}
//				System.err.println("bc");
			}
			if(changed>0) {
				int removed=0;
				for(int i=0; i<list.size(); i++) {
//					System.err.println("bd");
					Crispr p=list.get(i);
					if(p.trimmedConsensus>0) {cTrackerT.trimmedByConsensus++;}
					if(p.matches<0) {
						assert(p.matches==-1);
						countMatches(p, bases);

						int repeat=p.maxLength();
						int gap=p.gap();
						if(gap<minSpacer || gap>maxSpacer || repeat<minRepeat || repeat>maxRepeat
								|| p.minLength()<minTailRepeat) {
							removed++;
							list.set(i, null);
						}
					}
					if(i>0) {
						Crispr a=list.get(i-1);
						Crispr b=list.get(i);
					}
				}
//				System.err.println("be");
				if(removed>0) {
					Tools.condenseStrict(list);
				}
			}
//			System.err.println("bf");
		}
		
		int removeShort(byte[] bases, ArrayList<Crispr> list) {
			int removed=0;
			for(int i=0; i<list.size(); i++) {
				Crispr p=list.get(i);
//				if(p.matches<0) {
//					assert(p.matches==-1);
//					countMatches(p, bases);
//				}

				int repeat=p.maxLength();
				int gap=p.gap();
				if(gap<minSpacer || gap>maxSpacer || repeat<minRepeat || repeat>maxRepeat
						|| p.minLength()<minTailRepeat) {
					removed++;
					list.set(i, null);
				}
			}
			if(removed>0) {
				Tools.condenseStrict(list);
			}
			return removed;
		}
		
		public int consensus(Crispr left, Crispr right, byte[] bases) {
//			System.err.println("Performing consensus on:\n"+left+"\n"+right);
			if(left.b.equals(right.a)){return 0;}//Perfect!
			final int max=bases.length-1;
			final int overlap0=left.b.overlap(right.a);
//			System.err.println("Overlap="+overlap0);
			assert(left.a.a<=right.a.a) : "\n"+left+"\n"+right+"\n"+new String(bases)+"\n";
			
			if(overlap0<minOverlapConsensus && overlap0<minRepeat-2) {return 0;}
			if(left.a.a<=0 && right.b.b>=max) {return 0;}//Unlikely situation of a single repeat pair spanning the whole read
			if(left.minLength()<minTailRepeat || right.minLength()<minTailRepeat) {return 0;}
			
			assert(left.a.length()>=0) : new String(bases);
			assert(left.b.length()>=0) : new String(bases);
			assert(right.a.length()>=0) : new String(bases);
			assert(right.b.length()>=0) : new String(bases);
			
//			System.err.println("\ninitial: left="+toString(left,bases)+"\nright="+toString(right,bases));
			int trimmedLeft=0, trimmedRight=0, extended=0;
			
			{//New concise version
				
				//If not touching sides and they don't overlap enough, skip
				if(left.a.a>0 && right.b.b<max && 
						(overlap0<minOverlapConsensus || overlap0<minRepeat)) {return 0;}
				
				if(left.a.a<=0) {//If touching the left
					int dif=left.b.a-right.a.a;
					if(dif>0) {//Extend left side of left repeat
						left.b.a-=dif;
						extended+=dif;
						left.extendedConsensus+=dif;
					}
//					System.err.println("a: left="+toString(left,bases)+"\nright="+toString(right,bases));
				}else {
					int dif=Tools.min(right.minLength(), right.b.a-right.a.a);
					if(dif>0 && dif<=maxTrimConsensus) {
//						while(left.b.a>right.a.a) {//Shrink left side of right repeat
//							right.b.a++;
//							right.a.a++;
//							trimmedRight++;
//						}
						right.b.a+=dif;
						right.a.a+=dif;
						trimmedRight+=dif;
					}
//					System.err.println("b: left="+toString(left,bases)+"\nright="+toString(right,bases));
				}
				
				if(right.b.b>=max) {//touches right side
					int dif=left.b.b-right.a.b;
					if(dif>0) {//Extend right side of right repeat
						right.a.b+=dif;
						extended+=dif;
						right.extendedConsensus+=dif;
					}
//					System.err.println("c: left="+toString(left,bases)+"\nright="+toString(right,bases));
				}else {
					int dif=Tools.min(left.minLength(), left.b.b-right.a.b);
					if(dif>0 && dif<=maxTrimConsensus) {//Shrink right side of left repeat
//						System.err.println();
//						System.err.println(left.toString(bases)+"\n"+right.toString(bases));
//						while(left.b.b>right.a.b) {
//							left.b.b--;
//							left.a.b--;
//							trimmedLeft++;
//						}
						left.a.b-=dif;
						left.b.b-=dif;
						trimmedLeft+=dif;
//						if(trimmedLeft>0) {
//							System.err.println("trimmed="+trimmedLeft);
//							System.err.println(left.toString(bases)+"\n"+right.toString(bases));
//						}
					}
//					System.err.println("d: left="+toString(left,bases)+"\nright="+toString(right,bases));
				}

				if(Tools.absdif(left.b.a, right.a.a)<=maxTrimConsensus) {
					while(left.b.a<right.a.a && left.b.a<=left.b.b) {//Shrink left side of left repeat
						if(left.sameLength()) {left.a.a++;}
						left.b.a++;
						trimmedLeft++;
					}
				}
//				System.err.println("e: left="+toString(left,bases)+"\nright="+toString(right,bases));
				if(Tools.absdif(left.b.b, right.a.b)<=maxTrimConsensus) {
					while(left.b.b<right.a.b && right.a.b>=right.a.a) {//Shrink right side of right repeat
						if(right.sameLength()) {right.b.b--;}
						right.a.b--;
						trimmedRight++;
					}
				}
//				System.err.println("f: left="+toString(left,bases)+"\nright="+toString(right,bases));
			}
			assert(left.a.length()>=0) : new String(bases)+"\n"+left.toString(bases)+"\n"+right.toString(bases);
			assert(left.b.length()>=0) : new String(bases)+"\n"+left.toString(bases)+"\n"+right.toString(bases);
			assert(right.a.length()>=0) : new String(bases)+"\n"+left.toString(bases)+"\n"+right.toString(bases);
			assert(right.b.length()>=0) : new String(bases)+"\n"+left.toString(bases)+"\n"+right.toString(bases);
			
			left.trimmedConsensus+=trimmedLeft;
			right.trimmedConsensus+=trimmedRight;
			int ret=trimmedLeft+trimmedRight+extended;
			if(ret>0) {
				//countMatches(left, bases);
				//countMatches(right, bases);
				left.matches=right.matches=-1;
			}
			return ret;
		}
		
		//TODO: Brute force should probably start and stop a couple bp in from the repeat edges, 
		//or else have an extra few bp mismatch allowance to try trimming from the ends...  since the repeat
		//it is searching for may have extended into the spacer
		/** Looks for repeats that did not get found from kmer seeds. */
		int bruteForce(byte[] bases, ArrayList<Crispr> list) {
			if(list==null || list.isEmpty()) {return 0;}
			int al=0, ar=0, am=0;
			
			if(bruteForceLeft) {
				for(Crispr rp=list.get(0); rp!=null; ) {
					rp=bruteForceLeft(bases, rp);
					if(rp==null) {break;}
					list.add(rp);//Puts the new one at the end, thus requiring sorting
//					if(rp.a.a<=0 && rp.a.length()<rp.b.length()) {cTrackerT.partialTipRepeats++;}
					al++;
				}
				if(al>0) {Collections.sort(list);}
			}

			if(bruteForceRight) {
				for(Crispr rp=list.get(list.size()-1); rp!=null; ) {
					rp=bruteForceRight(bases, rp);
					if(rp==null) {break;}
					list.add(rp);//Puts the new one at the end
//					if(rp.b.b>=bases.length-1 && rp.a.length()>rp.b.length()) {cTrackerT.partialTipRepeats++;}
					ar++;
				}
				if(ar>0) {Collections.sort(list);}//Should not be needed...
			}
			
			if(list.size()<2) {return 0;}
			
			//TODO
			if(bruteForceMiddle) {
				Crispr prev=null;
				for(int i=0, max=list.size(); i<max; i++) {
					Crispr p=list.get(i);
					if(prev!=null) {
						Crispr mid=bruteForceMiddle(bases, prev, p, list);
						if(mid!=null) {
							list.add(i+1, mid);//Can be slow on long sequences
							p=mid;
							i++;
							am++;
						}
					}
					prev=p;
				}
				if(am>0) {Collections.sort(list);}//Probably not necessary
			}
			return al+ar+am;
		}
		
		/** Looks for a repeat to the left of the leftmost detected repeat. */
		Crispr bruteForceLeft(byte[] bases, Crispr rp0) {
			final Range r=rp0.a;
			if(r.a<minTailPeriod || r.b>=bases.length-1) {return null;}
			final int a0=r.a, b0=r.b;
			
			final int limit=Tools.max(0, minTailRepeat-1, a0-maxSpacer-1);
			final int mrmm=maxRepeatMismatches+mrmmPad;

//			System.err.println("\n"+new String(bases)+"\n"+rp0.toString(bases));
//			System.err.println("A: "+limit+", "+mrmm);
			
			for(int i=a0-minSpacer; i>=limit; i--) {
				int mm=alignLeft(bases, a0, b0, i, mrmm);
//				System.err.println("i="+i+", mm="+mm);
				if(mm<=mrmm && mm>=0) {
					final int b1=i;
					final int a1=Tools.max(0, b1-r.length()+1);
					final boolean edge=a1<=0 || b1>=bases.length-1;
					final int rlen=b1-a1+1;
					final int extramm=mm-maxRepeatMismatches;
//					System.err.println("B: a1="+a1+", b1="+b1+", edge="+edge+", extramm="+extramm);
					if(rlen-extramm>=minRepeat || (rlen-extramm>=minTailRepeat && edge)) {
						Crispr rp=new Crispr(a1, b1, a0, b0);
						int gap=rp.gap(), repeat=rp.a.length();

						assert(rp.b.length()>=rp.a.length());
						assert(rp.b.length()==r.length());
						assert(a1==0 || rp.a.length()==r.length());
						assert(a1==0 || rp.sameLength());
						assert(/*gap>=minSpacer &&*/ gap<=maxSpacer && repeat>=minTailRepeat && repeat<=maxRepeat);
						shrink(bases, rp);
						gap=rp.gap();
						repeat=rp.a.length();
						countMatches(rp, bases);
						assert(rp.mismatches<=mrmm);
//						assert(rp.mismatches!=0) : rp.toString(bases);
						if(gap>=minSpacer && gap<=maxSpacer && 
								(repeat>=minRepeat || (repeat>=minTailRepeat && edge)) && 
								repeat<=maxRepeat && (rp.mismatches<=maxRepeatMismatchesTail ||
								(repeat>=minRepeat && !edge && rp.mismatches<=maxRepeatMismatches))) {
							return rp;
						}else {
//							assert(false) : rp.toString(bases)+"\n"+ "gap="+gap+", repeat="+repeat+
//							", minRepeat="+minRepeat+", minTailRepeat="+minTailRepeat+
//							", minSpacer="+minSpacer+", maxSpacer="+maxSpacer+
//							", matches="+rp.matches+", mm="+rp.mismatches;
						}
					}
				}
			}
			return null;
		}
		
		/** Looks for a repeat to the right of the rightmost detected repeat. */
		Crispr bruteForceRight(byte[] bases, Crispr rp0) {
			final Range r=rp0.b;
			if(r.b+minTailPeriod>=bases.length || r.a<=0) {return null;}
			final int a0=r.a, b0=r.b;
			
			final int limit=Tools.min(bases.length-1, bases.length-minTailRepeat-1, b0+maxSpacer-1);
			final int mrmm=maxRepeatMismatches+mrmmPad;
			for(int i=b0+minSpacer; i<=limit; i++) {
//				System.err.println("alignRight("+a0+", "+b0+", "+i+"), limit="+limit+", minTailRepeat="+minTailRepeat+", rlen="+r.length()+", blen="+bases.length);
				int mm=alignRight(bases, a0, b0, i, mrmm);
				if(mm<=mrmm && mm>=0) {
					final int a1=i;
					final int b1=Tools.min(bases.length-1, a1+r.length()-1);
					final boolean edge=a1<=0 || b1>=bases.length-1;
					final int rlen=b1-a1+1;
					final int extramm=mm-maxRepeatMismatches;
					if(rlen-extramm>=minRepeat || (rlen-extramm>=minTailRepeat && edge)) {
						Crispr rp=new Crispr(a0, b0, a1, b1);
						int gap=rp.gap(), repeat=rp.a.length();
						assert(rp.b.length()<=rp.a.length());
						assert(rp.a.length()==r.length());
						assert(b1==bases.length-1 || rp.b.length()==r.length());
						assert(b1==bases.length-1 || rp.sameLength());
						assert(/*gap>=minSpacer &&*/ gap<=maxSpacer && repeat>=minTailRepeat && repeat<=maxRepeat) : +mm+", "+mrmm+", "+extramm+", "+rlen+", "+gap+", "+repeat+", "+minSpacer;
						shrink(bases, rp);
						gap=rp.gap();
						repeat=rp.b.length();
						countMatches(rp, bases);
						assert(rp.mismatches<=mrmm);
						if(gap>=minSpacer && gap<=maxSpacer && 
								(repeat>=minRepeat || (repeat>=minTailRepeat && edge)) && 
								repeat<=maxRepeat && (rp.mismatches<=maxRepeatMismatchesTail ||
								(repeat>=minRepeat && !edge && rp.mismatches<=maxRepeatMismatches))) {
							return rp;
						}
					}
				}
			}
			return null;
		}
		
		/** Looks for repeats between repeats. */
		Crispr bruteForceMiddle(byte[] bases, Crispr left, Crispr right, ArrayList<Crispr> list) {
			//These are not really important but I want to make sure left is to the left of right.
			//Probably they'll be violated sometimes and should be fixed or prevented.
			//...after testing, it looks like these are short tandem repeats or interleaved repeats.  They do happen.
//			assert(left.a.a<right.a.a) : list+"\n"+left.toString(bases)+"\n"+right.toString(bases)+"\n"+new String(bases)+"\n";
//			assert(left.a.b<right.a.b) : list+"\n"+left.toString(bases)+"\n"+right.toString(bases)+"\n"+new String(bases)+"\n";
//			assert(left.b.a<right.b.a) : list+"\n"+left.toString(bases)+"\n"+right.toString(bases)+"\n"+new String(bases)+"\n";
//			assert(left.b.b<right.b.b) : list+"\n"+left.toString(bases)+"\n"+right.toString(bases)+"\n"+new String(bases)+"\n";
			
			final int space=right.a.a-left.b.b-1;
			if(space<minSpacer-2) {return null;}//Too short
			if(space>2*maxSpacer+maxRepeat) {return null;}//Too long; TODO: could use max of left and right rlen
			
			if(space>=2*minSpacer+minRepeat) {
				
				//Look in the middle
				return null;
			}
			
			if(space<=maxSpacer) {
				//See if the repeats match
				return null;
			}
			return null;
		}
		
		int alignLeft(byte[] bases, int a2, int b2, int b1, int maxmm) {
			int matches=0, mismatches=0;
			for(int i=b1, j=b2; i>=0 && j>=a2 && mismatches<=maxmm; i--, j--) {
				int x=symbolToNumber[bases[i]], y=symbolToNumber[bases[j]];
				matches+=(x==y && x>=0 ? 1 : 0);
				mismatches+=(x==y && x>=0 ? 0 : 1);
			}
			if(mismatches>maxmm) {return -mismatches;}
			assert(matches+mismatches>=minTailRepeat) : matches+", "+mismatches+", "+minTailRepeat+", "+b1;
			return mismatches;
		}
		
		int alignRight(byte[] bases, int a1, int b1, int a2, int maxmm) {
//			System.err.println("alignRight("+new String(Arrays.copyOfRange(bases, a1, b1+1))+", "+a1+", "+b1+", "+a2+"); blen="+bases.length);
			int matches=0, mismatches=0;
			for(int i=a1, j=a2; i<=b1 && j<bases.length && mismatches<=maxmm; i++, j++) {
				int x=symbolToNumber[bases[i]], y=symbolToNumber[bases[j]];
				matches+=(x==y && x>=0 ? 1 : 0);
				mismatches+=(x==y && x>=0 ? 0 : 1);
			}
			if(mismatches>maxmm) {return -mismatches;}
			assert(matches+mismatches>=minTailRepeat) : matches+", "+mismatches+", "+minTailRepeat+", "+a2;
			return mismatches;
		}
		
//		/** Looks for a repeat to the left of the leftmost detected repeat. */
//		int bruteForceMiddle(byte[] bases, RangePair left, RangePair right, ArrayList<RangePair> list) {
//			//These are not really important but I want to make sure left is to the left of right.
//			//Probably they'll be violated sometimes and should be fixed or prevented.
//			assert(left.a.a<right.a.a);
//			assert(left.a.b<right.a.b);
//			assert(left.b.a<right.b.a);
//			assert(left.b.b<right.b.b);
//		}
		int countMismatches(Crispr p, byte[] bases) {
			countMatches(p, bases);
			return p.mismatches;
		}
		
		int countMatches(Crispr p, byte[] bases) {
//			System.err.println(p.toString(bases));
//			System.err.println("countMatches("+p.toString(bases)+"); blen="+bases.length);
			int matches=0, mismatches=0;
			if(p.internal(bases.length) && p.sameLength()) {
				for(int i=p.a.a, j=p.b.a; j<=p.b.b; i++, j++) {
					final byte x=bases[i], y=bases[j];//Branchless version
					final int m=((x==y) & (symbolToNumber[x]>=0) ? 1 : 0);
					matches+=m;
					mismatches+=(m^1);
				}
			}else if(p.a.a<1 && p.b.b>=bases.length-1) {//Spans range
				int m1=countMatchesLeft(p, bases);
				int m2=countMatchesRight(p, bases);
				matches=Tools.max(m1, m2);
				mismatches=Tools.min(p.a.length(), p.b.length())-matches;
				//Hopefully this will not cause problems with assertion errors.
			}else if(p.a.a<1) {
				
				//This should probably be prevented since their alignment is unknown.
				//However, it is hard to prevent and this assertion fired once during consensus.
//				assert(p.b.b<bases.length-1 || p.sameLength()) : "Pair spans entire sequence: "+
//						p.toString(bases)+"\n"+new String(bases)+"\n";
				
				for(int i=p.a.b, j=p.b.b; i>=0 && i>=p.a.a && j>=p.b.a; i--, j--) {
//					byte a=bases[i], b=bases[j];
//					if(a==b && AminoAcid.isFullyDefined(a)){matches++;}
//					else {mismatches++;}
					final byte x=bases[i], y=bases[j];//Branchless version
					final int m=((x==y) & (symbolToNumber[x]>=0) ? 1 : 0);
					matches+=m;
					mismatches+=(m^1);
				}
			}else {
				for(int i=p.a.a, j=p.b.a; j<bases.length && j<=p.b.b && i<=p.b.a; i++, j++) {
//					byte a=bases[i], b=bases[j];
//					if(a==b && AminoAcid.isFullyDefined(a)){matches++;}
//					else {mismatches++;}
					final byte x=bases[i], y=bases[j];//Branchless version
					final int m=((x==y) & (symbolToNumber[x]>=0) ? 1 : 0);
					matches+=m;
					mismatches+=(m^1);
				}
			}
			assert(mismatches<=maxRepeatMismatches+mrmmPad || ref!=null) : mismatches+", "+maxRepeatMismatches+", "+p.toString(bases)+"\n"+bases.length+": "+new String(bases);
			
			//This assertion fired for a couple tip inserts where the short side was 'N' or 'TN'
			//assert(matches>0 || p.minLength()<=1) : matches+", "+mismatches+", "+bases.length+", "+p.toString(bases)+"\n";
			p.matches=matches;
			p.mismatches=mismatches;
			return matches;
		}
		
		public final int countMatchesLeft(Crispr p, byte[] bases) {
			int matches=0;
			for(int i=p.a.b, j=p.b.b; i>=0 && i>=p.a.a && j>=p.b.a; i--, j--) {
				final byte x=bases[i], y=bases[j];//Branchless version
				final int m=((x==y) & (symbolToNumber[x]>=0) ? 1 : 0);
				matches+=m;
			}
			return matches;
		}
		
		public final int countMatchesRight(Crispr p, byte[] bases) {
			int matches=0;
			for(int i=p.a.a, j=p.b.a; j<bases.length && j<=p.b.b && i<=p.b.a; i++, j++) {
				final byte x=bases[i], y=bases[j];//Branchless version
				final int m=((x==y) & (symbolToNumber[x]>=0) ? 1 : 0);
				matches+=m;
			}
			return matches;
		}
		
		void clearMap() {
			kmerPositionMap.clear();
			allKmers.clear();
//			uniqueKmers.clear();
//			repeatKmers.clear();
		}

		private LongLongListHashMap kmerPositionMap;
		private LongList allKmers;
		private int[] acgtn;
		private PalindromeFinder palFinder;
		private CrisprTracker cTrackerT;
		private SeqMap refMapT;
		private EntropyTracker eTracker;
		private CellNet net;
		private float[] vec;

		private HashMap<SeqCount,SeqCountM> repeatMapT;
		private HashMap<SeqCount,SeqCountM> spacerMapT;
		
		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		/** Number of reads merged by this thread */
		protected long readsMergedT=0;
		
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
		/** Shared reject output stream */
		private final ConcurrentReadOutputStream rosu;
		/** Crispr output stream */
		private final ConcurrentReadOutputStream rosCrispr;
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	
	class Alignment implements Comparable<Alignment>{
		
		public Alignment() {
			clear();
		}
		
		void set(byte[] s, int p, int a0, int b0, int m, int mm, int c) {
			seq=s;
			pos=p;
			matches=m;
			mismatches=mm;
			score=matches-mismatches*3;
			count=c;
			a2=pos;
			b2=pos+seq.length-1;
			shift=Tools.absdif(a0, a2)+Tools.absdif(b0, b2);
			final int maxMM=Tools.max(maxRefMismatches, Math.round(s.length*maxRefMismatchFraction));
			valid=(matches>=minRefMatches && mismatches<=maxMM);//TODO: test gap
			perfect=(shift==0 && mismatches==0);
		}
		
		Alignment setFrom(Alignment o) {
			seq=o.seq;
			pos=o.pos;
			matches=o.matches;
			mismatches=o.mismatches;
			score=o.score;
			count=o.count;
			a2=o.a2;
			b2=o.b2;
			shift=o.shift;
			valid=o.valid;
			perfect=o.perfect;
			return this;
		}
		
		void clear() {
			valid=perfect=false;
			count=pos=matches=score=-1;
			mismatches=shift=10;
		}
		
		public int compareTo(Alignment b) {
			if(!valid) {return -1;}
			if(b==null) {return 1;}
			if(perfect!=b.perfect) {return perfect ? 1 : -1;}
			if(valid!=b.valid) {return valid ? 1 : -1;}
			if(score!=b.score) {return score-b.score;}
			if(count!=b.count) {return count-b.count;}
			return b.shift-shift;
		}
//		boolean pick=valid && (bestSeq==null || score>bestScore || (mismatches==0 && shift==0) ||
//				(score==bestScore && (count>bestCount || shift<bestShift)));
//		
//		if(shift==0 && mismatches<=breakMismatches) {break;}
//		if(mismatches==0 && a2>=0 && b2<bases.length && 
//				(sortRefCandidates || (a2<=a0 && b2>=b0 && (a2==a0 || b2==b0)))) {break;}
		
		public String toString() {
			return "sc="+score+", co="+count+", va="+valid+
					", pe="+perfect+", sh="+shift+", m="+matches+", mm="+mismatches;
		}
		
		byte[] seq;
		int count;
		int pos;
		int matches;
		int mismatches;
		int a2, b2;
		int score;
		int shift;
		boolean valid;
		boolean perfect;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;

	/** Primary reject output file path */
	private String outu1=null;
	/** Secondary reject output file path */
	private String outu2=null;

	private String ref=null;
	private String outCrispr=null;
	private String outCrisprHist=null;
	private String outPalindromeHist=null;
	private String outRepeat=null;
	private String outSpacer=null;
	private String outRef=null;
	private String outRefAndRepeats=null;
	private boolean forceMaps=true;
	private int minCountToPrint=0;
	
	private boolean printPals=true;
	private boolean printScore=true;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;

	private int kRepeat=13;//ideal is probably around 13-15
	private int rqhdist=0;
	private int maskMiddle=0;//middle mask of 1 or even 2 is good
	private final long midmaskRepeat;//This is an OR mask not an AND mask
	private int minSpacer=14;
	private int maxSpacer=60;
	private int minRepeat0=11;
	private int minRepeat=22;
	private int maxRepeat=56;
	private int minTailRepeat=9;
	private final int minPeriod;
	private final int minCrispr;
	private final int minTailPeriod;
	private final int minTailCrispr;
	private float minRGC=0.09f;
	private float maxRGC=0.89f;
	private int maxRepeatMismatches=3;
	private int maxRepeatMismatchesTail=1;
	private int mrmmPad=3;
	private int crisprOutPad=0;
	private int minOverlapConsensus=18;
	private int maxTrimConsensus=5;
	
	private float minEntropy=0.84f;
	private int entropyWindow=80;
	private int entropyK=4;
	private int maxPoly=7;
	private int maxNs=0;
	
	private int kRef=13;
	private int maskMiddleRef=1;
	private int minRefMatches=18;
	private int maxRefMismatches=5;
	private float maxRefMismatchFraction=0.2f;
	private int minRefCount=0;
	private int minRefOverlap=7;
	private int maxRefSkew=8;
	private int maxRefTrim=10;
	private float minRefOverlapFractionQ=.00f;//Has no effect until about 0.6, but not very useful
	private float minRefOverlapFractionR=.00f;
	private boolean sortRefCandidates=false;
	private boolean autoSortRefCandidates=true;
	private boolean incrementRef=false;
	private HashMap<SeqCountM, AtomicLong> refUseSet;

	private boolean revertAlignmentFailures=false;
	private boolean shrinkAlignmentFailures=false;
	private boolean discardAlignmentFailures=true;
	private boolean discardUnaligned=false;
	private boolean doubleAlign=true;
	private boolean doubleFetch=true;
	
	private int minRepeats=2;
	private boolean masked=false;
	private boolean consensus=true;
	private final boolean bruteForce;
	private boolean bruteForceLeft=true, bruteForceRight=true, bruteForceMiddle=false;
	private boolean refinement=false; //TODO: subsumption by contained sequences
	
	private boolean grow=true;
	private int growLookahead=5;//5 seems better than 6.  4 is probably too low and 7 is too high.

	private int minPal=5;
	private int minPalMatches=4;
	private int maxPalMismatches=2;
	private int minLoop=3;//4 excludes 0.04%
	private int maxLoop=26;
	private int minTail=0;
	private int maxTail=24;
	private int maxTailDif=21;
	private boolean forceSymmetry=false;
	private boolean requirePalindrome=false;
	private boolean annotate=false;
	
	private final CrisprTracker cTracker=new CrisprTracker();
	private final PalindromeTracker pTracker=new PalindromeTracker();
	private final PalindromeTracker pTrackerFull=new PalindromeTracker();
	//Only use palindromes from full-length sequences in the histogram.
	private boolean fullLengthPalStats=true;
	private boolean printStats=true;
	private boolean printHist=false;

	private boolean useQueues=false;
	private final ArrayBlockingQueue<SeqCountM> repeatQueue=new ArrayBlockingQueue<SeqCountM>(256);
	private final ArrayBlockingQueue<SeqCountM> spacerQueue=new ArrayBlockingQueue<SeqCountM>(256);
	
	//At one time I changed these to HashMap to fix a concurrency issue due to mutable keys
	//Changing the key to a String also worked
	//Ultimately I just ensured that the key and value were independent objects
	private HashMap<SeqCount,SeqCountM> repeatMap;
	private HashMap<SeqCount,SeqCountM> spacerMap;
	private SeqMap refMap;
	private int maxRefCount=0;
	private long refMapSum=0;

	static final byte[] symbolToNumber=AminoAcid.baseToNumber;
	static final byte[] symbolToComplementNumber=AminoAcid.baseToComplementNumber;
	

	private final ThreadLocal<Alignment[]> localAlignments=new ThreadLocal<Alignment[]>(){
        @Override protected Alignment[] initialValue() {
        	return new Alignment[] {new Alignment(), new Alignment(), new Alignment()};
        }
    };
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	
	protected long repeatsFound; //TODO: And etc - crisprs out and so forth, and palindromes.

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	private boolean merge=false;
	private long readsMerged=0;
	
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
	
	/** Primary reject output file */
	private final FileFormat ffoutu1;
	/** Secondary reject output file */
	private final FileFormat ffoutu2;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	
	public static synchronized void loadNet() {loadNet(netFile);}
	private static synchronized void loadNet(String fname) {
		if(!useNet) {
			net0=null;
			return;
		}
		if(fname==null) {
			netFile=null;
			net0=null;
			return;
		}
		if(netFile!=null && netFile.equals(fname) && net0!=null){
			return;
		}
		netFile=fname;
		net0=CellNetParser.load(netFile);
		netCutoff=(setCutoff ? netCutoff : net0.cutoff);
	}
	
	private static final String defaultNetFile=Data.findPath("?crispr.bbnet.gz", false);
	private static String netFile=defaultNetFile;
	private static CellNet net0;
	private static boolean setCutoff=false;
	private static float netCutoff=-1;
	private static boolean useNet=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** Print more verbose messages */
	public static final boolean verbose2=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	/** Reads are output in input order */
	private boolean ordered=false;
	
}
