package bbduk;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicLongArray;

import aligner.SideChannel3;
import cardinality.CardinalityTracker;
import dna.AminoAcid;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import hiseq.FlowcellCoordinate;
import jgi.BBMerge;
import jgi.BBMergeOverlapper;
import jgi.CalcTrueQuality;
import json.JsonObject;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.SamLine;
import structures.IntList;
import structures.ListNum;
import structures.Quantizer;
import structures.StringCount;
import tracker.EntropyTracker;
import tracker.PolymerTracker;
import tracker.ReadStats;
import var2.AnalyzeVars;
import var2.ScafMap;
import var2.Var;
import var2.VarMap;

/**
 * Handles read-processing part of BBDuk
 * @author Brian Bushnell
 * @date November 20, 2025
 *
 */
public class BBDukProcessor {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public BBDukProcessor(BBDukParser p, BBDukIndexAndLoader index_, ConcurrentReadInputStream cris_, 
		ConcurrentReadOutputStream ros_, ConcurrentReadOutputStream rosb_, ConcurrentReadOutputStream ross_){
		
		parser=p;
		index=index_;
		cris=cris_;
		ros=ros_;
		rosb=rosb_;
		ross=ross_;
		
		silent=p.silent;
		json=p.json;
		swift=p.swift;
		khistIn=p.khistIn;
		khistOut=p.khistOut;
		ALLOW_LOCAL_ARRAYS=p.ALLOW_LOCAL_ARRAYS;
		ref=p.ref;
		literal=p.literal;
		in1=p.in1;
		in2=p.in2;
		outstats=p.outstats;
		outrqc=p.outrqc;
		outrpkm=p.outrpkm;
		outrefstats=p.outrefstats;
		polymerStatsFile=p.polymerStatsFile;
		tossJunk=p.tossJunk;
		maxBasesOutm=p.maxBasesOutm;
		maxBasesOutu=p.maxBasesOutu;
		useShortKmers=p.useShortKmers;
		maskMiddle=p.maskMiddle;
		midMaskLen=p.midMaskLen;
		qHammingDistance=p.qHammingDistance;
		qHammingDistance2=p.qHammingDistance2;
		trimPad=p.trimPad;
		xMinLoc=p.xMinLoc;
		yMinLoc=p.yMinLoc;
		xMaxLoc=p.xMaxLoc;
		yMaxLoc=p.yMaxLoc;
		locationFilter=p.locationFilter;
		varMap=p.varMap;
		scafMap=p.scafMap;
		fixVariants=p.fixVariants;
		unfixVariants=p.unfixVariants;
		filterVars=p.filterVars;
		maxBadSubs=p.maxBadSubs;
		maxBadSubAlleleDepth=p.maxBadSubAlleleDepth;
		minBadSubReadDepth=p.minBadSubReadDepth;
		minBadSubEDist=p.minBadSubEDist;
		maxBadAlleleFraction=p.maxBadAlleleFraction;
		entropyCutoff=p.entropyCutoff;
		entropyHighpass=p.entropyHighpass;
		entropyMark=p.entropyMark;
		entropyMask=p.entropyMask;
		entropyTrim=p.entropyTrim;
		entropyMaskLowercase=p.entropyMaskLowercase;
		calcEntropy=p.calcEntropy;
		jsonStats=p.jsonStats;
		countPolymers=p.countPolymers;
		polymerChar1=p.polymerChar1;
		polymerChar2=p.polymerChar2;
		polymerLength=p.polymerLength;
		skipR1=p.skipR1;
		skipR2=p.skipR2;
		ecc=p.ecc;
		makeReadStats=p.makeReadStats;
		rcomp=p.rcomp;
		forbidNs=p.forbidNs;
		middleMask=p.middleMask;
		k=p.k;
		k2=p.k2;
		kbig=p.kbig;
		keff=p.keff;
		mink=p.mink;
		maxBadKmers0=p.maxBadKmers0;
		minKmerFraction=p.minKmerFraction;
		minCoveredFraction=p.minCoveredFraction;
		recalibrateQuality=p.recalibrateQuality;
		quantizeQuality=p.quantizeQuality;
		qtrimLeft=p.qtrimLeft;
		qtrimRight=p.qtrimRight;
		trimClip=p.trimClip;
		trimPolyA=p.trimPolyA;
		trimPolyGLeft=p.trimPolyGLeft;
		trimPolyGRight=p.trimPolyGRight;
		filterPolyG=p.filterPolyG;
		maxNonPoly=p.maxNonPoly;
		trimPolyCLeft=p.trimPolyCLeft;
		trimPolyCRight=p.trimPolyCRight;
		filterPolyC=p.filterPolyC;
		trimq=p.trimq;
		trimE=p.trimE;
		minAvgQuality=p.minAvgQuality;
		minBaseQuality=p.minBaseQuality;
		minAvgQualityBases=p.minAvgQualityBases;
		chastityFilter=p.chastityFilter;
		failBadBarcodes=p.failBadBarcodes;
		removeBadBarcodes=p.removeBadBarcodes;
		failIfNoBarcode=p.failIfNoBarcode;
		barcodes=p.barcodes;
		maxNs=p.maxNs;
		minConsecutiveBases=p.minConsecutiveBases;
		minBaseFrequency=p.minBaseFrequency;
		minReadLength=p.minReadLength;
		maxReadLength=p.maxReadLength;
		minLenFraction=p.minLenFraction;
		kfilter=p.kfilter;
		ktrimLeft=p.ktrimLeft;
		ktrimRight=p.ktrimRight;
		ktrimN=p.ktrimN;
		ktrimExclusive=p.ktrimExclusive;
		ksplit=p.ksplit;
		trimSymbol=p.trimSymbol;
		kmaskLowercase=p.kmaskLowercase;
		kmaskFullyCovered=p.kmaskFullyCovered;
		addTrimmedToBad=p.addTrimmedToBad;
		findBestMatch=p.findBestMatch;
		trimPairsEvenly=p.trimPairsEvenly;
		forceTrimLeft=p.forceTrimLeft;
		forceTrimRight=p.forceTrimRight;
		forceTrimRight2=p.forceTrimRight2;
		forceTrimModulo=p.forceTrimModulo;
		minGC=p.minGC;
		maxGC=p.maxGC;
		filterGC=p.filterGC;
		usePairGC=p.usePairGC;
		restrictLeft=p.restrictLeft;
		restrictRight=p.restrictRight;
		removePairsIfEitherBad=p.removePairsIfEitherBad;
		trimFailuresTo1bp=p.trimFailuresTo1bp;
		printNonZeroOnly=p.printNonZeroOnly;
		rename=p.rename;
		speed=p.speed;
		noAccel=p.noAccel;
		accel=p.accel;
		pairedToSingle=p.pairedToSingle;
		amino=p.amino;
		bitsPerBase=p.bitsPerBase;
		minlen=p.minlen;
		minminlen=p.minminlen;
		minlen2=p.minlen2;
		shift2=p.shift2;
		mask=p.mask;
		kmask=p.kmask;
		lengthMasks=p.lengthMasks;
		symbolToNumber=p.symbolToNumber;
		symbolToNumber0=p.symbolToNumber0;
		symbolToComplementNumber0=p.symbolToComplementNumber0;
		trimByOverlap=p.trimByOverlap;
		useQualityForOverlap=p.useQualityForOverlap;
//		strictOverlap=p.strictOverlap;
		minOverlap0=p.minOverlap0;
		minOverlap=p.minOverlap;
		minInsert0=p.minInsert0;
		minInsert=p.minInsert;
		maxRatio=p.maxRatio;
		ratioMargin=p.ratioMargin;
		ratioOffset=p.ratioOffset;
		efilterRatio=p.efilterRatio;
		efilterOffset=p.efilterOffset;
		pfilterRatio=p.pfilterRatio;
		meeFilter=p.meeFilter;
		histogramsBeforeProcessing=p.histogramsBeforeProcessing;

		MAKE_IHIST=p.MAKE_IHIST;
		outstream=BBDukParser.outstream;
		overwrite=BBDukParser.overwrite;
		append=BBDukParser.append;
		showSpeed=BBDukParser.showSpeed;
		THREADS=BBDukParser.workers;
		STATS_COLUMNS=BBDukParser.STATS_COLUMNS;
		REPLICATE_AMBIGUOUS=BBDukParser.REPLICATE_AMBIGUOUS;
		
		loglogIn=(p.loglog ? CardinalityTracker.makeTracker(p.parser) : null);
		loglogOut=(p.loglogOut ? CardinalityTracker.makeTracker(p.parser) : null);
		
		//Initialize polymer-tracking
		if(countPolymers){
			pTracker=new PolymerTracker();
		}else{
			pTracker=null;
		}
		
		{
			//2. Sync State - References and Arrays
			sidechannel=index.sidechannel;
			ref=index.ref;
			scaffoldNames=index.scaffoldNames;
			scaffoldLengths=index.scaffoldLengths;
			scaffoldReadCounts=index.scaffoldReadCounts; 
			scaffoldBaseCounts=index.scaffoldBaseCounts;
			scafMap=index.scafMap;
			varMap=index.varMap;

			//3. Sync Primitives and Status Flags
			ALLOW_LOCAL_ARRAYS=index.ALLOW_LOCAL_ARRAYS;
			storedKmers=index.storedKmers; // Critical for logic checks
			fixVariants=index.fixVariants; // Critical for variant logic
			errorState|=index.errorState;  // Merge error states
		}
		
		readstats=makeReadStats ? new ReadStats() : null;
		
		final int alen=(scaffoldNames==null ? 0 : scaffoldNames.size());
		
		if(findBestMatch){
			countArray=new int[alen];
			idList=new IntList();
			countList=new IntList();
		}else{
			countArray=null;
			idList=countList=null;
		}
		
		overlapVector=(trimByOverlap ? new int[5] : null);
		
		if(ALLOW_LOCAL_ARRAYS && alen>0 && alen<10000 && scaffoldReadCounts!=null && scaffoldBaseCounts!=null){
			scaffoldReadCountsT=new long[alen];
			scaffoldBaseCountsT=new long[alen];
		}else{
			scaffoldReadCountsT=scaffoldBaseCountsT=null;
		}
		
		if(calcEntropy){
			eTracker=new EntropyTracker(amino, Tools.max(0, entropyCutoff), entropyHighpass);
		}else{
			eTracker=null;
		}
		
		maxBasesOutm=(maxBasesOutm>0 ? Tools.max(1, maxBasesOutm/THREADS) : -1);
		maxBasesOutu=(maxBasesOutu>0 ? Tools.max(1, maxBasesOutu/THREADS) : -1);

		flowCoords=(locationFilter ? new FlowcellCoordinate() : null);
	}
	
	public BBDukProcessor(BBDukProcessor proc) {
		this(proc.parser, proc.index, proc.cris, proc.ros, proc.rosb, proc.ross);
	}
	
	@Override
	public BBDukProcessor clone() {
		return new BBDukProcessor(this);
	}
	
	public void add(BBDukProcessor pt) {
		readsIn+=pt.readsIn;
		basesIn+=pt.basesIn;
		readsOut+=pt.readsOutu;
		basesOut+=pt.basesOutu;
		readsKFiltered+=pt.readsKFiltered;
		basesKFiltered+=pt.basesKFiltered;
		readsQTrimmed+=pt.readsQTrimmed;
		basesQTrimmed+=pt.basesQTrimmed;
		readsFTrimmed+=pt.readsFTrimmed;
		basesFTrimmed+=pt.basesFTrimmed;
		readsKTrimmed+=pt.readsKTrimmed;
		basesKTrimmed+=pt.basesKTrimmed;
		readsTrimmedBySwift+=pt.readsTrimmedBySwift;
		basesTrimmedBySwift+=pt.basesTrimmedBySwift;
		readsTrimmedByOverlap+=pt.readsTrimmedByOverlap;
		basesTrimmedByOverlap+=pt.basesTrimmedByOverlap;
		badGcReads+=pt.badGcReads;
		badGcBases+=pt.badGcBases;
		badHeaderReads+=pt.badHeaderReads;
		badHeaderBases+=pt.badHeaderBases;
		readsQFiltered+=pt.readsQFiltered;
		basesQFiltered+=pt.basesQFiltered;
		readsNFiltered+=pt.readsNFiltered;
		basesNFiltered+=pt.basesNFiltered;
		readsEFiltered+=pt.readsEFiltered;
		basesEFiltered+=pt.basesEFiltered;
		readsPolyTrimmed+=pt.readsPolyTrimmed;
		basesPolyTrimmed+=pt.basesPolyTrimmed;
		
		if(pTracker!=null){
			pTracker.add(pt.pTracker);
		}
		if(pt.scaffoldReadCountsT!=null && scaffoldReadCounts!=null){
			for(int i=0; i<pt.scaffoldReadCountsT.length; i++){scaffoldReadCounts.addAndGet(i, pt.scaffoldReadCountsT[i]);}
//			pt.scaffoldReadCountsT=null;//TODO: final
		}
		if(pt.scaffoldBaseCountsT!=null && scaffoldBaseCounts!=null){
			for(int i=0; i<pt.scaffoldBaseCountsT.length; i++){scaffoldBaseCounts.addAndGet(i, pt.scaffoldBaseCountsT[i]);}
//			pt.scaffoldBaseCountsT=null;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void printOutput(long startTime){
		
		/* Write statistics to files */
		writeStats();
		writeRPKM();
		index.writeRefStats(in1, in2, readsIn);
		writeRqc();
		if(pTracker!=null && polymerStatsFile!=null){
			ReadWrite.writeString(pTracker.toHistogramCumulative(), polymerStatsFile);
		}
		
		/* Unload kmers and sequence data to save memory */
		if(RELEASE_TABLES){
			index.cleanup();
		}
		
		if(silent){return;}
		if(json){
			outstream.println(toJson(startTime));
			return;
		}
		
		outstream.println("\nInput:                  \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		
		if((ref!=null || literal!=null) && !(ktrimLeft || ktrimRight || ktrimN)){
			outstream.println("Contaminants:           \t"+readsKFiltered+" reads ("+toPercent(readsKFiltered, readsIn)+") \t"+
					basesKFiltered+" bases ("+toPercent(basesKFiltered, basesIn)+")");
			outstream.flush();
		}
		if(qtrimLeft || qtrimRight){
			outstream.println("QTrimmed:               \t"+readsQTrimmed+" reads ("+toPercent(readsQTrimmed, readsIn)+") \t"+
					basesQTrimmed+" bases ("+toPercent(basesQTrimmed, basesIn)+")");
		}
		if(trimPolyA>0 || trimPolyGLeft>0 || trimPolyGRight>0 || filterPolyG>0 || trimPolyCLeft>0 || trimPolyCRight>0 || filterPolyC>0){
			outstream.println("Polymer-trimmed:        \t"+readsPolyTrimmed+" reads ("+toPercent(readsPolyTrimmed, readsIn)+") \t"+
					basesPolyTrimmed+" bases ("+toPercent(basesPolyTrimmed, basesIn)+")");
		}
		if(forceTrimLeft>0 || forceTrimRight>0 || forceTrimRight2>0 || forceTrimModulo>0){
			outstream.println("FTrimmed:               \t"+readsFTrimmed+" reads ("+toPercent(readsFTrimmed, readsIn)+") \t"+
					basesFTrimmed+" bases ("+toPercent(basesFTrimmed, basesIn)+")");
		}
		if(ktrimLeft || ktrimRight || ktrimN){
			String x=(ktrimN ? "KMasked: " : "KTrimmed:");
			outstream.println(x+"               \t"+readsKTrimmed+" reads ("+toPercent(readsKTrimmed, readsIn)+") \t"+
					basesKTrimmed+" bases ("+toPercent(basesKTrimmed, basesIn)+")");
		}
		if(swift){
			outstream.println("Trimmed by Swift:       \t"+readsTrimmedBySwift+" reads ("+toPercent(readsTrimmedBySwift, readsIn)+") \t"+
					basesTrimmedBySwift+" bases ("+toPercent(basesTrimmedBySwift, basesIn)+")");
		}
		if(trimByOverlap){
			outstream.println("Trimmed by overlap:     \t"+readsTrimmedByOverlap+" reads ("+toPercent(readsTrimmedByOverlap, readsIn)+") \t"+
					basesTrimmedByOverlap+" bases ("+toPercent(basesTrimmedByOverlap, basesIn)+")");
		}
		if(filterGC){
			outstream.println("Filtered by GC:         \t"+badGcReads+" reads ("+toPercent(badGcReads, readsIn)+") \t"+
					badGcBases+" bases ("+toPercent(badGcBases, basesIn)+")");
		}
		if(locationFilter || chastityFilter || removeBadBarcodes){
			outstream.println("Filtered by header:     \t"+badHeaderReads+" reads ("+toPercent(badHeaderReads, readsIn)+") \t"+
					badHeaderBases+" bases ("+toPercent(badHeaderBases, basesIn)+")");
		}
		if(minAvgQuality>0 || minBaseQuality>0 || maxNs>=0 || minBaseFrequency>0 || chastityFilter || removeBadBarcodes){
			outstream.println("Low quality discards:   \t"+readsQFiltered+" reads ("+toPercent(readsQFiltered, readsIn)+") \t"+
					basesQFiltered+" bases ("+toPercent(basesQFiltered, basesIn)+")");
		}
		if(polymerChar1>=0 && polymerChar2>=0){
			outstream.println("Polymer Counts:         \t"+padRight(pTracker.getCountCumulative(polymerChar1, polymerLength)+" "+Character.toString((char)polymerChar1), 18)+"\t"+
					padRight(pTracker.getCountCumulative(polymerChar2, polymerLength)+" "+Character.toString((char)polymerChar2), 18)+"\t"+
					"("+Tools.format("%.4f", pTracker.calcRatioCumulative(polymerChar1, polymerChar2, polymerLength))+" ratio)");
		}
		if(sidechannel!=null){
			outstream.println(sidechannel.stats(readsIn, basesIn));
		}
		if(calcEntropy){
			String prefix;
			if(entropyTrim>0){
				prefix=("Entropy-trimmed:        \t");
			}else if(entropyMask){
				prefix=("Entropy-masked:         \t");
			}else{
				prefix=("Low entropy discards:   \t");
			}
			outstream.println(prefix+readsEFiltered+" reads ("+toPercent(readsEFiltered, readsIn)+") \t"+
					basesEFiltered+" bases ("+toPercent(basesEFiltered, basesIn)+")");
		}

		final long readsRemoved=readsIn-readsOut;
		final long basesRemoved=basesIn-basesOut;
		
		outstream.println("Total Removed:          \t"+readsRemoved+" reads ("+toPercent(readsRemoved, readsIn)+") \t"+
				basesRemoved+" bases ("+toPercent(basesRemoved, basesIn)+")");
		
		outstream.println("Result:                 \t"+readsOut+" reads ("+toPercent(readsOut, readsIn)+") \t"+
				basesOut+" bases ("+toPercent(basesOut, basesIn)+")");
		
		if(loglogIn!=null){
			outstream.println("Unique "+loglogIn.k+"-mers:         \t"+loglogIn.cardinality());
			if(khistIn!=null){
				loglogIn.printKhist(khistIn, overwrite, append, true, 2);
			}
		}
		if(loglogOut!=null){
			outstream.println("Unique "+loglogOut.k+"-mers out:     \t"+loglogOut.cardinality());
			if(khistOut!=null){
				loglogOut.printKhist(khistOut, overwrite, append, true, 2);
			}
		}
	}
	
	/**
	 * Formats ratio as percentage string with two decimal places.
	 * @param numerator Numerator value
	 * @param denominator Denominator value
	 * @return Percentage string (e.g., "45.32%")
	 */
	public static String toPercent(long numerator, long denominator){
		if(denominator<1){return "0.00%";}
		return Tools.format("%.2f%%",numerator*100.0/denominator);
	}
	
	/**
	 * Right-pads string with spaces to minimum length.
	 * @param s Input string
	 * @param minLen Minimum desired length
	 * @return Padded string
	 */
	private static String padRight(String s, int minLen){
		while(s.length()<minLen){s=s+" ";}
		return s;
	}
	
	/**
	 * Formats processing statistics as JSON string.
	 * @param startTime Processing start time in nanoseconds
	 * @return JSON-formatted statistics string
	 */
	String toJson(long startTime){

		jsonStats.add("k", k);
		jsonStats.add("mode", ktrimLeft ? "ktrimLeft" : ktrimRight ? "ktrimRight" : ktrimN ? "ktrimN" : "kFilter");
		jsonStats.add("readsIn", readsIn);
		jsonStats.add("basesIn", basesIn);
		
		if((ref!=null || literal!=null) && !(ktrimLeft || ktrimRight || ktrimN)){
			jsonStats.add("readsKFiltered", readsKFiltered);
			jsonStats.add("basesKFiltered", basesKFiltered);
		}
		if(qtrimLeft || qtrimRight){
			jsonStats.add("readsQTrimmed", readsQTrimmed);
			jsonStats.add("basesQTrimmed", basesQTrimmed);
		}
		if(trimPolyA>0 || trimPolyGLeft>0 || trimPolyGRight>0 || filterPolyG>0 || trimPolyCLeft>0 || trimPolyCRight>0 || filterPolyC>0){
			jsonStats.add("readsPolyTrimmed", readsPolyTrimmed);
			jsonStats.add("basesPolyTrimmed", basesPolyTrimmed);
		}
		if(forceTrimLeft>0 || forceTrimRight>0 || forceTrimRight2>0 || forceTrimModulo>0){
			jsonStats.add("readsFTrimmed", readsFTrimmed);
			jsonStats.add("basesFTrimmed", basesFTrimmed);
		}
		if(ktrimLeft || ktrimRight || ktrimN){
			String x=(ktrimN ? "KMasked: " : "KTrimmed:");
			jsonStats.add("reads"+x, readsKTrimmed);
			jsonStats.add("bases+x", basesKTrimmed);
		}
		if(swift){
			jsonStats.add("readsTrimmedBySwift", readsTrimmedBySwift);
			jsonStats.add("basesTrimmedBySwift", basesTrimmedBySwift);
		}
		if(trimByOverlap){
			jsonStats.add("readsTrimmedByOverlap", readsTrimmedByOverlap);
			jsonStats.add("basesTrimmedByOverlap", basesTrimmedByOverlap);
		}
		if(filterGC){
			jsonStats.add("badGcReads", badGcReads);
			jsonStats.add("badGcBases", badGcBases);
		}
		if(locationFilter || chastityFilter || removeBadBarcodes){
			jsonStats.add("badHeaderReads", badHeaderReads);
			jsonStats.add("badHeaderBases", badHeaderBases);
		}
		if(minAvgQuality>0 || minBaseQuality>0 || maxNs>=0 || minBaseFrequency>0 || chastityFilter || removeBadBarcodes){
			jsonStats.add("readsQFiltered", readsQFiltered);
			jsonStats.add("basesQFiltered", basesQFiltered);
		}
		if(polymerChar1>=0 && polymerChar2>=0){
			jsonStats.add("poly"+Character.toString((char)polymerChar1), pTracker.getCountCumulative(polymerChar1, polymerLength));
			jsonStats.add("poly"+Character.toString((char)polymerChar2), pTracker.getCountCumulative(polymerChar2, polymerLength));
			jsonStats.add("polyRatio", pTracker.calcRatioCumulative(polymerChar1, polymerChar2, polymerLength));
		}
		if(calcEntropy){
			String suffix;
			if(entropyTrim>0){
				suffix=("EntropyTrimmed");
			}else if(entropyMask){
				suffix=("EntropyMasked");
			}else{
				suffix=("EntropyFiltered");
			}
			jsonStats.add("reads"+suffix, readsEFiltered);
			jsonStats.add("bases"+suffix, basesEFiltered);
		}

		final long readsRemoved=readsIn-readsOut;
		final long basesRemoved=basesIn-basesOut;

		jsonStats.add("readsRemoved", readsRemoved);
		jsonStats.add("basesRemoved", basesRemoved);
		jsonStats.add("readsOut", readsOut);
		jsonStats.add("basesOut", basesOut);
		
		if(loglogIn!=null){
			jsonStats.add("uniqueKmersIn", loglogIn.cardinality());
		}
		if(loglogOut!=null){
			jsonStats.add("uniqueKmersOut", loglogOut.cardinality());
		}
		jsonStats.add("time", (System.nanoTime()-startTime)/1000000000.0);
		return jsonStats.toString();
	}
	
	/**
	 * Write statistics about how many reads matched each reference scaffold.
	 */
	void writeStats(){
		if(outstats==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outstats, overwrite, false, false);
		tsw.start();
		
		long rsum=0, bsum=0;
		
		/* Create StringCount list of scaffold names and hitcounts */
		ArrayList<StringCount> list=new ArrayList<StringCount>();
		for(int i=1; i<scaffoldNames.size(); i++){
			final long num1=scaffoldReadCounts.get(i), num2=scaffoldBaseCounts.get(i);
			if(num1>0 || !printNonZeroOnly){
				rsum+=num1;
				bsum+=num2;
				final String s=scaffoldNames.get(i);
				final int len=scaffoldLengths.get(i);
				final StringCount sn=new StringCount(s, len, num1, num2);
				list.add(sn);
			}
		}
		Shared.sort(list);
		final double rmult=100.0/(readsIn>0 ? readsIn : 1);
		final double bmult=100.0/(basesIn>0 ? basesIn : 1);
		
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		
		if(STATS_COLUMNS==3){
			tsw.print(Tools.format("#Total\t%d\n",readsIn));
			tsw.print(Tools.format("#Matched\t%d\t%.5f%%\n",rsum,rmult*rsum));
			tsw.print("#Name\tReads\tReadsPct\n");
			for(int i=0; i<list.size(); i++){
				StringCount sn=list.get(i);
				tsw.print(Tools.format("%s\t%d\t%.5f%%\n",sn.name,sn.reads,(sn.reads*rmult)));
			}
		}else{
			tsw.print(Tools.format("#Total\t%d\t%d\n",readsIn,basesIn));
			tsw.print(Tools.format("#Matched\t%d\t%.5f%%\n",rsum,rmult*rsum,bsum,bsum*bmult));
			tsw.print("#Name\tReads\tReadsPct\tBases\tBasesPct\n");
			for(int i=0; i<list.size(); i++){
				StringCount sn=list.get(i);
				tsw.print(Tools.format("%s\t%d\t%.5f%%\t%d\t%.5f%%\n",sn.name,sn.reads,(sn.reads*rmult),sn.bases,(sn.bases*bmult)));
			}
		}
		
		tsw.poisonAndWait();
	}
	
	/**
	 * Write RPKM statistics.
	 */
	void writeRPKM(){
		if(outrpkm==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outrpkm, overwrite, false, false);
		tsw.start();

		/* Count mapped reads */
		long mapped=0;
		for(int i=0; i<scaffoldReadCounts.length(); i++){
			mapped+=scaffoldReadCounts.get(i);
		}
		
		/* Print header */
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		tsw.print(Tools.format("#Reads\t%d\n",readsIn));
		tsw.print(Tools.format("#Mapped\t%d\n",mapped));
		tsw.print(Tools.format("#RefSequences\t%d\n",Tools.max(0, scaffoldNames.size()-1)));
		tsw.print("#Name\tLength\tBases\tCoverage\tReads\tRPKM\n");
		
		final float mult=1000000000f/Tools.max(1, mapped);
		
		/* Print data */
		for(int i=1; i<scaffoldNames.size(); i++){
			final long reads=scaffoldReadCounts.get(i);
			final long bases=scaffoldBaseCounts.get(i);
			final String s=scaffoldNames.get(i);
			final int len=scaffoldLengths.get(i);
			final double invlen=1.0/Tools.max(1, len);
			final double mult2=mult*invlen;
			if(reads>0 || !printNonZeroOnly){
				tsw.print(Tools.format("%s\t%d\t%d\t%.4f\t%d\t%.4f\n",s,len,bases,bases*invlen,reads,reads*mult2));
			}
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Write RQCFilter stats.
	 * @param time Elapsed time, nanoseconds
	 */
	void writeRqc(){
		if(outrqc==null){return;}
		addToRqcMap();
		if(outrqc.endsWith("hashmap")){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outrqc, overwrite, false, false);
		tsw.start();
		tsw.println(rqcString());
		tsw.poisonAndWait();
	}
	
	/** Formats RQC statistics map as string output.
	 * @return Formatted RQC statistics string */
	public static String rqcString(){
		if(RQC_MAP==null){return null;}
		StringBuilder sb=new StringBuilder();
		
		String[] keys=new String[] {"inputReads", "inputBases", "qtrimmedReads", "qtrimmedBases", 
			"qfilteredReads", "qfilteredBases", "ktrimmedReads", "ktrimmedBases", "kfilteredReads",
			"kfilteredBases", "outputReads", "outputBases"};
		
		for(String key : keys){
			Object value=RQC_MAP.get(key);
			if(value!=null){
				sb.append(key+"="+value+"\n");
			}
		}
		
		return sb.toString();
	}
	
	/** Populates RQC statistics map with processing counts.
	 * Adds input, filtered, trimmed, and output read/base counts. */
	void addToRqcMap(){
		putRqc("inputReads", readsIn, false, false);
		putRqc("inputBases", basesIn, false, false);
		if(qtrimLeft || qtrimRight){
			putRqc("qtrimmedReads", readsQTrimmed, false, true);
			putRqc("qtrimmedBases", basesQTrimmed, false, true);
		}
		putRqc("qfilteredReads", readsQFiltered, false, true);
		putRqc("qfilteredBases", basesQFiltered, false, true);
		
		if(ktrimLeft || ktrimRight || ktrimN){
			putRqc("ktrimmedReads", readsKTrimmed, true, true);
			putRqc("ktrimmedBases", basesKTrimmed, true, true);
		}else{
			putRqc("kfilteredReads", readsKFiltered, false, true);
			putRqc("kfilteredBases", basesKFiltered, false, true);
		}
		putRqc("outputReads", readsOut, true, false);
		putRqc("outputBases", basesOut, true, false);
	}
	
	/**
	 * Adds or updates entry in RQC statistics map.
	 *
	 * @param key Statistics key name
	 * @param value Count value to store
	 * @param evict Whether to replace existing value
	 * @param add Whether to add to existing value
	 */
	public static void putRqc(String key, Long value, boolean evict, boolean add){
		if(RQC_MAP==null){RQC_MAP=new HashMap<String,Long>();}
		Long old=RQC_MAP.get(key);
		if(evict || old==null){RQC_MAP.put(key, value);}
		else if(add){RQC_MAP.put(key, value+old);}
	}
	
	/** Calculates ratio of two specified polymer bases.
	 * @return Ratio of polymer base counts (base1/base2) */
	double getPolymerRatio(){
		return pTracker.calcRatioCumulative(polymerChar1, polymerChar2, polymerLength);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/
		
	
	public void process(){
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		ArrayList<Read> bad=(rosb==null ? null : new ArrayList<Read>(Shared.bufferLen()));
		ArrayList<Read> single=new ArrayList<Read>(Shared.bufferLen());
		ArrayList<Read> sideList=(sidechannel==null ? null : new ArrayList<Read>(Shared.bufferLen()));

		final boolean ktrimLeftOrRight=ktrimLeft || ktrimRight;
		final boolean ktrimTips=(ktrimLeft && ktrimRight);
		final boolean doKmerTrimming=storedKmers>0 && (ktrimLeft || ktrimRight || ktrimN || ksplit);
		final boolean doKmerFiltering=storedKmers>0 && !doKmerTrimming;

		//While there are more reads lists...
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			int removed=0;

			//For each read (or pair) in the list...
			for(int i=0; i<reads.size(); i++){
				final Read r1=reads.get(i);
				if(pairedToSingle && r1.mate!=null){
					r1.mate=null;
				}
				final Read r2=r1.mate;

				if(!r1.validated()){r1.validate(true);}
				if(r2!=null && !r2.validated()){r2.validate(true);}

				boolean remove=false;
				if(tossJunk){
					if(r1!=null && r1.junk()){
						setDiscarded(r1);
						remove=true;
					}
					if(r2!=null && r2.junk()){
						setDiscarded(r2);
						remove=true;
					}
				}

				if(isNotDiscarded(r1)){

					if(histogramsBeforeProcessing){addToHistograms(r1, r2);}

					if(loglogIn!=null){loglogIn.hash(r1);}

				}

				final int initialLength1=r1.length();
				final int initialLength2=r1.mateLength();
				final int initialPairLength=initialLength1+initialLength2;
				final int pairCount=r1.pairCount();

				final int minlen1=(int)Tools.max(initialLength1*minLenFraction, minReadLength);
				final int minlen2=(int)Tools.max(initialLength2*minLenFraction, minReadLength);

				if(verbose){outstream.println("Considering read "+r1.id+" "+new String(r1.bases));}

				readsIn+=pairCount;
				basesIn+=initialPairLength;

				if(!remove){//due to being junk
					if(chastityFilter){
						if(r1!=null && r1.failsChastity()){
							setDiscarded(r1);
							if(r2!=null){setDiscarded(r2);}
						}
					}

					if(locationFilter && isNotDiscarded(r1)){
						flowCoords.setFrom(r1.id);
						boolean discard=false;
						if(xMinLoc>-1 && flowCoords.x<xMinLoc){discard=true;}
						if(xMaxLoc>-1 && flowCoords.x>xMaxLoc){discard=true;}
						if(yMinLoc>-1 && flowCoords.y<yMinLoc){discard=true;}
						if(yMaxLoc>-1 && flowCoords.y>yMaxLoc){discard=true;}
						if(discard){
							setDiscarded(r1);
							if(r2!=null){setDiscarded(r2);}
						}
					}

					if(removeBadBarcodes){
						if(isNotDiscarded(r1) && r1.failsBarcode(barcodes, failIfNoBarcode)){
							if(failBadBarcodes){KillSwitch.kill("Invalid barcode detected: "+r1.id+"\nThis can be disabled with the flag barcodefilter=f");}
							setDiscarded(r1);
							if(r2!=null){setDiscarded(r2);}
						}
					}

					if(isDiscarded(r1)){
						badHeaderBases+=initialPairLength;
						badHeaderReads+=pairCount;
					}

					if(recalibrateQuality){
						if(isNotDiscarded(r1)){
							CalcTrueQuality.recalibrate(r1);
						}
						if(isNotDiscarded(r2)){
							CalcTrueQuality.recalibrate(r2);
						}
					}

					if(filterGC && (initialLength1>0 || initialLength2>0)){
						float gc1=(initialLength1>0 ? r1.gc() : -1);
						float gc2=(initialLength2>0 ? r2.gc() : gc1);
						if(gc1==-1){gc1=gc2;}
						if(usePairGC){
							final float gc;
							if(r2==null){
								gc=gc1;
							}else{
								gc=(gc1*initialLength1+gc2*initialLength2)/(initialPairLength);
							}
							gc1=gc2=gc;
						}
						if(isNotDiscarded(r1) && (gc1<minGC || gc1>maxGC)){
							setDiscarded(r1);
							badGcBases+=initialLength1;
							badGcReads++;
						}
						if(isNotDiscarded(r2) && (gc2<minGC || gc2>maxGC)){
							setDiscarded(r2);
							badGcBases+=initialLength2;
							badGcReads++;
						}
					}

					if(forceTrimLeft>0 || forceTrimRight>0 || forceTrimRight2>0 || forceTrimModulo>0){
						if(isNotDiscarded(r1)){
							final int len=r1.length();
							final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
							final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
							final int b1=forceTrimRight>0 ? forceTrimRight : len;
							final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
							final int b=Tools.min(b0, b1, b2);
							final int x=TrimRead.trimToPosition(r1, a, b, 1);
							basesFTrimmed+=x;
							readsFTrimmed+=(x>0 ? 1 : 0);
							if(r1.length()<minlen1){setDiscarded(r1);}
						}
						if(isNotDiscarded(r2)){
							final int len=r2.length();
							final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
							final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
							final int b1=forceTrimRight>0 ? forceTrimRight : len;
							final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
							final int b=Tools.min(b0, b1, b2);
							final int x=TrimRead.trimToPosition(r2, a, b, 1);
							basesFTrimmed+=x;
							readsFTrimmed+=(x>0 ? 1 : 0);
							if(r2.length()<minlen2){setDiscarded(r2);}
						}
					}

					if(filterVars){
						if(isNotDiscarded(r1)){
							boolean b=passesVariantFilter(r1);
							if(!b){setDiscarded(r1);}
						}
						if(isNotDiscarded(r2)){
							boolean b=passesVariantFilter(r2);
							if(!b){setDiscarded(r2);}
						}
					}

					if(isNotDiscarded(r1) && r1.length()<minlen1){setDiscarded(r1);}
					if(isNotDiscarded(r2) && r2.length()<minlen2){setDiscarded(r2);}

					if(removePairsIfEitherBad){remove=isDiscarded(r1) || isDiscarded(r2);}
					else{remove=isDiscarded(r1) && isNullOrDiscarded(r2);}
				}

				if(remove){
					if(r1!=null){
						basesQFiltered+=initialLength1;
						readsQFiltered++;
					}
					if(r2!=null){
						basesQFiltered+=initialLength2;
						readsQFiltered++;
					}
					if(bad!=null){bad.add(r1);}
				}else{

					if(ecc && r1!=null && r2!=null){BBMerge.findOverlapStrict(r1, r2, true);}

					//Process kmers
					if(doKmerTrimming){

						int rlen1=0, rlen2=0;
						int xsum=0;
						int rktsum=0;

						if(ktrimTips){
							if(r1!=null){
								int x=ktrimTips(r1);
								xsum+=x;
								rktsum+=(x>0 ? 1 : 0);
								rlen1=r1.length();
								if(rlen1<minlen1){setDiscarded(r1);}
							}
							if(r2!=null){
								int x=ktrimTips(r2);
								xsum+=x;
								rktsum+=(x>0 ? 1 : 0);
								rlen2=r2.length();
								if(rlen2<minlen2){setDiscarded(r2);}
							}
						}else if(ktrimLeftOrRight){
							if(r1!=null){
								int x=ktrim(r1);
								xsum+=x;
								rktsum+=(x>0 ? 1 : 0);
								rlen1=r1.length();
								if(rlen1<minlen1){setDiscarded(r1);}
							}
							if(r2!=null){
								int x=ktrim(r2);
								xsum+=x;
								rktsum+=(x>0 ? 1 : 0);
								rlen2=r2.length();
								if(rlen2<minlen2){setDiscarded(r2);}
							}
						}else if(ktrimN){
							if(r1!=null){
								int x=kmask(r1);
								xsum+=x;
								rktsum+=(x>0 ? 1 : 0);
								rlen1=r1.length();
								if(rlen1<minlen1){setDiscarded(r1);}
							}
							if(r2!=null){
								int x=kmask(r2);
								xsum+=x;
								rktsum+=(x>0 ? 1 : 0);
								rlen2=r2.length();
								if(rlen2<minlen2){setDiscarded(r2);}
							}
						}else if(ksplit){
							assert(r2==null);
							if(r1!=null){
								int oldLen=r1.pairLength();
								boolean b=ksplit(r1);
								int trimmed=oldLen-r1.pairLength();
								xsum+=trimmed;
								rktsum+=(trimmed>0 ? 1 : 0);
								rlen1=r1.length();
							}
						}

						if(ksplit){
							remove=(r1.mate!=null);
							if(remove && addTrimmedToBad && bad!=null){bad.add(r1);}
						}else if(shouldRemove(r1, r2)){
							if(!ktrimN){
								xsum+=(rlen1+rlen2);
								rktsum=pairCount;
							}
							remove=true;
							if(addTrimmedToBad && bad!=null){bad.add(r1);}
						}else if(ktrimRight && trimPairsEvenly && xsum>0 && r2!=null && r1.length()!=r2.length()){
							int x;
							if(r1.length()>r2.length()){
								x=TrimRead.trimToPosition(r1, 0, r2.length()-1, 1);
							}else{
								x=TrimRead.trimToPosition(r2, 0, r1.length()-1, 1);
							}
							if(rktsum<2){rktsum++;}
							xsum+=x;
							assert(r1.length()==r2.length()) : r1.length()+", "+r2.length();
						}
						basesKTrimmed+=xsum;
						readsKTrimmed+=rktsum;

					}else if(doKmerFiltering){
						//Do kmer matching

						if(minCoveredFraction>0){
							if(isNotDiscarded(r1)){
								final int minCoveredBases=(int)Math.ceil(minCoveredFraction*r1.length());
								final int covered=countCoveredBases(r1, minCoveredBases);
								if(covered>=minCoveredBases){setDiscarded(r1);}
							}
							if(isNotDiscarded(r2)){
								final int minCoveredBases=(int)Math.ceil(minCoveredFraction*r2.length());
								final int covered=countCoveredBases(r2, minCoveredBases);
								if(covered>=minCoveredBases){setDiscarded(r2);}
							}
						}else{

							final int maxBadKmersR1, maxBadKmersR2;
							if(minKmerFraction==0){
								maxBadKmersR1=maxBadKmersR2=maxBadKmers0;
							}else{
								final int vk1=r1.numValidKmers(keff), vk2=(r2==null ? 0 : r2.numValidKmers(keff));
								maxBadKmersR1=Tools.max(maxBadKmers0, (int)((vk1-1)*minKmerFraction));
								maxBadKmersR2=Tools.max(maxBadKmers0, (int)((vk2-1)*minKmerFraction));
							}

							if(!findBestMatch){
								final int a=(kbig<=k ? countSetKmers(r1, maxBadKmersR1) : countSetKmersBig(r1, maxBadKmersR1));
								final int b=(kbig<=k ? countSetKmers(r2, maxBadKmersR2) : countSetKmersBig(r2, maxBadKmersR2));

								if(r1!=null && a>maxBadKmersR1){setDiscarded(r1);}
								if(r2!=null && b>maxBadKmersR2){setDiscarded(r2);}

							}else{
								final int a=findBestMatch(r1, maxBadKmersR1);
								final int b=findBestMatch(r2, maxBadKmersR2);

								if(r1!=null && a>0){setDiscarded(r1);}
								if(r2!=null && b>0){setDiscarded(r2);}
							}
						}

						if(shouldRemove(r1, r2)){
							remove=true;
							if(r1!=null){
								readsKFiltered++;
								basesKFiltered+=initialLength1;
							}
							if(r2!=null){
								readsKFiltered++;
								basesKFiltered+=initialLength2;
							}
							if(bad!=null){bad.add(r1);}
						}

					}
				}

				if(!remove && trimByOverlap && r2!=null && expectedErrors(r1, r2)<meeFilter){

					if(aprob==null || aprob.length<r1.length()){aprob=new float[r1.length()];}
					if(bprob==null || bprob.length<r2.length()){bprob=new float[r2.length()];}

					//Do overlap trimming
					r2.reverseComplementFast();
					int bestInsert=BBMergeOverlapper.mateByOverlapRatio(r1, r2, aprob, bprob, overlapVector, minOverlap0, minOverlap,
						minInsert0, minInsert, maxRatio, 0.12f, ratioMargin, ratioOffset, 0.95f, 0.95f, useQualityForOverlap);

					if(bestInsert<minInsert){bestInsert=-1;}
					boolean ambig=(overlapVector[4]==1);
					final int bestBad=overlapVector[2];

					if(bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null && useQualityForOverlap){
						if(efilterRatio>0 && bestInsert>0 && !ambig){
							float bestExpected=BBMergeOverlapper.expectedMismatches(r1, r2, bestInsert);
							if((bestExpected+efilterOffset)*efilterRatio<bestBad){ambig=true;}
						}
						if(pfilterRatio>0 && bestInsert>0 && !ambig){
							float probability=BBMergeOverlapper.probability(r1, r2, bestInsert);
							if(probability<pfilterRatio){bestInsert=-1;}
						}
						if(meeFilter>=0 && bestInsert>0 && !ambig){
							float expected=BBMergeOverlapper.expectedMismatches(r1, r2, bestInsert);
							if(expected>meeFilter){bestInsert=-1;}
						}
					}

					r2.reverseComplementFast();

					if(bestInsert>0 && !ambig){
						if(bestInsert<r1.length()){
							if(verbose){outstream.println("Overlap right trimming r1 to "+0+", "+(bestInsert-1));}
							int x=TrimRead.trimToPosition(r1, 0, bestInsert-1, 1);
							if(verbose){outstream.println("Trimmed "+x+" bases: "+new String(r1.bases));}
							readsTrimmedByOverlap++;
							basesTrimmedByOverlap+=x;
						}
						if(bestInsert<r2.length()){
							if(verbose){outstream.println("Overlap right trimming r2 to "+0+", "+(bestInsert-1));}
							int x=TrimRead.trimToPosition(r2, 0, bestInsert-1, 1);
							if(verbose){outstream.println("Trimmed "+x+" bases: "+new String(r2.bases));}
							readsTrimmedByOverlap++;
							basesTrimmedByOverlap+=x;
						}
					}
				}

				if(!remove && swift){
					//Do Swift trimming

					int rlen1=0, rlen2=0;
					if(r1!=null){
						int x=trimSwift(r1);
						basesTrimmedBySwift+=x;
						readsTrimmedBySwift+=(x>0 ? 1 : 0);
						rlen1=r1.length();
						if(rlen1<minlen1){setDiscarded(r1);}
					}
					if(r2!=null){
						int x=trimSwift(r2);
						basesTrimmedBySwift+=x;
						readsTrimmedBySwift+=(x>0 ? 1 : 0);
						rlen2=r2.length();
						if(rlen2<minlen2){setDiscarded(r2);}
					}

					//Discard reads if too short
					if(shouldRemove(r1, r2)){
						basesTrimmedBySwift+=r1.pairLength();
						remove=true;
						if(addTrimmedToBad && bad!=null){bad.add(r1);}
					}
				}

				if(!remove && trimPolyA>0){
					//Do poly-A trimming

					int rlen1=0, rlen2=0;
					if(r1!=null){
						int x=trimPolyA(r1, trimPolyA);
						basesPolyTrimmed+=x;
						readsPolyTrimmed+=(x>0 ? 1 : 0);
						rlen1=r1.length();
						if(rlen1<minlen1){setDiscarded(r1);}
					}
					if(r2!=null){
						int x=trimPolyA(r2, trimPolyA);
						basesPolyTrimmed+=x;
						readsPolyTrimmed+=(x>0 ? 1 : 0);
						rlen2=r2.length();
						if(rlen2<minlen2){setDiscarded(r2);}
					}

					//Discard reads if too short
					if(shouldRemove(r1, r2)){
						basesPolyTrimmed+=r1.pairLength();
						remove=true;
						if(addTrimmedToBad && bad!=null){bad.add(r1);}
					}
				}

				if(!remove && (trimPolyGLeft>0 || trimPolyGRight>0 || filterPolyG>0)){
					//Do poly-G trimming

					int rlen1=0, rlen2=0;
					if(r1!=null){
						if(filterPolyG>0 && detectPolyLeft(r1, filterPolyG, maxNonPoly, (byte)'G')>=filterPolyG) {
							setDiscarded(r1);
							readsPolyTrimmed++;
						}else if(trimPolyGLeft>0 || trimPolyGRight>0){
							int x=trimPoly(r1, trimPolyGLeft, trimPolyGRight, maxNonPoly, (byte)'G');
							basesPolyTrimmed+=x;
							readsPolyTrimmed+=(x>0 ? 1 : 0);
							rlen1=r1.length();
							if(rlen1<minlen1){setDiscarded(r1);}
						}
					}
					if(r2!=null){
						if(filterPolyG>0 && detectPolyLeft(r2, filterPolyG, maxNonPoly, (byte)'G')>=filterPolyG) {
							setDiscarded(r2);
							readsPolyTrimmed++;
						}else if(trimPolyGLeft>0 || trimPolyGRight>0){
							int x=trimPoly(r2, trimPolyGLeft, trimPolyGRight, maxNonPoly, (byte)'G');
							basesPolyTrimmed+=x;
							readsPolyTrimmed+=(x>0 ? 1 : 0);
							rlen2=r2.length();
							if(rlen2<minlen2){setDiscarded(r2);}
						}
					}

					//Discard reads if too short
					if(shouldRemove(r1, r2)){
						basesPolyTrimmed+=r1.pairLength();
						remove=true;
						if(addTrimmedToBad && bad!=null){bad.add(r1);}
					}
				}

				if(!remove && (trimPolyCLeft>0 || trimPolyCRight>0 || filterPolyC>0)){
					//Do poly-C trimming

					int rlen1=0, rlen2=0;
					if(r1!=null){
						if(filterPolyC>0 && detectPolyLeft(r1, filterPolyC, maxNonPoly, (byte)'C')>=filterPolyC) {
							setDiscarded(r1);
							readsPolyTrimmed++;
						}else if(trimPolyCLeft>0 || trimPolyCRight>0){
							int x=trimPoly(r1, trimPolyCLeft, trimPolyCRight, maxNonPoly, (byte)'C');
							basesPolyTrimmed+=x;
							readsPolyTrimmed+=(x>0 ? 1 : 0);
							rlen1=r1.length();
							if(rlen1<minlen1){setDiscarded(r1);}
						}
					}
					if(r2!=null){
						if(filterPolyC>0 && detectPolyLeft(r1, filterPolyC, maxNonPoly, (byte)'C')>=filterPolyC) {
							setDiscarded(r2);
							readsPolyTrimmed++;
						}else if(trimPolyCLeft>0 || trimPolyCRight>0){
							int x=trimPoly(r2, trimPolyCLeft, trimPolyCRight, maxNonPoly, (byte)'C');
							basesPolyTrimmed+=x;
							readsPolyTrimmed+=(x>0 ? 1 : 0);
							rlen2=r2.length();
							if(rlen2<minlen2){setDiscarded(r2);}
						}
					}

					//Discard reads if too short
					if(shouldRemove(r1, r2)){
						basesPolyTrimmed+=r1.pairLength();
						remove=true;
						if(addTrimmedToBad && bad!=null){bad.add(r1);}
					}
				}

				if(!remove && (entropyMask || entropyTrim>0)){
					//Mask/trim entropy
					if(isNotDiscarded(r1)){
						int masked=(entropyTrim>0 ? trimLowEntropy(r1, null, eTracker) : maskLowEntropy(r1, null, eTracker));
						basesEFiltered+=masked;
						readsEFiltered+=(masked>0 ? 1 : 0);
					}
					if(isNotDiscarded(r2)){
						int masked=(entropyTrim>0 ? trimLowEntropy(r2, null, eTracker) : maskLowEntropy(r2, null, eTracker));
						basesEFiltered+=masked;
						readsEFiltered+=(masked>0 ? 1 : 0);
					}
				}

				if(entropyMark){
					markLowEntropy(r1, eTracker);
					markLowEntropy(r2, eTracker);
				}

				if(!remove){
					//Do quality trimming

					if(qtrimLeft || qtrimRight || trimClip){
						if(r1!=null){
							int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, trimE, 1, trimClip);
							basesQTrimmed+=x;
							readsQTrimmed+=(x>0 ? 1 : 0);

							//								assert(false) : trimClip+", "+x;
						}
						if(r2!=null){
							int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, trimE, 1, trimClip);
							basesQTrimmed+=x;
							readsQTrimmed+=(x>0 ? 1 : 0);
						}
					}

					if(isNotDiscarded(r1)){
						int len=r1.length();
						if(len<minlen1 || len>maxReadLength){setDiscarded(r1);}
					}
					if(isNotDiscarded(r2)){
						int len=r2.length();
						if(len<minlen2 || len>maxReadLength){setDiscarded(r2);}
					}

					//Discard reads if too short
					if(shouldRemove(r1, r2)){
						basesQTrimmed+=r1.pairLength();
						remove=true;
						if(addTrimmedToBad && bad!=null){bad.add(r1);}
					}

				}

				if(!remove){
					//Do quality filtering

					//Determine whether to discard the reads based on average quality
					if(minAvgQuality>0){
						if(r1!=null && r1.quality!=null && r1.avgQuality(false, minAvgQualityBases)<minAvgQuality){setDiscarded(r1);}
						if(r2!=null && r2.quality!=null && r2.avgQuality(false, minAvgQualityBases)<minAvgQuality){setDiscarded(r2);}
					}
					//Determine whether to discard the reads based on lowest quality base
					if(minBaseQuality>0){
						if(r1!=null && r1.quality!=null && r1.minQuality()<minBaseQuality){setDiscarded(r1);}
						if(r2!=null && r2.quality!=null && r2.minQuality()<minBaseQuality){setDiscarded(r2);}
					}
					//Determine whether to discard the reads based on the presence of Ns
					if(maxNs>=0){
						if(r1!=null && r1.countUndefined()>maxNs){
							readsNFiltered++;
							basesNFiltered+=r1.length();
							setDiscarded(r1);
						}
						if(r2!=null && r2.countUndefined()>maxNs){
							readsNFiltered++;
							basesNFiltered+=r2.length();
							setDiscarded(r2);
						}
					}
					//Determine whether to discard the reads based on a lack of useful kmers
					if(minConsecutiveBases>0){
						if(isNotDiscarded(r1) && !r1.hasMinConsecutiveBases(minConsecutiveBases)){setDiscarded(r1);}
						if(isNotDiscarded(r2) && !r2.hasMinConsecutiveBases(minConsecutiveBases)){setDiscarded(r2);}
					}
					//Determine whether to discard the reads based on minimum base frequency
					if(minBaseFrequency>0){
						if(r1!=null && r1.minBaseCount()<minBaseFrequency*r1.length()){setDiscarded(r1);}
						if(r2!=null && r2.minBaseCount()<minBaseFrequency*r2.length()){setDiscarded(r2);}
					}

					//Discard reads if too short
					if(shouldRemove(r1, r2)){
						basesQFiltered+=r1.pairLength();
						readsQFiltered+=pairCount;
						remove=true;
						if(addTrimmedToBad && bad!=null){bad.add(r1);}
					}
				}

				if(!remove && quantizeQuality) {
					Quantizer.quantize(r1);
					Quantizer.quantize(r2);
				}

				if(!remove && calcEntropy && entropyCutoff>=0 && !entropyMask && entropyTrim<1){
					//Test entropy
					if(isNotDiscarded(r1) && !eTracker.passes(r1.bases, true)){setDiscarded(r1);}
					if(isNotDiscarded(r2) && !eTracker.passes(r2.bases, true)){setDiscarded(r2);}

					if(shouldRemove(r1, r2)){
						basesEFiltered+=r1.pairLength();
						readsEFiltered+=pairCount;
						remove=true;
						if(bad!=null){bad.add(r1);}
					}
				}

				if(!remove && !histogramsBeforeProcessing){
					addToHistograms(r1, r2);
				}

				if(!remove && sidechannel!=null) {
					boolean mapped=sidechannel.map(r1, r2);
					if(mapped) {
						assert(r1.mapped() || r2.mapped());
						sideList.add(r1);
					}
				}

				if(ross!=null){
					if(isNotDiscarded(r1) && isNullOrDiscarded(r2)){
						Read clone=r1.clone();
						clone.mate=null;
						single.add(clone);
					}else if(r2!=null && isDiscarded(r1) && isNotDiscarded(r2)){
						Read clone=r2.clone();
						clone.mate=null;
						single.add(clone);
					}
				}

				if(remove && !trimFailuresTo1bp){
					//Evict read
					removed++;
					if(r2!=null){removed++;}
					reads.set(i, null);

					readsOutm+=pairCount;
					basesOutm+=r1.pairLength();
				}else{
					if(loglogOut!=null){loglogOut.hash(r1);}
					readsOutu+=pairCount;
					basesOutu+=r1.pairLength();
				}
			}

			//Send matched list to matched output stream
			if(rosb!=null){
				rosb.add(bad, ln.id);
				bad.clear();
			}

			//Send unmatched list to unmatched output stream
			if(ros!=null){
				ros.add((removed>0 ? Tools.condenseNew(reads) : reads), ln.id); //Creates a new list if old one became empty, to prevent shutting down the cris.
			}

			if(ross!=null){
				ross.add(single, ln.id);
				single.clear();
			}

			if(sidechannel!=null) {
				sidechannel.writeToMapped(sideList, ln.id);
				sideList.clear();
			}

			if(maxBasesOutm>=0 && basesOutm>=maxBasesOutm){break;}
			if(maxBasesOutu>=0 && basesOutu>=maxBasesOutu){break;}

			//Fetch a new read list
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		finishedSuccessfully=true;
	}

	private void setDiscarded(Read r){
		if(trimFailuresTo1bp){
			if(r.length()>1){TrimRead.trimByAmount(r, 0, r.length()-1, 1, false);}
		}else{
			r.setDiscarded(true);
		}
	}

	private boolean isDiscarded(Read r){
		if(r==null){return false;}
		if(r.discarded()){return true;}
		return trimFailuresTo1bp && r.length()==1;
	}

	private boolean isNullOrDiscarded(Read r){
		if(r==null){return true;}
		if(r.discarded()){return true;}
		return trimFailuresTo1bp && r.length()==1;
	}

	private boolean isNotDiscarded(Read r){
		if(r==null){return false;}
		if(r.discarded()){return false;}
		return !(trimFailuresTo1bp && r.length()==1);
	}

	private boolean shouldRemove(Read r1, Read r2){
		return (removePairsIfEitherBad && (isDiscarded(r1) || isDiscarded(r2))) || 
			(isDiscarded(r1) && isNullOrDiscarded(r2));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private void addToHistograms(Read r1, Read r2) {

		if(pTracker!=null){
			pTracker.addPair(r1);
		}

		if(fixVariants){
			AnalyzeVars.fixVars(r1, varMap, scafMap);
			AnalyzeVars.fixVars(r2, varMap, scafMap);
		}

		if(readstats!=null){
			readstats.addToHistograms(r1);

			if(MAKE_IHIST){
				SamLine sl1=r1.samline;
				if(sl1!=null && !r1.secondary() && sl1.pairnum()==0){
					readstats.addToInsertHistogram(sl1);
				}
			}
		}

		if(fixVariants && unfixVariants){
			AnalyzeVars.unfixVars(r1);
			AnalyzeVars.unfixVars(r2);
		}
	}


	/**
	 * Counts the number of kmer hits for a read.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return Number of hits
	 */
	private final int countSetKmers(final Read r, final int maxBadKmers){
		if(r==null || r.length()<k || storedKmers<1){return 0;}
		if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){return 0;}
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		int found=0;
		int len=0;

		final int start=(restrictRight<1 ? 0 : Tools.max(0, bases.length-restrictRight));
		final int stop=(restrictLeft<1 ? bases.length : Tools.min(bases.length, restrictLeft));

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=start; i<stop; i++){
			byte b=bases[i];
			long x=symbolToNumber0[b];
			long x2=symbolToComplementNumber0[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(forbidNs && !isFullyDefined(b)){len=0; rkmer=0;}else{len++;}
			if(verbose){
				String fwd=new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k));
				String rev=AminoAcid.reverseComplementBases(fwd);
				String fwd2=kmerToString(kmer, Tools.min(len, k));
				outstream.println("fwd="+fwd+", fwd2="+fwd2+", rev="+rev+", kmer="+kmer+", rkmer="+rkmer);
				outstream.println("b="+(char)b+", x="+x+", x2="+x2+", bitsPerBase="+bitsPerBase+", shift2="+shift2);
				if(!amino){
					assert(AminoAcid.stringToKmer(fwd)==kmer) : fwd+", "+AminoAcid.stringToKmer(fwd)+", "+kmer+", "+len;
					if(len>=k){
						assert(rcomp(kmer, Tools.min(len, k))==rkmer);
						assert(rcomp(rkmer, Tools.min(len, k))==kmer);
						assert(AminoAcid.kmerToString(kmer, Tools.min(len, k)).equals(fwd));
						assert(AminoAcid.kmerToString(rkmer, Tools.min(len, k)).equals(rev)) : AminoAcid.kmerToString(rkmer, Tools.min(len, k))+" != "+rev+" (rkmer)";
					}
					assert(fwd.equalsIgnoreCase(fwd2)) : fwd+", "+fwd2; //may be unsafe
				}
				outstream.println("Scanning6 i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+", bases="+fwd+", rbases="+rev);
			}
			if(len>=minlen2 && i>=minlen){
				final int id=index.getValue(kmer, rkmer, kmask, i, k, qHammingDistance);
				if(verbose){outstream.println("Testing kmer "+kmer+"; id="+id);}
				if(id>0){
					if(verbose){outstream.println("Found = "+(found+1)+"/"+maxBadKmers);}
					if(found==maxBadKmers){
						if(scaffoldReadCountsT!=null){
							scaffoldReadCountsT[id]++;
							scaffoldBaseCountsT[id]+=bases.length;
						}else{
							scaffoldReadCounts.addAndGet(id, 1);
							scaffoldBaseCounts.addAndGet(id, bases.length);
						}
						return (found=found+1);
						//Early exit, but prevents generation of histogram that goes over maxBadKmers+1.
					}
					found++;
				}
			}
		}
		return found;
	}


	/**
	 * Counts the number of kmer hits for a read.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return Number of hits
	 */
	private final int countCoveredBases(final Read r, final int minCoveredBases){
		if(r==null || r.length()<k || storedKmers<1){return 0;}
		if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){return 0;}
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		int found=0;
		int len=0;
		int lastFound=-1;
		boolean recorded=false;

		final int start=(restrictRight<1 ? 0 : Tools.max(0, bases.length-restrictRight));
		final int stop=(restrictLeft<1 ? bases.length : Tools.min(bases.length, restrictLeft));

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=start; i<stop; i++){
			byte b=bases[i];
			long x=symbolToNumber0[b];
			long x2=symbolToComplementNumber0[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(forbidNs && !isFullyDefined(b)){len=0; rkmer=0;}else{len++;}
			if(verbose){outstream.println("Scanning6b i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=minlen2 && i>=minlen){
				final int id=index.getValue(kmer, rkmer, kmask, i, k, qHammingDistance);
				if(verbose){outstream.println("Testing kmer "+kmer+"; id="+id);}
				if(id>0){

					int extra=Tools.min(k, i-lastFound);
					found+=extra;
					lastFound=i;

					if(verbose){outstream.println("Found = "+found+"/"+minCoveredBases);}
					if(found>=minCoveredBases){
						if(!recorded){
							if(scaffoldReadCountsT!=null){
								scaffoldReadCountsT[id]++;
								scaffoldBaseCountsT[id]+=bases.length;
							}else{
								scaffoldReadCounts.addAndGet(id, 1);
								scaffoldBaseCounts.addAndGet(id, bases.length);
							}
						}
						return found;
					}
				}
			}
		}
		return found;
	}

	/**
	 * Returns the id of the sequence with the most kmer matches to this read, or -1 if none are over maxBadKmers.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return id of best match
	 */
	private final int findBestMatch(final Read r, final int maxBadKmers){
		idList.size=0;
		if(r==null || r.length()<k || storedKmers<1){return -1;}
		if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){return -1;}
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		int len=0;
		int found=0;

		final int start=(restrictRight<1 ? 0 : Tools.max(0, bases.length-restrictRight));
		final int stop=(restrictLeft<1 ? bases.length : Tools.min(bases.length, restrictLeft));

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=start; i<stop; i++){
			byte b=bases[i];
			long x=symbolToNumber0[b];
			long x2=symbolToComplementNumber0[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(forbidNs && !isFullyDefined(b)){len=0; rkmer=0;}else{len++;}
			if(verbose){outstream.println("Scanning6 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=minlen2 && i>=minlen){
				if(verbose){outstream.println("Lookup kmer="+AminoAcid.kmerToString(kmer, k)+", rkmer="+AminoAcid.kmerToString(rkmer, k));}
				final int id=index.getValue(kmer, rkmer, kmask, i, k, qHammingDistance);
				if(id>0){
					countArray[id]++;
					if(countArray[id]==1){idList.add(id);}
					found++;
					if(verbose){outstream.println("Found = "+found+"/"+maxBadKmers);}
				}
			}
		}

		final int id, max;
		if(found>maxBadKmers){
			max=condenseLoose(countArray, idList, countList);
			int id0=-1;
			for(int i=0; i<countList.size; i++){
				if(countList.get(i)==max){
					id0=idList.get(i); break;
				}
			}
			if(rename){rename(r, idList, countList);}
			id=id0;
		}else{
			max=0;
			id=-1;
		}

		if(found>maxBadKmers){
			if(scaffoldReadCountsT!=null){
				scaffoldReadCountsT[id]++;
				scaffoldBaseCountsT[id]+=bases.length;
			}else{
				scaffoldReadCounts.addAndGet(id, 1);
				scaffoldBaseCounts.addAndGet(id, bases.length);
			}
		}
		return id;
	}

	/** Estimates kmer hit counts for kmers longer than k using consecutive matches
	 * @param r
	 * @param sets
	 * @return Number of sets of consecutive hits of exactly length kbig
	 */
	private final int countSetKmersBig(final Read r, final int maxBadKmers){
		if(r==null || r.length()<kbig || storedKmers<1){return 0;}
		if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){return 0;}
		assert(kbig>k);
		final int sub=kbig-k-1;
		assert(sub>=0) : kbig+", "+sub;
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		int found=0;
		int len=0;

		int bkStart=-1;
		int bkStop=-1;
		int id=-1, lastId=-1;

		final int start=(restrictRight<1 ? 0 : Tools.max(0, bases.length-restrictRight));
		final int stop=(restrictLeft<1 ? bases.length : Tools.min(bases.length, restrictLeft));

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=start; i<stop; i++){
			byte b=bases[i];
			long x=symbolToNumber0[b];
			long x2=symbolToComplementNumber0[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(forbidNs && !isFullyDefined(b)){len=0; rkmer=0;}else{len++;}
			if(verbose){outstream.println("Scanning7 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=minlen2 && i>=minlen){
				id=index.getValue(kmer, rkmer, kmask, i, k, qHammingDistance);
				if(verbose){outstream.println("Testing kmer "+kmer+"; id="+id);}
				if(id>0){
					lastId=id;
					if(bkStart==-1){bkStart=i;}
					bkStop=i;
				}else{
					if(bkStart>-1){
						int dif=bkStop-bkStart-sub;
						bkStop=bkStart=-1;
						if(dif>0){
							int old=found;
							found+=dif;
							if(found>maxBadKmers && old<=maxBadKmers){
								if(scaffoldReadCountsT!=null){
									scaffoldReadCountsT[lastId]++;
									scaffoldBaseCountsT[lastId]+=bases.length;
								}else{
									scaffoldReadCounts.addAndGet(lastId, 1);
									scaffoldBaseCounts.addAndGet(lastId, bases.length);
								}
								return found;
								//Early exit, but prevents generation of histogram that goes over maxBadKmers+1.
							}
						}
					}
				}
			}
		}

		// This catches the case where valid kmers extend to the end of the read
		if(bkStart>-1){
			int dif=bkStop-bkStart-sub;
			bkStop=bkStart=-1;
			if(dif>0){
				int old=found;
				found+=dif;
				if(found>maxBadKmers && old<=maxBadKmers){
					if(scaffoldReadCountsT!=null){
						scaffoldReadCountsT[lastId]++;
						scaffoldBaseCountsT[lastId]+=bases.length;
					}else{
						scaffoldReadCounts.addAndGet(lastId, 1);
						scaffoldBaseCounts.addAndGet(lastId, bases.length);
					}
				}
			}
		}
		return found;
	}

	private final int ktrim(final Read r){
		final int len=r.length();
		final int start=(restrictRight<1 ? 0 : Tools.max(0, len-restrictRight));
		final int stop=(restrictLeft<1 ? len : Tools.min(len, restrictLeft));
		return ktrim(r, start, stop);
	}

	private final int ktrimTips(final Read r){
		final int len=r.length();
		final int mid=len/2-(k-1)/2;
		int sum=0;
		if(ktrimRight){
			int start=Tools.max(0, (restrictRight<1 ? mid : len-restrictRight));
			sum+=ktrimTip(r, start, len, true, false);
		}
		if(ktrimLeft){
			int stop=Tools.min(r.length(), (restrictLeft<1 ? mid+k-1 : restrictLeft));
			sum+=ktrimTip(r, 0, stop, false, true);
		}
		return sum;
	}

	/**
	 * Trim a read to remove matching kmers and everything to their left or right, on one end.
	 * Allows both trimRight and trimLeft, in 2 passes.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return Number of bases trimmed
	 */
	private final int ktrimTip(final Read r, final int start, final int stop,
		final boolean right, final boolean left){
		assert(left || right);
		assert(!(left && right));
		if(r==null || r.length()<Tools.max(1, (useShortKmers ? Tools.min(k, mink) : k)) || storedKmers<1){return 0;}
		if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){return 0;}
		if(verbose){outstream.println("KTrimming read "+r.id);}
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		int found=0;
		int len=0;
		int id0=-1; //ID of first kmer found.

		int minLoc=999999999, minLocExclusive=999999999;
		int maxLoc=-1, maxLocExclusive=-1;

		//Scan for normal kmers
		for(int i=start; i<stop; i++){
			byte b=bases[i];
			long x=symbolToNumber0[b];
			long x2=symbolToComplementNumber0[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(forbidNs && !isFullyDefined(b)){len=0; rkmer=0;}else{len++;}
			if(verbose){outstream.println("Scanning3 i="+i+", kmer="+kmer+", rkmer="+rkmer+", len="+len+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=minlen2 && i>=minlen){
				final int id=index.getValue(kmer, rkmer, kmask, i, k, qHammingDistance);
				if(id>0){
					if(id0<0){id0=id;}
					minLoc=Tools.min(minLoc, i-k+1);
					assert(minLoc>=0);
					maxLoc=i;
					found++;
				}
			}
		}

		if(minLoc!=minLocExclusive){minLocExclusive=minLoc+k;}
		if(maxLoc!=maxLocExclusive){maxLocExclusive=maxLoc-k;}

		//If nothing was found, scan for short kmers.  Only used for trimming.
		if(useShortKmers && found==0){
			assert(!maskMiddle && middleMask==-1) : midMaskLen+", "+middleMask+", k="+", mink="+mink;

			//Look for short kmers on left side
			if(left){
				kmer=0;
				rkmer=0;
				len=0;
				final int lim=Tools.min(k, stop);
				for(int i=start; i<lim; i++){
					byte b=bases[i];
					long x=symbolToNumber0[b];
					long x2=symbolToComplementNumber0[b];
					kmer=((kmer<<bitsPerBase)|x)&mask;
					rkmer=rkmer|(x2<<(bitsPerBase*len));
					len++;
					if(verbose){outstream.println("Scanning4 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=mink){

						if(verbose){
							outstream.println("Looking for left kmer  "+kmerToString(kmer, len));
							outstream.println("Looking for left rkmer "+kmerToString(rkmer, len));
						}

						final int id=index.getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2);
						if(id>0){
							if(id0<0){id0=id;}
							if(verbose){outstream.println("Found "+kmer);}
							minLoc=0;
							minLocExclusive=Tools.min(minLocExclusive, i+1);
							maxLoc=Tools.max(maxLoc, i);
							maxLocExclusive=Tools.max(maxLocExclusive, 0);
							found++;
						}
					}
				}
			}

			//Look for short kmers on right side
			if(right){
				kmer=0;
				rkmer=0;
				len=0;
				final int lim=Tools.max(-1, stop-k);
				for(int i=stop-1; i>lim; i--){
					byte b=bases[i];
					long x=symbolToNumber0[b];
					long x2=symbolToComplementNumber0[b];
					kmer=kmer|(x<<(bitsPerBase*len));
					rkmer=((rkmer<<bitsPerBase)|x2)&mask;
					len++;
					if(verbose){outstream.println("Scanning5 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=mink){
						if(verbose){
							outstream.println("Looking for right kmer "+
								AminoAcid.kmerToString(kmer&~lengthMasks[len], len)+"; value="+toValue(kmer, rkmer, lengthMasks[len])+"; kmask="+lengthMasks[len]);
						}
						final int id=index.getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2);
						if(id>0){
							if(id0<0){id0=id;}
							if(verbose){outstream.println("Found "+kmer);}
							minLoc=i;
							minLocExclusive=Tools.min(minLocExclusive, bases.length);
							maxLoc=bases.length-1;
							maxLocExclusive=Tools.max(maxLocExclusive, i-1);
							found++;
						}
					}
				}
			}
		}


		if(verbose){outstream.println("found="+found+", minLoc="+minLoc+", maxLoc="+maxLoc+
			", minLocExclusive="+minLocExclusive+", maxLocExclusive="+maxLocExclusive);}

		if(found==0){return 0;}
		assert(found>0) : "Overflow in 'found' variable.";

		{//Increment counter for the scaffold whose kmer was first detected
			if(scaffoldReadCountsT!=null){
				scaffoldReadCountsT[id0]++;
				scaffoldBaseCountsT[id0]+=bases.length;
			}else{
				scaffoldReadCounts.addAndGet(id0, 1);
				scaffoldBaseCounts.addAndGet(id0, bases.length);
			}
		}

		if(trimPad!=0){
			maxLoc=Tools.mid(0, maxLoc+trimPad, bases.length);
			minLoc=Tools.mid(0, minLoc-trimPad, bases.length);
			maxLocExclusive=Tools.mid(0, maxLocExclusive+trimPad, bases.length);
			minLocExclusive=Tools.mid(0, minLocExclusive-trimPad, bases.length);
		}

		if(left){ //Trim from the read start to the rightmost kmer base
			if(verbose){outstream.println("Left trimming to "+(ktrimExclusive ? maxLocExclusive+1 : maxLoc+1)+", "+0);}
			int x=TrimRead.trimToPosition(r, ktrimExclusive ? maxLocExclusive+1 : maxLoc+1, bases.length-1, 1);
			if(verbose){outstream.println("Trimmed "+x+" bases: "+new String(r.bases));}
			return x;
		}else{ //Trim from the leftmost kmer base to the read stop
			assert(right);
			if(verbose){outstream.println("Right trimming to "+0+", "+(ktrimExclusive ? minLocExclusive-1 : minLoc-1));}
			int x=TrimRead.trimToPosition(r, 0, ktrimExclusive ? minLocExclusive-1 : minLoc-1, 1);
			if(verbose){outstream.println("Trimmed "+x+" bases: "+new String(r.bases));}
			return x;
		}
	}

	/**
	 * Trim a read to remove matching kmers and everything to their left or right.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return Number of bases trimmed
	 */
	private final int ktrim(final Read r, final int start, final int stop){
		assert(ktrimLeft || ktrimRight);
		if(r==null || r.length()<Tools.max(1, (useShortKmers ? Tools.min(k, mink) : k)) || storedKmers<1){return 0;}
		if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){return 0;}
		if(verbose){outstream.println("KTrimming read "+r.id);}
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		int found=0;
		int len=0;
		int id0=-1; //ID of first kmer found.

		int minLoc=999999999, minLocExclusive=999999999;
		int maxLoc=-1, maxLocExclusive=-1;

		//Scan for normal kmers
		for(int i=start; i<stop; i++){
			byte b=bases[i];
			long x=symbolToNumber0[b];
			long x2=symbolToComplementNumber0[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(forbidNs && !isFullyDefined(b)){len=0; rkmer=0;}else{len++;}
			if(verbose){outstream.println("Scanning3 i="+i+", kmer="+kmer+", rkmer="+rkmer+", len="+len+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=minlen2 && i>=minlen){
				final int id=index.getValue(kmer, rkmer, kmask, i, k, qHammingDistance);
				if(id>0){
					if(id0<0){id0=id;}
					minLoc=Tools.min(minLoc, i-k+1);
					assert(minLoc>=0);
					maxLoc=i;
					found++;
				}
			}
		}

		if(minLoc!=minLocExclusive){minLocExclusive=minLoc+k;}
		if(maxLoc!=maxLocExclusive){maxLocExclusive=maxLoc-k;}

		//If nothing was found, scan for short kmers.  Only used for trimming.
		if(useShortKmers && found==0){
			assert(!maskMiddle && middleMask==-1) : midMaskLen+", "+middleMask+", k="+", mink="+mink;

			//Look for short kmers on left side
			if(ktrimLeft){
				kmer=0;
				rkmer=0;
				len=0;
				final int lim=Tools.min(k, stop);
				for(int i=start; i<lim; i++){
					byte b=bases[i];
					long x=symbolToNumber0[b];
					long x2=symbolToComplementNumber0[b];
					kmer=((kmer<<bitsPerBase)|x)&mask;
					rkmer=rkmer|(x2<<(bitsPerBase*len));
					len++;
					if(verbose){outstream.println("Scanning4 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=mink){

						if(verbose){
							outstream.println("Looking for left kmer  "+kmerToString(kmer, len));
							outstream.println("Looking for left rkmer "+kmerToString(rkmer, len));
						}

						final int id=index.getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2);
						if(id>0){
							if(id0<0){id0=id;}
							if(verbose){outstream.println("Found "+kmer);}
							minLoc=0;
							minLocExclusive=Tools.min(minLocExclusive, i+1);
							maxLoc=Tools.max(maxLoc, i);
							maxLocExclusive=Tools.max(maxLocExclusive, 0);
							found++;
						}
					}
				}
			}

			//Look for short kmers on right side
			if(ktrimRight){
				kmer=0;
				rkmer=0;
				len=0;
				final int lim=Tools.max(-1, stop-k);
				for(int i=stop-1; i>lim; i--){
					byte b=bases[i];
					long x=symbolToNumber0[b];
					long x2=symbolToComplementNumber0[b];
					kmer=kmer|(x<<(bitsPerBase*len));
					rkmer=((rkmer<<bitsPerBase)|x2)&mask;
					len++;
					if(verbose){outstream.println("Scanning5 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=mink){
						if(verbose){
							outstream.println("Looking for right kmer "+
								AminoAcid.kmerToString(kmer&~lengthMasks[len], len)+"; value="+toValue(kmer, rkmer, lengthMasks[len])+"; kmask="+lengthMasks[len]);
						}
						final int id=index.getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2);
						if(id>0){
							if(id0<0){id0=id;}
							if(verbose){outstream.println("Found "+kmer);}
							minLoc=i;
							minLocExclusive=Tools.min(minLocExclusive, bases.length);
							maxLoc=bases.length-1;
							maxLocExclusive=Tools.max(maxLocExclusive, i-1);
							found++;
						}
					}
				}
			}
		}


		if(verbose){outstream.println("found="+found+", minLoc="+minLoc+", maxLoc="+maxLoc+", minLocExclusive="+minLocExclusive+", maxLocExclusive="+maxLocExclusive);}

		if(found==0){return 0;}
		assert(found>0) : "Overflow in 'found' variable.";

		{//Increment counter for the scaffold whose kmer was first detected
			if(scaffoldReadCountsT!=null){
				scaffoldReadCountsT[id0]++;
				scaffoldBaseCountsT[id0]+=bases.length;
			}else{
				scaffoldReadCounts.addAndGet(id0, 1);
				scaffoldBaseCounts.addAndGet(id0, bases.length);
			}
		}

		if(trimPad!=0){
			maxLoc=Tools.mid(0, maxLoc+trimPad, bases.length);
			minLoc=Tools.mid(0, minLoc-trimPad, bases.length);
			maxLocExclusive=Tools.mid(0, maxLocExclusive+trimPad, bases.length);
			minLocExclusive=Tools.mid(0, minLocExclusive-trimPad, bases.length);
		}

		if(ktrimLeft){ //Trim from the read start to the rightmost kmer base
			if(verbose){outstream.println("Left trimming to "+(ktrimExclusive ? maxLocExclusive+1 : maxLoc+1)+", "+0);}
			int x=TrimRead.trimToPosition(r, ktrimExclusive ? maxLocExclusive+1 : maxLoc+1, bases.length-1, 1);
			if(verbose){outstream.println("Trimmed "+x+" bases: "+new String(r.bases));}
			return x;
		}else{ //Trim from the leftmost kmer base to the read stop
			assert(ktrimRight);
			if(verbose){outstream.println("Right trimming to "+0+", "+(ktrimExclusive ? minLocExclusive-1 : minLoc-1));}
			int x=TrimRead.trimToPosition(r, 0, ktrimExclusive ? minLocExclusive-1 : minLoc-1, 1);
			if(verbose){outstream.println("Trimmed "+x+" bases: "+new String(r.bases));}
			return x;
		}
	}


	/**
	 * Mask a read to cover matching kmers.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return Number of bases masked
	 */
	private final int kmask(final Read r){
		assert(ktrimN);
		if(r==null || r.length()<Tools.max(1, (useShortKmers ? Tools.min(k, mink) : k)) || storedKmers<1){return 0;}
		if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){return 0;}
		if(verbose){outstream.println("KMasking read "+r.id);}
		final byte[] bases=r.bases, quals=r.quality;
		if(bases==null || bases.length<k){return 0;}
		long kmer=0;
		long rkmer=0;
		int found=0;
		int len=0;
		int id0=-1; //ID of first kmer found.

		final BitSet bs=new BitSet(bases.length+trimPad+1);
		if(kmaskFullyCovered){bs.set(0, bases.length);}

		final int minus=k-1-trimPad;
		final int plus=trimPad+1;

		final int start=(restrictRight<1 ? 0 : Tools.max(0, bases.length-restrictRight));
		final int stop=(restrictLeft<1 ? bases.length : Tools.min(bases.length, restrictLeft));

		//Scan for normal kmers
		for(int i=start; i<stop; i++){
			byte b=bases[i];
			long x=symbolToNumber0[b];
			long x2=symbolToComplementNumber0[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(forbidNs && !isFullyDefined(b)){len=0; rkmer=0;}else{len++;}
			if(verbose){outstream.println("Scanning3 i="+i+", kmer="+kmer+", rkmer="+rkmer+", len="+len+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}

			if(i>=minlen){
				final int id;
				if(len>=minlen2){
					id=index.getValue(kmer, rkmer, kmask, i, k, qHammingDistance);
				}else{
					id=-1;
				}
				if(id>0){
					if(id0<0){id0=id;}
					if(verbose){
						outstream.println("a: Found "+kmer);
						outstream.println("Setting "+Tools.max(0, i-minus)+", "+(i+plus));
						outstream.println("i="+i+", minus="+minus+", plus="+plus+", trimpad="+trimPad+", k="+k);
					}
					if(!kmaskFullyCovered){bs.set(Tools.max(0, i-minus), i+plus);}
					found++;
				}else if(kmaskFullyCovered){
					bs.clear(Tools.max(0, i-minus), i+plus);
				}
			}
		}

		//If nothing was found, scan for short kmers.
		if(useShortKmers){
			assert(!maskMiddle && middleMask==-1) : midMaskLen+", "+middleMask+", k="+", mink="+mink;

			//Look for short kmers on left side
			{
				kmer=0;
				rkmer=0;
				len=0;
				int len2=0;
				final int lim=Tools.min(k, stop);
				for(int i=start; i<lim; i++){
					byte b=bases[i];
					long x=symbolToNumber0[b];
					long x2=symbolToComplementNumber0[b];
					kmer=((kmer<<bitsPerBase)|x)&mask;
					rkmer=rkmer|(x2<<(bitsPerBase*len));
					len++;
					len2++;
					if(verbose){outstream.println("Scanning4 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}

					if(len2>=minminlen){
						if(verbose){
							outstream.println("Looking for left kmer  "+kmerToString(kmer, len));
							outstream.println("Looking for left rkmer "+kmerToString(rkmer, len));
						}
						final int id;
						if(len>=mink){
							id=index.getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2);
						}else{
							id=-1;
						}
						if(id>0){
							if(id0<0){id0=id;}
							if(verbose){
								outstream.println("b: Found "+kmer);
								outstream.println("Setting "+0+", "+Tools.min(bases.length, i+trimPad+1));
							}
							if(!kmaskFullyCovered){bs.set(0, Tools.min(bases.length, i+trimPad+1));}
							found++;
						}else if(kmaskFullyCovered){
							bs.clear(0, Tools.min(bases.length, i+trimPad+1));
						}
					}
				}
			}

			//Look for short kmers on right side
			{
				kmer=0;
				rkmer=0;
				len=0;
				int len2=0;
				final int lim=Tools.max(-1, stop-k);
				for(int i=stop-1; i>lim; i--){
					byte b=bases[i];
					long x=symbolToNumber0[b];
					long x2=symbolToComplementNumber0[b];
					kmer=kmer|(x<<(bitsPerBase*len));
					rkmer=((rkmer<<bitsPerBase)|x2)&mask;
					len++;
					len2++;
					if(verbose){outstream.println("Scanning5 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}

					if(len2>=minminlen){
						if(verbose){
							outstream.println("Looking for right kmer "+
								AminoAcid.kmerToString(kmer&~lengthMasks[len], len)+"; value="+toValue(kmer, rkmer, lengthMasks[len])+"; kmask="+lengthMasks[len]);
						}
						final int id;
						if(len>=mink){
							id=index.getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2);
						}else{
							id=-1;
						}
						if(id>0){
							if(id0<0){id0=id;}
							if(verbose){
								outstream.println("c: Found "+kmer);
								outstream.println("Setting "+Tools.max(0, i-trimPad)+", "+bases.length);
							}
							if(!kmaskFullyCovered){bs.set(Tools.max(0, i-trimPad), bases.length);}
							found++;
						}else if(kmaskFullyCovered){
							bs.clear(Tools.max(0, i-trimPad), bases.length);
						}
					}
				}
			}
		}


		if(verbose){outstream.println("found="+found+", bitset="+bs);}

		if(found==0){return 0;}
		assert(found>0) : "Overflow in 'found' variable.";

		{//Increment counter for the scaffold whose kmer was first detected
			if(scaffoldReadCountsT!=null){
				scaffoldReadCountsT[id0]++;
				scaffoldBaseCountsT[id0]+=bases.length;
			}else{
				scaffoldReadCounts.addAndGet(id0, 1);
				scaffoldBaseCounts.addAndGet(id0, bases.length);
			}
		}
		int cardinality=bs.cardinality();

		//Replace kmer hit zone with the trim symbol
		for(int i=0; i<bases.length; i++){
			if(bs.get(i)){
				if(kmaskLowercase){
					bases[i]=(byte)Tools.toLowerCase(bases[i]);
				}else{
					bases[i]=trimSymbol;
					if(quals!=null && trimSymbol=='N'){quals[i]=0;}
				}
			}
		}
		return cardinality;
	}


	/**
	 * Mask a read to cover matching kmers.
	 * @param r Read to process
	 * @param sets Kmer tables
	 * @return Number of bases masked
	 */
	private final boolean ksplit(final Read r){
		if(r==null || r.length()<Tools.max(1, (useShortKmers ? Tools.min(k, mink) : k)) || storedKmers<1){return false;}
		assert(r.mate==null) : "Kmer splitting should only be performed on unpaired reads.";
		assert(ksplit);
		if(verbose){outstream.println("KSplitting read "+r.id);}
		final byte[] bases=r.bases;
		if(bases==null || bases.length<k){return false;}
		long kmer=0;
		long rkmer=0;
		long found=0;
		int len=0;
		int id0=-1; //ID of first kmer found.
		int leftmost=Integer.MAX_VALUE, rightmost=-1;

		final int minus=k-1-trimPad;
		final int plus=trimPad;

		final int start=(restrictRight<1 ? 0 : Tools.max(0, bases.length-restrictRight));
		final int stop=(restrictLeft<1 ? bases.length : Tools.min(bases.length, restrictLeft));

		//Scan for normal kmers
		for(int i=start; i<stop; i++){
			byte b=bases[i];
			long x=symbolToNumber0[b];
			long x2=symbolToComplementNumber0[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(forbidNs && !isFullyDefined(b)){len=0; rkmer=0;}else{len++;}
			if(verbose){outstream.println("Scanning3 i="+i+", kmer="+kmer+", rkmer="+rkmer+", len="+len+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}

			if(i>=minlen){
				final int id;
				if(len>=minlen2){
					id=index.getValue(kmer, rkmer, kmask, i, k, qHammingDistance);
				}else{
					id=-1;
				}
				if(id>0){
					if(id0<0){id0=id;}
					if(verbose){
						outstream.println("a: Found "+kmer);
						outstream.println("Setting "+Tools.max(0, i-minus)+", "+(i+plus));
						outstream.println("i="+i+", minus="+minus+", plus="+plus+", trimpad="+trimPad+", k="+k);
					}
					leftmost=Tools.min(leftmost, Tools.max(0, i-minus));
					rightmost=Tools.max(rightmost, i+plus);
					found++;
				}
			}
		}

		//If nothing was found, scan for short kmers.
		if(useShortKmers && id0==-1){
			assert(!maskMiddle && middleMask==-1) : midMaskLen+", "+middleMask+", k="+", mink="+mink;

			//Look for short kmers on right side
			{
				kmer=0;
				rkmer=0;
				len=0;
				int len2=0;
				final int lim=Tools.max(-1, stop-k);
				for(int i=stop-1; i>lim; i--){
					byte b=bases[i];
					long x=symbolToNumber0[b];
					long x2=symbolToComplementNumber0[b];
					kmer=kmer|(x<<(bitsPerBase*len));
					rkmer=((rkmer<<bitsPerBase)|x2)&mask;
					len++;
					len2++;
					if(verbose){outstream.println("Scanning5 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}

					if(len2>=minminlen){
						if(verbose){
							outstream.println("Looking for right kmer "+
								AminoAcid.kmerToString(kmer&~lengthMasks[len], len)+"; value="+toValue(kmer, rkmer, lengthMasks[len])+"; kmask="+lengthMasks[len]);
						}
						final int id;
						if(len>=mink){
							id=index.getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2);
						}else{
							id=-1;
						}
						if(id>0){
							if(id0<0){id0=id;}
							if(verbose){
								outstream.println("b: Found "+kmer);
								outstream.println("Setting "+Tools.max(0, i-trimPad)+", "+bases.length);
							}
							leftmost=Tools.min(leftmost, Tools.max(0, i-trimPad));
							rightmost=bases.length-1;
							found++;
						}
					}
				}
			}

			//Look for short kmers on left side
			if(id0==-1){
				kmer=0;
				rkmer=0;
				len=0;
				int len2=0;
				final int lim=Tools.min(k, stop);
				for(int i=start; i<lim; i++){
					byte b=bases[i];
					long x=symbolToNumber0[b];
					long x2=symbolToComplementNumber0[b];
					kmer=((kmer<<bitsPerBase)|x)&mask;
					rkmer=rkmer|(x2<<(bitsPerBase*len));
					len++;
					len2++;
					if(verbose){outstream.println("Scanning4 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}

					if(len2>=minminlen){
						if(verbose){
							outstream.println("Looking for left kmer  "+kmerToString(kmer, len));
							outstream.println("Looking for left rkmer "+kmerToString(rkmer, len));
						}
						final int id;
						if(len>=mink){
							id=index.getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2);
						}else{
							id=-1;
						}
						if(id>0){
							if(id0<0){id0=id;}
							if(verbose){
								outstream.println("c: Found "+kmer);
								outstream.println("Setting "+0+", "+(i+trimPad));
							}
							leftmost=0;
							rightmost=Tools.max(rightmost, i+trimPad);
							found++;
						}
					}
				}
			}
		}


		if(verbose){outstream.println("found="+found);}

		if(found==0){return false;}

		{//Increment counter for the scaffold whose kmer was first detected
			if(scaffoldReadCountsT!=null){
				scaffoldReadCountsT[id0]++;
				scaffoldBaseCountsT[id0]+=bases.length;
			}else{
				scaffoldReadCounts.addAndGet(id0, 1);
				scaffoldBaseCounts.addAndGet(id0, bases.length);
			}
		}

		if(leftmost==0){
			TrimRead.trimToPosition(r, rightmost+1, bases.length-1, 1);
			return false;
		}else if(rightmost==bases.length-1){
			TrimRead.trimToPosition(r, 0, leftmost-1, 1);
			return false;
		}else{
			Read r2=r.subRead(rightmost+1, bases.length-1);
			TrimRead.trimToPosition(r, 0, leftmost-1, 1);
			r.mate=r2;
			r2.mate=r;
			r2.setPairnum(1);
			return true;
		}
	}

	/**
	 * @param r
	 * @param idList
	 * @param countList
	 */
	private void rename(Read r, IntList idList, IntList countList) {
		if(r==null || idList.size<1){return;}
		StringBuilder sb=new StringBuilder();
		if(r.id==null){sb.append(r.numericID);}
		else{sb.append(r.id);}
		for(int i=0; i<idList.size; i++){
			int id=idList.get(i);
			int count=countList.get(i);
			sb.append('\t');
			sb.append(scaffoldNames.get(id));
			sb.append('=');
			sb.append(count);
		}
		r.id=sb.toString();
	}

	/**
	 * Pack a list of counts from an array to an IntList.
	 * @param loose Counter array
	 * @param packed Unique values
	 * @param counts Counts of values
	 * @return
	 */
	private int condenseLoose(int[] loose, IntList packed, IntList counts){
		counts.size=0;
		if(packed.size<1){return 0;}

		int max=0;
		for(int i=0; i<packed.size; i++){
			final int p=packed.get(i);
			final int c=loose[p];
			counts.add(c);
			loose[p]=0;
			max=Tools.max(max, c);
		}
		return max;
	}

	private float expectedErrors(Read r1, Read r2){
		float a=(r1==null ? 0 : r1.expectedErrors(false, -1));
		float b=(r2==null ? 0 : r2.expectedErrors(false, -1));
		return Tools.max(a, b);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Entropy Methods       ----------------*/
	/*--------------------------------------------------------------*/

	private int maskLowEntropy(final Read r, BitSet bs, EntropyTracker et){
		final int window=et.windowBases();
		if(r==null || r.length()<window){return 0;}
		final byte[] bases=r.bases;
		if(bs==null){bs=new BitSet(r.length());}
		else{bs.clear();}

		et.clear();
		for(int i=0, min=window-1; i<bases.length; i++){
			et.add(bases[i]);
			if(i>=min && et.ns()<1 && !et.passes()){bs.set(et.leftPos(), et.rightPos()+1);}
		}

		return maskFromBitset(r, bs, entropyMaskLowercase);
	}

	private int trimLowEntropy(final Read r, BitSet bs, EntropyTracker et){
		final int window=et.windowBases();
		System.err.println("Trimming "+r.id+", len "+r.length()+", window "+window);
		if(r==null || r.length()<window){return 0;}
		final byte[] bases=r.bases;
		if(bs==null){bs=new BitSet(r.length());}
		else{bs.clear();}

		et.clear();
		for(int i=0, min=window-1; i<bases.length; i++){
			et.add(bases[i]);
			if(i>=min && et.ns()<1 && !et.passes()){bs.set(et.leftPos(), et.rightPos()+1);}
		}

		//Now, trim from bitset - which could be spun out into a function
		int left=0, right=0;
		final int len=bases.length;
		for(int i=0; i<len; i++){
			if(bs.get(i)){left++;}
			else {break;}
		}
		for(int i=len-1; i>=0; i--){
			if(bs.get(i)){right++;}
			else {break;}
		}
		if(left==0 && right==0){return 0;}
		return TrimRead.trimByAmount(r, left, right, 1);
	}

	private void markLowEntropy(final Read r, EntropyTracker et){
		final int window=et.windowBases();
		if(r==null || r.length()<window){return;}
		final byte[] bases=r.bases;

		float[] values=new float[r.length()];
		Arrays.fill(values, 1);

		et.clear();
		for(int i=0, min=window-1; i<bases.length; i++){
			et.add(bases[i]);
			if(i>=min && et.ns()<1){
				float e=et.calcEntropy();
				for(int j=et.leftPos(), max=et.rightPos(); j<=max; j++){
					values[j]=Tools.min(e, values[j]);
				}
			}
		}

		if(r.quality==null){
			r.quality=new byte[r.length()];
		}
		for(int i=0; i<values.length; i++){
			byte q=(byte)(values[i]*41);
			r.quality[i]=q;
		}
	}

	private int maskFromBitset(final Read r, final BitSet bs, final boolean lowercase){
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		int sum=0;
		if(!lowercase){
			for(int i=bs.nextSetBit(0); i>=0; i=bs.nextSetBit(i+1)){
				if(bases[i]!='N'){
					sum++;
					bases[i]='N';
					if(quals!=null){quals[i]=0;}
				}
			}
		}else{
			for(int i=bs.nextSetBit(0); i>=0; i=bs.nextSetBit(i+1)){
				if(!Tools.isLowerCase(bases[i])){
					if(bases[i]!='N'){sum++;}
					bases[i]=(byte)Tools.toLowerCase(bases[i]);//Don't change quality
				}
			}
		}
		return sum;
	}

	public final boolean passesVariantFilter(Read r){
		if(!r.mapped() || r.bases==null || r.samline==null || r.match==null){return true;}
		//TODO: Add Vars as well, like in FilterSam
		if(Read.countSubs(r.match)<=maxBadSubs){return true;}
		ArrayList<Var> list=AnalyzeVars.findUniqueSubs(r, r.samline, varMap, scafMap, maxBadSubAlleleDepth, maxBadAlleleFraction, minBadSubReadDepth, minBadSubEDist);
		return list==null || list.size()<=maxBadSubs;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses literal sequence argument into array of sequences.
	 * @param arg Comma-separated list of literal sequences
	 * @return Array of processed literal sequences
	 */
	public static final String[] processLiteralArg(String arg) {
		if(arg==null) {return null;}
		String[] split=arg.split(",");
		ArrayList<String> list=new ArrayList<String>(split.length);
		for(String b : split) {
			String c=processLiteralTerm(b);
			if(c!=null) {list.add(c);}
		}
		String[] ret=list.isEmpty() ? null : list.toArray(new String[0]);
		return ret;
	}

	/**
	 * Processes individual literal sequence term, expanding polymer shortcuts.
	 * @param b Single literal sequence or polymer specification
	 * @return Processed literal sequence
	 */
	public static final String processLiteralTerm(String b) {
		assert(b.length()>0) : "Invalid literal sequence: '"+b+"'";
		if(AminoAcid.isACGTN(b)) {
			return b;
		}
		b=b.toUpperCase().replaceAll("-", "");
		if(b.startsWith("POLY")) {
			b=b.replace("POLY", "");
			assert(AminoAcid.isACGTN(b)) : "Invalid literal sequence: '"+b+"'";
			StringBuilder sb=new StringBuilder(40);
			final int minlen=Tools.max(31+b.length(), b.length()*3);
			while(sb.length()<minlen) {sb.append(b);}
			System.err.println("Adding literal polymer "+sb);
			return sb.toString();
		}else {
			assert(false) : "Invalid literal sequence: '"+b+"'";
		}
		return null;
	}

	/**
	 * Transforms a kmer into a canonical value stored in the table.  Expected to be inlined.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @return Canonical value
	 */
	final long toValue(long kmer, long rkmer, long lengthMask){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : 
			"\n"+Long.toBinaryString(lengthMask)+
			"\n"+Long.toBinaryString(kmer)+
			"\n"+Long.toBinaryString(rkmer)+
			"\n"+Long.toBinaryString(rcomp(kmer, k));
		if(verbose){outstream.println("toValue("+AminoAcid.kmerToString(kmer, k)+", "+AminoAcid.kmerToString(rkmer, k)+", "+lengthMask+")");}
		final long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		if(verbose){outstream.println("value="+AminoAcid.kmerToString(value, k)+" = "+value);}
		final long ret=(value&middleMask)|lengthMask;
		if(verbose){outstream.println("ret="+AminoAcid.kmerToString(ret, k)+" = "+ret);}
		return ret;
	}

	/**
	 * Computes reverse complement of kmer.
	 * @param kmer Input kmer
	 * @param len Kmer length
	 * @return Reverse complement kmer
	 */
	final long rcomp(long kmer, int len){
		return amino ? kmer : AminoAcid.reverseComplementBinaryFast(kmer, len);
	}

	/**
	 * Determines if kmer passes speed filtering threshold.
	 * @param key Kmer hash value
	 * @return true if kmer should be processed
	 */
	final boolean passesSpeed(long key){
		return speed<1 || ((key&Long.MAX_VALUE)%17)>=speed;
	}

	/**
	 * Determines if kmer fails speed filtering threshold.
	 * @param key Kmer hash value
	 * @return true if kmer should be skipped
	 */
	final boolean failsSpeed(long key){
		return speed>0 && ((key&Long.MAX_VALUE)%17)<speed;
	}

	public static int trimPolyA(final Read r, final int minPoly){
		assert(minPoly>0);
		if(r==null || r.length()<minPoly){return 0;}

		int left=Tools.max(r.countLeft((byte)'A'), r.countLeft((byte)'T'));
		int right=Tools.max(r.countRight((byte)'A'), r.countRight((byte)'T'));

		if(left<minPoly){left=0;}
		if(right<minPoly){right=0;}
		int trimmed=0;
		if(left>0 || right>0){
			trimmed=TrimRead.trimByAmount(r, left, right, 1);
		}
		return trimmed;
	}

	public static int trimPoly(final Read r, final int minPolyLeft, final int minPolyRight, 
		int maxNonPoly, final byte c){
		assert(minPolyLeft>0 || minPolyRight>0);
		if(r==null){return 0;}

		int left=minPolyLeft>0 ? detectPolyLeft(r, minPolyLeft, maxNonPoly, c) : 0;
		int right=minPolyRight>0 ? detectPolyRight(r, minPolyRight, maxNonPoly, c) : 0;

		int trimmed=0;
		if(left>0 || right>0){
			trimmed=TrimRead.trimByAmount(r, left, right, 1);
		}
		return trimmed;
	}

	public static int detectPolyLeft(final Read r, final int minPoly, final int maxNonPoly, final byte c){
		assert(minPoly>0);
		final byte[] bases=r.bases;
		if(bases==null || bases.length<minPoly) {return 0;}

		int trimTo=-1;//Inclusive
		for(int i=0, polymer=0, nonpoly=0; i<bases.length && nonpoly<=maxNonPoly; i++) {
			final byte b=bases[i];
			if(b==c) {
				polymer++;
				if(polymer>=minPoly) {
					nonpoly=0;//Only reset when a valid homopolymer occurs
					trimTo=i;
				}
			}else{
				polymer=0;
				nonpoly++;
			}
		}
		return trimTo+1;
	}

	public static int detectPolyRight(final Read r, final int minPoly, final int maxNonPoly, final byte c){
		assert(minPoly>0);
		final byte[] bases=r.bases;
		if(bases==null || bases.length<minPoly) {return 0;}

		int trimTo=bases.length;//Inclusive
		for(int i=bases.length-1, polymer=0, nonpoly=0; i>=0 && nonpoly<=maxNonPoly; i--) {
			final byte b=bases[i];
			if(b==c) {
				polymer++;
				if(polymer>=minPoly) {
					nonpoly=0;
					trimTo=i;
				}
			}else{
				polymer=0;
				nonpoly++;
			}
		}
		return bases.length-trimTo;
	}

	private static int trimSwift(Read r){
		int left=0, right=0, trimmed=0;
		if(r.pairnum()==0){
			for(int i=r.length()-1; i>=0; i--){
				byte b=r.bases[i];
				if(b=='C' || b=='T' || b=='N'){right++;}
				else{break;}
			}

		}else{
			for(int i=0; i<r.length(); i++){
				byte b=r.bases[i];
				if(b=='G' || b=='A' || b=='N'){left++;}
				else{break;}
			}
		}
		if(left>0 || right>0){
			trimmed=TrimRead.trimByAmount(r, left, right, 1);
		}
		return trimmed;
	}

	/** For verbose / debugging output */
	final String kmerToString(long kmer, int k){
		return amino ? AminoAcid.kmerToStringAA(kmer, k) : AminoAcid.kmerToString(kmer, k);
	}

	/** Returns true if the symbol is not degenerate (e.g., 'N') for the alphabet in use. */
	final boolean isFullyDefined(byte symbol){
		return symbol>=0 && symbolToNumber[symbol]>=0;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Has this class encountered errors while processing? */
	public boolean errorState;

	/** Stores JSON output */
	private final JsonObject jsonStats;

	/** Total number of input reads processed */
	long readsIn;
	/** Total number of input bases processed */
	long basesIn;
	/** Total number of reads written to output */
	long readsOut;
	/** Total number of bases written to output */
	private long basesOut;

	private long readsQTrimmed;
	private long basesQTrimmed;
	private long readsFTrimmed;
	private long basesFTrimmed;
	private long readsQFiltered;
	private long basesQFiltered;
	private long readsEFiltered;
	private long basesEFiltered;
	private long readsNFiltered;
	private long basesNFiltered;

	private long readsPolyTrimmed;
	private long basesPolyTrimmed;

	/** Number of reads modified by kmer trimming */
	private long readsKTrimmed;
	/** Number of bases removed by kmer trimming */
	private long basesKTrimmed;
	/** Number of reads removed by kmer filtering */
	private long readsKFiltered;
	/** Number of bases removed by kmer filtering */
	private long basesKFiltered;

	private long badGcReads;
	private long badGcBases;

	private long badHeaderReads;
	private long badHeaderBases;

	private long readsTrimmedByOverlap;
	private long basesTrimmedByOverlap;

	private long readsTrimmedBySwift;
	private long basesTrimmedBySwift;

	/** Number of unique kmers actually stored in hash tables */
	private long storedKmers;

	/** scaffoldCounts[id] stores the number of reads with kmer matches to that scaffold */
	private AtomicLongArray scaffoldReadCounts;
	/** scaffoldBaseCounts[id] stores the number of bases with kmer matches to that scaffold */
	private AtomicLongArray scaffoldBaseCounts;
	/** scaffoldLengths[id] stores the length of that scaffold */
	private IntList scaffoldLengths;
	/** Set to false to force threads to share atomic counter arrays. */
	private boolean ALLOW_LOCAL_ARRAYS;

	/** Used instead of scaffoldReadCounts if ALLOW_LOCAL_ARRAYS */
	private final long[] scaffoldReadCountsT;
	/** Used instead of scaffoldBaseCounts if ALLOW_LOCAL_ARRAYS */
	private final long[] scaffoldBaseCountsT;


	long readsOutu=0;
	long basesOutu=0;

	long readsOutm=0;
	long basesOutm=0;

	boolean finishedSuccessfully=false;

	private final FlowcellCoordinate flowCoords;

	private final ReadStats readstats;
	private final int[] overlapVector;
	private final int[] countArray;

	private final IntList idList;
	private final IntList countList;

	private final EntropyTracker eTracker;

	private float[] aprob, bprob;

	/*--------------------------------------------------------------*/

	private final BBDukIndexAndLoader index;
	private final BBDukParser parser;
	/** Input read stream */
	private final ConcurrentReadInputStream cris;
	/** Output read streams */
	private final ConcurrentReadOutputStream ros, rosb, ross;

	/** A scaffold's name is stored at scaffoldNames.get(id).
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	private ArrayList<String> scaffoldNames;
	/** Array of reference files from which to load kmers */
	private String[] ref;
	/** Array of literal strings from which to load kmers */
	private final String[] literal;

	/*--------------------------------------------------------------*/
	/*----------------          Immutable           ----------------*/
	/*--------------------------------------------------------------*/

	//TODO: Document
	private final boolean silent;
	private final boolean json;
	
	private final boolean swift; //https://issues.jgi.doe.gov/browse/AUTOQC-2193

	/** For calculating kmer cardinality in input */
	private final CardinalityTracker loglogIn;
	/** For calculating kmer cardinality in output */
	private final CardinalityTracker loglogOut;
	/** Requires (and sets) cardinality tracking.  This is for input kmers. */
	private final String khistIn;
	/** Requires (and sets) cardinality tracking.  This is for output kmers. */
	private final String khistOut;

	/** Input reads */
	private final String in1, in2;
	/** Statistics output files */
	private final String outstats, outrqc, outrpkm, outrefstats, polymerStatsFile;

	final boolean tossJunk;

	/** Quit after this many bases written to outm */
	private long maxBasesOutm;
	/** Quit after this many bases written to outu */
	private long maxBasesOutu;
	
	/** Attempt to match kmers shorter than normal k on read ends when doing kTrimming. */
	private final boolean useShortKmers;
	/** Make the middle base in a kmer a wildcard to improve sensitivity */
	private final boolean maskMiddle;
	private final int midMaskLen;

	/** Search for query kmers with up to this many substitutions */
	private final int qHammingDistance;
	/** Search for short query kmers with up to this many substitutions */
	private final int qHammingDistance2;

	/** Trim this much extra around matched kmers */
	private final int trimPad;

	/*--------------------------------------------------------------*/
	/*----------------      Flowcell Filtering      ----------------*/
	/*--------------------------------------------------------------*/

	//TODO: Document
	private final int xMinLoc;
	private final int yMinLoc;
	private final int xMaxLoc;
	private final int yMaxLoc;
	private final boolean locationFilter;

	/*--------------------------------------------------------------*/
	/*----------------       Variant-Related        ----------------*/
	/*--------------------------------------------------------------*/

	//TODO: Document
	private VarMap varMap;
	private ScafMap scafMap;
	private boolean fixVariants;
	private final boolean unfixVariants;

	/** Filter reads with unsupported substitutions */
	private final boolean filterVars;
	/** Maximum allowed unsupported substitutions in a read */
	private final int maxBadSubs;
	/** Maximum variant depth for a variant to be considered unsupported */
	private final int maxBadSubAlleleDepth;
	/** Minimum read depth for a variant to be considered unsupported */
	private final int minBadSubReadDepth;
	//TODO
	private final int minBadSubEDist;
	//TODO
	private final float maxBadAlleleFraction;

	/*--------------------------------------------------------------*/
	/*----------------        Entropy Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Minimum entropy to be considered "complex", on a scale of 0-1 */
	private final float entropyCutoff;
	/** Mask entropy with a highpass filter */
	private final boolean entropyHighpass;

	/** Change the quality scores to be proportional to the entropy */
	private final boolean entropyMark;
	/** Mask low-entropy areas (e.g., with N) */
	private final boolean entropyMask;
	/** Trim only trailing or leading low-entropy areas, ignoring middle areas */
	private final int entropyTrim;
	/** Convert low-entropy areas to lower case */
	private final boolean entropyMaskLowercase;

	/** Perform entropy calculation */
	private final boolean calcEntropy;

	/*--------------------------------------------------------------*/
	/*----------------         Homopolymers         ----------------*/
	/*--------------------------------------------------------------*/

	private final boolean countPolymers;
	private final byte polymerChar1;
	private final byte polymerChar2;
	/** Minimum length to consider a homopolymer, for the purpose of statistics */
	private final int polymerLength;

	/** Tracks homopolymer statistics */
	private final PolymerTracker pTracker;

	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/

	/** Don't look for kmers in read 1 */
	private final boolean skipR1;
	/** Don't look for kmers in read 2 */
	private final boolean skipR2;
	/** Correct errors via read overlap */
	private final boolean ecc;
	/** True if a ReadStats object is being used for collecting data */
	private final boolean makeReadStats;

	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	/** Don't allow a read 'N' to match a reference 'A'.
	 * Reduces sensitivity when hdist>0 or edist>0.  Default: false. */
	private final boolean forbidNs;
	/** AND bitmask with 0's at the middle base */
	private final long middleMask;

	/** Normal kmer length */
	private final int k;
	/** k-1; used in some expressions */
	private final int k2;
	/** Emulated kmer greater than k */
	private final int kbig;
	/** Effective kmer size */
	private final int keff;
	/** Shortest kmer to use for trimming */
	private final int mink;
	/** A read may contain up to this many kmers before being considered a match.  Default: 0 */
	private final int maxBadKmers0;
	/** A read must share at least this fraction of its kmers to be considered a match.  Default: 0 */
	private final float minKmerFraction;
	/** Reference kmers must cover at least this fraction of read bases to be considered a match.  Default: 0 */
	private final float minCoveredFraction;

	/** Recalibrate quality scores using matrices */
	private final boolean recalibrateQuality;
	/** Quantize quality scores to reduce file size */
	boolean quantizeQuality=false;
	/** Quality-trim the left side */
	private final boolean qtrimLeft;
	/** Quality-trim the right side */
	private final boolean qtrimRight;
	/** Trim soft-clipped bases */
	private final boolean trimClip;
	/** Trim poly-A tails of at least this length */
	private final int trimPolyA;

	/** Trim poly-G prefixes of at least this length */
	private final int trimPolyGLeft;
	/** Trim poly-G tails of at least this length */
	private final int trimPolyGRight;
	/** Remove reads with poly-G prefixes of at least this length */
	private final int filterPolyG;
	/** Allow this many consecutive mismatching symbols in the homopolymer */
	int maxNonPoly=1;

	/** Trim poly-C prefixes of at least this length */
	private final int trimPolyCLeft;
	/** Trim poly-C tails of at least this length */
	private final int trimPolyCRight;
	/** Remove reads with poly-C prefixes of at least this length */
	private final int filterPolyC;

	/** Trim bases at this quality or below.  Default: 4 */
	private final float trimq;
	/** Error rate for trimming (derived from trimq) */
	private final float trimE;
	/** Throw away reads below this average quality after trimming.  Default: 0 */
	private final float minAvgQuality;
	/** Throw away reads with any base below this quality after trimming.  Default: 0 */
	private final byte minBaseQuality;
	/** If positive, calculate average quality from the first X bases only.  Default: 0 */
	private final int minAvgQualityBases;
	/** Throw away reads failing chastity filter (:Y: in read header) */
	private final boolean chastityFilter;
	/** Crash if a barcode is encountered that contains Ns or is not in the table */
	private final boolean failBadBarcodes;
	/** Remove reads with Ns in barcodes or that are not in the table */
	private final boolean removeBadBarcodes;
	/** Fail reads missing a barcode */
	private final boolean failIfNoBarcode;
	/** A set of valid barcodes; null if unused */
	private final HashSet<String> barcodes;
	/** Throw away reads containing more than this many Ns.  Default: -1 (disabled) */
	private final int maxNs;
	/** Throw away reads containing without at least this many consecutive called bases. */
	private final int minConsecutiveBases;
	/** Throw away reads containing fewer than this fraction of any particular base. */
	private final float minBaseFrequency;
	/** Throw away reads shorter than this after trimming.  Default: 10 */
	private final int minReadLength;
	/** Throw away reads longer than this after trimming.  Default: Integer.MAX_VALUE */
	private final int maxReadLength;
	/** Toss reads shorter than this fraction of initial length, after trimming */
	private final float minLenFraction;
	/** Filter reads by whether or not they have matching kmers */
	private final boolean kfilter;
	/** Trim matching kmers and all bases to the left */
	private final boolean ktrimLeft;
	/** Trim matching kmers and all bases to the right */
	boolean ktrimRight;
	/** Don't trim, but replace matching kmers with a symbol (default N) */
	private final boolean ktrimN;
	/** Exclude kmer itself when ktrimming */
	private final boolean ktrimExclusive;
	/** Split into two reads around the kmer */
	private final boolean ksplit;
	/** Replace bases covered by matched kmers with this symbol */
	private final byte trimSymbol;
	/** Convert kmer-masked bases to lowercase */
	private final boolean kmaskLowercase;
	/** Only mask fully-covered bases **/
	private final boolean kmaskFullyCovered;
	/** Output over-trimmed reads to outbad (outmatch).  If false, they are discarded. */
	private final boolean addTrimmedToBad;
	/** Find the sequence that shares the most kmer matches when filtering. */
	private final boolean findBestMatch;
	/** Trim pairs to the same length, when adapter-trimming */
	private final boolean trimPairsEvenly;
	/** Trim left bases of the read to this position (exclusive, 0-based) */
	private final int forceTrimLeft;
	/** Trim right bases of the read after this position (exclusive, 0-based) */
	private final int forceTrimRight;
	/** Trim this many rightmost bases of the read */
	private final int forceTrimRight2;
	/** Trim right bases of the read modulo this value.
	 * e.g. forceTrimModulo=50 would trim the last 3bp from a 153bp read. */
	private final int forceTrimModulo;

	/** Discard reads with GC below this. */
	private final float minGC;
	/** Discard reads with GC above this. */
	private final float maxGC;
	/** Discard reads outside of GC bounds. */
	private final boolean filterGC;
	/** Average GC for paired reads. */
	private final boolean usePairGC;

	/** If positive, only look for kmer matches in the leftmost X bases */
	private final int restrictLeft;
	/** If positive, only look for kmer matches the rightmost X bases */
	private final int restrictRight;

	/** Pairs go to outbad if either of them is bad, as opposed to requiring both to be bad.
	 * Default: true. */
	private final boolean removePairsIfEitherBad;

	/** Rather than discarding, trim failures to 1bp.
	 * Default: false. */
	private final boolean trimFailuresTo1bp;

	/** Print only statistics for scaffolds that matched at least one read
	 * Default: true. */
	private final boolean printNonZeroOnly;

	/** Rename reads to indicate what they matched.
	 * Default: false. */
	private final boolean rename;

	/** Fraction of kmers to skip, 0 to 16 out of 17 */
	private final int speed;

	/** noAccel is true if speed and qSkip are disabled, accel is the opposite. */
	private final boolean noAccel;
	private final boolean accel;

	private final boolean pairedToSingle;

	/*--------------------------------------------------------------*/
	/*-----------        Symbol-Specific Constants        ----------*/
	/*--------------------------------------------------------------*/

	/** True for amino acid data, false for nucleotide data */
	private final boolean amino;
	private final int bitsPerBase;

	private final int minlen;
	private final int minminlen;
	/** The length of half of a kmer outside the middle mask */
	private final int minlen2;
	//	private final int shift;
	private final int shift2;
	private final long mask;
	private final long kmask;

	/** x|kMasks[i] will set the bit to the left of the leftmost base */
	private final long[] lengthMasks;

	/** Symbol code; -1 for undefined */
	private final byte[] symbolToNumber;
	/** Symbol code; 0 for undefined */
	private final byte[] symbolToNumber0;
	/** Complementary symbol code; 0 for undefined */
	private final byte[] symbolToComplementNumber0;

	/*--------------------------------------------------------------*/
	/*----------------         BBMerge Flags        ----------------*/
	/*--------------------------------------------------------------*/

	/** Trim implied adapters based on overlap, for reads with insert size shorter than read length */
	private final boolean trimByOverlap;
	private final boolean useQualityForOverlap;
	//	private final boolean strictOverlap;

	private final int minOverlap0;
	private final int minOverlap;
	private final int minInsert0;
	private final int minInsert;

	private final float maxRatio;
	private final float ratioMargin;
	private final float ratioOffset;
	private final float efilterRatio;
	private final float efilterOffset;
	private final float pfilterRatio;
	private final float meeFilter;

	/*--------------------------------------------------------------*/
	/*----------------         Side Channel         ----------------*/
	/*--------------------------------------------------------------*/

	private SideChannel3 sidechannel;

	/*--------------------------------------------------------------*/
	/*----------------        Histogram Flags       ----------------*/
	/*--------------------------------------------------------------*/

	/** 
	 * Generate histograms from the reads before rather than after processing;
	 * default is true.  Khist is handled independently.
	 */
	private final boolean histogramsBeforeProcessing;

	private final boolean MAKE_IHIST;
	
	/** Number of ProcessThreads */
	private final int THREADS;

	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Verbose messages */
	public static final boolean verbose=false;
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Print speed statistics upon completion */
	public static boolean showSpeed=true;
	/** Number of columns for statistics output, 3 or 5 */
	public static int STATS_COLUMNS=3;
	/** Release memory used by kmer storage after processing reads */
	public static boolean RELEASE_TABLES=true;
	/** Make unambiguous copies of ref sequences with ambiguous bases */
	public static boolean REPLICATE_AMBIGUOUS=false;

	/** Stores some data for statistics when running RQCFilter; not used otherwise. */
	public static HashMap<String, Long> RQC_MAP=null;

}
