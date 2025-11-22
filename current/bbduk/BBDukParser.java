package bbduk;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import cardinality.CardinalityTracker;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import json.JsonObject;
import kmer.AbstractKmerTableSet;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Tools;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.SamLine;
import structures.Quantizer;
import tracker.EntropyTracker;
import tracker.ReadStats;
import var2.ScafMap;
import var2.VarMap;

/**
 * Handles parsing for BBDuk
 * @author Brian Bushnell
 * @date November 19, 2025
 *
 */
class BBDukParser {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	BBDukParser(String[] args, Class<?> c){
		args=preParse(args, c);
		runParseLoop(args);
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			samplerate=parser.samplerate;
			sampleseed=parser.sampleseed;
			recalibrateQuality=parser.recalibrateQuality;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			loglog=parser.loglog;
			loglogOut=parser.loglogOut;
			forceTrimModulo=parser.forceTrimModulo;
			forceTrimLeft=parser.forceTrimLeft;
			forceTrimRight=parser.forceTrimRight;
			forceTrimRight2=parser.forceTrimRight2;
			trimq=parser.trimq;
			trimE=parser.trimE();
			qtrimLeft=parser.qtrimLeft && trimE<1;
			qtrimRight=parser.qtrimRight && trimE<1;
			trimClip=parser.trimClip;
			trimPolyA=parser.trimPolyA;
			trimPolyGLeft=parser.trimPolyGLeft;
			trimPolyGRight=parser.trimPolyGRight;
			filterPolyG=parser.filterPolyG;
			trimPolyCLeft=parser.trimPolyCLeft;
			trimPolyCRight=parser.trimPolyCRight;
			filterPolyC=parser.filterPolyC;
			maxNonPoly=parser.maxNonPoly;
			minLenFraction=parser.minLenFraction;
			minAvgQuality=parser.minAvgQuality;
			minAvgQualityBases=parser.minAvgQualityBases;
			minBaseQuality=parser.minBaseQuality;
			chastityFilter=parser.chastityFilter;
			failBadBarcodes=parser.failBadBarcodes;
			removeBadBarcodes=parser.removeBadBarcodes;
			failIfNoBarcode=parser.failIfNoBarcode;
			barcodes=parser.barcodes;
			minReadLength=parser.minReadLength;
			maxReadLength=parser.maxReadLength;
			maxNs=parser.maxNs;
			minConsecutiveBases=parser.minConsecutiveBases;
			removePairsIfEitherBad=(!parser.requireBothBad) && (!trimFailuresTo1bp);
			tossJunk=parser.tossJunk;

			minGC=parser.minGC;
			maxGC=parser.maxGC;
			filterGC=(minGC>0 || maxGC<1);
			usePairGC=parser.usePairGC;
		}
		
		adjustFiles();
		
		if(WAYS!=7){
			System.err.println("WARNING! WAYS="+WAYS+".  Probably debug mode.");
		}
		
		if(khistIn!=null || khistOut!=null){
			if(khistIn!=null){parser.loglog=true;}
			if(khistOut!=null){parser.loglogOut=true;}
			parser.loglogbuckets=Tools.max(parser.loglogbuckets, 128000);
			CardinalityTracker.trackCounts=true;
		}
		
		if(!ordered && !silent){
			outstream.println("Set ORDERED to "+ordered);
			if(jsonStats!=null){jsonStats.add("ordered", ordered);}
		}
		if(silent || json){
			AbstractKmerTableSet.DISPLAY_PROGRESS=false; //TODO: Test to make sure this occurs for silent mode 
		}
		
		if(hammingDistance2==-1){hammingDistance2=hammingDistance;}
		if(qHammingDistance2==-1){qHammingDistance2=qHammingDistance;}
		if(editDistance2==-1){editDistance2=editDistance;}
		
		if(minoverlap_>=0){
			minOverlap=Tools.max(minoverlap_, 1);
			minOverlap0=Tools.min(minOverlap0, minOverlap);
		}
		
		if(mininsert_>=0){
			minInsert=Tools.max(mininsert_, 1);
			minInsert0=Tools.min(minInsert0, minInsert);
		}
		
		/* Set final variables; post-process and validate argument combinations */
		
		hammingDistance=Tools.max(editDistance, hammingDistance);
		hammingDistance2=Tools.max(editDistance2, hammingDistance2);
		minSkip=Tools.max(1, Tools.min(minSkip, maxSkip));
		maxSkip=Tools.max(minSkip, maxSkip);
		forbidNs=(forbidNs || hammingDistance<1);
		restrictLeft=Tools.max(restrictLeft, 0);
		restrictRight=Tools.max(restrictRight, 0);
		findBestMatch=(rename || findBestMatch);
		noAccel=(speed<1 && qSkip<2);
		accel=!noAccel;
		locationFilter=(xMinLoc>0 || yMinLoc>0 || xMaxLoc>-1 || yMaxLoc>-1);
		
		amino=Shared.AMINO_IN;
		rcomp=rcomp && !amino;
		
		//Set K
		maxSupportedK=(amino ? 12 : 31);
		if(!setk){k=(amino ? 11 : 27);}
		kbig=(k>maxSupportedK ? k : -1);
		k=Tools.min(k, maxSupportedK);
		
		if(strictOverlap){
			maxRatio=0.05f;
			ratioMargin=9f;
			ratioOffset=0.5f;
			efilterRatio=3.5f;
			efilterOffset=0.05f;
			pfilterRatio=0.001f;
			meeFilter=15f;
		}else{
			maxRatio=0.10f;
			ratioMargin=5f;
			ratioOffset=0.4f;
			efilterRatio=6f;
			efilterOffset=0.05f;
			pfilterRatio=0.00005f;
			meeFilter=999999999;
		}
		MAKE_IHIST=ReadStats.COLLECT_INSERT_STATS;
		
		{
			long usableMemory;
			long tableMemory;

			{
				long memory=Runtime.getRuntime().maxMemory();
				double xmsRatio=Shared.xmsRatio();
				usableMemory=(long)Tools.max(((memory-96000000-(20*400000 /* for atomic arrays */))*(xmsRatio>0.97 ? 0.82 : 0.72)), memory*0.45);
				tableMemory=(long)(usableMemory*.95);
			}

			if(initialSize<1){
				final long memOverWays=tableMemory/(12*WAYS);
				final double mem2=(prealloc ? preallocFraction : 1)*tableMemory;
				initialSize=(prealloc || memOverWays<initialSizeDefault ? (int)Tools.min(2142000000, (long)(mem2/(12*WAYS))) : initialSizeDefault);
				if(initialSize!=initialSizeDefault && !silent){
					outstream.println("Initial size set to "+initialSize);
				}
			}
		}
		
		if(ktrimLeft || ktrimRight || ktrimN || ksplit){
			if(kbig>k){
				outstream.println("***********************   WARNING   ***********************");
				outstream.println("WARNING: When kmer-trimming or masking, the maximum value of K is "+k+".");
				outstream.println("K has been reduced from "+kbig+" to "+k+".");
				outstream.println("***********************************************************");
				kbig=k;
			}
		}
		
		if((speed>0 || qSkip>1) && kbig>k){
			outstream.println("***********************   WARNING   ***********************");
			outstream.println("WARNING: When speed>0 or qskip>1, the maximum value of K is "+k+".");
			outstream.println("K has been reduced from "+kbig+" to "+k+".");
			outstream.println("***********************************************************");
			kbig=k;
		}
		
		if((speed>0 && qSkip>1) || (qSkip>1 && maxSkip>1) || (speed>0 && maxSkip>1)){
			outstream.println("WARNING: It is not recommended to use more than one of 'qskip', 'speed', and 'rskip/maxskip' together.");
			outstream.println("qskip="+qSkip+", speed="+speed+", maxskip="+maxSkip);
		}
		
		k2=k-1;
		keff=Tools.max(k, kbig);
		if(maskMiddle) {
			midMaskLen=(midMaskLen>0 ? midMaskLen : 2-(k&1));
		}else{
			midMaskLen=0;
		}
		if(kbig>k){
			minSkip=maxSkip=0;
			if(maskMiddle){
				outstream.println("maskMiddle was disabled because kbig>k");
				maskMiddle=false;
				midMaskLen=0;
			}
		}
		mink=Tools.min(mink, k);
		
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
			lengthMasks=new long[symbolArrayLen];
			setMasks=new long[symbols][symbolArrayLen];
			for(int i=0; i<symbolArrayLen; i++){
				clearMasks[i]=~(symbolMask<<(bitsPerBase*i));
				leftMasks[i]=((-1L)<<(bitsPerBase*i));
				rightMasks[i]=~((-1L)<<(bitsPerBase*i));
				lengthMasks[i]=((1L)<<(bitsPerBase*i));
				for(long j=0; j<symbols; j++){
					setMasks[(int)j][i]=(j<<(bitsPerBase*i));
				}
			}
			
			minlen=k-1;
			minminlen=mink-1;
			minlen2=(maskMiddle ? (k-midMaskLen)/2 : k);
			shift=bitsPerBase*k;
			shift2=shift-bitsPerBase;
			mask=(shift>63 ? -1L : ~((-1L)<<shift));
			kmask=lengthMasks[k];
		}
		
		minKmerFraction=Tools.max(minKmerFraction, 0);
		assert(minKmerFraction<=1) : "minKmerFraction must range from 0 to 1; value="+minKmerFraction;
		
		minCoveredFraction=Tools.max(minCoveredFraction, 0);
		assert(minCoveredFraction<=1) : "minCoveredFraction must range from 0 to 1; value="+minCoveredFraction;
		
		if(mink>0 && mink<k){useShortKmers=true;}
		if(useShortKmers){
			if(maskMiddle){
				outstream.println("maskMiddle was disabled because useShortKmers=true");
				maskMiddle=false;
				midMaskLen=0;
			}
		}
		
		kfilter=(ref!=null || literal!=null) && !(ktrimRight || ktrimLeft || ktrimN || ksplit);
		assert(findBestMatch==false || kfilter==false || kbig<=k) : "K must be less than 32 in 'findBestMatch' mode";
		
		assert(!useShortKmers || ktrimRight || ktrimLeft || ktrimN || ksplit) : "\nSetting mink ("+mink+") or useShortKmers also requires setting a ktrim mode, such as 'r', 'l', or 'n'\n";
		
		if(maskMiddle){
			assert(k>midMaskLen+1);
			int bits=midMaskLen*bitsPerBase;
//			int shift=(k-maskMiddle)&(~1);//Equivalent to (x/2)*2
			int shift=((k-midMaskLen)/2)*bitsPerBase; //old behavior before moving to variable width can be restored with +1: "((k-maskMiddle+1)/2)*bitsPerBase"
			middleMask=~((~((-1L)<<bits))<<shift);
//			mask<<=(shift);
		}else{
			middleMask=-1L;
		}
		
		makeReadStats=ReadStats.collectingStats();
		
		if(recalibrateQuality || true){
			SamLine.SET_FROM_OK=true;//TODO:  Should ban other operations
		}
		
		align=(align || alignRef!=null);
		
		//Initialize entropy
		calcEntropy=(entropyCutoff>=0 || entropyMark);
		if(calcEntropy){
			assert(EntropyTracker.defaultWindowBases>0 && (entropyMark || (entropyCutoff>=0 && entropyCutoff<=1)));
		}
		assert(calcEntropy || (!entropyMask && !entropyMark && entropyTrim==0)) : "Entropy masking/trimming operations require the entropy flag to be set.";
		
		if(polymerStatsFile!=null || (polymerChar1>=0 && polymerChar2>=0)){
			countPolymers=true;
		}
	}
	
	private boolean adjustFiles() {
		
		ref=modifyRefPath(ref, refNames);
		altref=modifyRefPath(altref, altRefNames);
		if(literal!=null){
			refNames.add("literal");
			if(!altRefNames.isEmpty()){altRefNames.add("literal");}
		}

		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(in1!=null && in1.contains("#") && !new File(in1).exists()){
			int pound=in1.lastIndexOf('#');
			String a=in1.substring(0, pound);
			String b=in1.substring(pound+1);
			in1=a+1+b;
			in2=a+2+b;
		}
		if(in2!=null && (in2.contains("=") || in2.equalsIgnoreCase("null"))){in2=null;}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED && !silent){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		if(qfin1!=null && qfin1.contains("#") && in2!=null && !new File(qfin1).exists()){
			int pound=qfin1.lastIndexOf('#');
			String a=qfin1.substring(0, pound);
			String b=qfin1.substring(pound+1);
			qfin1=a+1+b;
			qfin2=a+2+b;
		}
		
		if(out1!=null && out1.contains("#")){
			int pound=out1.lastIndexOf('#');
			String a=out1.substring(0, pound);
			String b=out1.substring(pound+1);
			out1=a+1+b;
			out2=a+2+b;
		}
		
		if(qfout1!=null && qfout1.contains("#") && in2!=null && !new File(qfout1).exists()){
			int pound=qfout1.lastIndexOf('#');
			String a=qfout1.substring(0, pound);
			String b=qfout1.substring(pound+1);
			qfout1=a+1+b;
			qfout2=a+2+b;
		}
		
		if(outb1!=null && outb1.contains("#")){
			int pound=outb1.lastIndexOf('#');
			String a=outb1.substring(0, pound);
			String b=outb1.substring(pound+1);
			outb1=a+1+b;
			outb2=a+2+b;
		}

		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		qfin1=Tools.fixExtension(qfin1);
		qfin2=Tools.fixExtension(qfin2);
		ref=Tools.fixExtension(ref);
		
		if((out2!=null || outb2!=null) && (in1!=null && in2==null)){
			if(!FASTQ.FORCE_INTERLEAVED){outstream.println("Forcing interleaved input because paired output was specified for a single input file.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=true;
		}

		if(!setOut){
			if(!silent && !json && alignOut==null){
				outstream.println("No output stream specified.  To write to stdout, please specify 'out=stdout.fq' or similar.");
			}
			out1=out2=null;
		}else if("stdout".equalsIgnoreCase(out1) || "standarddout".equalsIgnoreCase(out1)){
			out1="stdout.fq";
			out2=null;
		}
		return true;
	}
	
	private String[] preParse(String[] args, Class<?> c) {

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, c, true);
			args=pp.args;
			outstream=pp.outstream;
			jsonStats=pp.jsonObject;
			json=pp.json;
		}
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		parser=new Parser();
		parser.trimq=6;
		parser.minAvgQuality=0;
		parser.minReadLength=10;
		parser.maxReadLength=Integer.MAX_VALUE;
		parser.minLenFraction=0f;
		parser.requireBothBad=false;
		parser.maxNs=-1;
		parser.overwrite=overwrite;
		
		return args;
	}
	
	private void runParseLoop(String[] args) {
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("threads") || a.equals("t")){
				int x=("auto".equals(b) ? Shared.LOGICAL_PROCESSORS : Integer.parseInt(b));
				Shared.setThreads(x);
				THREADS=x;
			}else if(a.equals("workers") || a.equals("wt") || a.equals("w") || a.equals("workerthreads")){
				int x=("auto".equals(b) ? Shared.LOGICAL_PROCESSORS : Integer.parseInt(b));
				THREADS=x;
			}else if(a.equals("threadsin") || a.equals("tin")){
				threadsIn=Integer.parseInt(b);
			}else if(a.equals("threadsout") || a.equals("tout")){
				threadsOut=Integer.parseInt(b);
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfin1=b;
			}else if(a.equals("qfin2")){
				qfin2=b;
			}else if(a.equals("qfout") || a.equals("qfout1")){
				qfout1=b;
			}else if(a.equals("qfin2")){
				qfout2=b;
			}else if(a.equals("out") || a.equals("out1") || a.equals("outu") || a.equals("outu1") || a.equals("outnonmatch") ||
					a.equals("outnonmatch1") || a.equals("outunnmatch") || a.equals("outunmatch1") || a.equals("outunnmatched") || a.equals("outunmatched1")){
				out1=b;
				setOut=true;
			}else if(a.equals("out2") || a.equals("outu2") || a.equals("outnonmatch2") || a.equals("outunmatch2") ||
					a.equals("outnonmatched2") || a.equals("outunmatched2")){
				out2=b;
			}else if(a.equals("outb") || a.equals("outm") || a.equals("outb1") || a.equals("outm1") || a.equals("outbad") ||
					a.equals("outbad1") || a.equals("outmatch") || a.equals("outmatch1")){
				outb1=b;
				setOut=true;
			}else if(a.equals("outb2") || a.equals("outm2") || a.equals("outbad2") || a.equals("outmatch2")){
				outb2=b;
			}else if(a.equals("outs") || a.equals("outsingle")){
				outsingle=b;
			}else if(a.equals("stats") || a.equals("scafstats")){
				outstats=b;
			}else if(a.equals("polymerstats") || a.equals("polymerstatsfile") || a.equals("pstats") || a.equals("phist")){
				polymerStatsFile=b;
			}else if(a.equals("refstats")){
				outrefstats=b;
			}else if(a.equals("rpkm") || a.equals("fpkm") || a.equals("cov") || a.equals("coverage")){
				outrpkm=b;
			}else if(a.equals("sam") || a.equals("bam")){
				samFile=b;
			}else if(a.equals("duk") || a.equals("outduk")){
				outduk=b;
			}else if(a.equals("rqc")){
				outrqc=b;
			}else if(a.equals("ref") || a.equals("adapters")){
				ref=(b==null) ? null : (new File(b).exists() ? new String[] {b} : b.split(","));
			}else if(a.equals("altref")){
				altref=(b==null) ? null : (new File(b).exists() ? new String[] {b} : b.split(","));
			}else if(a.equals("samref") || a.equals("bamref")){
				samref=b;
			}else if(a.equals("literal")){
				literal=processLiteralArg(b);
			}else if(a.equals("forest")){
				useForest=Parse.parseBoolean(b);
				if(useForest){useTable=useArray=false;}
			}else if(a.equals("table")){
				useTable=Parse.parseBoolean(b);
				if(useTable){useForest=useArray=false;}
			}else if(a.equals("array")){
				useArray=Parse.parseBoolean(b);
				if(useArray){useTable=useForest=false;}
			}else if(a.equals("ways")){
//				WAYS=Integer.parseInt(b);
				assert(false) : "WAYS is final.";
			}else if(a.equals("ordered") || a.equals("ord")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("skipr1")){
				skipR1=Parse.parseBoolean(b);
			}else if(a.equals("skipr2")){
				skipR2=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				assert(b!=null) : "\nThe k key needs an integer value greater than 0, such as k=27\n";
				k=Integer.parseInt(b);
				setk=true;
			}else if(a.equals("mink") || a.equals("kmin")){
				mink=Integer.parseInt(b);
				assert(mink<0 || (mink>0 && mink<32)) : "kmin must be between 1 and 31; default is 4, negative numbers disable it.";
			}else if(a.equals("useshortkmers") || a.equals("shortkmers") || a.equals("usk")){
				useShortKmers=Parse.parseBoolean(b);
			}else if(a.equals("trimextra") || a.equals("trimpad") || a.equals("tp")){
				trimPad=Integer.parseInt(b);
			}else if(a.equals("hdist") || a.equals("hammingdistance")){
				hammingDistance=Integer.parseInt(b);
				assert(hammingDistance>=0 && hammingDistance<4) : "hamming distance must be between 0 and 3; default is 0.";
			}else if(a.equals("qhdist") || a.equals("queryhammingdistance")){
				qHammingDistance=Integer.parseInt(b);
				assert(qHammingDistance>=0 && qHammingDistance<4) : "hamming distance must be between 0 and 3; default is 0.";
			}else if(a.equals("edits") || a.equals("edist") || a.equals("editdistance")){
				editDistance=Integer.parseInt(b);
				assert(editDistance>=0 && editDistance<3) : "edit distance must be between 0 and 2; default is 0.\n" +
						"You can bypass this error message with the -da flag, but edist=3 at K=31 " +
						"requires 15,000,000x the time and memory for indexing compared to edist=0.";
			}else if(a.equals("hdist2") || a.equals("hammingdistance2")){
				hammingDistance2=Integer.parseInt(b);
				assert(hammingDistance2>=0 && hammingDistance2<4) : "hamming distance must be between 0 and 3; default is 0.";
			}else if(a.equals("qhdist2") || a.equals("queryhammingdistance2")){
				qHammingDistance2=Integer.parseInt(b);
				assert(qHammingDistance2>=0 && qHammingDistance2<4) : "hamming distance must be between 0 and 3; default is 0.";
			}else if(a.equals("edits2") || a.equals("edist2") || a.equals("editdistance2")){
				editDistance2=Integer.parseInt(b);
				assert(editDistance2>=0 && editDistance2<3) : "edit distance must be between 0 and 2; default is 0.";
			}else if(a.equals("maxskip") || a.equals("maxrskip") || a.equals("mxs")){
				maxSkip=Integer.parseInt(b);
			}else if(a.equals("minskip") || a.equals("minrskip") || a.equals("mns")){
				minSkip=Integer.parseInt(b);
			}else if(a.equals("skip") || a.equals("refskip") || a.equals("rskip")){
				minSkip=maxSkip=Integer.parseInt(b);
			}else if(a.equals("qskip")){
				qSkip=Integer.parseInt(b);
			}else if(a.equals("speed")){
				speed=Integer.parseInt(b);
				assert(speed>=0 && speed<=16) : "Speed range is 0 to 16.  Value: "+speed;
			}else if(a.equals("skipreads")){
				skipreads=Parse.parseKMG(b);
			}else if(a.equals("maxbadkmers") || a.equals("mbk")){
				maxBadKmers0=Integer.parseInt(b);
			}else if(a.equals("minhits") || a.equals("minkmerhits") || a.equals("mkh")){
				maxBadKmers0=Integer.parseInt(b)-1;
			}else if(a.equals("minkmerfraction") || a.equals("minfraction") || a.equals("mkf")){
				minKmerFraction=Float.parseFloat(b);
			}else if(a.equals("mincoveredfraction") || a.equals("mincovfraction") || a.equals("mcf")){
				minCoveredFraction=Float.parseFloat(b);
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
//				assert(WAYS>1) : "WAYS=1 is for debug mode.";
			}else if(a.equals("mm") || a.equals("maskmiddle")){
				if(b==null || Tools.startsWithLetter(b)) {
					maskMiddle=Parse.parseBoolean(b);
				}else{
					midMaskLen=Integer.parseInt(b);
					maskMiddle=midMaskLen>0;
				}
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("forbidns") || a.equals("forbidn") || a.equals("fn")){
				forbidNs=Parse.parseBoolean(b);
			}else if(a.equals("findbestmatch") || a.equals("fbm")){
				findBestMatch=Parse.parseBoolean(b);
			}else if(a.equals("kfilter")){
				boolean x=Parse.parseBoolean(b);
				if(x){ktrimLeft=ktrimRight=ktrimN=ksplit=false;}
			}else if(a.equals("ksplit")){
				boolean x=Parse.parseBoolean(b);
				if(x){ksplit=true; ktrimLeft=ktrimRight=ktrimN=false;}
				else{ksplit=false;}
			}else if(a.equals("ktrim")){
				if(b==null){b="";}
				if(b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr") || b.equalsIgnoreCase("tips")){
					ktrimLeft=ktrimRight=true;
					ktrimN=ksplit=false;
				}else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){
					ktrimLeft=true;
					ktrimRight=false;
					ktrimN=ksplit=false;
				}else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){
					ktrimLeft=false;
					ktrimRight=true;
					ktrimN=ksplit=false;
				}else if(b.equalsIgnoreCase("n")){
					ktrimLeft=ktrimRight=ksplit=false;
					ktrimN=true;
				}else if(b.length()==1 && !b.equalsIgnoreCase("t") && !b.equalsIgnoreCase("f")){
					ktrimLeft=ktrimRight=ksplit=false;
					ktrimN=true;
					trimSymbol=(byte)b.charAt(0);
				}else{
					assert(b!=null && (b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false"))) :
						"\nInvalid setting for ktrim - values must be f (false), l (left), r (right), rl (tips), or n.";
					ktrimRight=ktrimLeft=false;
				}
			}else if(a.equals("trimtips") || a.equals("ktrimtips")){
				if(b==null){b="";}
				if(b.length()>0){
					ktrimLeft=ktrimRight=true;
					ktrimN=ksplit=false;
					restrictLeft=restrictRight=Integer.parseInt(b);
				}
			}else if(a.equals("kmask") || a.equals("mask")){
				if("lc".equalsIgnoreCase(b) || "lowercase".equalsIgnoreCase(b)){
					kmaskLowercase=true;
					ktrimN=true;
					ktrimLeft=ktrimRight=ksplit=false;
				}else{
					if(Parse.parseBoolean(b)){b="N";}
					if(b!=null && b.length()==1){
						ktrimLeft=false;ktrimRight=false;ktrimN=true;
						trimSymbol=(byte)b.charAt(0);
					}else{
						boolean x=Parse.parseBoolean(b);
//						assert(!x) : "\nInvalid setting for kmask - values must be f (false), t (true), or a single character for replacement.";
						ktrimN=x;
					}
				}
			}else if(a.equals("kmaskfullycovered") || a.equals("maskfullycovered") || a.equals("mfc")){
				kmaskFullyCovered=Parse.parseBoolean(b);
			}else if(a.equals("ktrimright")){
				ktrimRight=Parse.parseBoolean(b);
				ktrimLeft=ktrimN=!(ktrimRight);
			}else if(a.equals("ktrimleft")){
				ktrimLeft=Parse.parseBoolean(b);
				ktrimRight=ktrimN=!(ktrimLeft);
			}else if(a.equals("ktrimn")){
				ktrimN=Parse.parseBoolean(b);
				ktrimLeft=ktrimRight=!(ktrimN);
			}else if(a.equals("ktrimexclusive")){
				ktrimExclusive=Parse.parseBoolean(b);
			}else if(a.equals("tbo") || a.equals("trimbyoverlap")){
				trimByOverlap=Parse.parseBoolean(b);
			}else if(a.equals("strictoverlap")){
				strictOverlap=Parse.parseBoolean(b);
			}else if(a.equals("usequality")){
				useQualityForOverlap=Parse.parseBoolean(b);
			}else if(a.equals("tpe") || a.equals("tbe") || a.equals("trimpairsevenly")){
				trimPairsEvenly=Parse.parseBoolean(b);
			}else if(a.equals("ottm") || a.equals("outputtrimmedtomatch")){
				addTrimmedToBad=Parse.parseBoolean(b);
			}else if(a.equals("minoverlap")){
				minoverlap_=Integer.parseInt(b);
			}else if(a.equals("mininsert")){
				mininsert_=Integer.parseInt(b);
			}else if(a.equals("prealloc") || a.equals("preallocate")){
				if(b==null || b.length()<1 || Character.isLetter(b.charAt(0))){
					prealloc=Parse.parseBoolean(b);
				}else{
					preallocFraction=Tools.max(0, Double.parseDouble(b));
					prealloc=(preallocFraction>0);
				}
			}else if(a.equals("restrictleft")){
				restrictLeft=Integer.parseInt(b);
			}else if(a.equals("restrictright")){
				restrictRight=Integer.parseInt(b);
			}else if(a.equals("statscolumns") || a.equals("columns") || a.equals("cols")){
				STATS_COLUMNS=Integer.parseInt(b);
				assert(STATS_COLUMNS==3 || STATS_COLUMNS==5) : "statscolumns bust be either 3 or 5. Invalid value: "+STATS_COLUMNS;
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				printNonZeroOnly=Parse.parseBoolean(b);
			}else if(a.equals("rename")){
				rename=Parse.parseBoolean(b);
			}else if(a.equals("refnames") || a.equals("userefnames")){
				useRefNames=Parse.parseBoolean(b);
			}else if(a.equals("initialsize")){
				initialSize=Parse.parseIntKMG(b);
			}else if(a.equals("dump")){
				dump=b;
			}else if(a.equals("minentropy") || a.equals("entropy") || a.equals("entropyfilter")){
				entropyCutoff=Float.parseFloat(b);
			}else if(a.equals("verifyentropy")){
				EntropyTracker.verify=Parse.parseBoolean(b);
			}else if(a.equals("entropymask") || a.equals("maskentropy")){
				if(b==null){
					entropyMask=true;
				}else if(b.equalsIgnoreCase("lc") || b.equalsIgnoreCase("lowercase")){
					entropyMask=true;
					entropyMaskLowercase=true;
				}else if(b.equalsIgnoreCase("filter")){
					entropyMask=false;
				}else{
					entropyMask=Parse.parseBoolean(b);
				}
			}else if(a.equals("entropytrim") || a.equals("trimentropy")){
				entropyTrim=parseEnd(b);
			}else if(a.equals("entropymark") || a.equals("markentropy")){
				entropyMark=Parse.parseBoolean(b);
			}
			
			else if(a.equals("countpolymers")){
				countPolymers=Parse.parseBoolean(b);
			}else if(a.equals("polybase1")){
				polymerChar1=(byte)b.charAt(0);
			}else if(a.equals("polybase2")){
				polymerChar2=(byte)b.charAt(0);
			}else if(a.equals("polymerratio") || a.equals("pratio")){
				assert(b!=null);
				b=b.toUpperCase();
				if(b.length()==2){
					polymerChar1=(byte)b.charAt(0);
					polymerChar2=(byte)b.charAt(1);
				}else if(b.length()==3){
					assert(b.charAt(1)==',');
					polymerChar1=(byte)b.charAt(0);
					polymerChar2=(byte)b.charAt(2);
				}else{
					assert(false) : "Format should be pratio=G,C";
				}
				assert(polymerChar1>=0 && polymerChar2>=0);
				assert(AminoAcid.baseToNumberACGTN[polymerChar1]>=0 && AminoAcid.baseToNumberACGTN[polymerChar2]>=0) : "Only ACGTN polymer tracking is possible: "+arg;
			}else if(a.equals("polymerlength") || a.equals("plen")){
				polymerLength=Integer.parseInt(b);
				assert(polymerLength>=1);
			}
			
			else if(a.equals("minbasefrequency")){
				minBaseFrequency=Float.parseFloat(b);
			}else if(a.equals("ecco") || a.equals("ecc")){
				ecc=Parse.parseBoolean(b);
			}else if(a.equals("copyundefined") || a.equals("cu")){
				REPLICATE_AMBIGUOUS=Parse.parseBoolean(b);
			}else if(a.equals("path")){
				Data.setPath(b);
			}else if(a.equals("maxbasesoutm")){
				maxBasesOutm=Parse.parseKMG(b);
			}else if(a.equals("maxbasesoutu") || a.equals("maxbasesout")){
				maxBasesOutu=Parse.parseKMG(b);
			}else if(a.equals("vars") || a.equals("variants") || a.equals("varfile") || a.equals("inv")){
				varFile=b;
			}else if(a.equals("vcf") || a.equals("vcffile")){
				vcfFile=b;
			}else if(a.equals("unfixvars") || a.equals("unfixvariants")){
				unfixVariants=Parse.parseBoolean(b);
			}else if(a.equals("histogramsbefore") || a.equals("histbefore")){
				histogramsBeforeProcessing=Parse.parseBoolean(b);
			}else if(a.equals("histogramsafter") || a.equals("histafter")){
				histogramsBeforeProcessing=!Parse.parseBoolean(b);
			}
			
			else if(a.equals("trimfailures") || a.equals("trimfailuresto1bp")){
				trimFailuresTo1bp=Parse.parseBoolean(b);
			}
			
			else if(a.equals("minx") || a.equals("xmin")){
				xMinLoc=Parse.parseIntKMG(b);
			}else if(a.equals("miny") || a.equals("ymin")){
				yMinLoc=Parse.parseIntKMG(b);
			}else if(a.equals("maxx") || a.equals("xmax")){
				xMaxLoc=Parse.parseIntKMG(b);
			}else if(a.equals("maxy") || a.equals("ymax")){
				yMaxLoc=Parse.parseIntKMG(b);
			}
			
			else if(a.equalsIgnoreCase("pairedToSingle")){
				pairedToSingle=Parse.parseBoolean(b);
			}
			
			else if(a.equals("filtersubs") || a.equals("filtervars")){
				filterVars=Parse.parseBoolean(b);
			}else if(a.equals("maxbadsubs") || a.equals("maxbadbars")){
				maxBadSubs=Integer.parseInt(b);
			}else if(a.equals("maxbadsubdepth") || a.equals("maxbadvardepth") || a.equals("maxbadsuballeledepth") || a.equals("maxbadvaralleledepth") || a.equals("mbsad")){
				maxBadSubAlleleDepth=Integer.parseInt(b);
			}else if(a.equals("minbadsubreaddepth") || a.equals("minbadvarreaddepth") || a.equals("mbsrd")){
				minBadSubReadDepth=Integer.parseInt(b);
			}
			
			else if(a.equals("khist") || a.equals("khistin")){
				khistIn=b;
			}else if(a.equals("khistout")){
				khistOut=b;
			}else if(a.equals("quantize") || a.equals("quantizesticky")){
				quantizeQuality=Quantizer.parse(arg, a, b);
			}
			
			else if(a.equals("json")){
				json=Parse.parseBoolean(b);
			}
			
			else if(a.equals("swift")){
				swift=Parse.parseBoolean(b);
			}
			
			else if(a.equals("alignout") || a.equals("sideout")){
				alignOut=b;
			}else if(a.equals("align")){
				align=Parse.parseBoolean(b);
			}else if(a.equals("alignref") || a.equals("sideref")){
				alignRef=b;
			}else if(a.equals("alignk") || a.equals("sidek") || a.equals("alignk1") || a.equals("sidek1")){
				alignK1=Integer.parseInt(b);
			}else if(a.equals("alignk2") || a.equals("sidek2")){
				alignK2=Integer.parseInt(b);
			}else if(a.equals("alignminid") || a.equals("alignminid1") || a.equals("sideminid") || a.equals("sideminid1")){
				alignMinid1=Float.parseFloat(b);
			}else if(a.equals("alignminid2") || a.equals("sideminid2")){
				alignMinid2=Float.parseFloat(b);
			}else if(a.equals("alignmm1") || a.equals("alignmidmask1") || a.equals("sidemm1") || a.equals("sidemidmask1")){
				alignMM1=Integer.parseInt(b);
			}else if(a.equals("alignmm2") || a.equals("alignmidmask2") || a.equals("sidemm2") || a.equals("sidemidmask2")){
				alignMM2=Integer.parseInt(b);
			}
			
			else if(i==0 && in1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				in1=args[i];
			}else if(i==1 && out1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				out1=args[i];
				setOut=true;
			}else if(i==2 && ref==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				ref=(new File(args[i]).exists() ? new String[] {args[i]} : args[i].split(","));
			}
			
			else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseHist(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseQualityAdjust(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(Parser.parseSam(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(parser.parseCommon(arg, a, b)){
				//do nothing
			}else if(parser.parseCardinality(arg, a, b)){
				//do nothing
			} 
			
			else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
	}
	
	/**
	 * Processes array of reference paths and adds modified paths to list.
	 * @param array Array of reference file paths
	 * @param list Output list to store processed paths
	 * @return Modified array of reference paths
	 */
	String[] modifyRefPath(String[] array, ArrayList<String> list){
		if(array==null){return array;}
		for(String s : array){
			String s2=modifyRefPath(s);
			list.add(s2);
		}
		return list.toArray(new String[0]);
	}
	
	/**
	 * Resolves reference path shortcuts to actual file paths.
	 * Handles special keywords like "phix", "adapters", "truseq", etc.
	 * to their corresponding resource files.
	 *
	 * @param s Reference path or keyword
	 * @return Resolved file path
	 */
	static String modifyRefPath(String s){
		if(s==null || Tools.isReadableFile(s)){
			//do nothing
		}else{
			if("phix".equalsIgnoreCase(s)){
				s=Data.findPath("?phix2.fa.gz");
			}else if("polya".equalsIgnoreCase(s) || "polyt".equalsIgnoreCase(s)){
				s=Data.findPath("?polyA.fa.gz");
			}else if("lambda".equalsIgnoreCase(s)){
				s=Data.findPath("?lambda.fa.gz");
			}else if("kapa".equalsIgnoreCase(s)){
				s=Data.findPath("?kapatags.L40.fa");
			}else if("pjet".equalsIgnoreCase(s)){
				s=Data.findPath("?pJET1.2.fa");
			}else if("mtst".equalsIgnoreCase(s)){
				s=Data.findPath("?mtst.fa");
			}else if("adapters".equalsIgnoreCase(s)){
				s=Data.findPath("?adapters.fa");
			}else if("phixadapters".equalsIgnoreCase(s)){
				s=Data.findPath("?phix_adapters.fa");
			}else if("pacbioadapter".equalsIgnoreCase(s) || "pacbioadapters".equalsIgnoreCase(s)){
				s=Data.findPath("?PacBioAdapter.fa");
			}else if("truseq".equalsIgnoreCase(s)){
				s=Data.findPath("?truseq.fa.gz");
			}else if("nextera".equalsIgnoreCase(s)){
				s=Data.findPath("?nextera.fa.gz");
			}else if("artifacts".equalsIgnoreCase(s)){
				s=Data.findPath("?sequencing_artifacts.fa.gz");
			}else if("crisprs".equalsIgnoreCase(s)){
				s=Data.findPath("?crisprs.fa.gz");
			}else if(s.startsWith("poly") && s.length()==5 && AminoAcid.baseToNumber[s.charAt(4)]>=0) {
				s=Data.findPath("?"+s.toLowerCase()+".fa");
			}else {
				assert(false) : "Can't find reference file "+s;
			}
		}
		return s;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses literal sequence argument into array of sequences.
	 * @param arg Comma-separated list of literal sequences
	 * @return Array of processed literal sequences
	 */
	static final String[] processLiteralArg(String arg) {
//		System.err.println("Caught "+arg);
		if(arg==null) {return null;}
		String[] split=arg.split(",");
		ArrayList<String> list=new ArrayList<String>(split.length);
		for(String b : split) {
			String c=processLiteralTerm(b);
			if(c!=null) {list.add(c);}
		}
		String[] ret=list.isEmpty() ? null : list.toArray(new String[0]);
//		System.err.println("Returning "+Arrays.toString(ret));
		return ret;
	}
	
	/**
	 * Processes individual literal sequence term, expanding polymer shortcuts.
	 * @param b Single literal sequence or polymer specification
	 * @return Processed literal sequence
	 */
	static final String processLiteralTerm(String b) {
//		System.err.println("Parsing "+b);
		assert(b.length()>0) : "Invalid literal sequence: '"+b+"'";
		if(AminoAcid.isACGTN(b)) {
//			System.err.println("Returning valid sequence "+b);
			return b;
		}
		b=b.toUpperCase().replaceAll("-", "");
		if(b.startsWith("POLY")) {
			b=b.replace("POLY", "");
//			System.err.println("Parsing poly "+b);
			assert(AminoAcid.isACGTN(b)) : "Invalid literal sequence: '"+b+"'";
			StringBuilder sb=new StringBuilder(40);
			final int minlen=Tools.max(31+b.length(), b.length()*3);
			while(sb.length()<minlen) {sb.append(b);}
			System.err.println("Adding literal polymer "+sb);
			return sb.toString();
		}else {
			assert(false) : "Invalid literal sequence: '"+b+"'";
		}
//		System.err.println("Returning null");
		return null;
	}
	
	static int parseEnd(String s){
		if(s==null){return RIGHTLEFT;}
		if(s.equalsIgnoreCase("right") || s.equalsIgnoreCase("r")){return RIGHT;}
		if(s.equalsIgnoreCase("left") || s.equalsIgnoreCase("l")){return LEFT;}
		if(s.equalsIgnoreCase("rl") || s.equalsIgnoreCase("lr") || s.equalsIgnoreCase("both") || s.equalsIgnoreCase("b")){return RIGHTLEFT;}
		return Parse.parseBoolean(s) ? RIGHTLEFT : 0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/* Initialize local variables with defaults */
	private boolean setk=false;
	private int mininsert_=-1;
	private int minoverlap_=-1;
	
	/*--------------------------------------------------------------*/
	
	Parser parser=null;

	//TODO: Document
	boolean silent=false;
	boolean json=false;
	
	boolean swift=false; //https://issues.jgi.doe.gov/browse/AUTOQC-2193
	
	/** Requires (and sets) cardinality tracking.  This is for input kmers. */
	String khistIn=null;
	/** Requires (and sets) cardinality tracking.  This is for output kmers. */
	String khistOut=null;
	
	boolean prealloc=false;
	/** Fraction of available memory preallocated to arrays */
	double preallocFraction=1.0;
	/** Initial size of data structures */
	int initialSize=-1;
	
	/** Names of reference files (refNames[0] is valid). */
	ArrayList<String> refNames=new ArrayList<String>();
	final ArrayList<String> altRefNames=new ArrayList<String>();
	/** Set to false to force threads to share atomic counter arrays. */
	boolean ALLOW_LOCAL_ARRAYS=true;
	/** Array of reference files from which to load kmers */
	String[] ref=null;
	/** Alternate reference to be used if main reference has no kmers */
	String[] altref=null;
	/** Array of literal strings from which to load kmers */
	String[] literal=null;
	/** Optional reference for sam file */
	String samref=null;
	
	boolean setOut=false;
	boolean setOutb=false;
	
	int threadsIn=-1;
	int threadsOut=-1;
	
	/** Input reads */
	String in1=null, in2=null;
	/** Input qual files */
	String qfin1=null, qfin2=null;
	/** Output qual files */
	String qfout1=null, qfout2=null;
	/** Output reads (unmatched and at least minlen) */
	String out1=null, out2=null;
	/** Output reads (matched or shorter than minlen) */
	String outb1=null, outb2=null;
	/** Output reads whose mate was discarded */
	String outsingle=null;
	/** Statistics output files */
	String outstats=null, outrqc=null, outrpkm=null, outrefstats=null, polymerStatsFile=null;
	@Deprecated
	/** duk-style statistics */
	String outduk=null;
	
	final boolean tossJunk;
	
	/** Dump kmers here. */
	String dump=null;

	/** Quit after this many bases written to outm */
	long maxBasesOutm=-1;
	/** Quit after this many bases written to outu */
	long maxBasesOutu=-1;
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	long maxReads=-1;
	/** Process this fraction of input reads. */
	float samplerate=1f;
	/** Set samplerate seed to this value. */
	long sampleseed=-1;
	
	/** Output reads in input order.  May reduce speed. */
	boolean ordered=true;
	/** Attempt to match kmers shorter than normal k on read ends when doing kTrimming. */
	boolean useShortKmers=false;
	/** Make the middle base in a kmer a wildcard to improve sensitivity */
	boolean maskMiddle=true;
	int midMaskLen=0;
	
	/** Store reference kmers with up to this many substitutions */
	int hammingDistance=0;
	/** Search for query kmers with up to this many substitutions */
	int qHammingDistance=0;
	/** Store reference kmers with up to this many edits (including indels) */
	int editDistance=0;
	/** Store short reference kmers with up to this many substitutions */
	int hammingDistance2=-1;
	/** Search for short query kmers with up to this many substitutions */
	int qHammingDistance2=-1;
	/** Store short reference kmers with up to this many edits (including indels) */
	int editDistance2=-1;
	/** Never skip more than this many consecutive kmers when hashing reference. */
	int maxSkip=1;
	/** Always skip at least this many consecutive kmers when hashing reference.
	 * 1 means every kmer is used, 2 means every other, etc. */
	int minSkip=1;
	
	/** Trim this much extra around matched kmers */
	int trimPad;
	
	/*--------------------------------------------------------------*/
	/*----------------      Flowcell Filtering      ----------------*/
	/*--------------------------------------------------------------*/

	//TODO: Document
	int xMinLoc=-1;
	int yMinLoc=-1;
	int xMaxLoc=-1;
	int yMaxLoc=-1;
	final boolean locationFilter;
	
	/*--------------------------------------------------------------*/
	/*----------------       Variant-Related        ----------------*/
	/*--------------------------------------------------------------*/

	//TODO: Document
	String varFile=null;
	String vcfFile=null;
	VarMap varMap=null;
//	HashSet<VarKey> varKeySet=null;
	ScafMap scafMap=null;
	boolean fixVariants=false;
	boolean unfixVariants=true;
	
	/** Optional file for quality score recalibration */
	String samFile=null;
	
	/** Filter reads with unsupported substitutions */
	boolean filterVars=false;
	/** Maximum allowed unsupported substitutions in a read */
	int maxBadSubs=2;
	/** Maximum variant depth for a variant to be considered unsupported */
	int maxBadSubAlleleDepth=1;
	/** Minimum read depth for a variant to be considered unsupported */
	int minBadSubReadDepth=2;
	//TODO
	int minBadSubEDist=0;
	//TODO
	float maxBadAlleleFraction=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Entropy Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Minimum entropy to be considered "complex", on a scale of 0-1 */
	float entropyCutoff=-1;
	/** Mask entropy with a highpass filter */
	boolean entropyHighpass=true;
	
	/** Change the quality scores to be proportional to the entropy */
	boolean entropyMark=false;
	/** Mask low-entropy areas (e.g., with N) */
	boolean entropyMask=false;
	/** Trim only trailing or leading low-entropy areas, ignoring middle areas */
	int entropyTrim=0;
	/** Convert low-entropy areas to lower case */
	boolean entropyMaskLowercase=false;
	
	/** Perform entropy calculation */
	final boolean calcEntropy;
	
	
	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Stores JSON output */
	JsonObject jsonStats;
	
	/*--------------------------------------------------------------*/
	/*----------------         Homopolymers         ----------------*/
	/*--------------------------------------------------------------*/
	
	boolean countPolymers=false;
	byte polymerChar1=-1;
	byte polymerChar2=-1;
	/** Minimum length to consider a homopolymer, for the purpose of statistics */
	int polymerLength=20;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Don't look for kmers in read 1 */
	boolean skipR1;
	/** Don't look for kmers in read 2 */
	boolean skipR2;
	/** Correct errors via read overlap */
	boolean ecc;
	/** True if a ReadStats object is being used for collecting data */
	final boolean makeReadStats;
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	boolean rcomp=true;
	/** Don't allow a read 'N' to match a reference 'A'.
	 * Reduces sensitivity when hdist>0 or edist>0.  Default: false. */
	boolean forbidNs;
	/** AND bitmask with 0's at the middle base */
	final long middleMask;
	/** Use HashForest data structure */
	boolean useForest;
	/** Use KmerTable data structure */
	boolean useTable;
	/** Use HashArray data structure (default) */
	boolean useArray=true;
	
	/** Normal kmer length */
	int k=31;
	/** k-1; used in some expressions */
	final int k2;
	/** Emulated kmer greater than k */
	int kbig=-1;
	/** Effective kmer size */
	final int keff;
	/** Shortest kmer to use for trimming */
	int mink=-1;
	/** A read may contain up to this many kmers before being considered a match.  Default: 0 */
	int maxBadKmers0=0;
	/** A read must share at least this fraction of its kmers to be considered a match.  Default: 0 */
	float minKmerFraction=0;
	/** Reference kmers must cover at least this fraction of read bases to be considered a match.  Default: 0 */
	float minCoveredFraction=0;
	
	/** Recalibrate quality scores using matrices */
	final boolean recalibrateQuality;
	/** Quantize quality scores to reduce file size */
	boolean quantizeQuality=false;
	/** Quality-trim the left side */
	final boolean qtrimLeft;
	/** Quality-trim the right side */
	final boolean qtrimRight;
	/** Trim soft-clipped bases */
	final boolean trimClip;
	/** Trim poly-A tails of at least this length */
	final int trimPolyA;
	
	/** Trim poly-G prefixes of at least this length */
	final int trimPolyGLeft;
	/** Trim poly-G tails of at least this length */
	final int trimPolyGRight;
	/** Remove reads with poly-G prefixes of at least this length */
	final int filterPolyG;
	/** Allow this many consecutive mismatching symbols in the homopolymer */
	int maxNonPoly=1;
	
	/** Trim poly-C prefixes of at least this length */
	final int trimPolyCLeft;
	/** Trim poly-C tails of at least this length */
	final int trimPolyCRight;
	/** Remove reads with poly-C prefixes of at least this length */
	final int filterPolyC;
	
	/** Trim bases at this quality or below.  Default: 4 */
	final float trimq;
	/** Error rate for trimming (derived from trimq) */
	final float trimE;
	/** Throw away reads below this average quality after trimming.  Default: 0 */
	final float minAvgQuality;
	/** Throw away reads with any base below this quality after trimming.  Default: 0 */
	final byte minBaseQuality;
	/** If positive, calculate average quality from the first X bases only.  Default: 0 */
	final int minAvgQualityBases;
	/** Throw away reads failing chastity filter (:Y: in read header) */
	final boolean chastityFilter;
	/** Crash if a barcode is encountered that contains Ns or is not in the table */
	final boolean failBadBarcodes;
	/** Remove reads with Ns in barcodes or that are not in the table */
	final boolean removeBadBarcodes;
	/** Fail reads missing a barcode */
	final boolean failIfNoBarcode;
	/** A set of valid barcodes; null if unused */
	final HashSet<String> barcodes;
	/** Throw away reads containing more than this many Ns.  Default: -1 (disabled) */
	final int maxNs;
	/** Throw away reads containing without at least this many consecutive called bases. */
	final int minConsecutiveBases;
	/** Throw away reads containing fewer than this fraction of any particular base. */
	float minBaseFrequency=0;
	/** Throw away reads shorter than this after trimming.  Default: 10 */
	final int minReadLength;
	/** Throw away reads longer than this after trimming.  Default: Integer.MAX_VALUE */
	final int maxReadLength;
	/** Toss reads shorter than this fraction of initial length, after trimming */
	final float minLenFraction;
	/** Filter reads by whether or not they have matching kmers */
	boolean kfilter;
	/** Trim matching kmers and all bases to the left */
	boolean ktrimLeft;
	/** Trim matching kmers and all bases to the right */
	boolean ktrimRight;
	/** Don't trim, but replace matching kmers with a symbol (default N) */
	boolean ktrimN;
	/** Exclude kmer itself when ktrimming */
	boolean ktrimExclusive;
	/** Split into two reads around the kmer */
	boolean ksplit;
	/** Replace bases covered by matched kmers with this symbol */
	byte trimSymbol='N';
	/** Convert kmer-masked bases to lowercase */
	boolean kmaskLowercase;
	/** Only mask fully-covered bases **/
	boolean kmaskFullyCovered;
	/** Output over-trimmed reads to outbad (outmatch).  If false, they are discarded. */
	boolean addTrimmedToBad=true;
	/** Find the sequence that shares the most kmer matches when filtering. */
	boolean findBestMatch;
	/** Trim pairs to the same length, when adapter-trimming */
	boolean trimPairsEvenly;
	/** Trim left bases of the read to this position (exclusive, 0-based) */
	final int forceTrimLeft;
	/** Trim right bases of the read after this position (exclusive, 0-based) */
	final int forceTrimRight;
	/** Trim this many rightmost bases of the read */
	final int forceTrimRight2;
	/** Trim right bases of the read modulo this value.
	 * e.g. forceTrimModulo=50 would trim the last 3bp from a 153bp read. */
	final int forceTrimModulo;
	
	/** Discard reads with GC below this. */
	final float minGC;
	/** Discard reads with GC above this. */
	final float maxGC;
	/** Discard reads outside of GC bounds. */
	final boolean filterGC;
	/** Average GC for paired reads. */
	final boolean usePairGC;
	
	/** If positive, only look for kmer matches in the leftmost X bases */
	int restrictLeft=0;
	/** If positive, only look for kmer matches the rightmost X bases */
	int restrictRight=0;
	
	/** Skip this many initial input reads */
	long skipreads=0;

	/** Pairs go to outbad if either of them is bad, as opposed to requiring both to be bad.
	 * Default: true. */
	final boolean removePairsIfEitherBad;

	/** Rather than discarding, trim failures to 1bp.
	 * Default: false. */
	boolean trimFailuresTo1bp;
	
	/** Print only statistics for scaffolds that matched at least one read
	 * Default: true. */
	boolean printNonZeroOnly=true;
	
	/** Rename reads to indicate what they matched.
	 * Default: false. */
	boolean rename;
	/** Use names of reference files instead of scaffolds.
	 * Default: false. */
	boolean useRefNames;
	
	/** Fraction of kmers to skip, 0 to 16 out of 17 */
	int speed=0;
	
	/** Skip this many kmers when examining the read.  Default 1.
	 * 1 means every kmer is used, 2 means every other, etc. */
	int qSkip=1;
	
	/** noAccel is true if speed and qSkip are disabled, accel is the opposite. */
	final boolean noAccel;
	final boolean accel;
	
	boolean pairedToSingle=false;
	final boolean loglog;
	final boolean loglogOut;
	
	/*--------------------------------------------------------------*/
	/*-----------        Symbol-Specific Constants        ----------*/
	/*--------------------------------------------------------------*/

	/** True for amino acid data, false for nucleotide data */
	final boolean amino;
	final int maxSupportedK;
	final int bitsPerBase;
	final int maxSymbol;
	final int symbols;
	final int symbolArrayLen;
	final int symbolSpace;
	final long symbolMask;
	
	final int minlen;
	final int minminlen;
	/** The length of half of a kmer outside the middle mask */
	final int minlen2;
	final int shift;
	final int shift2;
	final long mask;
	final long kmask;
	
	/** x&clearMasks[i] will clear base i */
	final long[] clearMasks;
	/** x|setMasks[j][i] will set position i to symbol j */
	final long[][] setMasks;
	/** x&leftMasks[i] will clear all bases to the right of i (exclusive) */
	final long[] leftMasks;
	/** x&rightMasks[i] will clear all bases to the left of i (inclusive) */
	final long[] rightMasks;
	/** x|kMasks[i] will set the bit to the left of the leftmost base */
	final long[] lengthMasks;
	
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
	/*----------------         BBMerge Flags        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Trim implied adapters based on overlap, for reads with insert size shorter than read length */
	boolean trimByOverlap;
	boolean useQualityForOverlap;
	boolean strictOverlap=true;
	
	int minOverlap0=7;
	int minOverlap=14;
	int minInsert0=16;
	int minInsert=40;
	
	final float maxRatio;
	final float ratioMargin;
	final float ratioOffset;
	final float efilterRatio;
	final float efilterOffset;
	final float pfilterRatio;
	final float meeFilter;
	
	/*--------------------------------------------------------------*/
	/*----------------         Side Channel         ----------------*/
	/*--------------------------------------------------------------*/
	
	boolean align=false;
	String alignOut=null;
	String alignRef=null;//="phix";
	float alignMinid1=0.66f;
	float alignMinid2=0.56f;
	int alignK1=17; //Phix is unique down to k=13
	int alignK2=13;
	int alignMM1=1;
	int alignMM2=1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Histogram Flags       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Generate histograms from the reads before rather than after processing;
	 * default is true.  Khist is handled independently.
	 */
	boolean histogramsBeforeProcessing=true;
	final boolean MAKE_IHIST;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of tables (and threads, during loading) */
	static final int WAYS=7; //123
	/** Default initial size of data structures */
	static final int initialSizeDefault=128000;

	/** Ends for some operations like entropytrim; could be migrated over to other operations */ 
	static final int RIGHT=1, LEFT=2, RIGHTLEFT=3;
	
	/** Print messages to this stream */
	public static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Print speed statistics upon completion */
	public static boolean showSpeed=true;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.threads();
	/** Number of columns for statistics output, 3 or 5 */
	public static int STATS_COLUMNS=3;
	/** Make unambiguous copies of ref sequences with ambiguous bases */
	public static boolean REPLICATE_AMBIGUOUS=false;
	/** Release memory used by kmer storage after processing reads */
	public static boolean RELEASE_TABLES=true;
	
}
