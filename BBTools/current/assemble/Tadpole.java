package assemble;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import jgi.BBMerge;
import kmer.AbstractKmerTableSet;
import kmer.HashBuffer;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import sort.ContigLengthComparator;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import structures.LongList;
import tracker.ReadStats;
import ukmer.Kmer;
import ukmer.KmerTableSetU;


/**
 * Short-kmer assembler based on KmerCountExact.
 * @author Brian Bushnell
 * @date May 15, 2015
 *
 */
public abstract class Tadpole extends ShaveObject{
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer(), t2=new Timer();
		t.start();
		t2.start();

		final Tadpole x=makeTadpole(args, true);
		t2.stop();
		outstream.println("Initialization Time:      \t"+t2);
		
		///And run it
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	public static Tadpole makeTadpole(String[] args, boolean setDefaults){
		final int k=preparseK(args);
		if(k>31 || FORCE_TADPOLE2){
			synchronized(Tadpole.class){
				AbstractKmerTableSet.MASK_CORE=false;
			}
			return new Tadpole2(args, true);
		}else{
			synchronized(Tadpole.class){
				AbstractKmerTableSet.MASK_CORE=true;
			}
			Tadpole tad=new Tadpole1(args, true);
//			if(!tad.setShave){tad.removeDeadEnds=true;}
//			if(!tad.setRinse){tad.removeBubbles=true;}
			return tad;
		}
	}
	
	public static final int preparseK(String[] args){
		int k=31;
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			while(a.charAt(0)=='-'){a=a.substring(1);}
			
			if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("tad2") || a.equals("tadpole2")){
				FORCE_TADPOLE2=Parse.parseBoolean(b);
			}
		}
		return Kmer.getMult(k)*Kmer.getK(k);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Tadpole(String[] args, boolean setDefaults){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), true);
			args=pp.args;
			outstream=pp.outstream;
		}
		kbig=preparseK(args);
		
		if(setDefaults){
			/* Set global defaults */
			ReadWrite.ZIPLEVEL=2;
			ReadWrite.USE_UNPIGZ=true;
			ReadWrite.USE_PIGZ=true;
			if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
				ByteFile.FORCE_MODE_BF2=true;
			}
			AbstractKmerTableSet.defaultMinprob=0.5;
		}
		
		/* Initialize local variables with defaults */
		Parser parser=new Parser();
		boolean ecc_=false, ecco_=false, merge_=false, testMerge_=true, vstrict_=false, setEcc_=false;
		boolean useOwnership_=false, setUseOwnership_=false;
		
		int prefilter=0;
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("in") || a.equals("in1") || a.equals("ine") || a.equals("ine1")){
				in1.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						in1.add(ss);
					}
				}
			}else if(a.equals("in2") || a.equals("ine2")){
				in2.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						in2.add(ss);
					}
				}
			}else if(a.equals("extra")){
				//Do nothing
			}else if(a.equals("out") || a.equals("out1") || a.equals("oute") || a.equals("oute1")){
				out1.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						out1.add(ss);
					}
				}
			}else if(a.equals("out2") || a.equals("oute2")){
				out2.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						out2.add(ss);
					}
				}
			}else if(a.equals("outd") || a.equals("outdiscard") || a.equals("outd1") || a.equals("outdiscard1")){
				outd1.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						outd1.add(ss);
					}
				}
			}else if(a.equals("outd2") || a.equals("outdiscard2")){
				outd2.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						outd2.add(ss);
					}
				}
			}else if(a.equals("dot") || a.equals("outdot")){
				outDot=b;
			}else if(a.equals("buildthreads") || a.equals("bt")){
				BUILD_THREADS=Integer.parseInt(b);
			}else if(a.equals("processcontigs")){
				processContigs=Parse.parseBoolean(b);
			}else if(a.equals("pop") || a.equalsIgnoreCase("popBubbles")){
				popBubbles=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("bubblePasses")){
				bubblePasses=Tools.max(1, Integer.parseInt(b));
			}else if(a.equalsIgnoreCase("popDirect")){
				BubblePopper.popDirect=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("popIndirect")){
				BubblePopper.popIndirect=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("debranch")){
				BubblePopper.debranch=Parse.parseBoolean(b);
			}else if(a.equals("etr") || a.equalsIgnoreCase("expandTandemRepeats")){
				expandTandemRepeats=Parse.parseBoolean(b);
				assert(false) : "TODO";
			}else if(a.equals("outkmers") || a.equals("outk") || a.equals("dump")){
				outKmers=b;
			}else if(a.equals("mincounttodump")){
				minToDump=Parse.parseIntKMG(b);
			}else if(a.equals("maxcounttodump")){
				maxToDump=Parse.parseIntKMG(b);
			}else if(a.equals("hist") || a.equals("khist")){
				outHist=b;
			}else if(a.equals("gchist")){
				gcHist=Parse.parseBoolean(b);
			}else if(a.equals("ihist") || a.equals("inserthistogram")){
				outInsert=b;
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("mode")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(Tools.isDigit(b.charAt(0))){
					processingMode=Parse.parseIntKMG(b);
				}else if(b.equalsIgnoreCase("contig")){
					processingMode=contigMode;
				}else if(b.equalsIgnoreCase("extend")){
					processingMode=extendMode;
				}else if(b.equalsIgnoreCase("correct") || b.equalsIgnoreCase("ecc") || b.equalsIgnoreCase("ecct")){
					processingMode=correctMode;
				}else if(b.equalsIgnoreCase("insert")){
					processingMode=insertMode;
				}else if(b.equalsIgnoreCase("discard") || b.equalsIgnoreCase("toss") || b.equalsIgnoreCase("filter")){
					processingMode=discardMode;
				}else{
					assert(false) : "Unknown mode "+b;
				}
			}else if(a.equals("ownership")){
				if("auto".equalsIgnoreCase(b)){
					setUseOwnership_=false;
				}else{
					useOwnership_=Parse.parseBoolean(b);
					setUseOwnership_=true;
				}
			}else if(a.equals("showstats") || a.equals("stats")){
				showStats=Parse.parseBoolean(b);
			}else if(a.equals("maxextension") || a.equals("maxe")){
				extendLeft=extendRight=Parse.parseIntKMG(b);
			}else if(a.equals("extendright") || a.equals("er")){
				extendRight=Parse.parseIntKMG(b);
			}else if(a.equals("extendleft") || a.equals("el")){
				extendLeft=Parse.parseIntKMG(b);
			}else if(a.equals("extensionrollback") || a.equals("extendrollback")){
				extensionRollback=Parse.parseIntKMG(b);
			}else if(a.equals("minextension") || a.equals("mine")){
				minExtension=Parse.parseIntKMG(b);
			}else if(a.equals("maxcontiglength") || a.equals("maxcontig") || a.equals("maxlength") || a.equals("maxlen") || a.equals("maxc")){
				maxContigLen=Parse.parseIntKMG(b);
				if(maxContigLen<0){maxContigLen=1000000000;}
			}else if(a.equals("mincontiglength") || a.equals("mincontiglen") || a.equals("mincontig") || a.equals("minc")){
				if("auto".equalsIgnoreCase(b)){
					minContigLen=-1;
				}else{
					minContigLen=Parse.parseIntKMG(b);
				}
			}else if(a.equals("mincoverage") || a.equals("mincov")){
				minCoverage=Float.parseFloat(b);
			}else if(a.equals("maxcoverage") || a.equals("maxcov")){
				if(b.equalsIgnoreCase("inf")){maxCoverage=Float.MAX_VALUE;}
				else{maxCoverage=Float.parseFloat(b);}
			}else if(a.equals("branchlower") || a.equals("branchlowerconst") || a.equals("blc")){
				branchLowerConst=Parse.parseIntKMG(b);
			}else if(a.equals("branchmult2") || a.equals("bm2")){
				branchMult2=Parse.parseIntKMG(b);
			}else if(a.equals("branchmult") || a.equals("branchmult1") || a.equals("bm1")){
				branchMult1=Parse.parseIntKMG(b);
			}else if(a.equals("mincount") || a.equals("mindepth") || a.equals("min")){
				minCountSeed=minCountExtend=Parse.parseIntKMG(b);
			}else if(a.equals("mindepthseed") || a.equals("mds") || a.equals("mincountseed") || a.equals("mcs")){
				minCountSeed=Parse.parseIntKMG(b);
			}else if(a.equals("mindepthextend") || a.equals("mde") || a.equals("mincountextend") || a.equals("mce")){
				minCountExtend=Parse.parseIntKMG(b);
			}else if(a.equals("mincountretain") || a.equals("mincr") || a.equals("mindepthretain") || a.equals("mindr")){
				kmerRangeMin=Parse.parseIntKMG(b);
			}else if(a.equals("maxcountretain") || a.equals("maxcr") || a.equals("maxdepthretain") || a.equals("maxdr")){
				kmerRangeMax=Parse.parseIntKMG(b);
			}else if(a.equals("contigpasses")){
				contigPasses=Parse.parseIntKMG(b);
			}else if(a.equals("contigpassmult")){
				contigPassMult=Double.parseDouble(b);
				assert(contigPassMult>=1) : "contigPassMult must be at least 1.";
			}else if(a.equals("threads") || a.equals("t")){
				Shared.setThreads(b);
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
//				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("verbose2")){
//				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
				verbose2=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("histcolumns")){
				histColumns=Parse.parseIntKMG(b);
			}else if(a.equals("histmax")){
				histMax=Parse.parseIntKMG(b);
			}else if(a.equals("histheader")){
				histHeader=Parse.parseBoolean(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				histZeros=!Parse.parseBoolean(b);
			}else if(a.equals("ilb") || a.equals("ignoreleftbranches") || a.equals("ignoreleftjunctions") || a.equals("ibb") || a.equals("ignorebackbranches")){
				extendThroughLeftJunctions=Parse.parseBoolean(b);
			}else if(a.equals("ibo") || a.equals("ignorebadowner")){
				IGNORE_BAD_OWNER=Parse.parseBoolean(b);
			}else if(a.equals("tad2") || a.equals("tadpole2")){
				FORCE_TADPOLE2=Parse.parseBoolean(b);
			}

			else if(a.equals("maskcore") || a.equals("coremask")){
				AbstractKmerTableSet.MASK_CORE=Kmer.MASK_CORE=Parse.parseBoolean(b);
			}else if(a.equals("fillfast") || a.equals("fastfill")){
				AbstractKmerTableSet.FAST_FILL=Parse.parseBoolean(b);
			}
			
			//Shaver
			else if(a.equals("shaverinse") || a.equals("shaveandrinse") || a.equals("wash") || a.equals("sr")){
				if(b==null || Character.isLetter(b.charAt(0))){
					removeDeadEnds=removeBubbles=Parse.parseBoolean(b);
				}else{
					maxShaveDepth=Integer.parseInt(b);
					removeDeadEnds=removeBubbles=(maxShaveDepth>0);
				}
				setShave=setRinse=true;
			}else if(a.equals("shave") || a.equals("removedeadends")){
				if(b==null || Character.isLetter(b.charAt(0))){
					removeDeadEnds=Parse.parseBoolean(b);
				}else{
					maxShaveDepth=Integer.parseInt(b);
					removeDeadEnds=(maxShaveDepth>0);
				}
				setShave=true;
			}else if(a.equals("rinse") || a.equals("shampoo") || a.equals("removebubbles")){
				removeBubbles=Parse.parseBoolean(b);
				setRinse=true;
			}else if(a.equals("maxshavedepth") || a.equals("shavedepth") || a.equals("msd")){
				maxShaveDepth=Integer.parseInt(b);
			}else if(a.equals("shavediscardlength") || a.equals("shavelength") || a.equals("discardlength") || a.equals("sdl")){
				shaveDiscardLen=Integer.parseInt(b);
			}else if(a.equals("shaveexploredistance") || a.equals("shaveexploredist") || a.equals("exploredist") || a.equals("sed")){
				shaveExploreDist=Integer.parseInt(b);
			}else if(a.equals("printeventcounts")){
				Shaver.printEventCounts=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("startFromHighCounts")){
				Shaver.startFromHighCounts=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("shaveFast") || a.equalsIgnoreCase("fastShave")){
				Shaver.shaveFast=Parse.parseBoolean(b);
			}
			
			
			//Junk removal
			else if(a.equals("tossjunk") || a.equals("tj")){
				tossJunk=Parse.parseBoolean(b);
			}else if(a.equals("tossdepth") || a.equals("discarddepth")){
				discardLowDepthReads=Integer.parseInt(b);
			}else if(a.equals("lowdepthfraction") || a.equals("ldf") || a.equals("lcf") || a.equals("lddf")){
				lowDepthDiscardFraction=Float.parseFloat(b);
			}else if(a.equals("requirebothbad") || a.equals("rbb")){
				requireBothBad=Parse.parseBoolean(b);
			}else if(a.equals("tossuncorrectable") || a.equals("tu")){
				discardUncorrectable=Parse.parseBoolean(b);
			}
			
			//Error Correction
			else if(a.equals("ecctail") || a.equals("eccright") || a.equals("tail")){
				ECC_TAIL=Parse.parseBoolean(b);
			}else if(a.equals("pincer") || a.equals("eccpincer")){
				ECC_PINCER=Parse.parseBoolean(b);
			}else if(a.equals("reassemble") || a.equals("eccreassemble")){
				ECC_REASSEMBLE=Parse.parseBoolean(b);
			}else if(a.equals("rollback") || a.equals("eccrollback")){
				ECC_ROLLBACK=Parse.parseBoolean(b);
			}else if(a.equals("bidirectional") || a.equals("requirebidirectional") || a.equals("eccbidirectional") || a.equals("rbi") || a.equals("rb")){
				ECC_REQUIRE_BIDIRECTIONAL=Parse.parseBoolean(b);
			}else if(a.equals("eccall") || a.equals("eccfull")){
				ECC_ALL=Parse.parseBoolean(b);
			}else if(a.equals("aecc") || a.equals("aec") || a.equals("aggressive")){
				ECC_AGGRESSIVE=Parse.parseBoolean(b);
				if(ECC_AGGRESSIVE){
					ECC_CONSERVATIVE=false;
					ecc_=setEcc_=true;
				}
			}else if(a.equals("cecc") || a.equals("cec") || a.equals("conservative")){
				ECC_CONSERVATIVE=Parse.parseBoolean(b);
				if(ECC_CONSERVATIVE){
					ECC_AGGRESSIVE=false;
					ecc_=setEcc_=true;
				}
			}else if(a.equals("ecc") || a.equals("ecct")){
				ecc_=Parse.parseBoolean(b);
				setEcc_=true;
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
			}else if(a.equals("ep") || a.equals("errorpath")){
				errorPath=Integer.parseInt(b);
			}else if(a.equals("em1") || a.equals("errormult1")){
				errorMult1=Float.parseFloat(b);
			}else if(a.equals("em2") || a.equals("errormult2")){
				errorMult2=Float.parseFloat(b);
			}else if(a.equals("emq") || a.equals("errormultqfactor")){
				errorMultQFactor=Float.parseFloat(b);
			}else if(a.equals("elc") || a.equals("errorlowerconst")){
				errorLowerConst=Integer.parseInt(b);
			}else if(a.equals("mcc") || a.equals("mincountcorrect")){
				minCountCorrect=Integer.parseInt(b);
			}else if(a.equals("psc") || a.equals("pathsimilarityconstant")){
				pathSimilarityConstant=Integer.parseInt(b);
			}else if(a.equals("psf") || a.equals("pathsimilarityfraction")){
				pathSimilarityFraction=Float.parseFloat(b);
			}else if(a.equals("eer") || a.equals("errorextensionreassemble")){
				errorExtensionReassemble=Integer.parseInt(b);
			}else if(a.equals("eep") || a.equals("errorextensionpincer")){
				errorExtensionPincer=Integer.parseInt(b);
			}else if(a.equals("eet") || a.equals("errorextensiontail")){
				errorExtensionTail=Integer.parseInt(b);
			}else if(a.equals("dz") || a.equals("deadzone")){
				deadZone=Integer.parseInt(b);
			}else if(a.equals("window") || a.equals("windowlen") || a.equals("w")){
				windowLen=Integer.parseInt(b);
			}else if(a.equals("windowcount") || a.equals("windowlimit") || a.equals("wc")){
				windowCount=Integer.parseInt(b);
			}else if(a.equals("windowcounthq") || a.equals("windowlimithq") || a.equals("wchq")){
				windowCountHQ=Integer.parseInt(b);
			}else if(a.equals("hqthresh") || a.equals("windowhqthresh") || a.equals("hqt")){
				windowHQThresh=(byte)Integer.parseInt(b);
			}else if(a.equals("qualsum") || a.equals("windowqualsum") || a.equals("qs")){
				windowQualSum=Integer.parseInt(b);
			}
			
			else if(a.equals("mbb") || a.equals("markbad") || a.equals("markbadbases")){
				if(b==null){b="1";}
				MARK_BAD_BASES=Integer.parseInt(b);
			}else if(a.equals("mdo") || a.equals("markdeltaonly")){
				MARK_DELTA_ONLY=Parse.parseBoolean(b);
			}else if(a.equals("meo") || a.equals("markerrorreadsonly")){
				MARK_ERROR_READS_ONLY=Parse.parseBoolean(b);
			}else if(a.equals("mq") || a.equals("markquality")){
				MARK_QUALITY=(byte)Integer.parseInt(b);
			}
			
			//Trimming
			else if(a.equals("trim") || a.equals("trimends")){
				if(b==null || Character.isLetter(b.charAt(0))){
					if("auto".equalsIgnoreCase(b)){trimEnds=-1;}
					else{trimEnds=Parse.parseBoolean(b) ? -1 : 0;}
				}else{
					trimEnds=Integer.parseInt(b);
				}
			}else if(a.equals("trimcircular") || a.equals("trimloop") || a.equals("trimlooploop")){
				trimCircular=Parse.parseBoolean(b);
			}
			
			else if(a.equalsIgnoreCase("sortbuffers") || a.equalsIgnoreCase("sortbuffer")){
				HashBuffer.SORT_BUFFERS=Parse.parseBoolean(b);
			}
			
			else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}
			
			else if(a.equals("prefilter")){
				if(b==null){prefilter=2;}
				else if(Character.isLetter(b.charAt(0))){prefilter=Parse.parseBoolean(b) ? 2 : 0;}
				else{prefilter=Integer.parseInt(b);}
			}else if(KmerTableSetU.isValidArgument(a)){
				//Do nothing
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		kmerRangeMin=Tools.max(prefilter+1, kmerRangeMin);
		
		if(outDot!=null || popBubbles){processContigs=true;}
		
		if(ECC_AGGRESSIVE){
			ECC_REQUIRE_BIDIRECTIONAL=false;
			errorExtensionReassemble=0;
			errorExtensionPincer=3;
			errorExtensionTail=4;
			errorMult1=Tools.max(4, errorMult1*0.6f);
			minCountCorrect=Tools.min(minCountCorrect, 3);
			deadZone=0;
			pathSimilarityFraction*=1.4f;
			windowCount=Tools.max(windowCount, 6);
			windowQualSum=Tools.min(windowQualSum, 65);
			windowLen=Tools.min(windowLen, 12);
			
			ECC_REQUIRE_BIDIRECTIONAL=false;
		}else if(ECC_CONSERVATIVE){
			errorExtensionReassemble+=2;
			errorExtensionPincer+=2;
			errorExtensionTail+=2;
			errorMult1=(errorMult1*1.2f);
			minCountCorrect=Tools.max(minCountCorrect, 4);
			deadZone=Tools.max(deadZone+1, 5);
			pathSimilarityFraction*=0.8f;
			windowLen=Tools.max(windowLen, 12);
			windowCount=Tools.min(windowCount, 4);
			windowQualSum=Tools.min(windowQualSum, 65);
			
			ECC_REQUIRE_BIDIRECTIONAL=true;
			ECC_ROLLBACK=true;
		}
		
//		assert(false) : ECC_REASSEMBLE;
		
		if(trimEnds<0){
			trimEnds=kbig/2;
		}
		if(minContigLen<0){
			minContigLen=Tools.max(124, 2*kbig);
		}
		
		if(verbose){
			BubblePopper.verbose=true;
//			assert(false) : "Verbose is disabled.";
//			AbstractKmerTableU.verbose=true;
		}
		THREADS=Shared.threads();
		
		assert(kmerRangeMax>=kmerRangeMin) : "kmerRangeMax must be at least kmerRangeMin: "+kmerRangeMax+", "+kmerRangeMin;
		
		if(processingMode<0){//unset
			if(ecc_ || discardUncorrectable){
				processingMode=correctMode;
				outstream.println("Switching to correct mode because ecc=t.");
			}else if(extendLeft>0 || extendRight>0){
				processingMode=extendMode;
				outstream.println("Switching to extend mode because an extend flag was set.");
			}else if(tossJunk || discardLowDepthReads>0){
				processingMode=discardMode;
				outstream.println("Switching to discard mode because a discard flag was set.");
			}else{
				processingMode=contigMode;
				//outstream.println("Operating in contig mode.");
			}
		}
		
		if(processingMode==extendMode || processingMode==correctMode || processingMode==discardMode){
			
//			{//Use in and out if ine and oute are not specified, in this mode.
//				if(ine1.isEmpty() && ine2.isEmpty()){
//					ine1.addAll(in1);
//					ine2.addAll(in2);
//				}
//				if(oute1.isEmpty() && oute2.isEmpty() && outContigs!=null){
//					oute1.add(outContigs);
//				}
//			}
			
			if(processingMode==extendMode){
				if(extendLeft==-1){extendLeft=100;}
				if(extendRight==-1){extendRight=100;}
			}else if(processingMode==correctMode){
//				extendLeft=extendRight=0;
				if(!setEcc_){ecc_=true;}
			}else if(processingMode==discardMode){
				extendLeft=extendRight=0;
				if(!setEcc_){ecc_=false;}
			}
		}
		
		{//Process parser fields
			Parser.parseQuality("","","");
			Parser.processQuality();
		}
		
		if(setUseOwnership_){
			useOwnership=useOwnership_;
		}else{
			useOwnership=(processingMode==contigMode || removeBubbles || removeDeadEnds);
		}
		
		final int extraBytesPerKmer;
		{
			int x=0;
			if(useOwnership){x+=4;}
			if(processingMode==correctMode || processingMode==discardMode){}
			else if(processingMode==contigMode || processingMode==extendMode){x+=1;}
			extraBytesPerKmer=x;
		}
		
		/* Set final variables; post-process and validate argument combinations */
		
		ecc=ecc_;
		ecco=ecco_;
		merge=merge_;
		testMerge=testMerge_;
		vstrict=vstrict_;
		
		/* Adjust I/O settings and filenames */
		
		assert(FastaReadInputStream.settingsOK());
		
		for(int i=0; i<in1.size(); i++){
			String s=in1.get(i);
			if(s!=null && s.contains("#") && !new File(s).exists()){
				int pound=s.lastIndexOf('#');
				String a=s.substring(0, pound);
				String b=s.substring(pound+1);
				in1.set(i, a+1+b);
				in2.add(a+2+b);
			}
		}
		
		for(int i=0; i<out1.size(); i++){
			String s=out1.get(i);
			if(s!=null && s.contains("#")){
				int pound=s.lastIndexOf('#');
				String a=s.substring(0, pound);
				String b=s.substring(pound+1);
				out1.set(i, a+1+b);
				out2.add(a+2+b);
			}
		}
		
		for(int i=0; i<outd1.size(); i++){
			String s=outd1.get(i);
			if(s!=null && s.contains("#")){
				int pound=s.lastIndexOf('#');
				String a=s.substring(0, pound);
				String b=s.substring(pound+1);
				outd1.set(i, a+1+b);
				outd2.add(a+2+b);
			}
		}

		nextTable=new AtomicInteger[contigPasses];
		nextVictims=new AtomicInteger[contigPasses];
		for(int i=0; i<contigPasses; i++){
			nextTable[i]=new AtomicInteger(0);
			nextVictims[i]=new AtomicInteger(0);
		}

		if(!Tools.testOutputFiles(overwrite, append, false, outKmers, outHist)){
			throw new RuntimeException("\nCan't write to some output files; overwrite="+overwrite+"\n");
		}
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			throw new RuntimeException("\nCan't write to some output files; overwrite="+overwrite+"\n");
		}		
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		assert(!in1.isEmpty()) : "Requires at least one input file.";
		if(!Tools.testInputFiles(true, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		assert(THREADS>0);
		outstream.println("Using "+THREADS+" threads.");
	}

	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public final void process(Timer t){
		
		/* Check for output file collisions */
		Tools.testOutputFiles(overwrite, append, false, outKmers, outHist);
		
		/* Count kmers */
		process2(processingMode);
		
		if(THREADS>1 && outHist!=null && outKmers!=null){
			Timer tout=new Timer();
			tout.start();
			Thread a=new DumpKmersThread();
			Thread b=new MakeKhistThread();
			a.start();
			b.start();
			while(a.getState()!=Thread.State.TERMINATED){
				try {
					a.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			while(b.getState()!=Thread.State.TERMINATED){
				try {
					b.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			tout.stop();
			outstream.println("Write Time:                 \t"+tout);
		}else{
			if(outHist!=null){
				makeKhist();
			}
			if(outKmers!=null){
				dumpKmersAsText();
			}
		}
		
		clearData();
		
		/* Stop timer and calculate speed statistics */
		t.stop();
		
		
		if(showSpeed){
			outstream.println("\nTotal Time:               \t"+t);
			
			if(processingMode==extendMode || processingMode==correctMode || processingMode==discardMode){
				outstream.println(Tools.readsBasesProcessed(t.elapsed, readsIn, basesIn, 8));
			}
		}
		
		{
			String outContigs=out1.isEmpty() ? null : out1.get(0);
			if(showStats && outContigs!=null && processingMode==contigMode && FileFormat.isFastaExt(ReadWrite.rawExtension(outContigs)) && !FileFormat.isStdio(outContigs)){
				outstream.println();
				jgi.AssemblyStats2.main(new String[] {"in="+outContigs, "printextended"});
			}
		}
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getSimpleName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	abstract void makeKhist();
	abstract void dumpKmersAsText();
	public abstract long loadKmers(Timer t);
	
	public final void clearData(){
		allContigs=null;
		tables().clear();
	}
	
	public final void process2(int mode){
		
		/* Start phase timer */
		Timer t=new Timer();
		
		/* Fill tables with kmers */
		outstream.println("\nLoading kmers.\n");
		loadKmers(t);
		
		t.stop();
//		outstream.println("Input:                      \t"+tables.readsIn+" reads \t\t"+tables.basesIn+" bases.");
//		outstream.println("Unique Kmers:               \t"+tables.kmersLoaded);
//		outstream.println("Load Time:                  \t"+t);
		
		
		t.start();
		
		if(kmerRangeMin>1 || kmerRangeMax<Integer.MAX_VALUE){
			AbstractRemoveThread.process(THREADS, kmerRangeMin, kmerRangeMax, tables(), true);
		}
		
		shaveAndRinse(t, removeDeadEnds, removeBubbles, true);
		
		if(mode==extendMode || mode==correctMode || mode==discardMode){
			outstream.println("\nExtending/error-correcting/discarding.\n");
			
			final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
			Read.VALIDATE_IN_CONSTRUCTOR=false;
			extendReads();
			Read.VALIDATE_IN_CONSTRUCTOR=vic;
			
			if(DISPLAY_PROGRESS){
				outstream.println("\nAfter extending reads:");
				Shared.printMemory();
				outstream.println();
			}
			
			t.stop();

			outstream.println("Input:                      \t"+readsIn+" reads \t\t"+basesIn+" bases.");
			outstream.println("Output:                     \t"+readsIn+" reads \t\t"+(basesIn+basesExtended)+" bases.");
			if(extendLeft>0 || extendRight>0){
				outstream.println("Bases extended:             \t"+basesExtended);
				outstream.println("Reads extended:             \t"+readsExtended+Tools.format(" \t(%.2f%%)", readsExtended*100.0/readsIn));
			}
			if(ecc){
				long partial=(readsCorrected-readsFullyCorrected);
				outstream.println("Errors detected:            \t"+(basesDetected+basesCorrectedEcco));
				{
					final long corrected=(basesCorrectedTail+basesCorrectedPincer+basesCorrectedReassemble+basesCorrectedEcco);
					StringBuilder sb=new StringBuilder();
					sb.append("Errors corrected:           \t"+Tools.padRight(corrected, 7)+" \t(");
					
					String comma="";
					if(ECC_PINCER){
						sb.append(comma).append(basesCorrectedPincer+" pincer");
						comma=", ";
					}
					if(ECC_TAIL || ECC_ALL){
						sb.append(comma).append(basesCorrectedTail+" tail");
						comma=", ";
					}
					if(ECC_REASSEMBLE){
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
						Tools.format(" \t(%.2f%%)", readsMerged*200.0/readsIn));}
				outstream.println("Reads with errors detected: \t"+Tools.padRight(readsDetected, 7)+
						Tools.format(" \t(%.2f%%)", readsDetected*100.0/readsIn));
				outstream.println("Reads fully corrected:      \t"+Tools.padRight(readsFullyCorrected, 7)+
						Tools.format(" \t(%.2f%% of detected)", readsFullyCorrected*100.0/readsDetected));
				outstream.println("Reads partly corrected:     \t"+Tools.padRight(partial, 7)+
						Tools.format(" \t(%.2f%% of detected)", partial*100.0/readsDetected));
				if(ECC_ROLLBACK || rollbacks>0){
					outstream.println("Rollbacks:                  \t"+Tools.padRight(rollbacks, 7)+
							Tools.format(" \t(%.2f%% of detected)", rollbacks*100.0/readsDetected));
				}
				
			}
			if(tossJunk || discardLowDepthReads>=0 || discardUncorrectable){
				outstream.println("Reads discarded:            \t"+readsDiscarded+Tools.format(" \t(%.2f%%)", readsDiscarded*100.0/readsIn));
				outstream.println("Bases discarded:            \t"+basesDiscarded+Tools.format(" \t(%.2f%%)", basesDiscarded*100.0/basesIn));
			}
			if(MARK_BAD_BASES>0){
				outstream.println("Reads marked:               \t"+readsMarked+Tools.format(" \t(%.2f%%)", readsMarked*100.0/readsIn));
				outstream.println("Bases marked:               \t"+basesMarked+Tools.format(" \t(%.2f%%)", basesMarked*100.0/basesIn));
			}
			
			outstream.println("Extend/error-correct time:  \t"+t);
		}else{
			/* Build contigs */
			outstream.println("\nBuilding contigs.\n");
			buildContigs(mode);
			
			if(DISPLAY_PROGRESS){
				outstream.println("\nAfter building contigs:");
				Shared.printMemory();
				outstream.println();
			}
			
			t.stop();
			
			if(readsIn>0){outstream.println("Input:                      \t"+readsIn+" reads \t\t"+basesIn+" bases.");}
			outstream.println("Bases generated:            \t"+basesBuilt);
			outstream.println("Contigs generated:          \t"+contigsBuilt);
			outstream.println("Longest contig:             \t"+longestContig);
			if(popBubbles){outstream.println("Bubbles popped:             \t"+bubblesPopped);}
			
			outstream.println("Contig-building time:       \t"+t);
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public final long shaveAndRinse(Timer t, boolean shave, boolean rinse, boolean print){
		long removed=0;
		if(shave || rinse){
			
			if(print){
				if(rinse && shave){
					outstream.println("\nRemoving dead ends and error bubbles.");
				}else if(shave){
					outstream.println("\nRemoving dead ends.");
				}else if(rinse){
					outstream.println("\nRemoving error bubbles.");
				}
			}

			removed=shave(shave, rinse);
			t.stop();

			if(print){
				outstream.println("Kmers removed:              \t"+removed);
				outstream.println("Removal time:               \t"+t);
			}

			t.start();
		}
		return removed;
	}
	
	abstract long shave(boolean shave, boolean rinse);
	abstract void initializeOwnership();
	
	/**
	 * Build contigs.
	 */
	private final void buildContigs(final int mode){
		
		if(mode==contigMode){
			allContigs=new ArrayList<Contig>();
			allInserts=null;
			
			if(useOwnership){
				Timer t=new Timer(outstream, true);
				t.start("Initializing ownership.");
				initializeOwnership();
				t.stop("Time: ");
			}
			
		}else if(mode==insertMode){
			allContigs=null;
			allInserts=new LongList();
		}else if(mode==extendMode){
			throw new RuntimeException("extendMode: TODO");
		}else{
			throw new RuntimeException("Unknown mode "+mode);
		}
		
		/* Create read input stream */
		final ConcurrentReadInputStream[] crisa=(mode==contigMode ? null : makeCrisArray(in1, in2));

		runBuildThreads(mode, crisa);
		
		/* Shut down I/O streams; capture error status */
		if(crisa!=null){
			for(ConcurrentReadInputStream cris : crisa){
				errorState|=ReadWrite.closeStreams(cris);
			}
		}
		
		if(allContigs!=null){
			ContigLengthComparator.comparator.setAscending(false);
			Shared.sort(allContigs, ContigLengthComparator.comparator);
			if(processContigs){
				processContigs();
			}
		}
		
		if(outInsert!=null){
			FileFormat ff=FileFormat.testOutput(outInsert, FileFormat.TEXT, 0, 0, true, overwrite, append, false);
			TextStreamWriter tsw=new TextStreamWriter(ff);
			tsw.start();
			for(int i=0; i<allInserts.size; i++){
				long count=allInserts.get(i);
				if(count>0 || histZeros){
					tsw.print(i+"\t"+count+"\n");
				}
			}
			errorState|=tsw.poisonAndWait();
		}
		
		String outContigs=out1.isEmpty() ? null : out1.get(0);
		if(outContigs!=null){
			FileFormat ff=FileFormat.testOutput(outContigs, FileFormat.FA, 0, 0, true, overwrite, append, false);
//			ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff, null, null, null, 4, null, false);
//			ros.start();
			ByteStreamWriter bsw=new ByteStreamWriter(ff);
			bsw.start();
			if(allContigs!=null){
				for(int i=0; i<allContigs.size(); i++){
					Contig r=allContigs.get(i);
					if(r.length()>=minContigLen){
						bsw.println(r);
					}
				}
			}
			errorState|=bsw.poisonAndWait();
		}
	}
	
	void runBuildThreads(int mode, ConcurrentReadInputStream[] crisa){
		Timer t=new Timer(outstream, true);
		
		if(BUILD_THREADS<1){BUILD_THREADS=THREADS;}
		
		/* Create ProcessThreads */
		ArrayList<AbstractBuildThread> alpt=new ArrayList<AbstractBuildThread>(BUILD_THREADS);
		for(int i=0; i<BUILD_THREADS; i++){alpt.add(makeBuildThread(i, mode, crisa));}
		for(AbstractBuildThread pt : alpt){pt.start();}
		
		/* Wait for threads to die, and gather statistics */
		for(AbstractBuildThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			for(Contig contig : pt.contigs){
				allContigs.add(contig);
				contigsBuilt++;
				basesBuilt+=contig.length();
				longestContig=Tools.max(longestContig, contig.length());
			}
			if(allInserts!=null){
				allInserts.incrementBy(pt.insertSizes);
			}
			
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			lowqReads+=pt.lowqReadsT;
			lowqBases+=pt.lowqBasesT;
		}
		t.stop("Time: ");
	}
	
	void processContigs(){
//		outstream.println("Initializing contigs.\n");
		Timer t=new Timer(outstream, true);
		outstream.println("Making contig graph.");
		initializeContigs(allContigs);
		runProcessContigThreads();
		outstream.println("Finished contig graph.");
		t.stop("Time: ");
		
		if(popBubbles){
			t.start();
			for(int popped=1, i=0; popped>0 && i<bubblePasses; i++){
				popped=popBubbles(BubblePopper.debranch && (i==bubblePasses-1));
			}
			t.stop("Time: ");
		}
		
		if(outDot!=null){
			outstream.println("Writing contig graph.");
			FileFormat ff=FileFormat.testOutput(outDot, FileFormat.TEXT, null, true, overwrite, append, false);
			ByteStreamWriter bsw=new ByteStreamWriter(ff);
			bsw.start();
			ByteBuilder bb=new ByteBuilder(1000);
			bb.append("digraph G {\n");
			for(Contig c : allContigs){
				bb.tab().append(c.id);
				bb.append(" [label=\"id=").append(c.id);
				bb.append("\\nlen=").append(c.bases.length);
				bb.append("\\ncov=").append(c.coverage, 1);
				bb.append("\\nleft=").append(codeStrings[c.leftCode]);
				bb.append("\\nright=").append(codeStrings[c.rightCode]).append("\"]").append('\n');
				if(c.leftEdges!=null){
					for(Edge e : c.leftEdges){
						bb.tab();
						e.toDot(bb);
					}
				}
				if(c.rightEdges!=null){
					for(Edge e : c.rightEdges){
						bb.tab();
						e.toDot(bb);
					}
				}
				bsw.print(bb);
				bb.clear();
			}
			bb.append("}\n");
			bsw.print(bb);
			bsw.poisonAndWait();
		}
	}
	
	HashMap<Integer, ArrayList<Edge>> destToEdgeMap(){
		HashMap<Integer, ArrayList<Edge>> destToEdgeMap=new HashMap<Integer, ArrayList<Edge>>();
		for(Contig c : allContigs){
			if(!c.used() && !c.associate()){
				if(c.leftEdges!=null){
					for(Edge e : c.leftEdges){
						Integer d=e.destination;
						ArrayList<Edge> list=destToEdgeMap.get(d);
						if(list==null){
							list=new ArrayList<Edge>(2);
							destToEdgeMap.put(d,  list);
						}
						list.add(e);
					}
				}
				if(c.rightEdges!=null){
					for(Edge e : c.rightEdges){
						Integer d=e.destination;
						ArrayList<Edge> list=destToEdgeMap.get(d);
						if(list==null){
							list=new ArrayList<Edge>(2);
							destToEdgeMap.put(d,  list);
						}
						list.add(e);
					}
				}
			}
		}
		return destToEdgeMap;
	}
	
	int popBubbles(boolean debranch){
		outstream.println("Popping bubbles; contigs="+allContigs.size());
		HashMap<Integer, ArrayList<Edge>> destToEdgeMap=destToEdgeMap();
		
		int bubblesPoppedThisPass=0;
		BubblePopper bp=new BubblePopper(allContigs, destToEdgeMap, kbig);
		
		if(debranch) {
			for(Contig c : allContigs){
				bp.debranch(c);
			}
		}
		
		for(Contig c : allContigs){
			if(!c.used() && (c.leftForwardBranch() || c.rightForwardBranch())) {
				bubblesPoppedThisPass+=bp.expand(c);
			}
		}
		
		ArrayList<Contig> temp=new ArrayList<Contig>(allContigs.size());
		
		basesBuilt=0;
		contigsBuilt=0;
		longestContig=0;
		
		for(Contig c : allContigs){
			if(c.used()){
//				assert(bp.validate(c));
				//Potentially do something about inbound edges...  there shouldn't be any, though.
			}else{
				bp.removeDeadEdges(c);//Should not be necessary but found cases where it is.
				
//				assert(bp.validate(c));
				temp.add(c);
				final int len=c.length();
				if(len>=minContigLen){
					basesBuilt+=len;
					contigsBuilt++;
					longestContig=Tools.max(len, longestContig);
				}
			}
		}
		Shared.sort(temp, ContigLengthComparator.comparator);
		
		int newID=0;
		destToEdgeMap=destToEdgeMap();
		for(Contig c : temp){
			c.renumber(newID, destToEdgeMap.get(c.id));
			newID++;
		}
		
		allContigs.clear();
		allContigs.addAll(temp);
		
//		destToEdgeMap=destToEdgeMap();
//		bp=new BubblePopper(allContigs, destToEdgeMap, kbig);
//		
//		for(int i=0; i<allContigs.size(); i++){
////			assert(i<allContigs.size()) : i+", "+allContigs.size();
//			Contig c=allContigs.get(i);
//			assert(c.id==i);
//			assert(bp.validate(c));
//		}
		
//		for(Contig c : temp){
//			if(!c.used && c.length()>minContigLen){
//				
//				allContigs.add(c);
//				
//				int len=c.length();
//				basesBuilt+=len;
//				contigsBuilt++;
//				longestContig=Tools.max(len, longestContig);
//			}
//		}
		bubblesPopped+=bubblesPoppedThisPass;
		
		outstream.println("Popped "+bubblesPoppedThisPass+" bubbles"+(bp.debranch ? "; removed "+bp.branchesRemoved+" branches." : "."));
		return bubblesPoppedThisPass;
	}
	
	void runProcessContigThreads(){
		/* Create ProcessThreads */
		AtomicInteger next=new AtomicInteger(0);
		ArrayList<AbstractProcessContigThread> alpt=new ArrayList<AbstractProcessContigThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(makeProcessContigThread(allContigs, next));}
		for(AbstractProcessContigThread pt : alpt){pt.start();}
		
		/* Wait for threads to die, and gather statistics */
		for(AbstractProcessContigThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			edgesMade+=pt.edgesMadeT;
		}
	}
	
	abstract void initializeContigs(ArrayList<Contig> contigs);
	
	abstract AbstractBuildThread makeBuildThread(int i, int mode, ConcurrentReadInputStream[] crisa);
	
	abstract AbstractProcessContigThread makeProcessContigThread(ArrayList<Contig> contigs, AtomicInteger next);
	
	/**
	 * Extend reads.
	 */
	private final void extendReads(){

		/* Create read input stream */
		final ConcurrentReadInputStream[] crisa=makeCrisArray(in1, in2);

		/* Create read input stream */
		final ConcurrentReadOutputStream[] rosa=makeCrosArray(out1, out2);

		/* Create read input stream */
		final ConcurrentReadOutputStream[] rosda=makeCrosArray(outd1, outd2);
		
		/* Create ProcessThreads */
		ArrayList<ExtendThread> alpt=new ArrayList<ExtendThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new ExtendThread(i, crisa, rosa, rosda));}
		for(ExtendThread pt : alpt){pt.start();}
		
		/* Wait for threads to die, and gather statistics */
		for(ExtendThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			lowqReads+=pt.lowqReadsT;
			lowqBases+=pt.lowqBasesT;
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
			readsDiscarded+=pt.readsDiscardedT;
			basesDiscarded+=pt.basesDiscardedT;

			readsMerged+=pt.readsMergedT;
			readsCorrectedEcco+=pt.readsCorrectedEccoT;
			basesCorrectedEcco+=pt.basesCorrectedEccoT;
		}
		
		/* Shut down I/O streams; capture error status */
		for(ConcurrentReadInputStream cris : crisa){
			errorState|=ReadWrite.closeStreams(cris);
		}
		/* Shut down I/O streams; capture error status */
		if(rosa!=null){
			for(ConcurrentReadOutputStream ros : rosa){
				errorState|=ReadWrite.closeStream(ros);
			}
		}
		if(rosda!=null){
			for(ConcurrentReadOutputStream ros : rosda){
				errorState|=ReadWrite.closeStream(ros);
			}
		}
	}
	
	private final ConcurrentReadInputStream[] makeCrisArray(ArrayList<String> list1, ArrayList<String> list2){
		final ConcurrentReadInputStream[] array;

		array=new ConcurrentReadInputStream[list1.size()];
		for(int i=0; i<list1.size(); i++){
			String a=list1.get(i);
			String b=(list2.size()>i ? list2.get(i): null);
			if(verbose){outstream.println("Creating cris for "+a);}

			final ConcurrentReadInputStream cris;
			{
				FileFormat ff1=FileFormat.testInput(a, FileFormat.FASTA, null, true, true);
				FileFormat ff2=(b==null ? null : FileFormat.testInput(b, FileFormat.FASTA, null, true, true));
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff1, ff2);
			}
			array[i]=cris;
		}
		return array;
	}
	
	private final static ConcurrentReadOutputStream[] makeCrosArray(ArrayList<String> list1, ArrayList<String> list2){
		final ConcurrentReadOutputStream[] array;

		array=new ConcurrentReadOutputStream[list1.size()];
		for(int i=0; i<list1.size(); i++){
			String a=list1.get(i);
			String b=(list2.size()>i ? list2.get(i): null);
			if(verbose){outstream.println("Creating cris for "+a);}

			final ConcurrentReadOutputStream cris;
			{
				final int buff=(!ordered ? 12 : Tools.max(32, 2*Shared.threads()));
				FileFormat ff1=FileFormat.testOutput(a, FileFormat.FASTQ, null, true, overwrite, append, ordered);
				FileFormat ff2=(b==null ? null : FileFormat.testOutput(b, FileFormat.FASTQ, null, true, overwrite, append, ordered));
//				assert(false) : list1+", "+list2;
				cris=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			}
			array[i]=cris;
		}
		return array;
	}
	
	/** Examines kmer counts around the merge borders to ensure the merge was not chimeric */ 
	public boolean mergeOK(Read merged, final int len1, final int len2, final BitSet bs, final IntList countList, final Kmer kmer, 
			final int width, final int thresh, final long mult){
		final int len=merged.length();
		final int overlap=Tools.min(len, len1+len2-len);
		final byte[] bases=merged.bases;
		if(len<len1+width+1 || len<len2+width+1 || len<kbig+width){return true;}
//		int valid=corrector.fillKmers(bases, kmers);
		final int a=len-len2-1; //Base to left of first boundary
		final int b=a+1; //Base to right of first boundary 
		final int c=len1-1;
		final int d=c+1;

		final int ak=a-kbig+1; //kmer to left of first boundary
		final int bk=b; //kmer to right of first boundary
		final int ck=c-kbig+1; //kmer to left of second boundary
		final int dk=d; //kmer to right of second boundary
		
		
		bs.clear();
		if(ak-width>=0 && ak+width<len){bs.set(ak-width, ak+width+1);} //probably from ak-width+1
		if(bk-width>=0 && ck+width<len){bs.set(bk-width, bk+width+1);} //Technically to bk+width
		if(ck-width>=0 && bk+width<len){bs.set(ck-width, ck+width+1);}
		if(dk-width>=0 && dk+width<len){bs.set(dk-width, dk+width+1);}
		tables().fillSpecificCounts(bases, countList, bs, kmer);
		
		if(ak-width>=0 && ak+width<len){
			int min=countList.get(ak);
			for(int i=ak-width+1; i<ak; i++){min=Tools.min(min, countList.get(i));}
			int min2=countList.get(ak+1);
			for(int i=ak+2; i<=ak+width; i++){min2=Tools.min(min2, countList.get(i));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*mult){return false;}
		}
		if(ck-width>=0 && ck+width<len){
//			if(dk>=len || countList.get(dk)>1){//Skip this operation if there is an error to the left
				int min=countList.get(ck);
				for(int i=ck-width+1; i<ck; i++){min=Tools.min(min, countList.get(i));}
				int min2=countList.get(ck+1);
				for(int i=ck+2; i<=ck+width; i++){min2=Tools.min(min2, countList.get(i));}
				assert(min>=0 && min2>=0);
				if(min>=thresh && min2<=1){return false;}
				if(min2>0 && min>min2*mult){return false;}
//			}
		}
		
		if(bk-width>=0 && bk+width<len){
//			if(ak<0 || countList.get(ak)>1){//Skip this operation if there is an error to the right
				int min=countList.get(bk);
				for(int i=bk+1; i<=bk+width; i++){min=Tools.min(min, countList.get(i));}
				int min2=countList.get(bk-1);
				for(int i=bk-width; i<bk-1; i++){min2=Tools.min(min2, countList.get(i));}
				assert(min>=0 && min2>=0);
				if(min>=thresh && min2<=1){return false;}
				if(min2>0 && min>min2*mult){return false;}
//			}
		}
		if(dk-width>=0 && dk+width<len){
			int min=countList.get(dk);
			for(int i=dk+1; i<=dk+width; i++){min=Tools.min(min, countList.get(i));}
			int min2=countList.get(dk-1);
			for(int i=dk-width; i<dk-1; i++){min2=Tools.min(min2, countList.get(i));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*mult){return false;}
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         ExtendThread         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Extends reads.
	 */
	private final class ExtendThread extends Thread{
		
		/**
		 * Constructor
		 * @param id_
		 * @param crisa_ Read input stream array
		 * @param rosa_
		 * @param rosda_
		 */
		public ExtendThread(int id_, ConcurrentReadInputStream[] crisa_, ConcurrentReadOutputStream[] rosa_, ConcurrentReadOutputStream[] rosda_){
			tid=id_;
			crisa=crisa_;
			rosa=rosa_;
			rosda=rosda_;
			leftCounts=extendThroughLeftJunctions ? null : KillSwitch.allocInt1D(4);
		}
		
		@Override
		public void run(){
			initializeThreadLocals();
			for(int i=0; i<crisa.length; i++){
				ConcurrentReadInputStream cris=crisa[i];
				ConcurrentReadOutputStream ros=(rosa!=null && rosa.length>i ? rosa[i] : null);
				ConcurrentReadOutputStream rosd=(rosda!=null && rosda.length>i ? rosda[i] : null);
				synchronized(crisa){
					if(!cris.started()){
						cris.start();
					}
				}
				if(ros!=null){
					synchronized(rosa){
						if(!ros.started()){
							ros.start();
						}
					}
				}
				if(rosd!=null){
					synchronized(rosda){
						if(!rosd.started()){
							rosd.start();
						}
					}
				}
				run(cris, ros, rosd);
			}
		}
		
		private void run(ConcurrentReadInputStream cris, ConcurrentReadOutputStream ros, ConcurrentReadOutputStream rosd){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//While there are more reads lists...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				final ArrayList<Read> listOut=new ArrayList<Read>(reads.size());
				final ArrayList<Read> listOutD=(rosd==null ? null : new ArrayList<Read>(reads.size()));
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					final Read r1=reads.get(i);
					final Read r2=r1.mate;
					
					processReadPair(r1, r2);
					if(r1.discarded() && (r2==null || r2.discarded())){
						readsDiscardedT+=r1.pairCount();
						basesDiscardedT+=r1.pairLength();
						if(listOutD!=null){listOutD.add(r1);}
					}else{
						listOut.add(r1);
					}
				}
				if(ros!=null){ros.add(listOut, ln.id);}
				if(rosd!=null){rosd.add(listOutD, ln.id);}
				
				//Fetch a new read list
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
		}
		
		final int findOverlap(Read r1, Read r2, boolean ecc){
			if(vstrict){
				return BBMerge.findOverlapVStrict(r1, r2, ecc);
			}else{
				return BBMerge.findOverlapStrict(r1, r2, ecc);
			}
		}
		
		private void processReadPair(final Read r10, final Read r20){
			Read r1=r10, r2=r20;
			if(verbose){outstream.println("Considering read "+r1.id+" "+new String(r1.bases));}
			final String r2id=r1.mateId();
			final int initialLength1=r1.length();
			final int initialLength2=r1.mateLength();
			
			readsInT++;
			basesInT+=r1.length();
			if(r2!=null){
				readsInT++;
				basesInT+=r2.length();
			}
			
			final Read r1_0=r1, r2_0=r2;
			if((ecco || merge) && r1!=null && r2!=null && !r1.discarded() && !r2.discarded()){
				final int insert=findOverlap(r1, r2, false);
				if(merge){
					if(insert>0){
						r2.reverseComplement();
						r1=r1.joinRead(insert);
						r2.reverseComplement();
						r2=null;
						if(testMerge && !mergeOK(r1, initialLength1, initialLength2, mergeOKBitsetT, countList, kmerT, testMergeWidth, testMergeThresh, testMergeMult)){
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
//					findOverlap(r1, r2, true);
					if(insert>0){
						r2.reverseComplement();
						Read merged=r1.joinRead(insert);
						if(!testMerge || mergeOK(merged, initialLength1, initialLength2, mergeOKBitsetT, countList, kmerT, testMergeWidth, testMergeThresh, testMergeMult)){
							int errors=BBMerge.errorCorrectWithInsert(r1, r2, insert);
							basesCorrectedEccoT+=errors;
							readsCorrectedEccoT+=(errors>0 ?1 : 0);
							readsMergedT++;
						}
						r2.reverseComplement();
					}
				}
			}
			
			processRead(r1);
			processRead(r2);
			
			//Unmerge
			if(merge && r2==null && r2id!=null && !trackerT.rollback){
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
				
				r10.bases=r1.bases;
				r10.quality=r1.quality;
				r10.flags=r1.flags;
				r20.bases=r2.bases;
				r20.quality=r2.quality;
				r20.flags=r2.flags;
			}
		}
		
		private void processRead(Read r){
			if(r==null){return;}
			if(!r.validated()){r.validate(true);}
			if(r.discarded()){
				lowqBasesT+=r.length();
				lowqReadsT++;
				return;
			}
			if(ecc || MARK_BAD_BASES>0){
				final int corrected=errorCorrect(r, leftCounts, rightCounts, kmerList, countList, countList2, builderT, builderT2, trackerT, bitsetT, kmerT, kmerT2);
				final int detected=trackerT.detected();
				final int correctedPincer=trackerT.correctedPincer;
				final int correctedTail=trackerT.correctedTail;
				final int correctedReassemble=trackerT.correctedReassemble();
				final int marked=trackerT.marked;
				assert(corrected==correctedPincer+correctedTail+correctedReassemble) : corrected+", "+trackerT;
				if(marked>0){
					readsMarkedT++;
					basesMarkedT+=marked;
				}
				if(trackerT.rollback){rollbacksT++;}
				if(detected>0){
					readsDetectedT++;
					basesDetectedT+=detected;
					if(corrected>0){
						readsCorrectedT++;
						basesCorrectedPincerT+=correctedPincer;
						basesCorrectedTailT+=correctedTail;
						basesCorrectedReassembleT+=correctedReassemble;
					}
					if(corrected==detected || (corrected>0 && countErrors(countList, r.quality)==0)){
						readsFullyCorrectedT++;
					}else if(discardUncorrectable){
						r.setDiscarded(true);
						if(r.mate!=null && !requireBothBad){
							r.mate.setDiscarded(true);
							return;
						}else if(r.mate==null){return;}
					}
				}
			}
			
			if(tossJunk && isJunk(r, rightCounts, kmerT)){
				r.setDiscarded(true);
				return;
			}
			
			if(discardLowDepthReads>=0 && hasKmersAtOrBelow(r, discardLowDepthReads, lowDepthDiscardFraction, kmerT)){
				r.setDiscarded(true);
				if(r.mate!=null && !requireBothBad){
					r.mate.setDiscarded(true);
					return;
				}else if(r.mate==null){return;}
			}
			
			int extensionRight=0, extensionLeft=0;
			if(extendRight>0){
				extensionRight=extendRead(r, builderT, leftCounts, rightCounts, extendRight, kmerT);
			}
			if(extendLeft>0){
				r.reverseComplement();
				extensionLeft=extendRead(r, builderT, leftCounts, rightCounts, extendLeft, kmerT);
				r.reverseComplement();
			}
			if(extensionRollback>0){
				int leftMod=0, rightMod=0;
				if(extensionLeft>0 && extensionLeft<extendLeft){
					leftMod=Tools.min(extensionLeft, (int)(r.numericID%(extensionRollback+1)));
					extensionLeft-=leftMod;
				}
				if(extensionRight>0 && extensionRight<extendRight){
					rightMod=Tools.min(extensionRight, (int)(r.numericID%(extensionRollback+1)));
					extensionRight-=rightMod;
				}
				if(leftMod>0 || rightMod>0){TrimRead.trimByAmount(r, leftMod, rightMod, 0);}
			}
			final int extension=extensionRight+extensionLeft;
			basesExtendedT+=extension;
			readsExtendedT+=(extension>0 ? 1 : 0);
		}
		
		/*--------------------------------------------------------------*/
		
		/** Input read stream */
		private final ConcurrentReadInputStream[] crisa;
		private final ConcurrentReadOutputStream[] rosa;
		private final ConcurrentReadOutputStream[] rosda;
		
		private final int[] leftCounts;
		private final int[] rightCounts=KillSwitch.allocInt1D(4);
		private final ErrorTracker trackerT=new ErrorTracker();
		private final ByteBuilder builderT=new ByteBuilder();
		private final ByteBuilder builderT2=new ByteBuilder();
		private final Kmer kmerT=new Kmer(kbig);
		private final Kmer kmerT2=new Kmer(kbig);
		private final BitSet bitsetT=new BitSet(300);
		private final BitSet mergeOKBitsetT=new BitSet(300);
		private final LongList kmerList=new LongList();
		private final IntList countList=new IntList();
		private final IntList countList2=new IntList();
		
		long readsInT=0;
		long basesInT=0;
		long lowqReadsT=0;
		long lowqBasesT=0;
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
		long readsDiscardedT=0;
		long basesDiscardedT=0;

		long readsMergedT=0;
		long readsCorrectedEccoT=0;
		long basesCorrectedEccoT=0;
		
		private final int tid;
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Extension Methods      ----------------*/
	/*--------------------------------------------------------------*/

	public abstract int extendRead(Read r, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int distance);
	public abstract int extendRead(Read r, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int distance, final Kmer kmer);
//	{
//		throw new RuntimeException("Must be overridden.");
//	}
	
	/**
	 * Extend these bases to the right by at most 'distance'.
	 * Stops at right junctions only.
	 * Does not claim ownership.
	 */
	public abstract int extendToRight2(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase);
	
	/**
	 * Extend these bases to the right by at most 'distance'.
	 * Stops at right junctions only.
	 * Does not claim ownership.
	 */
	public int extendToRight2(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase, Kmer kmer){
		throw new RuntimeException("Must be overridden.");
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Junk Detection        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Test a read to see if it could possibly assemble. */
	public abstract boolean isJunk(Read r);
	
	public abstract boolean isJunk(Read r, final int[] localLeftCounts, Kmer kmer);
	
	/** True if at least fraction of the reads kmers are at or below count. */
	public abstract boolean hasKmersAtOrBelow(Read r, int count, float fraction);
	
	/** True if at least fraction of the reads kmers are at or below count. */
	public abstract boolean hasKmersAtOrBelow(Read r, int count, float fraction, Kmer kmer);
	
	/*--------------------------------------------------------------*/
	/*----------------       Error Correction       ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int countErrors(IntList counts, byte[] quals){
		int possibleErrors=0;
		for(int i=1; i<counts.size; i++){
			final int a=counts.get(i-1), b=counts.get(i);
			boolean error;
			if(quals!=null){
				error=isErrorBidirectional(a, b, quals[i-1], quals[i+kbig-1]);
			}else{
				error=isErrorBidirectional(a, b, (byte)20, (byte)20);
			}
			if(error){
				possibleErrors++;
				i+=kbig;
			}
		}
		return possibleErrors;
	}
	
	public abstract int errorCorrect(Read r);
	
	public abstract int errorCorrect(Read r, final int[] leftCounts, final int[] rightCounts, LongList kmers, IntList counts, IntList counts2,
			final ByteBuilder bb, final ByteBuilder bb2, final ErrorTracker tracker, final BitSet bs, Kmer kmer, Kmer kmer2);
	
	public abstract int reassemble_inner(ByteBuilder bb, byte[] quals, final int[] rightCounts, IntList counts,
			final int extension, final Kmer kmer, final Kmer kmer2);
	
	public final int reassemble(final byte[] bases, final byte[] quals, final int[] rightCounts, final IntList counts, final IntList counts2,
			final ErrorTracker tracker, final int errorExtension, final ByteBuilder bb, final ByteBuilder bb2, final Kmer kmer, final Kmer regenKmer, BitSet bs){
		if(bases.length<kbig+1+deadZone){return 0;}
		final ByteBuilder fromLeft=new ByteBuilder(bases.length);
		final ByteBuilder fromRight=new ByteBuilder(bases.length);
		
		int detected0=tracker.detectedReassemble;
		int corrected=reassemble_pass(bases, quals, fromLeft, fromRight, rightCounts, counts, counts2, tracker, errorExtension, kmer, regenKmer, bs);
		
		int correctedIncr=corrected;
		int detectedIncr=tracker.detectedReassemble-detected0;
		int uncorrected=detectedIncr-correctedIncr;
		
		for(int passes=1; passes<6 && correctedIncr>0 && uncorrected>0; passes++){//Without a pass limit this could, in rare cases, make an infinite loop
			tracker.detectedReassemble-=uncorrected;
			detected0=tracker.detectedReassemble;
			correctedIncr=reassemble_pass(bases, quals, fromLeft, fromRight, rightCounts, counts, counts2, tracker, errorExtension, kmer, regenKmer, bs);
			
			corrected+=correctedIncr;
			detectedIncr=tracker.detectedReassemble-detected0;
			uncorrected=detectedIncr-correctedIncr;
		}
		
		return corrected;
	}
	
	public final int reassemble_pass(final byte[] bases, final byte[] quals, final ByteBuilder fromLeft, final ByteBuilder fromRight,
			final int[] rightCounts, final IntList counts, final IntList counts2, final ErrorTracker tracker, final int errorExtension,
			final Kmer kmer, final Kmer kmer2, final BitSet bs){
		if(bases.length<kbig+1+deadZone){return 0;}

		fromLeft.clear();
		fromRight.clear();
		for(byte b : bases){
			fromLeft.append(b);
			fromRight.append(b);
		}
		
		assert(counts.size>0) : counts+", "+bases.length;
		
		counts2.clear();
		counts2.addAll(counts);
		reassemble_inner(fromLeft, quals, rightCounts, counts2, errorExtension, kmer, kmer2);
		
		fromRight.reverseComplementInPlace();
		counts2.clear();
		counts2.addAll(counts);
		counts2.reverse();
		
		reassemble_inner(fromRight, quals, rightCounts, counts2, errorExtension, kmer, kmer2);
		fromRight.reverseComplementInPlace();
		
//		outstream.println(new String(fromRight));
//		outstream.println(copy);
//		outstream.println();

		int correctedInner=0;
		int correctedOuter=0;
		int detectedInner=0;
		int detectedOuter=0;
		boolean rollback=false;
		
		for(int i=0; i<bases.length; i++){
			byte a=bases[i];
			byte b=fromLeft.get(i);
			byte c=fromRight.get(i);
			if(a!=b || a!=c){
				if(b==c){detectedInner++;}
				else{
					detectedOuter++;
					if(a!=b && a!=c){
						assert(b!=c);
						rollback=true;
					}
				}
			}
			if(b==a){fromLeft.set(i, (byte)0);}
			if(c==a){fromRight.set(i, (byte)0);}
		}
		
		final int detected=detectedInner+detectedOuter;
		tracker.detectedReassemble+=detected;
		if(rollback || detected==0){return 0;}
		bs.clear();
		
		int clearedLeft=clearWindow2(fromLeft, quals, windowLen, windowCount, windowQualSum/*, windowCountHQ, windowHQThresh*/);
		fromRight.reverseInPlace();
		Tools.reverseInPlace(quals);
		int clearedRight=clearWindow2(fromRight, quals, windowLen, windowCount, windowQualSum/*, windowCountHQ, windowHQThresh*/);
		fromRight.reverseInPlace();
		Tools.reverseInPlace(quals);
		
		for(int i=0; i<bases.length; i++){
			byte a=bases[i];
			byte b=fromLeft.get(i);
			byte c=fromRight.get(i);
			byte d=a;
			if(b==0 && c==0){
				//do nothing
			}else if(b==c){
				d=b;
			}else if(b==0){
				d=c;
			}else if(c==0){
				d=b;
			}else if(b!=c){
//				if(AminoAcid.isFullyDefined(a)){
//					quals[i]=(byte)Tools.max(2, q-3);
//				}
			}
			
			if(ECC_REQUIRE_BIDIRECTIONAL && b!=c && i>=kbig && i<bases.length-kbig){d=a;}//Clause to force pincer mode in the middle
			
			if(d!=a){
				byte q=(quals==null ? 30 : quals[i]);
				if(b==c){
					correctedInner++;
					q=(byte)Tools.mid(q+qIncreasePincer, qMinPincer, qMaxPincer);
				}else{
					correctedOuter++;
					q=(byte)Tools.mid(q+qIncreaseTail, qMinTail, qMaxTail);
				}
				if(!rollback){
					bs.set(i);
					bases[i]=d;
					if(quals!=null){quals[i]=q;}
				}
			}
		}
		
		if(rollback && correctedInner+correctedOuter>0){
			tracker.rollback=true;
			return 0;
		}
		
		{
			tracker.correctedReassembleInner+=correctedInner;
			tracker.correctedReassembleOuter+=correctedOuter;
		}
		int corrected=correctedOuter+correctedInner;
		
		if(corrected>0){
			AbstractKmerTableSet tables=tables();
//			counts.clear();
//			tables.fillCounts(bases, counts, kmer);
			tables.regenerateCounts(bases, counts, kmer, bs);
			assert(counts.size>0);
		}
		
		return corrected;
	}
	
	private static int clearWindow2(final ByteBuilder bb, final byte[] quals, final int window,
			final int limit, final int qsumLimit/*, final int limitHQ, final byte hqThresh*/){
		final int len=bb.length;
		final byte[] array=bb.array;
		
		int cleared=0;
		int count=0, countHQ=0, qsum=0;
		for(int i=0, prev=-window; i<len; i++, prev++){
			byte b=array[i];
			
			if(b!=0 && (quals==null || quals[prev]>0)){
				count++;
				if(quals!=null){
					qsum+=quals[i];
//					if(quals!=null && quals[i]>=hqThresh){countHQ++;}
				}
				if(count>limit || qsum>qsumLimit /*|| countHQ>limitHQ*/){
					for(int j=Tools.max(0, i-window), lim=bb.length(); j<lim; j++){
						if(array[j]!=0){
							array[j]=0;
							cleared++;
						}
					}
					return cleared;
				}
			}
			if(prev>=0 && array[prev]>0 && (quals==null || quals[prev]>0)){
				count--;
				if(quals!=null){
					countHQ--;
//					 if(quals[prev]>=hqThresh){qsum-=quals[i];}
				}
			}
		}
		return cleared;
	}
	
	/** Changes to N any base covered strictly by kmers with count below minCount */
	public final int markBadBases(final byte[] bases, final byte[] quals, final IntList counts, final BitSet bs,
			final int minCount, boolean deltaOnly, final byte markQuality){
		if(counts.size<1){return 0;}
		
		bs.clear();
		assert(counts.size==bases.length-kbig+1) : counts.size+", "+bases.length;
		
//		boolean flag=true;
		for(int i=0; i<counts.size;){
			final int count=counts.get(i);
			if(count>=minCount){
				bs.set(i, i+kbig);
				i+=kbig;
			}else{
//				flag=false;
				i++;
			}
		}
		{//Last cycle
			final int i=counts.size-1;
			final int count=counts.get(i);
			if(count>=minCount){
				bs.set(i, i+kbig);
			}
		}
		
		final int card=bs.cardinality();
		final int toMark=bases.length-card;
		int marked=0;
		assert(card<=bases.length);
		
		int consecutiveBad=0;
		for(int i=0; i<bases.length; i++){
			if(bs.get(i)){
				consecutiveBad=0;
			}else{
				consecutiveBad++;
				boolean mark=((quals!=null && quals[i]>markQuality) || bases[i]!='N');
				if(mark && deltaOnly){
					mark=(consecutiveBad>=kbig) || bs.get(i+1) || (i>0 && bs.get(i-1));
				}
				if(mark){
					marked++;

					if(markQuality<1){
						bases[i]='N';
					}
					if(quals!=null){
						quals[i]=Tools.min(quals[i], (bases[i]=='N' ? 0 : markQuality));
					}
				}
				if(bases[i]=='N' || (quals!=null && quals[i]<=markQuality)){consecutiveBad=0;}
			}
		}
		
//		assert(toMark==0 && flag) : "toMark="+toMark+"card="+card+"len="+bases.length+"\n"+bs+"\n"+new String(bases)+"\n"+counts+"\nminCount="+minCount+"\n";
		
		return marked;
	}
	
	protected final boolean isSimilar(final int a, int loc1, int loc2, final IntList counts){
		loc1=Tools.max(loc1, 0);
		loc2=Tools.min(loc2, counts.size-1);
		for(int i=loc1; i<=loc2; i++){
			if(!isSimilar(a, counts.get(i))){return false;}
		}
		return true;
	}
	
	protected final boolean isSimilar(final int a, final int b){
		int min=Tools.min(a, b);
		int max=Tools.max(a, b);
		int dif=max-min;
		assert(dif>=0);
		return (dif<pathSimilarityConstant || dif<max*pathSimilarityFraction);
	}
	
	protected final boolean isError(final int a, int loc1, int loc2, final IntList counts){
		loc1=Tools.max(loc1, 0);
		loc2=Tools.min(loc2, counts.size-1);
		for(int i=loc1; i<=loc2; i++){
			if(!isError(a, counts.get(i))){return false;}
		}
		return true;
	}
	
	protected final boolean isErrorBidirectional(final int a, final int b, final byte qa, final byte qb){
		return (a>=b ? isError(a, b, qb) : isError(b, a, qa));
	}
	
	protected final boolean isError(final int high, final int low){
		float em1;
		if(errorPath==1){
			em1=errorMult1;
		}else if(errorPath==2){
//		if(low<minCountCorrect && high>=minCountCorrect && high>2*low){return true;}
//		final float em1=4;//Tools.mid(errorMult1, 4, low-1); //4;//(low<3 ? 4 : Tools.min(errorMult1, 2*low));
			em1=Tools.mid(errorMult1, 4, low*1.6f-3);
//		final float em1=4;
		}else{throw new RuntimeException(""+errorPath);}
		
		return (low*em1<high || (low<=errorLowerConst && high>=Tools.max(minCountCorrect, low*errorMult2)));
	}
	
	protected final boolean isError(final int high, final int low, final byte q){
		float em1;
		if(errorPath==1){
			em1=errorMult1*(1+q*errorMultQFactor);
		}else if(errorPath==2){
			if(low<minCountCorrect && high>=minCountCorrect && q<20 && high>2*low){return true;}
			//		final float em1=errorMult1*(1+q*errorMultQFactor);
			em1=Tools.mid(errorMult1, 4, low*(q<=10 ? 1.6f : 2f)-3); //final float em1=4;//(low<3 ? 4 : Tools.min(errorMult1, 2*low));
			em1=em1*(1+q*errorMultQFactor);
		}else{throw new RuntimeException(""+errorPath);}
		
		return (low*em1<high || (low<=errorLowerConst && high>=Tools.max(minCountCorrect, low*errorMult2)));
	}
	
	protected final boolean isSubstitution(int ca, int errorExtension, byte qb, IntList counts){
		final int cb=ca+1;
		final int aCount=counts.get(ca);
		final int bCount=counts.get(cb);
		if(isError(aCount, bCount, qb) && isSimilar(aCount, ca-errorExtension, ca-1, counts) &&
				isError(aCount, ca+2, ca+kbig, counts)){
			final int cc=ca+kbig;
			final int cd=cc+1;
			if(cd<counts.size){
				final int cCount=counts.get(cc);
				final int dCount=counts.get(cd);
				if(isError(aCount, dCount) || isError(dCount, cCount, qb)){
					return true;
				}
			}else{return true;}
		}
		return false;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	protected static final Kmer getKmer(byte[] bases, int loc, Kmer kmer){
		kmer.clear();
		for(int i=loc, lim=loc+kmer.kbig; i<lim; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){return null;}
			kmer.addRightNumeric(x);
		}
		assert(kmer.len==kmer.kbig);
		return kmer;
	}
	
	protected final boolean isJunction(int rightMax, int rightSecond, int leftMax, int leftSecond){
		if(isJunction(rightMax, rightSecond)){return true;}
		return isJunction(leftMax, leftSecond);
	}
	
	/** Technically, returns true if it is junction or max */
	protected final boolean isJunction(int max, int second){
		if(second<1 || second*branchMult1<max || (second<=branchLowerConst && max>=Tools.max(minCountExtend, second*branchMult2))){
			return false;
		}
		if(verbose){outstream.println("Breaking because second-highest was too high:\n" +
				"max="+max+", second="+second+", branchMult1="+branchMult1+"\n" +
				"branchLowerConst="+branchLowerConst+", minCountExtend="+minCountExtend+", branchMult2="+branchMult2+"\n" +
				(second*branchMult1<max)+", "+(second<=branchLowerConst)+", "+(max>=Tools.max(minCountExtend, second*branchMult2)));}
		return true;
	}
	
	float calcRatio(int[] counts){
		int a=0, b=0;
		for(int i=0; i<counts.length; i++){
			int x=counts[i];
			if(x>a){
				b=a;
				a=x;
			}else if(x>b){
				b=x;
			}
		}
		return (b<1 ? 99f : ((float)a)/b);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private final class DumpKmersThread extends Thread {
		
		DumpKmersThread(){}
		
		@Override
		public void run(){
			dumpKmersAsText();
		}
		
	}
	
	private final class MakeKhistThread extends Thread {
		
		MakeKhistThread(){}
		
		@Override
		public void run(){
			makeKhist();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public abstract AbstractKmerTableSet tables();
	
//	int ways; //MUST be set by subclass
	/** Big kmer length */
	final int kbig;
	public final int k(){return kbig;}

	private ArrayList<Contig> allContigs;
	private LongList allInserts;
	private long contigsBuilt=0;
	private long basesBuilt=0;
	private long longestContig=0;
	
	protected boolean extendThroughLeftJunctions=true;
	
	private boolean removeBubbles=false;
	private boolean removeDeadEnds=false;
	private boolean setShave=false;
	private boolean setRinse=false;
	protected int maxShaveDepth=1;
	protected int shaveDiscardLen=150;
	protected int shaveExploreDist=300;
	
	protected int kmerRangeMin=0;
	protected int kmerRangeMax=Integer.MAX_VALUE;
	
	protected int processingMode=-1;
	
	protected int extendLeft=-1;
	protected int extendRight=-1;
	protected int extensionRollback=3;
	
	/** Track kmer ownership */
	public final boolean useOwnership;
	
	public int maxContigLen=1000000000;
	public int minExtension=2;
	public int minContigLen=-1;
	public float minCoverage=1;
	public float maxCoverage=Float.MAX_VALUE;
	public boolean joinContigs;
	
	int trimEnds=0;
	boolean trimCircular=true;

	int minCountSeed=3;

	protected int minCountExtend=2;
	protected float branchMult1=20;
	protected float branchMult2=3;
	private int branchLowerConst=3;
	
	private int errorPath=1;
	private float errorMult1=16;
	private float errorMult2=2.6f;
	private float errorMultQFactor=0.002f;
	private int errorLowerConst=4;//3 seems fine
	private int minCountCorrect=3;//5 is more conservative...
	int minCountCorrect(){return minCountCorrect;}
	private int pathSimilarityConstant=3;
	private float pathSimilarityFraction=0.45f;//0.3
	protected int errorExtensionReassemble=5;//default 2; higher is more conservative
	protected int errorExtensionPincer=5;//default 5; higher is more conservative
	protected int errorExtensionTail=9;//default 9; higher is more conservative
	protected int deadZone=0;
	protected int windowLen=12;
	protected int windowCount=6;
	protected int windowQualSum=80;
	protected int windowCountHQ=2;
	protected byte windowHQThresh=24;
	
	protected byte qIncreasePincer=8;
	protected byte qMinPincer=24;
	protected byte qMaxPincer=32;
	
	protected byte qIncreaseTail=4;
	protected byte qMinTail=20;
	protected byte qMaxTail=28;
	
	
	public boolean showStats=true;
	
	/** Has this class encountered errors while processing? */
	public boolean errorState=false;
	
	/** Input reads */
	private ArrayList<String> in1=new ArrayList<String>(), in2=new ArrayList<String>();
	/** Output reads */
	private ArrayList<String> out1=new ArrayList<String>(), out2=new ArrayList<String>();
	/** Output discarded reads */
	private ArrayList<String> outd1=new ArrayList<String>(), outd2=new ArrayList<String>();
	/** Output graph */
	private String outDot=null;
	
//	/** Extra reads */
//	private ArrayList<String> extra=new ArrayList<String>();
	
//	/** Contig output file */
//	private String outContigs=null;
	/** Insert size histogram */
	private String outInsert=null;
	/** Kmer count output file */
	protected String outKmers=null;
	/** Histogram output file */
	protected String outHist=null;
	/** Add gc information to kmer histogram */
	protected boolean gcHist=false;

	/** Histogram columns */
	protected int histColumns=2;
	/** Histogram rows */
	protected int histMax=100000;
	/** Print a histogram header */
	protected boolean histHeader=true;
	/** Histogram show rows with 0 count */
	protected boolean histZeros=false;
	
	protected boolean smoothHist=false;
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	private long maxReads=-1;

	long edgesMade=0;
	long readsIn=0;
	long basesIn=0;
	long readsOut=0;
	long basesOut=0;
	long lowqReads=0;
	long lowqBases=0;
	long basesExtended=0;
	long readsExtended=0;
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
	long readsDiscarded=0;
	long basesDiscarded=0;
	
	long readsMerged=0;
	long readsCorrectedEcco=0;
	long basesCorrectedEcco=0;
	
	protected boolean ECC_PINCER=false;
	protected boolean ECC_TAIL=false;
	protected boolean ECC_ALL=false;
	protected boolean ECC_REASSEMBLE=true;
	protected boolean ECC_AGGRESSIVE=false;
	protected boolean ECC_CONSERVATIVE=false;
	protected boolean ECC_ROLLBACK=true;
	protected boolean ECC_REQUIRE_BIDIRECTIONAL=true;
	
	/** Mark bases as bad if they are completely covered by kmers with a count below this */
	protected int MARK_BAD_BASES=0;
	/** Only mark bad bases that are adjacent to good bases */
	protected boolean MARK_DELTA_ONLY=true;
	/** Only mark bad bases in reads that appear to have errors */
	protected boolean MARK_ERROR_READS_ONLY=true;
	/** Assign this quality score to marked bases */
	protected byte MARK_QUALITY=0;
	
	/** Discard reads that cannot be assembled */
	protected boolean tossJunk=false;
	
	/** Discard reads with kmers below this depth */
	protected int discardLowDepthReads=-1;
	
	/** Only discard reads in which at least this fraction of kmers are low depth */
	protected float lowDepthDiscardFraction=0f;
	
	/** Only discard reads if both in a pair fail */
	protected boolean requireBothBad=false;
	
	/** Discard reads with uncorrectable errors */
	protected boolean discardUncorrectable=false;
	
	/** Look for contig-contig edges */
	boolean processContigs=false;
	
	/** Pop simple bubbles in contig graph */
	boolean popBubbles=true;
	
	int bubblePasses=1;
	
	/** Expand simple tandem repeats in contig graph */
	boolean expandTandemRepeats=false;
	
	/*--------------------------------------------------------------*/
	/*----------------       ThreadLocal Temps      ----------------*/
	/*--------------------------------------------------------------*/
	
	protected final void initializeThreadLocals(){
		if(localLeftCounts.get()!=null){return;}
		localLeftCounts.set(new int[4]);
		localRightCounts.set(new int[4]);
		localExtraCounts.set(new int[4]);
		localLongList.set(new LongList());
		localIntList.set(new IntList());
		localIntList2.set(new IntList());
		localByteBuilder.set(new ByteBuilder());
		localByteBuilder2.set(new ByteBuilder());
		localBitSet.set(new BitSet(300));
		localKmer.set(new Kmer(kbig));
		localKmer2.set(new Kmer(kbig));
		localTracker.set(new ErrorTracker());
	}
	
	protected ThreadLocal<int[]> localLeftCounts=new ThreadLocal<int[]>();
	protected ThreadLocal<int[]> localRightCounts=new ThreadLocal<int[]>();
	protected ThreadLocal<int[]> localExtraCounts=new ThreadLocal<int[]>();
	protected ThreadLocal<LongList> localLongList=new ThreadLocal<LongList>();
	protected ThreadLocal<IntList> localIntList=new ThreadLocal<IntList>();
	protected ThreadLocal<IntList> localIntList2=new ThreadLocal<IntList>();
	protected ThreadLocal<ByteBuilder> localByteBuilder=new ThreadLocal<ByteBuilder>();
	protected ThreadLocal<ByteBuilder> localByteBuilder2=new ThreadLocal<ByteBuilder>();
	protected ThreadLocal<BitSet> localBitSet=new ThreadLocal<BitSet>();
	private ThreadLocal<Kmer> localKmer=new ThreadLocal<Kmer>();
	private ThreadLocal<Kmer> localKmer2=new ThreadLocal<Kmer>();
	protected ThreadLocal<ErrorTracker> localTracker=new ThreadLocal<ErrorTracker>();
	
	protected Kmer getLocalKmer(){
		Kmer local=localKmer.get();
		local.clearFast();
		return local;
	}
	
	protected Kmer getLocalKmer2(){
		Kmer local=localKmer2.get();
		local.clearFast();
		return local;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** min kmer count to dump to text */
	protected int minToDump=1;
	
	/** max kmer count to dump to text */
	protected int maxToDump=Integer.MAX_VALUE;
	
	/** Correct via kmers */
	final boolean ecc;
	
	/** Correct via overlap */
	final boolean ecco;
	
	/** Use stricter settings for merging */
	final boolean vstrict;
	
	/** Merge, correct, and unmerge */
	final boolean merge;
	
	/** Check bordering counts to see if the merge seems correct */
	final boolean testMerge;
	int testMergeWidth=4;
	long testMergeMult=80L;
	int testMergeThresh=3;
	
	/** For numbering contigs */
	final AtomicLong contigNum=new AtomicLong(0);
	
	int bubblesPopped=0;
	
	int contigPasses=16;
	double contigPassMult=1.7;
	
	/** For controlling access to tables for contig-building */
	final AtomicInteger nextTable[];
	
	/** For controlling access to victim buffers for contig-building */
	final AtomicInteger nextVictims[];
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Force output reads to stay in input order */
	public static boolean ordered=false;
	/** Print speed statistics upon completion */
	public static boolean showSpeed=true;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.threads();
	/** Number of BuildThreads */
	public static int BUILD_THREADS=-1;
	/** Do garbage collection prior to printing memory usage */
	private static final boolean GC_BEFORE_PRINT_MEMORY=false;

	static boolean IGNORE_BAD_OWNER=false;
	static boolean FORCE_TADPOLE2=false;
	
}
