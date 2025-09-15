package jgi;

import java.io.File;
import java.io.PrintStream;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicLongArray;

import bloom.KCountArray;
import bloom.KmerCountAbstract;
import bloom.ReadCounter;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import sort.ReadErrorComparator;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;
import ukmer.Kmer;


/**
 * AKA BBNorm.
 * Normalize depth by downsampling reads with high coverage.
 * Uses atomic arrays for bloom filter and an atomic histogram.
 * Succeeds KmerDownsampleAH.
 * Includes fast error correction and keep-count-based (rather than random) normalization.
 * 
 * @author Brian Bushnell
 * @date May 30, 2013
 *
 */
public class KmerNormalize {

	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		String in1=(args[0].indexOf("=")>0 ? null : args[0]);
		String in2=(in1!=null && args.length>1 ? args[1] : null);
		if(in2!=null && "null".equalsIgnoreCase(in2)){in2=null;}
		
		{
			if(in1!=null && !in1.contains(",")){
				File f=new File(in1);
				if(!f.exists() || !f.isFile()){
					in1=null;
//					throw new RuntimeException(reads1+" does not exist.");
				}
			}
			if(in2!=null && !in2.contains(",")){
				File f=new File(in2);
				if(!f.exists() || !f.isFile()){
					in2=null;
//					throw new RuntimeException(reads2+" does not exist.");
				}else if(in2.equalsIgnoreCase(in1)){
					throw new RuntimeException("Both input files are the same.");
				}
			}
		}
		
		KmerCountAbstract.minQuality=5;
		KmerCountAbstract.minProb=0.5f;
		
		Parser parser=new Parser();
		parser.trimq=TRIM_QUALITY;
		parser.minReadLength=MIN_LENGTH;
		
		int k=31;
		int cbits=32;
		int precbits=2;
		int cbits1=-1;
		int hashes=3;
//		int matrixbits=-1;
		long cells=-1;
		long maxReads=-1;
		int buildpasses=1;
		long tablereads=-1; //How many reads to process when building the hashtable
		int buildStepsize=4;
		
		String outKeep1=null;
		String outToss1=null;
		String outLow1=null, outMid1=null, outHigh1=null, outUnc1=null;

		String outKeep2=null;
		String outToss2=null;
		String outLow2=null, outMid2=null, outHigh2=null, outUnc2=null;
		
		int prehashes=-1;
		long precells=-1;
		String khistFile=null;
		String rhistFile=null;
		String peakFile=null;
		String khistFileOut=null;
		String rhistFileOut=null;
		String peakFileOut=null;
		int threads=-1;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		
		int minq=KmerCountAbstract.minQuality;
		KmerCountAbstract.CANONICAL=true;
		KmerCountAbstract.KEEP_DUPLICATE_KMERS=false;
		
		
		int targetDepthF=TARGET_DEPTH_F;
		int targetDepth1=TARGET_DEPTH_1;
		int maxDepth=MAX_DEPTH;
		int minDepth=MIN_DEPTH;
		int minKmersOverMinDepth=MIN_KMERS_OVER_MIN_DEPTH;
		float depthPercentile=DEPTH_PERCENTILE;

		int passes=2;
		boolean tossErrorReadsF=TOSS_ERROR_READS_F;
		boolean tossErrorReads1=TOSS_ERROR_READS_1;
		boolean discardBadOnlyF=DISCARD_BAD_ONLY_F;
		boolean discardBadOnly1=DISCARD_BAD_ONLY_1;
		boolean fixSpikes=FIX_SPIKES;
		float highPercentile=HIGH_PERCENTILE;
		float lowPercentile=LOW_PERCENTILE;
		int errorDetectRatio=ERROR_DETECT_RATIO;
		int hthresh=HTHRESH;
		int lthresh=LTHRESH;
		boolean countup=COUNTUP;
		boolean rbb=REQUIRE_BOTH_BAD;
		boolean setOverlap=false;
		
		
		boolean auto=true;
		
		List<String> extra=null;
		
		long memory=Runtime.getRuntime().maxMemory();
		long tmemory=Runtime.getRuntime().totalMemory();
//		assert(false) : memory+", "+tmemory;
		
		for(int i=(in1==null ? 0 : 1); i<args.length; i++){
			if(args[i]==null){args[i]="null";}
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseCommonStatic(arg, a, b)){
				if(a.equals("tbr")){//Handle conflated case
					tossErrorReads1=tossErrorReadsF=Parse.parseBoolean(b);
				}
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQualityAdjust(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("keepall")){
				KEEP_ALL=Parse.parseBoolean(b);
			}else if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("bits") ||a.equals("cbits") || a.equals("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.equals("bits1") ||a.equals("cbits1") || a.equals("cellbits1")){
				cbits1=Integer.parseInt(b);
			}else if(a.equals("prefilterbits") ||a.equals("prebits")){
				precbits=Integer.parseInt(b);
			}else if(a.equals("histlen") ||a.equals("histogramlen")){
				HIST_LEN_PRINT=Tools.min(Integer.MAX_VALUE, Long.parseLong(b)+1);
			}else if(a.equals("matrixbits")){
				int matrixbits=Integer.parseInt(b);
				assert(matrixbits<63);
				cells=1L<<matrixbits;
			}else if(a.equals("cells")){
				cells=Parse.parseKMG(b);
			}else if(a.equals("precells") || a.equals("prefiltercells")){
				precells=Parse.parseKMG(b);
				prefilter=prefilter || precells!=0;
			}else if(a.equals("prefiltersize") || a.equals("prefilterfraction")){
				prefilterFraction=Double.parseDouble(b);
				prefilter=prefilterFraction>0;
			}else if(a.equals("minq") || a.equals("minqual")){
				minq=Byte.parseByte(b);
			}else if(a.equals("zerobin")){
				ZERO_BIN=Parse.parseBoolean(b);
			}else if(a.equals("deterministic") || a.equals("dr") || a.equals("det")){
				boolean x=Parse.parseBoolean(b);
				DETERMINISTIC=x;
			}else if(a.equals("minprob")){
				KmerCountAbstract.minProb=Float.parseFloat(b);
				assert(KmerCountAbstract.minProb<1) : "minprob must be less than 1.  At 1, even kmers with 100% probablity of being error-free will be discarded.";
			}else if(a.equals("hashes")){
				hashes=Integer.parseInt(b);
			}else if(a.equals("prehashes") || a.equals("prefilterhashes")){
				prehashes=Integer.parseInt(b);
				prefilter=prefilter || prehashes!=0;
			}else if(a.equals("prefilter")){
				prefilter=Parse.parseBoolean(b);
			}else if(a.equals("countup")){
				countup=Parse.parseBoolean(b);
			}else if(a.equals("stepsize") || a.equals("buildstepsize")){
				buildStepsize=Integer.parseInt(b);
			}else if(a.equals("passes") || a.equals("p")){
				passes=Integer.parseInt(b);
				assert(passes>=1 && passes<=4) : "Passes should be in range 1~4.";
			}else if(a.equals("1pass") || a.equals("1p")){
				passes=1;
			}else if(a.equals("2pass") || a.equals("2p")){
				passes=2;
			}else if(a.equals("buildpasses")){
				buildpasses=Integer.parseInt(b);
			}else if(a.equals("printcoverage")){
				assert(false) : "This is not the program you are looking for.  Try KmerCoverage.  Move along.";
			}else if(a.equals("threads") || a.equals("t")){
				threads=Integer.parseInt(b);
			}else if(a.equals("rn") || a.equals("rename") || a.equals("renamereads")){
				renameReads=Parse.parseBoolean(b);
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("tablereads") || a.equals("buildreads")){
				tablereads=Parse.parseKMG(b);
			}else if(a.equals("out") || a.equals("out1") || a.equals("outk") || a.equals("outkeep") || a.equals("outgood")){
				outKeep1=b;
			}else if(a.equals("outt") || a.equals("outt1") || a.equals("outtoss") || a.equals("outoss") || a.equals("outbad")){
				outToss1=b;
			}else if(a.equals("outl") || a.equals("outl1") || a.equals("outlow") || a.equals("outlow1")){
				outLow1=b;
			}else if(a.equals("outm") || a.equals("outm1") || a.equals("outmid") || a.equals("outmid1") || a.equals("outmiddle")){
				outMid1=b;
			}else if(a.equals("outh") || a.equals("outh1") || a.equals("outhigh") || a.equals("outhigh1")){
				outHigh1=b;
			}else if(a.equals("outu") || a.equals("outu1") || a.equals("outuncorrected")){
				outUnc1=b;
			}else if(a.equals("out2") || a.equals("outk2") || a.equals("outkeep2") || a.equals("outgood2")){
				outKeep2=b;
			}else if(a.equals("outt2") || a.equals("outtoss2") || a.equals("outoss2") || a.equals("outbad2")){
				outToss2=b;
			}else if(a.equals("outl2") || a.equals("outlow2")){
				outLow2=b;
			}else if(a.equals("outm2") || a.equals("outmid2") || a.equals("outmiddle2")){
				outMid2=b;
			}else if(a.equals("outh2") || a.equals("outhigh2")){
				outHigh2=b;
			}else if(a.equals("outu2") || a.equals("outuncorrected2")){
				outUnc2=b;
			}else if(a.equals("lbd") || a.equals("lowbindepth") || a.equals("lowerlimit")){
				LOW_BIN_DEPTH=Integer.parseInt(b);
			}else if(a.equals("hbd") || a.equals("highbindepth") || a.equals("upperlimit")){
				HIGH_BIN_DEPTH=Integer.parseInt(b);
			}else if(a.equals("hist") || a.equals("histin") || a.equals("inhist") || a.equals("khist")){
				khistFile=b;
			}else if(a.equals("rhist")){
				rhistFile=b;
			}else if(a.equals("histout") || a.equals("outhist") || a.equals("hist2") || a.equals("khistout")){
				khistFileOut=b;
			}else if(a.equals("rhistout")){
				rhistFileOut=b;
			}else if(a.equals("peaks")){
				peakFile=b;
			}else if(a.equals("peaksout")){
				peakFileOut=b;
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered") || a.equals("ord")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("auto") || a.equals("automatic")){
				auto=Parse.parseBoolean(b);
			}else if(a.equals("canonical")){
				CANONICAL=KmerCountAbstract.CANONICAL=Parse.parseBoolean(b);
			}else if(a.equals("fixspikes") || a.equals("fs")){
				fixSpikes=Parse.parseBoolean(b);
			}else if(a.equals("printzerocoverage") || a.equals("pzc")){
				PRINT_ZERO_COVERAGE=Parse.parseBoolean(b);
			}else if(a.equals("removeduplicatekmers") || a.equals("rdk")){
				KmerCountAbstract.KEEP_DUPLICATE_KMERS=!Parse.parseBoolean(b);
			}else if(a.equals("target") || a.equals("targetdepth") || a.equals("tgt")){
				targetDepthF=Integer.parseInt(b);
			}else if(a.equals("target1") || a.equals("targetdepth1") || a.equals("tgt1")){
				targetDepth1=Integer.parseInt(b);
			}else if(a.equals("max") || a.equals("maxdepth")){
				maxDepth=Integer.parseInt(b);
			}else if(a.equals("min") || a.equals("mindepth")){
				minDepth=Integer.parseInt(b);
			}else if(a.equals("minkmers") || a.equals("minkmersovermindepth") || a.equals("mingoodkmersperread") || a.equals("mgkpr")){
				minKmersOverMinDepth=Tools.max(1, Integer.parseInt(b));
			}else if(a.equals("percentile") || a.equals("depthpercentile") || a.equals("dp")){
				depthPercentile=Float.parseFloat(b);
				if(depthPercentile>1 && depthPercentile<=100){depthPercentile/=100;}
				assert(depthPercentile>=0 && depthPercentile<=1) : "Depth percentile must be between 0 and 100.";
			}else if(a.equals("highdepthpercentile") || a.equals("highpercentile") || a.equals("hdp")){
				highPercentile=Float.parseFloat(b);
				if(highPercentile>1 && highPercentile<=100){highPercentile/=100;}
				assert(highPercentile>=0 && highPercentile<=1) : "Depth percentile must be between 0 and 100.";
			}else if(a.equals("lowdepthpercentile") || a.equals("lowpercentile") || a.equals("ldp")){
				lowPercentile=Float.parseFloat(b);
				if(lowPercentile>1 && lowPercentile<=100){lowPercentile/=100;}
				assert(lowPercentile>=0 && highPercentile<=1) : "Depth percentile must be between 0 and 100.";
			}else if(a.equals("targetbadpercentilelow") || a.equals("tbpl")){
				double d=Double.parseDouble(b);
				if(d>1 && d<=100){d/=100;}
				assert(d>=0) : "TARGET_BAD_PERCENT_LOW must be at least 0.";
				TARGET_BAD_PERCENT_LOW=d;
				TARGET_BAD_PERCENT_HIGH=Tools.max(TARGET_BAD_PERCENT_HIGH, TARGET_BAD_PERCENT_LOW);
			}else if(a.equals("targetbadpercentilehigh") || a.equals("tbph")){
				double d=Double.parseDouble(b);
				if(d>1 && d<=100){d/=100;}
				assert(d>=0 && lowPercentile<=1) : "TARGET_BAD_PERCENT_HIGH must be at least 0.";
				TARGET_BAD_PERCENT_HIGH=d;
				TARGET_BAD_PERCENT_LOW=Tools.min(TARGET_BAD_PERCENT_HIGH, TARGET_BAD_PERCENT_LOW);
			}else if(a.equals("errordetectratio") || a.equals("edr")){
				errorDetectRatio=Integer.parseInt(b);
			}else if(a.equals("errorcorrectratio") || a.equals("ecr")){
				ERROR_CORRECT_RATIO=Integer.parseInt(b);
			}else if(a.equals("highthresh") || a.equals("hthresh") || a.equals("ht")){
				hthresh=Integer.parseInt(b);
			}else if(a.equals("lowthresh") || a.equals("lthresh") || a.equals("lt")){
				lthresh=Integer.parseInt(b);
			}else if(a.equals("echighthresh") || a.equals("echthresh") || a.equals("echt")){
				EC_HTHRESH=Integer.parseInt(b);
			}else if(a.equals("eclowthresh") || a.equals("eclthresh") || a.equals("eclt")){
				EC_LTHRESH=Integer.parseInt(b);
			}else if(a.equals("markerrors") || a.equals("markonly") || a.equals("meo")){
				MARK_ERRORS_ONLY=Parse.parseBoolean(b);
			}else if(a.equals("markuncorrectableerrors") || a.equals("markuncorrectable") || a.equals("mue")){
				MARK_UNCORRECTABLE_ERRORS=Parse.parseBoolean(b);
			}else if(a.equals("tam") || a.equals("trimaftermarking")){
				TRIM_AFTER_MARKING=Parse.parseBoolean(b);
			}else if(a.equals("markwith1") || a.equals("markwithone") || a.equals("mw1")){
				MARK_WITH_1=Parse.parseBoolean(b);
//				TrimRead.PROB1=10;
			}else if(a.equals("aec") || a.equals("aecc") || a.equals("aggressiveerrorcorrection")){
				boolean x=Parse.parseBoolean(b);
				if(x){
					USE_ECC1=USE_ECCF=true;
					EC_HTHRESH=Tools.min(EC_HTHRESH, 16);
					EC_LTHRESH=Tools.max(EC_LTHRESH, 3);
					ERROR_CORRECT_RATIO=Tools.min(ERROR_CORRECT_RATIO, 100);
					MAX_ERRORS_TO_CORRECT=Tools.max(MAX_ERRORS_TO_CORRECT, 7);
					SUFFIX_LEN=Tools.min(SUFFIX_LEN, 3);
					PREFIX_LEN=Tools.min(PREFIX_LEN, 2);
				}
			}else if(a.equals("cec") || a.equals("cecc") || a.equals("conservativeerrorcorrection")){
				boolean x=Parse.parseBoolean(b);
				if(x){
					USE_ECC1=USE_ECCF=true;
					EC_HTHRESH=Tools.max(EC_HTHRESH, 30);
					EC_LTHRESH=Tools.min(EC_LTHRESH, 1);
					ERROR_CORRECT_RATIO=Tools.max(ERROR_CORRECT_RATIO, 170);
					MAX_ERRORS_TO_CORRECT=Tools.min(MAX_ERRORS_TO_CORRECT, 2);
					MAX_QUAL_TO_CORRECT=Tools.min(MAX_QUAL_TO_CORRECT, 25);
					SUFFIX_LEN=Tools.max(SUFFIX_LEN, 4);
					PREFIX_LEN=Tools.max(PREFIX_LEN, 4);
				}
			}else if(a.equals("tossbadreads") || a.equals("tosserrorreads") || a.equals("tbr") || a.equals("ter")){
				tossErrorReads1=tossErrorReadsF=Parse.parseBoolean(b);
			}else if(a.equals("tossbadreads2") || a.equals("tosserrorreads2") || a.equals("tbr2") || a.equals("ter2") ||
					a.equals("tossbadreadsf") || a.equals("tosserrorreadsf") || a.equals("tbrf") || a.equals("terf")){
				tossErrorReadsF=Parse.parseBoolean(b);
			}else if(a.equals("tossbadreads1") || a.equals("tosserrorreads1") || a.equals("tbr1") || a.equals("ter1")){
				tossErrorReads1=Parse.parseBoolean(b);
			}else if(a.equals("abrc") || a.equals("addbadreadscountup")){
				ADD_BAD_READS_COUNTUP=Parse.parseBoolean(b);
			}else if(a.equals("discardbadonly") || a.equals("dbo")){
				discardBadOnly1=discardBadOnlyF=Parse.parseBoolean(b);
			}else if(a.equals("discardbadonly1") || a.equals("dbo1")){
				discardBadOnly1=Parse.parseBoolean(b);
			}else if(a.equals("discardbadonlyf") || a.equals("dbof") || a.equals("discardbadonly2") || a.equals("dbo2")){
				discardBadOnlyF=Parse.parseBoolean(b);
			}else if(a.equals("requirebothbad") || a.equals("rbb")){
				rbb=Parse.parseBoolean(b);//Already caught by parser
			}else if(a.equals("saverarereads") || a.equals("srr")){
				SAVE_RARE_READS=Parse.parseBoolean(b);
			}else if(a.equals("eccbyoverlap") || a.equals("ecco") || a.equals("overlap")){
				if("auto".equalsIgnoreCase(b)){eccByOverlapAuto=true;}
				else{
					eccByOverlap=Parse.parseBoolean(b);
					eccByOverlapAuto=false;
				}
				setOverlap=true;
			}else if(a.equals("ecc")){
				USE_ECC1=USE_ECCF=Parse.parseBoolean(b);
			}else if(a.equals("ecc1")){
				USE_ECC1=Parse.parseBoolean(b);
			}else if(a.equals("ecc2") || a.equals("eccf")){
				USE_ECCF=Parse.parseBoolean(b);
			}else if(a.equals("ecclimit")){
				MAX_ERRORS_TO_CORRECT=Integer.parseInt(b);
			}else if(a.equals("eccmaxqual")){
				MAX_QUAL_TO_CORRECT=Integer.parseInt(b);
			}else if(a.equals("cfl")){
				CORRECT_FROM_LEFT=Parse.parseBoolean(b);
			}else if(a.equals("cfr")){
				CORRECT_FROM_RIGHT=Parse.parseBoolean(b);
			}else if(a.equals("sl") || a.equals("suflen") || a.equals("suffixlen")){
				SUFFIX_LEN=Integer.parseInt(b);
			}else if(a.equals("pl") || a.equals("prelen") || a.equals("prefixlen")){
				PREFIX_LEN=Integer.parseInt(b);
			}else if(a.equals("histcol") || a.equals("histcolumns") || a.equals("histogramcolumns")){
				HIST_COLUMNS=Integer.parseInt(b);
			}else if(a.equals("minheight")){
				minHeight=Long.parseLong(b);
			}else if(a.equals("minvolume")){
				minVolume=Long.parseLong(b);
			}else if(a.equals("minwidth")){
				minWidth=Integer.parseInt(b);
			}else if(a.equals("minpeak")){
				minPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeak")){
				maxPeak=Integer.parseInt(b);
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("maxpeakcount") || a.equals("maxpc") || a.equals("maxpeaks")){
				maxPeakCount=Integer.parseInt(b);
			}else if(a.equals("usetmpdir")){
				USE_TMPDIR=Parse.parseBoolean(b);
			}else if(a.equals("uselowerdepth") || a.equals("uld")){
				USE_LOWER_DEPTH=Parse.parseBoolean(b);
			}else if(a.equals("tmpdir")){
				TMPDIR=b;
				if(b!=null){
					b=b.trim();
					if(b.length()==0){b=null;}
					else{b=(b.replace('\\', '/')+"/").replaceAll("//", "/");}
				}
			}else if(a.equals("extra")){
				if(b!=null && !b.equalsIgnoreCase("null")){
					if(new File(b).exists()){
						extra=new ArrayList<String>();
						extra.add(b);
					}else{
						extra=Arrays.asList(b.split(","));
					}
				}
			}else{
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}
		
		{//Process parser fields
			Parser.processQuality(); //TODO: This may cause problems with multiple passes!  Best to not change the quality.
			
			TRIM_LEFT=parser.qtrimLeft;
			TRIM_RIGHT=parser.qtrimRight;
			TRIM_QUALITY=parser.trimq;
			trimE=parser.trimE();
			MIN_LENGTH=parser.minReadLength;
			rbb=parser.requireBothBad;
		}
		
		assert(passes<2 || (outLow1==null && outMid1==null && outHigh1==null && outUnc1==null)) :
			"\noutLow, outMid, outHigh, and outUnc don't work with multiple passes.  Set passes=1 or eliminate those output streams.";
		
		assert(in1!=null && !in1.equalsIgnoreCase("stdin") && !in1.toLowerCase().startsWith("stdin.")) :
			"\nThis program does not allow input from standard in,\nbecause it needs to read the input multiple times.\nOnly files are permitted.";
		
		if(MARK_ERRORS_ONLY){
			MAX_ERRORS_TO_CORRECT=Tools.max(MAX_ERRORS_TO_CORRECT, 9999);
			if(!USE_ECC1 && !USE_ECCF){USE_ECC1=true;}
		}
		
		if(!setOverlap && (USE_ECC1 || USE_ECCF)){
			eccByOverlapAuto=true;
		}
		
		if(in1!=null && in1.contains("#") && !new File(in1).exists()){
			int pound=in1.lastIndexOf('#');
			String a=in1.substring(0, pound);
			String b=in1.substring(pound+1);
			in1=a+1+b;
			in2=a+2+b;
		}
		if(in2!=null && (in2.contains("=") || in2.equalsIgnoreCase("null"))){in2=null;}
		final boolean paired=(FASTQ.FORCE_INTERLEAVED || in2!=null);
		if(in2!=null){
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		
		if(DETERMINISTIC){ordered=true;}
		
		boolean ok=Tools.testOutputFiles(overwrite, append, false, outKeep1, outToss1, outKeep2, outToss2, khistFile, khistFileOut, rhistFile, rhistFileOut, peakFile, peakFileOut);
		
		if(cbits>16 && passes>1){cbits=16;}
		
		maxDepth=Tools.max(maxDepth, targetDepthF);
		assert(targetDepthF>0);
		
		assert(FastaReadInputStream.settingsOK());
		if(k>31){CANONICAL=KmerCountAbstract.CANONICAL=false;}
		assert(CANONICAL==KmerCountAbstract.CANONICAL);
		
//		if(output!=null && reads1.contains(",")){
//			throw new RuntimeException("\nLists of input files can only be used with histogram output, not full output.\n" +
//					"Please set output=null or move the extra input files to 'extra=file1,file2,...fileN'");
//		}
		if(extra!=null){
			for(String s : extra){
				File f=new File(s);
				if(!f.exists() || !f.isFile()){throw new RuntimeException(s+" does not exist.");}
				assert(!s.equalsIgnoreCase(in1) && (in2==null || !s.equalsIgnoreCase(in2))) : "\nInput file "+s+" should not be included as an extra file.\n";
			}
		}
		
//		outstream.println("ForceInterleaved = "+FASTQ.FORCE_INTERLEAVED);
		
//		assert(false) : reads1+", "+reads2+", "+output;
//		if(FASTQ.FORCE_INTERLEAVED && in2==null){
//			outstream.println()
//		}
		
		if(threads<=0){
			THREADS=Shared.threads();
		}else{
			THREADS=threads;
		}
		
//		assert(false) : THREADS;
		
//		System.err.println("THREADS="+THREADS+", KmerCountAbstract.THREADS="+KmerCountAbstract.THREADS);
		
		long bases=0;
		qhist_total=new long[128];
		Timer t=new Timer();

		if(passes>1){
			String lastTemp1=null;
			String lastTemp2=null;

			if(passes>2){
//				System.out.println(">>>A");
				
				ERROR_DETECT_RATIO+=50;
				EC_HTHRESH=EC_HTHRESH*2+20;
				
				for(int pass=1; pass<passes-1; pass++){
					final String tempOutPrefix1=getTempPrefix(in1, outKeep1, pass, 1);
					final String tempOut1=getTempOut(outKeep1, tempOutPrefix1);
					final String tempOutToss1=(outToss1==null ? null : getTempOut(outToss1, tempOutPrefix1));

					final String tempOutPrefix2=(outKeep2==null ? null : getTempPrefix(in1, outKeep2, pass, 2));
					final String tempOut2=(outKeep2==null ? null : getTempOut(outKeep2, tempOutPrefix2));
					final String tempOutToss2=(outToss2==null ? null : getTempOut(outToss2, tempOutPrefix2));

					outstream.println("\n   ***********   Pass "+pass+"   **********   \n");

					int tgt=(targetDepth1<1 ? targetDepthF*4 : targetDepth1*2);
					int max=(tgt+tgt/4);
					
					int tgtBadLow=(int)Math.ceil(Tools.min(tgt, targetDepthF*TARGET_BAD_PERCENT_LOW*1.5));
					int tgtBadHigh=(int)Math.ceil(Tools.min(tgt, targetDepthF*TARGET_BAD_PERCENT_HIGH*1.5));
					
					CORRECT_ERRORS_THIS_PASS=USE_ECC1;
					TRIM_LEFT_THIS_PASS=(pass==1 && TRIM_LEFT);
					TRIM_RIGHT_THIS_PASS=(pass==1 && TRIM_RIGHT);
					bases+=runPass(auto, memory, (cbits1<1 ? cbits : cbits1), cells, precbits, precells, buildpasses, hashes, prehashes, k,
							maxReads, tablereads, minq, buildStepsize,
							(pass==1 ? in1 : lastTemp1), (pass==1 ? in2 : lastTemp2),
							tempOut1, tempOutToss1, null, null, null, null,
							tempOut2, tempOutToss2, null, null, null, null,
							(pass==1 ? khistFile : null), (pass==1 ? rhistFile : null), (pass==1 ? peakFile : null), (pass==1 ? extra : null),
							tgt, tgtBadLow, tgtBadHigh, max, Tools.min(minDepth, 2), Tools.min(minKmersOverMinDepth, 5),
							Tools.min(0.8f, Tools.max(0.4f, depthPercentile)*1.2f), false, rbb, true,
							highPercentile, 0, (errorDetectRatio>100 ? 100+(errorDetectRatio-100)/2 : errorDetectRatio), hthresh, lthresh, false, false, false);
					lastTemp1=tempOut1;
					lastTemp2=tempOut2;
					FASTQ.TEST_INTERLEAVED=true;
					FASTQ.FORCE_INTERLEAVED=(paired && outKeep2==null);
				}
				FASTQ.TEST_INTERLEAVED=true;
				FASTQ.FORCE_INTERLEAVED=(paired && outKeep2==null);
//				System.out.println(">>>C");
				
				ERROR_DETECT_RATIO-=50;
				EC_HTHRESH=(EC_HTHRESH-20)/2;

				for(int pass=passes-1; pass<passes; pass++){
					final String tempOutPrefix1=getTempPrefix(in1, outKeep1, pass, 1);
					final String tempOut1=getTempOut(outKeep1, tempOutPrefix1);
					final String tempOutToss1=(outToss1==null ? null : getTempOut(outToss1, tempOutPrefix1));

					final String tempOutPrefix2=(outKeep2==null ? null : getTempPrefix(in1, outKeep2, pass, 2));
					final String tempOut2=(outKeep2==null ? null : getTempOut(outKeep2, tempOutPrefix2));
					final String tempOutToss2=(outToss2==null ? null : getTempOut(outToss2, tempOutPrefix2));

					outstream.println("\n   ***********   Pass "+pass+"   **********   \n");

					int tgt=(targetDepth1<1 ? targetDepthF*2 : targetDepth1);
					int max=(tgt+tgt/4);
					
					int tgtBadLow=(int)Math.ceil(Tools.min(tgt, targetDepthF*TARGET_BAD_PERCENT_LOW));
					int tgtBadHigh=(int)Math.ceil(Tools.min(tgt, targetDepthF*TARGET_BAD_PERCENT_HIGH));
					
					CORRECT_ERRORS_THIS_PASS=USE_ECC1;
					TRIM_LEFT_THIS_PASS=(pass==1 && TRIM_LEFT);
					TRIM_RIGHT_THIS_PASS=(pass==1 && TRIM_RIGHT);
					bases+=runPass(auto, memory, (cbits1<1 ? cbits : cbits1), cells, precbits, precells, buildpasses, hashes, prehashes, k,
							maxReads, tablereads, minq, buildStepsize,
							(pass==1 ? in1 : lastTemp1), (pass==1 ? in2 : lastTemp2),
							tempOut1, tempOutToss1, null, null, null, null,
							tempOut2, tempOutToss2, null, null, null, null,
							(pass==1 ? khistFile : null), (pass==1 ? rhistFile : null), (pass==1 ? peakFile : null), (pass==1 ? extra : null),
							tgt, tgtBadLow, tgtBadHigh, max, Tools.min(minDepth, 3), minKmersOverMinDepth,
							Tools.min(0.8f, Tools.max(0.4f, depthPercentile)*1.2f), tossErrorReads1, rbb, discardBadOnly1,
							highPercentile, lowPercentile, errorDetectRatio, hthresh, lthresh, false, false, false);
					lastTemp1=tempOut1;
					lastTemp2=tempOut2;
				}
			}else{
//				System.out.println(">>>E");
				for(int pass=1; pass<passes; pass++){
					final String tempOutPrefix1=getTempPrefix(in1, outKeep1, pass, 1);
					final String tempOut1=getTempOut(outKeep1, tempOutPrefix1);
					final String tempOutToss1=(outToss1==null ? null : getTempOut(outToss1, tempOutPrefix1));

					final String tempOutPrefix2=(outKeep2==null ? null : getTempPrefix(in1, outKeep2, pass, 2));
					final String tempOut2=(outKeep2==null ? null : getTempOut(outKeep2, tempOutPrefix2));
					final String tempOutToss2=(outToss2==null ? null : getTempOut(outToss2, tempOutPrefix2));

					outstream.println("\n   ***********   Pass "+pass+"   **********   \n");

					int tgt=(targetDepth1<1 ? targetDepthF*4 : targetDepth1);
					int max=(tgt+tgt/4);
					
					int tgtBadLow=(int)Math.ceil(Tools.min(tgt, targetDepthF*TARGET_BAD_PERCENT_LOW));
					int tgtBadHigh=(int)Math.ceil(Tools.min(tgt, targetDepthF*TARGET_BAD_PERCENT_HIGH));
					
					CORRECT_ERRORS_THIS_PASS=USE_ECC1;
					TRIM_LEFT_THIS_PASS=(pass==1 && TRIM_LEFT);
					TRIM_RIGHT_THIS_PASS=(pass==1 && TRIM_RIGHT);
					bases+=runPass(auto, memory, (cbits1<1 ? cbits : cbits1), cells, precbits, precells, buildpasses, hashes, prehashes, k,
							maxReads, tablereads, minq, buildStepsize,
							(pass==1 ? in1 : lastTemp1), (pass==1 ? in2 : lastTemp2),
							tempOut1, tempOutToss1, null, null, null, null,
							tempOut2, tempOutToss2, null, null, null, null,
							(pass==1 ? khistFile : null), (pass==1 ? rhistFile : null), (pass==1 ? peakFile : null), (pass==1 ? extra : null),
							tgt, tgtBadLow, tgtBadHigh, max, Tools.min(minDepth, 3), minKmersOverMinDepth,
							Tools.min(0.8f, Tools.max(0.4f, depthPercentile)*1.2f), tossErrorReads1, rbb, discardBadOnly1,
							highPercentile, lowPercentile, errorDetectRatio, hthresh, lthresh, false, false, false);
					lastTemp1=tempOut1;
					lastTemp2=tempOut2;
					FASTQ.TEST_INTERLEAVED=true;
					FASTQ.FORCE_INTERLEAVED=(paired && outKeep2==null);
				}
//				System.out.println(">>>G");
			}
			
			outstream.println("\n   ***********   Pass "+(passes)+"   **********   \n");

			CORRECT_ERRORS_THIS_PASS=USE_ECCF;
			TRIM_LEFT_THIS_PASS=false;
			TRIM_RIGHT_THIS_PASS=false;
			bases+=runPass(auto, memory, cbits, cells, precbits, precells, buildpasses, hashes, prehashes, k,
					maxReads, tablereads, minq, buildStepsize,
					lastTemp1, lastTemp2,
					outKeep1, outToss1, outLow1, outMid1, outHigh1, outUnc1,
					outKeep2, outToss2, outLow2, outMid2, outHigh2, outUnc2,
					null, null, null, null,
					targetDepthF, targetDepthF, targetDepthF, maxDepth, minDepth, minKmersOverMinDepth, depthPercentile, tossErrorReadsF, rbb, discardBadOnlyF,
					highPercentile, lowPercentile, errorDetectRatio, hthresh, lthresh, fixSpikes, countup, renameReads);
		}else{
			
//			assert(false) : "=auto"+", m="+memory+", cb="+cbits+", cl="+cells+", pcb="+precbits+
//					"\n pcl="+precells+", bp="+buildpasses+", h="+hashes+", ph="+prehashes+", k="+k+
//					"\n ="+maxReads+", ="+tablereads+", ="+minq+", ="+buildStepsize+
//					"\n ="+in1+", ="+in2+
//					"\n ="+outKeep1+", ="+outToss1+", ="+outLow1+", ="+outMid1+", ="+outHigh1+", ="+outUnc1+
//					"\n ="+outKeep2+", ="+outToss2+", ="+outLow2+", ="+outMid2+", ="+outHigh2+", ="+outUnc2+
//					"\n ="+khistFile+", ="+rhistFile+", ="+peakFile+", ="+extra+
//					"\n ="+targetDepthF+", ="+targetDepthF+", ="+targetDepthF+", ="+maxDepth+", ="+minDepth+
//					"\n ="+minKmersOverMinDepth+", ="+depthPercentile+", ="+tossErrorReadsF+", ="+rbb+", ="+discardBadOnlyF+
//					"\n ="+highPercentile+", ="+lowPercentile+", ="+errorDetectRatio+", ="+hthresh+", ="+lthresh+
//					"\n ="+fixSpikes+", ="+countup+", ="+renameReads;
			
			CORRECT_ERRORS_THIS_PASS=(USE_ECC1 || USE_ECCF);
			TRIM_LEFT_THIS_PASS=(TRIM_LEFT);
			TRIM_RIGHT_THIS_PASS=(TRIM_RIGHT);
			bases+=runPass(auto, memory, cbits, cells, precbits, precells, buildpasses, hashes, prehashes, k,
					maxReads, tablereads, minq, buildStepsize,
					in1, in2,
					outKeep1, outToss1, outLow1, outMid1, outHigh1, outUnc1,
					outKeep2, outToss2, outLow2, outMid2, outHigh2, outUnc2,
					khistFile, rhistFile, peakFile, extra,
					targetDepthF, targetDepthF, targetDepthF, maxDepth, minDepth, minKmersOverMinDepth, depthPercentile, tossErrorReadsF, rbb, discardBadOnlyF,
					highPercentile, lowPercentile, errorDetectRatio, hthresh, lthresh, fixSpikes, countup, renameReads);
		}
		
		if(outKeep1!=null && (khistFileOut!=null || rhistFileOut!=null || peakFileOut!=null)){
			outstream.println("\n   ***********   Output Histogram Generation   **********   \n");
			FASTQ.TEST_INTERLEAVED=true;
			FASTQ.FORCE_INTERLEAVED=(paired && outKeep2==null);
			CORRECT_ERRORS_THIS_PASS=false;
			TRIM_LEFT_THIS_PASS=false;
			TRIM_RIGHT_THIS_PASS=false;
			bases+=runPass(auto, memory, cbits, cells, precbits, precells, buildpasses, hashes, prehashes, k,
					maxReads, tablereads, minq, buildStepsize,
					outKeep1, outKeep2,
					null, null, null, null, null, null,
					null, null, null, null, null, null,
					khistFileOut, rhistFileOut, peakFileOut, extra,
					99999999, 99999999, 99999999, 99999999, 0, 0, .5f, false, rbb, false,
					1, 0, 100, 10, 3, fixSpikes, false, false);
		}

		if(REMOVE_TEMP_FILES && temp_file_set!=null){
			outstream.println("\nRemoving temp files.");
			for(String s : temp_file_set){
				File f=new File(s);
				if(f.exists()){
//					System.out.println("Deleting "+s);
					boolean success=false;
					for(int i=0; i<100 && f.exists() && !success; i++){
						success=f.delete();
						f=new File(s);
					}
					if(f.exists() && !success){
//						System.err.println(f.canExecute());
//						System.err.println(f.canRead());
//						System.err.println(f.canWrite());
//						System.err.println(f.lastModified());
//						try {
//							java.nio.file.Files.delete(f.toPath());
//						} catch (IOException e) {
//							// TODO Auto-generated catch block
//							e.printStackTrace();
//						}
						System.err.println("Some temp files (prefixed TEMPFILE_BBNORM) could not be removed may need to be deleted manually.");
						f.deleteOnExit();
					}
				}
			}
		}
		
		t.stop();
		
		
		outstream.println("\nTotal time:      \t\t"+t+"   \t"+Tools.format("%.2f", bases*1000000.0/(t.elapsed))+" kb/sec");

		if(errorState){throw new RuntimeException("KmerNormalize terminated in an error state; the output may be corrupt.");}
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	private static String getTempPrefix(String inFname, String outFname, int pass, int pairnum){
		String tempOut=null, tempOutPrefix=null;
		for(int i=0; i<2000 && tempOut==null; i++){
			tempOutPrefix=(useTmpdir() ? Shared.tmpdir() : "")+"TEMPFILE_BBNORM_P"+pass+"_R"+pairnum+"_"+getSalt(inFname, i)+"_";
			tempOut=getTempOut(outFname, tempOutPrefix);
			if(new File(tempOut).exists()){
				tempOut=null;
				tempOutPrefix=null;
			}
		}
		if(tempOutPrefix==null){
			throw new RuntimeException("Can't generate a random temp file name.  Try deleting old temp files.");
		}
		return tempOutPrefix;
	}
	
	private static String getTempOut(String outFname, final String tempOutPrefix){
		assert(tempOutPrefix!=null);
		String tempOut=null;
		if(outFname==null || useTmpdir()){
			tempOut=tempOutPrefix+".fq.gz";
		}else{
			outFname=outFname.replace('\\', '/');
			int idx=outFname.lastIndexOf('/');
			if(idx<0){
				tempOut=tempOutPrefix+outFname;
			}else{
				tempOut=outFname.substring(0, idx+1)+tempOutPrefix+outFname.substring(idx+1);
			}
		}
		if(temp_file_set==null){temp_file_set=new HashSet<String>();}
		if(temp_file_set.contains(tempOut) || new File(tempOut).exists()){
			return getTempOut(outFname, tempOutPrefix+"_"+(100000*Math.random()));
		}
		temp_file_set.add(tempOut);
		return tempOut;
	}
	
	public static String getSalt(String fname, int attempt){
		return Long.toHexString(System.nanoTime()+attempt)+Long.toHexString(Long.rotateLeft(fname.hashCode(), 31)^System.currentTimeMillis());
	}
	
	private static boolean inMemorySort(ArrayList<Read> reads, String sorted, boolean reverse){
		try{
			Shared.sort(reads, ReadErrorComparator.comparator);
			if(reverse){Collections.reverse(reads);}
			TextStreamWriter tsw=new TextStreamWriter(sorted, overwrite, false, true);
			tsw.start();
//			assert(false) : "\nreads: "+reads.size()+"\n"+tsw+"\n";
			for(Read r : reads){
				tsw.println(r);
				if(r.mate!=null){tsw.println(r.mate);}
			}
			tsw.poisonAndWait();
		}catch(Throwable t){
			System.err.println("ERROR: "+t);
			return false;
		}
		return true;
	}
	
	private static long runPass(boolean auto, long memory, int cbits, long cells, int pcbits, long precells, int buildpasses, int hashes, int prehashes, int k,
			long maxReads, long tablereads, int minq, int buildStepsize,
			String in1, String in2,
			String outKeep1, String outToss1, String outLow1, String outMid1, String outHigh1, String outUnc1,
			String outKeep2, String outToss2, String outLow2, String outMid2, String outHigh2, String outUnc2,
			String khistFile, String rhistFile, String peakFile, List<String> extra,
			int targetDepth, int targetDepthBadLow, int targetDepthBadHigh, int maxDepth, int minDepth,
			int minKmersOverMinDepth, float depthPercentile, boolean tossErrorReads, boolean rbb, boolean discardBadOnly,
			float highPercentile, float lowPercentile, int errorDetectRatio, int hthresh, int lthresh, boolean fixSpikes, boolean countup,
			boolean rename){
		assert(in1!=null);
		TARGET_DEPTH=targetDepth;
		TARGET_DEPTH_BAD_LOW=targetDepthBadLow;
		TARGET_DEPTH_BAD_HIGH=targetDepthBadHigh;
		MAX_DEPTH=maxDepth;
		MIN_DEPTH=minDepth;
		MIN_KMERS_OVER_MIN_DEPTH=minKmersOverMinDepth;
		DEPTH_PERCENTILE=depthPercentile;
		RENAME_THIS_PASS=rename;

		COUNTUP=countup;
		if(COUNTUP){
			TARGET_DEPTH=(int)Math.round(TARGET_DEPTH*0.95);
		}
		TOSS_ERROR_READS=tossErrorReads;
//		REQUIRE_BOTH_BAD=(rbb);
		REQUIRE_BOTH_BAD=(rbb || COUNTUP);
		DISCARD_BAD_ONLY=discardBadOnly;
		HIGH_PERCENTILE=highPercentile;
		LOW_PERCENTILE=(COUNTUP ? LOW_PERCENTILE_COUNTUP : lowPercentile);
//		assert(!COUNTUP) : COUNTUP+", "+LOW_PERCENTILE_COUNTUP+", "+lowPercentile+", "+LOW_PERCENTILE;
		ERROR_DETECT_RATIO=errorDetectRatio;
		HTHRESH=hthresh;
		LTHRESH=lthresh;
		FIX_SPIKES=fixSpikes;
		
		{
			if(khistFile!=null || peakFile!=null){USE_KHISTOGRAM=true;}
			if(rhistFile!=null){USE_RHISTOGRAM=true;}
			
			final int maxCount=(int)(cbits>16 ? Shared.MAX_ARRAY_LEN-1 : (1L<<cbits)-1);
			assert(maxCount>0);

//			assert(false) : maxCount+", "+HIST_LEN_PRINT+", "+HIST_LEN+", "+THREAD_HIST_LEN;
				
			HIST_LEN_PRINT=Tools.max(2, Tools.min(HIST_LEN_PRINT, maxCount+1));
			assert(HIST_LEN_PRINT<=Integer.MAX_VALUE) : HIST_LEN_PRINT+", "+Integer.MAX_VALUE;
			HIST_LEN=(int)Tools.min(maxCount+1, Tools.max(HIST_LEN_PRINT, HIST_LEN));
			THREAD_HIST_LEN=Tools.min(THREAD_HIST_LEN, HIST_LEN);
			
			khistogram=new AtomicLongArray(HIST_LEN);
			if(USE_RHISTOGRAM && rhistFile!=null){
				rhistogram=new AtomicLongArray(HIST_LEN);
				bhistogram=new AtomicLongArray(HIST_LEN);
			}else{
				rhistogram=null;
				bhistogram=null;
			}
		}
		
		if(auto && cells==-1){
			final long usable=(long)Tools.max(((memory-96000000)*.73), memory*0.45);
			long mem=usable-(khistogram!=null ? (HIST_LEN*16L*(1)) : 0)-(rhistogram!=null ? (HIST_LEN*8L*(1)) : 0)-(bhistogram!=null ? (HIST_LEN*8L*(1)) : 0);
			if(buildpasses>1){mem/=2;}
			
			FILTERBYTES=(COUNTUP ? mem/2 : mem);
			cells=(FILTERBYTES*8)/cbits;
//
//			long tablebytes=((1L<<matrixbits)*cbits)/8;
//			if(tablebytes*3<usable){matrixbits++;}
//			outstream.println(tablebytes/1000000+", "+usable/1000000+", "+(tablebytes*3)/1000000);
			
		}else if(cells==-1){
			cells=1L<<34;
		}
		
		if(prefilter){
			if(precells<1){
				long totalbits=cells*cbits;
				long prebits=(long)(totalbits*prefilterFraction);
				precells=prebits/pcbits;
				cells=(totalbits-prebits+cbits-1)/cbits; //Steal memory from cell allocation
			}
			if(prehashes<1){
				prehashes=(hashes+1)/2;
			}
		}
		
		{
			outstream.println("\nSettings:");
			outstream.println("threads:          \t"+THREADS);
			outstream.println("k:                \t"+k);
			outstream.println("deterministic:    \t"+DETERMINISTIC);
			outstream.println("toss error reads: \t"+TOSS_ERROR_READS);
			outstream.println("passes:           \t"+buildpasses);
			outstream.println("bits per cell:    \t"+cbits);
//			outstream.println("matrixbits: \t"+matrixbits);
			outstream.println("cells:            \t"+Tools.toKMG(cells));
			outstream.println("hashes:           \t"+hashes);
			if(prefilter){
				outstream.println("prefilter bits:   \t"+pcbits);
//				outstream.println("matrixbits: \t"+matrixbits);
				outstream.println("prefilter cells:  \t"+(precells>0 && prehashes>0 ? Tools.toKMG(precells) : "?"));
				outstream.println("prefilter hashes: \t"+(precells>0 && prehashes>0 ? ""+prehashes : "?"));
			}
//			outstream.println("base min quality: \t"+KmerCountAbstract.minQuality);
			outstream.println("base min quality: \t"+minq);
			outstream.println("kmer min prob:    \t"+KmerCountAbstract.minProb);
			
			outstream.println();
			outstream.println("target depth:     \t"+TARGET_DEPTH);
			outstream.println("min depth:        \t"+MIN_DEPTH);
			outstream.println("max depth:        \t"+MAX_DEPTH);
			outstream.println("min good kmers:   \t"+MIN_KMERS_OVER_MIN_DEPTH);
			outstream.println("depth percentile: \t"+Tools.format("%.1f", 100*DEPTH_PERCENTILE));
			outstream.println("ignore dupe kmers:\t"+!KmerCountAbstract.KEEP_DUPLICATE_KMERS);
			outstream.println("fix spikes:       \t"+FIX_SPIKES);
			if((USE_KHISTOGRAM || USE_RHISTOGRAM) && HIST_LEN>0){
				outstream.println("histogram length: \t"+HIST_LEN);
			}
			if(khistFile!=null || rhistFile!=null){
				outstream.println("print zero cov:   \t"+PRINT_ZERO_COVERAGE);
			}
			
			outstream.println();
		}
		
		if(!prefilter && k<32 && cells>(1L<<(2*k))){cells=(1L<<(2*k));}
		assert(cells>0);
		
//		KmerCountAbstract.THREADS=Tools.max(THREADS/2, KmerCountAbstract.THREADS);  //Seems like 4 is actually optimal...
		
		FastaReadInputStream.MIN_READ_LEN=k;
		
		if(eccByOverlapAuto){
			eccByOverlapAuto=false;
			float overlapRatio=BBMerge.mergeableFraction(in1, in2, 1000000, 0.01f);
			eccByOverlap=(overlapRatio>0.25f);
			if(eccByOverlap){
				System.err.println("Enabled overlap correction ("+Tools.format("%.1f%% percent overlap)", 100*overlapRatio));
			}
		}
		
		Timer t=new Timer();
		Timer ht=new Timer();
		t.start();
		ht.start();
		KCountArray kca;
		KCountArray prefilterArray=null;
//		outstream.println();
		ReadCounter rc=new ReadCounter(k, true, eccByOverlap, false, Shared.AMINO_IN);
		rc.detectStepsize=buildStepsize;
		if(prefilter){
			prefilterArray=rc.makeKca(in1, in2, extra, pcbits, precells, prehashes, minq,
					tablereads, 1, 1, 1, null, 0);
			outstream.println("Made prefilter:   \t"+prefilterArray.toShortString(prehashes));
			double uf=prefilterArray.usedFraction();
			if(uf>0.6){
				outstream.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.95 ? "incredibly" : uf>0.9 ? "extremely" : uf>0.8 ? "very" :
					uf>0.7 ? "fairly" : "somewhat")+" full, which may reduce accuracy for kmers of depth under 3.  Ideal load is under 60% used." +
						"\nFor better accuracy, run on a node with more memory; quality-trim or error-correct reads; " +
						"or increase the values of the minprob flag to reduce spurious kmers.");
			}
		}
		
		kca=rc.makeKca(in1, in2, extra, cbits, cells, hashes, minq,
				tablereads, buildpasses, 2, 2, prefilterArray, (prefilterArray==null ? 0 : prefilterArray.maxValue));
		ht.stop();
		
		outstream.println("Made hash table:  \t"+kca.toShortString(hashes));
		double uf=kca.usedFraction();
		if(uf>0.6){
			outstream.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.95 ? "incredibly" : uf>0.9 ? "extremely" : uf>0.8 ? "very" :
				uf>0.7 ? "fairly" : "somewhat")+" full, which may reduce accuracy.  Ideal load is under 60% used." +
				"\nFor better accuracy, use the 'prefilter' flag; run on a node with more memory; quality-trim or error-correct reads; " +
					"or increase the values of the minprob flag to reduce spurious kmers.  In practice you should still get good normalization results " +
					"even with loads over 90%, but the histogram and statistics will be off.");
		}
		
		long estUnique;
		outstream.println();
		if(prefilterArray!=null){
			int lim1=prefilterArray.maxValue, lim2=prefilterArray.maxValue+1;
			double a=prefilterArray.estimateUniqueKmers(prehashes);
			double b=kca.estimateUniqueKmers(hashes, lim2);
			a=a-b;
			if(CANONICAL){
//				a=(a*KCountArray.canonMask)/(KCountArray.canonMask+1);
//				b=(b*KCountArray.canonMask)/(KCountArray.canonMask+1);
			}else{
				a/=2;
				b/=2;
			}
			estUnique=((long)((a+b)));
			outstream.println("Estimated kmers of depth 1-"+lim1+": \t"+(long)a);
			outstream.println("Estimated kmers of depth "+lim2+"+ : \t"+(long)b);
		}else{
//			double est=kca.cells*(1-Math.pow(1-Math.sqrt(kca.usedFraction()), 1.0/hashes));
//			double est=kca.cells*(1-Math.pow(1-kca.usedFraction(), 1.0/hashes));
			double est=kca.estimateUniqueKmers(hashes);
//			outstream.println("Used cells: "+kca.cellsUsed(1));
			if(CANONICAL){
//				est=(est*KCountArray.canonMask)/(KCountArray.canonMask+1);
			}else{
				est/=2;
			}
			estUnique=((long)((est)));
			
		}
		outstream.println("Estimated unique kmers:     \t"+estUnique);//+", or "+estUnique+" counting forward kmers only.");
//		outstream.println("(Includes forward and reverse kmers)");
		outstream.println();
		outstream.println("Table creation time:\t\t"+ht);//+"   \t"+Tools.format("%.2f", totalBases*1000000.0/(ht.elapsed))+" kb/sec");
		
		ListNum.setDeterministicRandom(DETERMINISTIC);
		
		long bases=0;
		if(COUNTUP){
			COUNTUP=false;
			
			int td0=TARGET_DEPTH, md0=MIN_DEPTH, mxd0=MAX_DEPTH, mkomd0=MIN_KMERS_OVER_MIN_DEPTH;
			TARGET_DEPTH=TARGET_DEPTH*4;
			MIN_DEPTH=MIN_DEPTH/2;
			MAX_DEPTH=MAX_DEPTH*4;
			MIN_KMERS_OVER_MIN_DEPTH=MIN_KMERS_OVER_MIN_DEPTH/2;
			
			int rnd=(int)(100+Math.random()*1000000);
			final String tempOutPrefix1=getTempPrefix(in1, outKeep1, rnd, 1);
			final String tempOut1=getTempOut(outKeep1, tempOutPrefix1);
			final String tempOutPrefix2=(outKeep2==null ? null : getTempPrefix(in1, outKeep2, rnd, 2));
			final String tempOut2=(outKeep2==null ? null : getTempOut(outKeep2, tempOutPrefix2));
			ArrayList<Read> storage=new ArrayList<Read>();
			
			if(in1!=null && in1.contains(",") && !new File(in1).exists()){
				String[] list1=in1.split(",");
				String[] list2=(in2==null ? null : in2.split(","));
				bases+=count(list1, list2, kca, k, maxReads, null, null, null, null, null, null,
						null, null, null, null, null, null, false, overwrite, null, null, null, estUnique, storage);
			}else{
				bases+=count(in1, in2, kca, k, maxReads, null, null, null, null, null, null,
						null, null, null, null, null, null, false, overwrite, null, null, null, estUnique, storage);
			}
			inMemorySort(storage, tempOut1, false);
			storage=null;
			in1=tempOut1;
			in2=null;
			
			TARGET_DEPTH=td0;
			MIN_DEPTH=md0;
			MAX_DEPTH=mxd0;
			MIN_KMERS_OVER_MIN_DEPTH=mkomd0;
			
			COUNTUP=true;
			
			
			if(in1!=null && in1.contains(",") && !new File(in1).exists()){
				String[] list1=in1.split(",");
				String[] list2=(in2==null ? null : in2.split(","));
				bases+=count(list1, list2, kca, k, maxReads, outKeep1, outToss1, outLow1, outMid1, outHigh1, outUnc1,
						outKeep2, outToss2, outLow2, outMid2, outHigh2, outUnc2, ordered, overwrite, khistFile, rhistFile, peakFile, estUnique, null);
			}else{
				bases+=count(in1, in2, kca, k, maxReads, outKeep1, outToss1, outLow1, outMid1, outHigh1, outUnc1,
						outKeep2, outToss2, outLow2, outMid2, outHigh2, outUnc2, ordered, overwrite, khistFile, rhistFile, peakFile, estUnique, null);
			}
			
		}else{


			if(in1!=null && in1.contains(",") && !new File(in1).exists()){
				String[] list1=in1.split(",");
				String[] list2=(in2==null ? null : in2.split(","));
				bases+=count(list1, list2, kca, k, maxReads, outKeep1, outToss1, outLow1, outMid1, outHigh1, outUnc1,
						outKeep2, outToss2, outLow2, outMid2, outHigh2, outUnc2, ordered, overwrite, khistFile, rhistFile, peakFile, estUnique, null);
			}else{
				bases+=count(in1, in2, kca, k, maxReads, outKeep1, outToss1, outLow1, outMid1, outHigh1, outUnc1,
						outKeep2, outToss2, outLow2, outMid2, outHigh2, outUnc2, ordered, overwrite, khistFile, rhistFile, peakFile, estUnique, null);
			}
		}
		
		if(ANALYZE_TOPOLOGY){printTopology();}
		
		t.stop();
//		outstream.println("\nTotal time:      \t\t"+t+"   \t"+Tools.format("%.2f", bases*1000000.0/(t.elapsed))+" kb/sec");
		return bases;
	}
	
	
	public static void printTopology(){
		long total=peaks.get()+spikes.get()+flats.get()+valleys.get()+slopes.get();
		double mult=100.0/total;
		
		long sp=spikes.get();
		long pe=peaks.get();
		long va=valleys.get();
		long sl=slopes.get();
		long fl=flats.get();
		double dsp=mult*sp;
		double dpe=mult*pe;
		double dva=mult*va;
		double dsl=mult*sl;
		double dfl=mult*fl;
		
		System.err.println("\nDepth Topology:\t");
		System.err.println("Spikes:     \t\t\t"+(dsp<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dsp,sp));
		System.err.println("Peaks:      \t\t\t"+(dpe<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dpe,pe));
		System.err.println("Valleys:    \t\t\t"+(dva<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dva,va));
		System.err.println("Slopes:     \t\t\t"+(dsl<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dsl,sl));
		System.err.println("Flats:      \t\t\t"+(dfl<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dfl,fl));
	}


	public static long count(String in1, String in2, KCountArray kca, int k, long maxReads,
			String outKeep1, String outToss1, String outLow1, String outMid1, String outHigh1, String outUnc1,
			String outKeep2, String outToss2, String outLow2, String outMid2, String outHigh2, String outUnc2,
			boolean ordered, boolean overwrite, String khistFile, String rhistFile, String peakFile, long estUnique, ArrayList<Read> storage) {
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		ConcurrentReadOutputStream rosKeep=null;
		if(outKeep1!=null){
			final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
			
			final String out=outKeep1;
			String out1=outKeep1.replaceFirst("#", "1");
			String out2=outKeep2;
			
			if(cris.paired() && out2==null){
				if(out.contains("#")){
					out2=out.replaceFirst("#", "2");
				}else{
					outstream.println("Writing interleaved.");
				}
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2)));
			
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			rosKeep=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
			rosKeep.start();
			outstream.println("Started output threads.");
		}
		
		ConcurrentReadOutputStream rosToss=null;
		if(outToss1!=null){
			final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
			
			final String out=outToss1;
			String out1=outToss1.replaceFirst("#", "1");
			String out2=outToss2;
			
			if(cris.paired() && out2==null){
				if(out.contains("#")){
					out2=out.replaceFirst("#", "2");
				}else{
					outstream.println("Writing interleaved.");
				}
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2)));
			
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			rosToss=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
			
			rosToss.start();
			outstream.println("Started output threads.");
		}
		
		ConcurrentReadOutputStream rosLow=null;
		if(outLow1!=null){
			final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
			
			final String out=outLow1;
			String out1=outLow1.replaceFirst("#", "1");
			String out2=outLow2;
			
			if(cris.paired() && out2==null){
				if(out.contains("#")){
					out2=out.replaceFirst("#", "2");
				}else{
					outstream.println("Writing interleaved.");
				}
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2)));
			
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			rosLow=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
			
			rosLow.start();
			outstream.println("Started output threads.");
		}
		
		ConcurrentReadOutputStream rosMid=null;
		if(outMid1!=null){
			final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
			
			final String out=outMid1;
			String out1=outMid1.replaceFirst("#", "1");
			String out2=outMid2;
			
			if(cris.paired() && out2==null){
				if(out.contains("#")){
					out2=out.replaceFirst("#", "2");
				}else{
					outstream.println("Writing interleaved.");
				}
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2)));
			
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			rosMid=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
			
			rosMid.start();
			outstream.println("Started output threads.");
		}
		
		ConcurrentReadOutputStream rosHigh=null;
		if(outHigh1!=null){
			final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
			
			final String out=outHigh1;
			String out1=outHigh1.replaceFirst("#", "1");
			String out2=outHigh2;
			
			if(cris.paired() && out2==null){
				if(out.contains("#")){
					out2=out.replaceFirst("#", "2");
				}else{
					outstream.println("Writing interleaved.");
				}
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2)));
			
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			rosHigh=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
			
			rosHigh.start();
			outstream.println("Started output threads.");
		}
		
		ConcurrentReadOutputStream rosUnc=null;
		if(outUnc1!=null){
			final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
			
			final String out=outUnc1;
			String out1=outUnc1.replaceFirst("#", "1");
			String out2=outUnc2;
			
			if(cris.paired() && out2==null){
				if(out.contains("#")){
					out2=out.replaceFirst("#", "2");
				}else{
					outstream.println("Writing interleaved.");
				}
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2)));
			
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			rosUnc=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
			
			rosUnc.start();
			outstream.println("Started output threads.");
		}
		
		long bases=downsample(cris, kca, k, maxReads, rosKeep, rosToss, rosLow, rosMid, rosHigh, rosUnc, khistFile, rhistFile, peakFile, overwrite, estUnique, storage);
		
		errorState|=ReadWrite.closeStreams(cris, rosKeep, rosToss, rosLow, rosMid, rosHigh, rosUnc);
		if(verbose){System.err.println("Closed streams");}
		
		return bases;
	}
	
	
	public static long count(String[] list1, String[] list2, KCountArray kca, int k, long maxReads,
			String outKeep1, String outToss1, String outLow1, String outMid1, String outHigh1, String outUnc1,
			String outKeep2, String outToss2, String outLow2, String outMid2, String outHigh2, String outUnc2,
			boolean ordered, boolean overwrite, String khistFile, String rhistFile, String peakFile, long estUnique, ArrayList<Read> storage) {
		
		ConcurrentReadOutputStream rosKeep=null, rosToss=null, rosLow=null, rosMid=null, rosHigh=null, rosUnc=null;
		String[] outKeep1a=null, outKeep2a=null;
		String[] outToss1a=null, outToss2a=null;
		String[] outLow1a=null, outLow2a=null;
		String[] outMid1a=null, outMid2a=null;
		String[] outHigh1a=null, outHigh2a=null;
		String[] outUnc1a=null, outUnc2a=null;
		

		final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
		if(outKeep1!=null){
			if(!new File(outKeep1).exists()){
				outKeep1a=outKeep1.split(",");
			}else{
				outKeep1a=new String[] {outKeep1};
			}
			if(outKeep2!=null){
				if(!new File(outKeep2).exists()){
					outKeep2a=outKeep2.split(",");
				}else{
					outKeep2a=new String[] {outKeep2};
				}
			}else{
				outKeep2a=new String[outKeep1a.length];
				for(int i=0; i<outKeep1a.length; i++){
					if(outKeep1a[i].contains("#")){
						outKeep2a[i]=outKeep1a[i].replaceFirst("#", "2");
						outKeep1a[i]=outKeep1a[i].replaceFirst("#", "1");
					}
				}
			}
		}
		if(outToss1!=null){
			if(!new File(outToss1).exists()){
				outToss1a=outToss1.split(",");
			}else{
				outToss1a=new String[] {outToss1};
			}
			if(outToss2!=null){
				if(!new File(outToss2).exists()){
					outToss2a=outToss2.split(",");
				}else{
					outToss2a=new String[] {outToss2};
				}
			}else{
				outToss2a=new String[outToss1a.length];
				for(int i=0; i<outToss1a.length; i++){
					if(outToss1a[i].contains("#")){
						outToss2a[i]=outToss1a[i].replaceFirst("#", "2");
						outToss1a[i]=outToss1a[i].replaceFirst("#", "1");
					}
				}
			}
		}
		if(outLow1!=null){
			if(!new File(outLow1).exists()){
				outLow1a=outLow1.split(",");
			}else{
				outLow1a=new String[] {outLow1};
			}
			if(outLow2!=null){
				if(!new File(outLow2).exists()){
					outLow2a=outLow2.split(",");
				}else{
					outLow2a=new String[] {outLow2};
				}
			}else{
				outLow2a=new String[outLow1a.length];
				for(int i=0; i<outLow1a.length; i++){
					if(outLow1a[i].contains("#")){
						outLow2a[i]=outLow1a[i].replaceFirst("#", "2");
						outLow1a[i]=outLow1a[i].replaceFirst("#", "1");
					}
				}
			}
		}
		if(outMid1!=null){
			if(!new File(outMid1).exists()){
				outMid1a=outMid1.split(",");
			}else{
				outMid1a=new String[] {outMid1};
			}
			if(outMid2!=null){
				if(!new File(outMid2).exists()){
					outMid2a=outMid2.split(",");
				}else{
					outMid2a=new String[] {outMid2};
				}
			}else{
				outMid2a=new String[outMid1a.length];
				for(int i=0; i<outMid1a.length; i++){
					if(outMid1a[i].contains("#")){
						outMid2a[i]=outMid1a[i].replaceFirst("#", "2");
						outMid1a[i]=outMid1a[i].replaceFirst("#", "1");
					}
				}
			}
		}
		if(outHigh1!=null){
			if(!new File(outHigh1).exists()){
				outHigh1a=outHigh1.split(",");
			}else{
				outHigh1a=new String[] {outHigh1};
			}
			if(outHigh2!=null){
				if(!new File(outHigh2).exists()){
					outHigh2a=outHigh2.split(",");
				}else{
					outHigh2a=new String[] {outHigh2};
				}
			}else{
				outHigh2a=new String[outHigh1a.length];
				for(int i=0; i<outHigh1a.length; i++){
					if(outHigh1a[i].contains("#")){
						outHigh2a[i]=outHigh1a[i].replaceFirst("#", "2");
						outHigh1a[i]=outHigh1a[i].replaceFirst("#", "1");
					}
				}
			}
		}
		if(outUnc1!=null){
			if(!new File(outUnc1).exists()){
				outUnc1a=outUnc1.split(",");
			}else{
				outUnc1a=new String[] {outUnc1};
			}
			if(outUnc2!=null){
				if(!new File(outUnc2).exists()){
					outUnc2a=outUnc2.split(",");
				}else{
					outUnc2a=new String[] {outUnc2};
				}
			}else{
				outUnc2a=new String[outUnc1a.length];
				for(int i=0; i<outUnc1a.length; i++){
					if(outUnc1a[i].contains("#")){
						outUnc2a[i]=outUnc1a[i].replaceFirst("#", "2");
						outUnc1a[i]=outUnc1a[i].replaceFirst("#", "1");
					}
				}
			}
		}
		
		long bases=0;
		
		for(int x=0; x<list1.length; x++){
			
			if(outKeep1a!=null){
				if(x==0 || outKeep1a.length>1){
					if(rosKeep!=null){
						rosKeep.close();
						rosKeep.join();
					}
					
					FileFormat ff1=FileFormat.testOutput(outKeep1a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					FileFormat ff2=FileFormat.testOutput(outKeep2a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					rosKeep=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
					
					rosKeep.start();
					outstream.println("Started output threads.");
				}else{
					rosKeep.resetNextListID();
				}
			}
			
			if(outToss1a!=null){
				if(x==0 || outToss1a.length>1){
					if(rosToss!=null){
						rosToss.close();
						rosToss.join();
					}
					
					FileFormat ff1=FileFormat.testOutput(outToss1a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					FileFormat ff2=FileFormat.testOutput(outToss2a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					rosToss=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
					
					rosToss.start();
					outstream.println("Started output threads.");
				}else{
					rosToss.resetNextListID();
				}
			}
			
			if(outLow1a!=null){
				if(x==0 || outLow1a.length>1){
					if(rosLow!=null){
						rosLow.close();
						rosLow.join();
					}
					
					FileFormat ff1=FileFormat.testOutput(outLow1a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					FileFormat ff2=FileFormat.testOutput(outLow2a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					rosLow=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
					
					rosLow.start();
					outstream.println("Started output threads.");
				}else{
					rosLow.resetNextListID();
				}
			}
			
			if(outMid1a!=null){
				if(x==0 || outMid1a.length>1){
					if(rosMid!=null){
						rosMid.close();
						rosMid.join();
					}
					
					FileFormat ff1=FileFormat.testOutput(outMid1a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					FileFormat ff2=FileFormat.testOutput(outMid2a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					rosMid=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
					
					rosMid.start();
					outstream.println("Started output threads.");
				}else{
					rosMid.resetNextListID();
				}
			}
			
			if(outHigh1a!=null){
				if(x==0 || outHigh1a.length>1){
					if(rosHigh!=null){
						rosHigh.close();
						rosHigh.join();
					}
					
					FileFormat ff1=FileFormat.testOutput(outHigh1a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					FileFormat ff2=FileFormat.testOutput(outHigh2a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					rosHigh=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
					
					rosHigh.start();
					outstream.println("Started output threads.");
				}else{
					rosHigh.resetNextListID();
				}
			}
			
			if(outUnc1a!=null){
				if(x==0 || outUnc1a.length>1){
					if(rosUnc!=null){
						rosUnc.close();
						rosUnc.join();
					}
					
					FileFormat ff1=FileFormat.testOutput(outUnc1a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					FileFormat ff2=FileFormat.testOutput(outUnc2a[x], FileFormat.FASTQ, null, true, overwrite, append, ordered);
					rosUnc=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
					
					rosUnc.start();
					outstream.println("Started output threads.");
				}else{
					rosUnc.resetNextListID();
				}
			}
				
			String in1=list1[x];
			String in2=(list2==null || list2.length<=x ? null : list2[x]);

			final ConcurrentReadInputStream cris;
			{
				FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
				FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
				if(verbose){System.err.println("Started cris");}
				cris.start(); //4567
			}
			boolean paired=cris.paired();
			if(verbose){System.err.println("Paired: "+paired);}

			bases+=downsample(cris, kca, k, maxReads, rosKeep, rosToss, rosLow, rosMid, rosHigh, rosUnc, khistFile, rhistFile, peakFile, overwrite, estUnique, storage);

			errorState|=ReadWrite.closeStream(cris);
			if(verbose){System.err.println("Closed stream");}
			
		}

		errorState|=ReadWrite.closeStreams(null, rosKeep, rosToss, rosLow, rosMid, rosHigh, rosUnc);

		return bases;
	}
	

	
	public static long downsample(ConcurrentReadInputStream cris, KCountArray kca, int k, long maxReads,
			ConcurrentReadOutputStream rosKeep, ConcurrentReadOutputStream rosToss, ConcurrentReadOutputStream rosLow, ConcurrentReadOutputStream rosMid, ConcurrentReadOutputStream rosHigh, ConcurrentReadOutputStream rosUnc,
			String khistFile, String rhistFile, String peakFile, boolean overwrite, long estUnique, ArrayList<Read> storage) {
		Timer tdetect=new Timer();
		tdetect.start();

		long totalBases=0;
		long totalReads=0;
		
		long readsKept=0;
		long readsTossed=0;
		long readsLowBin=0;
		long readsMidBin=0;
		long readsHighBin=0;
		long readsUncorrected=0;
		long basesKept=0;
		long basesTossed=0;
		long basesLowBin=0;
		long basesMidBin=0;
		long basesHighBin=0;
		long basesUncorrected=0;

		
		long errorReads=0;
		long errorPairs=0;
		long errorType1=0;
		long errorType2=0;
		long errorType3=0;
		
		long errorsDetected=0;
		long errorsMarked=0;
		long errorsCorrected=0;
		long basesTrimmed=0;
		
		{
			KCountArray kcaup=null;
			if(COUNTUP){
				final int bits;
				if(TARGET_DEPTH<=15){
					bits=4;
				}else if(TARGET_DEPTH<=255){
					bits=8;
				}else{
					bits=16;
				}

				long cells=(FILTERBYTES*8)/bits;
				int kbits=2*k;
				kcaup=KCountArray.makeNew(cells, bits, 3, null, 0);
			}

			ProcessThread[] pta=new ProcessThread[THREADS];
			for(int i=0; i<pta.length; i++){
				pta[i]=new ProcessThread(cris, kca, kcaup, k, rosKeep, rosToss, rosLow, rosMid, rosHigh, rosUnc, storage);
				pta[i].start();
			}

			kca=kcaup=null;

			for(int i=0; i<pta.length; i++){
				ProcessThread ct=pta[i];
				synchronized(ct){
					while(ct.getState()!=State.TERMINATED){
						try {
							ct.join(1000);
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					totalBases+=ct.totalBases;
					totalReads+=ct.totalReads;
					errorReads+=ct.errorReads;
					errorPairs+=ct.errorPairs;
					errorType1+=ct.errorType1;
					errorType2+=ct.errorType2;
					errorType3+=ct.errorType3;

					readsKept+=ct.readsKept;
					readsTossed+=ct.readsTossed;
					readsLowBin+=ct.readsLowBin;
					readsMidBin+=ct.readsMidBin;
					readsHighBin+=ct.readsHighBin;
					readsUncorrected+=ct.readsUncorrected;

					basesKept+=ct.basesKept;
					basesTossed+=ct.basesTossed;
					basesLowBin+=ct.basesLowBin;
					basesMidBin+=ct.basesMidBin;
					basesHighBin+=ct.basesHighBin;
					basesUncorrected+=ct.basesUncorrected;

					errorsDetected+=ct.errorsDetected;
					errorsMarked+=ct.errorsMarked;
					errorsCorrected+=ct.errorsCorrected;
					basesTrimmed+=ct.basesTrimmed;
					errorState|=ct.errorStateT;

					for(int j=0; j<ct.hist.length; j++){
						khistogram.addAndGet(j, ct.hist[j]);
					}
					for(int j=0; j<ct.qhist.length; j++){
						qhist_total[j]+=ct.qhist[j];
					}
				}
			}
		}
		
		if(!ZERO_BIN && khistogram!=null && khistogram.length()>1){
			khistogram.addAndGet(1, khistogram.get(0));
			khistogram.set(0, 0);
		}
		
//		outstream.println();
		tdetect.stop();
		outstream.println("Table read time: \t\t"+tdetect+"   \t"+Tools.format("%.2f", totalBases*1000000.0/(tdetect.elapsed))+" kb/sec");
		
		{
			String pad="";
			String s=""+totalReads;
			while(pad.length()+s.length()<9){pad+=" ";}
			outstream.println("Total reads in:  \t\t"+totalReads+pad+Tools.format("\t%.3f%% Kept", (readsKept*100.0/totalReads)));
			s=""+totalBases;
			while(pad.length()+s.length()<9){pad+=" ";}
			outstream.println("Total bases in:  \t\t"+totalBases+pad+Tools.format("\t%.3f%% Kept", (basesKept*100.0/totalBases)));

			if(rosLow!=null){
				s=""+readsLowBin;
				while(pad.length()+s.length()<9){pad+=" ";}
				outstream.println("Low bin reads:   \t\t"+readsLowBin+pad+Tools.format("\t%.3f%%", (readsLowBin*100.0/totalReads)));
				s=""+basesLowBin;
				while(pad.length()+s.length()<9){pad+=" ";}
				outstream.println("Low bin bases:   \t\t"+basesLowBin+pad+Tools.format("\t%.3f%%", (basesLowBin*100.0/totalBases)));
			}
			if(rosMid!=null){
				s=""+readsMidBin;
				while(pad.length()+s.length()<9){pad+=" ";}
				outstream.println("Mid bin reads:   \t\t"+readsMidBin+pad+Tools.format("\t%.3f%%", (readsMidBin*100.0/totalReads)));
				s=""+basesMidBin;
				while(pad.length()+s.length()<9){pad+=" ";}
				outstream.println("Mid bin bases:   \t\t"+basesMidBin+pad+Tools.format("\t%.3f%%", (basesMidBin*100.0/totalBases)));
			}
			if(rosHigh!=null){
				s=""+readsHighBin;
				while(pad.length()+s.length()<9){pad+=" ";}
				outstream.println("High bin reads:   \t\t"+readsHighBin+pad+Tools.format("\t%.3f%%", (readsHighBin*100.0/totalReads)));
				s=""+basesHighBin;
				while(pad.length()+s.length()<9){pad+=" ";}
				outstream.println("High bin bases:   \t\t"+basesHighBin+pad+Tools.format("\t%.3f%%", (basesHighBin*100.0/totalBases)));
			}
			
			s=""+errorReads;
			while(pad.length()+s.length()<9){pad+=" ";}
			outstream.println("Error reads in:  \t\t"+errorReads+pad+Tools.format("\t%.3f%%", (errorReads*100.0/totalReads)));
			if(cris.paired()){
				s=""+errorPairs;
				while(pad.length()+s.length()<9){pad+=" ";}
				outstream.println("Error pairs in:  \t\t"+errorPairs+pad+Tools.format("\t%.3f%%", (errorPairs*200.0/totalReads)));
			}
			s=""+errorType1;
			while(pad.length()+s.length()<9){pad+=" ";}
			outstream.println("Error type 1:    \t\t"+errorType1+pad+Tools.format("\t%.3f%%", (errorType1*100.0/totalReads)));
			s=""+errorType2;
			while(pad.length()+s.length()<9){pad+=" ";}
			outstream.println("Error type 2:    \t\t"+errorType2+pad+Tools.format("\t%.3f%%", (errorType2*100.0/totalReads)));
			s=""+errorType3;
			while(pad.length()+s.length()<9){pad+=" ";}
			outstream.println("Error type 3:    \t\t"+errorType3+pad+Tools.format("\t%.3f%%", (errorType3*100.0/totalReads)));


			if(TRIM_LEFT_THIS_PASS || TRIM_RIGHT_THIS_PASS){
				outstream.println("\nDuring Trimming:");
				s=""+(errorsDetected+errorsCorrected+errorsMarked);
				while(pad.length()+s.length()<9){pad+=" ";}
				outstream.println("Bases Trimmed:   \t\t"+(basesTrimmed));
			}
			
			if(CORRECT_ERRORS_THIS_PASS){
				outstream.println("\nDuring Error Correction:");
				s=""+(errorsDetected+errorsCorrected+errorsMarked);
				while(pad.length()+s.length()<9){pad+=" ";}
				outstream.println("Errors Suspected:\t\t"+(errorsDetected+errorsCorrected+errorsMarked));
				s=""+errorsCorrected;
				pad="";
				while(pad.length()+s.length()<9){pad+=" ";}
				outstream.println("Errors Corrected:\t\t"+errorsCorrected);
				s=""+errorsMarked;
				pad="";
				while(pad.length()+s.length()<9){pad+=" ";}
				outstream.println("Errors Marked:   \t\t"+errorsMarked+"\n");
			}
		}
		
//		outstream.println();
		if(khistogram!=null){
			
			if(peakFile!=null){
				CallPeaks.printClass=false;
				long[] array=Tools.toArray(khistogram);
				for(int i=0; i<array.length; i++){
					long x=array[i];
					long y=((x+i/2)/(i<1 ? 1 : i)); //x+i/2 rounds to compensate for colliding kmers being put in an overly high bin
					array[i]=y;
				}
				
				ArrayList<String> args=new ArrayList<String>();
				args.add("smoothradius=1");
				args.add("smoothprogressive=t");
				CallPeaks.printPeaks(array, null, peakFile, overwrite, minHeight, minVolume, minWidth, minPeak, maxPeak, maxPeakCount, k, ploidy, doLogScale, logWidth, args);
			}
			
			ByteStreamWriter bsw=null;
			if(USE_KHISTOGRAM && khistFile!=null){
				bsw=new ByteStreamWriter(khistFile, overwrite, false, false);
				bsw.start();

				if(HIST_COLUMNS==1){
					bsw.print("#tUnique_Kmers\n");
				}else if(HIST_COLUMNS==2){
					bsw.print("#Depth\tUnique_Kmers\n");
				}else if(HIST_COLUMNS==3){
					bsw.print("#Depth\tRaw_Count\tUnique_Kmers\n");
				}
				
			}
			int lim=(int)(HIST_LEN_PRINT-1);
			long remaining=Tools.sum(khistogram);
			long sumRaw1=0;
			long sumRaw2=0;
			long sum1=0;
			long sum2=0;
			long sumsquare=0;
			for(int i=0; i<lim; i++){
				long x=khistogram.get(i);
				long y=((x+i/2)/(i<1 ? 1 : i)); //x+i/2 rounds to compensate for colliding kmers being put in an overly high bin
//				long y=((x)/(i<1 ? 1 : i));
				sumRaw1+=x;
				sum1+=y;
				sumsquare+=(x*Tools.max(1, i));
				if(bsw!=null){
					if(PRINT_ZERO_COVERAGE /*|| x>0*/ || y>0 || HIST_COLUMNS==1){
						if(HIST_COLUMNS>1){
							bsw.print(i);
							bsw.print('\t');
						}
						if(HIST_COLUMNS==3){
							bsw.print(x);
							bsw.print('\t');
						}
						bsw.print(y);
						bsw.print('\n');
					}
				}
				if(sumRaw1>=remaining){break;} //Stop once there is no more coverage, even if PRINT_ZERO_COVERAGE is not set.
			}
			for(int i=lim; i<khistogram.length(); i++){
				long x=khistogram.get(i);
				sumRaw2+=x;
				long y=((x+i/2)/(i<1 ? 1 : i)); //x+i/2 rounds to compensate for colliding kmers being put in an overly high bin
//				long y=((x)/(i<1 ? 1 : i));
				sum2+=y;
			}
			if(bsw!=null){
				if(sumRaw2>0 || sum2>0){
					if(HIST_COLUMNS>1){
						bsw.print(lim);
						bsw.print('\t');
					}
					if(HIST_COLUMNS==3){
						bsw.print(sumRaw2);
						bsw.print('\t');
					}
					bsw.print(sum2);
					bsw.print('\n');
				}
				bsw.poisonAndWait();
				outstream.println("\nWrote histogram to "+khistFile);
			}
			
			long histCount=Tools.sum(khistogram); //Total number of kmers counted
			long halfCount=(histCount+1)/2;
			double histCountU=0; //Unique kmers counted
			long temp1=0;
			double temp2=0;
			int median_all=-1;
			int median_unique=-1;
			for(int i=0; i<khistogram.length(); i++){
				long x=khistogram.get(i);
				temp1+=x;
				if(temp1>=halfCount && median_all<0){median_all=i;}
//				histSum+=(x*(double)i);
				histCountU+=(x/(double)Tools.max(1, i));
			}
			double halfCount2=(histCountU)/2;
			for(int i=0; i<khistogram.length(); i++){
				long x=khistogram.get(i);
				temp2+=(x/Tools.max(i, 1.0));
				if(temp2>=halfCount2 && median_unique<0){
					median_unique=i;
					break;
				}
			}
			if(median_all<0){median_all=0;}
			double avg_all=sumsquare/(double)histCount;
			double avg_unique=histCount/histCountU;
			double stdev_unique=Tools.standardDeviationHistogramKmer(khistogram);
			double stdev_all=Tools.standardDeviationHistogram(khistogram);
			outstream.println("Total kmers counted:          \t"+(sumRaw1+sumRaw2));
			
			double uniqueC=((sum1+sum2)*100.0/(sumRaw1+sumRaw2));
			double uniqueE=((estUnique)*100.0/(sumRaw1+sumRaw2));
			double uniqueM=Tools.max(uniqueC, uniqueE);
			outstream.println("Total unique kmer count:      \t"+(sum1+sum2));
			if(CANONICAL){outstream.println("Includes forward kmers only.");}
			outstream.println("The unique kmer estimate can be more accurate than the unique count, if the tables are very full.");
			outstream.println("The most accurate value is the greater of the two.");
			outstream.println();
			
			outstream.println("Percent unique:               \t"+(uniqueM<10 ? " " : "")+Tools.format("%.2f%%", uniqueM));

			outstream.println("Depth average:                \t"+Tools.format("%.2f\t(unique kmers)", avg_unique));
			outstream.println("Depth median:                 \t"+Tools.format("%d\t(unique kmers)", median_unique));
			outstream.println("Depth standard deviation:     \t"+Tools.format("%.2f\t(unique kmers)", stdev_unique));
			outstream.println("Corrected depth average:      \t"+Tools.format("%.2f\t", Tools.observedToActualCoverage(avg_unique)));
			
			outstream.println("\nDepth average:                \t"+Tools.format("%.2f\t(all kmers)", avg_all));
			outstream.println("Depth median:                 \t"+Tools.format("%d\t(all kmers)", median_all));
			outstream.println("Depth standard deviation:     \t"+Tools.format("%.2f\t(all kmers)", stdev_all));
			
			double avgReadLen=totalBases*1.0/totalReads;
			double readDepth=median_all*(avgReadLen/(avgReadLen-k+1));
			
			outstream.println("\nApprox. read depth median:    \t"+Tools.format("%.2f", (readDepth)));
		}
		
		
		if(rhistogram!=null){
			TextStreamWriter tswh=null;
			StringBuilder sb=new StringBuilder(100);
			if(USE_RHISTOGRAM && rhistFile!=null){
				tswh=new TextStreamWriter(rhistFile, overwrite, false, false);
				tswh.start();
				tswh.print("#Depth\tReads\tBases\n");
			}
			int lim=(int)(HIST_LEN_PRINT-1);
			long remaining=Tools.sum(rhistogram);
			long sumReads1=0;
			long sumReads2=0;
			long sumBases1=0;
			long sumBases2=0;
			long sumSquareReads=0;
			long sumSquareBases=0;
			for(int i=0; i<lim; i++){
				final long x=rhistogram.get(i);
				final long y=bhistogram.get(i);
				sumReads1+=x;
				sumBases1+=y;
				sumSquareReads+=(x*Tools.max(1, i));
				sumSquareBases+=(y*Tools.max(1, i));
				if(tswh!=null){
					if(PRINT_ZERO_COVERAGE /*|| x>0*/ || y>0){
						sb.append(i).append('\t');
						sb.append(x).append('\t');
						sb.append(y).append('\n');
					}
					tswh.print(sb.toString());
					sb.setLength(0);
				}
				if(sumReads1>=remaining){break;} //Stop once there is no more coverage, even if PRINT_ZERO_COVERAGE is not set.
			}
			for(int i=lim; i<rhistogram.length(); i++){
				final long x=rhistogram.get(i);
				final long y=bhistogram.get(i);
				sumReads2+=x;
				sumBases2+=y;
			}
			if(tswh!=null){
				if(sumReads2>0 || sumBases2>0){
					sb.append(lim).append('\t');
					sb.append(sumReads2).append('\t');
					sb.append(sumBases2).append('\n');
				}
				tswh.print(sb.toString());
				tswh.poison();
				tswh.waitForFinish();
				outstream.println("\nWrote histogram to "+rhistFile);
			}

			long rhistCount=Tools.sum(rhistogram); //Total number of reads counted
			long bhistCount=Tools.sum(bhistogram); //Total number of bases counted
			int median_reads=-1;
			int median_bases=-1;
			{
				long halfCount=(rhistCount+1)/2;
				long temp=0;
				for(int i=0; i<rhistogram.length(); i++){
					long x=rhistogram.get(i);
					temp+=x;
					if(temp>=halfCount && median_reads<0){median_reads=i;}
				}
				if(median_reads<0){median_reads=0;}
			}
			{
				long halfCount=(bhistCount+1)/2;
				long temp=0;
				for(int i=0; i<bhistogram.length(); i++){
					long x=bhistogram.get(i);
					temp+=x;
					if(temp>=halfCount && median_bases<0){median_bases=i;}
				}
				if(median_bases<0){median_bases=0;}
			}
			double avg_reads=sumSquareReads/(double)rhistCount;
			double avg_bases=sumSquareBases/(double)bhistCount;
//			double read_stdev_unique=Tools.standardDeviationHistogramKmer(rhistogram_total);
			double read_stdev_all=Tools.standardDeviationHistogram(rhistogram);
			double base_stdev_all=Tools.standardDeviationHistogram(bhistogram);
			outstream.println("Total reads counted:          \t"+(sumReads1+sumReads2));
			
//			double uniqueC=((sumBases1+sumBases2)*100.0/(sumReads1+sumReads2));
//			double uniqueE=((estUnique)*100.0/(sumReads1+sumReads2));
//			double uniqueM=Tools.max(uniqueC, uniqueE);
			outstream.println("Total bases counted:          \t"+(sumBases1+sumBases2));

			outstream.println("Read depth average:           \t"+Tools.format("%.2f", avg_reads));
			outstream.println("Read depth median:            \t"+Tools.format("%d", median_reads));
			outstream.println("Read depth standard deviation:\t"+Tools.format("%.2f", read_stdev_all));
			
			outstream.println("\nBase depth average:           \t"+Tools.format("%.2f)", avg_bases));
			outstream.println("Base depth median:            \t"+Tools.format("%d", median_bases));
			outstream.println("Base depth standard deviation:\t"+Tools.format("%.2f", base_stdev_all));
		}
		
		if(errorState){throw new RuntimeException("BBNorm terminated in an error state; the output may be corrupt.");}
		
		return totalBases;
	}
	
	
	
	/**
	 * Locates and fixes spikes in a coverage profile (potentially) caused by false positives in a bloom filter.
	 * Theory:  If a high-count kmer is adjacent on both sides to low-count kmers, it may be a false positive.
	 * It could either be reduced to the max of the two flanking points or examined in more detail.
	 * @param cov An array of kmer counts for adjacent kmers in a read.
	 */
	private static void fixSpikes(int[] cov){
		
		for(int i=1; i<cov.length-1; i++){
			long a=Tools.max(1, cov[i-1]);
			int b=cov[i];
			long c=Tools.max(1, cov[i+1]);
			if(b>1 && b>a && b>c){
				//peak
				if((b>=2*a || b>a+2) && (b>=2*c || b>c+2)){
					//spike
					cov[i]=(int)Tools.max(a, c);
				}
			}
		}
	}
	
	private static void fixSpikes(int[] cov, long[] kmers, KCountArray kca, final int k){
		assert(k<32) : "this function not tested with k>31";
		if(cov.length<3){return;}
		if(cov[1]-cov[0]>1){
			cov[0]=kca.readPrecise(kmers[0], k, true);
		}
		if(cov[cov.length-1]-cov[cov.length-2]>1){
			cov[cov.length-1]=kca.readPrecise(kmers[cov.length-1], k, true);
		}
		
		for(int i=1; i<cov.length-1; i++){
			int b=cov[i];
			if(b>1){
				long a=Tools.max(1, cov[i-1]);
				long c=Tools.max(1, cov[i+1]);
				long key=kmers[i];

				if(b>a && b>c){
					//peak
					if(b<6 || b>a+1 || b>c+1){
						cov[i]=kca.readPreciseMin(key, k, true);
					}
					//				if((b>=2*a || b>a+2) && (b>=2*c || b>c+2)){
					//					//spike
					//					int b1=(int)((a+c)/2);
					//					int b2=kca.readLeft(key, k, CANONICAL);
					//					int b3=kca.readRight(key, k, CANONICAL);
					//					array[i]=Tools.min(b, b1, b2, b3);
					//				}
					//				else
					//				{
					////					array[i]=kca.readPreciseMin(key, k, CANONICAL);
					//				}
				}
				//			else
				//				if(Tools.max(ada, adc)>=Tools.max(2, Tools.min((int)a, b, (int)c)/4)){
				//					array[i]=kca.readPrecise(key, k, CANONICAL);
				//				}
				//			else
				//				if(b>a+1 || b>c+1){
				//					//steep
				//					array[i]=kca.readPrecise(key, k, CANONICAL);
				//				}
			}
		}
	}


	private static int correctErrors(Read r, int[] cov, long[] kmers, KCountArray kca, final int k,
			final int low, final int high, final int mult, int maxToCorrect, int maxQual, boolean kmersAlreadyValid, boolean coverageAlreadyValid, long[] qhist,
			final boolean markOnly, Kmer longkmer){
		assert(k<32) : "this function not tested with k>31";
		assert(maxToCorrect>0) : "Don't do error correction with a maximum of 0 errors; it's a waste of time.";
		if(maxToCorrect<1){return 0;}
		
		if(!kmersAlreadyValid){kmers=r.toKmers(k, 0, kmers, false, longkmer);}
		if(kmers==null || kmers.length<3){return -99;}
		
		if(!coverageAlreadyValid){cov=generateCoverage(kca, k, cov, kmers, true);}
		
		int disc=countDiscontinuities(cov, low, high, mult);
		
		if(disc==0){return 0;}
		
		byte[] copy=r.bases.clone();
		
		byte[] suffix=new byte[SUFFIX_LEN];
		
		int cfl=0, cfr=0;
		
		if(CORRECT_FROM_LEFT){
			cfl=correctErrorsFromLeft(r, cov, kmers, kca, k, low, high, mult, maxToCorrect, maxQual, suffix, qhist, markOnly, longkmer);
			if(cfl<0){
				//Failed correction.
				r.bases=copy;
				return cfl;
			}
			maxToCorrect-=cfl;
		}
		
		if(CORRECT_FROM_RIGHT && maxToCorrect>0){
			{//Optional block - allows saving of errors corrected from left even if correctErrorsFromRight fails.
				if(cfl>0){
					for(int i=0; i<r.length(); i++){
						copy[i]=r.bases[i];
					}
				}
			}
			cfr=correctErrorsFromRight(r, cov, kmers, kca, k, low, high, mult, maxToCorrect, maxQual, suffix, qhist, markOnly, longkmer);
			if(cfr<0){
				//Failed correction.
				r.bases=copy;
				return cfr;
			}
		}
		
		return cfl+cfr;
	}

	private static int markErrors(Read r, int[] cov, long[] kmers, KCountArray kca, final int k,
			final int low, final int high, final int mult, int maxToCorrect, boolean kmersAlreadyValid, boolean coverageAlreadyValid, long[] qhist, Kmer longkmer){
		assert(k<32) : "this function not tested with k>31";
		assert(maxToCorrect>0) : "Don't do error correction with a maximum of 0 errors; it's a waste of time.";
		if(maxToCorrect<1){return 0;}
		
		if(!kmersAlreadyValid){kmers=r.toKmers(k, 0, kmers, false, longkmer);}
		if(kmers==null || kmers.length<3){return 0;}
		
		if(!coverageAlreadyValid){cov=generateCoverage(kca, k, cov, kmers, true);}
		
		int disc=countDiscontinuities(cov, low, high, mult);
		
		if(disc==0){return 0;}
		
		int cfl=0, cfr=0;
		
		if(CORRECT_FROM_LEFT){
			cfl=markErrorsFromLeft(r, cov, k, low, high, mult, maxToCorrect, qhist);
			maxToCorrect-=cfl;
		}
		
		if(CORRECT_FROM_RIGHT){
			cfr=markErrorsFromRight(r, cov, k, low, high, mult, maxToCorrect, qhist);
		}
		
		int marked=cfl+cfr;
		final byte[] quals=r.quality;
		if(marked>0){
			int found=0;
			if(quals!=null){
				for(int i=0; i<quals.length; i++){
					byte q=quals[i];
					if(q<0){
						byte q2=(byte)(MARK_WITH_1 ? 1 : Tools.max(1, -(q/2+3)));
						quals[i]=q2;
						found++;
					}
				}
			}
			assert(found==marked);
		}else{
			if(quals!=null){
				for(int i=0; i<quals.length; i++){
					assert(quals[i]>=0);
				}
			}
			assert(marked==0) : marked;
		}
		
		return marked;
	}
	
	/** Returns number of discontinuities detected.  This is not the same as the number of errors,
	 * but the presence of discontinuities indicates the presence of errors.
	 * @param cov
	 * @param low
	 * @param high
	 * @param mult
	 * @return
	 */
	private static int countDiscontinuities(final int[] cov, final int low, final int high, final int mult){
		
		int found=0;
		
		for(int i=2; i<cov.length; i++){
			int a=Tools.min(cov[i-2], cov[i-1]);
			int b=cov[i];
			if(a>=high && (b<=low || a>=b*mult)){//error
				found++;
			}
		}
		
		for(int i=cov.length-3; i>=0; i--){
			int a=Tools.min(cov[i+2], cov[i+1]);
			int b=cov[i];
			if(a>=high && (b<=low || a>=b*mult)){//error
				found++;
			}
		}
		
		return found;
	}
	
	
	private static void regenerateKmersAndCoverage(final Read r, final long[] kmers, final int[] cov, final KCountArray kca, final int k, boolean makeCanonical,
			Kmer longkmer){
		assert(r!=null && kmers!=null && cov!=null && kca!=null);
		final byte[] bases=r.bases;
		if(bases==null || bases.length<k){return;}
		
		if(k>31){
			r.toKmers(k, 0, kmers, false, longkmer);
			generateCoverage(kca, k, cov, kmers, true);
			return;
		}
		
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		
		int len=0;
		long kmer=0;
		final int arraylen=bases.length-k+1;
		assert(kmers.length==arraylen && cov.length==arraylen);
		
		for(int i=0, j=1-k; i<bases.length; i++, j++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
//			assert(x>=0) : "This program does not allow degenerate bases other than N.  Invalid symbol: ASCII character "+b+" ("+(char)(b<33 ? ' ' : b)+")";
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;

				if(len>=k){
					long y=kmer;
					if(makeCanonical){
						y=KCountArray.makeCanonical2(y, k);
					}
					
					if(kmers[j]!=y){
						kmers[j]=y;
						cov[j]=kca.read(y, k, !makeCanonical);
					}
				}
			}
		}
	}
	
	
	private static int correctErrorsFromLeft(final Read r, final int[] cov, final long[] kmers, final KCountArray kca, final int k,
			final int low, final int high, final int mult, final int maxToCorrect, int maxQual, final byte[] suffix, final long[] qhist, boolean markOnly,
			Kmer longkmer){
		
		int found=0;
		int corrected=0;
		int uncorrected=0;
		final byte[] quals=r.quality;
		
		for(int i=PREFIX_LEN; i<cov.length; i++){
//			int a=Tools.min(cov[i-2], cov[i-1]);
			final int a=Tools.min(cov, i-PREFIX_LEN, i-1);
			final int b=cov[i];
			if(a>=high && (b<=low || a>=b*mult)){//error
				found++;
				final int loc=i+k-1;
				final byte q=(quals==null ? 10 : quals[loc]);
				if(qhist!=null){qhist[q]++;}
				
				if(markOnly){
					corrected++;
					if(quals==null){r.bases[loc]='N';}
//					else if(q>0){quals[loc]=(byte)Tools.max(1, q/2);}
					else if(q>0){quals[loc]=(byte)Tools.max(1, q/2-3);}
				}else{
					if(found>maxToCorrect || q>maxQual){return 0-found;}
					boolean success=correctErrorFromLeft(r, cov, kmers, kca, k, low, Tools.max(high, a/2), 2*a, mult, i, suffix);
					if(success){
						corrected++;
						//					r.toKmers(k, 0, kmers, false);
						//					generateCoverage(kca, k, cov, kmers, true);
						regenerateKmersAndCoverage(r, kmers, cov, kca, k, false, longkmer);
					}else{
						uncorrected++;
						break;
					}
				}
			}
		}
		
//		assert(false) : Arrays.toString(cov)+"\nlow="+low+", high="+high+", mult="+mult+", found="+found+", corrected="+corrected+", uncorrected="+uncorrected;
		
		
		return (uncorrected>0 ? 0-found : corrected);
	}
	
	
	private static int correctErrorsFromRight(final Read r, final int[] cov, final long[] kmers, final KCountArray kca, final int k,
			final int low, final int high, final int mult, final int maxToCorrect, int maxQual, final byte[] suffix, final long[] qhist, final boolean markOnly,
			Kmer longkmer){
		
		int found=0;
		int corrected=0;
		int uncorrected=0;
		final byte[] quals=r.quality;
		
		final int start=(markOnly ? Tools.min(cov.length-PREFIX_LEN-1, k-1) : cov.length-PREFIX_LEN-1);
		for(int i=start; i>=0; i--){
//			int a=Tools.min(cov[i+2], cov[i+1]);
			int a=Tools.min(cov, i+1, i+PREFIX_LEN);
			int b=cov[i];
			if(a>=high && (b<=low || a>=b*mult)){//error
				found++;
				final byte q=(quals==null ? 10 : quals[i]);
				if(qhist!=null){qhist[q]++;}
				
				if(markOnly){
					corrected++;
					if(quals==null){r.bases[i]='N';}
//					else if(q>0){quals[i]=(byte)Tools.max(1, q/2);}
					else if(q>0){quals[i]=(byte)Tools.max(1, q/2-3);}
				}else{
					if(found>maxToCorrect || q>maxQual){return 0-found;}
					boolean success=correctErrorFromRight(r, cov, kmers, kca, k, low, Tools.max(high, a/2), 2*a, mult, i, suffix);
					if(success){
						corrected++;
						//					r.toKmers(k, 0, kmers, false);
						//					generateCoverage(kca, k, cov, kmers, true);
						regenerateKmersAndCoverage(r, kmers, cov, kca, k, false, longkmer);
					}else{
						uncorrected++;
						break;
					}
				}
			}
		}
		
//		assert(false) : Arrays.toString(cov)+"\nlow="+low+", high="+high+", mult="+mult+", found="+found+", corrected="+corrected+", uncorrected="+uncorrected;
		
		
		return (uncorrected>0 ? 0-found : corrected);
	}
	
	
	private static int markErrorsFromLeft(final Read r, final int[] cov, final int k,
			final int low, final int high, final int mult, final int maxToCorrect, final long[] qhist){
		
		int found=0;
		final byte[] quals=r.quality, bases=r.bases;
		
		for(int i=PREFIX_LEN; i<cov.length; i++){
			final int a=Tools.min(cov, i-PREFIX_LEN, i-1);
			final int b=cov[i];
			if(a>=high && (b<=low || a>=b*mult)){//error
				final int loc=i+k-1;
				final byte q=(quals==null ? 10 : quals[loc]);
				
				if(q>0){
					found++;
					if(qhist!=null){qhist[q]++;}
					if(quals==null){bases[loc]='N';}
					else{quals[loc]=(byte)-q;}
				}
			}
		}
		return found;
	}
	
	
	private static int markErrorsFromRight(final Read r, final int[] cov, final int k,
			final int low, final int high, final int mult, final int maxToCorrect, final long[] qhist){
		
		int found=0;
		final byte[] quals=r.quality, bases=r.bases;
		
		final int start=cov.length-PREFIX_LEN-1;
		for(int i=start; i>=0; i--){
			int a=Tools.min(cov, i+1, i+PREFIX_LEN);
			int b=cov[i];
			if(a>=high && (b<=low || a>=b*mult)){//error
				final byte q=(quals==null ? 10 : quals[i]);
				
				if(q>0){
					found++;
					if(qhist!=null){qhist[q]++;}
					if(quals==null){bases[i]='N';}
					else{quals[i]=(byte)-q;}
				}
			}
		}
		
		return found;
	}
	
	private static boolean correctErrorFromLeft(final Read r, final int[] cov, final long[] kmers, final KCountArray kca, final int k,
			final int low, final int targetLowerBound, final int targetUpperBound, final int mult, final int loc, final byte[] suffix){
		
		for(int i=0, j=loc+k-1; i<suffix.length; i++, j++){
			if(j<r.length()){
				suffix[i]=r.bases[j];
			}else{
				suffix[i]='N';
			}
		}
		
//		if(r.numericID!=3500){return false;}
		
		
		long kmer=kmers[loc];
		final boolean defined=(AminoAcid.isFullyDefined(suffix[0]));
		
		//This block added to allow correction of no-calls
		if(!defined && loc>0){
			assert(kmer==-1L) : new String(suffix)+"\t"+kmer;
			if(kmers[loc-1]!=-1L){
				final int kbits=2*k;
				final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
				kmer=((kmers[loc-1]<<2)&mask);
			}
		}
//		int leftCov=Tools.min(cov[loc-1], cov[loc-2]);

//		assert(false) : "kmer = "+AminoAcid.kmerToString(kmers[0], k);
//		assert(false) : "suffix = "+new String(suffix);
		
//		assert(false) : new String(suffix)+"\t"+kmer;
		
		suffix[0]='A';
		final int a=testRightSuffix(kca, k, kmer, suffix);
		suffix[0]='C';
		final int c=testRightSuffix(kca, k, kmer, suffix);
		suffix[0]='G';
		final int g=testRightSuffix(kca, k, kmer, suffix);
		suffix[0]='T';
		final int t=testRightSuffix(kca, k, kmer, suffix);
		
		final int max=Tools.max(a, c, g, t);
		byte best='N';
		
//		assert(false) : "rid="+r.numericID+"\n"+Arrays.toString(cov)+"\n" +
//				new String(r.bases)+"\n" +
//				"loc="+loc+", "+new String(suffix)+"\n" +
//				"low="+low+", high="+high+", mult="+mult+", a="+a+", c="+c+", g="+g+", t="+t+", max="+max;
		
		if(max>=targetLowerBound && max<=targetUpperBound){
			//Found correct answer!
			final int max2;
			if(a==max){
				max2=Tools.max(c, g, t);
				best='A';
			}else if(c==max){
				max2=Tools.max(a, g, t);
				best='C';
			}else if(g==max){
				max2=Tools.max(a, c, t);
				best='G';
			}else if(t==max){
				max2=Tools.max(a, c, g);
				best='T';
			}else{
				max2=max;
				assert(false);
			}

//			assert(false) : max+", "+max2+", "+low+", "+(char)best;
			if(max2<=low || max2*mult<=max){
				final int bnum=loc+k-1;
				r.bases[bnum]=best;
				if(!defined && r.quality!=null){
					assert(r.quality[bnum]==0) : r;
					r.quality[bnum]=FIXED_N_QUAL;
				}
				return true;
			}
		}
		
//		assert(false) : max+", "+targetLowerBound+", "+targetUpperBound+", "+low+", "+(char)best;
		
		return false;
	}
	
	private static boolean correctErrorFromRight(final Read r, final int[] cov, final long[] kmers, final KCountArray kca, final int k,
			final int low, final int targetLowerBound, final int targetUpperBound, final int mult, final int loc, final byte[] suffix){
		
		for(int i=0, j=loc; i<suffix.length; i++, j--){
			if(j>=0){
				suffix[i]=r.bases[j];
			}else{
				suffix[i]='N';
			}
		}
//		if(r.numericID!=3500){return false;}
		
		long kmer=kmers[loc];
		final boolean defined=(AminoAcid.isFullyDefined(suffix[0]));
		
		//This block added to allow correction of no-calls
		if(!defined && loc+1<kmers.length){
			assert(kmer==-1L) : new String(suffix)+"\t"+kmer;
			if(kmers[loc+1]!=-1L){
				final int kbits=2*k;
				final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
				kmer=((kmers[loc+1]>>2)&mask);
			}
		}
//		int rightCov=Tools.min(cov[loc+1], cov[loc+2]);

//		assert(false) : "kmer = "+AminoAcid.kmerToString(kmers[0], k);
//		assert(false) : "suffix = "+new String(suffix);
		
		suffix[0]='A';
		final int a=testLeftSuffix(kca, k, kmer, suffix);
		suffix[0]='C';
		final int c=testLeftSuffix(kca, k, kmer, suffix);
		suffix[0]='G';
		final int g=testLeftSuffix(kca, k, kmer, suffix);
		suffix[0]='T';
		final int t=testLeftSuffix(kca, k, kmer, suffix);
		
		final int max=Tools.max(a, c, g, t);
		byte best='N';
		
//		assert(false) : "\nrid="+r.numericID+"\n"+Arrays.toString(cov)+"\n" +
//				new String(r.bases)+"\n"+
//				"kmer-2 = "+AminoAcid.kmerToString(kmers[loc-2], k)+"\n"+
//				"kmer-1 =  "+AminoAcid.kmerToString(kmers[loc-1], k)+"\n"+
//				"kmer   =   "+AminoAcid.kmerToString(kmer, k)+"\n"+
//				"kmer+1 =    "+AminoAcid.kmerToString(kmers[loc+1], k)+"\n"+
//				"kmer+2 =     "+AminoAcid.kmerToString(kmers[loc+2], k)+"\n"+
//				"count=("+kca.read(kmers[loc-2], k, true)+", "+kca.read(kmers[loc-1], k, true)+", "+
//				kca.read(kmer, k, true)+", "+kca.read(kmers[loc+1], k, true)+", "+kca.read(kmers[loc+2], k, true)+")\n"+
//				"loc="+loc+", suffix="+new String(suffix)+"\n" +
//				"low="+low+", high="+high+", mult="+mult+", a="+a+", c="+c+", g="+g+", t="+t+", max="+max;
		
		if(max>=targetLowerBound && max<=targetUpperBound){
			//Found correct answer!
			final int max2;
			if(a==max){
				max2=Tools.max(c, g, t);
				best='A';
			}else if(c==max){
				max2=Tools.max(a, g, t);
				best='C';
			}else if(g==max){
				max2=Tools.max(a, c, t);
				best='G';
			}else if(t==max){
				max2=Tools.max(a, c, g);
				best='T';
			}else{
				max2=max;
				assert(false);
			}
			
			if(max2<=low || max2*mult<=max){
				r.bases[loc]=best;
				if(!defined && r.quality!=null){
					assert(r.quality[loc]==0) : r;
					r.quality[loc]=FIXED_N_QUAL;
				}
				return true;
			}
		}
		return false;
	}
	
	private static int testRightSuffix(final KCountArray kca, final int k, final long kmer0, final byte[] suffix){
		assert(k<=31);
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		
		long kmer=kmer0>>2;
		int min=Integer.MAX_VALUE;

//		System.out.println("Processing suffix "+new String(suffix));
//		System.out.println("kmer = "+AminoAcid.kmerToString(kmer0, k));
//		System.out.println("cov = "+kca.read(kmer0, k, true));
		
		
		for(int i=0; i<suffix.length && min>0; i++){
			byte b=suffix[i];
			long x=AminoAcid.baseToNumber[b];
			if(x<0){
				//TODO: Find best next letter
				return 0;
			}
			
			kmer=((kmer<<2)|x)&mask;
			int cov=kca.read(kmer, k, true);
			min=Tools.min(min, cov);

//			System.out.println("kmer = "+AminoAcid.kmerToString(kmer, k));
//			System.out.println("cov = "+cov);
		}
//		System.out.println("returning "+min);
		
		assert(min<Integer.MAX_VALUE);
		return min;
	}
	
	private static int testLeftSuffix(final KCountArray kca, final int k, final long kmer0, final byte[] suffix){
		assert(k<=31);
		final int kbits=2*k;
		final int shift=kbits-2;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		
		long kmer=(kmer0<<2)&mask;
		int min=Integer.MAX_VALUE;

//		System.out.println("Processing suffix "+new String(suffix));
//		System.out.println("kmer = "+AminoAcid.kmerToString(kmer0, k));
//		System.out.println("cov = "+kca.read(kmer0, k, true));
		
		
		for(int i=0; i<suffix.length && min>0; i++){
			byte b=suffix[i];
			long x=AminoAcid.baseToNumber[b];
			if(x<0){
				//TODO: Find best next letter
				return 0;
			}

			kmer=((kmer>>2)|(x<<shift));
			int cov=kca.read(kmer, k, true);
			min=Tools.min(min, cov);

//			System.out.println("kmer = "+AminoAcid.kmerToString(kmer, k));
//			System.out.println("cov = "+cov);
		}
//		System.out.println("returning "+min);
		
		assert(min<Integer.MAX_VALUE);
		return min;
	}
	
	private static void analyzeSpikes(int[] array, int width){
		if(array.length<3){return;}
		int peakcount=0, valleycount=0, spikecount=0, flatcount=0, slopecount=0;
		for(int i=1; i<array.length-1; i++){
			long a=array[i-1];
			int b=array[i];
			long c=array[i+1];
			if(b>a && b>c){
				peakcount++;
				if((b>=2*a || b>a+2) && (b>=2*c || b>c+2)){
					spikecount++;
				}
			}else if(b<a && b<c){
				valleycount++;
			}else if(b==a && b==c){
				flatcount++;
			}else{
				slopecount++;
			}
		}
		if(peakcount>0){peaks.addAndGet(peakcount);}
		if(valleycount>0){valleys.addAndGet(valleycount);}
		if(spikecount>0){spikes.addAndGet(spikecount);}
		if(flatcount>0){flats.addAndGet(flatcount);}
		if(slopecount>0){slopes.addAndGet(slopecount);}
	}
	
	/**
	 * kmer array must be valid at this point
	 * @param r
	 * @param kca
	 * @param k
	 * @param out
	 * @param kmers
	 * @return Array of coverage per kmer
	 */
	public static int[] generateCoverage(Read r, KCountArray kca, final int k, int[] out, long[] kmers){
		
		assert(kmers!=null);
		if(kmers==null){return null;} //Read is too short
		
		out=generateCoverage(kca, k, out, kmers, k<=31);
		
		if(ANALYZE_TOPOLOGY){analyzeSpikes(out, 1);}
		return out;
	}
	
	/**
	 * kmer array must be valid at this point
	 * @param kca
	 * @param k
	 * @param out
	 * @param kmers
	 * @param makeCanonical
	 * @return Array of coverage per kmer
	 */
	public static int[] generateCoverage(KCountArray kca, int k, int[] out, long[] kmers, boolean makeCanonical){
		if(kmers==null){return null;}
		
		if(out==null || out.length!=kmers.length){out=new int[kmers.length];}
		Arrays.fill(out, -1);
		
		for(int i=0; i<kmers.length; i++){
			long kmer=kmers[i];
			if(kmer!=-1){
				int count=kca.read(kmer, k, makeCanonical);
				out[i]=count;
			}
		}
		
		if(FIX_SPIKES){fixSpikes(out, kmers, kca, k);}
		return out;
	}
	
	/** Returns {depth1, depth2, errors1, errors2} */
	public static int[] parseDepth(String s, int[] array){
		if(s==null || !s.startsWith("id=")){return null;}
		if(array==null){array=KillSwitch.allocInt1D(4);}
		String[] split=s.split("[, ]");
		Arrays.fill(array, -1);
//		assert(false) : s+"\n"+Arrays.toString(split);
		try {
			for(int i=1; i<split.length; i++){
				final String ss=split[i];
				if(ss.startsWith("d1=")){array[0]=Integer.parseInt(ss.substring(3));}
				else if(ss.startsWith("d2=")){array[2]=Integer.parseInt(ss.substring(3));}
				else if(ss.startsWith("e1=")){array[3]=Integer.parseInt(ss.substring(3));}
				else if(ss.startsWith("e2=")){array[4]=Integer.parseInt(ss.substring(3));}
			}
			return array;
		} catch (NumberFormatException e) {
			return null;
		}
	}
	
	
	private static class ProcessThread extends Thread{
		
		ProcessThread(ConcurrentReadInputStream cris_, KCountArray kca_, KCountArray kcaup_, int k_,
				ConcurrentReadOutputStream rosk_, ConcurrentReadOutputStream rost_, ConcurrentReadOutputStream rosl_, ConcurrentReadOutputStream rosm_, ConcurrentReadOutputStream rosh_, ConcurrentReadOutputStream rosu_,
				ArrayList<Read> storage_){
			cris=cris_;
			kca=kca_;
			kcaup=kcaup_;
			k=k_;
			rosk=rosk_;
			rost=rost_;
			rosl=rosl_;
			rosm=rosm_;
			rosh=rosh_;
			rosu=rosu_;
			storage=storage_;
		}
		
		@Override
		public void run(){
			errorStateT=true;
			randy=Shared.threadLocalRandom();
			if(COUNTUP){
				normalizeInThreadByCountup();
			}else{
				normalizeInThread();
			}
			errorStateT=false;
		}
		
		void normalizeInThread() {
			
			Kmer longkmer=(k<32 ? null : new Kmer(k));
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			final ArrayList<Read> keepList=(rosk==null ? null : new ArrayList<Read>(Shared.bufferLen()));
			final ArrayList<Read> tossList=(rost==null ? null : new ArrayList<Read>(Shared.bufferLen()));
			final ArrayList<Read> lowList=(rosl==null ? null : new ArrayList<Read>(Shared.bufferLen()));
			final ArrayList<Read> midList=(rosm==null ? null : new ArrayList<Read>(Shared.bufferLen()));
			final ArrayList<Read> highList=(rosh==null ? null : new ArrayList<Read>(Shared.bufferLen()));
			final ArrayList<Read> uncList=(rosu==null ? null : new ArrayList<Read>(Shared.bufferLen()));
			
			int[] cov1=null, cov2=null;
			long[] kmers1=null, kmers2=null;
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				for(int rnum=0; rnum<reads.size(); rnum++){
					Read r1=reads.get(rnum);
					Read r2=r1.mate;
					assert(r1!=r2);
					
					if(eccByOverlap && r1!=null && r2!=null){BBMerge.findOverlapStrict(r1, r2, true);}
					
					if(!TRIM_AFTER_MARKING && (TRIM_LEFT_THIS_PASS || TRIM_RIGHT_THIS_PASS)){
						if(r1!=null){basesTrimmed+=TrimRead.trimFast(r1, TRIM_LEFT_THIS_PASS, TRIM_RIGHT_THIS_PASS, TRIM_QUALITY, trimE, 1);}
						if(r2!=null){basesTrimmed+=TrimRead.trimFast(r2, TRIM_LEFT_THIS_PASS, TRIM_RIGHT_THIS_PASS, TRIM_QUALITY, trimE, 1);}
					}

					int depthAL1=-1, depthAL2=-1;
					int truedepth1=-1, truedepth2=-1;
					int mintruedepth1=-1, mintruedepth2=-1;
					
					int readcount=0;
					int basecount=0;
					
					int lowcount=0, totalcount=0, ec1=0, ec2=0;
					boolean error1=false, error2=false, uncorrectable1=false, uncorrectable2=false, marked1=false, marked2=false;
					
					if(r1!=null && r1.bases!=null){
						readcount++;
						basecount+=r1.length();
						if(r1.length()>=k){
							if(verbose){outstream.println();}
							kmers1=r1.toKmers(k, 0, kmers1, true, longkmer);
							cov1=generateCoverage(kca, k, cov1, kmers1, k<32);
							int[] cov=cov1.clone();
							sortCoverageAndIncrementHistogram(cov);
							
							if(cov!=null){
								final int covlast=cov.length-1;
								final int high=cov[(int)((covlast)*(1-HIGH_PERCENTILE))];
								final int low=cov[(int)((covlast)*(1-LOW_PERCENTILE))];
								mintruedepth1=low;
								int aboveLimit=covlast;
								int lc=0;
								final int mindepth=Tools.max(MIN_DEPTH, high/ERROR_DETECT_RATIO);
								truedepth1=cov[(int)((covlast)*(1-DEPTH_PERCENTILE))];
								while(aboveLimit>=0 && cov[aboveLimit]<mindepth){aboveLimit--;}
								if(aboveLimit+1>=MIN_KMERS_OVER_MIN_DEPTH || (aboveLimit>=0 && MIN_KMERS_OVER_MIN_DEPTH>cov.length)){
									depthAL1=cov[(int)(aboveLimit*(1-DEPTH_PERCENTILE))];
								}
								if(high<=LTHRESH || (high>=HTHRESH && low<=LTHRESH) || high>=low*ERROR_DETECT_RATIO){
									error1=true;
									if(high<=LTHRESH){errorType1++;}
									if(high>=HTHRESH && low<=LTHRESH){errorType2++;}
									if(high>=low*ERROR_DETECT_RATIO){errorType3++;}
								}
								
								totalcount+=cov.length;
								if(cov[0]<=LTHRESH){
									lc+=cov.length;
								}else if(high>=HTHRESH){
									int lim=Tools.min(LTHRESH, high/ERROR_DETECT_RATIO);
									for(int i=covlast; i>=0 && cov[i]<=lim; i--){lc++;}
								}
								lowcount+=lc;
								
								if(rhistogram!=null){
									int d=depthAL1>=0 ? depthAL1 : truedepth1>=0 ? truedepth1 : 0;
									d=Tools.min(d, HIST_LEN-1);
									rhistogram.incrementAndGet(d);
									bhistogram.addAndGet(d, r1.length());
								}
							}
						}
					}
					if(r2!=null && r2.bases!=null){
						readcount++;
						basecount+=r2.length();
						if(r2.length()>=k){
							if(verbose){outstream.println();}
							kmers2=r2.toKmers(k, 0, kmers2, true, longkmer);
							cov2=generateCoverage(kca, k, cov2, kmers2, k<32);
							int[] cov=cov2.clone();
							sortCoverageAndIncrementHistogram(cov);
							
							if(cov!=null){
								final int covlast=cov.length-1;
								final int high=cov[(int)((covlast)*(1-HIGH_PERCENTILE))];
								final int low=cov[(int)((covlast)*(1-LOW_PERCENTILE))];
								mintruedepth2=low;
								int aboveLimit=covlast;
								int lc=0;
								final int mindepth=Tools.max(MIN_DEPTH, high/ERROR_DETECT_RATIO);
								truedepth2=cov[(int)((covlast)*(1-DEPTH_PERCENTILE))];
								while(aboveLimit>=0 && cov[aboveLimit]<mindepth){aboveLimit--;}
								if(aboveLimit+1>=MIN_KMERS_OVER_MIN_DEPTH || (aboveLimit>=0 && MIN_KMERS_OVER_MIN_DEPTH>cov.length)){
									depthAL2=cov[(int)(aboveLimit*(1-DEPTH_PERCENTILE))];
								}
								if(high<=LTHRESH || (high>=HTHRESH && low<=LTHRESH) || high>=low*ERROR_DETECT_RATIO){
									error2=true;
									if(high<=LTHRESH){errorType1++;}
									if(high>=HTHRESH && low<=LTHRESH){errorType2++;}
									if(high>=low*ERROR_DETECT_RATIO){errorType3++;}
								}
								
								totalcount+=cov.length;
								if(cov[0]<=LTHRESH){
									lc+=cov.length;
								}else if(high>=HTHRESH){
									int lim=Tools.min(LTHRESH, high/ERROR_DETECT_RATIO);
									for(int i=covlast; i>=0 && cov[i]<=lim; i--){lc++;}
								}
								lowcount+=lc;
								cov2=cov;
								
								if(rhistogram!=null){
									int d=depthAL2>=0 ? depthAL2 : truedepth2>=0 ? truedepth2 : 0;
									d=Tools.min(d, HIST_LEN-1);
									rhistogram.incrementAndGet(d);
									bhistogram.addAndGet(d, r2.length());
								}
							}
						}
					}
					
					if(RENAME_THIS_PASS){
						if(r2==null){
							final String s="id="+r1.numericID+",d1="+depthAL1+(CORRECT_ERRORS_THIS_PASS ? ",e1="+(ec1<0 ? -ec1 : 0) : "");
							r1.id=s;
							if(EA){
								int[] quad=parseDepth(r1.id, null);
								assert(quad[0]==depthAL1);
								assert(quad[1]==-1);
								assert(quad[2]==(CORRECT_ERRORS_THIS_PASS ? (ec1<0 ? -ec1 : 0) : -1));
								assert(quad[3]==-1);
							}
						}else{
							final String s="id="+r1.numericID+",d1="+depthAL1+",d2="+depthAL2+(CORRECT_ERRORS_THIS_PASS ? ",e1="+(ec1<0 ? -ec1 : 0)+",e2="+(ec2<0 ? -ec2 : 0) : "");
							r1.id=s+" /1";
							r2.id=s+" /2";
							if(EA){
								int[] quad=parseDepth(r1.id, null);
								assert(quad[0]==depthAL1);
								assert(quad[1]==depthAL2);
								assert(quad[2]==(CORRECT_ERRORS_THIS_PASS ? (ec1<0 ? -ec1 : 0) : -1));
								assert(quad[3]==(CORRECT_ERRORS_THIS_PASS ? (ec2<0 ? -ec2 : 0) : -1));
							}
						}
					}
					
					r1.errors=lowcount;
					
					int maxDepth=MAX_DEPTH;
					int targetDepth=TARGET_DEPTH;
					
					if(lowcount>0){
//						targetDepth=(int)((TARGET_DEPTH_BAD_LOW*(long)lowcount+TARGET_DEPTH_BAD_HIGH*(totalcount-(long)lowcount))/totalcount);
						
						double fractionGood=(totalcount-lowcount)/(float)totalcount;
						targetDepth=(int)(TARGET_DEPTH_BAD_LOW+(TARGET_DEPTH_BAD_HIGH-TARGET_DEPTH_BAD_LOW)*(fractionGood*fractionGood));
						assert(TARGET_DEPTH_BAD_LOW<=TARGET_DEPTH_BAD_HIGH);
						assert(TARGET_DEPTH>=99999999 || (targetDepth>0 && targetDepth<=TARGET_DEPTH)) :
							targetDepth+", "+TARGET_DEPTH+", "+TARGET_DEPTH_BAD_LOW+", "+TARGET_DEPTH_BAD_HIGH+", "+lowcount+", "+totalcount;
						assert(TARGET_DEPTH>=99999999 || (targetDepth>=TARGET_DEPTH_BAD_LOW && targetDepth<=TARGET_DEPTH_BAD_HIGH)) :
							targetDepth+", "+TARGET_DEPTH+", "+TARGET_DEPTH_BAD_LOW+", "+TARGET_DEPTH_BAD_HIGH+", "+lowcount+", "+totalcount;
						maxDepth=targetDepth;
					}
					
					final int minAL=(depthAL1>=0 ? (depthAL2>=0 ? Tools.min(depthAL1, depthAL2) : depthAL1) : depthAL2);
					final int maxAL=Tools.max(depthAL1, depthAL2);
					final int minTrueDepth=(r2==null ? truedepth1 : Tools.min(truedepth1, truedepth2));
					final int maxTrueDepth=Tools.max(truedepth1, truedepth2);
					final int depthproxyAL=USE_LOWER_DEPTH ? minAL : maxAL;
					final int truedepthproxy=USE_LOWER_DEPTH ? minTrueDepth : maxTrueDepth;
					long coin=0;
					if(depthproxyAL>maxDepth && (error1 || error2 || !DISCARD_BAD_ONLY)){
						if(r1.rand<0){
							coin=randy.nextInt(depthproxyAL)+1;
						}else{
							coin=((long)(r1.rand*depthproxyAL))+1;
						}
					}
					
					totalReads+=readcount;
					totalBases+=basecount;
					
					boolean toss=(depthproxyAL<0 || coin>targetDepth || (r1!=null && r1.length()<MIN_LENGTH) || (r2!=null && r2.length()<MIN_LENGTH));
					if(TOSS_ERROR_READS && (error1 || error2)){
						if(SAVE_RARE_READS && depthproxyAL<=targetDepth && depthproxyAL>=HTHRESH){
							//do nothing
						}else if(!REQUIRE_BOTH_BAD || r2==null || (error1 && error2)){
							toss=true;
						}
					}
					
					if(TOSS_BY_LOW_TRUEDEPTH && !SAVE_RARE_READS && maxTrueDepth<MIN_DEPTH && (!REQUIRE_BOTH_BAD || (mintruedepth1<MIN_DEPTH && mintruedepth2<MIN_DEPTH))){
						toss=true;
					}
					
					if(KEEP_ALL){toss=false;}
					
//					if((r==null || verybad1) && (r2==null || verybad2)){toss=true;} //Always toss verybad reads.  Turned out to not be helpful.
					
					if(error1){errorReads++;}
					if(error2){errorReads++;}
					if(error1 || error2){errorPairs++;}
					
					if(toss){
						if(tossList!=null){tossList.add(r1);}
						readsTossed+=readcount;
						basesTossed+=basecount;
					}else{
						if(CORRECT_ERRORS_THIS_PASS){
							if(r1!=null && r1.length()>=k){
								if(MARK_ERRORS_ONLY){
									ec1=markErrors(r1, cov1, kmers1, kca, k, EC_LTHRESH, EC_HTHRESH, ERROR_CORRECT_RATIO, MAX_ERRORS_TO_CORRECT, true, true, qhist, longkmer);
									errorsMarked+=ec1;
									if(ec1>0){marked1=true;}
								}else{
									ec1=correctErrors(r1, cov1, kmers1, kca, k, EC_LTHRESH, EC_HTHRESH, ERROR_CORRECT_RATIO, MAX_ERRORS_TO_CORRECT, MAX_QUAL_TO_CORRECT, true, true, qhist, MARK_ERRORS_ONLY, longkmer);
									if(ec1>=0){
										errorsCorrected+=ec1;
									}else{
										uncorrectable1=true;
										if(MAX_ERRORS_TO_CORRECT>0){regenerateKmersAndCoverage(r1, kmers1, cov1, kca, k, false, longkmer);}
										if(MARK_UNCORRECTABLE_ERRORS){
											ec1=markErrors(r1, cov1, kmers1, kca, k, EC_LTHRESH, EC_HTHRESH, ERROR_CORRECT_RATIO, MAX_ERRORS_TO_CORRECT, true, true, qhist, longkmer);
											errorsMarked+=ec1;
											if(ec1>0){marked1=true;}
										}else{
											errorsDetected-=ec1;
										}
									}
								}
								if(TRIM_AFTER_MARKING){
									if(marked1 || TRIM_EVEN_IF_NO_ERRORS_DETECTED){
										basesTrimmed+=TrimRead.trimFast(r1, TRIM_LEFT_THIS_PASS, TRIM_RIGHT_THIS_PASS, TRIM_QUALITY, trimE, 1);
									}
								}
							}

							if(r2!=null && r2.length()>=k){
								if(MARK_ERRORS_ONLY){
									ec2=markErrors(r2, cov2, kmers2, kca, k, EC_LTHRESH, EC_HTHRESH, ERROR_CORRECT_RATIO, MAX_ERRORS_TO_CORRECT, true, true, qhist, longkmer);
									errorsMarked+=ec2;
									if(ec2>0){marked2=true;}
								}else{
									ec2=correctErrors(r2, cov2, kmers2, kca, k, EC_LTHRESH, EC_HTHRESH, ERROR_CORRECT_RATIO, MAX_ERRORS_TO_CORRECT, MAX_QUAL_TO_CORRECT, true, true, qhist, MARK_ERRORS_ONLY, longkmer);
									if(ec2>=0){
										errorsCorrected+=ec2;
									}else{
										uncorrectable2=true;
										if(MAX_ERRORS_TO_CORRECT>0){regenerateKmersAndCoverage(r2, kmers2, cov2, kca, k, false, longkmer);}
										if(MARK_UNCORRECTABLE_ERRORS){
											ec2=markErrors(r2, cov2, kmers2, kca, k, EC_LTHRESH, EC_HTHRESH, ERROR_CORRECT_RATIO, MAX_ERRORS_TO_CORRECT, true, true, qhist, longkmer);
											errorsMarked+=ec2;
											if(ec2>0){marked2=true;}
										}else{
											errorsDetected-=ec2;
										}
									}
								}
								if(TRIM_AFTER_MARKING){
									if(marked2 || TRIM_EVEN_IF_NO_ERRORS_DETECTED){
										basesTrimmed+=TrimRead.trimFast(r2, TRIM_LEFT_THIS_PASS, TRIM_RIGHT_THIS_PASS, TRIM_QUALITY, trimE, 1);
									}
								}
							}
						}

						if(keepList!=null){keepList.add(r1);}
						readsKept+=readcount;
						basesKept+=basecount;
					}
					
					if(depthAL1<LOW_BIN_DEPTH && depthAL2<LOW_BIN_DEPTH){
						readsLowBin+=readcount;
						basesLowBin+=basecount;
						if(lowList!=null){lowList.add(r1);}
//					}else if((depth1<0 || depth1>HIGH_BIN_DEPTH) && (depth2<0 || depth2>=HIGH_BIN_DEPTH)){
					}else if((depthAL1<LOW_BIN_DEPTH || depthAL1>HIGH_BIN_DEPTH) && (depthAL2<LOW_BIN_DEPTH || depthAL2>=HIGH_BIN_DEPTH)){
						readsHighBin+=readcount;
						basesHighBin+=basecount;
						if(highList!=null){highList.add(r1);}
					}else{
						assert((depthAL1>=LOW_BIN_DEPTH && depthAL1<=HIGH_BIN_DEPTH) || (depthAL2>=LOW_BIN_DEPTH && depthAL2<=HIGH_BIN_DEPTH)) :
							depthAL1+", "+depthAL2+", "+LOW_BIN_DEPTH+", "+HIGH_BIN_DEPTH;
						readsMidBin+=readcount;
						basesMidBin+=basecount;
						if(midList!=null){midList.add(r1);}
					}
					
					if(uncorrectable1 || uncorrectable2){
						readsUncorrected+=readcount;
						basesUncorrected+=basecount;
						if(uncList!=null){uncList.add(r1);}
					}
				}
				
				
				if(storage!=null){
					synchronized(storage){
						storage.addAll(keepList);
						if(ADD_BAD_READS_COUNTUP){storage.addAll(tossList);}
					}
				}
				

				if(rosk!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
//					System.err.println("Adding list "+ln.id+" of length "+reads.size());
					rosk.add(keepList, ln.id);
					keepList.clear();
				}
				if(rost!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
//					System.err.println("Adding list "+ln.id+" of length "+reads.size());
					rost.add(tossList, ln.id);
					tossList.clear();
				}
				
				if(rosl!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
//					System.err.println("Adding list "+ln.id+" of length "+reads.size());
					rosl.add(lowList, ln.id);
					lowList.clear();
				}
				if(rosm!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
//					System.err.println("Adding list "+ln.id+" of length "+reads.size());
					rosm.add(midList, ln.id);
					midList.clear();
				}
				if(rosh!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
//					System.err.println("Adding list "+ln.id+" of length "+reads.size());
					rosh.add(highList, ln.id);
					highList.clear();
				}
				if(rosu!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
//					System.err.println("Adding list "+ln.id+" of length "+reads.size());
					rosu.add(uncList, ln.id);
					uncList.clear();
				}
				
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln);
			if(verbose){System.err.println("Returned list");}
		}
		
		
		void normalizeInThreadByCountup() {
			
			Kmer longkmer=(k<32 ? null : new Kmer(k));
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			final ArrayList<Read> keepList=(rosk==null ? null : new ArrayList<Read>(Shared.bufferLen()));
			final ArrayList<Read> tossList=(rost==null ? null : new ArrayList<Read>(Shared.bufferLen()));
			
			int[] cov=null, covSorted=null, covup=null;
			long[] kmers1=null, kmers2=null, kmers3=null;
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				for(int rnum=0; rnum<reads.size(); rnum++){
					Read r=reads.get(rnum);
					Read r2=r.mate;
					assert(r!=r2);
					
					int readcount=0;
					int basecount=0;
					int errors=0, nonerrors=0;
					
					boolean k1valid=false, k2valid=false;
					
					if(r!=null && r.bases!=null){
						readcount++;
						basecount+=r.length();
						if(r.length()>=k){
							if(verbose){outstream.println();}
							kmers1=r.toKmers(k, 0, kmers1, true, longkmer);
							k1valid=true;
						}
					}
					if(r2!=null && r2.bases!=null){
						readcount++;
						basecount+=r2.length();
						if(r2.length()>=k){
							if(verbose){outstream.println();}
							kmers2=r2.toKmers(k, 0, kmers2, true, longkmer);
							k2valid=true;
						}
					}
					
					final int mergelen=(k1valid ? kmers1.length : 0)+(k2valid ? kmers2.length : 0);
					int valid=0, unique=0, desired=0, needed=0, badlyneeded=0;
					if(mergelen>0){
						if(kmers3==null || kmers3.length!=mergelen){kmers3=new long[mergelen];}
						int j=0;
						if(k1valid){
							for(int i=0; i<kmers1.length; i++, j++){kmers3[j]=kmers1[i];}
						}
						if(k2valid){
							for(int i=0; i<kmers2.length; i++, j++){kmers3[j]=kmers2[i];}
						}
						Arrays.sort(kmers3);

						if(cov==null || cov.length!=mergelen){cov=new int[mergelen];}
						if(covup==null || covup.length!=mergelen){covup=new int[mergelen];}
						for(int i=0; i<mergelen; i++){
							long kmer=kmers3[i];
							if(kmer==-1){
								cov[i]=-1;
								covup[i]=-1;
							}else if(IGNORE_DUPLICATE_KMERS_COUNTUP && i>0 && kmer==kmers3[i-1]){
								cov[i]=-1;
								covup[i]=-1;
								valid++;
							}else{
								cov[i]=kca.read(kmer);
								covup[i]=kcaup.read(kmer);
								valid++;
								unique++;
								if(cov[i]>=MIN_DEPTH){
									desired++;
									if(covup[i]<TARGET_DEPTH){
										needed++;
										if(covup[i]<(Tools.min(TARGET_DEPTH, cov[i])*3)/4){badlyneeded++;}
									}
								}
							}
						}
						final int invalid=cov.length-valid;
						
						if(covSorted==null || covSorted.length!=mergelen){covSorted=new int[mergelen];}
						for(int i=0; i<cov.length; i++){covSorted[i]=cov[i];}
						Arrays.sort(covSorted);
						
						int prev=-1;
						for(int i=0; i<covSorted.length; i++){
							int x=covSorted[i];
							if(prev>-1){
								if((x>=HTHRESH && prev<=LTHRESH) || x>=prev*ERROR_DETECT_RATIO){
									errors=covSorted.length-i;
									break;
								}else{nonerrors++;}
							}
							prev=x;
						}
					}

					int t1=Tools.max(8, (unique+5)/6);
					int t2=Tools.max(2, (unique+23)/24);
					
					boolean toss=!((needed>=t1 || badlyneeded>=t2) && (desired>=MIN_KMERS_OVER_MIN_DEPTH || unique<MIN_KMERS_OVER_MIN_DEPTH));
					if(TOSS_ERROR_READS && errors>8 && (needed<2*t1 && badlyneeded<2*t2)){toss=true;}
					if(TOSS_ERROR_READS && errors>unique/2 && (needed<3*t1 && badlyneeded<4*t2)){toss=true;}
//					assert(false) : "\n"+TOSS_ERROR_READS+", "+unique+", "+desired+", "+needed+", "+badlyneeded+", "+errors+", "+t1+", "+t2;
//					System.out.println("valid="+valid+", unique="+unique+", desired="+desired+", needed="+needed+", toss="+toss);
					if(KEEP_ALL){toss=false;}
					
					totalReads+=readcount;
					totalBases+=basecount;
					
					if(toss){
//						System.out.println("valid="+valid+", unique="+unique+", desired="+desired+", needed="+needed+", toss="+toss);
						if(tossList!=null){tossList.add(r);}
						readsTossed+=readcount;
						basesTossed+=basecount;
					}else{
//						System.out.println("valid="+valid+", unique="+unique+", desired="+desired+", needed="+needed+", toss="+toss);
//						System.out.println("valid="+valid+", unique="+unique+", desired="+desired+", needed="+needed+", toss="+toss+"\n"+Arrays.toString(cov)
//								+"\n"+Arrays.toString(covup)+"\n"+Arrays.toString(kmers3));
						for(int i=0; i<mergelen; i++){
							if(cov[i]>=MIN_DEPTH){
								long kmer=kmers3[i];
								kcaup.increment(kmer);
							}
						}
						if(keepList!=null){keepList.add(r);}
						readsKept+=readcount;
						basesKept+=basecount;
					}
					
					if(mergelen>0){
//						Arrays.sort(cov);
//						incrementHistogramSorted(cov);
						incrementHistogramSorted(covSorted);
					}
				}
				
				if(storage!=null){
					synchronized(storage){
						storage.addAll(keepList);
					}
				}
				
				if(rosk!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
//					System.err.println("Adding list "+ln.id+" of length "+reads.size());
					rosk.add(keepList, ln.id);
					keepList.clear();
				}
				if(rost!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
//					System.err.println("Adding list "+ln.id+" of length "+reads.size());
					rost.add(tossList, ln.id);
					tossList.clear();
				}

				assert(rosl==null) : "Low fraction out not supported by countup.";
				assert(rosm==null) : "Mid fraction out not supported by countup.";
				assert(rosh==null) : "High fraction out not supported by countup.";
				assert(rosu==null) : "TODO - Uncorrectable fraction out not supported by countup.";
				
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln);
			if(verbose){System.err.println("Returned list");}
		}
		
		private final int[] getSortedCoverageAndIncrementHistogram(Read r, int[] cov, long[] kmers,
				boolean kmersAlreadyValid, boolean kmersAlreadyCanonical, boolean coverageAlreadyValid, Kmer longkmer){
			assert(r!=null && r.bases!=null && r.length()>=k) : r;
			
			if(!coverageAlreadyValid){
				if(!kmersAlreadyValid){kmers=r.toKmers(k, 0, kmers, false, longkmer);}
				cov=generateCoverage(kca, k, cov, kmers, (!kmersAlreadyCanonical && k<32));
			}
			
			sortCoverageAndIncrementHistogram(cov);
			return cov;
		}
		
		private void sortCoverageAndIncrementHistogram(int[] cov){
			if(cov==null || cov.length==0){return;}
			Arrays.sort(cov);
			Tools.reverseInPlace(cov);
			incrementHistogramSorted(cov);
		}
		
		/** Handles coverage sorted in either direction */
		private final void incrementHistogramSorted(int[] cov){
			if(hist==null || cov==null || cov.length==0){return;}
			
//			outstream.println(Arrays.toString(cov));
			
			int last=cov[0];
			long sum=0;
//			long sum2=0;
			int i=0;
			while(i<cov.length && cov[i]<0){i++;}
			for(; i<cov.length; i++){
				int x=cov[i];
//				outstream.println("Processing "+x);
				if(x<0){break;}
				int y=Tools.min(x, HIST_LEN-1);
				if(y==last){sum++;}
				else if(sum>0){
//					outstream.println("Incrementing "+last+" by "+sum);
//					sum2+=sum;
					if(last<hist.length){hist[last]+=sum;}
					else{khistogram.addAndGet(last, sum);}
					sum=1;
				}
				last=y;
			}
//			outstream.println("Ended loop");
			if(sum>0){
//				outstream.println("Incrementing "+last+" by "+sum);
//				sum2+=sum;
				if(last<hist.length){hist[last]+=sum;}
				else{khistogram.addAndGet(last, sum);}
			}
//			assert(sum2==cov.length) : sum2+", "+cov.length+", "+last+", "+sum;
		}
		
		private final ConcurrentReadInputStream cris;
		/** Premade table holding counts of input kmers */
		private final KCountArray kca;
		/** Dynamic table holding counts of output kmers */
		private final KCountArray kcaup;
		/** kmer length */
		private final int k;
		/** Stream for kept reads */
		private final ConcurrentReadOutputStream rosk;
		/** Stream for tossed reads */
		private final ConcurrentReadOutputStream rost;
		/** Stream for low-count reads */
		private final ConcurrentReadOutputStream rosl;
		/** Stream for mid-count reads */
		private final ConcurrentReadOutputStream rosm;
		/** Stream for high-count reads */
		private final ConcurrentReadOutputStream rosh;
		/** Stream for reads with uncorrectable errors */
		private final ConcurrentReadOutputStream rosu;

		public final long[] hist=new long[THREAD_HIST_LEN];//(USE_HISTOGRAM ? new long[HIST_LEN] : null);
		public final long[] qhist=new long[128];
		
		private final ArrayList<Read> storage;
		
		private long totalBases=0;
		private long totalReads=0;
//		private final java.util.Random randy=new java.util.Random();
		private Random randy; //Note that Random does not support nextLong(long)

		public long readsKept=0;
		public long readsTossed=0;
		public long readsLowBin=0;
		public long readsMidBin=0;
		public long readsHighBin=0;
		public long readsUncorrected=0;
		public long basesKept=0;
		public long basesTossed=0;
		public long basesLowBin=0;
		public long basesMidBin=0;
		public long basesHighBin=0;
		public long basesUncorrected=0;

		public long errorReads=0;
		public long errorPairs=0;
		public long errorType1=0;
		public long errorType2=0;
		public long errorType3=0;

		public long errorsDetected=0;
		public long errorsCorrected=0;
		public long errorsMarked=0;
		public long basesTrimmed=0;
		
		boolean errorStateT=false;
	}

	public static boolean errorState(){return errorState;}
	public static boolean setErrorState(boolean b){return errorState=b;}
	
	public static PrintStream outstream=System.err;

	private static long minHeight=2;
	private static long minVolume=5;
	private static int minWidth=3;
	private static int minPeak=2;
	private static int maxPeak=Integer.MAX_VALUE;
	private static int maxPeakCount=10;
	private static int ploidy=-1;
	private static boolean doLogScale=false;
	private static double logWidth=0.05;

	public static int THREAD_HIST_LEN=1<<12;
	public static int HIST_LEN=1<<20;
	public static long HIST_LEN_PRINT=HIST_LEN;
	public static long HIST_COLUMNS=3;
	public static boolean USE_KHISTOGRAM=false;
	public static boolean USE_RHISTOGRAM=false;
	public static boolean PRINT_ZERO_COVERAGE=false;
	public static AtomicLongArray khistogram;
	public static AtomicLongArray rhistogram;
	public static AtomicLongArray bhistogram;
	public static long[] qhist_total;
	
	private static int THREADS=Shared.threads();
	private static boolean verbose=false;
	private static boolean errorState=false;
	
	private static boolean EA=Shared.EA();
	
	private static boolean eccByOverlap=false;
	private static boolean eccByOverlapAuto=false;
	
	/** High-depth reads will be downsampled to this level in the current pass */
	private static int TARGET_DEPTH=100;
	/** Error-containing reads will be downsampled to at least this level in the current pass */
	private static int TARGET_DEPTH_BAD_LOW=100;
	/** Error-containing reads will be downsampled to at most this level in the current pass */
	private static int TARGET_DEPTH_BAD_HIGH=100;
	/** High-depth reads will be downsampled to this level in the final pass */
	private static int TARGET_DEPTH_F=100;
	/** High-depth reads will be downsampled to this level in the first pass */
	private static int TARGET_DEPTH_1=-1;
	/** Reads under this depth will not be downsampled */
	private static int MAX_DEPTH=-1;
	/** Reads under this depth will be discarded, and kmers under this depth will be ignored */
	private static int MIN_DEPTH=5;
	/** Reads without this many kmers of at least min depth will be discarded */
	private static int MIN_KMERS_OVER_MIN_DEPTH=15;
	/** Position in sorted kmer depths array to use as proxy for overall read depth */
	private static float DEPTH_PERCENTILE=0.54f;
	
	/** Normalize based on depth of read with lower depth, instead of read with higher depth */
	public static boolean USE_LOWER_DEPTH=true;
	/** Throw out reads with depth at absolute depth percentile below mindepth */
	public static boolean TOSS_BY_LOW_TRUEDEPTH=true;
	/** Throw out reads containing errors in the current pass */
	public static boolean TOSS_ERROR_READS=false;
	/** Throw out reads containing errors in the final pass */
	public static boolean TOSS_ERROR_READS_F=false;
	/** Throw out reads containing errors in the first pass */
	public static boolean TOSS_ERROR_READS_1=false;
	/** Only downsample error reads on current pass (keep all error-free reads) */
	public static boolean DISCARD_BAD_ONLY=false;
	/** Only downsample error reads on first pass (keep all error-free reads) */
	public static boolean DISCARD_BAD_ONLY_F=false;
	/** Only downsample error reads on final pass (keep all error-free reads) */
	public static boolean DISCARD_BAD_ONLY_1=false;
	/** Require both reads in a pair to be bad before tossing the read */
	public static boolean REQUIRE_BOTH_BAD=false;
	/** Don't toss error reads with depth below max */
	public static boolean SAVE_RARE_READS=false;
	/** Position in sorted kmer depths array to use as proxy for high depth kmer */
	public static float HIGH_PERCENTILE=0.90f;
	/** Position in sorted kmer depths array to use as proxy for low depth kmer */
	public static float LOW_PERCENTILE=0.25f;
	/** Position in sorted kmer depths array to use as proxy for low depth kmer, during countup presort pass */
	public static float LOW_PERCENTILE_COUNTUP=0.20f;
	/** Set to true to keep error reads during countup presort pass */
	public static boolean ADD_BAD_READS_COUNTUP=false;
	
	/** Reads with a high/low ratio of at least this are considered error reads. */
	public static int ERROR_DETECT_RATIO=125;
	/** Threshold for high kmer in detection.  A high kmer at this or above is considered possibly non-error. */
	public static int HTHRESH=12;
	/** Threshold for low kmer in detection.  Kmers at this and below are always considered errors. */
	public static int LTHRESH=3;
	
	/** Reads with a high/low ratio of at least this are considered error reads. */
	public static int ERROR_CORRECT_RATIO=140;
	/** Threshold for high kmer in correction.  A high kmer at this or above considered possibly non-error. */
	public static int EC_HTHRESH=22;
	/** Threshold for low kmer in correction.  Kmers at this and below are considered errors if an adjacent kmer is at or above the high thresh. */
	public static int EC_LTHRESH=2;

	public static double TARGET_BAD_PERCENT_LOW=0.85;
	public static double TARGET_BAD_PERCENT_HIGH=1.5;

	private static long FILTERBYTES=-1;

	private static int SUFFIX_LEN=3;
	private static int PREFIX_LEN=3;
	
	private static boolean TRIM_LEFT_THIS_PASS=false;
	private static boolean TRIM_RIGHT_THIS_PASS=false;
	private static boolean RENAME_THIS_PASS=false;

	private static boolean CORRECT_ERRORS_THIS_PASS=false;
	private static boolean MARK_ERRORS_ONLY=false;
	private static boolean TRIM_AFTER_MARKING=false;
	private static boolean TRIM_EVEN_IF_NO_ERRORS_DETECTED=true;
	private static boolean MARK_WITH_1=false;
	private static boolean MARK_UNCORRECTABLE_ERRORS=false;
	private static boolean USE_ECC1=false;
	private static boolean USE_ECCF=false;
	private static boolean CORRECT_FROM_LEFT=true;
	private static boolean CORRECT_FROM_RIGHT=true;
	
	private static double prefilterFraction=0.35;
	
	private static int LOW_BIN_DEPTH=10;
	private static int HIGH_BIN_DEPTH=80;
	
	/** ECC_LIMIT */
	private static int MAX_ERRORS_TO_CORRECT=3;
	private static int MAX_QUAL_TO_CORRECT=127;
	

	public static boolean IGNORE_DUPLICATE_KMERS_COUNTUP=true;
	
	public static boolean CANONICAL=true;
	public static boolean ZERO_BIN=false;
	public static boolean FIX_SPIKES=false;
	public static boolean KEEP_ALL=false;
	public static boolean ordered=false;
	public static boolean overwrite=true;
	public static boolean append=false;
	public static boolean prefilter=false;
	public static boolean renameReads=false;
	public static boolean DETERMINISTIC=true;
	public static boolean COUNTUP=false;
	public static boolean ANALYZE_TOPOLOGY=false;
	/** Quality-trim left side of reads before further processing. */
	public static boolean TRIM_LEFT=false;
	/** Quality-trim right side of reads before further processing. */
	public static boolean TRIM_RIGHT=false;
	public static int MIN_LENGTH=1;
	/** Trim until 2 consecutive bases are encountered with at least this quality. */
	public static float TRIM_QUALITY=5;
	/** Error rate for trimming (derived from trimq) */
	private static float trimE;
	
	public static boolean REMOVE_TEMP_FILES=true;
	public static boolean USE_TMPDIR=true;
	public static String TMPDIR=Shared.tmpdir();
	public static boolean useTmpdir(){return USE_TMPDIR && TMPDIR!=null;}
	
	private static HashSet<String> temp_file_set=null;

	public static AtomicLong peaks=new AtomicLong();
	public static AtomicLong spikes=new AtomicLong();
	public static AtomicLong flats=new AtomicLong();
	public static AtomicLong valleys=new AtomicLong();
	public static AtomicLong slopes=new AtomicLong();
	
	public static final byte FIXED_N_QUAL=20;
	
}
