package parse;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import align2.AbstractMapThread;
import align2.IndexMaker4;
import align2.IndexMaker5;
import align2.QualityTools;
import aligner.SideChannel3;
import aligner.SideChannel4;
import bloom.KCountArray;
import bloom.KmerCountAbstract;
import cardinality.CardinalityTracker;
import cardinality.LogLog16;
import cardinality.LogLog2;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile4;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import jgi.BBMerge;
import jgi.CalcTrueQuality;
import kmer.AbstractKmerTable;
import shared.KillSwitch;
import shared.Shared;
import shared.SyncHeart;
import shared.Tools;
import shared.TrimRead;
import simd.Vector;
import sketch.SketchObject;
import stream.BamStreamer;
import stream.BamWriter;
import stream.ConcurrentDepot;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.FastaStreamer;
import stream.FastqStreamer;
import stream.FastqWriter;
import stream.Read;
import stream.ReadStreamByteWriter;
import stream.ReadStreamWriter;
import stream.SamLine;
import stream.SamStreamer;
import stream.SamReadInputStream;
import stream.SamWriter;
import stream.StreamerFactory;
import stream.bam.BgzfOutputStreamMT;
import stream.bam.BgzfSettings;
import structures.IntList;
import tax.TaxTree;
import tracker.EntropyTracker;
import tracker.ReadStats;
import var2.CallVariants;

/**
 * @author Brian Bushnell
 * @date Mar 21, 2014
 *
 */
public class Parser {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a new Parser instance with default settings */
	public Parser(){}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main parsing entry point that delegates to specialized parsing methods.
	 * Attempts to parse the given argument by trying different parameter categories
	 * in sequence until a match is found.
	 *
	 * @param arg The complete argument string (e.g., "qtrim=r")
	 * @param a The parameter name portion (e.g., "qtrim")
	 * @param b The parameter value portion (e.g., "r")
	 * @return true if the argument was successfully parsed, false otherwise
	 */
	public boolean parse(String arg, String a, String b){
		return (parseStatic(arg, a, b) || parseNonStatic(arg, a, b));
	}
	
	/**
	 * Entry point for parsing static fields.
	 *
	 * @param arg The complete argument string (e.g., "ziplevel=4")
	 * @param a The parameter name portion (e.g., "ziplevel")
	 * @param b The parameter value portion (e.g., "4")
	 * @return true if the argument was successfully parsed, false otherwise
	 */
	public static boolean parseStatic(String arg, String a, String b) {
		if(isJavaFlag(arg)){return true;}

		if(parseQuality(arg, a, b)){return true;}
		if(parseZip(arg, a, b)){return true;}
		if(parseSam(arg, a, b)){return true;}
		if(parseFasta(arg, a, b)){return true;}
		if(parseCommonStatic(arg, a, b)){return true;}
		if(parseHist(arg, a, b)){return true;}
		if(parseQualityAdjust(arg, a, b)){return true;}
		return false;
	}
	
	/**
	 * Entry point for parsing non-static fields.
	 *
	 * @param arg The complete argument string (e.g., "qtrim=r")
	 * @param a The parameter name portion (e.g., "qtrim")
	 * @param b The parameter value portion (e.g., "r")
	 * @return true if the argument was successfully parsed, false otherwise
	 */
	public boolean parseNonStatic(String arg, String a, String b){
		if(parseFiles(arg, a, b)){return true;}
		if(parseCommon(arg, a, b)){return true;}
		if(parseTrim(arg, a, b)){return true;}
		if(parseInterleaved(arg, a, b)){return true;}
		if(parseMapping(arg, a, b)){return true;}
		if(parseCardinality(arg, a, b)){return true;}
		return false;
	}

	/**
	 * Parses common parameters used across multiple BBTools programs.
	 * Handles parameters like reads limit, sample rate, overwrite permissions,
	 * and basic processing options.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public boolean parseCommon(String arg, String a, String b){
		if(a.equals("reads") || a.equals("maxreads")){
			maxReads=Parse.parseKMG(b);
		}else if(a.equals("samplerate")){
			samplerate=Float.parseFloat(b);
			assert(samplerate<=1f && samplerate>=0f) : "samplerate="+samplerate+"; should be between 0 and 1";
		}else if(a.equals("sampleseed")){
			sampleseed=Long.parseLong(b);
		}else if(a.equals("append") || a.equals("app")){
			append=ReadStats.append=Parse.parseBoolean(b);
		}else if(a.equals("overwrite") || a.equals("ow")){
			overwrite=Parse.parseBoolean(b);
		}else if(a.equals("testsize")){
			testsize=Parse.parseBoolean(b);
		}else if(a.equals("breaklen") || a.equals("breaklength")){
			breakLength=Integer.parseInt(b);
		}else if(a.equals("recalibrate") || a.equals("recalibratequality") || a.equals("recal")){
			recalibrateQuality=Parse.parseBoolean(b);
		}else if(a.equals("silent")){
			silent=Parse.parseBoolean(b);
		}else if(a.equals("wt") || a.equals("w") || a.equals("workers") || a.equals("workerthreads")){
			workers="auto".equalsIgnoreCase(b) ? -1 : Integer.parseInt(b);
		}else if(a.equals("threadsin") || a.equals("tin")){
			threadsIn="auto".equalsIgnoreCase(b) ? -1 : Integer.parseInt(b);
		}else if(a.equals("threadsout") || a.equals("tout")){
			threadsOut="auto".equalsIgnoreCase(b) ? -1 : Integer.parseInt(b);
		}else{
			return false;
		}
		return true;
	}

	/**
	 * Parses cardinality estimation parameters for LogLog algorithms.
	 * Handles k-mer cardinality tracking configuration including bucket counts,
	 * hash functions, and statistical methods.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public boolean parseCardinality(String arg, String a, String b){
		if(a.equals("cardinality") || a.equals("loglog")){
			if(b!=null && b.length()>0 && Tools.isDigit(b.charAt(0))){
				try {
					loglogk=Integer.parseInt(b);
					loglog=loglogk>0;
				} catch (NumberFormatException e) {
					loglog=Parse.parseBoolean(b);
				}
			}else{
				loglog=Parse.parseBoolean(b);
			}
		}else if(a.equals("cardinalityout") || a.equals("loglogout")){
			if(b!=null && b.length()>0 && Tools.isDigit(b.charAt(0))){
				try {
					loglogk=Integer.parseInt(b);
					loglogOut=loglogk>0;
				} catch (NumberFormatException e) {
					loglogOut=Parse.parseBoolean(b);
				}
			}else{
				loglogOut=Parse.parseBoolean(b);
			}
		}else if(a.equals("buckets") || a.equals("loglogbuckets")){
			loglogbuckets=Parse.parseIntKMG(b);
			assert(loglogbuckets>0);
			loglogbuckets=CardinalityTracker.powerOf2AtLeast(loglogbuckets);
		}else if(a.equals("loglogbits")){
			loglogbits=Integer.parseInt(b);
		}else if(a.equals("loglogk") || a.equals("cardinalityk") || a.equals("kcardinality")){
			loglogk=Integer.parseInt(b);
			loglog=loglogk>0;
		}else if(a.equals("loglogklist")){
			String[] split2=b.split(",");
			for(String k : split2){
				loglogKlist.add(Integer.parseInt(k));
			}
		}else if(a.equals("loglogseed")){
			loglogseed=Long.parseLong(b);
		}else if(a.equals("loglogminprob")){
			loglogMinprob=Float.parseFloat(b);
		}else if(a.equals("loglogtype")){
			loglogType=b;
		}else if(a.equals("loglogmean")){
			CardinalityTracker.USE_MEAN=true;
			CardinalityTracker.USE_MEDIAN=CardinalityTracker.USE_MWA=CardinalityTracker.USE_HMEAN=CardinalityTracker.USE_GMEAN=false;
		}else if(a.equals("loglogmedian")){
			CardinalityTracker.USE_MEDIAN=true;
			CardinalityTracker.USE_MEAN=CardinalityTracker.USE_MWA=CardinalityTracker.USE_HMEAN=CardinalityTracker.USE_GMEAN=false;
		}else if(a.equals("loglogmwa")){
			CardinalityTracker.USE_MWA=true;
			CardinalityTracker.USE_MEDIAN=CardinalityTracker.USE_MEAN=CardinalityTracker.USE_HMEAN=CardinalityTracker.USE_GMEAN=false;
		}else if(a.equals("logloghmean")){
			CardinalityTracker.USE_HMEAN=true;
			CardinalityTracker.USE_MEDIAN=CardinalityTracker.USE_MEAN=CardinalityTracker.USE_MWA=CardinalityTracker.USE_GMEAN=false;
		}else if(a.equals("logloggmean")){
			CardinalityTracker.USE_GMEAN=true;
			CardinalityTracker.USE_MEDIAN=CardinalityTracker.USE_MEAN=CardinalityTracker.USE_HMEAN=CardinalityTracker.USE_MWA=false;
		}else if(a.equals("loglogmantissa")){
			LogLog2.setMantissaBits(Integer.parseInt(b));
//			LogLog16.setMantissaBits(Integer.parseInt(b));
		}else if(a.equals("loglogcounts") || a.equals("loglogcount")){
			CardinalityTracker.trackCounts=Parse.parseBoolean(b);
		}else{
			return false;
		}
		
		return true;
	}
	
	/**
	 * Parses k-mer length parameter.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public boolean parseK(String arg, String a, String b){
		if(a.equalsIgnoreCase("k")){
			k=Integer.parseInt(b);
		}else{
			return false;
		}
		return true;
	}
	
	/**
	 * Parses interleaved file format parameters.
	 * Handles settings for paired-end read interleaving detection and enforcement,
	 * including auto-detection and manual override options.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public boolean parseInterleaved(String arg, String a, String b){
		if(a.equals("testinterleaved")){
			FASTQ.TEST_INTERLEAVED=Parse.parseBoolean(b);
			System.err.println("Set TEST_INTERLEAVED to "+FASTQ.TEST_INTERLEAVED);
			setInterleaved=true;
		}else if(a.equals("forceinterleaved")){
			FASTQ.FORCE_INTERLEAVED=Parse.parseBoolean(b);
			System.err.println("Set FORCE_INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			setInterleaved=true;
		}else if(a.equals("interleaved") || a.equals("int")){
			if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
			else{
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Parse.parseBoolean(b);
				System.err.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				setInterleaved=true;
			}
		}else if(a.equals("overrideinterleaved")){
			boolean x=Parse.parseBoolean(b);
			ReadStreamByteWriter.ignorePairAssertions=x;
			if(x){setInterleaved=true;}
		}else{
			return false;
		}
		return true;
	}
	
	/**
	 * Parses quality-based trimming parameters.
	 * Handles quality trimming direction (left/right/both), window-based trimming,
	 * optimal trimming algorithms, and polymer trimming options.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public boolean parseQTrim(String arg, String a, String b){
		if(a.equals("qtrim1")){
			if(b!=null && ("f".equalsIgnoreCase(b) || "false".equalsIgnoreCase(b))){qtrim1=false;}
			else{
				qtrim1=true;
				qtrim2=false;
			}
			a="qtrim";
		}else if(a.equals("qtrim2")){
			if(b!=null && ("f".equalsIgnoreCase(b) || "false".equalsIgnoreCase(b))){qtrim2=false;}
			else{
				qtrim2=true;
				qtrim1=false;
			}
			a="qtrim";
		}else if(a.equals("trimq2")){
			if(b!=null && ("f".equalsIgnoreCase(b) || "false".equalsIgnoreCase(b))){qtrim2=false;}
			else{
				qtrim2=true;
				qtrim1=false;
			}
			a="trimq";
		}
		
		if(a.equals("qtrim")/* || a.equals("trim")*/){
			if(b==null || b.length()==0){qtrimRight=qtrimLeft=true;}
			else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){qtrimLeft=true;qtrimRight=false;}
			else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){qtrimLeft=false;qtrimRight=true;}
			else if(b.equalsIgnoreCase("both") || b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr")){qtrimLeft=qtrimRight=true;}
			else if(b.equalsIgnoreCase("window") || b.equalsIgnoreCase("w") || b.startsWith("window,") || b.startsWith("w,")){
				qtrimLeft=false;
				qtrimRight=true;
				TrimRead.windowMode=true;
				TrimRead.optimalMode=false;
				String[] split=b.split(",");
				if(b.length()>1){
					TrimRead.windowLength=Integer.parseInt(split[1]);
				}
			}else if(Tools.isDigit(b.charAt(0))){
				parseTrimq(a, b);
				qtrimRight=true;
			}else{qtrimRight=qtrimLeft=Parse.parseBoolean(b);}
		}else if(a.equals("optitrim") || a.equals("otf") || a.equals("otm")){
			if(b!=null && (b.charAt(0)=='.' || Tools.isDigit(b.charAt(0)))){
				TrimRead.optimalMode=true;
				TrimRead.optimalBias=Float.parseFloat(b);
				assert(TrimRead.optimalBias>=0 && TrimRead.optimalBias<1);
			}else{
				TrimRead.optimalMode=Parse.parseBoolean(b);
			}
		}else if(a.equals("trimgoodinterval")){
			TrimRead.minGoodInterval=Integer.parseInt(b);
		}else if(a.equals("trimright") || a.equals("qtrimright")){
			qtrimRight=Parse.parseBoolean(b);
		}else if(a.equals("trimleft") || a.equals("qtrimleft")){
			qtrimLeft=Parse.parseBoolean(b);
		}else if(a.equals("trimq") || a.equals("trimquality") || a.equals("trimq2")){
			parseTrimq(a, b);
		}else if(a.equals("trimclip")){
			trimClip=Parse.parseBoolean(b);
		}else if(a.equals("trimpolya")){
			trimPolyA=parsePoly(b);
		}
		
		else if(a.equals("trimpolyg")){
			trimPolyGLeft=trimPolyGRight=parsePoly(b);
		}else if(a.equals("trimpolygleft")){
			trimPolyGLeft=parsePoly(b);
		}else if(a.equals("trimpolygright")){
			trimPolyGRight=parsePoly(b);
		}else if(a.equals("filterpolyg")){
			filterPolyG=parsePoly(b);
		}
		
		else if(a.equals("trimpolyc")){
			trimPolyCLeft=trimPolyCRight=parsePoly(b);
		}else if(a.equals("trimpolycleft")){
			trimPolyCLeft=parsePoly(b);
		}else if(a.equals("trimpolycright")){
			trimPolyCRight=parsePoly(b);
		}else if(a.equals("filterpolyc")){
			filterPolyC=parsePoly(b);
		}
		
		else if(a.equals("maxnonpoly")){
			maxNonPoly=parsePoly(b);
		}
		
		else{
			return false;
		}
		return true;
	}
	
	/**
	 * Parses polymer trimming parameter values.
	 * Converts string values to integer thresholds for polymer detection.
	 * @param b The parameter value string
	 * @return Integer threshold for polymer detection (0=disabled, >0=enabled)
	 */
	public static int parsePoly(String b){
		int r=2;
		if(b!=null){
			if(Tools.isDigit(b.charAt(0))){
				r=Integer.parseInt(b);
			}else{
				boolean x=Parse.parseBoolean(b);
				r=x ? 2 : 0;
			}
		}
		return r;
	}
	
	/**
	 * Parses general trimming and filtering parameters.
	 * Handles forced trimming, barcode filtering, read length filters,
	 * GC content filtering, and quality-based read filtering.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public boolean parseTrim(String arg, String a, String b){
		
		if(parseQTrim(arg, a, b)){
			//do nothing
		}else if(a.equals("forcetrimmod") || a.equals("forcemrimmodulo") || a.equals("ftm")){
			forceTrimModulo=Integer.parseInt(b);
		}else if(a.equals("ftl") || a.equals("forcetrimleft")){
			forceTrimLeft=Integer.parseInt(b);
		}else if(a.equals("ftr") || a.equals("forcetrimright")){
			forceTrimRight=Integer.parseInt(b);
		}else if(a.equals("ftr2") || a.equals("forcetrimright2")){
			forceTrimRight2=Integer.parseInt(b);
		}else if(a.equals("trimbadsequence")){
			trimBadSequence=Parse.parseBoolean(b);
		}else if(a.equals("chastityfilter") || a.equals("cf")){
			chastityFilter=Parse.parseBoolean(b);
		}else if(a.equals("failnobarcode")){
			failIfNoBarcode=Parse.parseBoolean(b);
		}else if(a.equals("badbarcodes") || a.equals("barcodefilter")){
			if(b!=null && (b.equalsIgnoreCase("crash") || b.equalsIgnoreCase("fail"))){
				failBadBarcodes=true;
				removeBadBarcodes=true;
			}else{
				removeBadBarcodes=Parse.parseBoolean(b);
				failBadBarcodes=false;
			}
		}else if(a.equals("barcodes") || a.equals("barcode")){
			if(b==null || b.length()<1){
				barcodes=null;
			}else{
				barcodes=new HashSet<String>();
				for(String s : b.split(",")){
					Tools.addNames(s, barcodes, false);
				}
			}
			if(barcodes!=null && barcodes.size()>0 && !failBadBarcodes && !removeBadBarcodes){
				removeBadBarcodes=true;
			}
		}else if(a.equals("requirebothbad") || a.equals("rbb")){
			requireBothBad=Parse.parseBoolean(b);
		}else if(a.equals("removeifeitherbad") || a.equals("rieb")){
			requireBothBad=!Parse.parseBoolean(b);
		}else if(a.equals("ml") || a.equals("minlen") || a.equals("minlength")){
			minReadLength=Parse.parseIntKMG(b);
		}else if(a.equals("maxlength") || a.equals("maxreadlength") || a.equals("maxreadlen") || a.equals("maxlen")){
			maxReadLength=Parse.parseIntKMG(b);
		}else if(a.equals("mingc")){
			minGC=Float.parseFloat(b);
//			if(minGC>0){filterGC=true;}
			assert(minGC>=0 && minGC<=1) : "mingc should be a decimal number between 0 and 1, inclusive.";
		}else if(a.equals("maxgc")){
			maxGC=Float.parseFloat(b);
//			if(maxGC<1){filterGC=true;}
			assert(minGC>=0 && minGC<=1) : "maxgc should be a decimal number between 0 and 1, inclusive.";
		}else if(a.equals("usepairgc") || a.equals("pairgc")){
			usePairGC=Parse.parseBoolean(b);
			ReadStats.usePairGC=usePairGC;
		}else if(a.equals("mlf") || a.equals("minlenfrac") || a.equals("minlenfraction") || a.equals("minlengthfraction")){
			minLenFraction=Float.parseFloat(b);
		}else if(a.equals("maxns")){
			maxNs=Integer.parseInt(b);
		}else if(a.equals("minconsecutivebases") || a.equals("mcb")){
			minConsecutiveBases=Integer.parseInt(b);
		}else if(a.equals("minavgquality") || a.equals("minaveragequality") || a.equals("maq")){
			if(b.indexOf(',')>-1){
				String[] split=b.split(",");
				assert(split.length==2) : "maq should be length 1 or 2 (at most 1 comma).\nFormat: maq=quality,bases; e.g. maq=10 or maq=10,20";
				minAvgQuality=Float.parseFloat(split[0]);
				minAvgQualityBases=Integer.parseInt(split[1]);
			}else{
				minAvgQuality=Float.parseFloat(b);
			}
		}else if(a.equals("minavgqualitybases") || a.equals("maqb")){
			minAvgQualityBases=Integer.parseInt(b);
		}else if(a.equals("minbasequality") || a.equals("mbq")){
			minBaseQuality=Byte.parseByte(b);
		}else if(a.equals("averagequalitybyprobability") || a.equals("aqbp")){
			Read.AVERAGE_QUALITY_BY_PROBABILITY=Parse.parseBoolean(b);
		}else if(a.equals("mintl") || a.equals("mintrimlen") || a.equals("mintrimlength")){
			minTrimLength=Integer.parseInt(b);
		}else if(a.equals("untrim")){
			untrim=Parse.parseBoolean(b);
		}else if(a.equals("tossjunk")){
			boolean x=Parse.parseBoolean(b);
			tossJunk=x;
			if(x){Read.JUNK_MODE=Read.FLAG_JUNK;}
		}else{
			return false;
		}
		return true;
	}
	
	/**
	 * Parses quality trimming threshold values.
	 * Supports single values or comma-separated arrays for position-specific thresholds.
	 * @param a The parameter name
	 * @param b The parameter value (single number or comma-separated list)
	 */
	private void parseTrimq(String a, String b){
		if(b.indexOf(',')>=0){
			String[] split=b.split(",");
			trimq2=new float[split.length];
			for(int i=0; i<split.length; i++){
				trimq2[i]=Float.parseFloat(split[i]);
			}
			trimq=trimq2.length<1 ? 0 : trimq2[0];
		}else{
			trimq=Float.parseFloat(b);
			trimq2=null;
		}
//		assert(false) : Arrays.toString(trimq2);
	}
	
	/**
	 * Parses input and output file path parameters.
	 * Handles standard input/output files, quality files, and file extensions.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value (file path)
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public boolean parseFiles(String arg, String a, String b){
		if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
			in1=b;
		}else if(a.equals("in2") || a.equals("input2")){
			in2=b;
		}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
			out1=b;
			setOut=true;
		}else if(a.equals("out2") || a.equals("output2")){
			out2=b;
			setOut=true;
		}else if(a.equals("qfin") || a.equals("qfin1")){
			qfin1=b;
		}else if(a.equals("qfout") || a.equals("qfout1")){
			qfout1=b;
			setOut=true;
		}else if(a.equals("qfin2")){
			qfin2=b;
		}else if(a.equals("qfout2")){
			qfout2=b;
			setOut=true;
		}else if(a.equals("extin")){
			extin=b;
		}else if(a.equals("extout")){
			extout=b;
		}else if(a.equals("outsingle") || a.equals("outs")){
			outsingle=b;
			setOut=true;
		}else{
			return false;
		}
		return true;
	}
	
	/**
	 * Parses alignment and mapping filter parameters.
	 * Handles identity filters, substitution/indel limits, and genome build settings.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public boolean parseMapping(String arg, String a, String b){
		if(a.equals("idfilter") || a.equals("identityfilter") || a.equals("minidfilter") || a.equals("minidentityfilter")){
			minIdFilter=Float.parseFloat(b);
			if(minIdFilter>1f){minIdFilter/=100;}
			assert(minIdFilter<=1f) : "idfilter should be between 0 and 1.";
		}else if(a.equals("maxidfilter") || a.equals("maxidentityfilter") || a.equals("maxid")){
			maxIdFilter=Float.parseFloat(b);
			if(maxIdFilter>1f){maxIdFilter/=100;}
			assert(maxIdFilter<=1f) : "idfilter should be between 0 and 1.";
		}else if(a.equals("subfilter")){
			subfilter=Integer.parseInt(b);
		}else if(a.equals("clipfilter")){
			clipfilter=Integer.parseInt(b);
		}else if(a.equals("nfilter")){
			nfilter=Integer.parseInt(b);
		}else if(a.equals("delfilter")){
			delfilter=Integer.parseInt(b);
		}else if(a.equals("insfilter")){
			insfilter=Integer.parseInt(b);
		}else if(a.equals("indelfilter")){
			indelfilter=Integer.parseInt(b);
		}else if(a.equals("dellenfilter")){
			dellenfilter=Integer.parseInt(b);
		}else if(a.equals("inslenfilter")){
			inslenfilter=Integer.parseInt(b);
		}else if(a.equals("editfilter")){
			editfilter=Integer.parseInt(b);
		}else if(a.equals("build") || a.equals("genome")){
			build=Integer.parseInt(b);
			Data.GENOME_BUILD=build;
		}else{
			return false;
		}
		return true;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Processes configuration file arguments and expands them into parameter arrays.
	 * Reads configuration files specified with "config=" parameter and inserts
	 * their contents into the argument list.
	 *
	 * @param args Original command-line arguments
	 * @return Expanded argument array with config file contents included
	 */
	static String[] parseConfig(String[] args){
		boolean found=false;
		for(String s : args){
			if(Tools.startsWithIgnoreCase(s, "config=")){
				found=true;
				break;
			}
		}
		if(!found){return args;}
		ArrayList<String> list=new ArrayList<String>();
		for(int i=0; i<args.length; i++){
			final String arg=(args[i]==null ? "null" : args[i]);
			final String[] split=arg.split("=");
			final String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			
			if(a.equals("config")){
				assert(b!=null) : "Bad parameter: "+arg;
				for(String bb : b.split(",")){
					try{
						TextFile tf=new TextFile(bb);
						for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
							String line2=line.trim();
							if(line2.length()>0 && !line2.startsWith("#")){
								list.add(line2);
							}
						}
						tf.close();
					}catch(Throwable t){
						throw new RuntimeException("Could not process config file "+b+"\nCaused by:\n"+t.toString()+"\n");
					}
				}
			}else if(arg!=null && !"null".equals(arg)){
				list.add(arg);
			}
		}
		return list.toArray(new String[list.size()]);
	}
	
	/**
	 * Parses global static parameters that affect all BBTools programs.
	 * Handles thread settings, memory options, file I/O modes, quality handling,
	 * and system-wide configuration parameters.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public static boolean parseCommonStatic(String arg, String a, String b){
		if(a.equals("null")){
			//Do nothing
		}else if(a.equals("monitor") || a.equals("killswitch")){
			if(Parse.isNumber(b)){
				String[] pair=b.split(",");
				if(pair.length==1){
					KillSwitch.launch(Double.parseDouble(pair[0]));
				}else{
					assert(pair.length==2) : "monitor takes one or two arguments, like this: monitor=600,0.002";
					KillSwitch.launch(Double.parseDouble(pair[0]), Double.parseDouble(pair[1]));
				}
			}else if(Parse.parseBoolean(b)){
				KillSwitch.launch();
			}
		}else if(a.equals("trd") || a.equals("trc") || a.equals("trimreaddescription") || a.equals("trimreaddescriptions")){
			Shared.TRIM_READ_COMMENTS=Parse.parseBoolean(b);
			if(!setTrimRname){Shared.TRIM_RNAME=Shared.TRIM_READ_COMMENTS;}
		}else if(a.equals("trimrefdescription") || a.equals("trimrefdescriptions") || a.equals("trimrname")){
			Shared.TRIM_RNAME=Parse.parseBoolean(b);
			setTrimRname=true;
		}else if(a.equals("tuc") || a.equals("touppercase")){
			Read.TO_UPPER_CASE=Parse.parseBoolean(b);
		}else if(a.equals("lctn") || a.equals("lowercaseton")){
			Read.LOWER_CASE_TO_N=Parse.parseBoolean(b);
		}else if(a.equals("changequality") || a.equals("cq")){
			Read.CHANGE_QUALITY=Parse.parseBoolean(b);
			BBMerge.changeQuality=Read.CHANGE_QUALITY;
		}else if(a.equals("tossbrokenreads") || a.equals("tbr")){
			boolean x=Parse.parseBoolean(b);
			Read.TOSS_BROKEN_QUALITY=x;
			ConcurrentReadInputStream.REMOVE_DISCARDED_READS=x;
		}else if(a.equals("nullifybrokenquality") || a.equals("nbq")){
			boolean x=Parse.parseBoolean(b);
			Read.NULLIFY_BROKEN_QUALITY=x;
		}else if(a.equals("dotdashxton")){
			boolean x=Parse.parseBoolean(b);
			Read.DOT_DASH_X_TO_N=x;
		}else if(a.equals("junk")){
			if("ignore".equalsIgnoreCase(b)){Read.JUNK_MODE=Read.IGNORE_JUNK;}
			else if("crash".equalsIgnoreCase(b) || "fail".equalsIgnoreCase(b)){Read.JUNK_MODE=Read.CRASH_JUNK;}
			else if("fix".equalsIgnoreCase(b)){Read.JUNK_MODE=Read.FIX_JUNK;}
			else if("flag".equalsIgnoreCase(b) || "discard".equalsIgnoreCase(b)){Read.JUNK_MODE=Read.FLAG_JUNK;}
			else if("iupacton".equalsIgnoreCase(b)){Read.JUNK_MODE=Read.FIX_JUNK_AND_IUPAC;}
			else{assert(false) : "Bad junk mode: "+arg;}
		}else if(a.equals("undefinedton") || a.equals("iupacton") || a.equals("itn")){
			Read.IUPAC_TO_N=Parse.parseBoolean(b);
		}else if(a.equals("ignorejunk")){
			boolean x=Parse.parseBoolean(b);
			if(x){Read.JUNK_MODE=Read.IGNORE_JUNK;}
			else if(Read.JUNK_MODE==Read.IGNORE_JUNK){Read.JUNK_MODE=Read.CRASH_JUNK;}
		}else if(a.equals("flagjunk")){
			boolean x=Parse.parseBoolean(b);
			if(x){Read.JUNK_MODE=Read.FLAG_JUNK;}
			else if(Read.JUNK_MODE==Read.FLAG_JUNK){Read.JUNK_MODE=Read.CRASH_JUNK;}
		}else if(a.equals("fixjunk")){
			boolean x=Parse.parseBoolean(b);
			if(x){Read.JUNK_MODE=Read.FIX_JUNK;}
			else if(Read.JUNK_MODE==Read.FIX_JUNK){Read.JUNK_MODE=Read.CRASH_JUNK;}
		}else if(a.equals("crashjunk") || a.equals("failjunk")){
			boolean x=Parse.parseBoolean(b);
			if(x){Read.JUNK_MODE=Read.CRASH_JUNK;}
			else if(Read.JUNK_MODE==Read.CRASH_JUNK){Read.JUNK_MODE=Read.IGNORE_JUNK;}
		}else if(a.equals("skipvalidation")){
			Read.SKIP_SLOW_VALIDATION=Parse.parseBoolean(b);
		}else if(a.equals("validate")){
			Read.SKIP_SLOW_VALIDATION=!Parse.parseBoolean(b);
		}else if(a.equals("validateinconstructor") || a.equals("vic")){
			Read.VALIDATE_IN_CONSTRUCTOR=Parse.parseBoolean(b);
		}else if(a.equals("validatebranchless")){
//			Read.VALIDATE_BRANCHLESS=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("bf1") || a.equalsIgnoreCase("bytefile1")){
			ByteFile.FORCE_MODE_BF1=Parse.parseBoolean(b);
			if(ByteFile.FORCE_MODE_BF1) {ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF4=false;}
		}else if(a.equals("bf1bufferlen")){
			ByteFile1.bufferlen=(int)Parse.parseKMGBinary(b);
		}else if(a.equalsIgnoreCase("bf2") || a.equalsIgnoreCase("bytefile2")){
			ByteFile.FORCE_MODE_BF2=Parse.parseBoolean(b);
			if(ByteFile.FORCE_MODE_BF2) {ByteFile.FORCE_MODE_BF1=ByteFile.FORCE_MODE_BF4=false;}
		}else if(a.equalsIgnoreCase("bf3") || a.equalsIgnoreCase("bytefile3")){
			ByteFile.FORCE_MODE_BF3=Parse.parseBoolean(b);
			if(ByteFile.FORCE_MODE_BF3) {ByteFile.FORCE_MODE_BF1=ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF4=false;}
		}else if(a.equalsIgnoreCase("bf4") || a.equalsIgnoreCase("bytefile4")){
			ByteFile.FORCE_MODE_BF4=Parse.parseBoolean(b);
			if(ByteFile.FORCE_MODE_BF4) {ByteFile.FORCE_MODE_BF1=ByteFile.FORCE_MODE_BF2=false;}
		}else if(a.equalsIgnoreCase("bf4threads") || a.equalsIgnoreCase("bfthreads")){
			ByteFile4.DEFAULT_THREADS=Integer.parseInt(b);
		}else if(a.equals("utot")){
			Read.U_TO_T=Parse.parseBoolean(b);
		}else if(a.equals("usejni") || a.equals("jni")){
			Shared.USE_JNI=Parse.parseBoolean(b);
		}else if(a.equals("usempi") || a.equals("mpi")){
			if(b!=null && Tools.isDigit(b.charAt(0))){
				Shared.MPI_NUM_RANKS=Integer.parseInt(b);
				Shared.USE_MPI=Shared.MPI_NUM_RANKS>0;
			}else{
				Shared.USE_MPI=Parse.parseBoolean(b);
			}
		}else if(a.equals("crismpi")){
			Shared.USE_CRISMPI=Parse.parseBoolean(b);
		}else if(a.equals("mpikeepall")){
			Shared.MPI_KEEP_ALL=Parse.parseBoolean(b);
		}else if(a.equals("readbufferlength") || a.equals("readbufferlen")){
			Shared.setBufferLen((int)Parse.parseKMG(b));
		}else if(a.equals("readbufferdata")){
			Shared.setBufferData(Parse.parseIntKMG(b));
		}else if(a.equals("readbuffers")){
			Shared.setBuffers(Integer.parseInt(b));
		}else if(a.equals("rbm") || a.equals("renamebymapping")){
			FASTQ.TAG_CUSTOM=Parse.parseBoolean(b);
		}else if(a.equals("don") || a.equals("deleteoldname")){
			FASTQ.DELETE_OLD_NAME=Parse.parseBoolean(b);
		}else if(a.equals("assertcigar")){
			ReadStreamWriter.ASSERT_CIGAR=Parse.parseBoolean(b);
		}else if(a.equals("verbosesamline")){
			SamLine.verbose=Parse.parseBoolean(b);
		}else if(a.equals("parsecustom") || a.equals("fastqparsecustom")){
			FASTQ.PARSE_CUSTOM=Parse.parseBoolean(b);
			System.err.println("Set FASTQ.PARSE_CUSTOM to "+FASTQ.PARSE_CUSTOM);
		}else if(a.equals("shrinkheaders")){
			FASTQ.SHRINK_HEADERS=Parse.parseBoolean(b);
		}else if(a.equals("fairqueues")){
			ConcurrentDepot.fair=Parse.parseBoolean(b);
		}else if(a.equals("fixheader") || a.equals("fixheaders")){
			Read.FIX_HEADER=Parse.parseBoolean(b);
		}else if(a.equals("allownullheader") || a.equals("allownullheaders")){
			Read.ALLOW_NULL_HEADER=Parse.parseBoolean(b);
		}else if(a.equals("aminoin") || a.equals("amino")){
			//TODO: ensure changes to this do not conflict with TranslateSixFrames "aain" flag.
			Shared.AMINO_IN=SketchObject.amino=Parse.parseBoolean(b);
		}else if(a.equals("amino8")){
			SketchObject.amino8=Parse.parseBoolean(b);
			if(SketchObject.amino8){
				Shared.AMINO_IN=SketchObject.amino=true;
				AminoAcid.AMINO_SHIFT=3;
			}
		}else if(a.equals("maxcalledquality")){
			int x=Tools.mid(1, Integer.parseInt(b), 93);
			Read.setMaxCalledQuality((byte)x);
		}else if(a.equals("mincalledquality")){
			int x=Tools.mid(0, Integer.parseInt(b), 93);
			Read.setMinCalledQuality((byte)x);
		}else if(a.equals("t") || a.equals("threads")){
			int old=Shared.threads();
			Shared.setThreads(b);
			if(!silent && printSetThreads && old!=Shared.threads()){
				System.err.println("Set threads to "+Shared.threads());
			}
		}else if(a.equals("recalpairnum") || a.equals("recalibratepairnum")){
			CalcTrueQuality.USE_PAIRNUM=Parse.parseBoolean(b);
		}else if(a.equals("taxpath")){
			if("auto".equalsIgnoreCase(b)){
				TaxTree.TAX_PATH=TaxTree.defaultTaxPath();
			}else{
				TaxTree.TAX_PATH=b.replaceAll("\\\\", "/");
			}
		}else if(a.equals("parallelsort") || a.equals("paralellsort")){
			boolean x=Parse.parseBoolean(b);
			Shared.setParallelSort(x);
		}else if(a.equals("gcbeforemem")){
			Shared.GC_BEFORE_PRINT_MEMORY=Parse.parseBoolean(b);
		}else if(a.equals("warnifnosequence")){
			FastaReadInputStream.WARN_IF_NO_SEQUENCE=Parse.parseBoolean(b);
		}else if(a.equals("warnfirsttimeonly")){
			FastaReadInputStream.WARN_FIRST_TIME_ONLY=Parse.parseBoolean(b);
		}else if(a.equals("silva")){
			TaxTree.SILVA_MODE=Parse.parseBoolean(b);
		}else if(a.equals("unite")){
			TaxTree.UNITE_MODE=Parse.parseBoolean(b);
		}else if(a.equals("imghq")){
			TaxTree.IMG_HQ=Parse.parseBoolean(b);
		}
		
		else if(a.equals("callins") ||  a.equals("callinss")){
			var2.Var.CALL_INS=Parse.parseBoolean(b);
		}else if(a.equals("calldel") || a.equals("calldels")){
			var2.Var.CALL_DEL=Parse.parseBoolean(b);
		}else if(a.equals("callsub") || a.equals("callsubs") || a.equals("callsnp") || a.equals("callsnps")){
			var2.Var.CALL_SUB=Parse.parseBoolean(b);
		}else if(a.equals("callindel") || a.equals("callindels")){
			var2.Var.CALL_INS=var2.Var.CALL_DEL=Parse.parseBoolean(b);
		}else if(a.equals("calljunct") || a.equals("calljunction") || a.equals("calljunctions")){
			var2.Var.CALL_JUNCTION=Parse.parseBoolean(b);
		}else if(a.equals("callnocall") || a.equals("callnocalls")){
			var2.Var.CALL_NOCALL=Parse.parseBoolean(b);
		}else if(a.equals("kmg") || a.equals("outputkmg")){
			Shared.OUTPUT_KMG=Parse.parseBoolean(b);
		}
		
		else if(a.equals("tmpdir")){
			Shared.setTmpdir(b);
		}
		
		else if(a.equals("comment")){
			Shared.comment=b;
		}
		
		else if(a.equals("fixextensions") || a.equals("fixextension") || a.equals("tryallextensions")){
			Shared.FIX_EXTENSIONS=Parse.parseBoolean(b);
		}
		
		else if(a.equals("2passresize") || a.equals("twopassresize")){
			AbstractKmerTable.TWO_PASS_RESIZE=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("forceJavaParseDouble")){
			Tools.FORCE_JAVA_PARSE_DOUBLE=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("simd")){
			if(b!=null && b.equalsIgnoreCase("auto")) {
				Shared.SIMD=(Vector.simd256);
			}else {
				Shared.SIMD=Parse.parseBoolean(b);
			}
		}else if(a.equalsIgnoreCase("simdsparse")){
			Vector.SIMD_MULT_SPARSE=Vector.SIMD_FMA_SPARSE=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("simdmultsparse")){
			Vector.SIMD_MULT_SPARSE=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("simdfmasparse")){
			Vector.SIMD_FMA_SPARSE=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("simdcopy")){
			Vector.SIMDCOPY=Parse.parseBoolean(b);
		}
//		else if(a.equalsIgnoreCase("simdminlen")) {//Disabled.
//			Vector.MINLEN=Integer.parseInt(b);
//		}
		
		else if(a.equalsIgnoreCase("awsServers") || a.equalsIgnoreCase("aws")){
			Shared.awsServers=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("nerscServers") || a.equalsIgnoreCase("nersc")){
			Shared.awsServers=!Parse.parseBoolean(b);
		}
		
		else if(a.equals("entropyk") || a.equals("ek")){
			EntropyTracker.defaultK=Integer.parseInt(b);
			EntropyTracker.setDefaultK=true;
		}else if(a.equals("entropywindow") || a.equals("ew")){
			EntropyTracker.defaultWindowBases=Integer.parseInt(b);
			EntropyTracker.setDefaultWindow=true;
		}else if(a.equalsIgnoreCase("protFull")){
			TaxTree.protFull=true;
		}
		
		else if(a.equalsIgnoreCase("lockedincrement") || a.equalsIgnoreCase("symmetricwrite")
				|| a.equalsIgnoreCase("symmetric") || a.equalsIgnoreCase("sw")){
			if("auto".equalsIgnoreCase(b)){
				KCountArray.LOCKED_INCREMENT=true;
				KCountArray.SET_LOCKED_INCREMENT=false;
			}else{
				KCountArray.LOCKED_INCREMENT=Parse.parseBoolean(b);
				KCountArray.SET_LOCKED_INCREMENT=true;
			}
		}
		
		else if(a.equalsIgnoreCase("buffer") || a.equalsIgnoreCase("buffered")){
			if(b!=null && Character.isDigit(b.charAt(0))){
				int x=Parse.parseIntKMG(b);
				KmerCountAbstract.BUFFERED=x>1;
				if(x>1){KmerCountAbstract.BUFFERLEN=x;}
			}else{
				KmerCountAbstract.BUFFERED=Parse.parseBoolean(b);
			}
		}
		
		else if(a.equalsIgnoreCase("pairreads") && b!=null){
			FASTQ.PAIR_READS=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("flipr2")){
			FASTQ.FLIP_R2=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("sidechannelstats")){
			SideChannel4.TRACK_STATS=SideChannel3.TRACK_STATS=Parse.parseBoolean(b);
		}
		
		else if(a.equals("lowmem") || a.equals("lowram") || a.equals("lowmemory")){
			boolean x=Parse.parseBoolean(b);
			shared.SyncHeart.setLowMemory(x);
			if(x){Shared.LOW_MEMORY=true;}
		}
		
//		else if(a.equalsIgnoreCase("sortserial")){
//			KmerCountAbstract.SORT_SERIAL=Parse.parseBoolean(b);
//		}
		
		else{
			return false;
		}
		return true;
	}
	
	/**
	 * Parses quality score encoding and handling parameters.
	 * Handles ASCII offset detection, Sanger vs Illumina encoding, and
	 * fake quality generation for FASTA files.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public static boolean parseQuality(String arg, String a, String b){
		parsedQuality=true; //For internal verification that this function was indeed called.
		if(a.equals("ignorebadquality") || a.equals("ibq")){
			FASTQ.IGNORE_BAD_QUALITY=Parse.parseBoolean(b);
			if(FASTQ.IGNORE_BAD_QUALITY){Read.CHANGE_QUALITY=false;}
		}else if(a.equals("ascii") || a.equals("asciioffset") || a.equals("quality") || a.equals("qual")){
			byte x;
			if(b.equalsIgnoreCase("sanger")){x=33;}
			else if(b.equalsIgnoreCase("illumina")){x=64;}
			else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;}
			else{x=(byte)Integer.parseInt(b);}
			qin=qout=x;
			FASTQ.SET_QIN=x>-1;
		}else if(a.equals("asciiin") || a.equals("qualityin") || a.equals("qualin") || a.equals("qin")){
			byte x;
			if(b.equalsIgnoreCase("sanger")){x=33;}
			else if(b.equalsIgnoreCase("illumina")){x=64;}
			else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY=true;}
			else{x=(byte)Integer.parseInt(b);}
			qin=x;
			FASTQ.SET_QIN=x>-1;
			FASTQ.DETECT_QUALITY=!FASTQ.SET_QIN;
		}else if(a.equals("asciiout") || a.equals("qualityout") || a.equals("qualout") || a.equals("qout")){
			byte x;
			if(b.equalsIgnoreCase("sanger")){x=33;}
			else if(b.equalsIgnoreCase("illumina")){x=64;}
			else if(b.equalsIgnoreCase("auto")){x=-1;FASTQ.DETECT_QUALITY_OUT=true;}
			else{x=(byte)Integer.parseInt(b);}
			qout=x;
		}else if(a.equals("fakequality") || a.equals("qfake")){
			Shared.FAKE_QUAL=Byte.parseByte(b);
		}else if(a.equals("fakefastaqual") || a.equals("fakefastaquality") || a.equals("ffq")){
			if(b==null || b.length()<1){b="f";}
			if(Character.isLetter(b.charAt(0))){
				FastaReadInputStream.FAKE_QUALITY=Parse.parseBoolean(b);
			}else{
				int x=Integer.parseInt(b);
				if(x<1){
					FastaReadInputStream.FAKE_QUALITY=false;
				}else{
					FastaReadInputStream.FAKE_QUALITY=true;
					Shared.FAKE_QUAL=(byte)Tools.min(x, 50);
				}
			}
		}else if(a.equals("qauto")){
			FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=true;
		}else{
			return false;
		}
		return true;
	}
	
	/**
	 * Checks if quality histogram files are all null.
	 * Used to determine whether quality statistics collection should be enabled.
	 * @return true if all quality histogram output files are null
	 */
	private static boolean qhistsNull(){
		return ReadStats.BQUAL_HIST_FILE==null && ReadStats.QUAL_HIST_FILE!=null && ReadStats.AVG_QUAL_HIST_FILE!=null && ReadStats.BQUAL_HIST_OVERALL_FILE!=null
				&& ReadStats.QUAL_COUNT_HIST_FILE==null;
	}
	
	/**
	 * Parses histogram output file parameters for various statistics.
	 * Handles quality histograms, length histograms, GC content histograms,
	 * and other statistical output options.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value (output file path)
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public static boolean parseHist(String arg, String a, String b){
		if(a.equals("qualityhistogram") || a.equals("qualityhist") || a.equals("qhist")){
			ReadStats.QUAL_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_QUALITY_STATS=!qhistsNull();
			if(ReadStats.COLLECT_QUALITY_STATS){System.err.println("Set quality histogram output to "+ReadStats.QUAL_HIST_FILE);}
		}else if(a.equals("basequalityhistogram") || a.equals("basequalityhist") || a.equals("bqhist")){
			ReadStats.BQUAL_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_QUALITY_STATS=!qhistsNull();
			if(ReadStats.BQUAL_HIST_FILE!=null){System.err.println("Set bquality histogram output to "+ReadStats.BQUAL_HIST_FILE);}
		}else if(a.equals("qualitycounthistogram") || a.equals("qualitycounthist") || a.equals("qchist") || a.equals("qdhist") || a.equals("qfhist")){
			ReadStats.QUAL_COUNT_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_QUALITY_STATS=!qhistsNull();
			if(ReadStats.QUAL_COUNT_HIST_FILE!=null){System.err.println("Set qcount histogram output to "+ReadStats.QUAL_COUNT_HIST_FILE);}
		}else if(a.equals("averagequalityhistogram") || a.equals("aqhist")){
			ReadStats.AVG_QUAL_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_QUALITY_STATS=!qhistsNull();
			if(ReadStats.COLLECT_QUALITY_STATS){System.err.println("Set average quality histogram output to "+ReadStats.AVG_QUAL_HIST_FILE);}
		}else if(a.equals("overallbasequalityhistogram") || a.equals("overallbasequalityhist") || a.equals("obqhist")){
			ReadStats.BQUAL_HIST_OVERALL_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_QUALITY_STATS=(ReadStats.BQUAL_HIST_FILE!=null || ReadStats.QUAL_HIST_FILE!=null || ReadStats.AVG_QUAL_HIST_FILE!=null || ReadStats.BQUAL_HIST_OVERALL_FILE!=null);
			if(ReadStats.COLLECT_QUALITY_STATS){System.err.println("Set quality histogram output to "+ReadStats.QUAL_HIST_FILE);}
		}else if(a.equals("matchhistogram") || a.equals("matchhist") || a.equals("mhist")){
			ReadStats.MATCH_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_MATCH_STATS=(ReadStats.MATCH_HIST_FILE!=null);
			if(ReadStats.COLLECT_MATCH_STATS){System.err.println("Set match histogram output to "+ReadStats.MATCH_HIST_FILE);}
		}else if(a.equals("inserthistogram") || a.equals("inserthist") || a.equals("ihist")){
			ReadStats.INSERT_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_INSERT_STATS=(ReadStats.INSERT_HIST_FILE!=null);
			if(ReadStats.COLLECT_INSERT_STATS){System.err.println("Set insert size histogram output to "+ReadStats.INSERT_HIST_FILE);}
		}else if(a.equals("basehistogram") || a.equals("basehist") || a.equals("bhist")){
			ReadStats.BASE_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_BASE_STATS=(ReadStats.BASE_HIST_FILE!=null);
			if(ReadStats.COLLECT_BASE_STATS){System.err.println("Set base content histogram output to "+ReadStats.BASE_HIST_FILE);}
		}else if(a.equals("qualityaccuracyhistogram") || a.equals("qahist")){
			ReadStats.QUAL_ACCURACY_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_QUALITY_ACCURACY=(ReadStats.QUAL_ACCURACY_FILE!=null);
			if(ReadStats.COLLECT_QUALITY_ACCURACY){
				ReadStats.COLLECT_QUALITY_STATS=true;
				ReadStats.COLLECT_MATCH_STATS=true;
				System.err.println("Set quality accuracy histogram output to "+ReadStats.QUAL_ACCURACY_FILE);
			}
		}else if(a.equals("indelhistogram") || a.equals("indelhist")){
			ReadStats.INDEL_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_INDEL_STATS=(ReadStats.INDEL_HIST_FILE!=null);
			if(ReadStats.COLLECT_INDEL_STATS){System.err.println("Set indel histogram output to "+ReadStats.INDEL_HIST_FILE);}
		}else if(a.equals("errorhistogram") || a.equals("ehist")){
			ReadStats.ERROR_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_ERROR_STATS=(ReadStats.ERROR_HIST_FILE!=null);
			if(ReadStats.COLLECT_ERROR_STATS){System.err.println("Set error histogram output to "+ReadStats.ERROR_HIST_FILE);}
		}else if(a.equals("lengthhistogram") || a.equals("lhist")){
			ReadStats.LENGTH_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_LENGTH_STATS=(ReadStats.LENGTH_HIST_FILE!=null);
			if(ReadStats.COLLECT_LENGTH_STATS){System.err.println("Set length histogram output to "+ReadStats.LENGTH_HIST_FILE);}
		}else if(a.equals("gchistogram") || a.equals("gchist")){
			ReadStats.GC_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_GC_STATS=(ReadStats.GC_HIST_FILE!=null);
			if(ReadStats.COLLECT_GC_STATS){System.err.println("Set GC histogram output to "+ReadStats.GC_HIST_FILE);}
		}else if(a.equals("gcbins") || a.equals("gchistbins")){
			if("auto".equalsIgnoreCase(b)){
				ReadStats.GC_BINS=4000;
				ReadStats.GC_BINS_AUTO=true;
			}else{
				ReadStats.GC_BINS=Integer.parseInt(b);
				ReadStats.GC_BINS_AUTO=false;
			}
		}else if(a.equals("gcchart") || a.equals("gcplot")){
			ReadStats.GC_PLOT_X=Parse.parseBoolean(b);
		}else if(a.equals("entropyhistogram") || a.equals("entropyhist") || a.equals("enhist") || a.equals("enthist")){
			ReadStats.ENTROPY_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_ENTROPY_STATS=(ReadStats.ENTROPY_HIST_FILE!=null);
			if(ReadStats.COLLECT_ENTROPY_STATS){System.err.println("Set entropy histogram output to "+ReadStats.ENTROPY_HIST_FILE);}
		}else if(a.equals("entropybins") || a.equals("entropyhistbins") || a.equals("entbins") || a.equals("enthistbins")){
			if("auto".equalsIgnoreCase(b)){
				ReadStats.ENTROPY_BINS=1000;
//				ReadStats.ENTROPY_BINS_AUTO=true;
			}else{
				ReadStats.ENTROPY_BINS=Integer.parseInt(b);
//				ReadStats.ENTROPY_BINS_AUTO=false;
			}
		}else if(a.equals("entropyns") || a.equals("entropyhistns")){
			ReadStats.allowEntropyNs=Parse.parseBoolean(b);
		}else if(a.equals("barcodestats") || a.equals("barcodecounts")){
			ReadStats.BARCODE_STATS_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_BARCODE_STATS=(ReadStats.BARCODE_STATS_FILE!=null);
		}else if(a.equals("timehistogram") || a.equals("thist")){
			ReadStats.TIME_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_TIME_STATS=(ReadStats.TIME_HIST_FILE!=null);
		}else if(a.equals("identityhistogram") || a.equals("idhist")){
			ReadStats.IDENTITY_HIST_FILE=(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none")) ? null : b;
			ReadStats.COLLECT_IDENTITY_STATS=(ReadStats.IDENTITY_HIST_FILE!=null);
			if(ReadStats.COLLECT_IDENTITY_STATS){System.err.println("Set identity histogram output to "+ReadStats.IDENTITY_HIST_FILE);}
		}else if(a.equals("idhistlen") || a.equals("idhistlength") || a.equals("idhistbins") || a.equals("idbins")){
			if("auto".equalsIgnoreCase(b)){
				ReadStats.ID_BINS=750;
				ReadStats.ID_BINS_AUTO=true;
			}else{
				ReadStats.ID_BINS=Integer.parseInt(b);
				ReadStats.ID_BINS_AUTO=false;
			}
		}
		
		else if(a.equals("maxhistlen")){
			ReadStats.MAXLEN=ReadStats.MAXINSERTLEN=ReadStats.MAXLENGTHLEN=Parse.parseIntKMG(b);
		}
		
		else if(a.equals("fixindels") || a.equals("ignorevcfindels")){
			CallVariants.fixIndels=Parse.parseBoolean(b);
		}
		
		else{
			return false;
		}
		return true;
	}

	/**
	 * Parses compression and decompression parameters.
	 * Handles gzip, bgzip, pigz, bzip2 settings including compression levels,
	 * thread counts, and tool preferences.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public static boolean parseZip(String arg, String a, String b){
		if(a.equals("ziplevel") || a.equals("zl")){
			int x=Integer.parseInt(b);
			if(x>=0){
				ReadWrite.ZIPLEVEL=Tools.min(x, 11);
				ReadWrite.ALLOW_ZIPLEVEL_CHANGE=false;
			}
		}else if(a.equals("bziplevel") || a.equals("bzl")){
			int x=Integer.parseInt(b);
			if(x>=0){
				ReadWrite.BZIPLEVEL=Tools.min(x, 9);
			}
		}else if(a.equals("allowziplevelchange")){
			ReadWrite.ALLOW_ZIPLEVEL_CHANGE=Parse.parseBoolean(b);
		}else if(a.equals("usegzip") || a.equals("gzip")){
			ReadWrite.USE_GZIP=Parse.parseBoolean(b);
		}else if(a.equals("usebgzip") || a.equals("bgzip")){
			
			if(b!=null && Tools.isDigit(b.charAt(0))){
				int zt=Integer.parseInt(b);
				if(zt<1){ReadWrite.USE_BGZIP=false;}
				else{
					ReadWrite.USE_BGZIP=true;
					ReadWrite.setZipThreads(zt);
				}
			}else{ReadWrite.USE_BGZIP=Parse.parseBoolean(b);}
			if(ReadWrite.USE_BGZIP){ReadWrite.PREFER_BGZIP=true;}
		}else if(a.equals("usebgzf") || a.equals("bgzf")){
			ReadWrite.USE_BGZF=Parse.parseBoolean(b);
		}else if(a.equals("forcepigz")){
			ReadWrite.FORCE_PIGZ=Parse.parseBoolean(b);
			if(ReadWrite.FORCE_PIGZ){ReadWrite.USE_PIGZ=true;}
		}else if(a.equals("forcebgzip")){
			ReadWrite.FORCE_BGZIP=Parse.parseBoolean(b);
			if(ReadWrite.FORCE_BGZIP){ReadWrite.USE_BGZIP=true;}
		}else if(a.equals("preferbgzip")){
			ReadWrite.PREFER_BGZIP=Parse.parseBoolean(b);
		}else if(a.equals("zipthreads")){
			ReadWrite.setZipThreads(Integer.parseInt(b));
		}else if(a.equals("usepigz") || a.equals("pigz")){
			if(b!=null && Tools.isDigit(b.charAt(0))){
				int zt=Integer.parseInt(b);
				if(zt<1){ReadWrite.USE_PIGZ=false;}
				else{
					ReadWrite.USE_PIGZ=true;
					ReadWrite.setZipThreads(zt);
				}
			}else{ReadWrite.USE_PIGZ=Parse.parseBoolean(b);}
		}else if(a.equals("zipthreaddivisor") || a.equals("ztd")){
			ReadWrite.setZipThreadMult(1/Float.parseFloat(b));
		}else if(a.equals("blocksize")){
			int x=Integer.parseInt(b);
			ReadWrite.PIGZ_BLOCKSIZE=Tools.mid(32, x, 1024);
		}else if(a.equals("pigziterations") || a.equals("pigziters")){
			int x=Integer.parseInt(b);
			ReadWrite.PIGZ_ITERATIONS=Tools.mid(32, x, 1024);
		}else if(a.equals("usegunzip") || a.equals("gunzip") || a.equals("ungzip")){
			ReadWrite.USE_GUNZIP=Parse.parseBoolean(b);
		}else if(a.equals("useunpigz") || a.equals("unpigz")){
			ReadWrite.USE_UNPIGZ=Parse.parseBoolean(b);
		}else if(a.equals("useunbgzip") || a.equals("unbgzip")){
			ReadWrite.USE_UNBGZIP=ReadWrite.PREFER_UNBGZIP=Parse.parseBoolean(b);
		}
		
		else if(a.equals("nativebgzip") || a.equals("nativebgzf")){
			ReadWrite.ALLOW_NATIVE_BGZF=Parse.parseBoolean(b);
			ReadWrite.PREFER_NATIVE_BGZF_IN=ReadWrite.PREFER_NATIVE_BGZF_OUT=Parse.parseBoolean(b);
		}else if(a.equals("usenativebgzip") || a.equals("usenativebgzf") ||
				a.equals("allownativebgzip") || a.equals("allownativebgzf")){
			ReadWrite.ALLOW_NATIVE_BGZF=Parse.parseBoolean(b);
		}else if(a.equals("nativebgzipin") || a.equals("nativebgzfin")){
			ReadWrite.PREFER_NATIVE_BGZF_IN=Parse.parseBoolean(b);
			ReadWrite.ALLOW_NATIVE_BGZF|=ReadWrite.PREFER_NATIVE_BGZF_IN;
		}else if(a.equals("nativebgzipout") || a.equals("nativebgzfout")){
			ReadWrite.PREFER_NATIVE_BGZF_OUT=Parse.parseBoolean(b);
			ReadWrite.ALLOW_NATIVE_BGZF|=ReadWrite.PREFER_NATIVE_BGZF_OUT;
		}else if(a.equals("prefernativebgzip") || a.equals("prefernativebgzf")){
			ReadWrite.PREFER_NATIVE_BGZF_IN=ReadWrite.PREFER_NATIVE_BGZF_OUT=Parse.parseBoolean(b);
			ReadWrite.ALLOW_NATIVE_BGZF|=ReadWrite.PREFER_NATIVE_BGZF_IN;
		}else if(a.equals("nativebgzipmt") || a.equals("nativebgzfmt") || a.equals("multithreadedbgzf")){
			BgzfSettings.USE_MULTITHREADED_BGZF=Parse.parseBoolean(b);
		}else if(a.equals("bgzfosmt2")){
			BgzfSettings.USE_BGZFOS_MT2=Parse.parseBoolean(b);
		}else if(a.equals("filteredbgzf")){
			BgzfOutputStreamMT.FILTERED_BGZF=Parse.parseBoolean(b);
		}else if(a.equals("bgzfthreadsin") || a.equals("bgzftin") || a.equals("bgzfreadthreads")){
			int x=Integer.parseInt(b);
			BgzfSettings.READ_THREADS=Tools.max(1, x>0 ? x : Shared.threads());
		}else if(a.equals("bgzfthreadsout") || a.equals("bgzftout") || a.equals("bgzfwritethreads")){
			int x=Integer.parseInt(b);
			BgzfSettings.WRITE_THREADS=Tools.max(1, x>0 ? x : Shared.threads());
		}
		
		else if(a.equals("preferunbgzip")){
			ReadWrite.PREFER_UNBGZIP=Parse.parseBoolean(b);
		}else if(a.equals("usebzip2") || a.equals("bzip2")){
			ReadWrite.USE_BZIP2=Parse.parseBoolean(b);
		}else if(a.equals("usepbzip2") || a.equals("pbzip2")){
			ReadWrite.USE_PBZIP2=Parse.parseBoolean(b);
		}else if(a.equals("uselbzip2") || a.equals("lbzip2")){
			ReadWrite.USE_LBZIP2=Parse.parseBoolean(b);
		}else{
			return false;
		}
		return true;
	}
	
	/**
	 * Parses SAM/BAM format parameters and output tags.
	 * Handles SAM version, tag generation options, read group settings,
	 * and alignment output formatting.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public static boolean parseSam(String arg, String a, String b){
		if(a.equals("samversion") || a.equals("samv") || a.equals("sam")){
			assert(b!=null) : "The sam flag requires a version number, e.g. 'sam=1.4'";
			SamLine.VERSION=Float.parseFloat(b);
		}else if(a.equals("sambamba")){
			Data.USE_SAMBAMBA=Parse.parseBoolean(b);
		}else if(a.equals("samtools")){
			Data.USE_SAMTOOLS=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printHeaderWait")) {
			SamReadInputStream.printHeaderWait=Parse.parseBoolean(b);
		}

		else if(a.equals("nativebam")){
			ReadWrite.ALLOW_NATIVE_BAM_OUT=ReadWrite.ALLOW_NATIVE_BAM_IN=Parse.parseBoolean(b);
			ReadWrite.PREFER_NATIVE_BAM_OUT=ReadWrite.PREFER_NATIVE_BAM_IN=Parse.parseBoolean(b);
		}else if(a.equals("usenativebam") || a.equals("allownativebam")){
			ReadWrite.ALLOW_NATIVE_BAM_OUT=ReadWrite.ALLOW_NATIVE_BAM_IN=Parse.parseBoolean(b);
		}else if(a.equals("nativebamout") || a.equals("usenativebamout")){
			ReadWrite.ALLOW_NATIVE_BAM_OUT=ReadWrite.PREFER_NATIVE_BAM_OUT=Parse.parseBoolean(b);
		}else if(a.equals("nativebamin") || a.equals("usenativebamin")){
			ReadWrite.ALLOW_NATIVE_BAM_IN=ReadWrite.PREFER_NATIVE_BAM_IN=Parse.parseBoolean(b);
		}else if(a.equals("prefernativebamout")){
			ReadWrite.PREFER_NATIVE_BAM_OUT=Parse.parseBoolean(b);
			ReadWrite.ALLOW_NATIVE_BAM_OUT|=ReadWrite.PREFER_NATIVE_BAM_OUT;
		}else if(a.equals("prefernativebamin")){
			ReadWrite.PREFER_NATIVE_BAM_IN=Parse.parseBoolean(b);
			ReadWrite.ALLOW_NATIVE_BAM_IN|=ReadWrite.PREFER_NATIVE_BAM_IN;
		}else if(a.equals("prefernativebam")){
			ReadWrite.PREFER_NATIVE_BAM_IN=ReadWrite.PREFER_NATIVE_BAM_OUT=Parse.parseBoolean(b);
			ReadWrite.ALLOW_NATIVE_BAM_IN|=ReadWrite.PREFER_NATIVE_BAM_IN;
			ReadWrite.ALLOW_NATIVE_BAM_OUT|=ReadWrite.PREFER_NATIVE_BAM_OUT;
		}else if(a.equals("userssw")){
			ReadWrite.USE_READ_STREAM_SAM_WRITER=Parse.parseBoolean(b);
		}

		else if(a.equals("attachedsamline") || a.equals("useattachedsamline")){
			ReadStreamWriter.USE_ATTACHED_SAMLINE=Parse.parseBoolean(b);
		}else if(a.equals("samtools")){
			Data.USE_SAMTOOLS=Parse.parseBoolean(b);
		}else if(a.equals("streamerthreads") || a.equals("ssthreads")){
			SamStreamer.DEFAULT_THREADS=Integer.parseInt(b);
		}else if(a.equals("bsthreads")){
			BamStreamer.DEFAULT_THREADS=Integer.parseInt(b);
		}else if(a.equals("fastqstreamerthreads") || a.equals("fqsthreads")){
			FastqStreamer.DEFAULT_THREADS=Integer.parseInt(b);
		}else if(a.equals("fastastreamerthreads") || a.equals("fasthreads")){
			FastaStreamer.DEFAULT_THREADS=Integer.parseInt(b);
		}else if(a.equals("samwriterthreads") || a.equals("swthreads")){
			SamWriter.DEFAULT_THREADS=Integer.parseInt(b);
		}else if(a.equals("bamwriterthreads") || a.equals("bwthreads")){
			BamWriter.DEFAULT_THREADS=Integer.parseInt(b);
		}else if(a.equals("fastqwriterthreads") || a.equals("fqwthreads")){
			FastqWriter.DEFAULT_THREADS=Integer.parseInt(b);
		}else if(a.equals("fastastreamer2")){
			StreamerFactory.FASTA_STREAMER_2=Parse.parseBoolean(b);
		}else if(a.equals("prefermd") || a.equals("prefermdtag")){
			SamLine.PREFER_MDTAG=Parse.parseBoolean(b);
		}else if(a.equals("notags")){
			SamLine.NO_TAGS=Parse.parseBoolean(b);
		}else if(a.equals("mdtag") || a.equals("md")){
			SamLine.MAKE_MD_TAG=Parse.parseBoolean(b);
		}else if(a.equals("idtag")){
			SamLine.MAKE_IDENTITY_TAG=Parse.parseBoolean(b);
		}else if(a.equals("mateqtag")){
			SamLine.MAKE_MATEQ_TAG=Parse.parseBoolean(b);
		}else if(a.equals("xmtag") || a.equals("xm")){
			SamLine.MAKE_XM_TAG=Parse.parseBoolean(b);
		}else if(a.equals("smtag")){
			SamLine.MAKE_SM_TAG=Parse.parseBoolean(b);
		}else if(a.equals("amtag")){
			SamLine.MAKE_AM_TAG=Parse.parseBoolean(b);
		}else if(a.equals("nmtag")){
			SamLine.MAKE_NM_TAG=Parse.parseBoolean(b);
		}else if(a.equals("xttag")){
			SamLine.MAKE_XT_TAG=Parse.parseBoolean(b);
		}else if(a.equals("stoptag")){
			SamLine.MAKE_STOP_TAG=Parse.parseBoolean(b);
		}else if(a.equals("lengthtag")){
			SamLine.MAKE_LENGTH_TAG=Parse.parseBoolean(b);
		}else if(a.equals("boundstag")){
			SamLine.MAKE_BOUNDS_TAG=Parse.parseBoolean(b);
		}else if(a.equals("scoretag")){
			SamLine.MAKE_SCORE_TAG=Parse.parseBoolean(b);
		}else if(a.equals("sortscaffolds")){
			SamLine.SORT_SCAFFOLDS=Parse.parseBoolean(b);
		}else if(a.equals("customtag")){
			SamLine.MAKE_CUSTOM_TAGS=Parse.parseBoolean(b);
		}else if(a.equals("nhtag")){
			SamLine.MAKE_NH_TAG=Parse.parseBoolean(b);
		}else if(a.equals("keepnames")){
			SamLine.KEEP_NAMES=Parse.parseBoolean(b);
		}else if(a.equals("saa") || a.equals("secondaryalignmentasterisks")){
			SamLine.SECONDARY_ALIGNMENT_ASTERISKS=Parse.parseBoolean(b);
		}else if(a.equals("inserttag")){
			SamLine.MAKE_INSERT_TAG=Parse.parseBoolean(b);
		}else if(a.equals("correctnesstag")){
			SamLine.MAKE_CORRECTNESS_TAG=Parse.parseBoolean(b);
		}else if(a.equals("intronlen") || a.equals("intronlength")){
			SamLine.INTRON_LIMIT=Integer.parseInt(b);
			SamLine.setintron=true;
		}else if(a.equals("suppressheader") || a.equals("noheader")){
			ReadStreamWriter.NO_HEADER=Parse.parseBoolean(b);
		}else if(a.equals("noheadersequences") || a.equals("nhs") || a.equals("suppressheadersequences")){
			ReadStreamWriter.NO_HEADER_SEQUENCES=Parse.parseBoolean(b);
		}else if(a.equals("tophat")){
			if(Parse.parseBoolean(b)){
				SamLine.MAKE_TOPHAT_TAGS=true;
				FastaReadInputStream.FAKE_QUALITY=true;
				Shared.FAKE_QUAL=40;
				SamLine.MAKE_MD_TAG=true;
			}
		}else if(a.equals("xstag") || a.equals("xs")){
			SamLine.MAKE_XS_TAG=true;
			if(b!=null){
				b=b.toLowerCase();
				if(b.startsWith("fr-")){b=b.substring(3);}
				if(b.equals("ss") || b.equals("secondstrand")){
					SamLine.XS_SECONDSTRAND=true;
				}else if(b.equals("fs") || b.equals("firststrand")){
					SamLine.XS_SECONDSTRAND=false;
				}else if(b.equals("us") || b.equals("unstranded")){
					SamLine.XS_SECONDSTRAND=false;
				}else{
					SamLine.MAKE_XS_TAG=Parse.parseBoolean(b);
				}
			}
			SamLine.setxs=true;
		}else if(a.equals("flipsam")){
			SamLine.FLIP_ON_LOAD=Parse.parseBoolean(b);
		}else if(parseReadgroup(arg, a, b)){
			//do nothing
		}else{
			return false;
		}
		return true;
	}

	/**
	 * Parses FASTA format-specific parameters.
	 * Handles read length splitting, minimum read lengths, and FASTA output formatting.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public static boolean parseFasta(String arg, String a, String b){
		if(a.equals("fastareadlen") || a.equals("fastareadlength")){
			FastaReadInputStream.TARGET_READ_LEN=Integer.parseInt(b);
			FastaReadInputStream.SPLIT_READS=(FastaReadInputStream.TARGET_READ_LEN>0);
		}else if(a.equals("fastaminread") || a.equals("fastaminlen") || a.equals("fastaminlength")){
			FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
		}else if(a.equals("forcesectionname")){
			FastaReadInputStream.FORCE_SECTION_NAME=Parse.parseBoolean(b);
		}else if(a.equals("fastawrap") || a.equals("wrap")){
			Shared.FASTA_WRAP=Parse.parseIntKMG(b);
		}else if(a.equals("fastadump")){
			AbstractKmerTable.FASTA_DUMP=Parse.parseBoolean(b);
		}else{
			return false;
		}
		return true;
	}
	
	/**
	 * Parses quality score recalibration parameters.
	 * Handles quality matrix loading, recalibration passes, and statistical methods
	 * for improving quality score accuracy.
	 *
	 * @param arg The complete argument string
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public static boolean parseQualityAdjust(String arg, String a, String b){
		int pass=0;
		if(a.endsWith("_p1") || a.endsWith("_p2")){
			pass=Integer.parseInt(a.substring(a.length()-1))-1;
			a=a.substring(0, a.length()-3);
		}
		
		if(a.equals("trackall")){
			CalcTrueQuality.TRACK_ALL=Parse.parseBoolean(b);
		}else if(a.equals("clearmatrices")){
			boolean x=Parse.parseBoolean(b);
			if(x){
				CalcTrueQuality.use_q102=new boolean[] {false, false};
				CalcTrueQuality.use_qap=new boolean[] {false, false};
				CalcTrueQuality.use_qbp=new boolean[] {false, false};
				CalcTrueQuality.use_qpt=new boolean[] {false, false};
				CalcTrueQuality.use_qbt=new boolean[] {false, false};
				CalcTrueQuality.use_q10=new boolean[] {false, false};
				CalcTrueQuality.use_q12=new boolean[] {false, false};
				CalcTrueQuality.use_qb12=new boolean[] {false, false};
				CalcTrueQuality.use_qb012=new boolean[] {false, false};
				CalcTrueQuality.use_qb123=new boolean[] {false, false};
				CalcTrueQuality.use_qb234=new boolean[] {false, false};
				CalcTrueQuality.use_q12b12=new boolean[] {false, false};
				CalcTrueQuality.use_qp=new boolean[] {false, false};
				CalcTrueQuality.use_q=new boolean[] {false, false};
			}
		}else if(a.equals("loadq102")){
			CalcTrueQuality.use_q102[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadqap")){
			CalcTrueQuality.use_qap[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadqbp")){
			CalcTrueQuality.use_qbp[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadqpt")){
			CalcTrueQuality.use_qpt[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadqbt")){
			CalcTrueQuality.use_qbt[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadq10")){
			CalcTrueQuality.use_q10[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadq12")){
			CalcTrueQuality.use_q12[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadqb12")){
			CalcTrueQuality.use_qb12[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadqb012")){
			CalcTrueQuality.use_qb012[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadqb123")){
			CalcTrueQuality.use_qb123[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadqb234")){
			CalcTrueQuality.use_qb234[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadq12b12")){
			CalcTrueQuality.use_q12b12[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadqp")){
			CalcTrueQuality.use_qp[pass]=Parse.parseBoolean(b);
		}else if(a.equals("loadq")){
			CalcTrueQuality.use_q[pass]=Parse.parseBoolean(b);
		}else if(a.equals("observationcutoff")){
			long x=Parse.parseIntKMG(b);
			CalcTrueQuality.OBSERVATION_CUTOFF[pass]=x;
		}else if(a.equals("recalpasses")){
			CalcTrueQuality.passes=Integer.parseInt(b);
		}else if(a.equals("recalqmax")){
			int x=Tools.mid(1, Integer.parseInt(b), 93);
			CalcTrueQuality.setQmax(x);
			Read.setMaxCalledQuality(Tools.max(x, Read.MAX_CALLED_QUALITY()));
		}else if(a.equals("recalqmin")){
			int x=Tools.mid(0, Integer.parseInt(b), 93);
			Read.setMinCalledQuality(Tools.min(x, Read.MIN_CALLED_QUALITY()));
		}else if(a.equals("recalwithposition") || a.equals("recalwithpos") || a.equals("recalusepos")){
			boolean x=Parse.parseBoolean(b);
			if(!x){
				Arrays.fill(CalcTrueQuality.use_qp, false);
				Arrays.fill(CalcTrueQuality.use_qbp, false);
				Arrays.fill(CalcTrueQuality.use_qap, false);
			}
		}else if(a.equals("qmatrixmode")){
			if("smr".equalsIgnoreCase(b) || "squaremeanroot".equalsIgnoreCase(b)){
				CalcTrueQuality.USE_SMR=true;
			}else if("weighted".equalsIgnoreCase(b) || "wa".equalsIgnoreCase(b) || 
					"weightedaverage".equalsIgnoreCase(b)){
				CalcTrueQuality.USE_SMR=false;
				CalcTrueQuality.USE_WEIGHTED_AVERAGE=true;
			}else if("average".equalsIgnoreCase(b) || "avg".equalsIgnoreCase(b)){
				CalcTrueQuality.USE_SMR=CalcTrueQuality.USE_WEIGHTED_AVERAGE=false;
				CalcTrueQuality.USE_AVERAGE=true;
			}else if("max".equalsIgnoreCase(b)){
				CalcTrueQuality.USE_SMR=false;
				CalcTrueQuality.USE_AVERAGE=CalcTrueQuality.USE_WEIGHTED_AVERAGE=false;
			}
		}else if(a.equals("recaltile") || a.equals("recaltiles") || a.equals("usetiles")){
			CalcTrueQuality.USE_TILES=Parse.parseBoolean(b);
			if(CalcTrueQuality.USE_TILES) {
				CalcTrueQuality.use_qpt[0]=true;
				CalcTrueQuality.use_qpt[1]=false;
				CalcTrueQuality.use_qbt[0]=false;
				CalcTrueQuality.use_qbt[1]=false;
			}else {
				CalcTrueQuality.use_qpt[0]=false;
				CalcTrueQuality.use_qpt[1]=false;
				CalcTrueQuality.use_qbt[0]=false;
				CalcTrueQuality.use_qbt[1]=false;
			}
		}else{
			return false;
		}
		return true;
	}

	/**
	 * Determines if an argument is a JVM flag rather than a program parameter.
	 * Recognizes memory settings, assertion flags, and other JVM options.
	 * @param arg The argument to check
	 * @return true if the argument is a JVM flag, false otherwise
	 */
	static boolean isJavaFlag(String arg){
		if(arg==null){return false;}
		if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.startsWith("-Xmn") || arg.startsWith("-xmx") || arg.startsWith("-xms") || arg.startsWith("-xmn")){
			return arg.length()>4 && Tools.isDigit(arg.charAt(4));
		}
		if(arg.startsWith("Xmx") || arg.startsWith("Xms") || arg.startsWith("Xmn") || arg.startsWith("xmx")){
			return arg.length()>3 && (Tools.isDigit(arg.charAt(3)) || arg.charAt(3)=='=');
		}
		if(arg.equals("-ea") || arg.equals("-da") || arg.equals("ea") || arg.equals("da")){
			return true;
		}
		if(arg.equals("ExitOnOutOfMemoryError") || arg.equals("exitonoutofmemoryerror") || arg.equals("eoom")){return true;}
		if(arg.equals("-ExitOnOutOfMemoryError") || arg.equals("-exitonoutofmemoryerror") || arg.equals("-eoom")){return true;}
		
		return false;
	}
	
	/** Return true if the user seems confused */
	static boolean parseHelp(String[] args, boolean autoExit){
		if(args==null || args.length==0 || (args.length==1 && args[0]==null)){
			if(autoExit){printHelp(1);}
			return true;
		}
		
		final String s=args[args.length-1].toLowerCase();
		
		if(s.equals("-version") || s.equals("--version") || (s.equals("version") && !new File(s).exists())){
			if(autoExit){printHelp(0);}
			return true;
		}else if(s.equals("-h") || s.equals("-help") || s.equals("--help")
				|| s.equals("?") || s.equals("-?") || (s.equals("help") && !new File(s).exists())){
			if(autoExit){printHelp(0);}
			return true;
		}
		return false;
	}
	
	/** Prints version information and help message, then exits.
	 * @param exitCode Exit code for System.exit() */
	public static void printHelp(int exitCode){
		System.err.println("BBTools version "+Shared.BBTOOLS_VERSION_STRING);
		System.err.println("For help, please run the shellscript with no parameters, or look in /docs/.");
		System.exit(exitCode);
	}
	
	/** Set SamLine Readgroup Strings */
	public static boolean parseReadgroup(String arg, String a, String b){
		if(a.equals("readgroup") || a.equals("readgroupid") || a.equals("rgid")){
			SamLine.READGROUP_ID=b;
			if(b!=null){SamLine.READGROUP_TAG="RG:Z:"+b;}
		}else if(a.equals("readgroupcn") || a.equals("rgcn")){
			SamLine.READGROUP_CN=b;
		}else if(a.equals("readgroupds") || a.equals("rgds")){
			SamLine.READGROUP_DS=b;
		}else if(a.equals("readgroupdt") || a.equals("rgdt")){
			SamLine.READGROUP_DT=b;
		}else if(a.equals("readgroupfo") || a.equals("rgfo")){
			SamLine.READGROUP_FO=b;
		}else if(a.equals("readgroupks") || a.equals("rgks")){
			SamLine.READGROUP_KS=b;
		}else if(a.equals("readgrouplb") || a.equals("rglb")){
			SamLine.READGROUP_LB=b;
		}else if(a.equals("readgrouppg") || a.equals("rgpg")){
			SamLine.READGROUP_PG=b;
		}else if(a.equals("readgrouppi") || a.equals("rgpi")){
			SamLine.READGROUP_PI=b;
		}else if(a.equals("readgrouppl") || a.equals("rgpl")){
			SamLine.READGROUP_PL=b;
		}else if(a.equals("readgrouppu") || a.equals("rgpu")){
			SamLine.READGROUP_PU=b;
		}else if(a.equals("readgroupsm") || a.equals("rgsm")){
			SamLine.READGROUP_SM=b;
		}else{
			return false;
		}
		return true;
	}
	
	/** Fix Readgroup Strings */
	public static void postparseReadgroup(String fname){
		if(fname==null){return;}
		fname=ReadWrite.stripToCore(fname);
		if("filename".equalsIgnoreCase(SamLine.READGROUP_ID)){
			SamLine.READGROUP_ID=fname;
			if(fname!=null){SamLine.READGROUP_TAG="RG:Z:"+fname;}
		}
		if("filename".equalsIgnoreCase(SamLine.READGROUP_CN)){
			SamLine.READGROUP_CN=fname;
		}
		if("filename".equalsIgnoreCase(SamLine.READGROUP_DS)){
			SamLine.READGROUP_DS=fname;
		}
		if("filename".equalsIgnoreCase(SamLine.READGROUP_DT)){
			SamLine.READGROUP_DT=fname;
		}
		if("filename".equalsIgnoreCase(SamLine.READGROUP_FO)){
			SamLine.READGROUP_FO=fname;
		}
		if("filename".equalsIgnoreCase(SamLine.READGROUP_KS)){
			SamLine.READGROUP_KS=fname;
		}
		if("filename".equalsIgnoreCase(SamLine.READGROUP_LB)){
			SamLine.READGROUP_LB=fname;
		}
		if("filename".equalsIgnoreCase(SamLine.READGROUP_PG)){
			SamLine.READGROUP_PG=fname;
		}
		if("filename".equalsIgnoreCase(SamLine.READGROUP_PI)){
			SamLine.READGROUP_PI=fname;
		}
		if("filename".equalsIgnoreCase(SamLine.READGROUP_PL)){
			SamLine.READGROUP_PL=fname;
		}
		if("filename".equalsIgnoreCase(SamLine.READGROUP_PU)){
			SamLine.READGROUP_PU=fname;
		}
		if("filename".equalsIgnoreCase(SamLine.READGROUP_SM)){
			SamLine.READGROUP_SM=fname;
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Converts quality trimming threshold to error probability.
	 * @return Error probability corresponding to the quality trimming threshold */
	public float trimE(){
		return (float)QualityTools.phredToProbError(trimq);
	}
	
	/**
	 * Converts quality trimming thresholds to error probabilities.
	 * Handles both single threshold and position-specific threshold arrays.
	 * @return Array of error probabilities for quality trimming
	 */
	public float[] trimE2(){
		return QualityTools.phredToProbError(trimq2==null ? new float[] {trimq} : trimq2);
	}

	/** Enable LogLog cardinality estimation */
	public boolean loglog=false;
	/** Output LogLog cardinality estimates to file */
	public boolean loglogOut=false;
	/** Number of buckets for LogLog cardinality estimation */
	public int loglogbuckets=2048;//1999
	/** Number of bits per bucket for LogLog estimation */
	public int loglogbits=8;
	/** K-mer length for LogLog cardinality estimation */
	public int loglogk=31;
	/** Random seed for LogLog hash functions */
	public long loglogseed=-1;
	/** Minimum probability threshold for LogLog counting */
	public float loglogMinprob=0;
	/** List of k-mer lengths for multi-k LogLog estimation */
	public IntList loglogKlist=new IntList();
	
	/** Enable quality score recalibration */
	public boolean recalibrateQuality=false;
	
	/** Force trimming to multiples of this length (-1 to disable) */
	public int forceTrimModulo=-1;
	/** Force trimming of bases from left end (-1 to disable) */
	public int forceTrimLeft=-1;
	/** Force trimming of bases from right end (-1 to disable) */
	public int forceTrimRight=-1;
	/** Force trimming of bases from right end of read 2 (-1 to disable) */
	public int forceTrimRight2=-1;
	/** Genome build number for coordinate system */
	public int build=1;

	/** Maximum number of reads to process (-1 for unlimited) */
	public long maxReads=-1;
	/** Fraction of reads to randomly sample (1.0 for all reads) */
	public float samplerate=1f;
	/** Random seed for read sampling (-1 for random seed) */
	public long sampleseed=-1;

	/** Enable quality trimming from left end of reads */
	public boolean qtrimLeft=false;
	/** Enable quality trimming from right end of reads */
	public boolean qtrimRight=false;
	/** Enable adapter clipping during trimming */
	public boolean trimClip=false;
	/** Minimum poly-A tail length to trim (0 to disable) */
	public int trimPolyA=0;
	
	/** Minimum poly-G length to trim from left end (0 to disable) */
	public int trimPolyGLeft=0;
	/** Minimum poly-G length to trim from right end (0 to disable) */
	public int trimPolyGRight=0;
	/** Minimum poly-G length to filter entire read (0 to disable) */
	public int filterPolyG=0;
	
	/** Minimum poly-C length to trim from left end (0 to disable) */
	public int trimPolyCLeft=0;
	/** Minimum poly-C length to trim from right end (0 to disable) */
	public int trimPolyCRight=0;
	/** Minimum poly-C length to filter entire read (0 to disable) */
	public int filterPolyC=0;
	/** Maximum non-polymer bases allowed in polymer detection */
	public int maxNonPoly=1;

	/** Apply quality trimming only to read 1 of pairs */
	public boolean qtrim1=false;
	/** Apply quality trimming only to read 2 of pairs */
	public boolean qtrim2=false;

	/** Quality threshold for trimming (Phred score) */
	public float trimq=6;
	/** Position-specific quality thresholds for trimming */
	public float[] trimq2=null;
	/** Minimum average quality score to retain read */
	public float minAvgQuality=0;
	/** Minimum individual base quality score required */
	public byte minBaseQuality=0;
	/** Number of bases to use for average quality calculation */
	public int minAvgQualityBases=0;
	/** Maximum ambiguous (N) bases allowed per read (-1 for unlimited) */
	public int maxNs=-1;
	/** Minimum consecutive non-N bases required */
	public int minConsecutiveBases=0;
	/** Minimum read length after trimming to retain */
	public int minReadLength=0;
	/** Maximum read length allowed (-1 for unlimited) */
	public int maxReadLength=-1;
	/** Minimum length after trimming (-1 to use minReadLength) */
	public int minTrimLength=-1;
	/** Minimum length as fraction of original read length */
	public float minLenFraction=0;
	/** Minimum GC content fraction to retain read */
	public float minGC=0;
	/** Maximum GC content fraction to retain read */
	public float maxGC=1;
	/** Use combined pair GC content for filtering rather than individual reads */
	public boolean usePairGC=true;
//	public boolean filterGC=false;
	/** Restore original read length by padding with Ns */
	public boolean untrim=false;
	/** Discard reads flagged as junk by quality checks */
	public boolean tossJunk=false;

	/** Minimum alignment identity to retain read (-1 to disable) */
	public float minIdFilter=-1;
	/** Maximum alignment identity to retain read */
	public float maxIdFilter=999999999;
	/** Maximum substitutions allowed in alignment (-1 to disable) */
	public int subfilter=-1;
	/** Maximum clipped bases allowed in alignment (-1 to disable) */
	public int clipfilter=-1;
	/** Maximum deletions allowed in alignment (-1 to disable) */
	public int delfilter=-1;
	/** Maximum insertions allowed in alignment (-1 to disable) */
	public int insfilter=-1;
	/** Maximum indels (insertions + deletions) allowed (-1 to disable) */
	public int indelfilter=-1;
	/** Maximum deletion length allowed in alignment (-1 to disable) */
	public int dellenfilter=-1;
	/** Maximum insertion length allowed in alignment (-1 to disable) */
	public int inslenfilter=-1;
	/** Maximum edit distance allowed in alignment (-1 to disable) */
	public int editfilter=-1;
	/** Maximum N bases allowed in alignment (-1 to disable) */
	public int nfilter=-1;
	
	/** Length threshold for breaking long reads into fragments */
	public int breakLength=0;
	/** Toss pair only if both reads are shorter than limit */
	public boolean requireBothBad=false;
	/** Trim sequences identified as adapters or contaminants */
	public boolean trimBadSequence=false;
	/** Apply Illumina chastity filtering to remove low-quality reads */
	public boolean chastityFilter=false;
	/** Remove reads with invalid or unrecognized barcodes */
	public boolean removeBadBarcodes=false;
	/** Fail program execution if bad barcodes are encountered */
	public boolean failBadBarcodes=false;
	/** Fail program execution if reads lack expected barcodes */
	public boolean failIfNoBarcode=false;
	
	/** Set of valid barcode sequences for filtering */
	public HashSet<String> barcodes=null;
	
	/** Permission to overwrite existing files */
	public boolean overwrite=true;
	
	/** Permission to append to existing files */
	public boolean append=false;
	/** Test and report file sizes without processing */
	public boolean testsize=false;
	
	/** Whether input file interleaving was explicitly set */
	public boolean setInterleaved=false;
	
	/** Primary input file path */
	public String in1=null;
	/** Secondary input file path for paired reads */
	public String in2=null;
	
	/** Primary quality file input path */
	public String qfin1=null;
	/** Secondary quality file input path */
	public String qfin2=null;

	/** Primary output file path */
	public String out1=null;
	/** Secondary output file path for paired reads */
	public String out2=null;
	/** Output file path for unpaired reads after filtering */
	public String outsingle=null;
	/** Whether output files were explicitly configured */
	public boolean setOut=false;

	/** Primary quality file output path */
	public String qfout1=null;
	/** Secondary quality file output path */
	public String qfout2=null;
	
	/** Default extension for input files */
	public String extin=null;
	/** Default extension for output files */
	public String extout=null;
	
	/** K-mer length for sequence analysis */
	public int k=31;

	public int workers() {return workers<1 ? Shared.threads() : workers;}
	public int workers=-1;
	public int threadsIn=-1;
	public int threadsOut=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Type of LogLog algorithm to use for cardinality estimation */
	public static String loglogType="LogLog2";
	/** Suppress progress and status messages */
	public static boolean silent=false;
	/** Print thread count changes to stderr */
	public static boolean printSetThreads=true;
	
	/** Whether read name trimming was explicitly set */
	private static boolean setTrimRname=false;
	/** Input quality score ASCII offset */
	private static byte qin=-1;
	/** Output quality score ASCII offset */
	private static byte qout=-1;
	/** Whether quality parameters have been parsed */
	private static boolean parsedQuality=false;
	
	/** Sets the quality score ASCII offset for both input and output.
	 * @param x The ASCII offset value (-1 for auto-detection) */
	public static void setQuality(int x){
		qin=(byte)x;
		parsedQuality=(x>-1);
	}
	/**
	 * Applies parsed quality settings to global FASTQ configuration.
	 * Sets ASCII offsets and enables/disables quality detection based on
	 * previously parsed parameters.
	 */
	public static void processQuality(){
//		assert(parsedQuality);
		if(!parsedQuality){return;}
		if(qin!=-1 && qout!=-1){
			FASTQ.ASCII_OFFSET=qin;
			FASTQ.ASCII_OFFSET_OUT=qout;
			FASTQ.DETECT_QUALITY=false;
		}else if(qin!=-1){
			FASTQ.ASCII_OFFSET=qin;
			FASTQ.DETECT_QUALITY=false;
		}else if(qout!=-1){
			FASTQ.ASCII_OFFSET_OUT=qout;
			FASTQ.DETECT_QUALITY_OUT=false;
		}
	}

	/**
	 * Validates multiple file formats for standard input/output compatibility.
	 * @param ffa Array of FileFormat objects to validate
	 * @return true if all formats are valid for stdio, false otherwise
	 */
	public boolean validateStdio(FileFormat... ffa) {
		boolean b=true;
		for(FileFormat ff:ffa){
			b=validateStdio(ff)&b;
		}
		return b;
	}

	/**
	 * Validates a file format for standard input/output compatibility.
	 * Checks that required format information is available when using stdin/stdout.
	 * @param ff The FileFormat to validate
	 * @return true if the format is valid for stdio, false otherwise
	 */
	public boolean validateStdio(FileFormat ff) {
		if(ff==null || !ff.stdio()){return true;}
		if(ff.fastq() && ff.stdin()){
			assert(setInterleaved) : "\nERROR: When piping fastq data from stdin, interleaving must be explicitly stated\n"
					+ "with the flag int=f for unpaired data or int=t for paired data.\n";
		}
		final int ext=ff.rawExtensionCode();
		if(ff.stdout() && ext==FileFormat.UNKNOWN){
			assert(false) : "\nERROR: When piping reads to stdout, the output format must be specified with an extension,\n"
					+ "such as stdout.fq or stdout.sam.gz.\n";
		}
		if(ff.stdin() && ext==FileFormat.UNKNOWN){
			assert(false) : "\nERROR: When piping reads from stdin, the input format should be specified with an extension,\n"
					+ "such as stdin.fq or stdin.sam.gz.\n";
		}
		return true;
	}
	
}
