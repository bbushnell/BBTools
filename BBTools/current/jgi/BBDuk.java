package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.atomic.AtomicLongArray;

import aligner.SideChannel3;
import cardinality.CardinalityTracker;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import hiseq.FlowcellCoordinate;
import json.JsonObject;
import kmer.AbstractKmerTable;
import kmer.AbstractKmerTableSet;
import kmer.ScheduleMaker;
import shared.KillSwitch;
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
import var2.VcfLoader;

/**
 * Separates, trims, or masks sequences based on matching kmers in a reference.
 * Supports Hamming and and edit distance.
 * Supports K 1-31 and emulated K>31.
 * @author Brian Bushnell
 * @date Aug 30, 2013
 *
 */
public class BBDuk {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
//		if(PreParser.isAmino(args)){
//			BBDukAA.main(args);//No longer needed since this supports AAs
//			return;
//		}
		
		//Create a new BBDuk instance
		BBDuk x=new BBDuk(args);
		
		//And run it
		x.process();
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
		assert(WAYS==7 && !verbose) : "Undo WAYS and verbose";
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public BBDuk(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), true);
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
		
		/* Initialize local variables with defaults */
		boolean setOut=false, setOutb=false;
		boolean ktrimRight_=false, ktrimLeft_=false, ktrimN_=false, ktrimExclusive_=false, ksplit_=false;
		boolean findBestMatch_=false;
		boolean addTrimmedToBad_=true;
		boolean rcomp_=true;
		boolean forbidNs_=false;
		boolean useForest_=false, useTable_=false, useArray_=true, prealloc_=false;
		int k_=27, kbig_=-1;
		int mink_=-1;
		int ways_=-1; //Currently disabled for speed
		int maxBadKmers_=0;
		long skipreads_=0;
		byte TRIM_SYMBOL_='N';
		boolean kmaskLowercase_=false;
		boolean kmaskFullyCovered_=false;
		boolean histogramsBeforeProcessing_=true;
		boolean trimFailuresTo1bp_=false;
		
		Parser parser=new Parser();
		parser.trimq=6;
		parser.minAvgQuality=0;
		parser.minReadLength=10;
		parser.maxReadLength=Integer.MAX_VALUE;
		parser.minLenFraction=0f;
		parser.requireBothBad=false;
		parser.maxNs=-1;
		parser.overwrite=overwrite;
		boolean trimByOverlap_=false, useQualityForOverlap_=false, strictOverlap_=true;
		boolean trimPairsEvenly_=false;
		boolean ordered_=false;
		int minoverlap_=-1, mininsert_=-1;
		int restrictLeft_=0, restrictRight_=0, speed_=0, qSkip_=1;
		boolean printNonZeroOnly_=true;
		boolean rename_=false, useRefNames_=false;
		boolean skipr1_=false, skipr2_=false;
		boolean ecc_=false;
		float minBaseFrequency_=0;
		float minKmerFraction_=0;
		float minCoveredFraction_=0;
		boolean setk=false;
		
		scaffoldNames.add(""); //Necessary so that the first real scaffold gets an id of 1, not zero
		scaffoldLengths.add(0);
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){
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
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(parser.parseCommon(arg, a, b)){
				//do nothing
			}else if(parser.parseCardinality(arg, a, b)){
				//do nothing
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
				useForest_=Parse.parseBoolean(b);
				if(useForest_){useTable_=useArray_=false;}
			}else if(a.equals("table")){
				useTable_=Parse.parseBoolean(b);
				if(useTable_){useForest_=useArray_=false;}
			}else if(a.equals("array")){
				useArray_=Parse.parseBoolean(b);
				if(useArray_){useTable_=useForest_=false;}
			}else if(a.equals("ways")){
				ways_=Integer.parseInt(b);
			}else if(a.equals("ordered") || a.equals("ord")){
				ordered_=Parse.parseBoolean(b);
			}else if(a.equals("skipr1")){
				skipr1_=Parse.parseBoolean(b);
			}else if(a.equals("skipr2")){
				skipr2_=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				assert(b!=null) : "\nThe k key needs an integer value greater than 0, such as k=27\n";
				k_=Integer.parseInt(b);
				setk=true;
			}else if(a.equals("mink") || a.equals("kmin")){
				mink_=Integer.parseInt(b);
				assert(mink_<0 || (mink_>0 && mink_<32)) : "kmin must be between 1 and 31; default is 4, negative numbers disable it.";
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
				qSkip_=Integer.parseInt(b);
			}else if(a.equals("speed")){
				speed_=Integer.parseInt(b);
				assert(speed_>=0 && speed_<=16) : "Speed range is 0 to 16.  Value: "+speed_;
			}else if(a.equals("skipreads")){
				skipreads_=Parse.parseKMG(b);
			}else if(a.equals("maxbadkmers") || a.equals("mbk")){
				maxBadKmers_=Integer.parseInt(b);
			}else if(a.equals("minhits") || a.equals("minkmerhits") || a.equals("mkh")){
				maxBadKmers_=Integer.parseInt(b)-1;
			}else if(a.equals("minkmerfraction") || a.equals("minfraction") || a.equals("mkf")){
				minKmerFraction_=Float.parseFloat(b);
			}else if(a.equals("mincoveredfraction") || a.equals("mincovfraction") || a.equals("mcf")){
				minCoveredFraction_=Float.parseFloat(b);
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
				assert(WAYS>1) : "WAYS=1 is for debug mode.";
			}else if(a.equals("mm") || a.equals("maskmiddle")){
				if(b==null || Tools.startsWithLetter(b)) {
					maskMiddle=Parse.parseBoolean(b);
				}else{
					midMaskLen=Integer.parseInt(b);
					maskMiddle=midMaskLen>0;
				}
			}else if(a.equals("rcomp")){
				rcomp_=Parse.parseBoolean(b);
			}else if(a.equals("forbidns") || a.equals("forbidn") || a.equals("fn")){
				forbidNs_=Parse.parseBoolean(b);
			}else if(a.equals("findbestmatch") || a.equals("fbm")){
				findBestMatch_=Parse.parseBoolean(b);
			}else if(a.equals("kfilter")){
				boolean x=Parse.parseBoolean(b);
				if(x){ktrimLeft_=ktrimRight_=ktrimN_=ksplit_=false;}
			}else if(a.equals("ksplit")){
				boolean x=Parse.parseBoolean(b);
				if(x){ksplit_=true; ktrimLeft_=ktrimRight_=ktrimN_=false;}
				else{ksplit_=false;}
			}else if(a.equals("ktrim")){
				if(b==null){b="";}
				if(b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr") || b.equalsIgnoreCase("tips")){
					ktrimLeft_=ktrimRight_=true;ktrimN_=ksplit_=false;
				}else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){
					ktrimLeft_=true;ktrimRight_=false;ktrimN_=ksplit_=false;
				}else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){
					ktrimLeft_=false;ktrimRight_=true;ktrimN_=ksplit_=false;
				}else if(b.equalsIgnoreCase("n")){
					ktrimLeft_=ktrimRight_=ksplit_=false;
					ktrimN_=true;
				}else if(b.length()==1 && !b.equalsIgnoreCase("t") && !b.equalsIgnoreCase("f")){
					ktrimLeft_=ktrimRight_=ksplit_=false;
					ktrimN_=true;
					TRIM_SYMBOL_=(byte)b.charAt(0);
				}else{
					assert(b!=null && (b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false"))) :
						"\nInvalid setting for ktrim - values must be f (false), l (left), r (right), rl (tips), or n.";
					ktrimRight_=ktrimLeft_=false;
				}
			}else if(a.equals("trimtips") || a.equals("ktrimtips")){
				if(b==null){b="";}
				if(b.length()>0){
					ktrimLeft_=ktrimRight_=true;ktrimN_=ksplit_=false;
					restrictLeft_=restrictRight_=Integer.parseInt(b);
				}
			}else if(a.equals("kmask") || a.equals("mask")){
				if("lc".equalsIgnoreCase(b) || "lowercase".equalsIgnoreCase(b)){
					kmaskLowercase_=true;
					ktrimN_=true;
					ktrimLeft_=ktrimRight_=ksplit_=false;
				}else{
					if(Parse.parseBoolean(b)){b="N";}
					if(b!=null && b.length()==1){
						ktrimLeft_=false;ktrimRight_=false;ktrimN_=true;
						TRIM_SYMBOL_=(byte)b.charAt(0);
					}else{
						boolean x=Parse.parseBoolean(b);
//						assert(!x) : "\nInvalid setting for kmask - values must be f (false), t (true), or a single character for replacement.";
						ktrimN_=x;
					}
				}
			}else if(a.equals("kmaskfullycovered") || a.equals("maskfullycovered") || a.equals("mfc")){
				kmaskFullyCovered_=Parse.parseBoolean(b);
			}else if(a.equals("ktrimright")){
				ktrimRight_=Parse.parseBoolean(b);
				ktrimLeft_=ktrimN_=!(ktrimRight_);
			}else if(a.equals("ktrimleft")){
				ktrimLeft_=Parse.parseBoolean(b);
				ktrimRight_=ktrimN_=!(ktrimLeft_);
			}else if(a.equals("ktrimn")){
				ktrimN_=Parse.parseBoolean(b);
				ktrimLeft_=ktrimRight_=!(ktrimN_);
			}else if(a.equals("ktrimexclusive")){
				ktrimExclusive_=Parse.parseBoolean(b);
			}else if(a.equals("tbo") || a.equals("trimbyoverlap")){
				trimByOverlap_=Parse.parseBoolean(b);
			}else if(a.equals("strictoverlap")){
				strictOverlap_=Parse.parseBoolean(b);
			}else if(a.equals("usequality")){
				useQualityForOverlap_=Parse.parseBoolean(b);
			}else if(a.equals("tpe") || a.equals("tbe") || a.equals("trimpairsevenly")){
				trimPairsEvenly_=Parse.parseBoolean(b);
			}else if(a.equals("ottm") || a.equals("outputtrimmedtomatch")){
				addTrimmedToBad_=Parse.parseBoolean(b);
			}else if(a.equals("minoverlap")){
				minoverlap_=Integer.parseInt(b);
			}else if(a.equals("mininsert")){
				mininsert_=Integer.parseInt(b);
			}else if(a.equals("prealloc") || a.equals("preallocate")){
				if(b==null || b.length()<1 || Character.isLetter(b.charAt(0))){
					prealloc_=Parse.parseBoolean(b);
				}else{
					preallocFraction=Tools.max(0, Double.parseDouble(b));
					prealloc_=(preallocFraction>0);
				}
			}else if(a.equals("restrictleft")){
				restrictLeft_=Integer.parseInt(b);
			}else if(a.equals("restrictright")){
				restrictRight_=Integer.parseInt(b);
			}else if(a.equals("statscolumns") || a.equals("columns") || a.equals("cols")){
				STATS_COLUMNS=Integer.parseInt(b);
				assert(STATS_COLUMNS==3 || STATS_COLUMNS==5) : "statscolumns bust be either 3 or 5. Invalid value: "+STATS_COLUMNS;
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				printNonZeroOnly_=Parse.parseBoolean(b);
			}else if(a.equals("rename")){
				rename_=Parse.parseBoolean(b);
			}else if(a.equals("refnames") || a.equals("userefnames")){
				useRefNames_=Parse.parseBoolean(b);
			}else if(a.equals("initialsize")){
				initialSize=Parse.parseIntKMG(b);
			}else if(a.equals("dump")){
				dump=b;
			}else if(a.equals("minentropy") || a.equals("entropy") || a.equals("entropyfilter")){
				entropyCutoff=Float.parseFloat(b);
			}else if(a.equals("verifyentropy")){
				EntropyTracker.verify=Parse.parseBoolean(b);
//			}else if(a.equals("entropytracker") || a.equals("usetracker")){
//				useEntropyTracker=Parse.parseBoolean(b);
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
				minBaseFrequency_=Float.parseFloat(b);
			}else if(a.equals("ecco") || a.equals("ecc")){
				ecc_=Parse.parseBoolean(b);
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
				histogramsBeforeProcessing_=Parse.parseBoolean(b);
			}else if(a.equals("histogramsafter") || a.equals("histafter")){
				histogramsBeforeProcessing_=!Parse.parseBoolean(b);
			}
			
			else if(a.equals("trimfailures") || a.equals("trimfailuresto1bp")){
				trimFailuresTo1bp_=Parse.parseBoolean(b);
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
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		if(WAYS!=7){
			System.err.println("WARNING! WAYS="+WAYS+".  Probably debug mode.");
		}
		
		if(khistIn!=null || khistOut!=null){
			if(khistIn!=null){parser.loglog=true;}
			if(khistOut!=null){parser.loglogOut=true;}
			parser.loglogbuckets=Tools.max(parser.loglogbuckets, 128000);
			CardinalityTracker.trackCounts=true;
		}
		
		if(ordered_ && !silent){
			outstream.println("Set ORDERED to "+ordered_);
			if(jsonStats!=null){jsonStats.add("ordered", ordered_);}
		}
		if(silent || json){
			AbstractKmerTableSet.DISPLAY_PROGRESS=false; //TODO: Test to make sure this occurs for silent mode 
		}
		
		if(hammingDistance2==-1){hammingDistance2=hammingDistance;}
		if(qHammingDistance2==-1){qHammingDistance2=qHammingDistance;}
		if(editDistance2==-1){editDistance2=editDistance;}
		minBaseFrequency=minBaseFrequency_;
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			samplerate=parser.samplerate;
			sampleseed=parser.sampleseed;
			recalibrateQuality=parser.recalibrateQuality;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
//			testsize=parser.testsize;
//			trimBadSequence=parser.trimBadSequence;
//			breakLength=parser.breakLength;
			
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
//			untrim=parser.untrim;
//			minTrimLength=(parser.minTrimLength>=0 ? parser.minTrimLength : minTrimLength);
//			requireBothBad=parser.requireBothBad;
			removePairsIfEitherBad=(!parser.requireBothBad) && (!trimFailuresTo1bp_);
			tossJunk=parser.tossJunk;

			minGC=parser.minGC;
			maxGC=parser.maxGC;
			filterGC=(minGC>0 || maxGC<1);
			usePairGC=parser.usePairGC;

			loglogIn=(parser.loglog ? CardinalityTracker.makeTracker(parser) : null);
			loglogOut=(parser.loglogOut ? CardinalityTracker.makeTracker(parser) : null);
			
			THREADS=Shared.threads();
			silent=Parser.silent;
			if(silent){
				DISPLAY_PROGRESS=false;
				showSpeed=false;
			}
		}

		ref=modifyRefPath(ref, refNames);
		altref=modifyRefPath(altref, altRefNames);
		if(literal!=null){
			refNames.add("literal");
			if(!altRefNames.isEmpty()){altRefNames.add("literal");}
		}
		refScafCounts=new int[refNames.size()];
		
		if(minoverlap_>=0){
			minOverlap=Tools.max(minoverlap_, 1);
			minOverlap0=Tools.min(minOverlap0, minOverlap);
		}
		
		if(mininsert_>=0){
			minInsert=Tools.max(mininsert_, 1);
			minInsert0=Tools.min(minInsert0, minInsert);
		}
		
		/* Set final variables; post-process and validate argument combinations */
		
		useForest=useForest_;
		useTable=useTable_;
		useArray=useArray_;
		hammingDistance=Tools.max(editDistance, hammingDistance);
		hammingDistance2=Tools.max(editDistance2, hammingDistance2);
		minSkip=Tools.max(1, Tools.min(minSkip, maxSkip));
		maxSkip=Tools.max(minSkip, maxSkip);
		addTrimmedToBad=addTrimmedToBad_;
		forbidNs=(forbidNs_ || hammingDistance<1);
		trimSymbol=TRIM_SYMBOL_;
		kmaskLowercase=kmaskLowercase_;
		kmaskFullyCovered=kmaskFullyCovered_;
		skipreads=skipreads_;
		trimByOverlap=trimByOverlap_;
		useQualityForOverlap=useQualityForOverlap_;
		strictOverlap=strictOverlap_;
		trimPairsEvenly=trimPairsEvenly_;
		ordered=ordered_;
		restrictLeft=Tools.max(restrictLeft_, 0);
		restrictRight=Tools.max(restrictRight_, 0);
		printNonZeroOnly=printNonZeroOnly_;
		rename=rename_;
		findBestMatch=(rename || findBestMatch_);
		useRefNames=useRefNames_;
		speed=speed_;
		qSkip=qSkip_;
		noAccel=(speed<1 && qSkip<2);
		accel=!noAccel;
		skipR1=skipr1_;
		skipR2=skipr2_;
		ecc=ecc_;
		locationFilter=(xMinLoc>0 || yMinLoc>0 || xMaxLoc>-1 || yMaxLoc>-1);
		trimFailuresTo1bp=trimFailuresTo1bp_;
		
		amino=Shared.AMINO_IN;
		rcomp=rcomp_ && !amino;
		
		//Set K
		maxSupportedK=(amino ? 12 : 31);
		if(!setk){k_=(amino ? 11 : 27);}
		kbig_=(k_>maxSupportedK ? k_ : -1);
		k_=Tools.min(k_, maxSupportedK);
		
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
		
		histogramsBeforeProcessing=histogramsBeforeProcessing_;
		MAKE_QUALITY_HISTOGRAM=ReadStats.COLLECT_QUALITY_STATS;
		MAKE_QUALITY_ACCURACY=ReadStats.COLLECT_QUALITY_ACCURACY;
		MAKE_MATCH_HISTOGRAM=ReadStats.COLLECT_MATCH_STATS;
		MAKE_BASE_HISTOGRAM=ReadStats.COLLECT_BASE_STATS;
		MAKE_EHIST=ReadStats.COLLECT_ERROR_STATS;
		MAKE_INDELHIST=ReadStats.COLLECT_INDEL_STATS;
		MAKE_LHIST=ReadStats.COLLECT_LENGTH_STATS;
		MAKE_GCHIST=ReadStats.COLLECT_GC_STATS;
		MAKE_ENTROPYHIST=ReadStats.COLLECT_ENTROPY_STATS;
		MAKE_IDHIST=ReadStats.COLLECT_IDENTITY_STATS;
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
				final double mem2=(prealloc_ ? preallocFraction : 1)*tableMemory;
				initialSize=(prealloc_ || memOverWays<initialSizeDefault ? (int)Tools.min(2142000000, (long)(mem2/(12*WAYS))) : initialSizeDefault);
				if(initialSize!=initialSizeDefault && !silent){
					outstream.println("Initial size set to "+initialSize);
				}
			}
		}
		
		if(ktrimLeft_ || ktrimRight_ || ktrimN_ || ksplit_){
			if(kbig_>k_){
				outstream.println("***********************   WARNING   ***********************");
				outstream.println("WARNING: When kmer-trimming or masking, the maximum value of K is "+k_+".");
				outstream.println("K has been reduced from "+kbig_+" to "+k_+".");
				outstream.println("***********************************************************");
				kbig_=k_;
			}
		}
		
		if((speed>0 || qSkip>1) && kbig_>k_){
			outstream.println("***********************   WARNING   ***********************");
			outstream.println("WARNING: When speed>0 or qskip>1, the maximum value of K is "+k_+".");
			outstream.println("K has been reduced from "+kbig_+" to "+k_+".");
			outstream.println("***********************************************************");
			kbig_=k_;
		}
		
		if((speed>0 && qSkip>1) || (qSkip>1 && maxSkip>1) || (speed>0 && maxSkip>1)){
			outstream.println("WARNING: It is not recommended to use more than one of 'qskip', 'speed', and 'rskip/maxskip' together.");
			outstream.println("qskip="+qSkip+", speed="+speed+", maxskip="+maxSkip);
		}
		
		k=k_;
		k2=k-1;
		kbig=kbig_;
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
		mink=Tools.min((mink_<1 ? 6 : mink_), k);
		maxBadKmers0=maxBadKmers_;
		
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
		
		minKmerFraction=Tools.max(minKmerFraction_, 0);
		assert(minKmerFraction<=1) : "minKmerFraction must range from 0 to 1; value="+minKmerFraction;
		
		minCoveredFraction=Tools.max(minCoveredFraction_, 0);
		assert(minCoveredFraction<=1) : "minCoveredFraction must range from 0 to 1; value="+minCoveredFraction;
		
		if(mink_>0 && mink_<k){useShortKmers=true;}
		if(useShortKmers){
			if(maskMiddle){
				outstream.println("maskMiddle was disabled because useShortKmers=true");
				maskMiddle=false;
				midMaskLen=0;
			}
		}
		
		ktrimRight=ktrimRight_;
		ktrimLeft=ktrimLeft_;
		ktrimN=ktrimN_;
		ksplit=ksplit_;
		ktrimExclusive=ktrimExclusive_;
		kfilter=(ref!=null || literal!=null) && !(ktrimRight || ktrimLeft || ktrimN || ksplit);
		assert(findBestMatch==false || kfilter==false || kbig<=k) : "K must be less than 32 in 'findBestMatch' mode";
		
		assert(!useShortKmers || ktrimRight || ktrimLeft || ktrimN || ksplit) : "\nSetting mink or useShortKmers also requires setting a ktrim mode, such as 'r', 'l', or 'n'\n";
		
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
//		assert(false) : Long.toBinaryString(middleMask)+", "+Long.toBinaryString((1L<<(2*k))-1);
//		middleMask=maskMiddle ? ~(symbolMask<<(bitsPerBase*(k/2))) : -1L;
		
		hitCounts=(outduk==null ? null : new long[HITCOUNT_LEN+1]);
		
		
		/* Adjust I/O settings and filenames */
		
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

		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		qfin1=Tools.fixExtension(qfin1);
		qfin2=Tools.fixExtension(qfin2);
		ref=Tools.fixExtension(ref);
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, qfout1, qfout2, outb1, outb2, outsingle, outstats, 
				outrpkm, outduk, outrqc, outrefstats, polymerStatsFile, khistIn, khistOut, alignOut)){
			throw new RuntimeException("\nCan't write to some output files; overwrite="+overwrite+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2, qfin1, qfin2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(!Tools.testInputFiles(true, true, ref)){
			throw new RuntimeException("\nCan't read to some reference files.\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, qfin1, qfin2, qfout1, qfout2,
				out1, out2, outb1, outb2, outsingle, outstats, outrpkm, outduk, outrqc, outrefstats, polymerStatsFile, khistIn, khistOut, alignOut)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		assert(THREADS>0) : "THREADS must be greater than 0.";

		assert(in1==null || in1.toLowerCase().startsWith("stdin") || in1.toLowerCase().startsWith("standardin") || new File(in1).exists()) : "Can't find "+in1;
		assert(in2==null || in2.toLowerCase().startsWith("stdin") || in2.toLowerCase().startsWith("standardin") || new File(in2).exists()) : "Can't find "+in2;
		
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
		
		final int defaultFormat=(ffin1==null ? FileFormat.FASTQ : ffin1.format());
		ffout1=FileFormat.testOutput(out1, defaultFormat, null, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, defaultFormat, null, true, overwrite, append, ordered);
		ffoutb1=FileFormat.testOutput(outb1, defaultFormat, null, true, overwrite, append, ordered);
		ffoutb2=FileFormat.testOutput(outb2, defaultFormat, null, true, overwrite, append, ordered);
		ffouts=FileFormat.testOutput(outsingle, defaultFormat, null, true, overwrite, append, ordered);

		parser.validateStdio(ffin1, ffout1, ffoutb1, ffouts);
		
		makeReadStats=ReadStats.collectingStats();
		
		//This block just causes problems when new features are added, so it's disabled.
		if(!((ref!=null || literal!=null) || qtrimLeft || qtrimRight || minAvgQuality>0 || minBaseQuality>0 || maxNs>=0 || trimByOverlap ||
				makeReadStats || entropyMask || entropyTrim>0 || entropyMark || entropyCutoff>=0 || filterVars ||
				forceTrimLeft>0 || forceTrimRight>0 || forceTrimRight2>0 || forceTrimModulo>0 || minBaseFrequency>0 || recalibrateQuality || 
				trimPolyA>0 || trimPolyGLeft>0 || trimPolyGRight>0 || filterPolyG>0 || trimPolyCLeft>0 || trimPolyCRight>0 || filterPolyC>0)){
//			outstream.println("NOTE: No reference files specified, no trimming mode, no min avg quality, no histograms - read sequences will not be changed.");
		}
		
		if(recalibrateQuality || true){
			SamLine.SET_FROM_OK=true;//TODO:  Should ban other operations
		}
				
		if(ref!=null){
			for(String s0 : ref){
				assert(s0!=null) : "Specified a null reference.";
				String s=s0.toLowerCase();
				assert(s==null || s.startsWith("stdin") || s.startsWith("standardin") || new File(s0).exists()) : "Can't find "+s0;
			}
		}
		
		align=(align || alignRef!=null);
		if(align) {
			sidechannel=new SideChannel3(alignRef, alignOut, null, alignK1, alignK2, 
					alignMinid1, alignMinid2, alignMM1, alignMM2, overwrite, ordered);
		}else {
			sidechannel=null;
		}
		
		//Initialize tables
		final int tableType=(useForest ? AbstractKmerTable.FOREST1D : useTable ? AbstractKmerTable.TABLE : useArray ? AbstractKmerTable.ARRAY1D : 0);
		ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, 12, prealloc_, (prealloc_ ? preallocFraction : 0.9));
		int[] schedule=scheduleMaker.makeSchedule();
		outstream.print("Allocating kmer table: \t");
		keySets=AbstractKmerTable.preallocate(WAYS, tableType, schedule, -1L);
		
		//Initialize entropy
		calcEntropy=(entropyCutoff>=0 || entropyMark);
		if(calcEntropy){
			assert(EntropyTracker.defaultWindowBases>0 && (entropyMark || (entropyCutoff>=0 && entropyCutoff<=1)));
		}
		assert(calcEntropy || (!entropyMask && !entropyMark && entropyTrim==0)) : "Entropy masking/trimming operations require the entropy flag to be set.";
		
		if(polymerStatsFile!=null || (polymerChar1>=0 && polymerChar2>=0)){
			countPolymers=true;
		}
		
		//Initialize polymer-tracking
		if(countPolymers){
//			assert(polymerChar1>=0 && AminoAcid.baseToNumberACGTN[polymerChar1]>=0);
//			assert(polymerChar2>=0 && AminoAcid.baseToNumberACGTN[polymerChar2]>=0);
			pTracker=new PolymerTracker();
		}else{
			pTracker=null;
		}
	}
	
	String[] modifyRefPath(String[] array, ArrayList<String> list){
		if(array==null){return array;}
		for(String s : array){
			String s2=modifyRefPath(s);
			list.add(s2);
		}
		return list.toArray(new String[0]);
	}
	
	public static String modifyRefPath(String s){
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
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public void process(){
		
		if(samref!=null){
			scafMap=ScafMap.loadReference(samref, true);
		}
		
		if(varFile!=null || vcfFile!=null || filterVars){
			if(scafMap==null){scafMap=ScafMap.loadSamHeader(in1);}
			assert(scafMap!=null && scafMap.size()>0) : "No scaffold names were loaded.";
			if(varFile!=null){
				outstream.println("Loading variants.");
				varMap=VcfLoader.loadVarFile(varFile, scafMap);
			}else if(vcfFile!=null){
				outstream.println("Loading variants.");
				varMap=VcfLoader.loadVcfFile(vcfFile, scafMap, false, false);
			}
//			if(varMap!=null && varMap.size()>0){
//				varKeySet=new HashSet<VarKey>();
//				for(Var v : varMap.toArray(false)){
//					if(v.type==Var.INS || v.type==Var.DEL){
//						varKeySet.add(VarKey.toVarKey(v));
//					}
//				}
//			}
			fixVariants=(makeReadStats && varMap!=null && varMap.size()>0 && scafMap!=null && scafMap.size()>0);
		}
		
		if(recalibrateQuality){
			if(samFile!=null){
				CalcTrueQuality.main2(new String[] {"in="+samFile, "showstats=f"});
			}
			CalcTrueQuality.initializeMatrices();
		}
		
		/* Check for output file collisions */
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outb1, outb2, outstats, outrpkm, outduk, outrqc, outrefstats)){
			throw new RuntimeException("One or more output files were duplicate or could not be written to.  Check the names or set the 'overwrite=true' flag.");
		}
		
		/* Start overall timer */
		Timer t=new Timer();
		
//		boolean dq0=FASTQ.DETECT_QUALITY;
//		boolean ti0=FASTQ.TEST_INTERLEAVED;
//		int rbl0=Shared.bufferLen();;
//		FASTQ.DETECT_QUALITY=false;
//		FASTQ.TEST_INTERLEAVED=false;
//		Shared.setBufferLen(16;
		
		process2(t.time1);
		
//		FASTQ.DETECT_QUALITY=dq0;
//		FASTQ.TEST_INTERLEAVED=ti0;
//		Shared.setBufferLen(rbl0;
		
		/* Stop timer and calculate speed statistics */
		t.stop();
		lastReadsOut=readsOut;
		
		if(showSpeed && !json){
			outstream.println();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		}
		
		if(outstream!=System.err && outstream!=System.out){outstream.close();}
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	public void process2(long startTime){
		
		/* Start phase timer */
		Timer t=new Timer();

		if(DISPLAY_PROGRESS && !json){
			outstream.println("Initial:");
			Shared.printMemory(outstream);
			outstream.println();
		}
		
		/* Fill tables with reference kmers */
		if((ref!=null && ref.length>0) || (literal!=null && literal.length>0)){
			final boolean oldTI=FASTQ.TEST_INTERLEAVED; //TODO: This needs to be changed to a non-static field, or somehow 'read mode' and 'ref mode' need to be distinguished.
			final boolean oldFI=FASTQ.FORCE_INTERLEAVED;
			final boolean oldPC=FASTQ.PARSE_CUSTOM;
			final boolean oldSplit=FastaReadInputStream.SPLIT_READS;
			final int oldML=FastaReadInputStream.MIN_READ_LEN;
			
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
			FASTQ.PARSE_CUSTOM=false;
			FastaReadInputStream.SPLIT_READS=false;
			FastaReadInputStream.MIN_READ_LEN=1;
			
			
			storedKmers=spawnLoadThreads();
			if(storedKmers<1 && altRefNames.size()>0){
				outstream.println("Ref had no kmers; using alt ref.");
				ref=altref;
				refNames=altRefNames;
				refScafCounts=new int[refNames.size()];
				scaffoldNames.clear();
				scaffoldLengths.clear();
				storedKmers=spawnLoadThreads();
			}
			
			FASTQ.TEST_INTERLEAVED=oldTI;
			FASTQ.FORCE_INTERLEAVED=oldFI;
			FASTQ.PARSE_CUSTOM=oldPC;
			FastaReadInputStream.SPLIT_READS=oldSplit;
			FastaReadInputStream.MIN_READ_LEN=oldML;
			
			if(useRefNames){toRefNames();}
			t.stop();
		}
		
		{
			long ram=freeMemory();
			ALLOW_LOCAL_ARRAYS=(scaffoldNames!=null && Tools.max(THREADS, 1)*3*8*scaffoldNames.size()<ram*5);
		}
		
		/* Dump kmers to text */
		if(dump!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(dump, overwrite, false, true);
			bsw.start();
			for(AbstractKmerTable set : keySets){
				set.dumpKmersAsBytes(bsw, k, 0, Integer.MAX_VALUE, null);
			}
			bsw.poisonAndWait();
		}
		
		if(storedKmers<1 && (ktrimRight || ktrimLeft || ktrimN || ksplit)){
			outstream.println("******  WARNING! A KMER OPERATION WAS CHOSEN BUT NO KMERS WERE LOADED.  ******");
			if(ref==null && literal==null){
				outstream.println("******  YOU NEED TO SPECIFY A REFERENCE FILE OR LITERAL SEQUENCE.       ******\n");
			}else{
				outstream.println("******  PLEASE ENSURE K IS LESS THAN OR EQUAL TO REF SEQUENCE LENGTHS.  ******\n");
			}
			if(ktrimRight && trimByOverlap){
				ktrimRight=false;
			}else{
				assert(false) : "You can bypass this assertion with the -da flag.";
			}
		}
		
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=THREADS<4;
		
		/* Do kmer matching of input reads */
		spawnProcessThreads(t);
		
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		/* Write legacy duk statistics (which requires tables) */
		writeDuk(System.nanoTime()-startTime);
		
		/* Unload kmers to save memory */
		if(RELEASE_TABLES){unloadKmers();}
		
		/* Write statistics to files */
		writeStats();
		writeRPKM();
		writeRefStats();
		writeRqc();
		if(pTracker!=null && polymerStatsFile!=null){
			ReadWrite.writeString(pTracker.toHistogramCumulative(), polymerStatsFile);
		}
		
		/* Unload sequence data to save memory */
		if(RELEASE_TABLES){unloadScaffolds();}
		
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
	
	private String toJson(long startTime){

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
	
	public static String toPercent(long numerator, long denominator){
		if(denominator<1){return "0.00%";}
		return Tools.format("%.2f%%",numerator*100.0/denominator);
	}
	
	private static String padRight(String s, int minLen){
		while(s.length()<minLen){s=s+" ";}
		return s;
	}
	
	/**
	 * Clear stored kmers.
	 */
	public void unloadKmers(){
		if(keySets!=null){
			for(int i=0; i<keySets.length; i++){keySets[i]=null;}
		}
	}
	
	/**
	 * Clear stored sequence data.
	 */
	public void unloadScaffolds(){
		if(scaffoldNames!=null && !scaffoldNames.isEmpty()){
			scaffoldNames.clear();
			scaffoldNames.trimToSize();
		}
		scaffoldReadCounts=null;
		scaffoldBaseCounts=null;
		hitCounts=null;
		scaffoldLengths=null;
	}
	
	/**
	 * Write statistics about how many reads matched each reference scaffold.
	 */
	private void writeStats(){
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
	private void writeRPKM(){
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
	 * Write statistics on a per-reference basis.
	 */
	private void writeRefStats(){
		if(outrefstats==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outrefstats, overwrite, false, false);
		tsw.start();
		
		/* Count mapped reads */
		long mapped=0;
		for(int i=0; i<scaffoldReadCounts.length(); i++){
			mapped+=scaffoldReadCounts.get(i);
		}
		
		final int numRefs=refNames.size();
		long[] refReadCounts=new long[numRefs];
		long[] refBaseCounts=new long[numRefs];
		long[] refLengths=new long[numRefs];
		
		for(int r=0, s=1; r<numRefs; r++){
			final int lim=s+refScafCounts[r];
			while(s<lim){
				refReadCounts[r]+=scaffoldReadCounts.get(s);
				refBaseCounts[r]+=scaffoldBaseCounts.get(s);
				refLengths[r]+=scaffoldLengths.get(s);
				s++;
			}
		}
		
		/* Print header */
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		tsw.print(Tools.format("#Reads\t%d\n",readsIn));
		tsw.print(Tools.format("#Mapped\t%d\n",mapped));
		tsw.print(Tools.format("#References\t%d\n",Tools.max(0, refNames.size())));
		tsw.print("#Name\tLength\tScaffolds\tBases\tCoverage\tReads\tRPKM\n");
		
		final float mult=1000000000f/Tools.max(1, mapped);
		
		/* Print data */
		for(int i=0; i<refNames.size(); i++){
			final long reads=refReadCounts[i];
			final long bases=refBaseCounts[i];
			final long len=refLengths[i];
			final int scafs=refScafCounts[i];
			final String name=ReadWrite.stripToCore(refNames.get(i));
			final double invlen=1.0/Tools.max(1, len);
			final double mult2=mult*invlen;
			if(reads>0 || !printNonZeroOnly){
				tsw.print(Tools.format("%s\t%d\t%d\t%d\t%.4f\t%d\t%.4f\n",name,len,scafs,bases,bases*invlen,reads,reads*mult2));
			}
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Write processing statistics in DUK's format.
	 * @param time Elapsed time, nanoseconds
	 */
	private void writeDuk(long time){
		if(outduk==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outduk, overwrite, false, false);
		tsw.start();
		tsw.println(dukString(time));
		tsw.poisonAndWait();
	}
	
	/**
	 * Write RQCFilter stats.
	 * @param time Elapsed time, nanoseconds
	 */
	private void writeRqc(){
		if(outrqc==null){return;}
		addToRqcMap();
		if(outrqc.endsWith("hashmap")){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outrqc, overwrite, false, false);
		tsw.start();
		tsw.println(rqcString());
		tsw.poisonAndWait();
	}
	
	public static String rqcString(){
		if(RQC_MAP==null){return null;}
		StringBuilder sb=new StringBuilder();
		
		String[] keys=new String[] {"inputReads", "inputBases", "qtrimmedReads", "qtrimmedBases", "qfilteredReads", "qfilteredBases",
				"ktrimmedReads", "ktrimmedBases", "kfilteredReads", "kfilteredBases", "outputReads", "outputBases"};
		
		for(String key : keys){
			Object value=RQC_MAP.get(key);
			if(value!=null){
				sb.append(key+"="+value+"\n");
			}
		}
		
		return sb.toString();
	}
	
	private void addToRqcMap(){
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
	
	public static void putRqc(String key, Long value, boolean evict, boolean add){
		if(RQC_MAP==null){RQC_MAP=new HashMap<String,Long>();}
		Long old=RQC_MAP.get(key);
		if(evict || old==null){RQC_MAP.put(key, value);}
		else if(add){RQC_MAP.put(key, value+old);}
	}
	
	/**
	 * Helper method; formats statistics to be duk-compatible
	 * @param time Elapsed time, nanoseconds
	 * @return duk output string
	 */
	private String dukString(long time){
		StringBuilder sb=new StringBuilder();
		sb.append("##INPUT PARAMETERS##\n");
		sb.append("#Reference file:	"+(ref==null || ref.length<1 ? null : ref.length==1 ? ref[0] : Arrays.toString(ref))+"\n");
		sb.append("#Query file:	"+in1+(in2==null ? "" : ","+in2)+"\n");
		sb.append("#Not matched reads file:	"+out1+(out2==null ? "" : ","+out2)+"\n");
		sb.append("#Matched reads file:	"+outb1+(outb2==null ? "" : ","+outb2)+"\n");
		sb.append("#Output file (duk):	"+outduk+"\n");
		sb.append("#Output file (stats):	"+outstats+"\n");
		sb.append("#Mer size:	"+k+"\n");
		long size=0;
		for(AbstractKmerTable x : keySets){size+=x.size();}
		sb.append("#Avg step size:	"+Tools.format("%.1f", refKmers/(double)(Tools.max(1, size)))+"\n");
		sb.append("#Cut off:	"+maxBadKmers0+"\n");
		sb.append("#Mask middle:	"+midMaskLen+"\n");
		sb.append("#Quality trim:	"+((qtrimLeft || qtrimRight) ? trimq : "false")+"\n");
		sb.append("\n");
		
		sb.append("##REFERENCE STAT##\n");
		sb.append("#Total Reads:	"+refReads+"\n");
		sb.append("#Total Bases:	"+refBases+"\n");
		sb.append("#Total kmers:	"+refKmers+"\n");
		sb.append("#Total stored kmers:	"+size+"\n");
		sb.append("\n");

		sb.append("## ELAPSED TIME##\n");
		sb.append("# Time:	"+Tools.format("%.2f", time/1000000000.0)+" seconds\n");
		sb.append("\n");

		sb.append("##QUERY FILE STAT##\n");
		sb.append("# Total number of reads:	"+readsIn+"\n");
		sb.append("# Total number of matched reads:	"+readsKFiltered+"\n");
		sb.append("# Match ratio:	"+Tools.format("%.6f", readsKFiltered*1.0/readsIn)+"\n");
		sb.append("\n");

		sb.append("##P-VALUE##\n");
		sb.append("#Avg number of Kmer for each read:	"+((basesIn/(Tools.max(readsIn, 1)))-k)+"\n");
//		sb.append("# P value for the given threshold 1 is 4.05231e-14\n"); //duk prints a P value; not sure what it means
		sb.append("\n");

		sb.append("## Histogram of kmer occurance for reads with at least one occurance ##\n");
		sb.append("#NumOcc\tNumReads\tPercentage\n");
		
		long sum=shared.Vector.sum(hitCounts);
		double mult=100.0/(sum<1 ? 1 : sum);
		for(int i=0; i<hitCounts.length; i++){
			long x=hitCounts[i];
			if(x>0){
				sb.append(i).append('\t').append(x).append('\t').append(Tools.format("%.4f",(x*mult))).append('\n');
			}
		}
		
		return sb.toString();
	}
	
	/**
	 * Fills the scaffold names array with reference names.
	 */
	private void toRefNames(){
		final int numRefs=refNames.size();
		for(int r=0, s=1; r<numRefs; r++){
			final int scafs=refScafCounts[r];
			final int lim=s+scafs;
			final String name=ReadWrite.stripToCore(refNames.get(r));
//			outstream.println("r="+r+", s="+s+", scafs="+scafs+", lim="+lim+", name="+name);
			while(s<lim){
//				outstream.println(r+", "+s+". Setting "+scaffoldNames.get(s)+" -> "+name);
				scaffoldNames.set(s, name);
				s++;
			}
		}
	}
	
	public double getPolymerRatio(){
		return pTracker.calcRatioCumulative(polymerChar1, polymerChar2, polymerLength);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	

	/**
	 * Fills tables with kmers from references, using multiple LoadThread.
	 * @return Number of kmers stored.
	 */
	private long spawnLoadThreads(){
		Timer t=new Timer();
		if((ref==null || ref.length<1) && (literal==null || literal.length<1)){return 0;}
		long added=0;
		
		/* Create load threads */
		LoadThread[] loaders=new LoadThread[WAYS];
		for(int i=0; i<loaders.length; i++){
			loaders[i]=new LoadThread(i);
			loaders[i].start();
		}
		
		/* For each reference file... */
		int refNum=0;
		if(ref!=null){
			for(String refname : ref){

				/* Start an input stream */
				FileFormat ff=FileFormat.testInput(refname, FileFormat.FASTA, null, false, true);
				ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1L, false, ff, null, null, null, Shared.USE_MPI, true);
				cris.start(); //4567
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				/* Iterate through read lists from the input stream */
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					{
						/* Assign a unique ID number to each scaffold */
						ArrayList<Read> reads2=new ArrayList<Read>(reads);
						for(Read r1 : reads2){
							final Read r2=r1.mate;
							final Integer id=scaffoldNames.size();
							refScafCounts[refNum]++;
							scaffoldNames.add(r1.id==null ? id.toString() : r1.id);
							int len=r1.length();
							r1.obj=id;
							if(r2!=null){
								r2.obj=id;
								len+=r2.length();
							}
							scaffoldLengths.add(len);
						}
						
						if(REPLICATE_AMBIGUOUS){
							reads2=Tools.replicateAmbiguous(reads2, Tools.min(k, mink));
						}

						/* Send a pointer to the read list to each LoadThread */
						for(LoadThread lt : loaders){
							boolean b=true;
							while(b){
								try {
									lt.queue.put(reads2);
									b=false;
								} catch (InterruptedException e) {
									//TODO:  This will hang due to still-running threads.
									throw new RuntimeException(e);
								}
							}
						}
					}

					/* Dispose of the old list and fetch a new one */
					cris.returnList(ln);
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				/* Cleanup */
				cris.returnList(ln);
				errorState|=ReadWrite.closeStream(cris);
				refNum++;
			}
		}

		/* If there are literal sequences to use as references */
		if(literal!=null){
			ArrayList<Read> list=new ArrayList<Read>(literal.length);
			if(verbose){outstream.println("Adding literals "+Arrays.toString(literal));}

			/* Assign a unique ID number to each literal sequence */
			for(int i=0; i<literal.length; i++){
				final Integer id=scaffoldNames.size();
				final Read r=new Read(literal[i].getBytes(), null, id);
				refScafCounts[refNum]++;
				scaffoldNames.add(id.toString());
				scaffoldLengths.add(r.length());
				r.obj=id;
				list.add(r);
			}
			
			if(REPLICATE_AMBIGUOUS){
				list=Tools.replicateAmbiguous(list, Tools.min(k, mink));
			}

			/* Send a pointer to the read list to each LoadThread */
			for(LoadThread lt : loaders){
				boolean b=true;
				while(b){
					try {
						lt.queue.put(list);
						b=false;
					} catch (InterruptedException e) {
						//TODO:  This will hang due to still-running threads.
						throw new RuntimeException(e);
					}
				}
			}
		}
		
		/* Signal loaders to terminate */
		for(LoadThread lt : loaders){
			boolean b=true;
			while(b){
				try {
					lt.queue.put(POISON);
					b=false;
				} catch (InterruptedException e) {
					//TODO:  This will hang due to still-running threads.
					throw new RuntimeException(e);
				}
			}
		}
		
		/* Wait for loaders to die, and gather statistics */
		boolean success=true;
		for(LoadThread lt : loaders){
			while(lt.getState()!=Thread.State.TERMINATED){
				try {
					lt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			added+=lt.addedT;
			refKmers+=lt.refKmersT;
			refBases+=lt.refBasesT;
			refReads+=lt.refReadsT;
//			modsum+=lt.modsumT;
			success&=lt.success;
		}
		if(!success){KillSwitch.kill("Failed loading ref kmers; aborting.");}
		
		//Correct statistics for number of threads, since each thread processes all reference data
		refKmers/=WAYS;
		refBases/=WAYS;
		refReads/=WAYS;
		
		scaffoldReadCounts=new AtomicLongArray(scaffoldNames.size());
		scaffoldBaseCounts=new AtomicLongArray(scaffoldNames.size());

		t.stop();
		if(DISPLAY_PROGRESS && !json){
			outstream.println("Added "+added+" kmers; time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
		
		if(verbose){
			TextStreamWriter tsw=new TextStreamWriter("stdout", false, false, false, FileFormat.TEXT);
			tsw.start();
			for(AbstractKmerTable table : keySets){
				table.dumpKmersAsText(tsw, k, 1, Integer.MAX_VALUE);
			}
			tsw.poisonAndWait();
		}
		
		return added;
	}

	/**
	 * Match reads against reference kmers, using multiple ProcessThread.
	 * @param t
	 */
	private void spawnProcessThreads(Timer t){
		t.start();
		
		/* Create read input stream */
		final ConcurrentReadInputStream cris;
		final boolean paired;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, ffin1.samOrBam(), ffin1, ffin2, qfin1, qfin2);
			cris.setSampleRate(samplerate, sampleseed);
			cris.start(); //4567
			paired=cris.paired();
			assert(ffin2!=null || !paired || (FASTQ.FORCE_INTERLEAVED || FASTQ.TEST_INTERLEAVED)) : 
				paired+", "+FASTQ.FORCE_INTERLEAVED+", "+FASTQ.TEST_INTERLEAVED+", "+ffin1.interleaved();
			if(!ffin1.samOrBam() && !silent){
				if(json){
					jsonStats.add("paired", paired);
				}else{
					outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));
				}
			}
		}
		
		/* Create read output streams */
		final ConcurrentReadOutputStream ros, rosb, ross;
		{
			final int buff=(!ordered ? 12 : Tools.max(32, 2*Shared.threads()));
			if(out1!=null){
				ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, true);
				ros.start();
			}else{ros=null;}
			if(outb1!=null){
				rosb=ConcurrentReadOutputStream.getStream(ffoutb1, ffoutb2, null, null, buff, null, true);
				rosb.start();
			}else{rosb=null;}
			if(outsingle!=null){
				ross=ConcurrentReadOutputStream.getStream(ffouts, null, null, null, buff, null, true);
				ross.start();
			}else{ross=null;}
			if(sidechannel!=null) {sidechannel.start();}
			if(ros!=null || rosb!=null || ross!=null || sidechannel!=null){
				t.stop();
				if(!silent && !json){outstream.println("Started output streams:\t"+t);}
				t.start();
			}
		}
		
		/* Optionally skip the first reads, since initial reads may have lower quality */
		if(skipreads>0){
			long skipped=0;

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(skipped<skipreads && reads!=null && reads.size()>0){
				skipped+=reads.size();
				
				if(rosb!=null){rosb.add(new ArrayList<Read>(1), ln.id);}
				if(ros!=null){ros.add(new ArrayList<Read>(1), ln.id);}
				if(ross!=null){ross.add(new ArrayList<Read>(1), ln.id);}
				
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);//if added for compiler benefit
			if(reads==null || reads.isEmpty()){
				ReadWrite.closeStreams(cris, ros, rosb, ross);
				if(sidechannel!=null) {sidechannel.shutdown();}
				outstream.println("Skipped all of the reads.");
				System.exit(0);
			}
		}
		
		/* Create ProcessThreads */
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new ProcessThread(cris, ros, rosb, ross, ALLOW_LOCAL_ARRAYS));}
		for(ProcessThread pt : alpt){pt.start();}
		
		/* Wait for threads to die, and gather statistics */
		for(ProcessThread pt : alpt){
			
			/* Wait for a thread to die */
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			
			/* Accumulate data from per-thread counters */
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			readsOut+=pt.readsOutuT;
			basesOut+=pt.basesOutuT;
			readsKFiltered+=pt.readsKFilteredT;
			basesKFiltered+=pt.basesKFilteredT;
			readsQTrimmed+=pt.readsQTrimmedT;
			basesQTrimmed+=pt.basesQTrimmedT;
			readsFTrimmed+=pt.readsFTrimmedT;
			basesFTrimmed+=pt.basesFTrimmedT;
			readsKTrimmed+=pt.readsKTrimmedT;
			basesKTrimmed+=pt.basesKTrimmedT;
			readsTrimmedBySwift+=pt.readsTrimmedBySwiftT;
			basesTrimmedBySwift+=pt.basesTrimmedBySwiftT;
			readsTrimmedByOverlap+=pt.readsTrimmedByOverlapT;
			basesTrimmedByOverlap+=pt.basesTrimmedByOverlapT;
			badGcReads+=pt.badGcReadsT;
			badGcBases+=pt.badGcBasesT;
			badHeaderReads+=pt.badHeaderReadsT;
			badHeaderBases+=pt.badHeaderBasesT;
			readsQFiltered+=pt.readsQFilteredT;
			basesQFiltered+=pt.basesQFilteredT;
			readsNFiltered+=pt.readsNFilteredT;
			basesNFiltered+=pt.basesNFilteredT;
			readsEFiltered+=pt.readsEFilteredT;
			basesEFiltered+=pt.basesEFilteredT;
			readsPolyTrimmed+=pt.readsPolyTrimmedT;
			basesPolyTrimmed+=pt.basesPolyTrimmedT;
			
			if(pTracker!=null){
				pTracker.add(pt.pTrackerT);
			}
			
			if(hitCounts!=null){
				for(int i=0; i<hitCounts.length; i++){hitCounts[i]+=pt.hitCountsT[i];}
				pt.hitCountsT=null;
			}
			if(pt.scaffoldReadCountsT!=null && scaffoldReadCounts!=null){
				for(int i=0; i<pt.scaffoldReadCountsT.length; i++){scaffoldReadCounts.addAndGet(i, pt.scaffoldReadCountsT[i]);}
				pt.scaffoldReadCountsT=null;
			}
			if(pt.scaffoldBaseCountsT!=null && scaffoldBaseCounts!=null){
				for(int i=0; i<pt.scaffoldBaseCountsT.length; i++){scaffoldBaseCounts.addAndGet(i, pt.scaffoldBaseCountsT[i]);}
				pt.scaffoldBaseCountsT=null;
			}
			errorState|=(!pt.finishedSuccessfully);
		}
		
		/* Shut down I/O streams; capture error status */
		{
			//Prevent a spurious error message in the event of a race condition when maxReads is set.
			boolean b=ReadWrite.closeStream(cris);
			if(maxReads<1 || maxReads==Long.MAX_VALUE || (maxReads!=readsIn && maxReads*2!=readsIn && samplerate<1)){errorState|=b;}
		}
		errorState|=ReadWrite.closeOutputStreams(ros, rosb, ross);
		if(sidechannel!=null) {errorState|=sidechannel.shutdown();}
		errorState|=ReadStats.writeAll();
		
		t.stop();
		if(showSpeed && !json){
			outstream.println("Processing time:   \t\t"+t);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads kmers into a table.  Each thread handles all kmers X such that X%WAYS==tnum.
	 */
	private class LoadThread extends Thread{
		
		public LoadThread(final int tnum_){
			tnum=tnum_;
			map=keySets[tnum];
		}
		
		/**
		 * Get the next list of reads (or scaffolds) from the queue.
		 * @return List of reads
		 */
		private ArrayList<Read> fetch(){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=queue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			return list;
		}
		
		@Override
		public void run(){
			ArrayList<Read> reads=fetch();
			while(reads!=POISON){
				for(Read r1 : reads){
					assert(r1.pairnum()==0);
					final Read r2=r1.mate;

					final int rblen=(r1==null ? 0 : r1.length());
					final int rblen2=r1.mateLength();
					
					addedT+=addToMap(r1, rblen>20000000 ? k : rblen>5000000 ? 11 : rblen>500000 ? 2 : 0);
					if(r2!=null){
						addedT+=addToMap(r2, rblen2>20000000 ? k : rblen2>5000000 ? 11 : rblen2>500000 ? 2 : 0);
					}
				}
				reads=fetch();
			}
			
			if(map.canRebalance() && map.size()>2L*map.arrayLength()){
				map.rebalance();
			}
			success=true;
		}

		/**
		 * Store the read's kmers in a table.
		 * @param r The current read to process
		 * @param skip Number of bases to skip between kmers
		 * @return Number of kmers stored
		 */
		private long addToMap(Read r, int skip){
			skip=Tools.max(minSkip, Tools.min(maxSkip, skip));
			final byte[] bases=r.bases;
			long kmer=0;
			long rkmer=0;
			long added=0;
			int len=0;
			
			if(bases!=null){
				refReadsT++;
				refBasesT+=bases.length;
			}
			if(bases==null || bases.length<k){return 0;}
			
			final int id=(Integer)r.obj;
			
			if(skip>1){ //Process while skipping some kmers
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=symbolToNumber0[b];
					long x2=symbolToComplementNumber0[b];
					kmer=((kmer<<bitsPerBase)|x)&mask;
					rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
					if(isFullyDefined(b)){len++;}else{len=0; rkmer=0;}
					if(verbose){outstream.println("Scanning1 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=k){
						refKmersT++;
						if(len%skip==0){
							final long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber[bases[i+1]]);
							added+=addToMap(kmer, rkmer, k, extraBase, id, kmask, hammingDistance, editDistance);
							if(useShortKmers){
								if(i==k2){added+=addToMapRightShift(kmer, rkmer, id);}
								if(i==bases.length-1){added+=addToMapLeftShift(kmer, rkmer, extraBase, id);}
							}
						}
					}
				}
			}else{ //Process all kmers
				for(int i=0; i<bases.length; i++){
					final byte b=bases[i];
					final long x=symbolToNumber0[b];
					final long x2=symbolToComplementNumber0[b];
//					assert(x!=x2) : x+", "+x2+", "+Character.toString((char)b)+"\n"+Arrays.toString(symbolToNumber0)+"\n"+Arrays.toString(symbolToComplementNumber);
					kmer=((kmer<<bitsPerBase)|x)&mask;
					//10000, 1111111111, 16, 16, 2, 10, 8
					rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
					if(isFullyDefined(b)){len++;}else{len=0; rkmer=0;}
					if(verbose){
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
					}
					if(len>=k){
//						assert(kmer==rcomp(rkmer, k)) : Long.toBinaryString(kmer)+", "+Long.toBinaryString(rkmer)+", "+Long.toBinaryString(mask)+", x="+x+", x2="+x2+", bits="+bitsPerBase+", s="+shift+", s2="+shift2+", b="+Character.toString((char)b);
						refKmersT++;
						final long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber[bases[i+1]]);
						final long atm=addToMap(kmer, rkmer, k, extraBase, id, kmask, hammingDistance, editDistance);
						added+=atm;
//						assert(false) : atm+", "+map.contains(toValue(kmer, rkmer, kmask));
						if(useShortKmers){
							if(i==k2){added+=addToMapRightShift(kmer, rkmer, id);}
							if(i==bases.length-1){added+=addToMapLeftShift(kmer, rkmer, extraBase, id);}
						}
					}
				}
			}
			return added;
		}
		

		/**
		 * Adds short kmers on the left end of the read.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @return Number of kmers stored
		 */
		private long addToMapLeftShift(long kmer, long rkmer, final long extraBase, final int id){
			if(verbose){outstream.println("addToMapLeftShift");}
			long added=0;
			for(int i=k-1; i>=mink; i--){
				kmer=kmer&rightMasks[i];
				rkmer=rkmer>>>bitsPerBase;
				long x=addToMap(kmer, rkmer, i, extraBase, id, lengthMasks[i], hammingDistance2, editDistance2);
				added+=x;
				if(verbose){
					if((toValue(kmer, rkmer, lengthMasks[i]))%WAYS==tnum){
						outstream.println("added="+x+"; i="+i+"; tnum="+tnum+"; Added left-shift kmer "+kmerToString(kmer&~lengthMasks[i], i)+"; value="+(toValue(kmer, rkmer, lengthMasks[i]))+"; kmer="+kmer+"; rkmer="+rkmer+"; kmask="+lengthMasks[i]+"; rightMasks[i+1]="+rightMasks[i+1]);
						outstream.println("i="+i+"; tnum="+tnum+"; Looking for left-shift kmer "+kmerToString(kmer&~lengthMasks[i], i));
						final long value=toValue(kmer, rkmer, lengthMasks[i]);
						if(map.contains(value)){outstream.println("Found "+value);}
					}
				}
			}
			return added;
		}
		

		/**
		 * Adds short kmers on the right end of the read.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param id Scaffold number
		 * @return Number of kmers stored
		 */
		private long addToMapRightShift(long kmer, long rkmer, final int id){
			if(verbose){outstream.println("addToMapRightShift");}
			long added=0;
			for(int i=k-1; i>=mink; i--){
				long extraBase=kmer&symbolMask;
				kmer=kmer>>>bitsPerBase;
				rkmer=rkmer&rightMasks[i];
//				assert(Long.numberOfLeadingZeros(kmer)>=2*(32-i)) : Long.numberOfLeadingZeros(kmer)+", "+i+", "+kmer+", "+kMasks[i];
//				assert(Long.numberOfLeadingZeros(rkmer)>=2*(32-i)) : Long.numberOfLeadingZeros(rkmer)+", "+i+", "+rkmer+", "+kMasks[i];
				long x=addToMap(kmer, rkmer, i, extraBase, id, lengthMasks[i], hammingDistance2, editDistance2);
				added+=x;
				if(verbose){
					if((toValue(kmer, rkmer, lengthMasks[i]))%WAYS==tnum){
						outstream.println("added="+x+"; i="+i+"; tnum="+tnum+"; Added right-shift kmer "+kmerToString(kmer&~lengthMasks[i], i)+"; value="+(toValue(kmer, rkmer, lengthMasks[i]))+"; kmer="+kmer+"; rkmer="+rkmer+"; kmask="+lengthMasks[i]+"; rightMasks[i+1]="+rightMasks[i+1]);
						outstream.println("i="+i+"; tnum="+tnum+"; Looking for right-shift kmer "+kmerToString(kmer&~lengthMasks[i], i));
						final long value=toValue(kmer, rkmer, lengthMasks[i]);
						if(map.contains(value)){outstream.println("Found "+value);}
					}
				}
			}
			return added;
		}
		
		
		/**
		 * Adds this kmer to the table, including any mutations implied by editDistance or hammingDistance.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param len Kmer length
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @param kmask0
		 * @return Number of kmers stored
		 */
		private long addToMap(final long kmer, final long rkmer, final int len, final long extraBase, final int id, final long kmask0, final int hdist, final int edist){
			
			assert(kmask0==lengthMasks[len]) : kmask0+", "+len+", "+lengthMasks[len]+", "+Long.numberOfTrailingZeros(kmask0)+", "+Long.numberOfTrailingZeros(lengthMasks[len]);
			
			if(verbose){outstream.println("addToMap_A; len="+len+"; kMasks[len]="+lengthMasks[len]);}
			assert((kmer&kmask0)==0);
			final long added;
			if(hdist==0){
				final long key=toValue(kmer, rkmer, kmask0);
				if(verbose){outstream.println("toValue ("+kmerToString(kmer, len)+", "+kmerToString(rkmer, len)+") = "+kmerToString(key, len)+" = "+key);}
				if(failsSpeed(key)){return 0;}
				if(key%WAYS!=tnum){return 0;}
				if(verbose){outstream.println("addToMap_B: "+kmerToString(key, len)+" ("+key+")");}
				added=map.setIfNotPresent(key, id);
			}else if(edist>0){
//				long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber2bases[i+1]]);
				added=mutate(kmer, rkmer, len, id, edist, extraBase);
			}else{
				added=mutate(kmer, rkmer, len, id, hdist, -1);
			}
			if(verbose){outstream.println("addToMap added "+added+" keys.");}
			return added;
		}
		
		/**
		 * Mutate and store this kmer through 'dist' recursions.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param id Scaffold number
		 * @param dist Number of mutations
		 * @param extraBase Base added to end in case of deletions
		 * @return Number of kmers stored
		 */
		private long mutate(final long kmer, final long rkmer, final int len, final int id, final int dist, final long extraBase){
			long added=0;
			
			final long key=toValue(kmer, rkmer, lengthMasks[len]);
			
//			if(dist==1){System.err.println(".\t.\t"+kmerToString(kmer, k)+" initial.");}//123
			
			if(verbose){outstream.println("mutate_A; len="+len+"; kmer="+kmer+"; rkmer="+rkmer+"; kMasks[len]="+lengthMasks[len]);}
			if(key%WAYS==tnum){
				if(verbose){outstream.println("mutate_B: "+kmerToString(kmer&~lengthMasks[len], len)+" = "+key);}
				int x=map.setIfNotPresent(key, id);
//				if(x>0){System.err.println(".\t.\t"+kmerToString(kmer, k)+" Added!");}//123
				if(verbose){outstream.println("mutate_B added "+x+" keys.");}
				added+=x;
				assert(map.contains(key));
			}
			
			if(dist>0){
				final int dist2=dist-1;
				
				//Sub
				for(int j=0; j<symbols; j++){
					for(int i=0; i<len; i++){
						final long temp=(kmer&clearMasks[i])|setMasks[j][i]; //TODO:  6/14/23, fixed incorrect description of setMasks that swapped i and j; may need changing in Seal and etc
//						System.err.println("cm:"+kmerToString(clearMasks[i], k));//123
//						System.err.println("sm:"+kmerToString(setMasks[j][i], k));//123
//						assert(Long.bitCount((temp^kmer))<=bitsPerBase) : //Warning: Slow assertion //123
//							"\n"+kmerToString(kmer, k)+"\n"+kmerToString(temp, k)+"\n"+i+","+j+","+k;
//						System.err.println(j+"\t"+i+"\t"+kmerToString(temp, k));//123
						if(temp!=kmer){
							long rtemp=rcomp(temp, len);
							added+=mutate(temp, rtemp, len, id, dist2, extraBase);
						}
					}
				}
				
				if(editDistance>0){
					//Del
					if(extraBase>=0 && extraBase<=maxSymbol){
						for(int i=1; i<len; i++){
							final long temp=(kmer&leftMasks[i])|((kmer<<bitsPerBase)&rightMasks[i])|extraBase;
							if(temp!=kmer){
								long rtemp=rcomp(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, -1);
							}
						}
					}

					//Ins
					final long eb2=kmer&symbolMask;
					for(int i=1; i<len; i++){
						final long temp0=(kmer&leftMasks[i])|((kmer&rightMasks[i])>>bitsPerBase);
						for(int j=0; j<symbols; j++){
							final long temp=temp0|setMasks[j][i-1];
							if(temp!=kmer){
								long rtemp=rcomp(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, eb2);
							}
						}
					}
				}

			}
			
//			if(dist==1){//123
//				System.err.println("Added "+added);
//				assert(false);
//			}
			
			return added;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Number of kmers stored by this thread */
		public long addedT=0;
		/** Number of items encountered by this thread */
		public long refKmersT=0, refReadsT=0, refBasesT=0;
		/** Thread number; used to determine which kmers to store */
		public final int tnum;
		/** Buffer of input read lists */
		public final ArrayBlockingQueue<ArrayList<Read>> queue=new ArrayBlockingQueue<ArrayList<Read>>(32);
//		/** Used to trick compiler */
//		public long modsumT=0; //123
		
		/** Destination for storing kmers */
		private final AbstractKmerTable map;
		
		/** Completed successfully */
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/

	/**
	 * Matches read kmers against reference kmers, performs binning and/or trimming, and writes output.
	 */
	private class ProcessThread extends Thread{
		
		/**
		 * Constructor
		 * @param cris_ Read input stream
		 * @param ros_ Unmatched read output stream (optional)
		 * @param rosb_ Matched read output stream (optional)
		 * @param ross_ Singleton read output stream (optional)
		 */
		public ProcessThread(ConcurrentReadInputStream cris_, ConcurrentReadOutputStream ros_, ConcurrentReadOutputStream rosb_, ConcurrentReadOutputStream ross_, boolean localArrays){
			cris=cris_;
			ros=ros_;
			rosb=rosb_;
			ross=ross_;
			
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
			
			hitCountsT=(hitCounts==null ? null : new long[hitCounts.length]);
			
			if(localArrays && alen>0 && alen<10000 && scaffoldReadCounts!=null && scaffoldBaseCounts!=null){
				scaffoldReadCountsT=new long[alen];
				scaffoldBaseCountsT=new long[alen];
			}else{
				scaffoldReadCountsT=scaffoldBaseCountsT=null;
			}
			
			if(calcEntropy){
				eTrackerT=new EntropyTracker(amino, Tools.max(0, entropyCutoff), entropyHighpass);
			}else{
				eTrackerT=null;
			}
			
			if(countPolymers){
				pTrackerT=new PolymerTracker();
			}else{
				pTrackerT=null;
			}
			
			maxBasesOutmT=(maxBasesOutm>0 ? Tools.max(1, maxBasesOutm/THREADS) : -1);
			maxBasesOutuT=(maxBasesOutu>0 ? Tools.max(1, maxBasesOutu/THREADS) : -1);

			flowCoords=(locationFilter ? new FlowcellCoordinate() : null);
		}
		
		@Override
		public void run(){
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

					readsInT+=pairCount;
					basesInT+=initialPairLength;

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
							badHeaderBasesT+=initialPairLength;
							badHeaderReadsT+=pairCount;
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
								badGcBasesT+=initialLength1;
								badGcReadsT++;
							}
							if(isNotDiscarded(r2) && (gc2<minGC || gc2>maxGC)){
								setDiscarded(r2);
								badGcBasesT+=initialLength2;
								badGcReadsT++;
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
								basesFTrimmedT+=x;
								readsFTrimmedT+=(x>0 ? 1 : 0);
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
								basesFTrimmedT+=x;
								readsFTrimmedT+=(x>0 ? 1 : 0);
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
							basesQFilteredT+=initialLength1;
							readsQFilteredT++;
						}
						if(r2!=null){
							basesQFilteredT+=initialLength2;
							readsQFilteredT++;
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
									int x=ktrimTips(r1, keySets);
									xsum+=x;
									rktsum+=(x>0 ? 1 : 0);
									rlen1=r1.length();
									if(rlen1<minlen1){setDiscarded(r1);}
								}
								if(r2!=null){
									int x=ktrimTips(r2, keySets);
									xsum+=x;
									rktsum+=(x>0 ? 1 : 0);
									rlen2=r2.length();
									if(rlen2<minlen2){setDiscarded(r2);}
								}
							}else if(ktrimLeftOrRight){
								if(r1!=null){
									int x=ktrim(r1, keySets);
									xsum+=x;
									rktsum+=(x>0 ? 1 : 0);
									rlen1=r1.length();
									if(rlen1<minlen1){setDiscarded(r1);}
								}
								if(r2!=null){
									int x=ktrim(r2, keySets);
									xsum+=x;
									rktsum+=(x>0 ? 1 : 0);
									rlen2=r2.length();
									if(rlen2<minlen2){setDiscarded(r2);}
								}
							}else if(ktrimN){
								if(r1!=null){
									int x=kmask(r1, keySets);
									xsum+=x;
									rktsum+=(x>0 ? 1 : 0);
									rlen1=r1.length();
									if(rlen1<minlen1){setDiscarded(r1);}
								}
								if(r2!=null){
									int x=kmask(r2, keySets);
									xsum+=x;
									rktsum+=(x>0 ? 1 : 0);
									rlen2=r2.length();
									if(rlen2<minlen2){setDiscarded(r2);}
								}
							}else if(ksplit){
								assert(r2==null);
								if(r1!=null){
									int oldLen=r1.pairLength();
									boolean b=ksplit(r1, keySets);
									int trimmed=oldLen-r1.pairLength();
									xsum+=trimmed;
									rktsum+=(trimmed>0 ? 1 : 0);
									rlen1=r1.length();
//									if(rlen1<minlen1){setDiscarded(r1);}
//									if(b && r1.mateLength()<minlen1){setDiscarded(r1);}
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
							basesKTrimmedT+=xsum;
							readsKTrimmedT+=rktsum;
							
						}else if(doKmerFiltering){
							//Do kmer matching
							
							if(minCoveredFraction>0){
								if(isNotDiscarded(r1)){
									final int minCoveredBases=(int)Math.ceil(minCoveredFraction*r1.length());
									final int covered=countCoveredBases(r1, keySets, minCoveredBases);
									if(covered>=minCoveredBases){setDiscarded(r1);}
								}
								if(isNotDiscarded(r2)){
									final int minCoveredBases=(int)Math.ceil(minCoveredFraction*r2.length());
									final int covered=countCoveredBases(r2, keySets, minCoveredBases);
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
									final int a=(kbig<=k ? countSetKmers(r1, keySets, maxBadKmersR1) : countSetKmersBig(r1, keySets, maxBadKmersR1));
									final int b=(kbig<=k ? countSetKmers(r2, keySets, maxBadKmersR2) : countSetKmersBig(r2, keySets, maxBadKmersR2));

									if(r1!=null && a>maxBadKmersR1){setDiscarded(r1);}
									if(r2!=null && b>maxBadKmersR2){setDiscarded(r2);}

								}else{
									final int a=findBestMatch(r1, keySets, maxBadKmersR1);
									final int b=findBestMatch(r2, keySets, maxBadKmersR2);

									if(r1!=null && a>0){setDiscarded(r1);}
									if(r2!=null && b>0){setDiscarded(r2);}
								}
							}
							
							if(shouldRemove(r1, r2)){
								remove=true;
								if(r1!=null){
									readsKFilteredT++;
									basesKFilteredT+=initialLength1;
								}
								if(r2!=null){
									readsKFilteredT++;
									basesKFilteredT+=initialLength2;
								}
								if(bad!=null){bad.add(r1);}
							}
							
						}
					}
					
//					assert(false) : remove+", "+trimByOverlap+", "+(r2!=null);
					
					if(!remove && trimByOverlap && r2!=null && expectedErrors(r1, r2)<meeFilter){
						
						if(aprob==null || aprob.length<r1.length()){aprob=new float[r1.length()];}
						if(bprob==null || bprob.length<r2.length()){bprob=new float[r2.length()];}
						
						//Do overlap trimming
						r2.reverseComplement();
//						int bestInsert=BBMergeOverlapper.mateByOverlap(r1, r2, aprob, bprob, overlapVector, minOverlap0, minOverlap,
//								overlapMargin, overlapMaxMismatches0, overlapMaxMismatches, overlapMinq);
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
						
						r2.reverseComplement();
						
						if(bestInsert>0 && !ambig){
							if(bestInsert<r1.length()){
								if(verbose){outstream.println("Overlap right trimming r1 to "+0+", "+(bestInsert-1));}
								int x=TrimRead.trimToPosition(r1, 0, bestInsert-1, 1);
								if(verbose){outstream.println("Trimmed "+x+" bases: "+new String(r1.bases));}
								readsTrimmedByOverlapT++;
								basesTrimmedByOverlapT+=x;
							}
							if(bestInsert<r2.length()){
								if(verbose){outstream.println("Overlap right trimming r2 to "+0+", "+(bestInsert-1));}
								int x=TrimRead.trimToPosition(r2, 0, bestInsert-1, 1);
								if(verbose){outstream.println("Trimmed "+x+" bases: "+new String(r2.bases));}
								readsTrimmedByOverlapT++;
								basesTrimmedByOverlapT+=x;
							}
						}
					}
					
					if(!remove && swift){
						//Do Swift trimming
						
						int rlen1=0, rlen2=0;
						if(r1!=null){
							int x=trimSwift(r1);
							basesTrimmedBySwiftT+=x;
							readsTrimmedBySwiftT+=(x>0 ? 1 : 0);
							rlen1=r1.length();
							if(rlen1<minlen1){setDiscarded(r1);}
						}
						if(r2!=null){
							int x=trimSwift(r2);
							basesTrimmedBySwiftT+=x;
							readsTrimmedBySwiftT+=(x>0 ? 1 : 0);
							rlen2=r2.length();
							if(rlen2<minlen2){setDiscarded(r2);}
						}
						
						//Discard reads if too short
						if(shouldRemove(r1, r2)){
							basesTrimmedBySwiftT+=r1.pairLength();
							remove=true;
							if(addTrimmedToBad && bad!=null){bad.add(r1);}
						}
					}
					
					if(!remove && trimPolyA>0){
						//Do poly-A trimming
						
						int rlen1=0, rlen2=0;
						if(r1!=null){
							int x=trimPolyA(r1, trimPolyA);
							basesPolyTrimmedT+=x;
							readsPolyTrimmedT+=(x>0 ? 1 : 0);
							rlen1=r1.length();
							if(rlen1<minlen1){setDiscarded(r1);}
						}
						if(r2!=null){
							int x=trimPolyA(r2, trimPolyA);
							basesPolyTrimmedT+=x;
							readsPolyTrimmedT+=(x>0 ? 1 : 0);
							rlen2=r2.length();
							if(rlen2<minlen2){setDiscarded(r2);}
						}
						
						//Discard reads if too short
						if(shouldRemove(r1, r2)){
							basesPolyTrimmedT+=r1.pairLength();
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
								readsPolyTrimmedT++;
							}else if(trimPolyGLeft>0 || trimPolyGRight>0){
								int x=trimPoly(r1, trimPolyGLeft, trimPolyGRight, maxNonPoly, (byte)'G');
								basesPolyTrimmedT+=x;
								readsPolyTrimmedT+=(x>0 ? 1 : 0);
								rlen1=r1.length();
								if(rlen1<minlen1){setDiscarded(r1);}
							}
						}
						if(r2!=null){
							if(filterPolyG>0 && detectPolyLeft(r2, filterPolyG, maxNonPoly, (byte)'G')>=filterPolyG) {
								setDiscarded(r2);
								readsPolyTrimmedT++;
							}else if(trimPolyGLeft>0 || trimPolyGRight>0){
								int x=trimPoly(r2, trimPolyGLeft, trimPolyGRight, maxNonPoly, (byte)'G');
								basesPolyTrimmedT+=x;
								readsPolyTrimmedT+=(x>0 ? 1 : 0);
								rlen2=r2.length();
								if(rlen2<minlen2){setDiscarded(r2);}
							}
						}
						
						//Discard reads if too short
						if(shouldRemove(r1, r2)){
							basesPolyTrimmedT+=r1.pairLength();
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
								readsPolyTrimmedT++;
							}else if(trimPolyCLeft>0 || trimPolyCRight>0){
								int x=trimPoly(r1, trimPolyCLeft, trimPolyCRight, maxNonPoly, (byte)'C');
								basesPolyTrimmedT+=x;
								readsPolyTrimmedT+=(x>0 ? 1 : 0);
								rlen1=r1.length();
								if(rlen1<minlen1){setDiscarded(r1);}
							}
						}
						if(r2!=null){
							if(filterPolyC>0 && detectPolyLeft(r1, filterPolyC, maxNonPoly, (byte)'C')>=filterPolyC) {
								setDiscarded(r2);
								readsPolyTrimmedT++;
							}else if(trimPolyCLeft>0 || trimPolyCRight>0){
								int x=trimPoly(r2, trimPolyCLeft, trimPolyCRight, maxNonPoly, (byte)'C');
								basesPolyTrimmedT+=x;
								readsPolyTrimmedT+=(x>0 ? 1 : 0);
								rlen2=r2.length();
								if(rlen2<minlen2){setDiscarded(r2);}
							}
						}
						
						//Discard reads if too short
						if(shouldRemove(r1, r2)){
							basesPolyTrimmedT+=r1.pairLength();
							remove=true;
							if(addTrimmedToBad && bad!=null){bad.add(r1);}
						}
					}
					
					if(!remove && (entropyMask || entropyTrim>0)){
						//Mask/trim entropy
						if(isNotDiscarded(r1)){
							int masked=(entropyTrim>0 ? trimLowEntropy(r1, null, eTrackerT) : maskLowEntropy(r1, null, eTrackerT));
							basesEFilteredT+=masked;
							readsEFilteredT+=(masked>0 ? 1 : 0);
						}
						if(isNotDiscarded(r2)){
							int masked=(entropyTrim>0 ? trimLowEntropy(r2, null, eTrackerT) : maskLowEntropy(r2, null, eTrackerT));
							basesEFilteredT+=masked;
							readsEFilteredT+=(masked>0 ? 1 : 0);
						}
					}
					
					if(entropyMark){
						markLowEntropy(r1, eTrackerT);
						markLowEntropy(r2, eTrackerT);
					}
					
					if(!remove){
						//Do quality trimming
						
						if(qtrimLeft || qtrimRight || trimClip){
							if(r1!=null){
								int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, trimE, 1, trimClip);
								basesQTrimmedT+=x;
								readsQTrimmedT+=(x>0 ? 1 : 0);

//								assert(false) : trimClip+", "+x;
							}
							if(r2!=null){
								int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, trimE, 1, trimClip);
								basesQTrimmedT+=x;
								readsQTrimmedT+=(x>0 ? 1 : 0);
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
							basesQTrimmedT+=r1.pairLength();
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
								readsNFilteredT++;
								basesNFilteredT+=r1.length();
								setDiscarded(r1);
							}
							if(r2!=null && r2.countUndefined()>maxNs){
								readsNFilteredT++;
								basesNFilteredT+=r2.length();
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
							basesQFilteredT+=r1.pairLength();
							readsQFilteredT+=pairCount;
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
						if(isNotDiscarded(r1) && !eTrackerT.passes(r1.bases, true)){setDiscarded(r1);}
						if(isNotDiscarded(r2) && !eTrackerT.passes(r2.bases, true)){setDiscarded(r2);}
						
						if(shouldRemove(r1, r2)){
							basesEFilteredT+=r1.pairLength();
							readsEFilteredT+=pairCount;
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
						
						readsOutmT+=pairCount;
						basesOutmT+=r1.pairLength();
					}else{
						if(loglogOut!=null){loglogOut.hash(r1);}
						readsOutuT+=pairCount;
						basesOutuT+=r1.pairLength();
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

				if(maxBasesOutmT>=0 && basesOutmT>=maxBasesOutmT){break;}
				if(maxBasesOutuT>=0 && basesOutuT>=maxBasesOutuT){break;}
				
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
			
			if(pTrackerT!=null){
				pTrackerT.addPair(r1);
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
		 * Transforms a kmer into all canonical values for a given Hamming distance.
		 * Returns the related id stored in the tables.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param lengthMask Bitmask with single '1' set to left of kmer
		 * @param qPos Position of kmer in query
		 * @param len kmer length
		 * @param qHDist Hamming distance
		 * @param sets Kmer hash tables
		 * @return Value stored in table, or -1
		 */
		private final int getValue(final long kmer, final long rkmer, final long lengthMask, final int qPos, final int len, final int qHDist, final AbstractKmerTable[] sets){
			assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
			if(verbose){outstream.println("getValue()");}
			int id=getValueInner(kmer, rkmer, lengthMask, len, qPos, sets);
			if(id<1 && qHDist>0){
				final int qHDist2=qHDist-1;
				
				//Sub
				for(int j=0; j<symbols && id<1; j++){
					for(int i=0; i<len && id<1; i++){
						final long temp=(kmer&clearMasks[i])|setMasks[j][i];
//						outstream.println(i+", "+j+", "+setMasks[j][i]+", "+qHDist);
						if(temp!=kmer){
							long rtemp=rcomp(temp, len);
//							assert(lengthMask==0 || (temp<lengthMask && rtemp<lengthMask)) : lengthMask+", "+temp+", "+rtemp+", "+kmer+", "+rkmer+
//							"\n"+len+", "+Long.numberOfTrailingZeros(lengthMask)+"\n"+
//									Long.toBinaryString(lengthMask|0x8000000000000000L)+"\n"+
//											Long.toBinaryString(temp|0x8000000000000000L)+"\n"+
//													Long.toBinaryString(rtemp|0x8000000000000000L);
							id=getValue(temp, rtemp, lengthMask, qPos, len, qHDist2, sets);
						}
					}
				}
			}
			return id;
		}
		
		/**
		 * Transforms a kmer into a canonical value stored in the table and search.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param lengthMask Bitmask with single '1' set to left of kmer
		 * @param qPos Position of kmer in query
		 * @param sets Kmer hash tables
		 * @return Value stored in table
		 */
		private final int getValueInner(final long kmer, final long rkmer, final long lengthMask, final int len, final int qPos, final AbstractKmerTable[] sets){
			assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
			if(qSkip>1 && (qPos%qSkip!=0)){return -1;}

			if(verbose){
				outstream.println("getValueInner(kmer="+AminoAcid.kmerToString(kmer, len)+", rkmer="+AminoAcid.kmerToString(rkmer, len)+", len="+len+", mask="+lengthMask+")");
				outstream.println("getValueInner(kmer="+kmer+", rkmer="+rkmer+", len="+len+", mask="+lengthMask+")");
			}
			final long max=(rcomp ? Tools.max(kmer, rkmer) : kmer);
			if(verbose){outstream.println("max="+AminoAcid.kmerToString(max, len)+" ("+max+")");}
			final long key=(max&middleMask)|lengthMask;
			if(verbose){outstream.println("key="+AminoAcid.kmerToString(key, len)+" ("+key+")");}
			if(passesSpeed(key)){
				if(verbose){outstream.println("Testing key "+kmerToString(key, len)+" ("+key+")");}
				AbstractKmerTable set=sets[(int)(key%WAYS)];
				final int id=set.getValue(key);
				if(verbose){outstream.println("getValueInner("+kmerToString(kmer, len)+", "+kmerToString(rkmer, len)+") > "+kmerToString(key, len)+" ("+key+") = "+id);}
				return id;
			}
			if(verbose){outstream.println("Invalid key.");}
			return -1;
		}
		
		
		/**
		 * Counts the number of kmer hits for a read.
		 * @param r Read to process
		 * @param sets Kmer tables
		 * @return Number of hits
		 */
		private final int countSetKmers(final Read r, final AbstractKmerTable[] sets, final int maxBadKmers){
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
					final int id=getValue(kmer, rkmer, kmask, i, k, qHammingDistance, sets);
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
							if(hitCounts==null){
								return (found=found+1);
							}//Early exit, but prevents generation of histogram that goes over maxBadKmers+1.
						}
						found++;
					}
				}
			}
			
			if(hitCountsT!=null){hitCountsT[Tools.min(found, HITCOUNT_LEN)]++;}
			return found;
		}
		
		
		/**
		 * Counts the number of kmer hits for a read.
		 * @param r Read to process
		 * @param sets Kmer tables
		 * @return Number of hits
		 */
		private final int countCoveredBases(final Read r, final AbstractKmerTable[] sets, final int minCoveredBases){
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
					final int id=getValue(kmer, rkmer, kmask, i, k, qHammingDistance, sets);
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
							if(hitCounts==null){
								return found;
							}
						}
					}
				}
			}
			
			if(hitCountsT!=null){hitCountsT[Tools.min(found, HITCOUNT_LEN)]++;}
			return found;
		}
		
		/**
		 * Returns the id of the sequence with the most kmer matches to this read, or -1 if none are over maxBadKmers.
		 * @param r Read to process
		 * @param sets Kmer tables
		 * @return id of best match
		 */
		private final int findBestMatch(final Read r, final AbstractKmerTable[] sets, final int maxBadKmers){
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
					final int id=getValue(kmer, rkmer, kmask, i, k, qHammingDistance, sets);
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
			
			if(hitCountsT!=null){hitCountsT[Tools.min(found, HITCOUNT_LEN)]++;}
			return id;
		}
		
		/** Estimates kmer hit counts for kmers longer than k using consecutive matches
		 * @param r
		 * @param sets
		 * @return Number of sets of consecutive hits of exactly length kbig
		 */
		private final int countSetKmersBig(final Read r, final AbstractKmerTable[] sets, final int maxBadKmers){
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
					id=getValue(kmer, rkmer, kmask, i, k, qHammingDistance, sets);
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
									if(hitCounts==null){
										return found;
									}//Early exit, but prevents generation of histogram that goes over maxBadKmers+1.
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
			
			if(hitCountsT!=null){hitCountsT[Tools.min(found, HITCOUNT_LEN)]++;}
			return found;
		}
		
		private final int ktrim(final Read r, final AbstractKmerTable[] sets){
			final int len=r.length();
			final int start=(restrictRight<1 ? 0 : Tools.max(0, len-restrictRight));
			final int stop=(restrictLeft<1 ? len : Tools.min(len, restrictLeft));
			return ktrim(r, sets, start, stop);
		}
		
		private final int ktrimTips(final Read r, final AbstractKmerTable[] sets){
			final int len=r.length();
			final int mid=len/2-(k-1)/2;
			int sum=0;
			if(ktrimRight){
				int start=Tools.max(0, (restrictRight<1 ? mid : len-restrictRight));
				sum+=ktrimTip(r, sets, start, len, true, false);
			}
			if(ktrimLeft){
				int stop=Tools.min(r.length(), (restrictLeft<1 ? mid+k-1 : restrictLeft));
				sum+=ktrimTip(r, sets, 0, stop, false, true);
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
		private final int ktrimTip(final Read r, final AbstractKmerTable[] sets, final int start, final int stop,
				final boolean right, final boolean left){
			assert(left || right);
			assert(!(left && right));
			if(r==null || r.length()<Tools.max(1, (useShortKmers ? Tools.min(k, mink) : k)) || storedKmers<1){return 0;}
			if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){return 0;}
			if(verbose){outstream.println("KTrimming read "+r.id);}
			final byte[] bases=r.bases, quals=r.quality;
			long kmer=0;
			long rkmer=0;
			int found=0;
			int len=0;
			int id0=-1; //ID of first kmer found.
			
			int minLoc=999999999, minLocExclusive=999999999;
			int maxLoc=-1, maxLocExclusive=-1;
			final int initialLength=r.length();
			
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
					final int id=getValue(kmer, rkmer, kmask, i, k, qHammingDistance, sets);
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
							
							final int id=getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2, sets);
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
							final int id=getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2, sets);
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
		private final int ktrim(final Read r, final AbstractKmerTable[] sets, final int start, final int stop){
			assert(ktrimLeft || ktrimRight);
			if(r==null || r.length()<Tools.max(1, (useShortKmers ? Tools.min(k, mink) : k)) || storedKmers<1){return 0;}
			if((skipR1 && r.pairnum()==0) || (skipR2 && r.pairnum()==1)){return 0;}
			if(verbose){outstream.println("KTrimming read "+r.id);}
			final byte[] bases=r.bases, quals=r.quality;
			long kmer=0;
			long rkmer=0;
			int found=0;
			int len=0;
			int id0=-1; //ID of first kmer found.
			
			int minLoc=999999999, minLocExclusive=999999999;
			int maxLoc=-1, maxLocExclusive=-1;
			final int initialLength=r.length();
			
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
					final int id=getValue(kmer, rkmer, kmask, i, k, qHammingDistance, sets);
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
							
							final int id=getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2, sets);
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
							final int id=getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2, sets);
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
		private final int kmask(final Read r, final AbstractKmerTable[] sets){
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
						id=getValue(kmer, rkmer, kmask, i, k, qHammingDistance, sets);
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
								id=getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2, sets);
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
								id=getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2, sets);
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
//			int y=r.countNocalls();
			int cardinality=bs.cardinality();
//			assert(cardinality>0);
			
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
//			assert(cardinality==r.countNocalls() || y>0) : cardinality+", "+r.countNocalls()+"\n"+r.length()+"\n"+bs+"\n"+r;//123
			return cardinality;
		}
		
		
		/**
		 * Mask a read to cover matching kmers.
		 * @param r Read to process
		 * @param sets Kmer tables
		 * @return Number of bases masked
		 */
		private final boolean ksplit(final Read r, final AbstractKmerTable[] sets){
			if(r==null || r.length()<Tools.max(1, (useShortKmers ? Tools.min(k, mink) : k)) || storedKmers<1){return false;}
			assert(r.mate==null) : "Kmer splitting should only be performed on unpaired reads.";
			assert(ksplit);
			if(verbose){outstream.println("KSplitting read "+r.id);}
			final byte[] bases=r.bases, quals=r.quality;
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
						id=getValue(kmer, rkmer, kmask, i, k, qHammingDistance, sets);
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
								id=getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2, sets);
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
								id=getValue(kmer, rkmer, lengthMasks[len], i, len, qHammingDistance2, sets);
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
//			System.err.println(new String(bases)+"\nleft="+left+", right="+right);
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
						 bases[i]=(byte)Tools.toLowerCase(bases[i]);
						 //Don't change quality
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
		
		/** Input read stream */
		private final ConcurrentReadInputStream cris;
		/** Output read streams */
		private final ConcurrentReadOutputStream ros, rosb, ross;
		
		private final FlowcellCoordinate flowCoords;
		
		private final ReadStats readstats;
		private final int[] overlapVector;
		private final int[] countArray;
		
		private final IntList idList;
		private final IntList countList;
		
		//These "*T" fields are used to store counts on a per-thread basis.
		
		long[] hitCountsT;
		long[] scaffoldReadCountsT;
		long[] scaffoldBaseCountsT;
		
		final EntropyTracker eTrackerT;
		final PolymerTracker pTrackerT;
		
		private float[] aprob, bprob;
		
		long readsInT=0;
		long basesInT=0;
		long readsOutuT=0;
		long basesOutuT=0;
		
		long readsOutmT=0;
		long basesOutmT=0;

		final long maxBasesOutmT;
		final long maxBasesOutuT;
		
		long readsQTrimmedT=0;
		long basesQTrimmedT=0;
		long readsFTrimmedT=0;
		long basesFTrimmedT=0;
		long readsQFilteredT=0;
		long basesQFilteredT=0;
		long readsNFilteredT=0;
		long basesNFilteredT=0;
		long readsEFilteredT=0;
		long basesEFilteredT=0;
		long readsPolyTrimmedT=0;
		long basesPolyTrimmedT=0;

		long readsKTrimmedT=0;
		long basesKTrimmedT=0;
		long readsKFilteredT=0;
		long basesKFilteredT=0;
		
		long readsTrimmedBySwiftT=0;
		long basesTrimmedBySwiftT=0;
		
		long readsTrimmedByOverlapT=0;
		long basesTrimmedByOverlapT=0;
		
		long badGcBasesT=0;
		long badGcReadsT=0;
		
		long badHeaderBasesT=0;
		long badHeaderReadsT=0;
		
		boolean finishedSuccessfully=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final String[] processLiteralArg(String arg) {
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
	
	public static final String processLiteralTerm(String b) {
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
	
	/** Current available memory */
	private static final long freeMemory(){
		Runtime rt=Runtime.getRuntime();
		return rt.freeMemory();
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
	
	final long rcomp(long kmer, int len){
		return amino ? kmer : AminoAcid.reverseComplementBinaryFast(kmer, len);
	}
	
	final boolean passesSpeed(long key){
		return speed<1 || ((key&Long.MAX_VALUE)%17)>=speed;
	}
	
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

//		int left=minPolyLeft>0 ? r.countLeft(c) : 0;
//		int right=minPolyRight>0 ? r.countRight(c) : 0;
//
//		if(left<minPolyLeft){left=0;}
//		if(right<minPolyRight){right=0;}

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
	
	private static int parseEnd(String s){
		if(s==null){return RIGHTLEFT;}
		if(s.equalsIgnoreCase("right") || s.equalsIgnoreCase("r")){return RIGHT;}
		if(s.equalsIgnoreCase("left") || s.equalsIgnoreCase("l")){return LEFT;}
		if(s.equalsIgnoreCase("rl") || s.equalsIgnoreCase("lr") || s.equalsIgnoreCase("both") || s.equalsIgnoreCase("b")){return RIGHTLEFT;}
		return Parse.parseBoolean(s) ? RIGHTLEFT : 0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	//TODO: Document
	boolean silent=false;
	boolean json=false;
	
	boolean swift=false; //https://issues.jgi.doe.gov/browse/AUTOQC-2193
	
	/** For calculating kmer cardinality in input */
	final CardinalityTracker loglogIn;
	/** For calculating kmer cardinality in output */
	final CardinalityTracker loglogOut;
	/** Requires (and sets) cardinality tracking.  This is for input kmers. */
	String khistIn=null;
	/** Requires (and sets) cardinality tracking.  This is for output kmers. */
	String khistOut=null;
	
	/** Has this class encountered errors while processing? */
	public boolean errorState=false;
	
	/** Fraction of available memory preallocated to arrays */
	private double preallocFraction=1.0;
	/** Initial size of data structures */
	private int initialSize=-1;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	final AbstractKmerTable[] keySets;
	/** A scaffold's name is stored at scaffoldNames.get(id).
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	final ArrayList<String> scaffoldNames=new ArrayList<String>();
	/** Names of reference files (refNames[0] is valid). */
	private ArrayList<String> refNames=new ArrayList<String>();
	private final ArrayList<String> altRefNames=new ArrayList<String>();
	/** Number of scaffolds per reference. */
	private int[] refScafCounts;
	/** scaffoldCounts[id] stores the number of reads with kmer matches to that scaffold */
	AtomicLongArray scaffoldReadCounts;
	/** scaffoldBaseCounts[id] stores the number of bases with kmer matches to that scaffold */
	AtomicLongArray scaffoldBaseCounts;
	/** Set to false to force threads to share atomic counter arrays. */
	private boolean ALLOW_LOCAL_ARRAYS=true;
	/** scaffoldLengths[id] stores the length of that scaffold */
	private IntList scaffoldLengths=new IntList();
	/** hitCounts[x] stores the number of reads with exactly x kmer matches */
	long[] hitCounts;
	/** Array of reference files from which to load kmers */
	private String[] ref=null;
	/** Alternate reference to be used if main reference has no kmers */
	private String[] altref=null;
	/** Array of literal strings from which to load kmers */
	private String[] literal=null;
	/** Optional reference for sam file */
	private String samref=null;

	/** Input reads */
	private String in1=null, in2=null;
	/** Input FileFormats */
	private final FileFormat ffin1, ffin2;
	/** Input qual files */
	private String qfin1=null, qfin2=null;
	/** Output qual files */
	private String qfout1=null, qfout2=null;
	/** Output reads (unmatched and at least minlen) */
	private String out1=null, out2=null;
	/** Output reads (matched or shorter than minlen) */
	private String outb1=null, outb2=null;
	/** Output FileFormats */
	private final FileFormat ffout1, ffout2, ffoutb1, ffoutb2, ffouts;
	/** Output reads whose mate was discarded */
	private String outsingle=null;
	/** Statistics output files */
	private String outstats=null, outrqc=null, outrpkm=null, outrefstats=null, polymerStatsFile=null;
	@Deprecated
	/** duk-style statistics */
	private String outduk=null;
	
	final boolean tossJunk;
	
	/** Dump kmers here. */
	private String dump=null;

	/** Quit after this many bases written to outm */
	long maxBasesOutm=-1;
	/** Quit after this many bases written to outu */
	long maxBasesOutu=-1;
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	private long maxReads=-1;
	/** Process this fraction of input reads. */
	private float samplerate=1f;
	/** Set samplerate seed to this value. */
	private long sampleseed=-1;
	
	/** Output reads in input order.  May reduce speed. */
	private final boolean ordered;
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
	private String varFile=null;
	private String vcfFile=null;
	private VarMap varMap=null;
//	private HashSet<VarKey> varKeySet=null;
	private ScafMap scafMap=null;
	private boolean fixVariants=false;
	private boolean unfixVariants=true;
	
	/** Optional file for quality score recalibration */
	private String samFile=null;
	
	/** Filter reads with unsupported substitutions */
	private boolean filterVars=false;
	/** Maximum allowed unsupported substitutions in a read */
	private int maxBadSubs=2;
	/** Maximum variant depth for a variant to be considered unsupported */
	private int maxBadSubAlleleDepth=1;
	/** Minimum read depth for a variant to be considered unsupported */
	private int minBadSubReadDepth=2;
	//TODO
	private int minBadSubEDist=0;
	//TODO
	private float maxBadAlleleFraction=0;
	
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
	
	long readsIn=0;
	long basesIn=0;
	long readsOut=0;
	long basesOut=0;
	
	long readsQTrimmed=0;
	long basesQTrimmed=0;
	long readsFTrimmed=0;
	long basesFTrimmed=0;
	long readsQFiltered=0;
	long basesQFiltered=0;
	long readsEFiltered=0;
	long basesEFiltered=0;
	long readsNFiltered=0;
	long basesNFiltered=0;

	long readsPolyTrimmed=0;
	long basesPolyTrimmed=0;
	
	long readsKTrimmed=0;
	long basesKTrimmed=0;
	long readsKFiltered=0;
	long basesKFiltered=0;
	
	long badGcReads;
	long badGcBases;

	long badHeaderReads=0;
	long badHeaderBases=0;

	long readsTrimmedByOverlap;
	long basesTrimmedByOverlap;
	
	long readsTrimmedBySwift;
	long basesTrimmedBySwift;
	
	long refReads=0;
	long refBases=0;
	long refKmers=0;
	
//	public long modsum=0; //123
	
	long storedKmers=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Homopolymers         ----------------*/
	/*--------------------------------------------------------------*/
	
	boolean countPolymers=false;
	byte polymerChar1=-1;
	byte polymerChar2=-1;
	/** Minimum length to consider a homopolymer, for the purpose of statistics */
	int polymerLength=20;
	
	/** Tracks homopolymer statistics */
	PolymerTracker pTracker;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Don't look for kmers in read 1 */
	final boolean skipR1;
	/** Don't look for kmers in read 2 */
	final boolean skipR2;
	/** Correct errors via read overlap */
	final boolean ecc;
	/** True if a ReadStats object is being used for collecting data */
	final boolean makeReadStats;
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	final boolean rcomp;
	/** Don't allow a read 'N' to match a reference 'A'.
	 * Reduces sensitivity when hdist>0 or edist>0.  Default: false. */
	final boolean forbidNs;
	/** AND bitmask with 0's at the middle base */
	final long middleMask;
	/** Use HashForest data structure */
	private final boolean useForest;
	/** Use KmerTable data structure */
	private final boolean useTable;
	/** Use HashArray data structure (default) */
	private final boolean useArray;
	
	/** Normal kmer length */
	final int k;
	/** k-1; used in some expressions */
	final int k2;
	/** Emulated kmer greater than k */
	final int kbig;
	/** Effective kmer size */
	final int keff;
	/** Shortest kmer to use for trimming */
	final int mink;
	/** A read may contain up to this many kmers before being considered a match.  Default: 0 */
	final int maxBadKmers0;
	/** A read must share at least this fraction of its kmers to be considered a match.  Default: 0 */
	final float minKmerFraction;
	/** Reference kmers must cover at least this fraction of read bases to be considered a match.  Default: 0 */
	final float minCoveredFraction;
	
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
	private final float trimE;
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
	final float minBaseFrequency;
	/** Throw away reads shorter than this after trimming.  Default: 10 */
	final int minReadLength;
	/** Throw away reads longer than this after trimming.  Default: Integer.MAX_VALUE */
	final int maxReadLength;
	/** Toss reads shorter than this fraction of initial length, after trimming */
	final float minLenFraction;
	/** Filter reads by whether or not they have matching kmers */
	private final boolean kfilter;
	/** Trim matching kmers and all bases to the left */
	final boolean ktrimLeft;
	/** Trim matching kmers and all bases to the right */
	boolean ktrimRight;
	/** Don't trim, but replace matching kmers with a symbol (default N) */
	final boolean ktrimN;
	/** Exclude kmer itself when ktrimming */
	final boolean ktrimExclusive;
	/** Split into two reads around the kmer */
	final boolean ksplit;
	/** Replace bases covered by matched kmers with this symbol */
	final byte trimSymbol;
	/** Convert kmer-masked bases to lowercase */
	final boolean kmaskLowercase;
	/** Only mask fully-covered bases **/
	final boolean kmaskFullyCovered;
	/** Output over-trimmed reads to outbad (outmatch).  If false, they are discarded. */
	final boolean addTrimmedToBad;
	/** Find the sequence that shares the most kmer matches when filtering. */
	final boolean findBestMatch;
	/** Trim pairs to the same length, when adapter-trimming */
	final boolean trimPairsEvenly;
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
	final int restrictLeft;
	/** If positive, only look for kmer matches the rightmost X bases */
	final int restrictRight;
	
	/** Skip this many initial input reads */
	private final long skipreads;

	/** Pairs go to outbad if either of them is bad, as opposed to requiring both to be bad.
	 * Default: true. */
	final boolean removePairsIfEitherBad;

	/** Rather than discarding, trim failures to 1bp.
	 * Default: false. */
	final boolean trimFailuresTo1bp;
	
	/** Print only statistics for scaffolds that matched at least one read
	 * Default: true. */
	private final boolean printNonZeroOnly;
	
	/** Rename reads to indicate what they matched.
	 * Default: false. */
	final boolean rename;
	/** Use names of reference files instead of scaffolds.
	 * Default: false. */
	private final boolean useRefNames;
	
	/** Fraction of kmers to skip, 0 to 16 out of 17 */
	final int speed;
	
	/** Skip this many kmers when examining the read.  Default 1.
	 * 1 means every kmer is used, 2 means every other, etc. */
	final int qSkip;
	
	/** noAccel is true if speed and qSkip are disabled, accel is the opposite. */
	final boolean noAccel;
	private final boolean accel;
	
	private boolean pairedToSingle=false;
	
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
	final boolean trimByOverlap;
	final boolean useQualityForOverlap;
	final boolean strictOverlap;
	
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
	final SideChannel3 sidechannel;
	
	/*--------------------------------------------------------------*/
	/*----------------        Histogram Flags       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Generate histograms from the reads before rather than after processing;
	 * default is true.  Khist is handled independently.
	 */
	final boolean histogramsBeforeProcessing;
	
	final boolean MAKE_QUALITY_ACCURACY;
	final boolean MAKE_QUALITY_HISTOGRAM;
	final boolean MAKE_MATCH_HISTOGRAM;
	final boolean MAKE_BASE_HISTOGRAM;
	
	final boolean MAKE_EHIST;
	final boolean MAKE_INDELHIST;
	final boolean MAKE_LHIST;
	final boolean MAKE_GCHIST;
	final boolean MAKE_ENTROPYHIST;
	final boolean MAKE_IDHIST;
	
	final boolean MAKE_IHIST;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of tables (and threads, during loading) */
	private static final int WAYS=7; //123
	/** Default initial size of data structures */
	private static final int initialSizeDefault=128000;
	/** Verbose messages */
	public static final boolean verbose=false; //123

	/** Ends for some operations like entropytrim; could be migrated over to other operations */ 
	private static final int RIGHT=1, LEFT=2, RIGHTLEFT=3;
	
	/** Number of reads output in the last run */
	public static long lastReadsOut;
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
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
	/** Indicates end of input stream */
	static final ArrayList<Read> POISON=new ArrayList<Read>(0);
	/** Number of columns for statistics output, 3 or 5 */
	public static int STATS_COLUMNS=3;
	/** Release memory used by kmer storage after processing reads */
	public static boolean RELEASE_TABLES=true;
	/** Max value of hitCount array */
	public static final int HITCOUNT_LEN=1000;
	/** Make unambiguous copies of ref sequences with ambiguous bases */
	public static boolean REPLICATE_AMBIGUOUS=false;
	
	/** Stores some data for statistics when running RQCFilter; not used otherwise. */
	public static HashMap<String, Long> RQC_MAP=null;
	
}
