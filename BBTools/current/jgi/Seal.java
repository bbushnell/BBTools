package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLongArray;

import barcode.BarcodeMappingStats;
import cardinality.CardinalityTracker;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import kmer.AbstractKmerTable;
import kmer.ScheduleMaker;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ArrayListSet;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.MultiCros;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;
import structures.IntList;
import structures.IntList3;
import structures.ListNum;
import structures.StringCount;
import tax.GiToTaxid;
import tax.TaxNode;
import tax.TaxTree;
import tracker.ReadStats;

/**
 * SEAL: Sequence Expression AnaLyzer
 * Derived from BBDuk.
 * Allows multiple values stored per kmer.
 * Intended for RNA-seq, coverage, and other reads-per-sequence quantification.
 * Also performs binning.
 * @author Brian Bushnell
 * @date November 10, 2014
 *
 */
public class Seal {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Create a new Seal instance
		Seal x=new Seal(args);
		
		///And run it
		x.process();
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Seal(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), true);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		SamLine.SET_FROM_OK=true;
		IntList3.defaultMode=IntList3.ASCENDING;
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		/* Initialize local variables with defaults */
		boolean rcomp_=true;
		boolean forbidNs_=false;
		boolean prealloc_=false;
		boolean useCountvector_=false;
		int tableType_=AbstractKmerTable.ARRAYHF;//Fine as long as numbers are assigned in order per thread; otherwise use ARRAYH.
		int k_=31;
		int ways_=-1; //Currently disabled
		int minKmerHits_=1;
		float minKmerFraction_=0;
		long skipreads_=0;
		
		Parser parser=new Parser();
		parser.trimq=6;
		parser.minAvgQuality=0;
		parser.minReadLength=10;
		parser.maxReadLength=Integer.MAX_VALUE;
		parser.minLenFraction=0f;
		parser.requireBothBad=false;
		parser.maxNs=-1;
		boolean ordered_=false;
		int restrictLeft_=0, restrictRight_=0, speed_=0, qSkip_=1;
		int ambigMode_=AMBIG_RANDOM;
		int matchMode_=MATCH_ALL;
		boolean keepPairsTogether_=true;
		boolean printNonZeroOnly_=true;
		boolean rename_=false, useRefNames_=false;
		boolean ecc_=false;
		int clearzone_=0;
		float clearzoneFraction_=0;
		
		{//Add a fake first scaffold so that the true first scaffold has ID 1, not 0.
			scaffoldNames.add("");
			scaffoldLengths.add(0);
			scaffoldKmers.add(0);
			scaffolds.add(null);
		}
		
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
			}else if(parser.parseCardinality(arg, a, b)){
				//do nothing
			}else if(parser.parseCommon(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfin1=b;
			}else if(a.equals("qfin2")){
				qfin2=b;
			}else if(a.equals("out") || a.equals("out1") || a.equals("outm") || a.equals("outm1") || a.equals("outmatched") || a.equals("outmatched1")){
				outm1=b;
			}else if(a.equals("out2") || a.equals("outm") || a.equals("outm2") || a.equals("outmatched") || a.equals("outmatched2")){
				outm2=b;
			}else if(a.equals("outu") || a.equals("outu1") || a.equals("outunmatched") || a.equals("outunmatched1")){
				outu1=b;
			}else if(a.equals("outu2") || a.equals("outunmatched") || a.equals("outunmatched2")){
				outu2=b;
			}else if(a.equals("outpattern") || a.equals("pattern") || a.equals("basename")){
				outpattern=b;
			}else if(a.equals("stats") || a.equals("scafstats")){
				outstats=b;
			}else if(a.equals("refstats")){
				outrefstats=b;
			}else if(a.equals("rpkm") || a.equals("fpkm") || a.equals("cov") || a.equals("coverage") || a.equals("covstats")){
				outrpkm=b;
			}else if(a.equals("tax") || a.equals("taxa") || a.equals("outtax")){
				outtax=b;
			}else if(a.equals("outbarcodes") || a.equals("outbc")){
				outbarcodes=b;
			}else if(a.equals("ref")){
				ref.clear();
				String[] b2=(b==null) ? null : (new File(b).exists() ? new String[] {b} : b.split(","));
				for(String b3 : b2){ref.add(b3);}
			}else if(a.equals("literal")){
				literal=(b==null) ? null : b.split(",");
//				assert(false) : b+", "+Arrays.toString(literal);
			}else if(a.equals("forest")){
				boolean x=Parse.parseBoolean(b);
				if(x){tableType_=AbstractKmerTable.FOREST2D;}
			}else if(a.equals("array") || a.equals("array2")){
				boolean x=Parse.parseBoolean(b);
				if(x){tableType_=AbstractKmerTable.ARRAY2D;}
			}else if(a.equals("array1")){
				boolean x=Parse.parseBoolean(b);
				if(x){tableType_=AbstractKmerTable.ARRAY1D;}
			}else if(a.equals("arrayh") || a.equals("hybrid")){
				boolean x=Parse.parseBoolean(b);
				if(x){tableType_=AbstractKmerTable.ARRAYH;}
			}else if(a.equals("arrayhf") || a.equals("hybridfast")){
				boolean x=Parse.parseBoolean(b);
				if(x){tableType_=AbstractKmerTable.ARRAYHF;}
			}else if(a.equals("ways")){
				ways_=Integer.parseInt(b);
			}else if(a.equals("ordered") || a.equals("ord")){
				ordered_=Parse.parseBoolean(b);
				outstream.println("Set ORDERED to "+ordered_);
			}else if(a.equals("k")){
				assert(b!=null) : "\nThe k key needs an integer value greater than 0, such as k=27\n";
				k_=Integer.parseInt(b);
				assert(k_>0 && k_<32) : "k must be at least 1; default is 31.";
			}else if(a.equals("hdist") || a.equals("hammingdistance")){
				hammingDistance=Integer.parseInt(b);
				assert(hammingDistance>=0 && hammingDistance<4) : "hamming distance must be between 0 and 3; default is 0.";
			}else if(a.equals("qhdist") || a.equals("queryhammingdistance")){
				qHammingDistance=Integer.parseInt(b);
				assert(qHammingDistance>=0 && qHammingDistance<4) : "hamming distance must be between 0 and 3; default is 0.";
			}else if(a.equals("edits") || a.equals("edist") || a.equals("editdistance")){
				editDistance=Integer.parseInt(b);
				assert(editDistance>=0 && editDistance<=3) : "edit distance must be between 0 and 2; default is 0.";
			}else if(a.equals("skip") || a.equals("refskip") || a.equals("rskip")){
				refSkip=Integer.parseInt(b);
			}else if(a.equals("qskip")){
				qSkip_=Integer.parseInt(b);
			}else if(a.equals("speed")){
				speed_=Integer.parseInt(b);
				assert(speed_>=0 && speed_<=16) : "Speed range is 0 to 16.  Value: "+speed_;
			}else if(a.equals("skipreads")){
				skipreads_=Parse.parseKMG(b);
			}else if(a.equals("minkmerhits") || a.equals("minhits") || a.equals("mh") || a.equals("mkh")){
				minKmerHits_=Integer.parseInt(b);
			}else if(a.equals("minkmerfraction") || a.equals("minfraction") || a.equals("mkf")){
				minKmerFraction_=Float.parseFloat(b);
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
			}else if(a.equals("ambiguous") || a.equals("ambig")){
				if(b==null){
					throw new RuntimeException(arg);
				}else if(b.equalsIgnoreCase("keep") || b.equalsIgnoreCase("best") || b.equalsIgnoreCase("first")){
					ambigMode_=AMBIG_FIRST;
				}else if(b.equalsIgnoreCase("all")){
					ambigMode_=AMBIG_ALL;
				}else if(b.equalsIgnoreCase("random") || b.equalsIgnoreCase("rand")){
					ambigMode_=AMBIG_RANDOM;
				}else if(b.equalsIgnoreCase("toss") || b.equalsIgnoreCase("discard") || b.equalsIgnoreCase("remove")){
					ambigMode_=AMBIG_TOSS;
				}else{
					throw new RuntimeException(arg);
				}
			}else if(a.equals("match") || a.equals("mode")){
				if(b==null){
					throw new RuntimeException(arg);
				}else if(b.equalsIgnoreCase("all") || b.equalsIgnoreCase("best")){
					matchMode_=MATCH_ALL;
				}else if(b.equalsIgnoreCase("first")){
					matchMode_=MATCH_FIRST;
				}else if(b.equalsIgnoreCase("unique") || b.equalsIgnoreCase("firstunique")){
					matchMode_=MATCH_UNIQUE;
				}else{
					throw new RuntimeException(arg);
				}
			}else if(a.equals("findbestmatch") || a.equals("fbm")){
				matchMode_=(Parse.parseBoolean(b) ? MATCH_ALL : MATCH_FIRST);
			}else if(a.equals("firstuniquematch") || a.equals("fum")){
				if(Parse.parseBoolean(b)){matchMode_=MATCH_UNIQUE;}
			}else if(a.equals("keeppairstogether") || a.equals("kpt")){
				keepPairsTogether_=Parse.parseBoolean(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				printNonZeroOnly_=Parse.parseBoolean(b);
			}else if(a.equals("rename")){
				rename_=Parse.parseBoolean(b);
			}else if(a.equals("renameduplicates")){
				RENAME_DUPLICATES=Parse.parseBoolean(b);
			}else if(a.equals("genesets") || a.equals("makesets") || a.equals("makegenesets")){
				MAKE_GENE_SETS=Parse.parseBoolean(b);
			}else if(a.equals("refnames") || a.equals("userefnames")){
				useRefNames_=Parse.parseBoolean(b);
			}else if(a.equals("addcount")){
				addCount=Parse.parseBoolean(b);
			}else if(a.equals("tophitonly")){
				topHitOnly=Parse.parseBoolean(b);
			}else if(a.equals("addunknown")){
				addUnknown=Parse.parseBoolean(b);
			}else if(a.equals("initialsize")){
				initialSize=Parse.parseIntKMG(b);
			}else if(a.equals("dump")){
				dump=b;
			}else if(a.equals("countvector")){
				useCountvector_=Parse.parseBoolean(b);
			}else if(a.equals("ecco") || a.equals("ecc")){
				ecc_=Parse.parseBoolean(b);
			}else if(a.equals("copyundefined") || a.equals("cu")){
				REPLICATE_AMBIGUOUS=Parse.parseBoolean(b);
			}else if(a.equals("bbsplit")){
				BBSPLIT_STYLE=Parse.parseBoolean(b);
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
				if("auto".equalsIgnoreCase(b)){giTableFile=TaxTree.defaultTableFile();}
			}else if(a.equals("taxnames") || a.equals("taxname")){
				taxNameFile=b;
			}else if(a.equals("taxnodes") || a.equals("taxnode")){
				taxNodeFile=b;
			}else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
			}else if(a.equals("mincount")){
				taxNodeCountLimit=Long.parseLong(b);
			}else if(a.equals("maxnodes")){
				taxNodeNumberLimit=Integer.parseInt(b);
			}else if(a.equals("minlevel")){
				taxNodeMinLevel=TaxTree.parseLevel(b);
			}else if(a.equals("maxlevel")){
				taxNodeMaxLevel=TaxTree.parseLevel(b);
			}else if(a.equals("clearzone") || a.equals("cz")){
				if(b.indexOf('.')>=0) {
					clearzoneFraction_=Float.parseFloat(b);
				}else {
					clearzone_=Integer.parseInt(b);
				}
			}else if(a.equals("clearzonefraction") || a.equals("czf")){
				clearzoneFraction_=Float.parseFloat(b);
			}
			
			else if(a.equals("processcontainedref")){
				processContainedRef=Parse.parseBoolean(b);
			}else if(a.equals("storerefbases")){
				storeRefBases=Parse.parseBoolean(b);
			}
			
			else if(a.equals("trackbarcodes")){
				trackBarcodes=Parse.parseBoolean(b);
			}
			
			else if(b==null && new File(arg).exists()){
				ref.add(arg);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}		
		
		if("auto".equals(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if(ref.isEmpty()){ref=null;}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			samplerate=parser.samplerate;
			sampleseed=parser.sampleseed;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			forceTrimModulo=parser.forceTrimModulo;
			forceTrimLeft=parser.forceTrimLeft;
			forceTrimRight=parser.forceTrimRight;
			forceTrimRight2=parser.forceTrimRight2;
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
			trimE=parser.trimE();
			minLenFraction=parser.minLenFraction;
			minAvgQuality=parser.minAvgQuality;
			minAvgQualityBases=parser.minAvgQualityBases;
			chastityFilter=parser.chastityFilter;
			minReadLength=parser.minReadLength;
			maxReadLength=parser.maxReadLength;
			maxNs=parser.maxNs;
			minConsecutiveBases=parser.minConsecutiveBases;
//			minGC=parser.minGC;
//			maxGC=parser.maxGC;
//			filterGC=parser.filterGC;
//			minTrimLength=(parser.minTrimLength>=0 ? parser.minTrimLength : minTrimLength);
//			requireBothBad=parser.requireBothBad;
			removePairsIfEitherBad=!parser.requireBothBad;

			loglog=(parser.loglog ? CardinalityTracker.makeTracker(parser) : null);
			loglogOut=(parser.loglogOut ? CardinalityTracker.makeTracker(parser) : null);
			
			THREADS=Shared.threads();
		}
		
		refNames.add(null);//Fake first reference also
		
		if(ref!=null){
			ArrayList<String> temp=new ArrayList<String>();
			for(String s : ref){
				if(s==null){
					assert(false) : "Null reference file.";
				}else if(new File(s).exists()){
					Tools.getFileOrFiles(s, temp, true, false, false, false);
				}else{
					String fname=null;
					if("phix".equalsIgnoreCase(s)){
						fname=Data.findPath("?phix2.fa.gz");
					}else if("polya".equalsIgnoreCase(s) || "polyt".equalsIgnoreCase(s)){
						fname=Data.findPath("?polyA.fa.gz");
					}else if("lambda".equalsIgnoreCase(s)){
						fname=Data.findPath("?lambda.fa.gz");
					}else if("kapa".equalsIgnoreCase(s)){
						fname=Data.findPath("?kapatags.L40.fa");
					}else if("pjet".equalsIgnoreCase(s)){
						fname=Data.findPath("?pJET1.2.fa");
					}else if("mtst".equalsIgnoreCase(s)){
						fname=Data.findPath("?mtst.fa");
					}else if("adapters".equalsIgnoreCase(s)){
						fname=Data.findPath("?adapters.fa");
					}else if("truseq".equalsIgnoreCase(s)){
						fname=Data.findPath("?truseq.fa.gz");
					}else if("nextera".equalsIgnoreCase(s)){
						fname=Data.findPath("?nextera.fa.gz");
					}else if("artifacts".equalsIgnoreCase(s)){
						fname=Data.findPath("?sequencing_artifacts.fa.gz");
					}else if("crisprs".equalsIgnoreCase(s)){
						s=Data.findPath("?crisprs.fa.gz");
					}else{
						assert(false) : "Can't find reference file "+s;
					}
					temp.add(fname);
				}
			}
			ref=temp;
			if(ref.size()<1){ref=null;}
			refNames.addAll(temp);
		}
		
//		if(ref!=null){
//			ArrayList<String> temp=new ArrayList<String>();
//			for(String s : ref){
//				Tools.getFileOrFiles(s, temp, true, false, false, false);
//			}
//			ref=temp.toArray(new String[0]);
//			if(ref.length<1){ref=null;}
//			refNames.addAll(temp);
//		}
		if(literal!=null){refNames.add("literal");}
		refScafCounts=new int[refNames.size()];
		
		if(prealloc_){
			outstream.println("Note - if this program runs out of memory, please disable the prealloc flag.");
		}
		
		/* Set final variables; post-process and validate argument combinations */
		
		tableType=tableType_;
		hammingDistance=Tools.max(editDistance, hammingDistance);
		refSkip=Tools.max(0, refSkip);
		rcomp=rcomp_;
		forbidNs=(forbidNs_ || hammingDistance<1);
		skipreads=skipreads_;
		ordered=ordered_;
		restrictLeft=Tools.max(restrictLeft_, 0);
		restrictRight=Tools.max(restrictRight_, 0);
		ambigMode=ambigMode_;
		matchMode=matchMode_;
		keepPairsTogether=keepPairsTogether_;
		printNonZeroOnly=printNonZeroOnly_;
		rename=rename_;
		useRefNames=useRefNames_;
		speed=speed_;
		qSkip=qSkip_;
		noAccel=(speed<1 && qSkip<2);
		clearzone=clearzone_;
		clearzoneFraction=clearzoneFraction_;
		parsecustom=FASTQ.PARSE_CUSTOM;
		ecc=ecc_;
		
		USE_TAXTREE=(taxNameFile!=null || taxNodeFile!=null || outtax!=null || taxTreeFile!=null);
		USE_COUNTVECTOR=useCountvector_;
		
		if(outbarcodes!=null){trackBarcodes=true;}
		if(trackBarcodes){
			barcodeData=new BarcodeMappingStats();
		}
		
		if((speed>0 && qSkip>1) || (qSkip>1 && refSkip>1) || (speed>0 && refSkip>1)){
			outstream.println("WARNING: It is not recommended to use more than one of qskip, speed, and rskip together.");
			outstream.println("qskip="+qSkip+", speed="+speed+", rskip="+refSkip);
		}
		
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
				final int factor=(tableType==AbstractKmerTable.ARRAY1D ? 12 : tableType==AbstractKmerTable.ARRAYH ? 22 : tableType==AbstractKmerTable.ARRAYHF ? 22 : 27);
				final long memOverWays=tableMemory/(factor*WAYS);
				final double mem2=(prealloc_ ? preallocFraction : 1)*tableMemory;
				initialSize=(prealloc_ || memOverWays<initialSizeDefault ? (int)Tools.min(2142000000, (long)(mem2/(factor*WAYS))) : initialSizeDefault);
				if(initialSize!=initialSizeDefault){
					outstream.println("Initial size set to "+initialSize);
				}
			}
		}
		
		k=k_;
		k2=k-1;
		if(maskMiddle) {
			midMaskLen=(midMaskLen>0 ? midMaskLen : 2-(k&1));
		}else{
			midMaskLen=0;
		}
		minKmerHits=minKmerHits_;
		minKmerFraction=Tools.max(minKmerFraction_, 0);
		assert(minKmerHits>=1) : "minKmerHits must be at least 1; value="+minKmerHits;
		assert(minKmerFraction<=1) : "minKmerFraction must range from 0 to 1; value="+minKmerFraction;
		
		kfilter=(ref!=null || literal!=null);
		assert(kfilter==false || (k>0 && k<32)) : "K must range from 1 to 31.";
		
		if(maskMiddle){
			final int bitsPerBase=2;
			assert(k>midMaskLen+1);
			int bits=midMaskLen*bitsPerBase;
//			int shift=(k-maskMiddle)&(~1);//Equivalent to (x/2)*2
			int shift=((k-midMaskLen)/2)*bitsPerBase; //old behavior before moving to variable width can be restored with +1: "((k-maskMiddle+1)/2)*bitsPerBase"
			middleMask=~((~((-1L)<<bits))<<shift);
//			mask<<=(shift);
		}else{
			middleMask=-1L;
		}
//		middleMask=maskMiddle ? ~(3L<<(2*(k/2))) : -1L;
		
		
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
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		if(qfin1!=null && qfin1.contains("#") && in2!=null && !new File(qfin1).exists()){
			int pound=qfin1.lastIndexOf('#');
			String a=qfin1.substring(0, pound);
			String b=qfin1.substring(pound+1);
			qfin1=a+1+b;
			qfin2=a+2+b;
		}
		
		if(outu1!=null && outu1.contains("#")){
			int pound=outu1.lastIndexOf('#');
			String a=outu1.substring(0, pound);
			String b=outu1.substring(pound+1);
			outu1=a+1+b;
			outu2=a+2+b;
		}
		
		if(outm1!=null && outm1.contains("#")){
			int pound=outm1.lastIndexOf('#');
			String a=outm1.substring(0, pound);
			String b=outm1.substring(pound+1);
			outm1=a+1+b;
			outm2=a+2+b;
		}
		
		if((outu2!=null || outm2!=null) && (in1!=null && in2==null)){
			if(!FASTQ.FORCE_INTERLEAVED){outstream.println("Forcing interleaved input because paired output was specified for a single input file.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=true;
		}

		if(!Tools.testOutputFiles(overwrite, append, false, outu1, outu2, outm1, outm2, outpattern, outstats, outrpkm, outrefstats)){
			throw new RuntimeException("\nCan't write to some output files; overwrite="+overwrite+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2, qfin1, qfin2, taxNameFile, taxNodeFile, giTableFile, taxTreeFile)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(!Tools.testInputFiles(true, true, ref)){
			throw new RuntimeException("\nCan't read to some reference files.\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, qfin1, qfin2, outu1, outu2, outm1, outm2, outpattern, outstats, outrpkm, outrefstats)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		assert(THREADS>0) : "THREADS must be greater than 0.";

		assert(in1==null || in1.toLowerCase().startsWith("stdin") || in1.toLowerCase().startsWith("standardin") || new File(in1).exists()) : "Can't find "+in1;
		assert(in2==null || in2.toLowerCase().startsWith("stdin") || in2.toLowerCase().startsWith("standardin") || new File(in2).exists()) : "Can't find "+in2;
		
		if(ref==null && literal==null){
			outstream.println("ERROR: No reference sequences specified.  Use the -da flag to run anyway.");
			assert(false) : "Please specify a reference.";
		}
				
		if(ref!=null){
			for(String s0 : ref){
				assert(s0!=null) : "Specified a null reference.";
				String s=s0.toLowerCase();
				assert(s==null || s.startsWith("stdin") || s.startsWith("standardin") || new File(s0).exists()) : "Can't find "+s0;
			}
		}
		
		if(MAKE_GENE_SETS) {
			setMap=new ConcurrentHashMap<String, StringCount>();
		}
		
		//Initialize tables
		ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, 14, prealloc_, (prealloc_ ? preallocFraction : 0.9));
		int[] schedule=scheduleMaker.makeSchedule();
		keySets=AbstractKmerTable.preallocate(WAYS, tableType, schedule, -1L);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Primary driver function of class.
	 */
	public void process(){
		
		/* Check for output file collisions */
		if(!Tools.testOutputFiles(overwrite, append, false, outu1, outu2, outm1, outm2, outstats, outrpkm, outrefstats)){
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
		lastReadsOut=readsUnmatched;
		
		
		if(showSpeed){
			outstream.println();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		}
		
		if(outstream!=System.err && outstream!=System.out){outstream.close();}
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException("Seal terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	public void process2(long startTime){
		
		/* Start phase timer */
		Timer t=new Timer();

		if(DISPLAY_PROGRESS){
			outstream.println("Initial:");
			Shared.printMemory();
			outstream.println();
		}
		
		/* Fill tables with reference kmers */
		{
			final boolean oldTI=FASTQ.TEST_INTERLEAVED; //TODO: This needs to be changed to a non-static field, or somehow 'read mode' and 'ref mode' need to be distinguished.
			final boolean oldFI=FASTQ.FORCE_INTERLEAVED;
			final boolean oldSplit=FastaReadInputStream.SPLIT_READS;
			final int oldML=FastaReadInputStream.MIN_READ_LEN;
			
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
			FastaReadInputStream.SPLIT_READS=false;
			FastaReadInputStream.MIN_READ_LEN=1;
			
			storedKmers=spawnLoadThreads();
			
			FASTQ.TEST_INTERLEAVED=oldTI;
			FASTQ.FORCE_INTERLEAVED=oldFI;
			FastaReadInputStream.SPLIT_READS=oldSplit;
			FastaReadInputStream.MIN_READ_LEN=oldML;
			
//			if(useRefNames){toRefNames();}
			t.stop();
		}
		
		/* Check memory */
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
		
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=THREADS<4;
		
		/* Do kmer matching of input reads */
		spawnProcessThreads(t);
		
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		/* Unload kmers to save memory */
		if(RELEASE_TABLES){
			unloadKmers();
		}
		
		if(USE_TAXTREE){
			if(giTableFile!=null){loadGiToTaxid();}
			if(USE_TAXTREE){tree=TaxTree.loadTaxTree(taxTreeFile, taxNameFile, taxNodeFile, null, DISPLAY_PROGRESS ? outstream : null, false, false);}
			addToTree();
		}
		
		/* Write statistics to files */
		writeStats();
		writeRPKM();
		writeBarcodeStats();
		addToRqcMap();
		
		if(!BBSPLIT_STYLE){
			writeRefStats();
		}else{
			writeRefStats_BBSplitStyle(readsIn);
		}
		writeTaxonomy();
		
		/* Unload sequence data to save memory */
		if(RELEASE_TABLES){
			unloadScaffolds();
			tree=null;
			GiToTaxid.unload();
		}
		
		outstream.println("\nInput:                  \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		
		if(ref!=null || literal!=null){
			outstream.println("Matched reads:          \t"+readsMatched+" reads ("+toPercent(readsMatched, readsIn)+") \t"+
					basesMatched+" bases ("+toPercent(basesMatched, basesIn)+")");
			outstream.println("Unmatched reads:        \t"+readsUnmatched+" reads ("+toPercent(readsUnmatched, readsIn)+") \t"+
					basesUnmatched+" bases ("+toPercent(basesUnmatched, basesIn)+")");
			outstream.flush();
		}
		if(qtrimLeft || qtrimRight){
			outstream.println("QTrimmed:               \t"+readsQTrimmed+" reads ("+toPercent(readsQTrimmed, readsIn)+") \t"+
					basesQTrimmed+" bases ("+toPercent(basesQTrimmed, basesIn)+")");
		}
		if(forceTrimLeft>0 || forceTrimRight>0 || forceTrimRight2>0 || forceTrimModulo>0){
			outstream.println("FTrimmed:               \t"+readsFTrimmed+" reads ("+toPercent(readsFTrimmed, readsIn)+") \t"+
					basesFTrimmed+" bases ("+toPercent(basesFTrimmed, basesIn)+")");
		}
		if(minAvgQuality>0 || maxNs>=0){
			outstream.println("Low quality discards:   \t"+readsQFiltered+" reads ("+toPercent(readsQFiltered, readsIn)+") \t"+
					basesQFiltered+" bases ("+toPercent(basesQFiltered, basesIn)+")");
		}
		if(loglog!=null){
			outstream.println("Unique "+loglog.k+"-mers:         \t"+loglog.cardinality());
		}
		if(loglogOut!=null){
			outstream.println("Unique "+loglogOut.k+"-mers out:     \t"+loglogOut.cardinality());
		}
		if(parsecustom){
			outstream.println();
			outstream.println("Correctly mapped:       \t"+correctReads+" reads ("+toPercent(correctReads, readsIn)+")");
			outstream.println("Incorrectly mapped:     \t"+incorrectReads+" reads ("+toPercent(incorrectReads, readsIn)+")");
		}
//		outstream.println("Result:                 \t"+readsMatched+" reads ("+toPercent(readsMatched*100.0/readsIn)+"%) \t"+
//				basesMatched+" bases ("+toPercent(basesMatched*100.0/basesIn)+"%)");
	}
	
	private static String toPercent(long numerator, long denominator){
		if(denominator<1){return "0.00%";}
		return Tools.format("%.2f%%",numerator*100.0/denominator);
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
		scaffoldFragCounts=null;
		scaffoldBaseCounts=null;
		scaffoldAmbigReadCounts=null;
		scaffoldLengths=null;
		scaffoldKmers=null;
		scaffolds=null;
	}
	
	private void addToRqcMap(){
		BBDuk.putRqc("inputReads", readsIn, false, false);
		BBDuk.putRqc("inputBases", basesIn, false, false);
		if(qtrimLeft || qtrimRight){
			BBDuk.putRqc("qtrimmedReads", readsQTrimmed, false, true);
			BBDuk.putRqc("qtrimmedBases", basesQTrimmed, false, true);
		}
		BBDuk.putRqc("qfilteredReads", readsQFiltered, false, true);
		BBDuk.putRqc("qfilteredBases", basesQFiltered, false, true);
		
		{//This is kind of a hack, to match BBDuk's syntax for RQCFilter.
			BBDuk.putRqc("kfilteredReads", readsMatched, false, true);
			BBDuk.putRqc("kfilteredBases", basesMatched, false, true);

			BBDuk.putRqc("outputReads", readsUnmatched, true, false);
			BBDuk.putRqc("outputBases", basesUnmatched, true, false);
		}
	}
	
	/**
	 * Write statistics about how many reads matched each reference scaffold.
	 */
	private void writeStats(){
		if(outstats==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outstats, overwrite, false, false);
		tsw.start();
		
		long rsum=0, bsum=0, asum=0;
		
		/* Create StringCount list of scaffold names and hitcounts */
		ArrayList<StringCount> list=new ArrayList<StringCount>();
		for(int i=1; i<scaffoldNames.size(); i++){
			final long num1=scaffoldReadCounts.get(i), num2=scaffoldBaseCounts.get(i), num3=scaffoldAmbigReadCounts.get(i);
			if(num1>0 || !printNonZeroOnly){
				rsum+=num1;
				bsum+=num2;
				asum+=num3;
				final String s=scaffoldNames.get(i);
				final int len=scaffoldLengths.get(i);
				final StringCount sn=new StringCount(s, len, num1, num2, num3);
				list.add(sn);
			}
		}
		Shared.sort(list);
		if(MAKE_GENE_SETS) {
			ArrayList<StringCount> list2=new ArrayList<StringCount>(setMap.values());
			Shared.sort(list2);
			list.addAll(list2);
		}
		final double rmult=100.0/(readsIn>0 ? readsIn : 1);
		final double bmult=100.0/(basesIn>0 ? basesIn : 1);
		
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		if(STATS_COLUMNS==3){
			tsw.print(Tools.format("#Total\t%d\n",readsIn));
//			tsw.print(Tools.format("#Matched\t%d\t%.5f%%\n",rsum,rmult*rsum)); //With ambig=all, gives over 100%
			tsw.print(Tools.format("#Matched\t%d\t%.5f%%\n",readsMatched,rmult*readsMatched));
			tsw.print("#Name\tReads\tReadsPct\n");
			for(int i=0; i<list.size(); i++){
				StringCount sn=list.get(i);
				tsw.print(Tools.format("%s\t%d\t%.5f%%\n",sn.name,sn.reads,(sn.reads*rmult)));
			}
		}else{
			tsw.print(Tools.format("#Total\t%d\t%d\n",readsIn,basesIn));
//			tsw.print(Tools.format("#Matched\t%d\t%.5f%%\n",rsum,rmult*rsum,bsum,bsum*bmult)); //With ambig=all, gives over 100%
			tsw.print(Tools.format("#Matched\t%d\t%.5f%%\t%d\n",readsMatched,rmult*readsMatched,basesMatched,basesMatched*bmult,asum));
			tsw.print("#Name\tReads\tReadsPct\tBases\tBasesPct\tAmbigReads\n");
			for(int i=0; i<list.size(); i++){
				StringCount sn=list.get(i);
				tsw.print(Tools.format("%s\t%d\t%.5f%%\t%d\t%.5f%%\t%d\n",sn.name,sn.reads,(sn.reads*rmult),sn.bases,(sn.bases*bmult),sn.ambigReads));
			}
		}
		tsw.poisonAndWait();
	}
	
	private void writeBarcodeStats(){
		if(barcodeData==null || outbarcodes==null) {return;}
		barcodeData.writeStats(outbarcodes, overwrite);
	}
	
	private void writeRPKM(){
		writeRPKM(outrpkm, in1, in2, readsIn, printNonZeroOnly,
				scaffoldNames, scaffoldLengths,
				scaffoldReadCounts, scaffoldFragCounts, scaffoldBaseCounts, scaffoldAmbigReadCounts);
	}
	
	/**
	 * Write RPKM statistics.
	 */
	public void writeRPKM(String out, String in1, String in2, long readsIn, boolean printNonZeroOnly,
			ArrayList<String> scaffoldNames, IntList scaffoldLengths, AtomicLongArray scaffoldReadCounts, 
			AtomicLongArray scaffoldFragCounts, AtomicLongArray scaffoldBaseCounts, AtomicLongArray scaffoldAmbigReadCounts){
		if(out==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(out, overwrite, false, false);
		tsw.start();
		
		/* Count mapped reads */
		long mappedReads=0;
		long mappedFrags=0;
		for(int i=0; i<scaffoldReadCounts.length(); i++){
			mappedReads+=scaffoldReadCounts.get(i);
			mappedFrags+=scaffoldFragCounts.get(i);
		}
		
		/* Print header */
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		tsw.print(Tools.format("#Reads\t%d\n",readsIn));
//		tsw.print(Tools.format("#Mapped\t%d\n",mappedReads));
		tsw.print(Tools.format("#Mapped\t%d\n",readsMatched));
		tsw.print(Tools.format("#RefSequences\t%d\n",Tools.max(0, scaffoldNames.size()-1)));
		tsw.print("#Name\tLength\tBases\tCoverage\tReads\tRPKM\tFrags\tFPKM\tAmbigReads\n");

		final float readMult=1000000000f/Tools.max(1, mappedReads);
		final float fragMult=1000000000f/Tools.max(1, mappedFrags);
		
		/* Print data */
		for(int i=1; i<scaffoldNames.size(); i++){
			final long reads=scaffoldReadCounts.get(i);
			final long frags=scaffoldFragCounts.get(i);
			final long bases=scaffoldBaseCounts.get(i);
			final long ambig=scaffoldAmbigReadCounts.get(i);
			final String s=scaffoldNames.get(i);
			final int len=scaffoldLengths.get(i);
			final double invlen=1.0/Tools.max(1, len);
			final double readMult2=readMult*invlen;
			final double fragMult2=fragMult*invlen;
			if(reads>0 || !printNonZeroOnly){
				tsw.print(Tools.format("%s\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\t%d\n",s,len,bases,bases*invlen,reads,reads*readMult2,frags,frags*fragMult2,ambig));
			}
		}
		if(MAKE_GENE_SETS) {
			ArrayList<StringCount> list=new ArrayList<StringCount>(setMap.values());
			Shared.sort(list);
			for(StringCount sc : list) {
				String s=sc.name;
				final long reads=sc.reads;
				final long frags=sc.reads;
				final long bases=sc.bases;
				final long ambig=sc.ambigReads;
				final int len=sc.length;
				final double invlen=1.0/Tools.max(1, len);
				final double readMult2=readMult*invlen;
				final double fragMult2=fragMult*invlen;
				if(reads>0 || !printNonZeroOnly){
					tsw.print(Tools.format("%s\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\t%d\n",s,len,bases,bases*invlen,reads,reads*readMult2,frags,frags*fragMult2,ambig));
				}
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
		final int totalScafs=scaffoldReadCounts.length();
		final int totalScafs2=(int)shared.Vector.sum(refScafCounts);
		assert(totalScafs<=totalScafs2+1) : totalScafs+", "+totalScafs2;
		long[] refReadCounts=new long[numRefs];
		long[] refFragCounts=new long[numRefs];
		long[] refBaseCounts=new long[numRefs];
		long[] refLengths=new long[numRefs];
		long[] ambigReadCounts=new long[numRefs];
		
		for(int refNum=1, scafNum=1; refNum<numRefs; refNum++){
			final int scafs=(useRefNames ? 1 : refScafCounts[refNum]);
			final int scafLimitThisRef=scafNum+scafs;
//			System.err.println("ref="+refNum+", scafs="+scafs+", scafnum="+scafNum+", limit="+scafLimitThisRef+", max="+totalScafs);
			assert(refScafCounts[refNum]>0) : refNum;
			assert(scafLimitThisRef<=scaffoldReadCounts.length()) : scafLimitThisRef+", "+scaffoldReadCounts.length()+", "+refNum+", "+refScafCounts[refNum]+", "+Arrays.toString(refScafCounts);
			while(scafNum<scafLimitThisRef){
				refReadCounts[refNum]+=scaffoldReadCounts.get(scafNum);
				refFragCounts[refNum]+=scaffoldFragCounts.get(scafNum);
				refBaseCounts[refNum]+=scaffoldBaseCounts.get(scafNum);
				refLengths[refNum]+=scaffoldLengths.get(scafNum);
				ambigReadCounts[refNum]+=scaffoldAmbigReadCounts.get(scafNum);
				scafNum++;
			}
		}
		
		/* Print header */
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		tsw.print(Tools.format("#Reads\t%d\n",readsIn));
		tsw.print(Tools.format("#Mapped\t%d\n",mapped));
		tsw.print(Tools.format("#References\t%d\n",refNames.size()-1));
		tsw.print("#Name\tLength\tScaffolds\tBases\tCoverage\tReads\tRPKM\tFrags\tFPKM\tAmbigReads\n");
		
		final float mult=1000000000f/Tools.max(1, mapped);
		
		/* Print data */
		for(int i=1; i<refNames.size(); i++){
			final long reads=refReadCounts[i];
			final long frags=refFragCounts[i];
			final long bases=refBaseCounts[i];
			final long len=refLengths[i];
			final long ambig=ambigReadCounts[i];
			final int scafs=refScafCounts[i];
			final String name=ReadWrite.stripToCore(refNames.get(i));
			final double invlen=1.0/Tools.max(1, len);
			final double mult2=mult*invlen;
			if(reads>0 || !printNonZeroOnly){
				tsw.print(Tools.format("%s\t%d\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\t%d\n",name,len,scafs,bases,bases*invlen,reads,reads*mult2,frags,frags*mult2,ambig));
			}
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Write statistics on a per-reference basis.
	 */
	private void writeRefStats_BBSplitStyle(long totalReads){
		if(outrefstats==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outrefstats, overwrite, false, false);
		tsw.start();
		
		final int numRefs=refNames.size();
		long[] refReadCounts=new long[numRefs];
		long[] refBaseCounts=new long[numRefs];
		long[] refAmbigCounts=new long[numRefs];
		
		for(int r=1, s=1; r<numRefs; r++){
			final int lim=s+(useRefNames ? 1 : refScafCounts[r]);
			while(s<lim){
				refReadCounts[r]+=scaffoldReadCounts.get(s);
				refBaseCounts[r]+=scaffoldBaseCounts.get(s);
				refAmbigCounts[r]+=scaffoldAmbigReadCounts.get(s);
				s++;
			}
		}
		
		/* Print header */
		tsw.print("#name\t%unambiguousReads\tunambiguousMB\t%ambiguousReads\tambiguousMB\tunambiguousReads\tambiguousReads\n");

		final float rmult=100f/Tools.max(1, totalReads);
		
		/* Print data */
		for(int i=1; i<refNames.size(); i++){
			final long reads=refReadCounts[i];
			final long bases=refBaseCounts[i];
			final float unambigMB=bases*0.000001f;
			
			final long ambigReads=refAmbigCounts[i];
			final long ambigBases=0; //TODO but not urgent
			final float ambigMB=ambigBases*0.000001f;
			
			final String name=ReadWrite.stripToCore(refNames.get(i));

			final double unambigReadP=rmult*reads;
			final double ambigReadP=rmult*ambigReads;
			if(reads>0 || !printNonZeroOnly){
				tsw.print(Tools.format("%s\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%d\n",name,unambigReadP,unambigMB,ambigReadP,ambigMB,reads,ambigReads));
			}
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Write taxonomic information.
	 */
	private void writeTaxonomy(){
		if(!USE_TAXTREE || outtax==null){return;}
		
		long mappedFrags=0;
		for(int i=0; i<scaffoldReadCounts.length(); i++){
			mappedFrags+=scaffoldFragCounts.get(i);
		}
		final double fragMult=100.0/Tools.max(1, fragsIn);
		
		final TextStreamWriter tsw=new TextStreamWriter(outtax, overwrite, false, false);
		tsw.start();
		
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		tsw.print(Tools.format("#Reads\t%d\n",fragsIn));
		tsw.print(Tools.format("#Mapped\t%d\n",mappedFrags));
		tsw.print(Tools.format("#Limits\t%d\t%d\t%d\t%d\n", taxNodeCountLimit, taxNodeNumberLimit, taxNodeMinLevel, taxNodeMaxLevel));
		tsw.print("#ID\tCount\tPercent\tLevel\tName\n");
		
//		assert(false) : taxNodeCountLimit+", "+taxNodeMinLevel+", "+taxNodeMaxLevel;
		ArrayList<TaxNode> nodes=tree.gatherNodesAtLeastLimit(taxNodeCountLimit, taxNodeMinLevel, taxNodeMaxLevel);
		
		for(int i=0, cap=Tools.min(nodes.size(), (taxNodeNumberLimit>0 ? taxNodeNumberLimit : Integer.MAX_VALUE)); i<cap; i++){
			TaxNode n=nodes.get(i);
			tsw.print(Tools.format("%d\t%d\t%.4f\t%s\t%s\n", n.id, n.countSum, n.countSum*fragMult, n.levelStringExtended(false), n.name));
		}
		
		tsw.poisonAndWait();
	}
	
//	/**
//	 * Fills the scaffold names array with reference names.
//	 */
//	private void toRefNames(){
//		final int numRefs=refNames.size();
//		for(int r=0, s=1; r<numRefs; r++){
//			final int scafs=refScafCounts[r];
//			final int lim=s+scafs;
//			final String name=ReadWrite.stripToCore(refNames.get(r));
////			outstream.println("r="+r+", s="+s+", scafs="+scafs+", lim="+lim+", name="+name);
//			while(s<lim){
////				outstream.println(r+", "+s+". Setting "+scaffoldNames.get(s)+" -> "+name);
//				scaffoldNames.set(s, name);
//				s++;
//			}
//		}
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static int numKmers(Read r1, Read r2, int k){
		int x=0;
		if(r1!=null){
			x+=Tools.max(r1.length()-k+1, 0);
		}
		if(r2!=null){
			x+=Tools.max(r2.length()-k+1, 0);
		}
		return x;
	}
	
	private void loadGiToTaxid(){
		Timer t=new Timer();
		outstream.println("Loading gi to taxa translation table.");
		GiToTaxid.initialize(giTableFile);
		t.stop();
		if(DISPLAY_PROGRESS){
			outstream.println("Time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
	}
	
//	private TaxTree loadTaxTree(){
//		assert(taxTreeFile!=null || (taxNameFile!=null && taxNodeFile!=null)) : "Must specify both taxname and taxnode files.";
//		Timer t=new Timer();
//		outstream.print("\nLoading tax tree; ");
//		final TaxTree tree;
//		if(taxTreeFile!=null){
//			tree=ReadWrite.read(TaxTree.class, taxTreeFile, true);
//		}else{
//			tree=new TaxTree(taxNameFile, taxNodeFile);
//		}
//		t.stop();
//		if(DISPLAY_PROGRESS){
//			outstream.println("time: \t"+t);
//			Shared.printMemory();
//			outstream.println();
//		}
//		return tree;
//	}
	
	private void addToTree(){
		for(int i=0; i<scaffoldFragCounts.length(); i++){
			long count=scaffoldFragCounts.get(i);
			if(count>0){
				String name=scaffoldNames.get(i);
				assert(name.startsWith("ncbi|") || name.startsWith("tid|") || (name.startsWith("gi|") && GiToTaxid.isInitialized())) :
					"\nFor taxonomy, all ref names must start with 'gi|' or 'ncbi|' or 'tid|'.\n" +
					"If the names start with 'gi', the gi= flag must be set.\n";
				int id=GiToTaxid.getID(name);
				if(id>-1){
					tree.incrementRaw(id, count);
				}
			}
		}
		tree.percolateUp();
	}
	
	/**
	 * Fills tables with kmers from references, using multiple LoadThread.
	 * @return Number of kmers stored.
	 */
	private long spawnLoadThreads(){
		Timer t=new Timer();
		if((ref==null || ref.size()<1) && (literal==null || literal.length<1)){return 0;}
		long added=0;
		
		final boolean oldParseCustom=FASTQ.PARSE_CUSTOM;
		FASTQ.PARSE_CUSTOM=false;
		
		/* Create load threads */
		LoadThread[] loaders=new LoadThread[WAYS];
		for(int i=0; i<loaders.length; i++){
			loaders[i]=new LoadThread(i);
			loaders[i].start();
		}
		
		boolean warned=false;
		
		/* For each reference file... */

		int refNum=1;
		if(ref!=null){
			
			HashMap<String, Integer> nameMap=new HashMap<String, Integer>();
			
			for(String refname : ref){

				/* Start an input stream */
				FileFormat ff=FileFormat.testInput(refname, FileFormat.FASTA, null, true, true);
				ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1L, false, ff, null, null, null, Shared.USE_MPI, true);
				cris.start(); //4567
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				final String core=ReadWrite.stripToCore(refname);
				if(useRefNames){
					assert(refNum==scaffoldNames.size());
					assert(!nameMap.containsKey(core)) : "Duplicate file name: "+core;
					Integer id=scaffoldNames.size();
					scaffoldNames.add(core);
					nameMap.put(core, id);
				}
				
				/* Iterate through read lists from the input stream */
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					{
						/* Assign a unique ID number to each scaffold */
						ArrayList<Read> reads2=new ArrayList<Read>(reads);
						for(Read r1 : reads2){
							final Read r2=r1.mate;
							if(useRefNames){
								r1.id=core;
								if(r2!=null){r2.id=core;}
							}else if(r1.id==null){r1.id=Integer.valueOf(scaffoldNames.size()).toString();}
							final Integer id;
							{
								String name=r1.id;
								Integer x=nameMap.get(r1.id);
								for(int suffix=2; !useRefNames && x!=null; suffix++){
									name=r1.id+"_"+suffix;
									if(!warned) {
										assert(RENAME_DUPLICATES) : ("Duplicate scaffold found: "+r1.id+"\n"+
												"Please ensure all scaffolds have unique names (using rename.sh), or add the flag 'renameduplicates'");
										System.err.println("WARNING: Duplicate scaffold is being renamed; future warnings will be suppressed."
												+ "\nTo force a crash when this occurs, add the flag 'renameduplicates=f'\nDuplicate name: "+r1.id);
									}
									warned=true;
									x=nameMap.get(name);
								}
								assert(x==null || useRefNames) : "Duplicate scaffold name "+r1.id+" found in scafnum "+x+" and "+scaffoldNames.size();
								if(x!=null){
									id=x;
								}else{
									id=scaffoldNames.size();
									scaffoldNames.add(name);
									Integer old=nameMap.put(name, id);
								}
							}
							if(useRefNames){assert(refNum==id);}
							
							refScafCounts[refNum]++;
							int len=r1.pairLength();
							r1.obj=id;
							if(r2!=null){r2.obj=id;}
							
							scaffoldLengths.increment(id, len);
						}
						
						if(REPLICATE_AMBIGUOUS){
							reads2=Tools.replicateAmbiguous(reads2, k);
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
		
//		int refNum=0;
//		if(ref!=null){
//			for(String refname : ref){
//
//				/* Start an input stream */
//				FileFormat ff=FileFormat.testInput(refname, FileFormat.FASTA, null, true, true);
//				ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1L, false, ff, null, null, null, Shared.USE_MPI, true);
//				cris.start(); //4567
//				ListNum<Read> ln=cris.nextList();
//				ArrayList<Read> reads=(ln!=null ? ln.list : null);
//
//				final String core=ReadWrite.stripToCore(refname);
//
//				/* Iterate through read lists from the input stream */
//				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//					{
//						/* Assign a unique ID number to each scaffold */
//						ArrayList<Read> reads2=new ArrayList<Read>(reads);
//						for(Read r1 : reads2){
//							final Read r2=r1.mate;
//							final Integer id=scaffoldNames.size();
//							refScafCounts[refNum]++;
//							scaffoldNames.add(r1.id==null ? id.toString() : r1.id);
//							int len=r1.length();
//							r1.obj=id;
//							if(r2!=null){
//								r2.obj=id;
//								len+=r2.length();
//							}
//							scaffoldLengths.add(len);
//						}
//
//						if(REPLICATE_AMBIGUOUS){
//							reads2=Tools.replicateAmbiguous(reads2, k);
//						}
//
//						/* Send a pointer to the read list to each LoadThread */
//						for(LoadThread lt : loaders){
//							boolean b=true;
//							while(b){
//								try {
//									lt.queue.put(reads2);
//									b=false;
//								} catch (InterruptedException e) {
//									//TODO:  This will hang due to still-running threads.
//									throw new RuntimeException(e);
//								}
//							}
//						}
//					}
//
//					/* Dispose of the old list and fetch a new one */
//					cris.returnList(ln);
//					ln=cris.nextList();
//					reads=(ln!=null ? ln.list : null);
//				}
//				/* Cleanup */
//				cris.returnList(ln);
//				errorState|=ReadWrite.closeStream(cris);
//				refNum++;
//			}
//		}

		/* If there are literal sequences to use as references */
		if(literal!=null){
			ArrayList<Read> list=new ArrayList<Read>(literal.length);
			if(verbose){outstream.println("Adding literals "+Arrays.toString(literal));}

			/* Assign a unique ID number to each literal sequence */
			if(useRefNames){
				final Integer id=scaffoldNames.size();
				scaffoldNames.add("literal");
				for(int i=0; i<literal.length; i++){
					final Read r=new Read(literal[i].getBytes(), null, id);
					refScafCounts[refNum]++;
					scaffoldLengths.increment(id, r.length());
					r.obj=id;
					list.add(r);
				}
			}else{
				for(int i=0; i<literal.length; i++){
					final int id=scaffoldNames.size();
					final Read r=new Read(literal[i].getBytes(), null, id);
					refScafCounts[refNum]++;
					scaffoldNames.add(""+id);
					scaffoldLengths.set(id, r.length());
					r.obj=id;
					list.add(r);
				}
			}
			
			if(REPLICATE_AMBIGUOUS){
				list=Tools.replicateAmbiguous(list, k);
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
			success&=lt.success;
		}
		if(!success){KillSwitch.kill("Failed loading ref kmers; aborting.");}
		
		//Correct statistics for number of threads, since each thread processes all reference data
		refKmers/=WAYS;
		refBases/=WAYS;
		refReads/=WAYS;
		
		scaffoldReadCounts=new AtomicLongArray(scaffoldNames.size());
		scaffoldFragCounts=new AtomicLongArray(scaffoldNames.size());
		scaffoldBaseCounts=new AtomicLongArray(scaffoldNames.size());
		scaffoldAmbigReadCounts=new AtomicLongArray(scaffoldNames.size());

		t.stop();
		if(DISPLAY_PROGRESS){
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
		
		FASTQ.PARSE_CUSTOM=oldParseCustom;
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
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, ff1.samOrBam(), ff1, ff2, qfin1, qfin2);
			cris.setSampleRate(samplerate, sampleseed);
			cris.start(); //4567
			paired=cris.paired();
			if(!ff1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		}
		
		/* Create read output streams */
		final ConcurrentReadOutputStream rosm, rosu;
		final MultiCros mcros;
		if(outu1!=null){
			final int buff=(!ordered ? 12 : Tools.max(32, 2*Shared.threads()));
			FileFormat ff1=FileFormat.testOutput(outu1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(outu2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			rosu=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, true);
			rosu.start();
		}else{rosu=null;}
		if(outm1!=null){
			final int buff=(!ordered ? 12 : Tools.max(32, 2*Shared.threads()));
			FileFormat ff1=FileFormat.testOutput(outm1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(outm2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			rosm=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, true);
			rosm.start();
		}else{rosm=null;}
		if(outpattern!=null){
			final int buff=(!ordered ? 12 : Tools.max(32, 2*Shared.threads()));
			mcros=new MultiCros(outpattern, null, ordered, overwrite, append, true, false, FileFormat.FASTQ, buff);
		}else{mcros=null;}
		
		if(rosu!=null || rosm!=null || mcros!=null){
			t.stop();
			outstream.println("Started output streams:\t"+t);
			t.start();
		}
		
		/* Optionally skip the first reads, since initial reads may have lower quality */
		if(skipreads>0){
			long skipped=0;

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(skipped<skipreads && reads!=null && reads.size()>0){
				skipped+=reads.size();
				
				if(rosm!=null){rosm.add(new ArrayList<Read>(1), ln.id);}
				if(rosu!=null){rosu.add(new ArrayList<Read>(1), ln.id);}
				
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
			if(reads==null || reads.isEmpty()){
				ReadWrite.closeStreams(cris, rosu, rosm);
				ReadWrite.closeStreams(mcros);
				outstream.println("Skipped all of the reads.");
				System.exit(0);
			}
		}
		
		/* Create ProcessThreads */
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new ProcessThread(cris, rosm, rosu, mcros, ALLOW_LOCAL_ARRAYS));}
		for(ProcessThread pt : alpt){pt.start();}
		
		/* Wait for threads to die, and gather statistics */
		for(ProcessThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			readsIn+=pt.readsInT;
			fragsIn+=pt.fragsInT;
			basesIn+=pt.basesInT;
			readsMatched+=pt.readsMatchedT;
			basesMatched+=pt.basesMatchedT;
			readsUnmatched+=pt.readsUnmatchedT;
			basesUnmatched+=pt.basesUnmatchedT;
			readsQTrimmed+=pt.readsQTrimmedT;
			basesQTrimmed+=pt.basesQTrimmedT;
			readsFTrimmed+=pt.readsFTrimmedT;
			basesFTrimmed+=pt.basesFTrimmedT;
			readsQFiltered+=pt.readsQFilteredT;
			basesQFiltered+=pt.basesQFilteredT;
			
			correctReads+=pt.correctT;
			incorrectReads+=pt.incorrectT;
			
			if(pt.scaffoldReadCountsT!=null && scaffoldReadCounts!=null){
				for(int i=0; i<pt.scaffoldReadCountsT.length; i++){scaffoldReadCounts.addAndGet(i, pt.scaffoldReadCountsT[i]);}
				pt.scaffoldReadCountsT=null;
			}
			if(pt.scaffoldBaseCountsT!=null && scaffoldBaseCounts!=null){
				for(int i=0; i<pt.scaffoldBaseCountsT.length; i++){scaffoldBaseCounts.addAndGet(i, pt.scaffoldBaseCountsT[i]);}
				pt.scaffoldBaseCountsT=null;
			}
			if(pt.scaffoldFragCountsT!=null && scaffoldFragCounts!=null){
				for(int i=0; i<pt.scaffoldFragCountsT.length; i++){scaffoldFragCounts.addAndGet(i, pt.scaffoldFragCountsT[i]);}
				pt.scaffoldFragCountsT=null;
			}
			if(pt.scaffoldAmbigReadCountsT!=null && scaffoldAmbigReadCounts!=null){
				for(int i=0; i<pt.scaffoldAmbigReadCountsT.length; i++){scaffoldAmbigReadCounts.addAndGet(i, pt.scaffoldAmbigReadCountsT[i]);}
				pt.scaffoldAmbigReadCountsT=null;
			}
			if(barcodeData!=null) {
				barcodeData.merge(pt.barcodeDataT);
			}
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris, rosu, rosm);
		errorState|=ReadWrite.closeStreams(mcros);
		errorState|=ReadStats.writeAll();
		
		t.stop();
		if(showSpeed){
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
					
					addedT+=addToMap(r1, refSkip);
					if(r2!=null){
						addedT+=addToMap(r2, refSkip);
					}
				}
				reads=fetch();
			}
			
//			if(AbstractKmerTable.TESTMODE){
//				for(int i=0; i<ll.size; i++){
//					assert(map.contains(ll.get(i), il.get(i)));
//					assert(!map.contains(ll.get(i), Integer.MAX_VALUE));
//				}
//				ll=null;
//				il=null;
//			}
			
			if(map.canRebalance() && map.size()>2L*map.arrayLength()){
				map.rebalance();
			}
			success=true;
		}

		/**
		 * @param r The current read to process
		 * @param skip Number of bases to skip between kmers
		 * @return Number of kmers stored
		 */
		private long addToMap(final Read r, final int skip){
			final byte[] bases=r.bases;
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			final long kmask=lengthMasks[k];
			long kmer=0;
			long rkmer=0;
			long added=0;
			int len=0;
			int totalKmers=0;
			
			if(tnum==0){
				if(storeRefBases){
					assert(r.mate==null);
					assert(scaffolds.size()==(Integer)r.obj) : scaffolds.size()+", "+(Integer)r.obj/*+"\n"+r.toFasta()*/;
					scaffolds.add(bases);
				}
				if(bases==null || bases.length<k){scaffoldKmers.add(0);}
			}
			
			if(bases!=null){
				refReadsT++;
				refBasesT+=bases.length;
			}
			if(bases==null || bases.length<k){return 0;}
			
			final int id=(Integer)r.obj;

			if(skip>1){ //Process while skipping some kmers
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=AminoAcid.baseToNumber[b];
					long x2=AminoAcid.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
					if(x<0){len=0; rkmer=0;}else{len++;}
					if(verbose){outstream.println("Scanning1 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=k){
						totalKmers++;
						if(len%skip==0){
							final long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
							added+=addToMap(kmer, rkmer, k, extraBase, id, kmask);
						}
					}
				}
			}else{ //Process all kmers
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=AminoAcid.baseToNumber[b];
					long x2=AminoAcid.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
					if(x<0){len=0; rkmer=0;}else{len++;}
					if(verbose){outstream.println("Scanning2 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=k){
						totalKmers++;
						final long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
						final long atm=addToMap(kmer, rkmer, k, extraBase, id, kmask);
						added+=atm;
					}
				}
			}
			refKmersT+=totalKmers;
			if(tnum==0){scaffoldKmers.add(totalKmers);}
			
			return added;
		}
		
		
		/**
		 * Adds this kmer to the table, including any mutations implied by editDistance or hammingDistance.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @param kmask0
		 * @return Number of kmers stored
		 */
		private long addToMap(final long kmer, final long rkmer, final int len, final long extraBase, final int id, final long kmask0){
			
			assert(kmask0==lengthMasks[len]) : kmask0+", "+len+", "+lengthMasks[len]+", "+Long.numberOfTrailingZeros(kmask0)+", "+Long.numberOfTrailingZeros(lengthMasks[len]);
			
			if(verbose){outstream.println("addToMap_A; len="+len+"; kMasks[len]="+lengthMasks[len]);}
			assert((kmer&kmask0)==0);
			final long added;
			if(hammingDistance==0){
				final long key=toValue(kmer, rkmer, kmask0);
				if(failsSpeed(key)){return 0;}
				if(key%WAYS!=tnum){return 0;}
				if(verbose){outstream.println("addToMap_B: "+AminoAcid.kmerToString(kmer&~lengthMasks[len], len)+" = "+key);}
//				int[] old=map.getValues(key, new int[1]);
				
//				int[] old=map.getValues(key, new int[1]); //123
				
				added=map.set(key, id);
//				assert(old==null || map.contains(key, old)); //123
//				assert(map.contains(key, id)); //123
//				ll.add(key);
//				il.add(id); assert(AbstractKmerTable.TESTMODE);
				
//				if(AbstractKmerTable.TESTMODE){
//					for(int i=0; i<ll.size; i++){
//						assert(map.contains(ll.get(i), il.get(i)));
//						assert(!map.contains(ll.get(i), Integer.MAX_VALUE));
//					}
//				}
				
			}else if(editDistance>0){
//				long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
				added=mutate(kmer, rkmer, len, id, editDistance, extraBase);
			}else{
				added=mutate(kmer, rkmer, len, id, hammingDistance, -1);
			}
			if(verbose){outstream.println("addToMap added "+added+" keys.");}
			return added;
		}

//		private LongList ll=new LongList();
//		private IntList il=new IntList();
		
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
			
			if(verbose){outstream.println("mutate_A; len="+len+"; kmer="+kmer+"; rkmer="+rkmer+"; kMasks[len]="+lengthMasks[len]);}
			if(key%WAYS==tnum){
				if(verbose){outstream.println("mutate_B: "+AminoAcid.kmerToString(kmer&~lengthMasks[len], len)+" = "+key);}
				int x=map.set(key, id);
				if(verbose){outstream.println("mutate_B added "+x+" keys.");}
				added+=x;
				assert(map.contains(key));
			}
			
			if(dist>0){
				final int dist2=dist-1;
				
				//Sub
				for(int j=0; j<4; j++){
					for(int i=0; i<len; i++){
						final long temp=(kmer&clearMasks[i])|setMasks[j][i];
						if(temp!=kmer){
							long rtemp=rcomp(temp, len);
							added+=mutate(temp, rtemp, len, id, dist2, extraBase);
						}
					}
				}
				
				if(editDistance>0){
					//Del
					if(extraBase>=0 && extraBase<=3){
						for(int i=1; i<len; i++){
							final long temp=(kmer&leftMasks[i])|((kmer<<2)&rightMasks[i])|extraBase;
							if(temp!=kmer){
								long rtemp=rcomp(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, -1);
							}
						}
					}

					//Ins
					final long eb2=kmer&3;
					for(int i=1; i<len; i++){
						final long temp0=(kmer&leftMasks[i])|((kmer&rightMasks[i])>>2);
						for(int j=0; j<4; j++){
							final long temp=temp0|setMasks[j][i-1];
							if(temp!=kmer){
								long rtemp=rcomp(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, eb2);
							}
						}
					}
				}

			}
			
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
		 * @param rosu_ Unmatched read output stream (optional)
		 * @param rosm_ Matched read output stream (optional)
		 */
		public ProcessThread(ConcurrentReadInputStream cris_, ConcurrentReadOutputStream rosm_, ConcurrentReadOutputStream rosu_,
				MultiCros mcros_, boolean localArrays){
			cris=cris_;
			rosm=rosm_;
			rosu=rosu_;
			mcros=mcros_;
			
			readstats=(ReadStats.collectingStats()) ? new ReadStats() : null;
			barcodeDataT=(trackBarcodes ? new BarcodeMappingStats() : null);
			
			final int alen=(scaffoldNames==null ? 0 : scaffoldNames.size());
			if(localArrays && alen>0 && alen<10000){
				scaffoldReadCountsT=new long[alen];
				scaffoldBaseCountsT=new long[alen];
				scaffoldFragCountsT=new long[alen];
				scaffoldAmbigReadCountsT=new long[alen];
			}else{
				scaffoldReadCountsT=scaffoldBaseCountsT=scaffoldFragCountsT=scaffoldAmbigReadCountsT=null;
			}
			
			if(USE_COUNTVECTOR){
				countVector=new IntList(1000);
				countArray=null;
			}else{
				countVector=null;
				countArray=new int[alen];
			}
		}
		
		@Override
		public void run(){
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			final ArrayList<Read> mlist=(rosm==null ? null : new ArrayList<Read>(Shared.bufferLen()));
			final ArrayList<Read> ulist=(rosu==null ? null : new ArrayList<Read>(Shared.bufferLen()));
			final ArrayListSet als=(outpattern==null ? null : new ArrayListSet(ordered));
			
			//While there are more reads lists...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					final Read r1=reads.get(i);
					final Read r2=r1.mate;
					
					if(!r1.validated()){r1.validate(true);}
					if(r2!=null && !r2.validated()){r2.validate(true);}
					
					if(readstats!=null){
						readstats.addToHistograms(r1);
					}

					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());

					final int minlen1=(int)Tools.max(initialLength1*minLenFraction, minReadLength);
					final int minlen2=(int)Tools.max(initialLength2*minLenFraction, minReadLength);
					
					if(loglog!=null){loglog.hash(r1);}
					
					if(verbose){outstream.println("Considering read "+r1.id+" "+new String(r1.bases));}
					
					fragsInT++;
					readsInT+=(1+r1.mateCount());
					basesInT+=r1.pairLength();
					
					boolean remove=false;
					
					if(chastityFilter){
						if(r1!=null && r1.failsChastity()){
							basesQFilteredT+=r1.length();
							readsQFilteredT++;
							r1.setDiscarded(true);
						}
						if(r2!=null && r2.failsChastity()){
							basesQFilteredT+=r2.length();
							readsQFilteredT++;
							r2.setDiscarded(true);
						}
					}
					
					if(forceTrimLeft>0 || forceTrimRight>0 || forceTrimRight2>0 || forceTrimModulo>0){
						if(r1!=null && !r1.discarded()){
							final int len=r1.length();
							final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
							final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
							final int b1=forceTrimRight>0 ? forceTrimRight : len;
							final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
							final int b=Tools.min(b0, b1, b2);
							final int x=TrimRead.trimToPosition(r1, a, b, 1);
							basesFTrimmedT+=x;
							readsFTrimmedT+=(x>0 ? 1 : 0);
							if(r1.length()<minlen1){r1.setDiscarded(true);}
						}
						if(r2!=null && !r2.discarded()){
							final int len=r2.length();
							final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
							final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
							final int b1=forceTrimRight>0 ? forceTrimRight : len;
							final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
							final int b=Tools.min(b0, b1, b2);
							final int x=TrimRead.trimToPosition(r2, a, b, 1);
							basesFTrimmedT+=x;
							readsFTrimmedT+=(x>0 ? 1 : 0);
							if(r2.length()<minlen2){r2.setDiscarded(true);}
						}
					}
					
					if(removePairsIfEitherBad){remove=r1.discarded() || (r2!=null && r2.discarded());}
					else{remove=r1.discarded() && (r2==null || r2.discarded());}
					

					
					if(!remove){
						//Do quality trimming
						
						int rlen1=0, rlen2=0;
						if(r1!=null){
							if(qtrimLeft || qtrimRight){
								int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, trimE, 1);
								basesQTrimmedT+=x;
								readsQTrimmedT+=(x>0 ? 1 : 0);
							}
							rlen1=r1.length();
							if(rlen1<minlen1 || rlen1>maxReadLength){
								r1.setDiscarded(true);
								if(verbose){outstream.println(r1.id+" discarded due to length.");}
							}
						}
						if(r2!=null){
							if(qtrimLeft || qtrimRight){
								int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, trimE, 1);
								basesQTrimmedT+=x;
								readsQTrimmedT+=(x>0 ? 1 : 0);
							}
							rlen2=r2.length();
							if(rlen2<minlen2 || rlen2>maxReadLength){
								r2.setDiscarded(true);
								if(verbose){outstream.println(r2.id+" discarded due to length.");}
							}
						}
						
						//Discard reads if too short
						if((removePairsIfEitherBad && (r1.discarded() || (r2!=null && r2.discarded()))) ||
								(r1.discarded() && (r2==null || r2.discarded()))){
							basesQFilteredT+=r1.pairLength();
							readsQTrimmedT+=r1.pairCount();
							remove=true;
						}
					}

					
					if(!remove){
						//Do quality filtering
						
						//Determine whether to discard the reads based on average quality
						if(minAvgQuality>0){
							if(r1!=null && r1.quality!=null && r1.avgQuality(false, minAvgQualityBases)<minAvgQuality){
								r1.setDiscarded(true);
								if(verbose){outstream.println(r1.id+" discarded due to low quality.");}
							}
							if(r2!=null && r2.quality!=null && r2.avgQuality(false, minAvgQualityBases)<minAvgQuality){
								r2.setDiscarded(true);
								if(verbose){outstream.println(r2.id+" discarded due to low quality.");}
							}
						}
						//Determine whether to discard the reads based on the presence of Ns
						if(maxNs>=0){
							if(r1!=null && r1.countUndefined()>maxNs){
								r1.setDiscarded(true);
								if(verbose){outstream.println(r1.id+" discarded due to Ns.");}
							}
							if(r2!=null && r2.countUndefined()>maxNs){
								r2.setDiscarded(true);
								if(verbose){outstream.println(r2.id+" discarded due to Ns.");}
							}
						}
						
						//Determine whether to discard the reads based on a lack of useful kmers
						if(minConsecutiveBases>0){
							if(r1!=null && !r1.discarded() && !r1.hasMinConsecutiveBases(minConsecutiveBases)){r1.setDiscarded(true);}
							if(r2!=null && !r2.discarded() && !r2.hasMinConsecutiveBases(minConsecutiveBases)){r2.setDiscarded(true);}
						}
						
						//Discard reads if too short
						if((removePairsIfEitherBad && (r1.discarded() || (r2!=null && r2.discarded()))) ||
								(r1.discarded() && (r2==null || r2.discarded()))){
							basesQFilteredT+=r1.pairLength();
							readsQFilteredT+=r1.pairCount();
							remove=true;
						}
					}
					
					final int sites, assigned;
					if(remove){
						if(r1!=null){
							basesQFilteredT+=r1.length();
							readsQFilteredT++;
						}
						if(r2!=null){
							basesQFilteredT+=r2.length();
							readsQFilteredT++;
						}
						sites=assigned=0;
					}else{
						
						if(ecc && r1!=null && r2!=null){BBMerge.findOverlapStrict(r1, r2, true);}
						
						//Do kmer matching
						if(keepPairsTogether){
							
							final int a, b, max;
							
							if(countArray==null){
								countVector.size=0;
								a=findBestMatch(r1, keySets, countVector);
								b=findBestMatch(r2, keySets, countVector);
								if(verbose){outstream.println("countVector: "+countVector);}
								max=condenseLoose(countVector, idList1, countList1);
							}else{
								idList1.size=0;
								a=findBestMatch(r1, keySets, countArray, idList1);
								b=findBestMatch(r2, keySets, countArray, idList1);

								max=condenseLoose(countArray, idList1, countList1);
							}
							
							if(verbose){
								outstream.println("idList1: "+idList1);
								outstream.println("countList1: "+countList1);
							}
							final int cz=clearzoneFraction<=0 ? clearzone : 
								Tools.max(clearzone, (int)Math.ceil(clearzoneFraction*(r1.numValidPairKmers(k))));
							filterTopScaffolds(r1, r2, idList1, countList1, finalList1, max, cz);
							if(verbose){
								outstream.println("idList1: "+idList1);
								outstream.println("countList1: "+countList1);
								outstream.println("finalList1: "+finalList1);
							}
							sites=finalList1.size;
							
							final int minhits=Tools.max(minKmerHits, (int)(minKmerFraction*numKmers(r1, r2, k)));
							if(max>=minhits){
								assigned=assignTogether(r1, r2, als);
							}else{
								readsUnmatchedT+=r1.pairCount();
								basesUnmatchedT+=r1.pairLength();
								assigned=0;
							}
							
							if(rename){
								rename(r1, idList1, countList1, max, assigned);
								rename(r2, idList1, countList1, max, assigned);
							}
							
						}else{
							final int max1, max2, a, b;
							{
								if(countArray==null){
									countVector.size=0;
									a=findBestMatch(r1, keySets, countVector);
									max1=condenseLoose(countVector, idList1, countList1);
								}else{
									idList1.size=0;
									a=findBestMatch(r1, keySets, countArray, idList1);
									max1=condenseLoose(countArray, idList1, countList1);
								}
//								if(rename){rename(r1, idList1, countList1, max1);}
								final int cz=clearzoneFraction<=0 ? clearzone : 
									Tools.max(clearzone, (int)Math.ceil(clearzoneFraction*(r1.numValidKmers(k))));
								filterTopScaffolds(r1, null, idList1, countList1, finalList1, max1, cz);
							}
							if(r2!=null){
								if(countArray==null){
									countVector.size=0;
									b=findBestMatch(r2, keySets, countVector);
									max2=condenseLoose(countVector, idList2, countList1);
								}else{
									idList2.size=0;
									b=findBestMatch(r2, keySets, countArray, idList2);
									max2=condenseLoose(countArray, idList2, countList2);
								}
//								if(rename){rename(r2, idList2, countList2, max2);}
								final int cz=clearzoneFraction<=0 ? clearzone : 
									Tools.max(clearzone, (int)Math.ceil(clearzoneFraction*(r2.numValidKmers(k))));
								filterTopScaffolds(r2, null, idList2, countList2, finalList2, max2, cz);
							}else{max2=0;}
							
							sites=finalList1.size+finalList2.size;
							
							assigned=assignIndependently(r1, r2, max1, max2, als);
							if(rename) {
								rename(r1, idList1, countList1, max1, assigned);
								if(r2!=null) {rename(r2, idList2, countList2, max2, assigned);}
							}
						}
					}
					
					if(remove || assigned<1){
						if(ulist!=null){
							if(loglogOut!=null){loglogOut.hash(r1);}
							ulist.add(r1);
						}
					}else{
						if(mlist!=null){mlist.add(r1);}
					}
					
				}
				
				//Send matched list to matched output stream
				if(rosu!=null){
					rosu.add(ulist, ln.id);
					ulist.clear();
				}
				
				//Send unmatched list to unmatched output stream
				if(rosm!=null){
					rosm.add(mlist, ln.id);
					mlist.clear();
				}
				
				if(mcros!=null){
					mcros.add(als, ln.id);
				}
				
				//Fetch a new read list
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
		}
		
		/*--------------------------------------------------------------*/
		/*----------------        Helper Methods        ----------------*/
		/*--------------------------------------------------------------*/
		
		/**
		 * @param r
		 * @param idList
		 * @param countList
		 */
		private void rename(Read r, IntList idList, IntList countList, final int max, int assignedSites) {
			if(r==null || (idList.size<1 && !addUnknown)){return;}
			StringBuilder sb=new StringBuilder();
			if(r.id==null){sb.append(r.numericID);}
			else{sb.append(r.id);}
			if(idList.size>0 && (assignedSites>0 || !topHitOnly)) {
				for(int i=0; i<idList.size; i++){
					int id=idList.get(i);
					int count=countList.get(i);
					if(count==max || !topHitOnly) {
						sb.append('\t');
						sb.append(scaffoldNames.get(id));
						if(addCount) {
							sb.append('=');
							sb.append(count);
						}
						if(topHitOnly) {break;}
					}
				}
			}else if(addUnknown) {
				sb.append('\t').append("unknown");
			}
			r.id=sb.toString();
		}
		
		private void addGeneSet(IntList finalList1, int start, int stop, int readSum, int lenSum) {
			if(MAKE_GENE_SETS && stop>start+1) {
				finalList1.sort();
				ByteBuilder bb=new ByteBuilder();
				int len=0;
				for(int j=start; j<stop; j++) {
					int id=finalList1.get(j);
					if(!bb.isEmpty()) {bb.comma();}
					bb.append(scaffoldNames.get(id));
					len+=scaffoldLengths.get(id);
				}
				String name=bb.toString();
				StringCount sc=setMap.get(name);
				if(sc==null) {
					synchronized(setMap) {
						sc=setMap.get(name);
						if(sc==null) {
							sc=new StringCount(name);
							sc.length=len;
							setMap.put(name, sc);
						}
					}
				}
				assert(sc!=null);
				synchronized(sc) {
					sc.reads+=readSum;
					sc.bases+=lenSum;
					sc.ambigReads+=readSum;
				}
			}
		}
		
		/**
		 * @param r1 Read 1
		 * @param r2 Read 2
		 * @return Number of sites assigned
		 */
		private int assignTogether(Read r1, Read r2, ArrayListSet als){
			final int sites=finalList1.size;
			final int lenSum=r1.length()+(r1.mateLength());
			final int readSum=1+(r2==null ? 0 : 1);
			final int start, stop;
			final String barcode=(trackBarcodes ? r1.barcode(true) : null);

			if(sites<2 || ambigMode==AMBIG_ALL){
				start=0;
				stop=sites;
			}else if(ambigMode==AMBIG_TOSS){
				start=stop=0;
			}else if(ambigMode==AMBIG_FIRST){
				finalList1.sort();
				start=0;
				stop=1;
			}else if(ambigMode==AMBIG_RANDOM){
				start=(int)(r1.numericID%sites);
				stop=start+1;
			}else{
				throw new RuntimeException("Unknown mode "+ambigMode);
			}
			
			if(MAKE_GENE_SETS && stop>start+1) {addGeneSet(finalList1, start, stop, readSum, lenSum);}
			
			for(int j=start; j<stop; j++){
				int id=finalList1.get(j);
				
				if(als!=null){
					als.add(r1, scaffoldNames.get(id));
				}
				
				if(parsecustom && j==start){
					String scafName=scaffoldNames.get(id);
					String rname=r1.parseCustomRname();
					if(scafName.equals(rname)){
						correctT+=(1+r1.mateCount());
					}else{
						incorrectT+=(1+r1.mateCount());
					}
				}
				
				if(trackBarcodes){barcodeDataT.incrementSourceMap(barcode, scaffoldNames.get(id), r1.pairCount());}
				if(scaffoldReadCountsT!=null){
					scaffoldReadCountsT[id]+=readSum;
					scaffoldBaseCountsT[id]+=lenSum;
					scaffoldFragCountsT[id]++;
					if(sites>1) {scaffoldAmbigReadCountsT[id]+=readSum;}
				}else{
					scaffoldReadCounts.addAndGet(id, readSum);
					scaffoldBaseCounts.addAndGet(id, lenSum);
					scaffoldFragCounts.addAndGet(id, 1);
					if(sites>1) {scaffoldAmbigReadCounts.addAndGet(id, readSum);}
				}
			}
			
			if(trackBarcodes){barcodeDataT.incrementCodeMap(barcode, r1.pairCount());}
			if(start<stop){
				readsMatchedT+=r1.pairCount();
				basesMatchedT+=r1.pairLength();
			}else{
				readsUnmatchedT+=r1.pairCount();
				basesUnmatchedT+=r1.pairLength();
				if(trackBarcodes){barcodeDataT.incrementSourceMap(barcode, sites>1 ? "AMBIG" : "UNKNOWN", r1.pairCount());}
			}
			
			return stop-start;
		}
		
		/**
		 * @param r1 Read 1
		 * @param r2 Read 2
		 * @param max1 Highest match count for read 1
		 * @param max2 Highest match count for read 2
		 * @return Number of sites assigned
		 */
		private int assignIndependently(Read r1, Read r2, int max1, int max2, ArrayListSet als){
			assert(als==null || r2==null) : "Pattern output does not work with keepPairsTogether=false and paired reads\n"+als;
			int assigned=0;
			final String barcode=(trackBarcodes ? r1.barcode(true) : null);
			if(trackBarcodes){barcodeDataT.incrementCodeMap(barcode, r1.pairCount());}
			
			if(max1>=Tools.max(minKmerHits, (int)(minKmerFraction*numKmers(r1, null, k)))){
				final int sites=finalList1.size;
				final int lenSum=r1.length();
				final int start, stop;
				
				if(sites<2 || ambigMode==AMBIG_ALL){
					start=0;
					stop=sites;
				}else if(ambigMode==AMBIG_TOSS){
					start=stop=0;
				}else if(ambigMode==AMBIG_FIRST){
					finalList1.sort();
					start=0;
					stop=1;
				}else if(ambigMode==AMBIG_RANDOM){
					start=(int)(r1.numericID%sites);
					stop=start+1;
				}else{
					throw new RuntimeException("Unknown mode "+ambigMode);
				}
				
				if(MAKE_GENE_SETS && stop>start+1) {addGeneSet(finalList1, start, stop, 1, lenSum);}
				
				for(int j=start; j<stop; j++){
					int id=finalList1.get(j);
					
					if(als!=null){
						als.add(r1, scaffoldNames.get(id));
					}
					
					if(parsecustom && j==start){
						String scafName=scaffoldNames.get(id);
						String rname=r1.parseCustomRname();
						if(scafName.equals(rname)){
							correctT++;
						}else{
							incorrectT++;
						}
					}

					if(trackBarcodes){barcodeDataT.incrementSourceMap(barcode, scaffoldNames.get(id), 1);}
					if(scaffoldReadCountsT!=null){
						scaffoldReadCountsT[id]++;
						scaffoldBaseCountsT[id]+=lenSum;
						if(max1>=max2){
							scaffoldFragCountsT[id]++;
						}
						if(sites>1) {scaffoldAmbigReadCountsT[id]++;}
					}else{
						scaffoldReadCounts.addAndGet(id, 1);
						scaffoldBaseCounts.addAndGet(id, lenSum);
						if(max1>=max2){
							scaffoldFragCounts.addAndGet(id, 1);
						}
						if(sites>1) {scaffoldAmbigReadCounts.incrementAndGet(id);}
					}
				}
				if(start<stop){
					readsMatchedT++;
					basesMatchedT+=r1.length();
					assigned+=(stop-start);
				}else{
					readsUnmatchedT++;
					basesUnmatchedT+=r1.length();
					if(trackBarcodes){barcodeDataT.incrementSourceMap(barcode, sites>1 ? "AMBIG" : "UNKNOWN", 1);}
				}
			}
			
			if(max2>=Tools.max(minKmerHits, (int)(minKmerFraction*numKmers(r2, null, k)))){
				final int sites=finalList2.size;
				final int lenSum=r2.length();
				final int start, stop;

				if(sites<2 || ambigMode==AMBIG_ALL){
					start=0;
					stop=sites;
				}else if(ambigMode==AMBIG_TOSS){
					start=stop=0;
				}else if(ambigMode==AMBIG_FIRST){
					finalList2.sort();
					start=0;
					stop=1;
				}else if(ambigMode==AMBIG_RANDOM){
					start=(int)(r2.numericID%sites);
					stop=start+1;
				}else{
					throw new RuntimeException("Unknown mode "+ambigMode);
				}
				
				if(MAKE_GENE_SETS && stop>start+1) {addGeneSet(finalList1, start, stop, 1, lenSum);}
				
				for(int j=start; j<stop; j++){
					int id=finalList2.get(j);
					
					if(als!=null){
						als.add(r2, scaffoldNames.get(id));
						throw new RuntimeException("Pattern output does not currently work with keepPairsTogether=false");
					}
					
					if(parsecustom && j==start){
						String scafName=scaffoldNames.get(id);
						String rname=r2.parseCustomRname();
						if(scafName.equals(rname)){
							correctT++;
						}else{
							incorrectT++;
						}
					}

					if(trackBarcodes){barcodeDataT.incrementSourceMap(barcode, scaffoldNames.get(id), 1);}
					if(scaffoldReadCountsT!=null){
						scaffoldReadCountsT[id]++;
						scaffoldBaseCountsT[id]+=lenSum;
						if(max2>max1){
							scaffoldFragCountsT[id]++;
						}
						if(sites>1) {scaffoldAmbigReadCountsT[id]++;}
					}else{
						scaffoldReadCounts.addAndGet(id, 1);
						scaffoldBaseCounts.addAndGet(id, lenSum);
						if(max2>max1){
							scaffoldFragCounts.addAndGet(id, 1);
						}
						if(sites>1) {scaffoldAmbigReadCounts.incrementAndGet(id);}
					}
				}
				if(start<stop){
					readsMatchedT++;
					basesMatchedT+=r2.length();
					assigned+=(stop-start);
				}else{
					readsUnmatchedT++;
					basesUnmatchedT+=r2.length();
					if(trackBarcodes){barcodeDataT.incrementSourceMap(barcode, sites>1 ? "AMBIG" : "UNKNOWN", 1);}
				}
			}
			
			return assigned;
		}
		
		/**
		 * Pack a list of nonunique values into a list of unique values and a list of their counts.
		 * @param loose Nonunique values
		 * @param packed Unique values
		 * @param counts Counts of values
		 * @return
		 */
		private int condenseLoose(IntList loose, IntList packed, IntList counts){
			packed.size=0;
			counts.size=0;
			if(loose.size<1){return 0;}
			loose.sort();
			int prev=-1;
			int max=0;
			int count=0;
			for(int i=0; i<loose.size; i++){
				int id=loose.get(i);
//				outstream.println("i="+i+", id="+id+", count="+count+", prev="+prev);
				if(id==prev){
					count++;
				}else{
					if(count>0){
						packed.add(prev);
						counts.add(count);
						max=Tools.max(count, max);
//						assert(false) : "i="+i+", "+id+", "+count+", "+packed+", "+counts;
					}
					prev=id;
					count=1;
				}
			}
			if(count>0){
				packed.add(prev);
				counts.add(count);
				max=Tools.max(count, max);
			}
			return max;
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
		
		private void filterTopScaffolds(Read r1, Read r2, IntList packed, IntList counts, IntList out, int max, int cz){
			out.size=0;
			if(packed.size<1){return;}
			if(processContainedRef){
				filterTopScaffolds_withContainedRef(r1, r2, packed, counts, out);
			}else{
				filterTopScaffolds_withClearzone(packed, counts, out, max, cz);
			}
		}
		
		private void filterTopScaffolds_withContainedRef(Read r1, Read r2, IntList packed, IntList counts, IntList out){
			for(int i=0; i<packed.size; i++){
				final int p=packed.get(i);
				final int c=counts.get(i);
				if(storeRefBases){
					if(Tools.containsForward(r1.bases, scaffolds.get(p), hammingDistance)>=0 ||
							(r2!=null && Tools.containsForward(r2.bases, scaffolds.get(p), hammingDistance)>=0)){
						out.add(p);
					}else if(rcomp && (Tools.containsReverse(r1.bases, scaffolds.get(p), hammingDistance)>=0 ||
							(r2!=null && Tools.containsReverse(r2.bases, scaffolds.get(p), hammingDistance)>=0))){
						out.add(p);
					}
				}else if(c>=scaffoldKmers.get(p)){
					out.add(p);
				}
			}
		}
		
		private void filterTopScaffolds_withClearzone(IntList packed, IntList counts, IntList out, int max, int cz){
			final int thresh=Tools.max(1, max-cz);
			for(int i=0; i<packed.size; i++){
				final int c=counts.get(i);
				assert((c>0) && c<=max) : c+"\n"+packed+"\n"+counts;
				if(c>=thresh){
					out.add(packed.get(i));
				}
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
		private final int[] getValues(final long kmer, final long rkmer, final long lengthMask, final int qPos, final int len, final int qHDist, final AbstractKmerTable[] sets){
			if(qHDist>0){return getValuesQHD(kmer, rkmer, lengthMask, qPos, len, qHDist, sets, qhList);}
			int[] ids=getValuesInner(kmer, rkmer, lengthMask, qPos, sets);
			if((ids==null || ids[0]<0) && qHDist>0){
				final int qHDist2=qHDist-1;
				
				//Sub
				for(int j=0; j<4; j++){
					for(int i=0; i<len; i++){
						final long temp=(kmer&clearMasks[i])|setMasks[j][i];
						if(temp!=kmer){
							long rtemp=rcomp(temp, len);
							ids=getValues(temp, rtemp, lengthMask, qPos, len, qHDist2, sets);
							if(ids!=null && ids[0]>-1){return ids;}
						}
					}
				}
			}
			return ids;
		}
		
		private final int[] getValuesQHD(final long kmer, final long rkmer, final long lengthMask, final int qPos, final int len, final int qHDist, final AbstractKmerTable[] sets, IntList list){
			assert(qHDist>0);
			list.clear();
			getValuesQHD_inner(kmer, rkmer, lengthMask, qPos, len, qHDist, sets, list);
			if(list.size>qhdistSizeLimit){
				list.sort();
				list.shrinkToUnique();
			}
			if(list.size>0){
				list.add(-1);//indicates end
				return list.array;
			}
			return null;
		}
			
		private final void getValuesQHD_inner(final long kmer, final long rkmer, final long lengthMask, final int qPos, final int len, final int qHDist, final AbstractKmerTable[] sets, IntList list){
			final int sizeLimit=10;
			int[] ids=getValuesInner(kmer, rkmer, lengthMask, qPos, sets);
			if(ids!=null){
				for(int x : ids){
					if(x<0){break;}
					if(list.size>sizeLimit || !list.contains(x)){list.add(x);}
				}
			}
			if(qHDist>0){
				final int qHDist2=qHDist-1;
				
				//Sub
				for(int j=0; j<4; j++){
					for(int i=0; i<len; i++){
						final long temp=(kmer&clearMasks[i])|setMasks[j][i];
						if(temp!=kmer){
							long rtemp=rcomp(temp, len);
							getValuesQHD_inner(temp, rtemp, lengthMask, qPos, len, qHDist2, sets, list);
						}
					}
				}
			}
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
		private final int[] getValuesInner(final long kmer, final long rkmer, final long lengthMask, final int qPos, final AbstractKmerTable[] sets){
			assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
			if(qSkip>1 && (qPos%qSkip!=0)){return null;}
			
			final long max=(rcomp ? Tools.max(kmer, rkmer) : kmer);
			final long key=(max&middleMask)|lengthMask;
			if(passesSpeed(key)){
				if(verbose){outstream.println("Testing key "+key);}
				AbstractKmerTable set=sets[(int)(key%WAYS)];
				if(verbose){outstream.println("Found set "+set.arrayLength()+", "+set.size());}
				final int[] ids=set.getValues(key, singleton);
				if(verbose){outstream.println("Found array "+(ids==null ? "null" : Arrays.toString(ids)));}
				return ids;
			}
			return null;
		}
		
		
		/**
		 * Returns number of matching kmers.
		 * @param r Read to process
		 * @param sets Kmer tables
		 * @return number of total kmer matches
		 */
		private final int findBestMatch(final Read r, final AbstractKmerTable sets[], IntList hits){
			if(r==null || r.bases==null || storedKmers<1){return 0;}
			final byte[] bases=r.bases;
			final int minlen=k-1;
			final int minlen2=(maskMiddle ? (k-midMaskLen)/2 : k);
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			final long kmask=lengthMasks[k];
			long kmer=0;
			long rkmer=0;
			int found=0;
			int len=0;
			
			if(bases==null || bases.length<k){return -1;}
			final int start=(restrictRight<1 ? 0 : Tools.max(0, bases.length-restrictRight));
			final int stop=(restrictLeft<1 ? bases.length : Tools.min(bases.length, restrictLeft));

			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			for(int i=start; i<stop; i++){
				byte b=bases[i];
				long x=AminoAcid.baseToNumber0[b];
				long x2=AminoAcid.baseToComplementNumber0[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(b=='N' && forbidNs){len=0; rkmer=0;}else{len++;}
				if(verbose){outstream.println("Scanning6 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=minlen2 && i>=minlen){
					final int[] ids=getValues(kmer, rkmer, kmask, i, k, qHammingDistance, sets);
					if(ids!=null && ids[0]>-1){
						for(int id : ids){
							if(id==-1){break;}
							hits.add(id);
						}
						if(verbose){outstream.println("Found = "+(found+1)+"/"+minKmerHits);}
						found++;
						//							assert(false) : (matchMode==MATCH_FIRST)+", "+(matchMode==MATCH_UNIQUE)+", "+ (ids.length==1 || ids[1]==-1);
						if(matchMode==MATCH_FIRST || (matchMode==MATCH_UNIQUE && (ids.length==1 || ids[1]==-1))){break;}
					}
				}
			}
			return found;
		}
		
		/**
		 * Returns number of matching kmers.
		 * @param r Read to process
		 * @param sets Kmer tables
		 * @return number of total kmer matches
		 */
		private final int findBestMatch(final Read r, final AbstractKmerTable sets[], int[] hits, IntList idList){
			if(r==null || r.bases==null || storedKmers<1){return 0;}
			final byte[] bases=r.bases;
			final int minlen=k-1;
			final int minlen2=(maskMiddle ? (k-midMaskLen)/2 : k);
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			final long kmask=lengthMasks[k];
			long kmer=0;
			long rkmer=0;
			int found=0;
			int len=0;
					
			if(bases==null || bases.length<k){return -1;}
			final int start=(restrictRight<1 ? 0 : Tools.max(0, bases.length-restrictRight));
			final int stop=(restrictLeft<1 ? bases.length : Tools.min(bases.length, restrictLeft));

			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			for(int i=start; i<stop; i++){
				byte b=bases[i];
				long x=AminoAcid.baseToNumber0[b];
				long x2=AminoAcid.baseToComplementNumber0[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(b=='N' && forbidNs){len=0; rkmer=0;}else{len++;}
				if(verbose){outstream.println("Scanning6 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=minlen2 && i>=minlen){
					final int[] ids=getValues(kmer, rkmer, kmask, i, k, qHammingDistance, sets);
					if(ids!=null && ids[0]>-1){
						for(int id : ids){
							if(id==-1){break;}
							hits[id]++;
							if(hits[id]==1){idList.add(id);}
						}
						if(verbose){outstream.println("Found = "+(found+1)+"/"+minKmerHits);}
						found++;
//						assert(false) : (matchMode==MATCH_FIRST)+", "+(matchMode==MATCH_UNIQUE)+", "+ (ids.length==1 || ids[1]==-1);
						if(matchMode==MATCH_FIRST || (matchMode==MATCH_UNIQUE && (ids.length==1 || ids[1]==-1))){break;}
					}
				}
			}
			return found;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Input read stream */
		private final ConcurrentReadInputStream cris;
		/** Output read streams */
		private final ConcurrentReadOutputStream rosm, rosu;
		/** Output pattern read streams */
		private final MultiCros mcros;
		
		private final ReadStats readstats;
		private final IntList countVector;
		private final int[] countArray;
		private final BarcodeMappingStats barcodeDataT;
		
		private final IntList idList1=new IntList(), idList2=new IntList();
		private final IntList countList1=new IntList(), countList2=new IntList();
		private final IntList finalList1=new IntList(), finalList2=new IntList();
		private final IntList qhList=new IntList();
		
		long[] scaffoldReadCountsT;
		long[] scaffoldBaseCountsT;
		long[] scaffoldFragCountsT;
		long[] scaffoldAmbigReadCountsT;
		final int[] singleton=new int[1];
		
		private long readsInT=0;
		private long fragsInT=0;
		private long basesInT=0;
		private long readsMatchedT=0;
		private long basesMatchedT=0;
		private long readsUnmatchedT=0;
		private long basesUnmatchedT=0;
		
		private long readsQTrimmedT=0;
		private long basesQTrimmedT=0;
		private long readsFTrimmedT=0;
		private long basesFTrimmedT=0;
		private long readsQFilteredT=0;
		private long basesQFilteredT=0;

		private long correctT=0;
		private long incorrectT=0;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
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
	private final long toValue(long kmer, long rkmer, long lengthMask){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return (value&middleMask)|lengthMask;
	}
	

	
	final long rcomp(long kmer, int len){
		return /*amino ? kmer :*/ AminoAcid.reverseComplementBinaryFast(kmer, len);
	}
	
	final boolean passesSpeed(long key){
		return speed<1 || ((key&Long.MAX_VALUE)%17)>=speed;
	}
	
	final boolean failsSpeed(long key){
		return speed>0 && ((key&Long.MAX_VALUE)%17)<speed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** For calculating kmer cardinality */
	final CardinalityTracker loglog;
	final CardinalityTracker loglogOut;
	
	/** Has this class encountered errors while processing? */
	public boolean errorState=false;
	
	/** Fraction of available memory preallocated to arrays */
	private double preallocFraction=0.5;
	/** Initial size of data structures */
	private int initialSize=-1;
	/** Default initial size of data structures */
	private static final int initialSizeDefault=128000; //123
	
	
//	asdf
	//This needs to be redone to better track ref (file) num and scaf num.
	
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	private final AbstractKmerTable[] keySets;
	/** A scaffold's name is stored at scaffoldNames.get(id).
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	private final ArrayList<String> scaffoldNames=new ArrayList<String>();
	/** Names of reference files (refNames[0] is valid). */
	private final ArrayList<String> refNames=new ArrayList<String>();
	/** Number of scaffolds per reference. */
	private final int[] refScafCounts;
	/** scaffoldCounts[id] stores the number of reads with kmer matches to that scaffold */
	private AtomicLongArray scaffoldReadCounts;
	/** scaffoldFragCounts[id] stores the number of fragments (reads or pairs) with kmer matches to that scaffold */
	private AtomicLongArray scaffoldFragCounts;
	/** scaffoldBaseCounts[id] stores the number of bases with kmer matches to that scaffold */
	private AtomicLongArray scaffoldBaseCounts;
	/** scaffoldAmbigReadCounts[id] stores the number of reads ambiguously matching that scaffold */
	private AtomicLongArray scaffoldAmbigReadCounts;
	/** Set to false to force threads to share atomic counter arrays. */
	private boolean ALLOW_LOCAL_ARRAYS=true;
	/** scaffoldLengths[id] stores the length of that scaffold */
	private IntList scaffoldLengths=new IntList();
	/** scaffoldLengths[id] stores the number of kmers in that scaffold (excluding mutants) */
	private IntList scaffoldKmers=new IntList();
	/** scaffolds[id] stores the number of kmers in that scaffold */
	private ArrayList<byte[]> scaffolds=new ArrayList<byte[]>();
	/** Array of reference files from which to load kmers */
	private ArrayList<String> ref=new ArrayList<String>();
	/** Array of literal strings from which to load kmers */
	private String[] literal=null;
	/** Taxonomic tree */
	private TaxTree tree;
	
	/** Input reads */
	private String in1=null, in2=null;
	/** Input qual files */
	private String qfin1=null, qfin2=null;
	/** Output reads (matched and at least minlen) */
	private String outm1=null, outm2=null;
	/** Output reads (unmatched or shorter than minlen) */
	private String outu1=null, outu2=null;
	/** Per-sequence or per-reference output pattern */
	private String outpattern=null;
	/** Statistics output files */
	private String outstats=null, outrpkm=null, outrefstats=null, outtax=null;
	/** Barcode output files */
	private String outbarcodes=null;
	/** NCBI file mapping gi numbers to taxa IDs (gi_taxid_nucl.dmp) */
	private String giTableFile=null;
	/** NCBI file of taxonomy names (names.dmp) */
	private String taxNameFile=null;
	/** NCBI file of taxonomic tree (nodes.dmp) */
	private String taxNodeFile=null;
	/** File containing a serialized TaxTree */
	private String taxTreeFile;
	
	/** Store reference sequences */
	private boolean storeRefBases=false;
	/** Only look for fully-contained reference sequences */
	private boolean processContainedRef=false;
	
	/** Dump kmers here. */
	private String dump=null;
	
	/** Maximum input reads (or pairs) to process.  Does not apply to references.  -1 means unlimited. */
	private long maxReads=-1;
	/** Process this fraction of input reads. */
	private float samplerate=1f;
	/** Set samplerate seed to this value. */
	private long sampleseed=-1;
	
	/** Output reads in input order.  May reduce speed. */
	private final boolean ordered;
	/** Make the middle base in a kmer a wildcard to improve sensitivity */
	private boolean maskMiddle=true;
	private int midMaskLen=0;
	
	/** Store reference kmers with up to this many substitutions */
	private int hammingDistance=0;
	/** Search for query kmers with up to this many substitutions */
	private int qHammingDistance=0;
	/** Store reference kmers with up to this many edits (including indels) */
	private int editDistance=0;
	/** Always skip this many kmers between used kmers when hashing reference. */
	private int refSkip=0;

	private long taxNodeCountLimit=1;
	private int taxNodeNumberLimit=-1;
	private int taxNodeMinLevel=0;
	private int taxNodeMaxLevel=TaxTree.parseLevel("domain");
	
	private boolean trackBarcodes=false;
	private BarcodeMappingStats barcodeData=null;
	
	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/
	
	long readsIn=0;
	long fragsIn=0;
	long basesIn=0;
	long readsMatched=0;
	long basesMatched=0;
	long readsUnmatched=0;
	long basesUnmatched=0;
	long readsAmbig=0;//TODO
	
	long readsQTrimmed=0;
	long basesQTrimmed=0;
	long readsFTrimmed=0;
	long basesFTrimmed=0;
	long readsQFiltered=0;
	long basesQFiltered=0;
	
	long refReads=0;
	long refBases=0;
	long refKmers=0;
	
	long correctReads=0;
	long incorrectReads=0;
	
	long storedKmers=0;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Correct errors via read overlap */
	private final boolean ecc;
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	/** Don't allow a read 'N' to match a reference 'A'.
	 * Reduces sensitivity when hdist>0 or edist>0.  Default: false. */
	private final boolean forbidNs;
	/** AND bitmask with 0's at the middle base */
	private final long middleMask;
	/** Data structure to use.  Default: ARRAYHF */
	private final int tableType;
	
	/** Normal kmer length */
	private final int k;
	/** k-1; used in some expressions */
	private final int k2;
	/** A read must share at least this many kmers to be considered a match.  Default: 1 */
	private final int minKmerHits;
	/** A read must share at least this fraction of its kmers to be considered a match.  Default: 0 */
	private final float minKmerFraction;
	/** Determines how to handle ambiguously-mapping reads */
	private final int ambigMode;
	/** Determines when to early-exit kmer matching */
	private final int matchMode;
	/** First and second must differ by more than this to be unambiguous. */
	private final int clearzone;
	/** First and second must differ by more than this fraction of read kmers to be unambiguous. */
	private final float clearzoneFraction;
	/** Calculate accuracy rate by parsing headers of synthetic reads */
	private final boolean parsecustom;
	
	/** Quality-trim the left side */
	private final boolean qtrimLeft;
	/** Quality-trim the right side */
	private final boolean qtrimRight;
	/** Trim bases at this quality or below.  Default: 4 */
	private final float trimq;
	/** Error rate for trimming (derived from trimq) */
	private final float trimE;
	/** Throw away reads below this average quality after trimming.  Default: 0 */
	private final float minAvgQuality;
	/** If positive, calculate average quality from the first X bases only.  Default: 0 */
	private final int minAvgQualityBases;
	/** Throw away reads failing chastity filter (:Y: in read header) */
	private final boolean chastityFilter;
	/** Throw away reads containing more than this many Ns.  Default: -1 (disabled) */
	private final int maxNs;
	/** Throw away reads containing without at least this many consecutive called bases. */
	private int minConsecutiveBases=0;
	/** Throw away reads shorter than this after trimming.  Default: 10 */
	private final int minReadLength;
	/** Throw away reads longer than this after trimming.  Default: Integer.MAX_VALUE */
	private final int maxReadLength;
	/** Toss reads shorter than this fraction of initial length, after trimming */
	private final float minLenFraction;
	/** Filter reads by whether or not they have matching kmers */
	private final boolean kfilter;
	/** Trim left bases of the read to this position (exclusive, 0-based) */
	private final int forceTrimLeft;
	/** Trim right bases of the read after this position (exclusive, 0-based) */
	private final int forceTrimRight;
	/** Trim this many rightmost bases of the read */
	private final int forceTrimRight2;
	/** Trim right bases of the read modulo this value.
	 * e.g. forceTrimModulo=50 would trim the last 3bp from a 153bp read. */
	private final int forceTrimModulo;
	
	/** If positive, only look for kmer matches in the leftmost X bases */
	private int restrictLeft;
	/** If positive, only look for kmer matches the rightmost X bases */
	private int restrictRight;
	
	/** Skip this many initial input reads */
	private final long skipreads;

	/** Pairs go to outbad if either of them is bad, as opposed to requiring both to be bad.
	 * Default: true. */
	private final boolean removePairsIfEitherBad;
	
	/** Print only statistics for scaffolds that matched at least one read
	 * Default: true. */
	private final boolean printNonZeroOnly;
	
	/** Rename reads to indicate what they matched.
	 * Default: false. */
	private final boolean rename;
	/** Use names of reference files instead of scaffolds.
	 * Default: false. */
	private final boolean useRefNames;
	
	/** Include the count when adding names */
	private boolean addCount=true;
	/** Only use the top hit when adding names */
	private boolean topHitOnly=false;
	/** Add 'unknown' to unmapped or ambiguous reads in rename mode */
	private boolean addUnknown=false;
	
	/** Fraction of kmers to skip, 0 to 15 out of 16 */
	private final int speed;
	
	/** Skip this many kmers when examining the read.  Default 1.
	 * 1 means every kmer is used, 2 means every other, etc. */
	private final int qSkip;
	
	/** True if speed and qSkip are disabled. */
	private final boolean noAccel;
	
	/** Pick a single scaffold per read pair, rather than per read */
	private final boolean keepPairsTogether;
	
	/** Store match IDs in an IntList rather than int array */
	private final boolean USE_COUNTVECTOR;
	
	/** Gather taxonomic information */
	private final boolean USE_TAXTREE;
	
	private ConcurrentHashMap<String, StringCount> setMap;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of tables (and threads, during loading) */
	private static final int WAYS=9; //123
	/** Verbose messages */
	public static final boolean verbose=false; //123

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
	private static final ArrayList<Read> POISON=new ArrayList<Read>(0);
	/** Number of columns for statistics output, 3 or 5 */
	public static int STATS_COLUMNS=5;
	/** Release memory used by kmer storage after processing reads */
	public static boolean RELEASE_TABLES=true;
	/** Max value of hitCount array */
	public static final int HITCOUNT_LEN=1000;
	/** Make unambiguous copies of ref sequences with ambiguous bases */
	public static boolean REPLICATE_AMBIGUOUS=false;
	/** Write refstats in similar style to BBSplit */
	public static boolean BBSPLIT_STYLE=false;
	/** Rename duplicate scaffold names instead of crashing. */
	public static boolean RENAME_DUPLICATES=true;
	/** Create sets of amiguously-mapping gene targets on the fly */
	public static boolean MAKE_GENE_SETS=false;
	
	/** x&clearMasks[i] will clear base i */
	private static final long[] clearMasks;
	/** x|setMasks[i][j] will set base i to j */
	private static final long[][] setMasks;
	/** x&leftMasks[i] will clear all bases to the right of i (exclusive) */
	private static final long[] leftMasks;
	/** x&rightMasks[i] will clear all bases to the left of i (inclusive) */
	private static final long[] rightMasks;
	/** x|kMasks[i] will set the bit to the left of the leftmost base */
	private static final long[] lengthMasks;
	
	private static final int qhdistSizeLimit=10;
	
	public static HashMap<String,String> RQC_MAP=null;

	public static final int AMBIG_ALL=1, AMBIG_FIRST=2, AMBIG_TOSS=3, AMBIG_RANDOM=4;
	public static final int MATCH_ALL=1, MATCH_FIRST=2, MATCH_UNIQUE=3;
	
	/*--------------------------------------------------------------*/
	/*----------------      Static Initializers     ----------------*/
	/*--------------------------------------------------------------*/
	
	static{
		clearMasks=new long[32];
		leftMasks=new long[32];
		rightMasks=new long[32];
		lengthMasks=new long[32];
		setMasks=new long[4][32];
		for(int i=0; i<32; i++){
			clearMasks[i]=~(3L<<(2*i));
		}
		for(int i=0; i<32; i++){
			leftMasks[i]=((-1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			rightMasks[i]=~((-1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			lengthMasks[i]=((1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			for(long j=0; j<4; j++){
				setMasks[(int)j][i]=(j<<(2*i));
			}
		}
	}
	
}
