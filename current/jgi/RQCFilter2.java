package jgi;

import java.io.File;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedHashSet;
import java.util.TimeZone;

import align2.BBMap;
import align2.BBSplitter;
import align2.RefToIndex;
import bloom.PolyFilter;
import cardinality.CardinalityTracker;
import clump.Clumpify;
import clump.KmerSort;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile1;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import hiseq.AnalyzeFlowCell;
import server.ServerTools;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import sketch.SendSketch;
import sketch.SketchObject;
import stream.FASTQ;
import stream.Read;
import stream.ReadInputStream;
import tax.FilterByTaxa;
import tax.GiToTaxid;
import tax.TaxFilter;
import tax.TaxNode;
import tax.TaxTree;
import tracker.ReadStats;

/**
 * Wrapper for several other programs to implement Rolling QC's filter stage.
 * Calls SendSketch, Clumpify, KmerCountExact, BBNorm, BBDuk, BBMap, BBMerge, and SplitNexteraLMP.
 * Trims adapters, removes contaminants, and does quality-trimming.
 * @author Brian Bushnell
 * @date Nov 26, 2013
 *
 */
public class RQCFilter2 {

	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Methods    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Program entrance from command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer
		Timer t=new Timer();
		
		//Create a filter instance
		RQCFilter2 x=new RQCFilter2(args);
		
		//Execute filtering.
		x.process();
		
		//Report time
		t.stop();
		outstream.println("\nOverall Time: \t"+t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	RQCFilter2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), true);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Parses some shared arguments
		Parser parser=new Parser();
		
		//Symbols to insert in output filename to denote operations performed; may be overriden from command line
		String symbols_=null;
		
		boolean doNextera_=false;
		FASTQ.DETECT_QUALITY_OUT=false;
		FASTQ.ASCII_OFFSET_OUT=33;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_BGZIP=ReadWrite.USE_UNBGZIP=ReadWrite.PREFER_BGZIP=true;
		ReadWrite.setZipThreads(Tools.max(Shared.threads()>1 ? 2 : 1, Shared.threads()>20 ? Shared.threads()/2 : Shared.threads()));
		TaxFilter.REQUIRE_PRESENT=false;
		SendSketch.suppressErrors=ServerTools.suppressErrors=true;
		ReadStats.GC_BINS_AUTO=true;
		
		boolean doMerge_=true;
		
		//Parse argument list
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("="); //Expect key=value pairs
			String a=split[0].toLowerCase(); //All keys are converted to lower case
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseCommonStatic(arg, a, b)){
				if(a.equalsIgnoreCase("taxpath")) {
					if(b==null) {
						taxTreeFile=null;
					}else {
						if(TaxTree.TAX_PATH==null) {taxTreeFile="";}
						else {taxTreeFile=TaxTree.TAX_PATH;}
						int lastchar=Tools.lastChar(taxTreeFile);
						if(taxTreeFile.length()>0 && lastchar!='/' && lastchar!='\\') {
							taxTreeFile+="/";
						}
						taxTreeFile+="tree.taxtree.gz";
					}
				}
			}else if(Parser.parseZip(arg, a, b)){
				if(a.equals("pigz")){
					pigz=b;
				}else if(a.equals("unpigz")){
					unpigz=b;
				}else if(a.equals("zl") || a.equals("ziplevel")){
					zl=b;
				}
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				primaryArgList.add(arg);
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
				out1=b;
			}else if(a.equals("out2") || a.equals("output2")){
				out2=b;
			}else if(a.equals("ref")){
				if(b!=null){
					if(!b.contains(",") || new File(b).exists()){
						bbdukFilterRefs.add(b);
					}else{
						String[] split2=b.split(",");
						for(String s2 : split2){
							bbdukFilterRefs.add(s2);
						}
					}
				}
			}else if(a.equals("spikeinref")){
				if(b!=null){
					if(!b.contains(",") || new File(b).exists()){
						spikeinRefs.add(b);
					}else{
						String[] split2=b.split(",");
						for(String s2 : split2){
							spikeinRefs.add(s2);
						}
					}
				}
			}else if(a.equals("artifactdb")){
				mainArtifactFile=b;
			}else if(a.equals("ribodb")){
				riboKmers=b;
			}else if(a.equals("phixref")){
				phixRef=b;
			}else if(a.equals("fragadapter")){
				fragAdapter=b;
			}else if(a.equals("rnaadapter")){
				rnaAdapter=b;
			}else if(a.equals("lfpelinker")){
				lfpeLinker=b;
			}else if(a.equals("cliplinker") || a.equals("jointseq")){
				clipLinker=b;
			}else if(a.equals("clrslinker")){
				clrsLinker=b;
			}else if(a.equals("bisulfite") || a.equals("bisulfate")){ //Allow for a common mispelling... ;)
				bisulfite=Parse.parseBoolean(b);
			}else if(a.equals("trimfragadapter") || a.equals("trimfragadapters")){
				fragAdapterFlag=Parse.parseBoolean(b);
			}else if(a.equals("swift") || a.equals("trimswift")){
				swiftFlag=Parse.parseBoolean(b);
			}else if(a.equals("trimrnaadapter") || a.equals("trimrnaadapters")){
				rnaAdapterFlag=Parse.parseBoolean(b);
			}else if(a.equals("removehuman") || a.equals("human")){
				humanFlag=Parse.parseBoolean(b);
			}else if(a.equals("removedog") || a.equals("dog")){
				dogFlag=Parse.parseBoolean(b);
			}else if(a.equals("removecat") || a.equals("cat")){
				catFlag=Parse.parseBoolean(b);
			}else if(a.equals("removemouse") || a.equals("mouse")){
				mouseFlag=Parse.parseBoolean(b);
			}else if(a.equals("catdoghumanmouse") || a.equals("mousecatdoghuman") || a.equals("")){
				mouseCatDogHumanFlag=Parse.parseBoolean(b);
			}else if(a.equals("keephumanreads") || a.equals("keephuman")){
				keepHumanReads=Parse.parseBoolean(b);
			}else if(a.equals("pairedtosingle")){
				pairedToSingle=Parse.parseBoolean(b);
			}else if(a.equals("aggressive")){
				aggressiveHumanFlag=aggressiveMicrobeFlag=Parse.parseBoolean(b);
			}else if(a.equals("aggressivehuman")){
				aggressiveHumanFlag=Parse.parseBoolean(b);
			}else if(a.equals("aggressivemicrobe")){
				aggressiveMicrobeFlag=Parse.parseBoolean(b);
			}else if(a.equals("removemicrobes") || a.equals("removecommonmicrobes") || a.equals("microbes")){
				removeCommonMicrobeFlag=Parse.parseBoolean(b);
			}else if(a.equals("trimuntrim")){
				trimUntrim=Parse.parseBoolean(b);
			}else if(a.equals("ribomap")){
				riboMapFlag=Parse.parseBoolean(b);
			}else if(a.equals("removechloroplast") || a.equals("chloroplast") || a.equals("chloro") || a.equals("chloromap")){
				chloroMapFlag=Parse.parseBoolean(b);
			}else if(a.equals("removemito") || a.equals("mito") || a.equals("mitomap")){
				mitoMapFlag=Parse.parseBoolean(b);
			}else if(a.equals("detectmicrobes") || a.equals("detectcommonmicrobes")){
				detectMicrobeFlag=Parse.parseBoolean(b);
			}else if(a.equals("detectmicrobes2") || a.equals("detectothermicrobes")){
				otherMicrobeFlag=Parse.parseBoolean(b);
			}else if(a.equals("removeribo") || a.equals("ribo")){
				riboKmerFlag=Parse.parseBoolean(b);
			}else if(a.equals("riboout") || a.equals("outribo")){
				riboOutFile=b;
			}else if(a.equals("fbtout") || a.equals("outfbt")){
				fbtOutFile=b;
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("ml") || a.equals("minlen") || a.equals("minlength")){
				minLen=Integer.parseInt(b);
			}else if(a.equals("mlf") || a.equals("minlenfrac") || a.equals("minlenfraction") || a.equals("minlengthfraction")){
				minLenFraction=Float.parseFloat(b);
			}else if(a.equals("libtype") || a.equals("library")){
				libType=toLibType(b);
			}else if(a.equals("path") || a.equals("outdir")){
				outDir=b;
			}else if(a.equals("symbols")){
				symbols_=b;
			}else if(a.equals("overallstats") || a.equals("stats")){
				rqcStatsName=b;
			}else if(a.equals("stats2")){
				rqcStatsName2=b;
			}else if(a.equals("scafstats")){
				scaffoldStatsName1=b;
			}else if(a.equals("scafstatskt") || a.equals("scafstatstrim")){
				scaffoldStatsName_kt=b;
			}else if(a.equals("refstats")){
				refStatsName=b;
			}else if(a.equals("kmerstats")){
				kmerStatsName1=b;
			}else if(a.equals("log")){
				logName=b;
			}else if(a.equals("donefile")){
				doneFile=b;
			}else if(a.equals("ihist")){
				ihistName=b;
			}else if(a.equals("merge") || a.equals("domerge")){
				doMerge_=Parse.parseBoolean(b);
			}else if(a.equals("khist") || a.equals("dokhist")){
				doKhist=Parse.parseBoolean(b);
			}else if(a.equals("filelist")){
				fileListName=b;
			}else if(a.equals("compress")){
				compress=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("polya") || a.equalsIgnoreCase("removepolya") || a.equalsIgnoreCase("filterpolya")){
				filterPolyA=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("trimpolya")){
				trimPolyA=Parse.parseBoolean(b);
			}else if(a.equals("trimpolyg")){
				trimPolyGLeft=trimPolyGRight=Parser.parsePoly(b);
			}else if(a.equals("trimpolygleft")){
				trimPolyGLeft=Parser.parsePoly(b);
			}else if(a.equals("trimpolygright")){
				trimPolyGRight=Parser.parsePoly(b);
			}else if(a.equals("filterpolyg")){
				filterPolyG=Parser.parsePoly(b);
			}else if(a.equals("phix") || a.equals("removephix")){
				phixFlag=Parse.parseBoolean(b);
			}else if(a.equals("lambda") || a.equals("removelambda")){
				lambdaFlag=Parse.parseBoolean(b);
			}else if(a.equals("pjet")){
				pjetFlag=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("sip")){
				sipFlag=Parse.parseBoolean(b);
			}else if(a.equals("mtst")){
				mtstFlag=Parse.parseBoolean(b);
			}else if(a.equals("kapa")){
				kapaFlag=Parse.parseBoolean(b);
			}else if(a.equals("jointseq")){
				jointSeq=b;
			}else if(a.equals("nextera") || a.equals("nexteralmp")){
				doNextera_=Parse.parseBoolean(b);
			}else if(a.equals("copyundefined") || a.equals("cu")){
				copyUndefined=Parse.parseBoolean(b);
			}else if(a.equals("ktrim")){
				ktrim=b;
			}else if(a.equals("ftl") || a.equals("forcetrimleft")){
				ftl=Integer.parseInt(b);
			}else if(a.equals("ftr") || a.equals("forcetrimright")){
				ftr=Integer.parseInt(b);
			}else if(a.equals("ftr2") || a.equals("forcetrimright2")){
				ftr2=Integer.parseInt(b);
			}else if(a.equals("mink")){
				mink=Integer.parseInt(b);
			}else if(a.equals("k")){
				assert(false) : "To specify kmer length, use filterk, trimk, mapk, or normalizek instead of just 'k'";
				filter_k=Integer.parseInt(b);
			}else if(a.equals("filterk")){
				filter_k=Integer.parseInt(b);
			}else if(a.equals("trimk")){
				trim_k=Integer.parseInt(b);
			}else if(a.equals("mapk")){
				map_k=Integer.parseInt(b);
			}else if(a.equals("normalizek") || a.equals("normk") || a.equals("ecck")){
				normalize_k=Integer.parseInt(b);
			}else if(a.equals("filterhdist")){
				hdist_filter=Integer.parseInt(b);
			}else if(a.equals("spikeink")){
				spikein_k=Integer.parseInt(b);
			}else if(a.equals("spikeinhdist")){
				hdist_spikein=Integer.parseInt(b);
			}else if(a.equals("filterqhdist")){
				qhdist_filter=Integer.parseInt(b);
			}else if(a.equals("trimhdist")){
				hdist_trim=Integer.parseInt(b);
			}else if(a.equals("trimhdist2")){
				hdist2_trim=Integer.parseInt(b);
			}else if(a.equals("ribohdist")){
				hdist_ribo=Integer.parseInt(b);
			}else if(a.equals("riboedist") || a.equals("riboedits")){
				edist_ribo=Integer.parseInt(b);
			}else if(a.equals("maq")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(b.indexOf(',')>-1){
					String[] x=b.split(",");
					assert(x.length==2) : "maq should be length 1 or 2 (at most 1 comma).\nFormat: maq=quality,bases; e.g. maq=10 or maq=10,20";
					minAvgQuality=Byte.parseByte(x[0]);
					minAvgQualityBases=Integer.parseInt(x[1]);
				}else{
					minAvgQuality=Byte.parseByte(b);
				}
			}else if(a.equals("forcetrimmod") || a.equals("forcemrimmodulo") || a.equals("ftm")){
				forceTrimModulo=Integer.parseInt(b);
			}else if(a.equals("trimq")){
				trimq=Byte.parseByte(b);
			}else if(a.equals("qtrim")){
				if(b==null){qtrim="r";}
				else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){qtrim="l";qtrimFlag=true;}
				else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){qtrim="r";qtrimFlag=true;}
				else if(b.equalsIgnoreCase("both") || b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr")){qtrim="lr";qtrimFlag=true;}
				else if(Tools.isDigit(b.charAt(0))){
					trimq=Byte.parseByte(b);
					qtrimFlag=(trimq>=0);
					qtrim=(qtrimFlag ? "lr" : "f");
				}else{
					qtrimFlag=Parse.parseBoolean(b);
					qtrim=""+qtrimFlag;
				}
			}else if(a.equals("optitrim") || a.equals("otf") || a.equals("otm")){
				if(b!=null && (b.charAt(0)=='.' || Tools.isDigit(b.charAt(0)))){
					TrimRead.optimalMode=true;
					TrimRead.optimalBias=Float.parseFloat(b);
					assert(TrimRead.optimalBias>=0 && TrimRead.optimalBias<1);
				}else{
					TrimRead.optimalMode=Parse.parseBoolean(b);
				}
			}else if(a.equals("maxns")){
				maxNs=Integer.parseInt(b);
			}else if(a.equals("usetmpdir")){
				writeTempToTmpdir=Parse.parseBoolean(b);
			}else if(a.equals("tmpdir")){
				tmpDir=b;
				writeTempToTmpdir=(b!=null);
			}else if(a.equals("delete") || a.equals("deletetemp")){
				deleteTemp=Parse.parseBoolean(b);
			}else if(a.equals("humanpath")){
				humanPath=b;
			}else if(a.equals("catpath")){
				catPath=b;
			}else if(a.equals("dogpath")){
				dogPath=b;
			}else if(a.equals("mousepath")){
				dogPath=b;
			}else if(a.equalsIgnoreCase("mouseCatDogHumanPath")){
				mouseCatDogHumanPath=b;
			}else if(a.equals("mapref") || a.equals("maprefs")){
				if(b==null){mappingRefs.clear();}
				else{
					for(String s : b.split(",")){
						mappingRefs.add(s);
					}
				}
			}else if(a.equals("chastityfilter") || a.equals("cf") || a.equals("chastity")){
				chastityfilter=b;
			}else if(a.equals("failnobarcode")){
				failnobarcode=b;
			}else if(a.equals("badbarcodes") || a.equals("barcodefilter")){
				barcodefilter=b;
			}else if(a.equals("barcodes") || a.equals("barcode")){
				barcodes=b;
			}else if(a.equals("extend")){
				extendFlag=Parse.parseBoolean(b);
				extendFlagAuto=false;
			}else if(a.equals("extendauto")){
				extendFlagAuto=Parse.parseBoolean(b);
			}else if(a.equals("taxlist") || a.equals("tax") || a.equals("taxa") || a.equals("taxid")){
				taxList=b;
			}else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
			}else if(a.equals("loadgitable")){
				loadGiTable=Parse.parseBoolean(b);
			}else if(a.equals("gitable")){
				giTable=b;
				loadGiTable=(b!=null);
			}else if(a.equals("taxlevel") || a.equals("tl") || a.equals("level") || a.equals("lv")){
				taxLevel=b;
			}else if(a.equals("microberef")){
				commonMicrobesRef=b;
			}else if(a.equals("microbepath")){
				commonMicrobesPath=b;
			}else if(a.equals("microbebuild")){
				commonMicrobesBuild=Integer.parseInt(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}
			
			else if(a.equals("sketch")){
				sketchFlag=Parse.parseBoolean(b);
			}else if(a.equals("silvalocal") || a.equals("localsilva")){
				SILVA_LOCAL=Parse.parseBoolean(b);
			}else if(a.equals("sketchdb") || a.equals("sketchserver")){
				sketchDB=b;
			}else if(a.equals("sketchreads")){
				sketchReads=Integer.parseInt(b);
			}else if(a.equals("sketchsamplerate")){
				sketchSamplerate=b;
				assert(b!=null);
			}else if(a.equals("sketchminprob")){
				sketchMinProb=b;
				assert(b!=null);
			}
			
			else if(a.equals("clump") || a.equals("clumpify")){
				doClump=Parse.parseBoolean(b);
			}else if(a.equals("recal") || a.equals("recalibrate")){
				recalibrate=Parse.parseBoolean(b);
			}else if(a.equals("quantize")){
				quantize=b;
			}else if(a.equals("phixsam")){
				phixSam=b;
			}else if(a.equals("tiledump")){
				tiledump=b;
			}else if(a.equals("polyfilter")){
				polyfilter=b;
			}else if(a.equals("fbt") || a.equals("filterbytile")){
				doFilterByTile=Parse.parseBoolean(b);
			}else if(a.equals("dedupe")){
				removeDuplicates=Parse.parseBoolean(b);
			}else if(a.equals("alldupes")){
				removeAllDuplicates=Parse.parseBoolean(b);
			}else if(a.equals("opticaldupes") || a.equals("optical")){
				removeOpticalDuplicates=Parse.parseBoolean(b);
			}else if(a.equals("edgedupes")){
				removeEdgeDuplicates=Parse.parseBoolean(b);
			}else if(a.equals("dpasses")){
				duplicatePasses=Integer.parseInt(b);
			}else if(a.equals("dsubs")){
				duplicateSubs=Integer.parseInt(b);
			}else if(a.equals("ddist")){
				duplicateDist=Integer.parseInt(b);
			}else if(a.equals("lowcomplexity")){
				lowComplexity=Parse.parseBoolean(b);
			}else if(a.equals("clumpifygroups") || a.equals("groups")){
				clumpifyGroups=Integer.parseInt(b);
			}else if(a.equals("clumpifytmpdir")){
				clumpifyTmpdir=Parse.parseBoolean(b);
			}
			
			else if(a.equals("entropy") || a.equals("efilter") || a.equals("minentropy") || a.equals("entropyfilter")){
				if(b==null){entropyFilter=true;}
				else if(Tools.isLetter(b.charAt(0))){
					entropyFilter=Parse.parseBoolean(b);
				}else{
					minEntropy=Double.parseDouble(b);
					entropyFilter=minEntropy>0;
				}
			}else if(a.equals("entropyk")){
				entropyk=Integer.parseInt(b);
			}else if(a.equals("entropywindow")){
				entropywindow=Integer.parseInt(b);
			}
			
			else if(a.equalsIgnoreCase("RQCFilterData") || a.equalsIgnoreCase("RQCFilterPath") || a.equalsIgnoreCase("data") || a.equalsIgnoreCase("datapath")){
				RQCFilterData=b;
			}
			
			else if(a.equals("discoveradapters") || a.equals("findadapters") || a.equals("detectadapters")){
				discoverAdaptersFlag=Parse.parseBoolean(b);
			}
			
			else if(a.equals("bloomfilter") || a.equals("bloom")){
				bloomFilter=Parse.parseBoolean(b);
			}else if(a.equals("bloomk") || a.equals("bloomfilterk")){
				bloomFilterK=Integer.parseInt(b);
			}else if(a.equals("bloomhashes") || a.equals("bloomfilterhashes")){
				bloomFilterHashes=Integer.parseInt(b);
			}else if(a.equals("bloomminhits") || a.equals("bloomfilterminhits")){
				bloomFilterMinHits=Integer.parseInt(b);
			}else if(a.equals("bloomfilterminreads") || a.equals("bloomminreads")){
				minReadsToBloomFilter=Parse.parseKMG(b);
			}else if(a.equals("bloomserial") || a.equals("serialbloom")){
				bloomSerial=Parse.parseBoolean(b);
			}
			
			else if(a.equals("dryrun")){
				dryrun=Parse.parseBoolean(b);
			}
			
			else if(a.equals("skipfilter")){
				skipFilter=Parse.parseBoolean(b);
			}
			
			//These no longer do anything
			else if(a.equals("dna")){
				dnaArtifactFlag=Parse.parseBoolean(b);
			}else if(a.equals("rna")){
				rnaArtifactFlag=Parse.parseBoolean(b);
				dnaArtifactFlag=!rnaArtifactFlag; //This line requested by Bryce.
			}
			
			else{
				//Uncaptured arguments are passed to BBDuk
				primaryArgList.add(arg);
			}
		}

		{
			if(RQCFilterData==null){RQCFilterData=".";}
			File f=new File(RQCFilterData);
			assert(f.exists() && f.isDirectory()) : "RQCFilterData path is not set correctly.";

			fixReferencePaths();
		}
		
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		
		//Maintain clumpified order in later stages
		if(doClump){ordered=true;}
		
		doNexteraLMP=doNextera_;
		
		if(removeCommonMicrobeFlag && detectMicrobeFlag){detectMicrobeFlag=false;}
		
//		assert(false) : rnaArtifactFlag+"\n"+primaryArgList+"\n"+libType+"\n"+outDir;
		
		if(writeTempToTmpdir){
			if(tmpDir==null){tmpDir=Shared.tmpdir();}
			if(tmpDir!=null){
				tmpDir=tmpDir.replace('\\', '/');
				if(tmpDir.length()>0 && !tmpDir.endsWith("/")){tmpDir+="/";}
			}
		}else{tmpDir=null;}
		
		if(hdist2_trim<0){hdist2_trim=hdist_trim;}
		
		//Pass overwrite flag to BBDuk
		primaryArgList.add("ow="+overwrite);
		
		if(outDir!=null){
			outDir=outDir.trim().replace('\\', '/');
			if(outDir.length()>0 && !outDir.endsWith("/")){outDir=outDir+"/";}
		}else{outDir="";}
		
		//TODO: Unify Name or File
		{//Prepend output directory to output files
			logName=appendOutDir(logName);
			reproduceName=appendOutDir(reproduceName);
			fileListName=appendOutDir(fileListName);
			ihistName=appendOutDir(ihistName);
			khistName=appendOutDir(khistName);
			peaksName=appendOutDir(peaksName);
			riboOutFile=appendOutDir(riboOutFile);
			fbtOutFile=appendOutDir(fbtOutFile);
			chloroOutFile=appendOutDir(chloroOutFile);
			humanOutFile=appendOutDir(humanOutFile);
			spikeinOutFile=appendOutDir(spikeinOutFile);
			synthOutFile1=appendOutDir(synthOutFile1);
			synthOutFile2=appendOutDir(synthOutFile2);
			microbeOutFile=appendOutDir(microbeOutFile);
			microbeStatsFile=appendOutDir(microbeStatsFile);
			microbesUsed=appendOutDir(microbesUsed);
			chloroStatsFile=appendOutDir(chloroStatsFile);
			sipStatsFile=appendOutDir(sipStatsFile);
			sipScafStats=appendOutDir(sipScafStats);

			cardinalityName=appendOutDir(cardinalityName);
			adaptersOutFile=appendOutDir(adaptersOutFile);
			
			
			phistName=appendOutDir(phistName);
			qhistName=appendOutDir(qhistName);
			bhistName=appendOutDir(bhistName);
			gchistName=appendOutDir(gchistName);
			doneFile=appendOutDir(doneFile);
		}
		
		{//Create unique output file names for first and second trimming passes
			if(rqcStatsName2!=null){
				rqcStatsName2=appendOutDir(rqcStatsName2);
			}
			if(rqcStatsName!=null){
				rqcStatsName_kt=appendOutDir("ktrim_"+rqcStatsName);
				rqcStatsName=appendOutDir(rqcStatsName);
			}
			if(kmerStatsName1!=null){
				kmerStatsName_kt=appendOutDir("ktrim_"+kmerStatsName1);
				kmerStatsName1=appendOutDir(kmerStatsName1);
			}
			if(kmerStatsName2!=null){
				kmerStatsName2=appendOutDir(kmerStatsName2);
			}
			if(scaffoldStatsNameSpikein!=null){
				scaffoldStatsNameSpikein=appendOutDir(scaffoldStatsNameSpikein);
			}
			if(scaffoldStatsName1!=null){
				scaffoldStatsName_kt=appendOutDir("ktrim_"+scaffoldStatsName1);
				scaffoldStatsName1=appendOutDir(scaffoldStatsName1);
			}
			if(scaffoldStatsName2!=null){
				scaffoldStatsName2=appendOutDir(scaffoldStatsName2);
			}
			if(refStatsName!=null){
				refStatsName=appendOutDir(refStatsName);
			}
		}
		
		doSpikein=(mtstFlag || kapaFlag || !spikeinRefs.isEmpty());
		
		//Determine execution path
		//This may be obsolete now
		if(libType==FRAG || ((libType==LFPE && lfpeLinker==null) || (libType==CLIP && clipLinker==null) || (libType==CLRS && clrsLinker==null))){
			doAdapterTrim=(fragAdapterFlag || rnaAdapterFlag);
			doFilter=!skipFilter;
		}else if(libType==LFPE){
			doAdapterTrim=true;
			doFilter=!skipFilter;
		}else if(libType==CLIP){
			doAdapterTrim=true;
			doFilter=!skipFilter;
		}else if(libType==CLRS){
			doAdapterTrim=true;
			doFilter=!skipFilter;
		}else{
			throw new RuntimeException("Unknown library type.");
		}
		
		if(phixSam==null) {recalibrate=false;}
		
		//Use combined index if all organisms are specified (common case)
		if(catFlag && dogFlag && humanFlag && mouseFlag){
			mouseCatDogHumanFlag=true;
		}
		
		//Turn off individual flags if a combined flag is used
		if(mouseCatDogHumanFlag){
			mouseFlag=false;
			catFlag=false;
			dogFlag=false;
			humanFlag=false;
		}
		
		if(dogFlag){mappingRefs.add("path="+dogPath);}
		if(catFlag){mappingRefs.add("path="+catPath);}
		if(mouseFlag){mappingRefs.add("path="+mousePath);}
		if(!doMerge_){ihistName=null;}
		doMerge=(ihistName!=null);
		
		//Set final field 'symbols'
		//This is for those magic letters in the filename
		symbols=(symbols_==null ? abbreviation() : symbols_);
		
		assert(in1!=null) : "No input file specified.";
		
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		
		//Create output filename from input filename if no output filename is specified
		if(out1==null){
			
			File f=new File(in1);
			String name=f.getName();
			rawName=ReadWrite.rawName(name);
			int dot=rawName.lastIndexOf('.');
			if(dot>-1){
				out1=rawName.substring(0, dot)+"."+symbols+rawName.substring(dot)+(compress ? ".gz" : "");
			}else{
				out1=rawName+"."+symbols+".fastq"+(compress ? ".gz" : "");
			}
		}else{
			File f=new File(out1);
			String name=f.getName();
			rawName=ReadWrite.rawName(name);
		}
		
		//Ensure temp files have random names
		tempSalt=KmerNormalize.getSalt(out1, 1);
		clumpPrefix="TEMP_CLUMP_"+tempSalt+"_";
		recalPrefix="TEMP_RECAL_"+tempSalt+"_";
		fbtPrefix="TEMP_FBT_"+tempSalt+"_";
		trimPrefix="TEMP_TRIM_"+tempSalt+"_";
		polyPrefix="TEMP_POLY_"+tempSalt+"_";
		humanPrefix="TEMP_HUMAN_"+tempSalt+"_";
		spikeinPrefix="TEMP_SPIKEIN_"+tempSalt+"_";
		filterPrefix1="TEMP_FILTER1_"+tempSalt+"_";
		filterPrefix2="TEMP_FILTER2_"+tempSalt+"_";
		sipPrefix="TEMP_SIP_"+tempSalt+"_";
		taxaPrefix="TEMP_TAXA_"+tempSalt+"_";
		microbePrefix="TEMP_MICROBE_"+tempSalt+"_";
		riboPrefix="TEMP_RIBO_"+tempSalt+"_";
		chloroPrefix="TEMP_CHLORO_"+tempSalt+"_";
		
		if(mappingRefs.size()>0){
			mappingPrefix=new String[mappingRefs.size()];
			for(int i=0; i<mappingRefs.size(); i++){
				mappingPrefix[i]="TEMP_MAP_"+tempSalt+"_"+i+"_";
			}
		}else{
			mappingPrefix=null;
		}
		
		if(reproduceName!=null){
			writeReproduceHeader(reproduceName, args, overwrite);
		}
	}

	private String appendOutDir(String s){
		if(s==null){return null;}
		if(s.startsWith(outDir)){return s;}
		return outDir+s;
	}

	private String removeOutDir(String s){
		if(s==null){return null;}
		if(outDir==null || outDir.length()==0){return s;}
		if(s.startsWith(outDir)){return s.substring(outDir.length());}
		return s;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Processing Methods      ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Primary method to fully execute the program.
	 */
	public void process(){
		
		if(doneFile!=null){
			try {
				File f=new File(doneFile);
				if(f.exists()){
					f.delete();
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		//Create output directory
		if(outDir!=null && outDir.length()>0){
			File f=new File(outDir);
			if(!f.exists()){
				f.mkdirs();
			}
		}
		
		//Create log file
		if(logName!=null){
			boolean b=Tools.canWrite(logName, overwrite);
			assert(b) : "Can't write to "+logName;
			log("bbtools filter start", false);
		}
		
		//Create file list file
		if(fileListName!=null){
			boolean b=Tools.canWrite(fileListName, overwrite);
			assert(b) : "Can't write to "+fileListName;
			
			StringBuilder sb=new StringBuilder();
			if(!doNexteraLMP){
				if(out1!=null){sb.append("filtered_fastq="+out1).append('\n');}
				if(out2!=null){sb.append("filtered_fastq_2="+out2).append('\n');}
			}
			
			String x=(outDir==null ? "" : outDir);
			int xlen=x.length();
			
			//Determine whether to append the output directory prefix in each case
			if(ihistName!=null){
				sb.append("ihist="+removeOutDir(ihistName)).append('\n');
			}
			if(doKhist){
				if(khistName!=null){
					sb.append("khist="+removeOutDir(khistName)).append('\n');
				}
				if(peaksName!=null){
					sb.append("peaks="+removeOutDir(peaksName)).append('\n');
				}
			}
			if(scaffoldStatsNameSpikein!=null && doSpikein){
				sb.append("spikeinstats="+removeOutDir(scaffoldStatsNameSpikein)).append('\n');
			}
			if(scaffoldStatsName1!=null){
				sb.append("scafstats1="+removeOutDir(scaffoldStatsName1)).append('\n');
			}
			if(scaffoldStatsName2!=null){
				sb.append("scafstats2="+removeOutDir(scaffoldStatsName2)).append('\n');
			}
			if(refStatsName!=null){
				sb.append("refstats="+removeOutDir(refStatsName)).append('\n');
			}
			if(riboKmerFlag && riboOutFile!=null){
				sb.append("ribo="+removeOutDir(riboOutFile)).append('\n');
			}
			if(doFilterByTile && fbtOutFile!=null){
				sb.append("filteredByTile="+removeOutDir(fbtOutFile)).append('\n');
			}
			
			if(sipFlag){
				if(sipStatsFile!=null){
					sb.append("sipSpikeinStats="+removeOutDir(sipStatsFile)).append('\n');
				}
				if(sipScafStats!=null){
					sb.append("sipScafStats="+removeOutDir(sipScafStats)).append('\n');
				}
			}
			
			if((chloroMapFlag || mitoMapFlag || riboMapFlag) && chloroOutFile!=null){
				sb.append("chloroMitoRiboReads="+removeOutDir(chloroOutFile)).append('\n');
			}
			if((chloroMapFlag || mitoMapFlag || riboMapFlag) && chloroStatsFile!=null){
				sb.append("chloroMitoRiboStats="+removeOutDir(chloroStatsFile)).append('\n');
			}
			
			if(removeCommonMicrobeFlag && microbeOutFile!=null){
				sb.append("chaffMicrobeReads="+removeOutDir(microbeOutFile)).append('\n');
				sb.append("microbeStats="+removeOutDir(microbeStatsFile)).append('\n');
			}else if(detectMicrobeFlag && microbeStatsFile!=null){
				sb.append("microbeStats="+removeOutDir(microbeStatsFile)).append('\n');
			}
			if(removeCommonMicrobeFlag && microbesUsed!=null){
				sb.append("microbesUsed="+removeOutDir(microbesUsed)).append('\n');
			}
			
			if(otherMicrobeFlag && microbeStats2File!=null){
				sb.append("microbeStats2="+removeOutDir(microbeStats2File)).append('\n');
			}
			
			if(discoverAdaptersFlag && adaptersOutFile!=null){
				sb.append("adaptersDetected="+removeOutDir(adaptersOutFile)).append('\n');
			}
			
			if(doFilter && synthOutFile1!=null){
				sb.append("chaffSynthReads1="+removeOutDir(synthOutFile1)).append('\n');
			}
			if(doFilter && synthOutFile2!=null){
				sb.append("chaffSynthReads2="+removeOutDir(synthOutFile2)).append('\n');
			}
			if((humanFlag || mouseCatDogHumanFlag) && humanOutFile!=null){
				sb.append("chaffHumanReads="+removeOutDir(humanOutFile)).append('\n');
			}
			
			if(sketchFlag && sketchName!=null){
				sb.append("sketchResult="+removeOutDir(sketchName)).append('\n');
			}

			if(doAdapterTrim){
				if(phistName!=null){
					sb.append("polymerHist="+removeOutDir(phistName)).append('\n');
				}
				if(qhistName!=null){
					sb.append("qualityHist="+removeOutDir(qhistName)).append('\n');
				}
				if(bhistName!=null){
					sb.append("baseFrequencyHist="+removeOutDir(bhistName)).append('\n');
				}
				if(gchistName!=null){
					sb.append("GCHist="+removeOutDir(gchistName)).append('\n');
				}
			}
			
			if(sb.length()>0){
				ReadWrite.writeString(sb, fileListName, false);
			}
		}
		
		Shared.GC_BEFORE_PRINT_MEMORY=true;
		Shared.printMemory();
		
		long lastReadsObserved=-1;
		{
			//Sketching has to be done first because after Clumpify, the entire file would need to be processed,
			//rather than just the first X reads
			if(sketchFlag){
				boolean success=runSketch(in1, sketchName, null);
				Shared.printMemory();
			}
			
			//Calculate number of total steps, to determine when to write to the output directory versus localdisk.
			int step=0;
			final int numSteps=(doClump ? 1 : 0)+(doFilterByTile ? 1 : 0)+(doSpikein ? 1 : 0)+(doFilter ? 2 : 0)+(doAdapterTrim ? 1 : 0)+(doNexteraLMP ? 1 : 0)+(riboKmerFlag ? 1 : 0)+
					((chloroMapFlag || mitoMapFlag || riboMapFlag) ? 1 : 0)+(sipFlag ? 1 : 0)+(recalibrate ? 1 : 0)+(polyfilter==null ? 0 : 1)+
					(removeCommonMicrobeFlag ? 1 : 0)+((humanFlag || mouseCatDogHumanFlag) ? 1 : 0)+mappingRefs.size();
			String inPrefix=null, outPrefix=null;
			
			if(discoverAdaptersFlag){
				discoverAdapters(in1, in2, adaptersOutFile);
				Data.unloadAll();
				Shared.printMemory();
			}
			
			//Clumpification
			if(doClump){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? clumpPrefix : null);
//				outstream.println("Clump. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					ReadWrite.ZIPLEVEL=Tools.min(ReadWrite.ZIPLEVEL, 4);
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				//in1z=bar.fq
				//out1z=bar.fq
				//inPrefix=TEMP_SKETCH_123455y546
				//outPrefix=TEMP_CLUMP_123455y546
				clumpify(in1z, in2z, out1z, out2z, inPrefix, outPrefix, step);
				lastReadsObserved=KmerSort.lastReadsOut;
				
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				if(inPrefix!=null){
					delete(inPrefix, out1z, out2z);//Why is this out, not in?
				}
				
				FASTQ.ASCII_OFFSET=33;
				FASTQ.DETECT_QUALITY=false;
				Data.unloadAll();
				Shared.printMemory();
			}
			
			if(recalibrate) {
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? recalPrefix : null);
//				outstream.println("Trim. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					ReadWrite.ZIPLEVEL=Tools.min(ReadWrite.ZIPLEVEL, 4);
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				calcTrueQuality(phixSam);
				recalibrate(in1z, in2z, out1z, out2z, inPrefix, outPrefix, step);
				lastReadsObserved=BBDuk.lastReadsOut;
				
				//TODO: This junk can go in a function
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				if(inPrefix!=null){
					delete(inPrefix, out1z, out2z);
				}
				
				FASTQ.ASCII_OFFSET=33;
				FASTQ.DETECT_QUALITY=false;
				ReadStats.clear();
				Data.unloadAll();
				Shared.printMemory();
			}
			
			//FilterByTile (Illumina only)
			if(doFilterByTile){
				step++;
				if(dryrun || lastReadsObserved<0 || lastReadsObserved>=400000 || tiledump!=null){
					inPrefix=outPrefix;
					outPrefix=(step<numSteps ? fbtPrefix : null);

					final String in1z, in2z, out1z, out2z;
					if(step==1){
						in1z=in1; in2z=in2;
					}else{
						in1z=stripDirs(out1); in2z=stripDirs(out2);
					}
					if(step>=numSteps){
						ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
						out1z=out1; out2z=out2;
					}else{
						ReadWrite.ZIPLEVEL=Tools.min(ReadWrite.ZIPLEVEL, 4);
						out1z=stripDirs(out1); out2z=stripDirs(out2);
					}

					filterByTile(in1z, in2z, out1z, out2z, fbtOutFile, tiledump, inPrefix, outPrefix, step);
					lastReadsObserved=AnalyzeFlowCell.lastReadsOut;
					
					if(in2!=null && out2==null){
						FASTQ.FORCE_INTERLEAVED=true;
						FASTQ.TEST_INTERLEAVED=false;
					}

					if(inPrefix!=null){
						delete(inPrefix, out1z, out2z);//Why is this out, not in?
					}

					FASTQ.ASCII_OFFSET=33;
					FASTQ.DETECT_QUALITY=false;
				}
				Data.unloadAll();
				Shared.printMemory();
			}
			
			//Adapter trimming
			if(doAdapterTrim){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? trimPrefix : null);
//				outstream.println("Trim. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					ReadWrite.ZIPLEVEL=Tools.min(ReadWrite.ZIPLEVEL, 4);
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				ktrim(in1z, in2z, out1z, out2z, inPrefix, outPrefix, step);
				lastReadsObserved=BBDuk.lastReadsOut;
				ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 6);
				
				//TODO: This junk can go in a function
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				if(inPrefix!=null){
					delete(inPrefix, out1z, out2z);
				}
				
				FASTQ.ASCII_OFFSET=33;
				FASTQ.DETECT_QUALITY=false;
				ReadStats.clear();
				Data.unloadAll();
				Shared.printMemory();
			}
			
			//Adapter trimming
			if(polyfilter!=null){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? polyPrefix : null);
//				outstream.println("Trim. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					ReadWrite.ZIPLEVEL=Tools.min(ReadWrite.ZIPLEVEL, 4);
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				polyfilter(in1z, in2z, out1z, out2z, inPrefix, outPrefix, polyfilter, step);
				lastReadsObserved=PolyFilter.lastReadsOut;
				ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 6);
				
				//TODO: This junk can go in a function
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				if(inPrefix!=null){
					delete(inPrefix, out1z, out2z);
				}
				
				FASTQ.ASCII_OFFSET=33;
				FASTQ.DETECT_QUALITY=false;
				Data.unloadAll();
				Shared.printMemory();
			}
			
			//Spikein filtering
			if(doSpikein){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? spikeinPrefix : null);
//				outstream.println("Filter. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				filterSpikein(in1z, in2z, out1z, out2z, spikeinOutFile, inPrefix, outPrefix, step);
				lastReadsObserved=Seal.lastReadsOut;
				
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				if(step>1){
					delete(inPrefix, out1z, out2z);
				}
				Data.unloadAll();
				Shared.printMemory();
			}
			
			//Synthetic contaminant filtering
			if(doFilter){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? filterPrefix1 : null);
//				outstream.println("Filter. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				filter1(in1z, in2z, out1z, out2z, synthOutFile1, inPrefix, outPrefix, step);
				lastReadsObserved=BBDuk.lastReadsOut;
				
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				if(step>1){
					delete(inPrefix, out1z, out2z);
				}
				Data.unloadAll();
				Shared.printMemory();
			}
			
			//Short synthetic contaminant filtering
			if(doFilter){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? filterPrefix2 : null);
//				outstream.println("Filter. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				filter2(in1z, in2z, out1z, out2z, synthOutFile2, inPrefix, outPrefix, step);
				lastReadsObserved=BBDuk.lastReadsOut;
				
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				if(step>1){
					delete(inPrefix, out1z, out2z);
				}
				Data.unloadAll();
				Shared.printMemory();
			}
			
			if(sipFlag){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? sipPrefix : null);
//				outstream.println("Filter. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				//Filter the file to only include specified organisms in the taxlist.
				final String ref=Data.findPath("?SIP_spikeins.fa.gz");
				
				{
					ArrayList<String> args=new ArrayList<String>();
					args.add("covstats="+sipStatsFile);
					final boolean addToOtherStats=(sipStatsFile==null);
					decontamByMapping(in1z, in2z, out1z, out2z, null, sipScafStats, inPrefix, outPrefix, ref, step, addToOtherStats, args);
					if(!addToOtherStats){
						filterstats.parseSip(sipScafStats);
						assert(filterstats.toString()!=null);
					}
				}
				
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				if(step>1){
					delete(inPrefix, out1z, out2z);
				}
				Data.unloadAll();
				Shared.printMemory();
			}
			
			//Non-specific ribosomal RNA removal via kmers
			if(riboKmerFlag){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? riboPrefix : null);
//				outstream.println("Filter. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				filterRibo(in1z, in2z, out1z, out2z, riboOutFile, inPrefix, outPrefix, step);
				lastReadsObserved=BBDuk.lastReadsOut;
				
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				if(step>1){
					delete(inPrefix, out1z, out2z);
				}
				Data.unloadAll();
				Shared.printMemory();
			}
			
			if(detectMicrobeFlag){
//				assert(false) : inPrefix+", "+outPrefix+", "+in1+", "+out1;
				inPrefix=outPrefix;
				final String in1z, in2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else if(step<numSteps){
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}else{
					inPrefix=null;
					in1z=outDir+stripDirs(out1);
					in2z=(out2==null ? null : outDir+stripDirs(out2));
				}
				
				String ref=commonMicrobesRef;//taxFilter(commonMicrobesRef, "taxa.fa.gz", false, false, true, overwrite, false);
//				assert(false) : in1z+" , "+inPrefix;
				detectCommonMicrobes(in1z, in2z, microbeStatsFile, inPrefix, ref, aggressiveMicrobeFlag, commonMicrobesBuild, commonMicrobesRef, commonMicrobesPath);
				Data.unloadAll();
				Shared.printMemory();
			}
			
			//Microbial contaminant removal
			if(removeCommonMicrobeFlag){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? microbePrefix : null);
//				outstream.println("Filter. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}

				String ref=taxFilter(commonMicrobesRef, "taxa.fa.gz", microbesUsed, taxLevel, false, false, true, overwrite, false);
//				outstream.println("in1z="+in1z+"\nout1z="+out1z+"\ninPrefix="+inPrefix+"\noutPrefix="+outPrefix);
				
				{//Microbial removal uses a kmer length of 13 for seeding
					int oldMapK=map_k;
					map_k=13;
					removeCommonMicrobes(in1z, in2z, out1z, out2z, microbeOutFile, microbeStatsFile, inPrefix, outPrefix, ref, step, aggressiveMicrobeFlag);
					map_k=oldMapK;
				}
				
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				if(step>1){
					delete(inPrefix, out1z, out2z);
				}
				Data.unloadAll();
				Shared.printMemory();
			}
			
			if(otherMicrobeFlag){
//				assert(false) : inPrefix+", "+outPrefix+", "+in1+", "+out1;
				inPrefix=outPrefix;
				final String in1z, in2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else if(step<numSteps){
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}else{
					inPrefix=null;
					in1z=outDir+stripDirs(out1);
					in2z=(out2==null ? null : outDir+stripDirs(out2));
				}
				
				String ref=otherMicrobesRef;//taxFilter(commonMicrobesRef, "taxa.fa.gz", false, false, true, overwrite, false);
//				assert(false) : in1z+" , "+inPrefix;
				detectCommonMicrobes(in1z, in2z, microbeStats2File, inPrefix, ref, aggressiveMicrobeFlag, otherMicrobesBuild, otherMicrobesRef, otherMicrobesPath);
				Data.unloadAll();
				Shared.printMemory();
			}
			
			//Organism-specific chloroplast, mito, and ribo contaminant removal via mapping
			//TODO: Describe ref construction
			if(chloroMapFlag || mitoMapFlag || riboMapFlag){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? chloroPrefix : null);
//				outstream.println("Filter. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				//Disable chloroMapFlag if the organism is not a plant
				if(chloroMapFlag && taxList!=null){
					boolean foundPlant=false;
					TaxTree tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, false);
					TaxNode plant=tree.getNodeByName("Viridiplantae");//Name of plant kingdom
					for(String s : taxList.split(",")){
						TaxNode tn=tree.parseNodeFromHeader(s, true);
						if(tn!=null && (tree.commonAncestor(plant, tn)==plant || tree.commonAncestor(plant, tn)==tn)){
							foundPlant=true;
							break;
						}
					}
					if(!foundPlant){
						outstream.println("Disabled chloroMapFlag because organism is not a plant.");
						chloroMapFlag=false;
					}
				}
				
				//A file containing relevant sequences from all available organisms
				final String ref0;
				if(chloroMapFlag){
					if(mitoMapFlag && riboMapFlag){ref0=chloroMitoRiboRef;}
					else if(mitoMapFlag){ref0=chloroMitoRef;}
					else if(riboMapFlag){ref0=chloroRiboRef;}
					else{ref0=chloroplastRef;}
				}else if(mitoMapFlag){
					if(riboMapFlag){ref0=mitoRiboRef;}
					else{ref0=mitoRef;}
				}else{
					ref0=riboRef;
				}
//				String ref=taxFilter(ref0, "chloroMitoRiboRef.fa.gz", true);
				
				//Filter the file to only include specified organisms in the taxlist.
				String ref;
				{//TODO: Work on this
					ArrayList<String> refList=new ArrayList<String>(3);
					if(chloroMapFlag){refList.add(chloroplastRef);}
					if(mitoMapFlag){refList.add(mitoRef);}
					if(riboMapFlag){refList.add(riboRef);}
					ref=taxFilterList(refList, "chloroMitoRiboRef.fa.gz", true, true, "species");
					if(ref!=null){
						File f=new File(ref);
						if(!f.exists() || f.length()<200){
							outstream.println("Can't find chloroplast/mito/ribo for taxa; using full dataset.");
							ref=ref0;
						}
					}else{ref=ref0;}
				}
				
				{
					final boolean addToOtherStats=(chloroStatsFile==null);
					decontamByMapping(in1z, in2z, out1z, out2z, chloroOutFile, chloroStatsFile, inPrefix, outPrefix, ref, step, addToOtherStats, null);
					if(!addToOtherStats){
						filterstats.parseChloro(chloroStatsFile);
						assert(filterstats.toString()!=null); //Must be here due to an ordering issue with checking
					}
				}
				
				//Delete the reference if it was just created.
				//TODO: Maybe this should be retained
				if(!ref0.equals(ref)){
					delete(null, ref);
				}
				Data.unloadAll();
				Shared.printMemory();
				
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				if(step>1){
					delete(inPrefix, out1z, out2z);
				}
			}
			
			//Human, cat, dog, and mouse removal
			if(humanFlag || mouseCatDogHumanFlag){
				step++;
				inPrefix=outPrefix;
				outPrefix=(step<numSteps ? humanPrefix : null);
//				outstream.println("Human. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				dehumanize(in1z, in2z, out1z, out2z, humanOutFile, inPrefix, outPrefix, step, mouseCatDogHumanFlag, aggressiveHumanFlag, lastReadsObserved);
				
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				Data.unloadAll(); //TODO: See if this needs to be added to other mapping steps
				Shared.printMemory();
				if(step>1){
					delete(inPrefix, out1z, out2z);
				}
			}
			
			//Removal of other assorted reference sequences by mapping
			//Cat, dog, and mouse individually go here
			if(mappingRefs.size()>0){
				for(int i=0; i<mappingRefs.size(); i++){
					step++;
					inPrefix=outPrefix;
					outPrefix=(step<numSteps ? mappingPrefix[i] : null);
					//				outstream.println("Human. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
					
					final String in1z, in2z, out1z, out2z;
					if(step==1){
						in1z=in1; in2z=in2;
					}else{
						in1z=stripDirs(out1); in2z=stripDirs(out2);
					}
					if(step>=numSteps){
						ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
						out1z=out1; out2z=out2;
					}else{
						out1z=stripDirs(out1); out2z=stripDirs(out2);
					}
					
					decontamByMapping(in1z, in2z, out1z, out2z, null, null, inPrefix, outPrefix, mappingRefs.get(i), step, true, null);
					
					if(in2!=null && out2==null){
						FASTQ.FORCE_INTERLEAVED=true;
						FASTQ.TEST_INTERLEAVED=false;
					}
					
					Data.unloadAll();
					Shared.printMemory();
					if(step>1){
						delete(inPrefix, out1z, out2z);
					}
				}
			}
			
			//Nextera LMP library processing
			//TODO: Rename to nexteraLMP
			if(doNexteraLMP){
				step++;
				inPrefix=outPrefix;
				outPrefix=null;
//				outstream.println("Nextera. step="+step+", in="+in1+", out="+out1+", inPrefix="+inPrefix+", outPrefix="+outPrefix);
				
				final String in1z, in2z, out1z, out2z;
				if(step==1){
					in1z=in1; in2z=in2;
				}else{
					in1z=stripDirs(out1); in2z=stripDirs(out2);
				}
				
				if(step>=numSteps){
					ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 9);
					out1z=out1; out2z=out2;
				}else{
					out1z=stripDirs(out1); out2z=stripDirs(out2);
				}
				
				//Insert size calculation
				//In this case it has to be done before SplitNextera
				if(doMerge){merge(in1z, in2z, inPrefix);}
				
				splitNextera(in1z, in2z, inPrefix, outPrefix, step);
				
				if(in2!=null && out2==null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
				}
				
				Data.unloadAll();
				Shared.printMemory();
				if(step>1){
					delete(inPrefix, out1z, out2z);
				}
			}else{
				if(doMerge){//Insert size calculation
					if(step==0){
						merge(in1, in2, null);
					}else{
						merge(out1, out2, null);
					}
				}
				if(doKhist){//TODO: This does not run in Nextera mode
					if(step==0){
						khist(in1, in2, null);
					}else{
						khist(out1, out2, null);
					}
				}
			}
		}
		
		if(doMerge){
			BBDuk.putRqc("outputReads", BBMerge.readsProcessedTotal*2, true, false);
			BBDuk.putRqc("outputBases", BBMerge.basesProcessedTotal, true, false);
		}
		
		//Write combined stats file (number of reads/bases present/removed in each stage)
		if(rqcStatsName!=null){
			final TextStreamWriter tsw=new TextStreamWriter(rqcStatsName, overwrite, false, false);
			tsw.start();
			tsw.println(BBDuk.rqcString());
			tsw.println("gcPolymerRatio="+Tools.format("%.6f", filterstats.gcPolymerRatio));
			tsw.poisonAndWait();
		}
		
		//Write the new combined stats file
		if(rqcStatsName2!=null){
			String s=filterstats.toString();
			ReadWrite.writeString(s, rqcStatsName2, false);
		}
		
		//Finish writing log file
		if(logName!=null){
			log("RQCFilter complete", true);
			if(logName.endsWith(".tmp")){ //Remove .tmp extension
				String old=logName;
				logName=logName.substring(0, logName.length()-4);
				new File(old).renameTo(new File(logName));
			}
		}
		
		if(doneFile!=null){
			Date d=new Date();
			ReadWrite.writeString(in1+"\n"+d+"\n"+d.getTime(), doneFile);
		}
	}
	
	private boolean discoverAdapters(String in1, String in2, String outa){
		log("findAdapters start", true);
		
		ArrayList<String> argList=new ArrayList<String>();
		
		{//Fill list with BBMerge arguments
			argList.add("overwrite="+overwrite);
			
			//Set read I/O files
			if(in1!=null){argList.add("in1="+in1);}
			if(in2!=null){argList.add("in2="+in2);}
			if(outa!=null){argList.add("outa="+outa);}
			argList.add("reads=1m");
		}
		
		String[] mergeargs=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbmerge.sh", mergeargs);
		}
		
//		assert(false) : ReadInputStream.toReads(outa, FileFormat.FA, -1);
		
		boolean success=false;
		if(!dryrun){//run BBMerge
			final boolean OLD_FORCE_INT=FASTQ.FORCE_INTERLEAVED, OLD_TEST_INT=FASTQ.TEST_INTERLEAVED;
			BBMerge merger=new BBMerge(mergeargs);
//			BBMerge merger=new BBMerge(new String[] {"overwrite=true", "in1=random.fq.gz", "outa=rqc4/adaptersDetected.fa"});
//			assert(false) : ReadInputStream.toReads(outa, FileFormat.FA, -1);
			try {
				merger.process();
				FASTQ.FORCE_INTERLEAVED=OLD_FORCE_INT;//These must be reset BEFORE reading the adapters
				FASTQ.TEST_INTERLEAVED=OLD_TEST_INT;
//				assert(false) : ReadInputStream.toReads(outa, FileFormat.FA, -1);
				long merged=merger.correctCountTotal+merger.incorrectCountTotal;
				if(outa!=null && merged>=5000){
					ArrayList<Read> list=ReadInputStream.toReads(outa, FileFormat.FA, -1);
					if(list!=null && list.size()>=2){
						success=true;
						for(int rnum=0; rnum<2 && success; rnum++){
							Read r=list.get(rnum);
							if(r.length()<trim_k){
//								assert(false) : r.length()+", "+trim_k; //123
								success=false; break;
							}
							for(int i=0; i<trim_k && success; i++){
								if(!AminoAcid.isFullyDefined(r.bases[i])){
//									assert(false) : r.length()+", "+trim_k+", "+i+", "+(char)r.bases[i]; //123
									success=false; break;
								}
							}
						}
					}else{
//						assert(false) : (list==null ? null : list.size()+"\n"+outa+"\n"+list+"\n");//123
					}
				}else{
//					assert(false) : merged+", "+outa; //123
				}
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
			FASTQ.FORCE_INTERLEAVED=OLD_FORCE_INT;
			FASTQ.TEST_INTERLEAVED=OLD_TEST_INT;
		}
		
		discoveredAdaptersValid=success;
//		assert(false) : discoveredAdaptersValid;//123
		log("findAdapters finish", true);
		return success;
	}
	
	
	/**
	 * Runs Clumpify for compression.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void clumpify(String in1, String in2, String out1, String out2, String inPrefix, String outPrefix, int stepNum){
		
		log("clumpify start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
		{//Fill list with Clumpify arguments
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			argList.add("reorder");
			
			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("out1="+outPre+out1);}
			if(out2!=null){argList.add("out2="+outPre+out2);}

			//Set dedupe params
			if(removeDuplicates){argList.add("dedupe");}
			if(removeAllDuplicates){argList.add("allduplicates");}
			if(removeOpticalDuplicates || removeEdgeDuplicates){argList.add("optical");}
			if(removeEdgeDuplicates){argList.add("spany"); argList.add("adjacent");}
			if(duplicateSubs>=0){argList.add("dsubs="+duplicateSubs);}
			if(duplicateDist>=0){argList.add("ddist="+duplicateDist);}
			if(duplicatePasses>0){argList.add("passes="+duplicatePasses);}
			if(lowComplexity){argList.add("lowcomplexity="+lowComplexity);}
			if(clumpifyGroups>0){argList.add("groups="+clumpifyGroups);}
			if(clumpifyTmpdir){argList.add("usetmpdir="+clumpifyTmpdir);}
		}
		
		String[] args=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "clumpify.sh", args);
		}
		
		if(!dryrun){//Run Clumpify
			Clumpify c=new Clumpify(args);
			try {
				c.process(new Timer());
				assert(KmerSort.lastReadsOut<=KmerSort.lastReadsIn) : KmerSort.lastReadsIn+", "+KmerSort.lastBasesIn+", "+KmerSort.lastReadsOut+", "+KmerSort.lastBasesOut;
				assert(KmerSort.lastBasesOut<=KmerSort.lastBasesIn) : KmerSort.lastReadsIn+", "+KmerSort.lastBasesIn+", "+KmerSort.lastReadsOut+", "+KmerSort.lastBasesOut;
				assert(KmerSort.lastBasesOut>=KmerSort.lastReadsOut) : KmerSort.lastReadsIn+", "+KmerSort.lastBasesIn+", "+KmerSort.lastReadsOut+", "+KmerSort.lastBasesOut;
				assert(KmerSort.lastBasesIn>=KmerSort.lastReadsIn) : KmerSort.lastReadsIn+", "+KmerSort.lastBasesIn+", "+KmerSort.lastReadsOut+", "+KmerSort.lastBasesOut;
				assert(KmerSort.lastReadsOut>=0);
				assert(KmerSort.lastBasesOut>=0);
				assert(KmerSort.lastReadsOut>0 || KmerSort.lastReadsIn==0);
				assert(KmerSort.lastBasesOut>0 || KmerSort.lastBasesIn==0);
				if(filterstats.readsIn<1){
					filterstats.readsIn=KmerSort.lastReadsIn;
					filterstats.basesIn=KmerSort.lastBasesIn;
				}
				filterstats.readsOut=KmerSort.lastReadsOut;
				filterstats.basesOut=KmerSort.lastBasesOut;
				filterstats.readsDuplicate=KmerSort.lastReadsIn-KmerSort.lastReadsOut;
				filterstats.basesDuplicate=KmerSort.lastBasesIn-KmerSort.lastBasesOut;
				
				assert(filterstats.toString()!=null); //123
				
//				assert(false) : KmerSort.lastReadsIn+", "+KmerSort.lastBasesIn+", "+KmerSort.lastReadsOut+", "+KmerSort.lastBasesOut;
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Optionally append files to file list here
		
		log("clumpify finish", true);
	}
	
	/**
	 * Runs CalcTrueQuality to generate recalibration matrices from a sam file.
	 * @param sam Path to a sam file
	 */
	private void calcTrueQuality(String sam){
		Data.setPath(".");

		log("calcTrueQuality start", true);
		
		String[] args=new String[] {"in="+sam, "ref=phix", "callvars", "usetiles"};
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "calctruequality.sh", args);
		}
		
		if(!dryrun){//Run CalcTrueQuality
			CalcTrueQuality ctq=new CalcTrueQuality(args);
			try {
				ctq.process();
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
			Shared.TRIM_READ_COMMENTS=Shared.TRIM_RNAME=false;
		}
		
		log("calcTrueQuality finish", true);
	}
		
	/**
	 * Runs BBDuk to recalibrate quality scores.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void recalibrate(String in1, String in2, String out1, String out2, 
			String inPrefix, String outPrefix, int stepNum){
		log("recalibrate start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
		{//Fill list with BBDuk arguments
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			if(ordered){argList.add("ordered");}
			if(overwrite) {argList.add("ow");}
			argList.add("usetiles");
			
			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("out1="+outPre+out1);}
			if(out2!=null){argList.add("out2="+outPre+out2);}
		}
		
		String[] args=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbduk.sh", args);
		}
		
		if(!dryrun){//Run BBDuk
			BBDuk duk=new BBDuk(args);
			try {
				duk.process();
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		CalcTrueQuality.unloadMatrices();
		log("recalibrate finish", true);
	}
	
	
	/**
	 * Runs FilterByTile to remove reads from low-quality areas of the flowcell.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param tiledump Dump from full-lane FilterByTile run (optional)
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void filterByTile(String in1, String in2, String out1, String out2, 
			String outbad, String tiledump, String inPrefix, String outPrefix, int stepNum){
		
		log("filterByTile start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
		{//Fill list with FilterByTile arguments
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			if(overwrite) {argList.add("ow");}
			
			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("out1="+outPre+out1);}
			if(out2!=null){argList.add("out2="+outPre+out2);}
			if(outbad!=null){argList.add("outb="+outbad);}
			if(tiledump!=null){argList.add("dumpin="+tiledump);}
		}
		
		String[] args=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "filterbytile.sh", args);
		}
		
		if(!dryrun){//Run FilterByTile
			AnalyzeFlowCell afc=new AnalyzeFlowCell(args);
			try {
				afc.process(new Timer());
				filterstats.readsLowQuality+=afc.readsDiscarded;
				filterstats.basesLowQuality+=afc.basesDiscarded;
				if(filterstats.readsIn<1){
					filterstats.readsIn=afc.readsProcessed;
					filterstats.basesIn=afc.basesProcessed;
				}
				filterstats.readsOut=afc.readsProcessed-afc.readsDiscarded;
				filterstats.basesOut=afc.basesProcessed-afc.basesDiscarded;
				
				assert(filterstats.toString()!=null); //123
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Optionally append files to file list here
		
		log("filterByTile finish", true);
	}
	
	
	/**
	 * Runs BBDuk to perform:
	 * Kmer trimming, short read removal.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void ktrim(String in1, String in2, String out1, String out2, String inPrefix, String outPrefix, int stepNum){
		
		log("ktrim start", true);
		ktrimFlag=true;
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
		{//Fill list with BBDuk arguments
			argList.add("ktrim="+(ktrim==null ? "f" : ktrim));
			if(ordered){argList.add("ordered");}
			if(minLen>0){argList.add("minlen="+minLen);}
			if(minLenFraction>0){argList.add("minlenfraction="+minLenFraction);}
			if((libType!=CLIP)){
				argList.add("mink="+mink);
				if(libType==FRAG && ("r".equalsIgnoreCase(ktrim) || "right".equalsIgnoreCase(ktrim))){
					if(tboFlag){argList.add("tbo");}
					if(tpeFlag){argList.add("tpe");}
				}
				argList.add("rcomp=f");
				argList.add("overwrite="+overwrite);
				argList.add("k="+trim_k);
				argList.add("hdist="+hdist_trim);
				if(hdist2_trim>=0){
					argList.add("hdist2="+hdist2_trim);
				}
				if(forceTrimModulo>0){
					argList.add("ftm="+forceTrimModulo);
				}

				argList.add("pratio=G,C");
				argList.add("plen=20");
				if(phistName!=null){argList.add("phist="+phistName);}
				if(qhistName!=null){argList.add("qhist="+qhistName);}
				if(bhistName!=null){argList.add("bhist="+bhistName);}
				if(gchistName!=null){argList.add("gchist="+gchistName);}
			}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			
			//Pass along uncaptured arguments
			for(String s : primaryArgList){argList.add(s);}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("out1="+outPre+out1);}
			if(out2!=null){argList.add("out2="+outPre+out2);}

//			if(rqcStatsName!=null){al.add("rqc="+rqcStatsName_kt);} //Old style for 2 log files
			if(rqcStatsName!=null){argList.add("rqc=hashmap");}
			if(kmerStatsName_kt!=null){argList.add("outduk="+kmerStatsName_kt);}
			if(scaffoldStatsName_kt!=null){argList.add("stats="+scaffoldStatsName_kt);}
			
			if(copyUndefined){argList.add("cu");}

			argList.add("loglog"); //Cardinality
			argList.add("loglogout");
		}
		
		{//Add BBDuk references
			ArrayList<String> refs=new ArrayList<String>();

			if(libType==FRAG){
				if(discoveredAdaptersValid){argList.add("adapters="+adaptersOutFile);}
				else{
					if(fragAdapterFlag){refs.add(fragAdapter);}
					if(rnaAdapterFlag){refs.add(rnaAdapter);}
				}
			}else if(libType==LFPE){
				refs.add(lfpeLinker);
			}else if(libType==CLIP){
//				refs.add(clipLinker);
				if(clipLinker!=null){
					argList.add("literal="+clipLinker);
					{//Special processing for literal strings of approx 4bp
						String[] split=clipLinker.split(",");
						int min=split[0].length();
						for(String s : split){min=Tools.min(min, s.length());}
						argList.add("k="+min);
						argList.add("mink=-1");
						argList.add("mm=f");
						argList.add("hdist=0");
						argList.add("edist=0");
						argList.add("ktrimexclusive=t");
					}
				}else{
					throw new RuntimeException("Null clip linker.");
				}
			}else if(libType==CLRS){
				refs.add(clrsLinker);
			}else{
				throw new RuntimeException("Unknown library type.");
			}
			
			StringBuilder refstring=new StringBuilder();
			for(String ref : refs){
				if(ref!=null){
					refstring.append(refstring.length()==0 ? "ref=" : ",");
					refstring.append(ref);
				}
			}
			
			if(refstring!=null && refstring.length()>0){
				argList.add(refstring.toString());
			}
		}
		
		String[] dukargs=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbduk.sh", dukargs);
		}
		
		if(!dryrun){//Run BBDuk
			BBDuk duk=new BBDuk(dukargs);
			try {
				duk.process();
				outstream.println(format("Adapter Sequence Removed:", duk.readsIn, duk.readsOut, duk.basesIn, duk.basesOut));
				log("#Input:\t"+duk.readsIn+" reads\t"+duk.basesIn+" bases\t"+duk.loglogIn.cardinality()+" kmers", true, false);
				log("#Remaining:\t"+duk.readsOut+" reads\t"+duk.basesOut+" bases\t"+duk.loglogOut.cardinality()+" kmers", true, false);
				if(duk.errorState){throw new Exception("BBDuk did not finish successfully; see above for the error message.");}

				filterstats.readsFTrimmed=duk.readsFTrimmed;
				filterstats.basesFTrimmed=duk.basesFTrimmed;
				
				filterstats.readsAdapter=duk.readsIn-duk.readsOut;
				filterstats.basesAdapter=duk.basesIn-duk.basesOut-duk.basesFTrimmed;
				if(filterstats.readsIn<1){
					filterstats.readsIn=duk.readsIn;
					filterstats.basesIn=duk.basesIn;
				}
				filterstats.readsOut=duk.readsOut;
				filterstats.basesOut=duk.basesOut;
				
				filterstats.gcPolymerRatio=duk.getPolymerRatio();
				
				assert(filterstats.toString()!=null); //123
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Optionally append files to file list here
		
		log("ktrim finish", true);
	}
	
	
	/**
	 * Runs PolyFilter to remove poly-G reads.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void polyfilter(String in1, String in2, String out1, String out2, 
			String inPrefix, String outPrefix, String polymers, int stepNum){
		log("polyfilter start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
		{//Fill list with BBDuk arguments
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			if(ordered){argList.add("ordered");}
			if(overwrite) {argList.add("ow");}
			argList.add("polymers="+polymers);
			
			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("out1="+outPre+out1);}
			if(out2!=null){argList.add("out2="+outPre+out2);}
		}
		
		String[] args=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "polyfilter.sh", args);
		}
		
		if(!dryrun){//Run BBDuk
			PolyFilter x=new PolyFilter(args);
			try {
				x.process(new Timer());
				filterstats.basesPolyG+=PolyFilter.lastBasesRemoved;
				filterstats.readsPolyG+=PolyFilter.lastReadsRemoved;
				filterstats.basesLowQuality+=PolyFilter.lastBasesRemoved;
				filterstats.readsLowQuality+=PolyFilter.lastReadsRemoved;
				if(filterstats.readsIn<1){
					filterstats.readsIn=x.readsProcessed;
					filterstats.basesIn=x.basesProcessed;
				}
				filterstats.readsOut=x.readsOut;
				filterstats.basesOut=x.basesOut;
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		log("polyfilter finish", true);
	}
	
	/**
	 * Runs Seal to perform:
	 * Spike-in removal and quantification.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void filterSpikein(String in1, String in2, String out1, String out2, String outbad, String inPrefix, String outPrefix,
			int stepNum){
		
		log("filterSpikein start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
//		outstream.println("inPre="+inPre+", outPre="+outPre+", outDir="+outDir+", tmpDir="+tmpDir); //123
		
		{//Fill list with BBDuk arguments
			if(ordered){argList.add("ordered");}
			argList.add("overwrite="+overwrite);
			argList.add("k="+spikein_k);
			argList.add("hdist="+hdist_spikein);
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			
//			//Pass along uncaptured arguments
//			for(String s : primaryArgList){argList.add(s);}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("outu1="+outPre+out1);}
			if(out2!=null){argList.add("outu2="+outPre+out2);}
			if(outbad!=null){argList.add("outm="+outbad);}

//			if(rqcStatsName!=null){al.add("rqc="+rqcStatsName);} //Old style for 2 log files
//			if(rqcStatsName!=null){argList.add("rqc=hashmap");}
			if(scaffoldStatsNameSpikein!=null){argList.add("stats="+scaffoldStatsNameSpikein);}
			
//			if(copyUndefined){argList.add("cu");}
//
			argList.add("loglog"); //Cardinality
			argList.add("loglogout");
		}
		
		{//Add references
			if(mtstFlag){spikeinRefs.add(mtstRef);}
			if(kapaFlag){spikeinRefs.add(kapaRef);}
			
			if(spikeinRefs.isEmpty()){assert(false) : "spikeinFlag is true but there are no spikein references.";}

			StringBuilder refstring=new StringBuilder();
			for(String ref : spikeinRefs){
				if(ref!=null){
					refstring.append(refstring.length()==0 ? "ref=" : ",");
					refstring.append(ref);
				}
			}

			if(refstring!=null && refstring.length()>0){
				argList.add(refstring.toString());
			}
		}
		
		String[] sealargs=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "seal.sh", sealargs);
		}
		
		if(!dryrun){//Run Seal
			Seal seal=new Seal(sealargs);
			try {
				seal.process();
				outstream.println(format("Spikein Sequence Removed:", seal.readsIn, seal.readsIn-seal.readsMatched, seal.basesIn, seal.basesIn-seal.basesMatched));
				log("#Input:\t"+seal.readsIn+" reads\t"+seal.basesIn+" bases\t"+seal.loglog.cardinality()+" kmers", true, false);
				log("#Remaining:\t"+seal.readsUnmatched+" reads\t"+seal.basesUnmatched+" bases\t"+seal.loglogOut.cardinality()+" kmers", true, false);
				
				filterstats.readsSpikin=seal.readsMatched;
				filterstats.basesSpikin=seal.basesMatched;
				if(filterstats.readsIn<1){
					filterstats.readsIn=seal.readsIn;
					filterstats.basesIn=seal.basesIn;
				}
				filterstats.readsOut=seal.readsIn-seal.readsMatched;
				filterstats.basesOut=seal.basesIn-seal.basesMatched;
				
				assert(filterstats.toString()!=null); //123
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Optionally append files to file list here
		
		log("filterSpikein finish", true);
	}
	
	/**
	 * Runs BBDuk to perform:
	 * Quality filtering, quality trimming, n removal, short read removal, artifact removal (via kmer filtering), phiX removal, lambda removal.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void filter1(String in1, String in2, String out1, String out2, String outbad, String inPrefix, String outPrefix,
			int stepNum){
		
		log("filter start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
//		outstream.println("inPre="+inPre+", outPre="+outPre+", outDir="+outDir+", tmpDir="+tmpDir); //123
		
		{//Fill list with BBDuk arguments
			if(swiftFlag){argList.add("swift");}
			if(minAvgQuality>-1){argList.add("maq="+minAvgQuality+","+minAvgQualityBases);}
			if(qtrim!=null){
				argList.add("trimq="+trimq);
				argList.add("qtrim="+qtrim);
			}
			if(ftl>=0){argList.add("ftl="+ftl);}
			if(ftr>=0){argList.add("ftr="+ftr);}
			if(ftr2>=0){argList.add("ftr2="+ftr2);}
			if(pairedToSingle){argList.add("pairedtosingle");}
			if(ordered){argList.add("ordered");}
			argList.add("overwrite="+overwrite);
			if(maxNs>=0){argList.add("maxns="+maxNs);}
			if(minLen>0){argList.add("minlen="+minLen);}
			if(minLenFraction>0){argList.add("minlenfraction="+minLenFraction);}
			argList.add("k="+filter_k);
			argList.add("hdist="+hdist_filter);
			if(qhdist_filter>0){argList.add("qhdist="+qhdist_filter);}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));

			if(chastityfilter!=null){argList.add("cf="+chastityfilter);}
			if(failnobarcode!=null){argList.add("failnobarcode="+failnobarcode);}
			if(barcodefilter!=null){argList.add("barcodefilter="+barcodefilter);}
			if(barcodes!=null){argList.add("barcodes="+barcodes);}
			
			if(entropyFilter && minEntropy>0){
				argList.add("minentropy="+minEntropy);
				argList.add("entropyk="+entropyk);
				argList.add("entropywindow="+entropywindow);
			}
			
			//Pass along uncaptured arguments
			for(String s : primaryArgList){argList.add(s);}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("out1="+outPre+out1);}
			if(out2!=null){argList.add("out2="+outPre+out2);}
			if(outbad!=null){argList.add("outm="+outbad);}

//			if(rqcStatsName!=null){al.add("rqc="+rqcStatsName);} //Old style for 2 log files
			if(rqcStatsName!=null){argList.add("rqc=hashmap");}
			if(kmerStatsName1!=null){argList.add("outduk="+kmerStatsName1);}
			if(scaffoldStatsName1!=null){argList.add("stats="+scaffoldStatsName1);}
			
			if(copyUndefined){argList.add("cu");}
			
			if(bisulfite){argList.add("ftr2=1");}
			if(trimPolyA){argList.add("trimpolya");}
			
			if(trimPolyGLeft>0){argList.add("trimpolygleft="+trimPolyGLeft);}
			if(trimPolyGRight>0){argList.add("trimpolygright="+trimPolyGRight);}
			if(filterPolyG>0){argList.add("filterpolyg="+filterPolyG);}
			argList.add("maxnonpoly=2");
			
			argList.add("loglog"); //Cardinality
			argList.add("loglogout");
		}
		
		{//Add BBDuk references
			bbdukFilterRefs.add(mainArtifactFile);
			if(!doNexteraLMP){bbdukFilterRefs.add(nexteraLinkerFile);}
			if(filterPolyA){bbdukFilterRefs.add(polyAFile);}

			if(phixFlag){bbdukFilterRefs.add(phixRef);}
			if(lambdaFlag){bbdukFilterRefs.add(lambdaRef);}
			if(pjetFlag){bbdukFilterRefs.add(pjetRef);}

			if(libType==FRAG){

			}else if(libType==LFPE){

			}else if(libType==CLIP){

			}else if(libType==CLRS){

			}else{
				throw new RuntimeException("Unknown library type.");
			}

			StringBuilder refstring=new StringBuilder();
			for(String ref : bbdukFilterRefs){
				if(ref!=null){
					refstring.append(refstring.length()==0 ? "ref=" : ",");
					refstring.append(ref);
				}
			}

			if(refstring!=null && refstring.length()>0){
				argList.add(refstring.toString());
			}
		}
		
		String[] dukargs=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbduk.sh", dukargs);
		}
		
		if(!dryrun){//Run BBDuk
			BBDuk duk=new BBDuk(dukargs);
			try {
				duk.process();
				outstream.println(format("Synthetic Contam Sequence Removed:", duk.readsIn, duk.readsOut, duk.basesIn, duk.basesOut));
				log("#Input:\t"+duk.readsIn+" reads\t"+duk.basesIn+" bases\t"+duk.loglogIn.cardinality()+" kmers", true, false);
				log("#Remaining:\t"+duk.readsOut+" reads\t"+duk.basesOut+" bases\t"+duk.loglogOut.cardinality()+" kmers", true, false);
				if(duk.errorState){throw new Exception("BBDuk did not finish successfully; see above for the error message.");}
				
				
//				System.err.println(filterstats.toString());

				long rRemoved=duk.readsIn-duk.readsOut;
				long bRemoved=duk.basesIn-duk.basesOut;

				filterstats.readsArtifact+=duk.readsKFiltered;
				filterstats.basesArtifact+=duk.basesKFiltered;
				filterstats.readsPolyG+=duk.readsPolyTrimmed;
				filterstats.basesPolyG+=duk.basesPolyTrimmed;
				filterstats.readsLowQuality+=(rRemoved-duk.readsKFiltered);
				filterstats.basesLowQuality+=(bRemoved-duk.basesKFiltered);
				filterstats.readsN+=duk.readsNFiltered;
				filterstats.basesN+=duk.basesNFiltered;
				if(filterstats.readsIn<1){
					filterstats.readsIn=duk.readsIn;
					filterstats.basesIn=duk.basesIn;
				}
				filterstats.readsOut=duk.readsOut;
				filterstats.basesOut=duk.basesOut;
				

//				System.err.println(filterstats.toString());
				
				assert(filterstats.toString()!=null); //123
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Optionally append files to file list here
		
		log("filter finish", true);
	}
	
	/**
	 * Runs BBDuk to perform:
	 * Short contaminant removal.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void filter2(String in1, String in2, String out1, String out2, String outbad, String inPrefix, String outPrefix,
			int stepNum){
		
		log("short filter start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
//		outstream.println("inPre="+inPre+", outPre="+outPre+", outDir="+outDir+", tmpDir="+tmpDir); //123
		
		{//Fill list with BBDuk arguments
			if(ordered){argList.add("ordered");}
			argList.add("overwrite="+overwrite);
			argList.add("k="+(doNexteraLMP ? 19 : 20));
			argList.add("hdist="+hdist_filter);
			if(qhdist_filter>0){argList.add("qhdist="+qhdist_filter);}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			
			if(recalibrate && quantize!=null) {
				argList.add("quantize="+quantize);
			}
			
			//Pass along uncaptured arguments
			for(String s : primaryArgList){argList.add(s);}
			
			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("out1="+outPre+out1);}
			if(out2!=null){argList.add("out2="+outPre+out2);}
			if(outbad!=null){argList.add("outm="+outbad);}
			
//			if(rqcStatsName!=null){argList.add("rqc=hashmap");} //TODO
			if(kmerStatsName2!=null){argList.add("outduk="+kmerStatsName2);}
			if(scaffoldStatsName2!=null){argList.add("stats="+scaffoldStatsName2);}
			
			if(copyUndefined){argList.add("cu");}

			argList.add("loglog"); //Cardinality
			argList.add("loglogout");
		}
		
		{//Add BBDuk references
			argList.add("ref="+shortArtifactFile+(doNexteraLMP ? "" : ","+shortArtifactFile));
		}
		
		String[] dukargs=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbduk.sh", dukargs);
		}
		
		if(!dryrun){//Run BBDuk
			BBDuk duk=new BBDuk(dukargs);
			try {
				duk.process();
				outstream.println(format("Short Synthetic Contam Sequence Removed:", duk.readsIn, duk.readsOut, duk.basesIn, duk.basesOut));
				log("#Input:\t"+duk.readsIn+" reads\t"+duk.basesIn+" bases\t"+duk.loglogIn.cardinality()+" kmers", true, false);
				log("#Remaining:\t"+duk.readsOut+" reads\t"+duk.basesOut+" bases\t"+duk.loglogOut.cardinality()+" kmers", true, false);
				if(duk.errorState){throw new Exception("BBDuk did not finish successfully; see above for the error message.");}
				
				long rRemoved=duk.readsIn-duk.readsOut;
				long bRemoved=duk.basesIn-duk.basesOut;

				filterstats.readsArtifact+=duk.readsKFiltered;
				filterstats.basesArtifact+=duk.basesKFiltered;
				filterstats.readsLowQuality+=(rRemoved-duk.readsKFiltered);
				filterstats.basesLowQuality+=(bRemoved-duk.basesKFiltered);
				filterstats.readsN+=duk.readsNFiltered;
				filterstats.basesN+=duk.basesNFiltered;
				if(filterstats.readsIn<1){
					filterstats.readsIn=duk.readsIn;
					filterstats.basesIn=duk.basesIn;
				}
				filterstats.readsOut=duk.readsOut;
				filterstats.basesOut=duk.basesOut;
				
				assert(filterstats.toString()!=null); //123
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Optionally append files to file list here
		
		log("short filter finish", true);
	}
	
	/**
	 * Runs BBDuk to perform:
	 * Ribosomal read removal.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param outRibo Output for ribosomal reads
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void filterRibo(String in1, String in2, String out1, String out2, String outRibo, String inPrefix, String outPrefix,
			int stepNum){
		
		log("filter ribo start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
//		outstream.println("inPre="+inPre+", outPre="+outPre+", outDir="+outDir+", tmpDir="+tmpDir); //123
		
		{//Fill list with BBDuk arguments
			if(ordered){argList.add("ordered");}
			argList.add("k=31");
			argList.add("ref="+riboKmers);
			if(hdist_ribo>0){argList.add("hdist="+hdist_ribo);}
			if(edist_ribo>0){argList.add("edist="+edist_ribo);}
			
			//Pass along uncaptured arguments
			for(String s : primaryArgList){argList.add(s);}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("out1="+outPre+out1);}
			if(out2!=null){argList.add("out2="+outPre+out2);}
			if(outRibo!=null){argList.add("outm="+outRibo);}

//			if(rqcStatsName!=null){al.add("rqc="+rqcStatsName);} //Old style for 2 log files
//			if(rqcStatsName!=null){argList.add("rqc=hashmap");}
//			if(kmerStatsName!=null){argList.add("outduk="+kmerStatsName);}
//			if(scaffoldStatsName!=null){argList.add("stats="+scaffoldStatsName);}
			argList.add("loglog");
			argList.add("loglogout");
		}
		
		String[] dukargs=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbduk.sh", dukargs);
		}
		
		if(!dryrun){//Run BBDuk
			BBDuk duk=new BBDuk(dukargs);
			try {
				duk.process();
				outstream.println(format("Ribosomal Sequence Removed:", duk.readsIn, duk.readsOut, duk.basesIn, duk.basesOut));
				log("#Input:\t"+duk.readsIn+" reads\t"+duk.basesIn+" bases\t"+duk.loglogIn.cardinality()+" kmers", true, false);
				log("#Remaining:\t"+duk.readsOut+" reads\t"+duk.basesOut+" bases\t"+duk.loglogOut.cardinality()+" kmers", true, false);
				if(duk.errorState){throw new Exception("BBDuk did not finish successfully; see above for the error message.");}
				
				filterstats.readsRiboKmer+=duk.readsKFiltered;
				filterstats.basesRiboKmer+=duk.basesKFiltered;
				if(filterstats.readsIn<1){
					filterstats.readsIn=duk.readsIn;
					filterstats.basesIn=duk.basesIn;
				}
				filterstats.readsOut=duk.readsOut;
				filterstats.basesOut=duk.basesOut;
				
				assert(filterstats.toString()!=null); //123
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Optionally append files to file list here
		
		log("filter ribo finish", true);
	}
	
	private static String toPercent(long numerator, long denominator){
		if(denominator<1){return "0.00%";}
		return Tools.format("%.2f%%", numerator*100.0/denominator);
	}
	
	private static String format(String prefix, long rin, long rout, long bin, long bout){
		long rrmvd=rin-rout;
		long brmvd=bin-bout;
		return prefix+"\t"+rrmvd+" reads ("+toPercent(rrmvd, rin)+")\t"+brmvd+" bases ("+toPercent(brmvd, bin)+")";
	}
	
	
	/**
	 * Runs SplitNexteraLMP.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void splitNextera(String in1, String in2, String inPrefix, String outPrefix, int stepNum){
		
		log("splitNextera start", true);
		splitNexteraFlag=true;
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
		final String lmpName, fragName, unknownName, singletonName;
		final String statsName=outPre+nexteraStats;
		
		int dot=rawName.lastIndexOf('.');
		if(dot>-1){
			lmpName=outPre+rawName.substring(0, dot)+"."+symbols+".lmp"+rawName.substring(dot)+(compress ? ".gz" : "");
			fragName=outPre+rawName.substring(0, dot)+"."+symbols+".frag"+rawName.substring(dot)+(compress ? ".gz" : "");
			unknownName=outPre+rawName.substring(0, dot)+"."+symbols+".unknown"+rawName.substring(dot)+(compress ? ".gz" : "");
			singletonName=outPre+rawName.substring(0, dot)+"."+symbols+".singleton"+rawName.substring(dot)+(compress ? ".gz" : "");
		}else{
			lmpName=outPre+rawName+"."+symbols+".lmp.fastq"+(compress ? ".gz" : "");
			fragName=outPre+rawName+"."+symbols+".frag.fastq"+(compress ? ".gz" : "");
			unknownName=outPre+rawName+"."+symbols+".unknown.fastq"+(compress ? ".gz" : "");
			singletonName=outPre+rawName+"."+symbols+".singleton.fastq"+(compress ? ".gz" : "");
		}
		
		{//Fill list with Nextera arguments
			argList.add("mask");
			argList.add("ow="+overwrite);
			if(minLen>0){argList.add("minlen="+minLen);}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));

			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}

			argList.add("out="+lmpName);
			argList.add("outu="+unknownName);
			argList.add("outf="+fragName);
			argList.add("outs="+singletonName);
			argList.add("stats="+statsName);
		}
		
		String[] splitargs=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "splitnextera.sh", splitargs);
		}
		
		if(!dryrun){//run SplitNexteraLMP
			SplitNexteraLMP split=new SplitNexteraLMP(splitargs);
			try {
				split.process();
				StringBuilder sb=new StringBuilder();
				sb.append("LMP:\t"+split.readsLmp()+" reads\t"+split.basesLmp()+" bases\n");
				sb.append("Frag:\t"+split.readsFrag()+" reads\t"+split.basesFrag()+" bases\n");
				sb.append("Unknown:\t"+split.readsUnk()+" reads\t"+split.basesUnk()+" bases\n");
				sb.append("Single:\t"+split.readsSingle()+" reads\t"+split.basesSingle()+" bases\n");
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		if(fileListName!=null){
			StringBuilder sb=new StringBuilder();
			sb.append("lmp="+lmpName).append('\n');
			sb.append("frag="+fragName).append('\n');
			sb.append("unknown="+unknownName).append('\n');
			sb.append("singleton="+singletonName).append('\n');
			sb.append("nexterastats="+statsName).append('\n');
			
			if(sb.length()>0){
				ReadWrite.writeString(sb, fileListName, true);
			}
		}
		
		log("splitNextera finish", true);
	}
	
	/**
	 * Runs BBMap to perform:
	 * Human contaminant removal.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void dehumanize(String in1, String in2, String out1, String out2, String outbad, String inPrefix, String outPrefix,
			int stepNum, boolean mouseCatDogHuman, boolean aggressive, long lastReadsObserved){
		
		log("dehumanize start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);

		final boolean useBBSplit=mouseCatDogHuman;
		{
			if(deterministic){argList.add("deterministic");}
			argList.add("ordered="+(ordered && false)); //Ordered is too slow here
			argList.add("k="+map_k);
			argList.add("idtag=t");
			argList.add("usemodulo");
			argList.add("printunmappedcount");
			argList.add("ow="+overwrite);
			if(trimUntrim){
				argList.add("qtrim=rl");
				argList.add("trimq=10");
				argList.add("untrim");
			}
			argList.add("kfilter=25");
			argList.add("maxsites=1");
			argList.add("tipsearch="+0);
			
			if(aggressive){
				argList.add("minratio=.75");
				argList.add("maxindel=8");
				argList.add("minhits="+1);
				argList.add("bw=26");
				argList.add("bwr=0.22");
				argList.add("build=2");
				argList.add("ef=0.02");
			}else{
				argList.add("minratio=.9");
				argList.add("maxindel=3");
				argList.add("minhits="+2);
				argList.add("bw=12");
				argList.add("bwr=0.16");
				argList.add("fast="+true);
				argList.add("maxsites2=10");
				argList.add("build=1");
				argList.add("ef=0.03");
			}
			
			final int genomes=Tools.max(1, (mouseCatDogHuman ? 4 : 0)+(humanFlag ? 1 : 0)+(dogFlag ? 1 : 0)+(catFlag ? 1 : 0)+(mouseFlag ? 1 : 0));
			if(bloomFilter && (dryrun || lastReadsObserved<0 || lastReadsObserved>=minReadsToBloomFilter) && Shared.memAvailableAdvanced()>=(genomes*5000000000L)){
				argList.add("bloomfilter");
				argList.add("bloomk="+bloomFilterK);
				argList.add("bloomhashes="+bloomFilterHashes);
				argList.add("bloomminhits="+bloomFilterMinHits);
				if(bloomSerial){argList.add("bloomserial");}
			}
			
			if(outbad!=null){argList.add("outm="+outbad);}
			
			if(mouseCatDogHuman){
				argList.add("path="+mouseCatDogHumanPath);
				if(refStatsName!=null && useBBSplit){argList.add("refstats="+refStatsName);}
			}else{
				if(humanRef==null){
					argList.add("path="+humanPath);
				}else{
					RefToIndex.NODISK=true;
					argList.add("ref="+humanRef);
				}
			}
			
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			if(useBBSplit){argList.add("forcereadonly");}
			
//			//Pass along uncaptured arguments
//			for(String s : primaryArgList){argList.add(s);}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(keepHumanReads){
				if(out1!=null){argList.add("out1="+outPre+out1);}
				if(out2!=null){argList.add("out2="+outPre+out2);}
			}else{
				if(out1!=null){argList.add("outu1="+outPre+out1);}
				if(out2!=null){argList.add("outu2="+outPre+out2);}
			}
			
		}
		
		String[] args=argList.toArray(new String[0]);
		
		if(!dryrun){//Run BBMap
			try {
				if(useBBSplit){
//					RefToIndex.FORCE_READ_ONLY=true;
					BBSplitter.main(args);
					RefToIndex.FORCE_READ_ONLY=false;
				}else{
					BBMap.main(args);
				}

//				System.err.println(filterstats.toString());

//				final long lastReadsUsed=BBMap.lastReadsUsed;
				final long lastReadsIn=BBMap.lastReadsIn;
				final long lastEitherMapped=BBMap.lastEitherMapped;
				final long lastReadsOut=lastReadsIn-lastEitherMapped;
				
//				final long lastBasesUsed=BBMap.lastBasesUsed;
				final long lastBasesIn=BBMap.lastBasesIn;
				final long lastEitherMappedBases=BBMap.lastEitherMappedBases;
				final long lastBasesOut=lastBasesIn-lastEitherMappedBases;

//				System.err.println(lastReadsUsed);
//				System.err.println(lastBasesUsed);
//				System.err.println();
//
//				System.err.println(lastReadsIn);
//				System.err.println(lastBasesIn);
//				System.err.println();
//				
//				System.err.println(lastEitherMapped);
//				System.err.println(lastEitherMappedBases);
//				System.err.println();
//
//				System.err.println(lastReadsOut);
//				System.err.println(lastBasesOut);
//				System.err.println();
				
				if(filterstats.readsIn>0){
					assert(lastReadsIn==filterstats.readsOut) : lastReadsIn+", "+filterstats.readsOut;
					assert(lastBasesIn==filterstats.basesOut) : lastBasesIn+", "+filterstats.basesOut;
				}
				
//				System.err.println(BBMap.lastBothUnmapped);
//				System.err.println(BBMap.lastBothUnmappedBases);
//				System.err.println();
				
//				System.err.println(BBMap.lastReadsPassedBloomFilter);
//				System.err.println(BBMap.lastBasesPassedBloomFilter);
//				System.err.println();
				
				if(refStatsName!=null && useBBSplit){
					filterstats.parseHuman(refStatsName);
				}else{
					filterstats.readsHuman=lastEitherMapped;
					filterstats.basesHuman=lastEitherMappedBases;
				}
				if(filterstats.readsIn<1){
					filterstats.readsIn=lastReadsIn;
					filterstats.basesIn=lastBasesIn;
				}
				filterstats.readsOut=lastReadsOut;
				filterstats.basesOut=lastBasesOut;

//				System.err.println(filterstats.toString());
				
				assert(filterstats.toString()!=null); //123
				
				outstream.println(format("Human Sequence Removed:", lastReadsIn, lastReadsOut,
						lastBasesIn, lastBasesOut));
				log("#Input:\t"+lastReadsIn+" reads\t"+lastBasesIn+" bases", true, false);
				log("#Remaining:\t"+lastReadsOut+" reads\t"+lastBasesOut+" bases", true, false);
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Clear the index
		Data.unloadAll();
		
		//Unset NODISK
		RefToIndex.NODISK=false;
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, useBBSplit ? "bbsplit.sh" : "bbmap.sh", args);
		}
		
		//Optionally append files to file list here
		
		log("dehumanize finish", true);
	}
	
	/** Fasta only */
	public static boolean printFileTaxonomy(String in, String results, String taxTreeFile, boolean ow, boolean app, boolean bestEffort){
		if(results==null){return false;}
		TextFile tf=new TextFile(in);
		TextStreamWriter tsw=new TextStreamWriter(results, ow, app, false);
		tsw.start();
		TaxTree tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, false);
		assert(tree!=null) : taxTreeFile+", "+in+", "+results;
		LinkedHashSet<TaxNode> nodes=new LinkedHashSet<TaxNode>();
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith(">")){
				TaxNode tn=tree.parseNodeFromHeader(line.substring(1), bestEffort);
				nodes.add(tn);
			}
		}
		tf.close();
		for(TaxNode tn : nodes){
			tsw.println(tn.id+"\t"+tn.levelStringExtended(false)+"\t"+tn.name);
		}
		return tsw.poisonAndWait();
	}
	
	/**
	 * Runs FilterByTaxa to remove sequences from a reference.
	 */
	private String taxFilter(String in, String out, String results, String level, boolean include, boolean bestEffort, boolean log, boolean ow, boolean app){
		if(taxList==null){
			printFileTaxonomy(in, results, taxTreeFile, ow, app, true);
			return in;
		}
		if(log){log("taxFilter start", true);}
		
		String temp=(tmpDir==null ? outDir : tmpDir)+taxaPrefix+out;
		ArrayList<String> argList=new ArrayList<String>();
		
		{
			argList.add("names="+taxList);
			argList.add("tree="+taxTreeFile);
			argList.add("level="+level);
			argList.add("in="+in);
			argList.add("out="+temp);
			argList.add("ow="+ow);
			argList.add("append="+app);
			argList.add("include="+include);
			argList.add("besteffort="+bestEffort);
			if(results!=null){argList.add("results="+results);}
		}
		
		if(loadGiTable){
			GiToTaxid.initialize(giTable);
		}
		
		String[] args=argList.toArray(new String[0]);
		final FilterByTaxa fbt;
		{
			final boolean OLD_FORCE_INT=FASTQ.FORCE_INTERLEAVED, OLD_TEST_INT=FASTQ.TEST_INTERLEAVED;
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
			fbt=new FilterByTaxa(args);
			fbt.process(new Timer());
			FASTQ.FORCE_INTERLEAVED=OLD_FORCE_INT;
			FASTQ.TEST_INTERLEAVED=OLD_TEST_INT;
			GiToTaxid.unload();
		}
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "filterbytaxa.sh", args);
		}
		
		if(log){log("taxFilter finish", true);}
		Shared.setBuffers();
		
//		outstream.println("*Returning "+(fbt.basesOut>0 ? temp : null));
		
		return fbt.basesOut>0 ? temp : null;
	}
	
	/**
	 * Runs FilterByTaxa to remove sequences from a reference.
	 */
	private String taxFilterList(ArrayList<String> in, String out, boolean include, boolean bestEffort, String initialLevel){
		String level=(bestEffort && include) ? initialLevel : taxLevel;
		
		if(in.size()<2){return taxFilter(in.get(0), out, null, level, include, bestEffort, true, overwrite, false);}
		
		log("taxFilterList start", true);
		
		for(int i=0; i<in.size(); i++){
			taxFilter(in.get(i), out, null, level, include, bestEffort, false, (i==0 && overwrite), i>0);
		}
		
		String temp=(tmpDir==null ? outDir : tmpDir)+taxaPrefix+out;
		
		return temp;
	}
	
	private void detectCommonMicrobes(final String in1, final String in2, final String scafstats, final String inPrefix, final String ref, boolean aggressive,
			int build, String cmRef, String cmPath){
		
		log("detectCommonMicrobes start", true);
		
		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		
		if(ref==null){
			String skipped="Tax filter removed all ref sequences; skipping microbe detection.";
			outstream.println(skipped);
			log(skipped, true);
			log("detectCommonMicrobes finish", true);
			return;
		}

		ArrayList<String> argList=new ArrayList<String>();
		{
			if(deterministic){argList.add("deterministic");}
			argList.add("k="+map_k);
			argList.add("idtag=t");
			argList.add("printunmappedcount");
			argList.add("ow="+overwrite);
			if(trimUntrim){
				argList.add("qtrim=rl");
				argList.add("trimq=10");
				argList.add("untrim");
			}
			argList.add("ef=0.001");
			if(cmPath!=null && cmRef.equals(ref) && cmRef.startsWith(cmPath)){
				RefToIndex.NODISK=false;
				argList.add("path="+cmPath);
			}else{
				RefToIndex.NODISK=true;
				argList.add("ref="+ref);
			}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			
			if(aggressive){
				argList.add("minid=.85");
				argList.add("maxindel=8");
				argList.add("minhits="+1);
				argList.add("bw=26");
				argList.add("bwr=0.22");
				argList.add("build=2");
				argList.add("tipsearch="+2);
			}else{
				argList.add("minid=.95");
				argList.add("idfilter=.95");
				argList.add("maxindel=3");
				argList.add("minhits="+2);
				argList.add("bw=12");
				argList.add("bwr=0.16");
				argList.add("fast="+true);
				argList.add("maxsites2=10");
				argList.add("build="+build);
				argList.add("tipsearch="+0);
			}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
//			if(outbad!=null){argList.add("outm="+outbad);}
			if(scafstats!=null){argList.add("scafstats="+scafstats);}
		}

		String[] args=argList.toArray(new String[0]);

		if(!dryrun){//Run BBMap
			try {
				BBMap.main(args);
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}

		//Clear the index
		Data.unloadAll();

		//Unset NODISK
		RefToIndex.NODISK=false;

		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbmap.sh", args);
		}

		if(ref!=null && !ref.equals(cmRef)){delete(null, ref);}

		log("detectCommonMicrobes finish", true);
	}
	
	/**
	 * Runs BBMap to perform:
	 * Microbial contaminant removal.
	 */
	private void removeCommonMicrobes(String in1, String in2, String out1, String out2, String outbad, String scafstats, String inPrefix, String outPrefix,
			final String ref, int stepNum, boolean aggressive){
		
		log("removeCommonMicrobes start", true);
		
		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
		if(ref==null){
			String skipped="Tax filter removed all ref sequences; skipping microbe removal.";
			outstream.println(skipped);
			log(skipped, true);

			try {
				if(in1!=null){
					File a=new File(inPre+in1);
					File b=new File(outPre+out1);
					outstream.println("Renaming "+a+" to "+b);
					assert(a.exists()) : a;
					assert(!b.exists() || overwrite) : b;
					a.renameTo(b);
					writeReproduceFile(reproduceName, "mv", new String[] {a.toString(), b.toString()});
				}
				if(in2!=null && out2!=null){
					File a=new File(inPre+in2);
					File b=new File(outPre+out2);
					outstream.println("Renaming "+a+" to "+b);
					assert(a.exists()) : a;
					assert(!b.exists() || overwrite) : b;
					a.renameTo(b);
					writeReproduceFile(reproduceName, "mv", new String[] {a.toString(), b.toString()});
				}
			} catch (Throwable e) {
				outstream.println(e.getMessage());
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
			log("removeCommonMicrobes finish", true);
			return;
		}

		ArrayList<String> argList=new ArrayList<String>();
		{
			if(deterministic){argList.add("deterministic");}
			if(ordered){argList.add("ordered");}
			argList.add("quickmatch");
			argList.add("k="+map_k);
			argList.add("idtag=t");
			argList.add("printunmappedcount");
			argList.add("ow="+overwrite);
			if(trimUntrim){
				argList.add("qtrim=rl");
				argList.add("trimq=10");
				argList.add("untrim");
			}
			argList.add("ef=0.001");
			if(commonMicrobesPath!=null && commonMicrobesRef.equals(ref) && commonMicrobesRef.startsWith(commonMicrobesPath)){
				RefToIndex.NODISK=false;
				argList.add("path="+commonMicrobesPath);
			}else{
				RefToIndex.NODISK=true;
				argList.add("ref="+ref);
			}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			
			if(aggressive){
				argList.add("minid=.85");
				argList.add("maxindel=8");
				argList.add("minhits="+1);
				argList.add("bw=26");
				argList.add("bwr=0.22");
				argList.add("build=2");
				argList.add("tipsearch="+2);
			}else{
				argList.add("minid=.95");
				argList.add("idfilter=.95");
				argList.add("maxindel=3");
				argList.add("minhits="+2);
				argList.add("bw=12");
				argList.add("bwr=0.16");
				argList.add("fast="+true);
				argList.add("maxsites2=10");
				argList.add("build="+commonMicrobesBuild);
				argList.add("tipsearch="+0);
			}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("outu1="+outPre+out1);}
			if(out2!=null){argList.add("outu2="+outPre+out2);}
			if(outbad!=null){argList.add("outm="+outbad);}
			if(scafstats!=null){argList.add("scafstats="+scafstats);}
			//			assert(false) : scafstats+", "+microbeStatsFile;
		}

		String[] args=argList.toArray(new String[0]);

		if(!dryrun){//Run BBMap
			try {
				BBMap.main(args);
				
				final long lastReadsIn=BBMap.lastReadsIn;
				final long lastEitherMapped=BBMap.lastEitherMapped;
				final long lastReadsOut=lastReadsIn-lastEitherMapped;
				
				final long lastBasesIn=BBMap.lastBasesIn;
				final long lastEitherMappedBases=BBMap.lastEitherMappedBases;
				final long lastBasesOut=lastBasesIn-lastEitherMappedBases;
				
				if(filterstats.readsIn>0){
					assert(lastReadsIn==filterstats.readsOut) : lastReadsIn+", "+filterstats.readsOut;
					assert(lastBasesIn==filterstats.basesOut) : lastBasesIn+", "+filterstats.basesOut;
				}
				
				filterstats.readsMicrobe=lastEitherMapped;
				filterstats.basesMicrobe=lastEitherMappedBases;
				if(filterstats.readsIn<1){
					filterstats.readsIn=lastReadsIn;
					filterstats.basesIn=lastBasesIn;
				}
				filterstats.readsOut=lastReadsOut;
				filterstats.basesOut=lastBasesOut;
				
				assert(filterstats.toString()!=null); //123
				
				outstream.println(format("Microbial Sequence Removed:", lastReadsIn, lastReadsOut,
						lastBasesIn, lastBasesOut));
				log("#Input:\t"+lastReadsIn+" reads\t"+lastReadsOut+" bases", true, false);
				log("#Remaining:\t"+lastBasesIn+" reads\t"+lastBasesOut+" bases", true, false);
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}

		//Clear the index
		Data.unloadAll();

		//Unset NODISK
		RefToIndex.NODISK=false;

		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbmap.sh", args);
		}

		if(ref!=null && !ref.equals(commonMicrobesRef)){delete(null, ref);}

		log("removeCommonMicrobes finish", true);
	}
	
	/**
	 * Runs BBMap to perform:
	 * Arbitrary contaminant removal.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param inPrefix Append this prefix to input filenames
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void decontamByMapping(String in1, String in2, String out1, String out2, String outbad, String scafstats, String inPrefix, String outPrefix,
			String ref, int stepNum, boolean addToOtherStats, ArrayList<String> otherArgs){
		
		log("decontamByMapping_"+ref+" start", true);
		assert(ref!=null) : "Reference was null.";
		
		ArrayList<String> argList=new ArrayList<String>();
		
		final String inPre=(inPrefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+inPrefix);
		final String outPre=(outPrefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+outPrefix);
		
		{
			if(deterministic){argList.add("deterministic");}
			if(ordered){argList.add("ordered");}
			argList.add("minratio=.9");
			argList.add("maxindel=3");
			argList.add("fast="+true);
			argList.add("minhits="+2);
			argList.add("tipsearch="+4);
			argList.add("bw=12");
			argList.add("bwr=0.16");
			argList.add("quickmatch");
			argList.add("k="+map_k);
			argList.add("idtag=t");
//			argList.add("usemodulo");
			argList.add("printunmappedcount");
			argList.add("ow="+overwrite);
			if(trimUntrim){
				argList.add("qtrim=rl");
				argList.add("trimq=10");
				argList.add("untrim");
			}
			argList.add("ef=0.03");
			if(ref.startsWith("path=")){
				argList.add(ref);
			}else{
				RefToIndex.NODISK=true;
				argList.add("ref="+ref);
			}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			
//			//Pass along uncaptured arguments
//			for(String s : primaryArgList){argList.add(s);}
			
			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			if(out1!=null){argList.add("outu1="+outPre+out1);}
			if(out2!=null){argList.add("outu2="+outPre+out2);}
			if(outbad!=null){argList.add("outm="+outbad);}
			if(scafstats!=null){argList.add("scafstats="+scafstats);}
//			assert(false) : scafstats+", "+microbeStatsFile;
			
			if(otherArgs!=null){
				argList.addAll(otherArgs);
			}
		}
		
		String[] args=argList.toArray(new String[0]);
		
		if(!dryrun){//Run BBMap
			try {
				BBMap.main(args);

				long lastReadsIn=BBMap.lastReadsIn;
				long lastEitherMapped=BBMap.lastEitherMapped;
				long lastReadsOut=lastReadsIn-lastEitherMapped;
				
				long lastBasesIn=BBMap.lastBasesIn;
				long lastEitherMappedBases=BBMap.lastEitherMappedBases;
				long lastBasesOut=lastBasesIn-lastEitherMappedBases;
				
				outstream.println(format("Other Contam Sequence Removed:", lastReadsIn, lastReadsOut,
						lastBasesIn, lastBasesOut));
				log("#Input:\t"+lastReadsIn+" reads\t"+lastBasesIn+" bases", true, false);
				log("#Remaining:\t"+lastReadsOut+" reads\t"+lastBasesOut+" bases", true, false);
				
//				System.err.println(filterstats.toString()); //123
//				
//				System.err.println(
//						lastReadsIn+", "+lastEitherMapped+", "+lastReadsOut+"\n"
//						+lastBasesIn+", "+lastEitherMappedBases+", "+lastBasesOut
//					);
				
				if(addToOtherStats){
					filterstats.readsOther+=lastEitherMapped;
					filterstats.basesOther+=lastEitherMappedBases;
				}
				if(filterstats.readsIn<1){
					filterstats.readsIn=lastReadsIn;
					filterstats.basesIn=lastBasesIn;
				}
				filterstats.readsOut=lastReadsOut;
				filterstats.basesOut=lastBasesOut;
				
				assert(!addToOtherStats || filterstats.toString()!=null); //123
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Clear the index
		Data.unloadAll();
		
		//Unset NODISK
		RefToIndex.NODISK=false;
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbmap.sh", args);
		}
		
		//Optionally append files to file list here
		
		log("decontamByMapping_"+ref+" finish", true);
	}
	
	private boolean runSketch(String in1, String outFname, String outPrefix){
		final String outPre=(outPrefix==null ? outDir : outDir+outPrefix);
		final String fname=outPre+outFname;
		
		boolean success=false;
		String[] split=new String[] {"refseq"};
		if(sketchDB!=null){split=sketchDB.split(",");}
		for(int i=0; i<split.length; i++){
			success=(runSketch_inner(in1, fname, split[i], i>0))|success;
		}
		return success;
	}
	
	private boolean runSketch_inner(String in, String fname, String db, boolean append){
		SketchObject.reset();
		ArrayList<String> argList=new ArrayList<String>();
		

		String name=in;
		try {
			if(taxList!=null){
				TaxTree tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, false);
				TaxNode tn=tree.parseNodeFromHeader(taxList.split(",")[0], true);
				if(tn!=null){name=tn.name;}
			}
		} catch (Exception e1) {}
		
		argList.add("in="+in);
		argList.add("out="+fname);
		argList.add("reads="+sketchReads);
		argList.add("name="+name);
		argList.add("fname="+in1);
		if(sketchMinProb!=null){argList.add("minprob="+sketchMinProb);}
		if(sketchSamplerate!=null){argList.add("samplerate="+sketchSamplerate);}
		if(sketchMerge){argList.add("merge");}
		argList.add("printname0=f");
		argList.add("records=20");
		argList.add("color=false");
		argList.add("depth");
		argList.add("depth2");
		argList.add("unique2");
		argList.add("unique3");
		argList.add("volume");
		argList.add("sortbyvolume");
		argList.add("contam2=genus");
		argList.add("minprob=0.6");
		argList.add("refhits");
		argList.add(db);
		
		//Special handling per server
		if("silva".equalsIgnoreCase(db)){
			argList.add("translate=f");
			if(SILVA_LOCAL){
				argList.add("local");
				argList.add("sizemult=1");
			}else{
				argList.add("sizemult=10");
			}
			argList.add("minkeycount=2");
		}else if("prokprot".equalsIgnoreCase(db) || "protein".equalsIgnoreCase(db)){
			argList.add("translate=t");
			argList.add("sizemult=2.2");
		}else if("refseq".equalsIgnoreCase(db)){
			argList.add("translate=f");
			argList.add("sizemult=2.5");
		}else if("nt".equalsIgnoreCase(db)){
			argList.add("translate=f");
			argList.add("sizemult=1.2");
		}
		
		if(append){argList.add("append");} //True because db name gets written first.
		else{argList.add("ow");}

		String[] argArray=argList.toArray(new String[0]);
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "sendsketch.sh", argArray);
		}
		
		try {
//			ReadWrite.writeString("DB:\t"+db, fname, append); //Already printed in query header
			
			if(!dryrun){SendSketch.main(argArray);}
			return true;
		} catch (Throwable e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return false;
		}
	}
	
	/**
	 * Runs BBMerge to generate an insert size histogram.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param prefix Append this prefix to input filenames
	 */
	private void merge(String in1, String in2, String prefix){
		
		log("merge start", true);
		
		ArrayList<String> argList=new ArrayList<String>();
		
		final String inPre=(prefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+prefix);
		
		if(extendFlagAuto){
			final long cardinality=CardinalityTracker.lastCardinality;
			final long capacity=kmerCapacity(20, true);
			outstream.println("cardinality="+cardinality+", capacity="+capacity);

			extendFlag=(cardinality>0 && cardinality*2<capacity);
		}
		
		{//Fill list with BBMerge arguments
			if(mergeStrictness!=null){argList.add(mergeStrictness);}
			argList.add("overwrite="+overwrite);
			
			//Set read I/O files
			if(in1!=null){argList.add("in1="+inPre+in1);}
			if(in2!=null){argList.add("in2="+inPre+in2);}
			
			if(ihistName!=null){argList.add("ihist="+ihistName);}
			if(cardinalityName!=null){argList.add("outc="+cardinalityName);}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			argList.add("zl="+(zl==null ? ""+ReadWrite.ZIPLEVEL : zl));
			if(discoveredAdaptersValid){argList.add("adapters="+adaptersOutFile);}
			else if(fragAdapter!=null){argList.add("adapters="+fragAdapter);}
			argList.add("mininsert=25");
			
			if(extendFlag){
				argList.add("ecct");
//				argList.add("extend2=20");
//				argList.add("iterations=10");
				argList.add("extend2=100");
				argList.add("rem");
				argList.add("k=62");
				argList.add("prefilter");
				argList.add("prealloc=0.9");
				System.gc();
			}
			argList.add("loglog");
		}
		
		String[] mergeargs=argList.toArray(new String[0]);
		
		if(reproduceName!=null){
			writeReproduceFile(reproduceName, "bbmerge.sh", mergeargs);
		}
		
		if(!dryrun){//run BBMerge
			final boolean OLD_FORCE_INT=FASTQ.FORCE_INTERLEAVED, OLD_TEST_INT=FASTQ.TEST_INTERLEAVED;
			BBMerge merger=new BBMerge(mergeargs);
			try {
				merger.process();
				log("#Input:\t"+(merger.readsProcessedTotal*2)+" reads\t"+merger.basesProcessedTotal+" bases\t"+merger.loglog.cardinality()+" kmers", true, false);
				log("#Remaining:\t"+(merger.readsProcessedTotal*2)+" reads\t"+merger.basesProcessedTotal+" bases\t"+merger.loglog.cardinality()+" kmers", true, false);
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
			FASTQ.FORCE_INTERLEAVED=OLD_FORCE_INT;
			FASTQ.TEST_INTERLEAVED=OLD_TEST_INT;
		}
		
		//Optionally append files to file list here
		
		log("merge finish", true);
	}
	
	
	/**
	 * Runs BBNorm or KmerCountExact to generate a kmer frequency histogram.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param prefix Append this prefix to input filenames
	 */
	private void khist(String in1, String in2, String prefix){
		
		log("khist start", true);
		
		ArrayList<String> argList=new ArrayList<String>();

		final String inPre=(prefix==null ? outDir : (tmpDir==null ? outDir : tmpDir)+prefix);
		
		final long cardinality=CardinalityTracker.lastCardinality;
		final long capacity=kmerCapacity(12, true);
		outstream.println("cardinality="+cardinality+", capacity="+capacity);
		
		if(cardinality<1 || cardinality*1.5>capacity){ //Too many kmers for exact counts; use BBNorm
			{//Fill list with BBNorm arguments
				argList.add("overwrite="+overwrite);

				//Set read I/O files
				if(in1!=null){argList.add("in1="+inPre+in1);}
				if(in2!=null){argList.add("in2="+inPre+in2);}

				if(khistName!=null){argList.add("khist="+khistName);}
				if(peaksName!=null){argList.add("peaks="+peaksName);}
				if(unpigz!=null){argList.add("unpigz="+unpigz);}
				argList.add("keepall");
				argList.add("prefilter");
				argList.add("passes=1");
				argList.add("bits=16");
				argList.add("minprob=0");
				argList.add("minqual=0");
				argList.add("histcolumns=2");
			}

			String[] khistargs=argList.toArray(new String[0]);

			if(reproduceName!=null){
				writeReproduceFile(reproduceName, "khist.sh", khistargs);
			}

			if(!dryrun){//run KmerNormalize
				try {
					KmerNormalize.main(khistargs);
				} catch (Exception e) {
					e.printStackTrace();
					log("failed", true);
					System.exit(1);
				}
			}
		}else{
			{//Fill list with KmerCountExact arguments
				argList.add("overwrite="+overwrite);
				
				//Set read I/O files
				if(in1!=null){argList.add("in1="+inPre+in1);}
				if(in2!=null){argList.add("in2="+inPre+in2);}

				if(khistName!=null){argList.add("khist="+khistName);}
				if(peaksName!=null){argList.add("peaks="+peaksName);}
				if(unpigz!=null){argList.add("unpigz="+unpigz);}
				argList.add("gchist");
				
				if(cardinality*4>capacity){
					argList.add("prealloc");
				}
			}

			String[] khistargs=argList.toArray(new String[0]);

			if(reproduceName!=null){
				writeReproduceFile(reproduceName, "kmercountexact.sh", khistargs);
			}

			if(!dryrun){//run KmerCountExact
				try {
					KmerCountExact.main(khistargs);
				} catch (Exception e) {
					e.printStackTrace();
					log("failed", true);
					System.exit(1);
				}
			}
		}
		
		//Optionally append files to file list here
		
		log("khist finish", true);
	}
	
	private static long kmerCapacity(int bytesPerKmer, boolean prealloc){
		System.gc();
		long memory=Runtime.getRuntime().maxMemory();
		double xmsRatio=Shared.xmsRatio();
		long usableMemory=(long)Tools.max(((memory-96000000)*(xmsRatio>0.97 ? 0.82 : 0.72)), memory*0.45);
		long tableMemory=(long)(usableMemory*.95);
		long estimatedKmerCapacity=(long)((tableMemory*1.0/bytesPerKmer)*(prealloc ? 0.9 : 0.6));
		return estimatedKmerCapacity;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Note - Seung-jin needs to be notified when new calls are added. */
	private void log(String message, boolean append){log(message, append, true);}
	
	/**
	 * Log a message in the log file
	 * @param message Message to log
	 * @param append True to append, false to overwrite
	 */
	private void log(String message, boolean append, boolean printTime){
		if(logName!=null){
			ReadWrite.writeString(message+(printTime ? ", "+timeString() : "")+"\n", logName, append);
		}
	}
	
	
	/**
	 * Delete all non-null filenames.
	 * @param prefix Append this prefix to filenames before attempting to delete them
	 * @param names Filenames to delete
	 */
	private void delete(String prefix, String...names){
		if(!deleteTemp){return;}
		log("delete temp files start", true);
		if(names!=null){
			final String pre=(prefix==null ? "" : (tmpDir==null ? outDir : tmpDir)+prefix);
			for(String s : names){
				if(s!=null){
					s=pre+s;
					if(verbose){outstream.println("Trying to delete "+s);}
					File f=new File(s);
					if(f.exists() || dryrun){
						if(!dryrun){f.delete();}
						writeReproduceFile(reproduceName, "rm", new String[] {s});
					}
				}
			}
		}
		log("delete temp files finish", true);
	}
	
	/**
	 * @return String of symbols indicating which processes were applied to the input reads
	 */
	private String abbreviation(){
		StringBuilder sb=new StringBuilder();
		
		if(mainArtifactFile!=null){sb.append("a");}
		
		if(maxNs>=0){sb.append("n");}
//		if(qtrim!=null && !qtrim.equalsIgnoreCase("f") && !qtrim.equalsIgnoreCase("false")){sb.append("q");}
		if(minAvgQuality>0){sb.append("q");}
		
		if(rnaArtifactFlag){sb.append("r");}
		if(dnaArtifactFlag){sb.append("d");}
		
		if(libType==CLIP){sb.append("c");}
		else if(libType==LFPE){sb.append("l");}
		else if(libType==CLRS){sb.append("s");}

		if(phixFlag){sb.append("p");}
		if(humanFlag || mouseCatDogHumanFlag){sb.append("h");}

//		if(ktrimFlag){sb.append("k");}
		
//		if(doTrim){sb.append("k");}
//		if(qtrimFlag){sb.append("t");}
		
		if(doAdapterTrim || qtrimFlag){sb.append("t");}
		
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * TODO:  Some machines are set to UTC rather than PST
	 * @return Timestamp in RQC's format
	 */
	public static String timeString(){
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
//		sdf.setTimeZone(TimeZone.getTimeZone("PST"));
		sdf.setTimeZone(TimeZone.getDefault());
		return sdf.format(new Date());
	}
	
	public static String stripDirs(String fname){
		return ReadWrite.stripPath(fname);
	}
	
	/**
	 * Set permissions on these files to 777
	 * @param names List of filenames
	 */
	private static void setPermissions(String...names){
		if(names==null){return;}
		for(String name : names){
			if(name!=null && name.trim().length()>0 && new File(name).exists()){
				ReadWrite.setPermissions(name, true, true, true, false);
			}
		}
	}
	
	/**
	 * Writes a single command to the reproduce file
	 * @param fname Filename to write, including path
	 * @param command Command to add to file
	 * @param args Arguments to the command
	 */
	private static void writeReproduceFile(String fname, String command, String[] args){
		StringBuilder sb=new StringBuilder();
		sb.append('\n');
		sb.append(command);
		if(args!=null){
			for(String s : args){
				sb.append(' ').append(s);
			}
		}
		sb.append('\n');
		ReadWrite.writeString(sb, fname, true);
	}
	
	/**
	 * Writes the header for the reproduce file
	 * @param fname Filename to write, including path
	 * @param command Command to add to file
	 * @param args Arguments to the command
	 * @param overwrite Permission to overwrite
	 */
	private static void writeReproduceHeader(String fname, String[] args, boolean overwrite){
		StringBuilder sb=new StringBuilder();
		boolean b=Tools.canWrite(fname, overwrite);
		assert(b) : "Can't write to "+fname;
		sb.append("#!/bin/bash\n");
		sb.append("#BBTools version "+Shared.BBTOOLS_VERSION_STRING+"\n");
		sb.append("#The steps below recapitulate the output of RQCFilter2 when run like this:\n");
		sb.append("#rqcfilter2.sh");
		if(args!=null){
			for(String s : args){
				sb.append(' ').append(s);
			}
		}
		sb.append('\n');
		sb.append("#Data dependencies are available at http://portal.nersc.gov/dna/microbial/assembly/bushnell/RQCFilterData.tar");
		sb.append('\n');
		sb.append('\n');
		ReadWrite.writeString(sb, fname, false);
	}
	
	/**
	 * @param s String representation of library type
	 * @return Numeric code for library type
	 */
	private static int toLibType(String s){
		if(s==null){return FRAG;}
		s=s.trim().toLowerCase();
		if(s.equals("lfpe")){return LFPE;}
		if(s.equals("clip")){return CLIP;}
		if(s.equals("clrs")){return CLRS;}
		if(s.equals("frag") || s.equals("fragment")){return FRAG;}
		throw new RuntimeException("Unknown library type "+s);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Screen output */
	private static PrintStream outstream=System.err;

	/** Synthetic contaminant filtering */
	private final boolean doFilter;
	/** Forcibly override synthetic contaminant filtering */
	private boolean skipFilter=false;
	/** Spikein filtering */
	private final boolean doSpikein;
	/** Clumpify */
	private boolean doClump=false;
	/** QScore recalibration */
	private boolean recalibrate=true;
	/** QScore quantizer */
	private String quantize="2";
	/** FilterByTile */
	private boolean doFilterByTile=true;
	/** Polyfilter */
	private String polyfilter="GC";
	/** Adapter-trimming */
	private final boolean doAdapterTrim;
	/** Run BBMerge for insert size calculation */
	private final boolean doMerge;
	/** Run KmerNormalize for kmer histogram generation */
	private boolean doKhist=false;
	/** Do NexteraLMP splitting */
	private final boolean doNexteraLMP;
	
	/** Symbols to insert in output filename to denote operations performed */
	private final String symbols;
	
	/** Name of raw input file, minus directory and file extension */
	private final String rawName;
	
	/** Type of library; controls processing methods and references to use */
	private int libType=FRAG;
	/** Legacy; affects output file name */
	private boolean rnaArtifactFlag=false;
	/** Legacy; affects output file name */
	private boolean dnaArtifactFlag=true;
	/** True if phix should be filtered out */
	private boolean phixFlag=true;
	/** True if Lambda should be filtered out by kmer-matching */
	private boolean lambdaFlag=true;
	/** True if pjet should be filtered out */
	private boolean pjetFlag=true;
	/** True if SIP spikeins should be filtered out */
	private boolean sipFlag=false;

	/** True if poly-A tails should be right-trimmed */
	private boolean trimPolyA=false;//TODO
	/** True if poly-A should be filtered out */
	private boolean filterPolyA=false;
	
	/** Positive if poly-G should be left-trimmed during filter1 */
	public int trimPolyGLeft=6;
	/** Positive if poly-G should be right-trimmed during filter1 */
	public int trimPolyGRight=6;
	/** Positive if poly-G prefixes should be filtered during filter1 */
	public int filterPolyG=0;

	/** True if mtst should be quantified and filtered out */
	private boolean mtstFlag=false;
	/** True if kapa should be quantified and filtered out */
	private boolean kapaFlag=true;
	
	/** Enables tbo during adapter trimming */
	private boolean tboFlag=true;
	/** Enables tpe during adapter trimming */
	private boolean tpeFlag=true;
	
	/** Unused */
	private String jointSeq=null;
	/** Toss reads shorter than this */
	private int minLen=25;
	/** Toss reads shorter than this fraction of initial length, after trimming */
	private float minLenFraction=0.333f;
	/** Trim bases at this quality or below */
	private float trimq=10;
	/** Throw away reads below this average quality after trimming.  Default: 5 */
	private float minAvgQuality=5;
	/** If positive, calculate the average quality from the first X bases. */
	private int minAvgQualityBases=0;
	/** Trim reads to be equal to 0 modulo this value.  Mainly for 151, 251, and 301bp runs. */
	private int forceTrimModulo=5;
	/** Quality-trimming mode */
	private String qtrim="f";//"rl";
	/** Kmer-trimming mode */
	private String ktrim="r";
	/** Kmer length to use for filtering */
	private int filter_k=31;
	/** Kmer length to use for trimming */
	private int trim_k=23;
	/** Kmer length to use for normalization and error-correction */
	private int normalize_k=31;
	/** Kmer length to use for mapping */
	private int map_k=14;
	/** Shortest kmer to use for trimming */
	private int mink=11;
	/** Throw away reads containing more than this many Ns.  Default: 0 (toss reads with any Ns) */
	private int maxNs=0;
	/** Use this Hamming distance when kmer filtering */
	private int hdist_filter=1;
	/** Use this query Hamming distance when kmer filtering */
	private int qhdist_filter=0;
	/** Use this Hamming distance when kmer trimming */
	private int hdist_trim=1;
	/** Use this Hamming distance when kmer trimming with short kmers */
	private int hdist2_trim=-1;
	/** Use this Hamming distance when kmer trimming with short kmers */
	private int hdist_ribo=0;
	/** Use this Hamming distance when kmer trimming with short kmers */
	private int edist_ribo=0;
	
	/** Kmer length to use for filtering spikeins */
	private int spikein_k=31;
	/** Use this Hamming distance when kmer filtering spikeins */
	private int hdist_spikein=0;
	
	/** Force-trim left during filter1 stage */
	private int ftl=-1;
	/** Force-trim right during filter1 stage */
	private int ftr=-1;
	/** Force-trim right 2 during filter1 stage */
	private int ftr2=-1;
	
	/** Merge strictness: strict, normal, loose, vloose */
	private String mergeStrictness="loose";
	
	/** Trim Truseq and Nextera adapters from right side of reads */
	private boolean fragAdapterFlag=true;
	/** Trim Truseq-RNA adapters from right side of reads */
	private boolean rnaAdapterFlag=false;
	/** Trim Swift low-complexity sequence */
	private boolean swiftFlag=false;
	
	/** Trim 1bp from right side after adapter-trimming/ */
	private boolean bisulfite=false;

	/** Performed quality-trimming on reads */
	private boolean qtrimFlag=false;
	/** Performed kmer-trimming on reads */
	private boolean ktrimFlag=false;
	/** Performed nextera splitting on reads */
	private boolean splitNexteraFlag=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Mapping Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Remove reads mapping to human with high identity */
	private boolean humanFlag=false;
	/** Remove reads mapping to dog with high identity */
	private boolean dogFlag=false;
	/** Remove reads mapping to cat with high identity */
	private boolean catFlag=false;
	/** Remove reads mapping to mouse with high identity */
	private boolean mouseFlag=false;
	/** Remove mouse, cat, dog, and human reads at the same time with BBSplit. */
	private boolean mouseCatDogHumanFlag=false;
	/** Perform cat, dog, mouse, and human removal aggressively, using unmasked genomes. */
	private boolean aggressiveHumanFlag=false;
	/** Report, but do not remove, cat/dog/mouse/human sequence */
	private boolean keepHumanReads=false;

	/** Perform microbe removal aggressively, using unmasked genomes. */
	private boolean aggressiveMicrobeFlag=false;
	/** Remove reads from common microbial contaminants with BBMap */
	private boolean removeCommonMicrobeFlag=false;
	/** Detect but do not remove reads from common microbial contaminants with BBMap */
	private boolean detectMicrobeFlag=false;
	/** Detect but do not remove reads from other microbial contaminants with BBMap */
	private boolean otherMicrobeFlag=false;

	/** Screen reads with a Bloom filter prior to mapping */
	private boolean bloomFilter=true;
	/** Bloom filter kmer length */
	private int bloomFilterK=29;
	/** Bloom filter hashes per kmer */
	private int bloomFilterHashes=1;
	/** Min consecutive hits to consider a read as matching */
	private int bloomFilterMinHits=6;
	/** Don't use Bloom filter for libraries under this size */
	private long minReadsToBloomFilter=4000000;
	/** Use the serialized Bloom filter */
	private boolean bloomSerial=true;
	
	private boolean trimUntrim=true;
	
	/*--------------------------------------------------------------*/
	/*----------------           BBMerge            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Use BBMerge to determine the adapter sequence, and use that instead of the default adapters */
	private boolean discoverAdaptersFlag=true;
	/** True if it looks like the discovered adapters are usable */
	private boolean discoveredAdaptersValid=false;
	
	/** Extend reads to merge longer inserts */
	private boolean extendFlag=false;
	/** Set extendFlag based on the number of unique kmers */
	private boolean extendFlagAuto=true;
	
	/*--------------------------------------------------------------*/
	/*----------------       Mito Ribo Chloro       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Remove ribosomal reads via kmer-matching */
	private boolean riboKmerFlag=false;
	/** Remove reads from specific ribosomal sequence with BBMap */
	private boolean riboMapFlag=false;
	/** Remove reads from specific chloroplast sequence with BBMap */
	private boolean chloroMapFlag=false;
	/** Remove reads from specific mitochondrial sequence with BBMap */
	private boolean mitoMapFlag=false;
	
	/*--------------------------------------------------------------*/
	/*----------------            Sketch            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Query the sketch server */
	private boolean sketchFlag=true;
	/** Use local flag for Silva */
	private boolean SILVA_LOCAL=true;
	/** Sketch server to query */
	private String sketchDB="nt,refseq,silva,prokprot";
	
	/** Estimate kmer cardinality */
	private boolean doCardinality=true;

	
	private boolean verbose=false;
	private boolean overwrite=true;
	private boolean compress=true;

	private boolean copyUndefined=false;
	
	/** Write temp files to $TMPDIR (localdisk) */
	private boolean writeTempToTmpdir=true;
	
	/** Captures the command line "pigz" flag */
	private String pigz="t";
	/** Captures the command line "unpigz" flag */
	private String unpigz="t";
	/** Captures the command line "zl" flag */
	private String zl;
	
	/** Mode for processing chastity flag in Illumina read names */
	private String chastityfilter="t";
	/** Consider the absence of a barcode to mean failure */
	private String failnobarcode=null;
	/** May be set to true, false, or crash to determine how to handle reads with no barcode */
	private String barcodefilter="f";
	/** An optional list of literal barcodes that are allowed */
	private String barcodes=null;
	
	/** Arguments to pass to BBDuk */
	private ArrayList<String> primaryArgList=new ArrayList<String>();
	/** References to pass to BBDuk for artifact removal */
	private ArrayList<String> bbdukFilterRefs=new ArrayList<String>();
	/** References to pass to Seal for spikein removal */
	private ArrayList<String> spikeinRefs=new ArrayList<String>();
	/** References to pass to BBMap for contaminant removal */
	private ArrayList<String> mappingRefs=new ArrayList<String>();
	
	/** List of taxa to NOT map against */
	private String taxList=null;
	/** Taxonomic level for filtering */
	private String taxLevel="order";
	/** Only needed if there are gi numbers in the references */
	private boolean loadGiTable=false;
	
	/** Number of reads used for SendSketch */
	private int sketchReads=2000000;
	private String sketchMinProb="0.2";
	private String sketchSamplerate="1.0";
	private boolean sketchMerge=true;
	
	private boolean entropyFilter=false;
	private double minEntropy=0.42;
	private int entropyk=2;
	private int entropywindow=40;
	
	/*--------------------------------------------------------------*/
	/*----------     Clumpify Deduplication Parameters    ----------*/
	/*--------------------------------------------------------------*/

	private boolean removeDuplicates=false;
	private boolean removeAllDuplicates=false;
	private boolean removeOpticalDuplicates=false;
	private boolean removeEdgeDuplicates=false;
	private int duplicatePasses=1;
	private int duplicateSubs=-1;
	private int duplicateDist=-1;
	private boolean lowComplexity=false;
	private int clumpifyGroups=-1;
	private boolean clumpifyTmpdir=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Read Data Files       ----------------*/
	/*--------------------------------------------------------------*/

	private final String tempSalt;

	private final String clumpPrefix;
	private final String recalPrefix;
	private final String fbtPrefix;
	private final String trimPrefix;
	private final String polyPrefix;
	private final String humanPrefix;
	private final String spikeinPrefix;
	private final String filterPrefix1;
	private final String filterPrefix2;
	private final String sipPrefix;
	private final String taxaPrefix;
	private final String microbePrefix;
	private final String riboPrefix;
	private final String chloroPrefix;
	private final String[] mappingPrefix;
	
	/** Directory in which to write all files */
	private String outDir="";
	
	/** Directory in which to write all temp files */
	private String tmpDir=Shared.tmpdir();
	
	/** Primary input reads file (required) */
	private String in1=null;
	/** Secondary input reads file */
	private String in2=null;
	/** Primary output reads file (required) */
	private String out1=null;
	/** Secondary output reads file */
	private String out2=null;
	
	private String phixSam=null;
	private String tiledump=null;

	private boolean deleteTemp=true;
	
	private boolean deterministic=true;
	
	private boolean ordered=false;
	
	private boolean dryrun=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Separated Reads       ----------------*/
	/*--------------------------------------------------------------*/

	private String riboOutFile="ribo.fq.gz";
	private String fbtOutFile="filteredByTile.fq.gz";
	private String chloroOutFile="chloro.fq.gz";
	private String humanOutFile="human.fq.gz";
	private String spikeinOutFile="spikein.fq.gz";
	private String synthOutFile1="synth1.fq.gz";
	private String synthOutFile2="synth2.fq.gz";
	private String microbeOutFile="microbes.fq.gz";
	
	/*--------------------------------------------------------------*/
	/*----------------           Log Files          ----------------*/
	/*--------------------------------------------------------------*/
	
	private RQCFilterStats filterstats=new RQCFilterStats();
	
	private String logName="status.log";
	private String reproduceName="reproduce.sh";
	private String fileListName="file-list.txt";

	private String rqcStatsName="filterStats.txt";
	private String rqcStatsName2="filterStats2.txt";
	private String kmerStatsName1="kmerStats1.txt";
	private String kmerStatsName2="kmerStats2.txt";
	private String scaffoldStatsName1="scaffoldStats1.txt";
	private String scaffoldStatsName2="scaffoldStats2.txt";
	private String scaffoldStatsNameSpikein="scaffoldStatsSpikein.txt";
	private String refStatsName="refStats.txt";
	private String microbeStatsFile="commonMicrobes.txt";
	private String microbeStats2File="otherMicrobes.txt";
	private String adaptersOutFile="adaptersDetected.fa";
	private String microbesUsed="microbesUsed.txt";
	private String chloroStatsFile="chloroStats.txt";
	private String sipStatsFile="sipSpikeinStats.txt";
	private String sipScafStats="sipScafStats.txt";
	private String nexteraStats="nexteraStats.txt";
	private String ihistName="ihist_merge.txt";
	private String khistName="khist.txt";
	private String peaksName="peaks.txt";
	private String phistName="phist.txt";
	private String qhistName="qhist.txt";
	private String bhistName="bhist.txt";
	private String gchistName="gchist.txt";
	private String sketchName="sketch.txt";
	
	private String cardinalityName="cardinality.txt";
	
	/** ktrim phase rqc stats file */
	private String rqcStatsName_kt;
	/** ktrim phase stats file */
	private String kmerStatsName_kt;
	/** ktrim phase scaffold stats file */
	private String scaffoldStatsName_kt;
	
	private String doneFile="RQCFilter_finished.txt";
	
	/*--------------------------------------------------------------*/
	/*----------------        Reference Files       ----------------*/
	/*--------------------------------------------------------------*/
	
	//TODO: Need to be backed up to tape
	private void fixReferencePaths(){

		if(RQCFilterData==null){RQCFilterData=".";}
		
		shortArtifactFile=fixPath(shortArtifactFile);
		mainArtifactFile=fixPath(mainArtifactFile);
		polyAFile=fixPath(polyAFile);
		nexteraLinkerFile=fixPath(nexteraLinkerFile);
		
		phixRef=fixPath(phixRef);
		lambdaRef=fixPath(lambdaRef);
		lfpeLinker=fixPath(lfpeLinker);
		clrsLinker=fixPath(clrsLinker);
		
		pjetRef=fixPath(pjetRef);
		riboKmers=fixPath(riboKmers);
		fragAdapter=fixPath(fragAdapter);
		rnaAdapter=fixPath(rnaAdapter);
		mtstRef=fixPath(mtstRef);
		kapaRef=fixPath(kapaRef);
		
		humanPath=fixPath(humanPath);
		catPath=fixPath(catPath);
		dogPath=fixPath(dogPath);
		mousePath=fixPath(mousePath);
		mouseCatDogHumanPath=fixPath(mouseCatDogHumanPath);
		humanRef=fixPath(humanRef);
		
		chloroplastRef=fixPath(chloroplastRef);
		mitoRef=fixPath(mitoRef);
		chloroMitoRef=fixPath(chloroMitoRef);
		chloroMitoRiboRef=fixPath(chloroMitoRiboRef);
		mitoRiboRef=fixPath(mitoRiboRef);
		chloroRiboRef=fixPath(chloroRiboRef);

		riboRef=fixPath(riboRef);
		
		commonMicrobesPath=fixPath(commonMicrobesPath);
		commonMicrobesRef=fixPath(commonMicrobesRef);

		otherMicrobesPath=fixPath(otherMicrobesPath);
		otherMicrobesRef=fixPath(otherMicrobesRef);
		
		taxTreeFile=fixPath(taxTreeFile);
	}
	
	private String fixPath(String s){
		if(s==null){return null;}
		s=s.replaceAll("RQCFILTER_PATH", RQCFilterData);
		File f=new File(s);
		assert(f.exists() && f.canRead()) : 
			"Cannot read file "+s+"\nPlease be sure RQCFilterData is set correctly.";
		return s;
	}
	
	private static String getRQCPath() {
		if(Shared.PERLMUTTER) {
			return "/global/cfs/cdirs/bbtools/RQCFilterData_Local";
		}else if(Shared.DORI) {
			return "/clusterfs/jgi/groups/gentech/genome_analysis/ref/RQCFilterData";
		}else if(Shared.AWS) {
			return "";
		}else {
			return "";
		}
	}
	
	private String RQCFilterData=getRQCPath();
	
	private String shortArtifactFile = "RQCFILTER_PATH/short.fa.gz";
	private String mainArtifactFile = "RQCFILTER_PATH/Illumina.artifacts.fa.gz";
	private String polyAFile = "RQCFILTER_PATH/polyA.fa.gz";
	private String nexteraLinkerFile = "RQCFILTER_PATH/nextera_LMP_linker.fa.gz";
	
	private String phixRef = "RQCFILTER_PATH/phix2.fa.gz";
	private String lambdaRef = "RQCFILTER_PATH/lambda.fa.gz";
	private String lfpeLinker = "RQCFILTER_PATH/lfpe.linker.fa.gz";
	private String clrsLinker = "RQCFILTER_PATH/crelox.fa.gz";
	private String clipLinker = clipLinkerDefault; //A literal string; "CATG" is supposed to be the normal linker.

//	private String sipRef = "RQCFILTER_PATH/pJET1.2.fa.gz";
	private String pjetRef = "RQCFILTER_PATH/pJET1.2.fa.gz";
	private String riboKmers = "RQCFILTER_PATH/riboKmers20fused.fa.gz";
	private String fragAdapter = "RQCFILTER_PATH/adapters2.fa.gz";
	private String rnaAdapter = "RQCFILTER_PATH/truseq_rna.fa.gz";
	private String mtstRef = "RQCFILTER_PATH/mtst.fa.gz";
	private String kapaRef = "RQCFILTER_PATH/kapatags.L40.fa.gz";
	
	private String humanPath = "RQCFILTER_PATH/human_genome/";
	private String catPath = "RQCFILTER_PATH/cat_genome/";
	private String dogPath = "RQCFILTER_PATH/dog_genome/";
	private String mousePath = "RQCFILTER_PATH/mouse_genome/";
	private String mouseCatDogHumanPath = "RQCFILTER_PATH/mousecatdoghuman/";
	private String humanRef = null;
	
	private String chloroplastRef = "RQCFILTER_PATH/plastid.fa.gz";
	private String mitoRef = "RQCFILTER_PATH/mitochondrion.fa.gz";
	private String chloroMitoRef = "RQCFILTER_PATH/chloro_mito.fa.gz";
	private String chloroMitoRiboRef = "RQCFILTER_PATH/chloro_mito_ribo.fa.gz";
	private String mitoRiboRef = "RQCFILTER_PATH/mito_ribo.fa.gz";
	private String chloroRiboRef = "RQCFILTER_PATH/chloro_ribo.fa.gz";

	private String riboRef = "RQCFILTER_PATH/SSU_LSU_sorted.fa.gz";
	
	private String commonMicrobesPath = "RQCFILTER_PATH/commonMicrobes/";
	private String commonMicrobesRef = "RQCFILTER_PATH/commonMicrobes/fusedERPBBmasked2.fa.gz";
	private int commonMicrobesBuild = 1;

	private String otherMicrobesPath = "RQCFILTER_PATH/otherMicrobes/"; //TODO
	private String otherMicrobesRef = "RQCFILTER_PATH/otherMicrobes/fusedEmasked.fa.gz"; //TODO
	private int otherMicrobesBuild = 1;
	
	private String taxTreeFile="RQCFILTER_PATH/tree.taxtree.gz";
	private String giTable=null;
	
	private boolean pairedToSingle=false;
	
	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Library type codes */
	private static final int FRAG=0, LFPE=1, CLIP=2, CLRS=3;
	private static final String clipLinkerDefault = "CATG";
	
}
