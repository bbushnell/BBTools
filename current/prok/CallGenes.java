package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import gff.CompareGff;
import gff.GffLine;
import jgi.BBMerge;
import json.JsonObject;
import map.LongHashSet;
import parse.Parse;
import bin.AdjustEntropy;
import clade.Clade;
import clade.CladeSearcher;
import clade.SendClade;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;

/**
 * This is the executable class for gene-calling.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 *
 */
public class CallGenes extends ProkObject {
	
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
		
		//Create an instance of this class
		CallGenes x=new CallGenes(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CallGenes(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, (args.length>40 ? null : getClass()), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;
			
			outGff=parser.out1;
			maxReads=parser.maxReads;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program
		
		ffoutGff=FileFormat.testOutput(outGff, FileFormat.GFF, null, true, overwrite, append, ordered);
		ffoutAmino=FileFormat.testOutput(outAmino, FileFormat.FA, null, true, overwrite, append, ordered);
		ffout16S=FileFormat.testOutput(out16S, FileFormat.FA, null, true, overwrite, append, ordered);
		ffout18S=FileFormat.testOutput(out18S, FileFormat.FA, null, true, overwrite, append, ordered);
		
		if(ffoutGff!=null){
			assert(!ffoutGff.isSequence()) : "\nout is for gff files.  To output sequence, please use outa.";
		}
		if(ffoutAmino!=null){
			assert(!ffoutAmino.gff()) : "\nouta is for sequence data.  To output gff, please use out.";
		}
		if(ffout16S!=null){
			assert(!ffout16S.gff()) : "\nout16S is for sequence data.  To output gff, please use out.";
		}
		if(ffout18S!=null){
			assert(!ffout18S.gff()) : "\nout18S is for sequence data.  To output gff, please use out.";
		}
		
		if(geneHistFile==null){geneHistBins=0;}
		else{
			assert(geneHistBins>1) : "geneHistBins="+geneHistBins+"; should be >1";
			assert(geneHistDiv>=1) : "geneHistDiv="+geneHistDiv+"; should be >=1";
		}
		geneHist=geneHistBins>1 ? new long[geneHistBins] : null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

//			outstream.println(arg+", "+a+", "+b);
			if(PGMTools.parseStatic(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("infna") || 
					a.equals("fnain") || a.equals("fna") || a.equals("ref")){
				assert(b!=null);
				Tools.addFiles(b, fnaList);
			}else if(b==null && new File(arg).exists() && FileFormat.isFastaFile(arg)){
				fnaList.add(arg);
			}else if(a.equals("pgm") || a.equals("gm") || a.equals("model")){
				assert(b!=null);
				if(b.equalsIgnoreCase("auto") || b.equalsIgnoreCase("default")){
					b=Data.findPath("?model.pgm");
					pgmList.add(b);
				}else{
					Tools.addFiles(b, pgmList);
				}
			}else if(b==null && new File(arg).exists() && FileFormat.isPgmFile(arg)){
				//[prok/CallGenes#001] FIXED: was pgmList.add(b), but b is null in this branch (guarded by b==null);
				//the fna branch above correctly adds arg. A bare-pgm-filename positional arg was adding null to pgmList.
				pgmList.add(arg);
			}else if(a.equals("outamino") || a.equals("aminoout") || a.equals("outa") || a.equals("outaa") || a.equals("aaout") || a.equals("amino")){
				outAmino=b;
			}else if(a.equalsIgnoreCase("out16s") || a.equalsIgnoreCase("16sout")){
				out16S=b;
			}else if(a.equalsIgnoreCase("out18s") || a.equalsIgnoreCase("18sout")){
				out18S=b;
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				//ReadWrite.verbose=verbose;
				GeneCaller.verbose=verbose;
			}else if(a.equalsIgnoreCase("ingff") || a.equalsIgnoreCase("gffin")){
				Tools.addFiles(b, inGffList);
			}
			
			else if(a.equals("json_out") || a.equalsIgnoreCase("json")){
				json_out=Parse.parseBoolean(b);
			}else if(a.equals("stats") || a.equalsIgnoreCase("outstats")){
				outStats=b;
			}else if(a.equals("hist") || a.equalsIgnoreCase("outhist") || a.equalsIgnoreCase("lengthhist") || a.equalsIgnoreCase("lhist") || a.equalsIgnoreCase("genehist")){
				geneHistFile=b;
			}else if(a.equals("bins")){
				geneHistBins=Integer.parseInt(b);
			}else if(a.equals("binlen") || a.equals("binlength") || a.equals("histdiv")){
				geneHistDiv=Integer.parseInt(b);
			}else if(a.equals("printzero") || a.equals("pz")){
				printZeroCountHist=Parse.parseBoolean(b);
			}
			
			else if(a.equals("merge")){
				merge=Parse.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco=Parse.parseBoolean(b);
			}else if(a.equals("2pass") || a.equals("twopass")){
				if(Parse.parseBoolean(b)) {passes=2;}
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
			}

			else if(a.equals("extended") || a.equals("extendedstats") || a.equals("verbosestats")){
				extendedStats=Parse.parseBoolean(b);
			}
			
			else if(a.equalsIgnoreCase("setbias16s")) {
				GeneCaller.biases[r16S]=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("setbias18s")) {
				GeneCaller.biases[r18S]=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("setbias23s")) {
				GeneCaller.biases[r23S]=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("setbias5s")) {
				GeneCaller.biases[r5S]=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("setbiastRNA")) {
				GeneCaller.biases[tRNA]=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("setbiasCDS")) {
				GeneCaller.biases[CDS]=Float.parseFloat(b);
			}
			
			else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}
			
			else if(a.equals("translate")){
				mode=TRANSLATE;
			}else if(a.equals("retranslate") || a.equals("detranslate")){
				mode=RETRANSLATE;
			}else if(a.equals("recode")){
				mode=RECODE;
			}
			
			else if(a.equalsIgnoreCase("minlen") || a.equals("minlength")){
				minLen=Integer.parseInt(b);
			}else if(a.equals("maxoverlapss") || a.equals("overlapss") || a.equals("overlapsamestrand") || a.equals("moss") || a.equalsIgnoreCase("maxOverlapSameStrand")){
				maxOverlapSameStrand=Integer.parseInt(b);
			}else if(a.equals("maxoverlapos") || a.equals("overlapos") || a.equals("overlapoppositestrand") || a.equals("moos") || a.equalsIgnoreCase("maxOverlapOppositeStrand")){
				maxOverlapOppositeStrand=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minStartScore")){
				minStartScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minStopScore")){
				minStopScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minInnerScore") || a.equalsIgnoreCase("minKmerScore")){
				minKmerScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minOrfScore") || a.equalsIgnoreCase("minScore")){
				minOrfScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minAvgScore")){
				minAvgScore=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("dumporfscores")){
				GeneCaller.dumpOrfScores=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("breakLimit")){
				GeneCaller.breakLimit=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("clearcutoffs") || a.equalsIgnoreCase("clearfilters")){
				GeneCaller.breakLimit=9999;
				minOrfScore=-9999;
				minAvgScore=-9999;
				minKmerScore=-9999;
				minStopScore=-9999;
				minStartScore=-9999;
			}

			else if(a.equalsIgnoreCase("e1")){
				Orf.e1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("e2")){
				Orf.e2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("e3")){
				Orf.e3=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("f1")){
				Orf.f1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("f2")){
				Orf.f2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("f3")){
				Orf.f3=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("p0")){
				GeneCaller.p0=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p1")){
				GeneCaller.p1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p2")){
				GeneCaller.p2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p3")){
				GeneCaller.p3=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p4")){
				GeneCaller.p4=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p5")){
				GeneCaller.p5=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("p6")){
				GeneCaller.p6=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("q1")){
				GeneCaller.q1=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q2")){
				GeneCaller.q2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q3")){
				GeneCaller.q3=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q4")){
				GeneCaller.q4=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("q5")){
				GeneCaller.q5=Float.parseFloat(b);
			}
			else if(a.equalsIgnoreCase("lookback")){
				GeneCaller.lookbackPlus=GeneCaller.lookbackMinus=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("lookbackplus")){
				GeneCaller.lookbackPlus=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("lookbackminus")){
				GeneCaller.lookbackMinus=Integer.parseInt(b);
			}
			
			else if(a.equalsIgnoreCase("compareto")){
				compareToGff=b;
			}

			else if(a.equalsIgnoreCase("orfnet") || a.equalsIgnoreCase("scorenet")){
				orfNetPath=b;
			}else if(a.equalsIgnoreCase("orfnetlow") || a.equalsIgnoreCase("netlow")){
				orfNetLowPath=b;
			}else if(a.equalsIgnoreCase("orfnetmid") || a.equalsIgnoreCase("netmid")){
				orfNetMidPath=b;
			}else if(a.equalsIgnoreCase("orfnethigh") || a.equalsIgnoreCase("nethigh")){
				orfNetHighPath=b;
			}else if(a.equalsIgnoreCase("orfnets") || a.equalsIgnoreCase("nets")){
				orfNetPaths=b.split(",");
			}else if(a.equalsIgnoreCase("gcmeans") || a.equalsIgnoreCase("gct") || a.equalsIgnoreCase("gcthresholds")){
				String[] parts=b.split(",");
				gcMeans=new float[parts.length];
				for(int j=0; j<parts.length; j++){gcMeans[j]=Float.parseFloat(parts[j]);}
			}else if(a.equalsIgnoreCase("stopnets")){
				stopNetPaths=b.split(",");
			}else if(a.equalsIgnoreCase("nnstrength")){
				Orf.NN_STRENGTH=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("nnthreshold")){
				Orf.NN_THRESHOLD=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("nnstopstrength")){
				Orf.NN_STOP_STRENGTH=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("nnstopthreshold")){
				Orf.NN_STOP_THRESHOLD=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("metanet")){
				metaNetPath=b;
			}else if(a.equalsIgnoreCase("metanets")){
				metaNetPaths=b.split(",");
			}else if(a.equalsIgnoreCase("nnmetastrength")){
				Orf.NN_META_STRENGTH=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("nnmetathreshold")){
				Orf.NN_META_THRESHOLD=Float.parseFloat(b);
			}

			else if(a.equalsIgnoreCase("taxonomy") || a.equalsIgnoreCase("tax")){
				useTaxonomy=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("local") || a.equalsIgnoreCase("localclade")){
				useLocalClade=(b==null || Parse.parseBoolean(b));
			}else if(a.equalsIgnoreCase("server") || a.equalsIgnoreCase("serverclade")){
				useLocalClade=!(b==null || Parse.parseBoolean(b));
			}else if(a.equalsIgnoreCase("percontig")){
				perContig=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("taxaddress") || a.equalsIgnoreCase("taxserver")){
				taxAddress=b;
			}else if(a.equalsIgnoreCase("trnaalign") || a.equalsIgnoreCase("aligntrna")){
				trnaAlign=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("trnalib")){
				GeneCaller.trnaLibrary=TrnaConsensusBuilder.loadLibrary(b);
				GeneCaller.trnaModelNames=TrnaConsensusBuilder.lastLoadedNames;
				trnaLibExplicit=true;
			}else if(a.equalsIgnoreCase("trnamodel")){
				GeneCaller.trnaModels=TrnaConsensusBuilder.loadModels(b);
				trnaLibExplicit=true;
			}else if(a.equalsIgnoreCase("idpass")){
				prok.TrnaCaller.ID_PASS=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("idborderline")){
				prok.TrnaCaller.ID_BORDERLINE=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("idborderlinelong")){
				prok.TrnaCaller.ID_BORDERLINE_LONG=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("hbmpass")){
				prok.TrnaCaller.HBM_PASS=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("indextopn")){
				prok.TrnaCaller.INDEX_TOP_N_OVERRIDE=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("indexminhits")){
				prok.TrnaCaller.INDEX_MINHITS_OVERRIDE=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("adaptiveminhits") || a.equalsIgnoreCase("adaptiveshortlist")){
				prok.TrnaCaller.ADAPTIVE_MINHITS=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("adaptfloor")){
				prok.TrnaCaller.ADAPT_FLOOR=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("adapttopfrac")){
				prok.TrnaCaller.ADAPT_TOPFRAC=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("adaptqfrac")){
				prok.TrnaCaller.ADAPT_QFRAC=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("shortliststats") || a.equalsIgnoreCase("slstats")){
				prok.TrnaCaller.SHORTLIST_STATS=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("refhist")){
				prok.TrnaCaller.REFHIST=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("indexk") || a.equalsIgnoreCase("shortlistk")){
				prok.TrnaCaller.INDEX_K=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("earlyexit")){
				prok.TrnaCaller.earlyExit=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("patience")){
				prok.TrnaCaller.earlyExitPatience=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("trnadebug")){
				prok.TrnaCaller.DEBUG=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("trnaintron")){
				prok.TrnaCaller.INTRON_PASS=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("trnaintronac") || a.equalsIgnoreCase("intronanticodon")){
				prok.TrnaCaller.INTRON_ANTICODON_RECOVERY=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("trnascavenge") || a.equalsIgnoreCase("scavenge")){
				prok.TrnaCaller.SCAVENGE=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("trnascavengeonly") || a.equalsIgnoreCase("scavengeonly")){
				prok.TrnaCaller.SCAVENGE_ONLY=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("mintrnakhits") || a.equalsIgnoreCase("mintrnakmerhits")){
				prok.TrnaCaller.MIN_TRNA_KHITS=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("scavpad") || a.equalsIgnoreCase("scavengepad")){
				prok.TrnaCaller.SCAV_PAD=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("scavcollapse") || a.equalsIgnoreCase("collapsefrac")){
				prok.TrnaCaller.SCAV_COLLAPSE_FRAC=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("scavmultilocus") || a.equalsIgnoreCase("trnamultilocus")){
				prok.TrnaCaller.SCAV_MULTILOCUS=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("scavlocusoverlap") || a.equalsIgnoreCase("trnalocusoverlap")){
				prok.TrnaCaller.SCAV_LOCUS_OVERLAP_FRAC=Float.parseFloat(b);
				if(!(prok.TrnaCaller.SCAV_LOCUS_OVERLAP_FRAC>=0 && prok.TrnaCaller.SCAV_LOCUS_OVERLAP_FRAC<=1)){
					throw new IllegalArgumentException("scavlocusoverlap must be in [0,1]: "+b);
				}
			}else if(a.equalsIgnoreCase("acextract") || a.equalsIgnoreCase("anticodonextract")){
				prok.TrnaCaller.extractAnticodons=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("acvalidate") || a.equalsIgnoreCase("acthresh")){
				prok.TrnaCaller.AC_VALIDATE=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("acmargin")){
				prok.TrnaCaller.AC_MARGIN=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("trimtrna") || a.equalsIgnoreCase("trnatrim")){
				prok.TrnaCaller.trimToAlignment=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("maxtrna") || a.equalsIgnoreCase("maxtrnalen")){
				prok.TrnaCaller.MAX_TRNA_OVERRIDE=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("trnacutoff1") || a.equalsIgnoreCase("trnaregion")){
				GeneCaller.cutoff1[ProkObject.tRNA]=Float.parseFloat(b);//region-open threshold
			}else if(a.equalsIgnoreCase("trnacutoff2") || a.equalsIgnoreCase("trnacand")){
				GeneCaller.cutoff2[ProkObject.tRNA]=Float.parseFloat(b);//composite candidate threshold
			}else if(a.equalsIgnoreCase("trnacutoff3") || a.equalsIgnoreCase("trnastart")){
				GeneCaller.cutoff3[ProkObject.tRNA]=Float.parseFloat(b);//start point-score threshold
			}else if(a.equalsIgnoreCase("trnacutoff4") || a.equalsIgnoreCase("trnastop")){
				GeneCaller.cutoff4[ProkObject.tRNA]=Float.parseFloat(b);//stop point-score threshold
			}else if(a.equalsIgnoreCase("trnacutoff5") || a.equalsIgnoreCase("trnainner")){
				GeneCaller.cutoff5[ProkObject.tRNA]=Float.parseFloat(b);//avg inner-score threshold
			}else if(a.equalsIgnoreCase("trnaboundarynet")){
				//DEFAULT-ON (Brian, Aug 22 2026, via Noire): the boundary refiner ships ACTIVE
				//by default (see the post-parse resolution block below) -- trnaboundarynet=f/
				//false is now the OPT-OUT switch, not "not given." Any other value is a path
				//override (matches the existing pgm=auto/default convention's spirit, inverted).
				if(b!=null && (b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false"))){
					trnaBoundaryDisabled=true;
				}else{
					trnaBoundaryNetPath=b;
				}
			}else if(a.equalsIgnoreCase("trnaboundary5net")){
				trnaBoundary5NetPath=b;
			}else if(a.equalsIgnoreCase("trnaboundary3net")){
				trnaBoundary3NetPath=b;
			}else if(a.equalsIgnoreCase("trnaboundarystarttable")){
				trnaBoundaryStartTablePath=b;
			}else if(a.equalsIgnoreCase("trnaboundarystoptable")){
				trnaBoundaryStopTablePath=b;
			}else if(a.equalsIgnoreCase("trnaboundarymargin")){
				prok.TrnaBoundaryScorer.MARGIN_THRESHOLD_START=Float.parseFloat(b);
				prok.TrnaBoundaryScorer.MARGIN_THRESHOLD_STOP=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("trnaboundarymarginstart")){
				prok.TrnaBoundaryScorer.MARGIN_THRESHOLD_START=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("trnaboundarymarginstop")){
				prok.TrnaBoundaryScorer.MARGIN_THRESHOLD_STOP=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("trnaboundarycrossenrichment")){
				prok.TrnaBoundaryVectorGen.INCLUDE_CROSS_ENRICHMENT=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("ncrna") || a.equalsIgnoreCase("generalncrna")){
				NCRNA_FAMILIES_ENABLED=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("ncrnaboundarynet")){
				NCRNA_BOUNDARY_NN_ENABLED=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("tmrna")){
				TMRNA_ENABLED=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("sixs") || a.equalsIgnoreCase("ssrs") || a.equalsIgnoreCase("6s")){
				SIXS_ENABLED=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("ncrnafamily")){
				NCRNA_FAMILY_FILTER=parseNcrnaFamily(b);
			}else if(a.equalsIgnoreCase("ncrnakmers")){
				NCRNA_KMERS_OVERRIDE=b;
			}else if(a.equalsIgnoreCase("rnasepkmers")){
				RNASEP_KMERS_OVERRIDE=b;
			}else if(a.equalsIgnoreCase("srpsmallkmers")){
				SRPSMALL_KMERS_OVERRIDE=b;
			}else if(a.equalsIgnoreCase("srplargekmers")){
				SRPLARGE_KMERS_OVERRIDE=b;
			}else if(a.equalsIgnoreCase("tmrnakmers")){
				TMRNA_KMERS_OVERRIDE=b;
			}else if(a.equalsIgnoreCase("sixsrf00013kmers")){
				SIXS_RF00013_KMERS_OVERRIDE=b;
			}else if(a.equalsIgnoreCase("sixsrf01685kmers")){
				SIXS_RF01685_KMERS_OVERRIDE=b;
			}else if(a.equalsIgnoreCase("tmrnaconsensus")){
				TMRNA_CONSENSUS_OVERRIDE=b;
			}else if(a.equalsIgnoreCase("tmrnamodels")){
				TMRNA_MODELS_OVERRIDE=b;
			}else if(a.equalsIgnoreCase("ncrnaidpass")){
				NCRNA_ID_PASS_OVERRIDE=parseUnitFloat(a, b);
			}else if(a.equalsIgnoreCase("ncrnaidborderline")){
				NCRNA_ID_BORDERLINE_OVERRIDE=parseUnitFloat(a, b);
			}else if(a.equalsIgnoreCase("ncrnahbmpass")){
				NCRNA_HBM_PASS_OVERRIDE=parseUnitFloat(a, b);
			}else if(a.equalsIgnoreCase("ncrnascorea")){
				NCRNA_SCORE_A_OVERRIDE=parseFiniteFloat(a, b);
			}else if(a.equalsIgnoreCase("ncrnascoreb")){
				NCRNA_SCORE_B_OVERRIDE=parseFiniteFloat(a, b);
			}else if(a.equalsIgnoreCase("rnasepscorea")){
				RNASEP_SCORE_A_OVERRIDE=parseFiniteFloat(a, b);
			}else if(a.equalsIgnoreCase("rnasepscoreb")){
				RNASEP_SCORE_B_OVERRIDE=parseFiniteFloat(a, b);
			}else if(a.equalsIgnoreCase("srpsmallscorea")){
				SRPSMALL_SCORE_A_OVERRIDE=parseFiniteFloat(a, b);
			}else if(a.equalsIgnoreCase("srpsmallscoreb")){
				SRPSMALL_SCORE_B_OVERRIDE=parseFiniteFloat(a, b);
			}else if(a.equalsIgnoreCase("srplargescorea")){
				SRPLARGE_SCORE_A_OVERRIDE=parseFiniteFloat(a, b);
			}else if(a.equalsIgnoreCase("srplargescoreb")){
				SRPLARGE_SCORE_B_OVERRIDE=parseFiniteFloat(a, b);
			}else if(a.equalsIgnoreCase("tmrnascorea")){
				TMRNA_SCORE_A_OVERRIDE=parseFiniteFloat(a, b);
			}else if(a.equalsIgnoreCase("tmrnascoreb")){
				TMRNA_SCORE_B_OVERRIDE=parseFiniteFloat(a, b);
			}else if(a.equalsIgnoreCase("ncrnacollapsefrac")){
				NCRNA_COLLAPSE_FRAC_OVERRIDE=parseUnitFloat(a, b);
			}else if(a.equalsIgnoreCase("ncrnawindowpad") || a.equalsIgnoreCase("ncrnapad")){
				NCRNA_WINDOW_PAD_OVERRIDE=parsePadOverride(a, b);
			}else if(a.equalsIgnoreCase("rnaseppad")){
				RNASEP_PAD_OVERRIDE=parsePadOverride("rnaseppad", b);
			}else if(a.equalsIgnoreCase("srpsmallpad")){
				SRPSMALL_PAD_OVERRIDE=parsePadOverride("srpsmallpad", b);
			}else if(a.equalsIgnoreCase("srplargepad")){
				SRPLARGE_PAD_OVERRIDE=parsePadOverride("srplargepad", b);
			}else if(a.equalsIgnoreCase("tmrnapad")){
				TMRNA_PAD_OVERRIDE=parsePadOverride("tmrnapad", b);
			}else if(a.equalsIgnoreCase("sixsrf00013pad")){
				SIXS_RF00013_PAD_OVERRIDE=parsePadOverride("sixsrf00013pad", b);
			}else if(a.equalsIgnoreCase("sixsrf01685pad")){
				SIXS_RF01685_PAD_OVERRIDE=parsePadOverride("sixsrf01685pad", b);
			}

			else if(ProkObject.parse(arg, a, b)){}

			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		if(pgmList.isEmpty()){
			String b=Data.findPath("?model.pgm");
			pgmList.add(b);
		}
		for(int i=0; i<pgmList.size(); i++){
			String s=pgmList.get(i);
			if(s.equalsIgnoreCase("auto") || s.equalsIgnoreCase("default")){
				pgmList.set(i, Data.findPath("?model.pgm"));
			}
		}
		
		if(trnaAlign && !trnaLibExplicit){
			loadTrnaResources();
		}else if(!trnaAlign){
			GeneCaller.trnaLibrary=null;
			GeneCaller.trnaModels=null;
			GeneCaller.trnaModelNames=null;
		}
		//Boundary-precision NN -- SHIPPED, DEFAULT ON (Brian, Aug 22 2026, via Noire): a
		//bare callgenes invocation (no tRNA boundary flags at all) now includes the refiner
		//automatically, matching the trna skill's "bare invocation IS the tuned config"
		//convention that already governs every other tRNA default. trnaboundarynet=f (or
		//=false) is the opt-out switch (parsed above). Explicit trnaboundarynet=/5net=/3net=
		//paths still override the shipped default when given, exactly as before.
		//TWO-NET RESOLUTION (Noire's spec, Aug 22 2026): trnaboundary5net=/trnaboundary3net=
		//each fall back to trnaboundarynet= when not individually given -- so trnaboundarynet=
		//alone still fully specifies both, while either 5net=/3net= alone requires
		//trnaboundarynet= as the fallback for the other side (asserted below, not silently
		//left half-specified).
		if(!trnaBoundaryDisabled){
			if(trnaBoundaryNetPath==null && trnaBoundary5NetPath==null && trnaBoundary3NetPath==null){
				trnaBoundaryNetPath=Data.findPath("?trna_boundary_net.bbnet");
			}
			if(trnaBoundaryStartTablePath==null){trnaBoundaryStartTablePath=Data.findPath("?trna_boundary_start_table.tsv");}
			if(trnaBoundaryStopTablePath==null){trnaBoundaryStopTablePath=Data.findPath("?trna_boundary_stop_table.tsv");}
		}
		final boolean boundaryNetGiven=(trnaBoundaryNetPath!=null || trnaBoundary5NetPath!=null || trnaBoundary3NetPath!=null);
		final String trnaBoundaryNet5Resolved=(trnaBoundary5NetPath!=null ? trnaBoundary5NetPath : trnaBoundaryNetPath);
		final String trnaBoundaryNet3Resolved=(trnaBoundary3NetPath!=null ? trnaBoundary3NetPath : trnaBoundaryNetPath);
		assert(boundaryNetGiven==(trnaBoundaryStartTablePath!=null) && boundaryNetGiven==(trnaBoundaryStopTablePath!=null))
			: "trnaboundarynet=/trnaboundary5net=/trnaboundary3net= must be given (at least one) together with "
			+"trnaboundarystarttable=/trnaboundarystoptable=, or none of them at all.";
		assert(!boundaryNetGiven || (trnaBoundaryNet5Resolved!=null && trnaBoundaryNet3Resolved!=null))
			: "trnaboundary5net= or trnaboundary3net= was given alone without trnaboundarynet= to supply the "
			+"other side -- give trnaboundarynet= as the shared fallback, or specify both 5net=/3net= explicitly.";
		if(boundaryNetGiven){
			prok.TrnaCaller.loadBoundaryNet(trnaBoundaryNet5Resolved, trnaBoundaryNet3Resolved,
				trnaBoundaryStartTablePath, trnaBoundaryStopTablePath);
		}
		//Gate A (forward-ported from Noire's tree, 2026-08-28, C3 merge): generic ncRNA family
		//loading is OFF by default -- Noire's original call here was unconditional
		//(`loadNcrnaResources();`); gating it behind NCRNA_FAMILIES_ENABLED is the entire
		//difference between "ncRNAs off by default" (Citan's explicit requirement) and the
		//shipped-default-on treatment tRNA's boundary NN got just above. When false, this line
		//never executes and GeneCaller.ncrnaFamilies stays empty -- zero cost, zero behavior
		//change from pre-merge master.
		//Citan's correction (2026-08-28): ncrnaboundarynet=t with ncrna=f is a nonsensical
		//combination (nothing to attach boundary refinement to) -- must fail loud as invalid
		//input, not silently no-op. Factored into a small testable method (see
		//CallGenesNcrnaGateTest); matches the existing tRNA flag-validation style (assert,
		//always live under BBTools' always-on -ea) rather than a special-cased exception.
		validateNcrnaGateCombo();
		validateNcrnaSweepOverrides();
		if(NCRNA_FAMILIES_ENABLED){loadNcrnaResources();}

		if(Shared.threads()<2){ordered=false;}
		assert(!fnaList.isEmpty()) : "At least 1 fasta file is required.";
		assert(!pgmList.isEmpty()) : "At least 1 pgm file is required.";
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		fnaList=Tools.fixExtension(fnaList);
		pgmList=Tools.fixExtension(pgmList);
		if(fnaList.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, outGff, outAmino, out16S, out18S, outStats, geneHistFile)){
			outstream.println((outGff==null)+", "+outGff);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "
					+outGff+", "+outAmino+", "+out16S+", "+out18S+", "+outStats+", "+geneHistFile+"\n");
		}
		
		//Ensure input files can be read
		ArrayList<String> foo=new ArrayList<String>();
		foo.addAll(fnaList);
		foo.addAll(pgmList);
		if(!Tools.testInputFiles(false, true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		foo.add(outGff);
		foo.add(outAmino);
		foo.add(out16S);
		foo.add(out18S);
		foo.add(outStats);
		foo.add(geneHistFile);
		if(!Tools.testForDuplicateFiles(true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create read streams and process all data */
	void process(Timer t){

		//Default/fallback model (respects the user pgm=/model= list).  With taxonomy=t, per-phylum PGMs are
		//now selected PER FILE inside the processing loop (classify -> cached per-phylum PGM), so ONE JVM
		//handles a whole multi-genome batch with correct per-genome taxonomy: the tRNA library, nets, and
		//taxonomy client load ONCE, and each phylum's PGM loads once (cached) -- no per-genome JVM restart.
		final GeneModel defaultPgm=PGMTools.loadAndMerge(pgmList);
		final boolean perFileTax=(useTaxonomy && !perContig);
		final java.util.HashMap<String,GeneModel> phylumCache=new java.util.HashMap<>();

		if(orfNetPath!=null){GeneCaller.loadOrfNet(orfNetPath);}
		if(orfNetLowPath!=null || orfNetMidPath!=null || orfNetHighPath!=null){
			GeneCaller.loadOrfNetsGC(orfNetLowPath, orfNetMidPath, orfNetHighPath);
		}
		if(gcMeans!=null){GeneCaller.gcMeans=gcMeans;}
		if(orfNetPaths!=null){
			GeneCaller.loadOrfNetsGCN(orfNetPaths, gcMeans);
		}
		if(stopNetPaths!=null){
			GeneCaller.loadStopNetsGCN(stopNetPaths);
		}
		if(metaNetPaths!=null){
			GeneCaller.loadMetaNetsGCN(metaNetPaths);
		}else if(metaNetPath!=null){
			GeneCaller.loadMetaNetSingle(metaNetPath);
		}

		if(call16S || call18S || call23S || calltRNA || call5S){
			loadLongKmers();
			loadConsensusSequenceFromFile(false, false);
		}
		
		ByteStreamWriter bsw=makeBSW(ffoutGff);
		if(bsw!=null){
			bsw.forcePrint("##gff-version 3\n");
		}

		ConcurrentReadOutputStream rosAmino=makeCros(ffoutAmino);
		ConcurrentReadOutputStream ros16S=makeCros(ffout16S);
		ConcurrentReadOutputStream ros18S=makeCros(ffout18S);
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Reset counters
		readsIn=genesOut=0;
		basesIn=bytesOut=0;
		
		for(int fnum=0; fnum<fnaList.size(); fnum++){
			final String fna=fnaList.get(fnum);
			String gffIn=(inGffList!=null && !inGffList.isEmpty()) ? inGffList.set(fnum, null) : null;
			//Per-file taxonomy: classify THIS genome, use its cached per-phylum PGM (else the default).
			final GeneModel pgm0=(perFileTax ? pgmForPhylumCached(classifyPhylum(fna), defaultPgm, phylumCache) : defaultPgm);
			if(perFileTax && bsw!=null){
				//GFF header provenance: what QuickClade detected for this file and how, so an
				//in-process caller reading the GFF alone (not the Java API) still gets it.
				bsw.forcePrint("##Domain "+(lastDetectedDomain!=null ? lastDetectedDomain : "unknown")+"\n");
				bsw.forcePrint("##Phylum "+(lastDetectedPhylum!=null ? lastDetectedPhylum : "unknown")+"\n");
				bsw.forcePrint("##ClassificationSource "+(lastClassificationSource!=null ? lastClassificationSource : "none")+"\n");
			}
			//Create a read input stream
			final GeneModel pgm=makeMultipassModel(pgm0, fna, gffIn, passes/*, maxReads*/);
			
			final ConcurrentReadInputStream cris=makeCris(fna);
			
			//Process the reads in separate threads
			spawnThreads(cris, bsw, rosAmino, ros16S, ros18S, pgm);
			
			//Close the input stream
			errorState|=ReadWrite.closeStream(cris);
		}
		
		//Close the input stream
		errorState|=ReadWrite.closeStreams(null, rosAmino, ros16S, ros18S);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the output stream
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		outstream.println(Tools.linesBytesOut(readsIn, basesIn, genesOut, bytesOut, 8, false));
		outstream.println();
		
		if(json_out){
			printStatsJson(outStats);
		}else{
			printStats(outStats);
		}
		
		if(geneHistFile!=null){
			printHist(geneHistFile);
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
		
		if(compareToGff!=null){
			if(compareToGff.equals("auto")){
				compareToGff=fnaList.get(0).replace(".fna", ".gff");
				compareToGff=compareToGff.replace(".fa", ".gff");
				compareToGff=compareToGff.replace(".fasta", ".gff");
			}
			CompareGff.main(new String[] {outGff, compareToGff});
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Taxonomy Classification     ----------------*/
	/*--------------------------------------------------------------*/

	/** Selects (and caches) the per-phylum PGM for one classified genome; falls back to defaultPgm when the
	 * phylum is unknown or has no shipped pgm_<phylum>.pgm.  Prints only on a cache miss (once per phylum). */
	private GeneModel pgmForPhylumCached(String phylum, GeneModel defaultPgm, java.util.HashMap<String,GeneModel> cache){
		if(phylum==null){return defaultPgm;}
		final GeneModel cached=cache.get(phylum);
		if(cached!=null){return cached;}
		final String path=Data.findPath("?pgm_"+phylum+".pgm", false);
		final GeneModel gm;
		if(path!=null && new File(path).exists()){
			final ArrayList<String> list=new ArrayList<>();
			list.add(path);
			gm=PGMTools.loadAndMerge(list);
			outstream.println("Taxonomy: using "+phylum+" PGM");
		}else{
			gm=defaultPgm;
			outstream.println("Taxonomy: "+phylum+" (no PGM available, using default)");
		}
		cache.put(phylum, gm);
		return gm;
	}

	/** Classifies one genome's domain+phylum via QuickClade (local database if local=t and
	 * available, else the network server), storing the full result (domain, phylum, and which
	 * path served it) in lastDetectedDomain/lastDetectedPhylum/lastClassificationSource for an
	 * in-process caller, and returning just the phylum for this method's original callers
	 * (per-phylum PGM selection). Never throws; a classification failure just returns null and
	 * leaves the per-file PGM selection to fall back to the default. */
	private String classifyPhylum(String fna){
		lastDetectedDomain=null; lastDetectedPhylum=null; lastClassificationSource=null;
		try{
			Clade.MAKE_DDLS=true;
			Clade.DDL_K=25;
			Clade.DDL_BUCKETS=32768;
			cardinality.DynamicDemiLog.setExponent(5);
			if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150){
				AdjustEntropy.load(4, 150);
			}
			ArrayList<Read> reads=ReadInputStream.toReads(fna, FileFormat.FASTA, -1);
			Clade clade=new Clade(0, 0, fna);
			for(Read r : reads){clade.add(r.bases, null);}
			clade.finish();

			String lineage=null;
			if(useLocalClade){
				lineage=CladeSearcher.classifyLocal(clade);
				if(lineage!=null){lastClassificationSource="local";}
				else{outstream.println("Warning: local=t requested but no local QuickClade "
					+"reference is available; falling back to the server.");}
			}
			if(lineage==null){
				//taxAddress's "refseq" default is a sentinel meaning "use SendClade's own default
				//address", not a literal server value -- pass null through in that case, exactly
				//as this call always did before local= existed.
				final String address=(taxAddress==null || taxAddress.equals("refseq")) ? null : taxAddress;
				final ArrayList<Clade> clades=new ArrayList<>();
				clades.add(clade);
				final String response=SendClade.sendClades(clades, address, true, 1, false, false, 1, false, 1);
				lastClassificationSource="server"+(address!=null ? ":"+address : "");
				//The response is tab-delimited with TWO different semicolon-joined columns: the
				//real lineage (d__Bacteria;k__...;p__...;...) and a separate per-rank confidence
				//column (d:100.0;k:99.7;p:99.7;...). Splitting the WHOLE line on ';' without first
				//isolating the right column corrupts extraction at tab boundaries -- real bug found
				//2026-08-29 (verified against a live server response): domain came back "unknown"
				//on every real classification because the d__ segment ran into the preceding
				//tab-separated field and never actually started with "d__" after trim(), while
				//phylum's p__ segment happened to land cleanly by luck. Isolate the lineage column
				//by its double-underscore rank markers (unique to it -- the confidence column uses
				//single colons, no "__") before splitting on ';'.
				if(response!=null){
					outer:
					for(String line : response.split("\n")){
						if(line.startsWith("#")){continue;}
						for(String field : line.split("\t")){
							if(field.contains("__")){lineage=field; break outer;}
						}
					}
				}
			}

			if(lineage!=null){
				lastDetectedDomain=findDomainInLineage(lineage);
				lastDetectedPhylum=findPhylumInLineage(lineage);
				if(lastDetectedPhylum!=null){return lastDetectedPhylum;}
			}
		}catch(Exception e){
			outstream.println("Warning: taxonomy classification failed: "+e.getMessage());
		}
		return null;
	}

	private static String findPhylumInLineage(String lineage){
		for(String part : lineage.split(";")){
			String trimmed=part.trim();
			if(trimmed.startsWith("p:")){return trimmed.substring(2);}
			else if(trimmed.startsWith("p__")){return trimmed.substring(3);}
		}
		return null;
	}

	/** Domain (superkingdom) counterpart to findPhylumInLineage -- current NCBI lineages carry
	 * a "d:"/"d__" rank (superkingdom was renamed domain; see clade/CladeIndex.java's own
	 * LCA_PREFIXES table, which documents the same convention for its lineage-LCA matching). */
	private static String findDomainInLineage(String lineage){
		for(String part : lineage.split(";")){
			String trimmed=part.trim();
			if(trimmed.startsWith("d:")){return trimmed.substring(2);}
			else if(trimmed.startsWith("d__")){return trimmed.substring(3);}
		}
		return null;
	}

	/** @return Domain detected by the most recent classifyPhylum() call (in-process caller access) */
	public String lastDetectedDomain(){return lastDetectedDomain;}
	/** @return Phylum detected by the most recent classifyPhylum() call (in-process caller access) */
	public String lastDetectedPhylum(){return lastDetectedPhylum;}
	/** @return "local" or "server[:address]", or null if classification was never attempted */
	public String lastClassificationSource(){return lastClassificationSource;}


	/**
	 * Writes gene calling statistics to specified file in human-readable format.
	 * Includes gene counts by type, coding fraction estimates, and detailed
	 * statistics for each gene category when extended stats are enabled.
	 * @param fname Output filename, or null to skip statistics output
	 */
	private void printStats(String fname){
		if(fname==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, false);
		bsw.start();
		
		if(callCDS){
			if(extendedStats || verbose) {
				bsw.println("Gene Starts Made:     \t "+Tools.padLeft(geneStartsMade, 12));
				bsw.println("Gene Stops Made:      \t "+Tools.padLeft(geneStopsMade, 12));
				bsw.println("Gene Starts Retained: \t "+Tools.padLeft(geneStartsRetained, 12));
				bsw.println("Gene Stops Retained:  \t "+Tools.padLeft(geneStopsRetained, 12));
			}
			bsw.print("Approx Coding Fraction:\t ");
			bsw.print(Tools.padLeft(Tools.format("%.4f", Tools.min(1.0, stCdsPass.lengthSum/(double)basesIn)), 12)).nl();
			bsw.println("CDS Out:              \t "+Tools.padLeft(geneStartsOut, 12));
		}
		if(call16S){bsw.println("16S Out:              \t "+Tools.padLeft(r16SOut, 12));}
		if(call18S){bsw.println("18S Out:              \t "+Tools.padLeft(r18SOut, 12));}
		if(call23S){bsw.println("23S Out:              \t "+Tools.padLeft(r23SOut, 12));}
		if(call5S){bsw.println("5S Out:               \t "+Tools.padLeft(r5SOut, 12));}
		if(calltRNA){bsw.println("tRNA Alignments:  \t "+Tools.padLeft(prok.TrnaCaller.alignmentCount(), 12));}
			if(calltRNA){bsw.println("tRNA Out:             \t "+Tools.padLeft(tRNAOut, 12));}
		if(!GeneCaller.ncrnaFamilies.isEmpty()){
			for(int i=0; i<GeneCaller.ncrnaFamilies.size(); i++){
				NcrnaFamily family=GeneCaller.ncrnaFamilies.get(i);
				bsw.println("ncRNA Alignments ("+family.name+"):\t "+Tools.padLeft(ncrnaAlignments[i], 12));
			}
			bsw.println("RNA Out:              \t "+Tools.padLeft(rnaOut, 12));
		}
		
		if(extendedStats || verbose) {
			if(callCDS){
				bsw.println();
				bsw.println("All ORF Stats:");
				bsw.print(stCds.toString());

				//			bsw.println();
				//			bsw.println("Retained ORF Stats:");
				//			bsw.print(stCds2.toString());

				bsw.println();
				bsw.println("Called ORF Stats:");
				stCdsPass.genomeSize=basesIn;
				bsw.print(stCdsPass.toString());
			}

			if(call16S && r16SOut>0){
				bsw.println();
				bsw.println("Called 16S Stats:");
				bsw.print(st16s.toString());
			}
			if(call23S && r23SOut>0){
				bsw.println();
				bsw.println("Called 23S Stats:");
				bsw.print(st23s.toString());
			}
			if(call5S && r5SOut>0){
				bsw.println();
				bsw.println("Called 5S Stats:");
				bsw.print(st5s.toString());
			}
			if(call18S && r18SOut>0){
				bsw.println();
				bsw.println("Called 18S Stats:");
				bsw.print(st18s.toString());
			}
			if(calltRNA && tRNAOut>0){
				bsw.println();
				bsw.println("Called tRNA Stats:");
				bsw.print(sttRNA.toString());
			}
		}
		bsw.poisonAndWait();
	}
	
	/**
	 * Writes gene calling statistics to specified file in JSON format.
	 * Provides machine-readable output suitable for downstream analysis pipelines.
	 * @param fname Output filename, or null to skip JSON output
	 */
	private void printStatsJson(String fname){
		if(fname==null){return;}
		
		JsonObject outer=new JsonObject();
		
		{
			JsonObject jo=new JsonObject();
			if(callCDS){
				jo.add("Gene Starts Made", geneStartsMade);
				jo.add("Gene Stops Made", geneStopsMade);
				jo.add("Gene Starts Retained", geneStartsRetained);
				jo.add("Gene Stops Retained", geneStopsRetained);
				jo.add("CDS Out", geneStartsOut);
			}
			if(call16S){jo.add("16S Out", r16SOut);}
			if(call18S){jo.add("18S Out", r18SOut);}
			if(call23S){jo.add("23S Out", r23SOut);}
			if(call5S){jo.add("5S Out", r5SOut);}
			if(calltRNA){jo.add("tRNA Out", tRNAOut);}
			if(!GeneCaller.ncrnaFamilies.isEmpty()){
				for(int i=0; i<GeneCaller.ncrnaFamilies.size(); i++){
					NcrnaFamily family=GeneCaller.ncrnaFamilies.get(i);
					jo.add("ncRNA Alignments ("+family.name+")", ncrnaAlignments[i]);
				}
				jo.add("RNA Out", rnaOut);
			}
			outer.add("Overall", jo);
		}
		
		if(callCDS){
			outer.add("All ORF Stats", stCds.toJson());
			outer.add("Retained ORF Stats", stCds2.toJson());
			stCdsPass.genomeSize=basesIn;
			outer.add("Called ORF Stats", stCdsPass.toJson());
		}

		//[prok/ScoreTracker#001] these toJson() calls are UNGUARDED (no >0 check) -- unlike printStats() above which guards each with "&& rXOut>0". An enabled-but-empty type -> count==0 -> 0.0/0=NaN -> a bare NaN token = invalid JSON. The cheap fix mirrors printStats' rXOut>0 guard, or fix non-finite centrally in ScoreTracker.toJson.
		if(call16S){outer.add("Called 16S Stats", st16s.toJson());}
		if(call18S){outer.add("Called 18S Stats", st18s.toJson());}
		if(call23S){outer.add("Called 23S Stats", st23s.toJson());}
		if(call5S){outer.add("Called 5S Stats", st5s.toJson());}
		if(calltRNA){outer.add("Called tRNA Stats", sttRNA.toJson());}
		

		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, false);
		bsw.start();
		bsw.println(outer.toText());
		bsw.poisonAndWait();
	}
	
	/**
	 * Writes gene length histogram to specified file with summary statistics.
	 * Includes mean, median, and standard deviation of gene lengths.
	 * @param fname Output filename for histogram data
	 */
	private void printHist(String fname){
		if(fname==null || geneHist==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, false);
		bsw.start();
		long sum=simd.Vector.sum(geneHist);
		double mean=Tools.averageHistogram(geneHist)*geneHistDiv;
		int median=Tools.medianHistogram(geneHist)*geneHistDiv;
		double std=Tools.standardDeviationHistogram(geneHist)*geneHistDiv;
		bsw.println("#Gene Length Histogram");
		bsw.print("#Genes:\t").println(sum);
		bsw.print("#Mean:\t").println(mean, 4);
		bsw.print("#Median:\t").println(median);
		bsw.print("#STDDev:\t").println(std, 4);
		bsw.print("#Length\tCount\n");
		long cum=0;
		for(int i=0; i<geneHist.length && cum<sum; i++){
			int len=i*geneHistDiv;
			long count=geneHist[i];
			cum+=count;
			if(count>0 || printZeroCountHist){
				bsw.print(len).tab().println(count);
			}
		}
		bsw.poisonAndWait();
	}
	
	/**
	 * Creates concurrent read input stream for specified FASTA file.
	 * Optimizes buffering and threading for high-throughput sequence processing.
	 * @param fname Input FASTA filename
	 * @return Configured ConcurrentReadInputStream
	 */
	private ConcurrentReadInputStream makeCris(String fname){
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ffin, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ByteStreamWriter bsw, 
			ConcurrentReadOutputStream rosAmino, ConcurrentReadOutputStream ros16S, ConcurrentReadOutputStream ros18S, GeneModel pgm){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, bsw, rosAmino, ros16S, ros18S, pgm, minLen, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for threads to finish
		waitForThreads(alpt);
		
		//Do anything necessary after processing
		
	}
	
	/**
	 * Waits for all processing threads to complete and aggregates results.
	 * Accumulates per-thread statistics including read counts, gene counts,
	 * and score tracking data across all worker threads.
	 * @param alpt List of ProcessThread instances to await
	 */
	private void waitForThreads(ArrayList<ProcessThread> alpt){
		
		//Wait for completion of all threads
		boolean success=true;
		for(ProcessThread pt : alpt){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			genesOut+=pt.genesOutT;
			bytesOut+=pt.bytesOutT;
			
			geneStopsMade+=pt.caller.geneStopsMade;
			geneStartsMade+=pt.caller.geneStartsMade;
			geneStartsRetained+=pt.caller.geneStartsRetained;
			geneStopsRetained+=pt.caller.geneStopsRetained;
			geneStartsOut+=pt.caller.geneStartsOut;

			r16SOut+=pt.caller.r16SOut;
			r18SOut+=pt.caller.r18SOut;
			r23SOut+=pt.caller.r23SOut;
			r5SOut+=pt.caller.r5SOut;
			tRNAOut+=pt.caller.tRNAOut;
			rnaOut+=pt.caller.rnaOut;
			long[] threadNcrnaAlignments=pt.caller.ncrnaAlignmentCounts();
			if(ncrnaAlignments==null){ncrnaAlignments=new long[threadNcrnaAlignments.length];}
			assert(ncrnaAlignments.length==threadNcrnaAlignments.length);
			for(int i=0; i<threadNcrnaAlignments.length; i++){
				ncrnaAlignments[i]+=threadNcrnaAlignments[i];
			}
			
			stCds.add(pt.caller.stCds);
			stCds2.add(pt.caller.stCds2);
//			stCdsPass.add(pt.caller.stCdsPass);
			
			for(int i=0; i<trackers.length; i++){
				trackers[i].add(pt.caller.trackers[i]);
			}
			
			if(geneHist!=null){Tools.add(geneHist, pt.geneHistT);}
			
			success&=pt.success;
		}

		//Track whether any threads failed. NOTE: this is the CORRECT monotonic pattern (set errorState=true on failure, NEVER clear) -- unlike the template Accumulator tools' "errorState&=!success" (prok/MergeRibo#001 / MergeRibo_Fast#001) which clears prior errors when a later pass succeeds. CallGenes loops spawnThreads per input file (process L403) but stays correct because it never clears.
		if(!success){errorState=true;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Multipass           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates refined gene model using iterative multi-pass approach.
	 * Each pass improves model accuracy by incorporating predictions from previous passes.
	 *
	 * @param pgm0 Initial gene model
	 * @param fna Input FASTA filename
	 * @param gff Optional input GFF file with known annotations
	 * @param passes Number of refinement iterations
	 * @return Refined gene model with improved accuracy
	 */
	public static GeneModel makeMultipassModel(GeneModel pgm0, String fna, String gff, int passes/*, long maxReads*/) {
		if(passes<2) {return pgm0;}
		ArrayList<Read> reads=ReadInputStream.toReads(fna, FileFormat.FASTA, -1/*maxReads*/);
		return makeMultipassModel(pgm0, reads, gff, passes);
	}
	
	/**
	 * Creates refined gene model using iterative multi-pass approach with pre-loaded reads.
	 * Progressively mixes original and derived models with different weights per pass.
	 *
	 * @param pgm0 Initial gene model
	 * @param reads Pre-loaded sequence reads
	 * @param gff Optional input GFF file with known annotations
	 * @param passes Number of refinement iterations
	 * @return Refined gene model with weighted combination of predictions
	 */
	public static GeneModel makeMultipassModel(GeneModel pgm0, ArrayList<Read> reads, String gff, int passes) {
		//Self-refinement: each pass re-calls the caller on the genome's OWN predictions (runOnePass) to derive a genome-specific pgm, then PGMTools.mix blends original+derived with pass-dependent weights (early passes trust the derived model less; the final pass weights CDS 0.50 + rRNA/tRNA ~0.12). Bootstraps a tuned model from an unannotated genome -- the praise-worthy core.
		GeneModel pgm=pgm0;
		for(int i=1; i<passes; i++) {
			pgm=runOnePass(reads, gff, pgm);
			if(i==1 && passes==2) {//only pass for 2-pass
				pgm=PGMTools.mix(0.50, 0.10, 0.10, true, pgm0, pgm);
			}else if(i==passes-1) {//final pass for 3+ passes
				pgm=PGMTools.mix(0.50, 0.12, 0.12, true, pgm0, pgm);
			}else if(i==1 && passes>2) {//first pass for 3+ passes
				pgm=PGMTools.mix(0.10, 0.01, 0.01, true, pgm0, pgm);
			}else {//middle pass for 4+ passes
				pgm=PGMTools.mix(0.20, 0.02, 0.02, true, pgm0, pgm);
			}
		}
		return pgm;
	}
	
	/** This needs a pgm OR a gff, not both */
	public static GeneModel runOnePass(String fna, String gff, GeneModel pgm0/*, long maxReads*/) {//TODO: Make this multithreaded
		ArrayList<Read> reads=ReadInputStream.toReads(fna, FileFormat.FASTA, -1/*maxReads*/);
		return runOnePass(reads, gff, pgm0);
	}
	
	/** This needs a pgm OR a gff, not both */
	public static GeneModel runOnePass(ArrayList<Read> reads, String gff, GeneModel pgm0) {//TODO: Make this multithreaded
		GeneCaller caller=new GeneCaller(minLen, maxOverlapSameStrand, maxOverlapOppositeStrand, 
				minStartScore, minStopScore, minKmerScore, minOrfScore, minAvgScore, pgm0);
		
		final ArrayList<GffLine> cds, rrna, trna;
		if(gff!=null && !"null".equalsIgnoreCase(gff)) {
			ArrayList<GffLine>[] allGffLines=GffLine.loadGffFileByType(gff, "CDS,rRNA,tRNA", true);
			cds=allGffLines[0];
			rrna=allGffLines[1];
			trna=allGffLines[2];
//			System.err.println("Loaded "+cds.size()+" CDSs");
		}else {
			ArrayList<Orf> orfs=caller.callGenes(reads);
			ArrayList<GffLine>[] allGffLines=Orf.toGffLinesByType(orfs, "CDS,rRNA,tRNA");
			cds=allGffLines[0];
			rrna=allGffLines[1];
			trna=allGffLines[2];
//			System.err.println("Called "+cds.size()+" CDSs");
		}
		
		GeneModel pgm=new GeneModel(true);
		pgm.process(reads, cds, rrna, trna);
		for(StatsContainer sc : pgm.allContainers) {
			sc.calculate();//Not sure if this is needed...
		}
		return pgm;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates ByteStreamWriter for specified file format, or null if no output needed.
	 * @param ff FileFormat specification, or null
	 * @return Configured ByteStreamWriter, or null
	 */
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/**
	 * Creates ConcurrentReadOutputStream with optimized buffering.
	 * Buffer size is adjusted based on threading requirements and ordering constraints.
	 * @param ff FileFormat specification, or null
	 * @return Configured output stream, or null
	 */
	private ConcurrentReadOutputStream makeCros(FileFormat ff){
		if(ff==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=(ordered ? Tools.mid(4, 64, (Shared.threads()*2)/3) : 4);

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff, null, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		/**
		 * Constructs ProcessThread with shared I/O streams and gene model.
		 * Initializes GeneCaller with configured scoring parameters and creates
		 * thread-local gene length histogram if enabled.
		 *
		 * @param cris_ Shared input stream for sequences
		 * @param bsw_ Shared output writer for GFF annotations
		 * @param rosAmino_ Output stream for amino acid sequences
		 * @param ros16S_ Output stream for 16S rRNA sequences
		 * @param ros18S_ Output stream for 18S rRNA sequences
		 * @param pgm_ Gene model for predictions
		 * @param minLen Minimum gene length threshold
		 * @param tid_ Thread identifier
		 */
		ProcessThread(final ConcurrentReadInputStream cris_, final ByteStreamWriter bsw_, 
				ConcurrentReadOutputStream rosAmino_, ConcurrentReadOutputStream ros16S_, ConcurrentReadOutputStream ros18S_, 
				GeneModel pgm_, final int minLen, final int tid_){
			cris=cris_;
			bsw=bsw_;
			rosAmino=rosAmino_;
			ros16S=ros16S_;
			ros18S=ros18S_;
			pgm=pgm_;
			tid=tid_;
			geneHistT=(geneHistBins>1 ? new long[geneHistBins] : null);
			caller=new GeneCaller(minLen, maxOverlapSameStrand, maxOverlapOppositeStrand, 
					minStartScore, minStopScore, minKmerScore, minOrfScore, minAvgScore, pgm);
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(ln);

				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		/**
		 * Processes single list of reads through gene calling pipeline.
		 * Handles read merging/error correction if enabled, calls genes for each read,
		 * and writes results to appropriate output streams with proper ordering.
		 * @param ln ListNum containing reads to process
		 */
		void processList(ListNum<Read> ln){
			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;

//			System.err.println(reads.size());
			
			ArrayList<Orf> orfList=new ArrayList<Orf>();
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				Read r1=reads.get(idx);
				Read r2=r1.mate;
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}
				if(r2!=null && !r2.validated()){r2.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r1.length();
				final int initialLength2=r1.mateLength();

				//Increment counters
				readsInT+=r1.pairCount();
				basesInT+=initialLength1+initialLength2;
				
				if(r2!=null){
					if(merge){
						final int insert=BBMerge.findOverlapStrict(r1, r2, false);
						if(insert>0){
							r2.reverseComplementFast();
							r1=r1.joinRead(insert);
							r2=null;
						}
					}else if(ecco){
						BBMerge.findOverlapStrict(r1, r2, true);
					}
				}
				
				{
					//Reads are processed in this block.
					{
						ArrayList<Orf> list=processRead(r1);
						if(list!=null){orfList.addAll(list);}
					}
					if(r2!=null){
						ArrayList<Orf> list=processRead(r2);
						if(list!=null){orfList.addAll(list);}
					}
				}
			}
			
			genesOutT+=orfList.size();
			ByteBuilder bb=new ByteBuilder();
			
			if(bsw!=null){
				if(bsw.ordered){
					for(Orf orf : orfList){
						orf.appendGff(bb);
						bb.nl();
					}
					bsw.add(bb, ln.id);
					bytesOutT+=bb.length();
				}else{
					for(Orf orf : orfList){
						orf.appendGff(bb);
						bb.nl();
//						if(bb.length()>=32000){
//							bytesOutT+=bb.length();
//							bsw.addJob(bb);
//							bb=new ByteBuilder();
//						}
					}
					if(bb.length()>0){
						bsw.addJob(bb);
						bytesOutT+=bb.length();
					}
				}
			}

			//Notify the input stream that the list was used
			cris.returnList(ln);
//			if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept, false if they should be discarded.
		 */
		ArrayList<Orf> processRead(final Read r){
			ArrayList<Orf> list=caller.callGenes(r, pgm, true);
			
			if(geneHistT!=null && list!=null){
				for(Orf o : list){
					int bin=Tools.min(geneHistT.length-1, o.length()/geneHistDiv);
					geneHistT[bin]++;
				}
			}
			
			//[prok/CallGenes#002 mechanism] these sequence streams add by r.numericID ONLY when a gene of the type is found; a read with none leaves a gap -> ordered output hangs. (processList's GFF path always adds per ln.id, so it's safe.)
			if(ros16S!=null){
				if(list!=null && !list.isEmpty()){
//					System.err.println("Looking for 16S.");
					ArrayList<Read> ssu=fetchType(r, list, r16S);
					if(ssu!=null && !ssu.isEmpty()){
//						System.err.println("Found "+ssu.size()+" 16S.");
						ros16S.add(ssu, r.numericID);
					}
				}
			}
			if(ros18S!=null){
				if(list!=null && !list.isEmpty()){
					ArrayList<Read> ssu=fetchType(r, list, r18S);
					if(ssu!=null && !ssu.isEmpty()){ros18S.add(ssu, r.numericID);}
				}
			}
			
			if(rosAmino!=null){
				if(mode==TRANSLATE){
					if(list!=null && !list.isEmpty()){
						ArrayList<Read> prots=translate(r, list);
						if(prots!=null){rosAmino.add(prots, r.numericID);}
					}
				}else if(mode==RETRANSLATE) {
					if(list!=null && !list.isEmpty()){
						ArrayList<Read> prots=translate(r, list);
						ArrayList<Read> ret=detranslate(prots);
						if(ret!=null){rosAmino.add(ret, r.numericID);}
					}
				}else if(mode==RECODE) {
					if(list!=null && !list.isEmpty()){
						Read recoded=recode(r, list);
						r.mate=null;
						ArrayList<Read> rec=new ArrayList<Read>(1);
						rec.add(recoded);
						if(rec!=null){rosAmino.add(rec, r.numericID);}
					}
				}else{
					assert(false) : mode;
				}
			}
			
			return list;
		}
		
		/** Number of reads processed by this thread */
		protected long readsInT=0;
		/** Number of bases processed by this thread */
		protected long basesInT=0;
		
		/** Number of genes called by this thread */
		protected long genesOutT=0;
		/** Number of bytes written by this thread */
		protected long bytesOutT=0;
		
		/** Thread-local gene length histogram array */
		final long[] geneHistT;

		/** Output stream for amino acid sequences */
		protected ConcurrentReadOutputStream rosAmino;
		/** Output stream for 16S rRNA sequences */
		protected ConcurrentReadOutputStream ros16S;
		/** Output stream for 18S rRNA sequences */
		protected ConcurrentReadOutputStream ros18S;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ByteStreamWriter bsw;
		/** Gene Model for annotation (not really needed) */
		private final GeneModel pgm;
		/** Gene Caller for annotation */
		final GeneCaller caller;
		/** Thread ID */
		final int tid;
	}
	
	/**
	 * Extracts sequences of specified gene type from ORF list.
	 * Processes both strands of the read, reverse complementing as needed.
	 *
	 * @param r Source read containing sequence data
	 * @param list List of ORFs to search
	 * @param type Gene type constant (r16S, r18S, tRNA, etc.)
	 * @return List of Read objects containing extracted sequences of specified type
	 */
	public static ArrayList<Read> fetchType(final Read r, final ArrayList<Orf> list, int type){
		if(list==null || list.isEmpty()){return null;}
		ArrayList<Read> ret=new ArrayList<Read>(list.size());
		for(int strand=0; strand<2; strand++){
			for(Orf orf : list){
				if(orf.strand==strand && orf.type==type){
					Read sequence=fetch(orf, r.bases, r.id);
					ret.add(sequence);
				}
			}
			r.reverseComplement();
		}
		return (ret.size()>0 ? ret : null);
	}
	
	/**
	 * Translates coding sequences (CDS) from ORF list to amino acids.
	 * Processes both strands and creates Read objects with amino acid sequences.
	 *
	 * @param r Source read containing nucleotide sequence
	 * @param list List of ORFs to translate
	 * @return List of Read objects containing amino acid translations
	 */
	public static ArrayList<Read> translate(final Read r, final ArrayList<Orf> list){
		if(list==null || list.isEmpty()){return null;}
		ArrayList<Read> prots=new ArrayList<Read>(list.size());
		for(int strand=0; strand<2; strand++){
			for(Orf orf : list){
				if(orf.strand==strand && orf.type==CDS){
					Read aa=translate(orf, r.bases, r.id);
					prots.add(aa);
				}
			}
			r.reverseComplement();
		}
		return prots.isEmpty() ? null : prots;
	}
	
	/**
	 * Recodes CDS regions in sequence using canonical codons.
	 * Replaces existing codons with canonical versions while preserving amino acids.
	 *
	 * @param r Read to recode in-place
	 * @param list List of ORFs defining coding regions
	 * @return Original read with recoded CDS regions
	 */
	public static Read recode(final Read r, final ArrayList<Orf> list){
		if(list==null || list.isEmpty()){return r;}
		for(int strand=0; strand<2; strand++){
			for(Orf orf : list){
				if(orf.strand==strand && orf.type==CDS){
					recode(orf, r.bases);
				}
			}
			r.reverseComplement();
		}
		return r;
	}
	
	/**
	 * Converts amino acid sequences back to nucleotide sequences using canonical codons.
	 * Each amino acid is replaced with its corresponding canonical codon triplet.
	 * @param prots List of protein (amino acid) sequences
	 * @return List of nucleotide sequences using canonical codons
	 */
	public static ArrayList<Read> detranslate(final ArrayList<Read> prots){
		if(prots==null || prots.isEmpty()){return null;}
		ArrayList<Read> nucs=new ArrayList<Read>(prots.size());
		for(int strand=0; strand<2; strand++){
			for(Read prot : prots){
				Read nuc=detranslate(prot);
				nucs.add(nuc);
			}
		}
		return nucs;
	}
	
	/**
	 * Translates single ORF to amino acid sequence.
	 * Handles strand orientation and creates Read with appropriate metadata.
	 *
	 * @param orf ORF defining coding region boundaries
	 * @param bases Source nucleotide sequence
	 * @param id Base identifier for naming translated sequence
	 * @return Read containing amino acid translation with coordinate metadata
	 */
	public static Read translate(Orf orf, byte[] bases, String id){
//		assert(orf.length()%3==0) : orf.length(); //Happens sometimes on genes that go off the end, perhaps
		if(orf.strand==1){orf.flip();}
		byte[] acids=AminoAcid.toAAs(bases, orf.start, orf.stop);
		if(orf.strand==1){orf.flip();}
		Read r=new Read(acids, null, id+"\t"+(Shared.strandCodes[orf.strand]+"\t"+orf.start+"-"+orf.stop), 0, Read.AAMASK);
//		assert((r.length()+1)*3==orf.length());
		return r;
	}
	
	/**
	 * Extracts nucleotide sequences for all CDS ORFs from both strands.
	 * Creates Read objects containing the raw nucleotide sequences of coding regions.
	 *
	 * @param r Source read containing sequence data
	 * @param list List of ORFs to extract
	 * @return List of Read objects containing CDS nucleotide sequences
	 */
	public static ArrayList<Read> fetch(final Read r, final ArrayList<Orf> list){
		if(list==null || list.isEmpty()){return null;}
		ArrayList<Read> genes=new ArrayList<Read>(list.size());
		for(int strand=0; strand<2; strand++){
			for(Orf orf : list){
				if(orf.strand==strand && orf.type==CDS){
					Read gene=fetch(orf, r.bases, r.id);
					genes.add(gene);
				}
			}
			r.reverseComplement();
		}
		return genes.isEmpty() ? null : genes;
	}
	
	/**
	 * Extracts nucleotide sequence for single ORF from Read.
	 * Handles reverse complement if needed and validates coordinate bounds.
	 *
	 * @param orf ORF defining sequence boundaries
	 * @param source Source Read containing sequence data
	 * @return Read containing extracted sequence with coordinate metadata
	 */
	public static Read fetch(Orf orf, Read source){
		assert(orf.start>=0 && orf.stop<source.length() && orf.length()>0) : 
			source.length()+"\n"+orf.length()+"\n"+orf;
		if(orf.strand==1){source.reverseComplement();}
		Read r=fetch(orf, source.bases, source.id);
		if(orf.strand==1){source.reverseComplement();}
		return r;
	}
	
	/**
	 * Extracts nucleotide sequence for single ORF from byte array.
	 * Creates new Read with subsequence and coordinate annotations.
	 *
	 * @param orf ORF defining sequence boundaries and strand
	 * @param bases Source nucleotide sequence array
	 * @param id Base identifier for naming extracted sequence
	 * @return Read containing extracted sequence with metadata
	 */
	public static Read fetch(Orf orf, byte[] bases, String id){
		assert(orf.start>=0 && orf.stop<bases.length) : bases.length+"\n"+orf;
		if(orf.strand==1){orf.flip();}
		byte[] sub=Arrays.copyOfRange(bases, orf.start, orf.stop+1);
		if(orf.strand==1){orf.flip();}
		Read r=new Read(sub, null, id+"\t"+(Shared.strandCodes[orf.strand]+"\t"+orf.start+"-"+orf.stop), 0, 0);
		assert(r.length()==orf.length()) : r.length()+", "+orf.length();
		return r;
	}
	
	/**
	 * Recodes ORF region in nucleotide array using canonical codons.
	 * Translates to amino acids then back to canonical codon sequences in-place.
	 * @param orf ORF defining coding region to recode
	 * @param bases Nucleotide array to modify in-place
	 */
	public static void recode(Orf orf, byte[] bases){
		if(orf.strand==1){orf.flip();}
		byte[] acids=AminoAcid.toAAs(bases, orf.start, orf.stop);
		for(int apos=0, bpos=orf.start; apos<acids.length; apos++){
			byte aa=acids[apos];
			int number=AminoAcid.acidToNumber[aa];
			String codon=(number>=0 ? AminoAcid.canonicalCodons[number] : "NNN");
			for(int i=0; i<3; i++, bpos++) {
				bases[bpos]=(byte)codon.charAt(i);
			}
		}
		if(orf.strand==1){orf.flip();}
	}
	
	/**
	 * Converts amino acid sequence to nucleotide sequence using canonical codons.
	 * Each amino acid becomes a three-nucleotide canonical codon sequence.
	 * @param prot Read containing amino acid sequence
	 * @return Read containing equivalent nucleotide sequence
	 */
	public static Read detranslate(Read prot){
		ByteBuilder bb=new ByteBuilder(prot.length()*3);
		for(byte aa : prot.bases){
			int number=AminoAcid.acidToNumber[aa];
			String codon=(number>=0 ? AminoAcid.canonicalCodons[number] : "NNN");
			bb.append(codon);
		}
		Read r=new Read(bb.array, null, prot.id, prot.numericID, 0);
		return r;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Factory method for creating GeneCaller with current static parameters.
	 * Uses configured scoring thresholds and overlap limits.
	 * @param pgm Gene model for probability calculations
	 * @return Configured GeneCaller instance
	 */
	public static GeneCaller makeGeneCaller(GeneModel pgm){
		GeneCaller caller=new GeneCaller(minLen, maxOverlapSameStrand, maxOverlapOppositeStrand,
				minStartScore, minStopScore, minKmerScore, minOrfScore, minAvgScore, pgm);
		return caller;
	}

	public static GeneCaller makeGeneCallerForPhylum(String phylum){
		GeneModel pgm=getPhylumPGM(phylum);
		return makeGeneCaller(pgm);
	}

	public static synchronized GeneModel getPhylumPGM(String phylum){
		if(defaultPGM==null){
			defaultPGM=GeneModelParser.loadModel(Data.findPath("?model.pgm"));
		}
		if(phylum==null){return defaultPGM;}
		GeneModel cached=pgmCache.get(phylum);
		if(cached!=null){return cached;}
		String path=Data.findPath("?pgm_"+phylum+".pgm", false);
		if(path!=null && new java.io.File(path).exists()){
			GeneModel pgm=GeneModelParser.loadModel(path);
			pgmCache.put(phylum, pgm);
			return pgm;
		}
		pgmCache.put(phylum, defaultPGM);
		return defaultPGM;
	}

	private static final java.util.HashMap<String, GeneModel> pgmCache=new java.util.HashMap<>();
	private static GeneModel defaultPGM;

	public static synchronized void loadTrnaResources(){
		if(GeneCaller.trnaLibrary!=null){return;}
		String libPath=Data.findPath("?trna_consensus.fa");
		String modelPath=Data.findPath("?trna_models.hbm");
		if(libPath!=null && new java.io.File(libPath).exists()){
			GeneCaller.trnaLibrary=TrnaConsensusBuilder.loadLibrary(libPath);
			GeneCaller.trnaModelNames=TrnaConsensusBuilder.lastLoadedNames;
		}
		if(modelPath!=null && new java.io.File(modelPath).exists()){
			GeneCaller.trnaModels=TrnaConsensusBuilder.loadModels(modelPath);
		}
	}

	/**
	 * Loads the generic ncRNA family bundles (RNase P, SRP small/large, and the
	 * separately gated experimental tmRNA family) into
	 * GeneCaller.ncrnaFamilies, one entry per family. Forward-ported from Noire's
	 * ncRNA-family-loading tree (2026-08-28, C3 merge) -- see that tree's original
	 * javadoc for the resource-naming rationale (hardcoded for now; srp_small/srp_large
	 * share one kmer set but have separate consensus libraries). ONE difference from
	 * Noire's original: the call site (in parse(), above) is now gated behind
	 * NCRNA_FAMILIES_ENABLED (Gate A) instead of unconditional -- Brian shipped ncRNAs
	 * OFF by default (Citan's explicit requirement), so this method itself stays
	 * callable directly (e.g. by NcrnaBoundaryVectorGen's training-vector generator, or
	 * eval scripts) but production callgenes invocations only reach it when ncrna=t.
	 * Idempotent via GeneCaller.ncrnaFamilies being non-empty.
	 *
	 * <p>C3 ADDITION (G11, 2026-08-28): each family also attempts to lazy-load its own
	 * boundary-precision-NN resources (net + start/stop tables + a matching HBM model
	 * file) when NCRNA_BOUNDARY_NN_ENABLED (Gate B) is on -- see addNcrnaFamily. Gate B
	 * is independent of and subordinate to Gate A: with Gate A off this method never
	 * runs at all (and Gate B being on without Gate A is rejected loudly at parse time,
	 * see the parse() assert); with Gate A on and Gate B off, families load with
	 * boundary refinement structurally disabled (null templates). With BOTH on, Gate B
	 * is STRICT, not tolerant (Citan's correction, 2026-08-28): every loaded family must
	 * have complete, index-aligned boundary resources, or addNcrnaFamily throws --
	 * a silent per-family skip would let a run report "boundary NN enabled" while only
	 * a subset of families actually got refined, an invalid mixed on/off comparison.
	 */
	public static synchronized void loadNcrnaResources(){
		//TODO: Probable bug - this production loader silently reuses an existing family set on
		//a second call, even if the requested ncrna/tmrna/sixs configuration differs. The
		//evaluation harness snapshots and restores this static state, but a persistent caller
		//should fail loud rather than emit calls using a stale resource configuration.
		if(!GeneCaller.ncrnaFamilies.isEmpty()){return;}
		LongHashSet rnasepKmers=loadEffectiveNcrnaKmerSet("rnasep", "rnasep_17mers.fa", RNASEP_KMERS_OVERRIDE, 17);
		addNcrnaFamily("rnasep", "rnasep_consensus.fa", "rnasep_models.hbm", rnasepKmers, 17, 200,
			resolveSweepPad("rnasep", RNASEP_PAD_OVERRIDE, 320),
			7, 60, true, 100f, 0.74f, 0.072f, 12,
			resolveFamilyScore("rnasep", RNASEP_SCORE_A_OVERRIDE, NCRNA_SCORE_A_OVERRIDE, 0f),
			resolveFamilyScore("rnasep", RNASEP_SCORE_B_OVERRIDE, NCRNA_SCORE_B_OVERRIDE, 1.6f),
			resolveSweepFloat("rnasep", NCRNA_ID_PASS_OVERRIDE, 0.65f), resolveSweepFloat("rnasep", NCRNA_ID_BORDERLINE_OVERRIDE, 0.58f),
			resolveSweepFloat("rnasep", NCRNA_HBM_PASS_OVERRIDE, 0.62f), resolveSweepFloat("rnasep", NCRNA_COLLAPSE_FRAC_OVERRIDE, 0.85f),
			boundaryStartOffsets("rnasep"), boundaryStopOffsets("rnasep"));
		LongHashSet srpSmallKmers=loadEffectiveNcrnaKmerSet("srp_small", "srp_17mers.fa", SRPSMALL_KMERS_OVERRIDE, 17);
		addNcrnaFamily("srp_small", "srp_small_consensus.fa", "srp_small_models.hbm", srpSmallKmers, 17, 80,
			resolveSweepPad("srp_small", SRPSMALL_PAD_OVERRIDE, 75),
			7, 60, true, 23f, 0.68f, 0.072f, 12,
			resolveFamilyScore("srp_small", SRPSMALL_SCORE_A_OVERRIDE, NCRNA_SCORE_A_OVERRIDE, 0f),
			resolveFamilyScore("srp_small", SRPSMALL_SCORE_B_OVERRIDE, NCRNA_SCORE_B_OVERRIDE, 1.8f),
			resolveSweepFloat("srp_small", NCRNA_ID_PASS_OVERRIDE, 0.70f), resolveSweepFloat("srp_small", NCRNA_ID_BORDERLINE_OVERRIDE, 0.64f),
			resolveSweepFloat("srp_small", NCRNA_HBM_PASS_OVERRIDE, 0.70f), resolveSweepFloat("srp_small", NCRNA_COLLAPSE_FRAC_OVERRIDE, 0.85f),
			boundaryStartOffsets("srp_small"), boundaryStopOffsets("srp_small"));
		LongHashSet srpLargeKmers=loadEffectiveNcrnaKmerSet("srp_large", "srp_17mers.fa", SRPLARGE_KMERS_OVERRIDE, 17);
		addNcrnaFamily("srp_large", "srp_large_consensus.fa", "srp_large_models.hbm", srpLargeKmers, 17, 150,
			resolveSweepPad("srp_large", SRPLARGE_PAD_OVERRIDE, 250),
			7, 60, true, 45f, 0.75f, 0.072f, 12,
			resolveFamilyScore("srp_large", SRPLARGE_SCORE_A_OVERRIDE, NCRNA_SCORE_A_OVERRIDE, 0f),
			resolveFamilyScore("srp_large", SRPLARGE_SCORE_B_OVERRIDE, NCRNA_SCORE_B_OVERRIDE, 1.0f),
			resolveSweepFloat("srp_large", NCRNA_ID_PASS_OVERRIDE, 0.70f), resolveSweepFloat("srp_large", NCRNA_ID_BORDERLINE_OVERRIDE, 0.64f),
			resolveSweepFloat("srp_large", NCRNA_HBM_PASS_OVERRIDE, 0.66f), resolveSweepFloat("srp_large", NCRNA_COLLAPSE_FRAC_OVERRIDE, 0.85f),
			boundaryStartOffsets("srp_large"), boundaryStopOffsets("srp_large"));
		if(TMRNA_ENABLED){
			final int before=GeneCaller.ncrnaFamilies.size();
			LongHashSet tmrnaKmers=loadEffectiveNcrnaKmerSet("tmrna", "tmrna_17mers.fa", TMRNA_KMERS_OVERRIDE, 17);
			addNcrnaFamily("tmrna", (TMRNA_CONSENSUS_OVERRIDE==null ? "tmrna_consensus.fa" : TMRNA_CONSENSUS_OVERRIDE),
				(TMRNA_MODELS_OVERRIDE==null ? "tmrna_models.hbm" : TMRNA_MODELS_OVERRIDE), tmrnaKmers, 17, 180,
				resolveSweepPad("tmrna", TMRNA_PAD_OVERRIDE, 370),
				7, 60, true, 11f, 0.48f, 0.072f, 12,
				resolveFamilyScore("tmrna", TMRNA_SCORE_A_OVERRIDE, NCRNA_SCORE_A_OVERRIDE, 5f),
				resolveFamilyScore("tmrna", TMRNA_SCORE_B_OVERRIDE, NCRNA_SCORE_B_OVERRIDE, 3.0f),
				resolveSweepFloat("tmrna", NCRNA_ID_PASS_OVERRIDE, 0.62f), resolveSweepFloat("tmrna", NCRNA_ID_BORDERLINE_OVERRIDE, 0.60f),
				resolveSweepFloat("tmrna", NCRNA_HBM_PASS_OVERRIDE, 0.62f), resolveSweepFloat("tmrna", NCRNA_COLLAPSE_FRAC_OVERRIDE, 0.85f),
				boundaryStartOffsets("tmrna"), boundaryStopOffsets("tmrna"));
			if(GeneCaller.ncrnaFamilies.size()!=before+1){
				throw new RuntimeException("tmrna=t requires tmrna_17mers.fa[.gz], tmrna_consensus.fa, and tmrna_models.hbm resources");
			}
			final NcrnaFamily f=GeneCaller.ncrnaFamilies.get(before);
			if(f.models==null || f.modelNames==null || f.models.length!=f.library.length || f.modelNames.length!=f.library.length){
				throw new RuntimeException("tmrna=t requires index-aligned tmrna_consensus.fa and tmrna_models.hbm resources");
			}
		}
		if(SIXS_ENABLED){
			final int before=GeneCaller.ncrnaFamilies.size();
			LongHashSet rf00013Kmers=loadEffectiveNcrnaKmerSet("sixs_rf00013", "sixs_rf00013_16mers.fa", SIXS_RF00013_KMERS_OVERRIDE, 16);
			addNcrnaFamily("sixs_rf00013", "sixs_rf00013_consensus.fa", "sixs_rf00013_models.hbm", rf00013Kmers, 16, 60,
				resolveSweepPad("sixs_rf00013", SIXS_RF00013_PAD_OVERRIDE, 250),
				7, 60, false, 0f, 0f, 0f, 1,
				resolveSweepFloat("sixs_rf00013", NCRNA_SCORE_A_OVERRIDE, 0f),
				resolveSweepFloat("sixs_rf00013", NCRNA_SCORE_B_OVERRIDE, 20f),
				resolveSweepFloat("sixs_rf00013", NCRNA_ID_PASS_OVERRIDE, 0.70f), resolveSweepFloat("sixs_rf00013", NCRNA_ID_BORDERLINE_OVERRIDE, 0.60f),
				resolveSweepFloat("sixs_rf00013", NCRNA_HBM_PASS_OVERRIDE, 0.64f), resolveSweepFloat("sixs_rf00013", NCRNA_COLLAPSE_FRAC_OVERRIDE, 0.85f),
				boundaryStartOffsets("sixs_rf00013"), boundaryStopOffsets("sixs_rf00013"));
			LongHashSet rf01685Kmers=loadEffectiveNcrnaKmerSet("sixs_rf01685", "sixs_rf01685_17mers.fa", SIXS_RF01685_KMERS_OVERRIDE, 17);
			addNcrnaFamily("sixs_rf01685", "sixs_rf01685_consensus.fa", "sixs_rf01685_models.hbm", rf01685Kmers, 17, 60,
				resolveSweepPad("sixs_rf01685", SIXS_RF01685_PAD_OVERRIDE, 125),
				7, 60, false, 0f, 0f, 0f, 1,
				resolveSweepFloat("sixs_rf01685", NCRNA_SCORE_A_OVERRIDE, 0f),
				resolveSweepFloat("sixs_rf01685", NCRNA_SCORE_B_OVERRIDE, 20f),
				resolveSweepFloat("sixs_rf01685", NCRNA_ID_PASS_OVERRIDE, 0.70f), resolveSweepFloat("sixs_rf01685", NCRNA_ID_BORDERLINE_OVERRIDE, 0.70f),
				resolveSweepFloat("sixs_rf01685", NCRNA_HBM_PASS_OVERRIDE, 0.75f), resolveSweepFloat("sixs_rf01685", NCRNA_COLLAPSE_FRAC_OVERRIDE, 0.85f),
				boundaryStartOffsets("sixs_rf01685"), boundaryStopOffsets("sixs_rf01685"));
			if(GeneCaller.ncrnaFamilies.size()!=before+2){
				throw new RuntimeException("sixs=t requires complete RF00013 and RF01685 consensus, HBM, and kmer resources");
			}
			for(int i=before; i<before+2; i++){
				final NcrnaFamily f=GeneCaller.ncrnaFamilies.get(i);
				if(f.models==null || f.modelNames==null || f.models.length!=f.library.length || f.modelNames.length!=f.library.length){
					throw new RuntimeException("sixs=t requires index-aligned consensus and HBM resources for "+f.name);
				}
			}
		}
	}

	/** Measured width-six candidate windows. Kept in one factory used by family loading;
	 * unknown families fail loud so a future pipeline cannot silently inherit tRNA geometry. */
	static int[] boundaryStartOffsets(String family){
		if(family.equals("rnasep") || family.equals("srp_small") || family.equals("srp_large") || family.equals("tmrna")){
			return new int[]{-3,-2,-1,0,1,2};
		}
		// No sixS boundary models are shipped yet.  These candidate positions keep the
		// family geometry explicit; Gate B remains strict and fails before use until its
		// family-specific net and tables are released.
		if(family.equals("sixs_rf00013") || family.equals("sixs_rf01685")){return new int[]{-3,-2,-1,0,1,2};}
		throw new IllegalArgumentException("No boundary-start offsets configured for ncRNA family: "+family);
	}

	static int[] boundaryStopOffsets(String family){
		if(family.equals("rnasep")){return new int[]{-1,0,1,2,3,4};}
		if(family.equals("srp_small") || family.equals("srp_large")){return new int[]{-2,-1,0,1,2,3};}
		if(family.equals("tmrna")){return new int[]{-3,-2,-1,0,1,2};}
		// No sixS boundary models are shipped yet.  These candidate positions keep the
		// family geometry explicit; Gate B remains strict and fails before use until its
		// family-specific net and tables are released.
		if(family.equals("sixs_rf00013") || family.equals("sixs_rf01685")){return new int[]{-3,-2,-1,0,1,2};}
		throw new IllegalArgumentException("No boundary-stop offsets configured for ncRNA family: "+family);
	}

	/** Parses and validates one padding-override CLI value (P1, 2026-08-28). Package-visible,
	 * factored out of parse() for direct testability without constructing a full CallGenes
	 * instance -- mirrors the validateNcrnaGateCombo/requireCompleteBoundaryResources pattern
	 * already used for the other C3 gates. */
	static int parsePadOverride(String flagName, String value){
		final int pad=Integer.parseInt(value);
		if(pad<0){throw new IllegalArgumentException(flagName+" must be >=0: "+value);}
		return pad;
	}

	/** Resolves one family's effective window pad: the override if set (>=0), else the
	 * family's shipped default. -1 is the "unset" sentinel (0 is a legitimate override
	 * value, so it cannot double as "unset"). Package-visible for direct testability. */
	static int resolvePad(int override, int shippedDefault){
		return override>=0 ? override : shippedDefault;
	}

	static int resolveSweepPad(String family, int familyOverride, int shippedDefault){
		return ncrnaSweepTarget(family) && NCRNA_WINDOW_PAD_OVERRIDE>=0
			? NCRNA_WINDOW_PAD_OVERRIDE : resolvePad(familyOverride, shippedDefault);
	}

	static float resolveSweepFloat(String family, float override, float shippedDefault){
		return ncrnaSweepTarget(family) && !Float.isNaN(override) ? override : shippedDefault;
	}

	/** Resolves an explicit per-family score before the legacy generic targeted-sweep
	 * override.  This permits all families' scores to be pinned in one run, avoiding
	 * cross-family arbitration confounding while retaining the old single-family CLI. */
	static float resolveFamilyScore(String family, float familyOverride,
			float genericOverride, float shippedDefault){
		return !Float.isNaN(familyOverride) ? familyOverride
			: resolveSweepFloat(family, genericOverride, shippedDefault);
	}

	/** Loads one family bundle from its own kmer-set resource file, with per-family index params. */
	/** Loads one family bundle, given an ALREADY-LOADED kmer set (for families that share one,
	 * e.g. srp_small/srp_large both seeding from srp_17mers.fa). Skips silently (no bundle added)
	 * if the library or kmer set resource can't be found -- see loadNcrnaResources' javadoc.
	 * C3 addition (G11, 2026-08-28): also attempts the family's boundary net/tables/models when
	 * NCRNA_BOUNDARY_NN_ENABLED -- resource names "<name>_boundary_net.bbnet",
	 * "<name>_boundary_start_table.tsv", "<name>_boundary_stop_table.tsv", plus the family's
	 * existing HBM model file, mirroring the tRNA boundary resource convention
	 * (?trna_boundary_net.bbnet etc.) one level down per family. STRICT, not tolerant, unlike
	 * the library/kmer-set resources above (Citan's correction, 2026-08-28): once Gate B is on,
	 * EVERY loaded family must have complete, index-aligned boundary resources or this method
	 * throws a RuntimeException naming the family, the missing/misaligned resource, and its
	 * path -- a silent per-family skip would let a benchmark report "boundary NN enabled" while
	 * actually only refining a subset of families, a methodologically invalid mixed on/off
	 * result. The library/kmer-set tolerance above is unrelated and unchanged: a family with no
	 * consensus library at all just doesn't load (visible as a missing family, not a silent
	 * partial-feature state). */
	private static void addNcrnaFamily(String name, String libResource, String modelResource,
			LongHashSet kmerSet, int kLong, int minLen, int windowPad,
			int indexK, int indexTopN, boolean adaptive,
			float adaptFloor, float adaptTopFrac, float adaptQFrac, int fixedMinHits,
			float scoreA, float scoreB, float idPass, float idBorderline,
			float hbmPass, float collapseFrac, int[] boundaryStartOffsets, int[] boundaryStopOffsets){
		if(kmerSet==null){return;}
		String libPath=findNcrnaResource(libResource);
		if(libPath==null || !new java.io.File(libPath).exists()){return;}
		byte[][] library=TrnaConsensusBuilder.loadLibrary(libPath);
		String[] modelNames=TrnaConsensusBuilder.lastLoadedNames;
		consensus.BaseGraph[] models=null;
		String modelPath=findNcrnaResource(modelResource);
		if(modelPath!=null && new java.io.File(modelPath).exists()){
			models=TrnaConsensusBuilder.loadModels(modelPath);
		}

		ml.CellNet boundaryNet=null;
		prok.TrnaBoundaryFeatures.NinemerTable boundaryStartTable=null, boundaryStopTable=null;
		int boundaryStartInside=-1, boundaryStartOutside=-1, boundaryStopInside=-1, boundaryStopOutside=-1;
		float boundaryMeanLen=0f;
		//Citan's correction (2026-08-28): Gate B must NEVER silently degrade. The original
		//version of this block tolerantly skipped boundary loading per-family when a resource
		//was missing -- reasonable for the LIBRARY/MODEL resources above (a family with no
		//consensus library genuinely can't be called at all, and that's already visible as a
		//missing family), but WRONG for boundary resources under Gate B specifically: a
		//benchmark run reports "boundary NN enabled" as a single flag, and a silent per-family
		//skip would produce a MIXED on/off result reported as fully enabled -- exactly the kind
		//of methodologically-invalid comparison Citan flagged. So once Gate B is on, EVERY
		//family that loads (via Gate A) MUST get complete boundary resources, or the whole run
		//fails loud identifying which family and which paths are missing -- not a per-family
		//silent fallback.
		if(NCRNA_BOUNDARY_NN_ENABLED){
			String netPath=Data.findPath("?"+name+"_boundary_net.bbnet", false);
			String startTablePath=Data.findPath("?"+name+"_boundary_start_table.tsv", false);
			String stopTablePath=Data.findPath("?"+name+"_boundary_stop_table.tsv", false);
			final boolean netOk=(netPath!=null && new java.io.File(netPath).exists());
			final boolean startOk=(startTablePath!=null && new java.io.File(startTablePath).exists());
			final boolean stopOk=(stopTablePath!=null && new java.io.File(stopTablePath).exists());
			requireCompleteBoundaryResources(name, netPath, startTablePath, stopTablePath, netOk, startOk, stopOk,
				models, library, modelNames, libPath, modelPath);
			boundaryNet=NcrnaBoundaryScorer.load(netPath);
			TrnaNinemerTableBuilder.LoadedTable lt1=TrnaNinemerTableBuilder.loadTable(startTablePath);
			TrnaNinemerTableBuilder.LoadedTable lt2=TrnaNinemerTableBuilder.loadTable(stopTablePath);
			assert(lt1.type==prok.TrnaBoundaryFeatures.BoundaryType.START) : name+"_boundary_start_table.tsv ("
				+startTablePath+") is not labeled a start-boundary table -- wrong file, or start/stop swapped.";
			assert(lt2.type==prok.TrnaBoundaryFeatures.BoundaryType.STOP) : name+"_boundary_stop_table.tsv ("
				+stopTablePath+") is not labeled a stop-boundary table -- wrong file, or start/stop swapped.";
			boundaryStartTable=lt1.table; boundaryStartInside=lt1.insideCount; boundaryStartOutside=lt1.outsideCount;
			boundaryStopTable=lt2.table; boundaryStopInside=lt2.insideCount; boundaryStopOutside=lt2.outsideCount;
			boundaryMeanLen=medianLength(library);
		}
		GeneCaller.ncrnaFamilies.add(new NcrnaFamily(name, library, models, modelNames, kmerSet, kLong, minLen, windowPad,
			indexK, indexTopN, adaptive, adaptFloor, adaptTopFrac, adaptQFrac, fixedMinHits,
			scoreA, scoreB, idPass, idBorderline, hbmPass, collapseFrac,
			boundaryNet, boundaryNet, boundaryStartTable, boundaryStopTable,
			boundaryStartInside, boundaryStartOutside, boundaryStopInside, boundaryStopOutside,
			boundaryMeanLen, boundaryStartOffsets, boundaryStopOffsets));
	}

	/** Resolves a literal development/sweep path when it exists or names a path, otherwise
	 * resolves the normal shipped resource name. This keeps production names unchanged while
	 * making a complete consensus/HBM A/B reproducible without swapping files between runs. */
	private static String findNcrnaResource(String pathOrResource){
		if(pathOrResource==null){return null;}
		final java.io.File f=new java.io.File(pathOrResource);
		if(f.exists() || pathOrResource.indexOf('/')>=0 || pathOrResource.indexOf('\\')>=0){return pathOrResource;}
		return Data.findPath("?"+pathOrResource, false);
	}

	/** Gate B's fail-loud resource-completeness+alignment check (Citan's corrections,
	 * 2026-08-28), factored out of addNcrnaFamily into its own package-visible, directly
	 * testable method: given the already-resolved net/table existence flags and the loaded
	 * library/models/modelNames, throws a RuntimeException naming the family and every
	 * failing check if boundary resources are incomplete OR the library/model pairing is
	 * uncheckable/misaligned; returns silently (no exception) if everything is complete and
	 * verified aligned. Package-visible (not public) purely for CallGenesNcrnaGateTest --
	 * not part of any public API contract.
	 *
	 * <p>Checks, in order: (1) net/start-table/stop-table files exist; (2) models and
	 * modelNames are both non-null and length-matched to library (a malformed resource that
	 * failed to parse any names/models is caught here, before touching index i); (3) neither
	 * side's first-whitespace-token-normalized names contain an empty token or a duplicate
	 * (Citan's uniqueness requirement -- a duplicate could hide a reorder behind a false
	 * per-index "match"); (4) per-index normalized-name equality between the FASTA header
	 * (modelNames[i]) and the HBM model name (models[i].name). First-token normalization
	 * (not raw full-string equality) specifically to stay immune to the global trd=/
	 * trimreaddescription flag, which truncates modelNames[i] (Read.id) but never touches
	 * models[i].name (the HBM parser doesn't go through Read parsing at all) -- see the call
	 * site's original comment history for the full trace, preserved here rather than
	 * duplicated at both the call site and this method. */
	static void requireCompleteBoundaryResources(String name, String netPath, String startTablePath,
			String stopTablePath, boolean netOk, boolean startOk, boolean stopOk,
			consensus.BaseGraph[] models, byte[][] library, String[] modelNames,
			String libPath, String modelPath){
		final boolean countsOk=(models!=null && models.length==library.length
			&& modelNames!=null && modelNames.length==library.length);
		String duplicateDetail=null;
		if(countsOk){
			duplicateDetail=findDuplicateOrEmptyToken("fasta", modelNames);
			if(duplicateDetail==null){duplicateDetail=findDuplicateOrEmptyToken("hbm", modelNamesOf(models));}
		}
		final int[] mismatchIndex={-1};
		if(countsOk && duplicateDetail==null){
			for(int i=0; i<library.length; i++){
				final String fastaTok=firstToken(modelNames[i]);
				final String hbmTok=firstToken(models[i].name);
				if(!fastaTok.equals(hbmTok)){mismatchIndex[0]=i; break;}
			}
		}
		final boolean modelsOk=(countsOk && duplicateDetail==null && mismatchIndex[0]<0);
		if(!netOk || !startOk || !stopOk || !modelsOk){
			final String indexDetail=(!countsOk ? "count mismatch"
				: duplicateDetail!=null ? "normalized-token uniqueness violated: "+duplicateDetail
				: mismatchIndex[0]>=0 ? "index "+mismatchIndex[0]+" name mismatch: fasta='"
					+modelNames[mismatchIndex[0]]+"' vs hbm='"+models[mismatchIndex[0]].name+"'"
				: "ok");
			throw new RuntimeException("ncrnaboundarynet=t is on but family '"+name+"' is missing required "
				+"boundary resources, or its library/model pairing is invalid -- refusing to run a mixed "
				+"on/off benchmark silently. "
				+"net("+(name+"_boundary_net.bbnet")+")="+netPath+" exists="+netOk+"; "
				+"startTable("+(name+"_boundary_start_table.tsv")+")="+startTablePath+" exists="+startOk+"; "
				+"stopTable("+(name+"_boundary_stop_table.tsv")+")="+stopTablePath+" exists="+stopOk+"; "
				+"models="+(models==null ? "null" : models.length+" entries")+" vs library="+library.length
				+" entries; per-index name alignment ("+libPath+" vs "+modelPath+"): "+indexDetail+". "
				+"Stage all four resources for this family (net, start table, stop table, an HBM model "
				+"file index-aligned with the consensus library), or turn ncrnaboundarynet=f off.");
		}
	}

	/** Resolves+loads one kmer-set resource via ProkObject's shared "?prefix_kmers.fa[.gz]"
	 * logic, or null (with a stderr note, matching loadLongKmersByType's own behavior) if
	 * the resource isn't staged. kmerResource already includes the "_<k>mers.fa" suffix, so
	 * this strips it back to the bare prefix ProkObject.loadLongKmersByType expects. */
	private static LongHashSet loadNcrnaKmerSet(String kmerResource, int kLong){
		final String suffix="_"+kLong+"mers.fa";
		if(!kmerResource.endsWith(suffix)){
			throw new RuntimeException("kmerResource '"+kmerResource+"' doesn't match expected suffix '"+suffix+"'");
		}
		final String prefix=kmerResource.substring(0, kmerResource.length()-suffix.length());
		return ProkObject.loadLongKmersByType(kLong, prefix);
	}

	/** Loads a sweep's explicit kmer FASTA, or the normal shipped resource when unset. */
	private static LongHashSet loadEffectiveNcrnaKmerSet(String family, String shippedResource,
			String familyOverride, int kLong){
		final String path=(ncrnaSweepTarget(family) && NCRNA_KMERS_OVERRIDE!=null
			? NCRNA_KMERS_OVERRIDE : familyOverride);
		if(path==null){return loadNcrnaKmerSet(shippedResource, kLong);}
		java.io.File f=new java.io.File(path);
		if(!f.isFile()){
			throw new IllegalArgumentException("ncRNA kmer file not found for "+family+": "+path);
		}
		LongHashSet set=ProkObject.loadLongKmers(path, kLong);
		System.err.println("Loaded "+set.size()+" "+family+" "+kLong+"-mers from explicit sweep file "+path);
		return set;
	}

	static String parseNcrnaFamily(String value){
		if(value==null){throw new IllegalArgumentException("ncrnafamily requires a value");}
		String s=value.toLowerCase().replace('-', '_');
		if(s.equals("rnase_p")){s="rnasep";}
		else if(s.equals("srpsmall")){s="srp_small";}
		else if(s.equals("srplarge")){s="srp_large";}
		else if(s.equals("tm_rna") || s.equals("ssra")){s="tmrna";}
		if(!s.equals("rnasep") && !s.equals("srp_small") && !s.equals("srp_large") && !s.equals("tmrna")
				&& !s.equals("sixs_rf00013") && !s.equals("sixs_rf01685")){
			throw new IllegalArgumentException("ncrnafamily must be rnasep, srp_small, srp_large, tmrna, sixs_rf00013, or sixs_rf01685: "+value);
		}
		return s;
	}

	static float defaultNcrnaIdPass(String family){
		if(family.equals("rnasep")){return 0.65f;}
		if(family.equals("srp_small") || family.equals("srp_large")){return 0.70f;}
		if(family.equals("tmrna")){return 0.62f;}
		if(family.equals("sixs_rf00013") || family.equals("sixs_rf01685")){return 0.70f;}
		throw new IllegalArgumentException("No default idpass for ncRNA family: "+family);
	}

	static float defaultNcrnaIdBorderline(String family){
		if(family.equals("rnasep")){return 0.58f;}
		if(family.equals("srp_small") || family.equals("srp_large")){return 0.64f;}
		if(family.equals("tmrna")){return 0.60f;}
		if(family.equals("sixs_rf00013")){return 0.60f;}
		if(family.equals("sixs_rf01685")){return 0.70f;}
		throw new IllegalArgumentException("No default idborderline for ncRNA family: "+family);
	}

	static float parseFiniteFloat(String flagName, String value){
		final float x=Float.parseFloat(value);
		if(!Float.isFinite(x)){throw new IllegalArgumentException(flagName+" must be finite: "+value);}
		return x;
	}

	static float parseUnitFloat(String flagName, String value){
		final float x=parseFiniteFloat(flagName, value);
		if(x<0 || x>1){throw new IllegalArgumentException(flagName+" must be in [0,1]: "+value);}
		return x;
	}

	static boolean ncrnaSweepTarget(String name){
		return NCRNA_FAMILY_FILTER!=null && NCRNA_FAMILY_FILTER.equals(name);
	}

	/** Rejects ambiguous global overrides and logically inverted identity thresholds. */
	static void validateNcrnaSweepOverrides(){
		final boolean anyGenericOverride=NCRNA_KMERS_OVERRIDE!=null || NCRNA_WINDOW_PAD_OVERRIDE>=0
			|| !Float.isNaN(NCRNA_ID_PASS_OVERRIDE) || !Float.isNaN(NCRNA_ID_BORDERLINE_OVERRIDE)
			|| !Float.isNaN(NCRNA_HBM_PASS_OVERRIDE) || !Float.isNaN(NCRNA_SCORE_A_OVERRIDE)
			|| !Float.isNaN(NCRNA_SCORE_B_OVERRIDE) || !Float.isNaN(NCRNA_COLLAPSE_FRAC_OVERRIDE);
		final boolean anyFamilyKmers=RNASEP_KMERS_OVERRIDE!=null || SRPSMALL_KMERS_OVERRIDE!=null
			|| SRPLARGE_KMERS_OVERRIDE!=null || TMRNA_KMERS_OVERRIDE!=null
			|| SIXS_RF00013_KMERS_OVERRIDE!=null || SIXS_RF01685_KMERS_OVERRIDE!=null;
		final boolean anyFamilyPad=RNASEP_PAD_OVERRIDE>=0 || SRPSMALL_PAD_OVERRIDE>=0 || SRPLARGE_PAD_OVERRIDE>=0;
		final boolean anyTmrnaModels=TMRNA_CONSENSUS_OVERRIDE!=null || TMRNA_MODELS_OVERRIDE!=null;
		final boolean anyTmrnaOverride=TMRNA_KMERS_OVERRIDE!=null || anyTmrnaModels || TMRNA_PAD_OVERRIDE>=0
			|| !Float.isNaN(TMRNA_SCORE_A_OVERRIDE) || !Float.isNaN(TMRNA_SCORE_B_OVERRIDE);
		final boolean anySixsOverride=SIXS_RF00013_KMERS_OVERRIDE!=null || SIXS_RF01685_KMERS_OVERRIDE!=null
			|| SIXS_RF00013_PAD_OVERRIDE>=0 || SIXS_RF01685_PAD_OVERRIDE>=0;
		final boolean anyFamilyScores=!Float.isNaN(RNASEP_SCORE_A_OVERRIDE)
			|| !Float.isNaN(RNASEP_SCORE_B_OVERRIDE) || !Float.isNaN(SRPSMALL_SCORE_A_OVERRIDE)
			|| !Float.isNaN(SRPSMALL_SCORE_B_OVERRIDE) || !Float.isNaN(SRPLARGE_SCORE_A_OVERRIDE)
			|| !Float.isNaN(SRPLARGE_SCORE_B_OVERRIDE) || !Float.isNaN(TMRNA_SCORE_A_OVERRIDE)
			|| !Float.isNaN(TMRNA_SCORE_B_OVERRIDE);
		if(NCRNA_FAMILY_FILTER!=null && !NCRNA_FAMILIES_ENABLED){
			throw new IllegalArgumentException("ncrnafamily requires ncrna=t");
		}
		if(NCRNA_FAMILY_FILTER!=null && NCRNA_FAMILY_FILTER.equals("tmrna") && !TMRNA_ENABLED){
			throw new IllegalArgumentException("ncrnafamily=tmrna requires tmrna=t");
		}
		if(NCRNA_FAMILY_FILTER!=null && NCRNA_FAMILY_FILTER.startsWith("sixs_") && !SIXS_ENABLED){
			throw new IllegalArgumentException("a sixS ncrnafamily requires sixs=t");
		}
		if(anyFamilyKmers && !NCRNA_FAMILIES_ENABLED){
			throw new IllegalArgumentException("family ncRNA kmer overrides require ncrna=t");
		}
		if(anyFamilyPad && !NCRNA_FAMILIES_ENABLED){
			throw new IllegalArgumentException("family ncRNA pad overrides require ncrna=t");
		}
		if(anyFamilyScores && !NCRNA_FAMILIES_ENABLED){
			throw new IllegalArgumentException("family ncRNA score overrides require ncrna=t");
		}
		if(anyTmrnaOverride && (!NCRNA_FAMILIES_ENABLED || !TMRNA_ENABLED)){
			throw new IllegalArgumentException("tmRNA-specific overrides require ncrna=t tmrna=t");
		}
		if(anySixsOverride && (!NCRNA_FAMILIES_ENABLED || !SIXS_ENABLED)){
			throw new IllegalArgumentException("sixS-specific overrides require ncrna=t sixs=t");
		}
		if(anyGenericOverride && NCRNA_FAMILY_FILTER==null){
			throw new IllegalArgumentException("ncRNA sweep overrides require exactly one ncrnafamily=");
		}
		if(NCRNA_FAMILY_FILTER!=null){
			final float defaultPass=defaultNcrnaIdPass(NCRNA_FAMILY_FILTER);
			final float defaultBorder=defaultNcrnaIdBorderline(NCRNA_FAMILY_FILTER);
			final float pass=resolveSweepFloat(NCRNA_FAMILY_FILTER, NCRNA_ID_PASS_OVERRIDE, defaultPass);
			final float borderline=resolveSweepFloat(NCRNA_FAMILY_FILTER, NCRNA_ID_BORDERLINE_OVERRIDE, defaultBorder);
			if(borderline>pass){
				throw new IllegalArgumentException("ncrnaidborderline ("+borderline+") must be <= ncrnaidpass ("+pass+")");
			}
		}
	}

	/** C3 addition (G11, 2026-08-28, Citan's per-index name-verification requirement): the
	 * first whitespace-delimited token of a header/name string -- see the caller's comment
	 * for why this (not raw full-string equality) is the robust comparison. Empty string
	 * (not null) for a null input, so callers can uniformly treat "" as the empty-token
	 * case to reject, without a separate null check. */
	private static String firstToken(String s){
		if(s==null){return "";}
		final int i=s.indexOf(' ');
		final int j=s.indexOf('\t');
		final int cut=(i<0 ? j : (j<0 ? i : Math.min(i, j)));
		return cut<0 ? s : s.substring(0, cut);
	}

	/** C3 addition (G11, 2026-08-28, Citan's uniqueness requirement): the model-name array
	 * from a BaseGraph[] (for uniformly reusing findDuplicateOrEmptyToken on either side of
	 * the fasta-vs-hbm comparison). Individual null entries map to null (firstToken treats
	 * null as ""), so a null BaseGraph in the array surfaces as an empty-token rejection
	 * rather than an NPE. */
	private static String[] modelNamesOf(consensus.BaseGraph[] models){
		String[] out=new String[models.length];
		for(int i=0; i<models.length; i++){out[i]=(models[i]==null ? null : models[i].name);}
		return out;
	}

	/** C3 addition (G11, 2026-08-28, Citan's uniqueness requirement): scans `names` (already
	 * length-matched to the library by the caller), normalizing each to firstToken, and
	 * returns a diagnostic string identifying the FIRST empty-token or duplicate-token
	 * violation found, or null if all normalized tokens are non-empty and unique. `side` is
	 * "fasta" or "hbm", folded into the message so a caller doesn't have to. Checked before
	 * any per-index equality comparison -- a duplicate on either side means first-token
	 * matching alone cannot distinguish a correct pairing from a reorder hiding behind it. */
	private static String findDuplicateOrEmptyToken(String side, String[] names){
		java.util.HashMap<String,Integer> seen=new java.util.HashMap<>();
		for(int i=0; i<names.length; i++){
			final String tok=firstToken(names[i]);
			if(tok.isEmpty()){
				return side+" index "+i+" has an empty/null normalized name (raw='"+names[i]+"')";
			}
			final Integer prior=seen.put(tok, i);
			if(prior!=null){
				return side+" indices "+prior+" and "+i+" both normalize to '"+tok
					+"' (raw='"+names[prior]+"' vs '"+names[i]+"') -- normalized names must be unique "
					+"for per-index alignment to be verifiable.";
			}
		}
		return null;
	}

	/** C3 addition (G11, 2026-08-28): family-mean length for the boundary NN's lengthRatio
	 * feature, matching NcrnaBoundaryVectorGen's own medianLength helper exactly (same
	 * formula, duplicated rather than shared since that method is private in a different
	 * class and this is a 4-line accessor, not worth widening NcrnaBoundaryVectorGen for). */
	private static float medianLength(byte[][] library){
		int[] lens=new int[library.length];
		for(int i=0; i<library.length; i++){lens[i]=library[i].length;}
		java.util.Arrays.sort(lens);
		return lens[lens.length/2];
	}

	/** Maximum number of reads to process, or -1 for no limit */
	private long maxReads=-1;
	/** Whether to attempt merging of overlapping paired reads */
	private boolean merge;
	/** Whether to perform error correction on overlapping paired reads */
	private boolean ecco;
	/** Number of iterative passes for model refinement */
	private int passes=1;
	
	/** Total number of input reads processed */
	private long readsIn=0;
	/** Total number of input bases processed */
	private long basesIn=0;
	/** Total number of genes called and output */
	private long genesOut=0;
	/** Total number of bytes written to output files */
	private long bytesOut=0;
	
	/** Minimum gene length in bases for calling genes */
	private static int minLen=80;//Actually a much higher value like 200 seems optimal compared to NCBI
	/** Maximum allowed overlap between genes on same strand */
	private static int maxOverlapSameStrand=80;
	/** Maximum allowed overlap between genes on opposite strands */
	private static int maxOverlapOppositeStrand=110;
	
	/* for kinnercds=6 */
//	private static float minStartScore=-0.10f;
//	private static float minStopScore=-0.5f;//Not useful; disabled
//	private static float minKmerScore=0.04f;//Does not seem useful.
//	private static float minOrfScore=40f; //Higher increases SNR dramatically but reduces TP rate
//	private static float minAvgScore=0.08f; //Not very effective

	/* for kinnercds=7 */
	/** Minimum score threshold for gene start sites */
	private static float minStartScore=-0.10f;
	/** Minimum score threshold for gene stop sites */
	private static float minStopScore=-0.5f;//Not useful; disabled
	/** Minimum score threshold for internal kmer content */
	private static float minKmerScore=0.02f;//Does not seem useful.
	/** Minimum total score threshold for complete ORF */
	private static float minOrfScore=50f; //Higher increases SNR dramatically but reduces TP rate
	/** Minimum average score per base for ORF retention */
	private static float minAvgScore=0.08f; //Not very effective
	
	/** Count of gene stop sites identified during processing */
	long geneStopsMade=0;
	/** Count of gene start sites identified during processing */
	long geneStartsMade=0;
	/** Count of gene start sites retained after filtering */
	long geneStartsRetained=0;
	/** Count of gene stop sites retained after filtering */
	long geneStopsRetained=0;
	/** Count of final gene start sites output */
	long geneStartsOut=0;

	/** Count of tRNA genes identified and output */
	long tRNAOut=0;
	/** Count of generic ncRNA (RNase P, SRP, ...) genes identified and output. Forward-ported
	 * from Noire's tree (2026-08-28, C3 merge; Citan's correction) -- GeneCaller already
	 * increments its own per-instance rnaOut (GeneCaller.java, orf.type==RNA branch), pre-
	 * dating this merge, but CallGenes never aggregated or reported it -- ncRNA calls were
	 * silently invisible in the summary/JSON output even when Gate A was on. Aggregated below
	 * alongside tRNAOut; reported only when a family actually loaded (gated on
	 * !GeneCaller.ncrnaFamilies.isEmpty(), matching tRNAOut's own calltRNA gate). */
	long rnaOut=0;
	/** Completed generic-ncRNA alignments, in GeneCaller.ncrnaFamilies order. */
	long[] ncrnaAlignments;
	/** Count of 16S rRNA genes identified and output */
	long r16SOut=0;
	/** Count of 23S rRNA genes identified and output */
	long r23SOut=0;
	/** Count of 5S rRNA genes identified and output */
	long r5SOut=0;
	/** Count of 18S rRNA genes identified and output */
	long r18SOut=0;
	
	/** Score tracker for all CDS ORFs before filtering */
	ScoreTracker stCds=new ScoreTracker(CDS);
	/** Score tracker for CDS ORFs after initial filtering */
	ScoreTracker stCds2=new ScoreTracker(CDS);
	/** Score tracker for final CDS ORFs that passed all filters */
	ScoreTracker stCdsPass=new ScoreTracker(CDS);
	/** Score tracker for tRNA genes */
	ScoreTracker sttRNA=new ScoreTracker(tRNA);
	/** Score tracker for 16S rRNA genes */
	ScoreTracker st16s=new ScoreTracker(r16S);
	/** Score tracker for 23S rRNA genes */
	ScoreTracker st23s=new ScoreTracker(r23S);
	/** Score tracker for 5S rRNA genes */
	ScoreTracker st5s=new ScoreTracker(r5S);
	/** Score tracker for 18S rRNA genes */
	ScoreTracker st18s=new ScoreTracker(r18S);
	
	/** Array of score trackers for all gene types for easy iteration */
	ScoreTracker[] trackers=new ScoreTracker[] {stCdsPass, sttRNA, st16s, st23s, st5s, st18s};
	
	/** Number of bins for gene length histogram */
	private int geneHistBins=1000;
	/** Divisor for gene length histogram binning (bases per bin) */
	private int geneHistDiv=21;
	/** Histogram array for tracking gene length distribution */
	private final long[] geneHist;
	/** Whether to include zero-count bins in histogram output */
	private boolean printZeroCountHist=false;
	
	/*--------------------------------------------------------------*/

	/** List of input FASTA filenames */
	private ArrayList<String> fnaList=new ArrayList<String>();
	/** List of probabilistic gene model filenames */
	private ArrayList<String> pgmList=new ArrayList<String>();
	/** List of input GFF annotation filenames */
	private ArrayList<String> inGffList=new ArrayList<String>();
	/** Output filename for GFF3 gene annotations */
	private String outGff=null;
	/** Output filename for amino acid translations */
	private String outAmino=null;
	/** Output filename for 16S rRNA sequences */
	private String out16S=null;
	/** Output filename for 18S rRNA sequences */
	private String out18S=null;
	/** GFF filename for comparison/validation of gene calling results */
	private String compareToGff=null;
	private String orfNetPath=null;
	private String orfNetLowPath=null;
	private String orfNetMidPath=null;
	private String orfNetHighPath=null;
	private String[] orfNetPaths=null;
	private String[] stopNetPaths=null;
	private String metaNetPath=null;
	private String[] metaNetPaths=null;
	private float[] gcMeans=null;
	//Default ON (Brian, 2026-08-16): classify each input via QuickClade and use the matching per-phylum PGM.
	//Intended for ISOLATES (one organism per file) -- a metagenome would be misclassified to a single phylum,
	//so set taxonomy=f for mixed samples.  Falls back to the general PGM (with a warning) if the clade server
	//is unreachable, so it never hard-fails offline.  Flag taxonomy=/tax= overrides.  Provisional pending re-sweep.
	private boolean useTaxonomy=true;
	/** Use the LOCAL QuickClade reference database (large multi-GB files) instead of the network
	 * server. Default false -- the server is the lightweight default, and local=t opts INTO the
	 * multi-GB local files if present (falls back to the server, with a warning, if they're not).
	 * This is the REVERSE of QuickClade's own CLI default (clade.CladeSearcher favors local,
	 * falling back to the server only when local deps are missing) -- CallGenes is meant to stay
	 * lightweight by default; the heavy local path is opt-in. (Brian, 2026-08-29) */
	private boolean useLocalClade=false;
	/** Domain detected by the most recent classifyPhylum() call, or null if none/unclassified */
	private String lastDetectedDomain=null;
	/** Phylum detected by the most recent classifyPhylum() call, or null if none/unclassified */
	private String lastDetectedPhylum=null;
	/** "local" or "server[:address]" -- which QuickClade path served the most recent classification,
	 * or null if classification was never attempted (taxonomy=f or percontig=t) */
	private String lastClassificationSource=null;
	private boolean perContig=false;
	private String taxAddress="refseq";
	private boolean trnaAlign=true;
	private boolean trnaLibExplicit=false;
	/** Boundary-precision NN paths (trnaboundarynet=/trnaboundarystarttable=/
	 * trnaboundarystoptable=). SHIPPED, DEFAULT ON (Brian, Aug 22 2026, via Noire): left null
	 * here means "not explicitly overridden," NOT "off" -- the post-parse resolution block in
	 * parse() fills these in from the shipped resources (?trna_boundary_net.bbnet etc.) unless
	 * trnaBoundaryDisabled is set. See prok.TrnaCaller.loadBoundaryNet. */
	private String trnaBoundaryNetPath=null, trnaBoundaryStartTablePath=null, trnaBoundaryStopTablePath=null;
	/** Opt-out switch (Brian, Aug 22 2026, via Noire): trnaboundarynet=f/false sets this,
	 * skipping the shipped-default resolution entirely -- the only way to get the pre-Aug-22
	 * no-refiner behavior back. */
	private boolean trnaBoundaryDisabled=false;
	/** Two-net dispatch (Noire's spec, Aug 22 2026, queued directive #2/#3 infrastructure):
	 * trnaboundary5net=/trnaboundary3net= override trnaBoundaryNetPath for just the 5' or 3'
	 * dispatch; whichever of the two is left null falls back to trnaBoundaryNetPath (backward
	 * compatible with the single-shared-net configuration). See the resolution in parse(). */
	private String trnaBoundary5NetPath=null, trnaBoundary3NetPath=null;

	/** Gate A (Noire's ncRNA-family-loading work, forward-ported 2026-08-28, C3 merge):
	 * generic ncRNA family loading (RNase P, SRP small/large) is OFF by default (Citan's
	 * explicit requirement -- Brian shipped ncRNAs off; this is NOT the same shipped-default-
	 * on treatment tRNA's boundary NN got above). ncrna=t / generalncrna=t turns it on.
	 * When false, loadNcrnaResources() is never called at all -- zero startup cost, and
	 * GeneCaller.ncrnaFamilies stays empty, so every ncRNA-gated code path (the
	 * !ncrnaFamilies.isEmpty() checks in GeneCaller) is a no-op, byte-identical to before
	 * this merge. The loader method itself stays public/callable directly (e.g. by
	 * NcrnaBoundaryVectorGen's training-vector generator, or eval scripts) -- this flag only
	 * gates PRODUCTION callgenes invocations. */
	static boolean NCRNA_FAMILIES_ENABLED=false;
	/** Experimental tmRNA bundle gate. Subordinate to ncrna=t and deliberately independent
	 * of the existing three-family gate so adding unfinished tmRNA resources cannot change an
	 * established ncrna=t run. */
	static boolean TMRNA_ENABLED=false;
	/** Experimental 6S/ssrS bundle gate.  Both RF00013 and RF01685 are loaded together,
	 * subordinate to ncrna=t, and OFF by default pending production-competition validation. */
	static boolean SIXS_ENABLED=false;
	/** Gate B (C3, Noire's spec plans/c3_ncrnaboundaryscorer_spec.md; G11, 2026-08-28):
	 * ncRNA boundary-precision-NN refinement, subordinate to Gate A -- ncrna=t alone gives
	 * plain scavenger calling with no boundary refinement; ncrnaboundarynet=t additionally
	 * attempts to lazy-load each family's boundary net/tables/model (see addNcrnaFamily).
	 * OFF by default, matching Gate A. REQUIRES Gate A: ncrnaboundarynet=t with ncrna=f is
	 * REJECTED LOUDLY at parse time (Citan's correction, 2026-08-28), not a silent no-op --
	 * see the assert in parse(). */
	static boolean NCRNA_BOUNDARY_NN_ENABLED=false;

	/** P1 (Brian via Citan, 2026-08-28): per-family scavenger candidate-window padding
	 * overrides. -1 means unset -- loadNcrnaResources() uses each family's shipped default
	 * (rnasep 320 / srp_small 75 / srp_large 250 / tmrna 370) when the corresponding flag is absent.
	 * A non-negative value passed via
	 * rnaseppad=/srpsmallpad=/srplargepad=/tmrnapad= (validated >=0 in parse()) overrides
	 * ONLY that family's windowPad argument to addNcrnaFamily; the other families are untouched.
	 * Deliberately narrow: does not touch NcrnaScavenger's window-construction logic itself,
	 * only the constant CallGenes has always passed into it. */
	static int RNASEP_PAD_OVERRIDE=-1;
	static int SRPSMALL_PAD_OVERRIDE=-1;
	static int SRPLARGE_PAD_OVERRIDE=-1;
	static int TMRNA_PAD_OVERRIDE=-1;
	static int SIXS_RF00013_PAD_OVERRIDE=-1;
	static int SIXS_RF01685_PAD_OVERRIDE=-1;

	/** Engineering-sweep controls. The generic overrides deliberately require one explicit
	 * family so a command cannot silently apply one value to biologically different callers. */
	static String NCRNA_FAMILY_FILTER=null;
	static String NCRNA_KMERS_OVERRIDE=null;
	static String RNASEP_KMERS_OVERRIDE=null;
	static String SRPSMALL_KMERS_OVERRIDE=null;
	static String SRPLARGE_KMERS_OVERRIDE=null;
	static String TMRNA_KMERS_OVERRIDE=null;
	static String SIXS_RF00013_KMERS_OVERRIDE=null;
	static String SIXS_RF01685_KMERS_OVERRIDE=null;
	static String TMRNA_CONSENSUS_OVERRIDE=null;
	static String TMRNA_MODELS_OVERRIDE=null;
	static int NCRNA_WINDOW_PAD_OVERRIDE=-1;
	static float NCRNA_ID_PASS_OVERRIDE=Float.NaN;
	static float NCRNA_ID_BORDERLINE_OVERRIDE=Float.NaN;
	static float NCRNA_HBM_PASS_OVERRIDE=Float.NaN;
	static float NCRNA_SCORE_A_OVERRIDE=Float.NaN;
	static float NCRNA_SCORE_B_OVERRIDE=Float.NaN;
	static float RNASEP_SCORE_A_OVERRIDE=Float.NaN;
	static float RNASEP_SCORE_B_OVERRIDE=Float.NaN;
	static float SRPSMALL_SCORE_A_OVERRIDE=Float.NaN;
	static float SRPSMALL_SCORE_B_OVERRIDE=Float.NaN;
	static float SRPLARGE_SCORE_A_OVERRIDE=Float.NaN;
	static float SRPLARGE_SCORE_B_OVERRIDE=Float.NaN;
	static float TMRNA_SCORE_A_OVERRIDE=Float.NaN;
	static float TMRNA_SCORE_B_OVERRIDE=Float.NaN;
	static float NCRNA_COLLAPSE_FRAC_OVERRIDE=Float.NaN;

	/** Gate A/B combination validator (Citan's correction, 2026-08-28), factored out of
	 * parse() into its own package-visible, directly testable method: ncrnaboundarynet=t
	 * with ncrna=f is a nonsensical combination (nothing to attach boundary refinement to),
	 * and must fail loud as invalid input rather than silently no-op. Reads the two static
	 * gate fields directly (not parameterized) since parse() has no other state to pass and
	 * a test can simply set the fields before calling. */
	static void validateNcrnaGateCombo(){
		assert(!NCRNA_BOUNDARY_NN_ENABLED || NCRNA_FAMILIES_ENABLED) : "ncrnaboundarynet=t requires "
			+"ncrna=t (or generalncrna=t) -- boundary refinement has no families to attach to "
			+"when generic ncRNA calling itself is off.";
		assert(!TMRNA_ENABLED || NCRNA_FAMILIES_ENABLED) : "tmrna=t requires ncrna=t (or generalncrna=t)";
		assert(!SIXS_ENABLED || NCRNA_FAMILIES_ENABLED) : "sixs=t requires ncrna=t (or generalncrna=t)";
	}

	/** Output filename for statistics summary */
	private String outStats="stderr";
	/** Output filename for gene length histogram */
	private String geneHistFile=null;
	/** Whether to output statistics in JSON format */
	private boolean json_out=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** FileFormat object for GFF output configuration */
	private final FileFormat ffoutGff;
	/** FileFormat object for amino acid output configuration */
	private final FileFormat ffoutAmino;
	/** FileFormat object for 16S rRNA output configuration */
	private final FileFormat ffout16S;
	/** FileFormat object for 18S rRNA output configuration */
	private final FileFormat ffout18S;
	
	/** Determines how sequence is processed if it will be output */
	int mode=TRANSLATE;
	
	/** Translate nucleotides to amino acids */
	private static final int TRANSLATE=1;
	/** Translate nucleotides to amino acids,
	 * then translate back to nucleotides */
	private static final int RETRANSLATE=2;
	/** Re-encode coding regions of nucleotide
	 * sequences as a canonical codons */
	private static final int RECODE=3;
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for progress messages and logging */
	private PrintStream outstream=System.err;
	/** Whether to print verbose progress and debugging information */
	public boolean verbose=false;
	/** Whether to include extended statistics in output */
	public boolean extendedStats=false;
	/** Flag indicating whether processing encountered errors */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	/**
	 * Whether to maintain input order in output (may cause hanging on some sequences)
	 */
	//TODO: Possible bug [prok/CallGenes#002] (MEDIUM; gated on ordered=true) - Brian-documented hang, mechanism localized: the GFF path (processList) calls bsw.add(bb, ln.id) for EVERY ln.id (even empty) -> safe. But the SEQUENCE streams ros16S/ros18S/rosAmino only add(..., r.numericID) when fetchType/translate return non-null/non-empty (processRead L~1017/1026/1034). A read producing no gene of that type leaves a GAP in the per-numericID ordered stream -> the ordered writer blocks forever waiting for the missing id (RefSeq mito: no CDS -> rosAmino gaps). Fix (Brian's recipe + SplitRibo's always-emit pattern): emit for each numericID even an empty list. Multi-site contract change -> flag-for-Brian, NOT patched.
	private boolean ordered=false; //this is OK sometimes, but sometimes hangs (e.g. on RefSeq mito), possibly if a sequence produces nothing.
	//To fix it, just ensure functions like translate always produce an array, even if it is empty.
	//TODO: Possible bug [prok/CallGenes#002] - real MEDIUM hang (gated on ordered=true) when a sequence yields no genes;
	//fix per the note above (have translate/fetchType return a non-null empty array). crash-loudly-not-hang.
	
}
