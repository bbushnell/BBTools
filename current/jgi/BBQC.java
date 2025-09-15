package jgi;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.TimeZone;

import align2.BBMap;
import dna.Data;
import fileIO.ByteFile1;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Tools;
import shared.TrimRead;
import stream.FASTQ;

/**
 * Wrapper for BBDukF, BBMap, and BBNorm to perform quality-control and artifact removal.
 * @author Brian Bushnell
 * @date Jan 20, 2013
 *
 */
public class BBQC {

	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Methods    ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Program entrance from command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		ReadWrite.USE_PIGZ=true;
		ReadWrite.USE_UNPIGZ=true;
		
		//Create a filter instance
		BBQC filter=new BBQC(args);
		
		///...and execute it.
		filter.process();
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	BBQC(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		//Symbols to insert in output filename to denote operations performed; may be overriden from command line
		String symbols_=null;//"filtered"
		int passes_=-1;
		
		//Parse argument list
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("="); //Expect key=value pairs
			String a=split[0].toLowerCase(); //All keys are converted to lower case
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){
				if(a.equals("pigz")){
					pigz=b;
				}else if(a.equals("unpigz")){
					unpigz=b;
				}else if(a.equals("zl") || a.equals("ziplevel")){
					zl=b;
				}
			}else if(Parser.parseCommonStatic(arg, a, b)){
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
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfin1=b;
			}else if(a.equals("qfout") || a.equals("qfout1")){
				qfout1=b;
			}else if(a.equals("qfin2")){
				qfin2=b;
			}else if(a.equals("qfout2")){
				qfout2=b;
			}else if(a.equals("ref")){
				if(b!=null){
					if(!b.contains(",") || new File(b).exists()){
						filterrefs.add(b);
					}else{
						String[] split2=b.split(",");
						for(String s2 : split2){
							filterrefs.add(s2);
						}
					}
				}
			}else if(a.equals("artifactdb")){
				mainArtifactFile=b;
			}else if(a.equals("rnadb")){
				artifactFileRna=b;
			}else if(a.equals("dnadb")){
				artifactFileDna=b;
			}else if(a.equals("phixref")){
				phixRef=b;
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("ml") || a.equals("minlen") || a.equals("minlength")){
				minLen=Integer.parseInt(b);
			}else if(a.equals("mlf") || a.equals("minlenfrac") || a.equals("minlenfraction") || a.equals("minlengthfraction")){
				minLenFraction=Float.parseFloat(b);
			}else if(a.equals("path") || a.equals("outdir")){
				outDir=b;
			}else if(a.equals("symbols")){
				symbols_=b;
			}else if(a.equals("overallstats") || a.equals("stats")){
				rqcStatsName=b;
			}else if(a.equals("scafstats")){
				scaffoldStatsName=b;
			}else if(a.equals("kmerstats")){
				kmerStatsName=b;
			}else if(a.equals("log")){
				logName=b;
			}else if(a.equals("filelist")){
				fileListName=b;
			}else if(a.equals("compress")){
				compress=Parse.parseBoolean(b);
			}else if(a.equals("rna")){
				rnaFlag=Parse.parseBoolean(b);
			}else if(a.equals("phix")){
				phixFlag=Parse.parseBoolean(b);
			}else if(a.equals("lambda")){
				lambdaFlag=Parse.parseBoolean(b);
			}else if(a.equals("pjet")){
				pjetFlag=Parse.parseBoolean(b);
			}else if(a.equals("ktrim")){
				ktrim=b;
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
			}else if(a.equals("trimhdist")){
				hdist_trim=Integer.parseInt(b);
			}else if(a.equals("trimhdist2")){
				hdist2_trim=Integer.parseInt(b);
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
			}else if(a.equals("human") || a.equals("removehuman")){
				removehuman=Parse.parseBoolean(b);
			}else if(a.equals("normalize") || a.equals("norm")){
				normalize=Parse.parseBoolean(b);
			}else if(a.equals("ecc")){
				ecc=Parse.parseBoolean(b);
			}else if(a.equals("aec") || a.equals("aecc")){
				aecc=Parse.parseBoolean(b);
				if(aecc){ecc=true;}
			}else if(a.equals("cecc")){
				cecc=Parse.parseBoolean(b);
				if(cecc){ecc=true;}
			}else if(a.equals("markerrorsonly") || a.equals("meo")){
				meo=Parse.parseBoolean(b);
			}else if(a.equals("tam")){
				tam=Parse.parseBoolean(b);
			}else if(a.equals("taf")){
				trimAfterFiltering=Parse.parseBoolean(b);
			}else if(a.equals("mue")){
				mue=Parse.parseBoolean(b);
			}else if(a.equals("mw1")){
				mw1=Parse.parseBoolean(b);
			}else if(a.equals("max") || a.equals("maxdepth")){
				maxdepth=Integer.parseInt(b);
			}else if(a.equals("min") || a.equals("mindepth")){
				mindepth=Integer.parseInt(b);
			}else if(a.equals("target") || a.equals("targetdepth")){
				target=Integer.parseInt(b);
			}else if(a.equals("prehashes")){
				prehashes=Integer.parseInt(b);
			}else if(a.equals("passes")){
				passes_=Integer.parseInt(b);
			}else if(a.equals("hashes")){
				hashes=Integer.parseInt(b);
			}else if(a.equals("bits")){
				bits=Integer.parseInt(b);
			}else if(a.equals("minratio")){
				minratio=Float.parseFloat(b);
			}else if(a.equals("maxindel")){
				maxindel=Integer.parseInt(b);
			}else if(a.equals("kfilter")){
				kfilter=Integer.parseInt(b);
			}else if(a.equals("hits") || a.equals("minhits")){
				minhits=Integer.parseInt(b);
			}else if(a.equals("fast")){
				fast=Parse.parseBoolean(b);
			}else if(a.equals("local")){
				local=Parse.parseBoolean(b);
			}else if(a.equals("mappath") || a.equals("indexpath")){
				indexPath=b;
			}else if(a.equals("mapref")){
				mapRef=b;
			}else if(a.equals("copyundefined") || a.equals("cu")){
				copyUndefined=Parse.parseBoolean(b);
			}else if(a.equals("qtrim")){
				if(b==null){qtrim="rl";}
				else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){qtrim="l";}
				else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){qtrim="r";}
				else if(b.equalsIgnoreCase("both") || b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr")){qtrim="lr";}
				else if(Tools.isDigit(b.charAt(0))){
					trimq=Byte.parseByte(b);
					qtrim=(trimq>=0 ? "lr" : "f");
				}else{qtrim=""+Parse.parseBoolean(b);}
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
			}else if(in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					in1=arg.replace("#", "1");
					in2=arg.replace("#", "2");
				}
			}else{
				//Uncaptured arguments are passed to BBDuk
				primaryArgList.add(arg);
			}
		}
		
		if(passes_>0){passes=passes_;}
		else if(!normalize){passes=1;}
		
		if(hdist2_trim<0){hdist2_trim=hdist_trim;}
		
		//Set final field 'symbols'
		symbols=(symbols_==null ? abbreviation() : symbols_);
		
		//Pass overwrite flag to BBDuk
		primaryArgList.add("ow="+overwrite);
		
		if(outDir!=null){
			outDir=outDir.trim().replace('\\', '/');
			if(outDir.length()>0 && !outDir.endsWith("/")){outDir=outDir+"/";}
		}else{outDir="";}
		
		{//Prepend output directory to output files
			if(logName!=null){logName=outDir+logName+".tmp";} //Add '.tmp' to log file
			if(fileListName!=null){fileListName=outDir+fileListName;}
		}

		{//Create unique output file names for second pass
			if(rqcStatsName!=null){
				rqcStatsName_kt=outDir+"ktrim_"+rqcStatsName;
				rqcStatsName=outDir+rqcStatsName;
			}
			if(kmerStatsName!=null){
				kmerStatsName_kt=outDir+"ktrim_"+kmerStatsName;
				kmerStatsName=outDir+kmerStatsName;
			}
			if(scaffoldStatsName!=null){
				scaffoldStatsName_kt=outDir+"ktrim_"+scaffoldStatsName;
				scaffoldStatsName=outDir+scaffoldStatsName;
			}
		}
		
		//Create output filename from input filename if no output filename is specified
		if(out1==null && in1!=null){
			File f=new File(in1);
			String name=f.getName();
			String raw=ReadWrite.rawName(name);
			int x=raw.lastIndexOf('.');
			if(x>-1){
				out1=raw.substring(0, x)+"."+symbols+raw.substring(x)+(compress ? ".gz" : "");
			}else{
				out1=raw+"."+symbols+".fastq"+(compress ? ".gz" : "");
			}
		}
		
		tempSalt=KmerNormalize.getSalt(out1, 0);
	}

	
	/*--------------------------------------------------------------*/
	/*----------------      Processing Methods      ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Primary method to fully execute the program.
	 */
	public void process(){
		
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
			log("start", false);
		}
		
		//Create file list file
		if(fileListName!=null){
			boolean b=Tools.canWrite(fileListName, overwrite);
			assert(b) : "Can't write to "+fileListName;
			
			StringBuilder sb=new StringBuilder();
			if(out1!=null){sb.append("filtered_fastq="+out1).append('\n');}
			if(qfout1!=null){sb.append("filtered_qual="+qfout1).append('\n');}
			if(out2!=null){sb.append("filtered_fastq_2="+out2).append('\n');}
			if(qfout2!=null){sb.append("filtered_qual_2="+qfout2).append('\n');}
			if(ihistName!=null){sb.append("ihist="+ihistName).append('\n');}
			if(scaffoldStatsName!=null){sb.append("scafstats="+scaffoldStatsName).append('\n');}
			
			if(sb.length()>0){
				ReadWrite.writeString(sb, fileListName, false);
			}
		}

		final String trimPrefix="TEMP_TRIM_"+tempSalt+"_";
		final String humanPrefix="TEMP_HUMAN_"+tempSalt+"_";
		final String filterPrefix="TEMP_FILTER_"+tempSalt+"_";
		
		int oldZL=ReadWrite.ZIPLEVEL;
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.ALLOW_ZIPLEVEL_CHANGE=false;

		final String in1s=stripDirs(in1), in2s=stripDirs(in2), qfin1s=stripDirs(qfin1), qfin2s=stripDirs(qfin2);
		final String out1s=stripDirs(out1), out2s=stripDirs(out2), qfout1s=stripDirs(qfout1), qfout2s=stripDirs(qfout2);
		
		ktrim(in1, in2, out1s, out2s, qfin1, qfin2, qfout1s, qfout2s, trimPrefix);
		
		if(in2!=null && out2==null){
			FASTQ.FORCE_INTERLEAVED=true;
			FASTQ.TEST_INTERLEAVED=false;
		}
		
		
		if(removehuman){
			filter(out1s, out2s, out1s, out2s, qfout1s, qfout2s, qfout1s, qfout2s, trimPrefix, filterPrefix, true, true, false);
			delete(trimPrefix,  out1s, out2s, qfout1s, qfout2s);
			if(normalize || ecc){
				dehumanize(out1s, out2s, out1s, out2s, qfout1s, qfout2s, filterPrefix, humanPrefix, true, true, false);
				delete(filterPrefix, out1s, out2s, qfout1s, qfout2s);
				Data.unloadAll();
				ReadWrite.ZIPLEVEL=oldZL;
				ReadWrite.ALLOW_ZIPLEVEL_CHANGE=true;
				normalize(out1s, out2s, out1, out2, qfout1s, qfout2s, qfout1, qfout2, humanPrefix, "", true, true, true);
				delete(humanPrefix, out1s, out2s, qfout1s, qfout2s);
			}else{
				ReadWrite.ZIPLEVEL=oldZL;
				ReadWrite.ALLOW_ZIPLEVEL_CHANGE=true;
				dehumanize(out1s, out2s, out1, out2, qfout1s, qfout2s, filterPrefix, "", true, true, true);
				delete(filterPrefix, out1s, out2s, qfout1s, qfout2s);
				Data.unloadAll();
			}
		}else{
			if(normalize || ecc){
				filter(out1s, out2s, out1s, out2s, qfout1s, qfout2s, qfout1s, qfout2s, trimPrefix, filterPrefix, true, true, false);
				delete(trimPrefix, out1s, out2s, qfout1s, qfout2s);
				normalize(out1s, out2s, out1, out2, qfout1s, qfout2s, qfout1, qfout2, filterPrefix, "", true, true, true);
				delete(filterPrefix, out1s, out2s, qfout1s, qfout2s);
			}else{
				filter(out1s, out2s, out1, out2, qfout1s, qfout2s, qfout1, qfout2, trimPrefix, "", true, true, true);
				delete(trimPrefix, out1s, out2s, qfout1s, qfout2s);
			}
		}
		
		//Write combined stats file (number of reads/bases present/removed in each stage)
		if(rqcStatsName!=null){
			final TextStreamWriter tsw=new TextStreamWriter(rqcStatsName, overwrite, false, false);
			tsw.start();
			tsw.println(BBDuk.rqcString());
			tsw.poisonAndWait();
		}
		
		//Finish writing log file
		if(logName!=null){
			log("complete", true);
			if(logName.endsWith(".tmp")){ //Remove .tmp extension
				String old=logName;
				logName=logName.substring(0, logName.length()-4);
				try {
					new File(old).renameTo(new File(logName));
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
	}
	
	
	/**
	 * Runs BBDuk to perform:
	 * Kmer trimming, short read removal.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param qfin1 Primary input qual file
	 * @param qfin2 Secondary input qual file
	 * @param qfout1 Primary output qual file
	 * @param qfout2 Secondary output qual file
	 * @param outPrefix Append this prefix to output filenames
	 */
	private void ktrim(String in1, String in2, String out1, String out2, String qfin1, String qfin2, String qfout1, String qfout2, String outPrefix){
		
		log("ktrim start", true);
		System.err.println("\nAdapter Trimming Phase Start");
		
		ArrayList<String> argList=new ArrayList<String>();
		
		{//Fill list with BBDuk arguments
			argList.add("ktrim="+(ktrim==null ? "f" : ktrim));
			if(minLen>0){argList.add("minlen="+minLen);}
			if(minLenFraction>0){argList.add("minlenfraction="+minLenFraction);}
			argList.add("mink="+mink);
			if("r".equalsIgnoreCase(ktrim) || "right".equalsIgnoreCase(ktrim)){
				argList.add("tbo");
				argList.add("tpe");
			}
			argList.add("k="+trim_k);
			argList.add("hdist="+hdist_trim);
			if(hdist2_trim>=0){
				argList.add("hdist2="+hdist2_trim);
			}
			if(forceTrimModulo>0){
				argList.add("ftm="+forceTrimModulo);
			}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			if(zl!=null){argList.add("zl="+zl);}
			
			//Pass along uncaptured arguments
			for(String s : primaryArgList){argList.add(s);}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+in1);}
			if(in2!=null){argList.add("in2="+in2);}
			if(out1!=null){argList.add("out1="+(tmpDir==null ? outDir : tmpDir)+outPrefix+out1);}
			if(out2!=null){argList.add("out2="+(tmpDir==null ? outDir : tmpDir)+outPrefix+out2);}
			if(qfin1!=null){argList.add("qfin1="+qfin1);}
			if(qfin2!=null){argList.add("qfin2="+qfin2);}
			if(qfout1!=null){argList.add("qfout1="+(tmpDir==null ? outDir : tmpDir)+outPrefix+qfout1);}
			if(qfout2!=null){argList.add("qfout2="+(tmpDir==null ? outDir : tmpDir)+outPrefix+qfout2);}

//			if(rqcStatsName!=null){al.add("rqc="+rqcStatsName_kt);} //Old style for 2 log files
			if(rqcStatsName!=null){argList.add("rqc=hashmap");}
			if(kmerStatsName!=null){argList.add("outduk="+kmerStatsName_kt);}
			if(scaffoldStatsName!=null){argList.add("stats="+scaffoldStatsName_kt);}
			
			if(copyUndefined){argList.add("cu");}
		}
		
		{//Add BBDuk references
			trimrefs.add(fragAdapters);

			StringBuilder refstring=new StringBuilder();
			for(String ref : trimrefs){
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
		
		{//run BBDuk
			BBDuk duk=new BBDuk(dukargs);
			try {
				duk.process();
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
	 * Runs BBDuk to perform:
	 * Quality filtering, quality trimming, n removal, short read removal, artifact removal (via kmer filtering), phiX removal.
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param qfin1 Primary input qual file
	 * @param qfin2 Secondary input qual file
	 * @param qfout1 Primary output qual file
	 * @param qfout2 Secondary output qual file
	 * @param inPrefix Append this prefix to input filenames
	 */
	private void filter(String in1, String in2, String out1, String out2, String qfin1, String qfin2, String qfout1, String qfout2,
			String inPrefix, String outPrefix, boolean prependIndir, boolean prependOutdir, boolean lastPhase){
		
		log("filter start", true);
		System.err.println("\nArtifact Filter/Quality Trim Phase Start");
		
		ArrayList<String> argList=new ArrayList<String>();
		
		{//Fill list with BBDuk arguments
			if(minAvgQuality>-1){argList.add("maq="+minAvgQuality+","+minAvgQualityBases);}
			if(maxNs>=0){argList.add("maxns="+maxNs);}
			if(minLen>0){argList.add("minlen="+minLen);}
			if(minLenFraction>0){argList.add("minlenfraction="+minLenFraction);}
			argList.add("k="+filter_k);
			argList.add("hdist="+hdist_filter);
			
			if(qtrim!=null && trimAfterFiltering){
				argList.add("trimq="+trimq);
				argList.add("qtrim="+qtrim);
			}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			if(zl!=null){argList.add("zl="+zl);}
			
			//Pass along uncaptured arguments
			for(String s : primaryArgList){argList.add(s);}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+in1);}
			if(in2!=null){argList.add("in2="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+in2);}
			if(out1!=null){argList.add("out1="+(prependOutdir ? (tmpDir==null || lastPhase ? outDir : tmpDir) : "")+outPrefix+out1);}
			if(out2!=null){argList.add("out2="+(prependOutdir ? (tmpDir==null || lastPhase ? outDir : tmpDir) : "")+outPrefix+out2);}
			if(qfin1!=null){argList.add("qfin1="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+qfin1);}
			if(qfin2!=null){argList.add("qfin2="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+qfin2);}
			if(qfout1!=null){argList.add("qfout1="+(prependOutdir ? (tmpDir==null || lastPhase ? outDir : tmpDir) : "")+outPrefix+qfout1);}
			if(qfout2!=null){argList.add("qfout2="+(prependOutdir ? (tmpDir==null || lastPhase ? outDir : tmpDir) : "")+outPrefix+qfout2);}
			
//			if(rqcStatsName!=null){al.add("rqc="+rqcStatsName);} //Old style for 2 log files
			if(rqcStatsName!=null){argList.add("rqc=hashmap");}
			if(kmerStatsName!=null){argList.add("outduk="+kmerStatsName);}
			if(scaffoldStatsName!=null){argList.add("stats="+scaffoldStatsName);}
		}
		
		{//Add BBDuk references
			filterrefs.add(mainArtifactFile);
			filterrefs.add(rnaFlag ? artifactFileRna : artifactFileDna);

			if(phixFlag){filterrefs.add(phixRef);}
			if(lambdaFlag){filterrefs.add(lambdaRef);}
			if(pjetFlag){filterrefs.add(pjetRef);}
			
			if(copyUndefined){argList.add("cu");}
			

			StringBuilder refstring=new StringBuilder();
			for(String ref : filterrefs){
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
		
		{//Run BBDuk
			BBDuk duk=new BBDuk(dukargs);
			try {
				duk.process();
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
	 * Runs BBMap to perform:
	 * Removal of reads that map to human with high identity (~88%).
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param qfin1 Primary input qual file
	 * @param qfin2 Secondary input qual file
	 * @param qfout1 Primary output qual file
	 * @param qfout2 Secondary output qual file
	 * @param inPrefix Append this prefix to input filenames
	 */
	private void dehumanize(String in1, String in2, String out1, String out2, String qfin1, String qfin2,
			String inPrefix, String outPrefix, boolean prependIndir, boolean prependOutdir, boolean lastPhase){
		
		log("dehumanize start", true);
		System.err.println("\nHuman Removal Phase Start");
		
		ArrayList<String> argList=new ArrayList<String>();
		
		{
			
			if(kfilter>map_k){argList.add("kfilter="+kfilter);}
			if(local){argList.add("local");}
			argList.add("minratio="+minratio);
			argList.add("maxindel="+maxindel);
			argList.add("fast="+fast);
			argList.add("minhits="+minhits);
			argList.add("tipsearch="+Tools.min(4, maxindel));
			argList.add("bw=18");
			argList.add("bwr=0.18");
			argList.add("quickmatch=f");
			argList.add("k="+map_k);
//			argList.add("cigar=f");
			argList.add("idtag=t");
			argList.add("sam=1.4");
			argList.add("usemodulo");
			argList.add("printunmappedcount");
			argList.add("ow="+overwrite);
			
			if(mapRef==null){
				argList.add("path="+indexPath);
			}else{
				argList.add("ref="+mapRef);
				argList.add("nodisk");
			}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			if(zl!=null){argList.add("zl="+zl);}
			
			//Pass along uncaptured arguments
			for(String s : primaryArgList){argList.add(s);}

			//Set read I/O files
			if(in1!=null){argList.add("in1="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+in1);}
			if(in2!=null){argList.add("in2="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+in2);}
			if(out1!=null){argList.add("outu1="+(prependOutdir ? (tmpDir==null || lastPhase ? outDir : tmpDir) : "")+outPrefix+out1);}
			if(out2!=null){argList.add("outu2="+(prependOutdir ? (tmpDir==null || lastPhase ? outDir : tmpDir) : "")+outPrefix+out2);}
			if(qfin1!=null){argList.add("qfin1="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+qfin1);}
			if(qfin2!=null){argList.add("qfin2="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+qfin2);}
			
		}
		
		String[] args=argList.toArray(new String[0]);
		
		{//Run BBMap
			try {
				BBMap.main(args);
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Optionally append files to file list here
		
		log("dehumanize finish", true);
	}
	
	/**
	 * Runs BBNorm to perform:
	 * Error correction, error marking, quality trimming, normalization
	 * 
	 * @param in1 Primary input reads file (required)
	 * @param in2 Secondary input reads file
	 * @param out1 Primary output reads file (required)
	 * @param out2 Secondary output reads file
	 * @param qfin1 Primary input qual file
	 * @param qfin2 Secondary input qual file
	 * @param qfout1 Primary output qual file
	 * @param qfout2 Secondary output qual file
	 * @param inPrefix Append this prefix to input filenames
	 */
	private void normalize(String in1, String in2, String out1, String out2, String qfin1, String qfin2, String qfout1, String qfout2,
			String inPrefix, String outPrefix, boolean prependIndir, boolean prependOutdir, boolean lastPhase){
		
		log("normalization start", true);
		System.err.println("\nNormalization/Error Correction Phase Start");
		
		ArrayList<String> argList=new ArrayList<String>();
		
		{//Fill list with BBDuk arguments
			if(qtrim!=null && !trimAfterFiltering){
				argList.add("trimq="+trimq);
				argList.add("qtrim="+qtrim);
			}
			if(minLen>0){argList.add("minlen="+minLen);}
			//if(minLenFraction>0){argList.add("minlenfraction="+minLenFraction);}

			argList.add("ecc="+ecc);
			if(aecc){argList.add("aec="+aecc);}
			if(cecc){argList.add("cecc="+cecc);}
			argList.add("meo="+meo);
			argList.add("tam="+tam);
			argList.add("mue="+mue);
			argList.add("mw1="+mw1);
			argList.add("prefilter=t");
			argList.add("prehashes="+prehashes);
			argList.add("hashes="+hashes);
			argList.add("bits="+bits);
			argList.add("k="+normalize_k);
			argList.add("passes="+passes);
			if(normalize){
				if(target>0){
					argList.add("target="+target);
					if(mindepth<0){mindepth=Tools.min(10, target/8);}
					if(maxdepth<0){maxdepth=Tools.max(target, (int)((target*17L)/16L));}
				}
				if(mindepth>=0){argList.add("min="+mindepth);}
				if(maxdepth>0){argList.add("max="+maxdepth);}
			}else{
				argList.add("keepall");
			}
			if(pigz!=null){argList.add("pigz="+pigz);}
			if(unpigz!=null){argList.add("unpigz="+unpigz);}
			if(zl!=null){argList.add("zl="+zl);}
			
			//Set read I/O files
			if(in1!=null){argList.add("in1="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+in1);}
			if(in2!=null){argList.add("in2="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+in2);}
//			if(out1!=null){argList.add("out="+outDir+out1);}
			if(out1!=null){argList.add("out1="+(prependOutdir ? (tmpDir==null || lastPhase ? outDir : tmpDir) : "")+outPrefix+out1);}
			if(out2!=null){argList.add("out2="+(prependOutdir ? (tmpDir==null || lastPhase ? outDir : tmpDir) : "")+outPrefix+out2);}
//			if(out2!=null){argList.add("out2="+outDir+out2);}
			if(qfin1!=null){argList.add("qfin1="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+qfin1);}
			if(qfin2!=null){argList.add("qfin2="+(prependIndir ? (tmpDir==null ? outDir : tmpDir) : "")+inPrefix+qfin2);}
//			if(qfout1!=null){argList.add("qfout1="+outDir+qfout1);}
//			if(qfout2!=null){argList.add("qfout2="+outDir+qfout2);}\
			
			
			if(kmerHistName!=null){argList.add("hist="+kmerHistName);}
		}
		
		String[] normargs=argList.toArray(new String[0]);
		
		{//Run BBNorm
			try {
				KmerNormalize.main(normargs);
			} catch (Exception e) {
				e.printStackTrace();
				log("failed", true);
				System.exit(1);
			}
		}
		
		//Optionally append files to file list here
		
		log("normalization finish", true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Log a message in the log file
	 * @param message Message to log
	 * @param append True to append, false to overwrite
	 */
	private void log(String message, boolean append){
		if(logName!=null){
			ReadWrite.writeString(message+", "+timeString()+"\n", logName, append);
		}
	}
	
	
	/**
	 * Delete all non-null filenames.
	 * @param prefix Append this prefix to filenames before attempting to delete them
	 * @param names Filenames to delete
	 */
	private void delete(String prefix, String...names){
		log("delete temp files start", true);
		if(names!=null){
			for(String s : names){
				if(s!=null){
					s=(tmpDir==null ? outDir : tmpDir)+prefix+s;
					if(verbose){System.err.println("Trying to delete "+s);}
					File f=new File(s);
					if(f.exists()){
						f.delete();
					}
				}
			}
		}
		log("delete temp files finish", true);
	}
	
	
	/**
	 * Delete all non-null filenames.
	 * @param prefix Append this prefix to filenames before attempting to delete them
	 * @param names Filenames to delete
	 */
	private void move(String prefix, String...names){
		log("delete temp files start", true);
		if(names!=null){
			for(String s : names){
				if(s!=null){
					s=(tmpDir==null ? outDir : tmpDir)+prefix+s;
					if(verbose){System.err.println("Trying to delete "+s);}
					File f=new File(s);
					if(f.exists()){
						f.delete();
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
		
		if(mainArtifactFile!=null || (rnaFlag ? artifactFileRna!=null : artifactFileDna!=null)){sb.append("a");}
		
		if(maxNs>=0){sb.append("n");}
//		if(qtrim!=null && !qtrim.equalsIgnoreCase("f") && !qtrim.equalsIgnoreCase("false")){sb.append("q");}
		if(minAvgQuality>0){sb.append("q");}
		
		if(rnaFlag){sb.append("r");}
		else{sb.append("d");}
		
		if(phixFlag){sb.append("p");}
		
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
	
	private static String stripDirs(String path){
		return ReadWrite.stripPath(path);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      BBNorm Parameters       ----------------*/
	/*--------------------------------------------------------------*/

	private boolean removehuman=true;
	private boolean normalize=false;
	private boolean ecc=false;
	private boolean aecc=false;
	private boolean cecc=false;
	private boolean meo=false;
	private boolean tam=false;
	private boolean trimAfterFiltering=true;
	private boolean mue=false;
	private boolean mw1=false;
	private int maxdepth=-1;
	private int mindepth=6;
	private int target=50;
	private int prehashes=3;
	private int passes=2;
	private int hashes=4;
	private int bits=16;
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Symbols to insert in output filename to denote operations performed */
	private final String symbols;
	
	/** True for rna artifacts, false for dna artifacts */
	private boolean rnaFlag=false;
	/** True if phix should be filtered out */
	private boolean phixFlag=true;
	/** True if lambda should be filtered out by kmer-match */
	private boolean lambdaFlag=true;
	/** True if pjet should be filtered out */
	private boolean pjetFlag=true;
	
	/** Unused */
	private boolean tboFlag=false;
	/** Unused */
	private boolean tpeFlag=false;
	
	/** Toss reads shorter than this */
	private int minLen=40;
	/** Toss reads shorter than this fraction of initial length, after trimming */
	private float minLenFraction=0.6f;
	/** Trim bases at this quality or below */
	private float trimq=12;
	/** Throw away reads below this average quality after trimming.  Default: 8 */
	private float minAvgQuality=8;
	/** If positive, calculate the average quality from the first X bases. */
	private int minAvgQualityBases=0;
	
	/** Trim reads to be equal to 0 modulo this value.  Mainly for 151, 251, and 301bp runs. */
	private int forceTrimModulo=5;
	/** Quality-trimming mode */
	private String qtrim="rl";
	/** Kmer-trimming mode */
	private String ktrim="r";
	/** Kmer length to use for filtering */
	private int filter_k=27;
	/** Kmer length to use for trimming */
	private int trim_k=23;
	/** Kmer length to use for normalization and error-correction */
	private int normalize_k=31;
	/** Kmer length to use for mapping */
	private int map_k=13;
	/** Shortest kmer to use for trimming */
	private int mink=11;
	/** Throw away reads containing more than this many Ns.  Default: 1 */
	private int maxNs=1;
	/** Use this Hamming distance when kmer filtering */
	private int hdist_filter=1;
	/** Use this Hamming distance when kmer trimming */
	private int hdist_trim=1;
	/** Use this Hamming distance when kmer trimming with short kmers */
	private int hdist2_trim=-1;
	
	/** Captures the command line "pigz" flag */
	private String pigz;
	/** Captures the command line "unpigz" flag */
	private String unpigz;
	/** Captures the command line "zl" flag */
	private String zl;
	
	private float minratio=0.84f;
	private int maxindel=6;
	private int kfilter=0;
	private int minhits=1;
	private boolean fast=true;
	private boolean local=true;
	private boolean copyUndefined=false;
	
	private boolean verbose=false;
	private boolean overwrite=true;
//	private boolean append=false;
	private boolean compress=true;
	
	/** Arguments to pass to BBDuk */
	private ArrayList<String> primaryArgList=new ArrayList<String>();
	/** References to pass to BBDuk for artifact removal */
	private ArrayList<String> trimrefs=new ArrayList<String>();
	/** References to pass to BBDuk for artifact removal */
	private ArrayList<String> filterrefs=new ArrayList<String>();
	
	/*--------------------------------------------------------------*/
	/*----------------        Read Data Files       ----------------*/
	/*--------------------------------------------------------------*/

	/** Directory in which to write all files */
	private String outDir="";
	
	/** Directory in which to write all temp files */
	private String tmpDir=Shared.tmpdir();
	
	private final String tempSalt;
	
	/** Primary input reads file (required) */
	private String in1=null;
	/** Secondary input reads file */
	private String in2=null;
	/** Primary output reads file (required) */
	private String out1=null;
	/** Secondary output reads file */
	private String out2=null;
	/** Primary input qual file */
	private String qfin1=null;
	/** Secondary input qual file */
	private String qfin2=null;
	/** Primary output qual file */
	private String qfout1=null;
	/** Secondary output qual file */
	private String qfout2=null;
	
	/*--------------------------------------------------------------*/
	/*----------------           Log Files          ----------------*/
	/*--------------------------------------------------------------*/
	
	private String logName="status.log";
	private String fileListName="file-list.txt";
	
	private String rqcStatsName="filterStats.txt";
	private String kmerStatsName="kmerStats.txt";
	private String scaffoldStatsName="scaffoldStats.txt";
	private String kmerHistName="khist.txt";
	
	private String ihistName=null;
	
	/** ktrim phase rqc stats file */
	private String rqcStatsName_kt;
	/** ktrim phase stats file */
	private String kmerStatsName_kt;
	/** ktrim phase scaffold stats file */
	private String scaffoldStatsName_kt;
	
	/*--------------------------------------------------------------*/
	/*----------------        Reference Files       ----------------*/
	/*--------------------------------------------------------------*/
	
	private String mainArtifactFile = "/global/dna/shared/rqc/ref_databases/qaqc/databases/illumina.artifacts/Illumina.artifacts.2013.12.no_DNA_RNA_spikeins.fa";
	private String artifactFileRna = "/global/dna/shared/rqc/ref_databases/qaqc/databases/illumina.artifacts/RNA_spikeins.artifacts.2012.10.NoPolyA.fa";
	private String artifactFileDna = "/global/dna/shared/rqc/ref_databases/qaqc/databases/illumina.artifacts/DNA_spikeins.artifacts.2012.10.fa";
	private String phixRef = "/global/dna/shared/rqc/ref_databases/qaqc/databases/phix174_ill.ref.fa";
	private String lambdaRef = "/global/dna/shared/rqc/ref_databases/qaqc/databases/lambda.fa.gz";
	private String pjetRef = "/global/dna/shared/rqc/ref_databases/qaqc/databases/pJET1.2.fasta";

	private String allArtifactsLatest = "/global/projectb/sandbox/rqc/qcdb/illumina.artifacts/Illumina.artifacts.fa";
	private String fragAdapters = "/global/cfs/cdirs/bbtools/data/adapters.fa";
	private String rnaAdapter = "/global/cfs/cdirs/bbtools/data/truseq_rna.fa.gz";
	private String indexPath = "/global/cfs/cdirs/bbtools/hg19/";
	private String mapRef = null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
}
