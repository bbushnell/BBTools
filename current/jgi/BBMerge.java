package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

import align2.QualityTools;
import assemble.Tadpole;
import bloom.BloomFilter;
import bloom.BloomFilterCorrector;
import bloom.KmerCountAbstract;
import cardinality.CardinalityTracker;
import dna.AminoAcid;
import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.KmerTableSet;
import ml.CellNet;
import ml.CellNetParser;
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
import stream.CustomHeader;
import stream.FASTQ;
import stream.Read;
import stream.ReadStreamWriter;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntList;
import structures.ListNum;
import structures.LongList;
import structures.Quantizer;
import tracker.AdapterTracker;
import tracker.ReadStats;
import ukmer.Kmer;

/**
 * @author Brian Bushnell
 * @date Aug 14, 2012
 *
 */
public class BBMerge {
	
	public static void main(String[] args){
//		boolean old=Shared.USE_JNI;
//		Shared.USE_JNI=false; //TODO: This is for RQCFilter.  Can be removed.
		
		MAKE_VECTOR=BBMergeOverlapper.MAKE_VECTOR=true;//This doesn't work in a static context, so...
		BBMerge x=new BBMerge(args);
		x.process();
		MAKE_VECTOR=BBMergeOverlapper.MAKE_VECTOR=false;
//		Shared.USE_JNI=old;
		Read.VALIDATE_IN_CONSTRUCTOR=true;
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	private static String[] preparse(String[] args){
		if(args==null){return new String[0];}
		int nulls=0;
		final int XSTRICT=1, USTRICT=2, VSTRICT=3, STRICT=4, NORMAL=5, LOOSE=6, VLOOSE=7, ULOOSE=8, XLOOSE=9, FAST=10;
		int mode=NORMAL;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("jni") || a.equals("usejni")){
				Shared.USE_JNI=Parse.parseBoolean(b);
			}else if(a.equals("showfullargs") || a.equalsIgnoreCase("showFullArgs")){
				showFullArgs=Parse.parseBoolean(b);
				args[i]=null;
				nulls++;
			}else if(a.equals("vstrict") || a.equals("verystrict")){
				boolean x=Parse.parseBoolean(b);
				if(x){mode=VSTRICT;}else if(mode==VSTRICT){mode=NORMAL;}
				args[i]=null;
				nulls++;
			}else if(a.equals("ustrict") || a.equals("ultrastrict")){
				boolean x=Parse.parseBoolean(b);
				if(x){mode=USTRICT;}else if(mode==USTRICT){mode=NORMAL;}
				args[i]=null;
				nulls++;
			}else if(a.equals("xstrict") || a.equals("hstrict") || a.equals("hyperstrict") || a.equals("maxstrict")){
				boolean x=Parse.parseBoolean(b);
				if(x){mode=XSTRICT;}else if(mode==XSTRICT){mode=NORMAL;}
				args[i]=null;
				nulls++;
			}else if(a.equals("strict")){
				boolean x=Parse.parseBoolean(b);
				if(x){mode=STRICT;}else if(mode==STRICT){mode=NORMAL;}
				args[i]=null;
				nulls++;
			}else if(a.equals("loose")){
				boolean x=Parse.parseBoolean(b);
				if(x){mode=LOOSE;}else if(mode==LOOSE){mode=NORMAL;}
				args[i]=null;
				nulls++;
			}else if(a.equals("vloose") || a.equals("veryloose")){
				boolean x=Parse.parseBoolean(b);
				if(x){mode=VLOOSE;}else if(mode==VLOOSE){mode=NORMAL;}
				args[i]=null;
				nulls++;
			}else if(a.equals("uloose") || a.equals("ultraloose")){
				boolean x=Parse.parseBoolean(b);
				if(x){mode=ULOOSE;}else if(mode==ULOOSE){mode=NORMAL;}
				args[i]=null;
				nulls++;
			}else if(a.equals("xloose") || a.equals("hloose") || a.equals("hyperloose") || a.equals("maxloose")){
				boolean x=Parse.parseBoolean(b);
				if(x){mode=XLOOSE;}else if(mode==XLOOSE){mode=NORMAL;}
				args[i]=null;
				nulls++;
			}else if(a.equals("fast")){
				boolean x=Parse.parseBoolean(b);
				if(x){mode=FAST;}else if(mode==FAST){mode=NORMAL;}
				args[i]=null;
				nulls++;
			}else if(a.equals("default")){
				boolean x=Parse.parseBoolean(b);
				if(x){mode=NORMAL;}
				args[i]=null;
				nulls++;
			}
		}
		
		if(mode==FAST){fast=true;}
		else if(mode==XSTRICT){xstrict=true;}
		else if(mode==USTRICT){ustrict=true;}
		else if(mode==VSTRICT){vstrict=true;}
		else if(mode==STRICT){strict=true;}
		else if(mode==LOOSE){loose=true;}
		else if(mode==VLOOSE){vloose=true;}
		else if(mode==ULOOSE){uloose=true;}
		else if(mode==XLOOSE){xloose=true;}
		
		if(nulls==0){return args;}
		ArrayList<String> args2=new ArrayList<String>(args.length-nulls+5);
		if(strict || vstrict || ustrict || xstrict){
			strict=true;
			loose=vloose=uloose=xloose=false;
			
			args2.add("maxbad=4");
			args2.add("margin=3");
			args2.add("minqo=8");
			args2.add("qualiters=2");
			
			if(xstrict){
				args2.add("ratiomode=t");
				args2.add("flatmode=t");
				args2.add("requireratiomatch=t");

				args2.add("minentropy=56");
				args2.add("minoverlap=14");
				args2.add("minoverlap0=3");
				
				args2.add("maxratio=0.055");
				args2.add("ratiomargin=12");
				args2.add("ratiooffset=0.65");
				args2.add("ratiominoverlapreduction=4");
				args2.add("efilter=2");
				args2.add("pfilter=0.25");
				args2.add("minsecondratio=0.24");
				args2.add("minapproxoverlap=18");
			}else if(ustrict){
				args2.add("ratiomode=t");
				args2.add("flatmode=t");
				args2.add("requireratiomatch=t");

				args2.add("minentropy=56");
				args2.add("minoverlap=14");
				args2.add("minoverlap0=3");
				
				args2.add("maxratio=0.045");
				args2.add("ratiomargin=12");
				args2.add("ratiooffset=0.5");
				args2.add("ratiominoverlapreduction=4");
				args2.add("efilter=2");
				args2.add("pfilter=0.03");
				args2.add("minsecondratio=0.20");
				args2.add("minapproxoverlap=20");
			}else if(vstrict){
				args2.add("ratiomode=t");
				args2.add("flatmode=f");

				args2.add("minentropy=52");
				args2.add("minoverlap=12");
				args2.add("minoverlap0=4");

				args2.add("maxratio=0.05");
				args2.add("ratiomargin=12");
				args2.add("ratiooffset=0.5");
				args2.add("ratiominoverlapreduction=4");
				args2.add("efilter=2");
				args2.add("pfilter=0.008");
				args2.add("minsecondratio=0.16");
				args2.add("minapproxoverlap=22");
			}else{
				args2.add("ratiomode=t");
				args2.add("flatmode=f");
				
				args2.add("minentropy=42");
				args2.add("minoverlap0=7");
				args2.add("minoverlap=11");
				
				args2.add("maxratio=0.075");
				args2.add("ratiomargin=7.5");
				args2.add("ratiooffset=0.55");
				args2.add("ratiominoverlapreduction=4");
				args2.add("efilter=4");
				args2.add("pfilter=0.0008");
				args2.add("minsecondratio=0.12");
				args2.add("minapproxoverlap=24");
			}
		}else if(loose || vloose || uloose || xloose){
			loose=true;
			strict=vstrict=ustrict=xstrict=false;
			args2.add("minoverlap=8");
			args2.add("minoverlap0=9");
			args2.add("qualiters=4");
			args2.add("mismatches=3");
			args2.add("margin=2");
			
			args2.add("ratiooffset=0.4");
			args2.add("minsecondratio=0.08");
			
			if(xloose){
				args2.add("owq=t");
				args2.add("ouq=t");
				args2.add("minentropy=22");
				args2.add("minoverlap=8");
				args2.add("minoverlap0=7");
				args2.add("maxratio=0.2");
				args2.add("mismatches=3");
				args2.add("ratiomargin=2");
				args2.add("flatmode=t");
				args2.add("pfilter=0.0000001");
				args2.add("efilter=8");
				args2.add("margin=2");
				args2.add("ratiominoverlapreduction=2");
				args2.add("minapproxoverlap=38");
			}else if(vloose || uloose){
				args2.add("owq=t");
				args2.add("ouq=t");
				if(uloose){
//					args2.add("maxratio=0.14");
//					args2.add("ratiomargin=2");
//					args2.add("flatmode=t");
//					args2.add("pfilter=0.0000001");
					
					
					args2.add("minoverlap=8");
					args2.add("minoverlap0=7");
					args2.add("mismatches=3");
					args2.add("margin=2");

					args2.add("ratiominoverlapreduction=2");
					args2.add("efilter=8");
					args2.add("maxratio=0.16");
					args2.add("ratiomargin=2.2");
					args2.add("pfilter=0.0000002");
					args2.add("minentropy=24");
					args2.add("minapproxoverlap=34");
				}else{
					args2.add("ratiominoverlapreduction=3");
					args2.add("maxratio=0.12");
					args2.add("ratiomargin=3");
					args2.add("pfilter=0.000004");
					args2.add("minentropy=28");
					args2.add("efilter=7.5");
					args2.add("ratiooffset=0.45");
					args2.add("minapproxoverlap=32");
				}
			}else{
				args2.add("maxratio=0.11");
				args2.add("ratiomargin=4.7");
				args2.add("ratiominoverlapreduction=2");
				args2.add("pfilter=0.00002");
				args2.add("efilter=8");
				args2.add("minentropy=30");
				args2.add("minapproxoverlap=30");
			}
		}else if(fast){
			args2.add("maxratio=0.08");
			args2.add("ratiomargin=2.5");
			args2.add("ratiominoverlapreduction=3");
			args2.add("pfilter=0.0002");
			args2.add("efilter=8");
			args2.add("minentropy=39");
			args2.add("mininsert0=50");
			args2.add("minsecondratio=0.08");
		}
		
		for(String s : args){
			if(s!=null){args2.add(s);}
		}
		return args2.toArray(new String[args2.size()]);
	}
	
	public BBMerge(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), true);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		{
			String[] args0=args;
			args=preparse(args);

			if(args0!=args && showFullArgs){
				outstream.println("Revised arguments: "+Arrays.toString(args)+"\n");
			}
		}
		
		Timer ttotal=new Timer();
		ttotal.start();
		
		in1=(args[0].indexOf('=')>0 ? null : args[0]);
		in2=(in1!=null && args.length>1 && args[1].indexOf('=')<0 ? args[1] : null);
		if(in2!=null && "null".equalsIgnoreCase(in2)){in2=null;}
		
		{
			if(in1!=null && !in1.contains(",") && !in1.startsWith("stdin.") && !in1.equals("stdin")){
				File f=new File(in1);
				if(!f.exists() || !f.isFile()){
					in1=null;
//					throw new RuntimeException(in1+" does not exist.");
				}
			}
			if(in2!=null && !in2.contains(",")){
				File f=new File(in2);
				if(!f.exists() || !f.isFile()){
					in2=null;
//					throw new RuntimeException(in2+" does not exist.");
				}else if(in1.equalsIgnoreCase(in2)){
					throw new RuntimeException("Both input files are the same.");
				}
			}
		}
		
		ReadWrite.setZipThreads(Tools.max(Shared.threads()>1 ? 2 : 1, Shared.threads()>20 ? Shared.threads()/2 : Shared.threads()));
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		Shared.setBufferLen(Tools.max(Shared.bufferLen(), 400));
		boolean setMix=false;
		
		boolean mm0set=false;
		
		Parser parser=new Parser();
		parser.trimq2=trimq;
		parser.minAvgQuality=minAvgQuality;
		parser.minReadLength=minReadLength;
		parser.maxReadLength=maxReadLength;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("parsecustom")){
				parseCustom=Parse.parseBoolean(b);
//			}else if(a.equals("tagcustom")){
//				TAG_CUSTOM=Parse.parseBoolean(b);
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseQualityAdjust(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(parser.parseCardinality(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("extra")){
				if(b!=null){
					for(String s : b.split(",")){
						extra.add(s);
					}
				}
			}else if(a.equals("useratio") || a.equals("ratio") || a.equals("ratiomode")){
				useRatioMode=Parse.parseBoolean(b);
			}else if(a.equals("useflatmode") || a.equals("flatmode") || a.equals("usenormalmode") || a.equals("normalmode")){
				useFlatMode=Parse.parseBoolean(b);
			}else if(a.equals("requireratiomatch") || a.equals("rrm")){
				requireRatioMatch=Parse.parseBoolean(b);
			}else if(a.equals("makevector") || a.equals("nn")){
				MAKE_VECTOR=BBMergeOverlapper.MAKE_VECTOR=Parse.parseBoolean(b);
			}else if(a.equals("net")){
				netFile=b;
			}else if(a.equals("netcutoff") || a.equals("cutoff")){
				netCutoff=Float.parseFloat(b);
				assert(netCutoff>=0 && netCutoff<=1) : netCutoff;
				setCutoff=true;
			}else if(a.equals("maxratio")){
				MAX_RATIO=Float.parseFloat(b);
//				useRatioMode=true;
			}else if(a.equals("maxmismatchesr") || a.equals("maxmismatches")){
				MAX_MISMATCHES_R=Integer.parseInt(b);
//				useRatioMode=true;
			}else if(a.equals("ratiomargin")){
				RATIO_MARGIN=Float.parseFloat(b);
//				useRatioMode=true;
			}else if(a.equals("ratiooffset")){
				RATIO_OFFSET=Float.parseFloat(b);
//				useRatioMode=true;
			}else if(a.equals("ratiominoverlapreduction")){
				MIN_OVERLAPPING_BASES_RATIO_REDUCTION=Integer.parseInt(b);
//				useRatioMode=true;
			}else if(a.equals("minentropy") || a.equals("entropy")){
				if(b!=null && Tools.isDigit(b.charAt(0))){
					minEntropyScore=Integer.parseInt(b);
				}else{
					useEntropy=Parse.parseBoolean(b);
				}
			}else if(a.equals("minoverlappingbases") || a.equals("minoverlapbases") || a.equals("minoverlap")){
				MIN_OVERLAPPING_BASES=Integer.parseInt(b);
			}else if(a.equals("minoverlappingbases0") || a.equals("minoverlapbases0") || a.equals("minoverlap0")){
				MIN_OVERLAPPING_BASES_0=Integer.parseInt(b);
			}else if(a.equals("minqo") || a.equals("minq")){
				MIN_QUALITY=(byte)Integer.parseInt(b);
			}else if(a.equals("maxq")){
				Read.MAX_MERGE_QUALITY=(byte)Integer.parseInt(b);
			}else if(a.equals("qualiters")){
				QUAL_ITERS=Tools.max(1, Integer.parseInt(b));
			}else if(a.equals("maxbadbases") || a.equals("maxbad") || a.equals("mismatches")){
				MAX_MISMATCHES=Integer.parseInt(b);
			}else if(a.equals("maxbadbases0") || a.equals("maxbad0") || a.equals("mismatches0")){
				MAX_MISMATCHES0=Integer.parseInt(b);
				mm0set=true;
			}else if(a.equals("margin")){
				MISMATCH_MARGIN=Integer.parseInt(b);
			}else if(a.equals("usemapping")){
				USE_MAPPING=Parse.parseBoolean(b);
			}else if(a.equals("bin")){
				bin=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=Shared.setThreads(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("outgood") || a.equals("outmerged") || a.equals("outm") || a.equals("out")){
				out1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outgood1") || a.equals("outmerged1") || a.equals("outm1") || a.equals("out1")){
				out1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outgood2") || a.equals("outmerged2") || a.equals("outm2") || a.equals("out2")){
				out2=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb") || a.equals("outu") || a.equals("outunmerged") || a.equals("outbad")){
				outb1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb1") || a.equals("outu1") || a.equals("outunmerged1") || a.equals("outbad1")){
				outb1=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outb2") || a.equals("outu2") || a.equals("outunmerged2") || a.equals("outbad2")){
				outb2=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outinsert") || a.startsWith("outi") || a.startsWith("outlength")){
				outinsert=(b==null || b.equals("null") ? null : b);
			}else if(a.startsWith("outhist") || a.equals("hist") || a.equals("histogram") || a.equals("ihist")){
				ihist=(b==null || b.equals("null") ? null : b);
			}else if(a.equals("outa") || a.equals("outadapter")){
				outAdapter=b;
				findAdapterSequence=(outAdapter!=null);
			}else if(a.equals("ignorephix") || a.equals("ignorephixadapters")){
				ignorePhixAdapters=AdapterTracker.ignorePhixAdapters=Parse.parseBoolean(b);
			}else if(a.equals("outc") || a.equals("outcardinality")){
				outCardinality=b;
//			}else if(a.equals("outputfailed")){
//				OUTPUT_FAILED=Parse.parseBoolean(b);outCardinality
			}else if(a.equals("trimpolya")){
				trimPolyA=AdapterTracker.trimPolyAorG=Parse.parseBoolean(b);
			}else if(a.equals("mix")){
				MIX_BAD_AND_GOOD=Parse.parseBoolean(b);
				setMix=true;
			}else if(a.equals("ooi") || a.equals("onlyoutputincorrect")){
				ONLY_OUTPUT_INCORRECT=Parse.parseBoolean(b);
			}else if(a.equals("ofp") || a.equals("outputfalsepositive")){
				OUTPUT_FALSE_POSITIVE=Parse.parseBoolean(b);
			}else if(a.equals("oofp") || a.equals("onlyoutputfalsepositive")){
				ONLY_OUTPUT_FALSE_POSITIVE=Parse.parseBoolean(b);
			}else if(a.equals("ofn") || a.equals("outputfalsenegative")){
				OUTPUT_FALSE_NEGATIVE=Parse.parseBoolean(b);
			}else if(a.equals("otp") || a.equals("outputtruepositive")){
				OUTPUT_TRUE_POSITIVE=Parse.parseBoolean(b);
			}else if(a.equals("otn") || a.equals("outputtruenegative")){
				OUTPUT_TRUE_NEGATIVE=Parse.parseBoolean(b);
			}else if(a.equals("train") || a.equals("training")){
				TRAINING=Parse.parseBoolean(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				NONZERO_ONLY=Parse.parseBoolean(b);
			}else if(a.equals("exemptshortoverlaps") || a.equals("exemptshort") || a.equals("ignoreshort")){
				EXEMPT_SHORT_OVERLAPS_FROM_STATS=Parse.parseBoolean(b);
			}else if(a.equals("showhiststats")){
				showHistStats=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				assert(false) : "verbose flag is static final; recompile to change it.";
//				verbose=Parse.parseBoolean(b);
			}else if(a.equals("trimnonoverlapping") || a.equals("tno")){
				trimNonOverlapping=Parse.parseBoolean(b);
			}else if(a.equals("join") || a.equals("merge")){
				join=Parse.parseBoolean(b);
				if(join){ecco=false;}
			}else if(a.equals("ecco") || a.equals("ecc") || a.equals("errorcorrect")){
				ecco=Parse.parseBoolean(b);
				if(ecco){join=false;}
			}else if(a.equals("tbo") || a.equals("trimbyoverlap")){
				trimByOverlap=Parse.parseBoolean(b);
			}else if(a.equals("useoverlap") || a.equals("usebases") || a.equals("matebyoverlap") || a.equals("matebybases")){
				MATE_BY_OVERLAP=Parse.parseBoolean(b);
			}
//			else if(a.startsWith("skipmated")){
//				SKIP_MATED_READS=Parse.parseBoolean(b);
//			}
			else if(a.equals("lowercase")){
				lowercaseAdapters=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("trimonfailure") || a.equals("tof")){
				if(b!=null && Tools.isDigit(b.charAt(0))){
					TRIM_ON_OVERLAP_FAILURE=Integer.parseInt(b);
				}else{
					TRIM_ON_OVERLAP_FAILURE=(Parse.parseBoolean(b) ? 1 : 0);
				}
			}else if(a.equals("overlapusingquality") || a.equals("ouq")){
				overlapUsingQuality=Parse.parseBoolean(b);
			}else if(a.equals("overlapwithoutquality") || a.equals("owoq") || a.equals("owuq") || a.equals("owq")){
				overlapWithoutQuality=Parse.parseBoolean(b);
			}else if(a.equals("maxexpectederrors") || a.equals("mee") || a.equals("meefilter")){
				maxExpectedErrors=Float.parseFloat(b);
			}else if(a.equals("mi") || a.equals("minins") || a.equals("mininsert")){
				minInsert=Integer.parseInt(b);
			}else if(a.equals("mi0") || a.equals("mininsert0")){
				minInsert0=Integer.parseInt(b);
			}else if(a.equals("minprob")){
				minProb=Float.parseFloat(b);
				assert(minProb<1) : "minprob must be less than 1.  At 1, even kmers with 100% probablity of being error-free will be discarded.";
			}else if(a.equals("prealloc")){
				prealloc=Parse.parseBoolean(b);
			}else if(a.equals("prefilter")){
				if(b==null){prefilter=2;}
				else if(Character.isLetter(b.charAt(0))){prefilter=Parse.parseBoolean(b) ? 2 : 0;}
				else{prefilter=Integer.parseInt(b);}
			}else if(a.equalsIgnoreCase("filterMemoryOverride") || a.equalsIgnoreCase("filterMemory") || 
					a.equalsIgnoreCase("prefilterMemory") || a.equalsIgnoreCase("filtermem")){
				filterMemoryOverride=Parse.parseKMG(b);
			}else if(a.equals("k")){
				kmerLength=Integer.parseInt(b);
			}else if(a.equals("tail")){
				eccTail=Parse.parseBoolean(b);
			}else if(a.equals("pincer")){
				eccPincer=Parse.parseBoolean(b);
			}else if(a.equals("reassemble")){
				eccReassemble=Parse.parseBoolean(b);
			}
			
			else if(a.equals("efilter")){
				if(b==null || Character.isLetter(b.charAt(0))){
					boolean x=Parse.parseBoolean(b);
					if(!x){efilterRatio=-1;}
				}else{
					efilterRatio=Float.parseFloat(b);
				}
				useEfilter=efilterRatio>=0;
			}else if(a.equals("pfilter")){
				if(b==null || Character.isLetter(b.charAt(0))){
					boolean x=Parse.parseBoolean(b);
					if(!x){pfilterRatio=0;}
				}else{
					pfilterRatio=Float.parseFloat(b);
				}
			}else if(a.equals("efilteroffset")){
				efilterOffset=Float.parseFloat(b);
			}else if(a.equals("kfilter")){
				if(b!=null && Tools.isDigit(b.charAt(0))){
					filterCutoff=Integer.parseInt(b);
					useKFilter=filterCutoff>0;
				}else{
					useKFilter=Parse.parseBoolean(b);
				}
			}else if(a.equals("usequality")){
				useQuality=Parse.parseBoolean(b);
			}else if(a.equals("ignorequality")){
				useQuality=!Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("samplerate")){
				samplerate=Float.parseFloat(b);
				assert(samplerate<=1f && samplerate>=0f) : "samplerate="+samplerate+"; should be between 0 and 1";
			}else if(a.equals("sampleseed")){
				sampleseed=Long.parseLong(b);
			}else if(a.equals("recalibrate") || a.equals("recalibratequality") || a.equals("recal")){
				recalibrateQuality=Parse.parseBoolean(b);
			}else if(a.equals("recalpairnum") || a.equals("recalibratepairnum")){
				CalcTrueQuality.USE_PAIRNUM=Parse.parseBoolean(b);
			}else if(a.equals("path")){
				Data.setPath(b);
			}else if(a.equals("iupacton") || a.equals("itn")){
				iupacToN=Parse.parseBoolean(b);
			}else if(a.equals("minsecondratio")){
				MIN_SECOND_RATIO=Float.parseFloat(b);
			}
			
			//Adapter parameters
			
			else if(a.equals("adapter") || a.equals("adapters")){
				adapterList1=adapterList2=getAdapterList(b);
			}else if(a.equals("adapter1") || a.equals("adapters1")){
				adapterList1=adapterList2=getAdapterList(b);
			}else if(a.equals("adapter2") || a.equals("adapters2")){
				adapterList1=adapterList2=getAdapterList(b);
			}
			
			//Extension parameters
			
			else if(a.equals("extendright") || a.equals("er") || a.equals("extend") || a.equals("extendright1") || a.equals("er1") || a.equals("extend1")){
				extendRight1=Parse.parseIntKMG(b);
			}else if(a.equals("extendright2") || a.equals("er2") || a.equals("extend2")){
				extendRight2=Parse.parseIntKMG(b);
			}else if(a.equals("extenditerations") || a.equals("iterations") || a.equals("ei") || a.equals("iters")){
				extendIterations=Tools.max(1, (int)Parse.parseKMG(b));
			}else if(a.equals("ecctadpole") || a.equals("ecct")){
				eccTadpole=Parse.parseBoolean(b);
			}else if(a.equals("eccbloom") || a.equals("eccb") || a.equals("ecccms")){
				eccBloom=Parse.parseBoolean(b);
			}else if(a.equals("exactkmercounts")){
				forceExactKmerCounts=Parse.parseBoolean(b);
			}else if(a.equals("bits") || a.equals("bloombits") || a.equals("cmsbits")){
				bloomBits=Integer.parseInt(b);
			}else if(a.equals("testmerge") || a.equals("testjunctions")){
				testMerge=Parse.parseBoolean(b);
			}else if(a.equals("testmergewidth")){
				testMergeWidth=Integer.parseInt(b);
			}else if(a.equals("testmergethresh")){
				testMergeThresh=Integer.parseInt(b);
			}else if(a.equals("testmergemult")){
				testMergeMult=Parse.parseKMG(b);
			}else if(a.equals("hashes") || a.equals("bloomhashes")){
				bloomHashes=Integer.parseInt(b);
			}else if(a.equals("shave") || a.equals("removedeadends")){
				shave=Parse.parseBoolean(b);
			}else if(a.equals("rinse") || a.equals("shampoo") || a.equals("removebubbles")){
				rinse=Parse.parseBoolean(b);
			}else if(a.equals("branchlower") || a.equals("branchlowerconst")){
				branchLowerConst=Parse.parseIntKMG(b);
			}else if(a.equals("branchmult2")){
				branchMult2=Parse.parseIntKMG(b);
			}else if(a.equals("branchmult1")){
				branchMult1=Parse.parseIntKMG(b);
			}else if(a.equals("mincount") || a.equals("mincov") || a.equals("mindepth") || a.equals("min")){
				minCountSeed=minCountExtend=Parse.parseIntKMG(b);
			}else if(a.equals("mindepthseed") || a.equals("mds") || a.equals("mincountseed") || a.equals("mcs")){
				minCountSeed=Parse.parseIntKMG(b);
			}else if(a.equals("mindepthextend") || a.equals("mde") || a.equals("mincountextend") || a.equals("mce")){
				minCountExtend=Parse.parseIntKMG(b);
			}else if(a.equals("ilb") || a.equals("ignoreleftbranches") || a.equals("ignoreleftjunctions") || a.equals("ibb") || a.equals("ignorebackbranches")){
				extendThroughLeftJunctions=Parse.parseBoolean(b);
			}else if(a.equals("rem") || a.equals("requireextensionmatch")){
				requireExtensionMatch=Parse.parseBoolean(b);
			}else if(a.equals("rsem") || a.equals("requirestrictextensionmatch")){
				requireStrictExtensionMatch=Parse.parseBoolean(b);
			}else if(a.equals("minapproxoverlaprem") || a.equals("minapproxoverlap")){
				minApproxOverlapRem=Parse.parseIntKMG(b);
			}
			
			else if(a.equals("quantize") || a.equals("quantizesticky")){
				quantizeQuality=Quantizer.parse(arg, a, b);
			}
			
			else if(a.equals("besteffort") || a.equals("forcemerge")){
				MAX_RATIO=0.5f;
				RATIO_MARGIN=1.1f;
				MAX_MISMATCHES_R=500;
				useEfilter=false;
				pfilterRatio=0;
			}
			
			else if(KmerTableSet.isValidArgument(a)){
				//Do nothing
			}
			
			else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(adapterList1!=null || adapterList2!=null){
			verifyAdapters=true;
		}
		
		if(requireStrictExtensionMatch){
			requireExtensionMatch=true;
		}

		if(requireExtensionMatch && extendRight2<1){
			outstream.println("Extend2 is defaulting to 50 because it was unset but r"+(requireStrictExtensionMatch ? "s" : "")+"em mode is being used.");
			extendRight2=50;
		}
		
		if(ecco && !setMix){
			outstream.println("Mergable and unmergable reads are all being sent to out because ecco=t.  To override this, set mix=f.");
			MIX_BAD_AND_GOOD=true;
		}
		
//		assert(!requireExtensionMatch || extendRight2>0) : "The requireExtensionMatch flag requires also setting the extend2 flag.\n"
//				+ "Suggested values are extend2=50 k=62.";
		
//		assert(false) : ecco;
		minInsert=Tools.max(minInsert, MIN_OVERLAPPING_BASES);
		if(minInsert0<1){
			minInsert0=(Tools.max((int)Math.ceil(minInsert*0.75), 5, MIN_OVERLAPPING_BASES_0));
			int cap=(loose ? 50 : 35);
			minInsert0=Tools.min(cap, minInsert0);
		}
		minInsert0=Tools.min(minInsert, minInsert0);
		
		if(MATE_BY_OVERLAP && !useFlatMode && !useRatioMode){
			outstream.println("\n*** WARNING! Both flat and ratio mode were disabled; using ratio mode. ***\n");
			useRatioMode=true;
		}
		
		loglog=(outCardinality==null && !parser.loglog ? null : CardinalityTracker.makeTracker(2048, 31, -1, 0));
		
		{//Process parser fields
			Parser.processQuality();
			
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=(parser.trimq2!=null ? parser.trimq2 : new float[] {parser.trimq});
			trimE=parser.trimE2();
			qtrim1=parser.qtrim1;
			qtrim2=(parser.qtrim2 || (parser.trimq2!=null && parser.trimq2.length>1));
			if(qtrim1==false && qtrim2==false){
				qtrim1=((qtrimLeft||qtrimRight)&&trimq[0]>=0);
			}
			minAvgQuality=parser.minAvgQuality;
			minAvgQualityBases=parser.minAvgQualityBases;
			minReadLength=Tools.max(1, parser.minReadLength);
			maxReadLength=(parser.maxReadLength<0 ? Integer.MAX_VALUE : parser.maxReadLength);
//			untrim=parser.untrim;
			
			forceTrimModulo=parser.forceTrimModulo;
			forceTrimLeft=parser.forceTrimLeft;
			forceTrimRight=parser.forceTrimRight;
			forceTrimRight2=parser.forceTrimRight2;
		}
//		parseCustom=FASTQ.PARSE_CUSTOM;
//		if(parseCustom){FASTQ.PARSE_CUSTOM_WARNING=false;}
		if(verbose){
//			assert(false) : "verbose flag is static final; recompile to change it.";
//			BBMergeOverlapper.verbose=true;
		}
		
		if(trimByOverlap){
			join=false;
		}
		
		if(!mm0set){
			MAX_MISMATCHES0=MAX_MISMATCHES+(loose ? 2 : 0);
		}
		
		if(MAX_MISMATCHES0<MAX_MISMATCHES){
			MAX_MISMATCHES0=MAX_MISMATCHES+(loose ? 2 : 0);
			outstream.println("MAX_MISMATCHES0 was set to "+MAX_MISMATCHES0+" to remain >=MAX_MISMATCHES");
		}
		
		if(MISMATCH_MARGIN>MAX_MISMATCHES){
			MISMATCH_MARGIN=MAX_MISMATCHES;
			outstream.println("MISMATCH_MARGIN was set to "+MISMATCH_MARGIN+" to remain >=MAX_MISMATCHES");
		}
		
		if(recalibrateQuality){CalcTrueQuality.initializeMatrices();}
		
		if(in2==null && in1!=null && in1.contains("#") && !new File(in1).exists()){
			in2=in1.replaceFirst("#", "2");
			in1=in1.replaceFirst("#", "1");
		}
		
		if(out2==null && out1!=null && out1.contains("#")){
			out2=out1.replaceFirst("#", "2");
			out1=out1.replaceFirst("#", "1");
		}
		
		if(outb2==null && outb1!=null && outb1.contains("#")){
			outb2=outb1.replaceFirst("#", "2");
			outb1=outb1.replaceFirst("#", "1");
		}
		
		if(extendRight1>0 || extendRight2>0 || useKFilter || eccTadpole || forceExactKmerCounts){
			
			final long mem=Shared.memAvailable();
			if(Shared.memAvailable()<2000000000L){
				System.err.println("\n***** WARNING *****\nUsing kmer counts uses a lot of memory, but only "+(mem/(1024*1024))+"MB is available.");
				System.err.println("If this process crashes, run bbmerge-auto.sh instead of bbmerge.sh, or set the -Xmx flag.\n");
			}
			
			ArrayList<String> list=new ArrayList<String>();
			{
				StringBuilder sb=new StringBuilder("in=");
				sb.append(in1);
				if(extra.size()>0){
					for(String s : extra){
						sb.append(',').append(s);
					}
				}
				list.add(sb.toString());
			}
			if(in2!=null){list.add("in2="+in2);}
			list.add("branchlower="+branchLowerConst);
			list.add("branchmult1="+branchMult1);
			list.add("branchmult2="+branchMult2);
			list.add("mincountseed="+minCountSeed);
			list.add("mincountextend="+minCountExtend);
			list.add("minprob="+minProb);
			list.add("k="+kmerLength);
			list.add("prealloc="+prealloc);
			list.add("prefilter="+prefilter);
			if(maxReads>0 && maxReads<Long.MAX_VALUE){
				list.add("reads="+maxReads);
			}

			list.add("ecctail="+eccTail);
			list.add("eccpincer="+eccPincer);
			list.add("eccreassemble="+eccReassemble);
			if(filterMemoryOverride>0){list.add("filtermemory="+filterMemoryOverride);}
			
			tadpole=Tadpole.makeTadpole(list.toArray(new String[0]), false);
			bloomFilter=null;
			corrector=null;
		}else if(eccBloom || testMerge){
			assert(kmerLength<32) : "Bloom filter kmer length must be in the range of 1-31.";
			Timer t=new Timer(outstream, true);
			KmerCountAbstract.CANONICAL=true;
			tadpole=null;
			bloomFilter=new BloomFilter(in1, in2, extra, kmerLength, kmerLength, bloomBits, bloomHashes, 1,
					true, false, false, 0.9f);
			corrector=new BloomFilterCorrector(bloomFilter, kmerLength, kmerLength);
			t.stop("Filter creation: \t\t");
			outstream.println(bloomFilter.filter.toShortString());
		}else{
			tadpole=null;
			bloomFilter=null;
			corrector=null;
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outb1, outb2, outinsert, ihist, outCardinality, outAdapter)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					out1+", "+out2+", "+outb1+", "+outb2+", "+outinsert+", "+ihist+"\n");
		}		
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outb1, outb2, outinsert, ihist)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		if(in2!=null){
			assert(!in1.equalsIgnoreCase(in2));
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
		}else{
			FASTQ.TEST_INTERLEAVED=true;
			FASTQ.FORCE_INTERLEAVED=true;
		}
		
		if(THREADS<1){THREADS=Shared.threads();}
		
		useMEEfilter=maxExpectedErrors>0;
		
		Read.VALIDATE_IN_CONSTRUCTOR=(THREADS<16);
		if(netFile==null) {net0=null;}
		else if(!netFile.equals(defaultNetFile)) {net0=CellNetParser.load(netFile);}
	}
	
	void process(){
		Timer ttotal=new Timer();
		ttotal.start();
		
		if(tadpole!=null){
			if(extendRight1>0 || extendRight2>0 || eccTadpole) {
				MAKE_VECTOR=BBMergeOverlapper.MAKE_VECTOR=false;//Currently incompatible...
			}
			
			Timer tload=new Timer();
			Tadpole.showSpeed=false;
			KmerTableSet.showSpeed=false;
			long kmers=tadpole.loadKmers(tload);
			tload.stop();
			outstream.println();
			
			if(shave || rinse){
				tload.start();
				long removed=tadpole.shaveAndRinse(tload, shave, rinse, true);
				tload.stop();
				outstream.println();
			}
			
//			outstream.println("Loaded "+kmers+" kmers in "+tload);
		}
		
		runPhase(join, maxReads, false);
		
		double stdev=0;
		if(histTotal!=null){
			stdev=Tools.standardDeviationHistogram(histTotal);
		}
		
		final long sum=correctCountTotal+incorrectCountTotal;
		final double divp=100d/readsProcessedTotal;
		final double div2=100d/sum;
		
		writeHistogram(ihist, sum*divp);
		
		if(outAdapter!=null){
			assert(findAdapterSequence);
			long count=atrack.writeAdapterConsensus(outAdapter);
			outstream.println("Adapters counted: "+count);
		}
		
		if(outCardinality!=null){
			ReadWrite.writeString(loglog.cardinality()+"\n", outCardinality);
		}
		
		ttotal.stop();
		outstream.println("Total time: "+ttotal+"\n");
		
		outstream.println("Pairs:               \t"+readsProcessedTotal);
		outstream.println("Joined:              \t"+sum+String.format((sum<10000 ? "       " : "   ")+"\t%.3f%%", sum*divp));
		outstream.println("Ambiguous:           \t"+ambiguousCountTotal+String.format((ambiguousCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", ambiguousCountTotal*divp));
		outstream.println("No Solution:         \t"+noSolutionCountTotal+String.format((noSolutionCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", noSolutionCountTotal*divp));
		if(minInsert>0){outstream.println("Too Short:           \t"+tooShortCountTotal+String.format((tooShortCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", tooShortCountTotal*divp));}
		if(maxReadLength<Integer.MAX_VALUE){outstream.println("Too Long:            \t"+tooLongCountTotal+String.format((tooLongCountTotal<10000 ? "       " : "   ")+"\t%.3f%%", tooLongCountTotal*divp));}
		
		if(extendRight1>0 || extendRight2>0){
			double dive=100d/extensionsAttempted;
			outstream.println("Fully Extended:      \t"+fullyExtendedTotal+String.format((fullyExtendedTotal<10000 ? "       " : "   ")+"\t%.3f%%", fullyExtendedTotal*dive));
			outstream.println("Partly Extended:     \t"+partlyExtendedTotal+String.format((partlyExtendedTotal<10000 ? "       " : "   ")+"\t%.3f%%", partlyExtendedTotal*dive));
			outstream.println("Not Extended:        \t"+notExtendedTotal+String.format((notExtendedTotal<10000 ? "       " : "   ")+"\t%.3f%%", notExtendedTotal*dive));
		}
		
		if(adapterList1!=null || adapterList2!=null){
			double diva=100d/(readsProcessedTotal*2);
			outstream.println("Adapters Expected:   \t"+adaptersExpected+String.format((adaptersExpected<10000 ? "       " : "   ")+"\t%.3f%%", adaptersExpected*diva));
			outstream.println("Adapters Found:      \t"+adaptersFound+String.format((adaptersFound<10000 ? "       " : "   ")+"\t%.3f%%", adaptersFound*diva));
		}
		
		if(ecco){
			outstream.println("Errors Corrected:    \t"+errorsCorrectedTotal);
		}
		
		if(loglog!=null){
			outstream.println("Unique "+loglog.k+"-mers:         \t"+loglog.cardinality());
		}
		
		if(parseCustom){
			outstream.println();
			outstream.println("Correct:             \t"+correctCountTotal+String.format((correctCountTotal<10000 ? "       " : "   ")+"\t%.4f%%", correctCountTotal*divp)+Tools.format("  \t%.4f%% of merged", correctCountTotal*div2));
			outstream.println("Incorrect:           \t"+incorrectCountTotal+String.format((incorrectCountTotal<10000 ? "       " : "   ")+"\t%.5f%%", incorrectCountTotal*divp)+Tools.format("  \t%.5f%% of merged", incorrectCountTotal*div2));
			double FPRate=incorrectCountTotal/(double)readsProcessedTotal;
			double FNRate=unmergedOverlappingTotal/(double)readsProcessedTotal;
			double snr=Tools.max(correctCountTotal, 0.001)/(Tools.max(incorrectCountTotal, 0.001));
			double snrDB=Tools.mid(-20, 80, 10*Math.log10(snr));
			outstream.println("SNR:                 \t"+Tools.format("%.3f dB", snrDB));
			outstream.println("FPR:                 \t"+Tools.format("%.4f%%", FPRate*100));
			outstream.println("FNR:                 \t"+Tools.format("%.4f%%", FNRate*100));
			outstream.println();
			outstream.println("Avg Insert Correct:  \t"+Tools.format("%.1f", (insertSumCorrectTotal)*1d/(correctCountTotal)));
			outstream.println("Avg Insert Incorrect:\t"+Tools.format("%.1f", (insertSumIncorrectTotal)*1d/(incorrectCountTotal)));
		}
		
		outstream.println("\nAvg Insert:          \t"+Tools.format("%.1f", (insertSumCorrectTotal+insertSumIncorrectTotal)*1d/(correctCountTotal+incorrectCountTotal)));
		outstream.println("Standard Deviation:  \t"+Tools.format("%.1f", stdev));
		outstream.println("Mode:                \t"+Tools.calcModeHistogram(histTotal));
		
		outstream.println();
		outstream.println("Insert range:        \t"+insertMinTotal+" - "+insertMaxTotal);
		outstream.println("90th percentile:     \t"+Tools.percentileHistogram(histTotal, .9));
		outstream.println("75th percentile:     \t"+Tools.percentileHistogram(histTotal, .75));
		outstream.println("50th percentile:     \t"+Tools.percentileHistogram(histTotal, .5));
		outstream.println("25th percentile:     \t"+Tools.percentileHistogram(histTotal, .25));
		outstream.println("10th percentile:     \t"+Tools.percentileHistogram(histTotal, .1));
	}
	
	public static void writeHistogram(String fname, double percentMerged){
		if(fname==null){return;}
		StringBuilder sb=new StringBuilder();

		if(showHistStats){
			sb.append("#Mean\t"+Tools.format("%.3f", Tools.averageHistogram(histTotal))+"\n");
			sb.append("#Median\t"+Tools.percentileHistogram(histTotal, 0.5)+"\n");
			sb.append("#Mode\t"+Tools.calcModeHistogram(histTotal)+"\n");
			sb.append("#STDev\t"+Tools.format("%.3f", Tools.standardDeviationHistogram(histTotal))+"\n");
			sb.append("#PercentOfPairs\t"+Tools.format("%.3f", percentMerged)+"\n");
		}
		sb.append("#InsertSize\tCount\n");
		for(int i=0; i<histTotal.length && i<=insertMaxTotal; i+=bin){
			int x=0;
			int y=0;
			for(int j=i; j<i+bin && j<histTotal.length; j++){
				x+=histTotal[j];
				y++;
			}
			x=(x+bin-1)/y;
			if(x>0 || !NONZERO_ONLY){
				sb.append(i+"\t"+x+"\n");
			}
		}
		ReadWrite.writeStringInThread(sb, fname);
	}
	
//	private static String toAdapterSequence(LongList[] lists, boolean trimPolyA){
//		StringBuilder adapter=new StringBuilder();
//		long max=0;
//		int lastBase=-1;
//		for(int i=0; true; i++){
//			long a=lists[0].get(i);
//			long c=lists[1].get(i);
//			long g=lists[2].get(i);
//			long t=lists[3].get(i);
//			long sum=(a+c+g+t);
//			max=Tools.max(max, sum);
//			if(sum==0 || (sum<10 && sum<=max/1000) || (max>100 && sum<8)){break;}
//			long thresh=(max>100 ? 4+(sum*2)/3 : (sum*2)/3);
//			if(a>thresh){
//				adapter.append('A');
//				lastBase=i;
//			}else if(c>thresh){
//				adapter.append('C');
//				lastBase=i;
//			}else if(g>thresh){
//				adapter.append('G');
//				lastBase=i;
//			}else if(t>thresh){
//				adapter.append('T');
//				lastBase=i;
//			}else{
//				adapter.append('N');
//			}
//		}
//		if(lastBase<0){return "N";}
//
//		String trimmed=trimPoly(adapter.toString(), 'N');
//		if(trimPolyA){
//			trimmed=trimPoly(trimmed, 'G');
//			trimmed=trimPoly(trimmed, 'A');
//		}
//		
//		if(lastBase>=0){
//			char A=(trimPolyA ? 'A' : 'N');
//			while(lastBase>=0 && (adapter.charAt(lastBase)=='N' || adapter.charAt(lastBase)==A)){lastBase--;}
//		}
//		
//		if(lastBase<0){return "N";}
//		return adapter.substring(0, lastBase+1);
//	}
//	
//	private static String trimPoly(String adapter, char trim){
//		int lastBase=-1;
//		for(int i=0; i<adapter.length(); i++){
//			char c=adapter.charAt(i);
//			if(AminoAcid.isFullyDefined(c)){
//				lastBase=i;
//			}
//		}
//		
//		int aCount=0;
//		int nCount=0;
//		int count=0;
//		if(lastBase>=0){
//			while(lastBase>=0){
//				char c=adapter.charAt(lastBase);
//				if(c=='N'){nCount++;}
//				else if(c==trim){aCount++;}
//				else{break;}
//				count++;
//				lastBase--;
//			}
//		}
//		
//		if(lastBase<0){return "N";}
//		if(count==nCount || (aCount>3)){
//			return adapter.substring(0, lastBase+1);
//		}
//		return adapter;
//	}

	
	public void runPhase(boolean join, long maxReads, boolean perfectonly){
		
		Timer talign=new Timer();
		
		ConcurrentReadOutputStream rosgood=null;
		ConcurrentReadOutputStream rosbad=null;
		ConcurrentReadOutputStream rosinsert=null;
		
		if(out1!=null){
			if(join==true){
				if(out2==null){outstream.println("Writing mergable reads merged.");}
				else{
					outstream.println("WARNING: 2 output files specified even though 'merge=true'.  out2 will be ignored.");
					out2=null;
				}
			}else{
				if(out2==null){outstream.println("Writing mergable reads interleaved.");}
				else{outstream.println("Writing mergable reads unmerged in two files.");}
			}
			
			final FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			final FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(16, 2*THREADS);
			rosgood=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			rosgood.start();
		}
		
		if(outb1!=null){

			final FileFormat ff1=FileFormat.testOutput(outb1, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			final FileFormat ff2=FileFormat.testOutput(outb2, FileFormat.FASTQ, null, true, overwrite, append, ordered);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";
			
			final int buff=Tools.max(16, 2*THREADS);
			rosbad=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
			rosbad.start();
		}
		
		if(outinsert!=null){
			final int buff=Tools.max(16, 2*THREADS);
			
			String out1=outinsert.replaceFirst("#", "1");

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			
			ReadStreamWriter.HEADER=header();
			final FileFormat ff=FileFormat.testOutput(out1, FileFormat.ATTACHMENT, ".info", true, overwrite, append, ordered);
			rosinsert=ConcurrentReadOutputStream.getStream(ff, null, null, null, buff, null, false);
			rosinsert.start();
		}
		
		
		if(rosgood!=null || rosbad!=null || rosinsert!=null){
			outstream.println("Started output threads.");
		}
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.setSampleRate(samplerate, sampleseed);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
//		assert(paired);//Fails on empty files.
		if(verbose){outstream.println("Paired: "+paired);}
		
		talign.start();
		
		
		MateThread[] pta=new MateThread[THREADS];
		for(int i=0; i<pta.length; i++){
			pta[i]=new MateThread(cris, rosgood, rosbad, rosinsert, join, trimByOverlap);
			pta[i].start();
		}

		resetCounters();
		
		for(int i=0; i<pta.length; i++){
			MateThread ct=pta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}

				readsProcessedTotal+=ct.pairsProcessed;
				basesProcessedTotal+=ct.basesProcessed;
				matedCountTotal+=ct.matedCount;
				correctCountTotal+=ct.correctCount;
				unmergedOverlappingTotal+=ct.unmergedOverlappingCount;
				ambiguousCountTotal+=ct.ambiguousCount;
				tooShortCountTotal+=ct.tooShortCount;
				tooLongCountTotal+=ct.tooLongCount;
				incorrectCountTotal+=ct.incorrectCount;
				noSolutionCountTotal+=ct.noSolutionCount;
				insertSumCorrectTotal+=ct.insertSumCorrect;
				insertSumIncorrectTotal+=ct.insertSumIncorrect;
				
				errorsCorrectedTotal+=ct.errorsCorrectedT;
				
				fullyExtendedTotal+=ct.fullyExtendedT;
				partlyExtendedTotal+=ct.partlyExtendedT;
				notExtendedTotal+=ct.notExtendedT;
				extensionsAttempted+=ct.extensionsAttemptedT;

				adaptersExpected+=ct.adaptersExpectedT;
				adaptersFound+=ct.adaptersFoundT;

				insertMinTotal=Tools.min(ct.insertMin, insertMinTotal);
				insertMaxTotal=Tools.max(ct.insertMax, insertMaxTotal);
				
//				outstream.println(ct.insertMin+", "+ct.insertMax);
				
				if(ct.hist!=null){
					for(int h=0; h<ct.hist.length; h++){
						histTotal[h]+=ct.hist[h];
					}
				}
				
				if(findAdapterSequence){
					atrack.merge(ct.atrackT);
				}
			}
		}
		
//		outstream.println("Finished reading");
		errorState|=ReadWrite.closeStreams(cris, rosgood, rosbad, rosinsert);
		
		talign.stop();
//		outstream.println("Align time: "+talign);
	}
	
	public static final float mergeableFraction(String fname1, String fname2, long numReads, float samplerate){
		long[] hist=makeInsertHistogram(fname1, fname2, numReads, samplerate);
		if(hist==null || hist.length<2){return 0;}
		long sum=shared.Vector.sum(hist);
		return sum<1 ? 0 : (sum-hist[0])/(float)sum;
	}
	
	public static final long[] makeInsertHistogram(String fname1, String fname2, long numReads, float samplerate){
		assert(fname1!=null);
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(fname1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(fname2, FileFormat.FASTQ, null, true, true);
			if(ff1.stdio()){return null;}
			assert(!ff1.stdio()) : "Standard in is not allowed as input when calculating insert size distributions for files.";
			cris=ConcurrentReadInputStream.getReadInputStream(numReads, true, ff1, ff2);
			cris.setSampleRate(samplerate, 1);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
			if(!cris.paired()){
				ReadWrite.closeStreams(cris);
				return null;
			}
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert(r.mate!=null);
		}

		LongList ll=new LongList(500);
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			for(Read r1 : reads){
				int x=findOverlapLoose(r1, r1.mate, false);
				if(x>0){ll.increment(x, 1);}
				else{ll.increment(0, 1);}
			}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);
		return ll.toArray();
	}

	/** Returns the insert size as calculated by overlap, or -1 */
	public static final int findOverlapUStrict(final Read r1, final Read r2, boolean ecc){
		final float maxRatio=0.045f;
		final float ratioMargin=12f;
		final float ratioOffset=0.5f;
		
		final float efilterRatio=2f;
		final float efilterOffset=0.03f;
		final float pfilterRatio=0.03f;
		
		final int minOverlap=14;
		final int minOverlap0=3;
		final int minInsert=35;
		final int minInsert0=20;
		final int entropy=56;
		final float minSecondRatio=0.16f;
		
//		final int ratioMinOverlapReduction=4;
		
		final int x=findOverlap(r1, r2, ecc,
				minOverlap, minOverlap0, minInsert, minInsert0, entropy,
				maxRatio, minSecondRatio, ratioMargin, ratioOffset,
				efilterRatio, efilterOffset, pfilterRatio);
		return x;
	}

	/** Returns the insert size as calculated by overlap, or -1 */
	public static final int findOverlapVStrict(final Read r1, final Read r2, boolean ecc){
		final float maxRatio=0.05f;
		final float ratioMargin=12f;
		final float ratioOffset=0.5f;
		
		final float efilterRatio=2f;
		final float efilterOffset=0.05f;
		final float pfilterRatio=0.008f;
		
		final int minOverlap=12;
		final int minOverlap0=4;
		final int minInsert=35;
		final int minInsert0=25;
		final int entropy=52;
		final float minSecondRatio=0.16f;
		
//		final int ratioMinOverlapReduction=4;
		
		final int x=findOverlap(r1, r2, ecc,
				minOverlap, minOverlap0, minInsert, minInsert0, entropy,
				maxRatio, minSecondRatio, ratioMargin, ratioOffset,
				efilterRatio, efilterOffset, pfilterRatio);
		return x;
	}
	
	public static final Read tryToMerge(final Read r1, final Read r2) {
		if(r2==null) {return null;}
		int insert=findOverlapStrict(r1, r2, false);
		if(insert<1){return null;}
		r2.reverseComplement();
		Read r3=r1.joinRead(insert);
		r2.reverseComplement();
		return r3;
	}

	/** Returns the insert size as calculated by overlap, or -1 */
	public static final int findOverlapStrict(final Read r1, final Read r2, boolean ecc){
		final float maxRatio=0.075f;
		final float ratioMargin=7.5f;
		final float ratioOffset=0.55f;
		
		final float efilterRatio=4f;
		final float efilterOffset=0.05f;
		final float pfilterRatio=0.0008f;
		
		final int minOverlap=11;
		final int minOverlap0=5;
		final int minInsert=35;
		final int minInsert0=25;
		final int entropy=42;
		final float minSecondRatio=0.12f;
		
//		final int ratioMinOverlapReduction=4;
		
		final int x=findOverlap(r1, r2, ecc,
				minOverlap, minOverlap0, minInsert, minInsert0, entropy,
				maxRatio, minSecondRatio, ratioMargin, ratioOffset,
				efilterRatio, efilterOffset, pfilterRatio);
		return x;
	}

	/** Returns the insert size as calculated by overlap, or -1 */
	public static final int findOverlapLoose(final Read r1, final Read r2, boolean ecc){
		final float maxRatio=0.11f;
		final float ratioMargin=4.7f;
		final float ratioOffset=0.45f;
		
		final float efilterRatio=8f;
		final float efilterOffset=0.55f;
		final float pfilterRatio=0.00002f;

		final int minOverlap=5;
		final int minOverlap0=6;
		final int minInsert=16;
		final int minInsert0=16;
		final int entropy=30;
		final float minSecondRatio=0.1f;
		
//		final int ratioMinOverlapReduction=2;
		
		final int x=findOverlap(r1, r2, ecc,
				minOverlap, minOverlap0, minInsert, minInsert0, entropy,
				maxRatio, minSecondRatio, ratioMargin, ratioOffset,
				efilterRatio, efilterOffset, pfilterRatio);
		return x;
	}
	
	/** Returns the insert size as calculated by overlap, or -1 */
	public static final int findOverlap(final Read r1, final Read r2, final boolean ecc,
			int minOverlap, final int minOverlap0, final int minInsert, final int minInsert0, final int entropy,
			final float maxRatio, final float minSecondRatio, final float ratioMargin, final float ratioOffset,
			final float efilterRatio, final float efilterOffset, final float pfilterRatio){
		
		assert(r1!=null && r2!=null);
		if(!r1.validated()){r1.validate(true);}
		if(!r2.validated()){r2.validate(true);}

		final int len1=r1.length(), len2=r2.length();
		final int minlen=Tools.min(len1, len2);
		
		if(minlen<MIN_OVERLAPPING_BASES || minlen<minInsert){
			return -1;
		}
		
		int[] rvector=localRvector.get();
		if(rvector==null){
			rvector=KillSwitch.allocInt1D(5);
			localRvector.set(rvector);
		}
		
		r2.reverseComplement();
		
		int bestInsert=-1;
		int bestBad=999999;
		boolean ambig, tooShort=false;
		
		if(USE_MAPPING && r1.chrom==r2.chrom && r1.start<r1.stop && r1.mapped() && r2.mapped()){
			bestBad=0;
			bestInsert=Read.insertSizeMapped(r1, r2, false);
			ambig=false;
		}else{
			if(entropy>0){
				int a=BBMergeOverlapper.calcMinOverlapByEntropy(r1.bases, 3, null, entropy);
				int b=BBMergeOverlapper.calcMinOverlapByEntropy(r2.bases, 3, null, entropy);
				minOverlap=Tools.max(MIN_OVERLAPPING_BASES, Tools.max(a, b));
			}else{minOverlap=MIN_OVERLAPPING_BASES;}
			if(verbose){outstream.println("minOverlap: "+minOverlap);}
			
			rvector[4]=0;

			int x=BBMergeOverlapper.mateByOverlapRatio(r1, r2, null, null, rvector, minOverlap0, minOverlap,
					minInsert0, minInsert, maxRatio, minSecondRatio, ratioMargin, ratioOffset, 0.95f, 0.95f, false);
			bestInsert=x;
			bestBad=rvector[2];
			ambig=(x>-1 ? rvector[4]==1 : false);
		}
		
		//TODO:  Crucial!  This line can vastly reduce merge rate, particularly if quality values are inaccurate.
		if(bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null){
			float bestExpected=BBMergeOverlapper.expectedMismatches(r1, r2, bestInsert);
			if((bestExpected+efilterOffset)*efilterRatio<bestBad){ambig=true;}
			if(verbose){outstream.println("Result after efilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
		}
		
		//TODO:  Crucial!  This line can vastly reduce merge rate, particularly if quality values are inaccurate.
		if(pfilterRatio>0 && bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null){
			float probability=BBMergeOverlapper.probability(r1, r2, bestInsert);
			if(probability<pfilterRatio){bestInsert=-1;}
			if(verbose){outstream.println("Result after pfilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
		}
		
		if(staticAdapterList!=null && bestInsert>0 && (bestInsert<len1 || bestInsert<len2)){
			boolean foundAdapter=verifyAdaptersStatic(r1, r2, bestInsert);
			if(!foundAdapter){
				bestInsert=RET_NO_SOLUTION;
				assert(r1.quality==null || (r1.quality.length==r1.bases.length && r2.quality.length==r2.bases.length));
				r2.reverseComplement();
				return RET_NO_SOLUTION;
			}
		}
		
		tooShort=(!ambig && bestInsert>0 && bestInsert<minInsert);
		
		if(ecc && bestInsert>-1 && !ambig && !tooShort){
			int errors=errorCorrectWithInsert(r1, r2, bestInsert);
		}
		
		if(r2!=null){r2.reverseComplement();}
		if(!ambig && bestInsert>-1){r1.setInsert(bestInsert);}
		
		return ambig ? -1 : bestInsert;
	}
	
	public static int errorCorrectWithInsert(Read r1, Read r2, int insert){
		if(insert>=r1.length()+r2.length()){return 0;}
		assert(insert>0);
		Read joined=r1.joinRead(insert);
		
		if(joined!=null && joined.length()>0){
			final int lenj=joined.length();
			final int lim1=Tools.min(joined.length(), r1.length());
			final int lim2=lenj-Tools.min(r2.length(), lenj);
			
			final byte[] old1=r1.bases, old2=r2.bases;
			r1.errors=r2.errors=0;
			
			r1.bases=KillSwitch.copyOfRange(joined.bases, 0, lim1);
			for(int i=0; i<r1.bases.length; i++){//count errors
				if(old1[i]!=r1.bases[i] && AminoAcid.isFullyDefined(r1.bases[i])){r1.errors++;}
			}
			if(r1.quality!=null){
				if(changeQuality){
					r1.quality=KillSwitch.copyOfRange(joined.quality, 0, lim1);
				}else{
					r1.quality=KillSwitch.copyOfRange(r1.quality, 0, lim1);
				}
			}

			r2.bases=KillSwitch.copyOfRange(joined.bases, lim2, lenj);
			for(int i=0, j=(old2.length-r2.bases.length); i<r2.bases.length; i++, j++){//count errors
				if(old2[j]!=r2.bases[i] && AminoAcid.isFullyDefined(r2.bases[i])){r2.errors++;}
			}
			if(r2.quality!=null){
				if(changeQuality){
					r2.quality=KillSwitch.copyOfRange(joined.quality, lim2, lenj);
				}else{
					r2.quality=KillSwitch.copyOfRange(r2.quality, 0, r2.bases.length);
				}
			}
		}else{
			return 0;
		}
		return r1.errors+r2.errors;
	}
	
	public static int countErrors(Read r1, Read r2, Read joined){
		if(joined==null){return 0;}
		final int insert=joined.length();
		if(insert<1 || insert>=r1.length()+r2.length()){return 0;}
		
		final int lenj=joined.length();

		int e1=0, e2=0;
		final byte[] bases1=r1.bases, bases2=r2.bases, basesj=joined.bases;
		for(int i=0; i<bases1.length && i<lenj; i++){//count errors
			if(bases1[i]!=basesj[i] && AminoAcid.isFullyDefined(basesj[i])){e1++;}
		}
		for(int i=0, j=lenj-bases2.length; i<bases2.length && j<lenj; i++, j++){//count errors
			if(j>=0 && bases2[i]!=basesj[j] && AminoAcid.isFullyDefined(basesj[j])){e2++;}
		}
		return e1+e2;
	}

	/** Assumes both reads are in original orientation */
	public static int countErrors(Read r1, Read r2, final int insert){
		return countErrors(r1.bases, r2.bases, insert);
	}
	
	/** Assumes both reads are in original orientation */
	public static int countErrors(byte[] bases1, byte[] bases2, final int insert){
		if(insert<1 || bases1==null || bases2==null){return 0;}
		if(insert>=bases1.length+bases2.length){return 0;}

		int offset1=(insert>=bases1.length ? insert-bases1.length : 0);
		int offset2=(insert>=bases2.length ? bases2.length-1 : insert-1);
		
		int errors=0;
//		int overlap=0;
		for(int i=offset1, j=offset2; i<bases1.length && j>=0; i++, j--) {
			byte a=bases1[i], b=bases2[j];
			errors+=(a==AminoAcid.baseToComplementExtended[b] ? 0 : 1);
//			overlap++;
		}
	
//		catch (Exception e) {
//			System.err.println(e);
//			assert(false) : insert+", "+overlap+", "+
//					errors+", "+bases1.length+", "+offset1+", "+offset2;
		//Should typically be true...  reads that mostly mismatch shouldn't merge
//		assert(errors*2<overlap+1) : insert+", "+overlap+", "+
//			errors+", "+bases1.length+", "+offset1+", "+offset2;
		return errors;
	}

	public static String header(){
		if(TAG_CUSTOM) {
			String i="##mo\tr1ee\tr2ee\talen\tblen\tbi\tbo\tbb\tbg\tbr\tbbi"+
					"\tsbi\tsbo\tsbb\tsbg\tsbr\tsbbi\trr\tbadr\tgoodr\tovr\tbe\tpr";
//			String o="truth\tgoal\tmerged";
//			String o="goal";
			String o="truth";
			final int inputs=Tools.tabPattern.split(i).length;
			final int outputs=Tools.tabPattern.split(o).length;
			return "#dims\t"+inputs+"\t"+outputs+"\n"+i+"\t"+o+"\n";
		}else {
			return "#id\tnumericID\tinsert\tstatus\tmismatches\n";
		}
	}
	
	private static void qtrim(Read r1, Read r2, int iter){
		if(false /*untrim*/){
			TrimRead.trim(r1, qtrimLeft, qtrimRight, trimq[iter], trimE[iter], 1);
			TrimRead.trim(r2, qtrimLeft, qtrimRight, trimq[iter], trimE[iter], 1);
		}else{
			TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq[iter], trimE[iter], 1);
			TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq[iter], trimE[iter], 1);
		}
	}
	
	public static final ArrayList<byte[]> getAdapterList(String name){
		if(name==null){return null;}
		ArrayList<byte[]> alb;
		if(name.equalsIgnoreCase("default") || (name.equalsIgnoreCase("adapters") && !new File(name).exists())){
			alb=new ArrayList<byte[]>(defaultAdapters.length);
			for(String s : defaultAdapters){alb.add(s.getBytes());}
		}else{
			alb=Tools.toAdapterList(name, maxAdapterLength);
		}
		return alb;
	}
	
	private static boolean verifyAdaptersStatic(Read r1, Read r2, int bestInsert){
		
		if(staticAdapterList==null){return true;}
		final int minAdapterOverlap=3;
		final int len1=r1.length(), len2=r2.length();
		if(Tools.max(len1, len2)-minAdapterOverlap<bestInsert){return true;}
		
		int good=0, bad=0, invalid=0;
		final float minAdapterRatio=0.6f;
		{
			final int adapterLen=len1-bestInsert;
			if(good<1 && adapterLen>=minAdapterOverlap){
				if(adapterLen>=minAdapterOverlap){
					for(byte[] adapter : staticAdapterList){
						int aiv=adapterIsValid(r1, adapter, bestInsert, minAdapterOverlap, minAdapterRatio);
						if(aiv==2){good++; break;}
						else if(aiv==1){invalid++;}
						else{bad++;}
					}
				}
			}
		}
		{
			final int adapterLen=len2-bestInsert;
			if(good<1 && adapterLen>=minAdapterOverlap){
				if(adapterLen>=minAdapterOverlap){
					r2.reverseComplement();
					for(byte[] adapter : staticAdapterList){
						int aiv=adapterIsValid(r2, adapter, bestInsert, minAdapterOverlap, minAdapterRatio);
						if(aiv==2){good++; break;}
						else if(aiv==1){invalid++;}
						else{bad++;}
					}
					r2.reverseComplement();
				}
			}
		}
		if(good>0){
			return true;
		}else if(bad>0){
			return false;
		}else if(invalid>0 && strict){
			return false;
		}
		return true;
	}
	
	/** Returns 0=bad, 1=unknown, 2=good */
	private static int adapterIsValid(final Read r, final byte[] adapter, final int insert,
			final int minAdapterOverlap, final float minAdapterRatio){
//		assert(false) : insert+", "+minAdapterOverlap;
//		outstream.println("Comparing adapter.");
		
		final byte[] bases=r.bases;
		final int limit=Tools.min(adapter.length, bases.length-insert);
		final int badLimit=Tools.max(2, limit-(int)(minAdapterRatio*limit));
		
//		assert(false) : limit+", "+minAdapterOverlap+", "+insert+", "+bases.length+", "+adapter.length;
		
//		int shortBadLimit=1+(int)Math.ceil(minAdapterOverlap*(1-minAdapterRatio));
//		outstream.println("limit="+limit);
		if(limit<minAdapterOverlap && limit<=badLimit){
//			assert(false);
			return 1;
		}

//		StringBuilder sb1=new StringBuilder();
//		StringBuilder sb2=new StringBuilder();
		
		int good=0, bad=0, nocall=0;
		for(int i=0, j=insert; i<limit; i++, j++){
			byte a=adapter[i], b=bases[j];
//			sb1.append((char)a);
//			sb2.append((char)b);
			if(a=='N' || b=='N'){
				nocall++;
			}else if(a==b){
				good++;
			}else{
				bad++;
				if(bad>badLimit){return 0;}
			}
		}
//		assert(false) : good+", "+bad+", "+nocall;
		final int overlap=bad+good;
		if(overlap<minAdapterOverlap){
			return 1;
		}
//		final int minGood=(int)minAdapterRatio*overlap;
		final int minGood=(int)Math.ceil(minAdapterRatio*overlap);
		final int ret=(good>=minGood ? 2 : 0);
//		assert(ret==2) : "\n"+sb1+"\n"+sb2+"\n"
//				+ "ratio="+ratio+", good="+good+", bad="+bad+", nocall="+nocall+", limit="+limit+", overlap="+overlap+"\n"
//				+ "(good/overlap)="+good+"/"+overlap+"="+(good/overlap);
//		outstream.println("Returning "+ret+" for ratio "+ratio);
		return ret;
	}
	
	private class MateThread extends Thread{
		
		
		public MateThread(ConcurrentReadInputStream cris_, ConcurrentReadOutputStream rosgood_, ConcurrentReadOutputStream rosbad_, ConcurrentReadOutputStream rosi_,
				boolean joinReads_, boolean trimByOverlap_) {
			cris=cris_;
			rosgood=rosgood_;
			rosbad=rosbad_;
			rosi=rosi_;
			joinReads=joinReads_;
			trimReadsByOverlap=trimByOverlap_;
			kmerT=(tadpole==null || tadpole.k()<32 ? null : new Kmer(tadpole.k()));
			net=(net0==null ? null : net0.copy(false));
			
			if(useEntropy){
				kmerCounts=new short[1<<(2*entropyK)];
			}else{
				kmerCounts=null;
			}
		}
		
		
		@Override
		public void run(){
			processReads();
		}
		
		private void processReads() {
			assert(USE_MAPPING || MATE_BY_OVERLAP);
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(r.mate!=null);
			}
			
			final byte[][] originals=((tadpole!=null || qtrimRight || qtrimLeft) &&
					(rosbad!=null || (rosgood!=null && (!join || MIX_BAD_AND_GOOD)))) ? new byte[4][] : null;
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				ArrayList<Read> listg=(rosgood==null /*&& rosi==null*/ ? null : new ArrayList<Read>(reads.size()));
				ArrayList<Read> listb=(rosbad==null ? null : new ArrayList<Read>(reads.size()));
				
				if(loglog!=null){
					for(Read r1 : reads){loglog.hash(r1);}
				}
				
				for(Read r1 : reads){
					int bestInsert=findOverlapInThread(r1, originals, listg, listb);
				}
				
				if(rosgood!=null){
					if(quantizeQuality) {Quantizer.quantize(listg);}
					rosgood.add(listg, ln.id);
				}
				if(rosi!=null){
					//This prints both merged and unmerged reads
					
					if(TRAINING) {
						for(Read r1 : reads){
							assert(false); //Maybe nothing needs doing here
						}
					}else {
//						for(Read r1 : reads){//Legacy outinsert support
//							StringBuilder sb=new StringBuilder(40);
//							sb.append(r1.id==null ? r1.numericID+"" : r1.id).append('\t');
//							sb.append(r1.numericID).append('\t');
//							final int bestInsert=r1.insert();
//							sb.append(bestInsert>=0 ? bestInsert : -1);
//							sb.append('\t');
//
//							if(bestInsert==RET_NO_SOLUTION){sb.append('F');}//Failed
//							else if(bestInsert==RET_AMBIG){sb.append('A');} //Ambiguous
//							else if(bestInsert==RET_SHORT){sb.append('S');} //Short
//							else{
//								if(r1.errors>0){sb.append('I');}//Imperfect
//								else{sb.append('P');}//Perfect
//								sb.append('\t');
//								sb.append(r1.errors);
//							}
//							r1.obj=sb;
//						}
					}
//					for(Read r : reads) {
//						ByteBuilder bb=(ByteBuilder) r.obj;
//						if(bb.length()<150 || bb.len
//					}
					rosi.add(reads, ln.id);
				}
				if(rosbad!=null){
					if(quantizeQuality) {Quantizer.quantize(listb);}
					rosbad.add(listb, ln.id);
				}
				
				//			outstream.println("returning list");
				cris.returnList(ln);
				//			outstream.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
				//			outstream.println("reads: "+(reads==null ? "null" : reads.size()));
			}
			cris.returnList(ln);
		}
		
		private int findOverlapInThread(final Read r1, final byte[][] originals, ArrayList<Read> listg, ArrayList<Read> listb){

			final Read r2=r1.mate;
			final int trueSize=parseCustom ? parseInsert(r1, r2) : -1;
			final int initialLen1=r1.length(), initialLen2=r2.length();
			final int minReadLen=Tools.min(initialLen1, initialLen2), lenSum=initialLen1+initialLen2;
			final int trueOverlap=Tools.min(trueSize, minReadLen, lenSum-trueSize);
			
			if(originals!=null){
				if(eccTadpole && false){
					originals[0]=r1.bases;
					originals[1]=r1.quality;
					originals[2]=r2.bases;
					originals[3]=r2.quality;
				}else{
					originals[0]=(r1.bases==null ? null : r1.bases.clone());
					originals[1]=(r1.quality==null ? null : r1.quality.clone());
					originals[2]=(r2.bases==null ? null : r2.bases.clone());
					originals[3]=(r2.quality==null ? null : r2.quality.clone());
				}
			}
			assert(r1.quality==null || (r1.quality.length==r1.bases.length && r2.quality.length==r2.bases.length));
			
			final int bestInsert=processReadPair(r1, r2);
			assert(r1.quality==null || (r1.quality.length==r1.bases.length && r2.quality.length==r2.bases.length));
			if(verbose){outstream.println("trueSize="+trueSize+", bestInsert="+bestInsert);}
			
			Read joined=null;
			
//			System.err.println(r2.obj);
			
			if(ecco){
				if(r1.length()>initialLen1){
					TrimRead.trimByAmount(r1, 0, r1.length()-initialLen1, 1);
				}
				if(r2.length()>initialLen2){
					TrimRead.trimByAmount(r2, 0, r2.length()-initialLen2, 1);
				}
			}
			
			if(bestInsert>0){//Decided to merge
				if(trueOverlap>0 && trueOverlap<=8 && EXEMPT_SHORT_OVERLAPS_FROM_STATS) {
					//Temporarily ignore for purpose of training
				}else 
					if(bestInsert==trueSize){
					correctCount++;
					insertSumCorrect+=bestInsert;
				}else{
					incorrectCount++;
					insertSumIncorrect+=bestInsert;
				}
				
				r1.setInsert(bestInsert);
				insertMin=Tools.min(bestInsert, insertMin);
				insertMax=Tools.max(bestInsert, insertMax);
				hist[Tools.min(bestInsert, hist.length-1)]++;
				if(joinReads){
					r2.reverseComplement();
					joined=r1.joinRead(bestInsert);
					r2.reverseComplement();
					assert(joined.length()==bestInsert);
					
					if(trimNonOverlapping){
						final int maxLen=Tools.max(initialLen1, initialLen2);
						if(bestInsert>=maxLen){//left trim originals, rl trim fused
							final int overlap=initialLen1+initialLen2-bestInsert;
							TrimRead.trimByAmount(joined, r1.length()-overlap, r2.length()-overlap, 1);
						}else{//right trim originals, no trim fused
							//Do nothing
						}
					}
					
				}else if(ecco){
					r2.reverseComplement();
					errorCorrectWithInsert(r1, r2, bestInsert);
					r2.reverseComplement();
					errorsCorrectedT+=(r1.errors+r2.errors);
					
					if(trimNonOverlapping){
						final int maxLen=Tools.max(initialLen1, initialLen2);
						if(bestInsert>=maxLen){//left trim originals, rl trim fused
							final int overlap=initialLen1+initialLen2-bestInsert;
							TrimRead.trimByAmount(r1, r1.length()-overlap, 0, 1);
							TrimRead.trimByAmount(r2, r2.length()-overlap, 0, 1);
						}else{//right trim originals, no trim fused
							final int overlap=bestInsert;
							if(r1.length()==overlap && r2.length()==overlap){
								//Do nothing, assuming they've already been trimmed.
							}else{
								TrimRead.trimByAmount(r1, 0, r1.length()-bestInsert, 1);
								TrimRead.trimByAmount(r2, 0, r2.length()-bestInsert, 1);
							}
						}
					}
					
				}
			}else{
				if(bestInsert==RET_AMBIG){ambiguousCount++;}
				else if(bestInsert==RET_SHORT){tooShortCount++;}
				else if(bestInsert==RET_LONG){tooLongCount++;}
				else if(bestInsert==RET_NO_SOLUTION){noSolutionCount++;}
				
				if(trueSize<initialLen1+initialLen2 && trueSize>0) {unmergedOverlappingCount++;}
			}
			
			r1.setInsert(bestInsert);
			
			if(findAdapterSequence && bestInsert>0 && bestInsert<r1.length()){
				storeAdapterSequence(r1, bestInsert);
				//r2.reverseComplement();
				storeAdapterSequence(r2, bestInsert);
				//outstream.println(new String(r2.bases));
				//r2.reverseComplement();
			}
			
			if(originals!=null && (!ecco || bestInsert<1)){
				r1.bases=originals[0];
				r1.quality=originals[1];
				r2.bases=originals[2];
				r2.quality=originals[3];
			}
			
			if(trimReadsByOverlap && bestInsert>0){
				int trimLim=bestInsert-1;
				if(trimLim<r1.length()){
					if(verbose){outstream.println("Overlap right trimming r1 to "+0+", "+(trimLim));}
					int x=TrimRead.trimToPosition(r1, 0, trimLim, 1);
					if(verbose){outstream.println("Trimmed "+x+" bases: "+new String(r1.bases));}
				}
				if(trimLim<r2.length()){
					if(verbose){outstream.println("Overlap right trimming r2 to "+0+", "+(trimLim));}
					int x=TrimRead.trimToPosition(r2, 0, trimLim, 1);
					if(verbose){outstream.println("Trimmed "+x+" bases: "+new String(r2.bases));}
				}
			}
			
			if(parseCustom) {

				final boolean fp=(bestInsert>0 && bestInsert!=trueSize);
				final boolean fn=(bestInsert<1 && trueOverlap>0);
				final boolean tp=(bestInsert>0 && bestInsert==trueSize);
				final boolean tn=(bestInsert<1 && trueOverlap<1);
				
				assert((fp ? 1 : 0)+(fn ? 1 : 0)+(tp ? 1 : 0)+(tn ? 1 : 0)==1) : 
					fp+", "+fn+", "+tp+", "+tn+", "+bestInsert+", "+trueOverlap+", "+trueSize+", "+(initialLen1+initialLen2);
				final boolean correct=(bestInsert==trueSize);
				if(MAKE_VECTOR) {
					FloatList fl=r2.fetchVector();
					if(fl.size==23) {fl.add((tp || fn) ? 1 : 0);}
//					assert(false) : fl;
				}
				if(TAG_CUSTOM) {
					ByteBuilder bb=r1.fetchBB();
					//				bb.append("_tr").append(correct ? 1 : 0);

//					if(correct) {//truth
//						bb.append(1);
//					}else {
//						bb.append(0);
//					}
					if(tp || fn) {//truth
						bb.append(1);
					}else {
						bb.append(0);
					}

					if(trueOverlap>0 && trueOverlap<8 && EXEMPT_SHORT_OVERLAPS_FROM_STATS) {
						r1.nullifyObj();
					}
				}

				if(ONLY_OUTPUT_INCORRECT && correct){
					if(TAG_CUSTOM) {r1.nullifyObj();}
					return bestInsert;
				}
				
//				false, 231, 58, 231, 289
				assert(fn==(trueSize<initialLen1+initialLen2 && trueSize>0 && bestInsert<=0)) : 
					fn+", "+bestInsert+", "+trueOverlap+", "+trueSize+", "+(initialLen1+initialLen2);
				
				if((OUTPUT_FALSE_POSITIVE && fp) || (OUTPUT_FALSE_NEGATIVE && fn) || 
						(OUTPUT_TRUE_POSITIVE && tp) || (OUTPUT_TRUE_NEGATIVE && tn)){
					//Print
				}else{
					if(TAG_CUSTOM) {r1.nullifyObj();}
					return bestInsert;
				}
			}
			
			if(bestInsert>0 || MIX_BAD_AND_GOOD){
				if(listg!=null){
					if(joined!=null){
						listg.add(joined);
					}else{
						listg.add(r1);
					}
				}
			}else if(listb!=null){
				listb.add(r1);
			}
			return bestInsert;
		}
		
		private final int preprocess(final Read r1, final Read r2, boolean qtrim){
			assert(r1!=null);
			if(!r1.validated()){r1.validate(true);}
			if(r2==null){return RET_BAD;}
			if(!r2.validated()){r2.validate(true);}
			
			if(iupacToN){
				if(r1!=null){r1.convertUndefinedTo((byte)'N');}
				if(r2!=null){r2.convertUndefinedTo((byte)'N');}
			}
			
			if(recalibrateQuality){
				CalcTrueQuality.recalibrate(r1);
				CalcTrueQuality.recalibrate(r2);
			}
			
			pairsProcessed++;
			basesProcessed+=r1.pairLength();
			
			if(forceTrimLeft>0 || forceTrimRight>0 || forceTrimModulo>0 || forceTrimRight2>0){
				if(r1!=null && !r1.discarded()){
					final int len=r1.length();
					final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
					final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
					final int b1=forceTrimRight>0 ? forceTrimRight : len;
					final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
					final int b=Tools.min(b0, b1, b2);
					final int x=TrimRead.trimToPosition(r1, a, b, 1);
				}
				if(r2!=null && !r2.discarded()){
					final int len=r2.length();
					final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
					final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
					final int b1=forceTrimRight>0 ? forceTrimRight : len;
					final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
					final int b=Tools.min(b0, b1, b2);
					final int x=TrimRead.trimToPosition(r2, a, b, 1);
				}
			}
			
			if(qtrim){qtrim(r1, r2, 0);}
			
			if(tadpole!=null && extendRight1>0){
				extendAndMerge(r1, r2, extendRight1, 1, false, extendRvec);
			}

			final int len1=r1.length(), len2=r2.length(); //Validated; these cannot be null.
			
			if(len1<minReadLength && len2<minReadLength){
				return RET_BAD;
			}else if(len1<2 || len2<2){
				return RET_AMBIG;
			}
			
			if(r1.quality!=null || r2.quality!=null){
				if(minAvgQuality>0){
					if(r1.avgQuality(false, minAvgQualityBases)<minAvgQuality || r2.avgQuality(false, minAvgQualityBases)<minAvgQuality){
						//Failed quality filter
						return RET_BAD;
					}
				}
				if(useMEEfilter && useQuality){
					int maxBasesToConsider=Tools.min(Tools.max(len1, len2), len1+len2-minInsert);
					if(r1.expectedTipErrors(false, maxBasesToConsider)>maxExpectedErrors || r2.expectedTipErrors(false, maxBasesToConsider)>maxExpectedErrors){
						//Failed MEEFilter
						return RET_BAD;
					}
				}
			}
			return 1;
		}
		
		private final int extendAndMerge(Read r1, Read r2, int amt, int iters, boolean merge, int[] extendRvec){
			assert(iters>0);
			assert(merge || iters==1);
			
			extendRvec[0]=extendRvec[1]=0;
			
			int sum1=0, sum2=0, attempted=0;
			int bestInsert=RET_AMBIG;
			for(int i=0; i<iters && (bestInsert==RET_AMBIG || bestInsert==RET_NO_SOLUTION); i++){
				
				int e1=(sum1==attempted ? extendRead(r1, amt) : 0);
				r2.reverseComplement();
				int e2=(sum2==attempted ? extendRead(r2, amt) : 0);
				r2.reverseComplement();
				
				attempted+=amt;
				sum1+=e1;
				sum2+=e2;
				
				if(merge){
					if(e1>0 || e2>0){
						bestInsert=processReadPair_inner(r1, r2);
					}else{
						break;
					}
				}
			}
			//Todo: un-extend.
			
			extensionsAttemptedT+=2;
			
			if(sum1==attempted){fullyExtendedT++;}
			else if(sum1>0){partlyExtendedT++;}
			else{notExtendedT++;}
			
			if(sum2==attempted){fullyExtendedT++;}
			else if(sum2>0){partlyExtendedT++;}
			else{notExtendedT++;}
			
			extendRvec[0]=sum1;
			extendRvec[1]=sum2;
			
			return bestInsert;
		}
		
		private final int extendRead(Read r, int amt){
			bbT.clear();
			bbT.append(r.bases);
			final int initialLen=r.length();
			final int extension=tadpole.extendToRight2(bbT, leftCounts, rightCounts, amt, false);
			
//			extensionsAttemptedT++;
//			if(extension==amt){
//				fullyExtendedT++;
//			}else if(extension>0){
//				partlyExtendedT++;
//			}else{
//				notExtendedT++;
//			}
			
			if(extension>0){
				r.bases=bbT.toBytes();
				if(r.quality!=null){
					r.quality=KillSwitch.copyOf(r.quality, r.bases.length);
					for(int i=initialLen; i<r.quality.length; i++){
						r.quality[i]=qfake;
					}
				}
			}
			return extension;
		}
		
		private final int mateByOverlap_ratioMode(Read r1, Read r2, int minOverlap){
			assert(useRatioMode);
			int min0=MIN_OVERLAPPING_BASES_0-MIN_OVERLAPPING_BASES_RATIO_REDUCTION;
			int min=minOverlap-MIN_OVERLAPPING_BASES_RATIO_REDUCTION;
			int x=-1;
			rvector[4]=0;

			float ratioMargin=RATIO_MARGIN;
			float maxRatio=MAX_RATIO;
			float minSecondRatio=MIN_SECOND_RATIO;

			boolean overlapped=false;
//			assert(false) : overlapUsingQuality+", "+(r1.quality!=null)+", "+(r2.quality!=null);
			if(overlapUsingQuality && r1.quality!=null && r2.quality!=null){
				overlapped=true;
				x=BBMergeOverlapper.mateByOverlapRatio(r1, r2, aprob, bprob, rvector,
						min0, min, minInsert0, minInsert, maxRatio, minSecondRatio, ratioMargin, RATIO_OFFSET, 0.95f, 0.95f, true);
				if(verbose){outstream.println("Result from ratiomode1:  \tinsert="+x+", bad="+rvector[2]+", ambig="+(rvector[4]==1));}
			}
			if(((!TAG_CUSTOM) || !overlapUsingQuality) && (!overlapped || (overlapWithoutQuality && (x<0 || rvector[4]==1)))){
				int x2=BBMergeOverlapper.mateByOverlapRatio(r1, r2, aprob, bprob, rvector, min0, min,
						minInsert0, minInsert, maxRatio, minSecondRatio, ratioMargin, RATIO_OFFSET, 0.95f, 0.95f, false);
				if(verbose){outstream.println("Result from ratiomode2:  \tinsert="+x+", bad="+rvector[2]+", ambig="+(rvector[4]==1));}
				x=x2;
//				if(x<0 || (x2>0 && rvector[4]==0)) { //Basically I want to keep the old info if it still didn't overlap...
			}
			return x;
		}
		
		private final int mateByOverlap_flatMode(Read r1, Read r2, int minOverlap){
			final int len1=r1.length(), len2=r2.length();
			boolean ambigNM=false;
			int bestInsertNM=-1;
			int bestBadNM=999999;
			
			assert(QUAL_ITERS>0);
			final int maxQualIters=(r1.quality==null || r2.quality==null ? 1 : QUAL_ITERS);
			final int maxTrims=(r1.quality==null || r2.quality==null ? 0 : TRIM_ON_OVERLAP_FAILURE);

			for(int i=0; i<maxQualIters && bestInsertNM<0 /*&& !ambigNM*/; i++){
				
				int x=BBMergeOverlapper.mateByOverlap(r1, r2, aprob, bprob, rvector, MIN_OVERLAPPING_BASES_0-i, minOverlap+i,
						minInsert0, MISMATCH_MARGIN, MAX_MISMATCHES0, MAX_MISMATCHES, (byte)(MIN_QUALITY-2*i));
				if(x>-1){
					bestInsertNM=x;
					bestBadNM=rvector[2];
					ambigNM=(rvector[4]==1);
					break;
				}
			}
			
			
			if(loose && bestInsertNM<0){//TODO check for estimated number of overlap errors
				int x=BBMergeOverlapper.mateByOverlap(r1, r2, aprob, bprob, rvector, MIN_OVERLAPPING_BASES_0-1, minOverlap+2,
						minInsert0, MISMATCH_MARGIN, MAX_MISMATCHES0+1, MAX_MISMATCHES+1, MIN_QUALITY-1);
				if(x>-1){
					bestInsertNM=x;
					bestBadNM=rvector[2];
					ambigNM=(rvector[4]==1);
				}
			}
			
			for(int trims=0, q=(int)trimq[0]; trims<maxTrims && !qtrim1 && bestInsertNM<0 /*&& !ambigNM*/; trims++, q+=8){
				float e=QualityTools.PROB_ERROR[q];
				int tr1=TrimRead.trimFast(r1, false, true, q, e, 1+len1*4/10); //r1.length());
				int tr2=TrimRead.trimFast(r2, true, false, q, e, 1+len2*4/10); //r2.length());
				if(tr1>0 || tr2>0){
					int x=BBMergeOverlapper.mateByOverlap(r1, r2, aprob, bprob, rvector, MIN_OVERLAPPING_BASES_0-1, minOverlap,
							minInsert0, MISMATCH_MARGIN, MAX_MISMATCHES0, MAX_MISMATCHES, MIN_QUALITY);
					if(x>-1){
						bestInsertNM=x;
						bestBadNM=rvector[2];
						ambigNM=(rvector[4]==1);
						trims=maxTrims;
					}
				}
			}
			if(verbose){outstream.println("Result from flatmode:  \tinsert="+bestInsertNM+", bad="+bestBadNM+", ambig="+ambigNM);}
			
			rvector[0]=bestInsertNM;
			rvector[2]=bestBadNM;
			rvector[4]=(ambigNM ? 1 : 0);
			return bestInsertNM;
		}
		
		private final int calcMinOverlapFromEntropy(final Read r1, final Read r2){
			if(!useEntropy){return MIN_OVERLAPPING_BASES;}
			final int minOverlap;
			if(loose){
				final int len1=r1.length(), len2=r2.length();
				int a=BBMergeOverlapper.calcMinOverlapByEntropy(r1.bases, entropyK, kmerCounts, minEntropyScore);
				int b=BBMergeOverlapper.calcMinOverlapByEntropy(r2.bases, entropyK, kmerCounts, minEntropyScore);
				float errorRate=r1.expectedErrors(false, len1)/len1+r2.expectedErrors(false, len2)/len2;
				minOverlap=(int)(Tools.max(MIN_OVERLAPPING_BASES, Tools.max(a, b))*(1+Tools.min(0.04f, errorRate)*4f));
			}else{
				int a=BBMergeOverlapper.calcMinOverlapByEntropyTail(r1.bases, entropyK, kmerCounts, minEntropyScore);
				int b=BBMergeOverlapper.calcMinOverlapByEntropyHead(r2.bases, entropyK, kmerCounts, minEntropyScore);
				minOverlap=Tools.max(MIN_OVERLAPPING_BASES, Tools.max(a, b));
			}
			return minOverlap;
		}
		
		private final int lookForAdapters(final Read r1, final Read r2){
			assert(lowercaseAdapters);
			if(!lowercaseAdapters){return -1;}
			if(!Tools.isLowerCase(r1.bases[r1.length()-1]) || !Tools.isLowerCase(r2.bases[0])){return -1;}
			
			final int lower1=r1.trailingLowerCase(), lower2=r2.leadingLowerCase();

			final int upper1=r1.length()-lower1, upper2=r2.length()-lower2;
			final int newlen=Tools.min(upper1, upper2);
			int good=0, bad=0;

			for(int i=0; i<newlen; i++){
				byte a=r1.bases[i];
				byte b=r2.bases[i+lower2];
				if(a!='N' && b!='N'){
					if(a==b){good++;}
					else{bad++;}
				}
			}
			if(bad*4<=good){
				rvector[0]=newlen;
				rvector[2]=bad;
				rvector[4]=0;
				return newlen;
			}
			return -1;
		}
		
		private final int mateByOverlap(Read r1, Read r2){
			final int len1=r1.length(), len2=r2.length();
			assert(r1.fetchObject()==null);
			final ByteBuilder bb=(TAG_CUSTOM ? new ByteBuilder(32) : null);
			if(bb!=null) {r1.setObj(bb);}
			if(MAKE_VECTOR && r2.obj==null) {
				if(r2.obj==null) {r2.setObj(new FloatList(24));}
				else{((FloatList)r2.obj).clear();}
			}
			final FloatList vector=(FloatList) r2.obj;
			
			final int minOverlap=calcMinOverlapFromEntropy(r1, r2);
			if(verbose){outstream.println("\nminOverlap: "+minOverlap);}
			if(TAG_CUSTOM || MAKE_VECTOR){
				int maxBasesToConsider=Tools.min(Tools.max(len1, len2), len1+len2-minInsert);
				float r1ee=r1.expectedTipErrors(false, maxBasesToConsider);
				float r2ee=r2.expectedTipErrors(false, maxBasesToConsider);
				if(TAG_CUSTOM) {
					bb.append(minOverlap*0.1f, 1).tab();
					bb.append(Tools.format("%.4f\t%.4f\t", r1ee, r2ee));
				}
				if(MAKE_VECTOR) {
					vector.add(minOverlap*0.1f); //mo
					vector.add(r1ee);
					vector.add(r2ee);
				}
			}
//			//TODO: Currently this is not used for anything.
//			final int bestInsertAD;
//			final int bestBadAD;
//			if(lowercaseAdapters){
//				bestInsertAD=lookForAdapters(r2, r2);
//				bestBadAD=(bestInsertAD>=0 ? rvector[2] : 0);
//			}
			
			if(aprob==null || aprob.length<Tools.max(len1, len2)){aprob=new float[Tools.max(len1, len2)];}
			if(bprob==null || bprob.length<Tools.max(len1, len2)){bprob=new float[Tools.max(len1, len2)];}
			
			final boolean ambigRM;
			final int bestBadRM, bestInsertRM;
			if(useRatioMode){
				bestInsertRM=mateByOverlap_ratioMode(r1, r2, minOverlap);
//				assert(!TAG_CUSTOM || bestInsertRM>0) : bestInsertRM;
				bestBadRM=rvector[2];
				ambigRM=(bestInsertRM>-1 ? rvector[4]==1 : false);
//				assert(false) : bestBadRM+", "+bestInsertRM+", "+Arrays.toString(rvector);
			}else{
				bestInsertRM=-1;
				bestBadRM=0;
				ambigRM=false;
			}
			if(verbose){System.err.println("A: "+ambigRM+", "+bestInsertRM);}
			final boolean ambigNM;
			final int bestInsertNM, bestBadNM;
			if(useFlatMode && ((!requireRatioMatch && (bestInsertRM<0 || ambigRM)) || 
					(requireRatioMatch && (bestInsertRM>0 && !ambigRM)))){
				bestInsertNM=mateByOverlap_flatMode(r1, r2, minOverlap);
				bestBadNM=rvector[2];
				ambigNM=(bestInsertNM>-1 ? rvector[4]==1 : false);
				if(verbose){System.err.println("A1: "+bestInsertNM);}
			}else{
				ambigNM=false;
				bestInsertNM=-1;
				bestBadNM=99999;
				if(verbose){System.err.println("A2: "+bestInsertNM);}
			}
			
			boolean ambig;
			int bestBad, bestInsert;
			if(requireRatioMatch && useFlatMode && useRatioMode){
				ambig=ambigRM || ambigNM;
				bestBad=bestBadRM;
				bestInsert=(bestInsertNM==bestInsertRM ? bestInsertNM : -1);

				if(verbose){outstream.println("Result after rrm:  \tinsert="+bestInsertNM+", bad="+bestBadNM+", ambig="+ambigNM);}
				if(verbose){System.err.println("A3: "+bestInsert);}
			}else if(useRatioMode && bestInsertRM>-1 && !ambigRM){
				ambig=ambigRM;
				bestBad=bestBadRM;
				bestInsert=bestInsertRM;
				if(verbose){System.err.println("A4: "+bestInsert);}
			}else{
				ambig=ambigNM;
				bestBad=bestBadNM;
				bestInsert=bestInsertNM;
				if(verbose){System.err.println("A5: "+bestInsert);}
			}
			if(verbose){System.err.println("B: "+ambig+", "+bestInsert+", "+bestInsertRM);}
			
			if((TAG_CUSTOM || MAKE_VECTOR) && bestInsert<0){bestInsert=bestInsertRM;}
			
			if(bestBad>MAX_MISMATCHES_R){ambig=true;}
			if(verbose){System.err.println("B2: "+ambig+", "+bestInsert+", "+bestBad+", "+MAX_MISMATCHES_R);}
			
			if((!TAG_CUSTOM && !MAKE_VECTOR) || ONLY_OUTPUT_FALSE_POSITIVE){
				if(ambig){return RET_AMBIG;}
				else if(bestInsert<0){return RET_NO_SOLUTION;}
			}
//			assert(false) : vector.size;
			//TODO:  Crucial!  This block can vastly reduce merge rate, particularly if quality values are inaccurate.
			if(useQuality && r1.quality!=null && r2.quality!=null){
				if(useEfilter && bestInsert>0 && (!ambig || TAG_CUSTOM || MAKE_VECTOR)){
					float bestExpected=BBMergeOverlapper.expectedMismatches(r1, r2, bestInsert);
					if((bestExpected+efilterOffset)*efilterRatio<bestBad){ambig=true;}
					
					if(verbose){outstream.println("Result after efilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
					if(TAG_CUSTOM){
//						bb.append("_be=f");
						bb.append(bestExpected, 4).tab();
					}
					if(MAKE_VECTOR){
						vector.add(bestExpected); //be
					}
				}
				if(verbose){System.err.println("C: "+ambig+", "+bestInsert);}

				if(pfilterRatio>0 && bestInsert>0 && (!ambig || TAG_CUSTOM || MAKE_VECTOR)){
					float probability=BBMergeOverlapper.probability(r1, r2, bestInsert);
					if(probability<pfilterRatio && ((!TAG_CUSTOM && !MAKE_VECTOR) || ONLY_OUTPUT_FALSE_POSITIVE)){bestInsert=-1;}
					if(verbose){outstream.println("Result after pfilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
					double p=Math.sqrt(probability)+0.0000015;
					if(TAG_CUSTOM){
						bb.append(p, 6).tab();
					}
					if(MAKE_VECTOR){
						vector.add((float)p);//pr
					}
				}else if(MAKE_VECTOR) {
					vector.add(0.1f);//For some reason this gets left off occasionally...
				}
				if(verbose){System.err.println("D: "+ambig+", "+bestInsert);}
			}
			if(verbose){System.err.println("E: "+ambig+", "+bestInsert);}
			
			//TODO: Uncomment.
//			if(TAG_CUSTOM && !ONLY_OUTPUT_FALSE_POSITIVE) {
//				assert(!TAG_CUSTOM || bestInsert>0) : bestInsert;
//				return bestInsert;
//			}
//			
//			if(ambig){return RET_AMBIG;}
			r1.errors=bestBad;
			
			if(MAKE_VECTOR && bestInsert>0 && r2.obj()!=null) {
				FloatList list=r2.fetchVector();
				if(list.size!=23) {
					System.err.println("\n"+r1.numericID+", "+r1.obj());
					System.err.println("*"+list.size+", "+parseCustom+", "+r1.numericID+", "+list);
					assert(false);
				}else {
//					System.err.println(list.size+", "+parseCustom+", "+r1.numericID+", "+list);
					//				assert((list.size()==23 && !parseCustom) || (list.size()==24 && parseCustom)) : 
					//					list.size+", "+parseCustom+", "+r1.numericID+", "+list;
					assert(list.size()==23) : 
						list.size+", "+parseCustom+", "+r1.numericID+", "+list;
					net.applyInput(list);
					float result=net.feedForward();
					//				System.err.println(list);
					//				System.err.println(r1.obj);
					
//					System.err.println("\n"+r1.numericID+", "+r1.obj());
//					System.err.println(r1.numericID+", "+list);
//					System.err.println(result);
					
//					if(parseCustom){
//						final int trueSize=parseCustom ? parseInsert(r1, r2) : -1;
//						final int initialLen1=r1.length(), initialLen2=r2.length();
//						final int minReadLen=Tools.min(initialLen1, initialLen2), lenSum=initialLen1+initialLen2;
//						final int trueOverlap=Tools.min(trueSize, minReadLen, lenSum-trueSize);
//
//						System.err.println(r1.numericID+", "+result+", "+netCutoff+", "+
//								(trueOverlap>0)+", "+trueOverlap);
//					}

					if(result<netCutoff){
						return RET_AMBIG;
					}
				}
			}
			
			return bestInsert>0 ? bestInsert : RET_NO_SOLUTION;
		}
		
		private final int parseInsert(Read r1, Read r2){
			int trueSize=-1;
			if(r1.id.startsWith("SYN")){
				CustomHeader h=new CustomHeader(r1.id, 0);
				trueSize=h.insert;
			}else if(r1.id.startsWith("insert=")){
				trueSize=GradeMergedReads.parseInsert(r1.id);
			}else{
				r1.setMapped(true);
				r2.setMapped(true);
				trueSize=Read.insertSizeMapped(r1, r2, false);
				assert(trueSize>0) : "Can't parse insert size.";
			}
			if(verbose){outstream.println("True Insert: "+trueSize);}
//			r1.setInsert(trueSize);
			return trueSize;
		}
		
		private boolean verifyAdapters(Read r1, Read r2, int bestInsert){
			final int minAdapterOverlap=3;
			final int len1=r1.length(), len2=r2.length();
			if(Tools.max(len1, len2)-minAdapterOverlap<bestInsert){return true;}
			
			int good=0, bad=0, invalid=0;
			final float minAdapterRatio=0.6f;
			{
				final int adapterLen=len1-bestInsert;
				if(adapterList1!=null && good<1 && adapterLen>=minAdapterOverlap){
					if(adapterLen>=minAdapterOverlap){
						adaptersExpectedT++;
						for(byte[] adapter : adapterList1){
							int aiv=adapterIsValid(r1, adapter, bestInsert, minAdapterOverlap, minAdapterRatio);
							if(aiv==2){good++; break;}
							else if(aiv==1){invalid++;}
							else{bad++;}
						}
					}
				}
			}
			{
				final int adapterLen=len2-bestInsert;
				if(adapterList2!=null && good<1 && adapterLen>=minAdapterOverlap){
					if(adapterLen>=minAdapterOverlap){
						adaptersExpectedT++;
						r2.reverseComplement();
						for(byte[] adapter : adapterList2){
							int aiv=adapterIsValid(r2, adapter, bestInsert, minAdapterOverlap, minAdapterRatio);
							if(aiv==2){good++; break;}
							else if(aiv==1){invalid++;}
							else{bad++;}
						}
						r2.reverseComplement();
					}
				}
			}
			if(good>0){
				adaptersFoundT+=good;
				return true;
			}else if(bad>0){
				return false;
			}else if(invalid>0 && strict){
				return false;
			}
			return true;
		}
		
		/**
		 * 
		 * @param r1 Read1
		 * @param r2 Read2
		 * @return A return code (RET_)
		 */
		private final int processReadPair(final Read r1, final Read r2){
			final int initialLen1=r1.length(), initialLen2=r2.length();
			{
				final int x=preprocess(r1, r2, (qtrim1 && !qtrim2));
				if(!TAG_CUSTOM && !MAKE_VECTOR && x<0){return x;}
			}
			
			r2.reverseComplement();
			
			byte[] qual1=r1.quality, qual2=r2.quality;
			if(!useQuality){//strip qualities
				r1.quality=r2.quality=null;
			}
			
			final int bestInsert0=processReadPair_inner(r1, r2);
			int bestInsert=bestInsert0;
			assert(!TAG_CUSTOM || ONLY_OUTPUT_FALSE_POSITIVE || MAKE_VECTOR || bestInsert>0) : bestInsert;
			if(qtrim2){
				for(int iter=0; iter<trimq.length && bestInsert<0; iter++){
					r1.quality=qual1;
					r2.quality=qual2;
					//				r2.reverseComplement();
					//				qtrim(r1, r2);
					//				r2.reverseComplement();

					TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq[iter], trimE[iter], 1);
					TrimRead.trimFast(r2, qtrimRight, qtrimLeft, trimq[iter], trimE[iter], 1);//Reversed because read is rcomped

					qual1=r1.quality;
					qual2=r2.quality;
					bestInsert=processReadPair_inner(r1, r2);
				}
			}
			
			//assert(false) : bestInsert+", "+r1.length()+", "+r2.length();
			final boolean foundAdapter;
//			final ByteBuilder bb=(TAG_CUSTOM ? r1.fetchBB() : null);
			if(verifyAdapters && bestInsert>0 && (extendRight1<1 || bestInsert<initialLen1 || bestInsert<initialLen2)){
				foundAdapter=verifyAdapters(r1, r2, bestInsert);
//				if(TAG_CUSTOM) {bb.append(foundAdapter ? "fa_1" : "fa_0");}
				if(!foundAdapter){
					if(!useQuality){//restore qualities
						restoreQualities(r1, r2, qual1, qual2);
					}
					assert(r1.quality==null || (r1.quality.length==r1.bases.length && r2.quality.length==r2.bases.length));
					r2.reverseComplement();
					return RET_NO_SOLUTION;
				}
			}else{
				foundAdapter=false;
//				if(TAG_CUSTOM) {bb.append("fa_0.5");}
			}
			
			//Kmer processing - correction and extension
			if(tadpole!=null){
				
				//Correction
				boolean corrected=false;
				if(eccTadpole /*&& !requireExtensionMatch*/ && (bestInsert==RET_AMBIG || bestInsert==RET_NO_SOLUTION)){
					int c1=tadpole.errorCorrect(r1);
					int c2=tadpole.errorCorrect(r2);
					corrected=true;
					if(c1>0 || c2>0){
						bestInsert=processReadPair_inner(r1, r2);
					}
				}
				
				//Extension
				if(extendRight2>0 && (requireExtensionMatch || bestInsert==RET_AMBIG || bestInsert==RET_NO_SOLUTION)
						&& !(requireStrictExtensionMatch && bestInsert<1)){
					final int lengthSum=r1.pairLength();
					final int approxMaxOverlappingInsert=lengthSum-minApproxOverlapRem;
					final int minExt=Tools.min(12, extendRight2*2);
					int bestInsertE=extendAndMerge(r1, r2, extendRight2, extendIterations, true, extendRvec);
					int extension1=extendRvec[0], extension2=extendRvec[1];
					int extension=extension1+extension2;
					
//					assert(r1.length()==150+extension1) : extension1+", "+r1.length();
					
					//If extension failed, try correction and extend again
					if(eccTadpole && requireExtensionMatch && !corrected && extension<=extendRight2/* *2 */ &&
							(extension==0 || (!strict && bestInsert<0 && bestInsertE<0) || (bestInsert>0 && bestInsertE==bestInsert))){
						int c1=0, c2=0;
						if(extension1==0){c1=tadpole.errorCorrect(r1);}
						if(extension2==0){c2=tadpole.errorCorrect(r2);}
						
						//If there were any corrections
						if(c1>0 || c2>0){
							//Un-extend r1
							if(extension1>0){TrimRead.trimToPosition(r1, 0, r1.length()-extension1-1, 1);}
							//Un-extend r2
							if(extension2>0){TrimRead.trimToPosition(r2, extension2, r2.length()-1, 1);}
							
							//Retry extension
							bestInsertE=extendAndMerge(r1, r2, extendRight2, extendIterations, true, extendRvec);
							extension1=extendRvec[0];
							extension2=extendRvec[1];
							extension=extension1+extension2;
						}
					}
					
//					if(extension==0 && eccTadpole && !corrected){
//						int c1=tadpole.errorCorrect(r1);
//						int c2=tadpole.errorCorrect(r2);
//						corrected=true;
//						if(c1>0 || c2>0){
//							bestInsertE=extendAndMerge(r1, r2, extendRight2, extendIterations, true, extendRvec);
//							extension=extendRvec[0]+extendRvec[1];
//						}
//					}
					
					if(requireExtensionMatch && extension>0){
						if(bestInsert!=bestInsertE){
							if(bestInsert>0){bestInsert=RET_AMBIG;}//There was an unextended overlap, not matching extended overlap
							else if(!requireStrictExtensionMatch && bestInsertE>0){//Extended had an overlap, but unextended did not
								if(bestInsertE>approxMaxOverlappingInsert && extension>=minExt){//Use only if the extended overlap was too long to see before
									bestInsert=bestInsertE;
								}
							}
						}
					}else if(bestInsertE>0){
						assert(bestInsert<0);
						bestInsert=bestInsertE;
					}
				}
				
				if(testMerge && bestInsert>0 && bestInsert>=r1.length()+testMergeWidth && bestInsert>=r2.length()+testMergeWidth){
					Read merged=r1.joinRead(bestInsert);
					if(!tadpole.mergeOK(merged, initialLen1, initialLen2, mergeOKBitsetT, countList, kmerT, testMergeWidth, testMergeThresh, testMergeMult)){
						bestInsert=RET_NO_SOLUTION;
					}
				}
			}else if(bloomFilter!=null){
				if(eccBloom && (bestInsert==RET_AMBIG || bestInsert==RET_NO_SOLUTION)){
					int c1=corrector.errorCorrect(r1);
					int c2=corrector.errorCorrect(r2);
					if(c1>0 || c2>0){
						bestInsert=processReadPair_inner(r1, r2);
					}
				}
				
				if(testMerge && bestInsert>0 && bestInsert>=r1.length()+testMergeWidth && bestInsert>=r2.length()+testMergeWidth){
					Read merged=r1.joinRead(bestInsert);
					if(!corrector.mergeOK(merged, initialLen1, initialLen2, kmers, testMergeWidth, testMergeThresh, testMergeMult)){bestInsert=RET_NO_SOLUTION;}
				}
			}
			
			if(verifyAdapters && !foundAdapter && bestInsert>0 && (extendRight1<1 || bestInsert<initialLen1 || bestInsert<initialLen2)){
				boolean x=verifyAdapters(r1, r2, bestInsert);
				if(!x){bestInsert=RET_NO_SOLUTION;}
			}
			
//			assert(bestInsert<0 || bestInsert>=50) : bestInsert+", "+foundAdapter;
			
			if(useKFilter && bestInsert>kmerLength){
				Read joined=r1.joinRead(bestInsert);
				if(useKFilter){
					int cov1=BBMergeOverlapper.minCoverage(r1, tadpole, kmerLength, filterCutoff+1);
					int cov2=(cov1>filterCutoff ? BBMergeOverlapper.minCoverage(r2, tadpole, kmerLength, filterCutoff+1) : cov1);
					if(cov1>filterCutoff && cov2>filterCutoff){
						int cov=BBMergeOverlapper.minCoverage(joined, tadpole, kmerLength, filterCutoff);
						if(cov<filterCutoff){bestInsert=RET_NO_SOLUTION;}
						if(verbose){outstream.println("Result after kfilter:  \tinsert="+bestInsert);}
					}
				}
			}
			
			if(!useQuality){//restore qualities
				restoreQualities(r1, r2, qual1, qual2);
			}
			r2.reverseComplement();
			assert(r1.quality==null || (r1.quality.length==r1.bases.length && r2.quality.length==r2.bases.length));
			return bestInsert;
		}
		
		private void restoreQualities(Read r1, Read r2, byte[] qual1, byte[] qual2){
			final int len1=r1.length(), len2=r2.length();
			r1.quality=qual1;
			r2.quality=qual2;
			final byte qf=Shared.FAKE_QUAL;
			if(qual1.length!=len1){
				r1.quality=Arrays.copyOf(qual1, len1);
				for(int i=len1; i<r1.quality.length; i++){
					r1.quality[i]=qf;
				}
			}
			if(qual2.length!=len2){
				r2.reverseComplement();
				r2.quality=Arrays.copyOf(qual2, len2);
				for(int i=len2; i<r2.quality.length; i++){
					r2.quality[i]=qf;
				}
				r2.reverseComplement();
			}
		}
		
		/**
		 * 
		 * @param r1 Read1
		 * @param r2 Read2
		 * @return A return code (RET_)
		 */
		private final int processReadPair_inner(final Read r1, final Read r2){
			int bestInsert=-1;
			boolean ambig;
			
			if(USE_MAPPING && r1.chrom==r2.chrom && r1.start<r1.stop && ((r1.mapped() || r1.synthetic()) 
					&& (r2.mapped() || r2.synthetic()))){
				bestInsert=r1.insert();
				ambig=false;
			}else{
				if(MATE_BY_OVERLAP){
					bestInsert=mateByOverlap(r1, r2);
					assert(!TAG_CUSTOM || ONLY_OUTPUT_FALSE_POSITIVE || MAKE_VECTOR || bestInsert>0) : bestInsert;
					ambig=(bestInsert==RET_AMBIG);
				}else{
					ambig=false;
					bestInsert=-1;
				}
			}
			
			if((TAG_CUSTOM) && !ONLY_OUTPUT_FALSE_POSITIVE) {return bestInsert;}
			if(ambig && !MAKE_VECTOR){return RET_AMBIG;}
			else if(bestInsert>0){
				if(bestInsert<minInsert){return RET_SHORT;}
				else if(bestInsert>maxReadLength){return RET_LONG;}
				return bestInsert;
			}
			else{return RET_NO_SOLUTION;}
		}
		
		private void storeAdapterSequence(Read r, int insert){
			atrackT.storeAdapterSequence(r, insert);
//			if(ignorePhixAdapters){
//				if(looksLikePhix(r, insert)){return;}
//			}
//			
//			LongList[] lists=adapterCountsT[r.pairnum()];
//			byte[] bases=r.bases;
//			
//			for(int i=insert, j=0; i<bases.length; i++, j++){
//				byte b=bases[i];
//				int num=AminoAcid.baseToNumber[b];
//				if(num>=0){
//					lists[num].increment(j);
//				}
//			}
		}
		
		private boolean looksLikePhix(Read r, int insert){
			return looksLikePhix(r.bases, insert) || looksLikePhix(r.mate.bases, insert);
		}
		
		private boolean looksLikePhix(byte[] bases, int insert){
			int len=bases.length-insert;
			if(len<phixPrefix.length){return false;}
			for(int i=insert, j=0; i<bases.length && j<phixPrefix.length; i++, j++){
				byte b=bases[i];
				if(b!='N' && b!=phixPrefix[j]){
					return false;
				}
			}
			outstream.println(new String(bases).substring(insert));
			outstream.println(new String(phixPrefix));
			return true;
		}
		
		/*--------------------------------------------------------------*/
		
//		final LongList[][] adapterCountsT=new LongList[2][4];
		final AdapterTracker atrackT=new AdapterTracker();

		private final BitSet mergeOKBitsetT=new BitSet(400);
		private final Kmer kmerT;
		private IntList countList=new IntList();
		private LongList kmers=new LongList();
		
		final byte qfake=Shared.FAKE_QUAL;
		
		private final int[] rvector=KillSwitch.allocInt1D(5);
		private final int[] extendRvec=KillSwitch.allocInt1D(2);

		private final int[] rightCounts=KillSwitch.allocInt1D(4);
		private final int[] leftCounts=(extendThroughLeftJunctions ? null : new int[4]);
		
		private final ByteBuilder bbT=new ByteBuilder();
		
		final long[] hist=new long[histlen];
		final short[] kmerCounts;
		
		private float[] aprob, bprob;

		long errorsCorrectedT=0;
		
		long pairsProcessed=0;
		long basesProcessed=0;
		long matedCount=0;
		long correctCount=0;
		long unmergedOverlappingCount=0;
		long ambiguousCount=0;
		long tooShortCount=0;
		long tooLongCount=0;
		long incorrectCount=0;
		long noSolutionCount=0;
		long insertSumCorrect=0;
		long insertSumIncorrect=0;
		int insertMax=0;
		int insertMin=999999999;

		long fullyExtendedT=0;
		long partlyExtendedT=0;
		long notExtendedT=0;
		long extensionsAttemptedT=0;
		
		long adaptersExpectedT=0;
		long adaptersFoundT=0;
		
		private final ConcurrentReadInputStream cris;
		private final ConcurrentReadOutputStream rosgood;
		private final ConcurrentReadOutputStream rosbad;
		private final ConcurrentReadOutputStream rosi;
		
		private final boolean joinReads;
		private final boolean trimReadsByOverlap;
		
		private final CellNet net;
	}
	
	/*--------------------------------------------------------------*/
	
	public static void resetCounters(){

		Arrays.fill(histTotal, 0);

		readsProcessedTotal=0;
		basesProcessedTotal=0;
		matedCountTotal=0;
		correctCountTotal=0;
		ambiguousCountTotal=0;
		tooShortCountTotal=0;
		tooLongCountTotal=0;
		incorrectCountTotal=0;
		noSolutionCountTotal=0;
		insertSumCorrectTotal=0;
		insertSumIncorrectTotal=0;
		fullyExtendedTotal=0;
		partlyExtendedTotal=0;
		notExtendedTotal=0;
		extensionsAttempted=0;
		adaptersExpected=0;
		adaptersFound=0;

		insertMinTotal=999999999;
		insertMaxTotal=0;
	}
	
	private static synchronized void loadNet(String fname) {
		if(fname==null) {
			netFile=null;
			net0=null;
			return;
		}
		if(netFile!=null && netFile.equals(fname) && net0!=null){
			return;
		}
		netFile=fname;
		net0=CellNetParser.load(netFile);
		netCutoff=(setCutoff ? netCutoff : net0.cutoff);
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1;
	private String in2;
	
	private ArrayList<String> extra=new ArrayList<String>();

	private String out1=null;
	private String out2=null;
	private String outb1=null;
	private String outb2=null;
	private String outinsert=null;
	private String ihist=null;
	private String outAdapter=null;
	private String outCardinality=null;

	/** List of R1 adapters to look for to verify a short-insert overlap */
	private ArrayList<byte[]> adapterList1;
	/** List of R2 adapters to look for to verify a short-insert overlap */
	private ArrayList<byte[]> adapterList2;
	/** Require short inserts to have adapter sequence at the expected location */
	private boolean verifyAdapters=false;
	
	final CardinalityTracker loglog;
	
	private long maxReads=-1;
	private boolean join=true;
	private boolean ecco=false;
	private boolean trimByOverlap=false;
	
	private float pfilterRatio=0.00004f;
	private float efilterRatio=6f;
	private float efilterOffset=0.05f;
	private boolean useEfilter=true;
	private boolean useMEEfilter=false;
	
	private boolean ordered=false;
	private boolean overlapUsingQuality=false;
	private boolean overlapWithoutQuality=true;
	private boolean useKFilter=false;
	private int filterCutoff=1;
	private int kmerLength=31;
	private boolean prealloc=false;
	private int prefilter=0;
	private long filterMemoryOverride=0;
	private boolean eccTail=false;
	private boolean eccPincer=false;
	private boolean eccReassemble=true;

	private boolean useEntropy=true;
	private int entropyK=3;
	private int minEntropyScore=39;//30 loose;//39 normal;//44 strict;
	
	private long sampleseed=-1;
	private float samplerate=1;
	
	private long errorsCorrectedTotal=0;
	
	private boolean findAdapterSequence=false;
	private boolean ignorePhixAdapters=false;
	
	private final AdapterTracker atrack=new AdapterTracker();
	
	private final Tadpole tadpole;
	private int extendRight1=0;
	private int extendRight2=0;
	private int extendIterations=1;
	private boolean eccTadpole=false;
	private boolean shave=false;
	private boolean rinse=false;

	private final BloomFilter bloomFilter;
	private final BloomFilterCorrector corrector;
	
	private boolean forceExactKmerCounts=false;
	private boolean testMerge=false;
	private boolean eccBloom=false;
	private int bloomBits=2;
	private int bloomHashes=3;
	int testMergeWidth=4;
	long testMergeMult=80L;
	int testMergeThresh=3;
	
	private boolean requireExtensionMatch=false;
	private boolean requireStrictExtensionMatch=false;
	private int minApproxOverlapRem=26;
	
	private boolean extendThroughLeftJunctions=true;
	private int minCountSeed=3, minCountExtend=2;
	private float branchMult1=20;
	private float branchMult2=3;
	private float minProb=0.5f;
	private int branchLowerConst=3;
	
	private boolean MATE_BY_OVERLAP=true;
	private boolean MIX_BAD_AND_GOOD=false;
	private boolean ONLY_OUTPUT_INCORRECT=false;
	
	private boolean OUTPUT_FALSE_POSITIVE=true;
	private boolean OUTPUT_FALSE_NEGATIVE=true;
	private boolean OUTPUT_TRUE_POSITIVE=true;
	private boolean OUTPUT_TRUE_NEGATIVE=true;
	
	private boolean quantizeQuality=false;
	
	/*--------------------------------------------------------------*/
	
	private static ThreadLocal<int[]> localRvector=new ThreadLocal<int[]>();
	
	static boolean errorState=false;

	private static boolean showFullArgs=true;

	/** Recalibrate quality scores using matrices */
	static boolean recalibrateQuality=false;
	static boolean useQuality=true;
	static boolean qtrimRight=false;
	static boolean qtrimLeft=false;
//	static boolean untrim=false;
	static float[] trimq=new float[] {6};
	static float[] trimE=QualityTools.phredToProbError(trimq);
	static float minAvgQuality=0;
	static int minAvgQualityBases=0;
	static float maxExpectedErrors=0;
	static int minReadLength=1;
	static int maxReadLength=-1;
	static int minInsert=15;
	static int minInsert0=-1;
	static boolean qtrim1=false;
	static boolean qtrim2=false;
	static int TRIM_ON_OVERLAP_FAILURE=1;
	static int QUAL_ITERS=3;
	static boolean parseCustom=false;
	static int maxAdapterLength=21;
	
	static boolean trimNonOverlapping=false;
	
	/** For adapter output */
	static boolean trimPolyA=true;
	
	static int forceTrimLeft;
	static int forceTrimRight;
	static int forceTrimRight2;
	/** Trim right bases of the read modulo this value.
	 * e.g. forceTrimModulo=50 would trim the last 3bp from a 153bp read. */
	static int forceTrimModulo;
	
	public static boolean strict=false;
	static boolean vstrict=false;
	static boolean ustrict=false;
	static boolean xstrict=false;
	static boolean loose=false;
	static boolean vloose=false;
	static boolean uloose=false;
	static boolean xloose=false;
	static boolean fast=false;
	
	/** If true, interpret lowercase bases as adapter sequence */
	static boolean lowercaseAdapters=false;
	
	private static final int histlen=2000;
	static long[] histTotal=new long[histlen];
	static int bin=1;

	static long readsProcessedTotal=0;
	static long basesProcessedTotal=0;
	static long matedCountTotal=0;
	static long correctCountTotal=0;
	static long unmergedOverlappingTotal=0;
	static long ambiguousCountTotal=0;
	static long tooShortCountTotal=0;
	static long tooLongCountTotal=0;
	static long incorrectCountTotal=0;
	static long noSolutionCountTotal=0;
	static long insertSumCorrectTotal=0;
	static long insertSumIncorrectTotal=0;
	static long fullyExtendedTotal=0;
	static long partlyExtendedTotal=0;
	static long notExtendedTotal=0;
	static long extensionsAttempted=0;
	static long adaptersExpected=0;
	static long adaptersFound=0;
	
	static int insertMinTotal=999999999;
	static int insertMaxTotal=0;
	
	private static int MIN_OVERLAPPING_BASES=11;
	private static int MIN_OVERLAPPING_BASES_0=8;
	private static int MISMATCH_MARGIN=2;
	private static int MIN_OVERLAPPING_BASES_RATIO_REDUCTION=3;
	
	/** Skip alignment and calculate insert from mapping info */
	protected static boolean USE_MAPPING=false;
	private static boolean NONZERO_ONLY=true;
	private static boolean showHistStats=true;
	static boolean TRAINING=false;
	
	static boolean useRatioMode=true;
	static boolean useFlatMode=false;
	static boolean requireRatioMatch=false;
	static int MAX_MISMATCHES_R=20;
	static float MAX_RATIO=0.09f;
	static float RATIO_MARGIN=5.5f;
	static float RATIO_OFFSET=0.55f;
	static float MIN_SECOND_RATIO=0.1f;
	
	public static int MAX_MISMATCHES=3;
	public static int MAX_MISMATCHES0=3;
	public static byte MIN_QUALITY=10;
	
	public static final int RET_NO_SOLUTION=-1;
	public static final int RET_AMBIG=-2;
	public static final int RET_BAD=-3;
	public static final int RET_SHORT=-4;
	public static final int RET_LONG=-5;
	
	private static boolean overwrite=true;
	private static boolean append=false;
	private static final boolean verbose=false;
	static final boolean TAG_CUSTOM=false;
	static boolean MAKE_VECTOR=false;
	private static String defaultNetFile=Data.findPath("?bbmerge.bbnet", false);
	static String netFile=defaultNetFile;
	static CellNet net0=CellNetParser.load(netFile, false);//False in case the net file is missing
	static float netCutoff=(net0==null ? 0.5f : net0.cutoff);
	private static boolean setCutoff=false;
	private static boolean ONLY_OUTPUT_FALSE_POSITIVE=false;
	private static boolean ONLY_OUTPUT_FALSE_NEGATIVE=false;

	//TODO: Remove after NN is trained
	private static boolean EXEMPT_SHORT_OVERLAPS_FROM_STATS=false;
	
	private static boolean iupacToN=false;
	
	public static boolean changeQuality=true;
	
	static PrintStream outstream=System.err;
	
	private static int THREADS=-1;
//	private static String version="9.02";
	
	private static final byte[] phixPrefix="AGATCGGAAGAGCG".getBytes();
	
	private static final String[] defaultAdapters={
			"AGATCGGAAGAGCACA", "AATGATACGGCGACCA",	"GATCGGAAGAGCACAC", "CTGTCTCTTATACACA",
			"GACGCTGCCGACGA", "CCGAGCCCACGAGAC", "CTGATGGCGCGAGGGA", "CTGAGCGGGCTGGCAA",
			"GATCGGAAGAGCGTCG", "GATCGTCGGACTGTAG", "CCTTGGCACCCGAGAA", "CCACGGGAACGTGGTG",
			"TGGAATTCTCGGGTGC", "TCGGACTGTAGAACTC", "AGATCGGAAGAGCGGT"
	};
	
	public static ArrayList<byte[]> staticAdapterList=getAdapterList("default");
	
}
