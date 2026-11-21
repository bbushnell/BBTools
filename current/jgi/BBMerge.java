package jgi;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

import align2.QualityTools;
import assemble.Tadpole;
import assemble.TargetKmerBridge;
import bloom.BloomFilter;
import bloom.BloomFilterCorrector;
import bloom.BloomFilterCorrector1;
import bloom.KmerCountAbstract;
import cardinality.CardinalityTracker;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.KmerTableSet;
import ml.CellNet;
import ml.CellNetParser;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import stream.CustomHeader;
import stream.FASTQ;
import stream.Read;
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
 * Merges overlapping paired-end reads and provides overlap/error-correction helpers.
 * Command-line processing uses shared static settings and delegates multi-kmer
 * requests to {@link BBMergeMulti}.
 * @author Brian Bushnell
 * @date Aug 14, 2012
 *
 */
public class BBMerge{
	
	/**
	 * Delegates multi-kmer requests to BBMergeMulti; otherwise constructs and runs BBMerge.
	 * The local path enables both feature-vector flags before construction. On normal
	 * return from processing, it clears those flags, enables Read constructor validation,
	 * and closes the redirected diagnostic stream.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		if(BBMergeMulti.isMultiK(args)){
			BBMergeMulti.main(args);
			return;
		}
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
	
	/**
	 * Applies shared preset, diagnostic and JNI settings, consuming preset and showfullargs
	 * tokens by nulling their input slots. JNI tokens remain for subsequent parsing.
	 * When tokens are consumed, prepends applicable preset defaults to the remaining
	 * explicit arguments in their original order. Presets also depend on existing shared mode flags.
	 * @param args Input argument array to mutate, or null for no arguments
	 * @return Empty array for null input, the original array if no tokens were consumed,
	 * or a new array containing preset defaults followed by the remaining arguments
	 */
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
	
	/**
	 * Parses arguments and configures this run and shared process settings.
	 * Initializes requested calibration and optional Tadpole or Bloom support before
	 * final input/output path checks; Bloom setup can count input reads.
	 * Loads a requested nondefault neural model after those checks.
	 * The configured merge/correction pass runs in {@link #process()}.
	 * @param args Command-line arguments handled by PreParser and preset expansion
	 */
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
			//JGI-092: Keep complete values; the native splitter preserves bare/empty options.
			String[] split=Parse.splitOnFirst(arg, '=');
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
			}else if(a.equals("noheader") || a.equals("suppressheader")){
				stream.ReadStreamWriter.NO_HEADER=Parse.parseBoolean(b);
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
				setKmerLength=true;
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
				//JGI-087: Each mate-specific option updates only its own list.
				adapterList1=getAdapterList(b);
			}else if(a.equals("adapter2") || a.equals("adapters2")){
				adapterList2=getAdapterList(b);
			}
			
			//Extension parameters
			
			else if(a.equals("extendright") || a.equals("er") || a.equals("extend") || a.equals("extendright1") || a.equals("er1") || a.equals("extend1")){
				extendRight1=Parse.parseIntKMG(b);
			}else if(a.equals("extendright2") || a.equals("er2") || a.equals("extend2")){
				extendRight2=Parse.parseIntKMG(b);
			}else if(a.equals("extenditerations") || a.equals("iterations") || a.equals("ei") || a.equals("iters")){
				extendIterations=Tools.max(1, (int)Parse.parseKMG(b));
			}else if(a.equals("hollow")){
				hollow=Parse.parseBoolean(b);
			}else if(a.equals("hollowbranchstop")){
				hollowBranchStop=Parse.parseBoolean(b);
			}else if(a.equals("hollowbranchratio")){
				hollowBranchRatio=Parse.parseIntKMG(b);
			}else if(a.equals("hollowbranchsupport")){
				hollowBranchSupport=Parse.parseIntKMG(b);
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
			}else if(a.equals("kmerbridge") || a.equals("targetbridge") || a.equals("bridge")){
				kmerBridge=Parse.parseBoolean(b);
				setKmerBridge=true;
			}else if(a.equals("bridgeonly") || a.equals("targetbridgeonly")){
				kmerBridge=Parse.parseBoolean(b);
				setKmerBridge=true;
				if(kmerBridge){MATE_BY_OVERLAP=false;}
			}else if(a.equals("bridgedistance") || a.equals("bridgedist") || a.equals("bridgemax")){
				bridgeDistance=Parse.parseIntKMG(b);
			}else if(a.equals("bridgetargetwindow") || a.equals("bridgewindow")){
				bridgeTargetWindow=Parse.parseIntKMG(b);
			}else if(a.equals("bridgeminhits")){
				bridgeMinHits=Parse.parseIntKMG(b);
			}else if(a.equals("bridgemaxoverlap")){
				bridgeMaxOverlap=Parse.parseIntKMG(b);
			}else if(a.equals("bridgereciprocal")){
				bridgeReciprocal=Parse.parseBoolean(b);
			}else if(a.equals("hashkmers") || a.equals("kmerhash") || a.equals("bridgehash") || a.equals("hashbridges")){
				hashKmerMode=parseHashMode(b);
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
		
		if(outinsert!=null){
			outinsert=outinsert.replaceFirst("#", "1");
			if(TRAINING){throw new IllegalArgumentException("BBMerge training export through outinsert is not implemented; normal outinsert writes five-column insert reports.");}
		}
		if(adapterList1!=null || adapterList2!=null){
			verifyAdapters=true;
		}
		
		if(requireStrictExtensionMatch){
			requireExtensionMatch=true;
		}
		if(setKmerLength && !setKmerBridge && kmerLength>0 && extendRight1<1 && extendRight2<1 &&
				!requireExtensionMatch && !useKFilter && !eccTadpole && !eccBloom && !forceExactKmerCounts && !testMerge){
			kmerBridge=true;
		}

		if(requireExtensionMatch && extendRight2<1){
			outstream.println("Extend2 is defaulting to 50 because it was unset but r"+(requireStrictExtensionMatch ? "s" : "")+"em mode is being used.");
			extendRight2=50;
		}
		if(kmerBridge && bridgeDistance<1){
			bridgeDistance=(extendRight2>0 ? 2*extendRight2+kmerLength-bridgeMaxOverlap : 500);
		}
		if(kmerBridge && (bridgeDistance<1 || bridgeTargetWindow!=0 || bridgeMinHits!=1)){
			throw new RuntimeException("Exact-tip kmer bridging requires distance>0, targetWindow=0, and minHits=1: distance="+
					bridgeDistance+", targetWindow="+bridgeTargetWindow+", minHits="+bridgeMinHits);
		}
		//Only the ordinary read-extension consumer has defined Hollow semantics.
		//Other consumers may reconstruct solid keys or interpret counts differently.
		if(hollowBranchRatio<1 || hollowBranchSupport<1){
			throw new IllegalArgumentException("hollowbranchratio and hollowbranchsupport must be positive.");
		}
		if(hollowBranchStop && (!hollow || kmerLength<67)){
			throw new IllegalArgumentException("hollowbranchstop=t requires hollow=t and k>=67.");
		}
		if(hollow && (kmerLength<=64 || extendRight1>0 || extendRight2<1
				|| kmerBridge || eccTadpole || eccBloom || ecco || useKFilter || testMerge
				|| shave || rinse || prefilter>0 || hashKmerMode==0 || hashKmerMode==2)){
			throw new IllegalArgumentException("BBMerge hollow=t requires single-k extend2/REM at k>64 with extend=0, "
					+"hashkmers=auto or pair, prefilter=0, shave=f, rinse=f, bridge=f, "
					+"ecct=f, eccb=f, ecco=f, kfilter=0 and testmerge=f.");
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
		
		loglog=(outCardinality==null && !parser.loglog ? null : CardinalityTracker.makeTracker(2048, 31, 0, 0));
		
		{//Process parser fields
			Parser.processQuality();
			
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=(parser.trimq2!=null ? parser.trimq2 : new float[]{parser.trimq});
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
		
		//JGI-093: Check the resolved report path before graph setup or any output worker starts.
		if(outinsert!=null && (!Tools.testOutputFiles(overwrite, append, false, outinsert) ||
				!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outb1, outb2, outinsert, ihist, outCardinality, outAdapter))){
			throw new IllegalArgumentException("Cannot write BBMerge insert report without overwriting an input or another output: "+outinsert);
		}
		if(extendRight1>0 || extendRight2>0 || useKFilter || eccTadpole || forceExactKmerCounts || kmerBridge){
			
			final long mem=Shared.memAvailable();
			if(Shared.memAvailable()<2000000000L){
				System.err.println("\n***** WARNING *****\nUsing kmer counts uses a lot of memory, but only "+(mem/(1024*1024))+"MB is available.");
				System.err.println("If this process crashes, run bbmerge-auto.sh instead of bbmerge.sh, or set the -Xmx flag.\n");
			}
			
			ArrayList<String> list=new ArrayList<String>();
			{
				StringBuilder sb=new StringBuilder("in=");
				sb.append(in1);
				list.add(sb.toString());
				if(extra.size()>0){
					StringBuilder eb=new StringBuilder("extra=");
					for(int i=0; i<extra.size(); i++){
						if(i>0){eb.append(',');}
						eb.append(extra.get(i));
					}
					list.add(eb.toString());
				}
			}
			if(in2!=null){list.add("in2="+in2);}
			list.add("branchlower="+branchLowerConst);
			list.add("branchmult1="+branchMult1);
			list.add("branchmult2="+branchMult2);
			list.add("mincountseed="+minCountSeed);
			list.add("mincountextend="+minCountExtend);
			list.add("minprob="+minProb);
			list.add("k="+kmerLength);
			if(hollow){
				//BBMerge only loads the table, but Tadpole validates its consumer mode.
				list.add("mode=extend");
				list.add("hollow=t");
				list.add("hollowbranchstop="+hollowBranchStop);
				list.add("hollowbranchratio="+hollowBranchRatio);
				list.add("hollowbranchsupport="+hollowBranchSupport);
			}
			final int hashMode=effectiveHashMode();
			list.add("hashkmers="+(hashMode==2 ? "fixed" : hashMode==1 ? "pair" : "f"));
			list.add("seedfromtable=f");
			//BBMerge loads counts but never calls Tadpole.process2 to build contigs.
			//Only optional Shaver cleaning allocates ownership arrays. Reserving
			//four bytes per entry otherwise needlessly lowers the resize ceiling.
			list.add("ownership="+(shave || rinse));
			list.add("prealloc="+(prealloc || hashMode==2));
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
			corrector=new BloomFilterCorrector1(bloomFilter, kmerLength, kmerLength);
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
		if(netFile==null){net0=null;}
		else if(!netFile.equals(defaultNetFile)){net0=CellNetParser.load(netFile);}
	}
	
	/**
	 * Loads configured Tadpole kmers, optionally cleans them, and runs the processing phase.
	 * After the phase returns, requests histogram output, writes requested adapter/cardinality
	 * results and prints shared statistics. Requested histogram writing is asynchronous.
	 */
	void process(){
		Timer ttotal=new Timer();
		ttotal.start();
		
		if(tadpole!=null){
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
		if(kmerBridge){
			outstream.println("Kmer Bridge Attempts: \t"+bridgeAttemptsTotal);
			outstream.println("Kmer Bridges:         \t"+bridgeSuccessTotal);
			outstream.println("Bridge No Path:       \t"+bridgeNoPathTotal);
			outstream.println("Bridge No Target:     \t"+bridgeNoTargetTotal);
			outstream.println("Bridge Ambiguous:     \t"+bridgeAmbiguousTotal);
			outstream.println("Bridge Cyclic:        \t"+bridgeCycleTotal);
			outstream.println("Bridge Nonreciprocal: \t"+bridgeNonreciprocalTotal);
			outstream.println("Bridge Rejected:      \t"+bridgeRejectedTotal);
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
	
	/**
	 * Formats the shared insert histogram, rounding each bin's mean frequency up.
	 * Uses shared bin width, maximum insert, statistics and zero-row settings. A bin reaching
	 * the array's end uses its actual entry count; optional statistics use the raw histogram.
	 * A null destination does nothing. Otherwise the histogram must be initialized, bin
	 * width positive, counts nonnegative and each grouped total must fit in a long.
	 * Starts asynchronous file writing and returns without waiting for completion.
	 * @param fname Output filename, or null to skip output
	 * @param percentMerged Caller-supplied percentage printed only with summary statistics
	 */
	public static void writeHistogram(String fname, double percentMerged){
		if(fname==null){return;}
		assert(bin>0) : "Histogram bin width must be positive to compute row averages: "+bin;
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
			//JGI-088: Preserve long frequencies when accumulating an output bin.
			long x=0;
			int y=0;
			for(int j=i; j<i+bin && j<histTotal.length; j++){
				final long count=histTotal[j];
				assert(count>=0 && count<=Long.MAX_VALUE-x) :
						"Histogram averages require nonnegative counts and a bin sum that fits in a long: index="+j+", count="+count+", sum="+x;
				x+=count;
				y++;
			}
			//JGI-088/091: Round by the actual bin size without overflowing an added numerator.
			x=x/y+(x%y==0 ? 0 : 1);
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

	
	/**
	 * Starts configured output writers and an input reader, then runs and joins workers.
	 * Uses shared settings and resets and aggregates run counters. Input limits are
	 * forwarded to the selected readers; they are not retained-sample targets.
	 * Interruption during a worker join does not cancel processing. The owner's
	 * interrupt status is restored after the stream-cleanup attempt, including
	 * when an existing worker failure is propagated.
	 * Earlier setup or result-aggregation failures can bypass stream cleanup.
	 *
	 * @param join Whether workers should join accepted overlapping mates
	 * @param maxReads Limit forwarded unchanged to each input reader; units depend on the reader
	 * @param perfectonly Currently unused
	 */
	public void runPhase(boolean join, long maxReads, boolean perfectonly){
		
		Timer talign=new Timer();
		
		Writer rosgood=null;
		Writer rosbad=null;
		ByteStreamWriter rosinsert=null;
		
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
			
			//FIXED [jgi/BBMerge two-file-ordered]: two-file paired output requires ordered=true so the new
			//WriterFactory/PairedWriter keeps R1/R2 mates synced across the two files (and to satisfy its assert).
			//Single-file (interleaved/merged) output keeps the user's 'ordered' preference.
			final boolean orderedG=(out2!=null) || ordered;
			final FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, orderedG);
			final FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, orderedG);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";

			final int buff=Tools.max(16, 2*THREADS);
			rosgood=WriterFactory.getStream(ff1, ff2, null, null, buff, null, false, 1);
			rosgood.start();
		}
		
		if(outb1!=null){

			//FIXED [jgi/BBMerge two-file-ordered]: same as the out1/out2 block above — outu1/outu2 (unmerged pairs in
			//two files) need ordered=true for the new PairedWriter. This is the path in the user's AssertionError report.
			final boolean orderedB=(outb2!=null) || ordered;
			final FileFormat ff1=FileFormat.testOutput(outb1, FileFormat.FASTQ, null, true, overwrite, append, orderedB);
			final FileFormat ff2=FileFormat.testOutput(outb2, FileFormat.FASTQ, null, true, overwrite, append, orderedB);
			assert(!ff1.samOrBam()) : "Sam files need reference info for the header.";

			final int buff=Tools.max(16, 2*THREADS);
			rosbad=WriterFactory.getStream(ff1, ff2, null, null, buff, null, false, 1);
			rosbad.start();
		}
		
		if(outinsert!=null){
			final FileFormat ff=FileFormat.testOutput(outinsert, FileFormat.TEXT, ".txt", true, overwrite && !append, append, ordered);
			final boolean printHeader=!(stream.ReadStreamWriter.NO_HEADER || (append && ff.exists()));
			rosinsert=new ByteStreamWriter(ff);
			//Header precedes startup; ordered data IDs remain the reader's dense IDs from zero.
			if(printHeader){rosinsert.forcePrint(header());}
			rosinsert.start();
		}
		
		
		if(rosgood!=null || rosbad!=null || rosinsert!=null){
			outstream.println("Started output threads.");
		}
		
		final Streamer cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=StreamerFactory.getReadInputStream(maxReads, true, ff1, ff2, 1);
			cris.setSampleRate(samplerate, sampleseed);
			if(verbose){outstream.println("Started cris");}
			cris.start();
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
		
		boolean interrupted=false;
		Throwable threadFailure=null;
		for(int i=0; i<pta.length; i++){
			MateThread ct=pta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try{
						ct.join(1000);
					}catch(InterruptedException e){
						//JGI-089: Keep joining with a clear flag; restore it only after stream cleanup.
						interrupted=true;
						e.printStackTrace();
					}
				}
				if(threadFailure==null && ct.failure!=null){threadFailure=ct.failure;}

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
				bridgeAttemptsTotal+=ct.bridgeAttemptsT;
				bridgeSuccessTotal+=ct.bridgeSuccessT;
				bridgeNoPathTotal+=ct.bridgeNoPathT;
				bridgeNoTargetTotal+=ct.bridgeNoTargetT;
				bridgeAmbiguousTotal+=ct.bridgeAmbiguousT;
				bridgeCycleTotal+=ct.bridgeCycleT;
				bridgeNonreciprocalTotal+=ct.bridgeNonreciprocalT;
				bridgeRejectedTotal+=ct.bridgeRejectedT;

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

				//Combine this thread's cardinality tracker into the shared loglog (thread-safe: merged at join).
				if(loglog!=null && ct.loglogT!=null){loglog.add(ct.loglogT);}
			}
		}
		
//		outstream.println("Finished reading");
		try{
			boolean reportError=false;
			try{
				errorState|=ReadWrite.closeStreams(cris, rosgood, rosbad);
			}finally{
				if(rosinsert!=null){
					try{reportError=rosinsert.poisonAndWait();}
					finally{errorState|=rosinsert.errorState;}
				}
			}

			talign.stop();
			if(reportError){throw new RuntimeException("Error writing BBMerge insert report: "+outinsert);}
			//FIXED [jgi/BBMerge worker failures]: uncaught MateThread exceptions previously
			//produced partial zero/NaN statistics and a successful process exit.
			if(threadFailure!=null){throw new RuntimeException("BBMerge worker thread failed.", threadFailure);}
		}finally{
			if(interrupted){Thread.currentThread().interrupt();}
		}
//		outstream.println("Align time: "+talign);
	}
	
	/**
	 * Estimates the positive-insert fraction using makeInsertHistogram's loose overlap search
	 * and current shared settings. Reads are not joined or error-corrected.
	 * @param fname1 Nonnull primary input path
	 * @param fname2 Optional mate input path
	 * @param numReads Limit forwarded unchanged to each input reader; units depend on the reader
	 * @param samplerate Sampling rate passed to the reader
	 * @return Positive-insert counts divided by all histogram counts, or zero when the
	 * histogram is absent, has fewer than two bins, or has a total below one
	 */
	public static final float mergeableFraction(String fname1, String fname2, long numReads, float samplerate){
		long[] hist=makeInsertHistogram(fname1, fname2, numReads, samplerate);
		if(hist==null || hist.length<2){return 0;}
		long sum=simd.Vector.sum(hist);
		return sum<1 ? 0 : (sum-hist[0])/(float)sum;
	}
	
	/**
	 * Counts loose overlap candidates from paired input using current shared settings,
	 * with error correction disabled and sampling seed 1. Positive insert lengths index
	 * their own bins; all nonpositive results enter bin zero. Counts are not normalized.
	 * Attempts stream closure on normal completion or an unpaired-input return; earlier
	 * read/search failures can bypass closure, and returned close-error flags are ignored.
	 * @param fname1 Nonnull primary input path
	 * @param fname2 Optional mate input path
	 * @param numReads Limit forwarded unchanged to each input reader; units depend on the reader
	 * @param samplerate Sampling rate passed to the reader
	 * @return Raw insert counts, or null for primary standard input or a reader reporting unpaired input
	 */
	public static final long[] makeInsertHistogram(String fname1, String fname2, long numReads, float samplerate){
		assert(fname1!=null);
		final Streamer cris;
		{
			FileFormat ff1=FileFormat.testInput(fname1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(fname2, FileFormat.FASTQ, null, true, true);
			if(ff1.stdio()){return null;}
			assert(!ff1.stdio()) : "Standard in is not allowed as input when calculating insert size distributions for files.";
			cris=StreamerFactory.getReadInputStream(numReads, true, ff1, ff2, 1);
			cris.setSampleRate(samplerate, 1);
			if(verbose){outstream.println("Started cris");}
			cris.start();
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
		while(ln!=null && reads!=null){//ln!=null prevents a compiler potential null access warning

			for(Read r1 : reads){
				int x=findOverlapLoose(r1, r1.mate, false);
				if(x>0){ll.increment(x, 1);}
				else{ll.increment(0, 1);}
			}
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		ReadWrite.closeStreams(cris);
		return ll.toArray();
	}

	/** Searches original-orientation mates with the ultra-strict overlap preset.
	 * Uses the validation, mutation and shared-setting rules of findOverlap.
	 * @param r1 Nonnull first read
	 * @param r2 Nonnull second read; must be r1.mate when ecc is true
	 * @param ecc Whether to correct an accepted candidate using the joined consensus
	 * @return Candidate insert length in bases, or a negative result on rejection */
	public static final int findOverlapUStrict(final Read r1, final Read r2, final boolean ecc){
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

	/** Searches original-orientation mates with the very-strict overlap preset.
	 * Uses the validation, mutation and shared-setting rules of findOverlap.
	 * @param r1 Nonnull first read
	 * @param r2 Nonnull second read; must be r1.mate when ecc is true
	 * @param ecc Whether to correct an accepted candidate using the joined consensus
	 * @return Candidate insert length in bases, or a negative result on rejection */
	public static final int findOverlapVStrict(final Read r1, final Read r2, final boolean ecc){
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
	
	/** Finds a strict overlap and joins linked mates through Read.joinRead.
	 * Inputs enter in original orientation. The search may validate reads and store
	 * the insert; it does not request ECC. R2 is temporarily reverse-complemented
	 * for joining and restored on normal return. Joining has Read.joinRead's metadata
	 * side effects; this is not a read-only overlap query or an exception-safe restore.
	 * @param r1 Nonnull first read, whose mate must be r2 for a successful join
	 * @param r2 Original-orientation mate; null returns null without searching
	 * @return Joined read for a positive candidate, or null when absent or rejected */
	public static final Read tryToMerge(final Read r1, final Read r2){
		if(r2==null){return null;}
		final int insert=findOverlapStrict(r1, r2, false);
		if(insert<1){return null;}
		r2.reverseComplementFast();
		final Read r3=r1.joinRead(insert);
		r2.reverseComplementFast();
		return r3;
	}

	/** Searches original-orientation mates with the strict overlap preset.
	 * Uses the validation, mutation and shared-setting rules of findOverlap.
	 * @param r1 Nonnull first read
	 * @param r2 Nonnull second read; must be r1.mate when ecc is true
	 * @param ecc Whether to correct an accepted candidate using the joined consensus
	 * @return Candidate insert length in bases, or a negative result on rejection */
	public static final int findOverlapStrict(final Read r1, final Read r2, final boolean ecc){
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

	/** Searches original-orientation mates with the loose overlap preset.
	 * Uses the validation, mutation and shared-setting rules of findOverlap.
	 * @param r1 Nonnull first read
	 * @param r2 Nonnull second read; must be r1.mate when ecc is true
	 * @param ecc Whether to correct an accepted candidate using the joined consensus
	 * @return Candidate insert length in bases, or a negative result on rejection */
	public static final int findOverlapLoose(final Read r1, final Read r2, final boolean ecc){
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
	
	/** Finds an insert from original-orientation mates using mapping or ratio overlap.
	 * Validates reads when needed, temporarily reverse-complements R2 and restores its
	 * orientation on normal return. Exceptions do not have a finally-based restore.
	 * Quality-dependent expected-mismatch/probability filters and the shared adapter
	 * list can reject a candidate. Ambiguity returns -1. An unambiguous nonnegative
	 * result is passed to r1.setInsert. For positive candidates, ECC additionally
	 * requires minInsert and may replace read arrays/error counters through the helper.
	 * A positive mapped candidate below minInsert can still be returned without ECC.
	 * Scratch is thread-local, but configuration and adapter lists are shared mutable
	 * state; callers must not treat this as an independently configured merger.
	 * @param r1 Nonnull first read
	 * @param r2 Nonnull second read; must be r1.mate when ecc is true
	 * @param ecc Whether to request correction for an eligible candidate
	 * @param minOverlap Incoming value, replaced by entropy/global thresholds on the overlap path
	 * @param minOverlap0 Overlap threshold in bases passed unchanged to the ratio search
	 * @param minInsert Search insert threshold and minimum length required for both reads
	 * @param minInsert0 Insert threshold in bases passed unchanged to the ratio search
	 * @param entropy Entropy target; nonpositive uses the global minimum overlap
	 * @param maxRatio Maximum mismatch ratio passed to the ratio search
	 * @param minSecondRatio Second-candidate ratio threshold passed to the ratio search
	 * @param ratioMargin Candidate-separation multiplier passed to the ratio search
	 * @param ratioOffset Mismatch-count offset passed to the ratio search
	 * @param efilterRatio Multiplier applied to the sum of expected mismatches and efilterOffset
	 * @param efilterOffset Offset added to expected mismatches before filtering
	 * @param pfilterRatio Positive probability threshold, or nonpositive to skip that filter
	 * @return Candidate insert length in bases, or a negative result on rejection */
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
		
		r2.reverseComplementFast();
		
		int bestInsert=-1;
		int bestBad=999999;
		boolean ambig, tooShort=false;
		
		if(USE_MAPPING && r1.chrom==r2.chrom && r1.start<r1.stop && r1.mapped() && r2.mapped()){
			bestBad=0;
			bestInsert=Read.insertSizeMapped(r1, r2, false);
			ambig=false;
		}else{
			if(entropy>0){
				final int a=BBMergeOverlapper.calcMinOverlapByEntropy(r1.bases, 3, null, entropy);
				final int b=BBMergeOverlapper.calcMinOverlapByEntropy(r2.bases, 3, null, entropy);
				minOverlap=Tools.max(MIN_OVERLAPPING_BASES, Tools.max(a, b));
			}else{minOverlap=MIN_OVERLAPPING_BASES;}
			if(verbose){outstream.println("minOverlap: "+minOverlap);}
			
			rvector[4]=0;

			final int x=BBMergeOverlapper.mateByOverlapRatio(r1, r2, null, null, rvector, minOverlap0, minOverlap,
					minInsert0, minInsert, maxRatio, minSecondRatio, ratioMargin, ratioOffset, 0.95f, 0.95f, false);
			bestInsert=x;
			bestBad=rvector[2];
			ambig=(x>-1 ? rvector[4]==1 : false);
		}
		
		//TODO:  Crucial!  This line can vastly reduce merge rate, particularly if quality values are inaccurate.
		if(bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null){
			final float bestExpected=BBMergeOverlapper.expectedMismatches(r1, r2, bestInsert);
			if((bestExpected+efilterOffset)*efilterRatio<bestBad){ambig=true;}
			if(verbose){outstream.println("Result after efilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
		}
		
		//TODO:  Crucial!  This line can vastly reduce merge rate, particularly if quality values are inaccurate.
		if(pfilterRatio>0 && bestInsert>0 && !ambig && r1.quality!=null && r2.quality!=null){
			final float probability=BBMergeOverlapper.probability(r1, r2, bestInsert);
			if(probability<pfilterRatio){bestInsert=-1;}
			if(verbose){outstream.println("Result after pfilter:  \tinsert="+bestInsert+", bad="+bestBad+", ambig="+ambig);}
		}
		
		if(staticAdapterList!=null && bestInsert>0 && (bestInsert<len1 || bestInsert<len2)){
			final boolean foundAdapter=verifyAdaptersStatic(r1, r2, bestInsert);
			if(!foundAdapter){
				bestInsert=RET_NO_SOLUTION;
				assert(r1.quality==null || (r1.quality.length==r1.bases.length && r2.quality.length==r2.bases.length));
				r2.reverseComplementFast();
				return RET_NO_SOLUTION;
			}
		}
		
		tooShort=(!ambig && bestInsert>0 && bestInsert<minInsert);
		
		if(ecc && bestInsert>-1 && !ambig && !tooShort){
			final int errors=errorCorrectWithInsert(r1, r2, bestInsert);
		}
		
		if(r2!=null){r2.reverseComplementFast();}
		if(!ambig && bestInsert>-1){r1.setInsert(bestInsert);}
		
		return ambig ? -1 : bestInsert;
	}
	
	/** Replaces linked, already-oriented mates with slices of their joined consensus.
	 * R2 must already face the same fragment direction as R1; this method does not
	 * reverse-complement either read. A nonempty join replaces base arrays, resets
	 * both error counters and counts changed calls whose new base is fully defined.
	 * Both mates must supply length-matched qualities or both omit them.
	 * Present quality arrays are copied from the consensus when changeQuality is true,
	 * otherwise from the original slices aligned to the retained bases, including R2's
	 * oriented suffix. Reads can be shortened to the insert.
	 * Joining also has Read.joinRead's metadata effects. Nonoverlap or an empty join
	 * returns zero without resetting the existing error counters.
	 * @param r1 Nonnull first read with r1.mate equal to r2
	 * @param r2 Nonnull mate already reverse-complemented into fragment orientation
	 * @param insert Positive fragment length in bases; no correction at or beyond the length sum
	 * @return Sum of the newly assigned error counters, or zero when no join is corrected */
	public static int errorCorrectWithInsert(final Read r1, final Read r2, final int insert){
		if(insert>=r1.length()+r2.length()){return 0;}
		assert(insert>0);
		final Read joined=r1.joinRead(insert);
		
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
					//Keep qualities aligned with R2's retained oriented suffix (JGI-095).
					//The discarded prefix contains bases outside the inferred fragment.
					r2.quality=KillSwitch.copyOfRange(r2.quality, old2.length-r2.bases.length, old2.length);
				}
			}
		}else{
			return 0;
		}
		return r1.errors+r2.errors;
	}
	
	/** Counts each read's differences from fully defined calls in an oriented consensus.
	 * R1 aligns at the consensus start and R2 at its end, both in the consensus
	 * direction. A shared coordinate may contribute once for each differing read;
	 * this is not the raw pairwise mismatch count. Inputs and error fields are unchanged.
	 * @param r1 First read with nonnull bases, already in consensus orientation
	 * @param r2 Second read with nonnull bases, already in consensus orientation
	 * @param joined Consensus; null or a nonpositive/nonoverlapping length returns zero
	 * @return Sum of read-to-consensus differences at fully defined consensus bases */
	public static int countErrors(final Read r1, final Read r2, final Read joined){
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

	/** Counts raw pairwise mismatches for original-orientation reads without modifying them.
	 * Delegates the coordinate and complement comparison to the array overload.
	 * @param r1 Nonnull first read; a null bases array yields zero
	 * @param r2 Nonnull mate in original orientation; a null bases array yields zero
	 * @param insert Fragment length in bases; nonpositive or nonoverlapping yields zero
	 * @return Number of unequal raw/complemented calls in the shared coordinates */
	public static int countErrors(final Read r1, final Read r2, final int insert){
		return countErrors(r1.bases, r2.bases, insert);
	}
	
	/** Counts raw pairwise mismatches at shared fragment coordinates without mutation.
	 * Both arrays enter in original orientation. R1 advances while R2 is traversed
	 * backward and compared through the extended-complement table. Qualities are not
	 * consulted and ambiguous bases have no separate skip rule.
	 * @param bases1 First read's encoded bases; null yields zero
	 * @param bases2 Mate's original-orientation encoded bases; null yields zero
	 * @param insert Fragment length in bases; nonpositive or at least the length sum yields zero
	 * @return Number of unequal raw/complemented calls across the covered coordinates */
	public static int countErrors(final byte[] bases1, final byte[] bases2, final int insert){
		if(insert<1 || bases1==null || bases2==null){return 0;}
		if(insert>=bases1.length+bases2.length){return 0;}

		//JGI-090: B covers fragment coordinates [insert-B.length, insert), so its length determines A's start.
		//At each shared coordinate, A's index plus B's raw index equals insert-1.
		int offset1=(insert>=bases2.length ? insert-bases2.length : 0);
		int offset2=(insert>=bases2.length ? bases2.length-1 : insert-1);
		
		int errors=0;
//		int overlap=0;
		for(int i=offset1, j=offset2; i<bases1.length && j>=0; i++, j--){
			final byte a=bases1[i], b=bases2[j];
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

	/**
	 * Returns header string for output files.
	 * Format depends on TAG_CUSTOM flag.
	 * @return Header string for data output
	 */
	public static String header(){
		if(TAG_CUSTOM){
			String i="##mo\tr1ee\tr2ee\talen\tblen\tbi\tbo\tbb\tbg\tbr\tbbi"+
					"\tsbi\tsbo\tsbb\tsbg\tsbr\tsbbi\trr\tbadr\tgoodr\tovr\tbe\tpr";
//			String o="truth\tgoal\tmerged";
//			String o="goal";
			String o="truth";
			final int inputs=Tools.tabPattern.split(i).length;
			final int outputs=Tools.tabPattern.split(o).length;
			return "#dims\t"+inputs+"\t"+outputs+"\n"+i+"\t"+o+"\n";
		}else{
			return "#id\tnumericID\tinsert\tstatus\tmismatches\n";
		}
	}

	/** Appends one normal report row before output correction/restoration mutates the pair.
	 * Counts raw/complemented current calls, ignoring qualities; rejected results have
	 * no insert/count. Does not mutate reads, scratch errors or feature attachments.
	 * @param report Worker-owned batch bytes, not yet transferred to the writer
	 * @param r1 First read with its mate in original sequencer orientation
	 * @param insert Signed final search result, before Read.setInsert clamps rejection codes
	 */
	private static void appendInsertReport(final ByteBuilder report, final Read r1, final int insert){
		assert(report!=null && r1!=null && r1.mate!=null) : "Insert report rows require a complete worker-owned pair and batch builder.";
		final int mismatches=(insert>0 ? countErrors(r1, r1.mate, insert) : -1);
		final char status=(insert>0 ? (mismatches==0 ? 'P' : 'I') : insert==RET_AMBIG ? 'A' : insert==RET_SHORT ? 'S' : 'F');
		if(r1.id==null){report.append(r1.numericID);}
		else{
			//ByteBuilder's String append truncates chars; retain ASCII speed and encode other names explicitly.
			for(int i=0; i<r1.id.length(); i++){
				final char c=r1.id.charAt(i);
				if(c<128){report.append(c);}
				else{report.append(r1.id.substring(i).getBytes(StandardCharsets.UTF_8));break;}
			}
		}
		report.tab().append(r1.numericID).tab().append(insert>0 ? insert : -1).tab().append(status).tab().append(mismatches).nl();
	}
	
	/** Trims both mates in their current orientation using shared side flags and
	 * the selected trimq/trimE entries. The active fast path may replace base/quality arrays
	 * without saving reversible trim state, requesting a minimum retained length of one.
	 * Calls are sequential; failure on the second mate does not undo the first trim.
	 * @param r1 Nonnull first read to trim
	 * @param r2 Nonnull second read to trim
	 * @param iter Index present in both shared trimming-parameter arrays
	 */
	private static void qtrim(Read r1, Read r2, int iter){
		if(false /*untrim*/){
			TrimRead.trim(r1, qtrimLeft, qtrimRight, trimq[iter], trimE[iter], 1);
			TrimRead.trim(r2, qtrimLeft, qtrimRight, trimq[iter], trimE[iter], 1);
		}else{
			TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq[iter], trimE[iter], 1);
			TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq[iter], trimE[iter], 1);
		}
	}
	
	/** Resolves an adapter specification into a caller-owned list and base arrays.
	 * The case-insensitive alias "default" always selects built-ins; "adapters" does
	 * so only when no path with that name exists. Other specifications go through
	 * Tools.toAdapterList, which reads files or comma-separated sequences/file paths,
	 * validates uppercase nucleotides, truncates to maxAdapterLength, converts degenerate
	 * bases to N, removes trailing Ns and deduplicates retained sequences. Built-ins
	 * bypass that validation/truncation path.
	 * @param name Built-in alias, file path, sequence list, or null
	 * @return New list and arrays, or null for a null specification/no retained sequences
	 */
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
	
	/** Applies the static adapter rejection policy; any good comparison wins.
	 * Tests tails using three informative bases and a 60% match threshold. With no
	 * good comparison, a bad one rejects; unknown results reject only in strict mode.
	 * A missing/empty list or no eligible tail passes without positive adapter evidence.
	 * R2 is temporarily reversed for its tail check and restored on normal completion;
	 * exceptional completion does not guarantee restoration. Callers own both reads.
	 * @param r1 First read in forward orientation
	 * @param r2 Mate already reverse-complemented into the first read's orientation
	 * @param bestInsert Nonnegative insert length; adapter tails begin at this index
	 * @return Whether the rejection policy accepts, including when no comparison is made
	 */
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
					r2.reverseComplementFast();
					for(byte[] adapter : staticAdapterList){
						int aiv=adapterIsValid(r2, adapter, bestInsert, minAdapterOverlap, minAdapterRatio);
						if(aiv==2){good++; break;}
						else if(aiv==1){invalid++;}
						else{bad++;}
					}
					r2.reverseComplementFast();
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
	
	/** Compares an adapter prefix with a read tail using exact, case-sensitive bytes.
	 * Only uppercase N in either sequence is excluded from informative-base counts.
	 * Stops early when mismatches exceed the length-based allowance; otherwise requires
	 * minAdapterOverlap informative bases and ceil(minAdapterRatio*informative) matches.
	 * Borrows both arrays without changing them; does not reverse-complement the read.
	 * @param r Read with nonnull bases in the adapter's sequencer orientation
	 * @param adapter Nonnull adapter bases
	 * @param insert Nonnegative read index at which the adapter is expected to start
	 * @param minAdapterOverlap Minimum informative overlap for a good result
	 * @param minAdapterRatio Required fraction of informative bases matching the adapter
	 * @return 0 for mismatch rejection, 1 for insufficient evidence, or 2 for a good match
	 */
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
	
	/** Per-phase worker consuming paired reads and submitting numbered output batches.
	 * Borrows shared input and optional writers; runPhase owns stream startup/cleanup,
	 * worker joins and aggregation of local counters, histograms and trackers.
	 * Owns mutable search scratch and optional processing state, including a network
	 * copy with separate row parameters/activation arrays but some shared metadata.
	 * Outer/static configuration and graph/Bloom/correction helpers remain shared;
	 * backend concurrency and helper safety are not established by this worker wrapper.
	 */
	private class MateThread extends Thread{
		
		
		/** Retains stream references and allocates this worker's processing state.
		 * Does not start this thread or any stream. Optional graph/cardinality/entropy
		 * state follows enclosing settings; the bridge borrows the enclosing Tadpole.
		 * net0.copy(false) copies network parameters without accumulated deltas, while
		 * command/statistics metadata can remain shared. Construction failures occur
		 * before run's failure capture and leave stream cleanup to the caller.
		 * @param cris_ Borrowed paired input, started by the owner before run
		 * @param rosgood_ Optional borrowed writer for accepted reads
		 * @param rosbad_ Optional borrowed writer for rejected reads
		 * @param rosi_ Optional borrowed byte writer receiving immutable normal report batches
		 * @param joinReads_ Whether to join accepted mates
		 * @param trimByOverlap_ Whether to trim accepted mates beyond their inferred insert
		 */
		public MateThread(Streamer cris_, Writer rosgood_, Writer rosbad_, ByteStreamWriter rosi_,
				boolean joinReads_, boolean trimByOverlap_){
			cris=cris_;
			rosgood=rosgood_;
			rosbad=rosbad_;
			rosi=rosi_;
			joinReads=joinReads_;
			trimReadsByOverlap=trimByOverlap_;
			kmerT=(tadpole==null || tadpole.k()<32 ? null : new Kmer(tadpole.k()));
			bridgeT=(kmerBridge ? new TargetKmerBridge(tadpole, bridgeDistance, bridgeReciprocal) : null);
			net=(net0==null ? null : net0.copy(false));
			//Unique per-thread tracker, same params (seed) as the shared loglog so they merge correctly.
			loglogT=(loglog==null ? null : CardinalityTracker.makeTracker(2048, 31, 0, 0));

			if(useEntropy){
				kmerCounts=new short[1<<(2*entropyK)];
			}else{
				kmerCounts=null;
			}
		}
		
		
		/** Runs the batch loop and records any Throwable for the owner to inspect after join.
		 * Hollow mode also prints the failure and invokes the process-fatal KillSwitch.
		 * Does not close streams, cancel peers or supply missing output batches on failure.
		 * Recording failure does not ensure cleanup or propagation: blocked input/output
		 * or missing ordered batches can prevent runPhase from reaching its final check.
		 */
		@Override
		public void run(){
			try{
				processReads();
			}catch(Throwable t){
				failure=t;
				//TODO: Probable bug - legacy ordered output can wait for this worker's
				//missing batch before runPhase reaches its final failure check.
				if(hollow){
					t.printStackTrace(outstream);
					KillSwitch.kill("Experimental Hollow BBMerge worker failed; incomplete outputs must not be used.");
				}
			}
		}
		
		/** Consumes nonnull paired Read roots; a null batch or null list ends the loop.
		 * Empty lists from nextList continue normal batch submission with their input IDs.
		 * Allocates the optional four-array originals holder once per call and reuses it;
		 * enabled good/bad lists are new per batch. Cardinality hashing precedes pair mutations.
		 * Submits good reads, immutable report bytes and bad reads in that order using the
		 * input ID, including empty report batches. Rows describe each selected pair before
		 * output mutation; good/bad qualities may be quantized before their submissions.
		 * Lists can share Read payloads; add may block or retain references, not deep-copy.
		 * Does not start/close streams, call returnList or poison writers; runPhase owns
		 * their lifecycle. Exceptions leave partial work for run's failure path.
		 */
		private void processReads(){
			assert(USE_MAPPING || MATE_BY_OVERLAP || kmerBridge);
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(r.mate!=null);
			}
			
			final byte[][] originals=((tadpole!=null || qtrimRight || qtrimLeft) &&
					(rosbad!=null || (rosgood!=null && (!join || MIX_BAD_AND_GOOD)))) ? new byte[4][] : null;
			while(ln!=null && reads!=null){//ln!=null prevents a compiler potential null access warning
				
				ArrayList<Read> listg=(rosgood==null /*&& rosi==null*/ ? null : new ArrayList<Read>(reads.size()));
				ArrayList<Read> listb=(rosbad==null ? null : new ArrayList<Read>(reads.size()));
				final ByteBuilder insertReport=(rosi==null ? null : new ByteBuilder());
				
				if(loglogT!=null){
					for(Read r1 : reads){loglogT.hash(r1);}
				}
				
				for(Read r1 : reads){
					int bestInsert=findOverlapInThread(r1, originals, listg, listb, insertReport);
				}
				
				if(rosgood!=null){
					if(quantizeQuality){Quantizer.quantize(listg);}
					rosgood.add(listg, ln.id);
				}
				if(rosi!=null){rosi.add(insertReport, ln.id);}//Transfers ownership, including empty batches.
				if(rosbad!=null){
					if(quantizeQuality){Quantizer.quantize(listb);}
					rosbad.add(listb, ln.id);
				}
				//			outstream.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
				//			outstream.println("reads: "+(reads==null ? "null" : reads.size()));
			}
		}
		
		/** Processes one linked pair, updates local statistics and selects output payloads.
		 * Pair search may trim/correct/extend reads; accepted pairs may then be joined or
		 * output-corrected. Optional originals preserve only base/quality arrays, not all
		 * Read state. Custom truth/output selection can omit a pair after statistics change.
		 * Rejected Read.insert values are clamped by setInsert, not preserved RET codes.
		 * @param r1 Nonnull input root with a nonnull mate
		 * @param originals Optional reusable four-slot base/quality snapshot holder
		 * @param listg Optional destination for accepted/joined or mixed output
		 * @param listb Optional destination for rejected roots
		 * @param insertReport Optional owned batch builder receiving the pre-output row
		 * @return Pair-search insert or RET code, not proof that a record was emitted
		 */
		private int findOverlapInThread(final Read r1, final byte[][] originals, ArrayList<Read> listg, ArrayList<Read> listb, final ByteBuilder insertReport){

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
			if(insertReport!=null){appendInsertReport(insertReport, r1, bestInsert);}
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
				if(trueOverlap>0 && trueOverlap<=8 && EXEMPT_SHORT_OVERLAPS_FROM_STATS){
					//Temporarily ignore for purpose of training
				}else if(bestInsert==trueSize){
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
					r2.reverseComplementFast();
					joined=r1.joinRead(bestInsert);
					r2.reverseComplementFast();
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
					r2.reverseComplementFast();
					errorCorrectWithInsert(r1, r2, bestInsert);
					r2.reverseComplementFast();
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
				
				if(trueSize<initialLen1+initialLen2 && trueSize>0){unmergedOverlappingCount++;}
			}
			
			r1.setInsert(bestInsert);
			
			if(findAdapterSequence && bestInsert>0 && bestInsert<r1.length()){
				storeAdapterSequence(r1, bestInsert);
				//r2.reverseComplementFast();
				storeAdapterSequence(r2, bestInsert);
				//outstream.println(new String(r2.bases));
				//r2.reverseComplementFast();
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
			
			if(parseCustom){

				final boolean fp=(bestInsert>0 && bestInsert!=trueSize);
				final boolean fn=(bestInsert<1 && trueOverlap>0);
				final boolean tp=(bestInsert>0 && bestInsert==trueSize);
				final boolean tn=(bestInsert<1 && trueOverlap<1);
				
				assert((fp ? 1 : 0)+(fn ? 1 : 0)+(tp ? 1 : 0)+(tn ? 1 : 0)==1) : 
					fp+", "+fn+", "+tp+", "+tn+", "+bestInsert+", "+trueOverlap+", "+trueSize+", "+(initialLen1+initialLen2);
				final boolean correct=(bestInsert==trueSize);
				if(MAKE_VECTOR){
					FloatList fl=r2.fetchVector();
					if(fl!=null && fl.size==23){fl.add((tp || fn) ? 1 : 0);}
//					assert(false) : fl;
				}
				if(TAG_CUSTOM){
					ByteBuilder bb=r1.fetchBB();
					//				bb.append("_tr").append(correct ? 1 : 0);

//					if(correct) {//truth
//						bb.append(1);
//					}else {
//						bb.append(0);
//					}
					if(tp || fn){//truth
						bb.append(1);
					}else{
						bb.append(0);
					}

					if(trueOverlap>0 && trueOverlap<8 && EXEMPT_SHORT_OVERLAPS_FROM_STATS){
						r1.nullifyObj();
					}
				}

				if(ONLY_OUTPUT_INCORRECT && correct){
					if(TAG_CUSTOM){r1.nullifyObj();}
					return bestInsert;
				}
				
//				false, 231, 58, 231, 289
				assert(fn==(trueSize<initialLen1+initialLen2 && trueSize>0 && bestInsert<=0)) : 
					fn+", "+bestInsert+", "+trueOverlap+", "+trueSize+", "+(initialLen1+initialLen2);
				
				if((OUTPUT_FALSE_POSITIVE && fp) || (OUTPUT_FALSE_NEGATIVE && fn) || 
						(OUTPUT_TRUE_POSITIVE && tp) || (OUTPUT_TRUE_NEGATIVE && tn)){
					//Print
				}else{
					if(TAG_CUSTOM){r1.nullifyObj();}
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
		
		/** Validates reads and applies configured normalization, trimming and prefilters.
		 * Counts pairs and their bases after validation but before trim/filter/extension.
		 * Null R2 is rejected here; outer callers may require it earlier. Changes are not
		 * rolled back on rejection. The legacy extend1 orientation limitation remains
		 * documented below; this helper does not establish safe extension behavior.
		 * @param r1 Nonnull first read
		 * @param r2 Mate, or null for this helper's RET_BAD path
		 * @param qtrim Whether to request quality trimming with parameter index zero
		 * @return 1 on passing prechecks, or a negative RET code
		 */
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
				//TODO: Probable bug - r2 is still in sequencer orientation here,
				//but extendAndMerge assumes the genomic-forward orientation supplied
				//by processReadPair's later extend2 call. This extends R2's wrong end.
				//Hollow rejects extend= until this legacy caller has a truth fixture.
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
		
		/** Extends each mate's sequencer-right end and optionally retries aligned search.
		 * R2 enters aligned with R1 and is temporarily reversed for extension. A mate
		 * continues only while it has received every requested base so far. Extensions
		 * remain on failure; no un-extension is implemented. Tallies two mate outcomes
		 * per call, not one extension attempt per loop iteration.
		 * @param r1 Owned first read in forward orientation
		 * @param r2 Owned mate already reverse-complemented into R1's orientation
		 * @param amt Requested additional bases per mate per iteration
		 * @param iters Positive retry count; must be one when merge is false
		 * @param merge Whether to search after a successful extension iteration
		 * @param extendRvec Writable array of at least two slots receiving appended bases per mate
		 * @return Last search result, or initial RET_AMBIG when no search is requested/performed
		 */
		private final int extendAndMerge(Read r1, Read r2, int amt, int iters, boolean merge, int[] extendRvec){
			assert(iters>0);
			assert(merge || iters==1);
			
			extendRvec[0]=extendRvec[1]=0;
			
			int sum1=0, sum2=0, attempted=0;
			int bestInsert=RET_AMBIG;
			for(int i=0; i<iters && (bestInsert==RET_AMBIG || bestInsert==RET_NO_SOLUTION); i++){
				
				int e1=(sum1==attempted ? extendRead(r1, amt) : 0);
				r2.reverseComplementFast();
				int e2=(sum2==attempted ? extendRead(r2, amt) : 0);
				r2.reverseComplementFast();
				
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
		
		/** Appends graph extension in the read's current forward orientation using local scratch.
		 * A positive extension replaces bases and copies present qualities, filling added
		 * quality positions with this worker's numeric qfake. Null qualities remain null.
		 * @param r Owned read with nonnull bases
		 * @param amt Requested extension limit in bases
		 * @return Number of appended bases
		 */
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
		
		/** Runs ratio search with optional quality weighting and an unweighted fallback.
		 * Uses current global ratio/insert limits and local probability/result scratch.
		 * Fallback may replace the weighted candidate. The helper's mismatch/ambiguity
		 * values remain in rvector[2]/[4]; a returned candidate can still be ambiguous.
		 * @param r1 Forward read
		 * @param r2 Mate already reverse-complemented into R1's orientation
		 * @param minOverlap Effective overlap support threshold before ratio reduction
		 * @return Candidate insert, or a negative search rejection
		 */
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
		
		/** Runs mismatch-based overlap retries and records candidate/mismatch/ambiguity scratch.
		 * Missing mate qualities restrict quality iteration and disable failure trimming.
		 * With qualities, configured retries may trim R1's right/R2's aligned left end;
		 * these mutations are not undone if search still fails. Loose mode has an extra retry.
		 * @param r1 Owned forward read
		 * @param r2 Owned mate already reverse-complemented into R1's orientation
		 * @param minOverlap Required overlap support before retry adjustments
		 * @return Candidate insert, or a negative search rejection
		 */
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
		
		/** Derives the search overlap threshold with reusable entropy-count scratch.
		 * Disabled entropy uses the global minimum. Normal mode examines R1's tail/R2's
		 * aligned head; loose mode examines both ends and applies its expected-error factor.
		 * @param r1 Forward read
		 * @param r2 Mate already reverse-complemented into R1's orientation
		 * @return Effective overlap support threshold in bases
		 */
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
		
		/** Dormant lowercase-tail heuristic; its only search call is commented out.
		 * Compares the remaining prefixes with uppercase N ignored, without mutating reads.
		 * Success writes rvector[0]/[2]/[4]; failure need not refresh those slots.
		 * @param r1 Nonempty forward read ending in lowercase adapter bases
		 * @param r2 Nonempty aligned mate beginning with lowercase adapter bases
		 * @return Shared prefix length, possibly zero, or -1 when the heuristic rejects
		 */
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
		
		/** Selects ratio/flat candidates, applies filters and optionally evaluates the network.
		 * Reuses/grows probability scratch and may populate read attachments. R1 must not
		 * carry an attachment; R2's attachment must be null or a FloatList. Candidate
		 * collection modes retain some ordinarily rejected results, and retries replace
		 * their feature vector. r1.errors is candidate mismatch scratch, not a promised
		 * final-output mismatch count. No training provenance or accuracy is established here.
		 * @param r1 Owned forward read
		 * @param r2 Owned mate already reverse-complemented into R1's orientation
		 * @return Positive candidate insert or a negative RET rejection/ambiguity code
		 */
		private final int mateByOverlap(Read r1, Read r2){
			final int len1=r1.length(), len2=r2.length();
			assert(r1.fetchObject()==null);
			final ByteBuilder bb=(TAG_CUSTOM ? new ByteBuilder(32) : null);
			if(bb!=null){r1.setObj(bb);}
			if(MAKE_VECTOR){
				if(r2.obj==null){r2.setObj(new FloatList(24));}
				else{((FloatList)r2.obj).clear();}//FIXED: Retries after correction or extension must replace the previous candidate vector.
			}
			final FloatList vector=(FloatList) r2.obj;
			
			final int minOverlap=calcMinOverlapFromEntropy(r1, r2);
			if(verbose){outstream.println("\nminOverlap: "+minOverlap);}
			if(TAG_CUSTOM || MAKE_VECTOR){
				int maxBasesToConsider=Tools.min(Tools.max(len1, len2), len1+len2-minInsert);
				float r1ee=r1.expectedTipErrors(false, maxBasesToConsider);
				float r2ee=r2.expectedTipErrors(false, maxBasesToConsider);
				if(TAG_CUSTOM){
					bb.append(minOverlap*0.1f, 1).tab();
					bb.append(Tools.format("%.4f\t%.4f\t", r1ee, r2ee));
				}
				if(MAKE_VECTOR){
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
				}else if(MAKE_VECTOR){
					vector.add(0.1f);//For some reason this gets left off occasionally...
				}
				if(verbose){System.err.println("D: "+ambig+", "+bestInsert);}
			}
			if(verbose){System.err.println("E: "+ambig+", "+bestInsert);}
			
			//JGI-086: Apply the same ambiguity policy after quality filtering as before it.
			//Custom/vector collection keeps its candidate; ordinary NN calls reach the network decision below.
			if(ambig && ((!TAG_CUSTOM && !MAKE_VECTOR) || ONLY_OUTPUT_FALSE_POSITIVE)){return RET_AMBIG;}
			r1.errors=bestBad;
			
			if(MAKE_VECTOR && bestInsert>0 && r2.obj()!=null){
				FloatList list=r2.fetchVector();
				if(list.size!=23){
					System.err.println("\n"+r1.numericID+", "+r1.obj());
					System.err.println("*"+list.size+", "+parseCustom+", "+r1.numericID+", "+list);
					assert(false);
				}else{
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
		
		/** Obtains the supplied truth insert from a recognized header or mapped coordinates.
		 * Mapping fallback sets both mapped flags before estimating and asserts a positive
		 * result. This is interpretation of supplied labels, not independent biological truth.
		 * @param r1 First read with a nonnull identifier
		 * @param r2 Linked mate used by the mapping fallback
		 * @return Header parser's insert or the asserted-positive mapping estimate
		 */
		private final int parseInsert(Read r1, Read r2){
			int trueSize=-1;
			if(r1.id.startsWith("SYN") || r1.id.startsWith("insert=") ||
					r1.id.startsWith("f_") || r1.id.indexOf("\tf_")>=0){
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
		
		/** Applies this worker's two adapter lists; a good comparison wins over bad/unknown.
		 * Without a good result, bad rejects and unknown rejects only in strict mode.
		 * Uses three informative bases and a 60% match threshold; if no comparison is made,
		 * the pair passes without positive adapter evidence. Temporarily reverses R2 for
		 * its tail check and restores it on normal completion, not guaranteed on exception.
		 * Counts eligible mate checks in adaptersExpectedT, including an empty list, and
		 * at most one success per call in adaptersFoundT; retries count again.
		 * @param r1 Owned first read in forward orientation
		 * @param r2 Owned mate already reverse-complemented into the first read's orientation
		 * @param bestInsert Nonnegative insert length; adapter tails begin at this index
		 * @return Whether the rejection policy accepts, including when no comparison is made
		 */
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
						r2.reverseComplementFast();
						for(byte[] adapter : adapterList2){
							int aiv=adapterIsValid(r2, adapter, bestInsert, minAdapterOverlap, minAdapterRatio);
							if(aiv==2){good++; break;}
							else if(aiv==1){invalid++;}
							else{bad++;}
						}
						r2.reverseComplementFast();
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
		
		/** Runs preprocessing, overlap search and configured correction/extension/filter retries.
		 * Accepts raw sequencer orientation, aligns R2 internally and restores its orientation
		 * on normal return. Base/quality arrays, errors and other Read state may change;
		 * exceptions do not guarantee orientation/quality restoration or rollback. This is
		 * search processing, not final output joining/submission. Collection modes retain
		 * their existing exceptions to ordinary rejection rules.
		 * @param r1 Owned input root with its linked mate
		 * @param r2 Nonnull mate in raw sequencer orientation
		 * @return Positive candidate insert or a negative RET code
		 */
		private final int processReadPair(final Read r1, final Read r2){
			final int initialLen1=r1.length(), initialLen2=r2.length();
			{
				final int x=preprocess(r1, r2, (qtrim1 && !qtrim2));
				if(!TAG_CUSTOM && !MAKE_VECTOR && x<0){return x;}
			}
			
			r2.reverseComplementFast();
			
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
					//				r2.reverseComplementFast();
					//				qtrim(r1, r2);
					//				r2.reverseComplementFast();

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
					r2.reverseComplementFast();
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

				//Target-directed graph bridge.  The reads are both in forward genomic orientation here.
				boolean bridged=false;
				if(kmerBridge && (bestInsert==RET_AMBIG || bestInsert==RET_NO_SOLUTION)){
					bridgeAttemptsT++;
					final int bridgeInsert=bridgeT.find(r1, r2);
					if(bridgeInsert>0){
						final int rawOverlap=r1.length()+r2.length()-bridgeInsert;
						if(bridgeInsert>=minInsert && bridgeInsert<=maxReadLength && rawOverlap<=bridgeMaxOverlap){
							bridgeT.apply(r1);
							bestInsert=bridgeInsert;
							bridged=true;
							bridgeSuccessT++;
						}else{bridgeRejectedT++;}
					}else{
						final int status=bridgeT.status();
						if(status==TargetKmerBridge.NO_PATH){bridgeNoPathT++;}
						else if(status==TargetKmerBridge.NO_TARGET){bridgeNoTargetT++;}
						else if(status==TargetKmerBridge.AMBIGUOUS){bridgeAmbiguousT++;}
						else if(status==TargetKmerBridge.CYCLE){bridgeCycleT++;}
						else if(status==TargetKmerBridge.NONRECIPROCAL){bridgeNonreciprocalT++;}
						else{throw new RuntimeException("Unknown kmer bridge status "+status);}
					}
				}
				
				//Extension
				if(!bridged && extendRight2>0 && (requireExtensionMatch || bestInsert==RET_AMBIG || bestInsert==RET_NO_SOLUTION)
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
			r2.reverseComplementFast();
			assert(r1.quality==null || (r1.quality.length==r1.bases.length && r2.quality.length==r2.bases.length));
			return bestInsert;
		}
		
		/** Restores saved quality references, copying them only when current lengths differ.
		 * Preserves null qualities. R2 enters aligned; its resize path reverses it around
		 * the copy and restores that orientation normally. The existing extension-padding
		 * TODO remains unresolved; this helper does not define a corrected padding policy.
		 * @param r1 Owned first read
		 * @param r2 Owned aligned mate
		 * @param qual1 Saved R1 numeric qualities, or null
		 * @param qual2 Saved R2 numeric qualities in the pair path's captured orientation, or null
		 */
		private void restoreQualities(Read r1, Read r2, byte[] qual1, byte[] qual2){
			//TODO: Probable bug - extension fill loops start at the new length,
			//so appended qualities stay zero; intended extension padding is unverified.
			//Absent input qualities are preserved without creating synthetic arrays.
			final int len1=r1.length(), len2=r2.length();
			r1.quality=qual1;
			r2.quality=qual2;
			final byte qf=Shared.FAKE_QUAL;
			if(qual1!=null && qual1.length!=len1){
				r1.quality=Arrays.copyOf(qual1, len1);
				for(int i=len1; i<r1.quality.length; i++){
					r1.quality[i]=qf;
				}
			}
			if(qual2!=null && qual2.length!=len2){
				r2.reverseComplementFast();
				r2.quality=Arrays.copyOf(qual2, len2);
				for(int i=len2; i<r2.quality.length; i++){
					r2.quality[i]=qf;
				}
				r2.reverseComplementFast();
			}
		}
		
		/** Uses a stored mapped insert when eligible, otherwise the configured overlap search.
		 * Applies insert-length bounds to normal results; custom collection retains its
		 * early-return behavior. May mutate search scratch, attachments and candidate errors.
		 * @param r1 Owned forward read
		 * @param r2 Owned mate already reverse-complemented into R1's orientation
		 * @return Candidate insert or a RET code, including short/long/ambiguous rejection
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
			
			if((TAG_CUSTOM) && !ONLY_OUTPUT_FALSE_POSITIVE){return bestInsert;}
			if(ambig && !MAKE_VECTOR){return RET_AMBIG;}
			else if(bestInsert>0){
				if(bestInsert<minInsert){return RET_SHORT;}
				else if(bestInsert>maxReadLength){return RET_LONG;}
				return bestInsert;
			}
			else{return RET_NO_SOLUTION;}
		}
		
		/** Delegates suffix counting to this worker's tracker and its current PhiX policy.
		 * The older local PhiX filter below is disabled, not an additional active check.
		 * @param r Read with bases, linked mate bases and pair number zero or one
		 * @param insert Nonnegative suffix boundary in read bases
		 */
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
		
		/** Dormant local PhiX diagnostic wrapper; active collection uses AdapterTracker.
		 * @param r Read with nonnull bases and mate bases
		 * @param insert Nonnegative suffix boundary in both mates
		 * @return Whether either suffix matches the local prefix heuristic
		 */
		private boolean looksLikePhix(Read r, int insert){
			return looksLikePhix(r.bases, insert) || looksLikePhix(r.mate.bases, insert);
		}
		
		/** Dormant local prefix diagnostic with uppercase N treated as a wildcard.
		 * A match prints the suffix and prefix to outstream; it is not a biological PhiX test.
		 * @param bases Nonnull sequence in the prefix's orientation
		 * @param insert Nonnegative suffix boundary
		 * @return Whether the suffix is long enough and matches the whole local prefix
		 */
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
		/** Local adapter suffix counts, merged by the owner after this worker joins. */
		final AdapterTracker atrackT=new AdapterTracker();

		/** Per-thread cardinality tracker; merged into the shared loglog after this thread joins.
		 * The shared loglog is not thread-safe, so each thread accumulates into its own tracker with
		 * the same params (same seed) and they are combined at the end. */
		final CardinalityTracker loglogT;

		/** Local scratch for graph merge validation. */
		private final BitSet mergeOKBitsetT=new BitSet(400);
		/** Local long-k graph scratch, absent without a graph or when k is below 32. */
		private final Kmer kmerT;
		/** Local bridge search state borrowing the enclosing graph, or null. */
		private final TargetKmerBridge bridgeT;
		/** Reused graph-validation counts. */
		private IntList countList=new IntList();
		/** Reused Bloom merge-validation k-mers. */
		private LongList kmers=new LongList();
		
		/** Numeric quality assigned by extendRead to newly appended quality positions. */
		final byte qfake=Shared.FAKE_QUAL;
		
		/** Overlap scratch: index 2 is candidate mismatch/rejection, 4 is ambiguity.
		 * Index 0 is written by flat/dormant-adapter paths, not every ratio search. */
		private final int[] rvector=KillSwitch.allocInt1D(5);
		/** Appended-base counts for R1/R2 from the latest extension call. */
		private final int[] extendRvec=KillSwitch.allocInt1D(2);

		/** Four-base extension-count scratch passed to the graph service. */
		private final int[] rightCounts=KillSwitch.allocInt1D(4);
		/** Optional left-count scratch; omitted when left junctions may be traversed. */
		private final int[] leftCounts=(extendThroughLeftJunctions ? null : new int[4]);
		
		/** Reused extension builder; successful extension installs its copied bytes. */
		private final ByteBuilder bbT=new ByteBuilder();
		
		/** Positive candidate-insert counts before output selection; the last bin includes larger lengths. */
		final long[] hist=new long[histlen];
		/** Reusable entropy k-mer counts, absent when entropy filtering is disabled. */
		final short[] kmerCounts;
		
		/** Per-base quality-probability scratch, allocated/grown by overlap search. */
		private float[] aprob, bprob;

		/** Corrections counted by the post-search output-correction path. */
		long errorsCorrectedT=0;
		
		/** Paired roots counted after validation and before trimming/filtering. */
		long pairsProcessed=0;
		/** Sum of paired input lengths at the same pretrim counting point. */
		long basesProcessed=0;
		/** Legacy merge counter with no increment in the current worker paths. */
		long matedCount=0;
		/** Accepted inserts agreeing with the supplied truth size, subject to exemptions. */
		long correctCount=0;
		/** Rejected pairs whose supplied truth insert indicates an overlap. */
		long unmergedOverlappingCount=0;
		/** Rejection tallies for the corresponding pair-search RET results. */
		long ambiguousCount=0;
		long tooShortCount=0;
		long tooLongCount=0;
		/** Accepted inserts differing from the supplied truth size, not an accuracy oracle. */
		long incorrectCount=0;
		/** Pair-search no-solution rejections. */
		long noSolutionCount=0;
		/** Insert-length sums in bases for the matching/mismatching truth tallies. */
		long insertSumCorrect=0;
		long insertSumIncorrect=0;
		/** Observed accepted insert extrema; minimum retains its sentinel until observed. */
		int insertMax=0;
		int insertMin=999999999;

		/** Mate outcomes from extension calls, not counts of individual extension iterations. */
		long fullyExtendedT=0;
		long partlyExtendedT=0;
		long notExtendedT=0;
		long extensionsAttemptedT=0;
		/** Per-pair bridge attempts and classified bridge outcomes. */
		long bridgeAttemptsT=0, bridgeSuccessT=0;
		long bridgeNoPathT=0, bridgeNoTargetT=0, bridgeAmbiguousT=0, bridgeCycleT=0;
		long bridgeNonreciprocalT=0, bridgeRejectedT=0;
		
		/** Eligible mate adapter checks and successful calls; repeated checks count again. */
		long adaptersExpectedT=0;
		long adaptersFoundT=0;
		
		/** Borrowed input and optional output streams; lifecycle belongs to runPhase. */
		private final Streamer cris;
		private final Writer rosgood;
		private final Writer rosbad;
		private final ByteStreamWriter rosi;
		
		/** Per-worker snapshots of the owner's output-processing choices. */
		private final boolean joinReads;
		private final boolean trimReadsByOverlap;
		/** Captured loop failure, inspected by the owner after joining this worker. */
		private Throwable failure=null;
		
		/** Network parameter/working-state copy; some metadata remains shared. */
		private final CellNet net;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Clears the shared aggregate histogram, counts and insert sums, including false
	 * negatives stored in unmergedOverlappingTotal. Restores the unobserved insert
	 * extrema to minimum 999999999 and maximum zero; configuration and worker-local
	 * counters are unchanged. Not synchronized: callers coordinate access to shared totals.
	 */
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
		unmergedOverlappingTotal=0;
		insertSumCorrectTotal=0;
		insertSumIncorrectTotal=0;
		fullyExtendedTotal=0;
		partlyExtendedTotal=0;
		notExtendedTotal=0;
		extensionsAttempted=0;
		bridgeAttemptsTotal=bridgeSuccessTotal=0;
		bridgeNoPathTotal=bridgeNoTargetTotal=bridgeAmbiguousTotal=bridgeCycleTotal=0;
		bridgeNonreciprocalTotal=bridgeRejectedTotal=0;
		adaptersExpected=0;
		adaptersFound=0;

		insertMinTotal=999999999;
		insertMaxTotal=0;
	}
	
	/** Dormant pathname-cached loader; BBMerge has no active call to this helper.
	 * Null clears the path/model but leaves the cutoff. An already loaded equal pathname
	 * skips parsing without checking file freshness. The path is assigned before parsing;
	 * errors or an absent-file null do not promise rollback or a clean disabled state.
	 * A loaded model supplies its cutoff unless an explicit cutoff was selected. Method
	 * synchronization does not make other shared model/configuration access synchronized.
	 * @param fname Model path, or null to clear this helper's path/model state
	 */
	private static synchronized void loadNet(String fname){
		if(fname==null){
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
	
	/** First input file path */
	private String in1;
	/** Second input file path */
	private String in2;
	
	/** Additional paths forwarded to graph/filter construction, not runPhase's read stream. */
	private ArrayList<String> extra=new ArrayList<String>();

	/** First output file for merged/good reads */
	private String out1=null;
	/** Second output file for merged/good reads */
	private String out2=null;
	/** First output file for unmerged/bad reads */
	private String outb1=null;
	/** Second output file for unmerged/bad reads */
	private String outb2=null;
	/** Resolved single-file normal insert-report path; null disables report work. */
	private String outinsert=null;
	/** Output file for insert size histogram */
	private String ihist=null;
	/** Output file for adapter sequences */
	private String outAdapter=null;
	/** Output file for cardinality estimation */
	private String outCardinality=null;

	/** List of R1 adapters to look for to verify a short-insert overlap */
	private ArrayList<byte[]> adapterList1;
	/** List of R2 adapters to look for to verify a short-insert overlap */
	private ArrayList<byte[]> adapterList2;
	/** Require short inserts to have adapter sequence at the expected location */
	private boolean verifyAdapters=false;
	
	/** Owner cardinality aggregate; workers hash through separate trackers and merge after join. */
	final CardinalityTracker loglog;
	
	/** Input limit forwarded unchanged to the selected reader; units depend on its contract. */
	private long maxReads=-1;
	/** Whether to join overlapping reads into single sequences */
	private boolean join=true;
	/** Correct accepted overlapping mates for unjoined output, independently of graph correction. */
	private boolean ecco=false;
	/** Whether to trim reads based on overlap detection */
	private boolean trimByOverlap=false;
	
	/** Positive overlap-probability rejection threshold; nonpositive disables this filter. */
	private float pfilterRatio=0.00004f;
	/** Multiplier on expected mismatches plus efilterOffset for the mismatch allowance. */
	private float efilterRatio=6f;
	/** Offset in expected-mismatch units, added before multiplication by efilterRatio. */
	private float efilterOffset=0.05f;
	/** Whether to use expected error filtering */
	private boolean useEfilter=true;
	/** Whether to use maximum expected errors filtering */
	private boolean useMEEfilter=false;
	
	/** Requested output ordering; two-file paired outputs force ordering in runPhase. */
	private boolean ordered=false;
	/** Whether ratio search first requests quality-weighted scoring when qualities are present. */
	private boolean overlapUsingQuality=false;
	/** Whether an unsuccessful/ambiguous weighted ratio search may retry unweighted. */
	private boolean overlapWithoutQuality=true;
	/** Enables graph coverage validation of candidate joins. */
	private boolean useKFilter=false;
	/** Coverage cutoff passed to the graph k-filter checks. */
	private int filterCutoff=1;
	/** Requested graph/Bloom k-mer length; constructor validation chooses the usable service. */
	private int kmerLength=31;
	/** Graph preallocation and prefilter settings forwarded to Tadpole. */
	private boolean prealloc=false;
	private int prefilter=0;
	/** Positive graph filter-memory override; zero leaves the delegated default. */
	private long filterMemoryOverride=0;
	/** Graph-correction strategy switches forwarded to Tadpole. */
	private boolean eccTail=false;
	private boolean eccPincer=false;
	private boolean eccReassemble=true;

	/** Enables sequence-complexity-derived minimum overlap using local worker count scratch. */
	private boolean useEntropy=true;
	/** Entropy k-mer width and target score; worker counts have 4^entropyK slots. */
	private int entropyK=3;
	private int minEntropyScore=39;//30 loose;//39 normal;//44 strict;
	
	/** Requested reader sampling seed/fraction; algorithms and negative-seed handling follow the reader. */
	private long sampleseed=-1;
	private float samplerate=1;
	
	/** Instance aggregate of workers' post-search output-correction counts. */
	private long errorsCorrectedTotal=0;
	
	/** Enables collecting adapter suffixes for consensus output. */
	private boolean findAdapterSequence=false;
	/** CLI mirror; active collection uses the corresponding AdapterTracker static setting. */
	private boolean ignorePhixAdapters=false;
	
	/** Owner adapter aggregate, merged from workers after join. */
	private final AdapterTracker atrack=new AdapterTracker();
	
	/** Shared preloaded graph/count service, or null; worker scratch is separate. */
	private final Tadpole tadpole;
	/** Requested pre-search/after-search extension bases per mate and retry count. */
	private int extendRight1=0;
	private int extendRight2=0;
	private int extendIterations=1;
	/** Enables graph error-correction retries after ambiguous/no-solution overlap results. */
	private boolean eccTadpole=false;
	/** Experimental endpoint-only counts, restricted to read extension. */
	private boolean hollow=false;
	/** Reject alternate-span conflicts before appending an extension base. */
	private boolean hollowBranchStop=false;
	/** Alternative-span branch thresholds forwarded to the Hollow graph service. */
	private int hollowBranchRatio=4, hollowBranchSupport=3;
	/** Optional graph-cleaning requests; these require explicit k-mer keys. */
	private boolean shave=false;
	private boolean rinse=false;
	/** Enables endpoint graph bridging after an ambiguous/no-solution overlap. */
	private boolean kmerBridge=false;
	/** Distinguish explicitly supplied k/bridge settings from constructor defaults. */
	private boolean setKmerLength=false, setKmerBridge=false;
	/** Maximum bridge walk distance in appended bases. */
	private int bridgeDistance=0;
	/** Reserved target-window request; the current endpoint bridge requires zero. */
	private int bridgeTargetWindow=0;
	/** Reserved minimum-hit request; the current endpoint bridge requires one. */
	private int bridgeMinHits=1;
	/** Maximum allowed length-sum minus inferred insert for an accepted bridge, in bases. */
	private int bridgeMaxOverlap=0;
	/** Whether the bridge must also agree when walked in the reverse direction. */
	private boolean bridgeReciprocal=true;
	/**
	 * Parses hash mode aliases: null/auto=-1, fixed/fingerprint/1=2, pair/two/2=1.
	 * Other values use Parse.parseBoolean, yielding pair mode (1) or explicit keys (0).
	 * @param b Requested mode or permissive boolean value
	 * @return -1 for automatic selection, 0 for explicit keys, 1 for hash pairs, or 2 for fixed fingerprints
	 */
	private static int parseHashMode(String b){
		if(b==null || b.equalsIgnoreCase("auto")){return -1;}
		if(b.equalsIgnoreCase("fixed") || b.equalsIgnoreCase("fingerprint") || b.equals("1")){return 2;}
		if(b.equalsIgnoreCase("pair") || b.equalsIgnoreCase("two") || b.equals("2")){return 1;}
		return Parse.parseBoolean(b) ? 1 : 0;
	}

	/**
	 * Resolves the configured mode without changing it. Hollow lookup forces hash pairs first.
	 * Otherwise graph cleaning requires explicit keys and rejects a positive requested mode;
	 * k below 64 forces explicit keys, then a nonnegative requested mode takes precedence.
	 * Automatic mode uses explicit keys at k=64 without preallocation, fixed fingerprints
	 * with preallocation, and hash pairs otherwise.
	 * @return 0 for explicit keys, 1 for hash pairs, or 2 for fixed fingerprints
	 */
	private int effectiveHashMode(){
		//Hollow lookup requires both full hashcodes, even when prealloc=t.
		if(hollow){return 1;}
		if(shave || rinse){
			if(hashKmerMode>0){throw new RuntimeException("hashkmers requires shave=f and rinse=f because graph cleaning enumerates kmers.");}
			return 0;
		}
		if(kmerLength<64){return 0;}
		if(hashKmerMode>=0){return hashKmerMode;}
		if(kmerLength==64 && !prealloc){return 0;}
		return prealloc ? 2 : 1;
	}

	/** -1=auto, 0=explicit, 1=resizable hash pair, 2=fixed xor2 fingerprint. */
	private int hashKmerMode=-1;

	/** Shared Bloom count/filter service, absent when the graph or neither service is selected. */
	private final BloomFilter bloomFilter;
	/** Shared Bloom correction/merge-check helper; this field is not worker-local state. */
	private final BloomFilterCorrector corrector;
	
	/** Requests a full graph/count table even without extension or correction. */
	private boolean forceExactKmerCounts=false;
	/** Validate sufficiently long candidate joins with the selected graph/Bloom service. */
	private boolean testMerge=false;
	/** Enables Bloom correction retries after ambiguous/no-solution overlap results. */
	private boolean eccBloom=false;
	/** Bloom construction bit depth and hash count. */
	private int bloomBits=2;
	private int bloomHashes=3;
	/** Width, multiplier and threshold forwarded to candidate merge validation. */
	int testMergeWidth=4;
	long testMergeMult=80L;
	int testMergeThresh=3;
	
	/** Reconcile original and extended insert candidates under the existing extension policy. */
	private boolean requireExtensionMatch=false;
	/** Requires a positive original candidate before the extend2 reconciliation attempt. */
	private boolean requireStrictExtensionMatch=false;
	/** Bases subtracted from the original pair-length sum for the approximate-overlap bound. */
	private int minApproxOverlapRem=26;
	
	/** Graph extension policy; true also omits worker left-count scratch. */
	private boolean extendThroughLeftJunctions=true;
	/** Minimum graph counts requested for seeding and extending. */
	private int minCountSeed=3, minCountExtend=2;
	/** Branch multipliers forwarded unchanged to graph construction. */
	private float branchMult1=20;
	private float branchMult2=3;
	/** Quality probability threshold forwarded to graph construction. */
	private float minProb=0.5f;
	/** Branch-count constant forwarded to graph construction. */
	private int branchLowerConst=3;
	
	/** Enables overlap search when the stored-mapping shortcut is unavailable. */
	private boolean MATE_BY_OVERLAP=true;
	/** Sends rejected roots to good output instead of the separate bad list. */
	private boolean MIX_BAD_AND_GOOD=false;
	/** Custom-label output selection; applied after statistics and candidate processing. */
	private boolean ONLY_OUTPUT_INCORRECT=false;
	
	/** Which supplied-truth categories survive custom-label output selection. */
	private boolean OUTPUT_FALSE_POSITIVE=true;
	private boolean OUTPUT_FALSE_NEGATIVE=true;
	private boolean OUTPUT_TRUE_POSITIVE=true;
	private boolean OUTPUT_TRUE_NEGATIVE=true;
	
	/** Quantize good/bad Read qualities immediately before their writer submissions. */
	private boolean quantizeQuality=false;
	
	/*--------------------------------------------------------------*/
	
	/** Per-calling-thread result scratch for the static overlap helper, not worker scratch. */
	private static ThreadLocal<int[]> localRvector=new ThreadLocal<int[]>();
	
	/** Shared error indication combined during cleanup; resetCounters does not clear it. */
	static boolean errorState=false;

	/** Print expanded arguments when the pre-parser replaces the original argument array. */
	private static boolean showFullArgs=true;

	/** Recalibrate quality scores using matrices */
	static boolean recalibrateQuality=false;
	/** When false, strips qualities after preprocessing and restores them on normal completion.
	 * qtrim2 retries restore saved qualities before trimming/searching again; those
	 * retry searches can therefore see qualities even when this setting is false.
	 */
	static boolean useQuality=true;
	/** Requested quality-trim sides; retries on aligned R2 exchange left and right. */
	static boolean qtrimRight=false;
	static boolean qtrimLeft=false;
//	static boolean untrim=false;
	/** Numeric Phred thresholds and corresponding error probabilities indexed by trim iteration. */
	static float[] trimq=new float[]{6};
	static float[] trimE=QualityTools.phredToProbError(trimq);
	/** Average numeric-quality prefilter threshold; nonpositive disables it. */
	static float minAvgQuality=0;
	/** Base-limit argument passed to Read.avgQuality for the average-quality precheck. */
	static int minAvgQualityBases=0;
	/** Expected-tip-error ceiling used when the MEE prefilter and qualities are enabled. */
	static float maxExpectedErrors=0;
	/** Both-short length threshold; a separate guard handles mates shorter than two bases. */
	static int minReadLength=1;
	/** Candidate insert upper bound inherited from parser max length; constructor resolves unlimited. */
	static int maxReadLength=-1;
	/** Candidate insert lower bound in bases, raised to the overlap minimum by construction. */
	static int minInsert=15;
	/** Detailed scan's insert start; a nonpositive request is resolved by construction. */
	static int minInsert0=-1;
	/** Initial versus retry trimming requests; qtrim2 suppresses the first qtrim pass. */
	static boolean qtrim1=false;
	static boolean qtrim2=false;
	/** Flat-search failure-trim retry limit; unavailable mate qualities disable these retries. */
	static int TRIM_ON_OVERLAP_FAILURE=1;
	/** Flat-search quality rounds; absent mate qualities restrict this to one round. */
	static int QUAL_ITERS=3;
	/** Interpret supplied header/mapping labels for truth statistics and category selection. */
	static boolean parseCustom=false;
	/** Maximum length passed to non-built-in adapter-list parsing. */
	static int maxAdapterLength=21;
	
	/** Requests accepted-output cropping in the join/output-correction branches. */
	static boolean trimNonOverlapping=false;
	
	/** Adapter-consensus trimming mirror, forwarded to AdapterTracker.trimPolyAorG. */
	static boolean trimPolyA=true;
	
	/** Positive left boundary to retain during forced trimming, in zero-based base coordinates. */
	static int forceTrimLeft;
	/** Positive inclusive right coordinate to retain; nonpositive values do not select it. */
	static int forceTrimRight;
	/** Positive number of rightmost bases to remove during forced trimming. */
	static int forceTrimRight2;
	/** Trim right bases of the read modulo this value.
	 * e.g. forceTrimModulo=50 would trim the last 3bp from a 153bp read. */
	static int forceTrimModulo;
	
	/** Strict preset selectors; the argument prepass also changes associated thresholds. */
	public static boolean strict=false;
	static boolean vstrict=false;
	static boolean ustrict=false;
	static boolean xstrict=false;
	/** Loose/speed preset selectors, not independent guarantees about accepted overlaps. */
	static boolean loose=false;
	static boolean vloose=false;
	static boolean uloose=false;
	static boolean xloose=false;
	static boolean fast=false;
	
	/** Retained switch for the dormant lowercase-adapter heuristic. */
	static boolean lowercaseAdapters=false;
	
	/** Number of insert-count bins; the final bin accumulates larger accepted inserts too. */
	private static final int histlen=2000;
	/** Shared sum of worker insert histograms, reset for phase aggregation. */
	static long[] histTotal=new long[histlen];
	/** Histogram reporting width in bases; does not resize histTotal. */
	static int bin=1;

	/** Shared phase totals of pretrim paired roots and paired base lengths. */
	static long readsProcessedTotal=0;
	static long basesProcessedTotal=0;
	/** Legacy sum of the currently unincremented worker matedCount field. */
	static long matedCountTotal=0;
	/** Accepted inserts matching supplied truth labels; not independent accuracy validation. */
	static long correctCountTotal=0;
	/** Rejected pairs with supplied truth overlap; resetCounters clears this false-negative tally. */
	static long unmergedOverlappingTotal=0;
	/** Aggregated pair-search rejection tallies. */
	static long ambiguousCountTotal=0;
	static long tooShortCountTotal=0;
	static long tooLongCountTotal=0;
	/** Accepted inserts differing from supplied truth labels. */
	static long incorrectCountTotal=0;
	/** Aggregated no-solution rejections. */
	static long noSolutionCountTotal=0;
	/** Accepted insert-length sums in bases, split by supplied truth agreement. */
	static long insertSumCorrectTotal=0;
	static long insertSumIncorrectTotal=0;
	/** Aggregated mate outcomes/attempts from extension calls, not per-iteration counts. */
	static long fullyExtendedTotal=0;
	static long partlyExtendedTotal=0;
	static long notExtendedTotal=0;
	static long extensionsAttempted=0;
	/** Aggregated bridge attempts and classified outcomes. */
	static long bridgeAttemptsTotal=0, bridgeSuccessTotal=0;
	static long bridgeNoPathTotal=0, bridgeNoTargetTotal=0, bridgeAmbiguousTotal=0, bridgeCycleTotal=0;
	static long bridgeNonreciprocalTotal=0, bridgeRejectedTotal=0;
	/** Aggregated eligible mate checks and successful adapter-check calls. */
	static long adaptersExpected=0;
	static long adaptersFound=0;
	
	/** Shared accepted insert extrema; reset uses the unobserved minimum sentinel and zero maximum. */
	static int insertMinTotal=999999999;
	static int insertMaxTotal=0;
	
	/** Minimum number of overlapping bases required */
	private static int MIN_OVERLAPPING_BASES=11;
	/** Alternative minimum overlap threshold */
	private static int MIN_OVERLAPPING_BASES_0=8;
	private static int MISMATCH_MARGIN=2;
	private static int MIN_OVERLAPPING_BASES_RATIO_REDUCTION=3;
	
	/** Enables eligible mapping shortcuts; the worker path uses the stored Read.insert value. */
	protected static boolean USE_MAPPING=false;
	/** Suppress histogram rows with zero reported frequency. */
	private static boolean NONZERO_ONLY=true;
	/** Include the histogram's summary-statistics header. */
	private static boolean showHistStats=true;
	/** Retained training selection; requesting outinsert rejects this unimplemented exporter early. */
	static boolean TRAINING=false;
	
	/** Whether to use ratio mode for overlap detection */
	static boolean useRatioMode=true;
	/** Whether to use flat mode for overlap detection */
	static boolean useFlatMode=false;
	/** Require matching candidate inserts when both ratio and flat searches are enabled. */
	static boolean requireRatioMatch=false;
	/** Raw candidate mismatch ceiling that marks the selected result ambiguous. */
	static int MAX_MISMATCHES_R=20;
	/** Candidate mismatch-ratio ceiling passed to ratio search. */
	static float MAX_RATIO=0.09f;
	/** Candidate-separation multiplier passed to ratio search. */
	static float RATIO_MARGIN=5.5f;
	/** Mismatch-count offset passed to ratio search. */
	static float RATIO_OFFSET=0.55f;
	/** Second-candidate ratio floor used by ratio-search ambiguity checks. */
	static float MIN_SECOND_RATIO=0.1f;
	
	/** Main mismatch limit passed to flat overlap search. */
	public static int MAX_MISMATCHES=3;
	/** Initial mismatch limit passed to flat overlap search. */
	public static int MAX_MISMATCHES0=3;
	/** Initial flat-search numeric Phred threshold; retry rounds may lower it. */
	public static byte MIN_QUALITY=10;
	
	/** Return code indicating no overlap solution found */
	public static final int RET_NO_SOLUTION=-1;
	/** Return code indicating ambiguous overlap */
	public static final int RET_AMBIG=-2;
	/** Return code indicating bad/unusable read pair */
	public static final int RET_BAD=-3;
	/** Return code indicating insert too short */
	public static final int RET_SHORT=-4;
	/** Return code indicating insert too long */
	public static final int RET_LONG=-5;
	
	/** Shared requested output-file policies applied when descriptors/writers are built. */
	private static boolean overwrite=true;
	private static boolean append=false;
	/** Compile-time diagnostic verbosity switch. */
	private static final boolean verbose=false;
	/** Compile-time custom-tag collection path; currently disabled. */
	static final boolean TAG_CUSTOM=false;
	/** Shared feature-collection/neural-decision switch, configured with the overlapper's flag. */
	static boolean MAKE_VECTOR=false;
	/** Resolved default network resource path used during class initialization. */
	private static String defaultNetFile=Data.findPath("?bbmerge.bbnet", false);
	/** Shared requested/cached model path; pathname equality does not establish file freshness. */
	static String netFile=defaultNetFile;
	/** Shared model template; worker copies have separate working/row arrays but shared metadata. */
	static CellNet net0=CellNetParser.load(netFile, false);//False disables absent-file null handling
	/** Shared neural acceptance threshold, initially from the model or 0.5 when absent. */
	static float netCutoff=(net0==null ? 0.5f : net0.cutoff);
	/** Whether a cutoff was explicitly selected, preventing loadNet from replacing it. */
	private static boolean setCutoff=false;
	/** Retained collection/filter control used in candidate rejection branches. */
	private static boolean ONLY_OUTPUT_FALSE_POSITIVE=false;
	/** Retained false-negative-only flag with no active worker consumer. */
	private static boolean ONLY_OUTPUT_FALSE_NEGATIVE=false;

	//TODO: Remove after NN is trained
	/** Existing short-truth-overlap exemptions; branch-specific cutoffs remain unchanged. */
	private static boolean EXEMPT_SHORT_OVERLAPS_FROM_STATS=false;
	
	/** Convert undefined bases to uppercase N during preprocessing. */
	private static boolean iupacToN=false;
	
	/** Whether output error correction copies consensus qualities into present quality arrays. */
	public static boolean changeQuality=true;
	
	/** Output stream for messages and statistics */
	static PrintStream outstream=System.err;
	
	/** Shared worker count; construction resolves an unset value from Shared.threads. */
	private static int THREADS=-1;
//	private static String version="9.02";
	
	/** Prefix used by the dormant local PhiX diagnostic, separate from AdapterTracker's check. */
	private static final byte[] phixPrefix="AGATCGGAAGAGCG".getBytes();
	
	/** Built-in adapter strings; getAdapterList creates new list/byte-array copies. */
	private static final String[] defaultAdapters={
			"AGATCGGAAGAGCACA", "AATGATACGGCGACCA",	"GATCGGAAGAGCACAC", "CTGTCTCTTATACACA",
			"GACGCTGCCGACGA", "CCGAGCCCACGAGAC", "CTGATGGCGCGAGGGA", "CTGAGCGGGCTGGCAA",
			"GATCGGAAGAGCGTCG", "GATCGTCGGACTGTAG", "CCTTGGCACCCGAGAA", "CCACGGGAACGTGGTG",
			"TGGAATTCTCGGGTGC", "TCGGACTGTAGAACTC", "AGATCGGAAGAGCGGT"
	};
	
	/** Shared mutable adapters for static overlap verification; null skips that check.
	 * Callers must keep the list/arrays stable while static verification uses them. */
	public static ArrayList<byte[]> staticAdapterList=getAdapterList("default");
	
}
