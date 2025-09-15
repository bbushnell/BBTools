package clump;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sort.SortByName;
import stream.FASTQ;
import stream.Read;
import structures.ByteBuilder;
import structures.Quantizer;

/**
 * @author Brian Bushnell
 * @date Nov 6, 2015
 *
 */
public class Clumpify {

	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 6);
		
		//Capture values of static variables that might be modified in case this is called by another class.
		final boolean oldCQ=Read.CHANGE_QUALITY;
		final boolean oldBgzip=ReadWrite.USE_BGZIP, oldPreferBgzip=ReadWrite.PREFER_BGZIP;
		
		BBMerge.changeQuality=Read.CHANGE_QUALITY=false;
		ReadWrite.USE_BGZIP=true;
		ReadWrite.PREFER_BGZIP=true;
		
		Clumpify x=new Clumpify(args);
		x.process(t);
		
		//Restore values of static variables.
//		Shared.setBuffers(oldCap);
//		ReadWrite.ZIPLEVEL=oldZl;
//		ReadWrite.USE_PIGZ=oldPigz;
		ReadWrite.USE_BGZIP=oldBgzip;
		ReadWrite.PREFER_BGZIP=oldPreferBgzip;
//		ReadWrite.USE_UNPIGZ=oldUnpigz;
//		ReadWrite.MAX_ZIP_THREADS=oldZipThreads;
		BBMerge.changeQuality=Read.CHANGE_QUALITY=oldCQ;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Clumpify(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), true);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		args2=new ArrayList<String>();
		args2.add("in1");
		args2.add("in2");
		args2.add("out1");
		args2.add("out2");
		args2.add("groups");
		args2.add("ecco=f");
		args2.add("rename=f");
		args2.add("shortname=f");
		args2.add("unpair=f");
		args2.add("repair=f");
		args2.add("namesort=f");
		args2.add("overwrite=t");
		
		String gString="auto";
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("out") || a.equals("out1")){
				out1=b;
			}else if(a.equals("out2")){
				out2=b;
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets") || a.equals("ways")){
				gString=b;
			}else if(a.equals("delete") || a.equals("deletetemp")){
				delete=Parse.parseBoolean(b);
			}else if(a.equals("deleteinput")){
				deleteInput=Parse.parseBoolean(b);
			}else if(a.equals("usetmpdir")){
				useTmpdir=Parse.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco=Parse.parseBoolean(b);
			}else if(a.equals("compresstemp") || a.equals("ct")){
				if(b!=null && b.equalsIgnoreCase("auto")){forceCompressTemp=forceRawTemp=false;}
				else{
					forceCompressTemp=Parse.parseBoolean(b);
					forceRawTemp=!forceCompressTemp;
				}
			}else if(a.equals("tmpdir")){
				Shared.setTmpdir(b);
			}else if(a.equals("rename") || a.equals("addname")){
				addName=Parse.parseBoolean(b);
			}else if(a.equals("shortname") || a.equals("shortnames")){
				shortName=b;
			}else if(a.equals("seed")){
				KmerComparator.defaultSeed=Long.parseLong(b);
			}else if(a.equals("hashes")){
				KmerComparator.setHashes(Integer.parseInt(b));
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
				args2.add(arg);
//			}else if(a.equals("k")){
//				k=Integer.parseInt(b);
//				args2.add(arg);
			}else if(a.equals("border")){
				KmerComparator.defaultBorder=Integer.parseInt(b);
			}

			else if(a.equals("unpair")){
				unpair=Parse.parseBoolean(b);
			}else if(a.equals("repair")){
				repair=Parse.parseBoolean(b);
			}else if(a.equals("namesort") || a.equals("sort")){
				namesort=Parse.parseBoolean(b);
			}else if(a.equals("overwrite")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("v1") || a.equals("kmersort1")){
				boolean x=Parse.parseBoolean(b);
				if(x){V2=V3=false;}
			}else if(a.equals("v2") || a.equals("kmersort2")){
				V2=Parse.parseBoolean(b);
				if(V2){V3=false;}
			}else if(a.equals("v3") || a.equals("kmersort3")){
				V3=Parse.parseBoolean(b);
				if(V3){V2=false;}
			}else if(a.equals("fetchthreads")){
				KmerSort3.fetchThreads=Integer.parseInt(b);
				assert(KmerSort3.fetchThreads>0) : KmerSort3.fetchThreads+"\nFetch threads must be at least 1.";
			}
			
			else if(a.equals("comparesequence")){
				KmerComparator.compareSequence=Parse.parseBoolean(b);
			}else if(a.equals("allowadjacenttiles") || a.equals("spantiles")){
				ReadKey.spanTilesX=ReadKey.spanTilesY=Parse.parseBoolean(b);
			}else if(a.equals("spanx") || a.equals("spantilesx")){
				ReadKey.spanTilesX=Parse.parseBoolean(b);
			}else if(a.equals("spany") || a.equals("spantilesy")){
				ReadKey.spanTilesY=Parse.parseBoolean(b);
			}else if(a.equals("spanadjacent") || a.equals("spanadjacentonly") || a.equals("adjacentonly") || a.equals("adjacent")){
				ReadKey.spanAdjacentOnly=Parse.parseBoolean(b);
			}
			
//			else if(a.equals("repair")){
//				repair=Parse.parseBoolean(b);
//			}else if(a.equals("namesort") || a.equals("sort")){
//				namesort=Parse.parseBoolean(b);
//			}
			
			else if(a.equals("interleaved") || a.equals("int")){
				if("auto".equalsIgnoreCase(b)){FASTQ.FORCE_INTERLEAVED=!(FASTQ.TEST_INTERLEAVED=true);}
				else{
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=Parse.parseBoolean(b);
					System.err.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}else if(a.equals("cq") || a.equals("changequality")){
				BBMerge.changeQuality=Read.CHANGE_QUALITY=Parse.parseBoolean(b);
			}else if(a.equals("quantize") || a.equals("quantizesticky")){
				quantizeQuality=Quantizer.parse(arg, a, b);
			}else if(a.equals("lowcomplexity")){
				lowComplexity=Parse.parseBoolean(b);
			}
			
			else if(Clump.parseStatic(arg, a, b)){
				//Do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//Do nothing
			}
			
			else{
				args2.add(arg);
			}
		}
		
		Clump.setXY();
		
		KmerSplit.quantizeQuality=KmerSort1.quantizeQuality=quantizeQuality;
		
		Parser.processQuality();
		
		assert(!unpair || !KmerComparator.mergeFirst) : "Unpair and mergefirst may not be used together.";
		
		if(in1==null){throw new RuntimeException("\nOne input file is required.\n");}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
//		assert(false) : ReadKey.spanTiles()+", "+ReadKey.spanTilesX+", "+ReadKey.spanTilesY+", "+Clump.sortX+", "+Clump.sortY;
		
		autoSetGroups(gString);
		
		if((in2!=null || out2!=null) && groups>1){FASTQ.FORCE_INTERLEAVED=true;} //Fix for crash with twin fasta files
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	public void process(Timer t){
		String[] args=args2.toArray(new String[0]);
		args[4]="groups="+groups;
		
		useSharedHeader=(FileFormat.hasSamOrBamExtension(in1) && out1!=null
				&& FileFormat.hasSamOrBamExtension(out1));
		
		if(groups==1){
			args[0]="in1="+in1;
			args[1]="in2="+in2;
			args[2]="out1="+out1;
			args[3]="out2="+out2;
			args[5]="ecco="+ecco;
			args[6]="rename="+addName;
			args[7]="shortname="+shortName;
			args[8]="unpair="+unpair;
			args[9]="repair="+repair;
			args[10]="namesort="+namesort;
			args[11]="ow="+overwrite;
			KmerSort1.main(args);
		}else{
			String pin1=in1, pin2=in2, temp;
			final int conservativePasses=Clump.conservativeFlag ? passes : Tools.max(1, passes/2);
			if(passes>1){Clump.setConservative(true);}
			long fileMem=-1;
			for(int pass=1; pass<=passes; pass++){
				if(/*passes>1 &&*/ (V2 || V3)){
//					System.err.println("Running pass with fileMem="+fileMem);
//					out=(pass==passes ? out1 : getTempFname("clumpify_p"+(pass+1)+"_temp%_"));
					temp=getTempFname("clumpify_p"+(pass+1)+"_temp%_");
					if(pass==passes){
						fileMem=runOnePass_v2(args, pass, pin1, pin2, out1, out2, fileMem);
					}else{
						fileMem=runOnePass_v2(args, pass, pin1, pin2, temp, null, fileMem);
					}
//					System.err.println("New fileMem="+fileMem);
				}else{
//					out=(pass==passes ? out1 : getTempFname("clumpify_temp_pass"+pass+"_"));
					temp=getTempFname("clumpify_temp_pass"+pass+"_");
					if(pass==passes){
						runOnePass(args, pass, pin1, pin2, out1, out2);
					}else{
						runOnePass(args, pass, pin1, pin2, temp, null);
					}
				}
				pin1=temp;
				pin2=null;
				KmerComparator.defaultBorder=Tools.max(0, KmerComparator.defaultBorder-1);
				KmerComparator.defaultSeed++;
				if(pass>=conservativePasses){Clump.setConservative(false);}
			}
		}
		
		if(deleteInput && !sharedErrorState && out1!=null && in1!=null){
			try {
				new File(in1).delete();
				if(in2!=null){new File(in2).delete();}
			} catch (Exception e) {
				System.err.println("WARNING: Failed to delete input files.");
			}
		}
		
		t.stop();
		System.err.println("Total time: \t"+t);
		
	}
	
	private void runOnePass(String[] args, int pass, String in1, String in2, String out1, String out2){
		assert(groups>1);
		if(pass>1){
			ecco=false;
			shortName="f";
			addName=false;
		}

		String temp=getTempFname("clumpify_p"+pass+"_temp%_");
		
		String temp2=temp.replace("%", "FINAL");
		final boolean externalSort=(pass==passes && (repair || namesort));

		args[0]="in1="+in1;
		args[1]="in2="+in2;
		args[2]="out="+temp;
		args[3]="out2="+null;
		args[5]="ecco="+ecco;
		args[6]="addname=f";
		args[7]="shortname="+shortName;
		args[8]="unpair="+unpair;
		args[9]="repair=f";
		args[10]="namesort=f";
		args[11]="ow="+overwrite;
		KmerSplit.maxZipLevel=2;
		KmerSplit.main(args);
		
		FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=false;
		FASTQ.ASCII_OFFSET=FASTQ.ASCII_OFFSET_OUT;

		args[0]="in="+temp;
		args[1]="in2="+null;
		args[2]="out="+(externalSort ? temp2 : out1);
		args[3]="out2="+(externalSort ? "null" : out2);
		args[5]="ecco=f";
		args[6]="addname="+addName;
		args[7]="shortname=f";
		args[8]="unpair=f";
		args[9]="repair="+(repair && externalSort);
		args[10]="namesort="+(namesort && externalSort);
		args[11]="ow="+overwrite;
		if(unpair){
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		KmerSort1.main(args);
		
		if(delete){
			for(int i=0; i<groups; i++){
				new File(temp.replaceFirst("%", ""+i)).delete();
			}
			if(pass>1){
				assert(in2==null);
				new File(in1).delete();
			}
		}
		
		if(externalSort){
			outstream.println();
			String[] sortArgs=new String[] {"in="+temp2, "out="+out1, "ow="+overwrite};
			if(out2!=null){sortArgs=new String[] {"in="+temp2, "out="+out1, "out2="+out2, "ow="+overwrite};}
			SortByName.main(sortArgs);
			if(delete){new File(temp2).delete();}
		}
	}
	
	private long runOnePass_v2(String[] args, int pass, String in1, String in2, String out1, String out2, long fileMem){
		assert(groups>1);
		if(pass>1){
			ecco=false;
			shortName="f";
			addName=false;
		}
		
		String temp=getTempFname("clumpify_p"+pass+"_temp%_");
		
//		String temp2=temp.replace("%", "FINAL");
		String namesorted=temp.replace("%", "namesorted_%");
		final boolean externalSort=(pass==passes && (repair || namesort));
		
		if(pass==1){
			args[0]="in1="+in1;
			args[1]="in2="+in2;
			args[2]="out="+temp;
			args[3]="out2="+null;
			args[5]="ecco="+ecco;
			args[6]="addname=f";
			args[7]="shortname="+shortName;
			args[8]="unpair="+unpair;
			args[9]="repair=f";
			args[10]="namesort=f";
			args[11]="ow="+overwrite;
			KmerSplit.maxZipLevel=2;
			KmerSplit.main(args);
			fileMem=KmerSplit.lastMemProcessed;
			
			FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=false;
			FASTQ.ASCII_OFFSET=FASTQ.ASCII_OFFSET_OUT;
		}
		
		args[0]="in1="+(pass==1 ? temp : in1);
		args[1]="in2="+null;
		args[2]="out="+(externalSort ? namesorted : out1);
		args[3]="out2="+(externalSort ? "null" : out2);
		args[5]="ecco=f";
		args[6]="addname="+addName;
		args[7]="shortname=f";
		args[8]="unpair=f";
		args[9]="repair="+(repair && externalSort);
		args[10]="namesort="+(namesort && externalSort);
		args[11]="ow="+overwrite;
		if(unpair){
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		if(externalSort){
			KmerSort.doHashAndSplit=false;
		}
		if(V3){
			KmerSort3.main(fileMem, pass, passes, args);
			if(fileMem<1){fileMem=KmerSort3.lastMemProcessed;}
		}else{KmerSort2.main(args);}
		
		if(delete){
			for(int i=0; i<groups; i++){
				new File((pass==1 ? temp : in1).replaceFirst("%", ""+i)).delete();
			}
		}
		
		if(externalSort){
			outstream.println();
			
			ArrayList<String> names=new ArrayList<String>();
			for(int i=0; i<groups; i++){
				names.add(namesorted.replaceFirst("%", ""+i));
			}
			ReadWrite.setZipThreads(Shared.threads());
			
			ReadWrite.USE_PIGZ=true;
			ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 6);
			FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
			FileFormat dest=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, false, false);
			FileFormat dest2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, false, false);
			SortByName.mergeAndDump(names, /*null, */dest, dest2, delete, useSharedHeader, false, outstream, 1000);
		}
		
//		if(externalSort){
//			outstream.println();
//			SortByName.main(new String[] {"in="+temp2, "out="+out, "ow="+overwrite});
//			if(delete){new File(temp2).delete();}
//		}
		return fileMem;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void autoSetGroups(String s) {
		if(s==null || s.equalsIgnoreCase("null")){return;}
		if(Tools.isDigit(s.charAt(0))){
			groups=Integer.parseInt(s);
			return;
		}
		assert(s.equalsIgnoreCase("auto")) : "Unknown groups setting: "+s;
		
		final long maxMem=Shared.memAvailable(1);
		FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, false, false);
		if(ff1==null || ff1.stdio()){return;}
		
//		outstream.println("in1="+in1+", overhead="+(0.5*(ReadKey.overhead+Clump.overhead)));
		
		double[] estimates=Tools.estimateFileMemory(in1, 1000, 0.5*(ReadKey.overhead+Clump.overhead), true, lowComplexity);
		if(in2!=null){
			double[] estimates2=Tools.estimateFileMemory(in2, 1000, 0.5*(ReadKey.overhead+Clump.overhead), true, lowComplexity);
			estimates[0]+=estimates2[0];
			estimates[1]+=estimates2[1];
			estimates[4]+=estimates2[4];
		}
		
//		outstream.println(Arrays.toString(estimates));
		
		double memEstimate=estimates==null ? 0 : estimates[0];
		double diskEstimate=estimates==null ? 0 : estimates[1];
		double readEstimate=estimates==null ? 0 : estimates[4];
		double worstCase=memEstimate*1.5;

//		outstream.println("Raw Disk Size Estimate: "+(long)(diskEstimate/(1024*1024))+" MB");
		outstream.println("Read Estimate:          "+(long)(readEstimate));
		outstream.println("Memory Estimate:        "+(long)(memEstimate/(1024*1024))+" MB");
		outstream.println("Memory Available:       "+(maxMem/(1024*1024))+" MB");
		
		if(maxMem>worstCase && readEstimate<Integer.MAX_VALUE){
			groups=1;
		}else{
			groups=Tools.max(11, (int)(3+(3*worstCase/maxMem)*(V3 ? KmerSort3.fetchThreads : 2)), (int)((2*readEstimate)/Integer.MAX_VALUE))|1;
		}
		outstream.println("Set groups to "+groups);
	}
	
	private String getTempFname(String core){
//		outstream.println(core);
		String temp;
		String path="", extension=".fq";
		if(out1!=null){
			core=ReadWrite.stripToCore(out1)+"_"+core;
			path=ReadWrite.getPath(out1);
			extension=ReadWrite.getExtension(out1);
		}
		
		if(useTmpdir && Shared.tmpdir()!=null){
			temp=Shared.tmpdir()+core+Long.toHexString((randy.nextLong()&Long.MAX_VALUE))+extension;
		}else{
			temp=path+core+Long.toHexString((randy.nextLong()&Long.MAX_VALUE))+extension;
		}
//		assert(false) : path+", "+temp+", "+core+", "+out1;
		
		String comp=ReadWrite.compressionType(temp);
		if(comp!=null){comp=".gz";} //Prevent bz2 temp files which cause a crash
		
		if(forceCompressTemp && comp==null){
			temp+=".gz";
		}else if(comp!=null && forceRawTemp){
			temp=temp.substring(0, temp.lastIndexOf('.'));
		}
		if(temp.endsWith(".bz2")){temp=temp.substring(0, temp.length()-4);} //Prevent bz2 temp files which cause a crash

//		outstream.println(temp);
		return temp;
	}
	
	public static void shrinkName(Read r) {
		if(r==null){return;}
		String s=r.id;
		if(s.contains("HISEQ")){s=s.replace("HISEQ", "H");}
		if(s.contains("MISEQ")){
			s=s.replace("MISEQ", "M");
		}
		if(s.contains(":000000000-")){
			s=s.replace(":000000000-", ":");
		}
		r.id=s;
	}
	
	public static void shortName(Read r) {
		ByteBuilder sb=new ByteBuilder(14);
		long x=r.numericID|1;
		
		while(x<1000000000L){
			x*=10;
			sb.append('0');
		}
		sb.append(r.numericID);
		
//		while(x<0x10000000L){
//			x*=16;
//			sb.append('0');
//		}
//		sb.append(Long.toHexString(r.numericID));
		
		sb.append(r.pairnum()==0 ? " 1:" : " 2:");
		r.id=sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private boolean lowComplexity=false;
	
	private boolean quantizeQuality=false;
	private Random randy=new Random();
	private int groups=31;
	private int passes=1;
	private boolean ecco=false;
	private boolean addName=false;
	private String shortName="f";
	private boolean useTmpdir=false;
	private boolean delete=true;
	private boolean deleteInput=false;
	private boolean useSharedHeader=false;
	private boolean forceCompressTemp=false;
	private boolean forceRawTemp=false;
	private boolean overwrite=true;

	private boolean unpair=false;
	private boolean repair=false;
	private boolean namesort=false;
	private boolean V2=false;
	private boolean V3=true;

	private String in1=null;
	private String in2=null;
	private String out1=null;
	private String out2=null;
	
	ArrayList<String> args2=new ArrayList<String>();
	private PrintStream outstream=System.err;

	public static boolean sharedErrorState=false;
	
}
