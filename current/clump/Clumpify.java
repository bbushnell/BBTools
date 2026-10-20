package clump;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Random;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sort.SortByName;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import structures.Quantizer;

/**
 * Orchestrates pivot-based read grouping, optional duplicate removal/correction,
 * and final ordering. One group uses KmerSort1; multiple groups use an initial
 * KmerSplit and the selected original sorter version, with temporary files
 * between passes.
 * Configuration is shared with the worker classes and changes process-wide
 * settings, so independent invocations must not run concurrently in one JVM.
 * This class owns intermediate-file deletion and optional original-input deletion.
 * @author Brian Bushnell
 * @date Nov 6, 2015
 *
 */
public class Clumpify{

	/**
	 * Configures compression/quality policy, runs the orchestrator, and restores
	 * the specifically saved settings on normal return. The historical global-state
	 * limitations below also apply when callers invoke the constructor directly.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 6);

		//TODO: Possible bug [clump/Clumpify#001] (concern ACKNOWLEDGED by Brian 2026-06-24, left as-is for now:
		//"they could be leaked and should be restored, but for now they don't seem to be causing a problem") -
		//this save/restore of statics is (a) INCOMPLETE: ZIPLEVEL (restore commented out below), FASTQ.DETECT_QUALITY/ASCII_OFFSET
		//(set in runOnePass), FASTQ.FORCE_INTERLEAVED, Read.VALIDATE_IN_CONSTRUCTOR, and KmerComparator
		//defaultSeed/defaultBorder are all mutated but NOT restored; and (b) BYPASSED by the primary library
		//caller - jgi/RQCFilter2+3 call `new Clumpify(args)` directly (not main()), so NONE of this runs for
		//them; only assemble/TadPipe (via Clumpify.main) gets the 3-static restore. Whether the leaked FASTQ/
		//ZIPLEVEL statics affect downstream pipeline stages depends on each caller re-initializing them
		//(stages re-parse args, which usually resets FASTQ); the clump-only statics don't leak outside clump.
		//Likely benign in practice but unverified - flagged rather than "fixed" (adding restores could change
		//pipeline behavior that may rely on e.g. ZIPLEVEL=6, and the full mutated-static set spans the workers).
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
	 * Parses orchestration options and forwards worker-specific arguments. Expands
	 * paired filename markers, rejects input/output name collisions, checks primary
	 * input readability, and chooses groups.
	 * Does not run the sorting passes. Several options change shared worker state.
	 * @param args Command line arguments
	 */
	public Clumpify(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), true);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//args2 is a POSITIONAL template passed to the worker mains (KmerSplit/KmerSort*): indices are fixed
		//and rewritten by-position in process()/runOnePass[_v2] - [0]in1 [1]in2 [2]out1 [3]out2 [4]groups
		//[5]ecco [6]rename/addname [7]shortname [8]unpair [9]repair [10]namesort [11]overwrite. Keep this
		//order in lockstep with the args[N]=... writes below, or a worker receives the wrong flag.
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
			}else if(a.equals("unpair")){
				unpair=Parse.parseBoolean(b);
			}else if(a.equals("repair")){
				repair=Parse.parseBoolean(b);
			}else if(a.equals("namesort") || a.equals("sort")){
				namesort=Parse.parseBoolean(b);
			}else if(a.equals("overwrite")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("phasetimers")){
				KmerSort.phaseTimers=Parse.parseBoolean(b);
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
			}else if(a.equals("comparesequence")){
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
			}else if(Clump.parseStatic(arg, a, b)){
				//Do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//Do nothing
			}else{
				args2.add(arg);
			}
		}
		
		if(passes<1){
			throw new IllegalArgumentException("passes must be at least 1; otherwise grouped processing "
					+"would produce no output before optional input deletion: passes="+passes);
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
		//FIXED [clump/Clumpify#002]: grouped workers see temporary input names, so
		//validate originals against final outputs here before overwrite or deleteInput.
		final String[] fileNames={in1, in2, out1, out2};
		for(int i=0; i<fileNames.length; i++){
			//Absent-file options are sentinels, not repeated physical filenames.
			if(Tools.isNullFileName(fileNames[i])){fileNames[i]=null;}
		}
		Tools.testForDuplicateFiles(true, fileNames);
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
//		assert(false) : ReadKey.spanTiles()+", "+ReadKey.spanTilesX+", "+ReadKey.spanTilesY+", "+Clump.sortX+", "+Clump.sortY;
		
		autoSetGroups(gString);
		
		//Paired output does not imply interleaved input when reads will be unpaired
		//and repaired by name. Preserve the requested input mode for that path.
		if(!unpair && (in2!=null || out2!=null) && groups>1){FASTQ.FORCE_INTERLEAVED=true;} //Fix for crash with twin fasta files
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs one group once, or the requested grouped passes, delegating stream ownership
	 * to the selected workers. Optionally deletes original inputs after processing when
	 * sharedErrorState is clear and an output name was supplied. Stops the timer. */
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
			try{
				new File(in1).delete();
				if(in2!=null){new File(in2).delete();}
			}catch(Exception e){
				System.err.println("WARNING: Failed to delete input files.");
			}
		}
		
		t.stop();
		System.err.println("Total time: \t"+t);
		
	}
	
	/**
	 * Executes a single clumpification pass using the original algorithm.
	 * First splits input into groups, then sorts each group by k-mer similarity.
	 *
	 * @param args Processing arguments array
	 * @param pass Current pass number
	 * @param in1 Primary input file
	 * @param in2 Secondary input file (may be null)
	 * @param out1 Primary output file
	 * @param out2 Secondary output file (may be null)
	 */
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
			final boolean pairAfterSort=unpair && out2!=null;
			final String sorted=pairAfterSort ? getTempFname("clumpify_repair_") : out1;
			String[] sortArgs=new String[] {"in="+temp2, "out="+sorted, "ow="+overwrite};
			if(pairAfterSort){
				//SortByName re-enables interleaving detection in its constructor. These
				//records are still independent until the entire name sort has completed.
				sortArgs=new String[] {"in="+temp2, "out="+sorted, "ow="+overwrite, "interleaved=f"};
			}
			if(out2!=null && !pairAfterSort){
				//SortByName supplies ordered descriptors for paired Writer outputs.
				sortArgs=new String[] {"in="+temp2, "out="+out1, "out2="+out2, "ow="+overwrite};
			}
			SortByName.main(sortArgs);
			if(pairAfterSort){
				writeSortedPairs(sorted, out1, out2);
				if(delete){new File(sorted).delete();}
			}
			if(delete){new File(temp2).delete();}
		}
	}
	
	/**
	 * Executes a grouped pass using version 2 or 3.
	 * Uses memory estimates from previous passes to optimize processing.
	 * Handles external sorting for final pass if repair or namesort is enabled.
	 *
	 * @param args Processing arguments array
	 * @param pass Current pass number
	 * @param in1 Primary input file
	 * @param in2 Secondary input file (may be null)
	 * @param out1 Primary output file
	 * @param out2 Secondary output file (may be null)
	 * @param fileMem Prior estimated read storage in bytes; negative when unavailable
	 * @return Estimated read storage used to size later grouped passes, not measured heap usage
	 */
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
			//Intact pairs remain interleaved in each sorted temporary file. Preserve their
			//reader configuration so the final paired Writer receives actual mate links.
			//Unpaired processing can place mates in different files; those must remain separate here.
			if(unpair){FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;}
			//The merge emits dense ascending batch IDs; paired Writer outputs require ordering.
			final boolean pairAfterSort=unpair && out2!=null;
			final String sorted=pairAfterSort ? getTempFname("clumpify_repair_") : out1;
			FileFormat dest=FileFormat.testOutput(sorted, FileFormat.FASTQ, null, true, overwrite, false, true);
			FileFormat dest2=FileFormat.testOutput(pairAfterSort ? null : out2, FileFormat.FASTQ, null, true, overwrite, false, true);
			final boolean mergeError=SortByName.mergeAndDump(names, dest, dest2, delete,
					useSharedHeader, false, outstream, 1000);
			sharedErrorState|=mergeError;
			if(mergeError){throw new RuntimeException("Name-sorted output merge failed; retaining the original input files.");}
			if(pairAfterSort){
				writeSortedPairs(sorted, out1, out2);
				if(delete){new File(sorted).delete();}
			}
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
	
	/**
	 * Restores adjacent matching names after the global sort and writes complete pairs.
	 * Holds one pending read across input batches. Uses the same name convention as
	 * KmerSort.nameSort, including identical names. Unmatched names are rejected because
	 * two mate outputs have no singleton destination. Owns and closes both streams;
	 * failure propagates before process can delete the original input files.
	 */
	private void writeSortedPairs(String sorted, String out1, String out2){
		assert(out1!=null && out2!=null) : "Global mate restoration needs both mate-output filenames";
		final boolean oldTest=FASTQ.TEST_INTERLEAVED, oldForce=FASTQ.FORCE_INTERLEAVED;
		ConcurrentReadInputStream cris=null;
		ConcurrentReadOutputStream ros=null;
		ListNum<Read> ln=null;
		boolean error=false;
		try{
			FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
			final FileFormat input=FileFormat.testInput(sorted, FileFormat.FASTQ, null, true, true);
			final FileFormat one=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, false, true);
			final FileFormat two=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, false, true);
			ClumpTools.validateFileNames(new FileFormat[] {input}, new FileFormat[] {one, two});
			cris=ConcurrentReadInputStream.getReadInputStream(-1, false, input, null, null, null);
			cris.start();
			ros=ConcurrentReadOutputStream.getStream(one, two, 4, null, useSharedHeader);
			ros.start();
			Read pending=null;
			long id=0;
			ln=cris.nextList();
			while(ln!=null && ln.list!=null && !ln.list.isEmpty()){
				final ArrayList<Read> pairs=new ArrayList<Read>((ln.list.size()+1)/2);
				for(Read r : ln.list){
					assert(r.mate==null) : "The sorted repair input must be read as individual records: "+r.id;
					if(pending==null){pending=r;}
					else if(pending.id!=null && r.id!=null && FASTQ.testPairNames(pending.id, r.id, true)){
						pending.mate=r; r.mate=pending;
						pending.setPairnum(0); r.setPairnum(1);
						pairs.add(pending);
						pending=null;
					}else{throw unmatchedRead(pending);}
				}
				if(!pairs.isEmpty()){ros.add(pairs, id++);}
				cris.returnList(ln);
				ln=null;
				ln=cris.nextList();
			}
			if(pending!=null){throw unmatchedRead(pending);}
		}catch(RuntimeException | Error e){
			sharedErrorState=true;
			throw e;
		}finally{
			try{
				if(ln!=null){cris.returnList(ln);}
			}finally{
				try{error|=ReadWrite.closeStream(cris);}
				finally{
					try{error|=ReadWrite.closeStream(ros);}
					finally{
						FASTQ.TEST_INTERLEAVED=oldTest; FASTQ.FORCE_INTERLEAVED=oldForce;
						sharedErrorState|=error;
					}
				}
			}
		}
		if(error){throw new RuntimeException("Paired output failed; retaining the original input files.");}
	}

	/** Explains why a singleton cannot be emitted into two mate files. */
	private static IllegalArgumentException unmatchedRead(Read read){
		return new IllegalArgumentException("Unmatched read '"+read.id+"' while writing two mate-output files; "
				+"name-sorted reads must have matching adjacent mates. Original inputs retained.");
	}

	/**
	 * Accepts an explicit group count or estimates one from sampled input and heap
	 * availability. Estimates guide partitioning but do not guarantee that skewed
	 * groups fit in memory. A standard-input descriptor retains the existing count.
	 * @param s Groups setting string ("auto" or numeric value)
	 */
	private void autoSetGroups(String s){
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
		
		final double overhead=0.5*(ReadKey.overhead+Clump.overhead);
		double[] estimates=Tools.estimateFileMemory(in1, 1000, overhead, true, lowComplexity);
		if(in2!=null){
			double[] estimates2=Tools.estimateFileMemory(in2, 1000, overhead, true, lowComplexity);
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
			final int fetchers=(V3 ? KmerSort3.fetchThreads : 2);
			groups=Tools.max(11, (int)(3+(3*worstCase/maxMem)*fetchers),
					(int)((2*readEstimate)/Integer.MAX_VALUE))|1;
		}
		outstream.println("Set groups to "+groups);
	}
	
	/**
	 * Constructs a temporary filename with a random suffix; does not create or reserve it.
	 * Uses output path and extension as template, optionally using system tmpdir.
	 * Applies requested temporary-compression policy and strips a final bz2 suffix.
	 *
	 * @param core Base name for the temporary file
	 * @return Complete temporary file path with random suffix
	 */
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
	
	/**
	 * Shortens read IDs by replacing common sequencer prefixes with abbreviations.
	 * Converts "HISEQ" to "H", "MISEQ" to "M", and removes the literal
	 * ":000000000-" prefix padding. A null read is ignored.
	 * Reduces file size by compacting verbose sequencer-generated names.
	 * @param r Read whose ID should be shortened
	 */
	public static void shrinkName(Read r){
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
	
	/**
	 * Replaces read ID with a compact numeric representation.
	 * Uses a nonnegative numeric ID padded to at least ten digits plus pair number.
	 * Creates IDs in format "0000000123 1:" or "0000000123 2:".
	 * @param r Read whose ID should be replaced with short form; null is ignored
	 */
	public static void shortName(Read r){
		if(r==null){return;}//Workers also pass the absent mate of an unpaired read.
		ByteBuilder sb=new ByteBuilder(14);
		//CLEVER [verified]: |1 guards against numericID==0 making x stay 0 and looping forever; it only flips
		//the low bit so it never changes the digit-count, so the zero-padding (to 10 digits) is unaffected.
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
	
	/** Selects the low-complexity input-memory estimate used for automatic grouping. */
	private boolean lowComplexity=false;
	
	/** Whether to quantize quality scores to reduce file size */
	private boolean quantizeQuality=false;
	/** Random suffix generator; filenames are not atomically reserved. */
	private Random randy=shared.Shared.random();
	/** Number of groups to split sequences into for processing */
	private int groups=31;
	/** Number of passes in the grouped path; the one-group path invokes KmerSort1 once. */
	private int passes=1;
	/** Whether to request strict overlap-based correction in the first split. */
	private boolean ecco=false;
	/** Whether to add names during processing */
	private boolean addName=false;
	/** Short name setting for read ID modification */
	private String shortName="f";
	/** Whether to use system temporary directory for intermediate files */
	private boolean useTmpdir=false;
	/** Whether to delete temporary files after processing */
	private boolean delete=true;
	/** Whether to delete input files after successful processing */
	private boolean deleteInput=false;
	/** Whether to use shared header for SAM/BAM file processing */
	private boolean useSharedHeader=false;
	/** Forces compression of temporary files regardless of output format */
	private boolean forceCompressTemp=false;
	/** Forces temporary files to be uncompressed regardless of output format */
	private boolean forceRawTemp=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;

	/** Whether to unpair paired-end reads into single-end format */
	private boolean unpair=false;
	/** Whether to repair broken paired-end read files */
	private boolean repair=false;
	/** Whether to sort output by read names */
	private boolean namesort=false;
	/** Whether to use version 2 of the clumpification algorithm */
	private boolean V2=false;
	/** Whether to use version 3 of the clumpification algorithm */
	private boolean V3=true;

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path for paired-end data */
	private String in2=null;
	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path for paired-end data */
	private String out2=null;
	
	/** Positional argument template for worker invocations; see constructor's index contract. */
	ArrayList<String> args2=new ArrayList<String>();
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;

	/** Error flag reported by workers and consulted before original-input deletion.
	 * This class does not reset the flag for subsequent invocations. */
	public static boolean sharedErrorState=false;
	
}
