package clump;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public class KmerSort1 extends KmerSort {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		final boolean pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		final float ztd=ReadWrite.ZIP_THREAD_MULT();
		final int mzt=ReadWrite.MAX_ZIP_THREADS();
		final int oldzl=ReadWrite.ZIPLEVEL;
		Timer t=new Timer();
		KmerSort1 x=new KmerSort1(args);
		x.process(t);
		ReadWrite.USE_PIGZ=pigz;
		ReadWrite.USE_UNPIGZ=unpigz;
		ReadWrite.setZipThreadMult(ztd);
		ReadWrite.setZipThreads(mzt);
		ReadWrite.ZIPLEVEL=oldzl;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerSort1(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=KmerComparator.verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
				assert(k>0 && k<32);
			}else if(a.equals("mincount") || a.equals("mincr")){
				minCount=Integer.parseInt(b);
			}else if(a.equals("rename") || a.equals("addname")){
				addName=Parse.parseBoolean(b);
			}else if(a.equals("shortname") || a.equals("shortnames")){
				if(b!=null && b.equals("shrink")){
					shrinkName=true;
				}else{
					shrinkName=false;
					shortName=Parse.parseBoolean(b);
				}
			}else if(a.equals("rcomp") || a.equals("reversecomplement")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco=Parse.parseBoolean(b);
			}else if(a.equals("condense") || a.equals("consensus") || a.equals("concensus")){//Note the last one is intentionally misspelled
				condense=Parse.parseBoolean(b);
			}else if(a.equals("correct") || a.equals("ecc")){
				correct=Parse.parseBoolean(b);
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
			}
			
			else if(a.equals("dedupe")){
				dedupe=Parse.parseBoolean(b);
			}else if(a.equals("markduplicates")){
				dedupe=Clump.markOnly=Parse.parseBoolean(b);
			}else if(a.equals("markall")){
				boolean x=Parse.parseBoolean(b);
				if(x){
					dedupe=Clump.markOnly=Clump.markAll=true;
				}else{
					Clump.markAll=false;
				}
			}
			
			else if(a.equals("prefilter")){
				KmerReduce.prefilter=Parse.parseBoolean(b);
			}else if(a.equals("entryfilter")){
				entryfilter=Parse.parseBoolean(b);
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets") || a.equals("ways")){
				groups=Integer.parseInt(b);
				splitInput=(groups>1);
			}else if(a.equals("seed")){
				KmerComparator.defaultSeed=Long.parseLong(b);
			}else if(a.equals("hashes")){
				KmerComparator.setHashes(Integer.parseInt(b));
			}else if(a.equals("border")){
				KmerComparator.defaultBorder=Integer.parseInt(b);
			}else if(a.equals("minprob")){
				KmerComparator.minProb=Float.parseFloat(b);
				
			}else if(a.equals("unpair")){
				unpair=Parse.parseBoolean(b);
			}else if(a.equals("repair")){
				repair=Parse.parseBoolean(b);
			}else if(a.equals("namesort") || a.equals("sort")){
				namesort=Parse.parseBoolean(b);
			}else if(a.equals("reorder") || a.equals("reorderclumps") || a.equals("reordermode")){
				reorderMode=REORDER_AUTO;
				if(b==null || b.equalsIgnoreCase("auto") || b.equalsIgnoreCase("a")){
					reorderMode=REORDER_AUTO;
				}else if(b.equalsIgnoreCase("unpaired") || b.equalsIgnoreCase("consensus") || b.equalsIgnoreCase("reorder") || b.equalsIgnoreCase("c")){
					reorderMode=REORDER_CONSENSUS;
				}else if(b.equalsIgnoreCase("pair") || b.equalsIgnoreCase("pairs") || b.equalsIgnoreCase("paired") || b.equalsIgnoreCase("p")){
					reorderMode=REORDER_PAIRED;
				}else{
					boolean x=Parse.parseBoolean(b);
					if(x){
						reorderMode=REORDER_AUTO;
					}else{
						reorderMode=REORDER_FALSE;
					}
				}
			}else if(a.equals("reorderpaired") || a.equals("reorderclumpspaired")){
				boolean x=Parse.parseBoolean(b);
				if(x){
					reorderMode=REORDER_PAIRED;
				}else{
					reorderMode=REORDER_FALSE;
				}
			}
			
			else if(a.equals("fetchthreads")){
				//Do nothing
			}else if(Clump.parseStatic(arg, a, b)){
				//Do nothing
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		Clump.renameConsensus=condense;
		if(dedupe){KmerComparator.compareSequence=true;}
		assert(!(reorderMode==REORDER_PAIRED && dedupe)) : "REORDER_PAIRED and dedupe are incompatible.";
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			out2=parser.out2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}

		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		
		if(groups>1 && in1.contains("%") && (splitInput || !new File(in1).exists())){
			ffin1=new FileFormat[groups];
			ffin2=new FileFormat[groups];
			for(int i=0; i<groups; i++){
				ffin1[i]=FileFormat.testInput(in1.replaceFirst("%", ""+i), FileFormat.FASTQ, extin, true, true);
				ffin2[i]=in2==null ? null : FileFormat.testInput(in2.replaceFirst("%", ""+i), FileFormat.FASTQ, extin, true, true);
			}
		}else{
			assert(!in1.contains("%") && groups==1) : "The % symbol must only be present in the input filename if groups>1.";
			ffin1=new FileFormat[1];
			ffin1[0]=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
			ffin2=new FileFormat[1];
			ffin2[0]=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
			groups=1;
		}
//		if(groups>1){ReadWrite.USE_UNPIGZ=false;} //Not needed since they are not concurrent
		
		if((reorderMode!=REORDER_FALSE) && (passes>1 || condense || correct || groups>1)){
			outstream.println("Clump reordering disabled because "+(passes>1 ? "passes>1" : condense ? " condense=t" : correct ? " ecc=t" : "groups>1"));
			reorderMode=REORDER_FALSE;
		}
		
		if(reorderMode==REORDER_PAIRED){
			if(!unpair || !repair){
				outstream.println("Unpair and repair enabled because clump reorder mode is set to paired.");
				unpair=true;
				repair=true;
			}
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	void process(Timer t){
		
		preprocess();

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=1; //This prevents more than 2 sets of reads from being in memory at once.
			assert(!out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, useSharedHeader);
			ros.start();
		}else{ros=null;}
		
		readsProcessed=basesProcessed=diskProcessed=memProcessed=0;
		
		//Process the read stream
		processInner(ros);
		
		printStats(t);
	}
	
	/** Collect and sort the reads */
	void processInner(final ConcurrentReadOutputStream ros){
		if(verbose){outstream.println("Making comparator.");}
		KmerComparator kc=new KmerComparator(k, addName, (rcomp || condense || correct));
		
		ClumpList.UNRCOMP=(!rcomp && !condense);
		Timer t=new Timer();
		
		final int conservativePasses=Clump.conservativeFlag ? passes : Tools.max(1, passes/2);
		if(groups==1 && passes>1){Clump.setConservative(true);}

		useSharedHeader=(ffin1[0].samOrBam() && ffout1!=null && ffout1.samOrBam());
		
		for(int group=0; group<groups; group++){
			if(verbose){outstream.println("Starting cris "+group+".");}
			
			final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads,
					useSharedHeader && groups==1, ffin1[group], ffin2[group], null, null);
			cris.start();
			if(reorderMode!=REORDER_FALSE){
				assert(groups==1) : "Too many groups for reorder: "+groups;
				if(reorderMode==REORDER_AUTO){
					if(cris.paired() && !dedupe){
						reorderMode=REORDER_PAIRED;
						unpair=repair=true;
					}else{
						reorderMode=REORDER_CONSENSUS;
					}
				}
			}
			
			if(verbose){t.start("Fetching reads.");}
			ArrayList<Read> reads=fetchReads1(cris, kc);
			quantizeQuality=false;
//			if(verbose){t.stop("Fetch time: ");}
			
			if(verbose){t.start("Sorting.");}
			Shared.sort(reads, kc);
			if(verbose){t.stop("Sort time: ");}
			
//			if(verbose){t.start("Counting clumps.");}
//			clumpsProcessed+=countClumps(reads);
//			if(verbose){t.stop("Count time: ");}
			
			if(verbose){t.start("Making clumps.");}
			readsProcessedThisPass=reads.size();
			
			ClumpList cl=new ClumpList(reads, k, reorderMode==REORDER_CONSENSUS);
			
			if(reorderMode!=REORDER_FALSE){
				reads.clear();
				if(reorderMode==REORDER_PAIRED){cl.reorderPaired();}
				else if(reorderMode==REORDER_CONSENSUS){cl.reorder();}
				else{assert(false) : reorderMode;}
				for(Clump c : cl){
					reads.addAll(c);
				}
			}
			
			clumpsProcessedThisPass=cl.size();
			clumpsProcessedTotal+=clumpsProcessedThisPass;
			if(verbose){t.stop("Clump time: ");}
			
			if(dedupe){
				reads.clear();
				if(verbose){t.start("Deduping.");}
				reads=processClumps(cl, ClumpList.DEDUPE);
				
				if(passes>1 && groups==1){
					
					FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=false;
					FASTQ.ASCII_OFFSET=FASTQ.ASCII_OFFSET_OUT;
					
					if(verbose){outstream.println("Pass 1.\n");}
					if(verbose){outstream.println("Reads:        \t"+readsProcessedThisPass);}
					outstream.println("Clumps:       \t"+clumpsProcessedThisPass);
					
					for(int pass=1; pass<passes; pass++){
						
						kc=new KmerComparator(k, kc.seed<0 ? -1 : kc.seed+1, kc.border-1, kc.hashes, false, kc.rcompReads);
						reads=runOnePass(reads, kc);

						if(verbose){outstream.println("Seed: "+kc.seed);}
						if(verbose){outstream.println("Pass "+(pass+1)+".");}
						outstream.println();
					}
				}
				
				if(verbose){t.stop("Dedupe time: ");}
			}else if(condense){
				reads.clear();
				if(verbose){t.start("Condensing.");}
				reads=processClumps(cl, ClumpList.CONDENSE);
				if(verbose){t.stop("Condense time: ");}
			}else if(correct){
				reads.clear();
				if(verbose){t.start("Correcting.");}
				reads=processClumps(cl, ClumpList.CORRECT);
				if(verbose){t.stop("Correct time: ");}
				
				if(verbose){outstream.println("Seed: "+kc.seed);}
				if(groups>1){
					if(verbose){outstream.println("Reads:        \t"+readsProcessedThisPass);}
					outstream.println("Clumps:       \t"+clumpsProcessedThisPass);
					if(correct){
						outstream.println("Corrections:  \t"+correctionsThisPass);
					}
					outstream.println();
				}
				
				if(passes>1 && groups==1){
					
					FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=false;
					FASTQ.ASCII_OFFSET=FASTQ.ASCII_OFFSET_OUT;
					
					if(verbose){outstream.println("Pass 1.");}
					if(verbose){outstream.println("Reads:        \t"+readsProcessedThisPass);}
					outstream.println("Clumps:       \t"+clumpsProcessedThisPass);
					if(correct){
						outstream.println("Corrections:  \t"+correctionsThisPass);
					}
					outstream.println();
					
					for(int pass=1; pass<passes; pass++){
						
						if(pass>=conservativePasses){Clump.setConservative(false);}
						
						kc=new KmerComparator(k, kc.seed<0 ? -1 : kc.seed+1, kc.border-1, kc.hashes, false, kc.rcompReads);
						reads=runOnePass(reads, kc);

						if(verbose){outstream.println("Seed: "+kc.seed);}
						if(verbose){outstream.println("Pass "+(pass+1)+".");}
						if(verbose){outstream.println("Reads:        \t"+readsProcessedThisPass);}
						outstream.println("Clumps:       \t"+clumpsProcessedThisPass);
						if(correct){
							outstream.println("Corrections:  \t"+correctionsThisPass);
						}
						outstream.println();
					}
				}
			}
			
			if(repair || namesort){
				if(groups>1){
					if(verbose){t.start("Name-sorting.");}
					reads=nameSort(reads, false);
					if(verbose){t.stop("Sort time: ");}
				}else{
					if(namesort){
						if(verbose){t.start("Name-sorting.");}
						reads=idSort(reads, repair);
						if(verbose){t.stop("Sort time: ");}
					}else{
						reads=read1Only(reads);
					}
				}
			}
			
			if(ros!=null){
				assert(reads.size()==readsProcessedThisPass || dedupe || condense ||
						(reads.size()*2==readsProcessedThisPass && repair)) :
					reads.size()+", "+readsProcessedThisPass;
				if(verbose){t.start("Writing.");}
				for(Read r : reads){
					readsOut+=r.pairCount();
					basesOut+=r.pairLength();
				}
				ros.add(reads, 0);
			}
		}
		
		if(ros!=null){
			if(verbose){outstream.println("Waiting for writing to complete.");}
			errorState=ReadWrite.closeStream(ros)|errorState;
			if(verbose){t.stop("Write time: ");}
		}
		
		if(verbose){outstream.println("Done!");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private final FileFormat ffin1[];
	private final FileFormat ffin2[];

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
}
