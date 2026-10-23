package clump2;

import java.io.File;
import java.util.ArrayList;

import bloom.KCountArray;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.LightweightWriterFactory;
import stream.MultiWriterPolicy;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import tracker.ReadStats;

/**
 * Multi-group sorter used after KmerSplit has written partitioned temp files.
 * Processes each group file into memory in sequence, clumps the reads, then
 * either writes them to one output set or re-hashes them for the next pass.
 *
 * @author Brian Bushnell
 * @date June 20, 2014
 */
public class KmerSort2 extends KmerSort{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Command-line entry point that saves and restores zip statics around one
	 * in-JVM sorter invocation.
	 *
	 * @param args Command-line arguments specifying input files, parameters, and processing options
	 */
	public static void main(String[] args){
		final boolean pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		final float ztd=ReadWrite.ZIP_THREAD_MULT();
		final int mzt=ReadWrite.MAX_ZIP_THREADS();
		final int oldzl=ReadWrite.ZIPLEVEL;
		Timer t=new Timer();
		KmerSort2 x=new KmerSort2(args);
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
	 * Parses one grouped sorter invocation, expands grouped input and output
	 * file names, and configures shared kmer, clump and read-normalization
	 * options before constructing per-group file formats.
	 *
	 * @param args Command-line arguments containing configuration parameters
	 */
	public KmerSort2(String[] args){
		
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
			}else if(a.equals("reorder") || a.equals("reorderclumps")){
				//reorder=Parse.parseBoolean(b);
			}else if(a.equals("reorderpaired") || a.equals("reorderclumpspaired")){
//				reorderpaired=Parse.parseBoolean(b);
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
		
		//Accept paired-input shorthand unless a literal input with '#' exists.
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		//Two explicit inputs carry mates in parallel files, not as interleaved records.
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		//Expand paired-output shorthand before creating final output descriptors.
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(out1==null && out2!=null){
			throw new RuntimeException("Error - cannot define out2 without defining out1.");
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		if(out1==null){ffout1=ffout2=null;}
		else{
			int g=out1.contains("%") ? groups : 1;
			ffout1=new FileFormat[g];
			ffout2=new FileFormat[g];
			if(g==1){
				ffout1[0]=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, out2!=null);
				ffout2[0]=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, out2!=null);
			}else{
				ReadWrite.ZIPLEVEL=2;
				ReadWrite.setZipThreadMult(Tools.min(0.5f, 2f/(g+1)));
				for(int i=0; i<g; i++){
					ffout1[i]=FileFormat.testOutput(out1.replaceFirst("%", ""+i), FileFormat.FASTQ, extout, (g<10), overwrite, append, out2!=null);
					ffout2[i]=(out2==null ? null : FileFormat.testOutput(out2.replaceFirst("%", ""+i), FileFormat.FASTQ, extout, (g<10), overwrite, append, true));
				}
			}
		}
		
		//Grouped callers pass a '%' template naming files split by KmerSplit.
		if(groups>1 && in1.contains("%") && (splitInput || !new File(in1).exists())){
			ffin1=new FileFormat[groups];
			ffin2=new FileFormat[groups];
			for(int i=0; i<groups; i++){
				ffin1[i]=FileFormat.testInput(in1.replaceFirst("%", ""+i), FileFormat.FASTQ, extin, true, true);
				ffin2[i]=in2==null ? null : FileFormat.testInput(in2.replaceFirst("%", ""+i), FileFormat.FASTQ, extin, true, true);
			}
		}else{
			//Ordinary single-group mode reads one input stream and writes one output set.
			assert(!in1.contains("%") && groups==1) : "The % symbol must only be present in the input filename if groups>1.";
			ffin1=new FileFormat[1];
			ffin1[0]=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
			ffin2=new FileFormat[1];
			ffin2[0]=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
			groups=1;
		}
		//FIXED [clump/KmerSort2#001]: reject in/out aliases by expanded names before startup.
		ClumpTools.validateFileNames(ffin1, ffin2, ffout1, ffout2);
//		if(groups>1){ReadWrite.USE_UNPIGZ=false;} //Not needed since they are not concurrent
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Builds optional pivot counts, opens one output stream per configured
	 * group or final destination, processes each input group, and prints stats.
	 *
	 * @param t Timer for tracking execution performance
	 */
	@Override
	void process(Timer t){
		
		preprocess();

		final Writer[] rosa=(ffout1==null ? null : new Writer[ffout1.length]);
		if(rosa!=null){
			if(rosa.length==1){
				rosa[0]=WriterFactory.getStream(ffout1[0], ffout2[0], null, null, 4, null, false, -1);
				rosa[0].start();
			}else{
				final int physicalOutputs=countPhysicalOutputs(ffout1, ffout2);
				final int sharedThreads=Shared.threads();
				final int requestedMode=MultiWriterPolicy.DEFAULT_MODE;
				final int writerMode=(requestedMode==MultiWriterPolicy.AUTO ?
						(physicalOutputs<=sharedThreads ? MultiWriterPolicy.ST : MultiWriterPolicy.ZT) : requestedMode);
				outstream.println("Clump shard writer mode: "+(writerMode==MultiWriterPolicy.ST ? "ST" : "ZT")+
						" for "+physicalOutputs+" physical outputs with t="+sharedThreads);
				final LightweightWriterFactory ztf=new LightweightWriterFactory(writerMode, 2);
				for(int i=0; i<rosa.length; i++){
					rosa[i]=ztf.getStream(ffout1[i], ffout2[i], null, null, null, false);
					rosa[i].start();
				}
			}
		}
		
		readsProcessed=basesProcessed=diskProcessed=memProcessed=0;
		
		//Process the read stream
		processInner(rosa);
		
		printStats(t);
	}

	/**
	 * Counts the current grouped temp descriptors into one table using the
	 * current comparator statics.
	 */
	@Override
	protected KCountArray countCurrentGroupedFilePivots(){
		return ClumpTools.getTable(ffin1, ffin2, k, minCount);
	}
	
	/**
	 * Processes groups sequentially in memory: fetch and hash reads, sort by
	 * pivot kmer, build clumps, apply the selected clump operation, optionally
	 * run same-list passes, and route the surviving reads to the output array.
	 *
	 * @param rosa Array of output streams for writing processed reads
	 */
	void processInner(final Writer rosa[]){
		if(verbose){outstream.println("Making comparator.");}
		KmerComparator kc=new KmerComparator(k, addName, (rcomp || condense || correct));
		
		ClumpList.UNRCOMP=(!rcomp && !condense);
		Timer t=new Timer();
		final long[] outIds=(rosa==null ? null : new long[rosa.length]);
		
//		final int conservativePasses=Clump.conservativeFlag ? passes : Tools.max(1, passes/2);
//		if(groups==1 && passes>1){Clump.setConservative(true);}

		useSharedHeader=(ffin1[0].samOrBam() && ffout1!=null && ffout1[0]!=null && ffout1[0].samOrBam());
		
		for(int group=0; group<groups; group++){
			if(verbose){outstream.println("Starting cris "+group+".");}
			
			final Streamer cris=StreamerFactory.getReadInputStream(maxReads, false, ffin1[group], ffin2[group], -1);
			cris.start();
			
//			if(verbose){t.start("Fetching reads.");}
			ArrayList<Read> reads=fetchReads1(cris, kc);
//			if(verbose){t.stop("Fetch time: ");}
			
			if(verbose){t.start("Sorting.");}
			Shared.sort(reads, kc);
			if(verbose){t.stop("Sort time: ");}
			
//			if(verbose){t.start("Counting clumps.");}
//			clumpsProcessed+=countClumps(reads);
//			if(verbose){t.stop("Count time: ");}
			
			if(verbose){t.start("Making clumps.");}
			readsProcessedThisPass=reads.size();
			ClumpList cl=new ClumpList(reads, k, false);
			clumpsProcessedThisPass=cl.size();
			clumpsProcessedTotal+=clumpsProcessedThisPass;
			if(verbose){t.stop("Clump time: ");}
			
			if(dedupe){
				reads.clear();
				if(verbose){t.start("Deduping.");}
				reads=processClumps(cl, ClumpList.DEDUPE);
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
						
//						if(pass>=conservativePasses){Clump.setConservative(false);}
						
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
						reads=nameSortById(reads, repair && unpair);
						if(verbose){t.stop("Sort time: ");}
					}else{
						reads=read1Only(reads);
					}
				}
			}
			
			for(Read r : reads){
				readsOut+=r.pairCount();
				basesOut+=r.pairLength();
			}
			
			//Two output models. NORMAL path (doHashAndSplit defaults true -> addToRos): each input group's reads are
			//re-hashed and SPLIT across all output files rosa[] (multi-group memory-efficient writeout). The else
			//path (one-group-per-output-file, streaming) is reachable when Clumpify.runOnePass_v2 sets
			//doHashAndSplit=false for externalSort before selecting KmerSort2 versus KmerSort3. It progressively
			//closes rosa[group-1] here AND the final loop below closes ALL rosa, so rosa[0..groups-2] are closed
			//twice; that relies on output close idempotence. [clump/KmerSort2#002, latent]
			if(doHashAndSplit || groups==0){
				addToRos(rosa, reads, t, kc, outIds);
			}else{
				if(group>0){
					Writer ros=rosa[group-1];
					errorState|=ReadWrite.closeOutputStreams(ros);
				}
				if(!reads.isEmpty()){rosa[group].add(reads, outIds[group]++);}
			}
			reads=null;
		}
		
		if(rosa!=null){
			if(verbose){outstream.println("Waiting for writing to complete.");}
			for(Writer ros : rosa){
				errorState=ReadWrite.closeOutputStreams(ros)|errorState;
			}
			if(verbose){t.stop("Write time: ");}
		}
		
		if(verbose){outstream.println("Done!");}
	}
	
	//CLEVER [verified]: the V2 multi-group writeout. With >1 output file, re-hash every read with a FRESH comparator
	//(seed+1, so the output-split hash is independent of the pivot-selection hash) and hashAndSplit them into div
	//buckets by hash(kmer)%div, one bucket per output file. This is what makes the groups>1 memory-efficient
	//pipeline work: each output file becomes a self-contained kmer-partition for a later pass. div==groups asserted.
	private void addToRos(Writer[] rosa, ArrayList<Read> list, Timer t, KmerComparator old, long[] outIds){
		if(rosa==null){return;}
		assert(rosa.length>0);
		if(rosa.length==1){
			if(verbose){t.start("Writing.");}
			if(!list.isEmpty()){rosa[0].add(list, outIds[0]++);}
			return;
		}
		KmerComparator kc=new KmerComparator(old.k, old.seed+1, old.border-1, old.hashes, false, false);
		final int div=rosa.length;
		assert(div==groups);
		@SuppressWarnings("unchecked")
		ArrayList<Read>[] array=new ArrayList[div];
		for(int i=0; i<array.length; i++){
			array[i]=new ArrayList<Read>();
		}
		if(verbose){t.start("Splitting.");}
		hashAndSplit(list, kc, array);
		if(verbose){t.stop("Split time: ");}
		if(verbose){t.start("Writing.");}
		for(int i=0; i<div; i++){
			if(!array[i].isEmpty()){rosa[i].add(array[i], outIds[i]++);}
			array[i]=null;
		}
		if(verbose){System.err.println("Sent writable reads.");}
	}
	
	/** Counts nonnull primary and mate shard outputs open at the same time. */
	private static int countPhysicalOutputs(final FileFormat[] ffout1, final FileFormat[] ffout2){
		int count=0;
		for(int i=0; i<ffout1.length; i++){
			if(ffout1[i]!=null){count++;}
			if(ffout2!=null && ffout2[i]!=null){count++;}
		}
		return count;
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

	/** Primary input descriptor per group; one element in unsplit mode. */
	private final FileFormat ffin1[];

	/** Optional mate input descriptor per group; null elements mean interleaved
	 * or unpaired input. */
	private final FileFormat ffin2[];

	/** Primary final or grouped output descriptor per output stream. */
	private final FileFormat ffout1[];

	/** Optional mate final or grouped output descriptor per output stream. */
	private final FileFormat ffout2[];
	
}
