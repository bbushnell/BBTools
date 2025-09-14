package clump;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.atomic.AtomicInteger;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
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
import structures.ListNum;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public class KmerSort3 extends KmerSort {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		assert(false);
		main(-1, 1, 1, args);
	}
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(long fileMem, int outerPassNum, int outerPasses, String[] args){
		final boolean pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		final float ztd=ReadWrite.ZIP_THREAD_MULT();
		final int mzt=ReadWrite.MAX_ZIP_THREADS();
		final int oldzl=ReadWrite.ZIPLEVEL;
		Timer t=new Timer();
		KmerSort3 x=new KmerSort3(fileMem, outerPassNum, outerPasses, args);
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
	public KmerSort3(long fileMem_, int outerPassNum_, int outerPasses_, String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		outerPassNum=outerPassNum_;
		outerPasses=outerPasses_;
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
			}else if(a.equals("fetchthreads")){
				fetchThreads=Integer.parseInt(b);
				assert(fetchThreads>0) : fetchThreads+"\nFetch threads must be at least 1.";
			}else if(a.equals("reorder") || a.equals("reorderclumps")){
				//reorder=Parse.parseBoolean(b);
			}else if(a.equals("reorderpaired") || a.equals("reorderclumpspaired")){
//				reorderpaired=Parse.parseBoolean(b);
			}
			
			else if(Clump.parseStatic(arg, a, b)){
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
		
		if(out1==null){ffout1=ffout2=null;}
		else{
			int g=out1.contains("%") ? groups : 1;
			ffout1=new FileFormat[g];
			ffout2=new FileFormat[g];
			if(g==1){
				ffout1[0]=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
				ffout2[0]=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
			}else{
				ReadWrite.ZIPLEVEL=2;
				ReadWrite.setZipThreadMult(Tools.min(0.5f, 2f/(g+1)));
				for(int i=0; i<g; i++){
					ffout1[i]=FileFormat.testOutput(out1.replaceFirst("%", ""+i), FileFormat.FASTQ, extout, (g<10), overwrite, append, false);
					ffout2[i]=(out2==null ? null : FileFormat.testOutput(out2.replaceFirst("%", ""+i), FileFormat.FASTQ, extout, (g<10), overwrite, append, false));
				}
			}
		}
		
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
		
		long sizeSum=0, expectedMemSum=0;
		for(FileFormat ff : ffin1){
			long x=new File(ff.name()).length();
			sizeSum+=x;
			if(ff.compressed()){x*=40;}else{x*=8;}
			expectedMemSum+=x;
		}
		for(FileFormat ff : ffin2){
			if(ff!=null){
				long x=new File(ff.name()).length();
				sizeSum+=x;
				if(ff.compressed()){x*=40;}else{x*=8;}
				expectedMemSum+=x;
			}
		}

		expectedSizePerGroup=(sizeSum+groups+1)/(Tools.max(groups, 1));
		expectedMemPerGroup=(expectedMemSum+groups+1)/(Tools.max(groups, 1));
		totalMem=Shared.memAvailable(1);
		fileSize=sizeSum;
		fileMem=fileMem_<1 ? 40*fileSize : fileMem_;
		memRatio=fileMem*1.0/Tools.max(1, fileSize);
		
//		if(groups>1){ReadWrite.USE_UNPIGZ=false;} //Not needed since they are not concurrent
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	void process(Timer t){
		
		preprocess();

		final ConcurrentReadOutputStream[] rosa=(ffout1==null ? null : new ConcurrentReadOutputStream[ffout1.length]);
		for(int i=0; rosa!=null && i<rosa.length; i++){
			final int buff=1;

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			rosa[i]=ConcurrentReadOutputStream.getStream(ffout1[i], ffout2[i], null, null, buff, null, false);
			rosa[i].start();
		}
		
		readsProcessed=basesProcessed=diskProcessed=memProcessed=0;
		
		//Process the read stream
		processInner(rosa);
		lastMemProcessed=memThisPass;
		
		printStats(t);
		
		if(outerPassNum<outerPasses){outstream.println();}
	}
	
	/** Collect and sort the reads */
	void processInner(final ConcurrentReadOutputStream rosa[]){
		if(verbose){outstream.println("Making comparator.");}
		KmerComparator kc=new KmerComparator(k, addName, (rcomp || condense || correct));
		
		ClumpList.UNRCOMP=(!rcomp && !condense);
		Timer t=new Timer();
		
//		final int conservativePasses=Clump.conservativeFlag ? passes : Tools.max(1, passes/2);
//		if(groups==1 && passes>1){Clump.setConservative(true);}

		useSharedHeader=(ffin1[0].samOrBam() && ffout1!=null && ffout1[0]!=null && ffout1[0].samOrBam());
		
		fetchThreads=Tools.min(groups, fetchThreads, ffin1.length);
		assert(fetchThreads>0);
		SynchronousQueue<ArrayList<Read>> listQ=new SynchronousQueue<ArrayList<Read>>();
		ArrayList<FetchThread3> alft=fetchReads(kc, fetchThreads, listQ, rosa);
		int poisonCount=0;
		
		for(int group=0; group<groups; group++){
			
			if(verbose){t.start("Fetching reads.");}
			ArrayList<Read> reads=null;
			
			//TODO: There appears to be something unsafe here which could lead to this loop being skipped.
			while(poisonCount<fetchThreads && reads==null){
				try {
					reads=listQ.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(reads==POISON){
					assert(reads.isEmpty());
					poisonCount++;
					if(verbose){System.err.println("Encountered poison; count="+poisonCount);}
					reads=null;
				}
			}
			if(verbose){t.stop("Fetched "+(reads==null ? 0 : reads.size())+" reads: ");}
			
			//Done by FetchThread3
//			if(verbose){t.start("Sorting.");}
//			Shared.sort(reads, kc);
//			if(verbose){t.stop("Sort time: ");}
			
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
						reads=idSort(reads, repair);
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
			
			if(doHashAndSplit || groups==0){
				addToRos(rosa, reads, t, kc);
			}else{
				if(group>0){
					ConcurrentReadOutputStream ros=rosa[group-1];
					errorState|=ReadWrite.closeStream(ros);
				}
				rosa[group].add(reads, 0);
			}
			reads=null;
		}
		
		if(verbose){outstream.println("Closing fetch threads.");}
		while(poisonCount<fetchThreads){
			ArrayList<Read> reads=null;
			try {
				reads=listQ.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(reads==POISON){
				poisonCount++;
				reads=null;
			}
		}
		long readsThisPass=closeFetchThread3s(alft);
		if(verbose){outstream.println("Closed fetch threads.");}
		
		if(rosa!=null){
			if(verbose){outstream.println("Waiting for writing to complete.");}
			for(ConcurrentReadOutputStream ros : rosa){
				errorState=ReadWrite.closeStream(ros)|errorState;
			}
			if(verbose){t.stop("Write time: ");}
		}
		
		if(verbose && outerPassNum==outerPasses){outstream.println("Done!");}
	}
	
	private void addToRos(ConcurrentReadOutputStream[] rosa, ArrayList<Read> list, Timer t, KmerComparator old){
		if(rosa==null){return;}
		assert(rosa.length>0);
		if(rosa.length==1){
			if(verbose){t.start("Writing.");}
			rosa[0].add(list, 0);
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
			rosa[i].add(array[i], 0);
			array[i]=null;
		}
		if(verbose){System.err.println("Sent writable reads.");}
	}
	
	public ArrayList<FetchThread3> fetchReads(final KmerComparator kc, final int fetchThreads, SynchronousQueue<ArrayList<Read>> listQ, ConcurrentReadOutputStream[] rosa){
		AtomicInteger nextGroup=new AtomicInteger(0);
		if(verbose){outstream.println("Making "+fetchThreads+" fetch thread"+(fetchThreads==1 ? "." : "s."));}
		ArrayList<FetchThread3> alft=new ArrayList<FetchThread3>(fetchThreads);
		for(int i=0; i<fetchThreads; i++){alft.add(new FetchThread3(kc, listQ, nextGroup, rosa));}
		
		if(verbose){outstream.println("Starting threads.");}
		for(FetchThread3 ft : alft){ft.start();}
		
		assert(alft.size()==fetchThreads);
		
		return alft;
	}
	

	private long closeFetchThread3s(ArrayList<FetchThread3> alft){
		readsThisPass=0;
		memThisPass=0;
		entryFilteredThisPass=0;
		/* Wait for threads to die */
		for(FetchThread3 ft : alft){

			/* Wait for a thread to die */
			while(ft.getState()!=Thread.State.TERMINATED){
				try {
					ft.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			entryFilteredThisPass+=ft.entryFilteredT;
			readsThisPass+=ft.readsProcessedT;
			basesProcessed+=ft.basesProcessedT;
			diskProcessed+=ft.diskProcessedT;
			memThisPass+=ft.memProcessedT;
			
			errorState|=ft.errorStateT;
		}
		entryFiltered+=entryFilteredThisPass;
		readsProcessed+=readsThisPass;
		memProcessed+=memThisPass;
		
		ecco=false;
		return readsThisPass;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class FetchThread3 extends Thread{
		
		FetchThread3(final KmerComparator kc_, SynchronousQueue<ArrayList<Read>> listQ_, AtomicInteger nextGroup_, ConcurrentReadOutputStream[] rosa_){
			kc=kc_;
			listQ=listQ_;
			nextGroup=nextGroup_;
			rosa=rosa_;
		}
		
		@Override
		public void run(){
			for(ArrayList<Read> reads=fetchNext(); reads!=null; reads=fetchNext()){
				boolean success=false;
				while(!success){
					try {
						listQ.put(reads);
						success=true;
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			if(verbose){System.err.println("No more reads to fetch.");}
			boolean success=false;
			while(!success){
				try {
					if(verbose){System.err.println("Adding poison.");}
					listQ.put(POISON);
					success=true;
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(verbose){System.err.println("A fetch thread finished.");}
		}
		
		private ArrayList<Read> fetchNext(){
			final int group=nextGroup.getAndIncrement();
			if(group>=groups){return null;}
			
//			assert(false) : ffin1[group]+", "+FASTQ.FORCE_INTERLEAVED+", "+FASTQ.TEST_INTERLEAVED;
			
			final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ffin1[group], ffin2[group], null, null);
			cris.start();
			
			//Check for file size imbalance
			if(!Clump.forceProcess){
				final long size;
				double expectedMem;
				{
					size=new File(ffin1[group].name()).length()+(ffin2[group]==null ? 0 : new File(ffin2[group].name()).length());
//					expectedMem=size*(ffin1[group].compressed() ? 40 : 8);
					expectedMem=(fileMem*(double)size)/fileSize;
				}
				
				if(size-1000>expectedSizePerGroup*3 || expectedMem*3>totalMem){
					
//					outstream.println("size="+size+", expectedSizePerGroup="+expectedSizePerGroup+", expectedMem="+expectedMem+", expectedSizePerGroup="+expectedMemPerGroup);
//					outstream.println("totalMem="+totalMem+", Shared.memFree()="+Shared.memFree()+", Shared.memUsed()="+Shared.memUsed());
//					outstream.println("fileSize="+fileSize+", "+"fileMem="+fileMem+", "+"memRatio="+memRatio+", ");
//					outstream.println((size-1000>expectedSizePerGroup*3)+", "+(expectedMem*3>totalMem));

					//TODO: This could also be based on the number of FetchThreads.
					if(expectedMem>0.3*totalMem && (fileMem<1 || fileMem>0.8*totalMem)){
						synchronized(this.getClass()){
							if(!memWarned){
								outstream.println("\n***  Warning  ***\n"
										+ "A temp file may be too big to store in memory, due to uneven splitting:");

								outstream.println("expectedMem="+((long)expectedMem)+", fileMem="+fileMem+", available="+totalMem);
							}
							if(entryfilter){
								if(!memWarned){
									outstream.println(
											"If the program crashes, increase memory or manually increase groups.\n");
								}
							}else if(dedupe & !Clump.opticalOnly){
								if(!memWarned){
									outstream.println(
											"entryfilter is automatically being set to true to reduce memory load while removing duplicates.\n"
													+"If the program crashes, increase memory or manually increase groups.\n");
								}
								entryfilter=true;
							}else if(repair || namesort){
								if(!memWarned){
									outstream.println("It cannot be streamed to output unaltered because "+(namesort ? " namesort=t" : "repair=t"));
									outstream.println("If this causes the program to crash, please re-run with more memory or groups.\n");
								}
							}else{
								if(!memWarned){
									outstream.println(
											"It will be streamed to output unaltered.\n"
													+ "To avoid this behavior, increase memory or increase groups.\n"
													+ "Set the flag forceprocess=t to disable this check.\n");
								}
								//							Timer t=new Timer();
								ArrayList<Read> list=streamNext_inner(cris);
								//							t.stop("Stream time: ");
								return list;
							}
							memWarned=true;
						}
					}
				}
			}
			
			return fetchNext_inner(cris);
		}
		
		private ArrayList<Read> streamNext_inner(ConcurrentReadInputStream cris){
			StreamToOutput sto=new StreamToOutput(cris, rosa, kc, (repair || namesort), false);
			errorStateT|=sto.process();
			readsProcessed+=sto.readsIn;
			basesProcessed+=sto.basesIn;
			readsOut+=sto.readsIn;
			basesOut+=sto.basesIn;
//			System.err.println(readsProcessed+", "+sto.readsIn);
			return new ArrayList<Read>();
		}
		
		private ArrayList<Read> fetchNext_inner(ConcurrentReadInputStream cris){
			
//			Timer t=new Timer();
//			if(verbose){t.start("Making hash threads.");}
			final int subthreads=Tools.mid(1, Shared.threads()/2, 8);
			ArrayList<FetchSubThread> alht=new ArrayList<FetchSubThread>(subthreads);
			long readsThisGroup=0, memThisGroup=0, entryFilteredThisGroup=0;
			for(int i=0; i<subthreads; i++){alht.add(new FetchSubThread(i, cris, kc, unpair));}
			
//			if(verbose){outstream.println("Starting threads.");}
			for(FetchSubThread ht : alht){ht.start();}
			
//			if(verbose){outstream.println("Waiting for threads.");}
			/* Wait for threads to die */
			for(FetchSubThread ht : alht){
				
				/* Wait for a thread to die */
				while(ht.getState()!=Thread.State.TERMINATED){
					try {
						ht.join();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
				entryFilteredThisGroup+=ht.entryFilteredT;
				readsThisGroup+=ht.readsProcessedST;
				basesProcessedT+=ht.basesProcessedST;
				diskProcessedT+=ht.diskProcessedST;
				memThisGroup+=ht.memProcessedST;
			}
			entryFilteredT+=entryFilteredThisGroup;
			readsProcessedT+=readsThisGroup;
			memProcessedT+=memThisGroup;

//			if(verbose){t.stop("Hash time: ");}
//			if(verbose){System.err.println("Closing input stream.");}
			errorStateT=ReadWrite.closeStream(cris)|errorStateT;
			
//			if(verbose){t.start("Combining thread output.");}
			long readsLeft=readsThisGroup-entryFilteredThisGroup;
			long slotsLeft=cris.paired() && !unpair ? readsLeft/2 : readsLeft;
			assert(slotsLeft<=Shared.MAX_ARRAY_LEN && readsProcessedT>=0);
			ArrayList<Read> list=new ArrayList<Read>((int)slotsLeft);
			for(int i=0; i<subthreads; i++){
				FetchSubThread ht=alht.set(i, null);
				list.addAll(ht.storage);
			}
			assert(list.size()==slotsLeft) : list.size()+", "+readsThisGroup+", "+slotsLeft+", "+cris.paired();
//			if(verbose){t.stop("Combine time: ");}
			
//			if(verbose){t.start("Sorting.");}
			Shared.sort(list, kc);
			
//			if(verbose){t.stop("Sort time: ");}
			return list;
		}
		
		final SynchronousQueue<ArrayList<Read>> listQ;
		final AtomicInteger nextGroup;
		final KmerComparator kc;
		final ConcurrentReadOutputStream[] rosa;

		protected long entryFilteredT=0;
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		protected long diskProcessedT=0;
		protected long memProcessedT=0;
		protected boolean errorStateT=false;
		
		
		private class FetchSubThread extends Thread{
			
			FetchSubThread(int id_, ConcurrentReadInputStream cris_, KmerComparator kc_, boolean unpair_){
				id=id_;
				cris=cris_;
				kcT=kc_;
				storage=new ArrayList<Read>();
				unpairT=unpair_;
				entryFilterTable=(entryfilter ? new HashMap<Long, Read>() : null);
			}
			
			@Override
			public void run(){
				ListNum<Read> ln=cris.nextList();
				final boolean paired=cris.paired();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					
					for(Read r : reads){
						if(!r.validated()){
							r.validate(true);
							if(r.mate!=null){r.mate.validate(true);}
						}
						readsProcessedST+=1+r.mateCount();
						basesProcessedST+=r.length()+r.mateLength();
						diskProcessedST+=r.countFastqBytes()+r.countMateFastqBytes();
						memProcessedST+=r.countBytes()+r.countMateBytes()+ReadKey.overhead;
						if(shrinkName){
							Clumpify.shrinkName(r);
							Clumpify.shrinkName(r.mate);
						}else if(shortName){
							Clumpify.shortName(r);
							Clumpify.shortName(r.mate);
						}
					}
					
					if(ecco){
						for(Read r : reads){
							Read r2=r.mate;
							assert(r.obj==null) : "TODO: Pivot should not have been generated yet, though it may be OK.";
							assert(r2!=null) : "ecco requires paired reads.";
							if(r2!=null){
								int x=BBMerge.findOverlapStrict(r, r2, true);
								if(x>=0){
									r.obj=null;
									r2.obj=null;
								}
							}
						}
					}
					
					if(entryFilterTable!=null){
						int removed=0;
						for(int i=0; i<reads.size(); i++){
							Read r=reads.get(i);
							final long key=Hasher.hashPair(r);
							final Long key2=Long.valueOf(key);
							final Read old=entryFilterTable.get(key2);
							if(old==null){
								entryFilterTable.put(key2, r);
							}else{
								boolean same=Hasher.equalsPaired(r, old);
								if(same){
									removed++;
									entryFilteredT+=r.pairCount();
									reads.set(i, null);
								}
							}
						}
						if(removed>0){Tools.condenseStrict(reads);}
					}
					
					ArrayList<Read> hashList=reads;
					if(paired && unpairT){
						hashList=new ArrayList<Read>(reads.size()*2);
						for(Read r1 : reads){
							Read r2=r1.mate;
							assert(r2!=null);
							hashList.add(r1);
							hashList.add(r2);
							if(groups>1 || !repair || namesort){
								r1.mate=null;
								r2.mate=null;
							}
						}
					}
					
					kcT.hash(hashList, table, minCount, true);
					storage.addAll(hashList);
					cris.returnList(ln.id, false);
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				if(ln!=null){
					cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
				}
				
				//Optimization for TimSort
				if(parallelSort){
					storage.sort(kcT);
//					Shared.sort(storage, kc); //Already threaded; this is not needed.
				}else{
					Collections.sort(storage, kcT);
				}
			}

			final int id;
			final ConcurrentReadInputStream cris;
			final KmerComparator kcT;
			final ArrayList<Read> storage;
			final boolean unpairT;
			final HashMap<Long, Read> entryFilterTable;
			public long entryFilteredT=0;

			protected long readsProcessedST=0;
			protected long basesProcessedST=0;
			protected long diskProcessedST=0;
			protected long memProcessedST=0;
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	protected static long lastMemProcessed=0;
	
	final long expectedSizePerGroup;
	private final long expectedMemPerGroup;
	final long totalMem;
	final long fileMem;
	final long fileSize;
	
	private final int outerPassNum;
	private final int outerPasses;
	
	private final double memRatio;
	
	/*--------------------------------------------------------------*/

	static final ArrayList<Read> POISON=new ArrayList<Read>();
	protected static int fetchThreads=2;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	final FileFormat ffin1[];
	final FileFormat ffin2[];

	private final FileFormat ffout1[];
	private final FileFormat ffout2[];
	
}
