package sort;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.PriorityQueue;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.CrisContainer;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.ListNum;
import tax.AccessionToTaxid;
import tax.GiToTaxid;
import tax.TaxTree;
import tracker.ReadStats;
import var2.ScafMap;

/**
 * Sorts reads by name, potentially from multiple input files.
 * 
 * @author Brian Bushnell
 * @date September 21, 2016
 *
 */
public class SortByName {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		final boolean oldFI=FASTQ.FORCE_INTERLEAVED, oldTI=FASTQ.TEST_INTERLEAVED;
		SortByName x=new SortByName(args);
		x.process(t);
		FASTQ.FORCE_INTERLEAVED=oldFI;
		FASTQ.TEST_INTERLEAVED=oldTI;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SortByName(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Set shared static variables
		Shared.capBuffers(4);
		Shared.setBufferData(1000000);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=true; //TEST_INTERLEAVED must explicitly be disabled with int=f to sort corrupt files.
		FASTQ.FORCE_INTERLEAVED=false;
		
		//Create a parser object
		Parser parser=new Parser();
		boolean ascending=true;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("verbose2")){
				assert(false) : "Verbose2 is disabled.";
//				verbose2=Parse.parseBoolean(b);
			}else if(a.equals("delete")){
				delete=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("deleteearly")){
				deleteEarly=Parse.parseBoolean(b);
			}else if(a.equals("allowtemp") || a.equals("usetemp")){
				allowTempFiles=Parse.parseBoolean(b);
			}else if(a.equals("memmult") || a.equals("mult")){
				memBlockMult=(float) Double.parseDouble(b);
			}else if(a.equals("memlimit")){
				memTotalMult=(float) Double.parseDouble(b);
			}else if(a.equals("ascending")){
				ascending=Parse.parseBoolean(b);
			}else if(a.equals("descending")){
				ascending=!Parse.parseBoolean(b);
			}else if(a.equals("length")){
				if(Parse.parseBoolean(b)){
					comparator=ReadLengthComparator.comparator;
				}
			}else if(a.equals("name")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorName.comparator;
				}
			}else if(a.equals("quality")){
				if(Parse.parseBoolean(b)){
					comparator=ReadQualityComparator.comparator;
				}
			}else if(a.equals("position")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorPosition.comparator;
				}
			}else if(a.equals("list") || a.equals("names")){
				comparator=new ReadComparatorList(b);
			}else if(a.equals("random") || a.equals("shuffle")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorRandom.comparator;
				}
			}else if(a.equals("taxa")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorTaxa.comparator;
				}
			}else if(a.equals("topo") || a.equals("topological") || a.equals("alphabetic") || a.equals("sequence") || a.equals("lexographic")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorTopological.comparator;
				}
			}else if(a.equals("flowcell")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorFlowcell.comparator;
				}
			}else if(a.equals("crispr")){
				if(Parse.parseBoolean(b)){
					ReadComparatorCrispr.loadNet();
					comparator=ReadComparatorCrispr.comparator;
				}
			}else if(a.equalsIgnoreCase("clump") || a.equalsIgnoreCase("clumpify")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorClump.comparator;
				}
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				if(b==null || "ignore".equalsIgnoreCase(b) || "skip".equalsIgnoreCase(b)){
					giTableFile=null;
					TaxTree.CRASH_IF_NO_GI_TABLE=false;
				}else{giTableFile=b;}
			}else if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
			}else if(a.equals("maxfiles") || a.equals("files")){
				maxFiles=Integer.parseInt(b);
				assert(maxFiles>1 && maxFiles<=10000000) : "Invalid value for maxfiles: "+maxFiles+"; range is 2-10000000";
				setMaxFiles=true;
			}else if(a.equals("genkmer")){
				genKmer=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		
		if(comparator==ReadComparatorTopological.comparator && genKmer){
			comparator=ReadComparatorTopological5Bit.comparator;
		}else{
			genKmer=false;
		}
		
		clump=(comparator==ReadComparatorClump.comparator);
		
		comparator.setAscending(ascending);
		SamLine.SET_FROM_OK=true;
		
		if(comparator==ReadComparatorRandom.comparator){
			ListNum.setDeterministicRandomSeed(-1);
			ListNum.setDeterministicRandom(true);
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			
			extin=parser.extin;
			extout=parser.extout;
			
			minlen=parser.minReadLength;
		}
		
//		assert(false) : setInterleaved+", "+FASTQ.FORCE_INTERLEAVED+", "+FASTQ.TEST_INTERLEAVED;
//		assert(!FASTQ.FORCE_INTERLEAVED) : setInterleaved+", "+FASTQ.FORCE_INTERLEAVED+", "+FASTQ.TEST_INTERLEAVED;
//		assert(!FASTQ.TEST_INTERLEAVED) : setInterleaved+", "+FASTQ.FORCE_INTERLEAVED+", "+FASTQ.TEST_INTERLEAVED;
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			if(in2==null && out2!=null){
				FASTQ.FORCE_INTERLEAVED=true;
				FASTQ.TEST_INTERLEAVED=false;
			}else if(in2!=null){
				FASTQ.FORCE_INTERLEAVED=false;
				FASTQ.TEST_INTERLEAVED=false;
			}
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		
		tempExt=Tools.getTempExt(ffin1, ffout1, extout);
		
		if(comparator==ReadComparatorPosition.comparator){
			if(ReadComparatorPosition.scafMap==null){
				ReadComparatorPosition.scafMap=ScafMap.loadSamHeader(in1);
			}
		}
		
		if((comparator==ReadComparatorTaxa.comparator)){
			if(giTableFile!=null){
				outstream.println("Loading gi table.");
				GiToTaxid.initialize(giTableFile);
			}
			if(accessionFile!=null){
				outstream.println("Loading accession table.");
				AccessionToTaxid.load(accessionFile);
			}
			if(taxTreeFile!=null){
				ReadComparatorTaxa.tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, false);
				assert(ReadComparatorTaxa.tree.nameMap()!=null);
			}else{
				throw new RuntimeException("No tree specified.");
			}
		}
		
		qualityObserved=!ffin1.fasta();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			useSharedHeader=(ffin1.samOrBam() && ffout1!=null && ffout1.samOrBam());
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
//		//Optionally create a read output stream
//		final ConcurrentReadOutputStream ros;
//		if(ffout1!=null){
//			final int buff=4;
//
//			if(cris.paired() && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
//				outstream.println("Writing interleaved.");
//			}
//
//			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
//			ros.start(); //Start the stream
//		}else{ros=null;}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the reads */
	public void processInner(final ConcurrentReadInputStream cris){
		//Do anything necessary prior to processing
		final int ziplevel0=ReadWrite.ZIPLEVEL;
		ReadWrite.ZIPLEVEL=Tools.mid(1, ReadWrite.ZIPLEVEL, 2);
		
		ArrayList<Read> storage=new ArrayList<Read>();
		
		final long maxMem=Shared.memAvailable(1);
		final long memLimit=(long)(maxMem*memTotalMult);//Crashed once at old 0.75 value
		final long currentLimit=(long)(maxMem*memBlockMult);
		final int readLimit=2000000000;
		long currentMem=0;
		long bytesDumped=0;
		long dumps=0;
//		IntList dumpCount=new IntList();
		AtomicLong outstandingMem=new AtomicLong();
		
		if(verbose){outstream.println("maxMem="+maxMem+", memLimit="+memLimit+
				", currentMem="+currentMem+", currentLimit="+currentLimit);}
		
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose2){outstream.println("Fetched "+reads.size()+" reads.");}
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;
					maxLengthObserved=Tools.max(maxLengthObserved, initialLength1, initialLength2);
					
					if(minlen<1 || initialLength1>=minlen || initialLength2>=minlen){
						if(genKmer){ReadComparatorTopological5Bit.genKmer(r1);}
						else if(clump){ReadComparatorClump.set(r1);}
						currentMem+=r1.countBytes()+(r2==null ? 0 : r2.countBytes());
						storage.add(r1);
					}
				}
				
				if(allowTempFiles && (currentMem>=currentLimit || storage.size()>=readLimit)){
					if(verbose){outstream.println("currentMem: "+currentMem+" >= "+currentLimit+", dumping. ");}
					outstandingMem.addAndGet(currentMem);
//					dumpCount.add(storage.size());
					sortAndDump(storage, currentMem, outstandingMem, null, null, false);
					storage=new ArrayList<Read>();
					bytesDumped+=currentMem;
					dumps++;
					currentMem=0;
					if(verbose){outstream.println("Waiting on memory; outstandingMem="+outstandingMem);}
					waitOnMemory(outstandingMem, memLimit);
					if(verbose){outstream.println("Done waiting; outstandingMem="+outstandingMem);}
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		outstream.println("Finished reading input.");
		
		outstandingMem.addAndGet(currentMem);
		if(dumps==0){
			ReadWrite.ZIPLEVEL=ziplevel0;
			outstream.println("Sorting.");
			if(out1!=null){
				sortAndDump(storage, currentMem, outstandingMem, out1, out2, useSharedHeader);
				storage=null;
				waitOnMemory(outstandingMem, 0);
			}else{
				Timer t=new Timer();
				Shared.sort(storage, comparator);//For timing/testing.
				t.stop("Sort Time: \t\t\t");
			}
		}else{
//			dumpCount.add(storage.size());
			sortAndDump(storage, currentMem, outstandingMem, null, null, false);
			storage=null;
			waitOnMemory(outstandingMem, 0);
			outstream.println("Merging "+(dumps+1)+" files.");
			ReadWrite.ZIPLEVEL=ziplevel0;
			mergeAndDump(outTemp, /*dumpCount, */useSharedHeader, false);
		}
		
	}

	private void waitOnMemory(AtomicLong outstandingMem, long target){
		if(outstandingMem.get()>target){
			if(verbose){outstream.println("Syncing; outstandingMem="+outstandingMem);}
			while(outstandingMem.get()>target){
				try {
					synchronized(outstandingMem){
						outstandingMem.wait(2000);
					}
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private final long maxSizeObserved(){//Includes some extra for a margin of error
		return Tools.max(maxSizeObservedStatic(), Tools.max(maxLengthObserved, 150)*(qualityObserved ? 2 : 1)+500);
	}
	final static long maxSizeObservedStatic(){//Includes some extra for a margin of error
		return Tools.max(maxLengthObservedStatic, 150)*(qualityObservedStatic ? 2 : 1)+500;
	}
	private static long mergeMemNeeded(long maxSizeObserved, int inFiles){
		return (long)(maxSizeObserved*2*(inFiles+1)*1.6);
	}
	private static int maxMergeFiles(long maxSizeObserved){
		final long memAvailable=Shared.memAvailableAdvanced();
		int files=(int)Tools.max(3, memAvailable/Tools.max(maxSizeObserved, 500))-1;
		return Tools.min(files, 10000);
	}
	private static void adjustBuffers(long maxSizeObserved, int inFiles, PrintStream outstream){
		if(Shared.numBuffers()<3){return;}
		long memNeeded=mergeMemNeeded(maxSizeObserved, inFiles);
		long memAvailable=Shared.memAvailableAdvanced();
		if(memNeeded>memAvailable*0.8){
			outstream.println("Reduced buffers prior to merging due to low memory.");
			Shared.capBufferLen(40);
			Shared.capBuffers(2);
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> mergeRecursive(final ArrayList<String> inList, boolean allowSubprocess){
		assert(maxFiles>1);
		
		ArrayList<String> currentList=inList;
		final int oldZL=ReadWrite.ZIPLEVEL;
		final boolean deleteAfterMerge=(delete && deleteEarly);
		final int limit=Tools.min(mergeFileLimit, maxFiles, maxMergeFiles(maxSizeObserved()));
		while(currentList.size()>maxFiles){
			ReadWrite.ZIPLEVEL=Tools.min(ReadWrite.ZIPLEVEL, 4);
			final int size=currentList.size();
			final int groups=(size+limit-1)/limit;
			assert(groups>0 && groups<size);
			ArrayList<ArrayList<String>> listList=new ArrayList<ArrayList<String>>();
			ArrayList<String> outList=new ArrayList<String>();
			for(int i=0; i<groups; i++){
				listList.add(new ArrayList<String>());
			}
			for(int i=0; i<size; i++){
				listList.get(i%groups).add(currentList.get(i));
			}
			for(ArrayList<String> subList : listList){
				String temp=getTempFile();
				FileFormat ffoutTemp=FileFormat.testOutput(temp, FileFormat.FASTQ, null, true, false, false, false);
				mergeRecursiveInner(subList, ffoutTemp, null, deleteAfterMerge, allowSubprocess);
				outList.add(temp);
			}
			if(delete && !deleteEarly){//delete them now that this pass is finished
				System.err.println("Deleting "+size+" temp files...");
				for(ArrayList<String> subList : listList){
					for(String fname : subList){new File(fname).delete();}
				}
				System.err.println(" done.");
			}
			currentList=outList;
		}
		ReadWrite.ZIPLEVEL=oldZL;
		return currentList;
	}
	
	public void mergeRecursiveInner(ArrayList<String> inList, FileFormat ffoutTemp1, FileFormat ffoutTemp2, boolean deleteAfterMerge, boolean allowSubprocess){
		errorState|=mergeAndDump(inList, /*null, */ffoutTemp1, ffoutTemp2, deleteAfterMerge, useSharedHeader, allowSubprocess, outstream, maxSizeObserved());
	}
	
	private String getTempFile(){
		String temp;
		File dir=new File(".");//(Shared.tmpdir()==null ? null : new File(Shared.tmpdir()));
		if(dir!=null && !dir.exists()){dir.mkdirs();}
		try {
			temp=File.createTempFile("sort_temp_", tempExt, dir).toString();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			KillSwitch.kill(e.getMessage());
			return null;
		}
		return temp;
	}
	
	private boolean mergeAndDump(ArrayList<String> fnames, /*IntList dumpCount, */boolean useHeader, boolean allowSubprocess) {
		final long maxSizeObserved=maxSizeObserved();
		final int limit=Tools.min(mergeFileLimit, maxMergeFiles(maxSizeObserved));
		if(fnames.size()>limit){
			outstream.println("Performing recursive merge to reduce open files.");
			fnames=mergeRecursive(fnames, allowSubprocess);
		}
		return mergeAndDump(fnames, /*dumpCount,*/ ffout1, ffout2, delete, useHeader, allowSubprocess, outstream, maxSizeObserved);
	}
	
	public static boolean mergeAndDump(ArrayList<String> fnames, FileFormat ffout1, FileFormat ffout2, 
			boolean deleteAfterMerge, boolean useHeader, boolean allowSubprocess, PrintStream outstream, final long maxSizeObserved) {
		
		final int oldBuffers=Shared.numBuffers();
		final int oldBufferLen=Shared.bufferLen();
		if(!FileFormat.isFastaFile(fnames.get(0))){
			synchronized(SortByName.class){
				qualityObservedStatic=true;
			}
		}
		
		adjustBuffers(maxSizeObserved, fnames.size(), outstream);
		
		System.err.println("Merging "+fnames);
//		outstream.println("zl="+ReadWrite.ZIPLEVEL);
//		outstream.println("ztd="+ReadWrite.ZIP_THREAD_DIVISOR);
//		outstream.println("mzt="+ReadWrite.MAX_ZIP_THREADS);
//		outstream.println("pigz="+ReadWrite.USE_PIGZ);
		
		ListNum.setDeterministicRandom(false);
		boolean errorState=false;
		final ConcurrentReadOutputStream ros;
		if(ffout1!=null){
			final int buff=1;
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, useHeader);
			ros.start(); //Start the stream
		}else{ros=null;}
		
		ArrayList<CrisContainer> cclist=new ArrayList<CrisContainer>(fnames.size());
		PriorityQueue<CrisContainer> q=new PriorityQueue<CrisContainer>(fnames.size());
		for(int i=0; i<fnames.size(); i++){
			String fname=fnames.get(i);
//			int size=(dumpCount==null ? -1 : dumpCount.get(i));
			CrisContainer cc=new CrisContainer(fname, /*size, */comparator, allowSubprocess);
			if(cc.peek()!=null){
				cclist.add(cc);
				q.add(cc);
			}
		}
		
		//TODO: Use the read counts, stored in dumpCount, to approximately restore random values for shuffling.
		mergeAndDump(q, ros, outstream);
		if(verbose){
			outstream.println("Finished processing "+fnames);
		}
		
		for(CrisContainer cc : cclist){
			errorState|=cc.close();
		}
		if(deleteAfterMerge){
			for(String fname : fnames){
				new File(fname).delete();
			}
		}
		if(ros!=null){errorState|=ReadWrite.closeStream(ros);}
		
		Shared.setBufferLen(oldBufferLen);
		Shared.setBuffers(oldBuffers);
		
		return errorState;
	}
	
	private static void mergeAndDump(final PriorityQueue<CrisContainer> q, final ConcurrentReadOutputStream ros, PrintStream outstream) {
		
//		for(CrisContainer cc : q){
//			assert(!cc.cris().paired()) : FASTQ.TEST_INTERLEAVED+", "+FASTQ.FORCE_INTERLEAVED+
//				", "+cc.cris().fname()+", "+cc.cris().getClass();
//		}
		
		long maxLen=0;
		
		final int limit=100000;
		ArrayList<Read> buffer=new ArrayList<Read>(2*limit);
		while(!q.isEmpty()){
			if(verbose2){outstream.println("q size: "+q.size());}
			for(int i=0; !q.isEmpty() && (buffer.size()<limit || i<1); i++){//For loop to force it to run at least once
				CrisContainer cc=q.poll();
				if(verbose2){outstream.println("Polled a cc.");}
				ArrayList<Read> list=cc.fetch();
				if(verbose2){outstream.println("Grabbed "+list.size()+" reads.");}
				buffer.addAll(list);
				if(verbose2){outstream.println("Buffer size: "+buffer.size());}
				if(cc.hasMore()){
					if(verbose2){outstream.println("Returned cc.");}
					q.add(cc);
				}else{
					if(verbose2){outstream.println("Discarded cc.");}
				}
			}
			Shared.sort(buffer, comparator);
			if(verbose2){outstream.println("Sorted buffer.");}
//			if(repair){
//				assert(false) : "TODO";
//			}
			
			for(int i=1; i<buffer.size(); i++){
				Read a=buffer.get(i-1);
				Read b=buffer.get(i);
				assert(comparator.compare(a, b)<=0) : "\n"+a+"\n"+b;
//				assert(a.mate==null) : "\n"+a+"\n"+b;
			}
			
			final Read peek=(q.isEmpty() ? null : q.peek().peek());
			final int maxIndex=(peek==null ? buffer.size() : indexOfLowestAbovePivot(buffer, peek));
			ArrayList<Read> list=new ArrayList<Read>(maxIndex);
			
			for(int index=0; index<maxIndex; index++){
				Read r=buffer.get(index);
				assert(peek==null || comparator.compare(peek, r)>0) : "\n"+peek+"\n"+r;
				list.add(r);
				maxLen=Tools.max(maxLen, r.length());
			}
			if(ros!=null){ros.add(list, 0);}
			
			ArrayList<Read> oldbuffer=buffer;
			buffer=new ArrayList<Read>(2*limit);
			for(int i=maxIndex, size=oldbuffer.size(); i<size; i++){
				buffer.add(oldbuffer.get(i));
			}
			if(verbose2){outstream.println("Buffer contains "+buffer.size()+" reads.");}
		}
		
		assert(buffer.isEmpty());
		synchronized(SortByName.class){
			maxLengthObservedStatic=Tools.max(maxLengthObservedStatic, maxLen);
		}
	}

//
//	final Read peek=(q.isEmpty() ? null : q.peek().peek());
//
//	final int maxIndex=(peek==null ? buffer.size() : indexOfLowestAbovePivot(buffer, peek));
//
//	for(int index=0; index<maxIndex; index++){
//		Read r=buffer.get(index);
//		list.add(r);
//		if(list.size()>=200){
//			if(ros!=null){ros.add(list, 0);}
//			list=new ArrayList<Read>(200);
//		}
//	}
//
//	if(!list.isEmpty()){
//		if(ros!=null){ros.add(list, 0);}
//		list=null;
//	}
////	if(verbose2){outstream.println("Stopped at index "+index+".");}
//
//	ArrayList<Read> oldbuffer=buffer;
//	buffer=new ArrayList<Read>(2*limit);
//	for(int i=maxIndex, size=oldbuffer.size(); i<size; i++){
//		buffer.add(oldbuffer.get(i));
//	}
//	if(verbose2){outstream.println("Buffer contains "+buffer.size()+" reads.");}
//}
	
	private static final int indexOfLowestAbovePivot(final ArrayList<Read> list, final Read pivot){
		final int size=list.size();
		final int maxIndex=binarySearch(list, pivot);
		if(maxIndex<0){return 0;}
		if(maxIndex>=size){return size;}
		Read r=list.get(maxIndex);
		final int x=comparator.compare(pivot, r);
		assert(x<=0) : x+"\n"+pivot.id+"\n"+r.id;
		final int ret=(x==0 ? maxIndex+1 : maxIndex);
		assert(ret==size || comparator.compare(pivot, list.get(ret))<0);
		assert(ret==0 || comparator.compare(pivot, list.get(ret-1))>=0);
		return ret;
	}
	
	/** Returns highest index of an equal value, or the insert position */
	private static final int binarySearch(final ArrayList<Read> list, final Read pivot){
		int a=0, b=list.size()-1;
		while(b>a){
			final int mid=(a+b)/2;
			final Read r=list.get(mid);
			final int x=comparator.compare(pivot, r);
			if(x<0){b=mid;}
			else if(x>0){a=mid+1;}
			else{return mid;}
		}
		
//		for(int i=1; i<list.size(); i++){
//			assert(comparator.compare(list.get(i-1), list.get(i))<=0) : i+"\n"+list.get(i-1)+"\n"+ list.get(i)+"\n";
//		}
//
//		assert(a>=list.size() || comparator.compare(pivot, list.get(a))<=0) :
//			comparator.compare(pivot, list.get(a))+"\n"+
//					a+", "+b+", "+list.size()+"\n"+pivot.id+"\n"+list.get(a).id;
//			comparator.compare(pivot, list.get(a))+", "+a+", "+b+", "+list.size()+"\n"+pivot.id+"\n"+list.get(a-1).id+"\n"+list.get(a).id+"\n"+list.get(a+1).id;
		while(a>=0 && a<list.size()-1){
			final int x=comparator.compare(pivot, list.get(a+1));
			assert(x<=0) : x+", "+a+", "+b+", "+list.size()+"\n"+pivot.id+"\n"+list.get(a).id;
			if(x<0){break;}
			a++;
		}
		
		if(a>=0 && a==list.size()-1){//Corner case; insert after last element was not handled
			if(comparator.compare(pivot, list.get(a))>0){a++;}
		}
		
		assert(a>=list.size() || comparator.compare(pivot, list.get(a))<=0) :
			comparator.compare(pivot, list.get(a))+"\n"+
					a+", "+b+", "+list.size()+"\n"+pivot.id+"\n"+list.get(a).id;
		
		return a;
	}
	
	private void sortAndDump(final ArrayList<Read> storage, final long currentMem, final AtomicLong outstandingMem, String fname1, String fname2, boolean useHeader) {
		String temp1=fname1;
		String temp2=fname2;
		if(temp1==null){
			synchronized(outTemp){
				if(verbose2){outstream.println("Synced to outTemp to make a temp file.");}
				temp1=getTempFile();
				if(verbose2){outstream.println("Temp file: "+temp1);}
				outTemp.add(temp1);
			}
		}

		if(verbose || true){outstream.println("Created a WriteThread for "+temp1+(temp2==null ? "" : ", "+temp2));}
		WriteThread wt=new WriteThread(storage, currentMem, outstandingMem, temp1, temp2, useHeader, outstream);
		wt.start();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static class WriteThread extends Thread{
		
		public WriteThread(final ArrayList<Read> storage_, final long currentMem_, final AtomicLong outstandingMem_, 
				String fname1_, String fname2_, boolean useHeader_, PrintStream outstream_){
			storage=storage_;
			currentMem=currentMem_;
			outstandingMem=outstandingMem_;
			fname1=fname1_;
			fname2=fname2_;
			useHeader=useHeader_;
			outstream=outstream_;
		}
		
		@Override
		public void run(){
			
			if(verbose){outstream.println("Started a WriteThread.");}
			final FileFormat ffout1=FileFormat.testOutput(fname1, FileFormat.FASTQ, null, true, false, false, false);
			final FileFormat ffout2=FileFormat.testOutput(fname2, FileFormat.FASTQ, null, true, false, false, false);
			final ConcurrentReadOutputStream ros;
			if(ffout1!=null){
				final int buff=4;
				ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, useHeader);
				ros.start(); //Start the stream
			}else{ros=null;}

			if(verbose){outstream.println("Started a ros.");}
			
			
			Shared.sort(storage, comparator);
			
			if(verbose){outstream.println("Sorted reads.");}
			
			ArrayList<Read> buffer=new ArrayList<Read>(200);
			long id=0;
			for(int i=0, lim=storage.size(); i<lim; i++){
				Read r=storage.set(i, null);
				buffer.add(r);
				if(buffer.size()>=200){
					if(ros!=null){ros.add(buffer, id);}
					id++;
					buffer=new ArrayList<Read>(200);
				}
			}
			if(ros!=null && buffer.size()>0){ros.add(buffer, id);}
			errorState|=ReadWrite.closeStream(ros);
			if(verbose){outstream.println("Closed ros.");}
			
			synchronized(outstandingMem){
				outstandingMem.addAndGet(-currentMem);
				if(verbose){outstream.println("Decremented outstandingMem: "+outstandingMem);}
				outstandingMem.notify();
				if(verbose){outstream.println("Notified outstandingMem.");}
			}
		}
		
		final ArrayList<Read> storage;
		final long currentMem;
		final AtomicLong outstandingMem;
		final String fname1;
		final String fname2;
		final boolean useHeader;
		boolean errorState=false;
		final PrintStream outstream;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;
	
	private ArrayList<String> outTemp=new ArrayList<String>();
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	private String tempExt=null;

	private String giTableFile=null;;
	private String taxTreeFile=null;
	private String accessionFile=null;
	
	/*--------------------------------------------------------------*/

	long maxLengthObserved=0;
	private static long maxLengthObservedStatic=0;
	private final boolean qualityObserved;
	private static boolean qualityObservedStatic;
	
	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;

	private boolean delete=true;
	private boolean deleteEarly=false;
	
	private boolean useSharedHeader=false;
	
	private boolean allowTempFiles=true;
	
	private int minlen=0;

	private float memBlockMult=0.30f; //Crashed once at old default 0.35
	private float memTotalMult=0.65f; //Crashed once at old default 0.75
	
	/** Max files to merge per pass */
	private int maxFiles=12;
	private boolean setMaxFiles=false;
	
	private static int mergeFileLimit=24;
	
	private boolean genKmer=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	
	private final boolean clump;
	
	private static ReadComparator comparator=ReadComparatorName.comparator;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** Print verbose messages */
	public static final boolean verbose2=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
