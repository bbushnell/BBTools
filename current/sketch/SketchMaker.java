package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicInteger;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import structures.LongList;
import tax.AccessionToTaxid;
import tax.GiToTaxid;
import tax.ImgRecord2;
import tax.TaxNode;
import tax.TaxTree;
import tracker.ReadStats;

/**
 * Creates MinHashSketches rapidly.
 * 
 * @author Brian Bushnell
 * @date July 6, 2016
 *
 */
public class SketchMaker extends SketchObject {
	
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
		
		final int mode=parseMode(args);
		if(mode==PER_FILE/* || mode==ONE_SKETCH || mode==PER_HEADER*/ || mode==PER_SEQUENCE){//ONE_SKETCH does not work for multiple input files.
			recallCompareSketch(args);
			return;
		}
		
		final int oldBufLen=Shared.bufferLen();
		
		//Create an instance of this class
		SketchMaker x=new SketchMaker(args);
		
		//Run the object
		x.process(t);
		
		Shared.setBufferLen(oldBufLen);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	private static void recallCompareSketch(String[] args){
		ArrayList<String> list=new ArrayList<String>(args.length+1);
		for(int i=0; i<args.length; i++){
			if(args[i].startsWith("out=")){
				args[i]=args[i].replaceFirst("out=", "outsketch=");
			}
			list.add(args[i]);
		}
		list.add("sketchonly");
		CompareSketch.main(list.toArray(new String[0]));
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SketchMaker(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		//Create a parser object
		Parser parser=new Parser();
		
		int minSizeKmers_=100;
		int files_=1;
		int mode_=ONE_SKETCH;
		hashNames=true;
		boolean setPrefilter=false;
		defaultParams.printDepth=defaultParams.printDepth2=defaultParams.printVolume=false;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("files")){
				files_=Integer.parseInt(b);
			}else if(a.equals("minsize")){
				minSizeKmers_=Parse.parseIntKMG(b);
			}else if(a.equals("prefilter")){
				prefilter=Parse.parseBoolean(b);
				setPrefilter=true;
			}
			
			else if(a.equals("name") || a.equals("taxname")){
				outTaxName=b;
			}else if(a.equals("name0")){
				outName0=b;
			}else if(a.equals("fname")){
				outFname=b;
			}else if(a.equals("taxid") || a.equals("tid")){
				outTaxID=Integer.parseInt(b);
			}else if(a.equals("spid")){
				outSpid=Integer.parseInt(b);
			}else if(a.equals("imgid")){
				outImgID=Integer.parseInt(b);
			}else if((a.startsWith("meta_") || a.startsWith("mt_")) && b!=null){
				if(outMeta==null){outMeta=new ArrayList<String>();}
				int underscore=a.indexOf('_', 0);
				outMeta.add(a.substring(underscore+1)+":"+b);
			}else if(a.equals("parsesubunit")){
				parseSubunit=Parse.parseBoolean(b);
			}
			
			else if(parseMode(arg, a, b)>-1){
				mode_=parseMode(arg, a, b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}
			
			else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
			}else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
			}else if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equalsIgnoreCase("img") || a.equals("imgfile") || a.equals("imgdump")){
				imgFile=b;
			}
			
			else if(a.equals("tossjunk")){
				tossJunk=Parse.parseBoolean(b);
			}
			
//			else if(a.equals("silva")){//Handled by parser
//				TaxTree.SILVA_MODE=Parse.parseBoolean(b);
//			}
			
			else if(a.equals("taxlevel") || a.equals("tl") || a.equals("level") || a.equals("lv")){
				taxLevel=TaxTree.parseLevel(b);
			}
			
			else if(parseSketchFlags(arg, a, b)){
				//do nothing
//				System.err.println("a: "+arg);
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
//				System.err.println("b: "+arg);
			}else if(defaultParams.parse(arg, a, b)){
				//do nothing
//				System.err.println("c: "+arg);
//				System.err.println(defaultParams.printDepth);
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		if("auto".equalsIgnoreCase(imgFile)){imgFile=TaxTree.defaultImgFile();}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		
		outMeta=SketchObject.fixMeta(outMeta);
		
		postParse();
		minSizeKmers=minSizeKmers_;
		mode=mode_;
		
		minSizeBases=minSizeKmers_+k-1;
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			
			extin=parser.extin;
		}
		files=(out1==null ? 0 : files_);

		if(!setPrefilter && !prefilter && (mode==PER_TAXA || mode==PER_IMG) && in1!=null && !in1.startsWith("stdin") && (AUTOSIZE || AUTOSIZE_LINEAR || targetSketchSize>200)){
			prefilter=true;
			System.err.println("Enabled prefilter due to running in per-"+(mode==PER_TAXA ? "taxa" : "IMG")+" mode; override with 'prefilter=f'.");
		}
		
		assert(mode!=ONE_SKETCH || files<2) : "Multiple output files are not allowed in single-sketch mode.";
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
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
		
		ffout=makeFFArray(out1, files, overwrite, append);
		if(ffout==null || ffout.length<1){
			System.err.println("WARNING: No output files were specified; no sketches will be written.\n");
		}
		
//		//Ensure output files can be written
//		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
//			outstream.println((out1==null)+", "+out1);
//			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
//		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2, taxTreeFile, giTableFile, imgFile, SSUMap.r16SFile, SSUMap.r18SFile)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, taxTreeFile, giTableFile, imgFile, SSUMap.r16SFile, SSUMap.r18SFile)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
//		assert(false) :  defaultParams.trackCounts();
		tool=new SketchTool(targetSketchSize, defaultParams);
		
		if(taxTreeFile!=null){setTaxtree(taxTreeFile, outstream);}
		
		if(giTableFile!=null){
			loadGiToTaxid();
		}
		if(accessionFile!=null){
			AccessionToTaxid.tree=taxtree;
			outstream.println("Loading accession table.");
			AccessionToTaxid.load(accessionFile);
			System.gc();
		}
		if(imgFile!=null){
			TaxTree.loadIMG(imgFile, true, outstream);
		}
		SSUMap.load(outstream);
		
		if(prefilter){
			if(mode==PER_TAXA){sizeList=sizeList(); sizeMap=null;}
			else if(mode==PER_IMG){sizeMap=sizeMap(); sizeList=null;}
			else{
				assert(false) : "Wrong mode for prefilter; should be taxa or img.";
				sizeList=null; sizeMap=null;
			}
		}else{
			sizeList=null; sizeMap=null;
		}
	}
	
	private static FileFormat[] makeFFArray(String fname0, int files, boolean overwrite, boolean append){
		if(files<1 || fname0==null){return null;}
		String[] fnames=new String[files];
		FileFormat[] ff=new FileFormat[files];
		for(int i=0; i<files; i++){
			String fname=fname0;
			if(files>1){
				assert(fname.indexOf('#')>-1) : "Output name requires # symbol for multiple files.";
				fname=fname.replaceFirst("#", ""+i);
			}
			fnames[i]=fname;
			ff[i]=FileFormat.testOutput(fname, FileFormat.TEXT, null, true, overwrite, append, false);
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, fnames)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+Arrays.toString(fnames)+"\n");
		}
		
		return ff;
	}
	
//	private static ByteStreamWriter[] makeTSWArray(FileFormat[] ff){
//		if(ff==null || ff.length==0){return null;}
//		ByteStreamWriter[] tsw=new ByteStreamWriter[ff.length];
//		for(int i=0; i<ff.length; i++){
//			tsw[i]=new ByteStreamWriter(ff[i]);
//			tsw[i].start();
//		}
//		return tsw;
//	}
	
	private static ByteStreamWriter[] makeTSWArray(FileFormat[] ff){
		if(ff==null || ff.length==0){return null;}
		ByteStreamWriter[] tsw=new ByteStreamWriter[ff.length];
		for(int i=0; i<ff.length; i++){
			tsw[i]=new ByteStreamWriter(ff[i]);
			tsw[i].start();
		}
		return tsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	//Makes a list of genome sizes (bases, not kmers) per taxa.
	private LongList sizeList(){
		Timer t=new Timer();
		Shared.GC_BEFORE_PRINT_MEMORY=true;
		t.start("Making taxa prefilter.");
		//For some reason this function is very slow, using only ~95% java and ~95% bgzip CPU.
		//But that's about the same speed as reformat...  BBDuk is only a little faster.
		
//		assert(false) : ReadWrite.USE_PIGZ+", "+ReadWrite.USE_UNPIGZ+", "+ReadWrite.MAX_ZIP_THREADS+", "+Shared.threads();
		
		LongList sizes=new LongList();
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
//			if(defaultParams.samplerate!=1){cris.setSampleRate(defaultParams.samplerate, sampleseed);} //Why would I ever want this here?
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();
		//Grab the actual read list from the ListNum
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		//As long as there is a nonempty read list...
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			//					if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);

				int taxID=-1;
				TaxNode tn=null;
				if(taxtree!=null){
					tn=taxtree.parseNodeFromHeader(r1.id, bestEffort);
//					System.err.println("A: "+bestEffort+", "+tn);//123
					while(tn!=null && tn.pid!=tn.id && tn.level<taxLevel){
						TaxNode temp=taxtree.getNode(tn.pid);
						if(temp==null || temp.level>=TaxTree.LIFE){break;}
						tn=temp;
					}
					if(tn!=null){taxID=tn.id;}
				}
				
				if(taxID>0){
					long a=r1.length();
					long b=r1.mateLength();
					if(a<k){a=0;}
					if(b<k){b=0;}
					sizes.increment(taxID, a+b);
				}
			}

			//Notify the input stream that the list was used
			cris.returnList(ln);
			//					if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

			//Fetch a new list
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		
		errorState|=ReadWrite.closeStream(cris);
//		outstream.println("Created prefilter.");
		t.stop("Created prefilter:");
		Shared.printMemory();
		System.err.println();
		
		return sizes;
	}
	
	//Makes a list of genome sizes (bases, not kmers) per img.
	private HashMap<Long, Long> sizeMap(){
		Timer t=new Timer();
		t.start("Making img prefilter.");
		
//		assert(false) : ReadWrite.USE_PIGZ+", "+ReadWrite.USE_UNPIGZ+", "+ReadWrite.MAX_ZIP_THREADS+", "+Shared.threads();
		
		HashMap<Long, Long> sizes=new HashMap<Long, Long>();
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			if(defaultParams.samplerate!=1){cris.setSampleRate(defaultParams.samplerate, sampleseed);}
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();
		//Grab the actual read list from the ListNum
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		//As long as there is a nonempty read list...
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			//					if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);

				long imgID=ImgRecord2.parseImgId(r1.id, true);
				assert(imgID>-1) : "IMG records must start with IMG number followed by a space: "+r1.id;
				
				if(imgID>0){
					long a=r1.length();
					long b=r1.mateLength();
					if(a<k){a=0;}
					if(b<k){b=0;}
					
					if(a+b>0){
						Long old=sizes.get(imgID);
						if(old==null){sizes.put(imgID, a+b);}
						else{sizes.put(imgID, a+b+old);}
					}
				}
			}

			//Notify the input stream that the list was used
			cris.returnList(ln);
			//					if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

			//Fetch a new list
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		
		errorState|=ReadWrite.closeStream(cris);
//		outstream.println("Created prefilter.");
		t.stop("Created prefilter:");
		Shared.printMemory();
		System.err.println();
		
		return sizes;
	}

	/** Create read streams and process all data */
	void process(Timer t){

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		if(mode==ONE_SKETCH && !forceDisableMultithreadedFastq && Shared.threads()>2 && ffin1.fastq()){
//			Shared.setBufferLen(2);
			singleSketchMT();
		}else{
			final int oldLen=Shared.bufferLen();
			Shared.capBufferLen(ffin1.fastq() ? 40 : 4);

			//Turn off read validation in the input threads to increase speed
			final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
			Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;

			//Create a read input stream
			final ConcurrentReadInputStream cris;
			{
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
				cris.start(); //Start the stream
				if(verbose){outstream.println("Started cris");}
			}

			//Process the reads in separate threads
			spawnThreads(cris);

			if(verbose){outstream.println("Finished; closing streams.");}

			//Write anything that was accumulated by ReadStats
			errorState|=ReadStats.writeAll();
			//Close the read streams
			errorState|=ReadWrite.closeStream(cris);

			//TODO: Write sketch

			//Reset read validation
			Read.VALIDATE_IN_CONSTRUCTOR=vic;
			Shared.setBufferLen(oldLen);
		}
		
		//Report timing and results
		t.stop();
		outstream.println("Wrote "+sketchesWritten+" of "+sketchesMade+" sketches.\n");
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private void singleSketchMT(){
		Timer t=new Timer();
		Sketch sketch=tool.processReadsMT(ffin1, ffin2, Shared.threads(), 
				maxReads, mode, defaultParams.samplerate, defaultParams.minEntropy, defaultParams.minProb, defaultParams.minQual, false);
		
		if(outTaxID>=0){sketch.taxID=outTaxID;}
		if(outTaxName!=null){sketch.setTaxName(outTaxName);}
		if(outFname!=null){sketch.setFname(outFname);}
		if(outName0!=null){sketch.setName0(outName0);}
		if(outSpid>=0){sketch.spid=outSpid;}
		if(outImgID>=0){sketch.imgID=outImgID;}
		sketch.setMeta(outMeta);
		
		//Accumulate per-thread statistics
		readsProcessed+=sketch.genomeSequences;
		basesProcessed+=sketch.genomeSizeBases;
		kmersProcessed+=sketch.genomeSizeKmers;

		sketchesMade+=1;
		
		t.stop("Finished sketching: ");
		Shared.printMemory();

		if(ffout!=null && ffout.length>0){
			sketch.addSSU();
			SketchTool.write(sketch, ffout[0]);
			sketchesWritten+=1;
		}
	}
	
	/** Spawn process threads */
	@SuppressWarnings("unchecked")
	private void spawnThreads(final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		Timer t=new Timer();
		
		//Determine how many threads may be used
		final int threads=Tools.mid(1, Shared.threads(), 14);//Probably capped for memory reasons.  Rarely hits 1200% anyway (was 12, bumped to 14).
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		
		if(mode==PER_TAXA || mode==PER_IMG){
			longMaps=new HashMap[MAP_WAYS];
			for(int i=0; i<longMaps.length; i++){
				longMaps[i]=new HashMap<Long, SketchHeap>();
			}
		}
		
		if(mode!=ONE_SKETCH){tsw=makeTSWArray(ffout);}
		
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		SketchHeap singleHeap=null;
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
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			kmersProcessed+=pt.smm.kmersProcessed;
			
			sketchesMade+=pt.sketchesMadeT;
			sketchesWritten+=pt.sketchesWrittenT;
			
//			System.err.println(pt.sketchesMadeT+" "+pt.sketchesWrittenT);
			
//			System.err.println("pt.readsProcessedT="+pt.readsProcessedT);
			if(mode==ONE_SKETCH){
				SketchHeap temp=pt.smm.heap;
				
				if(temp==null){
					//do nothing
				}else if(singleHeap==null){singleHeap=pt.smm.heap;}
				else{singleHeap.add(pt.smm.heap);}
				
				if(singleHeap!=null){
					if(outTaxID>=0){singleHeap.taxID=outTaxID;}
					if(outTaxName!=null){singleHeap.setTaxName(outTaxName);}
					if(outFname!=null){singleHeap.setFname(outFname);}
					if(outName0!=null){singleHeap.setName0(outName0);}
					if(outImgID>=0){singleHeap.imgID=outImgID;}
					
					singleHeap.genomeSizeBases=basesProcessed;
					singleHeap.genomeSizeKmers=kmersProcessed;
					singleHeap.genomeSequences=readsProcessed;
				}
			}
			success&=pt.success;
		}
		
		if(singleHeap!=null){
			singleHeap.setFname(ffin1.simpleName());
			if(singleHeap.name0()==null){
				singleHeap.setName0(ffin1.simpleName());
			}
		}
		
		t.stop("Finished sketching: ");
		Shared.printMemory();
		
		if(ffout!=null){
			if(mode==PER_TAXA || mode==PER_IMG){
				if(tsw==null){tsw=makeTSWArray(ffout);}
				success&=writeMap(longMaps);
			}else if(mode==ONE_SKETCH){
				Sketch sketch=new Sketch(singleHeap, true, tool.trackCounts, outMeta);
				if(outTaxID>=0){sketch.taxID=outTaxID;}
				if(outTaxName!=null){sketch.setTaxName(outTaxName);}
				if(outFname!=null){sketch.setFname(outFname);}
				if(outName0!=null){sketch.setName0(outName0);}
				if(outSpid>=0){sketch.spid=outSpid;}
				if(outImgID>=0){sketch.imgID=outImgID;}
				if(ffout!=null && ffout.length>0){
					sketch.addSSU();
					SketchTool.write(sketch, ffout[0]);
				}
				sketchesMade++;
				sketchesWritten++;
			}
		}
		
		if(tsw!=null){
			for(int i=0; i<tsw.length; i++){tsw[i].poisonAndWait();}
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private boolean writeMap(HashMap<Long, SketchHeap>[] maps){
		
		//Determine how many threads may be used
		final int threads=files;

		//Fill a list with WriteThreads
		ArrayList<WriteThread> alwt=new ArrayList<WriteThread>(threads);
		
		@SuppressWarnings("unchecked")
		ArrayDeque<SketchHeap>[] heaps=new ArrayDeque[threads];
		for(int i=0; i<threads; i++){
			heaps[i]=new ArrayDeque<SketchHeap>();
			WriteThread wt=new WriteThread(i, heaps[i]);
			alwt.add(wt);
		}
		
		for(int i=0; i<maps.length; i++){
			HashMap<Long, SketchHeap> longMap=maps[i];
			for(Entry<Long, SketchHeap> entry : longMap.entrySet()){
				//			set.remove(entry);  This will probably not work
				SketchHeap entryHeap=entry.getValue();
				sketchesMade++;
				if(entryHeap.size()>0 && entryHeap.genomeSizeKmers>=minSizeKmers){
					heaps[(entry.hashCode()&Integer.MAX_VALUE)%threads].add(entryHeap);
				}
			}
			//		intMap.clear();
			maps[i]=null;
		}

		//Start the threads
		for(WriteThread wt : alwt){wt.start();}
		
		//Wait for completion of all threads
		boolean success=true;
		for(WriteThread wt : alwt){

			//Wait until this thread has terminated
			while(wt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					wt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
//			sketchesMade+=wt.sketchesMadeT;
			sketchesWritten+=wt.sketchesWrittenT;
			success&=wt.success;
		}
		return success;
	}
	
	private class WriteThread extends Thread{
		
		WriteThread(int tnum_, ArrayDeque<SketchHeap> queue_){
			tnum=tnum_;
			queue=queue_;
		}
		
		@Override
		public void run(){
			success=false;
			for(SketchHeap polledHeap=queue.poll(); polledHeap!=null; polledHeap=queue.poll()){
				if(polledHeap.sketchSizeEstimate()>0){
					Sketch sketch=new Sketch(polledHeap, true, tool.trackCounts, outMeta);
					if(outTaxID>=0 && sketch.taxID<0){sketch.taxID=outTaxID;}
					if(outTaxName!=null && sketch.taxName()==null){sketch.setTaxName(outTaxName);}
					if(outFname!=null && sketch.fname()==null){sketch.setFname(outFname);}
					if(outName0!=null && sketch.name0()==null){sketch.setName0(outName0);}
					if(outSpid>=0 && sketch.spid<0){sketch.spid=outSpid;}
					if(outImgID>=0 && sketch.imgID<0){sketch.imgID=outImgID;}
					sketch.addSSU();
					SketchTool.write(sketch, tsw[tnum], bb);
					sketchesWrittenT++;
				}
			}
			bb=null;
			success=true;
			queue=null;
		}
		
		ArrayDeque<SketchHeap> queue;
		final int tnum;
		private ByteBuilder bb=new ByteBuilder();
//		long sketchesMadeT=0;
		long sketchesWrittenT=0;
		boolean success=false;
	}
	
//	private void writeOutput(ConcurrentHashMap<Integer, SketchHeap> map){
//		ByteStreamWriter tsw=new ByteStreamWriter(ffout);
//		tsw.start();
//		KeySetView<Integer, SketchHeap> y=map.keySet();
//		for(Integer x : map.keySet()){
//			SketchHeap smm.heap=map.get(x);
//			Sketch s=tool.toSketch(smm.heap);
//			tool.write(s, tsw);
//		}
//		tsw.poisonAndWait();
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Tax Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private void loadGiToTaxid(){
		Timer t=new Timer();
		outstream.println("Loading gi to taxa translation table.");
		GiToTaxid.initialize(giTableFile);
		t.stop();
		if(true){
			outstream.println("Time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final int tid_){
			cris=cris_;
			threadID=tid_;
			
			smm=new SketchMakerMini(tool, mode, defaultParams.minEntropy, defaultParams.minProb, defaultParams.minQual);
			localMap=(mode==PER_TAXA ? new HashMap<Long, SketchHeap>() : null);
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			bb=null;
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}
			
//			long cntr1=0, cntr2=0, cntr3=0, cntr4=0;

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;

					processReadPair(r1, r2);
				}

				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
			
			if(mode==PER_TAXA && localMap.size()>0){dumpLocalMap();}
			
//			System.out.println(cntr1+", "+cntr2+", "+cntr3+", "+cntr4);
		}
		
		void processReadPair(Read r1, Read r2){
			
			if(mode==PER_TAXA){
				assert(smm.heap==null || smm.heap.size()==0) : smm.heap.genomeSizeBases+", "+smm.heap;
				assert(smm.heap==null || smm.heap.genomeSizeBases==0) : smm.heap.genomeSizeBases+", "+smm.heap;
			}else if(mode==PER_SEQUENCE || mode==PER_IMG){
				assert(smm.heap==null || smm.heap.size()==0) : smm.heap.genomeSizeBases+", "+smm.heap;
				assert(smm.heap==null || smm.heap.genomeSizeBases==0) : smm.heap.genomeSizeBases+", "+smm.heap;
			}else{
				assert(smm.heap!=null && smm.heap.capacity()>=targetSketchSize) : targetSketchSize+", "+(smm.heap==null ? "null" : ""+smm.heap.capacity());
			}
			
			//Track the initial length for statistics
			final int initialLength1=r1.length();
			final int initialLength2=r1.mateLength();
			final String rid=r1.id;
			
			//Increment counters
			readsProcessedT+=r1.pairCount();
			basesProcessedT+=initialLength1+initialLength2;

			if(initialLength1<k && initialLength2<k){return;}
			final long expectedBases;
			final long imgID;
			{
				int taxID=-1;
				TaxNode tn=null;
				if(taxtree!=null && (mode==PER_TAXA || mode==PER_IMG || mode==PER_SEQUENCE || ((mode==ONE_SKETCH/* || mode==PER_HEADER*/) && smm.heap.taxID<0))){
					if(mode==PER_IMG){
						imgID=ImgRecord2.parseImgId(rid, true);
						tn=taxtree.imgToTaxNode(imgID);
						if(tn==null){tn=taxtree.parseNodeFromHeader(rid, bestEffort);}
//						assert(tn!=null || !rid.startsWith("tid")) : imgID+", "+taxID+", "+rid; //123
					}else{
						imgID=ImgRecord2.parseImgId(rid, false);
						tn=taxtree.parseNodeFromHeader(rid, bestEffort);
//						assert(tn!=null || !rid.startsWith("tid")) : imgID+", "+taxID+", "+rid; //123
					}
//					assert(false) : imgID +", "+rid;
					//					System.err.println("B: "+bestEffort+", "+tn);//123
					while(tn!=null && tn.pid!=tn.id && tn.level<taxLevel){
						TaxNode temp=taxtree.getNode(tn.pid);
						if(temp==null || temp.level>=TaxTree.LIFE){break;}
						tn=temp;
					}
//					assert(tn!=null) : imgID+", "+taxID+", "+rid; //123
					if(tn!=null){taxID=tn.id;}
					//					System.err.println("Node: "+rid+"\n->\n"+tn);
					
//					assert(taxID>0) : imgID+", "+taxID+", "+rid; //123
				}else{
					imgID=-1;
				}
				
				final long unitSizeBases;
				if(sizeList!=null){
					unitSizeBases=taxID<0 ? -1 : sizeList.get(taxID);
				}else if(sizeMap!=null){
					unitSizeBases=sizeMap.get(imgID);
				}else{
					unitSizeBases=-1;
				}
				
				
				if(mode==PER_TAXA){
					if(tossJunk && tn==null){return;}
					if(tn!=null){
						if(taxID==0 || (tn.level>taxLevel && tn.level>=TaxTree.PHYLUM)){return;}
						TaxNode parent=taxtree.getNode(tn.pid);
						if(parent.pid==parent.id){return;}
						if(prefilter && unitSizeBases>=0 && unitSizeBases<minSizeBases){return;}
						if(tossJunk){
							if(parent.level==TaxTree.LIFE){return;}
							if(tn.level==TaxTree.NO_RANK){return;}
						}
					}
				}
				
				if(mode==PER_SEQUENCE){
					expectedBases=initialLength1+initialLength2;
					if(expectedBases<minSizeBases){return;}
					int expectedSketchSize=toSketchSize(expectedBases, -1, -1, targetSketchSize);
					if(expectedSketchSize<minSketchSize){return;}
					if(smm.heap==null || smm.heap.capacity()<expectedSketchSize){
						smm.heap=new SketchHeap(expectedSketchSize, defaultParams.minKeyOccuranceCount, defaultParams.trackCounts());
					}
				}else if(mode==PER_TAXA){
					if(taxID>=0 && localMap.containsKey((long)taxID)){
						smm.heap=localMap.get((long)taxID);
						assert(smm.heap.taxID==taxID);
					}else if(sizeList==null){
						if(smm.heap==null){smm.heap=new SketchHeap(targetSketchSize, defaultParams.minKeyOccuranceCount, defaultParams.trackCounts());}
					}else{
						expectedBases=unitSizeBases>-1 ? unitSizeBases : initialLength1+initialLength2;
						if(expectedBases<minSizeBases){return;}
						int expectedSketchSize=toSketchSize(expectedBases, -1, -1, targetSketchSize);
						if(expectedSketchSize<minSketchSize){return;}
						if(smm.heap==null || smm.heap.capacity()!=expectedSketchSize){
							smm.heap=new SketchHeap(expectedSketchSize, defaultParams.minKeyOccuranceCount, defaultParams.trackCounts());
						}
//						assert(taxID!=96897) : taxID+", "+expectedBases+", "+expectedSketchSize+", "+unitSizeBases+", "+initialLength1+", "+sizeList.get(taxID);
					}
					assert(!localMap.containsKey((long)taxID) || localMap.get((long)taxID)==smm.heap) : taxID;
				}else if(mode==PER_IMG){
					if(sizeMap==null){
						if(smm.heap==null){smm.heap=new SketchHeap(targetSketchSize, defaultParams.minKeyOccuranceCount, defaultParams.trackCounts());}
					}else{
						expectedBases=unitSizeBases>-1 ? unitSizeBases : initialLength1+initialLength2;
						if(expectedBases<minSizeBases){return;}
						int expectedSketchSize=toSketchSize(expectedBases, -1, -1, targetSketchSize);
						if(expectedSketchSize<minSketchSize){return;}
						if(smm.heap==null || smm.heap.capacity()!=expectedSketchSize){
							smm.heap=new SketchHeap(expectedSketchSize, defaultParams.minKeyOccuranceCount, defaultParams.trackCounts());
						}
					}
				}
				
				assert(smm.heap!=null) : mode+", "+(sizeList==null);
				assert(taxID<0 || smm.heap.taxID<0 || smm.heap.taxID==taxID); //This is important.
				assert(imgID<0 || smm.heap.imgID<0 || smm.heap.imgID==imgID);

				if(smm.heap.taxID<0){smm.heap.taxID=taxID;}
				if(smm.heap.imgID<0){smm.heap.imgID=imgID;}
				if(smm.heap.name0()==null){
					smm.heap.setName0(rid);
				}
				if(smm.heap.taxName()==null && tn!=null){
					smm.heap.setTaxName(tn.name);
				}
				
				if(!bestEffort && tn==null){//Try to get a higher-level node
					TaxNode tn2=taxtree.parseNodeFromHeader(rid, true);
					if(tn2!=null){
						while(tn2!=null && tn2.pid!=tn2.id && tn2.level<taxLevel){
							TaxNode temp=taxtree.getNode(tn2.pid);
							if(temp==null || temp.level>=TaxTree.LIFE){break;}
							tn2=temp;
						}
						if(tn2.level<=taxLevel){
							taxID=tn2.id;
						}
						if(smm.heap.taxID<0){smm.heap.taxID=tn2.id;}
						if(smm.heap.taxName()==null && tn2!=null){
							smm.heap.setTaxName(tn2.name);
						}
					}
				}
				//assert(!localMap.containsKey((long)taxID) || localMap.get((long)taxID)==smm.heap) : taxID; //123
				
				assert(smm.heap.taxID<0 || smm.heap.taxName()!=null) : smm.heap.taxID+", "+smm.heap.taxName()+", "+smm.heap.name()+", "+tn;

				//assert(!localMap.containsKey((long)taxID) || localMap.get((long)taxID)==smm.heap) : taxID; //123
				if(initialLength1>=k){smm.processRead(r1);}
				//assert(!localMap.containsKey((long)taxID) || localMap.get((long)taxID)==smm.heap) : taxID; //123
				if(initialLength2>=k){smm.processRead(r2);}
				//assert(!localMap.containsKey((long)taxID) || localMap.get((long)taxID)==smm.heap) : taxID; //123
				
				
				if(mode==PER_SEQUENCE){
					manageHeap_perSequence();
				}else if(mode==PER_TAXA || mode==PER_IMG){
					//assert(!localMap.containsKey((long)taxID) || localMap.get((long)taxID)==smm.heap) : taxID; //123
					manageHeap_perTaxa(taxID, imgID, unitSizeBases);
					if(localMap.size()>20){dumpLocalMap();}
				}else if(mode==ONE_SKETCH/* || mode==PER_HEADER*/){
					//do nothing
				}else{
					assert(false) : mode;
				}
			}
		}
		
		private void manageHeap_perSequence(){
			assert(mode==PER_SEQUENCE);
			writeHeap(smm.heap);
		}
		
		private void manageHeap_perTaxa(final int taxID, final long imgID, final long unitSizeBases){
			//assert(!localMap.containsKey((long)taxID) || localMap.get((long)taxID)==smm.heap) : taxID; //123
			assert(mode==PER_TAXA || mode==PER_IMG);
			
			if(smm.heap.size()<=0 || ((taxID<0 && imgID<0) && smm.heap.genomeSizeKmers<minSizeKmers)){//Discard
				assert(!localMap.containsKey((long)taxID));
				smm.heap.clear(true);
				return;
			}

			//TODO:
			//I could, at this point, write to disk if smm.heap.genomeSize==taxSize.
			//Or if the taxID is unknown.

			final boolean known=(taxID>=0 || imgID>=0);
			final boolean unknown=!known;
			final boolean hasSize=(known && (sizeList!=null || sizeMap!=null));
			boolean finished=(unknown || (hasSize && smm.heap.genomeSizeBases>=unitSizeBases));
			
			//For some reason, this assertion fired for a single taxID (52271) halfway through sketching RefSeq.
			//Likely a transient hardware error.
			assert(!finished || smm.heap.genomeSizeBases==unitSizeBases) : finished+", "+unknown+", "+hasSize+", "+(sizeList==null)+"\n"
					+taxID+", "+unitSizeBases+", "+smm.heap.genomeSizeBases+", "+smm.heap.genomeSizeKmers;
//			if(finished && smm.heap.genomeSizeBases!=unitSizeBases){
//				System.err.println("Warning: tid="+taxID+", finished="+finished+", known="+known+
//						", smm.heap.genomeSizeBases="+smm.heap.genomeSizeBases+", unitSizeBases="+unitSizeBases);
//			}
			
			smm.heap.taxID=taxID;
			smm.heap.imgID=imgID;
			
			//assert(!localMap.containsKey((long)taxID) || localMap.get((long)taxID)==smm.heap) : taxID; //123
			
			final Long key;
			if(imgID>-1 && (mode==PER_IMG || taxID<1)){key=imgID;}
			else if(taxID>-1){key=(long)taxID;}
			else{key=Long.valueOf(nextUnknown.getAndIncrement());}
			
			//assert(!localMap.containsKey((long)taxID) || localMap.get((long)taxID)==smm.heap) : taxID; //123

			if(unknown || finished){
				writeHeap(smm.heap);
				smm.heap.clear(true);
				localMap.remove((long)taxID);
				return;
			}
			
			//At this point, the taxID is known and this heap does not constitute the whole taxSize, or the taxSize is unknown.
			if(!hasSize){
				final HashMap<Long, SketchHeap> map=longMaps[(int)(key&MAP_MASK)];
				final SketchHeap old;
				synchronized(map){
					old=map.get(key);
					if(old==null){
						map.put(key, smm.heap);
					}else{
						old.add(smm.heap);
					}
				}
				if(old==null){
					smm.heap=null;
				}else{
					smm.heap.clear(true);
					assert(!localMap.containsKey((long)taxID));
				}
				localMap.remove((long)taxID);
				return;
			}
			
			//assert(!localMap.containsKey((long)taxID) || localMap.get((long)taxID)==smm.heap) : taxID; //123
			
			{
				//At this point, the taxID is and taxSize are known, and this heap does not constitute the whole taxSize.
				final int expectedHeapSize=toSketchSize(unitSizeBases, -1, -1, targetSketchSize);
				assert(expectedHeapSize>=3) : expectedHeapSize;
				//			boolean writeHeap=false;
				//			boolean makeHeap=false;

				SketchHeap local=null;
				{
					local=localMap.get(key);
					assert(local==null || local.taxID==key.intValue());
					if(local==smm.heap){
						//do nothing
						smm.heap=null;
					}else if(local==null){
						if(expectedHeapSize==smm.heap.capacity()){//Store the current heap
							assert(smm.heap.taxID==key.intValue());
							assert(smm.heap.name()!=null);
							localMap.put(key, smm.heap);//Safe
							smm.heap=null;
						}else{//Store a precisely-sized heap
							SketchHeap temp=new SketchHeap(expectedHeapSize, defaultParams.minKeyOccuranceCount, defaultParams.trackCounts());
							temp.add(smm.heap);
							assert(temp.taxID==key.intValue());
							assert(temp.name()!=null);
							localMap.put(key, temp);//Looks safe
//							smm.heap=null; //Not sure which is better
							smm.heap.clear(true);
						}
					}else{
						assert(local.taxID==smm.heap.taxID);
						assert(local!=smm.heap);
//						assert(!localMap.containsKey((long)taxID) || localMap.get((long)taxID)==smm.heap) : taxID+", "+key; //123
//						assert(false) : taxID+", "+key;
						local.add(smm.heap); //This WAS the slow line.  It should no longer be executed.
						//					if(local.genomeSizeBases>=unitSizeBases){
						//						writeHeap=true;
						//						localMap.remove(key);
						//					}
//						smm.heap=null; //Not sure which is better
						smm.heap.clear(true);
					}
				}
			}
		}
		
		private void dumpLocalMap(){

			for(Entry<Long, SketchHeap> e : localMap.entrySet()){
				boolean writeHeap=false;
				final Long key=e.getKey();
				final SketchHeap localHeap=e.getValue();
				final HashMap<Long, SketchHeap> map=longMaps[(int)(key&MAP_MASK)];
				final SketchHeap old;

				final long unitSizeBases;
				if(sizeList!=null){
					unitSizeBases=localHeap.taxID<0 ? -1 : sizeList.get((int)localHeap.taxID);
				}else if(sizeMap!=null){
					unitSizeBases=sizeMap.get(localHeap.imgID);
				}else{
					unitSizeBases=-1;
				}
				final int expectedHeapSize=toSketchSize(unitSizeBases, -1, -1, targetSketchSize);

				synchronized(map){
					old=map.get(key);
					if(old==null){
						if(expectedHeapSize==localHeap.capacity()){//Store the current heap
							map.put(key, localHeap);
						}else{//Store a precisely-sized heap
//							assert(key.intValue()!=96897) : key+", "+unitSizeBases+", "+", "+sizeList.get(key.intValue());
							assert(expectedHeapSize>0) : expectedHeapSize+", "+key+", "+localHeap.taxID+", "+localHeap.name()+", "+unitSizeBases;
							SketchHeap temp=new SketchHeap(expectedHeapSize, defaultParams.minKeyOccuranceCount, defaultParams.trackCounts());
							temp.add(localHeap);
							map.put(key, temp);
						}
					}else{
						old.add(localHeap); //This was the slow line; should be faster now.
						if(old.genomeSizeBases>=unitSizeBases){
							writeHeap=true;
							map.remove(key);
						}
					}
				}

				if(writeHeap){
					assert(old!=null); //For compiler
					assert(old.genomeSizeBases>0) : unitSizeBases+", "+old.genomeSizeBases+", "+old.genomeSizeKmers;
					assert(old.genomeSizeBases==unitSizeBases) : unitSizeBases+", "+old.genomeSizeBases+", "+old.genomeSizeKmers+", "+old.size()+", "+old.taxID;
					writeHeap(old);
				}
			}
			localMap.clear();
		}
		
		private boolean writeHeap(SketchHeap heap){
			sketchesMadeT++;
//			assert(heap.size()>0) : heap.size(); //Not really necessary
			boolean written=false;
//			assert(heap.heap.size()==heap.set.size()) : heap.heap.size()+", "+heap.set.size();
			if(heap.size()>0 && heap.genomeSizeKmers>=minSizeKmers && heap.sketchSizeEstimate()>0){
				Sketch sketch=new Sketch(heap, true, tool.trackCounts, outMeta);
				if(outTaxID>=0 && sketch.taxID<0){sketch.taxID=outTaxID;}
				if(outTaxName!=null && sketch.taxName()==null){sketch.setTaxName(outTaxName);}
				if(outFname!=null && sketch.fname()==null){sketch.setFname(outFname);}
				if(outName0!=null && sketch.name0()==null){sketch.setName0(outName0);}
				if(outSpid>=0 && sketch.spid<0){sketch.spid=outSpid;}
				if(outImgID>=0 && sketch.imgID<0){sketch.imgID=outImgID;}
				if(parseSubunit && sketch.name0()!=null){
					if(outMeta!=null){
						sketch.meta=(ArrayList<String>)sketch.meta.clone();
					}else if(sketch.meta==null){
						if(sketch.name0().contains("SSU_")){
							sketch.addMeta("subunit:ssu");
						}else if(sketch.name0().contains("LSU_")){
							sketch.addMeta("subunit:lsu");
						}
					}
				}
				if(tsw!=null){
					final int choice=(sketch.hashCode()&Integer.MAX_VALUE)%files;
					sketch.addSSU();
					synchronized(tsw[choice]){
						SketchTool.write(sketch, tsw[choice], bb);
						sketchesWrittenT++;
						written=true;
					}
				}
			}else{
				heap.clear(true);
			}
			assert(heap.genomeSizeBases==0);
			assert(heap.genomeSequences==0);
			return written;
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		long sketchesMadeT=0;
		long sketchesWrittenT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Thread ID */
		final int threadID;
		private ByteBuilder bb=new ByteBuilder();
		
		final SketchMakerMini smm;
		final HashMap<Long, SketchHeap> localMap;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;

	/** Primary output file path */
	private String out1=null;
	
	/** Override input file extension */
	private String extin=null;
	
	private String giTableFile=null;
	private String taxTreeFile=null;
	private String accessionFile=null;
	private String imgFile=null;
	
	/*Override metadata */
	String outTaxName=null;
	String outFname=null;
	String outName0=null;
	int outTaxID=-1;
	long outSpid=-1;
	long outImgID=-1;
	ArrayList<String> outMeta=null;
	static boolean parseSubunit=false;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	/** Number of bases processed */
	protected long kmersProcessed=0;
	/** Number of sketches started */
	protected long sketchesMade=0;
	/** Number of sketches written */
	protected long sketchesWritten=0;
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	final LongList sizeList;
	final HashMap<Long, Long> sizeMap;

	private HashMap<Long, SketchHeap> longMaps[];
	private ByteStreamWriter tsw[];
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/** Primary output files */
	private final FileFormat ffout[];
	/** Number of output files */
	private final int files;
	
	final int mode;
	
	private final SketchTool tool;
	
	/** Don't make sketches from sequences smaller than this */
	final int minSizeBases;
	/** Don't make sketches from sequences smaller than this */
	final int minSizeKmers;
	
	private int taxLevel=1;
	private boolean prefilter=false;
	private boolean tossJunk=true;
	boolean bestEffort=true;
//	private final HashMap<Integer, byte[]> ssuMap=null;
	
	private final AtomicInteger nextUnknown=new AtomicInteger(minFakeID);

	private static final int MAP_WAYS=32;
	private static final int MAP_MASK=MAP_WAYS-1;
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
