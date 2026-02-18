package scalar;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.PriorityQueue;

import bin.DataLoader;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectDoubleMap;
import parse.LineParser1;
import parse.LineParserS1;
import parse.LineParserS4;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.IntList;
import structures.ListNum;
import structures.SuperLongList;
import tracker.KmerTracker;
import tracker.ReadStats;

/**
 * Partitions sequence reads into multiple output files for parallel processing.
 * Supports round-robin distribution or balanced distribution by base pairs.
 * Handles paired-end reads, PacBio subreads, and various sequence formats.
 * @author Brian Bushnell
 * @contributor Neptune
 * @date June 1, 2016
 */
public class PartitionReads {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Program entry point.
	 * Creates PartitionReads instance and processes input files according to command-line parameters.
	 * @param args Command-line arguments specifying input/output files and options
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		PartitionReads x=new PartitionReads(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor that parses command-line arguments and initializes processing parameters.
	 * Sets up input/output file formats, validates file paths, and configures partitioning options including number of ways, PacBio mode, and BP balancing.
	 * @param args Command-line arguments containing file paths and processing options
	 * @throws RuntimeException If required parameters are missing or files are invalid
	 */
	public PartitionReads(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		boolean setInterleaved=false; //Whether interleaved was explicitly set.

		//Set shared static variables
		Shared.setBufferLen(Tools.max(400, Shared.bufferLen()));
		Shared.capBuffers(2);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		SamLine.SET_FROM_OK=true;

		//Create a parser object
		Parser parser=new Parser();

		String qfout1=null, qfout2=null;

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ways")){
				ways=Integer.parseInt(b);
			}else if(a.equals("splitbybp") || a.equals("bybp") || a.equals("bp")){
				splitByBP=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("pacbio") || a.equals("subreads")){
				pacBioMode=Parse.parseBoolean(b);
			}else if(a.equals("partitionby") || a.equals("partition") || a.equals("splitby") || a.equals("mode")){
				mode=Tools.find(b==null ? null : b.toLowerCase(), modeNames);
				if(mode<0){throw new RuntimeException("Unknown partition mode: "+b+". Options: "+Arrays.toString(modeNames));}
			}else if(a.equals("cov") || a.equals("coverage") || a.equals("sam") || a.equals("bam")){
				if(FileFormat.isSamOrBamFile(b)){
					depthFile=b;
				}else{
					covFile=b;
				}
			}else if(a.equals("depth")){
				if(b==null || b.indexOf('.')<1){
					if(Parse.parseBoolean(b)){mode=DEPTH;}
				}else{
					if(FileFormat.isSamOrBamFile(b)){
						depthFile=b;
					}else{
						covFile=b;
					}
				}
			}else if(a.equals("cutoff") || a.equals("cutoffs")){
				String[] cutoffStrs=b.split(",");
				customCutoffs=new float[cutoffStrs.length];
				for(int j=0; j<cutoffStrs.length; j++){
					customCutoffs[j]=Float.parseFloat(cutoffStrs[j]);
				}
				//Validate cutoffs
				for(int j=0; j<customCutoffs.length; j++){
					if(customCutoffs[j]<0){
						throw new RuntimeException("Cutoff values cannot be negative: "+customCutoffs[j]);
					}
				}
				//Check for duplicates
				for(int j=0; j<customCutoffs.length-1; j++){
					for(int k=j+1; k<customCutoffs.length; k++){
						if(customCutoffs[j]==customCutoffs[k]){
							throw new RuntimeException("Duplicate cutoff values not allowed: "+customCutoffs[j]);
						}
					}
				}
				//Sort cutoffs ascending
				java.util.Arrays.sort(customCutoffs);
				int requestedWays=ways;
				ways=customCutoffs.length+1;  //Auto-set ways from cutoff count
				if(requestedWays>0 && requestedWays!=ways){
					outstream.println("Note: partitions="+requestedWays+" overridden by cutoff count; using "+ways+" partitions");
				}
			}else if(parser.in1==null && new File(arg).exists() && Tools.canRead(arg)){
				//First existing readable file becomes input
				FileFormat ff=FileFormat.testInput(arg, FileFormat.FASTA, null, false, false);
				if(ff!=null && (ff.fasta() || ff.fastq() || ff.samOrBam())){
					parser.in1=arg;
					outstream.println("Using "+arg+" as input.");
				}else{
					outstream.println("Unknown parameter "+args[i]);
					assert(false) : "Unknown parameter "+args[i];
				}
			}else if(parser.out1==null && (arg.contains("%") || arg.contains("#"))){
				//First file with % or # becomes output pattern
				parser.out1=arg;
				outstream.println("Using "+arg+" as output pattern.");
			}else if(new File(arg).exists() && Tools.canRead(arg)){
				//Detect file type for depth sources
				FileFormat ff=FileFormat.testInput(arg, FileFormat.SAM, null, false, false);
				if(ff!=null && ff.samOrBam()){
					depthFile=arg;
					outstream.println("Using "+arg+" for depth information.");
				}else if(looksLikeCovFile(arg)){
					covFile=arg;
					outstream.println("Using "+arg+" for depth information.");
				}else{
					outstream.println("Unknown file type: "+arg);
					assert(false) : "Unknown file type: "+arg;
				}
			}else if(Tools.find(a, modeNames)>=0){
				mode=Tools.find(a, modeNames);
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
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;

			extin=parser.extin;
			extout=parser.extout;
		}

		//Warn if no output pattern specified
		if(out1==null){
			outstream.println("*** Warning: No output pattern specified, no output will be generated. ***");
		}

		assert(ways>0) : "Ways must be at least 1.";

		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}

		//Do qfin # replacement
		if(qfin1!=null && qfin2==null && qfin1.indexOf('#')>-1 && !new File(qfin1).exists()){
			qfin2=qfin1.replace("#", "2");
			qfin1=qfin1.replace("#", "1");
		}

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}

		//Do output file # replacement
		if(qfout1!=null && qfout2==null && qfout1.indexOf('#')>-1){
			qfout2=qfout1.replace("#", "2");
			qfout1=qfout1.replace("#", "1");
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
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		assert(out1==null || out1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(out2==null || out2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout1==null || qfout1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout2==null || qfout2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";

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

		//Create output FileFormat objects
		if(out1==null){
			ffout1=ffout2=null;
			qfout1Array=qfout2Array=null;
		}else{
			ffout1=new FileFormat[ways];
			ffout2=new FileFormat[ways];
			qfout1Array=new String[ways];
			qfout2Array=new String[ways];
			for(int i=0; i<ways; i++){
				String temp1=out1==null ? null : out1.replaceFirst("%", ""+i);
				String temp2=out2==null ? null : out2.replaceFirst("%", ""+i);
				ffout1[i]=FileFormat.testOutput(temp1, FileFormat.FASTQ, extout, true, overwrite, append, false);
				ffout2[i]=FileFormat.testOutput(temp2, FileFormat.FASTQ, extout, true, overwrite, append, false);

				qfout1Array[i]=qfout1==null ? null : qfout1.replaceFirst("%", ""+i);
				qfout2Array[i]=qfout2==null ? null : qfout2.replaceFirst("%", ""+i);
			}
		}

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

		useSharedHeader=(ffin1.samOrBam() && ffout1!=null && ffout1.length>0 && ffout1[0]!=null && ffout1[0].samOrBam());

		//Set optimization flags based on mode
		needDimers=(mode==GC || mode==HH || mode==CAGA);
		needDepth=(mode==DEPTH);

		//Load depth map if in DEPTH mode with external files
		if(mode==DEPTH){
			if(covFile!=null){
				depthMap=loadCoverageFile(covFile);
			}else if(depthFile!=null){
				depthMap=loadDepthFromAlignment(depthFile);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates array of output writers based on partitioning mode.
	 * For single-dimension: creates 'ways' writers
	 * TODO: For multi-dimension: creates product of waysPerDimension
	 * @return Array of Writer objects, or null if no output specified
	 */
	private Writer[] getWriters(){
		if(ffout1==null || ffout1.length==0){return null;}

		final int numWriters=ways;  //TODO: calculate from waysPerDimension for multidimensional
		Writer[] ros=new Writer[numWriters];
		final int buff=1;

		for(int i=0; i<numWriters; i++){
			ros[i]=WriterFactory.getStream(ffout1[i], ffout2[i], qfout1Array[i], qfout2Array[i], buff, null, useSharedHeader, 1);
			ros[i].start();
		}
		return ros;
	}

	/**
	 * Main processing method that partitions reads into multiple output files.
	 * Creates input/output streams, processes reads using either round-robin or base pair balanced distribution, and reports timing statistics.
	 * @param t Timer for measuring execution time
	 * @throws RuntimeException If processing encounters errors
	 */
	void process(Timer t){

		//Create a read input stream
		final Streamer cris;
		{
			cris=StreamerFactory.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2, -1);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}

		//Create output writers
		final Writer[] ros=getWriters();

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Process the read stream
		if(splitByBP){mode=BP;}  //Legacy support for splitbybp flag

		if(mode>BP){  //Metric-based partitioning
			processInner_metric(cris, ros, mode);
		}else if(mode==BP){
			processInner_heap(cris, ros);
		}else{  //COUNT or default
			processInner(cris, ros);
		}

		if(verbose){outstream.println("Finished; closing streams.");}

		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros);

		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));

		//Print partition statistics
		if(ros!=null && ffout1!=null){
			printPartitionStatistics(ros, mode);
		}

		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Processes reads using round-robin distribution across output files.
	 * For PacBio mode, groups subreads by ZMW identifier into the same partition.
	 * Otherwise distributes reads sequentially across available output streams.
	 * @param cris Input stream for reading sequence data
	 * @param ros Array of output streams for writing partitioned reads
	 */
	void processInner(final Streamer cris, final Writer ros[]){

		//Do anything necessary prior to processing

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

			ArrayList<Read>[] outLists=new ArrayList[ways];
			for(int i=0; i<ways; i++){
				outLists[i]=new ArrayList<Read>();
			}

			int nextIndex=0;

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

				//Loop through each read in the list
				for(Read r1 : reads){
					final Read r2=r1.mate;

					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());

					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;

					if(pacBioMode){
						int zmw=Parse.parseZmw(r1.id);
						assert(zmw>=0) : "Invalid zmw in "+r1.id;
						nextIndex=(zmw%ways);
					}

					outLists[nextIndex].add(r1);
					nextIndex=(nextIndex+1)%ways;
				}

				//Output reads to the output stream
				for(int i=0; i<ways; i++){
					if(ros!=null){ros[i].add(outLists[i], ln.id);}
					outLists[i]=new ArrayList<Read>();
				}

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
		}

		//Do anything necessary after processing

	}

	/**
	 * Processes reads using base pair balanced distribution across output files.
	 * Uses a priority queue to track the number of base pairs in each partition and assigns new reads to the partition with the fewest base pairs.
	 * @param cris Input stream for reading sequence data
	 * @param ros Array of output streams for writing partitioned reads
	 */
	void processInner_heap(final Streamer cris, final Writer ros[]){

		//Do anything necessary prior to processing

		PriorityQueue<Partition> queue=new PriorityQueue<Partition>(ways);
		for(int i=0; i<ways; i++){
			queue.add(new Partition(i));
		}

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

			ArrayList<Read>[] outLists=new ArrayList[ways];
			for(int i=0; i<ways; i++){
				outLists[i]=new ArrayList<Read>();
			}

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

				//Loop through each read in the list
				for(Read r1 : reads){
					final Read r2=r1.mate;

					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());

					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;

					Partition p=queue.poll();
					outLists[p.id].add(r1);
					p.bp+=r1.pairLength();
					queue.add(p);
				}

				//Output reads to the output stream
				for(int i=0; i<ways; i++){
					if(ros!=null){ros[i].add(outLists[i], ln.id);}
					outLists[i]=new ArrayList<Read>();
				}
				if(verbose){outstream.println("Returned a list.");}

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
		}

		//Do anything necessary after processing

	}

	/**
	 * Builds histogram of metric distribution for balanced partitioning.
	 * @param cris Input stream (consumed during histogram building)
	 * @param mode Metric mode to calculate
	 * @return Histogram - long[1024] for bounded metrics (GC/HH/CAGA), SuperLongList for unbounded (LENGTH)
	 */
	private Object buildHistogram(final Streamer cris, int mode){
		final boolean unbounded=(mode==LENGTH);
		final int histSize=(mode==DEPTH ? 2048 : 1024);  //DEPTH needs more bins due to log quantization
		long[] histogram=unbounded ? null : new long[histSize];
		SuperLongList sll=unbounded ? new SuperLongList() : null;

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		while(ln!=null && reads!=null && reads.size()>0){
			if(verbose){outstream.println("Pass 1: Fetched "+reads.size()+" reads for histogram.");}

			for(Read r1 : reads){
				final Read r2=r1.mate;
				final int initialLength1=r1.length();
				final int initialLength2=(r1.mateLength());

				readsProcessed+=r1.pairCount();
				basesProcessed+=initialLength1+initialLength2;

				if(needDimers){
					dimers.clearAll();
					dimers.add(r1.bases);
					if(r2!=null){dimers.add(r2.bases);}
				}

				float metric=getMetricValue(r1, r2, dimers, mode);

				if(unbounded){
					sll.add((long)metric);
				}else{
					int bin;
					if(mode==DEPTH){
						bin=quantizeDepth(metric);  //Quantize real depth for histogram
						bin=Tools.min(bin, histSize-1);
					}else{
						bin=(int)(metric*1023);
						bin=Tools.min(bin, 1023);
					}
					histogram[bin]+=r1.pairLength();
				}
			}

			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}

		if(unbounded){
			sll.sort();  //Required before percentile operations!
			if(verbose){outstream.println("Pass 1 complete. Built SuperLongList histogram.");}
			return sll;
		}else{
			if(verbose){outstream.println("Pass 1 complete. Built histogram with "+histogram.length+" bins.");}
			return histogram;
		}
	}

	/**
	 * Partitions reads to output streams using pre-calculated boundaries.
	 * @param cris Input stream (creates new stream from beginning)
	 * @param ros Output writers
	 * @param mode Metric mode to calculate
	 * @param boundaries Partition boundaries from histogram analysis
	 */
	private void partitionReads(final Streamer cris, final Writer[] ros, int mode, float[] boundaries){
		final Streamer cris2=StreamerFactory.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2, -1);
		cris2.start();
		if(verbose){outstream.println("Started second stream for pass 2");}

		ListNum<Read> ln=cris2.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris2.paired());
		}

		ArrayList<Read>[] outLists=new ArrayList[ways];
		for(int i=0; i<ways; i++){
			outLists[i]=new ArrayList<Read>();
		}

		while(ln!=null && reads!=null && reads.size()>0){
			if(verbose){outstream.println("Pass 2: Fetched "+reads.size()+" reads for partitioning.");}

			for(Read r1 : reads){
				final Read r2=r1.mate;

				if(needDimers){
					dimers.clearAll();
					dimers.add(r1.bases);
					if(r2!=null){dimers.add(r2.bases);}
				}

				float metric=getMetricValue(r1, r2, dimers, mode);
				int partition=findPartitionIndex(metric, boundaries);

				outLists[partition].add(r1);
			}

			for(int i=0; i<ways; i++){
				if(ros!=null){ros[i].add(outLists[i], ln.id);}
				outLists[i]=new ArrayList<Read>();
			}

			ln=cris2.nextList();
			reads=(ln!=null ? ln.list : null);
		}

		errorState|=ReadWrite.closeStream(cris2);
	}

	/**
	 * Processes reads using metric-based balanced distribution across output files.
	 * Uses two-pass streaming: first pass builds histogram of metric distribution,
	 * second pass routes reads to partitions for balanced base pair distribution.
	 * @param cris Input stream for reading sequence data (consumed in pass 1)
	 * @param ros Array of output streams for writing partitioned reads
	 * @param mode Metric mode: GC, HH, CAGA, or LENGTH
	 */
	void processInner_metric(final Streamer cris, final Writer ros[], int mode){
		// Phase 1: Build histogram (returns long[] or SuperLongList)
		Object histogram=buildHistogram(cris, mode);

		// Phase 2: Calculate boundaries
		float[] boundaries;
		if(customCutoffs!=null){
			boundaries=customCutoffs;  //Use provided cutoffs
			if(verbose){outstream.println("Using custom cutoffs");}
		}else{
			boundaries=calculateBalancedBoundaries(histogram, ways, mode);  //Auto-calculate
		}

		if(verbose){
			outstream.println("Partition boundaries:");
			for(int i=0; i<boundaries.length; i++){
				outstream.println("  Partition "+i+": metric < "+String.format("%.4f", boundaries[i]));
			}
		}

		// Phase 3: Partition reads
		partitionReads(cris, ros, mode, boundaries);
	}
	void dimensionSplit(ArrayList<String> dimensionNames, IntList waysPerDimension, float overlap){
		final int dims=dimensionNames.size();
		//TODO:
		//1) Process the input file to make histograms of each dimension (total number of bp in each)
		//2) Determine optimal partition for each dimension so that each bucket will get a similar number of bp
		//3) Start the streams and process the data
	}

	/**
	 * Calculates the specified metric value for a read.
	 * @param r1 Primary read
	 * @param r2 Mate read (may be null)
	 * @param dimers KmerTracker for composition metrics (k=2); may be null for non-composition metrics
	 * @param mode Metric mode: GC, HH, CAGA, LENGTH
	 * @return Metric value (0.0-1.0 for composition metrics, unbounded for length)
	 */
	private float getMetricValue(Read r1, Read r2, KmerTracker dimers, int mode){
		switch(mode){
			case GC: return dimers.GC();
			case HH: return dimers.HH();
			case CAGA: return dimers.CAGA();
			case LENGTH: return r1.pairLength();
			case DEPTH:
				if(depthMap!=null){
					//Use external depth from coverage or alignment file - return REAL depth
					String name=Tools.trimToWhitespace(r1.id);
					double depth=depthMap.get(name);
					if(depth<0){
						//Not found in map - warn once and set to 0
						if(!warnedMissingCoverage){
							outstream.println("***Warning! "+name+" not found in coverage data; further warnings will be suppressed.***");
							warnedMissingCoverage=true;
						}
						depth=0;
					}
					return (float)depth;
				}else{
					//Parse depth from header - return REAL depth
					float depth=parseDepth(r1.id, lps, lpt);
					return (depth<0 ? 0 : depth);  //Handle not found
				}
			default: throw new RuntimeException("Unknown metric mode: "+mode);
		}
	}

	/**
	 * Quantizes depth value using log transform for compact histogram storage.
	 * @param depth Raw depth value
	 * @return Quantized depth level (0-10000 range, ~1% resolution)
	 */
	private static int quantizeDepth(float depth){
		if(depth<=0){return 0;}
		depth=Tools.min(depth, 65536f);  //Max depth cap
		float yf=((float)(Tools.log2(depth+0.0625f)+4));
		int level=(int)(yf*100);  //100 steps per log unit = ~1% resolution
		return level;
	}

	/**
	 * De-quantizes depth level back to real depth value.
	 * Inverse of quantizeDepth().
	 * @param level Quantized depth level
	 * @return Real depth value
	 */
	private static float dequantizeDepth(int level){
		if(level<=0){return 0;}
		float yf=level/100.0f;  //Reverse: level = yf*100
		float depth=(float)(Math.pow(2, yf-4) - 0.0625);  //Reverse: yf = log2(depth+0.0625)+4
		return Tools.max(0, depth);  //Ensure non-negative
	}

	/**
	 * Parses depth from read header in Spades or Tadpole format.
	 * @param name Read name/header
	 * @param lps LineParser for underscore-delimited format (Spades)
	 * @param lpt LineParser for comma/equals-delimited format (Tadpole)
	 * @return Depth value, or -1 if not found
	 */
	private static float parseDepth(String name, LineParserS1 lps, LineParserS4 lpt){
		if(name.startsWith("NODE_") && name.contains("_cov_")){//Spades
			lps.set(name);
			float depth=lps.parseFloat(5);
			return depth;
		}else if(name.startsWith("contig_") && name.contains(",cov=")){//Tadpole
			lpt.set(name);
			float depth=lpt.parseFloat(3);
			return depth;
		}else if(name.contains("_cov_")){//Generic
			lps.set(name);
			for(int i=0; i<lps.terms()-1; i++){
				if(lps.termEquals("cov", i)){
					i++;
					return lps.parseFloat(i);
				}
			}
		}
		return -1;
	}

	/**
	 * Calculates balanced partition boundaries from histogram to ensure equal base pairs per partition.
	 * @param histogramObj Histogram - long[] for bounded metrics, SuperLongList for unbounded
	 * @param ways Number of output partitions to create
	 * @return Array of boundary values; partition i contains reads with metric < boundaries[i]
	 */
	private float[] calculateBalancedBoundaries(Object histogramObj, int ways, int mode){
		if(histogramObj instanceof SuperLongList){
			//Unbounded metric (LENGTH, DEPTH) - use percentiles
			SuperLongList sll=(SuperLongList)histogramObj;
			float[] boundaries=new float[ways];
			for(int i=0; i<ways-1; i++){
				double fraction=(i+1)/(double)ways;
				boundaries[i]=sll.percentileValueBySum(fraction);
			}
			boundaries[ways-1]=Float.MAX_VALUE;  //Last boundary = infinity
			return boundaries;
		}else{
			//Bounded metric (GC, HH, CAGA) - use fixed bins
			long[] histogram=(long[])histogramObj;

			//Calculate total bp
			long total=0;
			for(long count : histogram){
				total+=count;
			}

			long targetPerPartition=total/ways;
			float[] boundaries=new float[ways];

			long accumulated=0;
			int partitionIndex=0;

			for(int i=0; i<histogram.length; i++){
				accumulated+=histogram[i];
				if(accumulated>=targetPerPartition*(partitionIndex+1) && partitionIndex<ways-1){
					boundaries[partitionIndex]=(mode==DEPTH ? dequantizeDepth(i+1) : (i+1)/1024.0f);
					partitionIndex++;
				}
			}
			boundaries[ways-1]=(mode==DEPTH ? Float.MAX_VALUE : 1.0f);  //Last boundary always 1.0

			return boundaries;
		}
	}

	/**
	 * Finds the appropriate partition index for a given metric value.
	 * @param metric Metric value in range 0.0-1.0
	 * @param boundaries Array of partition boundaries
	 * @return Partition index (0 to boundaries.length-1)
	 */
	private int findPartitionIndex(float metric, float[] boundaries){
		for(int i=0; i<boundaries.length-1; i++){
			if(metric<boundaries[i]){
				return i;
			}
		}
		return boundaries.length-1;
	}

	/**
	 * Checks if a file looks like a coverage file by examining the first line.
	 * Recognizes pileup (#ID) and covmaker (#Contigs) formats.
	 * @param fname Path to file to check
	 * @return True if file appears to be a coverage file
	 */
	private static boolean looksLikeCovFile(String fname){
		try{
			ByteFile bf=ByteFile.makeByteFile(fname, true);
			byte[] firstLine=bf.nextLine();
			bf.close();
			if(firstLine==null){return false;}
			String header=new String(firstLine);
			return header.startsWith("#ID") || header.startsWith("#Contigs") || header.startsWith("#Depths");
		}catch(Exception e){
			return false;
		}
	}

	/**
	 * Loads coverage file in either pileup or covmaker format.
	 * Detects format by checking first line header.
	 * @param fname Path to coverage file
	 * @return Map of contig name to depth value
	 */
	private static ObjectDoubleMap<String> loadCoverageFile(String fname){
		//Check format by reading first line
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		byte[] firstLine=bf.nextLine();
		bf.close();

		String header=new String(firstLine);
		if(header.startsWith("#ID")){
			//Pileup format
			return loadPileupCoverage(fname);
		}else if(header.startsWith("#Contigs") || header.startsWith("#Depths")){
			//Covmaker format - use DataLoader
			return DataLoader.loadCovFile2(fname);
		}else{
			throw new RuntimeException("Unknown coverage file format. Expected pileup (#ID) or covmaker (#Contigs) format.");
		}
	}

	/**
	 * Loads coverage from pileup.sh output format.
	 * Format: #ID header, then ID\tAvg_fold\t...
	 * @param fname Path to pileup file
	 * @return Map of contig name to average depth
	 */
	private static ObjectDoubleMap<String> loadPileupCoverage(String fname){
		System.err.print("Loading pileup coverage from "+fname+": ");
		Timer t=new Timer(System.err, false);
		LineParser1 lp=new LineParser1('\t');
		ByteFile bf=ByteFile.makeByteFile(fname, true);

		ObjectDoubleMap<String> map=new ObjectDoubleMap<String>();
		byte[] line=bf.nextLine();

		//Skip header
		while(line!=null && Tools.startsWith(line, '#')){
			line=bf.nextLine();
		}

		//Parse data lines
		int loaded=0;
		while(line!=null){
			lp.set(line);
			String name=Tools.trimToWhitespace(lp.parseString(0));
			float depth=lp.parseFloat(1);  //Avg_fold column
			map.put(name, depth);
			loaded++;
			line=bf.nextLine();
		}
		bf.close();
		t.stopAndPrint();
		System.err.println("Loaded "+loaded+" contigs from pileup file.");
		return map;
	}

	/**
	 * Calculates depth from SAM/BAM alignment file.
	 * Sums aligned base pairs per contig and divides by contig length.
	 * @param fname Path to SAM or BAM file
	 * @return Map of contig name to average depth
	 */
	private static ObjectDoubleMap<String> loadDepthFromAlignment(String fname){
		System.err.print("Loading depth from alignment file "+fname+": ");
		Timer t=new Timer(System.err, false);

		//First pass: get contig lengths from @SQ headers
		ObjectDoubleMap<String> lengthMap=new ObjectDoubleMap<String>();
		ObjectDoubleMap<String> bpMap=new ObjectDoubleMap<String>();

		ByteFile bf=ByteFile.makeByteFile(fname, true);
		byte[] line=bf.nextLine();

		//Parse @SQ headers for contig lengths
		while(line!=null && Tools.startsWith(line, '@')){
			String lineStr=new String(line);
			if(lineStr.startsWith("@SQ")){
				//Parse @SQ\tSN:name\tLN:length
				String[] parts=lineStr.split("\t");
				String name=null;
				int length=0;
				for(String part : parts){
					if(part.startsWith("SN:")){
						name=Tools.trimToWhitespace(part.substring(3));
					}else if(part.startsWith("LN:")){
						length=Integer.parseInt(part.substring(3));
					}
				}
				if(name!=null && length>0){
					lengthMap.put(name, length);
					bpMap.put(name, 0);  //Initialize bp count
				}
			}
			line=bf.nextLine();
		}

		//Second pass: count aligned bp per contig
		int readCount=0;
		while(line!=null){
			if(Tools.startsWith(line, '@')){
				line=bf.nextLine();
				continue;
			}

			//Parse SAM line: QNAME FLAG RNAME POS MAPQ CIGAR
			String lineStr=new String(line);
			String[] parts=lineStr.split("\t");
			if(parts.length>=6){
				int flag=Integer.parseInt(parts[1]);
				String rname=Tools.trimToWhitespace(parts[2]);
				String cigar=parts[5];

				//Skip unmapped reads (FLAG & 4) - invert check to avoid continue/infinite loop
				if((flag&4)==0){
					//Calculate aligned bases from CIGAR
					int alignedBases=calculateAlignedBases(cigar);
					if(alignedBases>0 && bpMap.contains(rname)){
						bpMap.put(rname, bpMap.get(rname)+alignedBases);
					}
				}
			}
			readCount++;
			if(readCount%100000==0){System.err.print(".");}
			line=bf.nextLine();
		}
		bf.close();

		//Calculate depth = bp/length
		ObjectDoubleMap<String> depthMap=new ObjectDoubleMap<String>();
		Object[] contigKeys=lengthMap.keys();
		for(Object obj : contigKeys){
			if(obj!=null){
				String contig=(String)obj;
				double length=lengthMap.get(contig);
				double bp=bpMap.get(contig);
				double depth=bp/length;
				depthMap.put(contig, depth);
			}
		}

		t.stopAndPrint();
		System.err.println("Loaded depth for "+depthMap.size()+" contigs from "+readCount+" alignments.");
		return depthMap;
	}

	/**
	 * Calculates aligned bases from CIGAR string.
	 * Counts M, =, X operations (matches/mismatches).
	 * @param cigar CIGAR string
	 * @return Number of aligned bases
	 */
	private static int calculateAlignedBases(String cigar){
		if(cigar.equals("*")){return 0;}
		int total=0;
		int num=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Character.isDigit(c)){
				num=num*10+(c-'0');
			}else{
				if(c=='M' || c=='=' || c=='X'){
					total+=num;
				}
				num=0;
			}
		}
		return total;
	}

	/**
	 * Prints partition statistics table showing reads and bases per output file.
	 * @param ros Array of output writers
	 * @param mode Partitioning mode for determining boundary display
	 */
	private void printPartitionStatistics(Writer[] ros, int mode){
		outstream.println("\nPartition Statistics:");
		outstream.println("File\tReads\tBases");

		for(int i=0; i<ways; i++){
			String filename=ffout1[i].name();
			long reads=ros[i].readsWritten();
			long bases=ros[i].basesWritten();

			outstream.println(filename+"\t"+reads+"\t"+bases);
		}
	}

	private static class Partition implements Comparable<Partition> {
		public Partition(final int id_) {id=id_;}
		/**
		 * Compares partitions by base pair count for priority queue ordering.
		 * Partitions with fewer base pairs have higher priority.
		 * Uses partition ID as tiebreaker for consistent ordering.
		 * @param o The partition to compare against
		 * @return Negative if this partition has fewer base pairs, positive if more, or ID difference if base pairs are equal
		 */
		@Override
		public int compareTo(Partition o) {
			long dif=bp-o.bp;
			return (dif>0 ? 1 : dif<0 ? -1 : id-o.id);
		}
		public final int id;
		public long bp=0;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;
	/** Secondary input file path for paired-end reads */
	private String in2=null;

	private String qfin1=null;
	private String qfin2=null;

	/**
	 * Primary output file path template with '%' placeholder for partition number
	 */
	private String out1=null;
	/**
	 * Secondary output file path template with '%' placeholder for partition number
	 */
	private String out2=null;

	private String[] qfout1Array=null;
	private String[] qfout2Array=null;

	private String extin=null;
	private String extout=null;

	/*--------------------------------------------------------------*/

	/** Number of reads processed so far */
	protected long readsProcessed=0;
	/** Number of bases processed so far */
	protected long basesProcessed=0;

	/** Maximum number of input reads to process; -1 means no limit */
	private long maxReads=-1;

	/** Number of output partitions to create */
	private int ways=-1;

	/** Keep PacBio subreads from the same ZMW together in the same partition */
	private boolean pacBioMode=false;

	/** Balance partitions by base pairs instead of read count */
	private boolean splitByBP=false;

	/** Partitioning mode: COUNT (round-robin), BP (balanced by base pairs), or metric-based (GC, HH, CAGA) */
	private int mode=COUNT;

	/** Coverage file path (pileup or covmaker format) for DEPTH mode */
	private String covFile=null;

	/** SAM/BAM file path for calculating depth in DEPTH mode */
	private String depthFile=null;

	/** Custom partition cutoffs (overrides automatic boundary calculation) */
	private float[] customCutoffs=null;

	/** Depth map loaded from coverage or alignment file (contig name -> depth) */
	private ObjectDoubleMap<String> depthMap=null;

	/** K-mer tracker for calculating composition metrics (GC, HH, CAGA); must use k=2 for these metrics */
	private final KmerTracker dimers=new KmerTracker(2, 0);

	/** Line parser for underscore-delimited headers (Spades format) */
	private final LineParserS1 lps=new LineParserS1('_');

	/** Line parser for comma/equals-delimited headers (Tadpole format) */
	private final LineParserS4 lpt=new LineParserS4("=,=,");

	/** Whether this mode requires dimer calculation (GC, HH, CAGA) */
	private final boolean needDimers;

	/** Whether this mode requires depth information (DEPTH) */
	private final boolean needDepth;

	/** Flag to warn once about missing coverage data */
	private boolean warnedMissingCoverage=false;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file format */
	private final FileFormat ffin1;
	/** Secondary input file format */
	private final FileFormat ffin2;

	/** Array of primary output file formats for each partition */
	private final FileFormat[] ffout1;
	/** Array of secondary output file formats for each partition */
	private final FileFormat[] ffout2;

	/*--------------------------------------------------------------*/
	/*----------------        Mode Constants        ----------------*/
	/*--------------------------------------------------------------*/

	static final int COUNT=0, BP=1, GC=2, HH=3, CAGA=4, LENGTH=5, DEPTH=6;
	static final String[] modeNames=new String[]{"count", "bp", "gc", "hh", "caga", "length", "depth"};

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private boolean useSharedHeader=false;
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging */
	public static boolean verbose=false;
	/** True if an error was encountered during processing */
	public boolean errorState=false;
	private boolean overwrite=true;
	/** Append to existing output files instead of overwriting */
	private boolean append=false;
	/** Flag for ordered processing (no effect on this singlethreaded program) */
	private final boolean ordered=false;

}
