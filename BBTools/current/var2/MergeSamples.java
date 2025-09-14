package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import ml.CellNet;
import ml.CellNetParser;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.ListNum;
import structures.StringPair;
import var2.CallVariants2.Sample;

/**
 * Merges VCF files from multiple samples into a unified variant call set.
 * Creates a complete variant list where every position that has a variant
 * in ANY input sample is represented in the output, with sample-specific
 * information preserved.
 * 
 * This is particularly useful for multi-sample variant calling workflows where
 * you want to ensure that interesting variants from any individual sample are
 * evaluated across all samples, even if they weren't initially called in every sample.
 * 
 * The merging process:
 * 1. Synchronously reads corresponding lines from all input VCF files
 * 2. Aggregates statistical evidence across samples for each position
 * 3. Preserves individual sample genotype information
 * 4. Combines header metadata appropriately
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date December 18, 2016
 */
public class MergeSamples {
	
	/**
	 * Command-line entry point for VCF file merging operations.
	 * Parses arguments, configures the merger, and manages I/O streams.
	 * The main processing is currently commented out (line 58).
	 * @param args Command line arguments including input files, output paths, and neural network parameters
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		MergeSamples x=new MergeSamples(args);
		//x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Default constructor for programmatic instantiation.
	 * Initializes threading infrastructure using system-detected thread count.
	 * Creates blocking queue with capacity threads+1 to allow producer to stay ahead.
	 */
	public MergeSamples(){
		threads=Shared.threads();
		inq=new ArrayBlockingQueue<ListNum<VCFLine[]>>(threads+1);
	}
	
	/**
	 * Primary constructor that processes command-line arguments and configures the merger.
	 * Handles file I/O setup, threading configuration, neural network loading,
	 * and parameter validation. Sets up compressed I/O with PIGZ support.
	 * @param args Command line arguments including input/output files, neural network settings, and processing options
	 * @throws RuntimeException if required input files are missing or output files cannot be created
	 */
	public MergeSamples(String[] args){
		
		{ //Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		// Neural network parameters
		String netFile=null;
		boolean autoCutoff=true;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}
			
			// Neural network parameters
			else if(a.equals("net") || a.equals("netfile")){
				netFile=b;
				useNet=(b!=null);
			}else if(a.equals("netcutoff")){
				if("auto".equalsIgnoreCase(b)){
					autoCutoff=true;
				}else{
					autoCutoff=false;
					netCutoff=Float.parseFloat(b);
				}
			}else if(a.equals("usenet") || a.equals("useann") || a.equals("usenn") || a.equals("nn")){
				useNet=Parse.parseBoolean(b);
			}else if(a.equals("netmode")){
				useNet=(b!=null);
				if(b!=null){FeatureVectorMaker.setMode(b);}
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{ //Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;
			out1=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		// Load neural network if specified
		if(netFile!=null && useNet){
			net0=CellNetParser.load(netFile);
			assert(net0!=null) : "Failed to load neural network: "+netFile;
			if(autoCutoff){netCutoff=net0.cutoff;}
			if(verbose){outstream.println("Loaded neural network: "+netFile+" (cutoff="+netCutoff+")");}
		}else{
			net0=null;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		threads=Shared.threads();
		inq=new ArrayBlockingQueue<ListNum<VCFLine[]>>(threads+1);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Core Methods          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * High-level interface for merging VCF files from Sample objects (CallVariants2 integration).
	 * Extracts sample names and VCF paths from Sample objects and delegates to mergeFiles().
	 * This method serves as the primary entry point for multi-sample variant calling workflows
	 * where samples are represented as structured objects rather than simple file paths.
	 * @param list ArrayList of Sample objects, each containing sample name and VCF file path
	 * @param scafMap ScafMap for scaffold name resolution and coordinate validation
	 * @param outVcf Output merged VCF filename (may be null for stdout)
	 * @param scoreHistFile Optional filename for quality score histogram output (may be null)
	 */
	public void mergeSamples(ArrayList<Sample> list, ScafMap scafMap, String outVcf, String scoreHistFile){
		map=scafMap;
		ArrayList<StringPair> vcfList=new ArrayList<StringPair>(list.size());
		for(Sample s : list){vcfList.add(new StringPair(s.name, s.vcfName));}
		mergeFiles(vcfList, outVcf, scoreHistFile);
	}
	
	/**
	 * Core merging implementation that synchronously processes multiple VCF files.
	 * Algorithm requirements and assumptions:
	 * - All input VCF files must contain identical genomic positions in identical order
	 * - Files are typically generated by CallVariants2 using the same reference genome
	 * - Each file represents variant calls for a different sample at the same loci
	 * - Uses subprocess mode only for ≤4 input files to avoid resource exhaustion
	 * @param list ArrayList of StringPair objects containing (sample_name, vcf_filename) mappings
	 * @param outVcf Output merged VCF filename (null writes to stdout)
	 * @param scoreHistFile Optional quality score histogram filename (null disables histogram)
	 */
	public void mergeFiles(ArrayList<StringPair> list, String outVcf, String scoreHistFile){
		System.err.println("Merging "+list);
		final int ways=list.size();
		ByteFile[] bfa=new ByteFile[ways];
		final boolean allowSubprocess=(ways<=4);
		for(int i=0; i<ways; i++){
			StringPair pair=list.get(i);
			FileFormat ff=FileFormat.testInput(pair.b, FileFormat.VCF, null, allowSubprocess, false);
			bfa[i]=ByteFile.makeByteFile(ff);
		}
		
		mergeMT(outVcf, bfa);

		if(scoreHistFile!=null){
			CVOutputWriter.writeScoreHist(scoreHistFile, scoreArray);
		}
	}
	
	/**
	 * Single-threaded merging implementation (legacy/debugging).
	 * Processes variant rows sequentially without parallelization.
	 * Sequential control flow makes debugging easier, but performance is limited
	 * by single-threaded processing and frequent buffer flushing at 32KB threshold.
	 * @param outVcf Output merged VCF filename (null writes to stdout)
	 * @param bfa Array of ByteFile objects representing input VCF files in sample order
	 */
	private void mergeST(String outVcf, ByteFile[] bfa){
		ByteStreamWriter bswVcf=null;
		if(outVcf!=null){
			bswVcf=new ByteStreamWriter(outVcf, true, false, true, FileFormat.VCF);
			bswVcf.start();
		}
		
		ByteBuilder bb=new ByteBuilder(34000);
		VCFLine[] row=processRow(bfa, bb);
		while(row!=null){
			if(row[0]!=null){
				VCFLine merged=merge(row);
				merged.toText(bb);
				bb.nl();
				if(bb.length>32000){ //Flush buffer periodically
					if(bswVcf!=null){bswVcf.print(bb);}
					bb=new ByteBuilder(34000);
				}
			}
			row=processRow(bfa, bb);
		}
		
		if(bswVcf!=null){
			if(bb.length>0){bswVcf.print(bb);}
			bswVcf.poisonAndWait();
		}
	}
	
	/**
	 * Multithreaded merging implementation using producer-consumer architecture.
	 * Main thread reads and batches input rows while worker threads perform merging.
	 * Processes header lines synchronously first, then dispatches data batches of 200 rows
	 * to worker threads. Uses ArrayBlockingQueue for thread coordination and maintains
	 * output ordering via ByteStreamWriter's internal sequencing.
	 * @param outVcf Output merged VCF filename (null disables file output)
	 * @param bfa Array of ByteFile objects representing synchronized input VCF files
	 */
	private void mergeMT(String outVcf, ByteFile[] bfa){
		ByteStreamWriter bswVcf=null;
		if(outVcf!=null){
			FileFormat ff=FileFormat.testOutput(outVcf, FileFormat.VCF, null, true, true, append, true);
			bswVcf=new ByteStreamWriter(ff);
			bswVcf.start();
		}
		
		ArrayList<MergeThread> alpt=spawnThreads(bswVcf);
		
		long nextID=0;
		ByteBuilder header=new ByteBuilder(34000);
		
		// Process header lines first
		VCFLine[] row=processRow(bfa, header);
		while(row!=null && row[0]==null){ //Header
			row=processRow(bfa, header);
		}
		if(bswVcf!=null){
			bswVcf.add(header, nextID);
			nextID++;
		}
		
		// Process data lines in batches
		ListNum<VCFLine[]> list=new ListNum<VCFLine[]>(new ArrayList<VCFLine[]>(200), nextID);
		while(row!=null){
			if(row[0]!=null){
				list.add(row);
				if(list.size()>=200){
					putList(list);
					nextID++;
					list=new ListNum<VCFLine[]>(new ArrayList<VCFLine[]>(200), nextID);
				}
			}
			row=processRow(bfa, header);
		}
		if(list.size()>0){
			putList(list);
			nextID++;
		}
		
		putList(POISON_LIST);
		
		waitForFinish(alpt);
		
		if(bswVcf!=null){bswVcf.poisonAndWait();}
	}
	
	/**
	 * Synchronously reads corresponding lines from all input VCF files.
	 * Critical synchronization method that ensures positional alignment across samples.
	 * Returns null if any file reaches EOF (assumes all files have identical length).
	 * Header lines (starting with '#') are processed separately and trigger header merging.
	 * Data lines are converted to VCFLine objects with position validation.
	 * @param bfa Array of ByteFile objects for synchronized reading
	 * @param bb ByteBuilder for accumulating merged header content
	 * @return Array of VCFLine objects (one per sample), or null if EOF reached
	 * @throws AssertionError if variant positions don't match across files
	 */
	VCFLine[] processRow(ByteFile[] bfa, ByteBuilder bb){
		byte[][] lines=new byte[bfa.length][];
		for(int i=0; i<bfa.length; i++){
			byte[] line=bfa[i].nextLine();
			if(line==null){return null;}
			lines[i]=line;
		}
		
		VCFLine[] row=new VCFLine[bfa.length];
		if(lines[0][0]=='#'){
			processHeader(lines, bb);
			return row;
		}
		for(int i=0; i<lines.length; i++){
			byte[] line=lines[i];
			row[i]=new VCFLine(line);
			if(i>0){assert(row[i].pos==row[0].pos) : "\n"+row[0]+"\n"+row[i];}
		}
		return row;
	}
	
	/**
	 * Intelligently merges header metadata from multiple VCF files.
	 * Merging strategy by header type:
	 * - Statistical fields (reads, pairs, quality averages): Sum across samples then average
	 * - Format definitions: Copy from first file (assumed identical across samples)
	 * - Sample columns (#CHROM line): Concatenate sample names from all files
	 * - Ploidy/rate calculations: Aggregate then recalculate derived values
	 * @param lines Array of header line bytes (one corresponding line from each input file)
	 * @param bb ByteBuilder for accumulating merged header output
	 */
	void processHeader(byte[][] lines, ByteBuilder bb){
		String[][] matrix=new String[lines.length][];
		for(int i=0; i<lines.length; i++){
			matrix[i]=new String(lines[i]).split("=");
		}
		
		if(matrix[0][0].equals("##ploidy")){
			ploidy=Integer.parseInt(matrix[0][1]);
			bb.append("##ploidy="+ploidy+"\n");
		}else if(matrix[0][0].equals("##reads")){
			for(String[] split : matrix){
				reads+=Long.parseLong(split[1]);
			}
			bb.append("##reads="+reads+"\n");
		}else if(matrix[0][0].equals("##pairedReads")){
			for(String[] split : matrix){
				pairedReads+=Long.parseLong(split[1]);
			}
			bb.append("##pairedReads="+pairedReads+"\n");
		}else if(matrix[0][0].equals("##properlyPairedReads")){
			for(String[] split : matrix){
				properlyPairedReads+=Long.parseLong(split[1]);
			}
			properPairRate=properlyPairedReads*1.0/(Tools.max(1, reads));
			bb.append("##properlyPairedReads="+properlyPairedReads+"\n");
			bb.append("##properPairRate="+Tools.format("%.4f\n", properPairRate));
		}else if(matrix[0][0].equals("##properPairRate")){
			//do nothing - recalculated above
		}else if(matrix[0][0].equals("##totalQualityAvg")){
			totalQualityAvg=0;
			for(String[] split : matrix){
				totalQualityAvg+=Float.parseFloat(split[1]);
			}
			totalQualityAvg/=lines.length;
			bb.append("##totalQualityAvg="+Tools.format("%.4f\n", totalQualityAvg));
		}else if(matrix[0][0].equals("##mapqAvg")){
			mapqAvg=0;
			for(String[] split : matrix){
				mapqAvg+=Float.parseFloat(split[1]);
			}
			mapqAvg/=lines.length;
			bb.append("##mapqAvg="+Tools.format("%.2f\n", mapqAvg));
		}else if(matrix[0][0].equals("##readLengthAvg")){
			readLengthAvg=0;
			for(String[] split : matrix){
				readLengthAvg+=Float.parseFloat(split[1]);
			}
			readLengthAvg/=lines.length;
			bb.append("##readLengthAvg="+Tools.format("%.2f\n", readLengthAvg));
		}else if(matrix[0][0].startsWith("#CHROM\tPOS\t")){
			// Combine sample columns from all input files
			bb.append(lines[0]);
			for(int i=1; i<lines.length; i++){
				String[] split=new String(lines[i]).split("\t");
				bb.tab().append(split[split.length-1]);
			}
			bb.nl();
		}else{
			// Copy other header lines as-is from first file
			bb.append(lines[0]);
			bb.nl();
		}
	}
	
	/**
	 * Merges variant calls from multiple samples at identical genomic positions.
	 * Multi-phase merging algorithm:
	 * 1. Converts each VCFLine to Var objects for statistical arithmetic
	 * 2. Accumulates read counts, base qualities, mapping qualities across samples
	 * 3. Combines coverage depth information using addCoverage()
	 * 4. Generates merged VCF line with neural network scoring if enabled
	 * 5. Preserves per-sample genotype columns in final output
	 * 6. Uses max individual quality score if better than merged score
	 * 7. Updates quality score histogram for downstream analysis
	 * @param row Array of VCFLine objects representing same genomic position across all samples
	 * @return Single merged VCFLine with combined statistical evidence and preserved sample data
	 */
	VCFLine merge(VCFLine[] row){
		Var sum=null;
		VCFLine best=null;
		
		// Find best individual call and aggregate statistics
		for(VCFLine line : row){
			if(best==null || line.qual>best.qual){best=line;}
			Var v=line.toVar();
			assert(v!=null);
			if(sum==null){sum=v;}
			else{
				sum.add(v); // Aggregate statistical evidence
				sum.addCoverage(v); // Combine coverage information
			}
		}
		assert(best!=null);
		assert(sum!=null) : row.length+", "+row[0];
		
		// Generate merged VCF line with combined statistics
		ByteBuilder bb=sum.toVCF(new ByteBuilder(), properPairRate, totalQualityAvg, mapqAvg, 
				readLengthAvg, ploidy, map, filter, net0, trimWhitespace);
		VCFLine merged=new VCFLine(bb.toBytes());
		
		// Preserve individual sample information
		merged.samples.clear();
		for(VCFLine line : row){
			merged.samples.addAll(line.samples);
		}
		
		// Use best individual quality if better than merged
		if(merged.qual<best.qual){
			merged.qual=best.qual;
			merged.filter=best.filter;
		}
		
		scoreArray[(int)merged.qual]++;
		return merged;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Threading Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Blocking queue retrieval for consumer threads.
	 * Repeatedly attempts to take a batch from the processing queue,
	 * handling InterruptedException by continuing the wait loop.
	 * @return ListNum containing batch of variant rows to process
	 */
	final ListNum<VCFLine[]> takeList(){
		ListNum<VCFLine[]> list=null;
		while(list==null){
			try {
				list=inq.take();
			}catch (InterruptedException e){
				e.printStackTrace();
			}
		}
		return list;
	}
	
	/**
	 * Blocking queue insertion for producer thread.
	 * Repeatedly attempts to add a batch to the processing queue,
	 * handling InterruptedException by retrying until successful.
	 * @param list ListNum containing batch of variant rows to queue for processing
	 */
	final void putList(ListNum<VCFLine[]> list){
		while(list!=null){
			try {
				inq.put(list);
				list=null;
			}catch (InterruptedException e){
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Creates and starts worker threads for parallel variant merging.
	 * Each thread processes batches independently and outputs via ByteStreamWriter
	 * for ordered result assembly. Thread count matches system CPU availability.
	 * @param bsw ByteStreamWriter for coordinated output ordering (may be null)
	 * @return ArrayList of started MergeThread instances
	 */
	private ArrayList<MergeThread> spawnThreads(ByteStreamWriter bsw){
		ArrayList<MergeThread> alpt=new ArrayList<MergeThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new MergeThread(bsw));
		}
		if(verbose){outstream.println("Spawned threads.");}
		
		for(MergeThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
		
		return alpt;
	}
	
	/**
	 * Blocks until all worker threads reach TERMINATED state.
	 * Uses Thread.join() with InterruptedException handling to ensure
	 * all processing completes before method returns.
	 * @param alpt ArrayList of MergeThread instances to wait for completion
	 */
	private void waitForFinish(ArrayList<MergeThread> alpt){
		boolean allSuccess=true;
		for(MergeThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				}catch (InterruptedException e){
					e.printStackTrace();
				}
			}
		}
	}
	
	/**
	 * Worker thread implementing consumer side of producer-consumer pattern.
	 * Continuously processes batches of variant rows from the shared queue
	 * until receiving poison pill termination signal. Each thread operates
	 * independently with its own ByteBuilder buffer for output generation.
	 * Thread-safe operation via BlockingQueue coordination.
	 */
	private class MergeThread extends Thread {

		/**
		 * Constructs worker thread with output writer reference.
		 * @param bsw_ ByteStreamWriter for ordered output (may be null for testing)
		 */
		MergeThread(ByteStreamWriter bsw_){
			bsw=bsw_;
		}

		/**
		 * Main worker thread loop implementing consumer pattern.
		 * Continuously takes batches from queue, processes them, and passes poison pill
		 * to next thread when termination is received. This ensures all threads
		 * eventually receive termination signal.
		 */
		@Override
		/** Main thread execution loop */
		public void run(){
			ListNum<VCFLine[]> list=takeList();
			while(list!=null && list!=POISON_LIST){
				processList(list);
				list=takeList();
			}
			putList(POISON_LIST);
		}

		/**
		 * Processes a complete batch of variant rows with output buffering.
		 * Merges each row in the batch and accumulates output in ByteBuilder
		 * before sending to ByteStreamWriter with preserved batch ID for ordering.
		 * @param list ListNum batch containing variant rows and sequence ID
		 */
		private void processList(ListNum<VCFLine[]> list){
			ByteBuilder bb=new ByteBuilder(4096);
			for(VCFLine[] row : list){
				mergeRow(row, bb);
			}
			if(bsw!=null){bsw.add(bb, list.id);}
		}
		
		/**
		 * Merges variant calls from multiple samples at single genomic position.
		 * Skips null rows (header lines), otherwise delegates to parent merge() method
		 * and appends formatted output with newline termination.
		 * @param row Array of VCFLine objects from same genomic position across samples
		 * @param bb ByteBuilder for accumulating formatted output text
		 */
		private void mergeRow(VCFLine[] row, ByteBuilder bb){
			if(row[0]!=null){
				VCFLine merged=merge(row);
				merged.toText(bb);
				bb.nl();
			}
		}
		
		/** ByteStreamWriter for maintaining output order across parallel threads */
		private final ByteStreamWriter bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Poison pill sentinel for terminating thread processing loops */
	final ListNum<VCFLine[]> POISON_LIST=new ListNum<VCFLine[]>(null, -1);
	/** Thread-safe queue for distributing variant row batches to worker threads */
	private final ArrayBlockingQueue<ListNum<VCFLine[]>> inq;
	/** Number of worker threads for parallel processing (matches CPU count) */
	private final int threads;
	
	/** Total read count accumulated from all input samples (legacy field) */
	long readsSum;
	/** Total read pair count accumulated from all input samples (legacy field) */
	long pairsSum;
	/** Chromosome ploidy level for variant calling (typically 1 for haploid, 2 for diploid) */
	int ploidy=1;
	
	/** Aggregate proper pair rate calculated as properlyPairedReads/max(1,reads) */
	double properPairRate;
	/** Mean base quality score averaged across all input samples */
	double totalQualityAvg;
	/** Mean mapping quality score averaged across all input samples */
	double mapqAvg;
	/** Mean sequencing read length averaged across all input samples */
	double readLengthAvg;
	
	/** Cumulative read count aggregated from all sample headers */
	long reads;
	/** Cumulative paired read count aggregated from all sample headers */
	long pairedReads;
	/** Cumulative properly paired read count aggregated from all sample headers */
	long properlyPairedReads;
	
	/** Variant filtering criteria (currently unused but passed to Var.toVCF()) */
	VarFilter filter;
	/** Scaffold name mapping for coordinate system resolution */
	ScafMap map;
	/** Flag to trim whitespace from scaffold names during processing */
	boolean trimWhitespace=true;
	
	/** Primary input VCF filename (first sample, used for validation) */
	private String in1=null;
	/** Primary merged output VCF filename (null writes to stdout) */
	private String out1=null;
	/** Optional output filename for variants that fail filtering criteria */
	private String outInvalid=null;
	
	/** Quality score histogram tracking distribution of final merged variant scores (0-199) */
	long[] scoreArray=new long[200];
	
	/** Total VCF lines processed across all input files (unused in current implementation) */
	private long linesProcessed=0;
	/** Count of valid variant lines processed (unused in current implementation) */
	private long linesValid=0;
	/** Total bytes read from input files (unused in current implementation) */
	private long bytesProcessed=0;
	/** Maximum number of lines to process before stopping (Long.MAX_VALUE = unlimited) */
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------     Neural Network Fields    ----------------*/
	/*--------------------------------------------------------------*/

	/** Master neural network model for variant quality scoring (loaded once, shared across threads) */
	private CellNet net0=null;
	/** Flag enabling neural network-based variant quality scoring */
	private boolean useNet=false;
	/** Probability threshold for neural network variant classification (0.0-1.0) */
	private float netCutoff=0.5f;
	
	/** PrintStream for diagnostic and progress messages (System.err by default) */
	private PrintStream outstream=System.err;
	/** Global flag controlling detailed progress and debug output */
	public static boolean verbose=false;
	/** Flag indicating whether processing encountered critical errors */
	public boolean errorState=false;
	/** Flag allowing overwrite of existing output files */
	private boolean overwrite=true;
	/** Flag enabling append mode for output files (vs overwrite) */
	private boolean append=false;
}