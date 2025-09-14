package var2;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import bloom.KCountArray7MTA;
import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import ml.CellNet;
import ml.CellNetParser;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.SamReadStreamer;
import stream.SamStreamer;
import stream.SamStreamerMF;
import structures.ListNum;

/**
* Calls variants from one or more SAM or BAM files using multithreaded processing.
* Supports prefiltering, realignment, quality trimming, and comprehensive variant analysis.
* Outputs results in VAR, VCF, and GFF formats with detailed statistics and histograms.
* 
* Key features:
* - Multithreaded variant calling with configurable thread pools
* - Optional prefiltering using Bloom filter-like structures for performance
* - Read realignment and quality-based trimming
* - Comprehensive variant statistics and filtering
* - Support for forced variants from input VCF files
* - Multiple output formats (VAR, VCF, GFF)
* 
* @author Brian Bushnell
* @contributor Isla
* @date November 4, 2016
*/
public class CallVariants {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * Automatically delegates to CallVariants2 if multisample mode is detected.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		//Execute CallVariants2 instead if multisample variant-calling is needed
		if(preparseMulti(args)){
			CallVariants2.main(args);
			return;
		}
		
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		CallVariants x=new CallVariants(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Parses command line arguments to detect multisample mode.
	 * Used to determine whether to delegate to CallVariants2.
	 * @param args Command line arguments
	 * @return True if multisample mode is requested
	 */
	private static boolean preparseMulti(String[] args){
		boolean multi=false;
		for(String arg : args){
			if(arg.contains("multi")){
				String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if(b==null || b.equalsIgnoreCase("null")){b=null;}
				while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
				
				if(a.equals("multi") || a.equals("multisample")){
					multi=Parse.parseBoolean(b);
				}
			}
		}
		return multi;
	}
	
	/**
	* Constructor that parses command line arguments and initializes all parameters.
	* Sets up input/output files, filtering parameters, threading options, and validation.
	* 
	* @param args Command line arguments containing file paths and processing options
	*/
	public CallVariants(String[] args){
		
		{ //Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		// Configure SAM parsing to optimize for variant calling
		// Only parse essential fields to improve performance
		SamLine.PARSE_0=false; // Don't parse read names by default
		SamLine.PARSE_8=false; // Don't parse next segment info
		SamLine.PARSE_OPTIONAL_MD_ONLY=true; // Only parse MD tag from optional fields
		
		SamLine.RNAME_AS_BYTES=false;
		
		//Set shared static variables for file I/O optimization
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		ReadWrite.USE_BGZIP=true;
		
		//Create a parser object for quality trimming
		Parser parser=new Parser();
		parser.qtrimLeft=qtrimLeft;
		parser.qtrimRight=qtrimRight;
		parser.trimq=trimq;
		Shared.TRIM_READ_COMMENTS=Shared.TRIM_RNAME=true;
		Read.IUPAC_TO_N=true; // Convert ambiguous bases to N
		
		// Configure default SAM filtering parameters
		// These settings focus on high-quality, properly mapped reads
		samFilter.includeUnmapped=false; // Skip unmapped reads
		samFilter.includeSupplimentary=false; // Skip supplementary alignments
		samFilter.includeDuplicate=false; // Skip duplicate reads
		samFilter.includeNonPrimary=false; // Skip secondary alignments
		samFilter.includeQfail=false; // Skip quality-failed reads
		samFilter.minMapq=4; // Minimum mapping quality threshold
		String atomic="auto"; // Atomic scaffold access mode
		
		// Neural network parameters
		String netFile=null;
		boolean autoCutoff=true;
		
		//Parse each argument and configure corresponding parameters
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			// Basic operation parameters
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("multi") || a.equals("multisample")){
				boolean multi=Parse.parseBoolean(b);
				assert(!multi) : "\nThis program does not support multi-sample variant calling.\n";
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
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
				if(b!=null) {FeatureVectorMaker.setMode(b);}
			}
			
			// Streaming and threading parameters
			else if(a.equals("ss") || a.equals("samstreamer") || a.equals("streamer")){
				if(b!=null && Tools.isDigit(b.charAt(0))){
					useStreamer=true;
					streamerThreads=Tools.max(1, Integer.parseInt(b));
				}else{
					useStreamer=Parse.parseBoolean(b);
				}
			}else if(a.equals("ssmf") || a.equals("samstreamermf") || a.equals("streamermf")){
				if(b!=null && Tools.isDigit(b.charAt(0))){
					SamStreamerMF.MAX_FILES=Integer.parseInt(b);
					useStreamerMF=SamStreamerMF.MAX_FILES>1;
					if(useStreamerMF){useStreamer=true;}
				}else{
					useStreamerMF=Parse.parseBoolean(b);
					if(useStreamerMF){
						SamStreamerMF.MAX_FILES=Tools.max(2, SamStreamerMF.MAX_FILES);
						useStreamer=true;
					}
				}
			}else if(a.equals("sslistsize")){
				SamStreamer.LIST_SIZE=Parse.parseIntKMG(b);
				assert(SamStreamer.LIST_SIZE>0);
			}
			
			// Analysis and output parameters
			else if(a.equals("cc") || a.equals("calccoverage") || a.equals("coverage")){
				calcCoverage=Parse.parseBoolean(b);
			}else if(a.equals("parsename")){
				SamLine.PARSE_0=Parse.parseBoolean(b);
			}
			
			// Variant-specific parameters
			else if(a.equalsIgnoreCase("noPassDotGenotype") || a.equalsIgnoreCase("noPassDot")){
				Var.noPassDotGenotype=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("minVarCopies")){
				Var.MIN_VAR_COPIES=Integer.parseInt(b);
			}else if(a.equals("extended")){
				Var.extendedText=Parse.parseBoolean(b);
			}else if(a.equals("useidentity")){
				Var.useIdentity=Parse.parseBoolean(b);
			}else if(a.equals("usehomopolymer") || a.equals("homopolymer")){
				Var.useHomopolymer=Parse.parseBoolean(b);
			}else if(a.equals("usepairing")){
				Var.usePairing=Parse.parseBoolean(b);
			}else if(a.equals("usebias")){
				Var.useBias=Parse.parseBoolean(b);
			}else if(a.equals("nscan") || a.equals("donscan")){
				Var.doNscan=Parse.parseBoolean(b);
			}else if(a.equals("useedist")){
				Var.useEdist=Parse.parseBoolean(b);
			}else if(a.equals("prefilter")){
				prefilter=Parse.parseBoolean(b);
			}
			
			// File path parameters
			else if(a.equals("ref")){
				ref=b;
				// Handle special case for PhiX reference
				if(ref!=null && !Tools.isReadableFile(ref)) {
					if(ref.equalsIgnoreCase("phix")) {
						ref=Data.findPath("?phix2.fa.gz");
					}
				}
			}else if(a.equals("vcf") || a.equals("vcfout") || a.equals("outvcf")){
				vcf=b;
			}else if(a.equals("invcf") || a.equals("vcfin") || a.equals("forced")){
				vcfin=b;
			}else if(a.equals("gff") || a.equals("gffout") || a.equals("outgff")){
				gffout=b;
			}
			
			// Histogram output files
			else if(a.equals("scorehist") || a.equals("shist")){
				scoreHistFile=b;
			}else if(a.equals("zygosityhist") || a.equals("ploidyhist") || a.equals("zhist") || a.equals("phist")){
				zygosityHistFile=b;
			}else if(a.equals("qualityhist") || a.equals("qualhist") || a.equals("qhist")){
				qualityHistFile=b;
			}
			
			// Processing parameters
			else if(a.equals("border")){
				border=Integer.parseInt(b);
			}else if(a.equals("sample") || a.equals("samplename")){
				sampleName=b;
			}
			
			// Memory and performance parameters
			else if(a.equals("ca3") || a.equals("32bit")){
				Scaffold.setCA3(Parse.parseBoolean(b));
			}else if(a.equals("atomic")){
				atomic=b;
			}else if(a.equals("strandedcov") || a.equals("trackstrand") || a.equals("stranded")){
				Scaffold.setTrackStrand(Parse.parseBoolean(b));
			}
			
			// Realignment parameters
			else if(a.equals("realign")){
				realign=Parse.parseBoolean(b);
			}else if(a.equals("unclip")){
				unclip=Parse.parseBoolean(b);
			}else if(a.equals("realignrows") || a.equals("rerows")){
				Realigner.defaultMaxrows=Integer.parseInt(b);
			}else if(a.equals("realigncols") || a.equals("recols")){
				Realigner.defaultColumns=Integer.parseInt(b);
			}else if(a.equals("realignpadding") || a.equals("repadding") || a.equals("padding")){
				Realigner.defaultPadding=Integer.parseInt(b);
			}else if(a.equals("msa")){
				Realigner.defaultMsaType=b;
			}else if(a.equals("vmtlimit")){
				vmtSizeLimit=Parse.parseIntKMG(b);
			}
			
			// Filter parsing (handled by filter objects)
			else if(samFilter.parse(arg, a, b)){
				//do nothing - handled by samFilter
			}
			
			// Additional analysis parameters
			else if(a.equalsIgnoreCase("countNearbyVars")){
				countNearbyVars=Parse.parseBoolean(b);
			}

			// Input file parameters
			else if(a.equals("in") || a.equals("in1") || a.equals("in2")){
				assert(b!=null) : "Bad parameter: "+arg;
				// Handle both single files and comma-separated lists
				if(new File(b).exists()){in.add(b);}
				else{
					for(String s : b.split(",")){in.add(s);}
				}
			}else if(a.equals("list")){
				// Read input file paths from a text file (one per line)
				for(String line : TextFile.toStringLines(b)){
					in.add(line);
				}
			}

			// Filter management
			else if(a.equals("clearfilters")){
				if(Parse.parseBoolean(b)){
					varFilter.clear();
					samFilter.clear();
				}
			}

			// Delegate remaining parameter parsing to filter objects and parser
			else if(varFilter.parse(a, b, arg)){
				//do nothing - handled by varFilter
			}else if(parser.parse(arg, a, b)){ //Parse standard flags in the parser
				//do nothing - handled by parser
			}

			// Handle file arguments without explicit flags
			else if(arg.indexOf('=')<0 && (new File(arg).exists() || arg.indexOf(',')>0)){
				if(new File(arg).exists()){
					// Automatically categorize files by type
					if(FileFormat.isSamOrBamFile(arg)){
						in.add(arg);
					}else if(FileFormat.isFastaFile(arg) && (ref==null || ref.equals(arg))){
						ref=arg;
					}else{
						assert(false) : "Unknown parameter "+arg;
						outstream.println("Warning: Unknown parameter "+arg);
					}
				}else{
					// Handle comma-separated file lists
					for(String s : arg.split(",")){
						if(FileFormat.isSamOrBamFile(s)){
							in.add(s);
						}else{
							assert(false) : "Unknown parameter "+arg+" part "+s;
							outstream.println("Warning: Unknown parameter "+arg+" part "+s);
						}
					}
				}
			}else{
				// Unknown parameter - warn but continue
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Warning: Unknown parameter "+args[i]);
			}
		}

		// Post-processing of parsed parameters

		// Configure atomic scaffold access based on thread count
		if("auto".equalsIgnoreCase(atomic)){Scaffold.setCA3A(Shared.threads()>8);}
		else{Scaffold.setCA3A(Parse.parseBoolean(atomic));}

		// Validate and set default ploidy
		if(ploidy<1){System.err.println("WARNING: ploidy not set; assuming ploidy=1."); ploidy=1;}
		samFilter.setSamtoolsFilter();

		// Configure streaming thread count within reasonable bounds
		streamerThreads=Tools.max(1, Tools.min(streamerThreads, Shared.threads()));
		assert(streamerThreads>0) : streamerThreads;

		{ //Process parser fields and extract standard parameters
			Parser.processQuality();

			maxReads=parser.maxReads;
			overwrite=parser.overwrite;
			append=parser.append;

			out=parser.out1;
			// Auto-detect VCF output from file extension
			if(vcf==null && out!=null){
				if(ReadWrite.rawExtension(out).equals("vcf")){
					vcf=out;
					out=null;
				}
			}

			extin=parser.extin;
			extout=parser.extout;

			// Quality trimming parameters
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
			trimE=parser.trimE();

			trimWhitespace=Shared.TRIM_READ_COMMENTS;
		}

		// Disable strand tracking if no VCF output is specified
		if(vcf==null){Scaffold.setTrackStrand(false);}

		// Initialize ploidy array for zygosity statistics
		ploidyArray=new long[ploidy+1];

		// Load neural network if specified
		if(netFile!=null && useNet){
			net0=CellNetParser.load(netFile);
			assert(net0!=null) : "Failed to load neural network: "+netFile;
			if(autoCutoff){netCutoff=net0.cutoff;}
			if(verbose){outstream.println("Loaded neural network: "+netFile+" (cutoff="+netCutoff+")");}
		}else{
			net0=null;
		}

		// Validate FASTA reader settings
		assert(FastaReadInputStream.settingsOK());

		//Ensure there is an input file
		if(in.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}

		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out, vcf, gffout)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}

		// Fix file extensions for compressed files
		in=Tools.fixExtension(in);
		ref=Tools.fixExtension(ref);

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		// Validate forced variant input files
		if(vcfin!=null && !Tools.testInputFiles(false, true, vcfin.split(","))){
			throw new RuntimeException("\nCan't read vcfin: "+vcfin+"\n");  
		}

		//Create output FileFormat objects
		ffout=FileFormat.testOutput(out, FileFormat.VAR, extout, true, overwrite, append, true);

		//Create input FileFormat objects for all input files
		for(String s : in){
			FileFormat ff=FileFormat.testInput(s, FileFormat.SAM, extin, true, false);
			ffin.add(ff);
		}

		// Set default sample name from first input file if not specified
		if(sampleName==null){
			sampleName=ReadWrite.stripToCore(ffin.get(0).name());
		}

		// Reference file is required for variant calling
		assert(ref!=null) : "Please specify a reference fasta.";
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads the reference genome file if not already loaded.
	 * Creates a ScafMap object containing all reference scaffolds with their sequences,
	 * lengths, and metadata. Configures the realigner to use the loaded scaffold mapping
	 * for alignment score calculations during read realignment.
	 * 
	 * This method is idempotent - safe to call multiple times without side effects.
	 * Uses ScafMap.loadReference() which handles FASTA parsing and scaffold indexing.
	 */
	private void loadReference(){
		if(loadedRef){return;} // Skip if already loaded
		assert(ref!=null);
		// Load reference and create scaffold mapping with SAM filtering
		scafMap=ScafMap.loadReference(ref, scafMap, samFilter, true);
		// Configure realigner to use the loaded scaffold map
		if(realign){Realigner.map=scafMap;}
		loadedRef=true;
	}
	
	/**
	 * Loads reference file or SAM headers to create scaffold map.
	 * Implements a two-tier approach for scaffold discovery:
	 * 1. Primary: Load from reference FASTA file using FastaReadInputStream
	 * 2. Fallback: Extract scaffold information from SAM/BAM headers using ScafMap.loadSamHeader()
	 * 
	 * The FASTA approach provides full sequence data for realignment and variant validation,
	 * while the SAM header approach only provides scaffold names and lengths for basic mapping.
	 * 
	 * @param t2 Timer for tracking load time and performance measurement
	 */
	private void loadScafMap(Timer t2) {
		if(ref!=null){
			// Load from reference FASTA file (preferred method)
			t2.start("Loading reference.");
			loadReference();
			t2.stop("Time: ");
		}else{
			// Extract scaffold information from SAM/BAM headers
			// This is a fallback when no reference FASTA is provided
			for(FileFormat ff : ffin){
				ScafMap.loadSamHeader(ff, scafMap);
			}
		}
	}
	
	/**
	 * Creates and populates a prefilter to reduce memory usage for low-frequency variants.
	 * Implements a two-pass algorithm using KCountArray7MTA (Bloom filter-like counter array):
	 * 
	 * Pass 1: Count variant occurrences using probabilistic counters with configurable bit width
	 * Pass 2: Only process variants that exceed minReads threshold in main processing
	 * 
	 * Algorithm details:
	 * 1. Calculates optimal counter bit width: cbits = 2^n where 2^cbits >= minReads
	 * 2. Allocates ~1/8th of available memory for counter array (1 bit per byte)
	 * 3. Creates KCountArray7MTA with calculated parameters and 2 hash functions
	 * 4. Chooses single-file vs multi-file processing based on thread count and dataset size
	 * 5. Adds forced variants from input VCF to ensure they pass filtering
	 * 
	 * Memory efficiency: Reduces main VarMap memory usage by filtering out low-confidence variants
	 * before full processing, critical for large-scale variant calling on limited-memory systems.
	 * 
	 * @param minReads Minimum number of reads required to pass prefilter
	 * @param vm Existing VarMap containing forced variants (may be null)
	 * @return Populated KCountArray7MTA prefilter, or null if insufficient memory
	 */
	private KCountArray7MTA prefilter(int minReads, VarMap vm){
		// Calculate bits needed per counter to track up to minReads occurrences
		int cbits=2;
		while((1L<<cbits)-1<minReads){
			cbits*=2;
		}
		
		// Allocate approximately 1/8th of available memory for prefilter
		long mem=Shared.memAvailable(4);
		long prebits=mem; //1 bit per byte; 1/8th of the memory

		long precells=prebits/cbits;
		if(precells<100000){ //Not enough memory - no point in prefiltering
			return null;
		}
		
		// Create Bloom filter-like counter array
		KCountArray7MTA kca=new KCountArray7MTA(precells, cbits, 2, null, minReads);
		
		// Choose processing method based on threading and file count
		if(!useStreamer || !useStreamerMF || ffin.size()<2 || Shared.threads()<5 || (maxReads>=0 && maxReads<Long.MAX_VALUE)){
			prefilter_SF(kca); // Single-file processing
		}else{
			prefilter_MF(kca); // Multi-file processing
		}
		
		// Add forced variants from input VCF to ensure they pass prefilter
		if(vm!=null && vm.size()>0){ //For forced vars from an input VCF
			for(Var v : vm){
				final long key=v.toKey();
				kca.incrementAndReturnUnincremented(key, minReads);
			}
		}
		
		kca.shutdown();
		return kca;
	}
	
	/**
	 * Performs prefiltering using single-file processing mode.
	 * Processes each input file sequentially using multithreaded read processing.
	 * Creates either SamReadStreamer or ConcurrentReadInputStream based on useStreamer setting,
	 * then spawns worker threads that extract variants and increment counters in the KCountArray7MTA.
	 * 
	 * This mode is used for smaller datasets or when memory/threading constraints prevent
	 * multi-file processing. Each file is processed completely before moving to the next.
	 * 
	 * @param kca The prefilter counter array to populate with variant counts
	 */
	private void prefilter_SF(final KCountArray7MTA kca){
		// Process each input file individually
		for(FileFormat ff : ffin){

			/** Optional SamReadStreamer for high throughput */
			final SamReadStreamer ss;
			/** Shared input stream */
			final ConcurrentReadInputStream cris;
			
			// Set up input stream (either streamer or standard concurrent reader)
			if(useStreamer){
				cris=null;
				ss=new SamReadStreamer(ff, streamerThreads, false, maxReads);
				ss.start();
				if(verbose){outstream.println("Started streamer");}
			}else{
				ss=null;
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff, null);
				cris.start(); //Start the stream
				if(verbose){outstream.println("Started cris");}
			}

			final int threads=Shared.threads();
			
			//Fill a list with ProcessThreads for parallel processing
			ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
			for(int i=0; i<threads; i++){
				alpt.add(new ProcessThread(cris, ss, null, i, kca, true, net0)); // true = prefilterOnly mode
			}
			
			//Start the threads
			for(ProcessThread pt : alpt){
				pt.start();
			}
			
			//Wait for completion of all threads and accumulate statistics
			boolean success=true;
			for(ProcessThread pt : alpt){
				
				//Wait until this thread has terminated
				while(pt.getState()!=Thread.State.TERMINATED){
					try {
						pt.join();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
				varsProcessed+=pt.varsProcessedT;
				
				//Accumulate per-thread statistics
				success&=pt.success;
			}
			
			//Track whether any threads failed
			if(!success){errorState=true;}
		}
	}
	
	/**
	 * Performs prefiltering using multi-file processing mode.
	 * Processes multiple input files simultaneously using SamStreamerMF for maximum I/O parallelism.
	 * Worker threads pull reads from all files concurrently, improving throughput for large datasets
	 * with multiple input files by reducing I/O bottlenecks.
	 * 
	 * This mode is optimal for high-throughput variant calling with many input files,
	 * sufficient threads (>4), and adequate memory. Provides better resource utilization
	 * than sequential file processing.
	 * 
	 * @param kca The prefilter counter array to populate with variant counts
	 */
	private void prefilter_MF(final KCountArray7MTA kca){
		// Create multi-file streamer for simultaneous file processing
		SamStreamerMF ssmf=new SamStreamerMF(ffin.toArray(new FileFormat[0]), streamerThreads, false, maxReads);
		ssmf.start();

		final int threads=Shared.threads();

		//Fill a list with ProcessThreads for parallel processing
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(null, null, ssmf, i, kca, true, net0)); // true = prefilterOnly mode
		}

		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}

		//Wait for completion of all threads and accumulate statistics
		boolean success=true;
		for(ProcessThread pt : alpt){

			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			varsProcessed+=pt.varsProcessedT;

			//Accumulate per-thread statistics
			success&=pt.success;
		}

		//Track whether any threads failed
		if(!success){errorState=true;}
	}
	
	/** 
	* Main processing method that orchestrates the complete variant calling pipeline.
	* Loads reference data, creates variant maps, processes input files, and generates output.
	* @param t Timer for overall execution timing
	* @return VarMap containing all discovered and filtered variants
	*/
	public VarMap process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Reset counters for this processing run
		readsProcessed=0;
		basesProcessed=0;
		trimmedBasesProcessed=0;
		
		Timer t2=new Timer();
		
		// Load reference or SAM headers to create scaffold map
		loadScafMap(t2);
		
		// Create and populate the variant map with all processing
		long[] types=makeVarMap(t2);
		
		// Write output files (Var, VCF, GFF)
		CVOutputWriter.writeOutput(varMap, varFilter, ffout, vcf, gffout,
				readsProcessed, pairedInSequencingReadsProcessed, properlyPairedReadsProcessed,
				trimmedBasesProcessed, ref, trimWhitespace, sampleName);
		
		// Write histogram files
		CVOutputWriter.writeHistograms(scoreHistFile, zygosityHistFile, qualityHistFile,
				scoreArray, ploidyArray, avgQualityArray, maxQualityArray);
		
		//Reset read validation to original state
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		// Print timing and results summary
		printResults(types, t);
		
		//Throw an exception if there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
		
		return varMap;
	}

	/**
	* Creates and populates a prefilter to reduce memory usage for low-frequency variants.
	* Uses a Bloom filter-like structure (KCountArray7MTA) to track variants that appear
	* fewer than minReads times, allowing them to be filtered out early to save memory.
	* 
	* @param minReads Minimum number of reads required to pass prefilter
	* @param vm Existing VarMap containing forced variants (may be null)
	* @return Populated KCountArray7MTA prefilter, or null if insufficient memory
	*/
	private long[] makeVarMap(Timer t2) {
		varMap=new VarMap(scafMap);
		
		// Load forced variants from input VCF if specified
		if(vcfin!=null){
			AnalyzeVars.loadForcedVCF(vcfin, scafMap, varMap, outstream);
		}
		
		// Set up prefilter to reduce memory usage for low-frequency variants
		final KCountArray7MTA kca;
		if(prefilter){
			t2.start("Loading the prefilter.");
			kca=prefilter(varFilter.minAlleleDepth, vcfin==null ? null : varMap);
			double used=(100.0*kca.cellsUsed())/kca.cells;
			outstream.println("Added "+varsProcessed+" events to prefilter; approximately "+(long)(kca.estimateUniqueKmers(2))+" were unique.");
			outstream.println(Tools.format("The prefilter is %.2f%% full.", used));
			varsProcessed=0; // Reset counter for main processing
			t2.stop("Time: ");
			outstream.println();
		}else{
			kca=null;
		}
		
		t2.start("Processing input files.");
		
		// Choose optimal processing strategy based on dataset characteristics
		if(Shared.threads()>4 && useStreamer && useStreamerMF){
			processInput_MF(ffin.toArray(new FileFormat[0]), kca); // Multi-file mode for large datasets
		}else{
			for(FileFormat ff : ffin){
				processInput_SF(ff, kca); // Single-file mode for typical use cases
			}
		}
		
		// Calculate summary statistics for variant map metadata
		final double properPairRate=properlyPairedReadsProcessed/(double)Tools.max(1, readsProcessed-readsDiscarded);
		final double pairedInSequencingRate=pairedInSequencingReadsProcessed/(double)Tools.max(1, readsProcessed-readsDiscarded);
		final double totalQualityAvg=totalQualitySum/(double)Tools.max(1, trimmedBasesProcessed);
		final double totalMapqAvg=totalMapqSum/(double)Tools.max(1, readsProcessed-readsDiscarded);
		
		// Store calculated statistics in variant map for output generation
		varMap.ploidy=ploidy;
		varMap.properPairRate=properPairRate;
		varMap.pairedInSequencingRate=pairedInSequencingRate;
		varMap.totalQualityAvg=totalQualityAvg;
		varMap.totalMapqAvg=totalMapqAvg;
		varMap.readLengthAvg=trimmedBasesProcessed/(double)Tools.max(1, readsProcessed-readsDiscarded);
		t2.stop("Time: ");
		Shared.printMemory();
		outstream.println();
		
		long initialCount=varMap.size(); // Store count before filtering
		
		t2.start("Processing variants.");
		final long[] types=processVariants(); // Apply scoring and filtering
		t2.stop("Time: ");
		outstream.println();
		
		// Optional nearby variant analysis for filtering artifact clusters
		if(countNearbyVars){
			t2.start("Counting nearby variants.");
			int x=varMap.countNearbyVars(varFilter);
			if(x>0 && varFilter.failNearby){
				Arrays.fill(types, 0); // Reset counts since we're removing variants
				for(Var v : varMap.toArray(false)){
					if(!v.forced() && v.nearbyVarCount>varFilter.maxNearbyCount){
						varMap.removeUnsynchronized(v); // Remove likely artifacts
					}else{
						//TODO: Recalculate passing statistics
						//Only relevant in failNearby mode... not really important.
					}
				}
			}
			t2.stop("Time: ");
			outstream.println();
		}
		return types;
	}

	/**
	* Prints comprehensive timing and results summary to output stream.
	* Includes variant type breakdown, statistics, and performance metrics.
	* @param types Array of variant counts by type from processVariants()
	* @param t Main timer for overall execution time
	*/
	private void printResults(long[] types, Timer t) {
		t.stop();
		
		// Calculate basic statistics for reporting
		long size=scafMap.lengthSum(); // Total reference length
		long initialCount=varMap.size(); // Current variant count after filtering
		long a=initialCount, b=varMap.size(), c=varsPrefiltered, d=varsProcessed;
		double amult=100.0/a; // Percentage multiplier for initial count
		double bmult=100.0/b; // Percentage multiplier for final count
		long homozygousCount=(ploidy<2 ? shared.Vector.sum(ploidyArray) : ploidyArray[ploidyArray.length-1]);
		double homozygousRate=homozygousCount*1.0/shared.Vector.sum(ploidyArray);
		
		// Report prefilter effectiveness if used
		if(prefilter){
			outstream.println(c+" of "+d+" events were screened by the prefilter ("+Tools.format("%.4f%%", c*100.0/d)+").");
		}
		outstream.println(b+" of "+a+" variants passed primary filters ("+Tools.format("%.4f%%", b*amult)+").");
		outstream.println();
		
		// Calculate detailed statistics by variant type
		final long sub=types[Var.SUB], del=types[Var.DEL], ins=types[Var.INS];
		final double smult=1.0/Tools.max(1, sub), dmult=1.0/Tools.max(1, del), imult=1.0/Tools.max(1, ins); // Avoid division by zero
		final long jun=types[Var.LJUNCT]+types[Var.RJUNCT]+types[Var.BJUNCT]; // Total junction variants
		
		// Calculate average statistics per variant type
		final double subAD=ADArray[0][Var.SUB]*smult, delAD=ADArray[0][Var.DEL]*dmult, insAD=ADArray[0][Var.INS]*imult; // Allele depth
		final double subRD=ADArray[1][Var.SUB]*smult, delRD=ADArray[1][Var.DEL]*dmult, insRD=ADArray[1][Var.INS]*imult; // Reference depth
		final double subAF=AFArray[Var.SUB]*smult, delAF=AFArray[Var.DEL]*dmult, insAF=AFArray[Var.INS]*imult; // Allele frequency
		final double subScore=Tools.sumHistogram(scoreArray[Var.SUB+1])*smult; // Average scores
		final double delScore=Tools.sumHistogram(scoreArray[Var.DEL+1])*dmult;
		final double insScore=Tools.sumHistogram(scoreArray[Var.INS+1])*imult;
		final double subQual=Tools.sumHistogram(avgQualityArray[Var.SUB+1])*smult; // Average qualities
		final double delQual=Tools.sumHistogram(avgQualityArray[Var.DEL+1])*dmult;
		final double insQual=Tools.sumHistogram(avgQualityArray[Var.INS+1])*imult;
		
		// Print detailed variant type statistics table
		outstream.println("Type           \tCount\tRate\tAD\tDepth\tAF\tScore\tQual");
		outstream.println("Substitutions: \t"+sub+Tools.format("\t%.1f%%\t%."+(subAD>1000 ? 0 : 1)+"f\t%."+(subRD>1000 ? 0 : 1)+"f\t%.3f\t%.1f\t%.1f", 
				sub*bmult, subAD, subRD, subAF, subScore, subQual));
		outstream.println("Deletions:     \t"+del+Tools.format("\t%.1f%%\t%."+(delAD>1000 ? 0 : 1)+"f\t%."+(delRD>1000 ? 0 : 1)+"f\t%.3f\t%.1f\t%.1f", 
				del*bmult, delAD, delRD, delAF, delScore, delQual));
		outstream.println("Insertions:    \t"+ins+Tools.format("\t%.1f%%\t%."+(insAD>1000 ? 0 : 1)+"f\t%."+(insRD>1000 ? 0 : 1)+"f\t%.3f\t%.1f\t%.1f", 
				ins*bmult, insAD, insRD, insAF, insScore, insQual));
		if(var2.Var.CALL_JUNCTION){
			outstream.println("Junctions:     \t"+jun+Tools.format("\t%.1f%%", jun*bmult));
		}
		outstream.println("Variation Rate:\t"+(b==0 ? 0 : 1)+"/"+(size/Tools.max(1,b))); // Variants per base
		outstream.println("Homozygous:    \t"+homozygousCount+Tools.format("\t%.1f%%", homozygousRate*100)+"\n");
		
		// Report realignment statistics if enabled
		if(realign){
			outstream.println("Realignments:  \t"+realignmentsAttempted);
			outstream.println("Successes:     \t"+realignmentsSucceeded);
			outstream.println("Improvements:  \t"+realignmentsImproved);
			outstream.println("Retained:      \t"+realignmentsRetained);
			outstream.println();
		}
		
		// Final performance summary
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
	}

	/** 
	* Processes input using single-file mode with multithreaded read processing.
	* Creates input streams and spawns worker threads for variant detection.
	* @param ff Input file format to process
	* @param kca Prefilter for memory efficiency (may be null)
	*/
	void processInput_SF(FileFormat ff, KCountArray7MTA kca){
		assert(ff.samOrBam());

		/** Optional SamReadStreamer for high throughput */
		final SamReadStreamer ss;
		/** Shared input stream */
		final ConcurrentReadInputStream cris;
		
		// Set up appropriate input stream based on configuration
		if(useStreamer){
			cris=null;
			ss=new SamReadStreamer(ff, streamerThreads, false, maxReads);
			ss.start();
			if(verbose){outstream.println("Started streamer");}
		}else{
			ss=null;
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff, null);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		
		//Process the reads in separate threads
		spawnThreads(cris, ss, null, kca);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
	}

	/** 
	 * Processes input using multi-file mode for high-throughput datasets.
	 * Uses SamStreamerMF to simultaneously read from multiple files, maximizing I/O parallelism
	 * and reducing processing time for large multi-file datasets.
	 * 
	 * This mode is optimal when:
	 * - Multiple input files are available
	 * - Sufficient threads are available (>4)
	 * - I/O bandwidth is the limiting factor
	 * 
	 * Worker threads pull reads from all files concurrently, improving resource utilization
	 * compared to sequential file processing.
	 * 
	 * @param ff Array of input file formats to process simultaneously
	 * @param kca Prefilter for memory efficiency (may be null)
	 */
	void processInput_MF(FileFormat[] ff, KCountArray7MTA kca){
		assert(useStreamer);
		assert(ff[0].samOrBam());

		// Create multi-file streamer for simultaneous file processing
		/** Optional SamStreamerMF for very high throughput */
		final SamStreamerMF ssmf;
		ssmf=new SamStreamerMF(ff, streamerThreads, false, maxReads);
		ssmf.start();
		if(verbose){outstream.println("Started streamer");}
		
		//Process the reads in separate threads
		spawnThreads(null, null, ssmf, kca);
		
		if(verbose){outstream.println("Finished; closing streams.");}
	}
	
	/**
	 * Processes all variants in the variant map using multithreaded scoring and filtering.
	 * Delegates to VarMap.processVariantsMT() which applies statistical filters,
	 * neural network scoring (if enabled), and generates histogram data for output.
	 * 
	 * @return Array of variant counts by type after filtering
	 */
	private long[] processVariants(){
		return varMap.processVariantsMT(varFilter, net0, scoreArray, ploidyArray, avgQualityArray, 
				maxQualityArray, ADArray, AFArray);
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final SamReadStreamer ss, final SamStreamerMF ssmf, final KCountArray7MTA kca){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ss, ssmf, i, kca, false, net0));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
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
			trimmedBasesProcessed+=pt.trimmedBasesProcessedT;
			readsDiscarded+=pt.readsDiscardedT;
			pairedInSequencingReadsProcessed+=pt.pairedInSequencingReadsProcessedT;
			properlyPairedReadsProcessed+=pt.properlyPairedReadsProcessedT;
			varsPrefiltered+=pt.prefilteredT;
			varsProcessed+=pt.varsProcessedT;
			totalQualitySum+=pt.totalQualitySumT;
			totalMapqSum+=pt.totalMapqSumT;
			success&=pt.success;
			if(pt.realigner!=null){
				realignmentsAttempted+=pt.realigner.realignmentsAttempted;
				realignmentsImproved+=pt.realigner.realignmentsImproved;
				realignmentsSucceeded+=pt.realigner.realignmentsSucceeded;
				realignmentsRetained+=pt.realigner.realignmentsRetained;
			}
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Dumps thread-local variants to the main variant map in a thread-safe manner.
	 * Transfers all variants from the thread's local HashMap to the global VarMap,
	 * merging duplicate variants and clearing the local map for continued processing.
	 * 
	 * @param mapT Thread-local variant map to dump
	 * @return Number of variants added to the main map
	 */
	private int dumpVars(HashMap<Var, Var> mapT){
		int added=varMap.dumpVars(mapT);
		assert(mapT.size()==0);
		return added;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final SamReadStreamer ss_, final SamStreamerMF ssmf_,
				final int tid_, final KCountArray7MTA kca_, final boolean prefilterOnly_, final CellNet net0_){
			cris=cris_;
			ss=ss_;
			ssmf=ssmf_;
			tid=tid_;
			kca=kca_;
			prefilterOnly=prefilterOnly_;
			net0=net0_;
			realigner=(realign ? new Realigner() : null);
		}
		
		@Override
		public void run(){
			
			// Initialize neural network for this thread
			if(net0!=null && useNet){
				net=net0.copy(false);
				netVec=new float[net.numInputs()];
				if(verbose){System.err.println("Thread "+tid+" initialized neural network");}
			}
			
			//Process the reads
			if(ss!=null){processInner_ss();}
			else if(cris!=null){processInner_cris();}
			else{processInner_ssmf();}
			
			//Do anything necessary after processing
			if(!varMapT.isEmpty()){
				dumpVars(varMapT);
			}
			assert(varMapT.isEmpty());
			
			//Indicate successful exit status
			success=true;
		}
		
		/**
		 * Processes read lists from ConcurrentReadInputStream in a loop.
		 * Continuously fetches read lists, processes each read individually,
		 * and manages proper list return to the input stream for memory efficiency.
		 */
		void processInner_cris(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){ //ln!=null prevents a compiler potential null access warning

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					assert(r.mate==null);
					
					if(!r.validated()){r.validate(true);}
					
					//Track the initial length for statistics
					final int initialLength=r.length();

					//Increment counters
					readsProcessedT++;
					basesProcessedT+=initialLength;
					
					boolean b=processRead(r);
					
					if(!b){
						readsDiscardedT++;
					}
				}

				//Notify the input stream that the list was used
				cris.returnList(ln);

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		/** Iterate through the reads */
		void processInner_ss(){
			
			//Grab the actual read list from the ListNum
			ListNum<Read> ln=ss.nextList();

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
				ArrayList<Read> reads=ln.list;

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					assert(r.mate==null);
					
					if(!r.validated()){r.validate(true);}
					
					//Track the initial length for statistics
					final int initialLength=r.length();

					//Increment counters
					readsProcessedT++;
					basesProcessedT+=initialLength;
					
					boolean b=processRead(r);
					if(!b){
						readsDiscardedT++;
					}
				}
				
				ln=ss.nextList();
			}
		}
		
		/** Iterate through the reads */
		void processInner_ssmf(){
			
			//Grab the actual read list from the ListNum
			ListNum<Read> ln=ssmf.nextList();

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
				ArrayList<Read> reads=ln.list;
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static accessmf

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					assert(r.mate==null);
					
					if(!r.validated()){r.validate(true);}
					
					//Track the initial length for statistics
					final int initialLength=r.length();

					//Increment counters
					readsProcessedT++;
					basesProcessedT+=initialLength;
					
					boolean b=processRead(r);
					if(!b){
						readsDiscardedT++;
					}
				}
				
				ln=ssmf.nextList();
			}
		}
		
		/**
		 * Processes a single read to extract variants and update statistics.
		 * Implements the core variant calling algorithm for individual reads:
		 * 
		 * Processing steps:
		 * 1. Validates read and applies SAM filtering criteria
		 * 2. Updates pairing and sequencing statistics
		 * 3. Performs read realignment if enabled
		 * 4. Calculates border and quality-based trimming amounts
		 * 5. Extracts variants using Var.toVars() with scaffold context
		 * 6. For prefilter mode: Increments variant counters in KCountArray7MTA
		 * 7. For main mode: Adds variants to thread-local map if they pass depth filter
		 * 
		 * Quality trimming algorithm: Combines border trimming (near scaffold ends)
		 * with Phred quality-based trimming using TrimRead.testOptimal().
		 * 
		 * Memory management: Uses thread-local variant map with periodic dumping
		 * to main VarMap when size exceeds vmtSizeLimit.
		 * 
		 * @param r Read to process (must be validated SAM read)
		 * @return True if read was processed successfully, false if discarded
		 */
		boolean processRead(final Read r){
			if(r.bases==null || r.length()<=1){return false;}
			final SamLine sl=r.samline;
			
			if(samFilter!=null && !samFilter.passesFilter(sl)){return false;}
			
			if(sl.properPair()){properlyPairedReadsProcessedT++;}
			if(sl.hasMate()){pairedInSequencingReadsProcessedT++;}
			final Scaffold scaf=scafMap.getScaffold(sl);
			final int scafnum=scaf.number;
			
			if(realign) {realigner.realign(r, sl, scaf, unclip);}
			
			//TODO ****** I am trimming of all left-clipping after this.
			//TODO Also, I'm adjusting the stop rather than start if on the revere strand!
//			assert(false) : new String(new String(r.match)+"\n"+r.start+"\n"+r.stop+"\n"+r.obj+"\n");
			
			int leftTrimAmount=border, rightTrimAmount=border;
			if(border>0){
				int skipTrimRange=Tools.max(10, border+5);
				if(r.start<skipTrimRange){
					if(r.strand()==Shared.PLUS){leftTrimAmount=0;}
					else{rightTrimAmount=0;}
				}
				if(r.stop>scaf.length-skipTrimRange){
					if(r.strand()==Shared.PLUS){rightTrimAmount=0;}
					else{leftTrimAmount=0;}
				}
			}
			if(qtrimLeft || qtrimRight){
				long packed=TrimRead.testOptimal(r.bases, r.quality, trimE);
				if(qtrimLeft){leftTrimAmount=Tools.max(leftTrimAmount, (int)((packed>>32)&0xFFFFFFFFL));}
				if(qtrimRight){rightTrimAmount=Tools.max(rightTrimAmount, (int)((packed)&0xFFFFFFFFL));}
			}
			
			int trimmed=(leftTrimAmount<1 && rightTrimAmount<1 ? 0 : TrimRead.trimReadWithMatch(r, sl, leftTrimAmount, rightTrimAmount, 0, scaf.length, false));
			if(trimmed<0){return false;} //In this case the whole read should be trimmed
			int extra=(qtrimLeft || qtrimRight) ? trimmed/2 : Tools.min(border, trimmed/2);
			ArrayList<Var> vars=Var.toVars(r, sl, callNs, scafnum);
			
			if(prefilterOnly){
				if(vars==null){return true;}
				for(Var v : vars){
					long key=v.toKey();
					kca.increment(key);
				}
			}else{
				trimmedBasesProcessedT+=r.length();
				totalQualitySumT+=shared.Vector.sum(r.quality);
				totalMapqSumT+=sl.mapq;
				if(calcCoverage){scaf.add(sl);}
				if(vars==null){return true;}

				for(Var v : vars){
//					assert(false) : v.alleleCount();
					int depth=Integer.MAX_VALUE;
					if(kca!=null){
						depth=kca.read(v.toKey());
					}
					if(depth>=varFilter.minAlleleDepth){
						v.endDistMax+=extra;
						v.endDistSum+=extra;

						Var old=varMapT.get(v);
						if(old==null){varMapT.put(v, v);}
						else{old.add(v);}
					}else{
						prefilteredT++;
					}
				}
				if(varMapT.size()>vmtSizeLimit){
					dumpVars(varMapT);
				}
			}
			varsProcessedT+=vars.size();
			return true;
		}

		private final KCountArray7MTA kca;
		private final boolean prefilterOnly;
		
		/** Number of vars blocked by the prefilter */
		protected long prefilteredT=0;
		/** Number of vars processed */
		protected long varsProcessedT=0;
		
		/** Sum of trimmed, mapped base qualities */
		protected long totalQualitySumT=0;
		/** Sum of mapqs */
		protected long totalMapqSumT=0;
		
		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		/** Number of trimmed, mapped bases processed. */
		protected long trimmedBasesProcessedT=0;
		/** Number of reads discarded by this thread */
		protected long readsDiscardedT=0;
		/** Number of paired reads processed by this thread, whether or not they mapped as pairs */
		protected long pairedInSequencingReadsProcessedT=0;
		/** Number of properly paired reads processed by this thread */
		protected long properlyPairedReadsProcessedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		HashMap<Var, Var> varMapT=new HashMap<Var, Var>();
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Optional SamReadStreamer for high throughput */
		private final SamReadStreamer ss;
		/** Optional SamStreamerMF for very high throughput */
		private final SamStreamerMF ssmf;
		/** For realigning reads */
		final Realigner realigner;
		
		private final CellNet net0;
		private CellNet net;
		private float[] netVec;
		
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** List of input SAM/BAM file paths */
	private ArrayList<String> in=new ArrayList<String>();

	/** Primary output file path for variant data */
	private String out=null;

	/** VCF format output file path */
	private String vcf=null;

	/** VCF input file path for forced variants */
	private String vcfin=null;

	/** GFF format output file path */
	private String gffout=null;

	/** GFF input file path (unused) */
	private String gffin=null;
	
	/** Output file for variant score histogram */
	private String scoreHistFile=null;
	/** Output file for zygosity/ploidy histogram */
	private String zygosityHistFile=null;
	/** Output file for base quality histogram */
	private String qualityHistFile=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	/** Reference FASTA file path */
	private String ref=null;
	
	/** Flag indicating whether reference has been loaded */
	private boolean loadedRef=false;

	/** Enable quality trimming from left end of reads */
	private boolean qtrimLeft=false;
	/** Enable quality trimming from right end of reads */
	private boolean qtrimRight=true;
	/** Quality threshold for trimming (Phred scale) */
	private float trimq=10;
	/** Quality threshold converted to error probability */
	private final float trimE;
	
	/*--------------------------------------------------------------*/

	/** Total number of reads processed */
	protected long readsProcessed=0;
	/** Total number of bases in processed reads */
	protected long basesProcessed=0;
	/** Total number of trimmed, mapped bases processed */
	protected long trimmedBasesProcessed=0;
	/** Number of reads discarded by filters */
	protected long readsDiscarded=0;
	/** Number of reads that were paired in sequencing */
	protected long pairedInSequencingReadsProcessed=0;
	/** Number of properly paired reads processed */
	protected long properlyPairedReadsProcessed=0;
	/** Number of variants filtered out by prefilter */
	protected long varsPrefiltered=0;
	/** Total number of variants encountered during processing */
	protected long varsProcessed=0;
	
	/** Sum of all base quality scores from trimmed, mapped bases */
	protected long totalQualitySum=0;
	/** Sum of all mapping quality scores */
	protected long totalMapqSum=0;
	
	/** Number of realignment attempts made */
	protected long realignmentsAttempted;
	/** Number of realignments that improved alignment score */
	protected long realignmentsImproved;
	/** Number of successful realignment operations */
	protected long realignmentsSucceeded;
	/** Number of realigned reads that were kept */
	protected long realignmentsRetained;

	/** Maximum reads to process (-1 for unlimited) */
	private long maxReads=-1;
	
	/** Scaffold mapping for reference sequences */
	public ScafMap scafMap;
	/** Container for all discovered variants */
	public VarMap varMap;
	
	/** Whether to calculate coverage statistics */
	public boolean calcCoverage=true;

	/** Expected ploidy level for variant calling */
	public int ploidy=-1;
	
	/** Number of bases to trim from read ends near scaffold boundaries */
	public int border=5;

	/** Enable read realignment for improved variant detection */
	public boolean realign=false;
	/** Remove soft clipping during realignment */
	public boolean unclip=false;
	
	/** Enable Bloom filter prefiltering to reduce memory usage */
	public boolean prefilter=false;
	
	/** Sample name for output files */
	public String sampleName=null;
	
	/** Enable analysis of nearby variant clusters for artifact detection */
	public boolean countNearbyVars=true;
	
	/*--------------------------------------------------------------*/
	/*----------------     Neural Network Fields    ----------------*/
	/*--------------------------------------------------------------*/

	/** Master neural network model (copied to each thread) */
	private CellNet net0=null;
	/** Whether to use neural network for variant filtering */
	private boolean useNet=false;
	/** Score threshold for neural network filtering */
	private float netCutoff=0.5f;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** List of input file format objects */
	private final ArrayList<FileFormat> ffin=new ArrayList<FileFormat>();
	
	/** Output file format object */
	private final FileFormat ffout;

	/** Variant filtering parameters and methods */
	public final VarFilter varFilter=new VarFilter();
	/** SAM/BAM filtering parameters and methods */
	public final SamFilter samFilter=new SamFilter();
	/** Histogram arrays for variant scores by type [type][score] */
	public final long[][] scoreArray=new long[8][200];
	/** Histogram array for ploidy/zygosity distribution */
	public final long[] ploidyArray;
	/** Histogram arrays for average base quality by type [type][quality] */
	public final long[][] avgQualityArray=new long[8][100];
	/** Histogram array for maximum base quality distribution */
	public final long[] maxQualityArray=new long[100];
	/** Allele and reference depth arrays [depth_type][variant_type] */
	public final long[][] ADArray=new long[2][7];
	/** Allele frequency totals by variant type */
	public final double[] AFArray=new double[7];
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Maximum size for thread-local variant maps before dumping */
	private static int vmtSizeLimit=10000;
	
	/** Whether to call variants at N bases in reference */
	static boolean callNs=false;
	/** Whether to trim whitespace from read names and comments */
	static boolean trimWhitespace=true;
	/** Whether to attempt fixing of indel variants in reads */
	public static boolean fixIndels=true;
	
	/** Enable high-performance SAM streaming */
	static boolean useStreamer=true;
	/** Enable multi-file streaming for large datasets */
	static boolean useStreamerMF=true;
	/** Number of threads for streaming operations */
	static int streamerThreads=SamStreamer.DEFAULT_THREADS;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages during processing */
	public static boolean verbose=false;
	/** True if an error was encountered during processing */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
