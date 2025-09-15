package synth;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import align2.QualityTools;
import dna.AminoAcid;
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
import structures.ByteBuilder;
import structures.ListNum;
import tax.TaxTree;
import tracker.ReadStats;

/**
*Generates synthetic metagenomic sequencing reads from reference genomes with realistic 
*coverage patterns, error profiles, and taxonomic diversity. Supports multiple sequencing 
*platforms (Illumina, PacBio, ONT) with platform-specific read length distributions, 
*error rates, and quality score modeling.
*
*Features include configurable coverage depth distributions (uniform, exponential, power-law),
*sine wave coverage modeling for realistic spatial bias, platform-specific error modeling
*with substitutions, indels, and homopolymer errors, taxonomic integration via TaxTree
*for proper species labeling, paired-end and single-end read generation with insert size
*variation, and multi-threaded processing for large reference datasets.
*
*Designed for benchmarking metagenomic analysis tools, generating training datasets,
*and validating computational pipelines with ground-truth synthetic data.
*
*@author Brian Bushnell
*@contributor Isla
*@contributor Janus
*@date Feb 8, 2025
 */
public class RandomReadsMG{

	/*
	 *TODO: Enhanced realism features (all optional, disabled by default)
	 *
	 *1. Adapter contamination for short inserts
	 *  -Tested!  Paired reads only.
	 *
	 *2. GC bias modeling  
	 *  -Implement coverage bias based on local GC content
	 *  -Bell curve with peak around 40-50% GC (typical Illumina bias)
	 *  -Parameter to control bias strength/shape
	 *
	 *3. Hexamer priming bias
	 *  -Model non-random "random" hexamer priming 
	 *  -Sequence-dependent coverage bias based on priming efficiency
	 *  -Use empirical hexamer preference data
	 *
	 *4. PCR duplicate simulation
	 *  -Tested for single-ended reads.
	 *
	 *5. Match string output option
	 *  -Add synthetic variant positions to read headers as match strings
	 *  -May significantly increase header length
	 *  -Useful for validation and benchmarking studies
	 */

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 *Program entry point for command-line execution. Instantiates RandomReadsMG,
	 *processes all input reference genomes to generate synthetic metagenomic reads,
	 *and reports timing statistics. The main workflow includes argument parsing,
	 *file validation, platform configuration, and multi-threaded read generation.
	 *
	 *@param args Command line arguments specifying input files, coverage parameters,
	 *            platform settings, and output destinations
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		RandomReadsMG x=new RandomReadsMG(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 *Constructs a RandomReadsMG instance with the specified command line arguments.
	 *Parses all parameters, validates file paths, configures platform-specific defaults,
	 *initializes output file formats, and loads the taxonomic tree for species labeling.
	 *The constructor performs complete setup but does not begin read generation.
	 *
	 *@param args Command line arguments containing input files, coverage settings,
	 *            platform parameters, and output file specifications
	 */
	public RandomReadsMG(String[] args){

		{ //Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;

		{ //Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();

			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			extout=parser.extout;
		}

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 

		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 *Parses command line arguments into internal configuration parameters.
	 *Handles file specifications, coverage depth settings, platform selection,
	 *error rate configuration, read length parameters, and custom depth mappings.
	 *Supports both explicit parameter=value syntax and positional file arguments.
	 *
	 *@param args Array of command line arguments to parse
	 *@return Configured Parser object with standard parameter settings
	 */
	private Parser parse(String[] args){

		//Create a parser object
		Parser parser=new Parser();

		//Set any necessary Parser defaults here
		//parser.foo=bar;

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("tree")){
				taxTreeFile=b;
			}else if(a.equals("in") || a.equals("ref")){
				Tools.getFileOrFiles(b, inputFiles, true, false, false, false);
			}else if(a.equals("depth") || a.equals("cov")){
				minDepth=maxDepth=Float.parseFloat(b);
			}else if(a.equals("mindepth") || a.equals("mincov")){
				minDepth=Float.parseFloat(b);
			}else if(a.equals("maxdepth") || a.equals("maxcov")){
				maxDepth=Float.parseFloat(b);
			}else if(a.equals("depthvariance") || a.equals("variance")){
				depthVariance=Float.parseFloat(b);
			}else if(a.equals("wavecov") || a.equals("sinewave")){
				waveCoverage=Parse.parseBoolean(b);
			}else if(a.equals("waves") || a.equals("sinewaves")){
				numSineWaves=Integer.parseInt(b);
				waveCoverage=numSineWaves>0;
			}else if(a.equals("maxamplitude") || a.equals("waveamp")){
				waveAmp=Float.parseFloat(b);
			}else if(a.equals("oribias")){
				oriBias=Float.parseFloat(b);
			}else if(a.equals("minprob") || a.equals("minwaveprob")){
				minWaveProb=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("minPeriod")){
				minPeriod=Parse.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("maxPeriod")){
				maxPeriod=Parse.parseIntKMG(b);

			}else if(a.equals("mode") || a.equals("depthmode")){
				depthMode=Tools.find(b.toUpperCase(), modes);
				assert(depthMode>=0) : depthMode;
			}else if(a.equals("platform")){
				String upper=b.toUpperCase();
				if(upper.equals("PB")){upper="PACBIO";}
				platform=Tools.find(upper, platforms);
				assert(platform>=0) : platform;
			}else if(a.equals("insert") || a.equals("avginsert")){
				avgInsert=Float.parseFloat(b);
			}else if(a.equals("paired") || a.equals("int") || a.equals("interleaved")){
				paired=Parse.parseBoolean(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("reads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("loud")){
				loud=Parse.parseBoolean(b);
			}else if(a.equals("hexamer") || a.equals("kprime") || a.equals("kmerprime") || 
					a.equalsIgnoreCase("randomPriming") || a.equalsIgnoreCase("randomkmer")){
				if(Tools.isNumeric(b)){
					randomPriming=true;
					RandomHexamer.setK(Integer.parseInt(b));
				}else{randomPriming=Parse.parseBoolean(b);}
			}else if(a.equals("minkmerprob") || a.equals("minkprob")){
				RandomHexamer.setMinProb(Float.parseFloat(b));
			}else if(a.equals("kmerpower") || a.equals("kpower") || a.equals("kmerexp") || a.equals("kexp")){
				RandomHexamer.setPower(Float.parseFloat(b));
			}

			else if(a.equalsIgnoreCase("addErrors")){
				addErrors=Parse.parseBoolean(b);
			}else if(a.equals("qscore") || a.equals("avgq") || a.equals("qavg") || a.equals("avgqual")){
				meanQScore=Integer.parseInt(b);
			}else if(a.equals("qrange")){
				qScoreRange=Integer.parseInt(b);
			}else if(a.equals("subrate") || a.equals("snprate")){
				subRate=Float.parseFloat(b);
			}else if(a.equals("indelrate")){
				indelRate=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("pacBioLengthSigma") || a.equals("pbsigma") || a.equals("pacbiosigma")){
				pacBioLengthSigma=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("ontLongTailFactor") || a.equals("tailfactor")){
				ontLongTailFactor=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("srate")){
				sRate=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("irate")){
				iRate=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("drate")){
				dRate=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("hrate")){
				hRate=Float.parseFloat(b);
			}

			else if(a.equals("len") || a.equals("length") || a.equals("readlen") || 
					a.equals("readlength") || a.equals("meanlen") || a.equals("avglen")){
				readlen=meanLength=Parse.parseIntKMG(b);
				setReadLength=true;
			}else if(a.equals("maxlen") || a.equals("maxlength")){
				maxLength=Parse.parseIntKMG(b);
				setMaxLength=true;
			}else if(a.equals("minlen") || a.equals("minlength")){
				minLength=Parse.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("maxConcurrentGenomes") || a.equals("concurrency") || 
					a.equalsIgnoreCase("mcg")){
				maxConcurrentGenomes=Tools.max(1, Integer.parseInt(b));
			}else if(a.equalsIgnoreCase("singleFileThreads") || a.equals("sfthreads") || 
					a.equalsIgnoreCase("sft")){
				singleFileThreads=Integer.parseInt(b);
				Shared.setThreads(Math.max(singleFileThreads, Shared.threads()));
			}else if(a.equalsIgnoreCase("threads") || a.equalsIgnoreCase("t")){
				singleFileThreads=Integer.parseInt(b);
				Shared.setThreads(singleFileThreads);
			}else if(a.equals("fragadapter") || a.equals("fragadapter1") || 
					a.equals("adapter") || a.equals("adapter1")){
				fragadapter1=(b==null ? null : b.getBytes());
			}else if(a.equals("fragadapter2") || a.equals("adapter2")){
				fragadapter2=(b==null ? null : b.getBytes());
			}else if(a.equals("addadapters")){
				addAdapters=Parse.parseBoolean(b);
			}else if(a.equals("pcr") || a.equals("pcrrate")){
				pcrRate=Float.parseFloat(b);
			}

			else if(b!=null && (a.startsWith("cov_") || a.startsWith("depth_"))){
				float f=Float.parseFloat(b);
				String name=split[0].substring(a.indexOf('_')+1);
				System.err.println("Setting custom depth "+f+" for "+name);
				depthMap.put(name, f);
			}else if(b!=null && Tools.isNumeric(b) && new File(split[0]).isFile()){
				float f=Float.parseFloat(b);
				inputFiles.add(split[0]);
				String name=ReadWrite.stripPath(split[0]);
				System.err.println("Setting custom depth "+f+" for "+name);
				depthMap.put(name, f);
			}else if(b==null && Tools.find(arg.toUpperCase(), modes)>=0){
				depthMode=Tools.find(arg.toUpperCase(), modes);
				assert(depthMode>=0) : depthMode;
			}else if(b==null && Tools.find(arg.toUpperCase(), platforms)>=0){
				platform=Tools.find(arg.toUpperCase(), platforms);
				assert(platform>=0) : platform;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){ //Parse standard flags in the parser
				//do nothing
			}else{
				File f=new File(arg);
				if(f.exists() && f.canRead()){
					Tools.getFileOrFiles(arg, inputFiles, true, false, false, false);
				}else{
					outstream.println("Unknown parameter "+args[i]);
					assert(false) : "Unknown parameter "+args[i];
				}
			}
		}

		return parser;
	}

	/**
	 *Replaces '#' placeholders in output file paths with '1' and '2' for paired-end output.
	 *Enables convenient specification of paired output files using a single path template.
	 *If out1 contains '#' and out2 is null, creates both out1 and out2 paths by replacing
	 *'#' with '1' and '2' respectively. Validates that out2 is not specified without out1.
	 */
	private void doPoundReplacement(){

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error-cannot define out2 without defining out1.");}
	}

	/**
	 *Validates that all input files exist and are readable, and all output files
	 *can be written without conflicts. Checks file permissions, prevents overwrites
	 *when overwrite=false, detects duplicate file specifications, and ensures
	 *input files are accessible before beginning processing.
	 *@throws RuntimeException if file validation fails
	 */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, inputFiles.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		assert(FastaReadInputStream.settingsOK());
	}

	/**
	 *Validates that all user-specified parameters are within acceptable ranges
	 *and required parameters are properly set. Performs bounds checking on depth
	 *ranges, error rates, read lengths, and platform-specific parameters to prevent
	 *invalid configurations that could cause runtime errors or unrealistic output.
	 *@return true if all parameters pass validation
	 */
	private boolean validateParams(){
		//		assert(false) : "TODO";
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 *Sets platform-specific default error rates and quality scores based on the selected
	 *sequencing technology. Illumina uses default rates (mostly error-free), while ONT
	 *and PacBio get realistic error profiles matching their characteristics.
	 */
	void setPlatformDefaults(){
		if(platform==ILLUMINA){return;}
		else if(platform==ONT){
			if(sRate<0){sRate=0.0025f;}
			if(iRate<0){iRate=0.0055f;}
			if(dRate<0){dRate=0.0045f;}
			if(hRate<0){hRate=0.02f;}
			Shared.FAKE_QUAL=20;
		}else if(platform==PACBIO){
			if(sRate<0){sRate=0.00015f;}
			if(iRate<0){iRate=0.000055f;}
			if(dRate<0){dRate=0.000045f;}
			if(hRate<0){hRate=0.000015f;}
			Shared.FAKE_QUAL=35;
		}else{
			throw new RuntimeException("Unknown Platform "+platform);
		}
	}

	/**
	 *Orchestrates the complete synthetic read generation process from input reference genomes.
	 *Sets platform-specific error rates, initializes random hexamer priming if enabled,
	 *creates output streams, spawns worker threads for parallel processing, and aggregates
	 *statistics. Handles thread synchronization, error state propagation, and resource cleanup.
	 *Reports detailed timing and throughput statistics upon completion.
	 *
	 *@param t Timer object for execution timing and reporting
	 */
	void process(Timer t){

		setPlatformDefaults();
		if(randomPriming){RandomHexamer.initialize(Shared.threadLocalRandom(seed));}

		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=true;

		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros=makeCros();

		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;

		//Process the reads in separate threads
		spawnThreads(inputFiles, ros);

		if(verbose){outstream.println("Finished; closing streams.");}

		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStream(ros);

		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;

		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(t.elapsed, readsOut, basesOut, 8));

		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 *Creates and initializes a ConcurrentReadInputStream for reading reference genome sequences.
	 *Configures the stream for multi-threaded access with automatic format detection and
	 *starts background reading threads. Used for processing individual reference files
	 *during synthetic read generation.
	 *
	 *@param ff FileFormat object specifying input file type and compression
	 *@return Initialized and started ConcurrentReadInputStream ready for reading
	 */
	private ConcurrentReadInputStream makeCris(FileFormat ff){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}

	/**
	 *Creates and initializes a ConcurrentReadOutputStream for writing generated synthetic reads.
	 *Configures paired-end or single-end output based on user settings, handles interleaved
	 *output when appropriate, and sets up quality score output streams. Returns null if
	 *no output file was specified (useful for benchmarking scenarios).
	 *@return Configured and started ConcurrentReadOutputStream, or null if no output specified
	 */
	private ConcurrentReadOutputStream makeCros(){
		if(ffout1==null){return null;}

		//Set output buffer size
		final int buff=4;

		//Notify user of output mode
		if(paired && out2==null){
			outstream.println("Writing interleaved.");
		}

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(
				ffout1, ffout2, qfout1, qfout2, buff, null, false);
		ros.start(); //Start the stream
		return ros;
	}

	/**
	 *Creates and launches worker threads for parallel processing of reference genome files.
	 *Handles special case of single-file multi-threading by duplicating file entries and
	 *adjusting depth calculations. Distributes work across available threads using atomic
	 *counters for load balancing, then waits for all threads to complete before returning.
	 *
	 *@param files Collection of reference genome file paths to process
	 *@param ros Output stream for generated synthetic reads
	 */
	private void spawnThreads(final Collection<String> files, final ConcurrentReadOutputStream ros){

		//Do anything necessary prior to processing
		ArrayList<String> flist=new ArrayList<String>(files);

		if(singleFileThreads>1 && Shared.threads()>1 && seed<0 && flist.size()==1){
			//This allows multithreaded processing of a single file.
			//Not really necessary though
			int t=singleFileThreads;
			for(int i=1; i<t; i++){flist.add(flist.get(0));}
			if(maxReads>0){maxReads=(maxReads+t-1)/t;}
			else{
				String name=ReadWrite.stripPath(flist.get(0));
				Float depth=depthMap.get(name);
				if(depth==null){
					depth=randomDepth(Shared.threadLocalRandom(seed))/t;
				}else{
					depth/=t;
				}
				depthMap.put(name, depth);
			}
			loud=false;
		}

		//Determine how many threads may be used
		final int threads=Tools.min(Shared.threads(), flist.size(), maxConcurrentGenomes);
		//		System.err.println("Using "+threads);
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		AtomicInteger atom=new AtomicInteger(0);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(flist, ros, i, atom));
		}

		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}

		//Wait for threads to finish
		waitForThreads(alpt);

		//Do anything necessary after processing

	}

	/**
	 *Blocks until all worker threads complete, then aggregates their statistics.
	 *Handles thread interruption gracefully and accumulates per-thread counters for
	 *reads processed, bases processed, reads generated, bases generated, and PCR duplicates.
	 *Sets global error state if any thread encountered errors during processing.
	 *@param alpt List of ProcessThread objects to wait for completion
	 */
	private void waitForThreads(ArrayList<ProcessThread> alpt){

		//Wait for completion of all threads
		boolean success=true;
		for(ProcessThread pt : alpt){

			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try{
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e){
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}

			//Accumulate per-thread statistics
			readsProcessed+=pt.readsInT;
			basesProcessed+=pt.basesInT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			pcrDupesOut+=pt.pcrDupesOutT;
			success&=pt.success;
		}

		//Track whether any threads failed
		if(!success){errorState=true;}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 *Creates a ConcurrentReadInputStream for a specific reference file.
	 *Automatically detects file format (FASTA/FASTQ), handles compression,
	 *and configures the stream for single-threaded reading with buffering.
	 *This overloaded version is used by ProcessThread workers for individual file access.
	 *
	 *@param fname Path to the reference genome file to read
	 *@return Initialized and started ConcurrentReadInputStream for the specified file
	 */
	private ConcurrentReadInputStream makeCris(String fname){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}

	/**
	 *Determines the coverage depth for a specific input file.
	 *First checks for custom depth settings, then uses the selected distribution mode.
	 *
	 *@param path Filename to check for custom depth setting
	 *@param taxID Taxonomy ID to check for custom depth setting
	 *@param fnum File number for reporting
	 *@param randy Random number generator
	 *@return The depth to use for this file
	 */
	float chooseDepthForFile(String path, int taxID, int fnum, Random randy){
		Float custom=null;
		String fname=ReadWrite.stripPath(path);
		if(depthMap.size()>0){
			custom=depthMap.get(path);
			if(taxID>0 && custom==null){custom=depthMap.get(Integer.toString(taxID));}
		}
		final float depth;
		if(custom!=null){depth=custom;}
		else{depth=randomDepth(randy);}
		if(loud){
			String dstring=(custom==null ? "" : " custom")+
					String.format("depth=%.2f", depth);
			String idstring=taxID>0 ? ("tid "+taxID) : ("name "+fname);
			System.err.println("File "+fnum+", "+idstring+": "+dstring);
		}
		return depth;
	}

	/**
	 *Generates a random coverage depth value based on the selected distribution mode.
	 *Delegates to the appropriate distribution function (MIN4, EXP, ROOT, or LINEAR)
	 *to create realistic metagenomic abundance patterns. Each mode produces different
	 *characteristics matching common metagenomic diversity patterns.
	 *
	 *@param randy Random number generator for consistent reproducible results
	 *@return Generated coverage depth value within the configured min/max range
	 */
	float randomDepth(Random randy){
		final float depth;
		if(depthMode==MIN4){depth=depthMin4(randy);}
		else if(depthMode==EXP){depth=depthExp(randy);}
		else if(depthMode==ROOT){depth=depthRoot(randy);}
		else if(depthMode==LINEAR){depth=depthLinear(randy);}
		else{throw new RuntimeException("Unknown mode "+depthMode);}
		return depth;
	}

	/**
	 *Generates a depth value using the minimum of 4 random values.
	 *This creates a distribution that is skewed toward lower values.
	 *@param randy Random number generator
	 *@return The generated depth value
	 */
	float depthMin4(Random randy){
		float minRoot=(float)Math.sqrt(minDepth);
		float range=(float)(Math.sqrt(maxDepth)-minRoot);
		final float rootDepth=minRoot+(Tools.min(randy.nextFloat(), randy.nextFloat(), 
				randy.nextFloat(), randy.nextFloat()))*range;
		final float fileDepth=rootDepth*rootDepth;
		return fileDepth;
	}

	/**
	 *Generates a depth value using a uniform linear distribution.
	 *@param randy Random number generator
	 *@return The generated depth value
	 */
	float depthLinear(Random randy){
		float range=(float)(maxDepth-minDepth);
		return randy.nextFloat()*range+minDepth;
	}

	/**
	 *Generates a depth value using a square root distribution.
	 *This creates a distribution that is moderately skewed toward lower values.
	 *@param randy Random number generator
	 *@return The generated depth value
	 */
	float depthRoot(Random randy){
		float range=(float)maxDepth-minDepth;
		float root=randy.nextFloat();
		return root*root*range+minDepth;
	}

	/**
	 *Generates a depth value using an exponential distribution.
	 *This creates a natural long-tailed distribution common in metagenomic samples.
	 *@param randy Random number generator
	 *@return The generated depth value
	 */
	float depthExp(Random randy){
		double lambda=1/Math.sqrt(minDepth*maxDepth);
		double depth=Tools.exponential(randy, lambda);
		while(depth<minDepth || depth>maxDepth){
			depth=Tools.exponential(randy, lambda);
		}
		return (float)depth;
	}

	/**
	 *Adds adapter contamination to a read starting at the specified position.
	 *Used when insert size is shorter than read length, causing sequencing
	 *to read through the insert and into the adapter sequence. Positions beyond
	 *the adapter are filled with 'G' bases (typical 2-dye Illumina behavior) with
	 *low quality scores.
	 *
	 *@param r The read to modify
	 *@param loc Position where adapter sequence begins
	 *@param adapter The adapter sequence to add
	 */
	public static void addFragAdapter(Read r, final int loc, final byte[] adapter){
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		final int initial=(bases==null ? 0 : bases.length);

		if(initial>0 && loc>=0 && loc<initial){
			final int lim=Tools.min(initial, adapter.length+loc);

			//Add adapter sequence
			for(int i=loc, j=0; i<lim; i++, j++){
				if(AminoAcid.isFullyDefined(bases[i])){
					bases[i]=adapter[j];
				}
			}

			//Fill remaining positions with G and low quality
			for(int i=lim; i<initial; i++){
				if(AminoAcid.isFullyDefined(bases[i])){
					bases[i]='G';
					if(quals!=null){quals[i]=(byte)Math.min(8, quals[i]);} // Low quality for off end
				}
			}
		}
	}

	/**
	 *Applies Illumina-specific sequencing errors based on quality scores. Generates
	 *random quality scores within the specified range, then introduces substitution
	 *errors probabilistically based on the quality score at each position. This models
	 *the relationship between base quality and error probability in Illumina data.
	 *
	 *@param r The read to modify
	 *@param meanQ Mean quality score for the read
	 *@param qRange Quality score range (±qRange around meanQ)
	 *@param randy Random number generator for randomized decisions
	 *@return The number of substitutions added to the read
	 */
	public static int mutateIllumina(Read r, int meanQ, int qRange, Random randy){
		if(r.quality==null){r.quality=new byte[r.length()];}
		final byte[] bases=r.bases, quals=r.quality;
		int fullRange=qRange*2+1;
		int baseQ=meanQ-qRange;
		int subs=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(AminoAcid.isFullyDefined(b)){
				int q=baseQ+randy.nextInt(fullRange);
				quals[i]=(byte)q;
				float prob=QualityTools.PROB_CORRECT[q];
				if(randy.nextFloat()>prob){
					int x=AminoAcid.baseToNumber[b];
					x=(x+(randy.nextInt(3)+1))&3;
					bases[i]=AminoAcid.numberToBase[x];
					subs++;
				}
			}else{
				quals[i]=0;
			}
		}
		return subs;
	}

	/**
	 *Adds substitution errors to a read at the specified rate, randomly changing bases
	 *without considering quality scores. Used for simple error modeling or in combination
	 *with platform-specific error functions.
	 *
	 *@param r The read to modify
	 *@param rate The probability (0.0-1.0) of introducing a substitution at each position
	 *@param randy Random number generator for randomized decisions
	 *@return The number of substitutions added to the read
	 */
	public static int addSubs(Read r, float rate, Random randy){
		final byte[] bases=r.bases;
		int subs=0;

		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(AminoAcid.isFullyDefined(b) && randy.nextFloat()<rate){
				//Make a substitution-choose a different base
				int x=AminoAcid.baseToNumber[b];
				x=(x+(randy.nextInt(3)+1))&3; // Add 1-3 to avoid same base
				bases[i]=AminoAcid.numberToBase[x];
				subs++;
			}
		}
		return subs;
	}

	/**
	 *Adds insertions and deletions to a read at the specified rate, ensuring the final read length
	 *matches the desired length.
	 *
	 *This method introduces indels randomly throughout the read while tracking available "padding"
	 *to ensure the final read is exactly the requested length. Insertions add random bases and
	 *increase available padding. Deletions are only performed when padding is available to prevent
	 *the read from becoming too short. Quality scores for inserted bases are randomly generated
	 *within the specified quality range.
	 *
	 *@param r The read to modify with indels
	 *@param rate The probability (0.0-1.0) of introducing an indel at each position
	 *@param desiredLength The target length for the modified read
	 *@param meanQ The mean quality score for inserted bases
	 *@param qRange The range around meanQ for quality score randomization
	 *@param randy Random number generator for randomized decisions
	 *@return The number of indels (insertions+deletions) added to the read
	 *
	 *@throws AssertionError If the resulting read length doesn't match the desired length
	 */
	public static int addIndels(Read r, float rate, int desiredLength, int meanQ, int qRange, Random randy){
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;

		int padding=r.bases.length-desiredLength;
		int indels=0;

		//Create arrays that can accommodate insertions
		ByteBuilder newBases=new ByteBuilder(desiredLength);
		ByteBuilder newQuals=quals==null ? null : new ByteBuilder(desiredLength);

		final int fullRange=qRange*2+1;
		final int baseQ=meanQ-qRange;

		for(int i=0; i<bases.length && newBases.length<desiredLength; i++){
			if(randy.nextFloat()<rate){
				//Make an indel-50% insertion, 50% deletion
				if(randy.nextBoolean()){
					//Insertion-add current base plus a random base
					newBases.append(bases[i]);
					if(newQuals!=null){newQuals.append(quals[i]);}

					//Insert a random base
					int x=randy.nextInt(4);
					newBases.append(AminoAcid.numberToBase[x]);
					if(newQuals!=null){
						int q=baseQ+randy.nextInt(fullRange);
						newQuals.append((byte)q);
					}
					indels++;
					padding++;
				}else if(padding>0){
					//Deletion-skip this base
					padding--;
					indels++;
					continue; // Skip to next base
				}
			}else{
				//No indel-keep base as is
				newBases.append(bases[i]);
				if(newQuals!=null){newQuals.append(quals[i]);}
			}
		}

		//Create right-sized result arrays
		assert(newBases.length==desiredLength && newBases.array.length==desiredLength);
		r.bases=newBases.array;
		if(newQuals!=null){r.quality=newQuals.array;}
		return indels;
	}

	/**
	 *Applies long-read sequencing errors including substitutions, insertions, deletions,
	 *and homopolymer-specific errors. Models the error characteristics of PacBio and
	 *ONT platforms with context-dependent error rates.
	 *
	 *@param r The read to modify
	 *@param sRate Substitution error rate
	 *@param iRate Insertion error rate  
	 *@param dRate Deletion error rate
	 *@param hRate Homopolymer error bonus rate
	 *@param randy Random number generator
	 *@return The number of changes made to the read
	 */
	public int mutateLongRead(Read r, float sRate, float iRate, float dRate, float hRate, Random randy){
		final float delProb=dRate/Math.max(0.000000000001f, (iRate+dRate));
		final float errProb=sRate+iRate+dRate;
		float bonus=0;
		byte prev=-1;
		int changes=0;
		ByteBuilder bb=new ByteBuilder(r.length()/8+10);
		final byte[] bases=r.bases;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			if(!AminoAcid.isFullyDefined(b)){bb.append(b);continue;}

			float f=randy.nextFloat();
			bonus=(b==prev ? bonus+hRate : 0);
			if(f>=errProb+bonus){
				bb.append(b);
			}else if(f<sRate){ //Substitution
				int x=AminoAcid.baseToNumber[b];
				x=(x+(randy.nextInt(3)+1))&3;
				bb.append(AminoAcid.numberToBase[x]);
				changes++;
			}else{ //Indel
				if(randy.nextFloat()<delProb){
					//Deletion, do nothing
				}else{ //Insertion
					bb.append(b);
					byte b2=(f>errProb ? b : AminoAcid.numberToBase[randy.nextInt(4)]);
					bb.append(b2); //This is a same-base insertion, sometimes.
				}
				bonus=0;
				changes++;
			}
		}
		r.bases=bb.toBytes();
		return changes;
	}

	/**
	 *Generates realistic PacBio HiFi read lengths using a log-normal distribution.
	 *PacBio HiFi reads follow a characteristic log-normal pattern with most reads
	 *clustering around the mean and a moderate tail of longer reads. The sigma parameter
	 *controls the distribution width, with typical values of 0.4-0.5 for HiFi data.
	 *Applies reasonable bounds to prevent extremely short or long reads.
	 *
	 *@param minLength Minimum allowable read length (hard lower bound)
	 *@param meanLength Target mean read length (typically 10000-15000 for HiFi)
	 *@param maxLength Maximum allowable read length (hard upper bound) 
	 *@param sigma Standard deviation in log space (0.4-0.5 typical for HiFi)
	 *@param randy Random number generator for consistent reproducible results
	 *@return Generated read length bounded within [minLength, max(meanLength*4, maxLength)]
	 */
	public static int generatePacBioHiFiLength(int minLength, int meanLength, int maxLength, double sigma, Random randy){
		//Log-normal distribution
		double mu=Math.log(meanLength)-0.5*sigma*sigma; //Bias correction for log-normal mean
		double logLength=mu+randy.nextGaussian()*sigma;
		int length=(int)Math.round(Math.exp(logLength));

		//Apply reasonable bounds
		return Math.max(minLength, Math.min(length, Math.max(meanLength*4, maxLength)));
	}

	/**
	 *Generates realistic ONT read lengths using a mixed distribution approach.
	 *ONT reads exhibit a bimodal pattern: most reads follow a log-normal core distribution,
	 *but a significant fraction form a heavy exponential tail of very long reads.
	 *The longTailFactor controls what fraction of reads come from the exponential tail
	 *versus the log-normal core, typically 0.1-0.3 for real ONT data.
	 *
	 *@param minLength Minimum allowable read length (hard lower bound)
	 *@param meanLength Approximate target mean for the core distribution
	 *@param maxLength Maximum allowable read length (hard upper bound)
	 *@param longTailFactor Probability of drawing from heavy tail (0.1-0.3 typical)
	 *@param randy Random number generator for consistent reproducible results
	 *@return Generated read length following ONT-like mixed distribution
	 */
	public static int generateONTLength(int minLength, int meanLength, int maxLength, double longTailFactor, Random randy){
		//Use mixed distribution approach
		if(randy.nextDouble()<longTailFactor){
			//Generate from heavy tail (exponential)
			double scale=meanLength*2; //Scale factor produces longer tail reads
			double x=-Math.log(randy.nextDouble())*scale;
			return (int)Math.min(x, maxLength);
		}else{
			//Generate from log-normal core distribution
			double sigma=0.5;
			double mu=Math.log(meanLength)-0.5*sigma*sigma;
			double logLength=mu+randy.nextGaussian()*sigma;
			return Math.max(minLength, (int)Math.min(Math.exp(logLength), maxLength));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 *Worker thread that processes reference genome files to generate synthetic reads.
	 *Each thread maintains its own random number generator for reproducible results
	 *and processes files from the shared queue using atomic coordination. Accumulates
	 *statistics locally and handles all aspects of read generation including coverage
	 *modeling, error introduction, and output formatting.
	 */
	private class ProcessThread extends Thread{

		/**
		 *Constructs a ProcessThread worker for parallel synthetic read generation.
		 *@param files_ Shared list of reference genome files to process
		 *@param ros_ Output stream for generated synthetic reads
		 *@param tid_ Thread ID for seeding random number generator
		 *@param nextFile_ Atomic counter for work distribution coordination
		 */
		ProcessThread(final ArrayList<String> files_, final ConcurrentReadOutputStream ros_, 
				int tid_, final AtomicInteger nextFile_){
			files=files_;
			ros=ros_;
			tid=tid_;
			nextFile=nextFile_;
		}

		/**
		 *Main thread execution method that processes assigned reference files.
		 *Initializes thread-local random number generator, processes files from
		 *the shared queue using atomic coordination, and sets success flag upon
		 *completion. All statistics are accumulated in thread-local variables.
		 */
		@Override
		public void run(){
			//Initialize thread-local random generator with deterministic seed
			randy=Shared.threadLocalRandom(seed>=0 ? seed+tid : -1);

			//Process files using atomic work distribution
			for(int i=nextFile.getAndIncrement(); i<files.size(); i=nextFile.getAndIncrement()){
				String fname=files.get(i);
				processFile(fname, i);
			}

			//Mark successful completion
			success=true;
		}

		/**
		 *Processes a single input file, generating synthetic reads at the appropriate depth.
		 *@param path File name to process
		 *@param fnum File number for tracking and reporting
		 */
		void processFile(String path, int fnum){
			//			System.err.println("Thread "+tid+" processing file "+fnum+"; next="+nextFile);
			final String fname=ReadWrite.stripPath(path);
			ConcurrentReadInputStream cris=makeCris(path);

			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			int taxID=TaxTree.parseHeaderStatic2(path, tree);
			if(taxID<0 && ln.size()>0){
				Read c0=ln.get(0);
				taxID=TaxTree.parseHeaderStatic2(c0.id, tree);
			}
			final float fileDepth=(maxReads>0 ? 0 : chooseDepthForFile(fname, taxID, fnum, randy));
			//			assert(taxID>0) : "Can't parse taxID from "+fname;

			if(waveCoverage){
				covModel=new CoverageModel(numSineWaves, waveAmp, oriBias, minWaveProb, randy, minPeriod, maxPeriod);
			}else{covModel=null;}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){

				for(Read c : ln){
					float depthRatio=1f;
					if(varyDepthPerContig){
						depthRatio=1f+(depthVariance*(randy.nextFloat()+randy.nextFloat()))-depthVariance;
					}
					float contigDepth=depthRatio*fileDepth;
					//					System.err.println("depthRatio="+depthRatio+"; contigDepth="+contigDepth);
					processContig(c, contigDepth, taxID, fnum, fname);
				}

				//Notify the input stream that the list was used
				cris.returnList(ln);

				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}

		}

		/**
		 *Generates synthetic reads from a single contig at the specified depth.
		 *@param contig Source contig to generate reads from
		 *@param depth Coverage depth to generate
		 *@param taxID Taxonomy ID for read headers
		 *@param fnum File number for read headers
		 */
		private void processContig(Read contig, float depth, int taxID, int fnum, String fname){
			final int basesPerRead=(paired ? 2*readlen : readlen);
			readsInT++;
			basesInT+=contig.length();
			if(platform==ILLUMINA){
				if(contig.length()<basesPerRead+10){return;}
				if(paired && contig.length()<avgInsert){return;}
			}else{
				// Long-read platforms (PacBio/ONT)
				int minContigLength=minLength+50; // Buffer for read placement
				if(contig.length()<minContigLength){return;}
			}

			lastStart=lastInsert=lastStrand-1; // Reset PCR duplicates
			long basesToGenerate=(maxReads>0 ? maxReads*basesPerRead : (long)(depth*contig.length()));
			long readsGenerated=0;
			long basesGenerated=0;
			ArrayList<Read> list=new ArrayList<Read>(200);
			float variance=varyDepthPerContig ? 0 : randy.nextFloat()*depthVariance;

			//			System.err.println("Generating "+basesToGenerate+" for depth-"+depth+" contig "+contig.id);

			for(long i=0; basesGenerated<basesToGenerate; i++){
				Read r=generateRead(contig, i, taxID, fnum, contig.numericID, variance, fname);
				if(r!=null){
					list.add(r);
					readsGenerated+=r.pairCount();
					basesGenerated+=r.pairLength();
				}
				if(list.size()>=200){
					if(ros!=null){ros.add(list, 0);}
					list=new ArrayList<Read>(200);
				}
			}
			if(list.size()>0){if(ros!=null){ros.add(list, 0);}}
			//			System.err.println("Generated "+basesGenerated+" for depth-"+depth+" contig "+contig.id);

			readsOutT+=readsGenerated;
			basesOutT+=basesGenerated;
		}

		/**
		 *Generates a single or paired read based on the paired flag.
		 *
		 *@param contig Source contig to generate reads from
		 *@param rnum Read number for identification
		 *@param taxID Taxonomy ID for read headers
		 *@param fnum File number for read headers
		 *@param cnum Contig number for read headers
		 *@param variance Variance value for depth variation
		 *@return The generated read or null if skipped
		 */
		private Read generateRead(Read contig, long rnum, int taxID, 
				int fnum, long cnum, float variance, String fname){
			if(platform==ILLUMINA){
				if(paired){return generateReadPair(contig, rnum, taxID, fnum, cnum, variance, fname);}
				else{return generateReadSingle(contig, rnum, taxID, fnum, cnum, variance, fname);}
			}else if(platform==PACBIO){return generateLongRead(contig, rnum, taxID, fnum, cnum, variance, fname);}
			else if(platform==ONT){return generateLongRead(contig, rnum, taxID, fnum, cnum, variance, fname);}
			else{
				throw new RuntimeException("Unknown Platform "+platform);
			}
		}

		/**
		 *Generates a single unpaired read from a contig.
		 *
		 *@param contig Source contig to generate reads from
		 *@param rnum Read number for identification
		 *@param taxID Taxonomy ID for read headers
		 *@param fnum File number for read headers
		 *@param cnum Contig number for read headers
		 *@param variance Variance value for depth variation
		 *@return The generated read or null if skipped
		 */
		private Read generateLongRead(Read contig, long rnum, int taxID, 
				int fnum, long cnum, float variance, String fname){
			final boolean novel=(pcrRate<=0 || lastInsert<1 || randy.nextFloat()>=pcrRate);
			int insert=lastInsert, start=lastStart, strand=lastStrand;

			if(novel){
				if(platform==PACBIO){
					insert=generatePacBioHiFiLength(minLength, meanLength, maxLength, pacBioLengthSigma, randy);
				}else{
					insert=generateONTLength(minLength, meanLength, maxLength, ontLongTailFactor, randy);
				}
				if(insert>=contig.length()){return null;}
				start=randy.nextInt(contig.length()-insert);
				strand=randy.nextInt(2);

				// Cache for potential duplicates
				lastInsert=insert;
				lastStart=start;
				lastStrand=strand;
			}else{pcrDupesOutT++;}

			int paddedLen=insert+(indelRate>0 ? 20 : 0);
			if(skip((start+paddedLen)/2, contig.length(), variance)
					|| start+paddedLen>contig.length() || start<0){return null;}

			byte[] bases=Arrays.copyOfRange(contig.bases, start, start+paddedLen);
			if(strand==1){AminoAcid.reverseComplementBasesInPlace(bases);}
			if(randomPriming && !RandomHexamer.keep(bases, randy)){return null;}
			String header=makeHeader(start, strand, paddedLen, taxID, fnum, cnum, 0, novel?0:1, fname);
			Read r=new Read(bases, null, header, rnum);
			if(addErrors){mutateLongRead(r, sRate, iRate, dRate, hRate, randy);}
			if(subRate>0){addSubs(r, subRate, randy);}
			if(indelRate>0){addIndels(r, indelRate, paddedLen, meanQScore, qScoreRange, randy);}
			return r;
		}

		/**
		 *Generates a single unpaired read from a contig.
		 *
		 *@param contig Source contig to generate reads from
		 *@param rnum Read number for identification
		 *@param taxID Taxonomy ID for read headers
		 *@param fnum File number for read headers
		 *@param cnum Contig number for read headers
		 *@param variance Variance value for depth variation
		 *@return The generated read or null if skipped
		 */
		private Read generateReadSingle(Read contig, long rnum, int taxID, 
				int fnum, long cnum, float variance, String fname){
			final int paddedLen=readlen+(indelRate>0 ? 5 : 0);
			final boolean novel=(pcrRate<=0 || lastStart<0 || randy.nextFloat()>=pcrRate);
			int start=lastStart, strand=lastStrand;

			if(novel){
				if(paddedLen>=contig.length()){return null;}
				start=randy.nextInt(contig.length()-paddedLen);
				strand=randy.nextInt(2);

				// Cache for potential duplicates
				lastStart=start;
				lastStrand=strand;
			}else{pcrDupesOutT++;}
			int insert=paddedLen;
			if(skip((start+insert)/2, contig.length(), variance) 
					|| start+paddedLen>=contig.length() || start<0){return null;}
			byte[] bases=Arrays.copyOfRange(contig.bases, start, start+paddedLen);
			if(strand==1){AminoAcid.reverseComplementBasesInPlace(bases);}
			if(randomPriming && !RandomHexamer.keep(bases, randy)){return null;}
			String header=makeHeader(start, strand, insert, taxID, fnum, cnum, 0, novel?0:1, fname);
			Read r=new Read(bases, null, header, rnum);
			if(addErrors){mutateIllumina(r, meanQScore, qScoreRange, randy);}
			if(subRate>0){addSubs(r, subRate, randy);}
			if(indelRate>0){addIndels(r, indelRate, readlen, meanQScore, qScoreRange, randy);}
			return r;
		}

		/**
		 *Generates a paired read from a contig with appropriate insert size and optional
		 *adapter contamination for short inserts.
		 *
		 *@param contig Source contig to generate reads from
		 *@param rnum Read number for identification
		 *@param taxID Taxonomy ID for read headers
		 *@param fnum File number for read headers
		 *@param cnum Contig number for read headers
		 *@param variance Variance value for depth variation
		 *@return The generated read pair or null if skipped
		 */
		private Read generateReadPair(Read contig, long rnum, int taxID, 
				int fnum, long cnum, float variance, String fname){
			final int paddedLen=readlen+(indelRate>0 ? 5 : 0);
			final boolean novel=(pcrRate<=0 || lastInsert<1 || randy.nextFloat()>=pcrRate);
			int insert=lastInsert, start1=lastStart, strand=lastStrand;
			if(novel){
				insert=-1;
				for(int i=0; i<9 && ((insert<paddedLen && !addAdapters) || insert<=0); i++){
					double g=randy.nextGaussian()*0.25f;
					insert=(int)((1+g)*avgInsert);
				}
				if(Math.max(insert,paddedLen)>=contig.length()){return null;}
				start1=randy.nextInt(contig.length()-Math.max(insert, paddedLen));
				strand=randy.nextInt(2);

				// Cache for potential duplicates
				lastInsert=insert;
				lastStart=start1;
				lastStrand=strand;
			}else{pcrDupesOutT++;}
			final int start2=start1+insert-paddedLen;
			if(skip((start1+insert)/2, contig.length(), variance)
					|| start1+paddedLen>=contig.length() || start1<0
					|| start2+paddedLen>=contig.length() || start2<0){return null;}

			byte[] bases1=Arrays.copyOfRange(contig.bases, start1, start1+paddedLen);
			byte[] bases2=Arrays.copyOfRange(contig.bases, start2, start2+paddedLen);
			AminoAcid.reverseComplementBasesInPlace(bases2);
			if(strand==1){
				byte[] temp=bases1;
				bases1=bases2;
				bases2=temp;
			}
			if(randomPriming && !RandomHexamer.keep(bases1, randy)){return null;}
			String header1=makeHeader(start1, strand, insert, taxID, fnum, cnum, 0, novel?0:1, fname);
			String header2=makeHeader(start1, strand, insert, taxID, fnum, cnum, 1, novel?0:1, fname);
			Read r1=new Read(bases1, null, header1, rnum);
			Read r2=new Read(bases2, null, header2, rnum);
			r2.setPairnum(1);
			r1.mate=r2;
			r2.mate=r1;

			// Add adapter contamination before other mutations
			if(addAdapters && insert<readlen){
				addFragAdapter(r1, insert, fragadapter1);
				addFragAdapter(r2, insert, fragadapter2);
			}

			if(addErrors){
				mutateIllumina(r1, meanQScore, qScoreRange, randy);
				mutateIllumina(r2, meanQScore, qScoreRange, randy);
			}
			if(subRate>0){
				addSubs(r1, subRate, randy);
				addSubs(r2, subRate, randy);
			}
			if(indelRate>0){
				addIndels(r1, indelRate, readlen, meanQScore, qScoreRange, randy);
				addIndels(r2, indelRate, readlen, meanQScore, qScoreRange, randy);
			}
			return r1;
		}

		/**
		 *Determines whether to skip generating a read at a particular position.
		 *This is used to create within-contig coverage variation.
		 *
		 *@param midpoint Position in the contig
		 *@param clen Length of the contig
		 *@param variance Variance parameter controlling the skip probability
		 *@return True if this position should be skipped
		 */
		private boolean skip(int midpoint, int clen, float variance){
			if(covModel!=null){
				return !covModel.shouldGenerateReadAt(midpoint, clen, randy);
			}else if(variance<=0){return false;}

			//Linear model only
			float maxSkipProb=1f-1f/(1f+variance);
			float relativePosition=midpoint/(float)clen;
			float skipProb=relativePosition*maxSkipProb;
			return randy.nextFloat()<skipProb;
		}

		/**
		 *Creates a read header with genome source information encoded.
		 *Format: f_[filenum]_c_[contignum]_s_[strand]_p_[position]_i_[insert]_d_[duplicate]_tid_[taxid] [pairnum]:
		 *Some fields are optional.
		 *@param start Starting position in the contig
		 *@param strand Read strand (0=forward, 1=reverse)
		 *@param insert Insert size for paired reads
		 *@param taxID Taxonomy ID
		 *@param fnum File number
		 *@param cnum Contig number
		 *@param pnum Pair number (0=first, 1=second)
		 *@param pcr PCR duplicate (0=original, 1+=duplicate)
		 *@return Formatted header string
		 */
		private String makeHeader(int start, int strand, int insert, int taxID, 
				int fnum, long cnum, int pnum, int pcr, String fname){
			bb.clear().append('f').under().append(fnum).under().append('c').under().append(cnum);
			bb.under().append('s').under().append(strand).under().append('p').under().append(start);
			bb.under().append('i').under().append(insert);
			if(pcrRate>0){bb.under().append('d').under().append(pcr);}
			if(taxID>0){bb.under().append("tid").under().append(taxID);}
			else{bb.under().append("name").under().append(fname);}
			bb.space().append(pnum+1).colon();
			return bb.toString();
		}

		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		/** Number of pcr duplicates produced by this thread */
		protected long pcrDupesOutT=0;
		/** Number of input contigs processed by this thread */
		protected long readsInT=0;
		/** Number of input bases processed by this thread */
		protected long basesInT=0;
		/** True only if this thread has completed successfully */
		boolean success=false;
		/** For generating PCR duplicates */
		private int lastStart=-1, lastInsert=-1, lastStrand=-1;
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Thread ID */
		final int tid;
		/** Atomic counter for coordinating file assignment across threads */
		final AtomicInteger nextFile;
		/** List of input reference files for this thread to process */
		private final ArrayList<String> files;
		/** Thread-local random number generator for reproducible results */
		private Random randy;
		/** Coverage model for sine wave spatial bias simulation */
		private CoverageModel covModel;
		private ByteBuilder bb=new ByteBuilder(128);
	}
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Set of input reference genome files to process */
	private LinkedHashSet<String> inputFiles=new LinkedHashSet<String>();
	/** Primary output file path for generated reads */
	private String out1=null;
	/** Secondary output file path for paired reads */
	private String out2=null;
	/** Primary quality output file path */
	private String qfout1=null;
	/** Secondary quality output file path */
	private String qfout2=null;
	/** Override output file extension */
	private String extout=null;
	/** Path to taxonomic tree file for species classification */
	private String taxTreeFile=null;
	/*--------------------------------------------------------------*/
	/** Number of reference sequences processed */
	protected long readsProcessed=0;
	/** Number of reference bases processed */
	protected long basesProcessed=0;
	/** Number of synthetic reads generated */
	protected long readsOut=0;
	/** Number of synthetic bases generated */
	protected long basesOut=0;
	/** Number of PCR duplicates generated */
	protected long pcrDupesOut=0;
	/** Enable verbose progress reporting */
	private boolean loud=true;
	/** Atomic counter for unique read ID assignment */
	private AtomicLong nextReadID=new AtomicLong(0);
	/** Minimum coverage depth for any reference sequence */
	private float minDepth=1;
	/** Maximum coverage depth for any reference sequence */
	private float maxDepth=256;
	/** Amount of within-contig coverage variation (0=uniform, 1=high variation) */
	private float depthVariance=0.5f;
	/** Enable sine wave coverage modeling for realistic spatial bias */
	private boolean waveCoverage=false;
	/** Number of overlapping sine waves for coverage modeling */
	private int numSineWaves=4;
	/** Maximum amplitude multiplier for sine wave coverage variation */
	private float waveAmp=0.7f;
	/** Orientation bias factor for coverage modeling */
	private float oriBias=0.25f;
	/** Minimum probability threshold for wave-based coverage */
	private float minWaveProb=0.1f;
	/** Minimum period length for sine wave coverage patterns */
	private int minPeriod=2000;
	/** Maximum period length for sine wave coverage patterns */
	private int maxPeriod=80000;
	/** Random seed for reproducible generation (-1 for random) */
	private long seed=-1;
	/** Enable per-contig depth variation within files */
	private boolean varyDepthPerContig=false;
	/** Custom depth settings for specific files or taxonomy IDs */
	private HashMap<String, Float> depthMap=new HashMap<String, Float>();
	/** Maximum number of reads to generate (-1 for depth-based) */
	private long maxReads=-1;
	/** Rate of substitution errors to add (0.0-1.0) */
	private float subRate=0;
	/** Rate of indel errors to add (0.0-1.0) */
	private float indelRate=0;
	/** Average insert size for paired-end reads */
	private float avgInsert=300;
	/** Length of individual reads to generate */
	private int readlen=150;
	/** Generate paired-end reads instead of single-end */
	private boolean paired=true;
	/** Add platform-specific sequencing errors */
	private boolean addErrors=false;
	/** Mean quality score for generated bases */
	private int meanQScore=25;
	/** Quality score range around the mean (±qScoreRange) */
	private int qScoreRange=0;
	/** Standard deviation for PacBio read length log-normal distribution */
	private float pacBioLengthSigma=0.5f;
	/** Long tail factor for ONT read length distribution */
	private float ontLongTailFactor=0.2f;
	/** Minimum read length for long-read platforms */
	private int minLength=1000;
	/** Mean read length for long-read platforms */
	private int meanLength=15000;
	/** Maximum read length for long-read platforms */
	private int maxLength=100000;
	/** Number of threads to use for single file processing */
	private int singleFileThreads=1;
	/** Maximum number of concurrent genomes to process */
	private int maxConcurrentGenomes=8192;
	/** Substitution error rate for long-read platforms (-1 for default) */
	private float sRate=-1;
	/** Insertion error rate for long-read platforms (-1 for default) */
	private float iRate=-1;
	/** Deletion error rate for long-read platforms (-1 for default) */
	private float dRate=-1;
	/** Homopolymer error bonus rate for long-read platforms (-1 for default) */
	private float hRate=-1;
	/** Available coverage depth distribution modes */
	static final String[] modes={"MIN4", "EXP", "ROOT", "LINEAR"};
	static final int MIN4=0, EXP=1, ROOT=2, LINEAR=3;
	/** Selected coverage depth distribution mode */
	int depthMode=MIN4;
	/** Available sequencing platform types */
	static final String[] platforms={"ILLUMINA", "ONT", "PACBIO"};
	static final int ILLUMINA=0, ONT=1, PACBIO=2;
	/** Selected sequencing platform */
	int platform=ILLUMINA;
	/** Flag indicating read length was explicitly set by user */
	boolean setReadLength=false;
	/** Flag indicating maximum length was explicitly set by user */
	boolean setMaxLength=false;
	/** Rate of PCR duplicates */
	private float pcrRate=0.0f;
	/** Use random kmer priming */
	private boolean randomPriming=false;
	/** Enable adapter contamination for short inserts */
	private boolean addAdapters=false;
	/** Forward read adapter sequence for contamination */
	private byte[] fragadapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGC".getBytes();
	/** Reverse read adapter sequence for contamination */
	private byte[] fragadapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT".getBytes();
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	/** Taxonomic tree for species classification and labeling */
	private final TaxTree tree;
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
