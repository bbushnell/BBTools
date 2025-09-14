package aligner;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.SIMDAlignByte;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import structures.IntListHashMap;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Performs high-throughput indel-free alignments using seed-and-extend or brute force strategies.
 * Implements k-mer indexing with rolling hash for query preprocessing and reference indexing.
 * Uses SIMD vectorization (AVX2/SSE) for diagonal alignment when sequences are short enough.
 * Supports multithreaded processing with work-stealing for reference sequence batches.
 * Algorithm: 1) Load all queries into memory with optional k-mer indexing,
 * 2) Stream reference sequences from disk in batches, 3) For each reference,
 * build k-mer index or use brute force, 4) Align all queries against reference,
 * 5) Output SAM format alignments with proper CIGAR strings and mapping quality.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date June 2, 2025
 */
public class IndelFreeAligner implements Accumulator<IndelFreeAligner.ProcessThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main entry point that initializes timer, creates aligner instance, and processes alignments.
	 * Follows standard BBTools pattern: timer start, object creation, processing, stream cleanup.
	 * @param args Command line arguments including reference file, query files, and alignment parameters
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		IndelFreeAligner x=new IndelFreeAligner(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	* Constructor.
	* @param args Command line arguments
	*/
	public IndelFreeAligner(String[] args){

		{ //Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{ //Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();

			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			in1=parser.in1;
			in2=parser.in2;
			extin=parser.extin;

			out1=parser.out1;
			extout=parser.extout;
		}

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 

		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses command line arguments using key=value format and configures alignment parameters.
	 * Processes specialized arguments like ref, subs, k, hits, and delegates standard arguments
	 * to Parser class. Validates parameter ranges and logical consistency during parsing.
	 * @param args Command line arguments in key=value format
	 * @return Configured Parser object with file paths and standard BBTools settings
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
			}else if(a.equals("ref")){
				refFile=b;
			}else if(a.equals("subs") || a.equals("maxsubs")){
				maxSubs=Integer.parseInt(b);
			}else if(a.equals("hits") || a.equals("minhits") || a.equals("seedhits")){
				minSeedHits=Math.max(1, Integer.parseInt(b));
			}else if(a.equals("minprob") || a.equals("minhitsprob")){
				minHitsProb=Float.parseFloat(b);
			}else if(a.equals("maxclip") || a.equals("clip")){
				Query.maxClip=Tools.max(0, Float.parseFloat(b));
				assert(Query.maxClip<1 || Query.maxClip==(int)Query.maxClip);
			}else if(a.equals("index")){
				indexQueries=Parse.parseBoolean(b);
			}else if(a.equals("prescan")){
				prescan=Parse.parseBoolean(b);
			}else if(a.equals("seedmap") || a.equals("map")){
				useSeedMap=Parse.parseBoolean(b);
			}else if(a.equals("seedlist") || a.equals("list")){
				useSeedMap=!Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
				assert(k<16) : "0<=k<16 : "+k;
				indexQueries=(k>0);
			}else if(a.equals("qstep") || a.equals("step") || a.equals("qskip")){
				qStep=Integer.parseInt(b);
				assert(qStep>0);
			}else if(a.equals("mm")){
				midMaskLen=(Tools.isNumeric(b) ? Integer.parseInt(b) :
					Parse.parseBoolean(b) ? 1 : 0);
			}else if(a.equals("blacklist") || a.equals("banhomopolymers")){
				Query.blacklistRepeatLength=(Tools.isNumeric(b) ? Integer.parseInt(b) :
					Parse.parseBoolean(b) ? 1 : 0);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){ //Parse standard flags in the parser
				//do nothing
			}else if(parser.out1==null && b==null && FileFormat.isSamOrBamFile(arg)){
				parser.out1=arg;
			}else if(parser.in1==null && b==null && FileFormat.isFastqFile(arg) && new File(arg).isFile()){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		return parser;
	}

	/**
	 * Expands paired-end filename template using # symbol replacement.
	 * Converts "reads#.fq" to separate "reads1.fq" and "reads2.fq" files when
	 * original file with # doesn't exist. Standard BBTools convention for paired files.
	 */
	private void doPoundReplacement(){
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}

		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}

	/**
	 * Automatically detects and adds compression extensions (.gz, .bz2) to file paths.
	 * Uses Tools.fixExtension() to ensure proper handling of compressed input/output files.
	 */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
	}

	/**
	* Performs comprehensive alignment testing all possible positions.
	* Uses SIMD optimization when available and appropriate.
	* @param query Query sequence bases
	* @param ref Reference sequence bases
	* @param maxSubs Maximum allowed substitutions
	* @param maxClips Maximum allowed clipped bases
	* @return List of all alignment positions with ≤maxSubs substitutions, or null if none found
	*/
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/**
	 * Optimizes static file I/O settings based on available thread count.
	 * Enables ByteFile mode BF2 (multithreaded) when >2 threads available for better
	 * I/O performance. Validates FastaReadInputStream configuration consistency.
	 */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		assert(FastaReadInputStream.settingsOK());
	}

	/**
	 * Validates alignment parameter constraints and logical consistency.
	 * Ensures: k-mer length 1-15 when indexing enabled, hit probability ≤1.0,
	 * middle mask length < k-1, substitution count ≥0. Critical for preventing
	 * array bounds violations and algorithmic correctness.
	 * @return True if all parameters satisfy mathematical and algorithmic constraints
	 */
	private boolean validateParams(){
		assert((k>=1 && k<=15) || !indexQueries);
		assert(minHitsProb<=1);
		assert(midMaskLen<k-1);
		assert(maxSubs>=0);
		return ((k>=1 && k<=15) || !indexQueries) && 
				(minHitsProb<=1) && (midMaskLen<k-1) && (maxSubs>=0);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Orchestrates complete alignment workflow using producer-consumer pattern.
	 * Algorithm: 1) Load all queries into memory with k-mer preprocessing,
	 * 2) Configure read validation and buffer sizes for optimal parallelization,
	 * 3) Create concurrent input stream for reference sequences,
	 * 4) Spawn worker threads using ThreadWaiter pattern,
	 * 5) Accumulate statistics and handle cleanup.
	 * Uses temporary buffer size reduction to improve parallelization for small references.
	 * @param t Timer for measuring total processing time and throughput calculation
	 */
	void process(Timer t){

		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;

		SamLine.RNAME_AS_BYTES=false;

		Query.setMode(k, midMaskLen, indexQueries);
		final ArrayList<Query> queries=fetchQueries(ffin1, ffin2);

		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		long oldBD=Shared.bufferData();
		Shared.setBufferData(Math.min(40000, oldBD));

		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris(refFile);

		//Optionally create a read output stream
		final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout1);

		//Process the reads in separate threads
		spawnThreads(cris, bsw, queries);

		if(verbose){outstream.println("Finished; closing streams.");}

		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);

		if(bsw!=null){bsw.poisonAndWait();}

		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		Shared.setBufferData(oldBD);

		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		outstream.println(Tools.things("Alignments", alignmentCount, 8));
		outstream.println(Tools.things("Seed Hits", seedHitCount, 8));

		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, 
			final ByteStreamWriter bsw, final ArrayList<Query> qList){

		//Do anything necessary prior to processing

		//Determine how many threads may be used
		final int threads=Shared.threads();

		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, bsw, qList, maxSubs, i));
		}

		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;

		//Do anything necessary after processing

	}

	/**
	 * Thread-safe accumulation of processing statistics from completed worker threads.
	 * Called by ThreadWaiter framework when each worker thread finishes. Synchronizes
	 * on thread object to safely aggregate read counts, base counts, alignment counts,
	 * and error states from all worker threads into main statistics counters.
	 * @param pt Completed ProcessThread containing per-thread statistics and success status
	 */
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt){
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			alignmentCount+=pt.alignmentsT;
			seedHitCount+=pt.seedHitsT;

			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			errorState|=(!pt.success);
		}
	}

	/**
	 * Indicates whether all worker threads completed successfully without errors.
	 * Used by ThreadWaiter framework to determine overall job success status.
	 * @return True if no worker thread encountered errors during alignment processing
	 */
	@Override
	public final boolean success(){return !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Initializes concurrent input stream for reference sequences with automatic format detection.
	 * Creates FileFormat object, configures ConcurrentReadInputStream for multi-threaded access,
	 * and starts the input stream immediately. Stream will provide reference sequence batches
	 * to worker threads using producer-consumer pattern.
	 * @param fname File path to reference sequences (FASTA/FASTQ, optionally compressed)
	 * @return Started ConcurrentReadInputStream ready for batch consumption by worker threads
	 */
	private ConcurrentReadInputStream makeCris(String fname){
		FileFormat ff=FileFormat.testInput(fname, null, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
		cris.start(); // Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}

	/**
	* Loads query sequences from input files and converts them to Query objects.
	* Pre-computes k-mer indices and alignment parameters for each query.
	* @param ff1 Primary input file format
	* @param ff2 Secondary input file format (may be null)
	* @return List of Query objects ready for alignment
	*/
	public ArrayList<Query> fetchQueries(FileFormat ff1, FileFormat ff2){
		Timer t=new Timer(outstream, false);
		ArrayList<Read> reads=ConcurrentReadInputStream.getReads(maxReads, false, ff1, ff2, null, null);
		ArrayList<Query> queries=new ArrayList<Query>(reads.size());
		if(indexQueries){
			Query.mhc=new MinHitsCalculator(k, maxSubs, midMaskLen, minHitsProb, Query.maxClip); // Initialize hit calculator
		}
		for(Read r : reads){ //TODO: Could be multithreaded.
			readsProcessed+=r.pairCount();
			basesProcessed+=r.pairLength();
			Query query=new Query(r.id, queries.size(), r.bases, r.quality);
			queries.add(query);
			if(r.mate==null){continue;}
			r=r.mate; // Process mate if present
			query=new Query(r.id, queries.size(), r.bases, r.quality);
			queries.add(query);
		}
		t.stop("Loaded "+queries.size()+" queries in ");
		return queries;
	}

	/**
	* Performs sparse alignment using seed hits to guide alignment positions.
	* More efficient than brute force when seed hits are available and selective.
	* @param query Query sequence bases
	* @param ref Reference sequence bases  
	* @param maxSubs Maximum allowed substitutions
	* @param maxClips Maximum allowed clipped bases
	* @param seedHits List of potential alignment start positions from seed matching
	* @return List of alignment start positions with ≤maxSubs substitutions, or null if none found
	*/
	public static IntList alignSparse(byte[] query, byte[] ref, int maxSubs, int maxClips, IntList seedHits){
		if(seedHits==null || seedHits.isEmpty()){
			return null; // No seed guidance available
		}

		IntList results=null;

		for(int i=0; i<seedHits.size; i++){
			int rStart=seedHits.array[i];
			int subs;

			//Choose appropriate alignment method based on position
			if(rStart<0){
				subs=alignClipped(query, ref, maxSubs, maxClips, rStart); // Left overhang
			}else if(rStart>ref.length-query.length){
				subs=alignClipped(query, ref, maxSubs, maxClips, rStart); // Right overhang
			}else{
				subs=align(query, ref, maxSubs, rStart); // Perfect fit within reference
			}

			if(subs<=maxSubs){
				if(results==null){results=new IntList(4);}
				results.add(rStart);
			}
		}

		return results;
	}

	/**
	* Performs comprehensive alignment testing all possible positions.
	* Uses SIMD optimization when available and appropriate.
	* @param query Query sequence bases
	* @param ref Reference sequence bases
	* @param maxSubs Maximum allowed substitutions
	* @param maxClips Maximum allowed clipped bases
	* @return List of all alignment positions with ≤maxSubs substitutions, or null if none found
	*/
	public static IntList alignAllPositions(byte[] query, byte[] ref, int maxSubs, int maxClips){
		if(Shared.SIMD && (query.length<256 || maxSubs<256)){
//			return SIMDAlignByte.alignDiagonal(query, ref, maxSubs); // Use vectorized alignment
			return SIMDAlignByte.alignDiagonal(query, ref, maxSubs, maxClips); // Use vectorized alignment
			//TODO: Pad short contigs to avoid scalar mode
		}
		IntList list=null;
		int rStart=-maxSubs;
		
		//Left overhang region (negative start positions)
		for(; rStart<0; rStart++){
			int subs=alignClipped(query, ref, maxSubs, maxClips, rStart);
			if(subs<=maxSubs){
				if(list==null){list=new IntList(4);}
				list.add(rStart);
			}
		}
		
		//Perfect fit region (query completely within reference)
		for(final int limit=ref.length-query.length; rStart<=limit; rStart++){
			int subs=align(query, ref, maxSubs, rStart);
			if(subs<=maxSubs){
				if(list==null){list=new IntList(4);}
				list.add(rStart);
			}
		}
		
		//Right overhang region (query extends past reference end)
		for(final int limit=ref.length-query.length+maxSubs; rStart<=limit; rStart++){
			int subs=alignClipped(query, ref, maxSubs, maxClips, rStart);
			if(subs<=maxSubs){
				if(list==null){list=new IntList(4);}
				list.add(rStart);
			}
		}
		return list;
	}

	/**
	* Aligns query to reference starting at specified position with no clipping.
	* Optimized for cases where query fits completely within reference bounds.
	* @param query Query sequence bases
	* @param ref Reference sequence bases
	* @param maxSubs Maximum allowed substitutions (for early termination)
	* @param rStart Starting position in reference (0-based)
	* @return Number of substitutions found
	*/
	static int align(byte[] query, byte[] ref, final int maxSubs, final int rStart){
		int subs=0;
		for(int i=0, j=rStart; i<query.length && subs<=maxSubs; i++, j++){
			final byte q=query[i], r=ref[j];
			final int incr=(q!=r || AminoAcid.baseToNumber[q]<0 ? 1 : 0); // Count mismatches and N's
			subs+=incr;
		}
		return subs;
	}

	/**
	 * Performs alignment with soft clipping support for query overhangs beyond reference boundaries.
	 * Algorithm: 1) Calculate left and right clipping amounts based on position,
	 * 2) Count excess clipping beyond maxClips threshold as substitutions,
	 * 3) Align overlapping region counting mismatches and ambiguous bases,
	 * 4) Return total penalty including substitutions and excess clipping.
	 * Used for alignments where query extends past reference start (negative rStart)
	 * or end (rStart + query.length > ref.length).
	 * @param query Query sequence bases to align with clipping
	 * @param ref Reference sequence bases to align against
	 * @param maxSubs Maximum substitutions allowed in aligned region
	 * @param maxClips Maximum soft-clipped bases allowed without penalty
	 * @param rStart Starting position in reference (may be negative for left overhang)
	 * @return Total alignment penalty: substitutions + max(0, excess_clipping)
	 */
	static int alignClipped(byte[] query, byte[] ref, int maxSubs, final int maxClips, 
			final int rStart){
		final int rStop1=rStart+query.length; // Position after final base
		final int leftClip=Math.max(0, -rStart), rightClip=Math.max(0, rStop1-ref.length);
		int clips=leftClip+rightClip;
		if(clips>=query.length){return query.length;} // Entirely clipped
		int subs=Math.max(0, clips-maxClips); // Excess clipping counts as substitutions
		int i=leftClip, j=rStart+leftClip; // Skip clipped bases
		
		//Align overlapping region
		for(final int limit=Math.min(rStop1, ref.length); j<limit && subs<=maxSubs; i++, j++){
			final byte q=query[i], r=ref[j];
			final int incr=(q!=r || AminoAcid.baseToNumber[q]<0 ? 1 : 0);
			subs+=incr;
		}
		return subs;
	}

	/**
	* Builds a k-mer index for a reference sequence.
	* Maps masked k-mers to their positions for efficient seed finding.
	* @param ref Reference sequence bases
	* @return Hash map from masked k-mers to lists of positions, or null if indexing disabled
	*/
	IntListHashMap buildReferenceIndex(byte[] ref){
		if(!indexQueries || k<=0){return null;}

		final int defined=Math.max(k-maxSubs, 2);
		final int kSpace=(1<<(2*defined));
		final long maxKmers=Math.min(kSpace, (ref.length-k+1)*2L);
		final int initialSize=(int)Math.min(4000000, ((maxKmers*3)/2));
		final IntListHashMap index=new IntListHashMap(initialSize);

		final int shift=2*k, shift2=shift-2, mask=~((-1)<<shift); // Bit manipulation constants
		int kmer=0, rkmer=0, len=0; // Rolling k-mer state

		for(int i=0; i<ref.length; i++){
			final byte b=ref[i];
			final int x=AminoAcid.baseToNumber[b], x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask; // Roll forward k-mer
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask; // Roll reverse k-mer

			if(x<0){len=0; rkmer=0;}else{len++;} // Reset on ambiguous base
			if(len>=k){
				//Apply wildcard mask and store both orientations
				int maskedKmer=(kmer&Query.midMask);
				int maskedRkmer=(rkmer&Query.midMask);
				index.put(maskedKmer, i-k+1); // Store k-mer start position
				if(maskedKmer!=maskedRkmer){ // Avoid duplicate entries for palindromes
					index.put(maskedRkmer, i-k+1);
				}
			}
		}
		return index;
	}

	/**
	 * Converts alignment hits into properly formatted SAM output with full metadata.
	 * Algorithm: 1) For each hit position, generate match string using toMatch(),
	 * 2) Create SamLine object with proper coordinates, flags, and sequence orientation,
	 * 3) Generate CIGAR string from match string accounting for clipping and reference bounds,
	 * 4) Calculate mapping quality based on alignment length and substitution count,
	 * 5) Add NM optional tag with edit distance, 6) Batch output for efficiency.
	 * Sets primary flag only for first hit, uses appropriate strand orientation,
	 * and converts 0-based internal coordinates to 1-based SAM coordinates.
	 * @param q Query sequence containing bases, quality scores, and metadata
	 * @param ref Reference sequence containing bases and sequence identifier
	 * @param hits List of valid alignment start positions (0-based coordinates)
	 * @param reverseStrand True if alignments use reverse complement query sequence
	 * @param count Running total of alignments for this query (affects primary/secondary flags)
	 * @param bsw Thread-safe output writer for batched SAM format data
	 * @return Number of alignment records generated and queued for output
	 */
	static int processHits(Query q, Read ref, IntList hits, boolean reverseStrand, int count,
			ByteStreamWriter bsw){
		if(hits==null || hits.size()==0){return 0;}
		ByteBuilder bb=new ByteBuilder();
		ByteBuilder match=new ByteBuilder(q.bases.length);

		for(int i=0; i<hits.size(); i++){
			count++;
			int start=hits.get(i);

			//Use appropriate query sequence for match calculation
			byte[] querySeq=reverseStrand ? q.rbases : q.bases;
			toMatch(querySeq, ref.bases, start, match.clear());

			SamLine sl=new SamLine();
			sl.pos=Math.max(start+1, 1); // Convert to 1-based SAM coordinates
			sl.qname=q.name;
			sl.setRname(ref.id);
			sl.seq=q.bases; // Always report original query sequence
			sl.qual=q.quals;
			sl.setPrimary(count==1); // Primary if first hit overall
			sl.setMapped(true);
			if(reverseStrand){sl.setStrand(Shared.MINUS);}
			sl.tlen=q.bases.length;
			sl.cigar=SamLine.toCigar14(match.toBytes(), start, start+q.bases.length-1, ref.length(), q.bases);
			int subs=sl.countSubs();
			sl.addOptionalTag("NM:i:"+subs); // Add edit distance
			sl.mapq=Tools.mid(0, (int)(40*(sl.length()*0.5-subs)/(sl.length()*0.5)), 40); // Calculate mapping quality
			sl.toBytes(bb).nl();
			if(bb.length()>=16384){ // Batch output for efficiency
				if(bsw!=null){bsw.addJob(bb);}
				bb=new ByteBuilder();
			}
		}
		if(bsw!=null && !bb.isEmpty()){bsw.addJob(bb);}
		return hits.size;
	}

	/**
	 * Creates position-by-position match string for CIGAR string generation.
	 * Algorithm: For each query position, determine if reference position is in-bounds,
	 * compare bases using fully-defined nucleotide check, and assign match codes:
	 * 'm'=perfect match, 'S'=substitution within reference, 'C'=clipped outside reference.
	 * Uses sentinel character '$' for out-of-bounds reference positions.
	 * @param query Query sequence bases to align
	 * @param ref Reference sequence bases to compare against
	 * @param rStart Starting position in reference sequence (may be negative)
	 * @param match ByteBuilder to append match string characters
	 */
	static void toMatch(byte[] query, byte[] ref, int rStart, ByteBuilder match){
		for(int i=0, j=rStart; i<query.length; i++, j++){
			boolean inbounds=(j>=0 && j<ref.length);
			byte q=query[i];
			byte r=(inbounds ? ref[j] : (byte)'$'); // Use sentinel for out-of-bounds
			boolean good=(q==r && AminoAcid.isFullyDefined(q));
			match.append(good ? 'm' : inbounds ? 'S' : 'C'); // m=match, S=substitution, C=clip
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	* Worker thread for processing reference sequences against query set.
	* This class is static to prevent accidental writing to shared variables.
	* Each thread processes a portion of the reference stream independently.
	*/
	class ProcessThread extends Thread {

		/**
		 * Initializes worker thread with shared resources and processing parameters.
		 * Thread stores references to shared input/output streams but maintains private
		 * statistics counters and processing state to avoid synchronization overhead.
		 * @param cris_ Thread-safe input stream providing reference sequence batches
		 * @param bsw_ Thread-safe output writer for SAM format alignment results
		 * @param qList Shared read-only list of preprocessed Query objects
		 * @param maxSubs_ Maximum substitutions threshold for alignment acceptance
		 * @param tid_ Unique thread identifier for debugging and coordination
		 */
		ProcessThread(final ConcurrentReadInputStream cris_, final ByteStreamWriter bsw_, 
				ArrayList<Query> qList, final int maxSubs_, final int tid_){
			cris=cris_;
			bsw=bsw_;
			queries=qList;
			maxSubs=maxSubs_;
			tid=tid_;
		}

		/**
		 * Main thread execution method implementing standard BBTools thread pattern.
		 * Calls processInner() to perform actual work, then synchronizes to set
		 * success flag indicating completion without errors. ThreadWaiter framework
		 * monitors this flag to determine when to call accumulate() method.
		 */
		@Override
		/**
		* Main thread execution method called by start().
		* Processes reference sequences and performs alignments.
		*/
		public void run(){
			synchronized(this){
				processInner(); // Process all assigned reference sequences
				success=true; // Indicate successful completion
			}
		}

		/**
		 * Implements producer-consumer pattern consuming reference sequence batches until exhaustion.
		 * Algorithm: 1) Request next batch from shared input stream,
		 * 2) Process all references in batch against all queries,
		 * 3) Return batch to stream signaling completion,
		 * 4) Repeat until stream empty. Proper batch return ensures input stream
		 * can track completion and terminate cleanly when all batches processed.
		 */
		void processInner(){
			ListNum<Read> ln=cris.nextList(); // Grab the first batch of reference sequences

			while(ln!=null && ln.size()>0){
				processList(ln); // Process this batch
				
				cris.returnList(ln); // Notify input stream that batch was processed
				ln=cris.nextList(); // Fetch next batch
			}

			//Clean up final batch
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		/**
		 * Processes one batch of reference sequences performing validation and alignment.
		 * Algorithm: 1) For each reference in batch, validate sequence if needed,
		 * 2) Track statistics for read counts and base counts,
		 * 3) Align all queries against this reference using processRefSequence().
		 * Validation occurs in worker threads to avoid blocking main thread.
		 * @param ln ListNum wrapper containing reference sequence batch with completion tracking
		 */
		void processList(ListNum<Read> ln){
			final ArrayList<Read> refList=ln.list;

			//Process each reference sequence in the batch
			for(int idx=0; idx<refList.size(); idx++){
				final Read ref=refList.get(idx);

				if(!ref.validated()){ref.validate(true);} // Validate in worker threads for speed

				//Track statistics
				final int initialLength1=ref.length();
				final int initialLength2=ref.mateLength();
				readsProcessedT+=ref.pairCount();
				basesProcessedT+=initialLength1+initialLength2;

				//Align all queries against this reference sequence
				processRefSequence(ref);
			}
		}

		/**
		 * Routes seed hit finding to appropriate algorithm based on configuration.
		 * Delegates to either map-based counting approach (better for high hit density)
		 * or list-based approach (better for low hit density) based on useSeedMap flag.
		 * @param q Query sequence containing k-mer arrays for seed matching
		 * @param refIndex Reference k-mer index for lookup of matching positions
		 * @param reverseStrand True to use reverse complement k-mers from query
		 * @param hitCounts Reusable IntHashMap for position hit counting (map approach only)
		 * @param rname Reference sequence name for debugging output
		 * @return List of alignment start positions meeting minimum seed hit threshold
		 */
		private IntList getSeedHits(Query q, IntListHashMap refIndex, 
				boolean reverseStrand, IntHashMap hitCounts, String rname){
			if(useSeedMap){return getSeedHitsMap(q, refIndex, reverseStrand, hitCounts);}
			else{return getSeedHitsList(q, refIndex, reverseStrand);}
		}

		/**
		* Finds seed hits using list-based approach.  Potentially slow with short kmers.
		* @param q Query sequence to search for
		* @param refIndex K-mer index of reference
		* @param reverseStrand True to use reverse complement query
		* @return List of alignment positions meeting minimum hit threshold
		*/
		private IntList getSeedHitsList(Query q, IntListHashMap refIndex, boolean reverseStrand){
			int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			if(queryKmers==null){return null;}
			final int minHits=Math.max(minSeedHits, q.minHits);
			if(prescan){
				int valid=prescan(q, refIndex, reverseStrand, minHits);
				if(valid<minHits){return null;}
			}
			IntList seedHits=new IntList(); //TODO: Test speed with and without lazy init

			//Check each query k-mer for matches in reference
			for(int i=0; i<queryKmers.length; i+=qStep){
				if(queryKmers[i]==-1){continue;} // Skip invalid k-mers

				IntList positions=refIndex.get(queryKmers[i]);
				if(positions!=null){
					for(int j=0; j<positions.size; j++){
						int refPos=positions.array[j];
						int alignStart=refPos-i; // Adjust for k-mer position in query
						seedHits.add(alignStart);
					}
				}
			}
			if(seedHits==null){return null;}

			seedHitsT+=seedHits.size();
			if(seedHits.size<minHits){return null;}
			
			//Remove duplicates and filter by minimum occurrence
			if(seedHits.size>1 || minHits>1){
				seedHits.sort();
				seedHits.condenseMinCopies(minHits);
			}
			alignmentsT+=seedHits.size();
			return seedHits.isEmpty() ? null : seedHits;
		}

		/**
		* Finds seed hits using list-based approach.  Potentially slow with short kmers.
		* @param q Query sequence to search for
		* @param refIndex K-mer index of reference
		* @param reverseStrand True to use reverse complement query
		* @return List of alignment positions meeting minimum hit threshold
		*/
		private IntList getSeedHitsMap(Query q, IntListHashMap refIndex, 
				boolean reverseStrand, IntHashMap hitCounts){
			final int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			if(queryKmers==null){return null;}
			final int minHits=Math.max(minSeedHits, q.minHits);
			//TODO: Early exit if there are not enough ref kmers.
			if(prescan){
				int valid=prescan(q, refIndex, reverseStrand, minHits);
				if(valid<minHits){return null;}
			}

			IntList seedHits=null;

			//Process query k-mers at specified step interval
			for(int i=0; i<queryKmers.length; i+=qStep){
				if(queryKmers[i]==-1){continue;} // Skip invalid k-mers

				IntList positions=refIndex.get(queryKmers[i]);
				if(positions!=null){
					seedHitsT+=positions.size;
					if(seedHits==null){ // Lazy allocation
						seedHits=new IntList();
						if(hitCounts==null){hitCounts=new IntHashMap();}
						else{hitCounts.clear();}
					}
					for(int j=0; j<positions.size; j++){
						int alignStart=positions.array[j]-i; // Calculate alignment start

						//Increment hit count and add to results when threshold met
						int newCount=hitCounts.increment(alignStart);
						if(newCount==minHits){seedHits.add(alignStart);}
					}
				}
			}

			alignmentsT+=(seedHits==null ? 0 : seedHits.size());
			return seedHits==null || seedHits.isEmpty() ? null : seedHits;
		}

		/**
		* Finds seed hits using map-based approach for efficient hit counting.
		* More efficient when many hits are expected per position.
		* @param q Query sequence to analyze
		* @param refIndex Reference index for k-mer lookup
		* @param reverseStrand Flag indicating reverse complement strand
		* @return Number of query kmers shared with the ref
		*/
		private int prescan(Query q, IntListHashMap refIndex, boolean reverseStrand, final int minHits){
			final int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			final int maxMisses=queryKmers.length-minHits; //TODO: Does not account for qStep
			if(queryKmers==null || maxMisses<0){return 0;}

			int misses=0;
			//Process query k-mers at specified step interval
			for(int i=0; i<queryKmers.length && misses<=maxMisses; i+=qStep){
				if(queryKmers[i]==-1){continue;} // Skip invalid k-mers
				boolean hit=refIndex.containsKey(queryKmers[i]);
				misses+=(hit ? 0 : 1);
			}

			return queryKmers.length-misses;
		}

		/**
		 * Dispatches reference sequence processing to optimal alignment strategy.
		 * Uses indexed seed-and-extend approach when k-mer indexing enabled for efficiency
		 * with selective seed hits, otherwise uses brute force testing all positions.
		 * @param ref Reference sequence to align all queries against
		 * @return Total alignment count across all queries for this reference
		 */
		long processRefSequence(final Read ref){
			return indexQueries ? processRefSequenceIndexed(ref) : processRefSequenceBrute(ref);
		}

		/**
		* Processes reference sequence using indexed seed-and-extend approach.
		* More efficient for longer references with selective seed hits.
		* @param ref Reference sequence to process
		* @return Total number of alignments found across all queries
		*/
		long processRefSequenceIndexed(final Read ref){
//			Timer t=new Timer();
//			t.start("Indexing ref.");
			IntListHashMap refIndex=buildReferenceIndex(ref.bases); // Build k-mer index for this reference
			IntHashMap seedMap=(useSeedMap ? new IntHashMap() : null);
//			t.stopAndStart("Time:");
//			System.err.println("refIndex: "+refIndex.size());
//			System.err.println("seedMap: "+(seedMap==null ? 0 : seedMap.size()));
//			Shared.printMemory();

			long sum=0;
			for(Query q : queries){
				int count=0;

				// Forward strand alignment
				IntList seedHits=getSeedHits(q, refIndex, false, seedMap, ref.name());
				IntList hits=alignSparse(q.bases, ref.bases, maxSubs, q.maxClips, seedHits);
				count+=processHits(q, ref, hits, false, 0, bsw);
//				System.err.println("seedHits+:"+(seedHits==null ? 0 : seedHits.size()));

				// Reverse strand alignment
				seedHits=getSeedHits(q, refIndex, true, seedMap, ref.name());
				hits=alignSparse(q.rbases, ref.bases, maxSubs, q.maxClips, seedHits);
				count+=processHits(q, ref, hits, true, count, bsw);
//				System.err.println("seedHits-:"+(seedHits==null ? 0 : seedHits.size()));

				readsOutT+=count;
				basesOutT+=count*q.bases.length;
				sum+=count;
			}
//			t.stop("Time:");
//			Shared.printMemory();
			return sum;
		}

		/**
		* Processes reference sequence using indexed seed-and-extend approach.
		* More efficient for longer references with selective seed hits.
		* @param ref Reference sequence to process
		* @return Total number of alignments found across all queries
		*/
		long processRefSequenceBrute(final Read ref){
			long sum=0;
			for(Query q : queries){
				int count=0;

				// Forward strand alignment - test all positions
				IntList hits=alignAllPositions(q.bases, ref.bases, maxSubs, q.maxClips);
				count+=processHits(q, ref, hits, false, 0, bsw);

				// Reverse strand alignment using pre-computed reverse complement query
				hits=alignAllPositions(q.rbases, ref.bases, maxSubs, q.maxClips);
				count+=processHits(q, ref, hits, true, count, bsw);

				readsOutT+=count;
				basesOutT+=count*q.bases.length;
				sum+=count;
			}
			alignmentsT+=(queries.size()*(long)ref.length()); // All positions tested
			return sum;
		}

		/** Reference sequences processed by this worker thread */
		protected long readsProcessedT=0;
		/** Reference bases processed by this worker thread */
		protected long basesProcessedT=0;
		/** Alignment operations performed by this worker thread */
		protected long alignmentsT=0;
		/** K-mer seed hits found by this worker thread during indexed alignment */
		protected long seedHitsT=0;

		/** Valid alignments generated by this worker thread */
		protected long readsOutT=0;
		/** Total bases in valid alignments generated by this worker thread */
		protected long basesOutT=0;

		/** Success flag indicating thread completed without errors */
		boolean success=false;

		/** Thread-safe input stream providing reference sequence batches */
		private final ConcurrentReadInputStream cris;
		/** Thread-safe output writer for SAM format alignment results */
		private final ByteStreamWriter bsw;
		/** Shared list of preprocessed Query objects to align against references */
		private final ArrayList<Query> queries;
		/** Maximum substitutions threshold for alignment acceptance */
		final int maxSubs;
		/** Unique identifier for this worker thread */
		final int tid;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path for query sequences (FASTA/FASTQ, optionally compressed) */
	private String in1=null;
	/** Secondary input file path for paired-end query sequences (may be null for single-end) */
	private String in2=null;

	/** Output file path for SAM format alignment results */
	private String out1=null;

	/** Extension override for input file format detection (.fq, .fa, .sam, etc.) */
	private String extin=null;
	/** Extension override for output file format (.sam, .bam, etc.) */
	private String extout=null;

	/** File path to reference sequences for alignment (FASTA/FASTQ format) */
	String refFile=null;
	/** Maximum substitutions allowed per alignment (includes mismatches and N's) */
	int maxSubs=5;
	/** K-mer length for seed matching and indexing (1-15, affects sensitivity vs speed) */
	int k=13;
	/** Length of middle region to mask in k-mers for fuzzy matching tolerance */
	int midMaskLen=1;
	/** Enable k-mer indexing for seed-and-extend alignment strategy */
	boolean indexQueries=true;
	/** Perform fast pre-screening to count shared k-mers before full seed matching */
	boolean prescan=true;
	/** Sampling interval for query k-mers (1=every k-mer, 2=every other, reduces sensitivity) */
	int qStep=1;

	/** Minimum seed hits required for alignment consideration (higher = more selective) */
	int minSeedHits=1;
	/** Probability threshold for calculating minimum hits based on query length */
	private float minHitsProb=0.9999f;
	/** Use hash map for seed hit counting instead of list-based approach */
	boolean useSeedMap=false;

	/*--------------------------------------------------------------*/

	/** Total reference sequences processed across all worker threads */
	protected long readsProcessed=0;
	/** Total reference bases processed across all worker threads */
	protected long basesProcessed=0;

	/** Total alignment operations attempted across all worker threads */
	protected long alignmentCount=0;
	/** Total k-mer seed hits found during indexed alignment across all threads */
	protected long seedHitCount=0;

	/** Total valid alignments written to output files */
	protected long readsOut=0;
	/** Total bases in valid alignments written to output files */
	protected long basesOut=0;

	/** Maximum reference sequences to process (-1 for unlimited processing) */
	private long maxReads=-1;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** File format configuration for primary query input stream */
	private final FileFormat ffin1;
	/** File format configuration for secondary query input stream (paired-end) */
	private final FileFormat ffin2;

	/** File format configuration for SAM/BAM output stream */
	private final FileFormat ffout1;

	/**
	 * Provides read-write lock for thread-safe access to shared resources.
	 * Required by Accumulator interface for coordinating worker thread synchronization.
	 * @return ReentrantReadWriteLock instance for managing concurrent access
	 */
	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	/** Thread synchronization lock for coordinating access to shared statistics and resources */
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Output stream for status messages and progress reporting */
	private PrintStream outstream=System.err;
	/** Enable verbose output including timing and debugging information */
	public static boolean verbose=false;
	/** Global error flag set when any worker thread encounters an error */
	public boolean errorState=false;
	/** Allow overwriting of existing output files */
	private boolean overwrite=true;
	/** Append results to existing output files instead of overwriting */
	private boolean append=false;

}
