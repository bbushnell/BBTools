package ifa;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMask;
import map.IntHashMap2;
import map.IntHashMap3;
import map.IntListHashMap3;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.SIMDAlignByte;
import simd.Vector;
import stream.Streamer;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamHeader;
import stream.SamHeaderWriter;
import stream.SamLine;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import structures.StringNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;
import tracker.ReadStats;

/**
 * Version 2 indel-free aligner using IntListHashMap3 reference indexes.
 * Materializes query buckets, then aligns both strands against reference batches
 * consumed by workers from one shared Streamer. Indexed mode uses masked two-bit
 * k-mer keys, a presence prescan and fresh seed lists; brute mode enumerates offsets.
 * Helpers conditionally delegate to SIMD. Optional fusion copies roots with N
 * padding and maps candidates back to original references using the query center.
 * <p>
 * Workers own indexes, seed scratch and counters. Query arrays are shared for
 * reading; an atomic query counter chooses the first counted alignment as primary.
 * Alignment output is headerless and unordered; an optional separate header writer
 * collects reference metadata. Query-root length filtering also gates mate admission.
 * <p>
 * Parsing/query setup changes shared I/O, Query, calculator and scoring settings.
 * Instances are not independent concurrent configurations. The normal process path
 * joins workers/writers and restores selected buffer/validation settings; earlier
 * exceptions can bypass restoration because it is not protected by finally.
 * @author Brian Bushnell
 * @contributor Isla, Amber
 * @date June 2, 2025
 */
public class IndelFreeAligner2 implements Accumulator<IndelFreeAligner2.ProcessThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs one CLI invocation and closes its diagnostic stream on normal completion. */
	public static void main(String[] args){
		Timer t=new Timer();
		IndelFreeAligner2 x=new IndelFreeAligner2(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/** Parses options/shared settings and prepares formats; brute mode uses one k=0 bucket. */
	public IndelFreeAligner2(String[] args){

		{ //Preparse block
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

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

		if(kArray==null || kArray.length<1 || (kArray.length==1 && kArray[0]==0)) {indexQueries=false;}
		if(indexQueries) {
			Arrays.sort(kArray);
			Tools.reverseInPlace(kArray);
		}else {
			kArray=new int[] {0};
		}
		kStep=Math.max(qStep, rStep);
		assert(qStep==1 || rStep==1) : "Don't use both qStep and rStep at once.";
		assert(kStep>=1) : "qStep and rStep must be at least 1: "+qStep+", "+rStep;
		assert(Integer.bitCount(rStep)==1) : "rStep must be a power of 2: "+rStep;

		Shared.BBMAP_CLASS=" "+this.getClass().getName();
		SamHeader.PN="IndelFreeAligner";
		validateParams();
		doPoundReplacement();
		fixExtensions();
		checkFileExistence();
		checkStatics();

		ffout1=FileFormat.testOutput(out1, FileFormat.SAM, extout, true, overwrite, append, false);
		ffheader=FileFormat.testOutput(headerOut, FileFormat.SAM, extout, true, overwrite, false, true);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses tool/common options and shared settings; k accepts scalars, lists and ranges. */
	private Parser parse(String[] args){
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			//FIXED IFA-014: preserve the full value after the first '=', including further separators.
			//Keep the existing null value for an absent or empty suffix (for example, ref=).
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=(split.length>1 && !split[1].isEmpty()) ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ref")){
				refFile=b;
			}else if(a.equals("subs") || a.equals("maxsubs") || a.equals("s")){
				maxSubs=Integer.parseInt(b);
			}else if(a.equals("ani") || a.equals("minani") || a.equals("identity") || a.equals("id") || a.equals("minid")){
				minid=Float.parseFloat(b);
				if(minid>1) {minid/=100f;}
			}else if(a.equals("hits") || a.equals("minhits") || a.equals("seedhits")){
				minSeedHits=Math.max(1, Integer.parseInt(b));
			}else if(a.equals("minprob") || a.equals("minhitsprob")){
				minHitsProb=Float.parseFloat(b);
			}else if(a.equals("iterations")){
				MinHitsCalculator3.iterations=MinHitsCalculator2.iterations=MinHitsCalculator.iterations=Parse.parseIntKMG(b);
			}else if(a.equals("maxclip") || a.equals("clip")){
				Query.maxClip=Tools.max(0, Float.parseFloat(b));
			}else if(a.equals("index")){
				indexQueries=Parse.parseBoolean(b);
			}else if(a.equals("brute") || a.equals("bruteforce")){
				indexQueries=!Parse.parseBoolean(b);
			}else if(a.equals("prescan")){
				prescan=Parse.parseBoolean(b);
			}else if(a.equals("seedmap") || a.equals("map")){
				useSeedMap=Parse.parseBoolean(b);
			}else if(a.equals("seedlist") || a.equals("list")){
				useSeedMap=!Parse.parseBoolean(b);
			}else if(a.equals("header") || a.equals("headerout") || a.equals("outheader") || a.equals("outh")){
				headerOut=b;
			}else if(a.equals("k")){
				if(b==null || b.equals("0") || b.equals("-1")){
					kArray=null;
				}else if(b.indexOf(',')>-1){
					int[] temp=Parse.parseIntArray(b, ",");
					kArray=temp;
				}else if(b.indexOf('-')>-1){
					int[] temp=Parse.parseIntArray(b, "-");
					Arrays.sort(temp);
					int range=temp[1]-temp[0]+1;
					kArray=new int[range];
					for(int x=0; x<range; x++) {kArray[x]=x+temp[0];}
				}else{
					kArray=new int[] {Integer.parseInt(b)};
				}
				indexQueries=(kArray!=null && kArray.length>0 && kArray[0]>0);
			}else if(a.equals("qstep") || a.equals("step") || a.equals("qskip")){
				qStep=Integer.parseInt(b);
			}else if(a.equals("rstep") || a.equals("rskip")){
				rStep=Integer.parseInt(b);
			}else if(a.equals("qlen") || a.equals("minqlen")){
				minQLen=Integer.parseInt(b);
			}else if(a.equals("rlen") || a.equals("minrlen")){
				minRLen=Integer.parseInt(b);
			}else if(a.equals("minlen")){
				minQLen=minRLen=Integer.parseInt(b);
			}else if(a.equals("mm")){
				midMaskLen=(Tools.isNumeric(b) ? Integer.parseInt(b) : Parse.parseBoolean(b) ? 1 : 0);
			}else if(a.equals("blacklist") || a.equals("banhomopolymers")){
				Query.blacklistRepeatLength=(Tools.isNumeric(b) ? Integer.parseInt(b) : Parse.parseBoolean(b) ? 1 : 0);
			}else if(a.equals("fuse")){
				fuse=Parse.parseBoolean(b);
			}else if(a.equals("padding")){
				padding=Integer.parseInt(b);
			}else if(a.equals("chunk") || a.equals("chunksize")){
				targetChunkSize=Parse.parseIntKMG(b);
			}else if(a.equals("entropymask") || a.equals("emask") || a.equals("mask")){
				entropyMask=Parse.parseBoolean(b);
			}else if(a.equals("entropywindow") || a.equals("ewindow") || a.equals("ew") || a.equals("window")){
				entropyWindow=Integer.parseInt(b);
			}else if(a.equals("entropycutoff") || a.equals("ecutoff") || a.equals("minentropy")){
				entropyCutoff=Float.parseFloat(b);
			}else if(a.equals("entropyk") || a.equals("ek") || a.equals("ke")){
				entropyK=Integer.parseInt(b);	
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(parser.out1==null && b==null && FileFormat.isSamOrBamFile(arg)){
				parser.out1=arg;
			}else if(parser.in1==null && b==null && FileFormat.isFastqFile(arg) && new File(arg).isFile()){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		return parser;
	}

	/** Expands a nonexistent query path containing # into mate paths; requires query input. */
	private void doPoundReplacement(){
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}

	/** Resolves alternate compressed extensions for query inputs. */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
	}

	/** Checks query inputs, outputs and duplicate paths; reference opening occurs later. */
	private void checkFileExistence(){
		if(!Tools.testOutputFiles(overwrite, append, false, out1, headerOut)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, headerOut)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/** Selects a shared byte-reader mode if unset and checks shared FASTA settings. */
	private static void checkStatics(){
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		assert(FastaReadInputStream.settingsOK());
	}

	/** Checks existing assertion-based constraints; not an exhaustive input validator. */
	private boolean validateParams(){
		for(int k : kArray){
			assert((k>=1 && k<=15) || !indexQueries);
			assert(midMaskLen<k-1 || !indexQueries);
		}
		assert(minHitsProb<=1);
		assert(maxSubs>=0);
		assert(qStep>=1) : "qStep must be at least 1: "+qStep;//qStep is an unconditional loop stride (i+=qStep); 0 would hang. rStep>=1 already ensured by bitCount(rStep)==1.
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads queries, starts reference/output services and joins workers.
	 * Then closes input, waits for alignment/header output and reports errors.
	 * Selected buffer/Read-validation settings are restored along this normal path
	 * before the final error-state exception; earlier failures can bypass restoration.
	 * @param t Caller-started timer, stopped after stream finalization
	 */
	void process(Timer t){
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		SamLine.RNAME_AS_BYTES=false;

		final ArrayList<ArrayList<Query>> queryBuckets=fetchQueries(ffin1, ffin2);

		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;

		int oldBD=Shared.bufferData();
		int oldBL=Shared.bufferLen();
		// If fusing, allow large chunks to maximize throughput. If not, keep small for parallelism.
		Shared.setBufferData(Math.min(targetChunkSize, (fuse ? targetChunkSize : 400000)));
		if(fuse) {Shared.setBufferLen(Math.max(oldBL, 2000));}

		final Streamer cris=makeCris(refFile);
		final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout1);
		final SamHeaderWriter shw=(ffheader==null ? null : new SamHeaderWriter(ffheader));

		spawnThreads(cris, bsw, shw, queryBuckets);

		if(verbose){outstream.println("Finished; closing streams.");}

		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(cris);

		//FIXED IFA-013: retain the output writer's flush/close/subprocess error flag for the final failure check.
		//Immediate write failures can abort in ByteStreamWriter; finalization failures may only set this flag.
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		if(shw!=null) {shw.poisonAndWait();}

		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		Shared.setBufferData(oldBD);
		Shared.setBufferLen(oldBL);

		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		outstream.println(Tools.things("Alignments", alignmentCount, 8));
		outstream.println(Tools.things("Seed Hits", seedHitCount, 8));

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts workers sharing services/buckets, waits for them and accumulates success. */
	private void spawnThreads(final Streamer cris, 
		final ByteStreamWriter bsw, final SamHeaderWriter shw, final ArrayList<ArrayList<Query>> queryBuckets){

		final int threads=Shared.threads();
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, bsw, shw, queryBuckets, maxSubs, minid, i));
		}
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState|=!success;
	}

	/** Merges a worker's counters and success while holding its monitor. */
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

	/** Reports errors accumulated so far; this is not a completion test. */
	@Override
	public final boolean success(){return !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates and starts the shared reference reader with maxReads and SAM-header retention. */
	private Streamer makeCris(String fname){
		FileFormat ff=FileFormat.testInput(fname, null, true);
		Streamer cris=StreamerFactory.getReadInputStream(maxReads, true, ff, null, -1);
		cris.start();
		if(verbose){outstream.println("Started cris");}
		return cris;
	}

	/** Unused scaffolding: builds and discards a root-length histogram; does not select k values. */
	private void analyzeQueries(ArrayList<Read> reads) {
		//TODO: Find qlen range and counts, optimal kmer length for each,
		//optionally autoselect k for ranges with SUFFICIENT members (or length) to make it useful.
		final int total=reads.size();
		if(total<1) {return;}
		IntHashMap3 lengthMap=new IntHashMap3();
		IntList lengthList=new IntList();
		for(Read r : reads) {
			final int length=r.length();
			int ret=lengthMap.increment(length);
			if(ret==1) {lengthList.add(length);}
		}
		lengthList.sort();
		//Now, decide which lengths to down-merge.
	}
	
	/**
	 * Materializes queries, installs shared calculators and prepares buckets.
	 * Counts every loaded root and mate before filtering. A root meeting minQLen
	 * admits itself and its mate; the mate has no separate length check here.
	 * @param ff1 Primary query input
	 * @param ff2 Optional separate mate input
	 * @return Owned Query buckets completed before alignment workers start
	 */
	public ArrayList<ArrayList<Query>> fetchQueries(FileFormat ff1, FileFormat ff2){
		Timer t=new Timer(outstream, false);
		ArrayList<Read> reads=StreamerFactory.getReads(maxReads, false, ff1, ff2, null, null);
		
		ArrayList<ArrayList<Query>> buckets=new ArrayList<>(kArray.length);
		for(int i=0; i<kArray.length; i++){
			buckets.add(new ArrayList<Query>(reads.size()/kArray.length));
		}

		MinHitsCalculator2[] calculators=(indexQueries ? new MinHitsCalculator2[kArray.length] : null);
		if(indexQueries){
			for(int i=0; i<kArray.length; i++){
				calculators[i]=new MinHitsCalculator2(kArray[i], maxSubs, minid, midMaskLen, 
					minHitsProb, Query.maxClip, Math.max(qStep, rStep));
			}
		}
		Query.setCalculators(calculators, minSeedHits);

		for(Read r : reads){ //TODO: Could be multithreaded.
			readsProcessed+=r.pairCount();
			basesProcessed+=r.pairLength();

			if(r.length()>=minQLen) {
				addReadToBucket(r, buckets);
				if(r.mate!=null){
					addReadToBucket(r.mate, buckets);
				}
			}
		}

		long totalLoaded=0;
		for(int i=0; i<buckets.size(); i++){
			totalLoaded+=buckets.get(i).size();
			if(verbose){outstream.println("K="+kArray[i]+": "+buckets.get(i).size()+" queries.");}
		}

		t.stop("Loaded "+totalLoaded+" queries in ");
		return buckets;
	}

	/** Wraps borrowed read bases/qualities in Query and selects its calculator bucket. */
	private void addReadToBucket(Read r, ArrayList<ArrayList<Query>> buckets){
		Query q=new Query(r.id, 0, r.bases, r.quality);
		if(q.calculatorIndex>=0){
			buckets.get(q.calculatorIndex).add(q);
		}else{
			// Handle queries that failed to match any K (too short?) 
			// or brute force mode. 
			// If brute force (calculatorIndex -1), put in bucket 0 so they get processed flatly
			buckets.get(0).add(q); 
		}
	}

	/**
	 * Scores borrowed zero-based offsets without modifying them. Overhangs use
	 * alignClipped; interior placements delegate to Vector.align. maxClips is the
	 * free overhang allowance; excess clips consume maxSubs with substitutions.
	 * @return New accepted-offset list in candidate order, or null if none pass
	 */
	public static IntList alignSparse(byte[] query, byte[] ref, int maxSubs, int maxClips, IntList seedHits){
		if(seedHits==null || seedHits.isEmpty()){return null;}
		IntList results=null;
		for(int i=0; i<seedHits.size; i++){
			int rStart=seedHits.array[i];
			int subs;
			if(rStart<0){
				subs=alignClipped(query, ref, maxSubs, maxClips, rStart);
			}else if(rStart>ref.length-query.length){
				subs=alignClipped(query, ref, maxSubs, maxClips, rStart);
			}else{
				subs=Vector.align(query, ref, maxSubs, rStart);
			}
			if(subs<=maxSubs){
				if(results==null){results=new IntList(4);}
				results.add(rStart);
			}
		}
//		System.err.println("as: "+results);
		return results;
	}
	
	/**
	 * Scores borrowed offsets in a fused reference. Interior placements delegate to
	 * Vector.alignFused; outer overhangs use alignClipped. Original-reference selection
	 * and final substitution/paid-clipping checks occur later in processHits.
	 * @return New accepted-offset list in candidate order, or null if none pass
	 */
	public static IntList alignSparseFused(byte[] query, byte[] ref, int maxSubs, int maxClips, IntList seedHits){
		if(seedHits==null || seedHits.isEmpty()){return null;}
		IntList results=null;
		final int rLen=ref.length;
		final int qLen=query.length;

		for(int i=0; i<seedHits.size; i++){
			int rStart=seedHits.array[i];
			int subs;

			// If we hit the absolute edge of the chunk buffer, use clipped alignment 
			// to avoid ArrayOutOfBounds.
			if(rStart<0 || rStart>rLen-qLen){
				subs=alignClipped(query, ref, maxSubs, maxClips, rStart);
			}else{
				// Internal: Safe to use SIMD, but might hit internal Padding (Ns).
				// Pass maxClips as the N-threshold.
				subs=Vector.alignFused(query, ref, maxSubs, maxClips, rStart);
			}

			if(subs<=maxSubs){
				if(results==null){results=new IntList(4);}
				results.add(rStart);
			}
		}
//		System.err.println("asf: "+results);
		return results;
	}

	/**
	 * Enumerates offsets for non-null arrays and nonnegative budgets. SIMD dispatch
	 * precedes the scalar path and has its own candidate enumeration. Scalar offsets
	 * require nonempty overlap; clips beyond maxClips consume maxSubs. Scalar results
	 * are ascending, or null if none pass. No scalar/SIMD equivalence is implied here.
	 */
	public static IntList alignAllPositions(byte[] query, byte[] ref, int maxSubs, int maxClips){
		if(Shared.SIMD && (query.length<256 || maxSubs<256)){
			return SIMDAlignByte.alignDiagonal(query, ref, maxSubs, maxClips);
		}
		if(query.length==0 || ref.length==0){return null;}
		IntList list=null;
		//FIXED IFA-007: include free and substitution-paid clipping, bounded by nonempty overlap.
		final long budget=(long)maxSubs+maxClips;
		final int last=(int)Math.min((long)ref.length-1, (long)ref.length-query.length+budget);
		int rStart=(int)Math.max(1L-query.length, -budget);
		for(; rStart<0 && rStart<=last; rStart++){
			int subs=alignClipped(query, ref, maxSubs, maxClips, rStart);
			if(subs<=maxSubs){
				if(list==null){list=new IntList(4);}
				list.add(rStart);
			}
		}
		for(final int limit=ref.length-query.length; rStart<=limit; rStart++){
			int subs=align(query, ref, maxSubs, rStart);
			if(subs<=maxSubs){
				if(list==null){list=new IntList(4);}
				list.add(rStart);
			}
		}
		for(; rStart<=last; rStart++){
			int subs=alignClipped(query, ref, maxSubs, maxClips, rStart);
			if(subs<=maxSubs){
				if(list==null){list=new IntList(4);}
				list.add(rStart);
			}
		}
		return list;
	}

	/**
	 * Scores a fully in-bounds placement, counting unequal or ambiguous query bases.
	 * Stops after exceeding maxSubs; a rejected score may be only a partial count.
	 */
	static int align(byte[] query, byte[] ref, final int maxSubs, final int rStart){
		int subs=0;
		for(int i=0, j=rStart; i<query.length && subs<=maxSubs; i++, j++){
			final byte q=query[i], r=ref[j];
			final int incr=(q!=r || AminoAcid.baseToNumber[q]<0 ? 1 : 0);
			subs+=incr;
		}
		return subs;
	}

	/**
	 * Scores overhangs, charging clips beyond maxClips. No overlap returns query.length;
	 * otherwise scoring stops after exceeding maxSubs. This pointwise helper does
	 * not enforce the scalar candidate-enumeration domain.
	 */
	static int alignClipped(byte[] query, byte[] ref, int maxSubs, final int maxClips, 
		final int rStart){
		final int rStop1=rStart+query.length;
		final int leftClip=Math.max(0, -rStart), rightClip=Math.max(0, rStop1-ref.length);
		int clips=leftClip+rightClip;
		if(clips>=query.length){return query.length;}
		int subs=Math.max(0, clips-maxClips);
		int i=leftClip, j=rStart+leftClip;
		for(final int limit=Math.min(rStop1, ref.length); j<limit && subs<=maxSubs; i++, j++){
			final byte q=query[i], r=ref[j];
			final int incr=(q!=r || AminoAcid.baseToNumber[q]<0 ? 1 : 0);
			subs+=incr;
		}
		return subs;
	}

	/**
	 * Builds owned masked forward-key lists of zero-based reference starts.
	 * Ambiguous bases invalidate windows even at masked positions; rStep samples
	 * starts at multiples of its power-of-two stride. Reads but does not retain or
	 * mutate ref. Short references produce an empty map; disabled/nonpositive k
	 * returns null. Indexed k is expected to be 1 through 15.
	 */
	IntListHashMap3 buildReferenceIndex(byte[] ref, int k){
		if(!indexQueries || k<=0){return null;}
		final int defined=Math.max(k-midMaskLen, 2);
		final int kSpace=(1<<(2*defined));
		final long maxKmers=Math.min(kSpace, (ref.length-k+1)*2L);
		//FIXED IFA-006: short references still need a valid empty map for prescan and seed consumers.
		final int initialSize=(int)Math.max(1, Math.min(4000000, ((maxKmers*3)/2)));
		final IntListHashMap3 index=new IntListHashMap3(initialSize, 0.7);

		final int shift=2*k, mask=~((-1)<<shift);
		final int stepMask=(rStep-1);
		final int stepTarget=((k-1)&stepMask);
		int kmer=0, len=0;

		for(int i=0; i<ref.length; i++){
			final byte b=ref[i];
			final int x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){len=0; kmer=0;}else{len++;}
			if(len>=k && ((i&stepMask)==stepTarget)){
				// Use static helper in Query for the mask
				int maskedKmer=(kmer&Query.makeMidMask(k, midMaskLen));
				index.put(maskedKmer, i-k+1);
			}
		}
		return index;
	}

//	/**
//	 * Translates fused coordinate to original contig and validates position.
//	 * @return SamLine if valid, null otherwise
//	 */
//	@Deprecated
//	static SamLine makeSamLine(Query q, Read ref, int start, 
//		ArrayList<Read> originalRefs, IntList starts){
//
//		int rStart=start;
//		Read realRef=ref;
//
//		if(originalRefs!=null){
//			// Translation Mode
//			int idx=Arrays.binarySearch(starts.array, 0, starts.size, start);
//			if(idx<0){idx=-idx-2;} // Calculate insertion point
//
//			if(idx<0 || idx>=originalRefs.size()){return null;} // Should be impossible
//
//			int refOffset=starts.get(idx);
//			realRef=originalRefs.get(idx);
//			rStart=start-refOffset;
//
//			// Validate bounds: check if alignment extends past end of real contig
//			// (Assuming N padding prevents cross-contig alignment, but safety first)
//			if(rStart<0 || rStart+q.length()>realRef.length()){return null;}
//		}
//
//		SamLine sl=new SamLine();
//		sl.pos=Math.max(rStart+1, 1);
//		sl.qname=q.name;
//		sl.setRname(realRef.id); // Use the original ID
//		sl.seq=q.bases;
//		sl.qual=q.quals;
//		sl.setPrimary(q.alignments.incrementAndGet()==1);
//		sl.setMapped(true);
//		sl.tlen=q.bases.length;
//		return sl;
//	}
	
	/*
	 * Historical fused-hit fix: use center-based interleaved range lookup
	 * to recover the original contig before checking local clipping and substitutions.
	 */
	/**
	 * Rechecks candidates against an original reference and builds SAM alignments.
	 * Fused candidates select a reference by query center in paired [start,stop)
	 * ranges, then convert to local offsets. originalRefs/ranges are unused in
	 * standard mode. A preliminary check rejects excess clipping above maxSubsQ;
	 * the combined substitution and excess-clipping cost must also fit maxSubsQ.
	 * NM counts substitutions only. The atomic query counter assigns
	 * primary status to the first counted alignment. Null bsw still builds/counts
	 * alignments and advances that counter. Borrowed hits are not modified.
	 * @return Accepted alignment count, not distinct query count
	 */
	static int processHits(Query q, Read ref, IntList hits, boolean reverseStrand,
		ByteStreamWriter bsw, ArrayList<Read> originalRefs, IntList ranges, int maxSubsQ){

		if(hits==null || hits.size()==0){return 0;}
		ByteBuilder bb=new ByteBuilder();
		ByteBuilder match=new ByteBuilder(q.bases.length);

		int added=0; 
		final int qLen=q.length();
		final int qHalf=qLen/2; // Center offset

		for(int i=0; i<hits.size(); i++){
			int start=hits.get(i);
			
			// 1. Recover the Real Reference
			Read realRef=ref;
			int rStart=start;

			if(originalRefs!=null){
				// SEARCH: Find which chunk contains the CENTER of the read.
				// ranges is [Start, Stop, Start, Stop...]
				int center = start + qHalf;
				int idx=Arrays.binarySearch(ranges.array, 0, ranges.size, center);
				if(idx<0){idx=-idx-2;}
				
				// VALIDATION:
				// idx must be nonnegative.
				// An EVEN index selects [Start, Stop).
				// An ODD index selects padding [Stop, NextStart) -> Invalid.
				if(idx < 0 || (idx&1)==1){ continue; }
				
				int refIdx = idx/2;
				if(refIdx >= originalRefs.size()){ continue; }

				realRef=originalRefs.get(refIdx);
				
				// Calculate start relative to the Chunk Start (ranges[idx])
				rStart=start-ranges.get(idx);
			}

			// 2. Calculate Clipping relative to Real Reference
			final int rLen=realRef.length();
			final int rStop=rStart+qLen;

			final int leftClip=Math.max(0, -rStart);
			final int rightClip=Math.max(0, rStop-rLen);
			final int totalClip=leftClip+rightClip;
			
			// 3. Calculate Clip Penalty
			int clipPenalty = Math.max(0, totalClip - q.maxClips);
			if(clipPenalty > maxSubsQ){continue;}

			// 4. Calculate Substitutions
			byte[] querySeq=reverseStrand ? q.rbases : q.bases;
			toMatch(querySeq, realRef.bases, rStart, match.clear());

			int subs=0;
			for(int m=0; m<match.length(); m++){
				final byte b=match.get(m);
				if(b=='S' || b=='N'){subs++;}
			}

			// 5. Final Filter
			if((subs + clipPenalty) > maxSubsQ){continue;}

			// 6. Build Valid SAM Line
			added++; 

			SamLine sl=new SamLine();
			sl.pos=Math.max(rStart+1, 1);
			sl.qname=q.name;
			sl.setRnameS(realRef.id);
			sl.setSeq(q.bases);
			sl.setQual(q.quals);
			sl.setPrimary(q.alignments.incrementAndGet()==1);
			sl.setMapped(true);
			if(reverseStrand){sl.setStrand(Shared.MINUS);}
			sl.tlen=qLen;
			sl.setCigar(SamLine.toCigar14(match.toBytes(), rStart, rStop-1, rLen, q.bases));
			sl.addOptionalTag("NM:i:"+subs);
			sl.mapq=Tools.mid(0, (int)(40*(sl.length()*0.5-subs)/(sl.length()*0.5)), 40);
			sl.toBytes(bb).nl();

			if(bb.length()>=16384){
				if(bsw!=null){bsw.addJob(bb);}
				bb=new ByteBuilder();
			}
		}
		if(bsw!=null && !bb.isEmpty()){bsw.addJob(bb);}
		return added;
	}

	/** Appends m for defined equal bases, S for in-bounds mismatches and C for overhangs. */
	static void toMatch(byte[] query, byte[] ref, int rStart, ByteBuilder match){
		for(int i=0, j=rStart; i<query.length; i++, j++){
			boolean inbounds=(j>=0 && j<ref.length);
			byte q=query[i];
			byte r=(inbounds ? ref[j] : (byte)'$');
			boolean good=(q==r && AminoAcid.isFullyDefined(q));
			match.append(good ? 'm' : inbounds ? 'S' : 'C');
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Entropy           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Masks bases in place and returns BBMask's count; an absent tracker is local to this call. */
	private int entropyMask(byte[] bases, EntropyTracker et) {
		if(bases==null || bases.length==0){return 0;}

		// 1. Create the shared mask (The "BitSet" point I missed!)
		BitSet bs=new BitSet(bases.length);

		// 2. Mark Low Entropy
		// We allocate the tracker once.
		if(et==null) {et=new EntropyTracker(entropyK, entropyWindow, false, entropyCutoff, true);}
		BBMask.maskLowEntropy(bases, bs, et);

//		// 3. Mark Repeats
//		// This sees the ORIGINAL bases, not Ns, so it finds repeats correctly.
//		BBMask.maskRepeats(bases, bs, 5, 20);

		// 4. Apply the mask in place; this is not a whole-array atomic update.
		// Uses the existing public method in BBMask
		return BBMask.maskBases(bases, bs, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Shared-reader consumer owning reference indexes/counters; borrows queries and writers. */
	class ProcessThread extends Thread {

		/** Borrows services/buckets and copies effective scoring settings for this worker. */
		ProcessThread(final Streamer cris_, final ByteStreamWriter bsw_, final SamHeaderWriter shw_, 
			ArrayList<ArrayList<Query>> qBuckets_, final int maxSubs_, final float minid_, final int tid_){
			cris=cris_;
			bsw=bsw_;
			shw=shw_;
			queryBuckets=qBuckets_;
			maxSubs=maxSubs_;
			minid=minid_;
			tid=tid_;
		}

		/** Initializes worker scratch and marks success only after normal processing completion. */
		@Override
		public void run(){
			synchronized(this){
				if(entropyMask) {et=new EntropyTracker(entropyK, entropyWindow, false, entropyCutoff, true);}
				processInner();
				success=true;
			}
		}

		/** Consumes shared-Streamer batches until its null terminal. */
		void processInner(){
			for(ListNum<Read> ln=cris.nextList(); ln!=null; ln=cris.nextList()){processList(ln);}
		}

		/**
		 * Removes short roots in place, validates survivors and counts roots plus mates.
		 * Queues optional header entries under the original batch ID, including empty
		 * header batches, before dispatching retained roots to standard/fused alignment.
		 */
		void processList(ListNum<Read> ln){
			final ArrayList<Read> refList=ln.list;
			final ArrayList<StringNum> alsn=(shw==null ? null : new ArrayList<StringNum>(refList.size()));
			int removed=0;
			// 1. Validation and Header Registration (Keep this!)
			for(int idx=0; idx<refList.size(); idx++){
				final Read ref=refList.get(idx);
				//FIXED IFA-010: condenseStrict removes nulls before either standard or fused alignment.
				if(ref.length()<minRLen){refList.set(idx, null); removed++; continue;}
				if(!ref.validated()){ref.validate(true);}
				if(alsn!=null) {alsn.add(new StringNum(ref.name(), ref.length()));}

				// Move stats counting here to ensure we count original inputs
				final int initialLength1=ref.length();
				final int initialLength2=ref.mateLength();
				readsProcessedT+=ref.pairCount();
				basesProcessedT+=initialLength1+initialLength2;
			}
			if(shw!=null) {shw.add(new ListNum<StringNum>(alsn, ln.id));}
			if(removed>0) {Tools.condenseStrict(refList);}

			// 2. Fork Logic: Fuse vs Standard
			if(fuse){processListFused(refList);}
			else{processListStandard(refList);}
		}

		/** Aligns reference roots independently; entropy masking may mutate their bases. */
		void processListStandard(ArrayList<Read> refList){
			for(int idx=0; idx<refList.size(); idx++){
				final Read ref=refList.get(idx);
				processRefSequence(ref, null, null);
			}
		}

		/**
		 * Copies roots with trailing N padding and records paired [start,stop) ranges.
		 * Falls back to standard processing when the computed total exceeds int.
		 */
		void processListFused(ArrayList<Read> refList){
			if(refList.isEmpty()){return;}

			// Calculate exact size needed
			long totalLen=0;
			for(Read r : refList){
				totalLen+=r.length()+padding;
			}
			if(totalLen>Integer.MAX_VALUE){
				// Fall back when the computed concatenation length exceeds the int-sized buffer.
				processListStandard(refList);
				return;
			}

			ByteBuilder bb=new ByteBuilder((int)totalLen);
			// Interleaved Ranges: [Start, Stop, Start, Stop...]
			IntList ranges=new IntList(refList.size()*2);

			for(Read r : refList){
				ranges.add(bb.length()); // Start (Inclusive)
				bb.append(r.bases);
				ranges.add(bb.length()); // Stop (Exclusive)
				for(int i=0; i<padding; i++){bb.append('N');}
			}

			Read fusedRef=new Read(bb.toBytes(), null, "fused_"+tid, 0);
			processRefSequence(fusedRef, refList, ranges);
		}

		/** Optionally masks supplied bases, then selects indexed or brute alignment. */
		long processRefSequence(final Read ref, ArrayList<Read> originalRefs, IntList ranges){
			if(entropyMask){entropyMask(ref.bases, et);}
			return indexQueries ? processRefSequenceIndexed(ref, originalRefs, ranges) : processRefSequenceBrute(ref, originalRefs, ranges);
		}

		/** Builds a reference index per nonempty k bucket and processes both query strands. */
		long processRefSequenceIndexed(final Read ref, ArrayList<Read> originalRefs, IntList ranges){
			long sum=0;
			final float subrate=1-minid;
			IntHashMap2 seedMap=(useSeedMap ? new IntHashMap2() : null);

			for(int i=0; i<queryBuckets.size(); i++){
				ArrayList<Query> queries=queryBuckets.get(i);
				if(queries.isEmpty()){continue;}

				int currentK=kArray[i];
				IntListHashMap3 refIndex=buildReferenceIndex(ref.bases, currentK);

				// BRANCH: Fused vs Standard
				if(originalRefs==null) {
					for(Query q : queries){
						int count=0;
						int maxSubsQ=Math.min(maxSubs, (int)(q.length()*subrate));

						IntList seedHits=getSeedHits(q, refIndex, false, seedMap, ref.name());
						IntList hits=alignSparse(q.bases, ref.bases, maxSubsQ, q.maxClips, seedHits);
						count+=processHits(q, ref, hits, false, bsw, originalRefs, ranges, maxSubsQ);

						// Repeat for Reverse Strand
						seedHits=getSeedHits(q, refIndex, true, seedMap, ref.name());
						hits=alignSparse(q.rbases, ref.bases, maxSubsQ, q.maxClips, seedHits);

						count+=processHits(q, ref, hits, true, bsw, originalRefs, ranges, maxSubsQ);
						readsOutT+=count;
						basesOutT+=count*q.bases.length;
						sum+=count;
					}
				}else {
					for(Query q : queries){
						int count=0;
						int maxSubsQ=Math.min(maxSubs, (int)(q.length()*subrate));

						IntList seedHits=getSeedHits(q, refIndex, false, seedMap, ref.name());
						IntList hits=alignSparseFused(q.bases, ref.bases, maxSubsQ, q.maxClips, seedHits);
						count+=processHits(q, ref, hits, false, bsw, originalRefs, ranges, maxSubsQ);

						// Repeat for Reverse Strand
						seedHits=getSeedHits(q, refIndex, true, seedMap, ref.name());
						hits=alignSparseFused(q.rbases, ref.bases, maxSubsQ, q.maxClips, seedHits);

						count+=processHits(q, ref, hits, true, bsw, originalRefs, ranges, maxSubsQ);
						readsOutT+=count;
						basesOutT+=count*q.bases.length;
						sum+=count;
					}
				}
			}
			return sum;
		}

		/** Enumerates both strands without seed indexes and returns/counts accepted alignments. */
		long processRefSequenceBrute(final Read ref, ArrayList<Read> originalRefs, IntList ranges){
			long sum=0;
			final float subrate=1-minid;
			for(ArrayList<Query> bucket : queryBuckets){
				for(Query q : bucket){
					int count=0;
					int maxSubsQ=Math.min(maxSubs, (int)(q.length()*subrate));

					IntList hits=alignAllPositions(q.bases, ref.bases, maxSubsQ, q.maxClips);
					count+=processHits(q, ref, hits, false, bsw, originalRefs, ranges, maxSubsQ);

					hits=alignAllPositions(q.rbases, ref.bases, maxSubsQ, q.maxClips);
					count+=processHits(q, ref, hits, true, bsw, originalRefs, ranges, maxSubsQ);

					readsOutT+=count;
					basesOutT+=count*q.bases.length;
					sum+=count;
				}
			}
			long totalQueries=0;
			for(ArrayList<Query> b : queryBuckets){totalQueries+=b.size();}
			alignmentsT+=(totalQueries*(long)ref.length());
			return sum;
		}

		/** Selects list/map aggregation; the retained rname argument is unused. */
		private IntList getSeedHits(Query q, IntListHashMap3 refIndex, 
			boolean reverseStrand, IntHashMap2 hitCounts, String rname){
			if(useSeedMap){return getSeedHitsMap(q, refIndex, reverseStrand, hitCounts);}
			else{return getSeedHitsList(q, refIndex, reverseStrand);}
		}

		/** Applies optional prescan, then returns a new sorted unique qualifying-offset list or null. */
		private IntList getSeedHitsList(Query q, IntListHashMap3 refIndex, boolean reverseStrand){
			int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			if(queryKmers==null){return null;}
			final int minHits=Math.max(minSeedHits, q.minHits);
			if(prescan){
				int valid=prescan(q, refIndex, reverseStrand, minHits);
				if(valid<minHits){return null;}
			}
			IntList seedHits=new IntList();
			for(int i=0; i<queryKmers.length; i+=qStep){
				if(queryKmers[i]==-1){continue;}
				IntList positions=refIndex.get(queryKmers[i]);
				if(positions!=null){
					for(int j=0; j<positions.size; j++){
						int refPos=positions.array[j];
						int alignStart=refPos-i;
						seedHits.add(alignStart);
					}
				}
			}
			if(seedHits==null){return null;}
			seedHitsT+=seedHits.size();
			if(seedHits.size<minHits){return null;}
			if(seedHits.size>1 || minHits>1){
				seedHits.sort();
				seedHits.condenseMinCopies(minHits);
			}
			alignmentsT+=seedHits.size();
			return seedHits.isEmpty() ? null : seedHits;
		}

		/**
		 * Applies optional prescan and counts offsets in a scratch map cleared on the
		 * first decoded hit. Returns a new threshold-arrival-order list, or null if empty.
		 */
		private IntList getSeedHitsMap(Query q, IntListHashMap3 refIndex, 
			boolean reverseStrand, IntHashMap2 hitCounts){
			final int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			if(queryKmers==null){return null;}
			final int minHits=Math.max(minSeedHits, q.minHits);
			if(prescan){
				int valid=prescan(q, refIndex, reverseStrand, minHits);
				if(valid<minHits){return null;}
			}
			IntList seedHits=null;
			for(int i=0; i<queryKmers.length; i+=qStep){
				if(queryKmers[i]==-1){continue;}
				IntList positions=refIndex.get(queryKmers[i]);
				if(positions!=null){
					seedHitsT+=positions.size;
					if(seedHits==null){
						seedHits=new IntList();
						if(hitCounts==null){hitCounts=new IntHashMap2();}
						else{hitCounts.clear();}
					}
					for(int j=0; j<positions.size; j++){
						int alignStart=positions.array[j]-i;
						int newCount=hitCounts.increment(alignStart);
						if(newCount==minHits){seedHits.add(alignStart);}
					}
				}
			}
			alignmentsT+=(seedHits==null ? 0 : seedHits.size());
			return seedHits==null || seedHits.isEmpty() ? null : seedHits;
		}

		/**
		 * Scans valid qStep-spaced keys until the adjusted miss budget is exceeded.
		 * Returns observed present keys, possibly from a partial scan; no miss-free
		 * floor-to-threshold shortcut is used in this version.
		 */
		private int prescan(Query q, IntListHashMap3 refIndex, boolean reverseStrand, final int minHits){
			final int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			final int maxMisses=q.maxMisses-(minHits-q.minHits);
			if(queryKmers==null || maxMisses<0){return 0;}
			int misses=0, total=0;
			for(int i=0; i<queryKmers.length && misses<=maxMisses; i+=qStep){
				if(queryKmers[i]==-1){continue;}
				total++;
				boolean hit=refIndex.containsKey(queryKmers[i]);
				misses+=(hit ? 0 : 1);
			}
			return total-misses;
		}

		/** Accepted reference roots plus mates; query-loading counts belong to the outer instance. */
		protected long readsProcessedT=0;
		/** Bases in accepted reference roots plus mates, before optional masking. */
		protected long basesProcessedT=0;
		/** Indexed candidates, or brute-mode query-count times reference length estimate. */
		protected long alignmentsT=0;
		/** Decoded seed occurrences before offset-threshold aggregation. */
		protected long seedHitsT=0;
		/** Accepted alignment count, including repeated alignments of one query. */
		protected long readsOutT=0;
		/** Full query lengths summed over accepted alignments, including clipped bases. */
		protected long basesOutT=0;
		/** Set only after normal processing completion; merged after joining the worker. */
		boolean success=false;

		/** Shared reference reader, owned and closed by process(). */
		private final Streamer cris;
		/** Shared alignment writer, or null when output is disabled. */
		private final ByteStreamWriter bsw;
		/** Optional separate shared header writer. */
		private final SamHeaderWriter shw;
		/** Prepared shared buckets; each query's alignment counter remains mutable. */
		private final ArrayList<ArrayList<Query>> queryBuckets;
		/** Per-worker copies of process-wide scoring options. */
		final int maxSubs;
		final float minid;
		/** Worker identifier used in temporary fused-reference names. */
		final int tid;
		/** Owned masking scratch, initialized only when entropy masking is enabled. */
		EntropyTracker et;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Query input paths; in2 is optional. */
	private String in1=null;
	private String in2=null;
	/** Alignment output and optional separate SAM-header output paths. */
	private String out1=null;
	private String headerOut=null;
	/** Query-input and alignment-output format overrides. */
	private String extin=null;
	private String extout=null;
	/** Reference input path, opened after query materialization. */
	String refFile=null;
	/** Process-wide scoring options, copied into workers. */
	static int maxSubs=5;
	static float minid=0.85f;

	/** Descending indexed k values after construction, or the sole zero value in brute mode. */
	int[] kArray=new int[] {8,9,10,12,14};

	/** Central masked bases in k-mer keys. */
	int midMaskLen=1;
	/** Enables query k-mer preparation and reference indexes. */
	boolean indexQueries=true;
	/** Enables the key-presence screen before position decoding. */
	boolean prescan=true;
	/** Query-start and power-of-two reference-start strides; at least one is one. */
	int qStep=1;
	int rStep=1;
	/** Maximum stride captured during construction. */
	final int kStep;

	/** Floor applied to each query's modeled seed threshold. */
	int minSeedHits=1;
	/** Probability forwarded to MinHitsCalculator2 during serial query setup. */
	private float minHitsProb=0.999f;
	/** Selects count-map aggregation rather than sorting decoded offsets. */
	boolean useSeedMap=false;

	/** Enables reference-root concatenation with trailing padding after each root. */
	boolean fuse=true;
	int padding=128;
	/** Reference-reader buffer-data request, additionally capped in standard mode. */
	int targetChunkSize=1000000;
	/** Root length filters; a passing query root also admits its mate. */
	int minQLen=1;
	int minRLen=1;

	/** Query-loading plus reference-processing totals, counting roots and mates. */
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	/** Sums of corresponding worker counters, including heuristic work estimates. */
	protected long alignmentCount=0;
	protected long seedHitCount=0;
	protected long readsOut=0;
	protected long basesOut=0;

	/** Limit independently forwarded to query/reference readers; -1 is unlimited. */
	private long maxReads=-1;

	/** Formats prepared by construction; makeCris prepares the reference format later. */
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	private final FileFormat ffout1;
	private final FileFormat ffheader;

	/** Lock supplied to the accumulator's caller; does not serialize arbitrary instance use. */
	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/** Diagnostic stream selected by PreParser. */
	private PrintStream outstream=System.err;
	/** Process-wide diagnostic verbosity. */
	public static boolean verbose=false;
	/** Accumulated worker/input/statistics/output-finalization error state. */
	public boolean errorState=false;
	/** Output policies also copied into shared ReadStats during construction. */
	private boolean overwrite=true;
	private boolean append=false;
	
	/** Masks standard reference bases or the temporary fused copy before alignment. */
	private boolean entropyMask=false;
	/** Window size for sliding window entropy/complexity analysis */
	private int entropyWindow=80;
	/** Entropy threshold below which regions are masked */
	private float entropyCutoff=0.70f;
	/** K-mer length used by the entropy tracker. */
	private int entropyK=4;
	/** Retained unused field; active repeat-key exclusion uses Query.blacklistRepeatLength. */
	private int repeatLen=4;
	
	
}
