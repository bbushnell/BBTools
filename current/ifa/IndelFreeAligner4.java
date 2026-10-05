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
 * Version 4 indel-free aligner using {@link PackedIndex4} reference indexes.
 * Materializes queries in descending-k buckets, then aligns both strands against
 * reference batches taken by workers from one shared Streamer. Indexed mode uses
 * masked two-bit k-mer keys, heuristic prescans and candidate substitution scoring;
 * brute mode enumerates offsets. Alignment helpers conditionally delegate to SIMD.
 * <p>
 * Workers own their reference indexes, seed scratch and counters. Query bases and
 * derived arrays are shared for reading; each query's atomic alignment counter
 * chooses the first counted alignment as primary. Output order is not stabilized.
 * Optional fusion copies a reference batch with N padding and converts accepted
 * positions back to an original reference using the query center.
 * <p>
 * Construction and query loading change process-wide I/O, Query and calculator
 * settings. Instances are not independent concurrent configurations. The normal
 * process path joins workers and writers and restores selected buffer/validation
 * settings; it does not provide a general finally-based restoration guarantee.
 * @author Brian Bushnell
 * @contributor Isla, Amber
 * @date June 3, 2025
 */
public class IndelFreeAligner4 implements Accumulator<IndelFreeAligner4.ProcessThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs one CLI invocation and closes its diagnostic stream on normal completion. */
	public static void main(String[] args){
		Timer t=new Timer();
		IndelFreeAligner4 x=new IndelFreeAligner4(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/**
	 * Parses CLI options, configures shared library settings and prepares file formats.
	 * Indexed k values are sorted descending; brute mode uses one bucket with k=0.
	 * @param args BBTools flag=value arguments, including query input and reference
	 */
	public IndelFreeAligner4(String[] args){

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

		if(kArray==null || kArray.length<1 || (kArray.length==1 && kArray[0]==0)){indexQueries=false;}
		if(indexQueries){
			Arrays.sort(kArray);
			Tools.reverseInPlace(kArray);
		}else{
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

	/** Parses tool options and delegates common options to Parser; also changes shared statics. */
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
				if(minid>1){minid/=100f;}
			}else if(a.equals("hits") || a.equals("minhits") || a.equals("seedhits")){
				minSeedHits=Math.max(1, Integer.parseInt(b));
			}else if(a.equals("minprob") || a.equals("minhitsprob")){
				minHitsProb=Float.parseFloat(b);
			}else if(a.equals("iterations")){
				MinHitsCalculator3.iterations=MinHitsCalculator2.iterations=Parse.parseIntKMG(b);
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
					for(int x=0; x<range; x++){kArray[x]=x+temp[0];}
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

	/** Resolves alternate compressed extensions for the two query input paths. */
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

	/** Checks the existing assertion-based parameter constraints; not an exhaustive input validator. */
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
	 * Loads queries, starts reference/output services and joins the alignment workers.
	 * Then closes input, waits for output/header completion and reports accumulated
	 * errors. Buffer sizes and Read validation mode are restored along this normal
	 * path, before the final error-state exception; earlier exceptions can bypass it.
	 * Processed counters include loaded query mates and accepted reference pairs.
	 * @param t Timer started by the caller, stopped after stream finalization
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
		Shared.setBufferData(Math.min(targetChunkSize, (fuse ? targetChunkSize : 400000)));
		if(fuse){Shared.setBufferLen(Math.max(oldBL, 2000));}

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
		if(shw!=null){shw.poisonAndWait();}

		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		Shared.setBufferData(oldBD);
		Shared.setBufferLen(oldBL);

		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		outstream.println(Tools.things("Alignments", alignmentCount, 8));
		outstream.println(Tools.things("Seed Hits ", seedHitCount, 8));
		outstream.println(Tools.things("Prescans  ", prescans, 8));
		outstream.println(Tools.things("Postscans ", postscans, 8));
		outstream.println(Tools.things("Postscan2 ", postscan2, 8));

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts workers sharing the services/buckets, waits for them and accumulates success. */
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
			prescans+=pt.prescansT;
			prescanFails+=pt.prescanFailsT;
			postscans+=pt.postscansT;
			postscan2+=pt.postscan2T;
			errorState|=(!pt.success);
		}
	}

	/** Reports the errors accumulated so far; this is not a completion test. */
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

	/**
	 * Materializes input reads, installs shared Query calculators and builds query buckets.
	 * Counts every loaded root and mate before filtering. A root meeting minQLen admits
	 * both itself and its mate; the mate has no separate length check here. Indexed
	 * queries select their calculator's bucket; unindexed queries use bucket zero.
	 * @param ff1 Primary query input
	 * @param ff2 Optional separate mate input
	 * @return Owned bucket lists of Query objects, completed before workers start
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

		for(Read r : reads){
			readsProcessed+=r.pairCount();
			basesProcessed+=r.pairLength();

			if(r.length()>=minQLen){
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

	/** Wraps one read's borrowed bases/qualities in Query and selects its calculator bucket. */
	private void addReadToBucket(Read r, ArrayList<ArrayList<Query>> buckets){
		Query q=new Query(r.id, 0, r.bases, r.quality);
		if(q.calculatorIndex>=0){
			buckets.get(q.calculatorIndex).add(q);
		}else{
			buckets.get(0).add(q); 
		}
	}

	/**
	 * Scores supplied zero-based offsets against one reference, preserving candidate order.
	 * Overhangs use alignClipped; interior offsets delegate to Vector.align.
	 * @param query Query bases for the selected strand
	 * @param ref Reference bases
	 * @param maxSubs Maximum substitution/paid-clipping cost
	 * @param maxClips Total free overhang allowance in bases
	 * @param seedHits Borrowed candidate offsets, possibly null; consumed without mutation
	 * @return New accepted-offset list, or null when no candidate passes
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
		return results;
	}
	
	/**
	 * Scores candidate offsets in a fused reference, preserving candidate order.
	 * Interior offsets delegate to Vector.alignFused; outer overhangs use alignClipped.
	 * Original-reference selection and final clipping checks occur later in processHits.
	 * @param query Query bases for the selected strand
	 * @param ref Fused reference bases with padding
	 * @param maxSubs Maximum substitution/paid-clipping cost
	 * @param maxClips Total free clipping allowance in bases
	 * @param seedHits Borrowed candidate offsets, possibly null; consumed without mutation
	 * @return New accepted-offset list, or null when no candidate passes
	 */
	public static IntList alignSparseFused(byte[] query, byte[] ref, int maxSubs, int maxClips, IntList seedHits){
		if(seedHits==null || seedHits.isEmpty()){return null;}
		IntList results=null;
		final int rLen=ref.length;
		final int qLen=query.length;

		for(int i=0; i<seedHits.size; i++){
			int rStart=seedHits.array[i];
			int subs;
			if(rStart<0 || rStart>rLen-qLen){
				subs=alignClipped(query, ref, maxSubs, maxClips, rStart);
			}else{
				subs=Vector.alignFused(query, ref, maxSubs, maxClips, rStart);
			}

			if(subs<=maxSubs){
				if(results==null){results=new IntList(4);}
				results.add(rStart);
			}
		}
		return results;
	}

	/**
	 * Enumerates offsets for non-null query/reference arrays and nonnegative budgets.
	 * SIMD dispatch precedes the scalar path and uses its own candidate enumeration.
	 * Scalar candidates require nonempty overlap; clipping beyond maxClips consumes
	 * maxSubs. Scalar accepted offsets are ascending, or null if none pass.
	 * @param query Query bases for the selected strand
	 * @param ref Reference bases
	 * @param maxSubs Maximum substitution/paid-clipping cost
	 * @param maxClips Total free overhang allowance in bases
	 * @return Accepted offsets from the selected implementation
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
	 * Stops after exceeding maxSubs, so a rejected score may be only a partial count.
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
	 * Scores an offset with left/right overhangs, charging clips beyond maxClips.
	 * A placement with no overlap returns query.length; otherwise scoring stops after
	 * exceeding maxSubs. This pointwise helper does not enforce the candidate domain.
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
	 * Rechecks candidates against an original reference and builds SAM alignments.
	 * Fused candidates select a reference by query-center membership in [start,stop)
	 * ranges, then convert to local offsets. Substitution and clipping costs are
	 * checked separately; NM counts substitutions only. The shared query counter
	 * assigns primary status in arrival order. Null output still builds/counts hits
	 * and advances that counter. Borrowed hits are consumed without modification.
	 * @param q Shared query, including original and reverse-complement bases
	 * @param ref Standard reference or temporary fused reference
	 * @param hits Candidate offsets, possibly null
	 * @param reverseStrand Whether to score reverse-complement query bases
	 * @param bsw Started shared output writer, or null to discard serialized output
	 * @param originalRefs Original references for fusion, or null in standard mode
	 * @param ranges Paired fused start/inclusive and stop/exclusive coordinates
	 * @param maxSubsQ Effective query-specific substitution/paid-clipping limit
	 * @return Number of accepted SAM alignments, not distinct queries
	 */
	static int processHits(Query q, Read ref, IntList hits, boolean reverseStrand,
		ByteStreamWriter bsw, ArrayList<Read> originalRefs, IntList ranges, int maxSubsQ){

		if(hits==null || hits.size()==0){return 0;}
		ByteBuilder bb=new ByteBuilder();
		ByteBuilder match=new ByteBuilder(q.bases.length);

		int added=0; 
		final int qLen=q.length();
		final int qHalf=qLen/2;

		for(int i=0; i<hits.size(); i++){
			int start=hits.get(i);
			Read realRef=ref;
			int rStart=start;

			if(originalRefs!=null){
				int center=start+qHalf;
				int idx=Arrays.binarySearch(ranges.array, 0, ranges.size, center);
				if(idx<0){idx=-idx-2;}
				
				if(idx<0 || (idx&1)==1){continue;}
				
				int refIdx=idx/2;
				if(refIdx>=originalRefs.size()){continue;}

				realRef=originalRefs.get(refIdx);
				rStart=start-ranges.get(idx);
			}

			final int rLen=realRef.length();
			final int rStop=rStart+qLen;

			final int leftClip=Math.max(0, -rStart);
			final int rightClip=Math.max(0, rStop-rLen);
			final int totalClip=leftClip+rightClip;
			
			int clipPenalty=Math.max(0, totalClip-q.maxClips);
			if(clipPenalty>maxSubsQ){continue;}

			byte[] querySeq=reverseStrand ? q.rbases : q.bases;
			toMatch(querySeq, realRef.bases, rStart, match.clear());

			int subs=0;
			for(int m=0; m<match.length(); m++){
				final byte b=match.get(m);
				//n comprehension: toMatch() only ever emits 'm' (match), 'S' (in-bounds mismatch), or 'C' (out-of-bounds clip) —
				//n never 'N'. So the 'N' arm is dead here (harmless); 'S' is the only thing counted. 'C' is deliberately NOT a sub:
				//n clip cost is accounted separately via clipPenalty above, so overhang isn't double-charged. Correct.
				if(b=='S' || b=='N'){subs++;}
			}

			if((subs+clipPenalty)>maxSubsQ){continue;}

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
	
	/** Masks low-complexity bases in place, returning the count reported by BBMask. */
	private int entropyMask(byte[] bases, EntropyTracker et){
		if(bases==null || bases.length==0){return 0;}
		BitSet bs=new BitSet(bases.length);
		if(et==null){et=new EntropyTracker(entropyK, entropyWindow, false, entropyCutoff, true);}
		BBMask.maskLowEntropy(bases, bs, et);
		return BBMask.maskBases(bases, bs, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** One shared-reader consumer with owned scratch/counters; borrows query buckets and writers. */
	class ProcessThread extends Thread {

		/** Borrows the services and query buckets; copies effective scoring settings for this worker. */
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

		/** Initializes scratch and marks success only after processInner returns normally. */
		@Override
		public void run(){
			synchronized(this){
				if(entropyMask){et=new EntropyTracker(entropyK, entropyWindow, false, entropyCutoff, true);}
				hitsList=new IntList(100); // Reused
				processInner();
				success=true;
			}
		}

		/** Consumes shared-reader batches until its null terminal. */
		void processInner(){
			for(ListNum<Read> ln=cris.nextList(); ln!=null; ln=cris.nextList()){processList(ln);}
		}

		/**
		 * Removes short reference roots in place, validates/counts survivors and queues
		 * header entries under the original batch ID, including an empty header batch.
		 * Dispatches the compacted list to standard or fused alignment.
		 */
		void processList(ListNum<Read> ln){
			final ArrayList<Read> refList=ln.list;
			final ArrayList<StringNum> alsn=(shw==null ? null : new ArrayList<StringNum>(refList.size()));
			int removed=0;
			for(int idx=0; idx<refList.size(); idx++){
				final Read ref=refList.get(idx);
				//FIXED IFA-010: condenseStrict removes nulls before either standard or fused alignment.
				if(ref.length()<minRLen){refList.set(idx, null); removed++; continue;}
				if(!ref.validated()){ref.validate(true);}
				if(alsn!=null){alsn.add(new StringNum(ref.name(), ref.length()));}

				final int initialLength1=ref.length();
				final int initialLength2=ref.mateLength();
				readsProcessedT+=ref.pairCount();
				basesProcessedT+=initialLength1+initialLength2;
			}
			if(shw!=null){shw.add(new ListNum<StringNum>(alsn, ln.id));}
			if(removed>0){Tools.condenseStrict(refList);}

			if(fuse){processListFused(refList);}
			else{processListStandard(refList);}
		}

		/** Aligns each reference root independently; entropy masking may mutate its bases. */
		void processListStandard(ArrayList<Read> refList){
			for(int idx=0; idx<refList.size(); idx++){
				final Read ref=refList.get(idx);
				processRefSequence(ref, null, null);
			}
		}

		/**
		 * Copies roots plus trailing N padding into one reference and records [start,stop)
		 * ranges. Falls back to standard processing when the computed length exceeds int.
		 */
		void processListFused(ArrayList<Read> refList){
			if(refList.isEmpty()){return;}

			long totalLen=0;
			for(Read r : refList){
				totalLen+=r.length()+padding;
			}
			if(totalLen>Integer.MAX_VALUE){
				processListStandard(refList);
				return;
			}

			ByteBuilder bb=new ByteBuilder((int)totalLen);
			IntList ranges=new IntList(refList.size()*2);

			for(Read r : refList){
				ranges.add(bb.length()); 
				bb.append(r.bases);
				ranges.add(bb.length()); 
				for(int i=0; i<padding; i++){bb.append('N');}
			}

			Read fusedRef=new Read(bb.toBytes(), null, "fused_"+tid, 0);
			processRefSequence(fusedRef, refList, ranges);
		}

		/** Optionally masks the supplied bases, then dispatches indexed or brute alignment. */
		long processRefSequence(final Read ref, ArrayList<Read> originalRefs, IntList ranges){
			if(entropyMask){entropyMask(ref.bases, et);}
			return indexQueries ? processRefSequenceIndexed(ref, originalRefs, ranges) : processRefSequenceBrute(ref, originalRefs, ranges);
		}

		/**
		 * Builds one packed reference index per nonempty k bucket and processes both strands.
		 * Each borrowed seed list is fully consumed before the next call reuses its buffer.
		 * Returns accepted alignment count and updates this worker's output counters.
		 */
		long processRefSequenceIndexed(final Read ref, ArrayList<Read> originalRefs, IntList ranges){
			long sum=0;
			final float subrate=1-minid;
			IntHashMap2 seedMap=(useSeedMap ? new IntHashMap2() : null);

			for(int i=0; i<queryBuckets.size(); i++){
				ArrayList<Query> queries=queryBuckets.get(i);
				if(queries.isEmpty()){continue;}

				int currentK=kArray[i];
				PackedIndex4 refIndex=new PackedIndex4(ref.bases, currentK, midMaskLen, rStep);

				for(Query q : queries){
					int count=0;
					int maxSubsQ=Math.min(maxSubs, (int)(q.length()*subrate));

					//n comprehension (reused-buffer safety, non-obvious): getSeedHits returns the per-thread REUSED IntList hitsList
					//n (cleared on entry), so the forward and reverse calls alias the SAME list. It's safe ONLY because alignSparse
					//n FULLY consumes seedHits into a fresh results IntList on the line right after each getSeedHits, BEFORE the next
					//n getSeedHits (reverse) clears/refills hitsList. Order is load-bearing: forward-get -> forward-consume -> reverse-
					//n get -> reverse-consume. hitsList is a per-ProcessThread field (run() inits it), so no cross-thread aliasing either.
					IntList seedHits=getSeedHits(q, refIndex, false, seedMap);
					if(originalRefs==null){
						IntList hits=alignSparse(q.bases, ref.bases, maxSubsQ, q.maxClips, seedHits);
						count+=processHits(q, ref, hits, false, bsw, originalRefs, ranges, maxSubsQ);

						seedHits=getSeedHits(q, refIndex, true, seedMap);
						hits=alignSparse(q.rbases, ref.bases, maxSubsQ, q.maxClips, seedHits);
						count+=processHits(q, ref, hits, true, bsw, originalRefs, ranges, maxSubsQ);
					}else{
						IntList hits=alignSparseFused(q.bases, ref.bases, maxSubsQ, q.maxClips, seedHits);
						count+=processHits(q, ref, hits, false, bsw, originalRefs, ranges, maxSubsQ);

						seedHits=getSeedHits(q, refIndex, true, seedMap);
						hits=alignSparseFused(q.rbases, ref.bases, maxSubsQ, q.maxClips, seedHits);
						count+=processHits(q, ref, hits, true, bsw, originalRefs, ranges, maxSubsQ);
					}
					readsOutT+=count;
					basesOutT+=count*q.bases.length;
					sum+=count;
				}
			}
			return sum;
		}

		/** Enumerates both strands without seed indexes; returns and counts accepted alignments. */
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
		
		/**
		 * Screens k-mer presence against the adjusted miss allowance in two interleaved passes.
		 * Returns a heuristic hit count: a miss-free pass returns at least minHits without
		 * enumerating remaining starts. This is not an exact seed count or sensitivity proof.
		 */
		private int prescan(Query q, PackedIndex4 refIndex, boolean reverseStrand, final int minHits){
			final int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			// Adjust maxMisses: if we require more hits than the query's internal baseline, we have less wiggle room.
			final int maxMisses=q.maxMisses-(minHits-q.minHits);
			
			if(queryKmers==null || maxMisses<0){return 0;}

			int misses=0, total=0, hits=0;
			// Direct access to the map is safe and fastest
			final map.IntHashMap2 map=refIndex.map;

//			final int step=2*qStep;
//			for(int start=0; start<=qStep && misses<=maxMisses; start+=qStep) {
//				for(int i=start; i<queryKmers.length && misses<=maxMisses; i+=step){
//					int kmer=queryKmers[i];
//					if(kmer==-1){continue;}
//
//					total++;//This could be built into query...
//					misses+=(map.containsKey(kmer) ? 0 : 1);
//					// Fail Fast: If we have already missed too many, success is impossible.
//				}
//				//Note:  If there are 0 misses after the first pass,
//				//second pass is unnecessary since kmers overlap.
//			}
			
			final int step=2*qStep;
			for(int start=0; start<=qStep && misses<=maxMisses && hits<minHits; start+=qStep) {
				for(int i=start; i<queryKmers.length && misses<=maxMisses && hits<minHits; i+=step){
					int kmer=queryKmers[i];
					if(kmer==-1){continue;}

					final int found=map.containsKeyBinary(kmer);
					hits+=found;
					misses+=(found^1);
				}
				//Note:  If there are 0 misses after the first pass,
				//second pass is unnecessary since kmers overlap.
				if(misses<1) {return Math.max(hits,  minHits);}
			}
//			return total-misses;
			return hits;
		}

		/**
		 * Applies optional prescan and decodes offsets meeting the effective hit threshold.
		 * Returns null on missing query k-mers or prescan rejection; otherwise returns the
		 * worker's reused hitsList, possibly empty. Consume it before another seed call.
		 */
		private IntList getSeedHits(Query q, PackedIndex4 refIndex, boolean reverseStrand, IntHashMap2 hitCounts){
			int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			if(queryKmers==null){return null;}
			final int minHits=Math.max(minSeedHits, q.minHits);
			
			if(prescan){
				prescansT++;
				int possibleHits=prescan(q, refIndex, reverseStrand, minHits);
				if(possibleHits<minHits){
					prescanFailsT++;
					return null;
				}
			}
			postscansT++;
			final IntList seeds;
			if(useSeedMap){seeds=getSeedHitsMap(queryKmers, refIndex, minHits, hitCounts);}
			else{seeds=getSeedHitsList(queryKmers, refIndex, minHits);}
			
			alignmentsT+=(seeds.size);
			postscan2T+=(seeds.size>0 ? 1 : 0);
			return seeds;
		}

		/**
		 * Counts decoded offsets in a cleared scratch map, adding each when it reaches minHits.
		 * Result is the reused hitsList in threshold-arrival order, not necessarily sorted.
		 */
		private IntList getSeedHitsMap(int[] queryKmers, PackedIndex4 refIndex, 
			int minHits, IntHashMap2 hitCounts){
			final IntList seedHits=hitsList; 
			seedHits.clear();

			if(hitCounts==null){hitCounts=new IntHashMap2();}
			else{hitCounts.clear();}

			final int[] positions=refIndex.positions;
			final IntHashMap2 map=refIndex.map;

			for(int i=0; i<queryKmers.length; i+=qStep){
				final int kmer=queryKmers[i];
				if(kmer==-1){continue;}

				// PackedIndex4 Logic:
				// -1       = Missing
				// < -1     = Singleton (RefPos | MIN_VALUE)
				// >= 0     = List Head Pointer
				//n studied praise: PackedIndex4 packs three states into one int per kmer with zero extra memory — val==-1 means
				//n ABSENT, val<-1 means a SINGLETON with its ref position stored in the low 31 bits (val&Integer.MAX_VALUE), and
				//n val>=0 means a HEAD POINTER into positions[], where each entry's sign bit is a STOP bit (entry<0 => last). So the
				//n common singleton case costs no positions[] slot at all, and the multi-hit list needs no length field. Verified the
				//n decode matches getSeedHitsList's identical scheme; the #001 guard is exactly what keeps the -1 state out of the else.
				final int val=map.get(kmer);

				if(val==-1){continue;}//FIXED [ifa/IndelFreeAligner4#001]: -1 = Missing (kmer absent from reference); skip. Same fix as getSeedHitsList.

				if(val<-1){
					// Singleton Case
					seedHitsT++;
					int refPos=val&Integer.MAX_VALUE;
					int alignStart=refPos-i;
					int newCount=hitCounts.increment(alignStart);
					if(newCount==minHits){seedHits.add(alignStart);}
				}else{
					assert(val>=0) : "Illegal value: "+kmer+","+val;
					// Multi-Hit Case
					int ptr=val;
					while(true){
						seedHitsT++;
						int entry=positions[ptr];
						int refPos=entry&Integer.MAX_VALUE;
						int alignStart=refPos-i;

						int newCount=hitCounts.increment(alignStart);
						if(newCount==minHits){seedHits.add(alignStart);}

						// Stop bit check
						if(entry<0){break;}
						ptr++;
					}
				}
			}
			return seedHits;
		}
		
		/**
		 * Decodes sampled k-mers, sorts offsets and retains those with enough copies.
		 * Missing keys are skipped; singleton and stop-bit position lists share the same
		 * reference-start minus query-start coordinate calculation.
		 * @param queryKmers K-mer keys by query start; -1 entries are skipped
		 * @param refIndex Borrowed packed reference index
		 * @param minHits Required occurrences of one offset
		 * @return Sorted unique offsets in reused hitsList; overwritten by the next seed call
		 */
		private IntList getSeedHitsList(int[] queryKmers, PackedIndex4 refIndex, int minHits){
			final IntList seedHits=hitsList; 
			seedHits.clear();
			final int[] positions=refIndex.positions;
			final IntHashMap2 map=refIndex.map;

			for(int i=0; i<queryKmers.length; i+=qStep){
				final int kmer=queryKmers[i];
				if(kmer==-1){continue;}
				
				// Fetch from map: 
				// -1       = Missing
				// < -1     = Singleton (RefPos | MIN_VALUE)
				// >= 0     = List Head (Index in positions array)
				final int val=map.get(kmer);

				//FIXED [ifa/IndelFreeAligner4#001]: -1 = Missing (kmer absent from the reference index); skip it.
				//It was falling into the else branch below and tripping assert(val>=0) on any real read with a
				//reference-absent kmer -> universal crash. Matches IndelFreeAligner3:882 / IndelFreeAligner2:923.
				if(val==-1){continue;}

				if(val<-1){
					// Singleton Case:
					// Decode RefPos from high bit (Mask off sign bit)
					int refPos=val&Integer.MAX_VALUE;
					int alignStart=refPos-i;
					seedHits.add(alignStart);
				}else{
					assert(val>=0) : "Illegal value: "+kmer+","+val;
					// Multi-Hit Case:
					// 'val' is the direct pointer to the positions array head.
					int ptr=val;

					while(true){
						int entry=positions[ptr];

						// Mask off the Stop Bit (sign bit) to get the real RefPos
						int refPos=entry&Integer.MAX_VALUE;
						int alignStart=refPos-i;

						// Add hit. Note: use add(), not addUnchecked(), as capacity is unknown.
						seedHits.add(alignStart);

						// Check Stop Bit: If entry is negative, this is the last item.
						if(entry<0){break;}
						ptr++;
					}
					//FIXED [ifa/IndelFreeAligner4#002]: deleted a duplicated byte-identical while-loop that was here.
					//After the loop above breaks (entry<0, before ptr++), ptr still points at the stop entry, so the
					//duplicate re-read positions[ptr] and re-added the last hit. Copy-paste error; absent from getSeedHitsMap.
				}
			}
			
			seedHitsT+=seedHits.size;
			if(seedHits.size>1 || (minHits>1 && seedHits.size>0)){
				seedHits.sort();
				seedHits.condenseMinCopies(minHits);
			}
			return seedHits;
		}

		/** Accepted reference roots plus mates; query-loading counts are added by the outer instance. */
		protected long readsProcessedT=0;
		/** Bases in accepted reference roots plus mates, before optional masking. */
		protected long basesProcessedT=0;
		/** Indexed candidate count, or brute-mode query-count times reference length estimate. */
		protected long alignmentsT=0;
		/** Decoded seed occurrences before threshold aggregation. */
		protected long seedHitsT=0;
		/** Accepted SAM alignment count, including repeated alignments of one query. */
		protected long readsOutT=0;
		/** Sum of full query lengths over accepted alignments, including clipped bases. */
		protected long basesOutT=0;
		/** Strand/bucket prescan calls. */
		protected long prescansT=0;
		/** Prescans returning below the required threshold. */
		protected long prescanFailsT=0;
		/** Seed-decoding calls that passed or skipped prescan. */
		protected long postscansT=0;
		/** Seed-decoding calls returning at least one candidate. */
		protected long postscan2T=0;
		/** Set only after normal completion; merged after joining this worker. */
		boolean success=false;

		/** Shared reference batch reader, owned and closed by process(). */
		private final Streamer cris;
		/** Shared alignment writer; null when alignment output is disabled. */
		private final ByteStreamWriter bsw;
		/** Shared separate header writer, or null. */
		private final SamHeaderWriter shw;
		/** Fully prepared shared buckets; per-query alignment counters remain mutable. */
		private final ArrayList<ArrayList<Query>> queryBuckets;
		/** Per-worker copies of the global scoring options. */
		final int maxSubs;
		final float minid;
		/** Worker identifier used in temporary fused reference names. */
		final int tid;
		/** Owned entropy scratch, initialized only when masking is enabled. */
		EntropyTracker et;
		/** Owned scratch returned by both seed decoders and cleared on each decoder entry. */
		IntList hitsList;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Query input paths; in2 is optional. */
	private String in1=null;
	private String in2=null;
	/** Alignment output and optional separate SAM header output paths. */
	private String out1=null;
	private String headerOut=null;
	/** Query-input and alignment-output format overrides. */
	private String extin=null;
	private String extout=null;
	/** Reference input path, opened after query loading. */
	String refFile=null;
	/** Process-wide scoring defaults/options, copied into each worker. */
	static int maxSubs=5;
	static float minid=0.85f;

	/** Descending k values after construction; the sole value is zero in brute mode. */
	int[] kArray=new int[] {8,9,10,12,14};

	/** Central masked bases in packed keys. */
	int midMaskLen=1;
	/** Enables query k-mer preparation and per-bucket reference indexes. */
	boolean indexQueries=true;
	/** Enables the presence-based heuristic before seed decoding. */
	boolean prescan=true;
	/** Query seed-start stride and power-of-two reference index stride; at least one is one. */
	int qStep=1;
	int rStep=1;
	/** Effective maximum stride captured during construction. */
	final int kStep;

	/** Floor applied to each query's modeled seed threshold. */
	int minSeedHits=1;
	/** Probability option forwarded to MinHitsCalculator2 during serial query setup. */
	private float minHitsProb=0.999f;
	/** Selects count-map aggregation instead of sorting decoded offsets. */
	boolean useSeedMap=false;

	/** Enables batch concatenation with trailing padding after every reference root. */
	boolean fuse=true;
	int padding=128;
	/** Reference-reader buffer-data request, additionally capped in standard mode. */
	int targetChunkSize=1000000;
	/** Root-read length filters; query mates are admitted with a passing root. */
	int minQLen=1;
	int minRLen=1;

	/** Loaded query roots/mates plus accepted reference roots/mates; not reference-only totals. */
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	/** Sums of the corresponding worker counters, including heuristic work estimates. */
	protected long alignmentCount=0;
	protected long seedHitCount=0;
	protected long readsOut=0;
	protected long basesOut=0;
	protected long prescans=0;
	protected long prescanFails=0;
	protected long postscans=0;
	protected long postscan2=0;

	/** Reader limit forwarded independently to query loading and reference streaming; -1 is unlimited. */
	private long maxReads=-1;

	/** Formats prepared during construction; reference format is prepared by makeCris. */
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	private final FileFormat ffout1;
	private final FileFormat ffheader;

	/** Lock supplied to the accumulator's caller; not a lock around arbitrary instance use. */
	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/** Diagnostic stream selected by PreParser. */
	private PrintStream outstream=System.err;
	/** Shared diagnostic verbosity, changed by CLI parsing. */
	public static boolean verbose=false;
	/** Accumulated worker/input/statistics/output-finalization error state. */
	public boolean errorState=false;
	/** Output policies passed to FileFormat; also copied to shared ReadStats during construction. */
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
	
}
