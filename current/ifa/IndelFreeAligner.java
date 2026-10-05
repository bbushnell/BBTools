package ifa;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.IntHashMap2;
import map.IntListHashMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.SIMDAlignByte;
import simd.Vector;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamHeader;
import stream.SamHeaderWriter;
import stream.SamLine;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import structures.StringNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Original indel-free aligner using IntListHashMap reference indexes.
 * Materializes query roots and mates into descending-k buckets, then aligns both
 * strands against reference roots consumed from a shared ConcurrentReadInputStream.
 * Workers return each consumed batch and its terminal container to that reader.
 * Indexed mode uses masked two-bit k-mer keys and substitution-scored candidates;
 * brute mode enumerates offsets. Alignment helpers conditionally delegate to SIMD.
 * This version has no reference fusion, entropy masking or explicit length filters.
 * <p>
 * Workers own their reference indexes, freshly allocated seed lists and counters.
 * Prepared query arrays are shared for reading; the atomic query alignment counter
 * assigns primary status to the first counted alignment. Output order is not
 * stabilized. Alignment output is headerless; an optional separate writer collects
 * reference names and lengths for the SAM header.
 * <p>
 * Scoring options belong to this instance, but construction and query loading also
 * change shared I/O, Query and calculator settings. Instances are therefore not
 * independent concurrent configurations. The normal process path joins workers
 * and writers and restores selected settings; earlier exceptions can bypass that
 * restoration because it is not protected by a finally block.
 * @author Brian Bushnell
 * @contributor Isla, Amber
 * @date June 2, 2025
 */
public class IndelFreeAligner implements Accumulator<IndelFreeAligner.ProcessThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs one CLI invocation and closes its diagnostic stream on normal completion. */
	public static void main(String[] args){
		Timer t=new Timer();
		IndelFreeAligner x=new IndelFreeAligner(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/**
	 * Parses options, configures shared libraries and prepares query/output formats.
	 * Indexed k values are sorted descending; brute mode uses one bucket with k=0.
	 * @param args BBTools flag=value arguments, including query input and reference
	 */
	public IndelFreeAligner(String[] args){

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

	/** Parses tool options and common Parser options; k accepts a scalar or comma-separated list. */
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
				MinHitsCalculator2.iterations=Parse.parseIntKMG(b);
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
				}else if(b.indexOf(',')>-1){ //else required: b is null in the branch above; a bare if would NPE on b.indexOf. Matches IFA4.
					int[] temp=Parse.parseIntArray(b, ",");
					kArray=temp;
				}else{
					kArray=new int[] {Integer.parseInt(b)};
				}
				indexQueries=(kArray!=null && kArray.length>0 && kArray[0]>0);
			}else if(a.equals("qstep") || a.equals("step") || a.equals("qskip")){
				qStep=Integer.parseInt(b);
			}else if(a.equals("rstep") || a.equals("rskip")){
				rStep=Integer.parseInt(b);
			}else if(a.equals("mm")){
				midMaskLen=(Tools.isNumeric(b) ? Integer.parseInt(b) : Parse.parseBoolean(b) ? 1 : 0);
			}else if(a.equals("blacklist") || a.equals("banhomopolymers")){
				Query.blacklistRepeatLength=(Tools.isNumeric(b) ? Integer.parseInt(b) : Parse.parseBoolean(b) ? 1 : 0);
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

	/** Resolves alternate compressed extensions for the query input paths. */
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
	 * Loads queries, starts reference/output services and joins the workers.
	 * Then closes input, waits for alignment/header output and reports accumulated
	 * errors. Shared buffer-data and Read validation settings are restored along this
	 * normal path before the final error-state exception; earlier failures may skip it.
	 * Processed counters include loaded query roots/mates and reference roots/mates.
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
		Shared.setBufferData(Math.min(400000, oldBD));

		final ConcurrentReadInputStream cris=makeCris(refFile);
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

	/** Starts workers sharing the services/buckets, waits for them and accumulates success. */
	private void spawnThreads(final ConcurrentReadInputStream cris, 
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
	private ConcurrentReadInputStream makeCris(String fname){
		FileFormat ff=FileFormat.testInput(fname, null, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
		cris.start();
		if(verbose){outstream.println("Started cris");}
		return cris;
	}

	/**
	 * Materializes query reads, installs shared Query calculators and builds buckets.
	 * Adds each root and mate without a length filter. Indexed queries select their
	 * calculator bucket; queries without a calculator index use bucket zero.
	 * @param ff1 Primary query input
	 * @param ff2 Optional separate mate input
	 * @return Owned bucket lists, prepared before any alignment worker starts
	 */
	public ArrayList<ArrayList<Query>> fetchQueries(FileFormat ff1, FileFormat ff2){
		Timer t=new Timer(outstream, false);
		ArrayList<Read> reads=ConcurrentReadInputStream.getReads(maxReads, false, ff1, ff2, null, null);

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

			addReadToBucket(r, buckets);
			if(r.mate!=null){
				addReadToBucket(r.mate, buckets);
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
	 * Scores supplied zero-based offsets, preserving candidate order.
	 * Overhangs use alignClipped; interior placements delegate to Vector.align.
	 * @param query Query bases for the selected strand
	 * @param ref Reference bases
	 * @param maxSubs Maximum substitution/paid-clipping cost
	 * @param maxClips Total free overhang allowance in bases
	 * @param seedHits Borrowed candidate offsets, possibly null; read without mutation
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
	 * Enumerates offsets for non-null arrays and nonnegative substitution/clip budgets.
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
	 * No overlap returns query.length; otherwise scoring stops after exceeding
	 * maxSubs. This pointwise helper does not enforce the candidate-emission domain.
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
	 * Builds forward masked k-mer keys with lists of zero-based reference starts.
	 * Ambiguous bases invalidate windows even at masked positions. rStep samples
	 * starts at multiples of its power-of-two stride. Short references yield a valid
	 * empty index; disabled indexing or nonpositive k returns null.
	 * @param ref Non-null reference bases, read without retention or mutation
	 * @param k Indexed k-mer length, expected to be 1 through 15
	 * @return Owned index, or null when indexing is disabled for this call
	 */
	IntListHashMap buildReferenceIndex(byte[] ref, int k){
		if(!indexQueries || k<=0){return null;}
		final int defined=Math.max(k-midMaskLen, 2);
		final int kSpace=(1<<(2*defined));
		final long maxKmers=Math.min(kSpace, (ref.length-k+1)*2L);
		//FIXED IFA-006: short references still need a valid empty map for prescan and seed consumers.
		final int initialSize=(int)Math.max(1, Math.min(4000000, ((maxKmers*3)/2)));
		final IntListHashMap index=new IntListHashMap(initialSize);

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

	/**
	 * Serializes all supplied offsets into SAM without another acceptance filter.
	 * Builds match strings and CIGARs, derives NM from SamLine.countSubs and advances
	 * the shared query counter for primary status. Null output still builds/counts
	 * alignments and advances that counter. Input arrays and hits are read only.
	 * @param q Shared query with original and reverse-complement bases
	 * @param ref Original reference root; this version has no fused-reference mapping
	 * @param hits Accepted candidate offsets, possibly null
	 * @param reverseStrand Whether to construct matches from reverse-complement bases
	 * @param bsw Started shared output writer, or null to discard serialized output
	 * @return Number of supplied offsets, or zero for null/empty hits
	 */
	static int processHits(Query q, Read ref, IntList hits, boolean reverseStrand,
		ByteStreamWriter bsw){
		if(hits==null || hits.size()==0){return 0;}
		ByteBuilder bb=new ByteBuilder();
		ByteBuilder match=new ByteBuilder(q.bases.length);

		for(int i=0; i<hits.size(); i++){
			int start=hits.get(i);
			byte[] querySeq=reverseStrand ? q.rbases : q.bases;
			toMatch(querySeq, ref.bases, start, match.clear());

			SamLine sl=new SamLine();
			sl.pos=Math.max(start+1, 1);
			sl.qname=q.name;
			sl.setRnameS(ref.id);
			sl.setSeq(q.bases);
			sl.setQual(q.quals);
			sl.setPrimary(q.alignments.incrementAndGet()==1);
			sl.setMapped(true);
			if(reverseStrand){sl.setStrand(Shared.MINUS);}
			sl.tlen=q.bases.length;
			sl.setCigar(SamLine.toCigar14(match.toBytes(), start, 
				start+q.bases.length-1, ref.length(), q.bases));
			int subs=sl.countSubs();
			sl.addOptionalTag("NM:i:"+subs);
			sl.mapq=Tools.mid(0, (int)(40*(sl.length()*0.5-subs)/(sl.length()*0.5)), 40);
			sl.toBytes(bb).nl();
			if(bb.length()>=16384){
				if(bsw!=null){bsw.addJob(bb);}
				bb=new ByteBuilder();
			}
		}
		if(bsw!=null && !bb.isEmpty()){bsw.addJob(bb);}
		return hits.size;
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
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Shared-reader consumer owning its indexes, seed lists and counters; borrows queries/writers. */
	class ProcessThread extends Thread {

		/** Borrows services/buckets and copies this invocation's scoring settings. */
		ProcessThread(final ConcurrentReadInputStream cris_, final ByteStreamWriter bsw_, final SamHeaderWriter shw_, 
			ArrayList<ArrayList<Query>> qBuckets_, final int maxSubs_, final float minid_, final int tid_){
			cris=cris_;
			bsw=bsw_;
			shw=shw_;
			queryBuckets=qBuckets_;
			maxSubs=maxSubs_;
			minid=minid_;
			tid=tid_;
		}

		/** Marks success only after normal completion of reference processing. */
		@Override
		public void run(){
			synchronized(this){
				processInner();
				success=true;
			}
		}

		/** Processes and returns nonempty batches, then returns the terminal container if present. */
		void processInner(){
			ListNum<Read> ln=cris.nextList();
			while(ln!=null && ln.size()>0){
				processList(ln);
				cris.returnList(ln);
				ln=cris.nextList();
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		/**
		 * Validates each reference root and counts it plus its mate; aligns roots only.
		 * Queues optional header entries under the input batch ID after processing.
		 */
		void processList(ListNum<Read> ln){
			final ArrayList<Read> refList=ln.list;
			final ArrayList<StringNum> alsn=(shw==null ? null : new ArrayList<StringNum>(refList.size()));

			for(int idx=0; idx<refList.size(); idx++){
				final Read ref=refList.get(idx);
				if(!ref.validated()){ref.validate(true);}
				if(alsn!=null) {alsn.add(new StringNum(ref.name(), ref.length()));}

				final int initialLength1=ref.length();
				final int initialLength2=ref.mateLength();
				readsProcessedT+=ref.pairCount();
				basesProcessedT+=initialLength1+initialLength2;

				processRefSequence(ref);
			}
			if(shw!=null) {shw.add(new ListNum<StringNum>(alsn, ln.id));}
		}

		/** Selects list or count-map seed aggregation; the retained rname argument is unused. */
		private IntList getSeedHits(Query q, IntListHashMap refIndex, 
			boolean reverseStrand, IntHashMap2 hitCounts, String rname){
			if(useSeedMap){return getSeedHitsMap(q, refIndex, reverseStrand, hitCounts);}
			else{return getSeedHitsList(q, refIndex, reverseStrand);}
		}

		/**
		 * Applies optional prescan, then sorts decoded offsets and filters by multiplicity.
		 * Returns a newly allocated sorted unique list, or null if no offset qualifies.
		 */
		private IntList getSeedHitsList(Query q, IntListHashMap refIndex, boolean reverseStrand){
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
		 * Applies optional prescan and adds each offset when its count reaches minHits.
		 * Clears a supplied scratch map on the first decoded hit, or allocates one if
		 * null. Returns a new list in threshold-arrival order, or null when none qualify.
		 */
		private IntList getSeedHitsMap(Query q, IntListHashMap refIndex, 
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
		 * Scans valid qStep-spaced query keys until the adjusted miss budget is exceeded.
		 * Returns observed present keys (total-minus-misses), possibly from a partial
		 * scan; unlike version4, no miss-free-pass floor is applied to this count.
		 */
		private int prescan(Query q, IntListHashMap refIndex, boolean reverseStrand, final int minHits){
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

		/** Dispatches one reference root to indexed or brute alignment. */
		long processRefSequence(final Read ref){
			return indexQueries ? processRefSequenceIndexed(ref) : processRefSequenceBrute(ref);
		}

		/** Builds one list-map index per nonempty k bucket and counts both strands' accepted hits. */
		long processRefSequenceIndexed(final Read ref){
			long sum=0;
			final float subrate=1-minid;
			IntHashMap2 seedMap=(useSeedMap ? new IntHashMap2() : null);

			for(int i=0; i<queryBuckets.size(); i++){
				ArrayList<Query> queries=queryBuckets.get(i);
				if(queries.isEmpty()){continue;}

				int currentK=kArray[i];
				IntListHashMap refIndex=buildReferenceIndex(ref.bases, currentK);

				for(Query q : queries){
					int count=0;
					int maxSubsQ=Math.min(maxSubs, (int)(q.length()*subrate));

					IntList seedHits=getSeedHits(q, refIndex, false, seedMap, ref.name());
					IntList hits=alignSparse(q.bases, ref.bases, maxSubsQ, q.maxClips, seedHits);
					count+=processHits(q, ref, hits, false, bsw);

					seedHits=getSeedHits(q, refIndex, true, seedMap, ref.name());
					hits=alignSparse(q.rbases, ref.bases, maxSubsQ, q.maxClips, seedHits);
					count+=processHits(q, ref, hits, true, bsw);

					readsOutT+=count;
					basesOutT+=count*q.bases.length;
					sum+=count;
				}
			}
			return sum;
		}

		/** Enumerates both strands without seed indexes and returns the accepted alignment count. */
		long processRefSequenceBrute(final Read ref){
			long sum=0;
			final float subrate=1-minid;
			for(ArrayList<Query> bucket : queryBuckets){
				for(Query q : bucket){
					int count=0;
					int maxSubsQ=Math.min(maxSubs, (int)(q.length()*subrate));

					IntList hits=alignAllPositions(q.bases, ref.bases, maxSubsQ, q.maxClips);
					count+=processHits(q, ref, hits, false, bsw);

					hits=alignAllPositions(q.rbases, ref.bases, maxSubsQ, q.maxClips);
					count+=processHits(q, ref, hits, true, bsw);

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

		/** Reference roots plus mates; query-loading counts are added by the outer instance. */
		protected long readsProcessedT=0;
		/** Bases in reference roots plus mates, although alignment processes roots only. */
		protected long basesProcessedT=0;
		/** Indexed candidate count, or brute-mode query-count times reference length estimate. */
		protected long alignmentsT=0;
		/** Decoded seed occurrences before offset-threshold aggregation. */
		protected long seedHitsT=0;
		/** Serialized alignment count, including repeated alignments of one query. */
		protected long readsOutT=0;
		/** Sum of full query lengths over accepted alignments, including clipped bases. */
		protected long basesOutT=0;
		/** Set only after normal processing completion; read by the accumulator after joining. */
		boolean success=false;

		/** Shared reference reader; batches are returned by processInner, reader closed by process. */
		private final ConcurrentReadInputStream cris;
		/** Shared alignment writer, or null when output is disabled. */
		private final ByteStreamWriter bsw;
		/** Optional separate shared header writer. */
		private final SamHeaderWriter shw;
		/** Prepared shared query buckets; each query's alignment counter remains mutable. */
		private final ArrayList<ArrayList<Query>> queryBuckets;
		/** Per-worker copies of this instance's scoring options. */
		final int maxSubs;
		final float minid;
		/** Retained worker identifier, currently unused by processing. */
		final int tid;
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
	/** Instance scoring limits; unlike later versions these fields are not static. */
	int maxSubs=5;
	float minid=0;

	/** Descending indexed k values after construction, or the sole value zero in brute mode. */
	int[] kArray=new int[] {10,12,14};

	/** Central masked bases in packed k-mer keys. */
	int midMaskLen=1;
	/** Enables query k-mer preparation and per-reference indexes. */
	boolean indexQueries=true;
	/** Enables the key-presence screen before decoding seed positions. */
	boolean prescan=true;
	/** Query-start and power-of-two reference-start strides; at least one is one. */
	int qStep=1;
	int rStep=1;
	/** Maximum stride captured during construction. */
	final int kStep;

	/** Floor on each query's modeled seed threshold. */
	int minSeedHits=1;
	/** Probability forwarded to MinHitsCalculator2 during serial query setup. */
	private float minHitsProb=0.999f;
	/** Selects count-map aggregation instead of sorting decoded offsets. */
	boolean useSeedMap=false;

	/** Query-loading and reference-processing totals, counting both roots and mates. */
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	/** Accumulated worker counters, including heuristic alignment-work estimates. */
	protected long alignmentCount=0;
	protected long seedHitCount=0;
	protected long readsOut=0;
	protected long basesOut=0;

	/** Limit independently forwarded to query and reference readers; -1 is unlimited. */
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
	/** Process-wide verbosity option. */
	public static boolean verbose=false;
	/** Accumulated worker/input/statistics/output-finalization error state. */
	public boolean errorState=false;
	/** FileFormat output policies, also copied into shared ReadStats settings. */
	private boolean overwrite=true;
	private boolean append=false;
}
