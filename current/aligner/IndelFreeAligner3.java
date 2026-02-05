package aligner;

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
 * Version 3: Uses PackedIndex (Hash-CSR) for reference indexing.
 * * Performs high-throughput indel-free alignments using seed-and-extend or brute force strategies.
 * Implements k-mer indexing with rolling hash for query preprocessing and reference indexing.
 * Uses SIMD vectorization (AVX2/SSE) for diagonal alignment when sequences are short enough.
 * Supports multithreaded processing with work-stealing for reference sequence batches.
 * @author Brian Bushnell
 * @contributor Isla, Amber
 * @date June 3, 2025
 */
public class IndelFreeAligner3 implements Accumulator<IndelFreeAligner3.ProcessThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		IndelFreeAligner3 x=new IndelFreeAligner3(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public IndelFreeAligner3(String[] args){

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

	private Parser parse(String[] args){
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
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

	private void doPoundReplacement(){
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}

	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
	}

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

	private static void checkStatics(){
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		assert(FastaReadInputStream.settingsOK());
	}

	private boolean validateParams(){
		for(int k : kArray){
			assert((k>=1 && k<=15) || !indexQueries);
			assert(midMaskLen<k-1 || !indexQueries);
		}
		assert(minHitsProb<=1);
		assert(maxSubs>=0);
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

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

		if(bsw!=null){bsw.poisonAndWait();}
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

	private void spawnThreads(final Streamer cris, 
		final ByteStreamWriter bsw, final SamHeaderWriter shw, final ArrayList<ArrayList<Query>> queryBuckets){

		final int threads=Shared.threads();
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, bsw, shw, queryBuckets, maxSubs, minid, i));
		}
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
	}

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

	@Override
	public final boolean success(){return !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private Streamer makeCris(String fname){
		FileFormat ff=FileFormat.testInput(fname, null, true);
		Streamer cris=StreamerFactory.getReadInputStream(maxReads, true, ff, null, -1);
		cris.start();
		if(verbose){outstream.println("Started cris");}
		return cris;
	}

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

	private void addReadToBucket(Read r, ArrayList<ArrayList<Query>> buckets){
		Query q=new Query(r.id, 0, r.bases, r.quality);
		if(q.calculatorIndex>=0){
			buckets.get(q.calculatorIndex).add(q);
		}else{
			buckets.get(0).add(q); 
		}
	}

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

	public static IntList alignAllPositions(byte[] query, byte[] ref, int maxSubs, int maxClips){
		if(Shared.SIMD && (query.length<256 || maxSubs<256)){
			return SIMDAlignByte.alignDiagonal(query, ref, maxSubs, maxClips);
		}
		IntList list=null;
		int rStart=-maxSubs;
		for(; rStart<0; rStart++){
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
		for(final int limit=ref.length-query.length+maxSubs; rStart<=limit; rStart++){
			int subs=alignClipped(query, ref, maxSubs, maxClips, rStart);
			if(subs<=maxSubs){
				if(list==null){list=new IntList(4);}
				list.add(rStart);
			}
		}
		return list;
	}

	static int align(byte[] query, byte[] ref, final int maxSubs, final int rStart){
		int subs=0;
		for(int i=0, j=rStart; i<query.length && subs<=maxSubs; i++, j++){
			final byte q=query[i], r=ref[j];
			final int incr=(q!=r || AminoAcid.baseToNumber[q]<0 ? 1 : 0);
			subs+=incr;
		}
		return subs;
	}

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
				if(b=='S' || b=='N'){subs++;}
			}

			if((subs+clipPenalty)>maxSubsQ){continue;}

			added++; 
			SamLine sl=new SamLine();
			sl.pos=Math.max(rStart+1, 1);
			sl.qname=q.name;
			sl.setRname(realRef.id);
			sl.seq=q.bases;
			sl.qual=q.quals;
			sl.setPrimary(q.alignments.incrementAndGet()==1);
			sl.setMapped(true);
			if(reverseStrand){sl.setStrand(Shared.MINUS);}
			sl.tlen=qLen;

			sl.cigar=SamLine.toCigar14(match.toBytes(), rStart, rStop-1, rLen, q.bases);

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

	class ProcessThread extends Thread {

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

		@Override
		public void run(){
			synchronized(this){
				if(entropyMask){et=new EntropyTracker(entropyK, entropyWindow, false, entropyCutoff, true);}
				hitsList=new IntList(100); // Reused
				processInner();
				success=true;
			}
		}

		void processInner(){
			for(ListNum<Read> ln=cris.nextList(); ln!=null; ln=cris.nextList()){processList(ln);}
		}

		void processList(ListNum<Read> ln){
			final ArrayList<Read> refList=ln.list;
			final ArrayList<StringNum> alsn=(shw==null ? null : new ArrayList<StringNum>(refList.size()));
			int removed=0;
			for(int idx=0; idx<refList.size(); idx++){
				final Read ref=refList.get(idx);
				if(ref.length()<minRLen){removed++; continue;}
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

		void processListStandard(ArrayList<Read> refList){
			for(int idx=0; idx<refList.size(); idx++){
				final Read ref=refList.get(idx);
				processRefSequence(ref, null, null);
			}
		}

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

		long processRefSequence(final Read ref, ArrayList<Read> originalRefs, IntList ranges){
			if(entropyMask){entropyMask(ref.bases, et);}
			return indexQueries ? processRefSequenceIndexed(ref, originalRefs, ranges) : processRefSequenceBrute(ref, originalRefs, ranges);
		}

		long processRefSequenceIndexed(final Read ref, ArrayList<Read> originalRefs, IntList ranges){
			long sum=0;
			final float subrate=1-minid;
			IntHashMap2 seedMap=(useSeedMap ? new IntHashMap2() : null);

			for(int i=0; i<queryBuckets.size(); i++){
				ArrayList<Query> queries=queryBuckets.get(i);
				if(queries.isEmpty()){continue;}

				int currentK=kArray[i];
				PackedIndex refIndex=new PackedIndex(ref.bases, currentK, midMaskLen, rStep);

				for(Query q : queries){
					int count=0;
					int maxSubsQ=Math.min(maxSubs, (int)(q.length()*subrate));

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
		
		private int prescan(Query q, PackedIndex refIndex, boolean reverseStrand, final int minHits){
			final int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			// Adjust maxMisses: if we require more hits than the query's internal baseline, we have less wiggle room.
			final int maxMisses=q.maxMisses-(minHits-q.minHits);
			
			if(queryKmers==null || maxMisses<0){return 0;}

			int misses=0, total=0;
			// Direct access to the map is safe and fastest
			final map.IntLongHashMap2 map=refIndex.map;

			final int step=2*qStep;
			for(int start=0; start<=qStep && misses<=maxMisses; start+=qStep) {
				for(int i=start; i<queryKmers.length && misses<=maxMisses; i+=step){
					int kmer=queryKmers[i];
					if(kmer==-1){continue;}

					total++;//This could be built into query...
					misses+=(map.containsKey(kmer) ? 0 : 1);
					// Fail Fast: If we have already missed too many, success is impossible.
				}
			}
			return total-misses;
		}

		private IntList getSeedHits(Query q, PackedIndex refIndex, boolean reverseStrand, IntHashMap2 hitCounts){
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

		private IntList getSeedHitsMap(int[] queryKmers, PackedIndex refIndex, 
				int minHits, IntHashMap2 hitCounts){
			//			IntList seedHits=new IntList();//hitsList; 
			IntList seedHits=hitsList; 
			seedHits.clear();

			if(hitCounts==null){hitCounts=new IntHashMap2();}
			else{hitCounts.clear();}

			for(int i=0; i<queryKmers.length; i+=qStep){
				if(queryKmers[i]==-1){continue;}
				long packed=refIndex.get(queryKmers[i]);
				if(packed==-1){continue;}

				int offset=(int)(packed>>>32);
				int count=(int)packed;
				seedHitsT+=count;
				int[] positions=refIndex.positions;
				
				if(count==1) {
					int alignStart=offset-i;
					int newCount=hitCounts.increment(alignStart);
					if(newCount==minHits){seedHits.add(alignStart);}
				}else {
					for(int j=0; j<count; j++){
						int refPos=positions[offset+j];
						int alignStart=refPos-i;
						int newCount=hitCounts.increment(alignStart);
						if(newCount==minHits){seedHits.add(alignStart);}
					}
				}
			}
			return seedHits;
		}

//		private IntList getSeedHitsList(int[] queryKmers, PackedIndex refIndex, int minHits) {
//			//			IntList seedHits=new IntList();//hitsList; 
//			final IntList seedHits=hitsList; 
//			seedHits.clear();
//
//			final int[] positions=refIndex.positions;
//			for(int i=0; i<queryKmers.length; i+=qStep){
//				if(queryKmers[i]==-1){continue;}
//				final long packed=refIndex.get(queryKmers[i]);
//				if(packed==-1){continue;}
//
//				int offset=(int)(packed>>>32);
//				int count=(int)packed;
//				seedHits.ensureCapacity(count);
//				for(int j=0; j<count; j++){
//					int refPos=positions[offset+j];
//					int alignStart=refPos-i;
//					seedHits.addUnchecked(alignStart);
//				}
//			}
//			seedHitsT+=seedHits.size;
//			if(seedHits.size>1 || (minHits>1 && seedHits.size>0)){
//				seedHits.sort();
//				seedHits.condenseMinCopies(minHits);
//			}
//			return seedHits;
//		}
		
		private IntList getSeedHitsList(int[] queryKmers, PackedIndex refIndex, int minHits) {
			final IntList seedHits=hitsList; 
			seedHits.clear();

			final int[] positions=refIndex.positions;
			
			for(int i=0; i<queryKmers.length; i+=qStep){
				if(queryKmers[i]==-1){continue;}
				
				final long packed=refIndex.get(queryKmers[i]);
				if(packed==-1){continue;} // Missing key

				// Unpack count (low 32 bits)
				final int count=(int)packed;
				
				if(count==1){
					// Singleton Optimization: RefPos is in high bits
					int refPos=(int)(packed>>>32);
					int alignStart=refPos-i;
					seedHits.add(alignStart);
				}else{
					// Multi-Hit: Offset is in high bits
					int offset=(int)(packed>>>32);
					seedHits.ensureCapacity(count);
					for(int j=0; j<count; j++){
						int refPos=positions[offset+j];
						int alignStart=refPos-i;
						seedHits.addUnchecked(alignStart);
					}
				}
			}
			seedHitsT+=seedHits.size;
			if(seedHits.size>1 || (minHits>1 && seedHits.size>0)){
				seedHits.sort();
				seedHits.condenseMinCopies(minHits);
			}
			return seedHits;
		}

		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		protected long alignmentsT=0;
		protected long seedHitsT=0;
		protected long readsOutT=0;
		protected long basesOutT=0;
		protected long prescansT=0;
		protected long prescanFailsT=0;
		protected long postscansT=0;
		protected long postscan2T=0;
		boolean success=false;

		private final Streamer cris;
		private final ByteStreamWriter bsw;
		private final SamHeaderWriter shw;
		private final ArrayList<ArrayList<Query>> queryBuckets;
		final int maxSubs;
		final float minid;
		final int tid;
		EntropyTracker et;
		IntList hitsList;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;
	private String in2=null;
	private String out1=null;
	private String headerOut=null;
	private String extin=null;
	private String extout=null;
	String refFile=null;
	static int maxSubs=5;
	static float minid=0.85f;

	int[] kArray=new int[] {8,9,10,12,14};

	int midMaskLen=1;
	boolean indexQueries=true;
	boolean prescan=true;
	int qStep=1;
	int rStep=1;
	final int kStep;

	int minSeedHits=1;
	private float minHitsProb=0.999f;
	boolean useSeedMap=false;

	boolean fuse=true;
	int padding=128;
	int targetChunkSize=1000000;
	int minQLen=1;
	int minRLen=1;

	protected long readsProcessed=0;
	protected long basesProcessed=0;
	protected long alignmentCount=0;
	protected long seedHitCount=0;
	protected long readsOut=0;
	protected long basesOut=0;
	protected long prescans=0;
	protected long prescanFails=0;
	protected long postscans=0;
	protected long postscan2=0;

	private long maxReads=-1;

	private final FileFormat ffin1;
	private final FileFormat ffin2;
	private final FileFormat ffout1;
	private final FileFormat ffheader;

	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
	private boolean entropyMask=false;
	/** Window size for sliding window entropy/complexity analysis */
	private int entropyWindow=80;
	/** Entropy threshold below which regions are masked */
	private float entropyCutoff=0.70f;
	private int entropyK=4;
	
}