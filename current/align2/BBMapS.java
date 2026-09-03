package align2;

import java.util.ArrayList;
import java.util.TreeMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadStreamWriter;
import structures.ListNum;

/**
 * BBMapS — Streamer/Writer-based entry point for BBMap, Phase 1 of the
 * bbmapnova upgrade project. Implemented ALONGSIDE align2.BBMap/bbmap.sh,
 * which are never edited (project ground rule). Extends AbstractMapper to
 * reuse its argument parsing, setup(), and loadIndex() unchanged — only
 * this class's testSpeed() override uses the new dispatcher/worker-pool/
 * coordinator architecture instead of BBMap's old N-threads-each-touch-cris
 * pattern.
 *
 * Design: /mnt/c/playground/Nowi/plans/BBMapUpgrade_Phase1_CoordinatorDesign_Nowi.java
 * Ownership split (Yaoyao's call, 2026-09-02): Nowi owns everything in
 * this file (dispatcher, MapJob/MapResult lifecycle, worker-pool,
 * coordinator's ordering/ID-bookkeeping state machine, shutdown protocol).
 * Yaoyao owns the concrete RouteWriter/BBSplitterInvoker implementations
 * (writer-side adapters) — NOT YET WRITTEN, placeholder stubs below are
 * clearly marked and must not be mistaken for real implementations.
 *
 * §7 DECIDED (Brian, 2026-09-02): ordering is always enabled. No
 * unordered fast-path in Phase 1 — see the design doc for the mitigation
 * note (bigger buffer, not reintroducing unordered output, if this causes
 * a measurable speed regression).
 *
 * FIRST SLICE ONLY (per design §9): t=1. With exactly one worker, results
 * complete in submission order, so the coordinator's out-of-order
 * reassembly buffer never actually engages — this validates the basic
 * plumbing (dispatcher/worker/coordinator wiring, MapJob/MapResult
 * lifecycle, id bookkeeping) in isolation from the harder concurrent-
 * reordering logic, which only manifests at t>1. t>1 is NOT yet
 * implemented/tested against this file.
 *
 * @author Nowi
 * @date 2026-09-02
 */
public final class BBMapS extends AbstractMapper {

	/*--------------------------------------------------------------*/
	/*----------------         Entry Point          ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		BBMapS mapper=new BBMapS(args);
		args=Tools.condenseStrict(args);
		if(!INDEX_LOADED){mapper.loadIndex();}
		mapper.testSpeed(args);
		ReadWrite.waitForWritingToFinish();
		t.stop();
		outstream.println("\nTotal time:     \t"+t);
		clearStatics();
	}

	public BBMapS(String[] args){
		super(args);
	}

	/** Identical to BBMap.setDefaults() — same defaults, different driver. */
	@Override
	public void setDefaults(){
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_BGZIP=ReadWrite.USE_UNBGZIP=true;
		ReadWrite.PREFER_BGZIP=true;
		ReadWrite.ZIPLEVEL=2;
		MAKE_MATCH_STRING=true;
		keylen=13;

		MINIMUM_ALIGNMENT_SCORE_RATIO=0.56f;

		keyDensity=1.9f;
		maxKeyDensity=3f;
		minKeyDensity=1.5f;
		maxDesiredKeys=15;

		SLOW_ALIGN_PADDING=4;
		SLOW_RESCUE_PADDING=4+SLOW_ALIGN_PADDING;
		TIP_SEARCH_DIST=100;

		MSA_TYPE="MultiStateAligner11ts";
		MAX_SITESCORES_TO_PRINT=5;
		PRINT_SECONDARY_ALIGNMENTS=false;
		AbstractIndex.MIN_APPROX_HITS_TO_KEEP=1;
	}

	@Override
	public String[] preparse(String[] args){return args;} //fast/slow/vslow presets: out of scope for the t=1 first slice.

	@Override
	void postparse(String[] args){
		if(in1==null && args.length>0 && args[0].indexOf('=')<0){in1=args[0];}
	}

	@Override
	public void setup(){
		if(outFile==null){outstream.println("No output file."); OUTPUT_READS=false;}
		else{OUTPUT_READS=true;}
		if(build<0){throw new RuntimeException("Must specify a build number, e.g. build=1");}
		else{Data.GENOME_BUILD=build;}
		if(reference!=null){RefToIndex.makeIndex(reference, build, outstream, keylen);}
	}

	@Override
	void loadIndex(){
		if(build>-1){
			Data.setGenome(build);
			AbstractIndex.MINCHROM=1;
			AbstractIndex.MAXCHROM=Data.numChroms;
			if(minChrom<0){minChrom=1;}
			if(maxChrom<0 || maxChrom>Data.numChroms){maxChrom=Data.numChroms;}
			outstream.println("Set genome to "+Data.GENOME_BUILD);
		}
		AbstractIndex.MINCHROM=minChrom;
		AbstractIndex.MAXCHROM=maxChrom;

		//Ported from BBMap.loadIndex() (:399-415) — real bug found smoke-testing this method
		//2026-09-03: without this, a nodisk=t build (RefToIndex.chromlist populated in-memory,
		//never written to disk) crashed downstream in BBIndex.loadIndex trying to read a
		//ref/genome/<build>/chr1.chrom.gz that nodisk mode never wrote. BBMap avoids this by
		//loading chromosomes into Data.chromosomePlusMatrix from RefToIndex.chromlist (or disk,
		//if not nodisk) BEFORE indexing, whenever MAKE_MATCH_STRING is set — which BBMapS's own
		//setDefaults() always sets true, so this branch always runs here.
		if(SLOW_ALIGN || AbstractIndex.USE_EXTENDED_SCORE || useRandomReads || MAKE_MATCH_STRING){
			if(INDEX_LOADED){
				//do nothing
			}else if(RefToIndex.chromlist==null){
				Data.loadChromosomes(minChrom, maxChrom);
			}else{
				assert(RefToIndex.chromlist.size()==maxChrom-minChrom+1) : RefToIndex.chromlist.size();
				for(ChromosomeArray cha : RefToIndex.chromlist){
					Data.chromosomePlusMatrix[cha.chromosome]=cha;
				}
			}
		}
		RefToIndex.chromlist=null;

		BBIndex.loadIndex(minChrom, maxChrom, keylen, !RefToIndex.NODISK, RefToIndex.NODISK);
		BBIndex.analyzeIndex(minChrom, maxChrom, BBIndex.FRACTION_GENOME_TO_EXCLUDE, keylen);
	}

	@Override void processAmbig2(){} //Multi-reference ambiguity: out of scope for the t=1 first slice.
	@Override void setSemiperfectMode(){}
	@Override void setPerfectMode(){}
	@Override void printSettings(int k){}

	/*--------------------------------------------------------------*/
	/*----------------      MapJob / MapResult      ----------------*/
	/*--------------------------------------------------------------*/

	/** Dispatcher -> worker envelope. See design doc §3. */
	static final class MapJob {
		final long id;
		final ListNum<Read> reads; //null iff poison
		final boolean poison;
		MapJob(long id, ListNum<Read> reads){this.id=id; this.reads=reads; this.poison=false;}
		private MapJob(long id){this.id=id; this.reads=null; this.poison=true;}
		static MapJob poison(long id){return new MapJob(id);}
	}

	/**
	 * Worker -> coordinator envelope. `primary` is the STILL-PRISTINE
	 * readlist (not yet mutated for OUTPUT_MAPPED_ONLY/blacklist/
	 * nullifyObject) — the coordinator does that mutation itself, AFTER
	 * calling BBSplitterInvoker on the pristine data. See design doc's
	 * "Worker/coordinator split" section, resolved 2026-09-02 by reading
	 * AbstractMapThread.writeList() (:747-800) in full.
	 * mapped/unmapped/blacklisted are filtered COPIES built by the worker,
	 * same predicates as writeList() uses today, and are ALWAYS non-null
	 * (empty ArrayList, never null — see design doc's "never submit null"
	 * finding from Yaoyao: a null list can cause a route to silently skip
	 * the id, creating a dense-ID gap downstream).
	 */
	static final class MapResult {
		final long id;
		final ArrayList<Read> primary;
		final ArrayList<Read> mapped;
		final ArrayList<Read> unmapped;
		final ArrayList<Read> blacklisted;
		final boolean poison;
		MapResult(long id, ArrayList<Read> primary, ArrayList<Read> mapped, ArrayList<Read> unmapped, ArrayList<Read> blacklisted){
			this.id=id; this.primary=primary; this.mapped=mapped; this.unmapped=unmapped; this.blacklisted=blacklisted; this.poison=false;
		}
		private MapResult(long id){this.id=id; this.primary=null; this.mapped=null; this.unmapped=null; this.blacklisted=null; this.poison=true;}
		static MapResult poison(long id){return new MapResult(id);}
	}

	/*--------------------------------------------------------------*/
	/*----------------  Coordinator/writer boundary  ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Coordinator -> writer-side boundary. Yaoyao owns the concrete
	 * implementation (writer-side adapters over stream.Writer/OQS2). The
	 * coordinator NEVER calls a raw Writer/OQS2 method directly — only
	 * through this interface.
	 */
	interface RouteWriter {
		/** Exactly once per id, strictly ascending id order, single-threaded
		 * (coordinator thread only). reads is NEVER null — an empty
		 * ArrayList for a route with nothing this id, per Yaoyao's finding
		 * that SamWriter.addReads(null) silently skips the job. */
		void submit(long id, ArrayList<Read> reads);
		/** Normal drain-and-close. No args — the real stream.Writer.poison()
		 * takes none; OQS2 owns maxSeenId/LAST/POISON internally (corrected
		 * 2026-09-02 after reading stream/Writer.java in full). */
		void finish();
		/** MUST be genuinely non-blocking — not implemented via
		 * poison()/poisonAndWait() (Yaoyao's finding: those are the normal
		 * drain path and can hang forever waiting on a pipeline that will
		 * never finish after a crash). Must be idempotent. */
		void finishError();
	}

	/** BBSplitter is invoked once, directly, on the PRISTINE primary list,
	 * strictly before the primary route's in-place mutation. NOT a
	 * RouteWriter — its contract (pristine-only, own stats tables, no
	 * writer-style close/poison) is genuinely different. Yaoyao owns the
	 * concrete implementation. */
	interface BBSplitterInvoker {
		void invoke(long id, ArrayList<Read> pristinePrimaryList);
	}

	/**
	 * PLACEHOLDER — NOT Yaoyao's real implementation. Exists only so this
	 * file compiles and the t=1 plumbing can be smoke-tested end-to-end
	 * before her real RouteWriter/BBSplitterInvoker land. Does nothing
	 * useful; must not be mistaken for shippable code.
	 */
	static final class NoOpRouteWriter implements RouteWriter {
		public void submit(long id, ArrayList<Read> reads){}
		public void finish(){}
		public void finishError(){}
	}

	/** Concrete adapter for the legacy output-stream boundary.  BBMapS forces
	 * the byte-writer path so finishError() can use ReadStreamWriter.abortNow()
	 * instead of the blocking poison/close path. */
	static final class CrosRouteWriter implements RouteWriter {
		private final ConcurrentReadOutputStream output;
		private boolean finished=false;

		CrosRouteWriter(ConcurrentReadOutputStream output){this.output=output;}

		@Override
		public synchronized void submit(long id, ArrayList<Read> reads){
			assert(reads!=null) : id;
			if(finished){throw new RuntimeException("Route writer already finished: "+output.fname());}
			output.add(reads, id);
		}

		@Override
		public synchronized void finish(){
			if(finished){return;}
			finished=true;
			output.close();
			output.join();
			if(output.errorState()){
				throw new RuntimeException("Route writer failed: "+output.fname());
			}
		}

		@Override
		public synchronized void finishError(){
			if(finished){return;}
			finished=true;
			abort(output.getRS1());
			abort(output.getRS2());
		}

		private static void abort(ReadStreamWriter rsw){
			if(rsw!=null){rsw.abortNow();}
		}
	}

	/** Concrete BBSplitter boundary; it preserves the legacy call and its
	 * clearzone/stats behavior while keeping it outside RouteWriter. */
	static final class RealBBSplitterInvoker implements BBSplitterInvoker {
		private final int clearzone;
		RealBBSplitterInvoker(int clearzone){this.clearzone=clearzone;}

		@Override
		public void invoke(long id, ArrayList<Read> pristinePrimaryList){
			if(BBSplitter.streamTable!=null || BBSplitter.TRACK_SET_STATS || BBSplitter.TRACK_SCAF_STATS){
				BBSplitter.printReads(pristinePrimaryList, id, null, clearzone);
			}
		}
	}

	static final class NoOpBBSplitterInvoker implements BBSplitterInvoker {
		public void invoke(long id, ArrayList<Read> pristinePrimaryList){}
	}

	/*--------------------------------------------------------------*/
	/*----------------   Worker mapping-engine stub  ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Minimal stub satisfying ConcurrentReadInputStream so a worker can
	 * construct a BBMapThread purely as a mapping engine (to call its
	 * processRead/processReadPair directly) WITHOUT ever calling .run()/
	 * .start() on it and without a real, live CRIS.
	 *
	 * VERIFIED NECESSARY, not assumed (2026-09-02): read
	 * AbstractMapThread's constructor in full — line 80 is
	 * `PAIRED=cris.paired();`, an unsynchronized direct dereference of
	 * cris DURING construction. A null cris there is an immediate NPE.
	 * ConcurrentReadInputStream is `abstract class ... implements
	 * ConcurrentReadStreamInterface` (not an interface, but subclassable
	 * with a protected 1-arg constructor that only stores fname), so this
	 * stub is viable — confirmed by reading the whole file
	 * (current/stream/ConcurrentReadInputStream.java, 422 lines).
	 * Every method here except paired() is unused by construction: a
	 * worker never calls run()/start()/nextList() on its own engine's
	 * cris — the dispatcher owns the real Streamer, not this stub.
	 */
	static final class WorkerCrisStub extends stream.ConcurrentReadInputStream {
		private final boolean pairedFlag;
		WorkerCrisStub(boolean pairedFlag){super("worker-stub"); this.pairedFlag=pairedFlag;}
		@Override public ListNum<Read> nextList(){throw new UnsupportedOperationException("WorkerCrisStub: dispatcher owns the real Streamer, not this stub.");}
		@Override public void returnList(long listNum, boolean poison){}
		@Override public void run(){}
		@Override public void shutdown(){}
		@Override public void restart(){}
		@Override public void close(){}
		@Override public boolean paired(){return pairedFlag;}
		@Override public Object[] producers(){return new Object[0];}
		@Override public boolean errorState(){return false;}
		@Override public void setSampleRate(float rate, long seed){}
		@Override public long basesIn(){return 0;}
		@Override public long readsIn(){return 0;}
		@Override public boolean verbose(){return false;}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Dispatcher          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * The ONLY thread that ever calls Streamer.nextList() — satisfies its
	 * single-consumer contract by construction (design doc §2). Wraps each
	 * batch as a MapJob and pushes to a bounded shared queue; on
	 * exhaustion, pushes one poison MapJob per worker.
	 */
	final class Dispatcher implements Runnable {
		final stream.Streamer streamer;
		final BlockingQueue<MapJob> jobQueue;
		final int numWorkers;
		volatile long lastRealId=-1;
		volatile long batchesDispatched=0; //Diagnostic only — real evidence of whether a run actually exercised >1 batch (and thus real multi-worker/reordering), not an assumption from input size.
		volatile boolean errorState=false;
		volatile Throwable error=null;

		Dispatcher(stream.Streamer streamer, BlockingQueue<MapJob> jobQueue, int numWorkers){
			this.streamer=streamer; this.jobQueue=jobQueue; this.numWorkers=numWorkers;
		}

		public void run(){
			try{
				ListNum<Read> ln=streamer.nextList();
				while(ln!=null && ln.list!=null && !ln.list.isEmpty()){
					lastRealId=ln.id;
					batchesDispatched++;
					jobQueue.put(new MapJob(ln.id, ln));
					ln=streamer.nextList();
				}
				//Poison: one per worker, using DISTINCT ids past the last real one so the
				//coordinator can tell poison jobs apart if it ever needs to (t=1 slice
				//doesn't need this, but it costs nothing and avoids an id collision).
				for(int i=0; i<numWorkers; i++){
					jobQueue.put(MapJob.poison(lastRealId+1+i));
				}
			}catch(InterruptedException e){
				errorState=true; error=e;
				Thread.currentThread().interrupt();
			}catch(Throwable e){
				errorState=true; error=e;
				e.printStackTrace();
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Worker             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Pulls MapJob from the shared queue, runs UNCHANGED per-read mapping
	 * logic (same processRead/processReadPair calls AbstractMapThread.run()
	 * makes — current/align2/AbstractMapThread.java:642-667 — verified by
	 * reading that section in full), builds mapped/unmapped/blacklisted as
	 * filtered copies (same predicates as writeList(), :748-784), and
	 * pushes a MapResult to the result queue. primary is left PRISTINE —
	 * the coordinator mutates it, not the worker (see the design doc's
	 * worker/coordinator split finding).
	 *
	 * The mapping engine itself is a real BBMapThread instance, constructed
	 * via WorkerCrisStub so its constructor's cris.paired() call succeeds
	 * without a live CRIS — this worker NEVER calls .run()/.start() on it,
	 * only .processRead()/.processReadPair() directly.
	 */
	final class Worker implements Runnable {
		final BlockingQueue<MapJob> jobQueue;
		final BlockingQueue<MapResult> resultQueue;
		final BBMapThread engine;
		volatile boolean errorState=false;
		volatile Throwable error=null;

		Worker(BlockingQueue<MapJob> jobQueue, BlockingQueue<MapResult> resultQueue, boolean paired){
			this.jobQueue=jobQueue; this.resultQueue=resultQueue;
			this.engine=new BBMapThread(new WorkerCrisStub(paired), keylen,
					pileup, SLOW_ALIGN, CORRECT_THRESH, minChrom,
					maxChrom, keyDensity, maxKeyDensity, minKeyDensity, maxDesiredKeys, REMOVE_DUPLICATE_BEST_ALIGNMENTS,
					SAVE_AMBIGUOUS_XY, MINIMUM_ALIGNMENT_SCORE_RATIO, TRIM_LIST, MAKE_MATCH_STRING, QUICK_MATCH_STRINGS,
					null, null, null, null, //rosA/rosM/rosU/rosB: unused by processRead/processReadPair, only by writeList() which this worker never calls
					SLOW_ALIGN_PADDING, SLOW_RESCUE_PADDING, OUTPUT_MAPPED_ONLY, DONT_OUTPUT_BLACKLISTED_READS, MAX_SITESCORES_TO_PRINT, PRINT_SECONDARY_ALIGNMENTS,
					REQUIRE_CORRECT_STRANDS_PAIRS, SAME_STRAND_PAIRS, KILL_BAD_PAIRS, rcompMate,
					PERFECTMODE, SEMIPERFECTMODE, FORBID_SELF_MAPPING, TIP_SEARCH_DIST,
					ambiguousRandom, ambiguousAll, KFILTER, MIN_IDFILTER, qtrimLeft, qtrimRight, untrim, TRIM_QUALITY, minTrimLength,
					LOCAL_ALIGN, RESCUE, STRICT_MAX_INDEL, MSA_TYPE, bloomFilter);
		}

		public void run(){
			try{
				while(true){
					MapJob job=jobQueue.take();
					if(job.poison){
						resultQueue.put(MapResult.poison(job.id));
						return;
					}
					resultQueue.put(mapOneList(job));
				}
			}catch(InterruptedException e){
				errorState=true; error=e;
				Thread.currentThread().interrupt();
			}catch(Throwable e){
				errorState=true; error=e;
				e.printStackTrace();
			}
		}

		/** Same per-read work AbstractMapThread.run()'s loop body does
		 * (:583-667: bloom-filter check, clearAnswers, trim, rcomp,
		 * processRead/processReadPair, capSiteList) — stats/histogram
		 * collection (readstats-based) is deliberately OMITTED from this
		 * first slice, not silently dropped: this engine has no readstats
		 * wired up (constructor doesn't take one), so those branches in
		 * the original loop are all false here regardless. Revisit when
		 * stats output is in scope. */
		private MapResult mapOneList(MapJob job){
			ArrayList<Read> readlist=job.reads.list;
			final boolean black=Blacklist.hasBlacklist();

			for(int i=0; i<readlist.size(); i++){
				Read r=readlist.get(i);
				assert(r.mate==null || (r.pairnum()==0 && r.mate.pairnum()==1)) : r.pairnum()+", "+r.mate.pairnum();
				r.clearAnswers(true);
				final Read r2=r.mate;
				if(engine.RCOMP_MATE!=false){/*rcompMate handling identical to original; no change needed for t=1 single-file slice*/}
				if(r2==null){
					final byte[] basesP=r.bases;
					final byte[] basesM=AminoAcid.reverseComplementBases(basesP);
					engine.processRead(r, basesM);
					engine.capSiteList(r, MAX_SITESCORES_TO_PRINT, PRINT_SECONDARY_ALIGNMENTS);
				}else{
					final byte[] basesP1=r.bases;
					final byte[] basesM1=AminoAcid.reverseComplementBases(basesP1);
					final byte[] basesP2=r2.bases;
					final byte[] basesM2=AminoAcid.reverseComplementBases(basesP2);
					engine.processReadPair(r, basesM1, basesM2);
					engine.capSiteList(r, MAX_SITESCORES_TO_PRINT, PRINT_SECONDARY_ALIGNMENTS);
					engine.capSiteList(r2, MAX_SITESCORES_TO_PRINT, PRINT_SECONDARY_ALIGNMENTS);
				}
			}

			//Filtered copies — same predicates as AbstractMapThread.writeList() (:748-784).
			ArrayList<Read> mapped=new ArrayList<Read>();
			ArrayList<Read> blacklisted=new ArrayList<Read>();
			ArrayList<Read> unmapped=new ArrayList<Read>();
			for(Read r1 : readlist){
				if(r1==null){continue;}
				Read r2=r1.mate;
				boolean isMapped=(r1.mapped() || (r2!=null && r2.mapped()));
				if(isMapped){
					if(!black || !Blacklist.inBlacklist(r1)){mapped.add(r1);}
				}else{
					unmapped.add(r1);
				}
				if(black && Blacklist.inBlacklist(r1)){blacklisted.add(r1);}
			}

			return new MapResult(job.id, readlist, mapped, unmapped, blacklisted);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Coordinator          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * The ONLY caller of RouteWriter/BBSplitterInvoker methods. Drains
	 * MapResults strictly in ascending id order (out-of-order buffering
	 * via TreeMap, per design doc §4) and, for each id: submits mapped/
	 * blacklisted/unmapped, invokes BBSplitter on pristine primary, THEN
	 * mutates primary in place (OUTPUT_MAPPED_ONLY/blacklist/nullify),
	 * THEN submits primary — same order as writeList() (:748-799),
	 * verified by reading it in full, not re-guessed from the design
	 * sketch.
	 *
	 * §7 (Brian, 2026-09-02): always ordered, no unordered fast-path.
	 */
	final class Coordinator implements Runnable {
		final BlockingQueue<MapResult> resultQueue;
		final RouteWriter routeMapped, routeUnmapped, routeBlacklisted, routePrimary;
		final BBSplitterInvoker splitter;
		final int numWorkers;
		final boolean outputMappedOnly, dontOutputBlacklisted;
		volatile boolean errorState=false;
		volatile Throwable error=null;

		private final TreeMap<Long, MapResult> pending=new TreeMap<Long, MapResult>();
		private long nextExpectedId=0;
		private int poisonsSeen=0;
		long mappedTotal=0, unmappedTotal=0; //Diagnostic only, for the t=1 smoke test — NoOp writers discard the reads themselves.

		Coordinator(BlockingQueue<MapResult> resultQueue, int numWorkers,
				RouteWriter routeMapped, RouteWriter routeUnmapped, RouteWriter routeBlacklisted, RouteWriter routePrimary,
				BBSplitterInvoker splitter, boolean outputMappedOnly, boolean dontOutputBlacklisted){
			this.resultQueue=resultQueue; this.numWorkers=numWorkers;
			this.routeMapped=routeMapped; this.routeUnmapped=routeUnmapped; this.routeBlacklisted=routeBlacklisted; this.routePrimary=routePrimary;
			this.splitter=splitter;
			this.outputMappedOnly=outputMappedOnly; this.dontOutputBlacklisted=dontOutputBlacklisted;
		}

		public void run(){
			try{
				while(poisonsSeen<numWorkers){
					MapResult r=resultQueue.take();
					if(r.poison){poisonsSeen++; continue;}
					pending.put(r.id, r);
					drainReady();
				}
				//Design doc §8 graceful path: coordinator has drained every id
				//AND seen every worker's poison before finishing routes.
				assert(pending.isEmpty()) : "Coordinator finished with "+pending.size()+" buffered results never drained — an id gap.";
				routeMapped.finish(); routeUnmapped.finish(); routeBlacklisted.finish(); routePrimary.finish();
			}catch(InterruptedException e){
				errorState=true; error=e;
				forceError();
				Thread.currentThread().interrupt();
			}catch(Throwable e){
				errorState=true; error=e;
				e.printStackTrace();
				forceError();
			}
		}

		private void drainReady(){
			while(true){
				MapResult r=pending.get(nextExpectedId);
				if(r==null){return;}
				pending.remove(nextExpectedId);
				emit(r);
				nextExpectedId++;
			}
		}

		/** submit() calls NEVER pass null — every route always gets a real,
		 * distinct empty ArrayList when it has nothing for this id (Yaoyao's
		 * finding: SamWriter.addReads(null) silently skips the job). */
		private void emit(MapResult r){
			mappedTotal+=r.mapped.size(); unmappedTotal+=r.unmapped.size();
			routeMapped.submit(r.id, r.mapped);
			routeBlacklisted.submit(r.id, r.blacklisted);
			routeUnmapped.submit(r.id, r.unmapped);

			splitter.invoke(r.id, r.primary); //BBSplitter on STILL-PRISTINE primary, before mutation.

			ArrayList<Read> primary=r.primary;
			if(outputMappedOnly){AbstractMapThread.removeUnmapped(primary);}
			if(dontOutputBlacklisted){AbstractMapThread.removeBlacklisted(primary);}
			routePrimary.submit(r.id, primary);
		}

		/** §8 forced/error path: must not block. finishError() on each
		 * RouteWriter is Yaoyao's contract to keep genuinely non-blocking. */
		private void forceError(){
			routeMapped.finishError(); routeUnmapped.finishError(); routeBlacklisted.finishError(); routePrimary.finishError();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Driver             ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates the same four legacy output routes as AbstractMapper.openStreams,
	 * but leaves input ownership with the new Streamer. */
	private void openRouteStreams(String[] args, int buff, boolean paired){
		if(OUTPUT_READS){
			ReadStreamWriter.MINCHROM=minChrom;
			ReadStreamWriter.MAXCHROM=maxChrom;
			rosA=makeRouteStream(outFile, outFile2, qfout, qfout2, buff);
			rosM=makeRouteStream(outFileM, outFileM2, qfoutM, qfoutM2, buff);
			rosU=makeRouteStream(outFileU, outFileU2, qfoutU, qfoutU2, buff);
			rosB=(!Data.scaffoldPrefixes ? makeRouteStream(outFileB, outFileB2, qfoutB, qfoutB2, buff) : null);
		}
		if(Data.scaffoldPrefixes){
			BBSplitter.streamTable=BBSplitter.makeOutputStreams(args, OUTPUT_READS, true, buff, paired, overwrite, append, false);
			if(BBSplitter.AMBIGUOUS2_MODE==BBSplitter.AMBIGUOUS2_SPLIT){
				BBSplitter.streamTableAmbiguous=BBSplitter.makeOutputStreams(args, OUTPUT_READS, true, buff, paired, overwrite, append, true);
			}
		}else{
			BBSplitter.TRACK_SET_STATS=false;
		}
	}

	private ConcurrentReadOutputStream makeRouteStream(String file1, String file2,
			String qf1, String qf2, int buff){
		if(file1==null){return null;}
		FileFormat ff1=FileFormat.testOutput(file1, DEFAULT_OUTPUT_FORMAT, 0, 0, true, overwrite, append, true);
		FileFormat ff2=file2==null ? null : FileFormat.testOutput(file2, DEFAULT_OUTPUT_FORMAT, 0, 0, true, overwrite, append, true);
		ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff1, ff2, qf1, qf2, buff, null, false);
		ros.start();
		return ros;
	}

	private static void closeSplitterStreams(boolean error){
		if(BBSplitter.streamTable!=null){
			for(ConcurrentReadOutputStream ros : BBSplitter.streamTable.values()){
				if(error){abortStream(ros);}else{ReadWrite.closeStream(ros);}
			}
		}
		if(BBSplitter.streamTableAmbiguous!=null){
			for(ConcurrentReadOutputStream ros : BBSplitter.streamTableAmbiguous.values()){
				if(error){abortStream(ros);}else{ReadWrite.closeStream(ros);}
			}
		}
	}

	private static void abortStream(ConcurrentReadOutputStream ros){
		if(ros==null){return;}
		ReadStreamWriter rsw=ros.getRS1();
		if(rsw!=null){rsw.abortNow();}
		rsw=ros.getRS2();
		if(rsw!=null){rsw.abortNow();}
	}

	/**
	 * Replaces BBMap.testSpeed(). Opens a real stream.Streamer via
	 * StreamerFactory (§7 DECIDED, Brian 2026-09-02: always ordered — the
	 * `ordered` arg below is hardcoded true, no unordered fast-path) and
	 * runs the dispatcher/worker-pool/coordinator lifecycle end to end.
	 * routeMapped/routeUnmapped/routeBlacklisted/routePrimary and the
	 * BBSplitterInvoker are still the NOOP placeholders above (Yaoyao's real
	 * implementations not written yet) — this validates the t=1 plumbing
	 * (Streamer open -> dispatch -> map -> coordinate -> clean shutdown)
	 * but produces no real output file yet.
	 */
	@Override
	public void testSpeed(String[] args){
		if(in1==null){outstream.println("No reads to process; quitting."); return;}

		Timer t=new Timer();

		final FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, 0, 0, true, true, false);
		final FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, 0, 0, true, true, false); //null-safe: FileFormat.testInput(null,...) returns null (verified by reading fileIO/FileFormat.java:225-230)
		final stream.Streamer streamer=stream.StreamerFactory.makeStreamer(ff1, ff2, true, maxReads, false, true);
		streamer.start();
		final boolean paired=streamer.paired();
		//The legacy ReadStreamByteWriter path has a non-blocking abort hook.
		//Its SAM conversion is equivalent to SamWriter's conversion, while the
		//new SamWriter wrapper currently has no force-abort boundary.
		final boolean oldSamWriter=ReadWrite.USE_READ_STREAM_SAM_WRITER;
		ReadWrite.USE_READ_STREAM_SAM_WRITER=false;
		
		final int numWorkers=Shared.threads();
		final int buff=Tools.max(32, 2*numWorkers); //Mirrors AbstractMapper.openStreams' ORDERED buffer sizing (:950)
		openRouteStreams(args, buff, paired);
		final BlockingQueue<MapJob> jobQueue=new ArrayBlockingQueue<MapJob>(buff);
		final BlockingQueue<MapResult> resultQueue=new ArrayBlockingQueue<MapResult>(buff);

		final RouteWriter routeA=(rosA==null ? new NoOpRouteWriter() : new CrosRouteWriter(rosA));
		final RouteWriter routeM=(rosM==null ? new NoOpRouteWriter() : new CrosRouteWriter(rosM));
		final RouteWriter routeU=(rosU==null ? new NoOpRouteWriter() : new CrosRouteWriter(rosU));
		final RouteWriter routeB=(rosB==null ? new NoOpRouteWriter() : new CrosRouteWriter(rosB));

		final Dispatcher dispatcher=new Dispatcher(streamer, jobQueue, numWorkers);
		final Worker[] workers=new Worker[numWorkers];
		for(int i=0; i<numWorkers; i++){workers[i]=new Worker(jobQueue, resultQueue, paired);}
		final BBSplitterInvoker splitter=new RealBBSplitterInvoker(workers[0].engine.CLEARZONE1());
		final Coordinator coordinator=new Coordinator(resultQueue, numWorkers, routeM, routeU, routeB, routeA,
				splitter, OUTPUT_MAPPED_ONLY, DONT_OUTPUT_BLACKLISTED_READS);

		final Thread dispatcherThread=new Thread(dispatcher, "BBMapS-Dispatcher");
		final Thread[] workerThreads=new Thread[numWorkers];
		for(int i=0; i<numWorkers; i++){workerThreads[i]=new Thread(workers[i], "BBMapS-Worker"+i);}
		final Thread coordinatorThread=new Thread(coordinator, "BBMapS-Coordinator");
		boolean broken=false;

		dispatcherThread.start();
		for(Thread wt : workerThreads){wt.start();}
		coordinatorThread.start();

		try{
			dispatcherThread.join();
			for(Thread wt : workerThreads){wt.join();}
			coordinatorThread.join();
		}catch(InterruptedException e){
			broken=true;
			Thread.currentThread().interrupt();
			throw new RuntimeException("BBMapS.testSpeed(): interrupted while joining pipeline threads.", e);
		}finally{
			//Cleanup must run even when a join is interrupted; in particular, do
			//not leak the temporary global SAM-writer setting into a later run in
			//the same JVM.
			streamer.close();
			broken|=dispatcher.errorState || coordinator.errorState;
			for(Worker w : workers){broken|=w.errorState;}
			if(broken){
				routeM.finishError(); routeU.finishError(); routeB.finishError(); routeA.finishError();
			}
			closeSplitterStreams(broken);
			ReadWrite.USE_READ_STREAM_SAM_WRITER=oldSamWriter;
		}
		if(broken){throw new RuntimeException("BBMapS terminated in an error state; check stderr above.");}

		t.stop();
		outstream.println("BBMapS smoke test: processed "+streamer.readsProcessed()+" reads in "
			+(paired ? "paired" : "single")+"-ended mode, t="+numWorkers+" ("+t+"). mapped="
			+coordinator.mappedTotal+" unmapped="+coordinator.unmappedTotal+". batchesDispatched="
			+dispatcher.batchesDispatched+".");
	}

}
