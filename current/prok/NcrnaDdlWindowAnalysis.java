package prok;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.Arrays;

import map.LongHashSet;
import map.ObjectSet;
import stream.Read;

/** One input worker's synchronous analysis of ordinary RF00174 windows.
 * The owner brackets its ONE ordinary call per Read with begin/end. This sink
 * verifies each model once on a separate target-only worker, then scores each
 * profile geometry once. It does not run a second genome reader or caller.
 * Seed arrays are cached only for the current Read and strand; native calling
 * reverse-complements the Read's bases in place. The window and current result
 * are borrowed only during Sink.window; consumers must finish using them there.
 * Actual ordinary attempts arrive later and retain the explicit event ordinal.
 * Repeated coordinate tuples and edge-retry attempts are not deduplicated.
 * This component admits no truth and awards no biological recovery credit.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlWindowAnalysis implements NcrnaWindowCaptureSink, NcrnaModelAttemptInstrumentSink{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Installs observation on this input's ordinary target worker. Profiles are
	 * shared read-only; all query, index and verifier scratch belongs to this
	 * instance. The native family controls, including a boundary-refiner factory,
	 * configure the separate worker through the existing initializer helper.
	 */
	public NcrnaDdlWindowAnalysis(final String input, final NcrnaFamily family,
			final NcrnaScavenger ordinary_, final NcrnaDdlProfile[][] profiles,
			final int[] ks, final Sink sink_){
		checkEnvironment();
		require(input!=null && !input.isEmpty() && family!=null && ordinary_!=null && sink_!=null
			&& TARGET.equals(family.name) && TARGET.equals(ordinary_.family)
			&& family.kLong==ConservedRnaSeedIndex.K && ordinary_.cmVerifier==null,
			"Window analysis requires one exact input and the ordinary RF00174 worker without CM verification");
		require(profiles!=null && ks!=null && profiles.length>0 && profiles.length==ks.length,
			"Declare every profile geometry and its matching k before workers start");
		ordinary_.validateOracleModels(family.modelNames);
		ordinary=ordinary_;
		sink=sink_;
		scoring=new NcrnaDdlRanks[profiles.length];
		geometryScoreCpu=new long[profiles.length];
		geometryRankCpu=new long[profiles.length];
		for(int geometry=0; geometry<profiles.length; geometry++){
			require(profiles[geometry]!=null && profiles[geometry].length==family.library.length,
				"Every geometry must preserve the entire original model order");
			scoring[geometry]=new NcrnaDdlRanks(profiles[geometry], ks[geometry]);
		}
		oracle=new NcrnaDdlOracle(input, oracleWorker(family), family.modelNames);
		seedIndex=new ConservedRnaSeedIndex(new LongHashSet[]{family.kmerSet}, new boolean[]{true});
		nativeIndex=new TrnaKmerIndex(family.library, ordinary.indexK, ordinary.adaptiveMinHits,
			ordinary.adaptFloor, ordinary.adaptTopFrac, ordinary.adaptQFrac, ordinary.indexMinHitsDefault);
		nativeScores=new float[family.library.length];
		shortlistRanks=new int[family.library.length];
		filtered=new boolean[family.library.length];
		bean=ManagementFactory.getThreadMXBean();
		owner=Thread.currentThread().getId();
		ordinary.setWindowCaptureSink(this);
		ordinary.setModelAttemptSink(this);
	}

	/** Must run before any output is opened. Nested oracle work would otherwise
	 * be charged to the ordinary caller's active ThreadLocal profiling span. */
	public static void checkEnvironment(){
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Window analysis requires -ea for requested-model identity checks");
		require(!NcrnaCpuProfile.enabled(),
			"Disable native ncRNA CPU profiling: synchronous oracle work needs separate analysis counters");
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported(), "Window analysis needs a supported current-thread CPU clock");
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
		require(clock.getCurrentThreadCpuTime()>=0, "The worker CPU clock must be available before analysis");
	}

	/** Constructs only this target family, using GeneCaller's native constructor
	 * arguments and shared family-control application. No observer is installed.
	 * The ordinary worker's external CM assignment is rejected by the constructor
	 * above rather than silently omitted.
	 */
	public static NcrnaScavenger oracleWorker(final NcrnaFamily f){
		require(f!=null && TARGET.equals(f.name), "The exhaustive worker supports RF00174 only");
		return worker(f, f.idPass, f.idBorderline, true);
	}

	/** Same family controls, selected profile order and unlimited stopping policy
	 * as the ordinary baseline. Only the explicit scalar identity pair changes. */
	public static NcrnaScavenger savedWindowWorker(final NcrnaFamily f,
			final NcrnaDdlProfileResource profiles, final float identity){
		require(f!=null && f.indexTopN==Integer.MAX_VALUE && f.modelPatience==Integer.MAX_VALUE
			&& profiles!=null && profiles.selected() && Float.isFinite(identity) && identity>=0 && identity<=1,
			"Saved acceptance replay requires the resolved unlimited family, selected profiles and an explicit valid identity");
		final NcrnaScavenger worker=worker(f, identity, identity, false);
		worker.useDdlProfiles(profiles, Integer.MAX_VALUE);
		worker.savedWindowFlank();
		return worker;
	}

	/** Same verifier with explicit discovery-seed bypass for shuffled candidates.
	 * Only the requested aligned-span floor changes between the two control arms.
	 */
	public static NcrnaScavenger shuffledWindowWorker(final NcrnaFamily f,
			final NcrnaDdlProfileResource profiles, final float identity, final boolean spanGate){
		require(f!=null && f.minLen==109 && f.indexTopN==Integer.MAX_VALUE && f.modelPatience==Integer.MAX_VALUE
			&& profiles!=null && profiles.selected() && Float.isFinite(identity) && identity>=0 && identity<=1,
			"Shuffled control requires the resolved unlimited109-base RF00174 baseline and explicit valid scalar identity");
		final NcrnaScavenger worker=worker(f, identity, identity, false, spanGate ? f.minLen : 1);
		worker.minKmerHits=0;
		worker.useDdlProfiles(profiles, Integer.MAX_VALUE);
		worker.savedWindowFlank();
		return worker;
	}

	private static NcrnaScavenger worker(final NcrnaFamily f, final float pass,
			final float borderline, final boolean nativeRanking){
		return worker(f, pass, borderline, nativeRanking, f.minLen);
	}

	private static NcrnaScavenger worker(final NcrnaFamily f, final float pass,
			final float borderline, final boolean nativeRanking, final int minimumLength){
		require(f!=null && TARGET.equals(f.name), "The exhaustive worker supports RF00174 only");
		final NcrnaScavenger worker=new NcrnaScavenger(f.library, f.models, f.modelNames,
			f.kmerSet, f.kLong, minimumLength, f.windowPad,
			f.indexK, f.indexTopN, f.adaptive,
			f.adaptFloor, f.adaptTopFrac, f.adaptQFrac, f.fixedMinHits,
			f.scoreA, f.scoreB, pass, borderline,
			f.boundary5NetTemplate, f.boundary3NetTemplate,
			f.boundaryStartTable, f.boundaryStopTable,
			f.boundaryStartInside, f.boundaryStartOutside, f.boundaryStopInside, f.boundaryStopOutside,
			f.boundaryMeanLen, f.boundaryStartOffsets, f.boundaryStopOffsets,
			f.boundaryMarginStart, f.boundaryMarginStop, nativeRanking);
		GeneCaller.applyNcrnaFamilyControls(worker, f);
		if(minimumLength!=f.minLen && worker.physicalEdgeMinLen>0){worker.physicalEdgeMinLen=minimumLength;}
		return worker;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Read / Event Scope     ----------------*/
	/*--------------------------------------------------------------*/

	/** Called immediately before the one ordinary call on an unpaired Read.
	 * Duplicate native GFF seqids are rejected within this input, even when the
	 * full descriptions differ. end() releases this component's Read/seed references;
	 * the native verifier's GC cache may retain its current contig.
	 */
	public void begin(final Read read){
		checkOwner();
		require(!closed && current==null && read!=null && read.bases!=null && read.mate==null,
			"Begin one unchanged unpaired original contig after ending its predecessor");
		final String token=NcrnaDdlOracle.seqid(read.id);
		require(seqids.add(token), "Repeated native GFF seqid would ambiguously join reference parents: "+token);
		current=read;
		seqid=token;
		contigOrdinal++;
		seedStrand=-1;
		seedHits=null;
		seedKeys=null;
	}

	/** Called after ordinary calling, before the Read is replaced. No oracle or
	 * attempt result is inferred from the returned final path; the owner retains
	 * ordinary GFF output for separate production grading.
	 */
	public void end(){
		checkOwner();
		require(current!=null && !capturing, "End exactly one completed ordinary contig call");
		finishEvent();
		current=null;
		seqid=null;
		seedStrand=-1;
		seedHits=null;
		seedKeys=null;
	}

	/** Disarms only the two observers installed by this instance. */
	public void finish(){
		checkOwner();
		require(current==null && !closed, "Finish this input once after ending every contig");
		ordinary.setWindowCaptureSink(null);
		ordinary.setModelAttemptSink(null);
		closed=true;
	}

	@Override
	public void capture(final String name, final int strand_, final int pass_,
			final int start_, final int stop_, final int source_, final byte[] bases, final byte[] window){
		checkOwner();
		require(!closed && !capturing && current!=null && name.equals(current.id) && bases==current.bases
			&& (strand_==0 || strand_==1) && (pass_==1 || pass_==2) && start_>=0
			&& stop_>=start_ && stop_<bases.length && window!=null && window.length==stop_-start_+1,
			"Capture must borrow the current original Read and exact native strand-oriented window");
		finishEvent();
		event=events++;
		strand=strand_;
		pass=pass_;
		start=start_;
		stop=stop_;
		source=source_;
		eventAttempts=0;
		final long seedStart=cpu();
		if(seedStrand!=strand){
			final ConservedRnaSeedIndex.ScanResult scan=seedIndex.scan(bases);
			seedHits=scan.hits(0);
			seedKeys=scan.keys(0);
			seedStrand=strand;
			seedScans++;
		}
		seedCpu+=cpu()-seedStart;
		final long nativeStart=cpu();
		Arrays.fill(shortlistRanks, 0);
		final int[] shortlist=nativeIndex.shortlist(window, ordinary.indexTopN, ordinary.indexScoreMargin);
		final boolean gate=!ordinary.strictIndexCutoff || nativeIndex.lastMaxShared()>=ordinary.indexMinHitsDefault;
		if(gate){
			for(int rank=0; rank<shortlist.length; rank++){shortlistRanks[shortlist[rank]]=rank+1;}
		}
		for(int model=0; model<nativeScores.length; model++){nativeScores[model]=nativeIndex.lastSharedCount(model);}
		nativeCpu+=cpu()-nativeStart;
		final long oracleStart=cpu();
		oracle.verify(name, bases, strand, pass, start, stop, source, seedHits, seedKeys);
		oracleCpu+=cpu()-oracleStart;
		final long rankStart=cpu();
		for(int geometry=0; geometry<scoring.length; geometry++){
			final long before=cpu();
			scoring[geometry].score(window);
			geometryScoreCpu[geometry]+=cpu()-before;
		}
		rankingCpu+=cpu()-rankStart;
		capturing=true;
		try{sink.window(this, window);}
		finally{capturing=false;}
	}

	@Override
	public void modelAttempt(final String name, final int strand_, final int pass_,
			final int start_, final int stop_, final int rank, final int model, final int shared,
			final float identity, final int alignedLength, final boolean identityPass){
		checkOwner();
		require(!closed && current!=null && !capturing && event>=0 && name.equals(current.id)
			&& strand_==strand && pass_==pass && start_==start && stop_==stop
			&& model>=0 && model<shortlistRanks.length && rank>0
			&& shortlistRanks[model]==rank && nativeScores[model]==shared,
			"An ordinary attempt must bind the latest explicit event and its independently counted native shortlist");
		eventAttempts++;
		attempts++;
		sink.attempt(this, rank, model, shared, identity, alignedLength, identityPass);
	}

	private void finishEvent(){
		if(event<0){return;}
		sink.endWindow(this);
		event=-1;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Borrowed Results      ----------------*/
	/*--------------------------------------------------------------*/

	/** Ranks only positive native-index support. Output is the existing four-int
	 * rank/tie/model tuple. Zero-support models remain unreachable, even when an
	 * exhaustive alignment accepts them. Cutoff eligibility is reported separately.
	 */
	public void nativeBest(final boolean[] acceptable, final int[] output){
		checkCapture();
		require(acceptable!=null && acceptable.length==filtered.length, "Native rank needs one oracle bit per original model");
		final long before=cpu();
		for(int i=0; i<filtered.length; i++){filtered[i]=acceptable[i] && nativeScores[i]>0;}
		NcrnaDdlRanks.best(nativeScores, filtered, output);
		final long elapsed=cpu()-before;
		nativeRankCpu+=elapsed;
		rankingCpu+=elapsed;
	}

	/** Scores are ready for every geometry; each query sketch was built once. */
	public void best(final int geometry, final int method, final boolean[] acceptable, final int[] output){
		checkCapture();
		require(geometry>=0 && geometry<scoring.length, "Unknown profile geometry");
		final long before=cpu();
		scoring[geometry].best(method, acceptable, output);
		final long elapsed=cpu()-before;
		geometryRankCpu[geometry]+=elapsed;
		rankingCpu+=elapsed;
	}

	public NcrnaDdlOracle oracle(){checkCapture();return oracle;}
	/** Complete positive-count order is reconstructed by descending count and
	 * ascending original model index. Zero-count models remain unreachable. */
	public int nativeSharedCount(final int model){
		checkCapture();
		require(model>=0 && model<nativeScores.length, "Native count provenance needs one original model index");
		return (int)nativeScores[model];
	}
	public int nativeShortlistRank(final int model){
		require(event>=0 && model>=0 && model<shortlistRanks.length, "Read native eligibility only for this active event/model");
		return shortlistRanks[model];
	}
	public int[] seedPositions(){checkCapture();return seedHits;}
	public long[] seedValues(){checkCapture();return seedKeys;}
	public String seqid(){return seqid;}
	public String header(){return current==null ? null : current.id;}
	public long contigOrdinal(){return contigOrdinal;}
	public long event(){return event;}
	public int strand(){return strand;}
	public int pass(){return pass;}
	public int start(){return start;}
	public int stop(){return stop;}
	public int source(){return source;}
	public int contigLength(){return current==null ? 0 : current.length();}
	public long eventAttempts(){return eventAttempts;}
	public long events(){return events;}
	public long attempts(){return attempts;}
	public long seedScans(){return seedScans;}
	public long seedCpu(){return seedCpu;}
	public long nativeCpu(){return nativeCpu;}
	public long oracleCpu(){return oracleCpu;}
	public long rankingCpu(){return rankingCpu;}
	/** Query sketch/raw-key construction and all three score vectors share this cost. */
	public long geometryScoreCpu(final int geometry){checkOwner();return geometryScoreCpu[geometry];}
	/** Locus-specific first-acceptable rank queries, excluding shared scoring/oracle work. */
	public long geometryRankCpu(final int geometry){checkOwner();return geometryRankCpu[geometry];}
	public long nativeRankCpu(){return nativeRankCpu;}

	private void checkCapture(){
		checkOwner();
		require(capturing, "Borrowed oracle/rank results are valid only inside the synchronous window callback");
	}
	private void checkOwner(){
		require(owner==Thread.currentThread().getId(), "One window-analysis instance belongs to one native input worker");
	}
	private long cpu(){
		final long value=bean.getCurrentThreadCpuTime();
		require(value>=0, "Current-thread CPU accounting became unavailable during analysis");
		return value;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/** The owner consumes results synchronously and emits rows through native
	 * writers. No borrowed window, oracle arrays or seed arrays may escape window().
	 * endWindow follows all actual attempts, including retries and zero attempts.
	 */
	public interface Sink{
		void window(NcrnaDdlWindowAnalysis view, byte[] window);
		void attempt(NcrnaDdlWindowAnalysis view, int rank, int model, int shared,
			float identity, int alignedLength, boolean identityPass);
		void endWindow(NcrnaDdlWindowAnalysis view);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final NcrnaScavenger ordinary;
	private final NcrnaDdlOracle oracle;
	private final NcrnaDdlRanks[] scoring;
	private final long[] geometryScoreCpu, geometryRankCpu;
	private final ConservedRnaSeedIndex seedIndex;
	private final TrnaKmerIndex nativeIndex;
	private final float[] nativeScores;
	private final int[] shortlistRanks;
	private final boolean[] filtered;
	private final ObjectSet<String> seqids=new ObjectSet<String>(String.class);
	private final Sink sink;
	private final ThreadMXBean bean;
	private final long owner;
	private Read current;
	private String seqid;
	private int[] seedHits;
	private long[] seedKeys;
	private int seedStrand=-1, strand, pass, start, stop, source;
	private long event=-1, contigOrdinal=-1, eventAttempts, events, attempts, seedScans;
	private long seedCpu, nativeCpu, oracleCpu, rankingCpu, nativeRankCpu;
	private boolean capturing, closed;
	private static final String TARGET="cobalamin_rf00174";
}
