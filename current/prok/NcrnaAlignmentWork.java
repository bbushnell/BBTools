package prok;

import consensus.BaseGraph;
import stream.Read;
import structures.LongList;
import idaligner.AlignmentStats;
import idaligner.ContigEdgeAlignment;
import idaligner.EndClippedAligner;
import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;

/** Worker-owned alignment work, attributed by the enclosing native family call.
 * Counts logical requests independently of physical cache fills. Model bases are
 * the actual model operand length, including an observed slice used by edge
 * refinement. CPU nanoseconds cover the alignment API call only,
 * not window copies, ranking, HBM scoring or cache lookup. No per-request objects.
 * The enclosing profiler installs one instance while executing a family and
 * clears it on exit. Independent workers never share a mutable counter.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaAlignmentWork{

	/** Selected-format primary operand order and fresh position hints, shared by
	 * production and saved-window measurements. This does not apply acceptance gates. */
	public static float selectedQuantum(final byte[] model, final byte[] window, final int[] positions){
		assert(model!=null && window!=null && positions!=null && positions.length>=4)
			: "Selected primary verification needs model query, genomic reference and four endpoint hints";
		java.util.Arrays.fill(positions, 0);
		return quantum(model, window, positions, model.length);
	}

	/** Real physical-edge kernel; explicit flags must come from original contig
	 * bounds, never from the ends of a detached saved window. */
	public static EndClippedAligner.Result contigEdges(final EndClippedAligner aligner,
			final byte[] model, final byte[] window, final boolean left, final boolean right,
			final int match, final int penalty){
		assert(aligner!=null && model!=null && window!=null && (left || right))
			: "Only an original physical contig edge permits model-end clipping";
		final NcrnaAlignmentWork work=current();
		final long start=work==null ? -1 : NcrnaCpuProfile.cpuNow();
		final EndClippedAligner.Result result=aligner.alignContigEdges(model, window, left, right, match, penalty);
		if(work!=null){work.executed(model.length, window.length, start);}
		return result;
	}

	/** Ordinary uncached Scrabble identity, retaining the caller's operand order. */
	public static float scrabble(final byte[] query, final byte[] reference,
			final int[] positions, final int modelLength){
		final NcrnaAlignmentWork work=current();
		final long start=work==null ? -1 : NcrnaCpuProfile.cpuNow();
		final float result=ScrabbleAligner.alignStatic(query, reference, positions);
		if(work!=null){work.executed(modelLength, candidateLength(query, reference, modelLength), start);}
		return result;
	}

	/** Ordinary uncached Quantum identity, including its incoming position hint. */
	public static float quantum(final byte[] query, final byte[] reference,
			final int[] positions, final int modelLength){
		final NcrnaAlignmentWork work=current();
		final long start=work==null ? -1 : NcrnaCpuProfile.cpuNow();
		final float result=QuantumAligner.alignStatic(query, reference, positions);
		if(work!=null){work.executed(modelLength, candidateLength(query, reference, modelLength), start);}
		return result;
	}

	/** Traceback is part of this actual alignment request, not another attempt. */
	public static float scrabbleTrace(final byte[] query, final byte[] reference,
			final AlignmentStats stats, final int modelLength){
		final NcrnaAlignmentWork work=current();
		final long start=work==null ? -1 : NcrnaCpuProfile.cpuNow();
		final float result=ScrabbleAligner.alignAndTraceStatic(query, reference, stats);
		if(work!=null){work.executed(modelLength, candidateLength(query, reference, modelLength), start);}
		return result;
	}

	/** Preserves the selected Quantum bandwidth multiplier and mutable stats. */
	public static float quantumTrace(final byte[] query, final byte[] reference,
			final AlignmentStats stats, final int modelLength, final int capMultiplier){
		final NcrnaAlignmentWork work=current();
		final long start=work==null ? -1 : NcrnaCpuProfile.cpuNow();
		final float result=capMultiplier==1 ? QuantumAligner.alignAndTraceStatic(query, reference, stats)
			: QuantumAligner.alignAndTraceStatic(query, reference, stats, capMultiplier);
		if(work!=null){work.executed(modelLength, candidateLength(query, reference, modelLength), start);}
		return result;
	}

	/** HBM likelihood evaluation is timed separately from its traceback alignment. */
	public static float hbm(final BaseGraph model, final Read read){
		final NcrnaAlignmentWork work=current();
		final long start=work==null ? -1 : NcrnaCpuProfile.cpuNow();
		final float result=model.score(read, false, true);
		if(work!=null){
			work.hbmCalls++;
			work.hbmCpuNanos+=NcrnaCpuProfile.cpuNow()-start;
		}
		return result;
	}

	/** Records an exact cache hit: logical work remains, physical kernel work does not. */
	public void reused(final int modelLength, final int queryLength){
		checkLengths(modelLength, queryLength);
		logical++;
		modelBases+=modelLength;
		queryBases+=queryLength;
		hits++;
		if(queryObserver!=null){queryObserver.query(modelLength, queryLength, false);}
	}

	/** Records a specialized kernel whose caller owns its timing boundary. */
	public void executed(final int modelLength, final int queryLength, final long start){
		checkLengths(modelLength, queryLength);
		final long elapsed=NcrnaCpuProfile.cpuNow()-start;
		assert(start>=0 && elapsed>=0) : shared.KillSwitch.assertDie("An installed family profile requires a supported monotonic worker CPU clock");
		completed(modelLength, queryLength, elapsed);
	}

	/** Specialized adapters may supply an already measured kernel CPU duration. */
	public void completed(final int modelLength, final int queryLength, final long elapsed){
		checkLengths(modelLength, queryLength);
		assert(elapsed>=0) : shared.KillSwitch.assertDie("Specialized alignment CPU must come from an enabled worker clock");
		logical++;
		physical++;
		modelBases+=modelLength;
		queryBases+=queryLength;
		physicalModelBases+=modelLength;
		cpuNanos+=elapsed;
		if(queryObserver!=null){queryObserver.query(modelLength, queryLength, true);}
	}

	/** One worker-owned observer borrows actual kernel operand lengths. It is
	 * removed at window exit; no request object or additional alignment is made. */
	public interface QueryObserver{void query(int modelLength, int queryLength, boolean physical);}
	void setQueryObserver(final QueryObserver observer){
		assert(observer!=null && (queryObserver==null || queryObserver==observer))
			: shared.KillSwitch.assertDie("A profiled window cannot replace another consumer's kernel observer");
		queryObserver=observer;
	}
	void clearQueryObserver(final QueryObserver observer){
		assert(queryObserver==null || queryObserver==observer)
			: shared.KillSwitch.assertDie("A window may clear only its own borrowed kernel observer");
		queryObserver=null;
	}

	/** One scheduled generic window, spanning primary, rescue and edge retry.
	 * Snapshots precede all work, so early rejects contribute a real zero bin.
	 * Nested alignments are requests inside this unit, never additional windows. */
	public void window(final long logicalBefore, final long physicalBefore){
		final long requests=logical-logicalBefore, executions=physical-physicalBefore;
		assert(logicalBefore>=0 && physicalBefore>=0 && requests>=executions && executions>=0
			&& requests<Integer.MAX_VALUE)
			: shared.KillSwitch.assertDie("One nonnested window must conserve monotonic logical/physical requests and fit an integer histogram bin");
		windows++;
		windowLogical+=requests;
		windowPhysical+=executions;
		logicalHistogram.increment((int)requests);
		physicalHistogram.increment((int)executions);
	}

	public int windowBins(){return Math.max(logicalHistogram.size, physicalHistogram.size);}

	/** One experimental verifier acceptance before candidate resolution and DP.
	 * Repeated windows/retries remain events, never unique loci or retained calls.
	 * Native support is unavailable when the sketch is the only ranking index. */
	public void proposal(final int nativeCount, final boolean fromSketch){
		assert(nativeCount>=-1 && nativeCount<=65535 && (nativeCount>=0 || fromSketch))
			: shared.KillSwitch.assertDie("An accepted experimental proposal must retain its winning stream and measured native count, or sketch-only NA");
		final int stream=nativeCount<0 ? 2 : fromSketch ? 1 : 0;
		final int support=nativeCount<0 ? 2 : nativeCount<2 ? 0 : 1;
		proposalCounts[stream*3+support]++;
	}

	/** Fixed stream/support cells avoid per-acceptance objects and boxed keys. */
	public long proposalCount(final int stream, final int support){
		assert(stream>=0 && stream<3 && support>=0 && support<3)
			: "Proposal report indices are native/sketch/sketch-only by below-floor/at-floor/NA";
		return proposalCounts[stream*3+support];
	}

	/** Attribution-only alphabet census, outside the shortlist CPU timer.
	 * Counts query invocations, including a repeated physical-edge query, not loci.
	 * N is supported; other undefined symbols expose V1's unnormalized-input limit. */
	public void rankingQuery(final byte[] sequence){
		assert(sequence!=null) : "Ranking alphabet attribution requires the actual query bytes before encoding";
		rankingQueries++;
		for(final byte base:sequence){
			final int upper=base&0xdf;
			if(upper!='A' && upper!='C' && upper!='G' && upper!='T' && upper!='N'){
				undefinedRankingQueries++;
				return;
			}
		}
	}
	public long windowCount(final boolean requests, final int bin){
		assert(bin>=0) : "A window-work histogram is indexed by nonnegative alignment-request counts";
		return (requests ? logicalHistogram : physicalHistogram).get(bin);
	}

	/** Nearest-rank count percentile over actual windows, including zero work. */
	public int windowPercentile(final boolean requests, final int percent){
		assert(windows>0 && percent>0 && percent<=100)
			: "Window percentiles require at least one observed window and an integer percentage1..100";
		final long target=(windows/100)*percent+((windows%100)*percent+99)/100;
		final LongList histogram=requests ? logicalHistogram : physicalHistogram;
		long cumulative=0;
		for(int bin=0; bin<histogram.size; bin++){
			cumulative+=histogram.get(bin);
			if(cumulative>=target){return bin;}
		}
		throw new IllegalStateException("Window histogram lost observations before percentile calculation");
	}

	/** Merge only after joining the owning worker; shared writers never mutate it. */
	public void add(final NcrnaAlignmentWork other){
		assert(other!=null && other!=this) : "Joined family counters merge once into a separate aggregate";
		logical+=other.logical;
		physical+=other.physical;
		hits+=other.hits;
		modelBases+=other.modelBases;
		queryBases+=other.queryBases;
		physicalModelBases+=other.physicalModelBases;
		cpuNanos+=other.cpuNanos;
		hbmCalls+=other.hbmCalls;
		hbmCpuNanos+=other.hbmCpuNanos;
		windows+=other.windows;
		windowLogical+=other.windowLogical;
		windowPhysical+=other.windowPhysical;
		rankingQueries+=other.rankingQueries;
		undefinedRankingQueries+=other.undefinedRankingQueries;
		for(int i=0; i<proposalCounts.length; i++){
			proposalCounts[i]=Math.addExact(proposalCounts[i], other.proposalCounts[i]);
		}
		logicalHistogram.incrementBy(other.logicalHistogram);
		physicalHistogram.incrementBy(other.physicalHistogram);
	}

	private static int candidateLength(final byte[] query, final byte[] reference, final int modelLength){
		assert(modelLength==query.length || modelLength==reference.length)
			: shared.KillSwitch.assertDie("The caller must identify the actual model operand, not infer its role from shorter/longer ordering");
		return modelLength==query.length ? reference.length : query.length;
	}

	private static void checkLengths(final int modelLength, final int queryLength){
		assert(modelLength>0 && queryLength>0) : shared.KillSwitch.assertDie("A measured kernel requires actual nonempty model and candidate operands");
	}

	/** One shared stateless observer; the active worker supplies counter ownership. */
	public static ContigEdgeAlignment.TraceAlignment refinementObserver(){
		return current()==null ? null : REFINEMENT;
	}

	private static final ContigEdgeAlignment.TraceAlignment REFINEMENT=(model, candidate, stats) ->
		scrabbleTrace(model, candidate, stats, model.length);

	/** Null outside an explicitly enabled native family profiling scope. */
	public static NcrnaAlignmentWork current(){return NcrnaCpuProfile.enabled() ? ACTIVE.get() : null;}

	static void enter(final NcrnaAlignmentWork work){
		assert(work!=null && ACTIVE.get()==null) : shared.KillSwitch.assertDie("Family profiling scopes cannot nest or overlap on one worker");
		ACTIVE.set(work);
	}

	static void leave(){ACTIVE.set(null);}

	public long logical, physical, hits, modelBases, queryBases, physicalModelBases, cpuNanos, hbmCalls, hbmCpuNanos;
	public long windows, windowLogical, windowPhysical;
	public long rankingQueries, undefinedRankingQueries;
	private final long[] proposalCounts=new long[9];
	private QueryObserver queryObserver;
	private final LongList logicalHistogram=new LongList(16), physicalHistogram=new LongList(16);
	private static final ThreadLocal<NcrnaAlignmentWork> ACTIVE=new ThreadLocal<NcrnaAlignmentWork>();
}
