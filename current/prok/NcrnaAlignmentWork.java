package prok;

import consensus.BaseGraph;
import stream.Read;
import idaligner.AlignmentStats;
import idaligner.ContigEdgeAlignment;
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
	private static final ThreadLocal<NcrnaAlignmentWork> ACTIVE=new ThreadLocal<NcrnaAlignmentWork>();
}
