package prok;

import java.util.Arrays;

import idaligner.AlignmentStats;

/** First-attempt evidence shared by two ranked walks of one exact query/route.
 * A repeated candidate still participates in the second walk's own patience,
 * using its first result rather than executing the model alignment again.
 * begin/end delimit one verifier phase; primary and clipped-edge results never
 * mix, even when their query bytes are equal. No hashing or sequence copies.
 * Primitive score statistics are copied; fresh trace results are retained only
 * until both walks finish. All storage belongs to one caller worker.
 * @author Brian, Raiden
 */
public final class NcrnaPairedWalk{

	public NcrnaPairedWalk(final int models){
		if(models<1 || models>Integer.MAX_VALUE/7){
			throw new IllegalArgumentException("Paired alignment storage needs a positive model count fitting seven statistics per model");
		}
		seen=new boolean[models];
		identities=new float[models];
		positions=new int[models*7];
		traces=new AlignmentStats[models];
	}

	/** The caller owns this immutable query until end; each new phase starts empty. */
	public void begin(final byte[] query, final int route_){
		if(sequence!=null || query==null || query.length<1 || route_<QUANTUM || route_>TRACE){
			throw new IllegalStateException("Begin one nonnested paired phase with its actual positive-length query and verifier route");
		}
		sequence=query;
		route=route_;
		Arrays.fill(seen, false);
	}
	public boolean active(){return sequence!=null;}
	public boolean contains(final int model, final byte[] query, final int requestedRoute){
		validate(model, query, requestedRoute);
		return seen[model];
	}
	public float identity(final int model){
		assert(sequence!=null && model>=0 && model<seen.length && seen[model])
			: "Only a completed first attempt can supply evidence to another ranked position";
		return identities[model];
	}
	public void savePositions(final int model, final byte[] query, final int requestedRoute,
			final float identity, final int[] values){
		validate(model, query, requestedRoute);
		if(route==TRACE || values==null || values.length!=7 || seen[model] || !Float.isFinite(identity)){
			throw new IllegalStateException("A first primary attempt supplies all seven native statistics and one finite identity exactly once");
		}
		System.arraycopy(values, 0, positions, model*7, 7);
		identities[model]=identity;
		seen[model]=true;
	}
	public void restorePositions(final int model, final int[] values){
		assert(route!=TRACE && sequence!=null && model>=0 && model<seen.length && seen[model]
			&& values!=null && values.length==7)
			: "A repeated primary candidate must restore all seven original statistics before score and span checks";
		System.arraycopy(positions, model*7, values, 0, 7);
	}
	public void saveTrace(final int model, final byte[] query, final float identity, final AlignmentStats trace){
		validate(model, query, TRACE);
		if(seen[model] || !Float.isFinite(identity)){
			throw new IllegalStateException("A trace result belongs to exactly one first model attempt in this paired phase");
		}
		identities[model]=identity;
		traces[model]=trace;
		seen[model]=true;
	}
	/** Null is a retained failed trace, distinguished from absence by contains(). */
	public AlignmentStats trace(final int model){
		assert(route==TRACE && sequence!=null && model>=0 && model<seen.length && seen[model])
			: "A repeated trace request may read only its retained first result in this exact phase";
		return traces[model];
	}
	public void end(){
		assert(sequence!=null) : "Every active paired phase must end before a different query or kernel route begins";
		sequence=null;
		Arrays.fill(traces, null);
	}
	private void validate(final int model, final byte[] query, final int requestedRoute){
		if(sequence==null || sequence!=query || route!=requestedRoute || model<0 || model>=seen.length){
			throw new IllegalArgumentException("Paired evidence cannot cross query identity, verifier route or bound model axis");
		}
	}

	private final boolean[] seen;
	private final float[] identities;
	private final int[] positions;
	private final AlignmentStats[] traces;
	private byte[] sequence;
	private int route;
	public static final int QUANTUM=1, SCRABBLE=2, TRACE=3;
}
