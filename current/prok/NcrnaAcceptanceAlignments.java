package prok;

import java.util.Arrays;

/** Worker-owned ordinary PRIMARY measurements for one borrowed genomic query.
 * Owns no aligner, ranking engine, family resources or CM policy. Measured
 * coordinates are local to the genomic query, regardless of kernel operands.
 * All storage is allocated once; recording attempts allocates no objects.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaAcceptanceAlignments{

	public NcrnaAcceptanceAlignments(final int modelCount){
		require(modelCount>0, "Acceptance reuse requires a nonempty fixed global model axis");
		states=new byte[modelCount];
		identities=new float[modelCount];
		starts=new int[modelCount];
		stops=new int[modelCount];
		rankedPresent=new boolean[modelCount];
		clear();
	}

	/** Begins one ordinary PRIMARY window. The query is borrowed and must remain
	 * unchanged until clear. Ranked IDs and all seven hint integers are copied;
	 * unused capacity of the caller's ranked buffer is never inspected.
	 */
	public void begin(final byte[] query, final int senseOrigin, final boolean quantum,
			final int[] rankedModels, final int liveRankedCount, final int[] initialHint){
		require(!active(), "Clear the preceding PRIMARY window before borrowing another query");
		require(query!=null && query.length>0 && senseOrigin>=0
			&& (long)senseOrigin+query.length-1<=Integer.MAX_VALUE,
			"The borrowed genomic query must occupy a valid sense-coordinate interval");
		require(rankedModels!=null && liveRankedCount>=0 && liveRankedCount<=rankedModels.length
			&& liveRankedCount<=states.length, "Use the live unique-model rank count, never backing-array capacity");
		requireHint(initialHint);
		final int count=Math.min(MAX_RANKED, liveRankedCount);
		Arrays.fill(rankedPresent, false);
		for(int rank=0; rank<count; rank++){
			final int model=rankedModels[rank];
			require(model>=0 && model<states.length && !rankedPresent[model],
				"Copied continuation ranks must contain unique in-range global model IDs");
			rankedPresent[model]=true;
		}
		System.arraycopy(rankedModels, 0, ranked, 0, count);
		System.arraycopy(initialHint, 0, hint, 0, HINT_LENGTH);
		rankCount=count;
		origin=senseOrigin;
		quantumRoute=quantum;
		borrowed=query;
	}

	public boolean active(){return borrowed!=null;}

	/** Detects accidental use by another query, coordinate frame or kernel route.
	 * Matching bytes in a different array are not the same borrowed window.
	 */
	public void requireWindow(final byte[] query, final int senseOrigin, final boolean quantum){
		requireActive();
		require(query==borrowed && senseOrigin==origin && quantum==quantumRoute,
			"Acceptance reuse must bind the exact active PRIMARY query, sense origin and ordinary route");
	}

	public byte[] query(){
		requireActive();
		return borrowed;
	}
	public int senseOrigin(){
		requireActive();
		return origin;
	}
	public boolean quantum(){
		requireActive();
		return quantumRoute;
	}
	public int rankedCount(){
		requireActive();
		return rankCount;
	}
	public int rankedModel(final int rank0){
		requireActive();
		require(rank0>=0 && rank0<rankCount, "Continuation cannot read spare rank capacity or ranks beyond the copied top100");
		return ranked[rank0];
	}
	public int state(final int model){
		requireModel(model);
		return states[model];
	}

	/** A prealignment length rejection consumes no kernel work or hint update. */
	public void recordLengthRejected(final int model){
		requireUnvisited(model);
		states[model]=LENGTH_REJECTED;
	}

	/** Snapshots an actual alignment, including failed-span attempts and models
	 * beyond the copied ranked list. The independent final outgoing hint always
	 * follows the last real attempt; cached per-model reads never roll it back.
	 * The caller converts Scrabble model-reference coordinates into query-local
	 * genomic coordinates before calling this method. A failed/empty alignment
	 * may retain (-1, -1) only when spanValid is false; no endpoint is invented.
	 */
	public void recordAligned(final int model, final float identity, final int localStart,
			final int localStop, final boolean spanValid, final int[] outgoingHint){
		requireUnvisited(model);
		require(Float.isFinite(identity) && identity>=0 && identity<=1,
			"Only a finite measured ordinary identity can be reused");
		final boolean absent=localStart==-1 && localStop==-1;
		require(!spanValid && absent || localStart>=0 && localStop>=localStart && localStop<borrowed.length,
			"Endpoints must be query-local; only an invalid-span result may explicitly lack both coordinates (-1, -1)");
		requireHint(outgoingHint);
		identities[model]=identity;
		starts[model]=localStart;
		stops[model]=localStop;
		states[model]=(byte)(spanValid ? ALIGNED_VALID : ALIGNED_INVALID_SPAN);
		System.arraycopy(outgoingHint, 0, hint, 0, HINT_LENGTH);
	}

	/** Invalid-span measurements retain their real identity/endpoints for reuse
	 * as failures. Their endpoints may both be -1 when no genomic span exists;
	 * consumers must require ALIGNED_VALID before making a CM anchor.
	 */
	public float identity(final int model){
		requireMeasured(model);
		return identities[model];
	}
	public int localStart(final int model){
		requireMeasured(model);
		return starts[model];
	}
	public int localStop(final int model){
		requireMeasured(model);
		return stops[model];
	}

	/** Copies the outgoing PRIMARY sequential state, not a per-model result.
	 * This seven-int vector is separate from normalized genomic-query endpoints.
	 */
	public void restoreHint(final int[] destination){
		requireActive();
		requireHint(destination);
		System.arraycopy(hint, 0, destination, 0, HINT_LENGTH);
	}

	/** Idempotent cleanup for finally blocks, including failed/empty windows.
	 * Drop the borrowed query and invalidate per-model states and ranked access.
	 * Stale payloads are inaccessible through the checked API: a measured record
	 * overwrites identity/endpoints, and begin overwrites the full hint and live
	 * rank prefix and resets rank-validation scratch before exposing a new window.
	 */
	public void clear(){
		borrowed=null;
		origin=-1;
		quantumRoute=false;
		rankCount=0;
		Arrays.fill(states, (byte)UNVISITED);
	}

	private void requireActive(){
		require(active(), "PRIMARY window has ended or has not begun; cached state is unavailable");
	}
	private void requireModel(final int model){
		requireActive();
		require(model>=0 && model<states.length, "Global model ID must belong to this worker's fixed axis");
	}
	private void requireUnvisited(final int model){
		requireModel(model);
		require(states[model]==UNVISITED, "A PRIMARY model attempt can be recorded only once per borrowed window");
	}
	private void requireMeasured(final int model){
		requireModel(model);
		require(states[model]==ALIGNED_INVALID_SPAN || states[model]==ALIGNED_VALID,
			"Unvisited and length-rejected models have no measured identity or endpoints");
	}
	private static void requireHint(final int[] vector){
		require(vector!=null && vector.length==HINT_LENGTH, "Ordinary sequential state must retain exactly seven integers");
	}
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	public static final int UNVISITED=0, LENGTH_REJECTED=1, ALIGNED_INVALID_SPAN=2, ALIGNED_VALID=3;
	public static final int HINT_LENGTH=7, MAX_RANKED=100;
	private final byte[] states;
	private final float[] identities;
	private final int[] starts, stops;
	private final boolean[] rankedPresent;
	private final int[] ranked=new int[MAX_RANKED], hint=new int[HINT_LENGTH];
	private byte[] borrowed;
	private int origin=-1, rankCount;
	private boolean quantumRoute;
}
