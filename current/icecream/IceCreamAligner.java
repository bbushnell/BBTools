package icecream;

import aligner.AlignmentResult;
import shared.Shared;

/**
 * Forward overlap-alignment interface used by IceCreamFinder workers.
 * Implementations retain mutable iteration counters; give each worker its own
 * instance. Both methods align the complete query against an inclusive reference
 * window. The Short method uses query-width dynamic-programming arrays; its name
 * does not specify numeric precision, which the JNI wrapper selects separately.
 * Caller-provided bounds and a nonempty query are assumed rather than validated.
 * @author Brian Bushnell, Sucrose
 */
public abstract class IceCreamAligner{

	/**
	 * Selects the current backend. JNI selection ignores bits and initializes the
	 * native implementation; otherwise only the 32-bit Java implementation exists.
	 * @param bits Requested Java precision when JNI is disabled
	 * @return New backend instance with its own counters
	 * @throws RuntimeException If JNI is disabled and bits is not 32
	 */
	public static IceCreamAligner makeAligner(int bits){
		if(Shared.USE_JNI){
			return new IceCreamAlignerJNI();
		}else if(bits==32){
			return new IceCreamAlignerJava();
		}else{
			throw new RuntimeException(""+bits);
		}
	}

	/**
	 * Aligns with reference-window-width working arrays and row-based early exits.
	 * The reference window must be at least as long as the query for initialization.
	 * Supplied implementations rescale the raw score before applying minRatio.
	 * @param query Nonempty complete query, read without mutation
	 * @param ref Reference bases, read without mutation
	 * @param rstart Inclusive window start, nonnegative
	 * @param rstop Inclusive window stop, less than ref.length and not before rstart
	 * @param minScore Aggressive early-exit threshold on a query row's best raw score
	 * @param minRatio Required ratio of rescaled score to query length
	 * @return Alignment with absolute reference endpoint and legacy maxQpos=query.length,
	 * or null after an early exit or a failing ratio check
	 */
	public abstract AlignmentResult alignForward(final byte[] query, final byte[] ref, final int rstart, final int rstop, final int minScore,
			final float minRatio);

	/**
	 * Aligns with query-width working arrays, selecting the best full-query endpoint
	 * over reference positions. Current Java and native short algorithms do not use
	 * minScore; retain that distinction from alignForward's early-exit behavior.
	 * @param query Nonempty complete query, read without mutation
	 * @param ref Reference bases, read without mutation
	 * @param rstart Inclusive window start, nonnegative
	 * @param rstop Inclusive window stop, less than ref.length and not before rstart
	 * @param minScore Compatibility argument, unused by the supplied short algorithms
	 * @param minRatio Required ratio of rescaled score to query length
	 * @return Alignment with absolute reference endpoint and legacy maxQpos=query.length-1,
	 * or null when the rescaled score ratio is below minRatio
	 */
	public abstract AlignmentResult alignForwardShort(final byte[] query, final byte[] ref, final int rstart, final int rstop, final int minScore,
			final float minRatio);

//	public static final int pointsMatch = 1;
//	public static final int pointsSub = -1;
//	public static final int pointsDel = -2;
//	public static final int pointsIns = -2;

	/**
	 * Returns the backend's accumulated full-alignment cell counter.
	 * Early-exit accounting differs between the Java and JNI implementations.
	 * @return Recorded cell count; not a timing measurement
	 */
	abstract long iters();
	/** Returns the backend's accumulated short-alignment cell counter.
	 * @return Recorded cell count; not a timing measurement
	 */
	abstract long itersShort();

}
