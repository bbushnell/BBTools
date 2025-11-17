package icecream;

import aligner.AlignmentResult;
import shared.Shared;

/**
 * Abstract base class for sequence alignment algorithms with configurable
 * implementations. Provides factory method to select alignment strategy (JNI or
 * Java) and defines abstract alignment methods for forward alignments with score
 * and ratio thresholds.
 *
 * @author Brian Bushnell
 */
public abstract class IceCreamAligner {

	/**
	 * Factory method that creates an appropriate aligner implementation based on
	 * configuration. Returns JNI implementation if USE_JNI is enabled, otherwise
	 * returns Java implementation for 32-bit alignment.
	 *
	 * @param bits Bit configuration for alignment (currently only 32-bit supported
	 * for Java implementation)
	 * @return IceCreamAlignerJNI if JNI enabled, IceCreamAlignerJava for 32-bit,
	 * throws RuntimeException for unsupported bit configurations
	 */
	public static IceCreamAligner makeAligner(int bits){
		if(Shared.USE_JNI){
			return new IceCreamAlignerJNI();
		}else if(bits==32){
			return new IceCreamAlignerJava();
		}
		else{
			throw new RuntimeException(""+bits);
		}
	}

	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @param rstop
	 * @param minScore Quit early if score drops below this
	 * @param minRatio Don't return results if max score is less than this fraction of max possible score
	 * @return
	 */
	public abstract AlignmentResult alignForward(final byte[] query, final byte[] ref, final int rstart, final int rstop, final int minScore,
			final float minRatio);

	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @param rstop
	 * @param minScore Quit early if score drops below this
	 * @param minRatio Don't return results if max score is less than this fraction of max possible score
	 * @return
	 */
	public abstract AlignmentResult alignForwardShort(final byte[] query, final byte[] ref, final int rstart, final int rstop, final int minScore,
			final float minRatio);

//	public static final int pointsMatch = 1;
//	public static final int pointsSub = -1;
//	public static final int pointsDel = -2;
//	public static final int pointsIns = -2;

	/**
	 * Returns the number of iterations performed by the standard alignment algorithm.
	 * Used for performance monitoring and optimization analysis.
	 * @return Number of alignment iterations executed
	 */
	abstract long iters();
	/**
	 * Returns the number of iterations performed by the short alignment algorithm.
	 * Used for performance monitoring and optimization analysis.
	 * @return Number of short alignment iterations executed
	 */
	abstract long itersShort();

}