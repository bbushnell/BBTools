package icecream;

import aligner.AlignmentResult;
import shared.Shared;

public abstract class IceCreamAligner {

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

	abstract long iters();
	abstract long itersShort();

}