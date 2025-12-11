package align2;

import java.util.Arrays;

import shared.Shared;
import shared.Tools;

/**
 * Abstract base class for banded sequence alignment algorithms.
 * Provides constrained alignment within a diagonal band to optimize speed and memory usage.
 * Supports both forward and reverse complement alignments with configurable width limits.
 *
 * @author Brian Bushnell
 * @date Aug 5, 2013
 */
public abstract class BandedAligner {
	
	/**
	 * Constructs a banded aligner with specified maximum band width.
	 * Ensures width is odd and at least 3 for proper diagonal band calculation.
	 * @param width_ Maximum width of the alignment band
	 */
	public BandedAligner(int width_){
		maxWidth=Tools.max(width_, 3)|1;
		assert(maxWidth>=3) : "width<3 : "+width_+" -> "+maxWidth;
		assert(big>maxWidth/2);
	}
	
	public static final BandedAligner makeBandedAligner(int width_){
		//TODO: Remove the false condition when BandedAlignerJNI yields identical results to BandedAlignerConcrete.
		BandedAligner ba=((Shared.USE_JNI && false) ? new BandedAlignerJNI(width_) : new BandedAlignerConcrete(width_));
		return ba;
	}
	
	/**
	 * Performs progressive quadruple alignment with increasing edit distance thresholds.
	 * Starts with minEdits limit and progressively increases by factor of 4 until maxEdits.
	 * Tests all four orientations (forward, reverse, forward RC, reverse RC) at each threshold.
	 *
	 * @param query Query sequence to align
	 * @param ref Reference sequence to align against
	 * @param minEdits Minimum edit distance to start testing
	 * @param maxEdits Maximum edit distance allowed
	 * @param exact Whether to require exact alignment within edit limit
	 * @return Best edit distance found across all orientations
	 */
	public final int alignQuadrupleProgressive(final byte[] query, final byte[] ref, int minEdits, int maxEdits, final boolean exact){
		maxEdits=Tools.min(maxEdits, Tools.max(query.length, ref.length));
		minEdits=Tools.min(minEdits, maxEdits);
		//System.err.println("maxEdits="+maxEdits+", "+minEdits);
		for(long i=minEdits, me=-1; me<maxEdits; i=i*4){
			me=Tools.min(i, maxEdits);
			if(me*2>maxEdits){me=maxEdits;}
			int edits=alignQuadruple(query, ref, (int)me, exact);
//			System.err.println("i="+i+", me="+me+", minEdits="+minEdits+", maxEdits="+maxEdits+", edits="+edits);
			if(edits<me){return edits;}
		}
		return maxEdits;
	}
	
	public final int alignQuadruple(final byte[] query, final byte[] ref, final int maxEdits, final boolean exact){
		final int a=alignForward(query, ref, 0, 0, maxEdits, exact);
		final int b=alignReverse(query, ref, query.length-1, ref.length-1, maxEdits, exact);
		final int me2=Tools.min(maxEdits, Tools.max(a, b));
		if(me2==0){return 0;}
		final int c=alignForwardRC(query, ref, query.length-1, 0, me2, exact);
		final int d=alignReverseRC(query, ref, 0, ref.length-1, me2, exact);
//		System.err.println("a="+a+", b="+b+", c="+c+", d="+d);
		return Tools.min(Tools.max(a, b), Tools.max(c, d));
	}
	
	public final int alignDouble(final byte[] query, final byte[] ref, final int maxEdits, final boolean exact){
		final int a=alignForward(query, ref, 0, 0, maxEdits, exact);
		if(a==0){return 0;}
		final int c=alignForwardRC(query, ref, query.length-1, 0, a, exact);
		return Tools.min(a, c);
	}
	
	public abstract int alignForward(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	public abstract int alignForwardRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	public abstract int alignReverse(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	public abstract int alignReverseRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	/**
	 * Fills array interior elements with large sentinel values.
	 * Preserves first and last elements while setting middle values to 'big'.
	 * Used to initialize alignment arrays before dynamic programming.
	 * @param array Array to initialize with sentinel values
	 */
	protected void fillBig(int[] array){
		final int lim=array.length-1;
		for(int i=1; i<lim; i++){array[i]=big;}
	}
	
	/**
	 * Calculates alignment score from the last alignment result.
	 * Score represents alignment quality: higher values indicate better alignments.
	 * @return Alignment score based on lastRow minus lastEdits plus one
	 */
	public final int score(){
		return lastRow-lastEdits+1;
	}
	
	/**
	 * Finds the position of minimum value in alignment array.
	 * Searches outward from center to find best alignment position within band.
	 *
	 * @param array Alignment scores array
	 * @param halfWidth Half-width of the search band
	 * @return Offset from center of the best alignment position
	 */
	protected int lastOffset(int[] array, int halfWidth){
		final int center=halfWidth+1;
		int minLoc=center;
		for(int i=1; i<=halfWidth; i++){
			if(array[center+i]<array[minLoc]){minLoc=center+i;}
			if(array[center-i]<array[minLoc]){minLoc=center-i;}
		}
		return center-minLoc;
	}
	
	/**
	 * Old version of off-center penalty function.
	 * Adds linear penalty to alignment scores based on distance from center.
	 * Deprecated in favor of penalizeOffCenter which uses max instead of addition.
	 *
	 * @param array Alignment scores array to modify
	 * @param halfWidth Half-width of the penalty band
	 * @return Minimum penalized score
	 */
	protected int penalizeOffCenter_old(int[] array, int halfWidth){
		if(verbose){
			System.err.println("penalizeOffCenter_old("+Arrays.toString(array)+", "+halfWidth);
		}
		final int center=halfWidth+1;
		int edits=array[center];
		for(int i=1; i<=halfWidth; i++){
			array[center+i]=Tools.min(big, array[center+i]+i);
			edits=Tools.min(edits, array[center+i]);
			array[center-i]=Tools.min(big, array[center-i]+i);
			edits=Tools.min(edits, array[center-i]);
		}
		if(verbose){
			System.err.println("returned "+edits);
		}
		return edits;
	}
	
	/**
	 * Applies penalty for alignments away from center diagonal.
	 * Uses max function to ensure minimum penalty based on distance from center.
	 * Prevents alignments that are heavily biased toward indels over matches.
	 *
	 * @param array Alignment scores array to modify
	 * @param halfWidth Half-width of the penalty band
	 * @return Minimum penalized score after applying off-center penalties
	 */
	protected int penalizeOffCenter(int[] array, int halfWidth){
		if(verbose){
			System.err.println("penalizeOffCenter("+Arrays.toString(array)+", "+halfWidth);
		}
		final int center=halfWidth+1;
		int edits=array[center];
		for(int i=1; i<=halfWidth; i++){
			array[center+i]=Tools.min(big, Tools.max(i, array[center+i]));
			edits=Tools.min(edits, array[center+i]);
			array[center-i]=Tools.min(big, Tools.max(i, array[center-i]));
			edits=Tools.min(edits, array[center-i]);
		}
		if(verbose){
			System.err.println("returned "+edits);
		}
		return edits;
	}
	
	/** Final row position reached in the last alignment operation */
	public int lastRow;
	/** Final edit distance calculated in the last alignment operation */
	public int lastEdits;

	/** Position offset of best alignment relative to center of band */
	protected int lastOffset;
	
	public int lastRefLoc;
	public int lastQueryLoc;
	
	public final int maxWidth;
	
	public static final int big=99999999;
	public static boolean verbose=false;
	/** Whether to apply penalties for alignments far from the center diagonal */
	public static boolean penalizeOffCenter=true;
	
}
