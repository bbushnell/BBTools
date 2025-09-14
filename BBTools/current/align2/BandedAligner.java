package align2;

import java.util.Arrays;

import shared.Shared;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Aug 5, 2013
 *
 */
public abstract class BandedAligner {
	
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
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public abstract int alignForward(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public abstract int alignForwardRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public abstract int alignReverse(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	public abstract int alignReverseRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact);
	
	protected void fillBig(int[] array){
		final int lim=array.length-1;
		for(int i=1; i<lim; i++){array[i]=big;}
	}
	
	/** Score is lastRow-edits */
	public final int score(){
		return lastRow-lastEdits+1;
	}
	
	/** Position of min value in array (meaning the best alignment) relative to the middle of the array. */
	protected int lastOffset(int[] array, int halfWidth){
		final int center=halfWidth+1;
		int minLoc=center;
		for(int i=1; i<=halfWidth; i++){
			if(array[center+i]<array[minLoc]){minLoc=center+i;}
			if(array[center-i]<array[minLoc]){minLoc=center-i;}
		}
		return center-minLoc;
	}
	
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
	
	/** Final row aligned in last alignment. */
	public int lastRow;
	/** Final edits value in last alignment. */
	public int lastEdits;

	/** Position of min value in array (meaning the best alignment) relative to the middle of the array.
	 * Positive value is to the right (ref sequence longer than query), negative value left (ref shorter than query) */
	protected int lastOffset;
	
	public int lastRefLoc;
	public int lastQueryLoc;
	
	public final int maxWidth;
	
	public static final int big=99999999;
	public static boolean verbose=false;
	/** Penalizes non-length-neutral alignments.
	 * This causes query-to-ref alignment to yield same score as ref-to-query alignment, which is useful for assertions.  */
	public static boolean penalizeOffCenter=true;
	
}
