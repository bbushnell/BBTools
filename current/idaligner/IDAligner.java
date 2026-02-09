package idaligner;

import java.util.Arrays;

/**
 * Interface for aligners that can calculate pairwise sequence identity.
 * Provides methods for aligning sequences and returning identity scores
 * as floating-point values between 0.0 and 1.0.
 *
 * @author Brian Bushnell
 * @contributor Isla, Collei
 * @date April 24, 2025
 */
public interface IDAligner {
	
	/** Returns the name identifier of this aligner implementation.
	 * @return Aligner name */
	public String name();
	
	/**
	 * Performs sequence alignment and calculates pairwise identity.
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @return Identity score between 0.0 and 1.0
	 */
	public float align(byte[] q, byte[] r);
	
	/**
	 * Performs sequence alignment with optional position tracking.
	 * If posVector is null, sequences may be swapped to make query shorter.
	 *
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @param posVector Optional int[2] array for returning {rStart, rStop}
	 * of the optimal alignment
	 * @return Identity score between 0.0 and 1.0
	 */
	public float align(byte[] q, byte[] r, int[] posVector);
	
	/**
	 * Lightweight wrapper for aligning to a window of the reference sequence.
	 * Extracts reference region and adjusts coordinates in result.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param pos Optional int[2] for returning {rStart, rStop} of optimal alignment
	 * @param rStart Alignment window start position
	 * @param rStop Alignment window end position
	 * @return Identity(0.0-1.0)
	 */
	public default float align(final byte[] query, final byte[] ref, 
			final int[] pos, int rStart, int rStop){
		rStart=Math.max(rStart, 0);
		rStop=Math.min(rStop, ref.length-1);
		final int rlen=rStop-rStart+1;
		final byte[] region=(rlen==ref.length ? ref : Arrays.copyOfRange(ref, rStart, rStop));
		final float id=align(query, region, pos);
		if(pos!=null){
			assert(pos[1]>0) : id+", "+Arrays.toString(pos)+", "+rStart;
			pos[0]+=rStart;
			pos[1]+=rStart;
		}
		return id;
	}
	
	/**
	 * Performs sequence alignment with minimum score threshold.
	 * If posVector is null, sequences may be swapped to make query shorter.
	 *
	 * @param q Query sequence as byte array
	 * @param r Reference sequence as byte array
	 * @param posVector Optional int[2] array for returning {rStart, rStop}
	 * of the optimal alignment
	 * @param minScore Legacy field to allow early exit in some implementations
	 * @return Identity score between 0.0 and 1.0
	 */
	public float align(byte[] q, byte[] r, int[] posVector, int minScore);
	
	/** Returns the number of alignment loops performed.
	 * @return Loop count as a long value */
	public long loops();
	
	/**
	 * Sets the loop counter for alignment operations.
	 * Typically used for bookkeeping or resetting between runs.
	 * @param i Loop count value to set
	 */
	public void setLoops(long i);
	
	/*--------------------------------------------------------------*/
	/*----------------        AlignmentStats        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Lightweight wrapper for aligning to a window of the reference sequence.
	 * Extracts reference region and adjusts coordinates in result.
	 *
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param stats Optional container for returning alignment information
	 * @param rStart Alignment window start position
	 * @param rStop Alignment window end position
	 * @return Identity(0.0-1.0)
	 */
	public default float align(final byte[] query, final byte[] ref, 
			final AlignmentStats stats, int rStart, int rStop){
		rStart=Math.max(rStart, 0);
		rStop=Math.min(rStop, ref.length-1);
		final int rlen=rStop-rStart+1;
		final byte[] region=(rlen==ref.length ? ref : Arrays.copyOfRange(ref, rStart, rStop));
		final float id=align(query, region, stats);
		if(stats!=null){
			assert(stats.rStop>0) : id+", "+stats.toString()+", "+rStart;
			stats.rStart+=rStart;
			stats.rStop+=rStart;
		}
		return id;
	}
	
	/**
	 * Retrofit swapping AlignmentStats in for posVector.
	 * Performs sequence alignment within a specified reference window.
	 * If stats is null, sequences may be swapped to make query shorter.
	 *
	 * @param q Query sequence
	 * @param r Ref sequence
	 * @param stats Optional AlignmentStats object for returning information
	 * @return Identity (0-1)
	 */
	public default float align(byte[] q, byte[] r, AlignmentStats stats) {
		int[] pos=null;
		if(stats!=null) {
			stats.clear();
			pos=new int[4];
			assert(!stats.doTrace) : "Class "+getClass()+" needs to override this method.";
		}
		float id=align(q, r, pos);
		if(stats!=null) {
			stats.setFromPos(pos, id);
			stats.qLen=q.length;
			stats.rLen=r.length;
			stats.solve();
		}
		return id;
	}
	
}
