package aligner;

import prok.GeneCaller;
import stream.Read;

/**
 * Alignment wrapper for a Read that uses SingleStateAlignerFlat2 for alignment
 * while caching identity, match string, and start/stop coordinates for reuse
 * and sorting. Intended for high-level alignment workflows.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class Alignment implements Comparable<Alignment>{
	
	public Alignment(Read r_){
		r=r_;
	}
	
	/**
	 * Compares alignments by identity score first, then by read length so higher
	 * quality alignments sort ahead of lower ones.
	 * @param o The alignment to compare against
	 * @return Negative if this alignment is lower quality, positive if higher, zero if equal
	 */
	@Override
	public int compareTo(Alignment o) {
		return id>o.id ? 1 : id<o.id ? -1 : r.length()>o.r.length() ? 1 : r.length()<o.r.length() ? -1 : 0;
	}
	
	/**
	 * Aligns the wrapped read against the reference sequence using SingleStateAlignerFlat2,
	 * caching identity, match string, and position coordinates in this object.
	 * @param ref Reference sequence to align against
	 * @return Identity score between 0.0 and 1.0
	 */
	public float align(byte[] ref){
		id=align(r, ref);
		match=r.match;
		start=r.start;
		stop=r.stop;
		return id;
	}
	
	/**
	 * Performs the full alignment pipeline using SingleStateAlignerFlat2: fillUnlimited,
	 * score calculation, and traceback to compute identity and update the read's
	 * position coordinates.
	 *
	 * @param r Read object to align (position coordinates will be updated)
	 * @param ref Reference sequence to align against
	 * @return Identity score calculated from the match string using Read.identity
	 */
	public static final float align(Read r, byte[] ref){
		SingleStateAlignerFlat2 ssa=GeneCaller.getSSA();
		final int a=0, b=ref.length-1;
		int[] max=ssa.fillUnlimited(r.bases, ref, a, b, 0);
		if(max==null){return 0;}
		
		final int rows=max[0];
		final int maxCol=max[1];
		final int maxState=max[2];
		
		//returns {score, bestRefStart, bestRefStop} 
		//padded: {score, bestRefStart, bestRefStop, padLeft, padRight};
		int[] score=ssa.score(r.bases, ref, a, b, rows, maxCol, maxState);
		int rstart=score[1];
		int rstop=score[2];
		r.start=rstart;
		r.stop=rstop;
		
		byte[] match=ssa.traceback(r.bases, ref, a, b, rows, maxCol, maxState);
		float id=Read.identity(match);
		return id;
	}
	
	public final Read r;
	public float id=-1;
	public byte[] match;
	public int start;
	public int stop;
	
}
