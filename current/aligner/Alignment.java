package aligner;

import prok.GeneCaller;
import stream.Read;

/**
 * Creates an Alignment wrapper for a Read object integrating with SingleStateAlignerFlat2
 * backend for high-level alignment operations with caching and sorting capabilities.
 * Stores alignment results including identity score, match string, and position coordinates.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class Alignment implements Comparable<Alignment>{
	
	/** Creates an Alignment wrapper for the specified Read object.
	 * @param r_ The Read object to wrap for alignment operations */
	public Alignment(Read r_){
		r=r_;
	}
	
	@Override
	public int compareTo(Alignment o) {
		return id>o.id ? 1 : id<o.id ? -1 : r.length()>o.r.length() ? 1 : r.length()<o.r.length() ? -1 : 0;
	}
	
	/**
	 * Aligns the wrapped read against the reference sequence and caches results.
	 * Uses SingleStateAlignerFlat2 for alignment computation and stores identity,
	 * match string, and position coordinates in instance variables.
	 *
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
	 * Static alignment method performing complete alignment pipeline using SingleStateAlignerFlat2.
	 * Executes fillUnlimited matrix computation, score boundary calculation, and traceback
	 * generation to compute precise identity and update Read position coordinates.
	 *
	 * @param r Read object to align (position coordinates will be updated)
	 * @param ref Reference sequence to align against
	 * @return Identity score calculated from match string using Read.identity
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
	
	/** The Read object being aligned */
	public final Read r;
	/** Identity score from alignment, -1 if not yet computed */
	public float id=-1;
	/** Match string encoding alignment operations (matches, mismatches, gaps) */
	public byte[] match;
	/** Start position of alignment on reference sequence */
	public int start;
	/** Stop position of alignment on reference sequence */
	public int stop;
	
}
