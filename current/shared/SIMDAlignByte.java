package shared;

import java.util.Arrays;

import dna.AminoAcid;
import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.FloatVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.ShortVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import structures.IntList;

/** 
 * Holds SIMD methods for alignment.
 * @author Brian Bushnell
 * @author Isla
 * @contributor Zephy
 * @date May 12, 2025
 *
 */
public class SIMDAlignByte {

	//Example from https://medium.com/@Styp/java-18-vector-api-do-we-get-free-speed-up-c4510eda50d2
	@SuppressWarnings("restriction")
	private static final VectorSpecies<Float> FSPECIES=FloatVector.SPECIES_256;//FloatVector.SPECIES_PREFERRED; //This needs to be final or performance drops.
	private static final int FWIDTH=FSPECIES.length();
	//		private static final int boundMask=~(FWIDTH-1);

	@SuppressWarnings("restriction")
	private static final VectorSpecies<Byte> BSPECIES=ByteVector.SPECIES_256;
	private static final int BWIDTH=BSPECIES.length();

	@SuppressWarnings("restriction")
	private static final VectorSpecies<Integer> ISPECIES=IntVector.SPECIES_256;
	private static final int IWIDTH=ISPECIES.length();

	@SuppressWarnings("restriction")
	private static final VectorSpecies<Short> SSPECIES=ShortVector.SPECIES_256;
	private static final int SWIDTH=SSPECIES.length();

	@SuppressWarnings("restriction")
	private static final VectorSpecies<Double> DSPECIES=DoubleVector.SPECIES_256;
	private static final int DWIDTH=DSPECIES.length();

	@SuppressWarnings("restriction")
	private static final VectorSpecies<Long> LSPECIES=LongVector.SPECIES_256;
	private static final int LWIDTH=LSPECIES.length();

	private static final byte MATCH=1, N_SCORE=0, SUB=-1, INS=-1, DEL=-1;

	public static void alignBandVector(byte q, byte[] ref, int bandStart, int bandEnd, 
			byte[] prev, byte[] curr) {

		// Calculate range size
		final int rangeSize = bandEnd - bandStart + 1;

		// Broadcast constants to vectors
		final ByteVector vq = ByteVector.broadcast(BSPECIES, q);
		final ByteVector vMATCH = ByteVector.broadcast(BSPECIES, MATCH);
		final ByteVector vN_SCORE = ByteVector.broadcast(BSPECIES, N_SCORE);
		final ByteVector vSUB = ByteVector.broadcast(BSPECIES, SUB);
		final ByteVector vINS = vSUB;
		final ByteVector v15 = ByteVector.broadcast(BSPECIES, (byte)15);

		// If range is at least one vector wide
		if(rangeSize >= BWIDTH) {
			// Process all complete vector blocks except the final one
			for(int j=bandStart; j<bandEnd-(BWIDTH-2); j+=BWIDTH) {
				// Load reference values
				ByteVector refVec = ByteVector.fromArray(BSPECIES, ref, j-1);

				// Match condition calculation
				VectorMask<Byte> isMatchMask = vq.compare(VectorOperators.EQ, refVec);

				// Check for N values
				VectorMask<Byte> hasNMask = vq.or(refVec).compare(VectorOperators.GE, v15);

				// Combined blend operations
				ByteVector scoreAddVec = vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

				// Load previous scores
				ByteVector prevJ1Vec = ByteVector.fromArray(BSPECIES, prev, j-1);
				ByteVector prevJVec = ByteVector.fromArray(BSPECIES, prev, j);

				// Calculate scores and find max
				ByteVector diagScoreVec = prevJ1Vec.add(scoreAddVec);
				ByteVector upScoreVec = prevJVec.add(vINS);
				ByteVector maxVec = diagScoreVec.max(upScoreVec);

				// Store results
				maxVec.intoArray(curr, j);
			}

			// Process the final block with a vector, aligned to include the last element
			final int lastVectorStart = Math.max(bandStart, bandEnd+1-BWIDTH);

			// Load reference values for final vector
			ByteVector refVec = ByteVector.fromArray(BSPECIES, ref, lastVectorStart-1);

			// Match condition calculation
			VectorMask<Byte> isMatchMask = vq.compare(VectorOperators.EQ, refVec);

			// Check for N values
			VectorMask<Byte> hasNMask = vq.or(refVec).compare(VectorOperators.GE, v15);

			// Combined blend operations
			ByteVector scoreAddVec = vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

			// Load previous scores
			ByteVector prevJ1Vec = ByteVector.fromArray(BSPECIES, prev, lastVectorStart-1);
			ByteVector prevJVec = ByteVector.fromArray(BSPECIES, prev, lastVectorStart);

			// Calculate scores and find max
			ByteVector diagScoreVec = prevJ1Vec.add(scoreAddVec);
			ByteVector upScoreVec = prevJVec.add(vINS);
			ByteVector maxVec = diagScoreVec.max(upScoreVec);

			// Store results
			maxVec.intoArray(curr, lastVectorStart);
		} else {
			// Range too small for vector processing, use scalar code
			for(int j=bandStart; j<=bandEnd; j++) {
				final byte r = ref[j-1];
				final boolean isMatch = (q==r);
				final boolean hasN = ((q|r)>=15);
				final byte scoreAdd = isMatch ? MATCH : (hasN ? N_SCORE : SUB);

				final byte pj1 = prev[j-1], pj = prev[j];
				final byte diagScore = (byte)(pj1+scoreAdd);
				final byte upScore = (byte)(pj+INS);

				curr[j] = (byte)Math.max(diagScore, upScore);
			}
		}
	}

	/**
	 * Counts substitutions (mismatches) between a query sequence and multiple diagonals 
	 * of a reference sequence using SIMD vectorization.
	 * 
	 * This function compares the query sequence against BWIDTH (32) different starting 
	 * positions in the reference sequence simultaneously, counting mismatches for each.
	 * It's designed for use in sequence alignment algorithms where multiple diagonal
	 * alignments need to be evaluated.
	 *
	 * @param query The query sequence as a byte array
	 * @param ref The reference sequence as a byte array  
	 * @param pos Output parameter - will contain the lane index (0-31) with minimum mismatches
	 * @return The minimum number of mismatches found across all lanes, or qLen if all lanes overflow
	 */
	public static final int countSubs(byte[] query, byte[] ref, int[] pos, int maxSubs){
		final int qLen=query.length;  // Cache query length
		final int rLen=ref.length;    // Cache reference length

		// Need at least BWIDTH positions to make SIMD worthwhile
		if(rLen<BWIDTH || qLen==0 || pos[0]==0 || countSubs(query, ref, pos, maxSubs)<10){//This could actually decide based on a recursive call to maxSubs/2
			// Fall back to scalar implementation for small inputs
			int subs=0;  // Mismatch counter
			// Compare up to the shorter of the two sequences
			for(int i=0, minlen=Math.min(qLen, rLen); i<minlen && subs<maxSubs; i++){
				subs+=(query[i]==ref[i] ? 0 : 1);  // Increment if mismatch found
			}
			pos[0]=0;  // Only one "lane" in scalar mode
			return subs;
		}
		pos[0]=0;

		// Initialize scores vector with -128 (allows counting up to 127 mismatches)
		ByteVector scores=ByteVector.broadcast(BSPECIES, (byte)-128);
		// Vector of ones for incrementing scores
		ByteVector ones=ByteVector.broadcast(BSPECIES, (byte)1);
		// Vector of -1 for overflow detection
		ByteVector minusOne=ByteVector.broadcast(BSPECIES, (byte)-1);
		// Mask tracking which lanes have overflowed
		VectorMask<Byte> overflowMask=scores.maskAll(false);

		// Process up to where all BWIDTH diagonals fit within reference
		int limit=Math.min(qLen, rLen-BWIDTH+1);
		int start=0;  // Current position in query

		// Main processing loop - processes 64 query positions at a time
		while(start<limit){
			int processed = 0;  // Track how many positions we actually process
			// Inner loop processes up to 64 positions before checking overflow
			for(int i=0; i<64 && (start+i)<limit; i++){
				// Broadcast current query byte to all lanes
				ByteVector q=ByteVector.broadcast(BSPECIES, query[start+i]);
				// Load BWIDTH reference bytes starting at different offsets
				ByteVector r=ByteVector.fromArray(BSPECIES, ref, start+i);
				// Compare query byte with reference bytes, creating mismatch mask
				VectorMask<Byte> mismatchMask=q.compare(VectorOperators.NE, r);
				// Increment scores where mismatches occurred
				scores=scores.add(ones, mismatchMask);
				processed++;  // Increment counter
			}

			// Check for new overflows (scores > -1 means they wrapped around)
			VectorMask<Byte> newOverflows=scores.compare(VectorOperators.GT, minusOne);
			// Update cumulative overflow mask
			overflowMask=overflowMask.or(newOverflows);

			// If all lanes overflowed, return query length
			if(overflowMask.trueCount()==BWIDTH){
				return qLen;
			}

			start+=processed;  // Only increment by positions actually processed
		}
		
		// Handle remaining query positions
		if(start<qLen){
		    int remaining=qLen-start;
		    // Only need BWIDTH + remaining - 1 bytes for the slice
		    byte[] refSlice = new byte[BWIDTH + remaining - 1];
		    
		    // Pre-fill the slice with reference data or 'X' for out-of-bounds
		    for(int j=0; j<refSlice.length; j++){
		        int refPos = start+j;
		        refSlice[j] = (refPos < rLen) ? ref[refPos] : (byte)'X';
		    }

		    // Process remaining positions in chunks of 64
		    for(int outer=0; outer<remaining; outer+=64){
		        for(int i=outer; i<Math.min(outer+64, remaining); i++){
		            ByteVector q=ByteVector.broadcast(BSPECIES, query[start+i]);
		            ByteVector r=ByteVector.fromArray(BSPECIES, refSlice, i);

		            VectorMask<Byte> mismatchMask=q.compare(VectorOperators.NE, r);
		            scores=scores.add(ones, mismatchMask);
		        }

		        // Check for overflows after each chunk of 64
		        VectorMask<Byte> newOverflows=scores.compare(VectorOperators.GT, minusOne);
		        overflowMask=overflowMask.or(newOverflows);
		        if(overflowMask.trueCount()==BWIDTH){
		            return qLen;
		        }
		    }
		}

		// Find minimum score among valid (non-overflowed) lanes
		VectorMask<Byte> validLanes=overflowMask.not();  // Invert mask to get valid lanes
		// Reduce to find minimum score in valid lanes
		long minScoreLong=scores.reduceLanesToLong(VectorOperators.MIN, validLanes);
		byte minScore=(byte)minScoreLong;  // Convert back to byte

		// Create vector with lane indices [0, 1, 2, ..., 31]
		ByteVector indices = ByteVector.broadcast(BSPECIES, (byte)0).addIndex(1);

		// Find which lanes have the minimum score
		VectorMask<Byte> minLanes=scores.compare(VectorOperators.EQ, minScore);
		// AND with valid lanes to exclude overflowed lanes
		VectorMask<Byte> validMinLanes=minLanes.and(validLanes);

		// Find the position of the first minimum lane
		// Blend in 127 for lanes that aren't valid minimums
		ByteVector maskedIndices = indices.blend((byte)127, validMinLanes.not());

		// Find minimum index among valid minimum lanes
		long firstMinIndexLong = maskedIndices.reduceLanesToLong(VectorOperators.MIN);
		// Extract lane index as int
		int lowLane = (int)(firstMinIndexLong & 0xFF);

		pos[0]=lowLane;  // Store the lane index in output parameter
		// Return the score, capping at qLen if overflow
		return minScore+128;
	}
	
	/**
	 * Aligns a sequence to a reference, with no indels.
	 * Returns all positions containing fewer than maxSubs mismatches.
	 * Appears broken due to lack of residual handling and AOOB exceptions.
	 */
	@Deprecated
	public static final IntList alignDiagonal(byte[] query, byte[] ref, int maxSubs){
		assert(maxSubs<256) : "This method can only handle up to 255 substitutions.";
	    final int qLen=query.length;
	    final int rLen=ref.length;
	    IntList results=null;
	    
	    // Only return things under limit
	    ByteVector limitV=ByteVector.broadcast(BSPECIES, (byte)(-128+maxSubs));
	    // Vector of ones for incrementing scores
	    ByteVector oneV=ByteVector.broadcast(BSPECIES, (byte)1);
	    
	    // Outer loop - slide window across reference
	    for(int start=0; start<=rLen-qLen; start+=BWIDTH){
	        // Initialize scores for this window
	        ByteVector scoreV=ByteVector.broadcast(BSPECIES, (byte)-128);
	        
	        // How many lanes are valid for this window
	        int validLanes = Math.min(BWIDTH, rLen-qLen-start+1);
	        
	        // Inner loop - process query
	        VectorMask<Byte> underLimitMask=null;
	        for(int i=0; i<qLen; i++){
	            ByteVector q=ByteVector.broadcast(BSPECIES, query[i]);
	            ByteVector r=ByteVector.fromArray(BSPECIES, ref, start+i);
	            VectorMask<Byte> mismatchMask=q.compare(VectorOperators.NE, r);
	            scoreV=scoreV.add(oneV, mismatchMask);
	            
	            // Early exit if all valid lanes exceeded limit
	            // This is the key optimization - one compare, one mask operation
	            underLimitMask=scoreV.compare(VectorOperators.LE, limitV);
	            if(!underLimitMask.anyTrue()){
	                break; // All lanes have too many mismatches
	            }
	        }
	        
	        // Collect results - positions with score <= limitV
	        VectorMask<Byte> successMask=scoreV.compare(VectorOperators.LE, limitV);
	        
	        // Add successful positions to results
	        if(underLimitMask!=null && underLimitMask.anyTrue()) {
	        	if(results==null) {results=new IntList();}
	        	for(int lane=0; lane<validLanes; lane++){
	        		if(successMask.laneIsSet(lane)){
	        			results.add(start + lane);
	        		}
	        	}
	        }
	    }
	    
	    return results;
	}
	
	/**
	 * Aligns a sequence to a reference, with no indels, supporting clipped alignments.
	 * Returns all positions containing fewer than maxSubs mismatches.
	 * Uses scalar loops for boundary conditions and SIMD for bulk processing.
	 * @TODO Pad short contigs to avoid spending time in scalar mode
	 */
	public static final IntList alignDiagonal(byte[] query, byte[] ref, int maxSubs, int maxClips){
		assert(maxSubs<256) : "This method can only handle up to 255 substitutions.";
		final int qLen=query.length;
		final int rLen=ref.length;
		IntList results=null;
//		System.err.println("\talignDiagonal(ql="+query.length+", rl="+ref.length+", ms="+maxSubs+", mc="+maxClips);
		
		// Calculate alignment range
		final int startPos=-maxClips;
		final int endPos=rLen-qLen+maxClips;
		
		// Vector constants
		ByteVector limitV=ByteVector.broadcast(BSPECIES, (byte)(-128+maxSubs));
		ByteVector oneV=ByteVector.broadcast(BSPECIES, (byte)1);
		
		int pos=startPos;
		
		// Scalar pre-loop for negative positions (left clips)
		while(pos<0 && pos<=endPos){
			int subs=alignClippedScalar(query, ref, maxSubs, maxClips, pos);
			if(subs<=maxSubs){
				if(results==null){results=new IntList();}
				results.add(pos);
			}
//			System.err.println("A: pos="+pos+", subs="+subs);
			pos++;
		}
		
		// SIMD main loop - only process when we have enough reference data
		final int maxSimdPos=rLen-qLen-BWIDTH+1;
		while(pos<=endPos && pos<=maxSimdPos){
			int validLanes=Math.min(BWIDTH, endPos-pos+1);
			if(validLanes<=0){break;}
			
			ByteVector scoreV=ByteVector.broadcast(BSPECIES, (byte)-128);
			VectorMask<Byte> underLimitMask=null;
			
			for(int i=0; i<qLen; i++){
				ByteVector q=ByteVector.broadcast(BSPECIES, query[i]);
				ByteVector r=ByteVector.fromArray(BSPECIES, ref, pos+i);
				VectorMask<Byte> mismatchMask=q.compare(VectorOperators.NE, r);
				scoreV=scoreV.add(oneV, mismatchMask);
				
				underLimitMask=scoreV.compare(VectorOperators.LE, limitV);
				if(!underLimitMask.anyTrue()){
					break;
				}
			}
			
			if(underLimitMask!=null && underLimitMask.anyTrue()){
				VectorMask<Byte> successMask=scoreV.compare(VectorOperators.LE, limitV);
				if(results==null){results=new IntList();}
				for(int lane=0; lane<validLanes; lane++){
					if(successMask.laneIsSet(lane)){
						results.add(pos+lane);
					}
//					System.err.println("B: pos="+(pos+lane)+", subs="+(int)scoreV.lane(lane));
				}
			}
			pos+=BWIDTH;
		}
		
		// Scalar post-loop for remaining positions
		while(pos<=endPos){
			int subs=alignClippedScalar(query, ref, maxSubs, maxClips, pos);
			if(subs<=maxSubs){
				if(results==null){results=new IntList();}
				results.add(pos);
			}
//			System.err.println("C: pos="+pos+", subs="+subs);
			pos++;
		}
		
//		System.err.println("Returning "+results);
		return results;
	}

	// Helper method for scalar alignment without clipping
	private static int alignScalar(byte[] query, byte[] ref, int maxSubs, int rStart) {
	    int subs = 0;
	    for(int i = 0, j = rStart; i < query.length && subs <= maxSubs; i++, j++) {
	        final byte q = query[i], r = ref[j];
	        if(q != r || AminoAcid.baseToNumber[q] < 0) {
	            subs++;
	        }
	    }
	    return subs;
	}

	// Helper method for scalar clipped alignment with proper clipping penalties
	private static int alignClippedScalar(byte[] query, byte[] ref, int maxSubs, int maxClips, int rStart){
		final int rStop1=rStart+query.length;
		final int leftClip=Math.max(0, -rStart);
		final int rightClip=Math.max(0, rStop1-ref.length);
		int clips=leftClip+rightClip;
		if(clips>=query.length){return query.length;}
		int subs=Math.max(0, clips-maxClips); // Excess clipping counts as substitutions
		int i=leftClip, j=rStart+leftClip;
		
		for(final int limit=Math.min(rStop1, ref.length); j<limit && subs<=maxSubs; i++, j++){
			final byte q=query[i], r=ref[j];
			final int incr=(q!=r || AminoAcid.baseToNumber[q]<0 ? 1 : 0);
			subs+=incr;
		}
		return subs;
	}

	//	public static void processDeletionsTailVector(byte[] curr, int bandStart, int bandEnd) {
	//	    final int maxDist = BWIDTH/2;
	//	    
	//	    // Initialize with scalar code up to the point where vector access is safe
	//	    byte leftCell = curr[bandStart-1];
	//	    for(int j = bandStart; j < bandStart+maxDist; j++) {
	//	        if (j <= bandEnd) {
	//	            final byte maxDiagUp = curr[j];
	//	            final byte leftScore = (byte)(leftCell + DEL);
	//	            leftCell = (byte)Math.max(maxDiagUp, leftScore);
	//	            curr[j] = leftCell;
	//	        }
	//	    }
	//	    
	//	    // Main vector processing loop - only when we have safe access
	//	    for(int j = bandStart+maxDist; j <= bandEnd-(BWIDTH-1); j += BWIDTH) {
	//	        // Load the current vector
	//	        ByteVector currVec = ByteVector.fromArray(BSPECIES, curr, j);
	//	        
	//	        // For each vector position, apply shell sort pattern of decreasing distances
	//	        for(int dist = maxDist; dist > 0; dist /= 2) {
	//	            // Load values from distance positions back
	//	            ByteVector leftVec = ByteVector.fromArray(BSPECIES, curr, j-dist);
	//	            
	//	            // Add appropriate penalty (DEL*dist)
	//	            ByteVector penaltyVec = ByteVector.broadcast(BSPECIES, (byte)(DEL*dist));
	//	            ByteVector leftScoreVec = leftVec.add(penaltyVec);
	//	            
	//	            // Update current vector with maximums
	//	            currVec = currVec.max(leftScoreVec);
	//	        }
	//	        
	//	        // Store updated vector
	//	        currVec.intoArray(curr, j);
	//	    }
	//	    
	//	    // Handle remaining elements with scalar code
	//	    for(int j = bandEnd-(bandEnd-bandStart+1)%BWIDTH+1; j <= bandEnd; j++) {
	//	        // Apply the same shell sort approach to remaining elements
	//	        for(int dist = maxDist; dist > 0; dist /= 2) {
	//	            if(j-dist >= bandStart) {
	//	                byte leftValue = curr[j-dist];
	//	                byte leftScore = (byte)(leftValue + DEL*dist);
	//	                curr[j] = (byte)Math.max(curr[j], leftScore);
	//	            }
	//	        }
	//	    }
	//	}

	public static void processDeletionsTailVector(byte[] curr, int bandStart, int bandEnd) {
		final int maxDist = BWIDTH/2;

		// Initialize with scalar code for the boundary conditions
		byte leftCell = curr[bandStart-1];
		int lastProcessed = bandStart - 1;  // Track what we've processed

		for(int j = bandStart; j < bandStart+maxDist && j <= bandEnd; j++) {
			final byte maxDiagUp = curr[j];
			final byte leftScore = (byte)(leftCell + DEL);
			leftCell = (byte)Math.max(maxDiagUp, leftScore);
			curr[j] = leftCell;
			lastProcessed = j;
		}

		// Main vector processing loop
		for(int j = bandStart; j <= bandEnd-(BWIDTH-1); j += BWIDTH) {
			// Load the current vector
			ByteVector currVec = ByteVector.fromArray(BSPECIES, curr, j);

			// For each vector position, apply shell sort pattern of decreasing distances
			for(int dist = maxDist; dist > 0; dist /= 2) {
				if(j-dist >= bandStart-1) {
					ByteVector leftVec = ByteVector.fromArray(BSPECIES, curr, j-dist);
					ByteVector penaltyVec = ByteVector.broadcast(BSPECIES, (byte)(DEL*dist));
					ByteVector leftScoreVec = leftVec.add(penaltyVec);
					currVec = currVec.max(leftScoreVec);
				}
			}

			// Store updated vector
			currVec.intoArray(curr, j);
			lastProcessed = j + BWIDTH - 1;  // We processed up through here
		}

		// Handle remaining elements with scalar code
		for(int j = lastProcessed + 1; j <= bandEnd; j++) {
			for(int dist = maxDist; dist > 0; dist /= 2) {
				if(j-dist >= bandStart-1) {
					byte leftValue = curr[j-dist];
					byte leftScore = (byte)(leftValue + DEL*dist);
					curr[j] = (byte)Math.max(curr[j], leftScore);
				}
			}
		}
	}

}
