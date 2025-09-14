package shared;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.FloatVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.ShortVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/** 
 * Holds SIMD methods for alignment.
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 12, 2025
 *
 */
public class SIMDAlign {

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

	private static final long MATCH=1L<<42, N_SCORE=0, SUB=(-1L)<<42, INS=(-1L)<<42;
	private static final long DEL=((-1L)<<42)|((1L)<<21);
	private static final long SCORE_MASK=((-1L)<<42);
	private static final boolean debug=false;
	
	private static final LongVector vMATCH=LongVector.broadcast(LSPECIES, MATCH);
	private static final LongVector vN_SCORE=LongVector.broadcast(LSPECIES, N_SCORE);
	private static final LongVector vSUB=LongVector.broadcast(LSPECIES, SUB);
	private static final LongVector vINS=vSUB;//LongVector.broadcast(LSPECIES, INS);
	private static final LongVector v15=LongVector.broadcast(LSPECIES, 15);

	/** Designed to eliminate scalar loop by reprocessing some final elements */
	public static void alignBandVector(long q, long[] ref, int bandStart, int bandEnd, 
			long[] prev, long[] curr){

		// Calculate range size
		final int rangeSize = bandEnd - bandStart + 1;

		// Broadcast constants to vectors
		final LongVector vq=LongVector.broadcast(LSPECIES, q);

		// If range is at least one vector wide
		if(rangeSize>=LWIDTH){
			// Process all complete vector blocks except the final one
			for(int j=bandStart; j<bandEnd-(LWIDTH-2); j+=LWIDTH){
				// Load reference values
				LongVector refVec=LongVector.fromArray(LSPECIES, ref, j-1);

				// Match condition calculation
				VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);

				// Check for N values
				VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);

				// Combined blend operations
				LongVector scoreAddVec=vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

				// Load previous scores
				LongVector prevJ1Vec=LongVector.fromArray(LSPECIES, prev, j-1);
				LongVector prevJVec=LongVector.fromArray(LSPECIES, prev, j);

				// Calculate scores and find max
				LongVector diagScoreVec=prevJ1Vec.add(scoreAddVec);
				LongVector upScoreVec=prevJVec.add(vINS);
				LongVector maxVec=diagScoreVec.max(upScoreVec);

				// Store results
				maxVec.intoArray(curr, j);
			}

			// Process the final block with a vector, aligned to include the last element
			final int lastVectorStart=Math.max(bandStart, bandEnd+1-LWIDTH);

			// Load reference values for final vector
			LongVector refVec=LongVector.fromArray(LSPECIES, ref, lastVectorStart-1);

			// Match condition calculation
			VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);

			// Check for N values
			VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);

			// Combined blend operations
			LongVector scoreAddVec=vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

			// Load previous scores
			LongVector prevJ1Vec=LongVector.fromArray(LSPECIES, prev, lastVectorStart-1);
			LongVector prevJVec=LongVector.fromArray(LSPECIES, prev, lastVectorStart);

			// Calculate scores and find max
			LongVector diagScoreVec=prevJ1Vec.add(scoreAddVec);
			LongVector upScoreVec=prevJVec.add(vINS);
			LongVector maxVec=diagScoreVec.max(upScoreVec);

			// Store results
			maxVec.intoArray(curr, lastVectorStart);
		}else{
			// Range too small for vector processing, use scalar code
			for(int j=bandStart; j<=bandEnd; j++){
				final long r=ref[j-1];
				final boolean isMatch=((q&r)>0 && q==r);
				final boolean hasN=(q==15 || r==31);
				final long scoreAdd=isMatch?MATCH:(hasN?N_SCORE:SUB);

				final long pj1=prev[j-1], pj=prev[j];
				final long diagScore=pj1+scoreAdd;
				final long upScore=pj+INS;

				curr[j]=Math.max(diagScore, upScore);
			}
		}
	}



	/** Designed to eliminate scalar loop by reprocessing some final elements 
	 * Handles deletions 
	 * TODOL Fails on Glocal+5 with input "AAAA" "AAGAA" */
	public static void alignBandVectorDel(long q, long[] ref, int bandStart, int bandEnd, 
			long[] prev, long[] curr){
		// Calculate range size
		final int rangeSize = bandEnd - bandStart + 1;

		// Broadcast constants to vectors
		final LongVector vq=LongVector.broadcast(LSPECIES, q);
		final LongVector vMATCH=LongVector.broadcast(LSPECIES, MATCH);
		final LongVector vN_SCORE=LongVector.broadcast(LSPECIES, N_SCORE);
		final LongVector vSUB=LongVector.broadcast(LSPECIES, SUB);
		final LongVector vINS=LongVector.broadcast(LSPECIES, INS);
		final LongVector v15=LongVector.broadcast(LSPECIES, 15);

		// If range is at least one vector wide
		if(rangeSize>=LWIDTH){
			// Process all complete vector blocks except the final one
			for(int j=bandStart; j<bandEnd-(LWIDTH-2); j+=LWIDTH){
				// Load reference values
				LongVector refVec=LongVector.fromArray(LSPECIES, ref, j-1);

				// Match condition calculation
				VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);

				// Check for N values
				VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);

				// Combined blend operations
				LongVector scoreAddVec=vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

				// Load previous scores
				LongVector prevJ1Vec=LongVector.fromArray(LSPECIES, prev, j-1);
				LongVector prevJVec=LongVector.fromArray(LSPECIES, prev, j);

				// Calculate scores and find max
				LongVector diagScoreVec=prevJ1Vec.add(scoreAddVec);
				LongVector upScoreVec=prevJVec.add(vINS);
				LongVector maxVec=diagScoreVec.max(upScoreVec);

				// Store results
				maxVec.intoArray(curr, j);

				//	            if(SKIP_MASK) {
				long left=curr[j-1];
				curr[j]=left=Math.max(left+DEL, curr[j]);
				curr[j+1]=left=Math.max(left+DEL, curr[j+1]);
				curr[j+2]=left=Math.max(left+DEL, curr[j+2]);
				curr[j+3]=left=Math.max(left+DEL, curr[j+3]);
				//	            }else {
				//	            	long left=curr[j-1];
				//	            	left+=DEL;
				//	            	curr[j]=left=((left&MASK)>curr[j] ? left : curr[j]);
				//	            	left+=DEL;
				//	            	curr[j+1]=left=((left&MASK)>curr[j+1] ? left : curr[j+1]);
				//	            	left+=DEL;
				//	            	curr[j+2]=left=((left&MASK)>curr[j+2] ? left : curr[j+2]);
				//	            	left+=DEL;
				//	            	curr[j+3]=left=((left&MASK)>curr[j+3] ? left : curr[j+3]);
				//	            }
			}

			// Process the final block with a vector, aligned to include the last element
			final int j=Math.max(bandStart, bandEnd+1-LWIDTH);

			// Load reference values for final vector
			LongVector refVec=LongVector.fromArray(LSPECIES, ref, j-1);

			// Match condition calculation
			VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);

			// Check for N values
			VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);

			// Combined blend operations
			LongVector scoreAddVec=vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

			// Load previous scores
			LongVector prevJ1Vec=LongVector.fromArray(LSPECIES, prev, j-1);
			LongVector prevJVec=LongVector.fromArray(LSPECIES, prev, j);

			// Calculate scores and find max
			LongVector diagScoreVec=prevJ1Vec.add(scoreAddVec);
			LongVector upScoreVec=prevJVec.add(vINS);
			LongVector maxVec=diagScoreVec.max(upScoreVec);

			// Store results
			maxVec.intoArray(curr, j);

			//	        print(curr, "?");

			//            if(SKIP_MASK) {
			long left=curr[j-1];
			curr[j]=left=Math.max(left+DEL, curr[j]);
			curr[j+1]=left=Math.max(left+DEL, curr[j+1]);
			curr[j+2]=left=Math.max(left+DEL, curr[j+2]);
			curr[j+3]=left=Math.max(left+DEL, curr[j+3]);
			//            }else {
			//            	long left=curr[j-1];
			//            	left+=DEL;
			//            	curr[j]=left=((left&MASK)>curr[j] ? left : curr[j]);
			//            	left+=DEL;
			//            	curr[j+1]=left=((left&MASK)>curr[j+1] ? left : curr[j+1]);
			//            	left+=DEL;
			//            	curr[j+2]=left=((left&MASK)>curr[j+2] ? left : curr[j+2]);
			//            	left+=DEL;
			//            	curr[j+3]=left=((left&MASK)>curr[j+3] ? left : curr[j+3]);
			//            }
		}else{
			//	    	System.err.println("Scalar");
			// Range too small for vector processing, use scalar code
			for(int j=bandStart; j<=bandEnd; j++){
				final long r=ref[j-1];
				final boolean isMatch=((q&r)>0 && q==r);
				final boolean hasN=(q==15 || r==31);
				final long scoreAdd=isMatch?MATCH:(hasN?N_SCORE:SUB);

				final long pj1 = prev[j-1], pj = prev[j], cj1=curr[j-1];
				final long diagScore = pj1 + scoreAdd;
				final long upScore = pj + INS;
				final long leftScore = cj1+DEL;

				long score = Math.max(leftScore, Math.max(diagScore, upScore));
				curr[j] = score;
			}
		}
	}

	/** 
	 * Vector implementation that also tracks and returns the position of maximum score
	 * @return Position of the maximum score in the processed band
	 */
	public static int alignBandVectorAndReturnMaxPos(long q, long[] ref, int bandStart, int bandEnd, 
			long[] prev, long[] curr, 
			long MATCH, long N_SCORE, long SUB, long INS, long SCORE_MASK) {

		// Calculate range size
		final int rangeSize = bandEnd - bandStart + 1;

		// Initialize tracking vectors for max score and position
		LongVector maxScoreVec = LongVector.broadcast(LSPECIES, Long.MIN_VALUE);
		LongVector maxPosVec = LongVector.broadcast(LSPECIES, -1);

		// Create initial position vector (0,1,2,3...)
		long[] posInit = new long[LWIDTH];
		for(int i=0; i<LWIDTH; i++) posInit[i] = i;
		LongVector posTemplate = LongVector.fromArray(LSPECIES, posInit, 0);

		// Broadcast constants to vectors
		final LongVector vq = LongVector.broadcast(LSPECIES, q);
		final LongVector vMATCH = LongVector.broadcast(LSPECIES, MATCH);
		final LongVector vN_SCORE = LongVector.broadcast(LSPECIES, N_SCORE);
		final LongVector vSUB = LongVector.broadcast(LSPECIES, SUB);
		final LongVector vINS = LongVector.broadcast(LSPECIES, INS);
		final LongVector v15 = LongVector.broadcast(LSPECIES, 15);
		final LongVector vSCORE_MASK = LongVector.broadcast(LSPECIES, SCORE_MASK);

		// If range is at least one vector wide
		if(rangeSize >= LWIDTH) {
			// Process all complete vector blocks except the final one
			for(int j=bandStart; j<bandEnd-(LWIDTH-2); j+=LWIDTH) {
				// Load reference values
				LongVector refVec = LongVector.fromArray(LSPECIES, ref, j-1);

				// Match condition calculation
				VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);

				// Check for N values
				VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);

				// Combined blend operations
				LongVector scoreAddVec = vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

				// Load previous scores
				LongVector prevJ1Vec = LongVector.fromArray(LSPECIES, prev, j-1);
				LongVector prevJVec = LongVector.fromArray(LSPECIES, prev, j);

				// Calculate scores and find max
				LongVector diagScoreVec = prevJ1Vec.add(scoreAddVec);
				LongVector upScoreVec = prevJVec.add(vINS);
				LongVector maxVec = diagScoreVec.max(upScoreVec);

				// Store results
				maxVec.intoArray(curr, j);

				// Extract score part for comparison (apply mask)
				LongVector scorePart = maxVec.and(vSCORE_MASK);

				// Calculate position vector for this chunk (j, j+1, j+2...)
				LongVector posVec = posTemplate.add(j);

				// Find elements where current scores > max scores
				VectorMask<Long> betterMask = scorePart.compare(VectorOperators.GT, maxScoreVec);

				// Update tracking vectors where better scores found
				maxScoreVec = maxScoreVec.blend(scorePart, betterMask);
				maxPosVec = maxPosVec.blend(posVec, betterMask);
			}

			// Process the final block with a vector, aligned to include the last element
			final int lastVectorStart = Math.max(bandStart, bandEnd+1-LWIDTH);

			// Load reference values for final vector
			LongVector refVec = LongVector.fromArray(LSPECIES, ref, lastVectorStart-1);

			// Match condition calculation
			VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);

			// Check for N values
			VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);

			// Combined blend operations
			LongVector scoreAddVec = vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

			// Load previous scores
			LongVector prevJ1Vec = LongVector.fromArray(LSPECIES, prev, lastVectorStart-1);
			LongVector prevJVec = LongVector.fromArray(LSPECIES, prev, lastVectorStart);

			// Calculate scores and find max
			LongVector diagScoreVec = prevJ1Vec.add(scoreAddVec);
			LongVector upScoreVec = prevJVec.add(vINS);
			LongVector maxVec = diagScoreVec.max(upScoreVec);

			// Store results
			maxVec.intoArray(curr, lastVectorStart);

			// Extract score part for comparison (apply mask)
			LongVector scorePart = maxVec.and(vSCORE_MASK);

			// Calculate position vector for this chunk
			LongVector posVec = posTemplate.add(lastVectorStart);

			// Find elements where current scores > max scores
			VectorMask<Long> betterMask = scorePart.compare(VectorOperators.GT, maxScoreVec);

			// Update tracking vectors where better scores found
			maxScoreVec = maxScoreVec.blend(scorePart, betterMask);
			maxPosVec = maxPosVec.blend(posVec, betterMask);
		} else {
			// Range too small for vector processing, use scalar code
			long maxScore = Long.MIN_VALUE;
			int maxPos = -1;

			for(int j=bandStart; j<=bandEnd; j++) {
				final long r = ref[j-1];
				final boolean isMatch = ((q&r)>0 && q==r);
				final boolean hasN = (q==15 || r==31);
				final long scoreAdd = isMatch ? MATCH : (hasN ? N_SCORE : SUB);

				final long pj1 = prev[j-1], pj = prev[j];
				final long diagScore = pj1 + scoreAdd;
				final long upScore = pj + INS;

				long score = Math.max(diagScore, upScore);
				curr[j] = score;

				// Track maximum score
				final boolean better=((score&SCORE_MASK)>maxScore);
				maxScore=better ? score : maxScore;
				maxPos=better ? j : maxPos;
			}

			return maxPos;
		}

		// Now find position of max score in the tracking vectors
		long[] maxScores = new long[LWIDTH];
		long[] maxPositions = new long[LWIDTH];
		maxScoreVec.intoArray(maxScores, 0);
		maxPosVec.intoArray(maxPositions, 0);

		// Find max element in reduced vector
		long globalMax = Long.MIN_VALUE;
		int maxPos = -1;

		for(int i=0; i<LWIDTH; i++) {
			if(maxScores[i] > globalMax) {
				globalMax = maxScores[i];
				maxPos = (int)maxPositions[i];
			}
		}

		return maxPos;
	}



	/** 
	 * Vector implementation that also tracks and returns the position of maximum score
	 * @return Position of the maximum score in the processed band
	 */
	public static int alignBandVectorAndReturnMaxPos(long q, long[] ref, int bandStart, int bandEnd, 
			long[] prev, long[] curr) {

		// Calculate range size
		final int rangeSize = bandEnd - bandStart + 1;

		// Initialize tracking vectors for max score and position
		LongVector maxScoreVec = LongVector.broadcast(LSPECIES, Long.MIN_VALUE);
		LongVector maxPosVec = LongVector.broadcast(LSPECIES, -1);

		// Create initial position vector (0,1,2,3...)
		long[] posInit = new long[LWIDTH];
		for(int i=0; i<LWIDTH; i++) posInit[i] = i;
		LongVector posTemplate = LongVector.fromArray(LSPECIES, posInit, 0);

		// Broadcast constants to vectors
		final LongVector vq = LongVector.broadcast(LSPECIES, q);
		final LongVector vMATCH = LongVector.broadcast(LSPECIES, MATCH);
		final LongVector vN_SCORE = LongVector.broadcast(LSPECIES, N_SCORE);
		final LongVector vSUB = LongVector.broadcast(LSPECIES, SUB);
		final LongVector vINS = vSUB;
		final LongVector v15 = LongVector.broadcast(LSPECIES, 15);

		// If range is at least one vector wide
		if(rangeSize >= LWIDTH) {
			// Process all complete vector blocks except the final one
			for(int j=bandStart; j<bandEnd-(LWIDTH-2); j+=LWIDTH) {
				// Load reference values
				LongVector refVec = LongVector.fromArray(LSPECIES, ref, j-1);

				// Match condition calculation
				VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);

				// Check for N values
				VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);

				// Combined blend operations
				LongVector scoreAddVec = vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

				// Load previous scores
				LongVector prevJ1Vec = LongVector.fromArray(LSPECIES, prev, j-1);
				LongVector prevJVec = LongVector.fromArray(LSPECIES, prev, j);

				// Calculate scores and find max
				LongVector diagScoreVec = prevJ1Vec.add(scoreAddVec);
				LongVector upScoreVec = prevJVec.add(vINS);
				LongVector maxVec = diagScoreVec.max(upScoreVec);

				// Store results
				maxVec.intoArray(curr, j);

				// Extract score part for comparison (apply mask)
				//	            LongVector scorePart = maxVec.and(vSCORE_MASK);

				// Calculate position vector for this chunk (j, j+1, j+2...)
				LongVector posVec = posTemplate.add(j);

				// Find elements where current scores > max scores
				VectorMask<Long> betterMask = maxVec.compare(VectorOperators.GT, maxScoreVec);

				// Update tracking vectors where better scores found
				maxScoreVec = maxScoreVec.blend(maxVec, betterMask);
				maxPosVec = maxPosVec.blend(posVec, betterMask);
			}

			// Process the final block with a vector, aligned to include the last element
			final int lastVectorStart = Math.max(bandStart, bandEnd+1-LWIDTH);

			// Load reference values for final vector
			LongVector refVec = LongVector.fromArray(LSPECIES, ref, lastVectorStart-1);

			// Match condition calculation
			VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);

			// Check for N values
			VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);

			// Combined blend operations
			LongVector scoreAddVec = vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

			// Load previous scores
			LongVector prevJ1Vec = LongVector.fromArray(LSPECIES, prev, lastVectorStart-1);
			LongVector prevJVec = LongVector.fromArray(LSPECIES, prev, lastVectorStart);

			// Calculate scores and find max
			LongVector diagScoreVec = prevJ1Vec.add(scoreAddVec);
			LongVector upScoreVec = prevJVec.add(vINS);
			LongVector maxVec = diagScoreVec.max(upScoreVec);

			// Store results
			maxVec.intoArray(curr, lastVectorStart);

			// Extract score part for comparison (apply mask)
			//	        LongVector scorePart = maxVec.and(vSCORE_MASK);

			// Calculate position vector for this chunk
			LongVector posVec = posTemplate.add(lastVectorStart);

			// Find elements where current scores > max scores
			VectorMask<Long> betterMask = maxVec.compare(VectorOperators.GT, maxScoreVec);

			// Update tracking vectors where better scores found
			maxScoreVec = maxScoreVec.blend(maxVec, betterMask);
			maxPosVec = maxPosVec.blend(posVec, betterMask);
		} else {
			// Range too small for vector processing, use scalar code
			long maxScore = Long.MIN_VALUE;
			int maxPos = -1;

			for(int j=bandStart; j<=bandEnd; j++) {
				final long r = ref[j-1];
				final boolean isMatch = ((q&r)>0 && q==r);
				final boolean hasN = (q==15 || r==31);
				final long scoreAdd = isMatch ? MATCH : (hasN ? N_SCORE : SUB);

				final long pj1 = prev[j-1], pj = prev[j];
				final long diagScore = pj1 + scoreAdd;
				final long upScore = pj + INS;

				long score = Math.max(diagScore, upScore);
				curr[j] = score;

				// Track maximum score
				// Track best score in row
				final boolean better=(score>maxScore);
				maxScore=better ? score : maxScore;
				maxPos=better ? j : maxPos;
			}

			return maxPos;
		}

		if(true) {// Vector version of finding max position in final vectors
			// Similar speed to scalar for AVX-256
			long maxValue=maxScoreVec.reduceLanes(VectorOperators.MAX);
			// Create a mask for positions equal to the max
			VectorMask<Long> maxMask=maxScoreVec.compare(VectorOperators.EQ, maxValue);
			// Find the first position where mask is true (more efficient than looping)
			int firstMaxPos=maxMask.firstTrue();
			// Get the corresponding position from the position vector
			return (int)maxPosVec.lane(firstMaxPos);
		} else {// Scalar version of finding max position in final vectors
			long globalMax=Long.MIN_VALUE;
			int maxPos=-1;
			for(int i=0; i<LWIDTH; i++) {
				long score=maxScoreVec.lane(i);
				if(score>globalMax) {
					globalMax=score;
					maxPos=(int)maxPosVec.lane(i);
				}
			}
			return maxPos;
		}
	}



	/** 
	 * Vector implementation that also tracks and returns the position of maximum score
	 * Also handles deletions
	 * @return Position of the maximum score in the processed band
	 */
	public static int alignBandVectorAndReturnMaxPos2(long q, long[] ref, int bandStart, int bandEnd, 
			long[] prev, long[] curr) {
		// Calculate range size
		final int rangeSize = bandEnd - bandStart + 1;

		// Initialize tracking vectors for max score and position
		LongVector maxScoreVec = LongVector.broadcast(LSPECIES, Long.MIN_VALUE);
		LongVector maxPosVec = LongVector.broadcast(LSPECIES, -1);

		// Create initial position vector (0,1,2,3...)
		long[] posInit = new long[LWIDTH];
		for(int i=0; i<LWIDTH; i++) posInit[i] = i;
		LongVector posTemplate = LongVector.fromArray(LSPECIES, posInit, 0);

		// Broadcast constants to vectors
		final LongVector vq = LongVector.broadcast(LSPECIES, q);
		final LongVector vMATCH = LongVector.broadcast(LSPECIES, MATCH);
		final LongVector vN_SCORE = LongVector.broadcast(LSPECIES, N_SCORE);
		final LongVector vSUB = LongVector.broadcast(LSPECIES, SUB);
		final LongVector vINS = vSUB;
		final LongVector v15 = LongVector.broadcast(LSPECIES, 15);

		// If range is at least one vector wide
		if(rangeSize >= LWIDTH) {
			// Process all complete vector blocks except the final one
			for(int j=bandStart; j<bandEnd-(LWIDTH-2); j+=LWIDTH) {
				// Load reference values
				LongVector refVec = LongVector.fromArray(LSPECIES, ref, j-1);

				// Match condition calculation
				VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);

				// Check for N values
				VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);

				// Combined blend operations
				LongVector scoreAddVec = vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

				// Load previous scores
				LongVector prevJ1Vec = LongVector.fromArray(LSPECIES, prev, j-1);
				LongVector prevJVec = LongVector.fromArray(LSPECIES, prev, j);

				// Calculate scores and find max
				LongVector diagScoreVec = prevJ1Vec.add(scoreAddVec);
				LongVector upScoreVec = prevJVec.add(vINS);
				LongVector maxVec = diagScoreVec.max(upScoreVec);

				// Store results
				maxVec.intoArray(curr, j);

				// Extract score part for comparison (apply mask)
				//	            LongVector scorePart = maxVec.and(vSCORE_MASK);

				// Calculate position vector for this chunk (j, j+1, j+2...)
				LongVector posVec = posTemplate.add(j);

				// Find elements where current scores > max scores
				VectorMask<Long> betterMask = maxVec.compare(VectorOperators.GT, maxScoreVec);

				// Update tracking vectors where better scores found
				maxScoreVec = maxScoreVec.blend(maxVec, betterMask);
				maxPosVec = maxPosVec.blend(posVec, betterMask);

				//	            if(SKIP_MASK) {
				long left=curr[j-1];
				curr[j]=left=Math.max(left+DEL, curr[j]);
				curr[j+1]=left=Math.max(left+DEL, curr[j+1]);
				curr[j+2]=left=Math.max(left+DEL, curr[j+2]);
				curr[j+3]=left=Math.max(left+DEL, curr[j+3]);
				//	            }else {
				//	            	long left=curr[j-1];
				//	            	left+=DEL;
				//	            	curr[j]=left=((left&MASK)>curr[j] ? left : curr[j]);
				//	            	left+=DEL;
				//	            	curr[j+1]=left=((left&MASK)>curr[j+1] ? left : curr[j+1]);
				//	            	left+=DEL;
				//	            	curr[j+2]=left=((left&MASK)>curr[j+2] ? left : curr[j+2]);
				//	            	left+=DEL;
				//	            	curr[j+3]=left=((left&MASK)>curr[j+3] ? left : curr[j+3]);
				//	            }
			}

			// Process the final block with a vector, aligned to include the last element
			final int j = Math.max(bandStart, bandEnd+1-LWIDTH);

			// Load reference values for final vector
			LongVector refVec = LongVector.fromArray(LSPECIES, ref, j-1);

			// Match condition calculation
			VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);

			// Check for N values
			VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);

			// Combined blend operations
			LongVector scoreAddVec = vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);

			// Load previous scores
			LongVector prevJ1Vec = LongVector.fromArray(LSPECIES, prev, j-1);
			LongVector prevJVec = LongVector.fromArray(LSPECIES, prev, j);

			// Calculate scores and find max
			LongVector diagScoreVec = prevJ1Vec.add(scoreAddVec);
			LongVector upScoreVec = prevJVec.add(vINS);
			LongVector maxVec = diagScoreVec.max(upScoreVec);

			// Store results
			maxVec.intoArray(curr, j);

			// Extract score part for comparison (apply mask)
			//	        LongVector scorePart = maxVec.and(vSCORE_MASK);

			// Calculate position vector for this chunk
			LongVector posVec = posTemplate.add(j);

			// Find elements where current scores > max scores
			VectorMask<Long> betterMask = maxVec.compare(VectorOperators.GT, maxScoreVec);

			// Update tracking vectors where better scores found
			maxScoreVec = maxScoreVec.blend(maxVec, betterMask);
			maxPosVec = maxPosVec.blend(posVec, betterMask);

			//          if(SKIP_MASK) {
			long left=curr[j-1];
			curr[j]=left=Math.max(left+DEL, curr[j]);
			curr[j+1]=left=Math.max(left+DEL, curr[j+1]);
			curr[j+2]=left=Math.max(left+DEL, curr[j+2]);
			curr[j+3]=left=Math.max(left+DEL, curr[j+3]);
			//          }else {
			//          	long left=curr[j-1];
			//          	left+=DEL;
			//          	curr[j]=left=((left&MASK)>curr[j] ? left : curr[j]);
			//          	left+=DEL;
			//          	curr[j+1]=left=((left&MASK)>curr[j+1] ? left : curr[j+1]);
			//          	left+=DEL;
			//          	curr[j+2]=left=((left&MASK)>curr[j+2] ? left : curr[j+2]);
			//          	left+=DEL;
			//          	curr[j+3]=left=((left&MASK)>curr[j+3] ? left : curr[j+3]);
			//          }
		} else {
			// Range too small for vector processing, use scalar code
			long maxScore = Long.MIN_VALUE;
			int maxPos = -1;

			for(int j=bandStart; j<=bandEnd; j++) {
				final long r = ref[j-1];
				final boolean isMatch = ((q&r)>0 && q==r);
				final boolean hasN = (q==15 || r==31);
				final long scoreAdd = isMatch ? MATCH : (hasN ? N_SCORE : SUB);

				final long pj1 = prev[j-1], pj = prev[j], cj1=curr[j-1];
				final long diagScore = pj1 + scoreAdd;
				final long upScore = pj + INS;
				final long leftScore = cj1+DEL;

				long score = Math.max(leftScore, Math.max(diagScore, upScore));
				curr[j] = score;

				// Track maximum score
				// Track best score in row
				final boolean better=(score>maxScore);
				maxScore=better ? score : maxScore;
				maxPos=better ? j : maxPos;
			}

			return maxPos;
		}

		if(true) {// Vector version of finding max position in final vectors
			// Similar speed to scalar for AVX-256
			long maxValue=maxScoreVec.reduceLanes(VectorOperators.MAX);
			// Create a mask for positions equal to the max
			VectorMask<Long> maxMask=maxScoreVec.compare(VectorOperators.EQ, maxValue);
			// Find the first position where mask is true (more efficient than looping)
			int firstMaxPos=maxMask.firstTrue();
			// Get the corresponding position from the position vector
			return (int)maxPosVec.lane(firstMaxPos);
		} else {// Scalar version of finding max position in final vectors
			long globalMax=Long.MIN_VALUE;
			int maxPos=-1;
			for(int i=0; i<LWIDTH; i++) {
				long score=maxScoreVec.lane(i);
				if(score>globalMax) {
					globalMax=score;
					maxPos=(int)maxPosVec.lane(i);
				}
			}
			return maxPos;
		}
	}



	/** Vector implementation of inner loop that tracks the maximum score.
	 * Max score is independent of the deletion tail loop since it only reduces scores.
	 * @return Position of the maximum score in the processed band */
	public static int alignBandVectorAndReturnMaxPosConcise(
			long q, long[] ref, int bandStart, int bandEnd, long[] prev, long[] curr, 
			long MATCH, long N_SCORE, long SUB, long INS) {
		final int rangeSize=bandEnd-bandStart+1;// Calculate range size

		if(rangeSize<LWIDTH) {
			return alignBandVectorAndReturnMaxPosScalar(
					q, ref, bandStart, bandEnd, prev, curr, MATCH, N_SCORE, SUB, INS);}
		// Initialize tracking vectors for max score and position
		LongVector maxScoreVec=LongVector.broadcast(LSPECIES, Long.MIN_VALUE);
		LongVector maxPosVec=LongVector.broadcast(LSPECIES, -1);
		long[] posInit=new long[LWIDTH];// Initialize position vector to (0,1,2,3...)
		for(int i=0; i<LWIDTH; i++) {posInit[i]=i;}
		LongVector posTemplate=LongVector.fromArray(LSPECIES, posInit, 0);
		// Broadcast constants to vectors
		final LongVector vq=LongVector.broadcast(LSPECIES, q);
		final LongVector vMATCH=LongVector.broadcast(LSPECIES, MATCH);
		final LongVector vN_SCORE=LongVector.broadcast(LSPECIES, N_SCORE);
		final LongVector vSUB=LongVector.broadcast(LSPECIES, SUB);
		final LongVector vINS=LongVector.broadcast(LSPECIES, INS);
		final LongVector v15=LongVector.broadcast(LSPECIES, 15);
		if(rangeSize>=LWIDTH) {// If range is at least one vector wide
			// Process all complete vector blocks except the final one
			for(int j=bandStart; j<bandEnd-(LWIDTH-2); j+=LWIDTH) {
				LongVector refVec=LongVector.fromArray(LSPECIES, ref, j-1);// Load reference
				// Match condition calculation
				VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);
				// Check for N values
				VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);
				// Combined blend operations
				LongVector scoreAddVec=vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);
				// Load previous scores
				LongVector prevJ1Vec=LongVector.fromArray(LSPECIES, prev, j-1);
				LongVector prevJVec=LongVector.fromArray(LSPECIES, prev, j);
				// Calculate scores and find max
				LongVector diagScoreVec=prevJ1Vec.add(scoreAddVec);
				LongVector upScoreVec=prevJVec.add(vINS);
				LongVector maxVec=diagScoreVec.max(upScoreVec);
				maxVec.intoArray(curr, j);// Store results
				// Calculate position vector for this chunk (j, j+1, j+2...)
				LongVector posVec=posTemplate.add(j);
				// Find elements where current scores > max scores
				VectorMask<Long> betterMask=maxVec.compare(VectorOperators.GT, maxScoreVec);
				// Update tracking vectors where better scores found
				maxScoreVec=maxScoreVec.blend(maxVec, betterMask);
				maxPosVec=maxPosVec.blend(posVec, betterMask);
			}
			// Process the final block with a vector, aligned to include the last element
			// It is safe to reprocess some elements because the answer is unchanged
			final int lastVectorStart=Math.max(bandStart, bandEnd+1-LWIDTH);
			// Load reference values for final vector
			LongVector refVec=LongVector.fromArray(LSPECIES, ref, lastVectorStart-1);
			// Match condition calculation
			VectorMask<Long> isMatchMask=vq.compare(VectorOperators.EQ, refVec);
			// Check for N values
			VectorMask<Long> hasNMask=vq.or(refVec).compare(VectorOperators.GE, v15);
			// Combined blend operations
			LongVector scoreAddVec=vSUB.blend(vMATCH, isMatchMask).blend(vN_SCORE, hasNMask);
			// Load previous scores
			LongVector prevJ1Vec=LongVector.fromArray(LSPECIES, prev, lastVectorStart-1);
			LongVector prevJVec=LongVector.fromArray(LSPECIES, prev, lastVectorStart);
			// Calculate scores and find max
			LongVector diagScoreVec=prevJ1Vec.add(scoreAddVec);
			LongVector upScoreVec=prevJVec.add(vINS);
			LongVector maxVec=diagScoreVec.max(upScoreVec);
			maxVec.intoArray(curr, lastVectorStart);// Store results
			LongVector posVec=posTemplate.add(lastVectorStart);// Calculate position vector
			// Find elements where current scores > max scores
			VectorMask<Long> betterMask=maxVec.compare(VectorOperators.GT, maxScoreVec);
			// Update tracking vectors where better scores found
			maxScoreVec=maxScoreVec.blend(maxVec, betterMask);
			maxPosVec=maxPosVec.blend(posVec, betterMask);
		} else {
			// Range too small for vector processing, use scalar code
			long maxScore=Long.MIN_VALUE;
			int maxPos=-1;
			for(int j=bandStart; j<=bandEnd; j++) {
				final long r=ref[j-1];
				final boolean isMatch=((q&r)>0 && q==r);
				final boolean hasN=(q==15 || r==31);
				final long scoreAdd=isMatch?MATCH:(hasN?N_SCORE:SUB);
				final long pj1=prev[j-1], pj=prev[j];
				final long diagScore=pj1+scoreAdd;
				final long upScore=pj+INS;
				long score=Math.max(diagScore, upScore);
				curr[j]=score;
				final boolean better=(score>maxScore);// Track maximum score
				maxScore=better?score:maxScore;
				maxPos=better?j:maxPos;
			}
			return maxPos;
		}
		// Vector version of finding max position in final vectors
		// Similar speed to scalar for AVX-256
		long maxValue=maxScoreVec.reduceLanes(VectorOperators.MAX);
		// Create a mask for positions equal to the max
		VectorMask<Long> maxMask=maxScoreVec.compare(VectorOperators.EQ, maxValue);
		// Find the first position where mask is true (more efficient than looping)
		int firstMaxPos=maxMask.firstTrue();
		// Get the corresponding position from the position vector
		return (int)maxPosVec.lane(firstMaxPos);
	}

	/** Scalar implementation of inner loop for band shorter than LWIDTH.
	 * @return Position of the maximum score in the processed band */
	public static int alignBandVectorAndReturnMaxPosScalar(
			long q, long[] ref, int bandStart, int bandEnd, long[] prev, long[] curr, 
			long MATCH, long N_SCORE, long SUB, long INS) {
		final int rangeSize=bandEnd-bandStart+1;// Calculate range size
		assert(rangeSize<LWIDTH);
		long maxScore=Long.MIN_VALUE;
		int maxPos=-1;
		for(int j=bandStart; j<=bandEnd; j++) {
			final long r=ref[j-1];
			final boolean isMatch=((q&r)>0 && q==r);
			final boolean hasN=(q==15 || r==31);
			final long scoreAdd=isMatch?MATCH:(hasN?N_SCORE:SUB);
			final long pj1=prev[j-1], pj=prev[j];
			final long diagScore=pj1+scoreAdd;
			final long upScore=pj+INS;
			long score=Math.max(diagScore, upScore);
			curr[j]=score;
			final boolean better=(score>maxScore);// Track maximum score
			maxScore=better?score:maxScore;
			maxPos=better?j:maxPos;
		}
		return maxPos;
	}

	public static void alignBandVectorInt(int q, int[] ref, int bandStart, int bandEnd,
			int[] prev, int[] curr,
			int MATCH, int N_SCORE, int SUB, int INS){

		// Calculate vector-aligned bounds
		int alignedStart=bandStart;
		int alignedEnd=bandStart+((bandEnd-bandStart+1)/IWIDTH)*IWIDTH;

		// Broadcast constants to vectors
		IntVector vq=IntVector.broadcast(ISPECIES, q);
		IntVector vMATCH=IntVector.broadcast(ISPECIES, MATCH);
		IntVector vN_SCORE=IntVector.broadcast(ISPECIES, N_SCORE);
		IntVector vSUB=IntVector.broadcast(ISPECIES, SUB);
		IntVector vINS=IntVector.broadcast(ISPECIES, INS);

		// Process in IWIDTH chunks
		for(int j=alignedStart; j<alignedEnd; j+=IWIDTH){
			// Load reference values
			IntVector refVec=IntVector.fromArray(ISPECIES, ref, j-1);

			// Calculate match conditions with 4-bit encoding
			IntVector andResult=vq.lanewise(VectorOperators.AND, refVec);
			VectorMask<Integer> isCompatible=andResult.compare(VectorOperators.GT, 0);
			VectorMask<Integer> isEqual=vq.compare(VectorOperators.EQ, refVec);
			VectorMask<Integer> isMatchMask=isEqual.and(isCompatible);

			// Check for N values
			VectorMask<Integer> qIsN=vq.compare(VectorOperators.EQ, 15);
			VectorMask<Integer> refIsN=refVec.compare(VectorOperators.EQ, 31);
			VectorMask<Integer> hasNMask=qIsN.or(refIsN);

			// Start with SUB score for all positions
			IntVector scoreAddVec=vSUB;

			// Apply MATCH score where matches occur (this overwrites SUB scores)
			scoreAddVec=scoreAddVec.blend(vMATCH, isMatchMask);

			// Apply N_SCORE where Ns occur (this overwrites any previous scores)
			scoreAddVec=scoreAddVec.blend(vN_SCORE, hasNMask);

			// Load previous scores
			IntVector prevJ1Vec=IntVector.fromArray(ISPECIES, prev, j-1);
			IntVector prevJVec=IntVector.fromArray(ISPECIES, prev, j);

			// Calculate diagonal and up scores
			IntVector diagScoreVec=prevJ1Vec.add(scoreAddVec);
			IntVector upScoreVec=prevJVec.add(vINS);

			// Find max score between diagonal and up
			IntVector maxVec=diagScoreVec.max(upScoreVec);

			// Store results
			maxVec.intoArray(curr, j);
		}

		// Handle remaining elements (scalar code)
		for(int j=alignedEnd; j<=bandEnd; j++){
			final int r=ref[j-1];
			final boolean isMatch=((q&r)>0 && q==r);
			final boolean hasN=(q==15 || r==31);
			final int scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

			final int pj1=prev[j-1], pj=prev[j];
			final int diagScore=pj1+scoreAdd;
			final int upScore=pj+INS;

			curr[j]=Math.max(diagScore, upScore);
		}
	}

	/**
	 * Pure vector implementation of deletion tail loop based on exact pseudocode
	 * TODO: Slow
	 */
	public static void processDeletionsTailVectorPure(
			long[] curr, int bandStart, int bandEnd) {
		//	    print(curr, "Initial");

		// Calculate vector-aligned bounds
		final int limit = bandStart + ((bandEnd - bandStart + 1)/LWIDTH)*LWIDTH;

		// Broadcast constant
		final LongVector delIncrementVec = LongVector.broadcast(LSPECIES, DEL);

		// Process full LWIDTH chunks with exactly two passes as specified

		//	    if(limit>bandStart) {System.err.println("Vector processing "+bandStart+" to "+(limit-1));}
		for(int i=bandStart; i<limit; i+=LWIDTH) {
			for(int pass=0, dist=LWIDTH/2; dist>0 && i-dist>=0; pass++, dist/=2) {
				int k = i-dist;
				// Load cells from distance 'dist' away
				LongVector a = LongVector.fromArray(LSPECIES, curr, k);
				LongVector b = LongVector.fromArray(LSPECIES, curr, i);
				// Adjust DEL_INCREMENT for the distance
				LongVector adjustedDelIncrement = LongVector.broadcast(LSPECIES, DEL * dist);
				LongVector aPlus = a.add(adjustedDelIncrement);
				b = b.max(aPlus);
				b.intoArray(curr, i);
			}

		}

		// Handle remaining elements with scalar code
		long leftCell = curr[limit-1];
		//	    if(bandEnd>=limit) {System.err.println("Scalar processing "+limit+" to "+bandEnd);}
		for(int j=limit; j<=bandEnd; j++) {
			final long maxDiagUp = curr[j];
			final long leftScore = leftCell + DEL;
			leftCell = (maxDiagUp >= leftScore) ? maxDiagUp : leftScore;
			curr[j] = leftCell;
		}
		//	    print(curr, "Final");
		//	    System.err.println();
	}

	/**
	 * Pure vector implementation of deletion tail loop based on exact pseudocode
	 * Inner loop unrolled
	 */
	public static void processDeletionsTailVectorUnrolled(
			long[] curr, int bandStart, int bandEnd, long DEL_INCREMENT) {
		//	    print(curr, "Initial");

		// Calculate vector-aligned bounds
		final int limit = bandStart + ((bandEnd - bandStart + 1)/LWIDTH)*LWIDTH;

		// Broadcast constant
		final LongVector delIncrementVec1 = LongVector.broadcast(LSPECIES, DEL_INCREMENT);
		final LongVector delIncrementVec2 = LongVector.broadcast(LSPECIES, DEL_INCREMENT*2);

		// Process full LWIDTH chunks with exactly two passes as specified

		//	    if(limit>bandStart) {System.err.println("Vector processing "+bandStart+" to "+(limit-1));}
		for(int i=bandStart; i<limit; i+=LWIDTH) {
			//Dist=2
			int k = Math.max(0, i-2);
			// Load cells from distance 'dist' away
			LongVector a = LongVector.fromArray(LSPECIES, curr, k);
			LongVector b = LongVector.fromArray(LSPECIES, curr, i);
			// Adjust DEL_INCREMENT for the distance
			LongVector aPlus = a.add(delIncrementVec2);
			b = b.max(aPlus);
			b.intoArray(curr, i);

			//Dist = 1
			k = Math.max(0, i-1);
			// Load cells from distance 'dist' away
			a = LongVector.fromArray(LSPECIES, curr, k);
			b = LongVector.fromArray(LSPECIES, curr, i);
			// Adjust DEL_INCREMENT for the distance
			aPlus = a.add(delIncrementVec1);
			b = b.max(aPlus);
			b.intoArray(curr, i);
		}

		// Handle remaining elements with scalar code
		long leftCell = curr[limit-1];
		//	    if(bandEnd>=limit) {System.err.println("Scalar processing "+limit+" to "+bandEnd);}
		for(int j=limit; j<=bandEnd; j++) {
			final long maxDiagUp = curr[j];
			final long leftScore = leftCell + DEL_INCREMENT;
			leftCell = (maxDiagUp >= leftScore) ? maxDiagUp : leftScore;
			curr[j] = leftCell;
		}
		//	    print(curr, "Final");
		//	    System.err.println();
	}

	private static final void print(long[] curr, String name) {
		System.err.print(name+" Score:\t");
		for(int i=0; i<curr.length; i++) {System.err.print((curr[i]>>42)+" ");}
		System.err.print("\n"+name+" Dels: \t");
		for(int i=0; i<curr.length; i++) {System.err.print((((curr[i]>>21)&0xFFFF)+" "));}
		System.err.println();
	}

	public static void processInnerDiagonalVectorized(
			int innerMinI, int innerMaxI,
			long[] revQuery, long[] ref, int jm1Start,
			long[] diag_km2, long[] diag_km1, long[] diag_k,
			int km2DiagOffset, int km1DiagOffset,
			int qLen, long MATCH, long N_SCORE, long SUB,
			long INS, long DEL_INCREMENT, long SCORE_MASK) {

		// For now, just implement a scalar version that exactly matches the original
		for(int i=innerMinI, offset=0, jm1=jm1Start; i<=innerMaxI; i++, offset++, jm1++) {
			long q = revQuery[i];
			long r = ref[jm1];

			boolean isMatch = (q == r);
			boolean hasN = (q | r) > 15;
			long scoreAdd = isMatch ? MATCH : (hasN ? N_SCORE : SUB);

			// Get adjacent cell values using incrementing indices
			int upDiagIdx = offset + km2DiagOffset;
			int leftIdx = offset + km1DiagOffset - 1;

			long diagValue = diag_km2[upDiagIdx];
			long upValue = diag_km1[upDiagIdx];

			// Handle first column special case
			long leftValue;
			if (leftIdx < 0) {
				int currentRow = qLen - i;
				leftValue = currentRow * INS;
			} else {
				leftValue = diag_km1[leftIdx];
			}

			// Calculate scores directly
			long diagScore = diagValue + scoreAdd;
			long upScore = upValue + INS;
			long leftScore = leftValue + DEL_INCREMENT;

			// Find max using conditional expressions - exactly as original
			final long maxDiagUp = Math.max(diagScore, upScore);
			final long maxValue = ((maxDiagUp & SCORE_MASK) >= leftScore) ? maxDiagUp : leftScore;

			diag_k[offset] = maxValue;

			if(debug) {
				int matrixRow = qLen - i;
				System.err.println("Scalar Setting (" + matrixRow + "," + (jm1+1) + ") to " + maxValue);
			}
		}
	}

	@SuppressWarnings("restriction")
	public static void processCrossCutDiagonalSIMD(
			long[] revQuery, long[] ref, int k, 
			int innerMinCol, int innerMaxCol, int qLen,
			long[] diag_km2, long[] diag_km1, long[] diag_k) {

		if(debug) {System.err.print(".");}

		// Calculate range size
		final int rangeSize = innerMaxCol - innerMinCol + 1;

		// Create constant vectors once
		final LongVector matchVector=LongVector.broadcast(LSPECIES, MATCH);
		final LongVector subVector=LongVector.broadcast(LSPECIES, SUB);
		final LongVector nScoreVector=LongVector.broadcast(LSPECIES, N_SCORE);
		final LongVector insVector=LongVector.broadcast(LSPECIES, INS);
		final LongVector delIncrementVector=LongVector.broadcast(LSPECIES, DEL);
		final LongVector fifteenVector=LongVector.broadcast(LSPECIES, 15);
		final LongVector scoreMaskVector=LongVector.broadcast(LSPECIES, SCORE_MASK);

		// Check if range is large enough for vector processing
		if(rangeSize >= LWIDTH) {
			// Process all complete vector blocks except potentially the last one
			for(int col=innerMinCol; col <= innerMaxCol - LWIDTH; col+=LWIDTH) {
				// Calculate base index for revQuery
				final int queryBaseIdx=qLen+col-k;

				// Load LWIDTH elements from each array with their respective offsets
				final LongVector qVector=LongVector.fromArray(LSPECIES, revQuery, queryBaseIdx);
				final LongVector rVector=LongVector.fromArray(LSPECIES, ref, col-1);
				final LongVector diagVector=LongVector.fromArray(LSPECIES, diag_km2, col-2);
				final LongVector upVector=LongVector.fromArray(LSPECIES, diag_km1, col-1);
				final LongVector leftVector=LongVector.fromArray(LSPECIES, diag_km1, col-2);

				// Check for matches and Ns
				final VectorMask<Long> isMatchMask=qVector.compare(VectorOperators.EQ, rVector);
				final VectorMask<Long> hasNMask=qVector.or(rVector).compare(VectorOperators.GT, fifteenVector);

				// Calculate score to add - start with SUB
				LongVector scoreAddVector=subVector;

				// Apply the masks using vector blend operations
				scoreAddVector=scoreAddVector.blend(nScoreVector, hasNMask);
				scoreAddVector=scoreAddVector.blend(matchVector, isMatchMask);

				// Calculate scores
				final LongVector diagScoreVector=diagVector.add(scoreAddVector);
				final LongVector upScoreVector=upVector.add(insVector);
				final LongVector leftScoreVector=leftVector.add(delIncrementVector);

				// Find maximum scores
				final LongVector maxDiagUpVector=diagScoreVector.max(upScoreVector);

				// Half-masked comparison 
				final LongVector maskedMaxVector=maxDiagUpVector.and(scoreMaskVector);
				final VectorMask<Long> comparisonMask=maskedMaxVector.compare(
						VectorOperators.GE, leftScoreVector);

				// Select final max
				final LongVector resultVector=leftScoreVector.blend(maxDiagUpVector, comparisonMask);

				// Store results directly
				resultVector.intoArray(diag_k, col-1);

				if(debug) {
					for(int i=0; i<LWIDTH; i++) {
						int currentCol=col+i;
						int row=k-currentCol;
						if(row >= 2 && row <= qLen && currentCol>=2 && currentCol<=innerMaxCol) {
							System.err.println("Cell ("+row+","+currentCol+") calculation (SIMD):");
							System.err.println("*Setting ("+row+","+currentCol+") (diag_k["+(currentCol-1)+"]) to "+
									resultVector.lane(i)+" ("+(resultVector.lane(i)>>42)+")");
						}
					}
				}
			}

			// Process the last elements with one final vector operation
			// Position the vector to exactly include the last element
			final int lastVectorStart = Math.max(innerMinCol, innerMaxCol + 1 - LWIDTH);

			// Only do final vector if we haven't already processed these elements
			if(lastVectorStart > innerMaxCol - LWIDTH) {
				// Calculate base index for revQuery for last vector
				final int queryBaseIdx=qLen+lastVectorStart-k;

				// Load LWIDTH elements from each array with their respective offsets
				final LongVector qVector=LongVector.fromArray(LSPECIES, revQuery, queryBaseIdx);
				final LongVector rVector=LongVector.fromArray(LSPECIES, ref, lastVectorStart-1);
				final LongVector diagVector=LongVector.fromArray(LSPECIES, diag_km2, lastVectorStart-2);
				final LongVector upVector=LongVector.fromArray(LSPECIES, diag_km1, lastVectorStart-1);
				final LongVector leftVector=LongVector.fromArray(LSPECIES, diag_km1, lastVectorStart-2);

				// Check for matches and Ns
				final VectorMask<Long> isMatchMask=qVector.compare(VectorOperators.EQ, rVector);
				final LongVector orVector=qVector.or(rVector);
				final VectorMask<Long> hasNMask=orVector.compare(VectorOperators.GT, fifteenVector);

				// Calculate score to add - start with SUB
				LongVector scoreAddVector=subVector;

				// Apply the masks using vector blend operations
				scoreAddVector=scoreAddVector.blend(nScoreVector, hasNMask);
				scoreAddVector=scoreAddVector.blend(matchVector, isMatchMask);

				// Calculate scores
				final LongVector diagScoreVector=diagVector.add(scoreAddVector);
				final LongVector upScoreVector=upVector.add(insVector);
				final LongVector leftScoreVector=leftVector.add(delIncrementVector);

				// Find maximum scores
				final LongVector maxDiagUpVector=diagScoreVector.max(upScoreVector);

				// Half-masked comparison 
				final LongVector maskedMaxVector=maxDiagUpVector.and(scoreMaskVector);
				final VectorMask<Long> comparisonMask=maskedMaxVector.compare(
						VectorOperators.GE, leftScoreVector);

				// Select final max
				final LongVector resultVector=leftScoreVector.blend(maxDiagUpVector, comparisonMask);

				// Store results directly
				resultVector.intoArray(diag_k, lastVectorStart-1);

				if(debug) {
					for(int i=0; i<LWIDTH; i++) {
						int currentCol=lastVectorStart+i;
						if(currentCol <= innerMaxCol) { // Only debug output for columns in our range
							int row=k-currentCol;
							if(row >= 2 && row <= qLen && currentCol>=2) {
								System.err.println("Cell ("+row+","+currentCol+") calculation (Final SIMD):");
								System.err.println("*Setting ("+row+","+currentCol+") (diag_k["+(currentCol-1)+"]) to "+
										resultVector.lane(i)+" ("+(resultVector.lane(i)>>42)+")");
							}
						}
					}
				}
			}
		} else {
			assert(innerMaxCol-innerMinCol<=4);
			// Range too small for vector processing, use scalar code
			for(int col=innerMinCol; col<=innerMaxCol; col++) {
				int row=k-col;

				if(row < 2 || row > qLen) continue; // Skip invalid rows

				final long q=revQuery[qLen-row];
				final long r=ref[col-1];

				final boolean isMatch=(q==r);
				final boolean hasN=(q|r)>15;
				long scoreAdd=isMatch?MATCH:(hasN?N_SCORE:SUB);

				long diagValue=diag_km2[col-2];
				long upValue=diag_km1[col-1];
				long leftValue=diag_km1[col-2];

				long diagScore=diagValue+scoreAdd;
				long upScore=upValue+INS;
				long leftScore=leftValue+DEL;

				final long maxDiagUp=Math.max(diagScore, upScore);
				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore?maxDiagUp:leftScore;

				diag_k[col-1]=maxValue;

				if(debug) {
					System.err.println("Cell ("+row+","+col+") calculation (scalar):");
					System.err.println("*Setting ("+row+","+col+") (diag_k["+(col-1)+"]) to "+maxValue+" ("+(maxValue>>42)+")");
				}
			}
		}
	}
	
	@SuppressWarnings("restriction")
	public static void processCrossCutDiagonalSIMD_old(
			long[] revQuery, long[] ref, int k, 
			int innerMinCol, int innerMaxCol, int qLen,
			long[] diag_km2, long[] diag_km1, long[] diag_k,
			long MATCH, long SUB, long INS, long DEL_INCREMENT, 
			long N_SCORE, long SCORE_MASK) {

		if(debug) {System.err.print(".");}

		// Calculate range size
		final int rangeSize = innerMaxCol - innerMinCol + 1;

		// Create constant vectors once
		final LongVector matchVector=LongVector.broadcast(LSPECIES, MATCH);
		final LongVector subVector=LongVector.broadcast(LSPECIES, SUB);
		final LongVector nScoreVector=LongVector.broadcast(LSPECIES, N_SCORE);
		final LongVector insVector=LongVector.broadcast(LSPECIES, INS);
		final LongVector delIncrementVector=LongVector.broadcast(LSPECIES, DEL_INCREMENT);
		final LongVector fifteenVector=LongVector.broadcast(LSPECIES, 15);
		final LongVector scoreMaskVector=LongVector.broadcast(LSPECIES, SCORE_MASK);

		// Check if range is large enough for vector processing
		if(rangeSize >= LWIDTH) {
			// Process all complete vector blocks except potentially the last one
			for(int col=innerMinCol; col <= innerMaxCol - LWIDTH; col+=LWIDTH) {
				// Calculate base index for revQuery
				final int queryBaseIdx=qLen+col-k;

				// Load LWIDTH elements from each array with their respective offsets
				final LongVector qVector=LongVector.fromArray(LSPECIES, revQuery, queryBaseIdx);
				final LongVector rVector=LongVector.fromArray(LSPECIES, ref, col-1);
				final LongVector diagVector=LongVector.fromArray(LSPECIES, diag_km2, col-2);
				final LongVector upVector=LongVector.fromArray(LSPECIES, diag_km1, col-1);
				final LongVector leftVector=LongVector.fromArray(LSPECIES, diag_km1, col-2);

				// Check for matches and Ns
				final VectorMask<Long> isMatchMask=qVector.compare(VectorOperators.EQ, rVector);
				final LongVector orVector=qVector.or(rVector);
				final VectorMask<Long> hasNMask=orVector.compare(VectorOperators.GT, fifteenVector);

				// Calculate score to add - start with SUB
				LongVector scoreAddVector=subVector;

				// Apply the masks using vector blend operations
				scoreAddVector=scoreAddVector.blend(nScoreVector, hasNMask);
				scoreAddVector=scoreAddVector.blend(matchVector, isMatchMask);

				// Calculate scores
				final LongVector diagScoreVector=diagVector.add(scoreAddVector);
				final LongVector upScoreVector=upVector.add(insVector);
				final LongVector leftScoreVector=leftVector.add(delIncrementVector);

				// Find maximum scores
				final LongVector maxDiagUpVector=diagScoreVector.max(upScoreVector);

				// Half-masked comparison 
				final LongVector maskedMaxVector=maxDiagUpVector.and(scoreMaskVector);
				final VectorMask<Long> comparisonMask=maskedMaxVector.compare(
						VectorOperators.GE, leftScoreVector);

				// Select final max
				final LongVector resultVector=leftScoreVector.blend(maxDiagUpVector, comparisonMask);

				// Store results directly
				resultVector.intoArray(diag_k, col-1);

				if(debug) {
					for(int i=0; i<LWIDTH; i++) {
						int currentCol=col+i;
						int row=k-currentCol;
						if(row >= 2 && row <= qLen && currentCol>=2 && currentCol<=innerMaxCol) {
							System.err.println("Cell ("+row+","+currentCol+") calculation (SIMD):");
							System.err.println("*Setting ("+row+","+currentCol+") (diag_k["+(currentCol-1)+"]) to "+
									resultVector.lane(i)+" ("+(resultVector.lane(i)>>42)+")");
						}
					}
				}
			}

			// Process the last elements with one final vector operation
			// Position the vector to exactly include the last element
			final int lastVectorStart = Math.max(innerMinCol, innerMaxCol + 1 - LWIDTH);

			// Only do final vector if we haven't already processed these elements
			if(lastVectorStart > innerMaxCol - LWIDTH) {
				// Calculate base index for revQuery for last vector
				final int queryBaseIdx=qLen+lastVectorStart-k;

				// Load LWIDTH elements from each array with their respective offsets
				final LongVector qVector=LongVector.fromArray(LSPECIES, revQuery, queryBaseIdx);
				final LongVector rVector=LongVector.fromArray(LSPECIES, ref, lastVectorStart-1);
				final LongVector diagVector=LongVector.fromArray(LSPECIES, diag_km2, lastVectorStart-2);
				final LongVector upVector=LongVector.fromArray(LSPECIES, diag_km1, lastVectorStart-1);
				final LongVector leftVector=LongVector.fromArray(LSPECIES, diag_km1, lastVectorStart-2);

				// Check for matches and Ns
				final VectorMask<Long> isMatchMask=qVector.compare(VectorOperators.EQ, rVector);
				final LongVector orVector=qVector.or(rVector);
				final VectorMask<Long> hasNMask=orVector.compare(VectorOperators.GT, fifteenVector);

				// Calculate score to add - start with SUB
				LongVector scoreAddVector=subVector;

				// Apply the masks using vector blend operations
				scoreAddVector=scoreAddVector.blend(nScoreVector, hasNMask);
				scoreAddVector=scoreAddVector.blend(matchVector, isMatchMask);

				// Calculate scores
				final LongVector diagScoreVector=diagVector.add(scoreAddVector);
				final LongVector upScoreVector=upVector.add(insVector);
				final LongVector leftScoreVector=leftVector.add(delIncrementVector);

				// Find maximum scores
				final LongVector maxDiagUpVector=diagScoreVector.max(upScoreVector);

				// Half-masked comparison 
				final LongVector maskedMaxVector=maxDiagUpVector.and(scoreMaskVector);
				final VectorMask<Long> comparisonMask=maskedMaxVector.compare(
						VectorOperators.GE, leftScoreVector);

				// Select final max
				final LongVector resultVector=leftScoreVector.blend(maxDiagUpVector, comparisonMask);

				// Store results directly
				resultVector.intoArray(diag_k, lastVectorStart-1);

				if(debug) {
					for(int i=0; i<LWIDTH; i++) {
						int currentCol=lastVectorStart+i;
						if(currentCol <= innerMaxCol) { // Only debug output for columns in our range
							int row=k-currentCol;
							if(row >= 2 && row <= qLen && currentCol>=2) {
								System.err.println("Cell ("+row+","+currentCol+") calculation (Final SIMD):");
								System.err.println("*Setting ("+row+","+currentCol+") (diag_k["+(currentCol-1)+"]) to "+
										resultVector.lane(i)+" ("+(resultVector.lane(i)>>42)+")");
							}
						}
					}
				}
			}
		} else {
			assert(innerMaxCol-innerMinCol<=4);
			// Range too small for vector processing, use scalar code
			for(int col=innerMinCol; col<=innerMaxCol; col++) {
				int row=k-col;

				if(row < 2 || row > qLen) continue; // Skip invalid rows

				final long q=revQuery[qLen-row];
				final long r=ref[col-1];

				final boolean isMatch=(q==r);
				final boolean hasN=(q|r)>15;
				long scoreAdd=isMatch?MATCH:(hasN?N_SCORE:SUB);

				long diagValue=diag_km2[col-2];
				long upValue=diag_km1[col-1];
				long leftValue=diag_km1[col-2];

				long diagScore=diagValue+scoreAdd;
				long upScore=upValue+INS;
				long leftScore=leftValue+DEL_INCREMENT;

				final long maxDiagUp=Math.max(diagScore, upScore);
				final long maxValue=(maxDiagUp&SCORE_MASK)>=leftScore?maxDiagUp:leftScore;

				diag_k[col-1]=maxValue;

				if(debug) {
					System.err.println("Cell ("+row+","+col+") calculation (scalar):");
					System.err.println("*Setting ("+row+","+col+") (diag_k["+(col-1)+"]) to "+maxValue+" ("+(maxValue>>42)+")");
				}
			}
		}
	}

}
