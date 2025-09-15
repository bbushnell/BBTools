package var2;

import shared.Tools;
import structures.ByteBuilder;

/**
 * Performs soft-clipping operations on sequence alignments by identifying
 * poorly aligned terminal regions and converting them to clipped sequences.
 * Uses a scoring system to determine optimal clipping boundaries.
 * 
 * @author Brian Bushnell
 * @contributor Isla Winglet
 */
public class SoftClipper {

	/**
	 * Applies soft clipping to a match string by identifying poorly aligned regions
	 * at the ends and converting them to clipped bases. Uses dynamic scoring to
	 * find the optimal alignment region and clips everything outside it.
	 * 
	 * @param match Original match string representing alignment operations
	 * @param minClipLength Minimum number of bases required to perform clipping
	 * @param allowMutation Whether to modify the original match array or create a copy
	 * @param oldStart Original start position of the alignment
	 * @param oldStop Original stop position of the alignment
	 * @param startStopRvec Return vector for adjusted start/stop positions [start, stop]
	 * @return Modified match string with soft clipping applied, or original if no clipping needed
	 */
	public static byte[] softClipMatch(byte[] match, int minClipLength, boolean allowMutation, 
			final int oldStart, final int oldStop, final int[] startStopRvec){

		// Scoring parameters for alignment quality assessment
		final int matchScore=100;      // Score for perfect matches
		final int subScore=-200;       // Penalty for substitutions
		final int subScore2=-100;      // Reduced penalty for consecutive substitutions
		final int insScore=-200;       // Penalty for insertions
		final int delScore=-200;       // Penalty for deletions
		final int delScore2=-10;       // Reduced penalty for consecutive deletions
		final int clipScore=-1;        // Small penalty for clipping
		final int nScore=1;            // Score for N bases (ambiguous matches)

		int insCount=0;
		int delCount=0;
		
		// Track scoring to find optimal alignment region
		long score=0;
		long maxScore=0;
		int maxPos=-1;                 // Position of maximum score (end of optimal region)
		int maxStart=-1;               // Start position of optimal scoring region
		int currentStart=-1;           // Current region start position
		byte current='?';              // Previous match character for consecutive penalty logic
		
		// Scan through match string to find highest-scoring contiguous region
		for(int mpos=0; mpos<match.length; mpos++){
			final byte m=match[mpos];
			
			if(m=='m' || m=='N' || m=='R'){
				// Good alignment - start new region if score was zero
				if(score==0){currentStart=mpos;}
				
				score=score+(m=='m' ? matchScore : nScore);
				
				// Track the best scoring region
				if(score>maxScore){
					maxScore=score;
					maxPos=mpos;
					maxStart=currentStart;
				}
			}else{
				// Poor alignment - apply penalties
				if(m=='S' || m=='s'){
					// Substitution penalty (reduced for consecutive subs)
					score=score+(m==current ? subScore2 : subScore);
				}else if(m=='D'){
					// Deletion penalty (reduced for consecutive dels)
					score=score+(m==current ? delScore2 : delScore);
					delCount++;
				}else if(m=='I' || m=='X' || m=='Y'){
					// Insertion penalty
					score=score+insScore;
					insCount++;
				}else if(m=='C'){
					// Already clipped penalty
					score=score+clipScore;
				}
				// Don't allow negative scores (start fresh)
				score=Tools.max(0, score);
			}
			current=m;
		}
		
		// If no good alignment region found, return original
		if(maxScore<1){return match;}
		
		// Calculate clipping lengths
		final int leftClipM=maxStart;                    // Match positions to clip on left
		final int rightClipM=(match.length-maxPos-1);   // Match positions to clip on right
		int leftClip=0, rightClip=0;                    // Actual base positions to clip
		
		// Count bases to be clipped (excluding deletions which don't consume read bases)
		for(int i=0; i<match.length; i++){
			byte m=match[i];
			if(i<maxStart){
				leftClip+=(m=='D' ? 0 : 1);
			}else if(i>maxPos){
				rightClip+=(m=='D' ? 0 : 1);
			}
		}
		
		// Only clip if we meet minimum length requirements
		if(leftClip<minClipLength && rightClip<minClipLength){return match;}
		
		int start=oldStart, stop=oldStop;
		
		// Simple case: no deletions, just replace ends with 'C' (clipped)
		if(delCount==0){
			final byte[] array=allowMutation ? match : match.clone();
			for(int i=0; i<leftClip; i++){array[i]='C';}
			for(int i=0, j=array.length-1; i<rightClip; i++, j--){array[j]='C';}
			startStopRvec[0]=start;
			startStopRvec[1]=stop;
			return array;
		}
		
		// Complex case: deletions present, need to rebuild match string and adjust coordinates
		ByteBuilder bb=new ByteBuilder(match.length);
		
		// Handle left clipping with coordinate adjustment
		if(leftClip>=minClipLength){
			for(int mpos=0, processed=0; mpos<match.length; mpos++){
				byte m=match[mpos];
				if(mpos>=leftClipM){
					bb.append(m);
				}else{
					if(m=='D'){
						// Deletion: advance reference position
						start++;
					}else if(m=='I'){
						// Insertion: retreat reference position, mark as clipped
						start--;
						bb.append('C');
						processed++;
					}else{
						// Other operations: mark as clipped
						bb.append('C');
						processed++;
					}
				}
			}
		}else{
			bb.append(match);
		}
		
		// Handle right clipping with coordinate adjustment (process in reverse)
		if(rightClip>=minClipLength){
			bb.reverseInPlace();
			byte[] temp=bb.toBytes();
			bb.clear();
			for(int mpos=0, processed=0; mpos<temp.length; mpos++){
				byte m=temp[mpos];
				if(mpos>=rightClipM){
					bb.append(m);
				}else{
					if(m=='D'){
						// Deletion: retreat reference end position
						stop--;
					}else if(m=='I'){
						// Insertion: advance reference end position, mark as clipped
						stop++;
						bb.append('C');
						processed++;
					}else{
						// Other operations: mark as clipped
						bb.append('C');
						processed++;
					}
				}
			}
			bb.reverseInPlace();
		}
		
		// Return adjusted coordinates
		startStopRvec[0]=start;
		startStopRvec[1]=stop;
		return bb.toBytes();
	}
}