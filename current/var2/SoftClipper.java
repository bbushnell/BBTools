package var2;

import shared.Tools;
import structures.ByteBuilder;

/**
 * Scores an expanded alignment-operation string and converts selected terminal
 * operations to soft clips. Used by VarHelper's junction caller after Var.toVars
 * expands the match string. This heuristic scans existing operations; it does
 * not realign bases or validate the operation alphabet.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @contributor Sucrose
 */
public class SoftClipper{

	/**
	 * Finds the first strictly best-scoring run encountered, resetting negative
	 * running scores to zero. Matches score 100; N/R score 1; S/s cost 200 initially
	 * and 100 when the identical symbol repeats; D costs 200 initially and 10 on
	 * repetition; I/X/Y cost 200 and existing C costs 1.
	 * Terminal clip lengths count non-D operations, not all match-array entries.
	 * Returns early when no positive score exists or both clip lengths are below
	 * minClipLength; those returns leave the input and output coordinates untouched.
	 *
	 * Without any D, reaching the clipping path converts both selected ends to C
	 * even if only one end reaches the threshold. Each clipped I decrements start
	 * on the left or increments stop on the right; other symbols leave them unchanged.
	 * This path mutates match only when allowMutation is true, otherwise clones it.
	 * With any D, a new array is built regardless of allowMutation and each end must
	 * separately meet the threshold. Clipped D operations are removed; clipped I
	 * operations become C. On the left these adjust start by +1/-1 respectively;
	 * on the right they adjust stop by -1/+1. Other clipped operations become C
	 * without coordinate adjustment.
	 *
	 * @param match Nonnull expanded match operations, one byte per operation
	 * @param minClipLength Threshold in read bases, with branch-specific end handling above
	 * @param allowMutation Permits in-place edits only on the deletion-free clipping path
	 * @param oldStart Original zero-based reference start
	 * @param oldStop Original zero-based inclusive reference stop
	 * @param startStopRvec Output [start, stop]; needs at least two entries on a clipping path,
	 * and is not initialized on early return
	 * @return Original array on early return or in-place clipping; otherwise a new array
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

		// Track the highest running score encountered; ties retain the earlier run.
		long score=0;
		long maxScore=0;
		int maxPos=-1;                 // End position of the best scoring run encountered
		int maxStart=-1;               // Start position of that run
		int currentStart=-1;           // Current region start position
		byte current='?';              // Previous match character for consecutive penalty logic

		// Scan expanded operations; this does not decode run-length counts.
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

		// At least one end must meet the threshold before either clipping path runs.
		if(leftClip<minClipLength && rightClip<minClipLength){return match;}

		int start=oldStart, stop=oldStop;

		// No deletions: replace both selected ends, including a shorter opposite end.
		// Preserve the junction position when an insertion becomes a clip (VAR2-007).
		if(delCount==0){
			final byte[] array=allowMutation ? match : match.clone();
			for(int i=0; i<leftClip; i++){
				if(array[i]=='I'){start--;}
				array[i]='C';
			}
			for(int i=0, j=array.length-1; i<rightClip; i++, j--){
				if(array[j]=='I'){stop++;}
				array[j]='C';
			}
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
