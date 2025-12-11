package aligner;

import stream.Read;

/**
 * Container for sequence alignment results and metadata.
 * Stores alignment scores, positions, and specialized flags for downstream processing.
 * Used by alignment algorithms to return comprehensive alignment information including
 * peak scores, sequence coordinates, quality metrics, and structural variant detection.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class AlignmentResult {
	public AlignmentResult(int maxScore_, int maxQpos_, int maxRpos_, int qLen_, int rLen_, int rStart_, int rStop_, float ratio_) {
		maxScore=maxScore_;
		maxQpos=maxQpos_;
		maxRpos=maxRpos_;
		qLen=qLen_;
		rLen=rLen_;
		rStart=rStart_;
		rStop=rStop_;
		ratio=ratio_;
	}
	public int maxScore;
	public int maxQpos;
	public int maxRpos;
	public int qLen;
	public int rLen;
	public int rStart;
	public int rStop;
	public boolean left;
	public int junctionLoc;
	public float ratio;
	/** True if this alignment exhibits ice cream cone geometry pattern */
	public boolean icecream=false;
	public boolean ambiguous=false;
	public Read alignedRead;
}