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
	/**
	 * Creates an alignment result with core metrics and boundary coordinates.
	 *
	 * @param maxScore_ Maximum alignment score achieved
	 * @param maxQpos_ Query position at maximum score
	 * @param maxRpos_ Reference position at maximum score
	 * @param qLen_ Total query sequence length
	 * @param rLen_ Total reference sequence length
	 * @param rStart_ Reference alignment start position (inclusive)
	 * @param rStop_ Reference alignment stop position (exclusive)
	 * @param ratio_ Alignment quality ratio for filtering
	 */
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
	/** Maximum alignment score achieved during alignment */
	public int maxScore;
	/** Query sequence position where maximum score was achieved */
	public int maxQpos;
	/** Reference sequence position where maximum score was achieved */
	public int maxRpos;
	/** Total length of the query sequence */
	public int qLen;
	/** Total length of the reference sequence */
	public int rLen;
	/** Reference alignment start position (inclusive) */
	public int rStart;
	/** Reference alignment stop position (exclusive) */
	public int rStop;
	/** Orientation flag for split alignments and structural variant detection */
	public boolean left;
	/** Junction location for split alignments and chimeric read processing */
	public int junctionLoc;
	/** Alignment quality ratio used for filtering decisions and quality control */
	public float ratio;
	/** True if this is an ice cream cone */
	public boolean icecream=false;
	/** Flag indicating alignment ambiguity requiring additional analysis */
	public boolean ambiguous=false;
	/**
	 * Reference to the aligned Read object maintaining connection with source data
	 */
	public Read alignedRead;
}