package scalar;

/**
 * Represents a single row/interval of scalar compositional data.
 * Used for row-based access to ScalarData.
 * Fields are mutable to allow object reuse.
 *
 * @author Neptune
 * @date Feb 15, 2026
 */
public class ScalarInterval {

	/** Sequence/contig name */
	public String name;

	/** Actual contig/interval length in bases */
	public float length;

	/** GC content (0-1) */
	public float gc;

	/** Homopolymer entropy (0-1) */
	public float hh;

	/** CAGA skew metric (0-1) */
	public float caga;

	/** Read depth/coverage */
	public float depth;

	/** Start position within contig (0 for whole contigs) */
	public int start;

	/** Primary taxonomy ID (0 if unavailable) */
	public int taxID;

	/** Secondary taxonomy ID (0 if unavailable) */
	public int taxID2;

	/** Default constructor */
	public ScalarInterval() {}

	/**
	 * Populate this interval with values from ScalarData at given index.
	 * @param data Source ScalarData
	 * @param index Row index
	 */
	public void set(ScalarData data, int index) {
		this.name = (data.names != null ? data.names.get(index) : null);
		this.length = (data.length != null ? data.length.get(index) : 0);
		this.gc = data.gc.get(index);
		this.hh = data.hh.get(index);
		this.caga = data.caga.get(index);
		this.depth = (data.depth != null ? data.depth.get(index) : 0);
		this.start = (data.start != null ? data.start.get(index) : 0);
		this.taxID = (data.taxIDs != null ? data.taxIDs.get(index) : 0);
		this.taxID2 = (data.taxIDs2 != null ? data.taxIDs2.get(index) : 0);
	}

}
