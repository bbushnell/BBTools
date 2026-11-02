package scalar;

/**
 * Represents a single row/interval of scalar compositional data.
 * Used for row-based access to ScalarData.
 * Fields are mutable to allow object reuse. Instances are not thread-safe.
 * Coordinate/length interpretation comes from the ScalarData producer.
 *
 * @author Neptune
 * @date Feb 15, 2026
 */
public class ScalarInterval{

	/** Sequence/contig name */
	public String name;

	/** Producer-supplied length; sequence sampling currently counts newly seen valid bases. */
	public float length;

	/** GC content (0-1) */
	public float gc;

	/** GC-normalized homopolymer dimer metric (0-1), not Shannon entropy. */
	public float hh;

	/** GC-normalized CAGA dimer metric (0-1). */
	public float caga;

	/** Read depth/coverage */
	public float depth;

	/** Zero-based inclusive metric-span start; -1 when unavailable. */
	public int start;
	/** Exclusive metric-span end; -1 when unavailable or absent from an older table. */
	public int end=-1;

	/** Primary taxonomy ID; nonpositive values indicate no assignment. */
	public int taxID;

	/** Secondary taxonomy ID (0 if unavailable) */
	public int taxID2;

	/** Default constructor */
	public ScalarInterval(){}

	/**
	 * Populate this interval with values from ScalarData at given index.
	 * Core GC/HH/CAGA lists must be populated. Optional lists may be null; otherwise
	 * they must contain the requested row. The table itself is not retained.
	 * @param data Source ScalarData
	 * @param index Row index
	 */
	public void set(ScalarData data, int index){
		this.name=(data.names!=null ? data.names.get(index) : null);
		this.length=(data.length!=null ? data.length.get(index) : 0);
		//scalar/ScalarInterval#001: the three core metric columns are required by contract.
		this.gc=data.gc.get(index);
		this.hh=data.hh.get(index);
		this.caga=data.caga.get(index);
		this.depth=(data.depth!=null ? data.depth.get(index) : 0);
		this.start=(data.start!=null ? data.start.get(index) : 0);
		this.end=(data.end!=null && data.end.size()>0 ? data.end.get(index) : -1);
		this.taxID=(data.taxIDs!=null ? data.taxIDs.get(index) : 0);
		this.taxID2=(data.taxIDs2!=null ? data.taxIDs2.get(index) : 0);
	}

}
