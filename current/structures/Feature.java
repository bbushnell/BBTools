package structures;

/**
 * Interface for genomic features with coordinate and metadata information.
 * Provides essential methods for accessing feature boundaries, strand orientation,
 * sequence identifier, and feature name. Used throughout BBTools for representing
 * genomic annotations like genes, exons, and other biological features.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public interface Feature {
	
	/** Gets the start position of this genomic feature.
	 * @return Start coordinate (typically 0-based or 1-based depending on context) */
	public int start();
	/** Gets the stop position of this genomic feature.
	 * @return Stop coordinate (inclusive or exclusive depending on context) */
	public int stop();
	/** Gets the strand orientation of this genomic feature.
	 * @return Strand identifier (typically +1 for forward, -1 for reverse, 0 for unknown) */
	public int strand();
	/** Gets the sequence identifier for this genomic feature.
	 * @return Sequence ID (chromosome name, contig ID, or other sequence identifier) */
	public String seqid();
	/** Gets the name of this genomic feature.
	 * @return Feature name or identifier */
	public String name();
	
}
