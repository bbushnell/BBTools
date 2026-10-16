package structures;

/**
 * Genomic feature coordinates and metadata used by range and sequence consumers.
 * Endpoints are zero-based and inclusive, in the sequence frame exposed by the
 * implementation. Strand codes use 0 for forward and 1 for reverse; additional
 * codes depend on the implementation, such as GffLine's 2 for '?' and 3 for '.'.
 *
 * @author Brian Bushnell
 * @contributor Shinobu
 * @date 2013
 */
public interface Feature{

	/**
	 * Gets the first position of this feature.
	 * @return Zero-based inclusive start coordinate.
	 */
	public int start();
	/**
	 * Gets the last position of this feature.
	 * @return Zero-based inclusive stop coordinate.
	 */
	public int stop();
	/**
	 * Gets the implementation's stored strand code.
	 * @return 0 for forward, 1 for reverse, or an implementation-specific code.
	 */
	public int strand();
	/**
	 * Gets the identifier used to locate the containing sequence.
	 * @return Sequence identifier supplied by the implementation.
	 */
	public String seqid();
	/**
	 * Gets optional descriptive metadata for this feature.
	 * @return Implementation-specific name or metadata, possibly null.
	 */
	public String name();

}
