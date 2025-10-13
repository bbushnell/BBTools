package clade;

/**
 * Context object for connection-specific parameters in CladeServer.
 * Ensures thread-safe handling of concurrent requests by giving
 * each connection its own isolated parameter set.
 *
 * @author Chloe
 * @date September 16, 2025
 */
public class CladeContext {

	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	/** Default constructor with default values */
	public CladeContext() {
		// Default values
	}

	/** Copy constructor for cloning contexts */
	public CladeContext(CladeContext other) {
		this.format = other.format;
		this.hits = other.hits;
		this.heap = other.heap;
		this.printQTID = other.printQTID;
		this.banSelf = other.banSelf;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Output format (HUMAN or ONELINE) */
	public int format = CladeSearcher.HUMAN;

	/** Number of hits to return per query */
	public int hits = 1;

	/** Heap size for comparisons */
	public int heap = 1;

	/** Print query TaxID */
	public boolean printQTID = false;

	/** Ban self-matches */
	public boolean banSelf = false;
}