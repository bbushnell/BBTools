package stream;

/** Byte-range operations for read identifiers, independent of decoding and shared settings.
 * @author Shinobu
 * @date September 29, 2026
 */
public final class ReadHeader{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Static-only utility. */
	private ReadHeader(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Finds the exclusive end of an identifier without copying or decoding its bytes.
	 * Trimming uses the legacy FASTA Character.isWhitespace(byte) predicate; leading
	 * whitespace produces an empty range and signed high bytes are not delimiters.
	 * This method does not remove a marker or validate header syntax. Callers retain
	 * their existing charset and decide whether an empty range becomes an empty ID.
	 * @param header Nonnull header bytes; checked by assertion
	 * @param start First identifier byte, from zero through header.length inclusive;
	 *        callers skip any format marker before passing this index
	 * @param trimDescription Whether to stop at the first whitespace byte
	 * @return First whitespace index at or after start when trimming, otherwise header.length
	 */
	public static int end(final byte[] header, final int start, final boolean trimDescription){
		assert(header!=null && start>=0 && start<=header.length) :
			"Identifier decoding requires a valid header range; start="+start+
			", length="+(header==null ? -1 : header.length);
		if(trimDescription){
			for(int i=start; i<header.length; i++){
				if(Character.isWhitespace(header[i])){return i;}
			}
		}
		return header.length;
	}

}
