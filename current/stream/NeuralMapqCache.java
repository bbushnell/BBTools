package stream;

/**
 * Stores loose and strict neural MAPQs in two seven-bit components of a Read short.
 * Each component encodes MAPQ+1: zero is absent, and values 1–64 represent MAPQs 0–63.
 * Loose occupies bits 0–6 and strict bits 7–13. A setter preserves all other bits;
 * clear resets the entire short. Getters assume values written through this API and
 * do not validate arbitrary contents of the public Read field.
 * The caller owns synchronization of access to each Read; component updates are
 * separate read-modify-write operations and do not publish an atomic pair of values.
 */
public final class NeuralMapqCache{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Prevents construction of this static utility. */
	private NeuralMapqCache(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Alias for setLoose, with the same range and nonnull requirements. */
	public static void set(final Read read, final int mapq){setLoose(read, mapq);}

	/** Sets the loose value while preserving the strict component and other bits.
	 * @param read Read owning the packed cache
	 * @param mapq Mapping quality in [0,63]
	 * @throws IllegalArgumentException if read is null or mapq is outside [0,63]
	 */
	public static void setLoose(final Read read, final int mapq){
		validate(read, mapq, "Loose");
		read.neuralMapqs=(short)((read.neuralMapqs&~LOOSE_MASK)|(mapq+1));
	}

	/** Sets the strict value while preserving the loose component and other bits.
	 * @param read Read owning the packed cache
	 * @param mapq Mapping quality in [0,63]
	 * @throws IllegalArgumentException if read is null or mapq is outside [0,63]
	 */
	public static void setStrict(final Read read, final int mapq){
		validate(read, mapq, "Strict");
		read.neuralMapqs=(short)((read.neuralMapqs&~STRICT_MASK)|((mapq+1)<<STRICT_SHIFT));
	}

	/** Alias for getLoose; returns -1 for a null Read or absent loose value. */
	public static int get(final Read read){return getLoose(read);}

	/** Returns the loose MAPQ, or -1 for a null Read or absent loose component. */
	public static int getLoose(final Read read){return decode(read, 0);}

	/** Returns the strict MAPQ, or -1 for a null Read or absent strict component. */
	public static int getStrict(final Read read){return decode(read, STRICT_SHIFT);}

	/** Clears the entire cache field, including both components; null is ignored. */
	public static void clear(final Read read){if(read!=null){read.neuralMapqs=0;}}

	/** Decodes one seven-bit component without validating the stored encoding.
	 * @param read Read owning the cache, or null
	 * @param shift Component offset: zero for loose, STRICT_SHIFT for strict
	 * @return Encoded value minus one, or -1 for a null Read
	 */
	private static int decode(final Read read, final int shift){
		if(read==null){return -1;}
		final int encoded=(read.neuralMapqs>>>shift)&COMPONENT_MASK;
		return encoded-1;
	}

	/** Rejects a null destination or an out-of-range value before modifying the cache.
	 * @param read Required destination Read
	 * @param mapq Mapping quality in [0,63]
	 * @param label Component name used in the diagnostic
	 */
	private static void validate(final Read read, final int mapq, final String label){
		if(read==null || mapq<0 || mapq>63){throw new IllegalArgumentException(label+" neural MAPQ outside [0,63]: "+mapq);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Seven-bit extraction/loose masks, strict bit offset and shifted strict mask. */
	private static final int COMPONENT_MASK=127, LOOSE_MASK=COMPONENT_MASK, STRICT_SHIFT=7, STRICT_MASK=COMPONENT_MASK<<STRICT_SHIFT;
}
