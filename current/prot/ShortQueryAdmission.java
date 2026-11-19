package prot;

/** Prespecified exact-K5 short-query admission floor, anchored at six hits for99 residues.
 * Length means normalized amino-acid residues, not occupied hashes or surviving k-mers.
 * The caller must select this policy explicitly; existing fixed-floor searches do not use it.
 * @author Keqing
 */
public final class ShortQueryAdmission{

	private ShortQueryAdmission(){}

	/** Integer ceiling of6*(length-4)/95 for the released40..99 interval.
	 * This validates a length rather than filtering or clipping a sequence. */
	public static int rejectFloor(final int length){
		if(length<40 || length>99){
			throw new IllegalArgumentException("Length-scaled exact-K5 admission requires40..99 normalized residues; received "+length);
		}
		final int floor=(6*(length-4)+94)/95;
		assert(floor>=3 && floor<=6) : "The released40..99 formula must produce only3..6 hits; length="+length;
		return floor;
	}
}
