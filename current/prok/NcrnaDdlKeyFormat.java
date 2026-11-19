package prok;

import cardinality.DynamicDemiLog;

/** Explicit absolute DDL key fields for the recipe study.
 * Exponent reduction uses the same bit mask for profiles and queries. Length
 * tags belong in the hash input; no origin metadata accompanies stored keys.
 * No global DDL setting changes.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlKeyFormat{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlKeyFormat(final int exponentBits_, final int mantissaBits_){
		if(exponentBits_<1 || exponentBits_>6 || mantissaBits_<1
				|| exponentBits_+mantissaBits_>16){
			throw new IllegalArgumentException("Absolute key fields must fit a nonzero char, with one to six exponent bits");
		}
		exponentBits=exponentBits_;
		mantissaBits=mantissaBits_;
		maxExponent=(1<<exponentBits)-1;
		maxKey=(1<<(exponentBits+mantissaBits))-1;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Encoding            ----------------*/
	/*--------------------------------------------------------------*/

	/** Encodes an already mixed 64-bit hash, before any bucket maximum or floor.
	 * The exponent wraps by masking, including for an all-zero hash.
	 * The shared encoder preserves the native insertion's shift semantics.
	 */
	public int encode(final long hash){
		final int key=DynamicDemiLog.encodeKey(hash, exponentBits, mantissaBits);
		assert(key>0 && key<=maxKey)
			: "Absolute key fields must preserve zero as empty and stay inside their declared posting keyspace";
		return key;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	public final int exponentBits, mantissaBits, maxExponent, maxKey;
}
