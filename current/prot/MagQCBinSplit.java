package prot;

/**
 * Stateless train/validation routing shared by independent D39 generators.
 * The input is a canonical SHA-256 source-selection or clean-complete genome key, not a model
 * name, sampling seed, row number or target organism alone. All cooperating
 * jobs must use the same modulus and pinned source/table indexing.
 *
 * @author Yoimiya
 */
final class MagQCBinSplit {

	private MagQCBinSplit(){}

	/** One bucket in 26 matches the starting 1,000,000/40,000 subnet row ratio. */
	static final int DEFAULT_MODULUS=26;
	static final String VERSION="source_index_and_complete_tid_sha256_v1";

	/** Rejects unusable routing fractions before a generator opens any output. */
	static void validate(int modulus){
		if(modulus<2 || modulus>1000){throw new IllegalArgumentException("splitmod must be in [2,1000]: "+modulus);}
	}

	/**
	 * Bucket zero is reserved for validation. A uniform 63-bit digest prefix
	 * gives a reproducible partition with negligible modulo bias. Requested
	 * output counts remain exact: generators retry candidates in the other split.
	 */
	static boolean validation(byte[] fingerprint,int modulus){
		assert(fingerprint!=null && fingerprint.length==32) : "D39 source fingerprints must be complete SHA-256 digests";
		assert(modulus>=2 && modulus<=1000) : "Constructor validation bounds the shared split modulus";
		long value=fingerprint[0]&127;
		for(int i=1; i<8; i++){value=(value<<8)|(fingerprint[i]&255);}
		return value%modulus==0;
	}
}
