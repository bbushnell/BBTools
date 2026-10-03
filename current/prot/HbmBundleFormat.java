package prot;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.zip.CRC32;

/**
 * Constants and shared helpers for the deterministic HBM model-bundle binary format.
 *
 * <p>Implements the accepted byte-level contract (design
 * {@code hbm_bundle_format_proposal_v5_final_ump45_20260915.md}, sha256
 * {@code 9d749188…}). All multi-byte integers are big-endian with no alignment padding.</p>
 *
 * <p><b>Java-8-safe (BBTools "no lazy sysadmin left behind").</b> The per-block checksum uses
 * {@link java.util.zip.CRC32} (standard CRC-32, Java 1.1+), NOT {@code CRC32C} (Java-9+); strong
 * integrity is SHA-256 via {@link java.security.MessageDigest}.</p>
 *
 * <p><b>Crash-loud in production, not assertion-gated.</b> Public helpers validate their arguments
 * with explicit {@code if}/throw (not {@code assert}), so a bad call fails the same way with
 * {@code -ea} or {@code -da}; the magic is exposed only through copy-returning access so no caller
 * can mutate a process-wide array.</p>
 *
 * @author UMP45
 */
public final class HbmBundleFormat {

	private HbmBundleFormat(){}

	/*--------------------------------------------------------------*/
	/*----------------      Magic and version       ----------------*/
	/*--------------------------------------------------------------*/

	private static final byte[] MAGIC={(byte)'M', (byte)'Q', (byte)'H', (byte)'B'};//"MQHB"; private → immutable
	/** Magic length in bytes. */
	public static final int MAGIC_LEN=4;
	/** A fresh copy of the magic marker (offset 0, and the trailer end sentinel). */
	public static byte[] magic(){return MAGIC.clone();}
	/** True iff {@code data[off..off+4)} equals the magic. Range/null checked. */
	public static boolean magicMatches(final byte[] data, final int off){
		checkRange(data, off, MAGIC_LEN, "magicMatches");
		for(int i=0; i<MAGIC_LEN; i++){if(data[off+i]!=MAGIC[i]){return false;}}
		return true;
	}
	/** Only format version this loader/builder speaks. */
	public static final int FORMAT_VERSION=1;

	/*--------------------------------------------------------------*/
	/*----------------      Fixed-header offsets     ----------------*/
	/*--------------------------------------------------------------*/
	// §1.1: literal offsets, no alignment padding.
	public static final int OFF_MAGIC=0;
	public static final int OFF_FORMAT_VERSION=4;        // u16
	public static final int OFF_NAA=6;                   // u16
	public static final int OFF_PAD=8;                   // u16
	public static final int OFF_N_FAMILIES=10;           // i32
	public static final int OFF_WEIGHT_BY_IDENTITY=14;   // u8
	public static final int OFF_BASE_WEIGHT=15;          // i32
	public static final int OFF_MAF_SUB=19;              // f32
	public static final int OFF_MAF_DEL=23;              // f32
	public static final int OFF_MAF_INS=27;              // f32
	public static final int OFF_MIN_DEPTH=31;            // i32
	public static final int OFF_IDENTITY_CEILING=35;     // f32
	public static final int OFF_TRIM_DEPTH_FRACTION=39;  // f32
	public static final int OFF_PROVENANCE=43;           // PROVENANCE_COUNT * SHA256_LEN

	public static final int PROVENANCE_COUNT=16;
	public static final int SHA256_LEN=32;
	/** Fixed-header byte count = 43 + 16*32 = 555. */
	public static final int FIXED_HEADER_LEN=OFF_PROVENANCE+PROVENANCE_COUNT*SHA256_LEN;
	public static final int DIRECTORY_START=FIXED_HEADER_LEN;

	/*--------------------------------------------------------------*/
	/*----------------   Provenance field indices    ----------------*/
	/*--------------------------------------------------------------*/
	// §1.2 fixed order. 0..7 RUNTIME-SEMANTIC, 8..9 BUILD-AUDIT, 10..11 INPUT runtime-compat, 12..15 INPUT audit.
	public static final int PROV_FORMAT_SPEC=0, PROV_LOADER_SOURCE=1, PROV_LOADER_CLASS=2, PROV_AAGRAPH=3,
		PROV_AAGRAPHNODE=4, PROV_AAGRAPHSCORER=5, PROV_GLOCALAMINOLINEAR=6, PROV_BLOSUM62=7,
		PROV_BUILDER_SOURCE=8, PROV_BUILDER_CLASS=9, PROV_ROSTER=10, PROV_CONSENSUS_REF=11,
		PROV_SOURCE_CORPUS=12, PROV_CLUSTER_MEMBERSHIP=13, PROV_MEMBER_POLICY=14, PROV_MEMBER_MANIFEST=15;
	/** Runtime-semantic provenance indices 0..7 (checked at load vs the trusted channel, §4). */
	public static final int RUNTIME_SEMANTIC_COUNT=8;

	/*--------------------------------------------------------------*/
	/*----------------      Canonical build knobs    ----------------*/
	/*--------------------------------------------------------------*/
	/** Residue-histogram length; source-grounded. AAGraphNode.NAA == Blosum62.X_CODE+1 == 22. */
	public static final int NAA=AAGraphNode.NAA;
	public static final int EXPECT_PAD=0;
	public static final boolean KNOB_WEIGHT_BY_IDENTITY=false;
	public static final int KNOB_BASE_WEIGHT=1;
	public static final float KNOB_MAF_SUB=0.25f;
	public static final float KNOB_MAF_DEL=0.5f;
	public static final float KNOB_MAF_INS=0.5f;
	public static final int KNOB_MIN_DEPTH=1;
	public static final float KNOB_IDENTITY_CEILING=40f;
	public static final float KNOB_TRIM_DEPTH_FRACTION=0f;

	/*--------------------------------------------------------------*/
	/*----------------      Resource caps (§2)       ----------------*/
	/*--------------------------------------------------------------*/
	public static final int CAP_L=100000;
	public static final int CAP_CHAIN=100000;
	public static final int CAP_REPID=512;
	public static final long CAP_TOTAL_BYTES=4L*1024*1024*1024;
	public static final long CAP_TOTAL_NODES=2000000000L;

	/*--------------------------------------------------------------*/
	/*----------------        Record sizes           ----------------*/
	/*--------------------------------------------------------------*/
	public static final int COUNT_ARRAY_BYTES=NAA*4;
	public static final int DIR_ENTRY_FIXED_BYTES=24;
	public static final int INS_RECORD_BYTES=4+COUNT_ARRAY_BYTES;

	/*--------------------------------------------------------------*/
	/*----------------          Helpers              ----------------*/
	/*--------------------------------------------------------------*/

	/** Overflow-safe range check for a byte array (crash-loud, not assertion-gated). */
	private static void checkRange(final byte[] data, final int off, final int len, final String who){
		if(data==null){throw new IllegalArgumentException(who+": null data");}
		if(off<0 || len<0 || off>data.length || len>data.length-off){
			throw new IllegalArgumentException(who+": bad range off="+off+" len="+len+" size="+data.length);
		}
	}

	/** Standard CRC-32 (Java-8-safe) of a byte range, as the raw unsigned 32-bit value in a long (0..2^32-1). */
	public static long crc32(final byte[] data, final int off, final int len){
		checkRange(data, off, len, "crc32");
		final CRC32 crc=new CRC32();
		crc.update(data, off, len);
		return crc.getValue();
	}

	/** SHA-256 of a byte range as 32 raw bytes. A range error throws {@link IllegalArgumentException};
	 *  only a genuinely absent SHA-256 provider is reported as unavailable. */
	public static byte[] sha256(final byte[] data, final int off, final int len){
		checkRange(data, off, len, "sha256");
		final MessageDigest md;
		try{md=MessageDigest.getInstance("SHA-256");}
		catch(NoSuchAlgorithmException e){throw new RuntimeException("SHA-256 unavailable", e);}
		md.update(data, off, len);
		return md.digest();
	}

	/** Lowercase 64-hex of a 32-byte SHA-256 digest. Length is checked with an explicit throw (not assert). */
	public static String toHexLower(final byte[] digest){
		if(digest==null || digest.length!=SHA256_LEN){
			throw new IllegalArgumentException("expected 32-byte digest, got "+(digest==null?"null":String.valueOf(digest.length)));
		}
		final char[] hex="0123456789abcdef".toCharArray();
		final char[] out=new char[digest.length*2];
		for(int i=0; i<digest.length; i++){
			final int b=digest[i]&0xff;
			out[i*2]=hex[b>>>4];
			out[i*2+1]=hex[b&0xf];
		}
		return new String(out);
	}

	/** True iff the 32-byte field at {@code off} is all zero (the all-zero loader-hash placeholder,
	 *  §1.2 — rejected in a built bundle). Null/range checked. */
	public static boolean isAllZero(final byte[] sha, final int off){
		checkRange(sha, off, SHA256_LEN, "isAllZero");
		for(int i=0; i<SHA256_LEN; i++){if(sha[off+i]!=0){return false;}}
		return true;
	}
}
