package prot;

import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Locale;

/**
 * Runtime format identity and checksum primitives for dense HBM text v2.
 * Family checksums cover their f/body rows; the root covers all headers and
 * ordered family-checksum rows. Each digest consumes UTF-8 line bytes plus LF.
 * Sixteen stored semantic pins preserve the native runtime/build bindings.
 * Offline conversion lives in HbmDenseTextPacker, not this runtime dependency.
 * @author Collei
 */
public final class HbmDenseTextBundle {

	static final String FORMAT="hbm_dense_v2";
	static final String CONTRACT="mqhb_v1_canonical_knobs_pad0_naa22_weight_equals_count";

	/** A v1 dump cannot enter the production loader merely by being compressed. */
	static boolean matches(final Path path){
		final String name=path.toString().toLowerCase(Locale.ROOT);
		return name.endsWith(".hbmt") || name.endsWith(".hbmt.gz");
	}

	/** CRLF and LF input have identical semantics after ByteFile normalization. */
	static void update(final MessageDigest md, final byte[] line){
		if(md!=null){md.update(line); md.update((byte)'\n');}
	}
	static MessageDigest digest(){
		try{return MessageDigest.getInstance("SHA-256");}
		catch(NoSuchAlgorithmException e){throw new IllegalStateException("Required SHA-256 is unavailable", e);}
	}
	private HbmDenseTextBundle(){}
}
