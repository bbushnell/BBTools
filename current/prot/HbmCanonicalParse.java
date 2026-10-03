package prot;

import java.nio.ByteBuffer;
import java.nio.charset.CharacterCodingException;
import java.nio.charset.CharsetDecoder;
import java.nio.charset.CodingErrorAction;
import java.nio.charset.StandardCharsets;

/**
 * Shared canonical/strict parsing helpers for Increment 3B slice-1 text artifacts (roster,
 * cluster rows, receipts, manifests): strict UTF-8 decode (rejects, never silently replaces,
 * invalid byte sequences), canonical nonnegative decimal integers (no sign, no leading zero
 * except the literal "0", overflow-checked), and canonical lowercase-hex SHA-256 validation.
 * Standalone -- mirrors patterns already accepted elsewhere in this project
 * ({@code HbmMemberManifest}'s {@code utf8StrictDecode}/{@code checkHexLower}/
 * {@code parseIntStrict}) without depending on that class.
 *
 * @author Eru
 */
public final class HbmCanonicalParse {

	private HbmCanonicalParse(){}

	/** Strict UTF-8 decode: rejects malformed/unmappable byte sequences (throws), never
	 *  silently substitutes a replacement character. */
	public static String utf8StrictDecode(final byte[] raw, final int off, final int len){
		final CharsetDecoder dec=StandardCharsets.UTF_8.newDecoder()
			.onMalformedInput(CodingErrorAction.REPORT).onUnmappableCharacter(CodingErrorAction.REPORT);
		try{return dec.decode(ByteBuffer.wrap(raw, off, len)).toString();}
		catch(CharacterCodingException e){throw new RuntimeException("UTF8_INVALID: malformed or unmappable UTF-8 bytes: "+e.getMessage());}
	}

	public static String utf8StrictDecode(final byte[] raw){return utf8StrictDecode(raw, 0, raw.length);}

	/** Canonical nonnegative decimal {@code long}: digits only, no leading '+'/'-', no
	 *  leading zero unless the value is exactly "0". Throws on any other form or overflow. */
	public static long parseCanonicalNonNegLong(final String s, final String fieldName){
		if(s==null || s.isEmpty()){throw new RuntimeException("CANON_INT: empty "+fieldName);}
		for(int i=0; i<s.length(); i++){
			final char c=s.charAt(i);
			if(c<'0' || c>'9'){throw new RuntimeException("CANON_INT: non-digit character in "+fieldName+": '"+s+"'");}
		}
		if(s.length()>1 && s.charAt(0)=='0'){throw new RuntimeException("CANON_INT: leading zero in "+fieldName+": '"+s+"'");}
		try{return Long.parseLong(s);}
		catch(NumberFormatException e){throw new RuntimeException("CANON_INT: overflow or unparseable "+fieldName+": '"+s+"'");}
	}

	/** Canonical lowercase 64-hex SHA-256 string check (throws if not exactly that shape). */
	public static void checkLowercaseHexSha256(final String s, final String fieldName){
		if(s==null || s.length()!=64){throw new RuntimeException("HEX_SHA: "+fieldName+" must be 64 lowercase hex characters, got '"+s+"'");}
		for(int i=0; i<64; i++){
			final char c=s.charAt(i);
			final boolean ok=(c>='0'&&c<='9')||(c>='a'&&c<='f');
			if(!ok){throw new RuntimeException("HEX_SHA: "+fieldName+" contains a non-lowercase-hex character: '"+s+"'");}
		}
	}

	/** Canonical lowercase 20-hex textual suffix of SHA-256 (80 retained bits). */
	public static void checkLowercaseHexSha80(final String s, final String fieldName){
		if(s==null || s.length()!=DigestSuffix.HEX_LENGTH){
			throw new RuntimeException("HEX_SHA80: "+fieldName+" must be 20 lowercase hex characters, got '"+s+"'");
		}
		for(int i=0; i<DigestSuffix.HEX_LENGTH; i++){
			final char c=s.charAt(i);
			final boolean ok=(c>='0'&&c<='9')||(c>='a'&&c<='f');
			if(!ok){throw new RuntimeException("HEX_SHA80: "+fieldName+" contains a non-lowercase-hex character: '"+s+"'");}
		}
	}

	/** Overflow-checked byte-index scan for a single-byte delimiter. */
	public static int indexOf(final byte[] a, final byte b, final int fromIndex){
		for(int i=fromIndex; i<a.length; i++){if(a[i]==b){return i;}}
		return -1;
	}
}
