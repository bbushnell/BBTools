package prot;

/**
 * Shared reversible identifier codec (converged Ady/UMP45, Yoimiya assignment 2026-09-09): real
 * contig/shred identifiers can contain any of the byte structural delimiters this project uses
 * to build tab-separated files, semicolon-joined lists, comma-joined lists, and pipe-separated
 * tokens (UMP45's real-data counts: 547,044 comma-bearing ids, 626 semicolon-bearing ids, plus
 * pipe-bearing prefix ids). {@link #encode(String)}/{@link #decode(String)} let every consumer
 * (MagQCBinManifest's selection tokens, ReferenceCdsSurvivalSidecar's containing_shred_ids,
 * ReferenceCdsSidecarIndex's selected_source_span_ids, ReferenceCdsSurvivalLabelReader's
 * comma-split ids) share ONE escaping rule instead of diverging.
 *
 * <p>Exactly seven ASCII characters are ever escaped -- every structural delimiter used anywhere
 * in this project's identifier-bearing formats, plus the escape character itself:
 * {@code \t , ; | \r \n %}. Every other character, including every non-ASCII/Unicode character
 * and plain space, passes through {@link #encode(String)} completely untouched -- this operates
 * on Java {@code char}s, not UTF-8 bytes, so no character can ever collide with one of the seven
 * reserved ASCII values by accident, and no Unicode text is ever re-encoded. Escaping is
 * {@code %} followed by two uppercase hex digits of the character's value (0-127 for all seven
 * reserved characters, matching plain percent-encoding).
 *
 * <p><b>Versioning (per-file schema-version bump, not a header flag):</b> this class does not
 * decide when a file's identifiers are encoded -- each consumer format gates that itself via its
 * own existing {@code #schema_version} mechanism (MagQCBinManifest v1→v2; the sidecar/index
 * formats v1→v2 on UMP45's side), so a legacy v1 file is read exactly as before, unchanged, and
 * only a v2 file has its identifiers decoded. {@link #decode(String)} is strict either way: a
 * bare trailing {@code %} or a {@code %} not followed by exactly two valid hex digits is a
 * malformed-input error, never silently passed through or dropped -- crash-loud on corruption is
 * the point of versioning this explicitly.
 *
 * @author Ady
 */
public final class IdentifierCodec {
	private IdentifierCodec(){ }

	/** Names this codec's escaping rule, for documentation/log purposes only -- versioning
	 *  itself is done by each consumer format's own schema-version field, not by this string. */
	public static final String ENCODING="percent_v1";

	/** The seven reserved characters, indexed by ASCII value (0-127) for O(1) lookup; true
	 *  means "this character is escaped by encode() and is a structural delimiter somewhere
	 *  in this project's identifier-bearing formats." */
	private static final boolean[] RESERVED=buildReserved();
	private static boolean[] buildReserved(){
		final boolean[] r=new boolean[128];
		r['%']=true; r[',']=true; r[';']=true; r['|']=true;
		r['\t']=true; r['\n']=true; r['\r']=true;
		return r;
	}

	private static final char[] HEX="0123456789ABCDEF".toCharArray();

	/** True iff {@code id} contains at least one of the seven reserved characters, i.e.
	 *  {@code encode(id) != id}. Zero-cost fast path for the common case (UMP45's real-data
	 *  measurement: the large majority of ids encode to themselves). */
	public static boolean needsEncoding(String id){
		if(id==null){throw new IllegalArgumentException("null id");}
		for(int i=0; i<id.length(); i++){
			final char c=id.charAt(i);
			if(c<128 && RESERVED[c]){return true;}
		}
		return false;
	}

	/** Escapes every reserved character in {@code id} as {@code %XX} (uppercase hex); every
	 *  other character, including all non-ASCII/Unicode characters and plain space, is copied
	 *  through unchanged. Never returns null for a non-null input; the empty string encodes to
	 *  itself. */
	public static String encode(String id){
		if(id==null){throw new IllegalArgumentException("null id");}
		if(!needsEncoding(id)){return id;}//fast path: no allocation beyond the scan already done
		final StringBuilder sb=new StringBuilder(id.length()+8);
		for(int i=0; i<id.length(); i++){
			final char c=id.charAt(i);
			if(c<128 && RESERVED[c]){
				sb.append('%').append(HEX[(c>>4)&0xF]).append(HEX[c&0xF]);
			}else{
				sb.append(c);
			}
		}
		return sb.toString();
	}

	/** True iff {@code decode(s)} would succeed without throwing -- every {@code %} in {@code s}
	 *  is followed by exactly two valid hex digits. Does not check whether {@code s} contains
	 *  any actual escape sequence (see {@link #needsEncoding} for that on the DECODED side);
	 *  this validates well-formedness only. */
	public static boolean isEncoded(String s){
		if(s==null){throw new IllegalArgumentException("null id");}
		return firstMalformedPercent(s)<0;
	}

	/** Returns the index of the first malformed {@code %} escape in {@code s} (a trailing bare
	 *  {@code %}, or a {@code %} not followed by two hex digits), or -1 if every {@code %} in
	 *  {@code s} starts a well-formed two-hex-digit escape. */
	private static int firstMalformedPercent(String s){
		for(int i=0; i<s.length(); i++){
			if(s.charAt(i)=='%'){
				if(i+2>=s.length() || !isHex(s.charAt(i+1)) || !isHex(s.charAt(i+2))){return i;}
				i+=2;//skip the two hex digits; the loop's i++ advances past the '%' itself
			}
		}
		return -1;
	}
	private static boolean isHex(char c){
		return (c>='0' && c<='9') || (c>='A' && c<='F') || (c>='a' && c<='f');
	}

	/** Reverses {@link #encode(String)}: each well-formed {@code %XX} becomes the character of
	 *  that hex value; every other character is copied through unchanged. Strict: throws
	 *  {@link IllegalArgumentException} naming the exact position of a bare trailing {@code %}
	 *  or a {@code %} not followed by two valid hex digits -- malformed input is never silently
	 *  passed through or truncated. Accepts (and decodes) any valid {@code %XX}, not only the
	 *  seven bytes {@link #encode} itself would have produced, matching standard
	 *  percent-decoding semantics. */
	public static String decode(String enc){
		if(enc==null){throw new IllegalArgumentException("null id");}
		final int bad=firstMalformedPercent(enc);
		if(bad>=0){
			throw new IllegalArgumentException("Malformed "+ENCODING+" escape at position "+bad
				+" in: "+enc);
		}
		if(enc.indexOf('%')<0){return enc;}//fast path: no escapes at all
		final StringBuilder sb=new StringBuilder(enc.length());
		for(int i=0; i<enc.length(); i++){
			final char c=enc.charAt(i);
			if(c=='%'){
				final int hi=hexVal(enc.charAt(i+1)), lo=hexVal(enc.charAt(i+2));
				sb.append((char)((hi<<4)|lo));
				i+=2;
			}else{
				sb.append(c);
			}
		}
		return sb.toString();
	}
	private static int hexVal(char c){
		if(c>='0' && c<='9'){return c-'0';}
		if(c>='A' && c<='F'){return c-'A'+10;}
		return c-'a'+10;//isHex already restricted callers to [0-9A-Fa-f]
	}
}
