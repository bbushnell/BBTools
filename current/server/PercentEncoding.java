package server;

import java.nio.charset.StandardCharsets;
import java.util.BitSet;
import java.util.HashMap;

/**
 * Handles targeted percent-encoding and decoding for server URL fragments.
 * Escapes selected reserved/common ASCII symbols while leaving non-ASCII text unchanged.
 * Decodes percent-encoded ASCII and UTF-8 bytes while preserving malformed percent syntax.
 * @author Brian Bushnell
 */
public class PercentEncoding {
	
/*--------------------------------------------------------------*/
/*----------------        Static Methods        ----------------*/
/*--------------------------------------------------------------*/

	/**
	 * Checks if a string contains any special symbols that require percent-encoding.
	 * Special symbols include both reserved URI characters and common problematic characters.
	 * @param s The string to check for special symbols
	 * @return true if the string contains any special symbols, false otherwise
	 */
	public static boolean containsSpecialSymbol(String s){
		if(s==null){return false;}
		for(int i=0, max=s.length(); i<max; i++){
			char c=s.charAt(i);
			if(isSpecial.get(c)){return true;}
		}
		return false;
	}
	
	/**
	 * Checks if a string contains any common symbols that require percent-encoding.
	 * Common symbols are a subset of special symbols including newlines, spaces, quotes, etc.
	 * @param s The string to check for common symbols
	 * @return true if the string contains any common symbols, false otherwise
	 */
	public static boolean containsCommonSymbol(String s){
		if(s==null){return false;}
		for(int i=0, max=s.length(); i<max; i++){
			char c=s.charAt(i);
			if(isCommon.get(c)){return true;}
		}
		return false;
	}
	
	/**
	 * Converts all special symbols in a string to their percent-encoded equivalents.
	 * Leaves non-ASCII characters unchanged.
	 * Returns the original string unchanged if no special symbols are found.
	 * @param s The string to encode
	 * @return The percent-encoded string with special symbols converted to %XX format
	 */
	public static String symbolToCode(String s){
		if(!containsSpecialSymbol(s)){return s;}
		StringBuilder sb=new StringBuilder(s.length());
		for(int i=0, max=s.length(); i<max; i++){
			char c=s.charAt(i);
			//FIXED [server/PercentEncoding#003]: symbolToCodeArray is ASCII-sized, so the unguarded lookup
			//crashed on any non-ASCII char that appeared with an encodable ASCII symbol.
			String code=(c<128 ? symbolToCodeArray[c] : null);
			if(code!=null){
				sb.append(code);
			}else{
				sb.append(c);
			}
		}
		return sb.toString();
	}
	
	/**
	 * Converts only common symbols in a string to their percent-encoded equivalents.
	 * Leaves non-ASCII characters unchanged.
	 * Returns the original string unchanged if no common symbols are found.
	 * @param s The string to encode
	 * @return The percent-encoded string with common symbols converted to %XX format
	 */
	public static String commonSymbolToCode(String s){
		if(!containsCommonSymbol(s)){return s;}
		StringBuilder sb=new StringBuilder(s.length());
		for(int i=0, max=s.length(); i<max; i++){
			char c=s.charAt(i);
			if(isCommon.get(c)){
				String code=symbolToCodeArray[c];
				assert(code!=null) : "isCommon and symbolToCodeArray both derive from commonSymbol; missing code for "+(int)c;
				sb.append(code);
			}else{
				sb.append(c);
			}
		}
		return sb.toString();
	}
	
	/**
	 * Parses a percent-encoded character sequence (%XX) into its numeric value.
	 *
	 * @param s The string containing the percent-encoded sequence
	 * @param start The index of the percent sign in the string
	 * @return The decoded byte value, or -1 unless start points to %XX with two hex digits
	 */
	private static int parseCode(String s, int start){
		if(s==null || start+2>=s.length()){return -1;}
		if(s.charAt(start)!='%'){return -1;}
		assert(s.charAt(start)=='%');
		int sum=0;
		for(int i=start+1; i<=start+2; i++){
			sum=sum<<4;
			final char c=s.charAt(i);
			if(c>='0' && c<='9'){
				sum=sum+(c-'0');
			//FIXED [server/PercentEncoding#001]: was 'c<'F'' which excluded hex digit F, so %2F (/) and %3F (?) failed to decode (parseCode returned -1), breaking the encode/decode round-trip on the TaxClient->TaxServer path and dna/Data URL decode.  Now 'c<='F''.
			}else if(c>='A' && c<='F'){
				sum=sum+(10+c-'A');
			//FIXED [server/PercentEncoding#002]: also accept lowercase a-f (RFC 3986: decoders should accept both cases).  BBTools encodes uppercase, but external clients may send lowercase.
			}else if(c>='a' && c<='f'){
				sum=sum+(10+c-'a');
			}else{
				return -1;
			}
		}
		return sum;
	}
	
	/**
	 * Converts targeted percent-encoded byte sequences back to symbols.
	 * Decodes ASCII escapes directly, decodes contiguous percent-byte runs that
	 * begin with a high byte as UTF-8,
	 * and leaves malformed UTF-8 replacement to the Java UTF-8 decoder.
	 * Preserves malformed percent syntax as literal characters.
	 * Returns the original string unchanged if no percent-encoded sequences are found.
	 *
	 * @param s The string containing percent-encoded sequences
	 * @return The decoded string with %XX sequences converted back to original characters
	 */
	public static String codeToSymbol(String s){
		int idx=s.indexOf('%');
		if(idx<0){return s;}
		
		StringBuilder sb=new StringBuilder(s.length());
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='%'){
				int sym=parseCode(s, i);
				if(sym<0){sb.append(c);}
				else if(sym<128){
					sb.append((char)sym);
					i+=2;//Skip next 2 characters
				}else{
					int start=i;
					int len=0;
					while(i<s.length() && s.charAt(i)=='%' && (sym=parseCode(s, i))>=0){
						len++;
						i+=3;
					}
					byte[] bytes=new byte[len];
					for(int j=0; j<len; j++){
						bytes[j]=(byte)parseCode(s, start+j*3);
					}
					sb.append(new String(bytes, StandardCharsets.UTF_8));
					i--;//Unskip the first invalid triplet or non-percent character
				}
			}else{sb.append(c);}
		}
		return (sb.length()==s.length() ? s : sb.toString());
	}
	
	/**
	 * Creates a mapping from percent-encoded strings to their original symbols.
	 * Combines both reserved and common symbol mappings into a single HashMap.
	 * @return HashMap mapping percent codes to original symbols
	 */
	private static HashMap<String, String> makeCodeToSymbolMap(){
		HashMap<String, String> map=new HashMap<String, String>(129);
		assert(reservedSymbol.length==reservedCode.length);
		assert(commonSymbol.length==commonCode.length);
		for(int i=0; i<reservedSymbol.length; i++){
			map.put(reservedCode[i], reservedSymbol[i]);
		}
		for(int i=0; i<commonSymbol.length; i++){
			map.put(commonCode[i], commonSymbol[i]);
		}
		return map;
	}
	
	/**
	 * Creates a mapping from original symbols to their percent-encoded equivalents.
	 * Combines both reserved and common symbol mappings into a single HashMap.
	 * @return HashMap mapping original symbols to percent codes
	 */
	private static HashMap<String, String> makeSymbolToCodeMap(){
		HashMap<String, String> map=new HashMap<String, String>(257);
		assert(reservedSymbol.length==reservedCode.length);
		assert(commonSymbol.length==commonCode.length);
		for(int i=0; i<reservedSymbol.length; i++){
			map.put(reservedSymbol[i], reservedCode[i]);
		}
		for(int i=0; i<commonSymbol.length; i++){
			map.put(commonSymbol[i], commonCode[i]);
		}
		return map;
	}
	
	/**
	 * Creates a fast-lookup array mapping ASCII characters to percent codes.
	 * Uses character values as array indices for O(1) lookup performance.
	 * @return Array where array[c] contains the percent code for character c, or null
	 */
	private static String[] makeSymbolToCodeArray(){
		final String[] array=new String[128];
		for(int i=0; i<reservedSymbol.length; i++){
			String s=reservedSymbol[i];
			String c=reservedCode[i];
			array[s.charAt(0)]=c;
		}
		for(int i=0; i<commonSymbol.length; i++){
			String s=commonSymbol[i];
			String c=commonCode[i];
			array[s.charAt(0)]=c;
		}
		return array;
	}
	
	/**
	 * Creates a BitSet for fast character membership testing.
	 * Sets bits corresponding to the first character of each string in the input arrays.
	 * @param matrix Variable number of string arrays containing characters to include
	 * @return BitSet with bits set for all characters found in the input arrays
	 */
	private static final BitSet makeBitSet(String[]...matrix){
		BitSet bs=new BitSet(128);
		for(String[] array : matrix){
			for(String s : array){
				char c=s.charAt(0);
				bs.set(c);
			}
		}
		return bs;
	}

/*--------------------------------------------------------------*/
/*----------------          Constants           ----------------*/
/*--------------------------------------------------------------*/

	//See https://en.wikipedia.org/wiki/Percent-encoding
	/**
	 * Reserved URI characters that require percent-encoding according to RFC 3986
	 */
	public static final String[] reservedSymbol=new String[] {
		"!", "#", "$", "&", "'", "(", ")", "*", "+", ",", "/", ":", ";", "=", "?", "@", "[", "]"
	};
	
	/** Percent-encoded equivalents for reserved URI characters */
	public static final String[] reservedCode=new String[] {
		"%21", "%23", "%24", "%26", "%27", "%28", "%29", "%2A", "%2B", "%2C", "%2F", "%3A", "%3B", "%3D", "%3F", "%40", "%5B", "%5D"
	};
	
	/** Common problematic characters that often need percent-encoding */
	public static final String[] commonSymbol=new String[] {
		"\n", " ", "\"", "%", "<", ">", "\\", "|",
	};
	
	/** Percent-encoded equivalents for common problematic characters */
	public static final String[] commonCode=new String[] {
		"%0A", "%20", "%22", "%25", "%3C", "%3E", "%5C", "%7C"
	};

	/** BitSet for fast lookup of special characters requiring encoding */
	private static final BitSet isSpecial=makeBitSet(reservedSymbol, commonSymbol);
	/** BitSet for fast lookup of common characters requiring encoding */
	private static final BitSet isCommon=makeBitSet(commonSymbol);

//	public static final HashMap<String, String> codeToSymbolMap=makeCodeToSymbolMap();
//	public static final HashMap<String, String> symbolToCodeMap=makeSymbolToCodeMap();
	/** Array for fast lookup of percent codes by character value */
	public static final String[] symbolToCodeArray=makeSymbolToCodeArray();
	
	/** Don't print caught exceptions */
	public static boolean suppressErrors=false;
	
}
