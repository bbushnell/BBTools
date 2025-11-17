package server;

import java.util.BitSet;
import java.util.HashMap;

import structures.ByteBuilder;

/**
 * Handles URL percent-encoding and decoding for special characters.
 * Provides encoding/decoding of reserved and common symbols according to RFC 3986.
 * Supports both full encoding of special symbols and targeted encoding of common symbols.
 * @author Brian Bushnell
 */
public class PercentEncoding {
	
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
			if(isSpecial.get(c)){
//				System.err.print("b");
				return true;}
		}
//		System.err.print("c");
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
	 * Returns the original string unchanged if no special symbols are found.
	 * @param s The string to encode
	 * @return The percent-encoded string with special symbols converted to %XX format
	 */
	public static String symbolToCode(String s){
//		System.err.print("a");
		if(!containsSpecialSymbol(s)){return s;}
//		System.err.print("d");
		ByteBuilder bb=new ByteBuilder();
		for(int i=0, max=s.length(); i<max; i++){
			char c=s.charAt(i);
			String code=symbolToCodeArray[c];
			if(code!=null){
//				System.err.print("e("+code+")");
				bb.append(code);
			}else{
//				System.err.print("f");
				bb.append(c);
			}
		}
//		System.err.println("g");
//		System.err.println(bb);
		return bb.toString();
	}
	
	/**
	 * Converts only common symbols in a string to their percent-encoded equivalents.
	 * Returns the original string unchanged if no common symbols are found.
	 * @param s The string to encode
	 * @return The percent-encoded string with common symbols converted to %XX format
	 */
	public static String commonSymbolToCode(String s){
		if(!containsCommonSymbol(s)){return s;}
		ByteBuilder bb=new ByteBuilder();
		for(int i=0, max=s.length(); i<max; i++){
			char c=s.charAt(i);
			if(isCommon.get(c)){
				String code=symbolToCodeArray[c];
				assert(code!=null);
				bb.append(code);
			}else{
				bb.append(c);
			}
		}
		return bb.toString();
	}
	
	/**
	 * Parses a percent-encoded character sequence (%XX) into its numeric value.
	 * Expects a percent sign followed by two hexadecimal digits.
	 *
	 * @param s The string containing the percent-encoded sequence
	 * @param start The index of the percent sign in the string
	 * @return The decoded numeric value, or -1 if parsing fails
	 */
	private static int parseCode(String s, int start){
		if(s==null || start+2>=s.length()){return -1;}
		assert(s.charAt(start)=='%');
		int sum=0;
		for(int i=start+1; i<=start+2; i++){
			sum=sum<<4;
			final char c=s.charAt(i);
			if(c>='0' && c<='9'){
				sum=sum+(c-'0');
			}else if(c>='A' && c<'F'){
				sum=sum+(10+c-'A');
			}else{
				return -1;
			}
		}
		return sum;
	}
	
	/**
	 * Converts percent-encoded sequences in a string back to their original symbols.
	 * Decodes all valid %XX sequences found in the string.
	 * Returns the original string unchanged if no percent-encoded sequences are found.
	 *
	 * @param s The string containing percent-encoded sequences
	 * @return The decoded string with %XX sequences converted back to original characters
	 */
	public static String codeToSymbol(String s){
		int idx=s.indexOf('%');
		if(idx<0){return s;}
		
		ByteBuilder bb=new ByteBuilder(s.length());
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='%'){
				int sym=parseCode(s, i);
				if(sym<0){bb.append(c);}
				else{
					bb.append((char)sym);
					i+=2;//Skip next 2 characters
				}
			}else{bb.append(c);}
		}
		return (bb.length()==s.length() ? s : bb.toString());
	}
	
	/**
	 * Creates a mapping from percent-encoded strings to their original symbols.
	 * Combines both reserved and common symbol mappings into a single HashMap.
	 * @return HashMap mapping percent codes to original symbols
	 */
	private static HashMap<String, String> makeCodeToSymbolMap() {
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
	private static HashMap<String, String> makeSymbolToCodeMap() {
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
	private static String[] makeSymbolToCodeArray() {
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
