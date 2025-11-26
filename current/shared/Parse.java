package shared;

import structures.ByteBuilder;
import structures.LongList;

public final class Parse {
	

	/**
	 * Parses a string with optional KMG suffix into an integer value.
	 * Supports suffixes like k (1000), m (million), g/b (billion), etc.
	 *
	 * @param b String to parse, may include KMG suffix
	 * @return Parsed integer value
	 * @throws AssertionError if the resulting value exceeds integer range
	 */
	public static int parseIntKMG(String b){
		long x=parseKMG(b);
		assert(x<=Integer.MAX_VALUE && x>Integer.MIN_VALUE) : "Value "+x+" is out of range for integers: "+b;
		return (int)x;
	}
	
	/**
	 * Parses a string with optional KMG suffix into a long value.
	 * Supports decimal numbers and special values like "inf", "max", "huge".
	 * KMG suffixes: k=1000, m=million, g/b=billion, t=trillion, p/q=quadrillion,
	 * e=quintillion, c/h=hundred, d=ten.
	 *
	 * @param b0 String to parse, may include KMG suffix or special values
	 * @return Parsed long value, or Long.MAX_VALUE for infinity values
	 * @throws RuntimeException if suffix is unrecognized
	 */
	public static long parseKMG(final String b0){
		String b=b0;
		if(b==null){return 0;}
		assert(b.length()>0);
		final char c=Tools.toLowerCase(b.charAt(b.length()-1));
		final boolean dot=b.indexOf('.')>=0;
		if(!dot && !Tools.isLetter(c)){return Long.parseLong(b);}
//		if(!Tools.isLetter(c) && !dot){return Long.parseLong(b);}
		
		if(b.equalsIgnoreCase("big") || b.equalsIgnoreCase("inf") || b.equalsIgnoreCase("infinity") || b.equalsIgnoreCase("max") || b.equalsIgnoreCase("huge")){
			return Long.MAX_VALUE;
		}
		
		long mult=1;
		if(Tools.isLetter(c)){
			if(c=='k'){mult=1000;}
			else if(c=='m'){mult=1000000;}
			else if(c=='g' || c=='b'){mult=1000000000;}
			else if(c=='t'){mult=1000000000000L;}
			else if(c=='p' || c=='q'){mult=1000000000000000L;}
			else if(c=='e'){mult=1000000000000000000L;}
//			else if(c=='z'){mult=1000000000000000000000L;}//Out of range
			else if(c=='c' || c=='h'){mult=100;}
			else if(c=='d'){mult=10;}
			else{throw new RuntimeException(b);}
			b=b.substring(0, b.length()-1);
		}
		assert(!Tools.endsWithLetter(b)) : "Too many letters at the end of "+b0;
		
		//Calculate product, check for overflow, and return
		if(!dot){
			long m=Long.parseLong(b);
			long p=m*mult;
			assert(p>=m) : p+", "+m+", "+b;
			return p;
		}else{
			double m=Double.parseDouble(b);
			long p=(long)(m*mult);
			assert(p>=m) : p+", "+m+", "+b;
			return p;
		}
	}
	
	/**
	 * Parses a string with optional KMG suffix into a double value.
	 * Similar to parseKMG but returns double precision result.
	 * Additional support for 'f' suffix which is ignored.
	 *
	 * @param b0 String to parse, may include KMG suffix or special values
	 * @return Parsed double value, or Long.MAX_VALUE for infinity values
	 * @throws RuntimeException if suffix is unrecognized
	 */
	public static double parseDoubleKMG(final String b0){
		String b=b0;
		if(b==null){return 0;}
		assert(b.length()>0);
		final char c=Tools.toLowerCase(b.charAt(b.length()-1));
		final boolean dot=b.indexOf('.')>=0;
		if(!dot && !Tools.isLetter(c)){return Long.parseLong(b);}
//		if(!Tools.isLetter(c) && !dot){return Long.parseLong(b);}
		
		if(b.equalsIgnoreCase("big") || b.equalsIgnoreCase("inf") || b.equalsIgnoreCase("infinity") || b.equalsIgnoreCase("max") || b.equalsIgnoreCase("huge")){
			return Long.MAX_VALUE;
		}
		
		long mult=1;
		if(Tools.isLetter(c)){
			if(c=='k'){mult=1000;}
			else if(c=='m'){mult=1000000;}
			else if(c=='g' || c=='b'){mult=1000000000;}
			else if(c=='t'){mult=1000000000000L;}
			else if(c=='p' || c=='q'){mult=1000000000000000L;}
			else if(c=='e'){mult=1000000000000000000L;}
//			else if(c=='z'){mult=1000000000000000000000L;}//Out of range
			else if(c=='c' || c=='h'){mult=100;}
			else if(c=='d'){mult=10;}
			else if(c=='f'){/* Ignore */}
			else{throw new RuntimeException(b);}
			b=b.substring(0, b.length()-1);
		}
		assert(!Tools.endsWithLetter(b)) : "Too many letters at the end of "+b0;
		
		//Calculate product, check for overflow, and return
		if(!dot){
			double m=Double.parseDouble(b);
			double p=m*mult;
			assert(m<0 || p>=m) : p+", "+m+", "+b;
			return p;
		}else{
			double m=Double.parseDouble(b);
			double p=m*mult;
			assert(m<0 || p>=m) : p+", "+m+", "+b;
			return p;
		}
	}
	
	/**
	 * Parses a string with binary KMG suffixes (powers of 1024).
	 * Supports k=1024, m=1024^2, g/b=1024^3, t=1024^4.
	 *
	 * @param b String to parse with binary suffix
	 * @return Parsed long value using binary multipliers
	 * @throws RuntimeException if suffix is unrecognized
	 */
	public static long parseKMGBinary(String b){
		if(b==null){return 0;}
		char c=Tools.toLowerCase(b.charAt(b.length()-1));
		boolean dot=b.indexOf('.')>=0;
		if(!Tools.isLetter(c) && !dot){return Long.parseLong(b);}
		
		long mult=1;
		if(Tools.isLetter(c)){
			if(c=='k'){mult=1024;}
			else if(c=='m'){mult=1024*1024;}
			else if(c=='g' || c=='b'){mult=1024*1024*1024;}
			else if(c=='t'){mult=1024L*1024L*1024L*1024L;}
			else{throw new RuntimeException(b);}
			b=b.substring(0, b.length()-1);
		}
		
		if(!dot){return Long.parseLong(b)*mult;}
		
		return (long)(Double.parseDouble(b)*mult);
	}
	
	/**
	 * Tests if a string represents a number by checking the first character.
	 * @param s String to test
	 * @return true if string starts with digit, '.', or '-'
	 */
	public static boolean isNumber(String s){
		if(s==null || s.length()==0){return false;}
		char c=s.charAt(0);
		return Tools.isDigit(c) || c=='.' || c=='-';
	}
	
	/**
	 * Tests if a string represents a boolean value.
	 * @param s String to test
	 * @return true if string is "t", "f", "true", or "false" (case insensitive)
	 */
	public static boolean isBoolean(String s){
		if(s==null || s.length()==0){return false;}
		return "t".equalsIgnoreCase(s) || "f".equalsIgnoreCase(s) || "true".equalsIgnoreCase(s) || "false".equalsIgnoreCase(s);
	}
	
	/**
	 * Parse this argument.  More liberal than Boolean.parseBoolean.
	 * Null, t, true, or 1 all yield true.
	 * Everything else, including the String "null", is false.
	 * @param s Argument to parse
	 * @return boolean form
	 */
	public static boolean parseBoolean(String s){
		if(s==null || s.length()<1){return true;}
		if(s.length()==1){
			char c=Tools.toLowerCase(s.charAt(0));
			return c=='t' || c=='1';
		}
		if(s.equalsIgnoreCase("null") || s.equalsIgnoreCase("none")){return false;}
		return Boolean.parseBoolean(s);
	}
	
	/**
	 * Parses yes/no strings into boolean values.
	 * Accepts "y"/"yes" for true, "n"/"no" for false.
	 * Special case: "unknown" returns false for IMG database compatibility.
	 *
	 * @param s String to parse
	 * @return true for yes values, false for no values
	 * @throws RuntimeException if string is not recognized yes/no format
	 */
	public static boolean parseYesNo(String s){
		if(s==null || s.length()<1){return true;}
		if(s.length()==1){
			char c=Tools.toLowerCase(s.charAt(0));
			if(c=='y'){return true;}
			if(c=='n'){return false;}
			throw new RuntimeException(s);
		}
		
		if(s.equalsIgnoreCase("yes")){return true;}
		if(s.equalsIgnoreCase("no")){return false;}
		if(s.equalsIgnoreCase("unknown")){return false;} //Special case for IMG database
		
		throw new RuntimeException(s);
	}
	
	public static boolean eic(String a, String b) {//False when both null
		return a==null ? false : a.equalsIgnoreCase(b);
	}
	
	public static boolean equalsIgnoreCase(String a, String b) {//false when both null
		return a==null ? false : a.equalsIgnoreCase(b);
	}
	
	public static boolean equalsIgnoreCase(byte[] a, String b) {//false when both null
		if(a==null || b==null || a.length!=b.length()) {return false;}
		for(int i=0; i<a.length; i++) {
			char c=b.charAt(i);
			if(a[i]!=c) {return false;}
		}
		return true;
	}
	
	public static boolean equalsIgnoreCase(String a, byte[] b) {return equalsIgnoreCase(b, a);}
	
	/**
	 * Parses a delimited string into a float array.
	 * @param s String containing delimited float values
	 * @param regex Delimiter pattern for splitting
	 * @return Array of parsed float values, or null if input is null
	 */
	public static float[] parseFloatArray(String s, String regex){return parseFloatArray(s, regex, null);}
	
	/**
	 * Parses a delimited string into a float array with wildcard support.
	 *
	 * @param s String containing delimited float values
	 * @param regex Delimiter pattern for splitting
	 * @param wildcard String that represents -1 value
	 * @return Array of parsed float values, or null if input is null
	 */
	public static float[] parseFloatArray(String s, String regex, String wildcard){
		if(s==null || "null".equals(s)){return null;}
		String[] split=s.split(regex);
		float[] array=new float[split.length];
		for(int i=0; i<split.length; i++){
			array[i]=(split[i].equals(wildcard)) ? -1 : Float.parseFloat(split[i]);
		}
		return array;
	}
	
	/**
	 * Parses a delimited string into an integer array.
	 * @param s String containing delimited integer values
	 * @param regex Delimiter pattern for splitting
	 * @return Array of parsed integer values, or null if input is null
	 */
	public static int[] parseIntArray(String s, String regex){return parseIntArray(s, regex, null);}
	
	/**
	 * Parses a delimited string into an integer array with wildcard support.
	 *
	 * @param s String containing delimited integer values
	 * @param regex Delimiter pattern for splitting
	 * @param wildcard String that represents -1 value
	 * @return Array of parsed integer values, or null if input is null
	 */
	public static int[] parseIntArray(String s, String regex, String wildcard){
		if(s==null || "null".equals(s)){return null;}
		String[] split=s.split(regex);
		int[] array=new int[split.length];
		for(int i=0; i<split.length; i++){
			array[i]=(split[i].equals(wildcard)) ? -1 : Integer.parseInt(split[i]);
		}
		return array;
	}
	
	/**
	 * Parses a delimited string into a byte array.
	 * @param s String containing delimited byte values
	 * @param regex Delimiter pattern for splitting
	 * @return Array of parsed byte values, or null if input is null
	 */
	public static byte[] parseByteArray(String s, String regex){
		if(s==null){return null;}
		String[] split=s.split(regex);
		byte[] array=new byte[split.length];
		for(int i=0; i<split.length; i++){
			array[i]=Byte.parseByte(split[i]);
		}
		return array;
	}
	
	public static int parseIntHexDecOctBin(final String s){
		if(s==null || s.length()<1){return 0;}
		int radix=10;
		if(s.length()>1 && s.charAt(1)=='0'){
			final char c=s.charAt(1);
			if(c=='x' || c=='X'){radix=16;}
			else if(c=='b' || c=='B'){radix=2;}
			else if(c=='o' || c=='O'){radix=8;}
		}
		return Integer.parseInt(s, radix);
	}
	
	/**
	 * @param array Text
	 * @param a Index of first digit
	 * @param b Index after last digit (e.g., array.length)
	 * @return Parsed number
	 */
	public static float parseFloat(byte[] array, int a, int b){
		return (float)parseDouble(array, a, b);
	}
	
	/**
	 * @param array Text
	 * @param a Index of first digit
	 * @param b Index after last digit (e.g., array.length)
	 * @return Parsed number
	 */
	public static float parseFloat(String array, int a, int b){
		return (float)parseDouble(array, a, b);
	}
	
	/**
	 * Parses a float value following a search term in a string.
	 *
	 * @param s String to search in
	 * @param term Search term to find
	 * @param delimiter Character that terminates the number
	 * @return Parsed float value
	 */
	public static float parseFloat(String s, String term, char delimiter) {
		return (float)parseDouble(s, term, delimiter);
	}

	/**
	 * Parses a double value following a search term in a string.
	 * Finds the term, then parses the number until delimiter is reached.
	 *
	 * @param s String to search in
	 * @param term Search term to find
	 * @param delimiter Character that terminates the number
	 * @return Parsed double value
	 * @throws AssertionError if term is not found
	 */
	public static double parseDouble(String s, String term, char delimiter) {
		int idx=s.indexOf(term);
		assert(idx>=0) : "No "+term+" in String: '"+s+"'";
		idx+=term.length();
		int idx2=Tools.indexOf(s, delimiter, idx);
		idx2=(idx2<0 ? s.length() : idx2);
		double d=Parse.parseDouble(s, idx, idx2);
		return d;
	}

	/**
	 * Extracts a substring following a search term.
	 *
	 * @param s String to search in
	 * @param term Search term to find
	 * @param delimiter Character that terminates the substring
	 * @return Extracted substring, or null if term not found
	 */
	public static String parseString(String s, String term, char delimiter) {
		int idx=s.indexOf(term);
//		assert(idx>=0) : "No "+term+" in String: '"+s+"'";
		if(idx<0) {return null;}
		idx+=term.length();
		int idx2=Tools.indexOf(s, delimiter, idx);
		idx2=(idx2<0 ? s.length() : idx2);
		return s.substring(idx, idx2);
	}

	/**
	 * @param array Text
	 * @param a Index of first digit
	 * @param b Index after last digit (e.g., array.length)
	 * @return Parsed number
	 */
	public static double parseDoubleSlow(byte[] array, int a, int b){
		String s=new String(array, a, b-a);
		return Double.parseDouble(s);
	}

	/**
	 * Parses a double from a byte array starting at given position.
	 * @param array Byte array containing ASCII digits
	 * @param start Starting index
	 * @return Parsed double value
	 */
	public static double parseDouble(final byte[] array, final int start){
		return parseDouble(array, start, array.length);
	}
	
	/**
	 * @param array Text
	 * @param a0 Index of first digit
	 * @param b Index after last digit (e.g., array.length)
	 * @return Parsed number
	 */
	public static double parseDouble(final byte[] array, final int a0, final int b){
		if(Tools.FORCE_JAVA_PARSE_DOUBLE){
			return Double.parseDouble(new String(array, a0, b-a0));
		}
		if(b-a0>1 && b-a0<5) {
			final byte x=array[a0];
			if(!Tools.numericMap[x]){
				if(x == 'N'){return Double.NaN;}
				if(x == 'I'){return Double.POSITIVE_INFINITY;}
			}
			if(x == '-' && array[a0+1] == 'I') 
				return Double.NEGATIVE_INFINITY;
		}
		int a=a0;
		assert(b>a);
		long upper=0;
		final byte z='0';
		long mult=1;
		if(array[a]=='-'){mult=-1; a++;}
		
		for(; a<b; a++){
			final byte c=array[a];
			if(c=='.'){break;}
			final int x=(c-z);
			assert(x<10 && x>=0) : x+" = "+(char)c+"\narray="+new String(array)+", start="+a+", stop="+b;
			upper=(upper*10)+x;
		}
		
		long lower=0;
		int places=0;
		for(a++; a<b; a++){
			final byte c=array[a];
			final int x=(c-z);
			assert(x<10 && x>=0) : x+" = "+(char)c+"\narray='"+new String(array)+"', start="+a+", stop="+b+", len="+array.length+" -> '"+new String(array, a, b-a)+"'"+
				"\nThis function does not support exponents; if the input has an exponent, add the flag 'forceJavaParseDouble'.";
			lower=(lower*10)+x;
			places++;
		}
		
		double d=mult*(upper+lower*ByteBuilder.decimalInvMult[places]);
//		assert(d==parseDoubleSlow(array, a0, b)) : d+", "+parseDoubleSlow(array, a0, b);
		return d;
	}
	
	/**
	 * @param array Text
	 * @param a0 Index of first digit
	 * @param b Index after last digit (e.g., array.length)
	 * @return Parsed number
	 */
	public static double parseDouble(final String array, final int a0, final int b){
		if(Tools.FORCE_JAVA_PARSE_DOUBLE){
			return Double.parseDouble(array.substring(a0, b));
		}
		int a=a0;
		assert(b>a);
		long upper=0;
		final char z='0';
		long mult=1;
		if(array.charAt(a)=='-'){mult=-1; a++;}
		
		for(; a<b; a++){
			final char c=array.charAt(a);
			if(c=='.'){break;}
			final int x=(c-z);
			assert(x<10 && x>=0) : x+" = "+(char)c+"\narray="+new String(array)+", start="+a+", stop="+b;
			upper=(upper*10)+x;
		}
		
		long lower=0;
		int places=0;
		for(a++; a<b; a++){
			final char c=array.charAt(a);
			final int x=(c-z);
			assert(x<10 && x>=0) : x+" = "+(char)c+"\narray='"+new String(array)+"', start="+a+", stop="+b+", len="+array.length()+" -> '"+array.substring(a, b)+"'"+
				"\nThis function does not support exponents; if the input has an exponent, add the flag 'forceJavaParseDouble'.";
			lower=(lower*10)+x;
			places++;
		}
		
		double d=mult*(upper+lower*ByteBuilder.decimalInvMult[places]);
//		assert(d==parseDoubleSlow(array, a0, b)) : d+", "+parseDoubleSlow(array, a0, b);
		return d;
	}

	/**
	 * Parses an integer from a byte array starting at given position.
	 * @param array Byte array containing ASCII digits
	 * @param start Starting index
	 * @return Parsed integer value
	 */
	public static int parseInt(byte[] array, int start){
		return parseInt(array, start, array.length);
	}
	
//	/**
//	 * @param array Text
//	 * @param a Index of first digit
//	 * @param b Index after last digit (e.g., array.length)
//	 * @return Parsed number
//	 */
//	public static int parseInt(byte[] array, int a, int b){
//		assert(b>a);
//		int r=0;
//		final byte z='0';
//		int mult=1;
//		if(array[a]=='-'){mult=-1; a++;}
//		for(; a<b; a++){
//			int x=(array[a]-z);
//			assert(x<10 && x>=0) : x+" = "+(char)array[a]+"\narray="+new String(array)+", start="+a+", stop="+b;
//			r=(r*10)+x;
//		}
//		return r*mult;
//	}
	
	/** 
	 * Returns the int representation of a number represented in ASCII text, from position a to b.
	 * This function is much faster than creating a substring and calling Integer.parseInt()
	 * Throws Assertions rather than Exceptions for invalid input.
	 * This function does NOT detect overflows, e.g., values over 2^31-1 (Integer.MAX_VALUE).
	 * This function has no side-effects.
	 * @param array byte array containing the text to parse.
	 * @param a Index of the first digit of the number.
	 * @param b Index after the last digit (e.g., array.length).
	 * @return int representation of the parsed number.
	 * @throws Assertions rather than Exceptions for invalid input.
	 * 
	 * @TODO Correctly represent Integer.MIN_VALUE
	 * @TODO Detect overflow.
	 */
	public static int parseInt(byte[] array, int a, int b){
		assert(b>a) : "The start position of the text to parse must come before the stop position: "+
			a+","+b+","+new String(array);
		int r=0; //Initialize the return value to 0.

		//z holds the ASCII code for 0, which is subtracted from other ASCII codes
		//to yield the int value of a character.  For example, '7'-'0'=7,
		//because ASCII '7'=55, while ASCII '0'=48, and 55-48=7. 
		final byte z='0';

		//mult is 1 for positive numbers, or -1 for negative numbers.
		//It will be multiplied by the unsigned result to yield the final signed result.
		int mult=1;
		
		//If the term starts with a minus sign, set the multiplier to -1 and increment the position.
		if(array[a]=='-'){mult=-1; a++;}
		
		//Iterate through every position, incrementing a, up to b (exclusive).
		for(; a<b; a++){
			//x is the numeric value of the character at position a.
			//In other words, if array[a]='7',
			//x would be 7, not the ASCII code for '7' (which is 55).
			int x=(array[a]-z);
			
			//Assert that x is in the range of 0-9; otherwise, the character was not a digit.
			//The ASCII code will be printed here because in some cases the character could be
			//a control character (like carriage return or vertical tab or bell) which is unprintable.
			//But if possible the character will be printed to, as well as the position,
			//and the entire String from which the number is to be parsed.
			assert(x<10 && x>=0) : "Non-digit character with ASCII code "+(int)array[a]+" was encountered.\n"
					+"x="+x+"; char="+(char)array[a]+"\narray="+new String(array)+", start="+a+", stop="+b;
			
			//Multiply the old value by 10, then add the new 1's digit.
			//This is because the text is assumed to be base-10,
			//so each subsequent character will represent 1/10th the significance of the previous character.
			r=(r*10)+x;
		}
		
		//Change the unsigned value into a signed result, and return it.
		return r*mult;
	}
	
	/**
	 * @param array Text
	 * @param a Index of first digit
	 * @param b Index after last digit (e.g., array.length)
	 * @return Parsed number
	 */
	public static int parseInt(String array, int a, int b){
//		assert(false) : Character.toString(array.charAt(a));
		assert(b>a);
		int r=0;
		final byte z='0';
		int mult=1;
		if(array.charAt(a)=='-'){mult=-1; a++;}
		for(; a<b; a++){
			int x=(array.charAt(a)-z);
			assert(x<10 && x>=0) : x+" = "+array.charAt(a)+"\narray="+new String(array)+", start="+a+", stop="+b;
			r=(r*10)+x;
		}
		return r*mult;
	}
	
	/**
	 * @param array Text
	 * @param a Index of first digit
	 * @return Parsed number
	 */
	public static int parseInt(String array, int a){
		int r=0;
		final byte z='0';
		int mult=1;
		if(array.charAt(a)=='-'){mult=-1; a++;}
		for(; a<array.length(); a++){
			int x=(array.charAt(a)-z);
			if(x<0 || x>9) {break;}//End of number
			r=(r*10)+x;
		}
		return r*mult;
	}
	
	/**
	 * Parses a double from a string starting at given position.
	 * Automatically detects end of number.
	 *
	 * @param array String containing digits
	 * @param a Starting index
	 * @return Parsed double value
	 */
	public static double parseDouble(String array, int a){
		int b=a;
		while(b<array.length()) {
			char c=array.charAt(b);
			if((c<'0' || c>'9') && c!='.' && c!='-') {break;}
			b++;
		}
		return parseDouble(array, a, b);
	}
	
	/**
	 * Parses a float from a string starting at given position.
	 * @param array String containing digits
	 * @param a Starting index
	 * @return Parsed float value
	 */
	public static float parseFloat(String array, int a){
		return (float)parseDouble(array, a);
	}
	
	/**
	 * Parses a long from entire byte array.
	 * @param array Byte array containing ASCII digits
	 * @return Parsed long value
	 */
	public static long parseLong(byte[] array){return parseLong(array, 0, array.length);}
	
	/**
	 * Parses a long from byte array starting at given position.
	 * @param array Byte array containing ASCII digits
	 * @param start Starting index
	 * @return Parsed long value
	 */
	public static long parseLong(byte[] array, int start){return parseLong(array, start, array.length);}
	
	/**
	 * @param array Text
	 * @param a Index of first digit
	 * @param b Index after last digit (e.g., array.length)
	 * @return Parsed number
	 */
	public static long parseLong(byte[] array, int a, int b){
		assert(b>a);
		long r=0;
		final byte z='0';
		long mult=1;
		if(array[a]=='-'){mult=-1; a++;}
		for(; a<b; a++){
			int x=(array[a]-z);
			assert(x<10 && x>=0) : x+" = "+(char)array[a]+"\narray="+new String(array)+", start="+a+", stop="+b;
			r=(r*10)+x;
		}
		return r*mult;
	}
	
	/**
	 * Parses a long using 6-bit encoding with ASCII offset of 48.
	 * Each byte is treated as a 6-bit value after subtracting 48.
	 *
	 * @param array Byte array containing encoded values
	 * @param a Starting index (inclusive)
	 * @param b Ending index (exclusive)
	 * @return Parsed long value using 6-bit shifts
	 */
	public static long parseLongA48(byte[] array, int a, int b){
		if(array.length==0){return 0;}
		long x=0;
		for(; a<b; a++) {
			x<<=6;
			x|=(((long)array[a])-48);
		}
		return x;
	}
	
	/**
	 * @param array Text
	 * @param a Index of first digit
	 * @param b Index after last digit (e.g., array.length)
	 * @return Parsed number
	 */
	public static long parseLong(String array, int a, int b){
		assert(b>a);
		long r=0;
		final byte z='0';
		long mult=1;
		if(array.charAt(a)=='-'){mult=-1; a++;}
		for(; a<b; a++){
			int x=(array.charAt(a)-z);
			assert(x<10 && x>=0) : x+" = "+array.charAt(a)+"\narray="+new String(array)+", start="+a+", stop="+b;
			r=(r*10)+x;
		}
		return r*mult;
	}


	//Note: clen is optional, but allows poorly-formatted input like trailing whitespace
	//Without clen ",,," would become {0,0,0,0} 
	/**
	 * Parses a string containing comma-separated numbers into a long array.
	 * Handles malformed input like trailing whitespace.
	 * @param sub String containing comma-separated numbers
	 * @return Array of parsed long values, or null if input is null/empty
	 */
	public static long[] parseLongArray(String sub) {
		if(sub==null || sub.length()<1){return null;}
		long current=0;
//		int clen=0;
		LongList list=new LongList(min(8, 1+sub.length()/2));
		for(int i=0, len=sub.length(); i<len; i++){
//			System.err.println();
			int c=sub.charAt(i)-'0';
			if(c<0 || c>9){
//				System.err.println('A');
				//assert(clen>0);
				list.add(current);
				current=0;
//				clen=0;
			}else{
//				System.err.println('B');
				current=(current*10)+c;
//				clen++;
			}
//			System.err.println("i="+i+", c="+c+", current="+current+", list="+list);
		}
//		if(clen>0){
			list.add(current);
//		}
//		assert(false) : "\n'"+sub+"'\n"+Arrays.toString(list.toArray());
		return list.toArray();
	}
	
	/**
	 * Extracts ZMW (Zero Mode Waveguide) ID from PacBio read identifier.
	 * Expected format: m54283_190403_183820/4194374/919_2614
	 * @param id PacBio read identifier string
	 * @return ZMW ID number, or -1 if format is invalid
	 */
	public static int parseZmw(String id){
		//Example: m54283_190403_183820/4194374/919_2614
		//Run ID is m54283_190403_183820
		//zmw ID is 4194374.
		//Read start/stop coordinates are 919_2614
		int under=id.indexOf('_');
		int slash=id.indexOf('/');
		if(under<0 || slash<0){return -1;}
		String[] split=id.split("/");
		String z=split[1];
		return Integer.parseInt(z);
	}
	
	/**
	 * Splits a string on the first occurrence of a character.
	 * Returns array with one or two elements depending on whether character is found.
	 *
	 * @param s String to split
	 * @param c Character to split on
	 * @return Array containing parts before and after first occurrence of character
	 */
	public static String[] splitOnFirst(String s, char c) {
		final int idx=s.indexOf(c);
		if(idx<0) {return new String[] {s};}
		final String a=s.substring(0, idx);
		if(idx>=s.length()-1) {return new String[] {a};}
		return new String[] {a, s.substring(idx+1)};
	}
	
	/**
	 * Converts a symbol name to its character representation.
	 * Handles escape sequences by removing leading backslashes.
	 * @param b Symbol name string
	 * @return Character representation, or 0 if string is multi-character after processing
	 */
	public static char parseSymbolToCharacter(String b){
		b=parseSymbol(b);
		while(b.length()>1 && b.charAt(0)=='\\'){
			b=b.substring(1);
		}
//		System.err.println("Returning "+Character.toString(b.charAt(0)));
		return b.length()>1 ? 0 : b.charAt(0);
	}
	
	/**
	 * Converts symbolic names to their string representations.
	 * Handles convenience names like "space", "tab", "comma", etc.
	 * Also handles Java regex metacharacters with proper escaping.
	 *
	 * @param b Symbol name to convert
	 * @return String representation of the symbol, or original string if no match
	 */
	public static String parseSymbol(String b){
		if(b==null || b.length()<2){return b;}
		
		//Convenience characters
		if(b.equalsIgnoreCase("space")){
			return " ";
		}else if(b.equalsIgnoreCase("tab")){
			return "\t";
		}else if(b.equalsIgnoreCase("whitespace")){
			return "\\s+";
		}else if(b.equalsIgnoreCase("pound")){
			return "#";
		}else if(b.equalsIgnoreCase("greaterthan")){
			return ">";
		}else if(b.equalsIgnoreCase("lessthan")){
			return "<";
		}else if(b.equalsIgnoreCase("equals")){
			return "=";
		}else if(b.equalsIgnoreCase("colon")){
			return ":";
		}else if(b.equalsIgnoreCase("semicolon")){
			return ";";
		}else if(b.equalsIgnoreCase("comma")){
			return ",";
		}else if(b.equalsIgnoreCase("bang")){
			return "!";
		}else if(b.equalsIgnoreCase("and") || b.equalsIgnoreCase("ampersand")){
			return "&";
		}else if(b.equalsIgnoreCase("quote") || b.equalsIgnoreCase("doublequote")){
			return "\"";
		}else if(b.equalsIgnoreCase("singlequote") || b.equalsIgnoreCase("apostrophe")){
			return "'";
		}else if(b.equalsIgnoreCase("underscore")){
			return "_";
		}
		
		//Java meta characters
		if(b.equalsIgnoreCase("backslash")){
			return "\\\\";
		}else if(b.equalsIgnoreCase("hat") || b.equalsIgnoreCase("caret")){
			return "\\^";
		}else if(b.equalsIgnoreCase("dollar")){
			return "\\$";
		}else if(b.equalsIgnoreCase("dot")){
			return "\\.";
		}else if(b.equalsIgnoreCase("pipe") || b.equalsIgnoreCase("or")){
			return "\\|";
		}else if(b.equalsIgnoreCase("questionmark")){
			return "\\?";
		}else if(b.equalsIgnoreCase("star") || b.equalsIgnoreCase("asterisk")){
			return "\\*";
		}else if(b.equalsIgnoreCase("plus")){
			return "\\+";
		}else if(b.equalsIgnoreCase("openparen")){
			return "\\(";
		}else if(b.equalsIgnoreCase("closeparen")){
			return "\\)";
		}else if(b.equalsIgnoreCase("opensquare")){
			return "\\[";
		}else if(b.equalsIgnoreCase("closesquare")){
			return "]";
		}else if(b.equalsIgnoreCase("opencurly")){
			return "\\{";
		}else if(b.equalsIgnoreCase("closecurly")){
			return "}";
		}
		
		//No matches, return the literal
		return b;
	}
	
	/**
	 * Creates a character remapping array from a string specification.
	 * String must have even length with character pairs (from->to).
	 *
	 * @param b Remap specification string, or "f"/"false" to disable
	 * @return Byte array for character remapping, or null if disabled
	 * @throws AssertionError if string length is odd
	 */
	public static byte[] parseRemap(String b){
		final byte[] remap;
		if(b==null || ("f".equalsIgnoreCase(b) || "false".equalsIgnoreCase(b))){
			remap=null;
		}else{
			assert((b.length()&1)==0) : "Length of remap argument must be even.  No whitespace is allowed.";
			
			remap=new byte[128];
			for(int j=0; j<remap.length; j++){remap[j]=(byte)j;}
			for(int j=0; j<b.length(); j+=2){
				char x=b.charAt(j), y=b.charAt(j+1);
				remap[x]=(byte)y;
			}
		}
		return remap;
	}
	
	/** Returns the smaller of two integers */
	public static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger of two integers */
	public static final int max(int x, int y){return x>y ? x : y;}

}
