package json;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

/**
 * Reusable UTF-8 parser for JSON objects and arrays, with historical token extensions.
 * Parsed integers use Long and fractional/exponent values use finite Double.
 * Structural errors set an observable flag and may return partial data; check
 * errorState() or validate() before trusting input from an external source.
 * Each thread must own its parser. Byte-array input is borrowed and must remain
 * unchanged while parsing; returned objects and strings belong to the caller.
 * @author Brian Bushnell
 */
public class JsonParser{
	
	/*--------------------------------------------------------------*/
	/*----------------          Entry Point         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Test method demonstrating JSON parsing capabilities.
	 * Prints a sample containing a deliberately invalid bare token, its partial
	 * parse, and a separate valid array. This is a demonstration, not a test oracle.
	 * @param args Command-line arguments (unused)
	 */
	public static void main(String[] args){
		String s;
		
		s="{\n"
			+"   \"33154\": {\n"
			+"      \"name\": \"Opisthokonta\",\n"
			+"      \"tax_id\": 33154,\n"
			+"      \"level\": \"no rank\",\n"
			+"      \"no rank\": {\n"
			+"         \"name\": \"Opisthokonta\",\n"
			+"         \"tax_id\": 33154\n"
			+"      },\n"
			+"      \"foo\": {\n"
			+"         \"bar\": \"bam\",\n"
			+"         \"sam\": \"cram\"\n"
			+"      },\n"
			+"      \"foo2\": {\n"
			+"         \"true\": false\n"
			+"      },\n"
			+"      \"foo3\": {\n"
			+"         \"null\": null\n"
			+"      },\n"
			+"      \"foo4\": {\n"
			+"         \"null\": invalid\n"
			+"      },\n"
			+"      \"superkingdom\": {\n"
			+"         \"name\": \"Eukaryota\",\n"
			+"         \"tax_id\": 2759,\n"
			+"         \"number1\": 2759,\n"
			+"         \"number2\": -2759,\n"
			+"         \"number3\": .2759,\n"
			+"         \"number4\": 2.759,\n"
			+"         \"number5\": -2.759,\n"
			+"         \"number6\": -2.759e17,\n"
			+"         \"number7\": -2.759e-1,\n"
			+"         \"number8\": -2.759E-1,\n"
			+"         \"number9\": -2E-1,\n"
			+"         \"slash\": \"hello \\\"world\\\"\",\n"
			+"         \"slash\": \"hello world\",\n"
			+"         \"complex\": [\"hello world\", 1, {\"tax_id\": 2759}, [3, 4, 5]]\n"
			+"      }\n"
			+"   }\n"
			+"}";
		
		System.out.println("Original:\n"+s);
		JsonParser jp=new JsonParser(s);
		JsonObject j=jp.parseJsonObject();
		System.out.println("Original:\n"+s);
		System.out.println("Regenerated:\n"+j);
		
		s="[\"complex\", 1, {\"b\": 2}, [3, 4, 5]]";
		
		System.out.println("Original:\n"+s);
		jp.set(s.getBytes(StandardCharsets.UTF_8));
		Object[] array=jp.parseJsonArray();
		System.out.println("Original:\n"+s);
		System.out.println("Regenerated:\n"+JsonObject.toString(array));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a parser with no input; set(byte[]) supplies or replaces it. */
	public JsonParser(){}
	
	/** Copies String input to UTF-8 bytes.
	 * @param s Non-null input text */
	public JsonParser(String s){
		set(s.getBytes(StandardCharsets.UTF_8));
	}
	
	/** Borrows an input buffer without copying it.
	 * @param s UTF-8 bytes, or null for no input */
	public JsonParser(byte[] s){
		set(s);
	}
	
	/** Parses an object using a new parser; null/empty input returns null.
	 * Use an instance when the caller needs its error flag or explicit validation.
	 * @param s Input text, or null
	 * @return Parsed object, possibly partial on malformed input, or null */
	public static JsonObject parseJsonObjectStatic(String s){
		if(s==null || s.length()<1){return null;}
		return new JsonParser(s).parseJsonObject();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Replaces borrowed input and clears position and error state.
	 * @param s New UTF-8 input, or null
	 * @return This parser for chaining */
	public JsonParser set(byte[] s){
		text=s;
		pos=0;
		errorState=false;
		return this;
	}
	
	/** Resets from a String and parses an object.
	 * @param s Non-null UTF-8-compatible text
	 * @return Object or partial object; null if input is empty or not an object */
	public JsonObject parseJsonObject(String s){
		set(s.getBytes(StandardCharsets.UTF_8));
		return parseJsonObject();
	}
	
	/** Resets from borrowed bytes and parses an object.
	 * @param s UTF-8 input, or null
	 * @return Object or partial object; null if input is absent or not an object */
	public JsonObject parseJsonObject(byte[] s){
		set(s);
		return parseJsonObject();
	}
	
	/** Resets from a String and parses an array.
	 * @param s Non-null UTF-8-compatible text
	 * @return Values or partial values; null if input is empty or not an array */
	public Object[] parseJsonArray(String s){
		set(s.getBytes(StandardCharsets.UTF_8));
		return parseJsonArray();
	}
	
	/** Resets from borrowed bytes and parses an array.
	 * @param s UTF-8 input, or null
	 * @return Values or partial values; null if input is absent or not an array */
	public Object[] parseJsonArray(byte[] s){
		set(s);
		return parseJsonArray();
	}
	
	/** Parses an object at the current position after optional whitespace.
	 * Leaves position on the closing brace; use set to restart parsing. Does not
	 * reject trailing input; validate checks a complete document from its start.
	 * @return Object or partial object; null if input is absent or has the wrong root */
	public JsonObject parseJsonObject(){
		if(text==null || text.length<1){return null;}
		// Fixed [json/JsonParser#001]: skip leading whitespace and flag a wrong root,
		// instead of depending on assertions to reject non-object input.
		while(pos<text.length && isWhitespace(text[pos])){pos++;}
		if(pos>=text.length || text[pos]!='{'){errorState=true; return null;}
		JsonObject o=makeObject();
		return o;
	}
	
	/** Parses an array at the current position after optional whitespace.
	 * Leaves position on the closing bracket; use set to restart parsing.
	 * @return Values or partial values; null if input is absent or has the wrong root */
	public Object[] parseJsonArray(){
		if(text==null || text.length<1){return null;}
		// [json/JsonParser#001] FIXED: skip leading whitespace; fail gracefully (errorState) on non-array input
		while(pos<text.length && isWhitespace(text[pos])){pos++;}
		if(pos>=text.length || text[pos]!='['){errorState=true; return null;}
		Object[] array=makeArray();
		return array;
	}
	
	/** Checks one complete value from the start of the saved input without consuming it.
	 * Whitespace is allowed around the value; empty input and trailing data are invalid.
	 * Requires standard JSON literals, number syntax and string escapes, unlike
	 * the permissive parsing methods. Long/finite-Double representation limits
	 * still apply. Updates the error flag; use set to clear it before another parse.
	 * @return True when a complete value was read without syntax or conversion errors */
	public boolean validate(){
		final int savedPos=pos;
		pos=0;
		errorState=false;
		strictValidation=true;
		try{
			if(text==null || text.length==0){errorState=true; return false;}
			skipWhitespace();
			readValue();
			if(!errorState){
				pos++;
				skipWhitespace();
				if(pos!=text.length){errorState=true;}
			}
			return !errorState;
		}catch(RuntimeException e){errorState=true; return false;}
		finally{pos=savedPos; strictValidation=false;}
	}

	/** Reports whether the most recent parse or validation encountered an error.
	 * @return Current error flag; set(byte[]) clears it */
	public boolean errorState(){return errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Converts a bare token into a typed value: null, Boolean (via
	 * parseBoolean, accepting t/f/true/false), Double (if it contains '.'/'e'/'E'),
	 * else Long. On a parse failure sets errorState and returns the raw string.
	 * @param s Token text without surrounding whitespace
	 * @return The parsed value (Long/Double/Boolean/null), or the raw String on error
	 */
	private Object tokenToObject(String s){
		if(strictValidation && !isStandardToken(s)){errorState=true; return s;}
		final char firstLetter=s.length()>0 ? s.charAt(0) : 0;
		Object value;
		try{
			if(Character.isLetter(firstLetter)){
				if(verbose){outstream.println("Letter");}
				if(s.equalsIgnoreCase("null")){
					value=null;
				}else{
					value=parseBoolean(s);
				}
			}else{
				if(verbose){outstream.println("Number");}
				if(s.indexOf('.')>=0 || s.indexOf('e')>=0 || s.indexOf('E')>=0){
					final double number=Double.parseDouble(s);
					if(!Double.isFinite(number)){throw new NumberFormatException("Non-finite JSON number: "+s);}
					value=number;
				}else{
					value=Long.parseLong(s);
				}
			}
		}catch(Exception e){
			//This handles an incorrectly formatted input file
			errorState=true;
			value=s;
		}
		return value;
	}
	
	/** Checks standard literal or number syntax without converting the token.
	 * Numbers require an integer part, no leading zeros, and digits after a dot
	 * or exponent marker. Conversion limits are enforced separately.
	 * @param s Non-null bare token
	 * @return True for lowercase true/false/null or a JSON number */
	private static boolean isStandardToken(final String s){
		if(s.equals("true") || s.equals("false") || s.equals("null")){return true;}
		final int length=s.length();
		int i=0;
		if(i<length && s.charAt(i)=='-'){i++;}
		if(i>=length){return false;}
		final char first=s.charAt(i++);
		if(first>='1' && first<='9'){
			while(i<length && isDigit(s.charAt(i))){i++;}
		}else if(first!='0'){return false;}
		if(i<length && s.charAt(i)=='.'){
			final int start=++i;
			while(i<length && isDigit(s.charAt(i))){i++;}
			if(i==start){return false;}
		}
		if(i<length && (s.charAt(i)=='e' || s.charAt(i)=='E')){
			i++;
			if(i<length && (s.charAt(i)=='+' || s.charAt(i)=='-')){i++;}
			final int start=i;
			while(i<length && isDigit(s.charAt(i))){i++;}
			if(i==start){return false;}
		}
		return i==length;
	}

	/** Recognizes an ASCII decimal digit.
	 * @param c Candidate character
	 * @return True for 0 through 9 */
	private static boolean isDigit(final char c){return c>='0' && c<='9';}

	/** Recognizes the historical case-insensitive true/false and t/f tokens.
	 * @param s Candidate token
	 * @return Parsed boolean
	 * @throws Exception If the token is not recognized */
	private static boolean parseBoolean(String s) throws Exception{
		if(s==null){throw INVALID_JSON;}
		if(s.equalsIgnoreCase("true") || s.equalsIgnoreCase("t")){return true;}
		if(s.equalsIgnoreCase("false") || s.equalsIgnoreCase("f")){return false;}
		throw INVALID_JSON;
	}

	/** Recognizes the four JSON whitespace bytes.
	 * @param b Candidate byte
	 * @return True for space, tab, newline or carriage return */
	private static boolean isWhitespace(byte b){return b==' ' || b=='\t' || b=='\n' || b=='\r';}

	/**
	 * Decodes a single-character JSON escape (the byte after a backslash) to its literal
	 * byte. Recognizes " \\ / n r t b f; parsing returns an unknown escape character
	 * unchanged after dropping its backslash, while validation also records an error.
	 * Unicode escapes are handled separately by readQuotedString.
	 * @param b The byte following the backslash
	 * @return The decoded byte
	 */
	private byte unescape(byte b){
		switch(b){
			case '"': return '"';
			case '\\': return '\\';
			case '/': return '/';
			case 'n': return '\n';
			case 'r': return '\r';
			case 't': return '\t';
			case 'b': return '\b';
			case 'f': return '\f';
			default: if(strictValidation){errorState=true;} return b;
		}
	}

	/** Reads quoted UTF-8 text, retaining escaped UTF-16 units without byte truncation.
	 * Leaves pos at the closing quote, or at end of input on an unterminated string.
	 * Unknown single-character escapes are accepted only outside validation.
	 * @return Decoded text; malformed input sets errorState and may yield partial text */
	private String readQuotedString(){
		assert(text[pos]=='"') : "Quoted-string decoding requires the opening quote at pos="+pos;
		int start=++pos;
		StringBuilder out=null;
		for(; pos<text.length; pos++){
			final byte b=text[pos];
			if(b=='"'){
				final String tail=decodeUtf8(start, pos);
				return out==null ? tail : out.append(tail).toString();
			}
			if(b=='\\'){
				if(out==null){out=new StringBuilder();}
				out.append(decodeUtf8(start, pos));
				if(++pos>=text.length){errorState=true; return out.toString();}
				if(text[pos]=='u'){
					int value=-1;
					if(pos+4<text.length){
						final int a=hexDigit(text[pos+1]), c=hexDigit(text[pos+2]);
						final int d=hexDigit(text[pos+3]), e=hexDigit(text[pos+4]);
						if((a|c|d|e)>=0){value=(a<<12)|(c<<8)|(d<<4)|e;}
					}
					if(value>=0){out.append((char)value); pos+=4;}
					else{errorState=true; out.append('u');}
				}else{out.append((char)(unescape(text[pos])&255));}
				start=pos+1;
			}else if(b>=0 && b<0x20){errorState=true;}
		}
		errorState=true;
		final String tail=decodeUtf8(start, pos);
		return out==null ? tail : out.append(tail).toString();
	}

	/** Decodes a raw span and checks ambiguous replacement characters by byte round trip.
	 * A replacement character triggers a byte round trip, distinguishing valid U+FFFD
	 * text from replacement of malformed input. Only that uncommon case allocates bytes.
	 * @param from First byte, inclusive
	 * @param to First byte after the span
	 * @return Decoded span; malformed bytes set errorState */
	private String decodeUtf8(final int from, final int to){
		final String decoded=new String(text, from, to-from, StandardCharsets.UTF_8);
		if(decoded.indexOf('\ufffd')>=0){
			final byte[] encoded=decoded.getBytes(StandardCharsets.UTF_8);
			if(encoded.length!=to-from){errorState=true;}
			else{
				for(int i=0; i<encoded.length; i++){
					if(encoded[i]!=text[from+i]){errorState=true; break;}
				}
			}
		}
		return decoded;
	}

	/** Decodes one hexadecimal digit.
	 * @param b Candidate ASCII byte
	 * @return Value 0-15, or -1 if not a hexadecimal digit */
	private static int hexDigit(byte b){
		if(b>='0' && b<='9'){return b-'0';}
		if(b>='a' && b<='f'){return b-'a'+10;}
		if(b>='A' && b<='F'){return b-'A'+10;}
		return -1;
	}
	
	/** Reads an object and leaves pos on its closing brace.
	 * Preserves the closing-delimiter convention of [json/JsonParser#004].
	 * Missing separators, keys or terminators set errorState rather than relying
	 * on assertions, continuing the error handling of #001 and #003.
	 * @return Parsed members, possibly partial on malformed input */
	private JsonObject makeObject(){
		assert(text[pos]=='{') : "Object parser requires its opening brace at pos="+pos;
		final JsonObject current=new JsonObject();
		pos++;
		skipWhitespace();
		if(pos<text.length && text[pos]=='}'){return current;}
		while(pos<text.length){
			if(text[pos]!='"'){errorState=true; return current;}
			final String key=readQuotedString();
			if(errorState){return current;}
			pos++;
			skipWhitespace();
			if(pos>=text.length || text[pos]!=':'){errorState=true; return current;}
			pos++;
			skipWhitespace();
			final Object value=readValue();
			current.add(key, value);
			if(errorState){return current;}
			pos++;
			skipWhitespace();
			if(pos<text.length && text[pos]=='}'){return current;}
			if(pos>=text.length || text[pos]!=','){errorState=true; return current;}
			pos++;
			skipWhitespace();
		}
		errorState=true;
		return current;
	}

	/** Reads an array and leaves pos on its closing bracket.
	 * Elements retain their parsed types; structural failures set errorState.
	 * @return Parsed values, possibly partial on malformed input */
	private Object[] makeArray(){
		assert(text[pos]=='[') : "Array parser requires its opening bracket at pos="+pos;
		final ArrayList<Object> current=new ArrayList<Object>();
		pos++;
		skipWhitespace();
		if(pos<text.length && text[pos]==']'){return current.toArray();}
		while(pos<text.length){
			current.add(readValue());
			if(errorState){return current.toArray();}
			pos++;
			skipWhitespace();
			if(pos<text.length && text[pos]==']'){return current.toArray();}
			if(pos>=text.length || text[pos]!=','){errorState=true; return current.toArray();}
			pos++;
			skipWhitespace();
		}
		errorState=true;
		return current.toArray();
	}

	/** Reads one value, leaving pos on its final byte or closing delimiter.
	 * Bare tokens retain the historical token extensions described by tokenToObject.
	 * @return Parsed value; malformed values may return null or their raw token */
	private Object readValue(){
		if(pos>=text.length){errorState=true; return null;}
		final byte first=text[pos];
		if(first=='"'){return readQuotedString();}
		if(first=='{'){return makeObject();}
		if(first=='['){return makeArray();}
		final int start=pos;
		while(pos<text.length){
			final byte b=text[pos];
			if(isWhitespace(b) || b==',' || b==']' || b=='}'){break;}
			pos++;
		}
		if(pos==start){errorState=true; return null;}
		final String token=new String(text, start, pos-start, StandardCharsets.US_ASCII);
		pos--;
		return tokenToObject(token);
	}

	/** Skips only the four JSON whitespace bytes, without crossing end of input. */
	private void skipWhitespace(){
		while(pos<text.length && isWhitespace(text[pos])){pos++;}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Borrowed UTF-8 input, or null before input is supplied. */
	byte[] text;
	/** Current byte index; completed container parses leave it on the closing delimiter. */
	int pos=0;
	/** Records syntax, encoding or unsupported numeric-range errors. */
	boolean errorState;
	/** Enabled only inside validate; finally restores permissive parsing mode. */
	private boolean strictValidation;
	
	/** Compile-time token tracing switch. */
	private static final boolean verbose=false;
	/** Destination for optional token tracing. */
	private static final PrintStream outstream=System.err;
	/** Reused internal token-rejection marker, caught without exposing its stack trace. */
	private static final Exception INVALID_JSON=new Exception("Invalid Json");
	
}
