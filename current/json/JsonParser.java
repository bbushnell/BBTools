package json;

import java.io.PrintStream;
import java.util.ArrayList;

import structures.ByteBuilder;

/**
 * Simple JSON parser for reading JSON objects and arrays from text.
 * Parses JSON text into JsonObject instances or Object arrays.
 * Thread-safe when each thread uses its own instance.
 * @author Brian Bushnell
 */
public class JsonParser {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Test method demonstrating JSON parsing capabilities.
	 * Parses sample JSON objects and arrays, printing original and regenerated forms.
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
		
//		s="{\"complex\": [\"a\", 1, {\"b\": 2}, [3, 4, 5]]\n}";
		
		System.out.println("Original:\n"+s);
		JsonParser jp=new JsonParser(s);
		JsonObject j=jp.parseJsonObject();
		System.out.println("Original:\n"+s);
		System.out.println("Regenerated:\n"+j);
		
		s="[\"complex\", 1, {\"b\": 2}, [3, 4, 5]]";
		
		System.out.println("Original:\n"+s);
		jp.set(s.getBytes());
		Object[] array=jp.parseJsonArray();
		System.out.println("Original:\n"+s);
		System.out.println("Regenerated:\n"+JsonObject.toString(array));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public JsonParser(){}
	
	public JsonParser(String s){
		set(s.getBytes());
	}
	
	public JsonParser(byte[] s){
		set(s);
	}
	
	public static JsonObject parseJsonObjectStatic(String s){
		if(s==null || s.length()<1) {return null;}
		return new JsonParser(s).parseJsonObject();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public JsonParser set(byte[] s){
		text=s;
		pos=0;
		errorState=false;
		return this;
	}
	
	public JsonObject parseJsonObject(String s){
		set(s.getBytes());
		return parseJsonObject();
	}
	
	public JsonObject parseJsonObject(byte[] s){
		set(s);
		return parseJsonObject();
	}
	
	public Object[] parseJsonArray(String s){
		set(s.getBytes());
		return parseJsonArray();
	}
	
	public Object[] parseJsonArray(byte[] s){
		set(s);
		return parseJsonArray();
	}
	
	public JsonObject parseJsonObject(){
		if(text==null || text.length<1){return null;}
		// [json/JsonParser#001] FIXED: skip leading whitespace; fail gracefully (errorState) on non-object input instead of asserting under -ea
		while(pos<text.length && isWhitespace(text[pos])){pos++;}
		if(pos>=text.length || text[pos]!='{'){errorState=true; return null;}
		JsonObject o=makeObject();
		return o;
	}
	
	public Object[] parseJsonArray(){
		if(text==null || text.length<1){return null;}
		// [json/JsonParser#001] FIXED: skip leading whitespace; fail gracefully (errorState) on non-array input
		while(pos<text.length && isWhitespace(text[pos])){pos++;}
		if(pos>=text.length || text[pos]!='['){errorState=true; return null;}
		Object[] array=makeArray();
		return array;
	}
	
	public boolean validate(){
		if(text==null || text.length<1){return true;}
		try {
			if(text[0]=='['){
				Object[] array=parseJsonArray();
				return !errorState;
			}else if(text[0]=='{'){
				JsonObject o=parseJsonObject();
				return !errorState;
			}
		} catch (Throwable e) {}
		return false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Converts an accumulated token buffer into a typed value: null, Boolean (via
	 * parseBoolean, accepting t/f/true/false), Double (if it contains '.'/'e'/'E'),
	 * else Long. On a parse failure sets errorState and returns the raw string.
	 * @param bb The token buffer (cleared by this call)
	 * @return The parsed value (Long/Double/Boolean/null), or the raw String on error
	 */
	private Object bufferToObject(ByteBuilder bb){
		String s=bb.toString();
		bb.clear();
		final char firstLetter=s.length()>0 ? s.charAt(0) : 0;
		Object value;
		try {
			if(Character.isLetter(firstLetter)){
				if(verbose){outstream.println("Letter");}
				if(s.equalsIgnoreCase("null")){
					value=null;
				}else{
//					value=Boolean.parseBoolean(s);
					value=parseBoolean(s);
				}
			}else{
				if(verbose){outstream.println("Number");}
				if(s.indexOf('.')>=0 || s.indexOf('e')>=0 || s.indexOf('E')>=0){
					value=Double.parseDouble(s);
				}else{
					value=Long.parseLong(s);
				}
			}
		} catch (Exception e) {
			//This handles an incorrectly formatted input file
			errorState=true;
			value=s;
		}
		return value;
	}
	
	private static boolean parseBoolean(String s) throws Exception{
		if(s==null){throw INVALID_JSON;}
		if(s.equalsIgnoreCase("true") || s.equalsIgnoreCase("t")){return true;}
		if(s.equalsIgnoreCase("false") || s.equalsIgnoreCase("f")){return false;}
		throw INVALID_JSON;
	}

	/** True for the four JSON whitespace bytes (space, tab, newline, carriage return). */
	private static boolean isWhitespace(byte b){return b==' ' || b=='\t' || b=='\n' || b=='\r';}

	/**
	 * Decodes a single-character JSON escape (the byte after a backslash) to its literal
	 * byte. Recognizes " \\ / n r t b f; an unrecognized escape is returned as-is (lenient).
	 * \\uXXXX is handled separately by appendUnicodeEscape.
	 * @param b The byte following the backslash
	 * @return The decoded byte
	 */
	private static byte unescape(byte b){
		switch(b){
			case '"': return '"';
			case '\\': return '\\';
			case '/': return '/';
			case 'n': return '\n';
			case 'r': return '\r';
			case 't': return '\t';
			case 'b': return '\b';
			case 'f': return '\f';
			default: return b;
		}
	}

	/**
	 * Decodes a \\uXXXX escape. On entry pos is at 'u'; reads the four hex digits at
	 * pos+1..pos+4, appends the decoded character, and returns the index of the last hex
	 * digit (so the caller's pos++ advances past it). On malformed hex, sets errorState,
	 * appends 'u', and returns pos unchanged.
	 * @param bb Buffer to append the decoded character to
	 * @param pos Index of the 'u'
	 * @return New pos (index of the last consumed hex digit), or pos on error
	 */
	private int appendUnicodeEscape(ByteBuilder bb, int pos){
		if(pos+4<text.length){
			final int d1=hexDigit(text[pos+1]), d2=hexDigit(text[pos+2]);
			final int d3=hexDigit(text[pos+3]), d4=hexDigit(text[pos+4]);
			if((d1|d2|d3|d4)>=0){
				bb.append((char)((d1<<12)|(d2<<8)|(d3<<4)|d4));
				return pos+4;
			}
		}
		errorState=true;
		bb.append('u');
		return pos;
	}

	/** Returns the value 0-15 of a hex-digit byte, or -1 if not a hex digit. */
	private static int hexDigit(byte b){
		if(b>='0' && b<='9'){return b-'0';}
		if(b>='a' && b<='f'){return b-'a'+10;}
		if(b>='A' && b<='F'){return b-'A'+10;}
		return -1;
	}
	
	/**
	 * Parses one JSON object beginning at the current '{' and advances pos past the
	 * matching '}'. Recurses into nested objects (makeObject) and arrays (makeArray).
	 * Both value and structural parse errors are graceful: errorState is set (see #001).
	 * @return The parsed JsonObject. The outermost object is returned directly; nested
	 *         objects are attached under their key in the enclosing object.
	 */
	private JsonObject makeObject(){
		assert(text[pos]=='{');
		pos++;
		
		if(verbose){outstream.println("Entering makeObject.");}
		
		JsonObject current=new JsonObject();
		ByteBuilder bb=new ByteBuilder();
		boolean quoteMode=false;
		boolean slashMode=false;
		String key=null;
		
		for(; pos<text.length; pos++){
			final byte b=text[pos];
//			if(verbose){outstream.println(pos+"=\t"+(char)b);
			
			if(quoteMode){
				// [json/JsonParser#002] FIXED: decode escapes (inverse of JsonObject.appendQuoted)
				if(slashMode){
					slashMode=false;
					if(b=='u'){pos=appendUnicodeEscape(bb, pos);}
					else{bb.append(unescape(b));}
				}else if(b=='\\'){
					slashMode=true;//consume backslash, do not append
				}else if(b=='"'){
					if(verbose){outstream.println(">Quote; quote mode="+quoteMode+", key="+key+", buffer="+bb);}
					String s=bb.toString();
					bb.clear();
					if(key==null){
						key=s;
						if(verbose){outstream.println("Set key to \""+key+"\"");}
					}else{
						current.add(key, s);
						if(verbose){outstream.println("Added \""+key+"\": \""+s+"\"");}
						key=null;
					}
					quoteMode=!quoteMode;
				}else{
					bb.append(b);
				}
			}else if(b=='"'){
				if(verbose){outstream.println(">Quote; quote mode="+quoteMode+", key="+key+", buffer="+bb);}
				quoteMode=!quoteMode;
			}else if(b==','){
				if(verbose){outstream.println(">Comma; key="+key+", buffer=\""+bb+"\""/*+"\n"+new String(text, 0, pos)*/);}
				if(key!=null){//number or boolean
					final Object value=bufferToObject(bb);
					current.add(key, value);
					key=null;
					if(verbose){outstream.println("Added "+value+"; current=\n"+current+"\n");}
				}
			}else if(b==':'){
				if(verbose){outstream.println(">Colon");}
				if(key==null){errorState=true;}// [json/JsonParser#001] FIXED: was assert(key!=null)
			}else if(b=='{'){
				if(verbose){outstream.println(">{, key="+key+", A) current object is:\n"+current+"\n");}
				JsonObject j=makeObject();
				if(key==null){//outermost?
					if(verbose){outstream.println("Returning.");}
					return j;
				}else{
					current.add(key, j);
					if(verbose){outstream.println("Added new object:\n"+j+"\n");}
					key=null;
				}
			}else if(b=='}'){
				if(verbose){outstream.println(">}, key="+key+", B) current object is:\n"+current+"\n");}
				if(key!=null){//number or boolean
					final Object value=bufferToObject(bb);
					current.add(key, value);
					key=null;
					if(verbose){outstream.println("Added "+value+"; current=\n"+current+"\n");}
				}
				// [json/JsonParser#004] FIXED: removed pos++ so makeObject leaves pos AT '}', symmetric with makeArray's ']' case (which has no pos++); the caller's loop pos++ then advances uniformly. The old double-advance skipped the char after a nested object, crashing (assert(false)) when an object was the last element of an array, and mis-exiting object-last-in-object via the fallthrough.
				return current;
			}else if(b=='['){
				if(verbose){outstream.println(">[, C) current object is:\n"+current+"\n");}
				Object[] array=makeArray();
				if(key==null || bb.length()!=0){errorState=true;}// [json/JsonParser#001] FIXED: was asserts
				current.add(key, array);
				key=null;
			}else if(b==']'){
				if(verbose){outstream.println(">], D) current object is:\n"+current+"\n");}
				errorState=true;// [json/JsonParser#001] FIXED: was assert(false) (']' inside object = malformed)
			}else if(b==' ' || b=='\t' || b=='\n' || b=='\r'){
				if(verbose){outstream.println(">Other");}
				//ignore
			}else{
				if(verbose){outstream.println(">NormalMode, buffer="+bb);}
				bb.append(b);
			}
		}
		errorState=true;// [json/JsonParser#003] FIXED: reached only on unterminated input (no closing '}')
		return current;
	}

	/**
	 * Parses one JSON array beginning at the current '[' and advances pos past the
	 * matching ']'. Recurses into nested objects/arrays.
	 * @return The array elements as an Object[] (String, Long, Double, Boolean, null,
	 *         JsonObject, or nested Object[]).
	 */
	private Object[] makeArray(){
		assert(text[pos]=='[');
		pos++;
		
		if(verbose){outstream.println("Entering makeArray.");}
		
		ArrayList<Object> current=new ArrayList<Object>();
		ByteBuilder bb=new ByteBuilder();
		boolean quoteMode=false;
		boolean slashMode=false;
		
		for(; pos<text.length; pos++){
			final byte b=text[pos];
			
			if(quoteMode){
				// [json/JsonParser#002] FIXED: decode escapes (inverse of JsonObject.appendQuoted)
				if(slashMode){
					slashMode=false;
					if(b=='u'){pos=appendUnicodeEscape(bb, pos);}
					else{bb.append(unescape(b));}
				}else if(b=='\\'){
					slashMode=true;//consume backslash, do not append
				}else if(b=='"'){
					if(verbose){outstream.println(">Quote; quote mode="+quoteMode+", buffer="+bb);}
					String s=bb.toString();
					bb.clear();
					current.add(s);
					quoteMode=!quoteMode;
				}else{
					bb.append(b);
				}
			}else if(b=='"'){
				if(verbose){outstream.println(">Quote; quote mode="+quoteMode+", buffer="+bb);}
				quoteMode=!quoteMode;
			}else if(b==','){
				if(verbose){outstream.println(">Comma; buffer=\""+bb+"\"");}
				if(bb.length()>0){
					final Object value=bufferToObject(bb);
					current.add(value);
					if(verbose){outstream.println("Added "+value+"; current=\n"+current+"\n");}
				}
			}else if(b==':'){
				if(verbose){outstream.println(">Colon");}
				errorState=true;// [json/JsonParser#001] FIXED: was assert(false) (':' inside array = malformed)
			}else if(b=='{'){
				if(verbose){outstream.println(">{, E) current object is:\n"+current+"\n");}
				JsonObject j=makeObject();
				current.add(j);
				if(verbose){outstream.println("Added new object:\n"+j+"\n");}
			}else if(b=='}'){
				if(verbose){outstream.println(">}, current array is:\n"+current+"\n");}
				errorState=true;// [json/JsonParser#001] FIXED: was assert(false) ('}' inside array = malformed)
			}else if(b=='['){
				if(verbose){outstream.println(">[, F) current object is:\n"+current+"\n");}
				Object[] array=makeArray();
				current.add(array);
			}else if(b==']'){
				if(verbose){outstream.println(">], G) current object is:\n"+current+"\n");}
				if(bb.length()>0){
					final Object value=bufferToObject(bb);
					current.add(value);
					if(verbose){outstream.println("Added "+value+"; current=\n"+current+"\n");}
				}
				if(verbose){outstream.println("Returning "+current+"; text="+new String(text, 0, pos));}
				return current.toArray(); 
			}else if(b==' ' || b=='\t' || b=='\n' || b=='\r'){
				if(verbose){outstream.println(">Other");}
				//ignore
			}else{
				if(verbose){outstream.println(">NormalMode, buffer="+bb);}
				bb.append(b);
			}
		}
		errorState=true;// [json/JsonParser#003] FIXED: reached only on unterminated input (no closing ']')
		return current.toArray();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	byte[] text;
	int pos=0;
	boolean errorState;
	
	private static final boolean verbose=false;
	private static final PrintStream outstream=System.err;
	private static final Exception INVALID_JSON=new Exception("Invalid Json");
	
}
