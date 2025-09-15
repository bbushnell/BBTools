package json;

import java.io.PrintStream;
import java.util.ArrayList;

import structures.ByteBuilder;

/**
 * How to use this class:
 * 1) Create one instance per thread 
 * 2) set() some Json text
 * 3) Call either parseJsonObject or parseJsonArray
 * @author Brian Bushnell
 *
 */
public class JsonParser {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** For testing */
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
		assert(text[0]=='{') : text[0]+"\n"+new String(text);
		JsonObject o=makeObject();
		return o;
	}
	
	public Object[] parseJsonArray(){
		if(text==null || text.length<1){return null;}
		assert(text[0]=='[') : text[0]+"\n"+new String(text);
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
	
	/** This handles cases when the value is not surrounded by quotes. */
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
	
	/** Not strictly correct, but allows t and f */
	private static boolean parseBoolean(String s) throws Exception{
		if(s==null){throw INVALID_JSON;}
		if(s.equalsIgnoreCase("true") || s.equalsIgnoreCase("t")){return true;}
		if(s.equalsIgnoreCase("false") || s.equalsIgnoreCase("f")){return false;}
		throw INVALID_JSON;
	}
	
	/** Create a JsonObject from { to the next } */
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
				if(slashMode){
					if(verbose){outstream.println(">SlashModeEnd, buffer="+bb);}
					bb.append(b);
					slashMode=false;
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
					if(verbose){outstream.println(">QuoteMode, buffer="+bb);}
					if(b=='\\'){
						if(verbose){outstream.println(">SlashMode, buffer="+bb);}
						slashMode=true;
					}
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
				assert(key!=null);
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
				pos++;
				return current; 
			}else if(b=='['){
				if(verbose){outstream.println(">[, C) current object is:\n"+current+"\n");}
				Object[] array=makeArray();
				assert(key!=null) : "Should be in makeArray.";
				current.add(key, array);
				key=null;
				assert(bb.length()==0);
			}else if(b==']'){
				if(verbose){outstream.println(">], D) current object is:\n"+current+"\n");}
				assert(false);
			}else if(b==' ' || b=='\t' || b=='\n' || b=='\r'){
				if(verbose){outstream.println(">Other");}
				//ignore
			}else{
				if(verbose){outstream.println(">NormalMode, buffer="+bb);}
				bb.append(b);
			}
		}
		return current;
	}
	
	/** Create an array from [ to the next ] */
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
				if(slashMode){
					if(verbose){outstream.println(">SlashModeEnd, buffer="+bb);}
					bb.append(b);
					slashMode=false;
				}else if(b=='"'){
					if(verbose){outstream.println(">Quote; quote mode="+quoteMode+", buffer="+bb);}
					String s=bb.toString();
					bb.clear();
					current.add(s);
					quoteMode=!quoteMode;
				}else{
					if(verbose){outstream.println(">QuoteMode, buffer="+bb);}
					if(b=='\\'){
						if(verbose){outstream.println(">SlashMode, buffer="+bb);}
						slashMode=true;
					}
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
				assert(false);
			}else if(b=='{'){
				if(verbose){outstream.println(">{, E) current object is:\n"+current+"\n");}
				JsonObject j=makeObject();
				current.add(j);
				if(verbose){outstream.println("Added new object:\n"+j+"\n");}
			}else if(b=='}'){
				if(verbose){outstream.println(">}, current array is:\n"+current+"\n");}
				assert(false);
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
		return current.toArray();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	byte[] text;
	int pos=0;
	boolean errorState;
	
	/** Always false except when testing */
	private static final boolean verbose=false;
	private static final PrintStream outstream=System.err;
	private static final Exception INVALID_JSON=new Exception("Invalid Json");
	
}
