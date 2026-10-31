package json;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import shared.Tools;
import structures.ByteBuilder;

/** Mutable JSON object with separate insertion-ordered scalar and child maps.
 * Scalar/array entries are emitted before nested-object entries; insertion order
 * is preserved within each map, not across both maps. Values are borrowed, not
 * copied. Instances require external coordination when mutated or serialized by
 * multiple threads, and nested values must not form cycles.
 * @author Brian Bushnell */
public class JsonObject{

	/**
	 * Test method demonstrating JsonObject usage and nested object creation.
	 * Creates a hierarchy of JsonObjects and demonstrates various add operations,
	 * including object nesting and value replacement.
	 * @param args Command line arguments (unused)
	 */
	public static void main(String[] args){
		JsonObject bob=new JsonObject("name", "bob");
		JsonObject joe=new JsonObject("name", "joe");
		JsonObject sue=new JsonObject("name", "sue");
		JsonObject dan=new JsonObject("name", "dan");
		bob.add("joe", joe, true);
		bob.add("sue", sue, true);
		joe.add("dan", dan, true);
		bob.add("a", 1, true);
		bob.add("b", 2, true);
		bob.add("c", "3", true);
		bob.add("a", "4", true);
		dan.add("e", 5, true);
		dan.add("f", "6", true);
		sue.add("g", "7", true);

		System.out.println("dan:\n"+dan);
		System.out.println("sue:\n"+sue);
		System.out.println("joe:\n"+joe);
		System.out.println("bob:\n"+bob);
		
		ArrayList<JsonObject> list=new ArrayList<JsonObject>();
		list.add(joe);
		list.add(sue);
		list.add(dan);
		System.out.println("list:\n"+toString(list));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an empty JsonObject with no initial key-value pairs. */
	public JsonObject(){}
	
	/** Creates an object borrowing one initial value.
	 * @param key Member name; null is serialized as the name "null"
	 * @param value Initial value, including null */
	public JsonObject(String key, Object value){
		add(key, value, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds a formatted numeric value with specified decimal places as a literal.
	 * The value will be output in JSON without quotes, formatted to the specified precision.
	 *
	 * @param key0 The key to associate with the value
	 * @param value The numeric value to add
	 * @param decimals Number of decimal places to include in the formatted output
	 */
	public void addLiteral(String key0, double value, int decimals){
		add(key0, new JsonLiteral(value, decimals));
	}

	/**
	 * Adds a string value as a literal (without quotes in JSON output).
	 * WARNING: This method should be used with caution as it can produce
	 * incorrectly formatted JSON files if the string contains special characters.
	 *
	 * @param key0 The key to associate with the value
	 * @param value The string value to add as a literal
	 */
	public void addLiteral(String key0, String value){
		add(key0, new JsonLiteral(value));
	}

	/** Adds or replaces a member, removing any older child with the same name.
	 * @param key0 Member name
	 * @param value Borrowed value; exact JsonObject values use the child map */
	public void add(String key0, Object value){add(key0, value, true);}
	
	/**
	 * Adds a key-value pair to this JsonObject, renaming the key if it already exists.
	 * If the key exists, appends " 2", " 3", etc. to make it unique.
	 * @param key0 The base key to add
	 * @param value The value to associate with the key
	 */
	public void addAndRename(String key0, Object value){add(key0, value, false);}

	/** Stores a scalar/array value or dispatches an exact JsonObject value to the child map.
	 * @param key0 Requested member name
	 * @param value Borrowed member value
	 * @param replace True to replace; false to find a name unused in either map */
	private void add(String key0, Object value, boolean replace){
		if(value!=null && value.getClass()==JsonObject.class){
			add(key0, (JsonObject)value, replace);
			return;
		}
		int x=2;
		String key=key0;
		if(omap==null){omap=new LinkedHashMap<String, Object>(8);}
		while(!replace && containsKey(key)){
			key=key0+" "+x;
			x++;
		}
		// One logical key cannot retain an older value in the other storage map.
		if(replace && jmap!=null){jmap.remove(key);}
		omap.put(key, value);
	}

	/** Adds or replaces a nested child; a null child is stored as a scalar JSON null.
	 * @param key0 Member name
	 * @param value Borrowed child, or null */
	public void add(String key0, JsonObject value){add(key0, value, true);}
	
	/** Adds a nested child under a name unused in either storage map.
	 * @param key0 Base member name; suffixes begin with " 2"
	 * @param value Borrowed child, or null */
	public void addAndRename(String key0, JsonObject value){add(key0, value, false);}

	/** Stores a child with replacement or cross-map renaming.
	 * @param key0 Requested member name
	 * @param value Borrowed child; null delegates to scalar storage
	 * @param replace True to replace an existing member of either type */
	private void add(final String key0, JsonObject value, boolean replace){
		if(value==null){add(key0, (Object)null, replace); return;}
		int x=2;
		String key=key0;
		if(jmap==null){jmap=new LinkedHashMap<String, JsonObject>(8);}
		while(!replace && containsKey(key)){
			key=key0+" "+x;
			x++;
		}
		if(replace && omap!=null){omap.remove(key);}
		jmap.put(key, value);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Serialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Converts a list of JsonObjects to a comma-and-newline-separated sequence of objects.
	 * NOTE: this does NOT wrap the output in '[' ']', so the result is not itself a valid
	 * JSON array. Callers that need an array must add the brackets (no current caller does;
	 * used only by main()). Contrast toString(Object[]), which does wrap in brackets.
	 * @param list The list of JsonObjects to convert
	 * @return Comma+newline-separated JsonObjects with no enclosing brackets
	 */
	public static String toString(ArrayList<JsonObject> list){
		ByteBuilder sb=new ByteBuilder();
		int commas=list.size()-1;
		for(JsonObject j : list){
			j.append(0, sb, false);
			if(commas>0){
				sb.append(",\n");
			}
			commas--;
		}
		return sb.toString();
	}
	
	/** Renders this object into a new ASCII JSON buffer.
	 * @return Newly allocated formatted text */
	public ByteBuilder toText(){
		return toText(null, 0, false);
	}
	
	/** Appends this object, retaining any existing destination contents.
	 * @param sb Destination, or null to allocate a new builder
	 * @param level Nonnegative indentation level
	 * @param inArray True to use compact formatting within an array
	 * @return The supplied or newly allocated builder */
	public ByteBuilder toText(ByteBuilder sb, int level, boolean inArray){
		if(sb==null){sb=new ByteBuilder();}
		append(level, sb, inArray);
		return sb;
	}
	
	/**
	 * Converts this JsonObject to a string with a named wrapper.
	 * Creates a JSON object with the specified name containing this object.
	 * @param name The name to use as the wrapper key
	 * @return JSON string with named wrapper
	 */
	public String toString(String name){
		ByteBuilder sb=new ByteBuilder();
		sb.append('{').append('\n');
		for(int i=0; i<padmult; i++){sb.append(' ');}
		appendQuoted(sb, name);
		sb.append(':').append(' ');
		toText(sb, 1, false);
		sb.append('\n').append('}');
		return sb.toString();
	}
	
	/** Renders an array with enclosing brackets; null is rendered as an empty array.
	 * @param array Borrowed values, or null
	 * @return Formatted array text */
	public static String toString(Object[] array){
		ByteBuilder sb=new ByteBuilder();
		appendArray(sb, array, 0);
		return sb.toString();
	}
	
	/** Renders this object with normal indentation.
	 * @return Newly allocated JSON text */
	@Override
	public String toString(){
		return toText(null, 0, false).toString();
	}
	
	/** Renders this object followed by a newline.
	 * @return Newly allocated JSON text ending in a newline */
	public String toStringln(){
		return toText(null, 0, false).nl().toString();
	}
	
	/**
	 * Appends this JsonObject's formatted representation to a ByteBuilder.
	 * Handles indentation, comma placement, and nested structure formatting.
	 *
	 * @param level Current indentation level for pretty-printing
	 * @param sb ByteBuilder to append the formatted JSON to
	 * @param inArray Whether this object is being rendered inside an array (affects formatting)
	 */
	public void append(int level, ByteBuilder sb, boolean inArray){
		int pad=padmult*level;
		int pad2=padmult*(level+1);
		
		sb.append('{');
		if(!inArray){sb.append('\n');}
		
		int commas=(omap==null ? 0 : omap.size())+(jmap==null ? 0 : jmap.size())-1;
		
		if(omap!=null){
			for(Entry<String, Object> e : omap.entrySet()){
				String key=e.getKey();
				Object value=e.getValue();
				if(!inArray){for(int i=0; i<pad2; i++){sb.append(' ');}}
				
				appendEntry(sb, key, value, level, inArray);

				if(commas>0){sb.append(',');}
				if(!inArray){sb.append('\n');}
				commas--;
			}
		}
		
		if(jmap!=null){
			for(Entry<String, JsonObject> e : jmap.entrySet()){
				String key=e.getKey();
				JsonObject value=e.getValue();
				if(!inArray){for(int i=0; i<pad2; i++){sb.append(' ');}}
				appendKey(sb, key);
				
				value.append(level+(inArray ? 0 : 1), sb, inArray);
				if(commas>0){sb.append(',');}
				if(!inArray){sb.append('\n');}
				commas--;
			}
		}
		
		if(!inArray){for(int i=0; i<pad; i++){sb.append(' ');}}
		sb.append('}');
	}
	
	/**
	 * Appends a key-value entry to the ByteBuilder in JSON format.
	 * Helper method that combines key and value formatting.
	 *
	 * @param sb ByteBuilder to append to
	 * @param key The key to format and append
	 * @param value The value to format and append
	 * @param level Current indentation level
	 * @param inArray Whether formatting is within an array context
	 */
	private static void appendEntry(ByteBuilder sb, String key, Object value, int level, boolean inArray){
		appendKey(sb, key);
		appendValue(sb, value, level, inArray);
	}
	
	/**
	 * Appends a JSON key to the ByteBuilder with proper formatting.
	 * Adds quotes around the key and the colon separator.
	 * @param sb ByteBuilder to append to
	 * @param key The key string to format and append
	 */
	private static void appendKey(ByteBuilder sb, String key){
		appendQuoted(sb, key);
		sb.append(':').append(' ');
	}

	/** Hexadecimal digits for escaped UTF-16 units; initialized before format switches. */
	private static final char[] HEX="0123456789abcdef".toCharArray();

	/**
	 * Appends a string as a quoted JSON string, escaping per RFC 8259: '"' and '\\'
	 * are backslash-escaped, common control chars use short escapes (\\n \\r \\t \\b \\f),
	 * and other control or non-ASCII UTF-16 units use hexadecimal escapes. Clean
	 * ASCII runs are appended in bulk. Escaping non-ASCII units avoids truncation
	 * by the deliberately byte-oriented ByteBuilder and preserves surrogate pairs.
	 * @param sb ByteBuilder to append to
	 * @param s String to quote and escape (null is rendered as the literal characters "null")
	 */
	private static void appendQuoted(ByteBuilder sb, String s){
		if(s==null){s="null";}
		sb.append('"');
		final int len=s.length();
		int start=0;
		for(int i=0; i<len; i++){
			final char c=s.charAt(i);
			if(c=='"' || c=='\\' || c<0x20 || c>=128){
				if(i>start){sb.append(s, start, i);}
				switch(c){
					case '"': sb.append('\\', '"'); break;
					case '\\': sb.append('\\', '\\'); break;
					case '\n': sb.append('\\', 'n'); break;
					case '\r': sb.append('\\', 'r'); break;
					case '\t': sb.append('\\', 't'); break;
					case '\b': sb.append('\\', 'b'); break;
					case '\f': sb.append('\\', 'f'); break;
					default: sb.append('\\', 'u').append(HEX[(c>>>12)&15], HEX[(c>>>8)&15],
							HEX[(c>>>4)&15], HEX[c&15]);
				}
				start=i+1;
			}
		}
		if(start<len){sb.append(s, start, len);}//flush trailing clean run
		sb.append('"');
	}
	
	/**
	 * Appends a value to the ByteBuilder with appropriate JSON formatting.
	 * Handles type detection and applies correct formatting rules for each type.
	 * Supports strings, numbers, booleans, null, arrays, collections, and nested objects.
	 * Non-finite Float/Double values are emitted as JSON null.
	 *
	 * @param sb ByteBuilder to append to
	 * @param value The value to format (may be null)
	 * @param level Current indentation level for nested structures
	 * @param inArray Whether this value is being formatted within an array
	 */
	private static void appendValue(ByteBuilder sb, Object value, int level, boolean inArray){
		final Class<?> c=(value==null ? null : value.getClass());
		if(c==null || value==null){
			sb.append("null");
		}else if(c==String.class){
			// [json/JsonObject#001] FIXED: escape string values (paired with JsonParser#002 un-escaper)
			appendQuoted(sb, value.toString());
		}else if(c==JsonLiteral.class){
			sb.append(((JsonLiteral)value).toString());
		}else if(c==Double.class){
			final double number=((Double)value).doubleValue();
			if(!Double.isFinite(number)){sb.append("null");}
			else if(restrictDecimals>=0){appendRestrictedNumber(sb, number, restrictDecimals);}
			else{sb.append(value.toString());}
		}else if(c==Float.class){
			final float number=((Float)value).floatValue();
			if(!Float.isFinite(number)){sb.append("null");}
			else if(restrictDecimals>=0){appendRestrictedNumber(sb, number, restrictDecimals);}
			else{sb.append(value.toString());}
		}else if(c==JsonObject.class){
			((JsonObject)value).append(level+(inArray ? 0 : 1), sb, inArray);
		}else if(value instanceof Object[]){
			appendArray(sb, (Object[])value, level);
		}else if(c.isArray()){// [json/JsonObject#003] FIXED: primitive array via reflection (was ClassCastException)
			appendPrimitiveArray(sb, value, level);
		}else if(c==Boolean.class || value instanceof Number){ //long, int, boolean
			sb.append(value.toString());
		}else if(value instanceof Collection){
			appendCollection(sb, (Collection<?>)value, level);
		}else{ //Default behavior for unhandled classes
			appendQuoted(sb, value.toString());
		}
	}
	
	/** Formats a floating value without overflowing ByteBuilder's long digit accumulator.
	 * Integral values retain the historical compact form. Ordinary magnitudes with
	 * up to 18 fraction digits use ByteBuilder; other cases use the ROOT formatter.
	 * @param sb Destination retaining its prefix
	 * @param value Floating value
	 * @param decimals Nonnegative requested fraction digits */
	private static void appendRestrictedNumber(final ByteBuilder sb, final double value,
			final int decimals){
		assert(decimals>=0) : "Restricted JSON formatting is selected only for nonnegative precision";
		if(value>=Long.MIN_VALUE && value<0x1.0p63 && value==(long)value){
			sb.append((long)value);
		}else if(decimals<=18 && Math.abs(value)<0x1.0p63){
			sb.append(value, decimals);
		}else{sb.append(Tools.format("%."+decimals+"f", value));}
	}

	/**
	 * Appends an Object array to the ByteBuilder in JSON array format.
	 * Elements are separated by commas and spaces for readability.
	 *
	 * @param sb ByteBuilder to append to
	 * @param array The Object array to format (may be null)
	 * @param level Current indentation level for nested elements
	 */
	private static void appendArray(ByteBuilder sb, Object[] array, int level){
		int commas=(array==null ? 0 : array.length)-1;
		sb.append('[');
		if(array!=null){
			for(Object value : array){
				appendValue(sb, value, level, noNewlinesInArrays);
				if(commas>0){sb.append(',').append(' ');}
				commas--;
			}
		}
		sb.append(']');
	}

	/**
	 * Appends a primitive array (int[], long[], double[], etc.) in JSON array format,
	 * using reflection so each element is boxed and routed through appendValue.
	 * @param sb ByteBuilder to append to
	 * @param array A non-null primitive array
	 * @param level Current indentation level for nested elements
	 */
	private static void appendPrimitiveArray(ByteBuilder sb, Object array, int level){
		final int n=java.lang.reflect.Array.getLength(array);
		sb.append('[');
		for(int i=0; i<n; i++){
			appendValue(sb, java.lang.reflect.Array.get(array, i), level, noNewlinesInArrays);
			if(i<n-1){sb.append(',').append(' ');}
		}
		sb.append(']');
	}

	/**
	 * Appends a Collection to the ByteBuilder in JSON array format.
	 * Similar to appendArray but works with any Collection type.
	 *
	 * @param sb ByteBuilder to append to
	 * @param stuff The Collection to format (may be null)
	 * @param level Current indentation level for nested elements
	 */
	private static void appendCollection(ByteBuilder sb, Collection<?> stuff, int level){
		int commas=(stuff==null ? 0 : stuff.size())-1;
		sb.append('[');
		if(stuff!=null){
			for(Object value : stuff){
				appendValue(sb, value, level, noNewlinesInArrays);
				if(commas>0){sb.append(',').append(' ');}
				commas--;
			}
		}
		sb.append(']');
	}

	/*--------------------------------------------------------------*/
	/*----------------            Access            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns a String member without coercing other value types.
	 * @param key Scalar member name
	 * @return Stored string, or null for an absent/null scalar member */
	public String getString(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		assert(o.getClass()==String.class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (String)o;
	}

	/** Returns a Long member, including integers produced by JsonParser.
	 * @param key Scalar member name; non-null values must have type Long
	 * @return Stored Long, or null for an absent/null scalar member */
	public Long getLong(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		assert(o.getClass()==Long.class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (Long)o;
	}

	/** Returns an Integer member or narrows an in-range parsed Long.
	 * @param key Scalar member name; Long values must fit a signed int
	 * @return Stored/coerced integer, or null for an absent/null scalar member */
	public Integer getInt(String key){
		// [json/JsonObject#002] FIXED: null-map guard + Long->int coercion (mirrors getDouble/getLong)
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		if(o.getClass()==Long.class){//JsonParser stores ints as Long; coerce
			long x=((Long)o).longValue();
			assert(x>=Integer.MIN_VALUE && x<=Integer.MAX_VALUE) : "Value out of int range: "+x;
			return (int)x;
		}
		assert(o.getClass()==Integer.class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (Integer)o;
	}
	
	/**
	 * Checks whether this JsonObject contains the specified key.
	 * Searches both the object map (omap) and JsonObject map (jmap).
	 * @param key The key to search for
	 * @return true if the key exists in either map, false otherwise
	 */
	public boolean containsKey(String key){
		if(omap!=null && omap.containsKey(key)){return true;}
		if(jmap!=null && jmap.containsKey(key)){return true;}
		return false;
	}

	/**
	 * Retrieves a Double value associated with the specified key.
	 * Automatically converts Long values to Double if needed.
	 *
	 * @param key The key to look up
	 * @return The Double value, or null if key doesn't exist or omap is null
	 * @throws AssertionError if the value exists but is not a Double or Long
	 */
	public Double getDouble(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		if(o.getClass()==Long.class){
			return ((Long)o).doubleValue();
		}
		assert(o.getClass()==Double.class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (Double)o;
	}

	/**
	 * Retrieves a Number value associated with the specified key.
	 * Accepts Double, Long, Integer, or Float values.
	 *
	 * @param key The key to look up
	 * @return The Number value, or null if key doesn't exist or omap is null
	 * @throws AssertionError if the value exists but is not a supported Number type
	 */
	public Number getNumber(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		Class<?> c=o.getClass();
		assert(c==Double.class || c==Long.class || c==Integer.class || c==Float.class) : "Wrong class: "+c+"\n"+o;
		return (Number)o;
	}

	/** Returns an Object[] member without converting primitive arrays or collections.
	 * @param key Scalar-map member name; non-null values must have type Object[]
	 * @return Stored array, or null for an absent/null member */
	public Object[] getArray(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		assert(o.getClass()==Object[].class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (Object[])o;
	}

	/** Returns a borrowed nested child.
	 * @param key Child member name
	 * @return Stored child, or null if absent */
	public JsonObject getJson(String key){
		if(jmap==null){return null;}
		return jmap.get(key);
	}

	/** Removes only a nested-child member.
	 * @param key Child member name
	 * @return Removed child, or null if absent */
	public JsonObject removeJson(String key){
		if(jmap==null){return null;}
		return jmap.remove(key);
	}

	/** Removes only a scalar/array member.
	 * @param key Scalar member name
	 * @return Removed value, or null if absent or mapped to null */
	public Object removeObject(String key){
		if(omap==null){return null;}
		return omap.remove(key);
	}

	/** Clears all nested JsonObjects by setting jmap to null */
	public void clearJson(){
		jmap=null;
	}

	/** Clears all objects by setting omap to null */
	public void clearOmap(){
		omap=null;
	}
	
	/** Copies child references in their insertion order.
	 * @return New array of borrowed children, or null when the child map is unallocated */
	public Object[] toJmapArray(){
		if(jmap==null){return null;}
		Object[] array=new Object[jmapSize()];
		int i=0;
		for(Entry<String, JsonObject> e : jmap.entrySet()){
			array[i]=e.getValue();
			i++;
		}
		return array;
	}
	
	/** Counts stored nested-child members.
	 * @return Child-map size, or zero if unallocated */
	public int jmapSize(){return jmap==null ? 0 : jmap.size();}
	
	/** Counts scalar/array members, including members mapped to null.
	 * @return Scalar-map size, or zero if unallocated */
	public int omapSize(){return omap==null ? 0 : omap.size();}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Insertion-ordered scalar/array values, including null; lazily allocated.
	 * An allocated map may be empty. Direct edits must not duplicate child-map keys. */
	public LinkedHashMap<String, Object> omap;

	/** Insertion-ordered non-null children. Direct mutations must not duplicate omap keys. */
	public LinkedHashMap<String, JsonObject> jmap;

	/** Global floating formatting precision; negative keeps each value's normal text. */
	private static int restrictDecimals=-1;

	/** Sets global Float/Double formatting precision before concurrent serialization.
	 * JsonLiteral values keep their own preformatted precision.
	 * Integral values remain compact; exceptional ranges/precision use ROOT formatting.
	 * @param d Negative for default formatting, otherwise requested fraction digits */
	public static synchronized void setDecimals(int d){
		restrictDecimals=d;
	}
	
	/** Spaces per indentation level. */
	public static final int padmult=3;
	
	/** Global compact-object switch for arrays; configure before concurrent use. */
	public static boolean noNewlinesInArrays=false;
	
}
