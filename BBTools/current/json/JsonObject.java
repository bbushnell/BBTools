package json;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import structures.ByteBuilder;

/**
 * Custom JSON object implementation for BBTools output formatting. Preserves key insertion
 * order and supports nested objects, arrays, and literal values for numeric formatting.
 * 
 * @author Brian Bushnell
 * @date Pre-2020 (before Java built-in JSON support)
 * @documentation Eru
 */
public class JsonObject {

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
		bob.add("a",1, true);
		bob.add("b",2, true);
		bob.add("c","3", true);
		bob.add("a","4", true);
		dan.add("e",5, true);
		dan.add("f","6", true);
		sue.add("g","7", true);

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
	
	/**
	 * Creates an empty JsonObject with no initial key-value pairs.
	 */
	public JsonObject(){}
	
	/**
	 * Creates a JsonObject with a single initial key-value pair.
	 * @param key The initial key to add
	 * @param value The initial value to associate with the key
	 */
	public JsonObject(String key, Object value){
		add(key, value, true);
	}
	
//	public JsonObject(String name_){
//		name=name_;
//	}
//	
//	public JsonObject(String name_, String key, Object value){
//		name=name_;
//		add(key, value);
//	}

	/**
	 * Adds a formatted numeric value with specified decimal places as a literal.
	 * The value will be output in JSON without quotes, formatted to the specified precision.
	 * @param key0 The key to associate with the value
	 * @param value The numeric value to add
	 * @param decimals Number of decimal places to include in the formatted output
	 */
	public void addLiteral(String key0, double value, int decimals){
		if(omap==null){omap=new LinkedHashMap<String, Object>(8);}
		omap.put(key0, new JsonLiteral(value, decimals));
	}

	/**
	 * Adds a string value as a literal (without quotes in JSON output).
	 * WARNING: This method should be used with caution as it can produce 
	 * incorrectly formatted JSON files if the string contains special characters.
	 * @param key0 The key to associate with the value
	 * @param value The string value to add as a literal
	 */
	public void addLiteral(String key0, String value){
		if(omap==null){omap=new LinkedHashMap<String, Object>(8);}
		omap.put(key0, new JsonLiteral(value));
	}

	/**
	 * Adds a key-value pair to this JsonObject, replacing any existing value.
	 * @param key0 The key to add
	 * @param value The value to associate with the key
	 */
	public void add(String key0, Object value){add(key0, value, true);}
	
	/**
	 * Adds a key-value pair to this JsonObject, renaming the key if it already exists.
	 * If the key exists, appends " 2", " 3", etc. to make it unique.
	 * @param key0 The base key to add
	 * @param value The value to associate with the key
	 */
	public void addAndRename(String key0, Object value){add(key0, value, false);}

	/**
	 * Internal method for adding objects with replacement or renaming behavior.
	 * @param key0 The base key name
	 * @param value The value to add
	 * @param replace If true, replaces existing values; if false, renames key to avoid conflicts
	 */
	private void add(String key0, Object value, boolean replace){
		if(value!=null && value.getClass()==JsonObject.class){
			add(key0, (JsonObject)value, replace);
			return;
		}
		int x=2;
		String key=key0;
		if(omap==null){omap=new LinkedHashMap<String, Object>(8);}
		while(!replace && omap.containsKey(key)){
			key=key0+" "+x;
			x++;
		}
		omap.put(key, value);
	}

	/**
	 * Adds a JsonObject as a nested object, replacing any existing value.
	 * @param key0 The key to add
	 * @param value The JsonObject to nest under this key
	 */
	public void add(String key0, JsonObject value){add(key0, value, true);}
	
	/**
	 * Adds a JsonObject as a nested object, renaming the key if it already exists.
	 * @param key0 The base key to add
	 * @param value The JsonObject to nest under this key
	 */
	public void addAndRename(String key0, JsonObject value){add(key0, value, false);}

	/**
	 * Internal method for adding JsonObjects with replacement or renaming behavior.
	 * @param key0 The base key name
	 * @param value The JsonObject to add
	 * @param replace If true, replaces existing values; if false, renames key to avoid conflicts
	 */
	private void add(final String key0, JsonObject value, boolean replace){
		int x=2;
		String key=key0;
		if(jmap==null){jmap=new LinkedHashMap<String, JsonObject>(8);}
		while(!replace && jmap.containsKey(key)){
			key=key0+" "+x;
			x++;
		}
		jmap.put(key, value);
	}
	
	/**
	 * Converts a list of JsonObjects to a JSON array string representation.
	 * Objects are separated by commas and newlines for readability.
	 * @param list The list of JsonObjects to convert
	 * @return JSON array string representation
	 */
	public static String toString(ArrayList<JsonObject> list) {
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
	
	/**
	 * Converts this JsonObject to formatted text using default parameters.
	 * @return ByteBuilder containing the formatted JSON text
	 */
	public ByteBuilder toText(){
		return toText(null, 0, false);
	}
	
	/**
	 * Converts this JsonObject to formatted text with specified parameters.
	 * @param sb ByteBuilder to append to (created if null)
	 * @param level Indentation level for pretty-printing
	 * @param inArray Whether this object is inside an array (affects formatting)
	 * @return ByteBuilder containing the formatted JSON text
	 */
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
		sb.append('"').append(name).append('"').append(':').append(' ');
		toText(sb, 1, false);
		sb.append('\n').append('}');
		return sb.toString();
	}
	
	/**
	 * Converts an array of objects to JSON array string representation.
	 * @param array The array to convert to JSON
	 * @return JSON array string
	 */
	public static String toString(Object[] array){
		ByteBuilder sb=new ByteBuilder();
		appendArray(sb, array, 0);
		return sb.toString();
	}
	
	/**
	 * Returns the JSON string representation of this object.
	 * @return Formatted JSON string
	 */
	@Override
	public String toString(){
		return toText(null, 0, false).toString();
	}
	
	/**
	 * Returns the JSON string representation with a trailing newline.
	 * @return Formatted JSON string with newline
	 */
	public String toStringln(){
		return toText(null, 0, false).nl().toString();
	}
	
	/**
	 * Appends this JsonObject's formatted representation to a ByteBuilder.
	 * Handles indentation, comma placement, and nested structure formatting.
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
		sb.append('"').append(key).append("\": ");
	}
	
	/**
	 * Appends a value to the ByteBuilder with appropriate JSON formatting.
	 * Handles type detection and applies correct formatting rules for each type.
	 * Supports strings, numbers, booleans, null, arrays, collections, and nested objects.
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
			sb.append("\"").append(value.toString()).append('"');
		}else if(c==JsonLiteral.class){
			sb.append(((JsonLiteral)value).toString());
		}else if(c==Double.class && restictDecimals>=0){
			sb.append(((Double)value).doubleValue(), restictDecimals);
		}else if(c==Float.class && restictDecimals>=0){
			sb.append(((Float)value).floatValue(), restictDecimals);
		}else if(c==JsonObject.class){
			((JsonObject)value).append(level+(inArray ? 0 : 1), sb, inArray);
		}else if(c.isArray()){
			appendArray(sb, (Object[])value, level);
		}else if(c==Boolean.class || value instanceof Number){ //long, int, boolean
			sb.append(value.toString());
		}else if(value instanceof Collection){
			appendCollection(sb, (Collection<?>)value, level);
		}else{ //Default behavior for unhandled classes
			sb.append("\"").append(value.toString()).append('"');
		}
	}
	
	/**
	 * Appends an Object array to the ByteBuilder in JSON array format.
	 * Elements are separated by commas and spaces for readability.
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
	 * Appends a Collection to the ByteBuilder in JSON array format.
	 * Similar to appendArray but works with any Collection type.
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

	/**
	 * Retrieves a String value associated with the specified key.
	 * @param key The key to look up
	 * @return The String value, or null if key doesn't exist or omap is null
	 * @throws AssertionError if the value exists but is not a String
	 */
	public String getString(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		assert(o.getClass()==String.class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (String)o;
	}

	/**
	 * Retrieves a Long value associated with the specified key.
	 * @param key The key to look up
	 * @return The Long value, or null if key doesn't exist or omap is null
	 * @throws AssertionError if the value exists but is not a Long
	 */
	public Long getLong(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		assert(o.getClass()==Long.class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (Long)o;
	}

	/**
	 * Retrieves an Integer value associated with the specified key.
	 * @param key The key to look up
	 * @return The Integer value, or null if key doesn't exist
	 * @throws AssertionError if omap is null or the value exists but is not an Integer
	 */
	public Integer getInt(String key){
		assert(omap!=null);
		Object o=omap.get(key);
//		assert(o!=null);
		assert(o==null || o.getClass()==Integer.class) : "Wrong class: "+o.getClass()+"\n"+o;
//		long x=((Long)o).longValue();
//		assert(x>=Integer.MIN_VALUE && x<=Integer.MAX_VALUE);
//		return (int)x;
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

//	public Double getDouble(String key){
//		if(smap==null){return null;}
//		Object o=smap.get(key);
//		if(o==null){return null;}
//		assert(o.getClass()==Double.class) : "Wrong class: "+o.getClass()+"\n"+o;
//		return (Double)o;
//	}

	/**
	 * Retrieves a Double value associated with the specified key.
	 * Automatically converts Long values to Double if needed.
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

	/**
	 * Retrieves an Object array associated with the specified key.
	 * @param key The key to look up
	 * @return The Object array, or null if key doesn't exist or omap is null
	 * @throws AssertionError if the value exists but is not an Object array
	 */
	public Object[] getArray(String key){
		if(omap==null){return null;}
		Object o=omap.get(key);
		if(o==null){return null;}
		assert(o.getClass()==Object[].class) : "Wrong class: "+o.getClass()+"\n"+o;
		return (Object[])o;
	}

	/**
	 * Retrieves a nested JsonObject associated with the specified key.
	 * @param key The key to look up
	 * @return The JsonObject, or null if key doesn't exist or jmap is null
	 */
	public JsonObject getJson(String key){
		if(jmap==null){return null;}
		return jmap.get(key);
	}

	/**
	 * Removes and returns a nested JsonObject associated with the specified key.
	 * @param key The key to remove
	 * @return The removed JsonObject, or null if key doesn't exist or jmap is null
	 */
	public JsonObject removeJson(String key){
		if(jmap==null){return null;}
		return jmap.remove(key);
	}

	/**
	 * Removes and returns an Object associated with the specified key.
	 * @param key The key to remove
	 * @return The removed Object, or null if key doesn't exist or omap is null
	 */
	public Object removeObject(String key){
		if(omap==null){return null;}
		return omap.remove(key);
	}

	/**
	 * Clears all nested JsonObjects by setting jmap to null.
	 */
	public void clearJson(){
		jmap=null;
	}

	/**
	 * Clears all objects by setting omap to null.
	 */
	public void clearOmap(){
		omap=null;
	}
	
	/**
	 * Converts the nested JsonObject map to an array of JsonObjects.
	 * @return Array of JsonObjects from jmap, or null if jmap is null
	 */
	public Object[] toJmapArray() {
		if(jmap==null){return null;}
		Object[] array=new Object[jmapSize()];
		int i=0;
		for(Entry<String, JsonObject> e : jmap.entrySet()){
			array[i]=e.getValue();
			i++;
		}
		return array;
	}
	
	/**
	 * Returns the number of entries in the nested JsonObject map.
	 * @return Size of jmap, or 0 if jmap is null
	 */
	public int jmapSize(){return jmap==null ? 0 : jmap.size();}
	
	/**
	 * Returns the number of entries in the object map.
	 * @return Size of omap, or 0 if omap is null
	 */
	public int omapSize(){return omap==null ? 0 : omap.size();}
	
//	public String name;
	/** Map for storing key-value pairs where values are primitives, arrays, or other objects */
	public LinkedHashMap<String, Object> omap;
	
	/** Map for storing key-value pairs where values are nested JsonObjects */
	public LinkedHashMap<String, JsonObject> jmap;

	/** Global setting for restricting decimal places in numeric output (-1 = no restriction) */
	private static int restictDecimals=-1;
	
	/** Format string for decimal output based on restictDecimals setting */
	private static String decimalFormat="%."+restictDecimals+"f";
	
	/**
	 * Sets the global decimal precision for numeric output.
	 * @param d Number of decimal places to display (-1 for no restriction)
	 */
	public static synchronized void setDecimals(int d){
		if(d!=restictDecimals){
			restictDecimals=d;
			decimalFormat="%."+restictDecimals+"f";
		}
	}
	
	/** Multiplier for indentation padding (number of spaces per level) */
	public static final int padmult=3;
	
	/** If true, arrays are formatted without newlines for more compact output */
	public static boolean noNewlinesInArrays=false;
	
}
