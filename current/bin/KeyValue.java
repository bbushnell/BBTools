package bin;

import java.util.ArrayList;
import java.util.Collections;

import structures.IntHashMap;

//Sorts by B descending then A ascending
/**
 * Key-value pair container with custom sorting: value descending, then key ascending.
 * Used for converting IntHashMap entries into a sortable list representation.
 * @author Brian Bushnell
 * @date June 3, 2025
 */
class KeyValue implements Comparable<KeyValue> {
	
	/**
	 * Creates a key-value pair.
	 * @param a_ The key value
	 * @param b_ The associated value
	 */
	KeyValue(int a_, int b_){key=a_; value=b_;}
	
	/**
	 * Converts an IntHashMap to a sorted list of KeyValue pairs, filtering out invalid entries.
	 * Sorts by value descending and key ascending; returns null for null or empty maps.
	 * @param map The IntHashMap to convert
	 * @return Sorted ArrayList of KeyValue pairs, or null if map is null or empty
	 */
	static ArrayList<KeyValue> toList(IntHashMap map){
		if(map==null || map.isEmpty()) {return null;}
		ArrayList<KeyValue> list=new ArrayList<KeyValue>(map.size());
		int[] keys=map.keys();
		int[] values=map.values();
		for(int i=0; i<keys.length; i++) {
			if(keys[i]!=map.invalid()) {
				list.add(new KeyValue(keys[i], values[i]));
			}
		}
		Collections.sort(list);
		return list;
	}
	
	/**
	 * Compares KeyValue objects for sorting: primary by value descending, secondary by key ascending.
	 * @param o The KeyValue to compare against
	 * @return Negative if this should come first, positive if other should come first, 0 if equal
	 */
	@Override
	public int compareTo(KeyValue o) {
		if(value!=o.value) {return value>o.value ? -1 : 1;}
		return key-o.key;
	}
	
	/** Value component used for primary descending sorting. */
	int key, value;
	
}
