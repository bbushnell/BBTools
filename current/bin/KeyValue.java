package bin;

import java.util.ArrayList;
import java.util.Collections;

import structures.IntHashMap;

//Sorts by B descending then A ascending
/**
 * A key-value pair container that implements custom sorting behavior.
 * Sorts by value descending, then by key ascending for tie-breaking.
 * Used for converting IntHashMap entries to sortable list format.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
class KeyValue implements Comparable<KeyValue> {
	
	/**
	 * Constructs a key-value pair.
	 * @param a_ The key value
	 * @param b_ The associated value
	 */
	KeyValue(int a_, int b_){key=a_; value=b_;}
	
	/**
	 * Converts an IntHashMap to a sorted list of KeyValue pairs.
	 * Filters out invalid entries and sorts by value descending, key ascending.
	 * @param map The IntHashMap to convert
	 * @return Sorted ArrayList of KeyValue pairs, or null if map is null/empty
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
	
	@Override
	public int compareTo(KeyValue o) {
		if(value!=o.value) {return value>o.value ? -1 : 1;}
		return key-o.key;
	}
	
	/** Value component used for primary descending sorting. */
	/** Key component used for secondary ascending sorting. */
	int key, value;
	
}
