package bin;

import java.util.ArrayList;
import java.util.Collections;

import structures.IntHashMap;

//Sorts by B descending then A ascending
class KeyValue implements Comparable<KeyValue> {
	
	KeyValue(int a_, int b_){key=a_; value=b_;}
	
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
	
	int key, value;
	
}
