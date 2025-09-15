package sketch;

import kmer.AbstractKmerTable;
import structures.LongList;

public class Whitelist {

	public static void initialize(AbstractKmerTable[] tableArray){
		assert(keySets==null);
		keySets=tableArray;
	}
	
	public static void apply(Sketch s){
		assert(exists());
		LongList list=new LongList(s.keys.length);
		for(long key : s.keys){
			if(contains(key)){
				list.add(key);
			}
		}
		if(list.size()!=s.keys.length){
			s.keys=list.toArray();
		}
	}
	
	/** Hashed value from an actual sketch */
	public static boolean contains(long key){
		if(keySets==null){return true;}
		int way=(int)(key%ways);
		return keySets[way].getValue(key)>0;
	}
	
	/** Raw hashed value which has not yet been subtracted from Long.MAX_VALUE */
	public static boolean containsRaw(long key){
		return contains(Long.MAX_VALUE-key);
	}
	
	public static boolean exists(){
		return keySets!=null;
	}
	
	/** Hold codes.  A code X such that X%WAYS=Y will be stored in keySets[Y] */
	private static AbstractKmerTable[] keySets;
	private static final int ways=31;
	
}
