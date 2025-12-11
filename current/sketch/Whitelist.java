package sketch;

import kmer.AbstractKmerTable;
import structures.LongList;

/**
 * Manages whitelisting of k-mer sketch keys for filtering operations.
 * Distributes keys across multiple k-mer tables for fast lookup and can prune sketches to whitelisted keys.
 * @author Brian Bushnell
 */
public class Whitelist {

	/**
	 * Initializes the whitelist with an array of k-mer tables, distributed by hash.
	 * Asserts that initialization only happens once.
	 * @param tableArray Array of k-mer tables to back the whitelist
	 */
	public static void initialize(AbstractKmerTable[] tableArray){
		assert(keySets==null);
		keySets=tableArray;
	}
	
	/**
	 * Filters a sketch in place, retaining only keys present in the whitelist.
	 * Rebuilds the key array if any keys are removed.
	 * @param s Sketch to filter using the whitelist
	 */
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
	
	/**
	 * Checks if a hashed sketch key exists in the whitelist, using hash-based distribution across tables.
	 * Returns true if no whitelist is active.
	 * @param key Hashed value from an actual sketch
	 * @return true if the key is whitelisted or no whitelist exists
	 */
	public static boolean contains(long key){
		if(keySets==null){return true;}
		int way=(int)(key%ways);
		return keySets[way].getValue(key)>0;
	}
	
	/**
	 * Checks if a raw hashed key exists in the whitelist by converting it to sketch format first.
	 * @param key Raw hashed value which has not yet been subtracted from Long.MAX_VALUE
	 * @return true if the converted key is whitelisted
	 */
	public static boolean containsRaw(long key){
		return contains(Long.MAX_VALUE-key);
	}
	
	/** Returns whether whitelist tables have been initialized.
	 * @return true if a whitelist is active; false otherwise */
	public static boolean exists(){
		return keySets!=null;
	}
	
	/** Whitelist tables distributed by hashing key%ways. */
	private static AbstractKmerTable[] keySets;
	private static final int ways=31;
	
}
