package sketch;

import kmer.AbstractKmerTable;
import structures.LongList;

/**
 * Manages whitelisting of k-mer keys for sketch filtering operations.
 * Uses hash-distributed storage across multiple k-mer tables for efficient lookup.
 * Provides filtering capability to retain only whitelisted keys from sketches.
 * @author Brian Bushnell
 */
public class Whitelist {

	/**
	 * Initializes the whitelist with an array of k-mer tables.
	 * Tables are distributed by hash to balance load across ways.
	 * @param tableArray Array of k-mer tables to use for whitelist storage
	 */
	public static void initialize(AbstractKmerTable[] tableArray){
		assert(keySets==null);
		keySets=tableArray;
	}
	
	/**
	 * Filters a sketch to retain only whitelisted keys.
	 * Creates new key array containing only keys found in the whitelist.
	 * Modifies the sketch in-place if filtering removes any keys.
	 * @param s The sketch to filter using the whitelist
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
	
	/** Checks if a whitelist has been initialized.
	 * @return true if whitelist tables have been set, false otherwise */
	public static boolean exists(){
		return keySets!=null;
	}
	
	/** Hold codes.  A code X such that X%WAYS=Y will be stored in keySets[Y] */
	private static AbstractKmerTable[] keySets;
	/** Number of ways to distribute keys across k-mer tables */
	private static final int ways=31;
	
}
