package structures;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Tools;

/**
 * A specialized hash map that maps integer keys to IntList values.
 * Uses open addressing with linear probing for collision resolution.
 * Optimized for scenarios where each key maps to a collection of integers.
 * 
 * @author Brian Bushnell
 * @contributor Isla Winglet
 * @date June 2, 2025
 */
public final class IntListHashMap implements Serializable {
	
	/** Serialization version identifier */
	private static final long serialVersionUID = 1L;
	
	/** Default constructor with initial capacity of 256 */
	public IntListHashMap(){this(256);}
	
	/** Constructor with specified initial size and default load factor */
	public IntListHashMap(int initialSize){this(initialSize, 0.7f);}
	
	/**
	 * Constructor with specified initial size and load factor.
	 * @param initialSize Initial capacity of the hash table
	 * @param loadFactor_ Maximum load factor before resizing (0.0-1.0)
	 */
	public IntListHashMap(int initialSize, float loadFactor_){
		invalid=randy.nextInt()|MINMASK; // Ensure negative invalid key
		assert(invalid<0);
		assert(initialSize>0);
		assert(loadFactor_>0 && loadFactor_<1);
		loadFactor=Tools.mid(0.25f, loadFactor_, 0.90f); // Clamp load factor
		keyList=new IntList(); // Track all valid keys
		resize(initialSize);
	}
	
	/** Removes all key-value pairs from the map */
	public void clear(){
		if(size<1){return;}
		for(int i=0; i<keyList.size; i++){
			int key=keyList.array[i];
			int cell=findCell(key);
			if(cell>=0){keys[cell]=invalid; values[cell]=null;} // Clear found cells
		}
		keyList.clear();
		size=0;
	}
	
	/** Returns true if the map contains the specified key */
	public final boolean contains(int key){return findCell(key)>=0;}
	
	/** Returns true if the map contains the specified key */
	public final boolean containsKey(int key){return findCell(key)>=0;}
	
	/**
	 * Returns the IntList associated with the specified key.
	 * @param key The key to look up
	 * @return The IntList value, or null if key not found
	 */
	public IntList get(int key){
		int cell=findCell(key);
		return cell<0 ? null : values[cell];
	}
	
	/**
	 * Returns the IntList for the key, creating a new one if necessary.
	 * @param key The key to look up or create
	 * @return The IntList value (never null)
	 */
	public IntList getOrCreate(int key){
		IntList list=get(key);
		if(list==null){list=new IntList(2); put(key, list);} // Create new list
		return list;
	}
	
	/**
	 * Associates the specified IntList with the specified key.
	 * @param key The key to store
	 * @param value The IntList value to associate
	 * @return The previous IntList value, or null
	 */
	public IntList put(int key, IntList value){
		if(key==invalid){resetInvalid();} // Handle collision with invalid key
		final int cell=findCellOrEmpty(key);
		final IntList oldV=values[cell];
		values[cell]=value;
		if(keys[cell]==invalid){ // New key
			keys[cell]=key;
			keyList.add(key); // Track in key list
			size++;
			if(size>sizeLimit){resize();} // Resize if needed
		}
		return oldV;
	}
	
	/** Adds a single integer value to the IntList for the specified key */
	public void put(int key, int value){getOrCreate(key).add(value);}
	
	/** Copies all key-value pairs from another IntListHashMap */
	public void putAll(IntListHashMap map){
		for(int i=0; i<map.keyList.size; i++){
			int key=map.keyList.array[i];
			IntList list=map.get(key);
			if(list!=null){put(key, list.copy());} // Deep copy lists
		}
	}
	
	/**
	 * Removes the specified key and its associated value.
	 * @param key The key to remove
	 * @return True if the key was found and removed
	 */
	public boolean remove(int key){
		if(key==invalid){return false;}
		final int cell=findCell(key);
		if(cell<0){return false;}
		assert(keys[cell]==key);
		keys[cell]=invalid; values[cell]=null; size--;
		
		// Remove from keyList efficiently
		for(int i=0; i<keyList.size; i++){
			if(keyList.array[i]==key){
				keyList.array[i]=keyList.array[keyList.size-1]; // Swap with last
				keyList.size--; break;
			}
		}
		rehashFrom(cell); // Maintain hash table integrity
		return true;
	}
	
	/** Rehashes entries following a deletion to maintain clustering */
	private void rehashFrom(int initial){
		if(size<1){return;}
		final int limit=keys.length;
		// Rehash entries after deletion point
		for(int cell=initial+1; cell<limit; cell++){
			final int key=keys[cell];
			if(key==invalid){return;} // Stop at first empty cell
			rehashCell(cell);
		}
		// Wrap around to beginning
		for(int cell=0; cell<initial; cell++){
			final int key=keys[cell];
			if(key==invalid){return;}
			rehashCell(cell);
		}
	}
	
	/** Attempts to move an entry closer to its ideal position */
	private boolean rehashCell(final int cell){
		final int key=keys[cell]; final IntList value=values[cell];
		assert(key!=invalid);
		if(key==invalid){resetInvalid();}
		final int dest=findCellOrEmpty(key);
		if(cell==dest){return false;} // Already in correct position
		assert(keys[dest]==invalid);
		keys[cell]=invalid; values[cell]=null; // Clear old position
		keys[dest]=key; values[dest]=value; // Set new position
		return true;
	}
	
	/** Generates a new invalid key when collision occurs */
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;} // Find unused negative
		assert(x<0);
		invalid=x;
		for(int i=0; i<keys.length; i++){
			if(keys[i]==old){keys[i]=invalid;} // Update old invalid markers
		}
	}
	
	/** Finds the cell containing the specified key, or -1 if not found */
	int findCell(final int key){
		if(key==invalid){return -1;}
		final int limit=keys.length, initial=(int)((key&MASK)%modulus);
		// Linear probe from initial position
		for(int cell=initial; cell<limit; cell++){
			final int x=keys[cell];
			if(x==key){return cell;} if(x==invalid){return -1;}
		}
		// Wrap around to beginning
		for(int cell=0; cell<initial; cell++){
			final int x=keys[cell];
			if(x==key){return cell;} if(x==invalid){return -1;}
		}
		return -1;
	}
	
	/** Finds cell containing key or first empty cell for insertion */
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		final int limit=keys.length, initial=(int)((key&Integer.MAX_VALUE)%modulus);
		// Linear probe for key or empty cell
		for(int cell=initial; cell<limit; cell++){
			final int x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		// Wrap around
		for(int cell=0; cell<initial; cell++){
			final int x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}
	
	/** Doubles the hash table size when load factor exceeded */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}
	
	/** Resizes hash table to specified size */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		long newPrime=Primes.primeAtLeast(size2); // Use prime for better distribution
		if(newPrime+extra>Integer.MAX_VALUE){
			newPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra); // Avoid overflow
		}
		assert(newPrime>modulus) : "Overflow: "+size+", "+size2+", "+modulus+", "+newPrime;
		modulus=(int)newPrime;
		
		final int size3=(int)(newPrime+extra);
		sizeLimit=(int)(modulus*loadFactor);
		final int[] oldK=keys; final IntList[] oldV=values;
		keys=KillSwitch.allocInt1D(size3); values=new IntList[size3];
		Arrays.fill(keys, invalid); // Initialize with invalid markers
		
		if(size<1){return;}
		
		size=0; keyList.clear(); // Reset for rehashing
		for(int i=0; i<oldK.length; i++){
			final int k=oldK[i]; final IntList v=oldV[i];
			if(k!=invalid){put(k, v);} // Rehash all valid entries
		}
	}
	
	/** Returns array of all keys in the map */
	public int[] toArray(){return keyList.toArray();}
	
	/** Returns the internal keys array (for advanced use) */
	public int[] keys(){return keys;}
	
	/** Returns the internal values array (for advanced use) */
	public IntList[] values(){return values;}
	
	/** Returns the current invalid key marker */
	public int invalid(){return invalid;}
	
	/** Returns the number of key-value pairs in the map */
	public int size(){return size;}
	
	/** Returns true if the map contains no key-value pairs */
	public boolean isEmpty(){return size==0;}
	
	/** Hash table for keys */
	private int[] keys;
	/** Hash table for values */
	private IntList[] values;
	/** List of all valid keys for iteration */
	private IntList keyList;
	/** Current number of key-value pairs */
	private int size=0;
	/** Sentinel value for empty cells (always negative) */
	private int invalid;
	/** Prime modulus for hash calculation */
	private int modulus;
	/** Maximum size before resizing */
	private int sizeLimit;
	/** Load factor threshold for resizing */
	private final float loadFactor;
	
	/** Mask for positive hash values */
	static final int MASK=Integer.MAX_VALUE;
	/** Mask to ensure negative invalid keys */
	static final int MINMASK=Integer.MIN_VALUE;
	/** Extra space buffer for hash table */
	private static final int extra=10;
	/** Random number generator for invalid key generation */
	private static final Random randy=new Random(1);
}