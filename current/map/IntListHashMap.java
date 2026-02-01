package map;

import java.io.Serializable;
import java.util.Arrays;
import shared.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Tools;
import structures.IntList;

/**
 * A specialized hash map that maps integer keys to IntList values.
 * Uses open addressing with linear probing for collision resolution.
 * Optimized for scenarios where each key maps to a collection of integers.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date June 2, 2025
 */
public final class IntListHashMap implements Serializable {
	
	private static final long serialVersionUID = 1L;
	
	/** Creates an IntListHashMap with default initial capacity and load factor. */
	public IntListHashMap(){this(256);}
	
	/** Creates an IntListHashMap with specified initial capacity and default load factor.
	 * @param initialSize Initial hash table capacity */
	public IntListHashMap(int initialSize){this(initialSize, 0.7f);}
	
	/**
	 * Creates an IntListHashMap with specified capacity and load factor, initializing
	 * the backing arrays and key tracker.
	 * @param initialSize Initial capacity for the hash table
	 * @param loadFactor_ Maximum load factor before resizing
	 */
	public IntListHashMap(int initialSize, double loadFactor_){
		invalid=randy.nextInt()|MINMASK; // Ensure negative invalid key
		assert(invalid<0);
		assert(initialSize>0);
		assert(loadFactor_>0 && loadFactor_<1);
		loadFactor=Tools.mid(0.25f, (float)loadFactor_, 0.90f); // Clamp load factor
		keyList=new IntList(); // Track all valid keys
		resize(initialSize);
	}
	
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
	
	public final boolean contains(int key){return findCell(key)>=0;}
	
	public final boolean containsKey(int key){return findCell(key)>=0;}
	
	public IntList get(int key){
		int cell=findCell(key);
		return cell<0 ? null : values[cell];
	}
	
	public IntList getOrCreate(int key){
		IntList list=get(key);
		if(list==null){list=new IntList(2); put(key, list);} // Create new list
		return list;
	}
	
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
	
	public void put(int key, int value){getOrCreate(key).add(value);}
	
	public void putAll(IntListHashMap map){
		for(int i=0; i<map.keyList.size; i++){
			int key=map.keyList.array[i];
			IntList list=map.get(key);
			if(list!=null){put(key, list.copy());} // Deep copy lists
		}
	}
	
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
	
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}
	
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
	
	public int[] toArray(){return keyList.toArray();}
	
	public int[] keys(){return keys;}
	
	public IntList[] values(){return values;}
	
	public int invalid(){return invalid;}
	
	public int size(){return size;}
	
	public boolean isEmpty(){return size==0;}
	
	private int[] keys;
	private IntList[] values;
	private IntList keyList;
	private int size=0;
	private int invalid;
	private int modulus;
	private int sizeLimit;
	private final float loadFactor;
	
	static final int MASK=Integer.MAX_VALUE;
	static final int MINMASK=Integer.MIN_VALUE;
	private static final int extra=10;
	private static final Random randy=shared.Shared.random(1);
}