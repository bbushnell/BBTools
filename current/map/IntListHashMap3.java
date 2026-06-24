package map;

import java.io.Serializable;
import java.util.Arrays;
import shared.Random;

import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.IntList;

/**
 * A specialized hash map that maps integer keys to IntList values.
 * Uses open addressing with linear probing for collision resolution.
 * Optimized for scenarios where each key maps to a collection of integers.
 * Uses 'hash32plus' mixer and Power-of-2 sizing for performance,
 * similar to IntHashMap3.
 *
 * @author Brian Bushnell
 * @contributor Isla, Amber
 * @date June 2, 2025
 */
public final class IntListHashMap3 implements Serializable {
	
	private static final long serialVersionUID = 1L;
	
	/** Creates an IntListHashMap3 with default initial capacity and load factor. */
	public IntListHashMap3(){this(256);}
	
	/** Creates an IntListHashMap3 with specified initial capacity and default load factor.
	 * @param initialSize Initial hash table capacity */
	public IntListHashMap3(int initialSize){this(initialSize, 0.7f);}
	
	/**
	 * Creates an IntListHashMap3 with specified capacity and load factor, initializing
	 * the backing arrays and key tracker.
	 * @param initialSize Initial capacity for the hash table
	 * @param loadFactor_ Maximum load factor before resizing
	 */
	public IntListHashMap3(int initialSize, double loadFactor_){
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
	
	public void putAll(IntListHashMap3 map){
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
	
	/** After a removal blanks one cell, re-probe the rest of that cluster so no key is stranded behind
	 * the new gap (identical algorithm to IntListHashMap; verified: rehashCell only moves keys backward
	 * to dest<=cell, the just-vacated cell is never misread as the stop sentinel). */
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
	
	//CLEVER [verified]: keys may be any int, so a real key can equal the negative sentinel; on that
	//insert this regenerates the sentinel and rewrites every old marker, keeping "empty" distinct from
	//data without reserving a key value (same trick as IntHashSet / the IntListHashMap twin).
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
	
	/** Thomas Wang / Bob Jenkins integer hash mixer, adapted for Java.
	 * @return a non-negative int strictly below 0x7FFFF800 (&lt; Shared.SAFE_ARRAY_LEN).
	 * CLEVER [verified]: the final ternary folds the top range [0x7FFFF800, 0x7FFFFFFF] back down via
	 * (key-0x7FFFF800)*64, so the output is always a safe array index -> hash(key)&mask is in range for
	 * ANY mask, even mask=-1 (same capped-hash guarantee as Tools.hash32plus; mask never needs capping). */
	private static final int hash(int key){
		key=~key+(key<<15);
		key=key^(key>>>12);
		key=key+(key<<2);
		key=key^(key>>>4);
		key=key*2057;
		key=key^(key>>>16);
		key=key&(0x7FFFFFFF);
		return key<0x7FFFF800 ? key : (key-0x7FFFF800)*64;
	}

	/** Probe for an existing key; -1 if absent (or if key==the current invalid sentinel).
	 * initial=hash(key)&mask is IDENTICAL to findCellOrEmpty's, so a key placed on insert is always
	 * found here (probe-consistency verified). Stops at the first invalid cell (open-addressing). */
	int findCell(final int key){
		if(key==invalid){return -1;}
		final int limit=keys.length;

		// Use mask for fast modulo (Power of 2)
		final int initial=hash(key) & mask;
		
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
	
	//CLEVER [verified]: initial=hash(key)&mask, byte-identical to findCell's. The two MUST agree or an
	//inserted key becomes unfindable; they do. Returns key's cell or the first empty cell.
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		final int limit=keys.length;
		final int initial=hash(key) & mask;
		
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
		//grow 2x: double the LOGICAL pow2 capacity (keys.length-extra==pow2), not keys.length (=pow2+extra).
		//keys.length*2L overshot the power-of-2 boundary: highestOneBit(2*pow2+20)=2*pow2, then the
		//size3<size2 round-up doubled AGAIN -> 4*pow2. Fixed to true 2x growth.
		//TODO: Possible bug [map/IntListHashMap3#001] - shared pow2-family over-allocation, anchor
		//[map/IntLongHashMap2#002] (Brian-greenlit family fix). Output unchanged; ~2x less memory.
		resize((keys.length-extra)*2L);
	}
	
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		
		// Round up to Power of 2
		long size3=Long.highestOneBit(size2);
		if(size3<size2){size3<<=1;}
		
		mask=(int)(size3-1); // Mask is length-1
		size3=Math.min(size3+extra, Shared.SAFE_ARRAY_LEN);
		
		// Load factor logic adapted for new sizing
		final float loadFactor2=(size3<Shared.SAFE_ARRAY_LEN ? loadFactor : Tools.max(loadFactor, 0.85f));
		sizeLimit=(int)((size3-extra)*loadFactor2);
		
		final int[] oldK=keys; final IntList[] oldV=values;
		keys=KillSwitch.allocInt1D((int)size3); values=new IntList[(int)size3];
		Arrays.fill(keys, invalid); // Initialize with invalid markers
		
		if(size<1){return;}
		
		size=0; keyList.clear(); // Reset for rehashing
		for(int i=0; i<oldK.length; i++){
			final int k=oldK[i]; final IntList v=oldV[i];
			if(k!=invalid){put(k, v);} // Rehash all valid entries
		}
	}
	
	/** Valid keys only, no sentinels (snapshot of keyList). Use this, not keys(), to iterate keys. */
	public int[] toArray(){return keyList.toArray();}

	/** Raw backing array - contains the invalid() sentinel in empty cells. Callers must filter on
	 * invalid() (or use toArray()); index alignment with values() is internal. No caller reads this today. */
	public int[] keys(){return keys;}

	/** Raw backing array - empty cells are null and live cells align by index with keys(). No caller reads this today. */
	public IntList[] values(){return values;}
	
	public int invalid(){return invalid;}
	
	public int size(){return size;}
	
	public boolean isEmpty(){return size==0;}
	
	private int[] keys;
	private IntList[] values;
	private IntList keyList;
	private int size=0;
	private int invalid;
	private int mask; // Replaces modulus
	private int sizeLimit;
	private final float loadFactor;
	
	static final int MASK=Integer.MAX_VALUE;
	static final int MINMASK=Integer.MIN_VALUE;
	private static final int extra=10;
	private static final Random randy=shared.Shared.random(1);
}