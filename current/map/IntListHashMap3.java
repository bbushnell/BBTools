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
 * Thread-safety: Not thread-safe. External synchronization required for concurrent access.
 * Stored IntList values are borrowed: get() returns a live mutable list, and
 * put(int, IntList) stores the supplied list by reference.
 *
 * @author Brian Bushnell
 * @contributor Isla, Amber
 * @date June 2, 2025
 */
public final class IntListHashMap3 implements Serializable{

	/** Serialization version for this map layout. */
	private static final long serialVersionUID=1L;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an IntListHashMap3 with default initial capacity and load factor. */
	public IntListHashMap3(){this(256);}

	/** Creates an IntListHashMap3 with specified initial capacity and default load factor.
	 * @param initialSize Positive minimum logical table size */
	public IntListHashMap3(int initialSize){this(initialSize, 0.7f);}

	/**
	 * Creates an IntListHashMap3 with specified capacity and load factor, initializing
	 * the backing arrays and key tracker.
	 * @param initialSize Positive minimum logical table size
	 * @param loadFactor_ Requested fraction in (0,1), narrowed to float and clamped to [0.25,0.90]
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

	/*--------------------------------------------------------------*/
	/*----------------         Public Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Releases all list references without modifying the lists or shrinking capacity.
	 * Dense maps clear whole backing arrays. Sparse maps resolve all tracked keys
	 * before deleting any cell, reusing keyList as temporary cell-index storage. */
	public void clear(){
		if(size<1){return;}
		if(size>0.25*sizeLimit){
			Arrays.fill(keys, invalid);
			Arrays.fill(values, null);
		}else{
			final int limit=keyList.size;
			for(int i=0; i<limit; i++){
				final int key=keyList.array[i];
				final int cell=findCell(key);
				assert(cell>=0) : "keyList tracks every live key before IntListHashMap3.clear deletion; "+
					"key="+key+", index="+i+", tracked="+keyList.size+", size="+size;
				keyList.array[i]=cell;
			}
			for(int i=0; i<limit; i++){
				final int cell=keyList.array[i];
				keys[cell]=invalid; values[cell]=null; // Clear resolved cells
			}
		}
		keyList.clear();
		size=0;
	}

	/** Tests key occupancy, independently of whether its stored list is null.
	 * @param key Key to query
	 * @return True when key is present */
	public final boolean contains(int key){return findCell(key)>=0;}

	/** Alias for contains(int).
	 * @param key Key to query
	 * @return True when key is present */
	public final boolean containsKey(int key){return findCell(key)>=0;}

	/** Borrows the currently stored mutable list.
	 * @param key Key to query
	 * @return Stored list reference, or null for either absence or a stored null */
	public IntList get(int key){
		int cell=findCell(key);
		return cell<0 ? null : values[cell];
	}

	/** Supplies a mutable list for key, allocating and storing one if absent or null.
	 * @param key Key whose value is needed
	 * @return Borrowed stored list, non-null on return */
	public IntList getOrCreate(int key){
		IntList list=get(key);
		if(list==null){list=new IntList(2); put(key, list);} // Create new list
		return list;
	}

	/** Stores the supplied list reference without cloning or clearing either list.
	 * At the physical capacity ceiling, a new entry exceeding 95% occupancy is stored
	 * before failure is reported; existing-key replacements do not request growth.
	 * @param key Key to insert or replace; every int is supported
	 * @param value List reference to retain, or null
	 * @return Previous list, or null for absence or a previously stored null
	 * @throws IllegalStateException If a new entry exceeds the terminal 95% limit */
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

	/** Appends to the stored list, creating one for an absent or stored-null mapping.
	 * Existing external aliases to the list observe this append.
	 * @param key Key to insert or extend
	 * @param value Integer to append, including duplicates */
	public void put(int key, int value){getOrCreate(key).add(value);}

	/** Replaces mappings with independent copies of the source's non-null lists.
	 * Source-null values and destination-only mappings are skipped. A distinct source
	 * is preserved; a self-copy replaces each non-null list with its copy.
	 * @param map Non-null source stable during the call unless it is this map itself */
	public void putAll(IntListHashMap3 map){
		for(int i=0; i<map.keyList.size; i++){
			int key=map.keyList.array[i];
			IntList list=map.get(key);
			if(list!=null){put(key, list.copy());} // Deep copy lists
		}
	}

	/** Removes a mapping without clearing its former list; repairs the probe cluster.
	 * Scans the tracked-key list and swaps its last key into the removed position.
	 * @param key Key to remove
	 * @return True if the key was present, including when its stored value was null */
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

	/*--------------------------------------------------------------*/
	/*----------------         Private Methods      ----------------*/
	/*--------------------------------------------------------------*/

	/** After a removal blanks one cell, re-probe the rest of that cluster so no key is stranded behind
	 * the new gap (identical algorithm to IntListHashMap; rehashCell moves keys earlier in their probe
	 * sequence, possibly across the array end, and never past the cursor).
	 * @param initial Cell just emptied by removal */
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

	/** Moves a live key/list reference to an earlier available slot without changing size.
	 * @param cell Occupied cell index
	 * @return True when the entry moved and its old cell was cleared */
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

	/** Chooses an unused negative sentinel and rewrites empty keys, preserving live entries. */
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

	/** Mixes an integer, folding the top output range below the physical array ceiling.
	 * Based on the Thomas Wang/Bob Jenkins integer mixer; the cap keeps initial
	 * probes in range even when the logical power-of-two mask exceeds the physical length.
	 * @param key Integer to mix
	 * @return Nonnegative mixed hash strictly below 0x7FFFF800 */
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
	 * initial=hash(key)&amp;mask is IDENTICAL to findCellOrEmpty's, so a key placed on insert is always
	 * found here. Stops at the first invalid cell.
	 * @param key Key to query
	 * @return Matching cell index, or -1 when absent */
	int findCell(final int key){
		if(key==invalid){return -1;}
		final int limit=keys.length;

		// Use mask for fast modulo (Power of 2)
		final int initial=hash(key)&mask;

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

	/** Finds the key or first empty cell using the same probe sequence as lookup.
	 * @param key Key unequal to the current empty-cell sentinel
	 * @return Matching or empty cell index
	 * @throws RuntimeException If no matching or empty cell exists */
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		final int limit=keys.length;
		final int initial=hash(key)&mask;

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

	/** Requests a larger table after the size threshold is reached. */
	private final void resize(){
		assert(size>=sizeLimit);
		//Exclude slack before doubling to avoid fourfold power-of-two rounding.
		//Fixed [map/IntListHashMap3#001], shared with [map/IntLongHashMap2#002].
		resize((keys.length-extra)*2L);
	}

	/** Rebuilds arrays at a rounded power-of-two logical size plus slack, capped physically.
	 * Retains value-list references and rebuilds keyList in old-table order.
	 * At the physical ceiling, permits up to 95% occupied cells without reallocating.
	 * @param size2 Positive requested logical size; must exceed entry count when growing
	 * @throws IllegalStateException If storage cannot grow and occupancy exceeds 95% */
	private final void resize(final long size2){
		assert(size2>0) : "A positive logical capacity is needed for the hash mask: "+size2;

		// Round up to Power of 2
		long size3=Long.highestOneBit(size2);
		if(size3<size2){size3<<=1;}

		final int newMask=(int)(size3-1);
		size3=Math.min(size3+extra, Shared.SAFE_ARRAY_LEN);
		if(keys!=null && size3<=keys.length){
			//Retain the table and probe mask when growth is capped; leave empty cells for probing.
			sizeLimit=(int)(keys.length*95L/100);
			if(size>sizeLimit){
				throw new IllegalStateException("IntListHashMap3 cannot grow beyond "+keys.length+
					" cells; 95% occupancy limit="+sizeLimit+", entries="+size);
			}
			return;
		}
		assert(size2>size) : "Growing storage must fit all existing map entries: "+size+", "+size2;
		mask=newMask;

		// Load factor logic adapted for new sizing
		sizeLimit=(size3<Shared.SAFE_ARRAY_LEN ? (int)((size3-extra)*loadFactor) : (int)(size3*95/100));

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

	/*--------------------------------------------------------------*/
	/*----------------           Accessors          ----------------*/
	/*--------------------------------------------------------------*/

	/** Copies valid keys in tracked-key order, which removal and resizing can change.
	 * @return Newly allocated dense key array without sentinels */
	public int[] toArray(){return keyList.toArray();}

	/** Borrows mutable key storage; modifying it can break hashing and key tracking.
	 * Resizing replaces the array; clear() mutates it in place.
	 * @return Map-owned key array with invalid() in empty cells */
	public int[] keys(){return keys;}

	/** Borrows value storage parallel to keys(); null alone does not establish absence.
	 * Resizing replaces the array; clear() releases its references but does not clear lists.
	 * @return Map-owned array of borrowed list references */
	public IntList[] values(){return values;}

	/** Returns the current empty-cell marker, which changes if inserted as a real key.
	 * @return Negative sentinel absent from the key set */
	public int invalid(){return invalid;}

	/** Returns the number of mapped keys, including stored-null mappings.
	 * @return Entry count */
	public int size(){return size;}

	/** Tests whether there are no mapped keys.
	 * @return True when size() is zero */
	public boolean isEmpty(){return size==0;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Key storage, with invalid in empty cells. */
	private int[] keys;
	/** List references parallel to keys; empty cells contain null. */
	private IntList[] values;
	/** Each live key once; also temporary cell-index scratch during sparse clear. */
	private IntList keyList;
	/** Number of occupied key cells. */
	private int size=0;
	/** Negative empty-cell marker, replaced when that key is inserted. */
	private int invalid;
	/** Logical power-of-two size minus one, retained when physical length is capped. */
	private int mask;
	/** Entry count above which insertion requests growth. */
	private int sizeLimit;
	/** Requested load fraction narrowed to float and clamped to [0.25,0.90]. */
	private final float loadFactor;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Mask for the nonnegative portion of a key. */
	static final int MASK=Integer.MAX_VALUE;
	/** Sign-bit mask used when choosing an empty-cell sentinel. */
	static final int MINMASK=Integer.MIN_VALUE;
	/** Extra probe cells beyond the logical table size. */
	private static final int extra=10;
	/** Seeded generator for replacement empty-cell markers. */
	private static final Random randy=shared.Shared.random(1);
}
