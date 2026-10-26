package map;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import shared.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.LongList;

/**
 * Hash map that maps long keys to lists of long values.
 * Uses open addressing with linear probing for collision resolution and is optimized for primitive longs to avoid boxing overhead.
 * Value lists returned by get() and values() are live internal lists and may be cleared or reused after remove() or clear().
 * This class is not thread-safe; use external synchronization around any overlapping reads and writes.
 * @author Brian Bushnell
 * @date August 29, 2023
 */
public final class LongLongListHashMap{

	/**
	 * Legacy ad hoc stress benchmark against {@link HashMap}.
	 * @param args ignored
	 */
	public static void main(String[] args){
		Random randy2=Shared.threadLocalRandom();
		LongLongListHashMap map=new LongLongListHashMap(20, 0.7f);
		HashMap<Long, Integer> map2=new HashMap<Long, Integer>(20, 0.7f);
		ArrayList<Long> list=new ArrayList<Long>();
		ArrayList<Long> list2=new ArrayList<Long>();
//		ArrayList<Integer> vals=new ArrayList<Integer>();
		for(long i=0; i<1000; i++){
			assert(!map.contains(i));
			assert(!map2.containsKey(i));
			list.add(Long.valueOf(i));
		}
		for(int i=0; i<1000; i++){
			long r=randy2.nextLong();
			list2.add(r);
		}

		for(long x : list){
			map.put(x, (int)(2*x));
			map2.put(x, (int)(2*x));
			assert(map.get(x).contains((int)(2*x)));
			assert(map2.get(x)==((int)(2*x)));
		}

		for(long x : list){
			assert(map.get(x).contains((int)(2*x)));
			assert(map2.get(x)==((int)(2*x)));
			map.remove(x);
			map2.remove(x);
			assert(!map.contains(x));
			assert(!map2.containsKey(x));
		}
		assert(map.isEmpty());
		assert(map2.isEmpty());

		for(long x : list2){
			map.put(x, (int)(2*x));
			map2.put(x, (int)(2*x));
			assert(map.get(x).contains((int)(2*x)));
			assert(map2.get(x)==((int)(2*x)));
		}

		for(long x : list2){
			assert(map.get(x).contains((int)(2*x)));
			assert(map2.get(x)==((int)(2*x)));
			map.remove(x);
			map2.remove(x);
			assert(!map.contains(x));
			assert(!map2.containsKey(x));
		}
		assert(map.isEmpty());
		assert(map2.isEmpty());

		int count=2000000;
		int runs=16;
		LongList ll=new LongList(count);
		for(int i=0; i<count; i++){ll.add(randy2.nextLong());}

		Shared.printMemory();
		Timer t=new Timer();
		for(int k=0; k<2; k++){
			System.err.println("LongLongListHashMap:");
			t.start();
			for(int i=0; i<runs; i++){
//				for(long x : ll.array){
//					map.add(x);
//				}
				final long[] y=ll.array;
				for(int z=0; z<count; z++){
					final long key=y[z];
					map.put(key, key);
					map.contains(key);
					map.remove(key);
					map.put(key, key);
				}
//				for(long x : ll.array){
//					map.remove(x);
//				}
//				map.clear();
//				assert(map.isEmpty());
//				System.err.println("Finished run "+i);
			}
			t.stop();
			System.err.println(t);
			Shared.printMemory();

			System.err.println("HashMap:");
			t.start();
			for(int i=0; i<runs; i++){
				for(long x : ll.array){
					map2.put(x, (int)x);
					map2.containsKey(x);
					map2.remove(x);
					map2.put(x, (int)x);
				}
//				assert(map2.isEmpty());
//				System.err.println("Finished run "+i);
			}
			t.stop();
			System.err.println(t);
			Shared.printMemory();
		}
		t.stop();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs a map with the default requested capacity and load factor.
	 */
	public LongLongListHashMap(){
		this(256);
	}

	/**
	 * Constructs a map with the default load factor.
	 * @param initialSize requested minimum prime modulus before adding slack cells;
	 * asserted positive
	 */
	public LongLongListHashMap(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Constructs a map with a caller-supplied initial modulus and load factor.
	 * @param initialSize requested minimum prime modulus before adding slack cells;
	 * asserted positive
	 * @param loadFactor_ requested resize threshold fraction; asserted in (0,1)
	 * before clamping to [0.25, 0.90]
	 */
	public LongLongListHashMap(int initialSize, float loadFactor_){
		invalid=randy.nextLong()|MINMASK;
		assert(invalid<0);
		assert(initialSize>0);
		assert(loadFactor_>0 && loadFactor_<1);
		loadFactor=Tools.mid(0.25f, loadFactor_, 0.90f);
		resize(initialSize);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Removes all mappings. Large value lists are released; smaller lists may be retained empty for reuse. */
	public void clear(){
		if(size<1){return;}
		for(int i=0; i<values.length /*&& size>0*/; i++){//TODO: Enable the early exit
			if(keys[i]!=invalid){
				keys[i]=invalid;
				values[i].clear();
				if(values[i].capacity()>10){values[i]=null;}//Release large lists; retain small lists as empty recycle slots.
				size--;
			}
			assert(values[i]==null || values[i].isEmpty());
		}
		assert(size==0);
		size=0;
//		assert(verify()); //123
	}

	/**
	 * Tests whether key is present in the table.
	 * @param key key to search for
	 * @return true if key is mapped
	 */
	public boolean contains(long key){
//		assert(verify()); //123
		return key==invalid ? false : findCell(key)>=0;
	}

	/**
	 * Alias for contains(long), matching java.util.Map naming.
	 * @param key key to search for
	 * @return true if key is mapped
	 */
	public boolean containsKey(long key){
		return contains(key);
	}

	/**
	 * Returns the live value list mapped from key.
	 * @param key key to search for
	 * @return Borrowed value list, or null when absent. Growth retains its identity;
	 * removal, clear or putCopyFree replacement may end its association with this key
	 */
	public LongList get(long key){
//		assert(verify()); //123
		final int cell=findCell(key);
		return cell<0 ? null : values[cell];
	}

	/**
	 * Adds a value to the list mapped by the specified key.
	 * Creates a new list if the key is not already present.
	 * @param key Key to map the value to
	 * @param value Value to add to the key's list
	 * @return true if the key was newly added, false if key already existed
	 */
	public boolean put(long key, long value){
//		assert(verify()); //123
		if(key==invalid){resetInvalid();}
		int cell=findCellOrEmpty(key);
		if(keys[cell]==invalid){ //new key
			keys[cell]=key;
			//reuse a recycled empty list if one is parked here (clear()/remove() keep small lists); else alloc
			if(values[cell]==null){
				values[cell]=new LongList(2);
			}
			values[cell].add(value);
			size++;
			if(size>sizeLimit){resize();}
			return true;
		}
		values[cell].add(value); //existing key: append to its list
		return false;
	}

	/**
	 * Adds all values from the specified list to the key's mapped list.
	 * Creates a new live internal list if the key is not already present and
	 * copies the caller-owned list into it; an empty input creates a present key
	 * with an empty list.
	 * @param key Key to map the values to
	 * @param value List of values to copy into the key's list
	 * @return true if the key was newly added, false if key already existed
	 */
	public boolean put(long key, LongList value){
//		assert(verify()); //123
		if(key==invalid){resetInvalid();}
		int cell=findCellOrEmpty(key);
		if(keys[cell]==invalid){
			keys[cell]=key;
			if(values[cell]==null){
				values[cell]=new LongList(2);
			}
			values[cell].addAll(value);
			size++;
			if(size>sizeLimit){resize();}
			return true;
		}
		values[cell].addAll(value);
		return false;
	}

	/** Adopts a list without copying elements, replacing rather than appending.
	 * The caller transfers exclusive ownership: do not independently reuse or modify
	 * the list, or share it with another key/map. Later append, remove and clear may
	 * mutate or clear it; a small cleared list may be recycled for a different key.
	 * Replacement detaches the old list without clearing it. Empty lists are permitted
	 * for compatibility with bulk put, but verify() rejects occupied empty lists.
	 * @param key Key whose mapping to store or replace
	 * @param value Non-null list to adopt; its identity is preserved across map growth
	 * @return True if this key was absent, false if its mapping was replaced
	 * @throws NullPointerException If value is null, before modifying this map */
	public final boolean putCopyFree(long key, final LongList value){
		if(value==null){throw new NullPointerException("Cannot adopt a null value list for key "+key);}
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		assert(keys[cell]==invalid || keys[cell]==key) :
			"findCellOrEmpty must select a matching or empty cell before transferring list ownership: "+key+", "+cell;
		values[cell]=value;
		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
			return true;
		}
		return false;
	}

	/**
	 * Removes the specified key and its associated list from the map.
	 * Clears the list; large backing arrays are released while small empty lists
	 * may be retained for reuse.
	 * @param key Key to remove
	 * @return true if the key was present and removed, false otherwise
	 */
	public boolean remove(long key){
//		assert(verify()); //123
		if(key==invalid){return false;}
		final int cell=findCell(key);
		if(cell<0){return false;}
		assert(keys[cell]==key);
		keys[cell]=invalid;
		values[cell].clear();
		if(values[cell].capacity()>=10){values[cell]=null;}
		size--;

		rehashFrom(cell);
//		assert(verify()); //123
		return true;
	}

	/**
	 * Returns the number of mapped keys.
	 * @return mapped-key count
	 */
	public int size(){return size;}

	/**
	 * Returns true when no keys are mapped.
	 * @return true if empty
	 */
	public boolean isEmpty(){return size==0;}

	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns a string representation showing all keys in this map.
	 * @return list-view string
	 */
	@Override
	public String toString(){
		return toStringListView();
	}

	/**
	 * Returns a list of occupied cells with cell index, key and value list.
	 * @return set-view string
	 */
	public String toStringSetView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<keys.length; i++){
			if(keys[i]!=invalid){
				sb.append(comma+"("+i+", "+keys[i]+", "+values[i]+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}

	/**
	 * Returns a list of mapped keys.
	 * @return list-view string
	 */
	public String toStringListView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<keys.length; i++){
			if(keys[i]!=invalid){
				sb.append(comma+keys[i]);
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}

	/**
	 * Returns a newly allocated array of all mapped keys in backing-array order.
	 * @return mapped keys
	 */
	public long[] toArray(){
		long[] x=KillSwitch.allocLong1D(size);
		int i=0;
		for(long key : keys){
			if(key!=invalid){
				x[i]=key;
				i++;
			}
		}
		return x;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** repOK-style dev invariant checker (assert(verify()) pattern; call sites currently //123-disabled).
	 * Requires a non-empty value list for every occupied key. It looks up
	 * findCell(key) and requires it to equal the stored index i, requires empty cells to have null/empty
	 * lists, and checks numValues==numFound==size.
	 * @return true when size, empty-cell and non-empty-list invariants hold */
	public boolean verify(){
		//Note: This method rejects empty lists; an occupied key must map to actual values.
		if(keys==null){return true;}
		int numValues=0;
		int numFound=0;
		for(int i=0; i<keys.length; i++){
			final long key=keys[i];
			final LongList value=values[i];

			if(key==invalid){
				if(value!=null && !value.isEmpty()){
					assert(false) : i+", "+key+", "+value;
					return false;
				}
			}else{
				numValues++;
				if(value==null || value.isEmpty()){
					assert(false) : i+", "+key+", "+value;
					return false;
				}
				final int cell=findCell(key);
				if(i==cell){
					numFound++;
				}else{
					assert(false) : i+", "+key+", "+value+", "+cell+"\n"+((cell>=0) ? keys[cell]+", "+values[cell]+"\n" : "");
					return false;
				}
			}
		}
		boolean pass=(numValues==numFound && numValues==size);
		assert(pass) : numValues+", "+numFound+", "+size;
		return pass;
	}

	/**
	 * Rehashes the contiguous probe chain after a removed cell.
	 * @param initial just-emptied cell that begins the possible chain break
	 */
	private void rehashFrom(int initial){
		if(size<1){return;}
		final int limit=keys.length;
		for(int cell=initial+1; cell<limit; cell++){
			final long x=keys[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final long x=keys[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
	}

	//CLEVER [verified]: unlike IntListHashMap (which nulls the source on a move), this SWAPS the value
	//slots: destCell is an empty cell whose values[destCell] is a recycled empty/null LongList, so moving
	//it back to the now-empty sourceCell preserves the recycled-list pool instead of orphaning it. After
	//the swap, occupied cells still hold non-empty lists and empty cells still hold null/empty lists.
	/**
	 * Moves the key and its value list to its earliest available probe-chain slot.
	 * @param sourceCell occupied cell to test
	 * @return true if the key moved
	 */
	private boolean rehashCell(final int sourceCell){
		final long key=keys[sourceCell];
		final LongList value=values[sourceCell];
		assert(key!=invalid);
		if(key==invalid){resetInvalid();}
		final int destCell=findCellOrEmpty(key);
		if(sourceCell==destCell){return false;}
		assert(keys[destCell]==invalid); //destCell empty -> its list is recycled (null or empty)
		keys[sourceCell]=invalid;
		values[sourceCell]=values[destCell]; //recycled empty/null list moves to the vacated cell
		keys[destCell]=key;
		values[destCell]=value;
		return true;
	}

	/** Chooses a new negative empty-cell sentinel and rewrites old empty cells. */
	private void resetInvalid(){
		final long old=invalid;
		long x=invalid;
		while(x==old || contains(x)){x=randy.nextLong()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<keys.length; i++){
			if(keys[i]==old){
				assert(values[i]==null || values[i].isEmpty());
				keys[i]=invalid;
			}
		}
	}

	/**
	 * Finds a key by linear probing from its hashed home cell.
	 * @param key key to find
	 * @return occupied cell index, or -1 if key is absent
	 */
	private int findCell(final long key){
		if(key==invalid){return -1;}

		final int limit=keys.length, initial=(int)((key&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final long x=keys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		for(int cell=0; cell<initial; cell++){
			final long x=keys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		return -1;
	}

	//CLEVER [verified]: initial index here ((key&MASK)%modulus) is byte-identical to findCell's, so a key
	//placed on insert is always found later (probe-consistency). MASK==Long.MAX_VALUE strips the sign bit
	//so negative keys map into [0,modulus); a lone divergence between the two would be the bug. None here.
	/**
	 * Finds a key or the first empty cell in its probe chain.
	 * @param key key to find or insert; must not equal invalid()
	 * @return matching or empty cell index
	 */
	private int findCellOrEmpty(final long key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";

		final int limit=keys.length, initial=(int)((key&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final long x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final long x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}

	//Prime-modulus growth ~2x via primeAtLeast(keys.length*2+1) -- NOT the pow2 family's 4x bug (no
	//highestOneBit round-up here). Move already-owned lists into the new empty table;
	//ordinary public put still copies caller lists.
	/** Expands the table and rehashes mapped keys into the new arrays. */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}

	/**
	 * Rehashes all mapped keys into a prime-sized table, moving owned lists without copying elements.
	 * @param size2 requested minimum modulus before adding slack cells
	 */
	private final void resize(final long size2){
//		assert(verify()); //123
		assert(size2>size) : size+", "+size2;
		long newPrime=Primes.primeAtLeast(size2);
		if(newPrime+extra>Integer.MAX_VALUE){
			newPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra);
		}
		assert(newPrime>modulus) : "Overflow: "+size+", "+size2+", "+modulus+", "+newPrime;
		modulus=(int)newPrime;

		final int size3=(int)(newPrime+extra);
		sizeLimit=(int)(modulus*loadFactor);
		final long[] oldKeys=keys;
		final LongList[] oldValues=values;
		keys=KillSwitch.allocLong1D(size3);
		values=KillSwitch.allocObject1D(size3, LongList.class);
		Arrays.fill(keys, invalid);

//		System.err.println("Resizing "+(old==null ? "null" : ""+old.length)+" to "+size3);

		if(size<1){return;}

		size=0;
		for(int i=0; i<oldKeys.length; i++){
			long key=oldKeys[i];
			if(key!=invalid){
				putCopyFree(key, oldValues[i]);
			}
		}
//		assert(verify()); //123
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed backing array - empty cells hold the invalid() sentinel; filter on invalid() (or use
	 * toArray()). Index-aligned with values().
	 * @return borrowed backing key array */
	public long[] keys(){return keys;}

	/** Borrowed backing array - empty cells are null OR a retained empty LongList; live cells align by index
	 * with keys(). Lists are borrowed and may be cleared or reused after remove() or clear().
	 * @return borrowed backing value-list array */
	public LongList[] values(){return values;}

	/**
	 * Returns the current negative empty-cell sentinel.
	 * @return empty-cell sentinel
	 */
	public long invalid(){return invalid;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Keys array; empty cells hold invalid. */
	private long[] keys;
	/** Value-list array parallel to keys; empty cells hold null or an empty recycled list. */
	private LongList[] values;
//	private LongList keyList;//TODO: Store keys here for quick clearing, perhaps
	/** Number of occupied cells. */
	private int size=0;
	/** Sentinel value representing empty cells in the hash table */
	private long invalid;
	/** Prime modulus used for the initial probe cell. */
	private int modulus;
	/** Resize threshold in occupied cells. */
	private int sizeLimit;
	/** Effective load factor used to compute sizeLimit. */
	private final float loadFactor;

	/** Deterministic generator for negative empty-cell sentinels. */
	private static final Random randy=shared.Shared.random(1);
	/** Mask that clears the sign bit before the modulus operation. */
	private static final long MASK=Long.MAX_VALUE;
	/** Mask that sets the sign bit when choosing empty-cell sentinels. */
	private static final long MINMASK=Long.MIN_VALUE;

	/** Extra cells after the prime modulus to simplify wrapped probing. */
	private static final int extra=10;

}
