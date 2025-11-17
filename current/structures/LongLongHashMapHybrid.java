package structures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * This class should allow mapping a long to one or more values.
 * For a single value, it will act as a LongLongHashMap.
 * For multiple values, it will act more like a LongLongListHashMap.
 * The primary value stored when there are multiple values will be the
 * (-index-OFFSET) of the list's index.
 * As such, this does NOT support negative values, though it could be
 * modified to support most negative values, by making OFFSET large.
 * However, that makes the logic of determining whether a key is present
 * from the return value more confusing.
 * @author Brian Bushnell
 * @date November 13, 2024
 *
 */
public final class LongLongHashMapHybrid{
	
	/** Test method demonstrating functionality and performance benchmarks.
	 * @param args Command-line arguments (unused) */
	public static void main(String[] args){
		Random randy2=Shared.threadLocalRandom();
		LongLongHashMapHybrid map=new LongLongHashMapHybrid(20, 0.7f);
		HashMap<Long, Long> map2=new HashMap<Long, Long>(20, 0.7f);
		ArrayList<Long> list=new ArrayList<Long>();
		ArrayList<Long> list2=new ArrayList<Long>();
//		ArrayList<Integer> vals=new ArrayList<Integer>();
		for(long i=0; i<1000; i++){
			assert(!map.contains(i));
			assert(!map2.containsKey(i));
			list.add(Long.valueOf(i));
		}
		for(int i=0; i<1000; i++){
			long r=randy2.nextLong(Long.MAX_VALUE/4);
			list2.add(r);
		}
		
		for(long x : list){
			map.put(x, (2*x));
			map2.put(x, (2*x));
			assert(map.get(x)==((2*x)));
			assert(map2.get(x)==((2*x)));
		}
		
		for(long x : list){
			assert(map.get(x)==((2*x)));
			assert(map2.get(x)==((2*x)));
			map.remove(x);
			map2.remove(x);
			assert(!map.contains(x));
			assert(!map2.containsKey(x));
		}
		assert(map.isEmpty());
		assert(map2.isEmpty());
		
		for(long x : list2){
			map.put(x, (2*x));
			map2.put(x, (2*x));
			assert(map.get(x)==((2*x)));
			assert(map2.get(x)==((2*x)));
		}
		
		for(long x : list2){
			assert(map.get(x)==((2*x)));
			assert(map2.get(x)==((2*x)));
			map.remove(x);
			map2.remove(x);
			assert(!map.contains(x));
			assert(!map2.containsKey(x));
		}
		assert(map.isEmpty());
		assert(map2.isEmpty());
		
		int count=4000000;
		int runs=32;
		LongList ll=new LongList(count);
		for(int i=0; i<count; i++){ll.add(randy2.nextLong());}

		Shared.printMemory();
		Timer t=new Timer();
		for(int k=0; k<2; k++){
			System.err.println("LongLongHashMapHbrid:");
			t.start();
			for(int i=0; i<runs; i++){
//				for(long x : ll.array){
//					map.add(x);
//				}
				final long[] y=ll.array;
				for(int z=0; z<count; z++){
					final long key=y[z];
					final long value=key&Long.MAX_VALUE;
					map.put(key, value);
					assert(map.contains(key));
					map.remove(key);
					assert(!map.contains(key));
					map.put(key, value);
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
			
//			System.err.println("HashMap:");
//			t.start();
//			for(int i=0; i<runs; i++){
//				for(long x : ll.array){
//					map2.add(x);
//				}
//				for(long x : ll.array){
//					map2.remove(x);
//				}
//				assert(map2.isEmpty());
////				System.err.println("Finished run "+i);
//			}
//			t.stop();
//			System.err.println(t);
//			Shared.printMemory();
		}
		t.stop();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a hash map with default initial size of 256 and load factor 0.7.
	 */
	public LongLongHashMapHybrid(){
		this(256);
	}
	
	/** Constructs a hash map with specified initial size and default load factor 0.7.
	 * @param initialSize Initial capacity of the hash map */
	public LongLongHashMapHybrid(int initialSize){
		this(initialSize, 0.7f);
	}
	
	/**
	 * Constructs a hash map with specified initial size and load factor.
	 * Generates a random invalid key marker for empty cells.
	 * @param initialSize Initial capacity of the hash map
	 * @param loadFactor_ Load factor between 0.25 and 0.90
	 */
	public LongLongHashMapHybrid(int initialSize, float loadFactor_){
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
	
	/** Removes all key-value mappings from the hash map.
	 * Resets all keys to invalid marker and values to zero. */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
		size=0;
//		assert(verify()); //123
	}
	
	/**
	 * Tests if the specified key exists in the hash map.
	 * @param key The key to search for
	 * @return true if the key is present, false otherwise
	 */
	public boolean contains(long key){
//		assert(verify()); //123
		return key==invalid ? false : findCell(key)>=0;
	}
	
	/**
	 * Tests if the specified key exists in the hash map.
	 * Alias for contains method.
	 * @param key The key to search for
	 * @return true if the key is present, false otherwise
	 */
	public boolean containsKey(long key){
		return contains(key);
	}
	
	/**
	 * Returns the value or list index for the specified key.
	 * For single values, returns the value directly.
	 * For multiple values, returns negative index (-index-OFFSET) to the list.
	 *
	 * @param key The key to look up
	 * @return The value, list index code, or NOTPRESENT if key not found
	 */
	public long get(long key){
//		assert(verify()); //123
		long value=NOTPRESENT;
		if(key!=invalid){
			int cell=findCell(key);
			if(cell>=0){value=values[cell];}
		}
		return value;
	}
	
	/**
	 * Retrieves the multi-value list from a negative list index code.
	 * @param code Negative index code from get() method
	 * @return The LongList containing multiple values
	 */
	public LongList getListFromCode(long code){
		assert(code<=-OFFSET) : code;
		return multivalues.get((int)(-code-OFFSET));
	}
	
	/**
	 * Fills the buffer with all values for the specified key.
	 * For single values, adds the value to buffer.
	 * For multiple values, returns the existing list directly.
	 *
	 * @param key The key to look up
	 * @param buffer Buffer to fill with values (cleared first)
	 * @return The buffer (for single values) or the multi-value list
	 */
	public LongList getOrFill(long key, LongList buffer){
//		assert(buffer.isEmpty());
		buffer.clear();
		long code=get(key);
		if(code==NOTPRESENT) {return buffer;}
		if(code>=0) {buffer.add(code); return buffer;}
		return multivalues.get((int)(-code-OFFSET));
	}
	
	/**
	 * Fills the buffer with all values for the specified key.
	 * Always uses the provided buffer, copying from multi-value lists if needed.
	 *
	 * @param key The key to look up
	 * @param buffer Buffer to fill with values (cleared first)
	 * @return The original get() return value or NOTPRESENT if not found
	 */
	public long fill(long key, LongList buffer){
//		assert(buffer.isEmpty());
		buffer.clear();
		long code=get(key);
		if(code==NOTPRESENT) {return code;}
		if(code>=0) {buffer.add(code);}
		else {buffer.addAll(multivalues.get((int)(-code-OFFSET)));}
		return code;
	}

	/**
	 * Map this key to value.
	 * @param key
	 * @param value
	 * @return true if the value was added, false if it was already contained.
	 */
	public boolean put(long key, long value){
		assert(value>=0) : "Unsupported negative value "+value;
		return putInner(key, value);
	}
	
	/**
	 * Map this key to value.
	 * @param key
	 * @param value
	 * @return true if the value was added, false if it was already contained.
	 */
	public boolean putInner(long key, long value){
//		assert(verify()); //123
		if(key==invalid){resetInvalid();}
		int cell=findCellOrEmpty(key);
		if(keys[cell]==invalid){
			keys[cell]=key;
			values[cell]=value;
			size++;
			if(size>sizeLimit){resize();}
//			assert(verify()); //123
			return true;
		}
		assert(keys[cell]==key);
		final long vCell=values[cell];
		if(vCell==value) {return false;}
		if(vCell>=0) {
			int listIndex=multivalues.size();
			LongList list=new LongList(4);
			multivalues.add(list);
			list.add(values[cell]);
			list.add(value);
			values[cell]=-listIndex-OFFSET;
			return true;
		}
		int listIndex=(int)(-vCell-OFFSET);
		LongList list=multivalues.get(listIndex);
//		if(list.contains(value)) {return false;}//Slow and not really needed for indexing.
		list.add(value);
//		assert(verify()); //123
		return true;
	}
	
	/**
	 * Remove this key from the map.
	 * @param key
	 * @return Old value.
	 */
	public long remove(long key){
		//This operation is difficult when the key has multiple values
		//It can be done, but will leave a hole (null list)
		//No reason to do it when used for indexing anyway
//		throw new RuntimeException("Unimplemented");
//		assert(verify());
		if(key==invalid){return NOTPRESENT;}
		final int cell=findCell(key);
		if(cell<0){return NOTPRESENT;}
		assert(keys[cell]==key || keys[cell]<NOTPRESENT);
		keys[cell]=invalid;
		final long value=values[cell];
		values[cell]=0;
		size--;
		if(value<0) {
			int idx=(int)(-value-OFFSET);
			multivalues.set(idx, null);
			if(idx+1==multivalues.size()) {multivalues.remove(idx);}
		}
		
		rehashFrom(cell);
//		assert(verify());
		return value;
	}
	
	/** Returns the number of key-value mappings in the hash map */
	public int size(){return size;}
	
	/** Returns true if the hash map contains no key-value mappings */
	public boolean isEmpty(){return size==0;}
	
	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	/**
	 * Returns detailed string showing array indices, keys, and values.
	 * Format: [(index, key, value), ...] for debugging purposes.
	 * @return Detailed string representation for debugging
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
	 * Returns compact string representation showing only keys.
	 * Format: [key1, key2, key3, ...] similar to ArrayList.
	 * @return Compact string representation of keys
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
	 * Returns an array containing all keys in the hash map.
	 * Order is not guaranteed.
	 * @return Array of all keys
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
	
	/**
	 * Returns an array of keys whose values meet the threshold.
	 * Only includes keys with values greater than or equal to thresh.
	 * @param thresh Minimum value threshold for inclusion
	 * @return Array of keys with values >= thresh
	 */
	public long[] toArray(long thresh){
		int len=0;
//		assert(verify());
		for(int i=0; i<values.length; i++){
			assert((values[i]==0)==(keys[i]==invalid)) : i+", "+values[i]+", "+keys[i]+", "+invalid+"\n"+toStringSetView();
			assert((keys[i]<0)==((keys[i]==invalid))) : toStringSetView();
			if(values[i]>=thresh){
				assert(keys[i]>=0) : "\nNegative key ("+keys[i]+", "+values[i]+", "+i+") for thresh "+thresh+":\n"+toStringSetView();
				len++;
			}
		}
		long[] x=KillSwitch.allocLong1D(len);
		for(int i=0, j=0; j<len; i++){
			if(values[i]>=thresh){
				x[j]=keys[i];
				assert(keys[i]>=0) : "\nNegative key ("+keys[i]+", "+values[i]+", "+i+") for thresh "+thresh+":\n"+toStringSetView();
				j++;
			}
		}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Validates internal consistency of the hash map structure.
	 * Checks that all keys hash to correct locations and counts match.
	 * Used for debugging and testing purposes.
	 * @return true if structure is valid, false otherwise
	 */
	public boolean verify(){
		if(keys==null){return true;}
		int numValues=0;
		int numFound=0;
		for(int i=0; i<keys.length; i++){
			final long key=keys[i];
			final long value=values[i];
			
			if(key==invalid){
				if(value!=0){
					assert(false) : i+", "+key+", "+value;
					return false;
				}
			}else{
				numValues++;
				if(value<1){
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
	 * Rehashes all entries after the initial position to maintain hash table integrity.
	 * Called after removing an entry to fill gaps and preserve lookup correctness.
	 * @param initial Starting position for rehashing
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
	
	private boolean rehashCell(final int cell){
		final long key=keys[cell];
		final long value=values[cell];
		assert(key!=invalid);
		if(key==invalid){resetInvalid();}
		final int dest=findCellOrEmpty(key);
		if(cell==dest){return false;}
		assert(keys[dest]==invalid);
		keys[cell]=invalid;
		values[cell]=0;
		keys[dest]=key;
		values[dest]=value;
		return true;
	}
	
	/** Generates a new invalid key marker when the current one collides.
	 * Updates all cells using the old invalid marker to the new one. */
	private void resetInvalid(){
		final long old=invalid;
		long x=invalid;
		while(x==old || contains(x)){x=randy.nextLong()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<keys.length; i++){
			if(keys[i]==old){
				assert(values[i]==0);
				keys[i]=invalid;
			}
		}
	}
	
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
	
	/** Resizes the hash table to approximately double the current capacity.
	 * Called when the size limit is reached. */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}
	
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
		final long[] oldValues=values;
		keys=KillSwitch.allocLong1D(size3);
		values=KillSwitch.allocLong1D(size3);
		Arrays.fill(keys, invalid);
		
//		System.err.println("Resizing "+(old==null ? "null" : ""+old.length)+" to "+size3);
		
		if(size<1){return;}
		
		size=0;
		for(int i=0; i<oldKeys.length; i++){
			long key=oldKeys[i];
			if(key!=invalid){
				put(key, oldValues[i]);
			}
		}
//		assert(verify()); //123
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Multivalue          ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns direct reference to the internal keys array */
	public long[] keys() {return keys;}

	/** Returns direct reference to the internal values array */
	public long[] values() {return values;}

	/** Returns the current invalid key marker used for empty cells */
	public long invalid() {return invalid;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Storage for keys with multiple values */
	private ArrayList<LongList> multivalues=new ArrayList<LongList>(4);
	
	/** Hash table array storing keys */
	private long[] keys;
	/** Hash table array storing values or list indices */
	private long[] values;
	/** Number of key-value mappings currently stored */
	private int size=0;
	/** Value for empty cells */
	private long invalid;
	/** Prime number used for hash table size and modular arithmetic */
	private int modulus;
	/** Maximum size before triggering resize based on load factor */
	private int sizeLimit;
	/** Load factor determining when to resize the hash table */
	private final float loadFactor;
	
	/** Random number generator for creating invalid key markers */
	private static final Random randy=new Random(1);
	/** Bit mask for positive long values used in hashing */
	private static final long MASK=Long.MAX_VALUE;
	/** Bit mask ensuring invalid markers are negative */
	private static final long MINMASK=Long.MIN_VALUE;
	/** Return value indicating a key was not found in the map */
	public static final long NOTPRESENT=-1;
	/** Offset added to list indices to create negative reference codes */
	private static final long OFFSET=2;
	
	/** Extra capacity added during resize to reduce collision probability */
	private static final int extra=10;
	
}
