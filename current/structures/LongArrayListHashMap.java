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
 * @author Brian Bushnell
 * @date September 15, 2023
 *
 */
public class LongArrayListHashMap <X> {
	
	/**
	 * Test harness comparing performance against standard HashMap.
	 * Runs benchmarks with sequential and random keys to validate correctness and speed.
	 * @param args Command line arguments (unused)
	 */
	public static void main(String[] args){
		Random randy2=Shared.threadLocalRandom();
		LongArrayListHashMap<Long> map=new LongArrayListHashMap<Long>(20, 0.7f);
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
			map.put(x, (2*x));
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
			map.put(x, (2*x));
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
		ArrayList<Long> ll=new ArrayList<Long>(count);
		for(int i=0; i<count; i++){ll.add(randy2.nextLong());}

		Shared.printMemory();
		Timer t=new Timer();
		for(int k=0; k<2; k++){
			System.err.println("LongHashMap:");
			t.start();
			for(int i=0; i<runs; i++){
//				for(long x : ll.array){
//					map.add(x);
//				}
				for(int z=0; z<count; z++){
					final Long key=ll.get(z);
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
				for(Long x : ll){
					map2.put(x, x.intValue());
					map2.containsKey(x);
					map2.remove(x);
					map2.put(x, x.intValue());
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
	
	/** Creates a new LongArrayListHashMap with default initial size of 256 */
	public LongArrayListHashMap(){
		this(256);
	}
	
	/**
	 * Creates a new LongArrayListHashMap with specified initial size.
	 * Uses default load factor of 0.7.
	 * @param initialSize Initial capacity of the hash map
	 */
	public LongArrayListHashMap(int initialSize){
		this(initialSize, 0.7f);
	}
	
	/**
	 * Creates a new LongArrayListHashMap with specified size and load factor.
	 * Load factor is clamped between 0.25 and 0.90 for stability.
	 * @param initialSize Initial capacity of the hash map
	 * @param loadFactor_ Target load factor for resize operations
	 */
	public LongArrayListHashMap(int initialSize, float loadFactor_){
		synchronized(this) {
			invalid=randy.nextLong()|MINMASK;
			assert(invalid<0);
			assert(initialSize>0);
			assert(loadFactor_>0 && loadFactor_<1);
			loadFactor=Tools.mid(0.25f, loadFactor_, 0.90f);
			resize(initialSize);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Removes all key-value mappings from the map.
	 * Recycles ArrayList instances longer than 10 elements to reduce memory usage.
	 * Resets size to zero and marks all cells as invalid.
	 */
	public void clear(){
		if(size<1){return;}
		for(int i=0; i<values.length /*&& size>0*/; i++) {//TODO: Enable the early exit
			if(keys[i]!=invalid) {
				keys[i]=invalid;
				int cap=values[i].size();
				values[i].clear();
				if(cap>10) {values[i]=null;}//Recycle long lists
				size--;
			}
			assert(values[i]==null || values[i].isEmpty());
		}
		assert(size==0);
		size=0;
//		assert(verify()); //123
	}
	
	/**
	 * Tests whether the specified key is present in the map.
	 * @param key The key to test for presence
	 * @return true if the key is present, false otherwise
	 */
	public boolean contains(long key){
//		assert(verify()); //123
		return key==invalid ? false : findCell(key)>=0;
	}
	
	/**
	 * Tests whether the specified key is present in the map.
	 * Alias for contains(long key) for HashMap compatibility.
	 * @param key The key to test for presence
	 * @return true if the key is present, false otherwise
	 */
	public boolean containsKey(long key){
		return contains(key);
	}
	
	/**
	 * Returns the ArrayList of values associated with the specified key.
	 * @param key The key whose associated values are to be returned
	 * @return ArrayList of values for the key, or null if key is not present
	 */
	public ArrayList<X> get(long key){
//		assert(verify()); //123
		final int cell=findCell(key);
		return cell<0 ? null : values[cell];
	}
	
	/**
	 * Map this key to value.
	 * @param key
	 * @param value
	 * @return true if the key was added, false if it was already contained.
	 */
	public boolean put(long key, X value){
//		assert(verify()); //123
		if(key==invalid){resetInvalid();}
		int cell=findCellOrEmpty(key);
		if(keys[cell]==invalid){
			keys[cell]=key;
			if(values[cell]==null) {
				values[cell]=new ArrayList<X>(2);
			}
			values[cell].add(value);
			size++;
			if(size>sizeLimit){resize();}
			return true;
		}
		values[cell].add(value);
		return false;
	}
	
	/**
	 * Map this key to value.
	 * @param key
	 * @param value
	 * @return true if the key was added, false if it was already contained.
	 */
	public boolean put(long key, ArrayList<X> value){
//		assert(verify()); //123
		if(key==invalid){resetInvalid();}
		int cell=findCellOrEmpty(key);
		if(keys[cell]==invalid){
			keys[cell]=key;
			if(values[cell]==null) {
				values[cell]=new ArrayList<X>(2);
			}
			values[cell].addAll(value);
			size++;
			if(size>sizeLimit){resize();}
			return true;
		}
		values[cell].addAll(value);
		return false;
	}
	
	/**
	 * Remove this key from the map.
	 * @param key
	 * @return true if it was present
	 */
	public boolean remove(long key){
//		assert(verify()); //123
		if(key==invalid){return false;}
		final int cell=findCell(key);
		if(cell<0){return false;}
		assert(keys[cell]==key);
		keys[cell]=invalid;
		int cap=values[cell].size();
		values[cell].clear();
		if(cap>=10) {values[cell]=null;}//Recycle long lists
		size--;
		
		rehashFrom(cell);
//		assert(verify()); //123
		return true;
	}
	
	/** Returns the number of key-value mappings in this map */
	public int size(){return size;}
	
	/** Returns true if this map contains no key-value mappings */
	public boolean isEmpty(){return size==0;}
	
	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	/**
	 * Returns a detailed string representation showing cell indices, keys, and values.
	 * Format: [(index, key, ArrayList), ...]
	 * Useful for debugging hash table structure.
	 * @return Detailed string representation of internal structure
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
	 * Returns a simple string representation showing only the keys.
	 * Format: [key1, key2, key3, ...]
	 * @return String representation showing only keys
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
	 * Returns an array containing all keys in this map.
	 * The array length equals the number of mappings.
	 * @return Array of all keys present in the map
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
	
	/**
	 * Validates the internal consistency of the hash table structure.
	 * Checks that all keys can be found at their expected positions
	 * and that value counts match the size field.
	 * @return true if the hash table is consistent, false otherwise
	 */
	public boolean verify(){
		if(keys==null){return true;}
		int numValues=0;
		int numFound=0;
		for(int i=0; i<keys.length; i++){
			final long key=keys[i];
			final ArrayList<X> value=values[i];
			
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
	 * Rehashes all entries starting from the specified position.
	 * Used after removal to maintain hash table integrity by moving
	 * displaced entries back to their proper positions.
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
	
	private boolean rehashCell(final int sourceCell){
		final long key=keys[sourceCell];
		final ArrayList<X> value=values[sourceCell];
		assert(key!=invalid);
		if(key==invalid){resetInvalid();}
		final int destCell=findCellOrEmpty(key);
		if(sourceCell==destCell){return false;}
		assert(keys[destCell]==invalid);
		keys[sourceCell]=invalid;
		values[sourceCell]=values[destCell];
		keys[destCell]=key;
		values[destCell]=value;
		return true;
	}
	
	/**
	 * Generates a new invalid marker value when the current one conflicts with data.
	 * Updates all empty cells to use the new invalid marker.
	 * Ensures the invalid marker is always negative and unique.
	 */
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
	
	/**
	 * Doubles the hash table size plus one when load factor is exceeded.
	 * Synchronized to prevent concurrent modifications during resize.
	 * Calls resize with keys.length*2L+1 as argument.
	 */
	private synchronized final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}
	
	/**
	 * Resizes the hash table to accommodate the specified size.
	 * Finds the next prime number for the new modulus to improve hash distribution.
	 * Rehashes all existing entries into the new table.
	 * @param size2 Target size for the new hash table
	 */
	@SuppressWarnings("unchecked")
	private synchronized final void resize(final long size2){
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
		final ArrayList<X>[] oldValues=values;
		keys=KillSwitch.allocLong1D(size3);
		values=KillSwitch.allocObject1D(size3, ArrayList.class);
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
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the internal keys array (includes invalid markers) */
	public long[] keys() {return keys;}

	/** Returns the internal values array (includes null entries) */
	public ArrayList<X>[] values() {return values;}

	/** Returns the current invalid marker value used for empty cells */
	public long invalid() {return invalid;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Array storing hash table keys, with invalid markers for empty cells */
	private long[] keys;
	/** Nested array of ArrayLists storing values corresponding to keys */
	private ArrayList<X>[] values;
//	private ArrayList<X> keyList;//TODO: Store keys here for quick clearing, perhaps
	/** Current number of key-value mappings in the map */
	private int size=0;
	/** Value for empty cells */
	private long invalid;
	/** Prime number used for hash calculation and array size */
	private int modulus;
	/** Maximum number of entries before triggering resize operation */
	private int sizeLimit;
	/** Target load factor for determining when to resize the hash table */
	private final float loadFactor;
	
	/** Random number generator for creating invalid marker values */
	private static final Random randy=new Random(1);
	/** Bit mask for ensuring positive hash values (Long.MAX_VALUE) */
	private static final long MASK=Long.MAX_VALUE;
	/** Bit mask for ensuring invalid markers are negative (Long.MIN_VALUE) */
	private static final long MINMASK=Long.MIN_VALUE;
	
	/** Extra slots added to hash table size for improved collision handling */
	private static final int extra=10;
	
}
