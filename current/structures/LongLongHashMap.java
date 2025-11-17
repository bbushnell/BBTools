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
 * @date December 12, 2017
 *
 */
public final class LongLongHashMap{
	
	/**
	 * Test harness and benchmark for LongLongHashMap performance.
	 * Compares against Java's HashMap for correctness and timing.
	 * @param args Command-line arguments (unused)
	 */
	public static void main(String[] args){
		Random randy2=Shared.threadLocalRandom();
		LongLongHashMap map=new LongLongHashMap(20, 0.7f);
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
			assert(map.get(x)==((int)(2*x)));
			assert(map2.get(x)==((int)(2*x)));
		}
		
		for(long x : list){
			assert(map.get(x)==((int)(2*x)));
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
			assert(map.get(x)==((int)(2*x)));
			assert(map2.get(x)==((int)(2*x)));
		}
		
		for(long x : list2){
			assert(map.get(x)==((int)(2*x)));
			assert(map2.get(x)==((int)(2*x)));
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
			System.err.println("LongLongHashMap:");
			t.start();
			for(int i=0; i<runs; i++){
//				for(long x : ll.array){
//					map.add(x);
//				}
				final long[] y=ll.array;
				for(int z=0; z<count; z++){
					final long key=y[z];
					map.add(key);
					map.contains(key);
					map.remove(key);
					map.add(key);
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
	
	/** Creates a LongLongHashMap with default initial capacity of 256. */
	public LongLongHashMap(){
		this(256);
	}
	
	/** Creates a LongLongHashMap with specified initial capacity and default load factor.
	 * @param initialSize Initial capacity for the hash table */
	public LongLongHashMap(int initialSize){
		this(initialSize, 0.7f);
	}
	
	/**
	 * Creates a LongLongHashMap with specified capacity and load factor.
	 * Load factor is clamped between 0.25 and 0.90 for performance.
	 * @param initialSize Initial capacity for the hash table
	 * @param loadFactor_ Target load factor before resizing (0.25-0.90)
	 */
	public LongLongHashMap(int initialSize, float loadFactor_){
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
	
	/** Removes all key-value mappings from the map.
	 * Resets all keys to the invalid sentinel value and values to zero. */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
		size=0;
//		assert(verify()); //123
	}
	
	/**
	 * Returns true if this map contains the specified key.
	 * @param key Key whose presence is to be tested
	 * @return true if this map contains the key, false otherwise
	 */
	public boolean contains(long key){
//		assert(verify()); //123
		return key==invalid ? false : findCell(key)>=0;
	}
	
	/**
	 * Returns true if this map contains the specified key.
	 * Alias for contains() method for consistency with Java Collections.
	 * @param key Key whose presence is to be tested
	 * @return true if this map contains the key, false otherwise
	 */
	public boolean containsKey(long key){
		return contains(key);
	}
	
	/**
	 * Returns the value to which the specified key is mapped.
	 * Returns -1 if the key is not present in the map.
	 * @param key The key whose associated value is to be returned
	 * @return The value to which the key is mapped, or -1 if not present
	 */
	public long get(long key){
//		assert(verify()); //123
		long value=-1;
		if(key!=invalid){
			int cell=findCell(key);
			if(cell>=0){value=values[cell];}
		}
		return value;//TODO: Should return invalid...
	}
	
	/**
	 * Increment this key's value by 1.
	 * @param key
	 * @return New value
	 */
	public long add(long key){
		return increment(key, 1);
	}
	
	/**
	 * Increment this key's value by incr.
	 * @param key
	 * @param incr
	 * @return New value
	 */
	public long increment(long key, long incr){
//		assert(verify()); //123
		if(key==invalid){resetInvalid();}
		int cell=findCellOrEmpty(key);
		if(keys[cell]==invalid){
			keys[cell]=key;
			values[cell]=incr;
			size++;
//			assert(verify()); //123
			if(size>sizeLimit){resize();}
//			assert(verify()); //123
			return incr;
		}else{
			values[cell]+=incr;
//			assert(verify()); //123
			return values[cell];
		}
	}
	
	/**
	 * Map this key to value.
	 * @param key
	 * @param value
	 * @return true if the key was added, false if it was already contained.
	 */
	public boolean put(long key, long value){
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
//		assert(verify()); //123
		return false;
	}
	
	/**
	 * Remove this key from the map.
	 * @param key
	 * @return Old value.
	 */
	public long remove(long key){
//		assert(verify()); //123
		if(key==invalid){return -1;}
		final int cell=findCell(key);
		if(cell<0){return -1;}
		assert(keys[cell]==key);
		keys[cell]=invalid;
		final long value=values[cell];
		values[cell]=0;
		size--;
		
		rehashFrom(cell);
//		assert(verify()); //123
		return value;
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
	 * Returns a detailed string representation showing array indices, keys, and values.
	 * Useful for debugging hash table internal structure.
	 * @return Detailed string showing (index, key, value) tuples
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
	
	/** Returns a compact string representation showing only the keys.
	 * @return String representation as comma-separated list of keys */
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
	 * The order is not guaranteed to be consistent.
	 * @return Array of all keys currently in the map
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
	 * Returns an array containing keys whose values meet or exceed the threshold.
	 * Includes validation to ensure returned keys are non-negative.
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
	 * Validates the internal consistency of the hash table.
	 * Checks that all occupied cells have valid keys and positive values,
	 * empty cells have zero values, and all keys can be found at their expected positions.
	 * @return true if the hash table is internally consistent, false otherwise
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
	 * Rehashes all entries starting from the specified position.
	 * Used after removal to maintain hash table clustering properties.
	 * Processes entries in circular order to avoid gaps.
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
	
	/**
	 * Generates a new sentinel value for marking empty cells.
	 * Called when a collision occurs with the current invalid value.
	 * Updates all empty cells to use the new sentinel value.
	 */
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
	
	/** Doubles the hash table size and rehashes all entries.
	 * Called automatically when the size limit is exceeded. */
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
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the internal key array (includes empty cells marked with invalid) */
	public long[] keys() {return keys;}

	/** Returns the internal value array (includes zeros for empty cells) */
	public long[] values() {return values;}

	/** Returns the current sentinel value used to mark empty cells */
	public long invalid() {return invalid;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Array storing hash table keys, with empty cells marked by invalid value */
	private long[] keys;
	/** Array storing hash table values, with empty cells containing zero */
	private long[] values;
	/** Number of occupied key-value pairs currently in the map */
	private int size=0;
	/** Value for empty cells */
	private long invalid;
	/** Prime number used as modulus for hash function calculations */
	private int modulus;
	/** Maximum entries before triggering automatic resize operation */
	private int sizeLimit;
	/** Target ratio of occupied cells to total capacity */
	private final float loadFactor;
	
	/** Random number generator for creating sentinel invalid values */
	private static final Random randy=new Random(1);
	/** Bit mask (Long.MAX_VALUE) used in hash calculations */
	private static final long MASK=Long.MAX_VALUE;
	/** Bit mask (Long.MIN_VALUE) ensuring generated invalid values are negative */
	private static final long MINMASK=Long.MIN_VALUE;
	
	/** Extra capacity added beyond the prime modulus for hash table size */
	private static final int extra=10;
	
}
