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
 * Primitive long-to-int hash map backed by parallel key and value arrays.
 * Uses open addressing with linear probing and a negative sentinel key to
 * mark empty cells. Not synchronized; callers need external coordination
 * for concurrent mutation. put(long,int) inserts only absent keys; existing
 * mappings are unchanged.
 *
 * @author Brian Bushnell
 * @date August 3, 2017
 */
public final class LongHashMap{

	/**
	 * Legacy ad hoc stress benchmark against {@link HashMap}.
	 * @param args ignored
	 */
	public static void main(String[] args){
		Random randy2=Shared.threadLocalRandom();
		LongHashMap map=new LongHashMap(20, 0.7f);
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
			System.err.println("LongHashMap:");
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

	/** Constructs a map with the default requested capacity and load factor. */
	public LongHashMap(){
		this(256);
	}

	/**
	 * Constructs a map with the default load factor.
	 * @param initialSize requested minimum prime modulus before adding slack cells;
	 * asserted positive
	 */
	public LongHashMap(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Constructs a map with a caller-supplied initial modulus and load factor.
	 * @param initialSize requested minimum prime modulus before adding slack cells;
	 * asserted positive
	 * @param loadFactor_ requested resize threshold fraction; asserted in (0,1)
	 * before clamping to [0.25, 0.90]
	 */
	public LongHashMap(int initialSize, float loadFactor_){
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

	/** Removes all mappings and resets every cell to empty. */
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
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
	 * Returns the value mapped from key.
	 * @param key key to search for
	 * @return stored value, or -1 when key is absent
	 */
	public int get(long key){
//		assert(verify()); //123
		int value=-1;
		if(key!=invalid){
			int cell=findCell(key);
			if(cell>=0){value=values[cell];}
		}
		return value;
	}

	/**
	 * Increments a key's value by 1.
	 * @param key key to increment
	 * @return new value after wrapping int addition, or 1 for a newly inserted key
	 */
	public int add(long key){
		return increment(key, 1);
	}

	/**
	 * Increments a key's value by incr.
	 * @param key key to increment
	 * @param incr amount to add; also the inserted value for a new key
	 * @return new value after wrapping int addition, or incr for a newly inserted key
	 */
	public int increment(long key, int incr){
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
	 * Maps key to value when key is absent.
	 * @param key key to map
	 * @param value value to store for a new key
	 * @return true if key was newly added, false if key was already present
	 */
	public boolean put(long key, int value){
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
	 * Removes a key and rehashes following cells to preserve probe chains.
	 * @param key key to remove
	 * @return previous value, or -1 when key was absent
	 */
	public int remove(long key){
//		assert(verify()); //123
		if(key==invalid){return -1;}
		final int cell=findCell(key);
		if(cell<0){return -1;}
		assert(keys[cell]==key);
		keys[cell]=invalid;
		final int value=values[cell];
		values[cell]=0;
		size--;

		rehashFrom(cell);
//		assert(verify()); //123
		return value;
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
	 * Returns string representation showing only the keys as a list.
	 * @return list-view string
	 */
	@Override
	public String toString(){
		return toStringListView();
	}

	/**
	 * Returns a list of occupied cells with cell index, key and value.
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

	/**
	 * Counts keys with values at or above thresh and returns their minimum key.
	 * Asserts that all mapped keys are non-negative and all occupied values are
	 * nonzero; thresh must be positive to exclude empty cells.
	 * @param thresh positive minimum value to count
	 * @return {@code [minimum key, count]}, or {@code [invalid(), 0]} if no
	 * values qualify
	 */
	public long[] getMin(int thresh){
		int found=0;
		long min=Long.MAX_VALUE;
		for(int i=0; i<values.length; i++){
			assert((values[i]==0)==(keys[i]==invalid)) : i+", "+values[i]+", "+keys[i]+", "+invalid+"\n"+toStringSetView();
			assert((keys[i]<0)==((keys[i]==invalid))) : toStringSetView();
			if(values[i]>=thresh){
				assert(keys[i]>=0) : "\nNegative key ("+keys[i]+", "+values[i]+", "+i+") for thresh "+thresh+":\n"+toStringSetView();
				found++;
				min=Tools.min(min, keys[i]);
			}
		}
		return new long[] {found>0 ? min : invalid, found};
	}

	/**
	 * Returns keys with values at or above thresh in backing-array order.
	 * Asserts that all mapped keys are non-negative and all occupied values are
	 * nonzero.
	 * @param thresh minimum value to include, clamped to at least 1
	 * @return qualifying keys
	 */
	public long[] toArray(int thresh){
		thresh=Tools.max(thresh, 1);
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
	 * Development invariant check for positive-counting map uses.
	 * @return true when size, empty-cell and positive-value invariants hold
	 */
	public boolean verify(){
		//repOK-style dev checker (all callers are commented `assert(verify()); //123` -> off in production, must work
		//when called). Uses findCell(key) -> CORRECT (index-vs-key family slip N/A; matches the Long* audit). Encodes
		//the positive-counting invariant: present key => value>=1; empty cell => value 0 (also asserted with keys>=0 in
		//getMin/toArray(thresh)). By-design for positive-count callers; put()/increment() don't enforce it.
		if(keys==null){return true;}
		int numValues=0;
		int numFound=0;
		for(int i=0; i<keys.length; i++){
			final long key=keys[i];
			final int value=values[i];

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

	/**
	 * Moves the key in cell to its earliest available probe-chain slot.
	 * @param cell occupied cell to test
	 * @return true if the key moved
	 */
	private boolean rehashCell(final int cell){
		final long key=keys[cell];
		final int value=values[cell];
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

	/** Chooses a new negative empty-cell sentinel and rewrites old empty cells. */
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

	/** Expands the table and rehashes mapped keys into the new arrays. */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}

	/**
	 * Rehashes all mapped keys into a prime-sized table.
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
		final int[] oldValues=values;
		keys=KillSwitch.allocLong1D(size3);
		values=KillSwitch.allocInt1D(size3);
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

	/**
	 * Returns the borrowed backing key array; empty cells contain invalid().
	 * Mutating the borrowed array can break invariants, and resize() can
	 * replace it while clear() mutates it in place.
	 * @return borrowed backing key array
	 */
	public long[] keys(){return keys;}

	/**
	 * Returns the borrowed backing value array parallel to keys().
	 * Mutating the borrowed array can break invariants, and resize() can
	 * replace it while clear() mutates it in place.
	 * @return borrowed backing value array
	 */
	public int[] values(){return values;}

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
	/** Values parallel to keys; empty cells hold zero. */
	private int[] values;
	/** Number of occupied cells. */
	private int size=0;
	/** Randomly generated marker value representing empty cells. */
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
