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
 * Hybrid hash map for mapping long keys to single or multiple long values.
 * Acts as a LongLongHashMap for single values and LongLongListHashMap for multi-value keys.
 * Stores multi-value entries as negative codes (-(index+OFFSET)) pointing to separate lists.
 * Supports non-negative values only. A repeated inline value is ignored; once a key is promoted
 * to a side list, subsequent values are appended even when they duplicate existing list entries.
 * Instances are unsynchronized; mutations require exclusive access. Returned side lists
 * and raw arrays are borrowed live storage, not independent snapshots.
 * @author Brian Bushnell
 * @date November 13, 2024
 */
public final class LongLongHashMapHybrid{

	/**
	 * Legacy ad hoc stress benchmark against {@link HashMap}.
	 * @param args ignored
	 */
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
			System.err.println("LongLongHashMapHybrid:");
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

	/** Constructs a map with the default initial modulus request and load factor. */
	public LongLongHashMapHybrid(){
		this(256);
	}

	/**
	 * Constructs a map with the default load factor.
	 * @param initialSize requested minimum prime modulus before adding slack cells;
	 * asserted positive
	 */
	public LongLongHashMapHybrid(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Constructs a map with a caller-supplied initial modulus and load factor.
	 * @param initialSize requested minimum prime modulus before adding slack cells;
	 * asserted positive
	 * @param loadFactor_ requested resize threshold fraction; asserted in (0,1)
	 * before clamping to [0.25, 0.90]
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

	/** Removes all key/value pairs and clears any side-list storage. */
	public void clear(){
		//Clear side lists even when size==0: removing only non-tail multi-value keys can leave null holes.
		multivalues.clear();
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
	 * Returns the inline value or multi-value code mapped from key.
	 * @param key key to search for
	 * @return inline value, negative side-list code, or NOTPRESENT
	 */
	public long get(long key){
//		assert(verify()); //123
		//returns: single value (>=0), a negative multi-value code (<=-2, decode via getOrFill/fill), or
		//NOTPRESENT(-1) on miss. The three ranges are disjoint (codes<=-2, -1=absent, singles>=0) so -1 is
		//unambiguous despite codes also being negative.
		long value=NOTPRESENT;
		if(key!=invalid){
			int cell=findCell(key);
			if(cell>=0){value=values[cell];}
		}
		return value;
	}

	/**
	 * Decodes a valid live multi-value code into its borrowed side list.
	 * @param code negative code from get(long), at most {@code -OFFSET};
	 *        {@code NOTPRESENT} is not a side-list code
	 * @return borrowed live side list for code
	 */
	public LongList getListFromCode(long code){
		assert(code<=-OFFSET) : code;
		return multivalues.get((int)(-code-OFFSET));
	}

	/**
	 * Returns key's values, clearing buffer first. The buffer must not alias any
	 * owned side list.
	 * @param key key to search for
	 * @param buffer caller-owned scratch list
	 * @return buffer for absent or inline keys, or the live side list for multi-value keys
	 */
	public LongList getOrFill(long key, LongList buffer){
//		assert(buffer.isEmpty());
		buffer.clear();
		long code=get(key);
		if(code==NOTPRESENT){return buffer;}
		if(code>=0){buffer.add(code); return buffer;}
		return multivalues.get((int)(-code-OFFSET));
	}

	/**
	 * Clears buffer and copies key's values into it. The buffer must not alias
	 * any owned side list.
	 * @param key key to search for
	 * @param buffer caller-owned destination list
	 * @return inline value, negative side-list code, or NOTPRESENT
	 */
	public long fill(long key, LongList buffer){
//		assert(buffer.isEmpty());
		buffer.clear();
		long code=get(key);
		if(code==NOTPRESENT){return code;}
		if(code>=0){buffer.add(code);}
		else{buffer.addAll(multivalues.get((int)(-code-OFFSET)));}
		return code;
	}

	/**
	 * Maps the specified key to the specified value.
	 * Validates that value is non-negative before storing.
	 * @param key The key to map
	 * @param value The non-negative value to associate with key
	 * @return true if the value was added, false only if it matched the existing inline single value
	 */
	public boolean put(long key, long value){
		//public contract: single values must be >=0. Internally values[] encodes multi-value keys as NEGATIVE
		//codes -(listIndex+OFFSET) (<=-2), so the sign bit distinguishes single (>=0) from multi (<=-2); -1 is
		//reserved as NOTPRESENT. (resize must bypass this assert to copy codes -> it uses putInner; see #001.)
		assert(value>=0) : "Unsupported negative value "+value;
		return putInner(key, value);
	}

	/**
	 * Internal method that performs the actual key-value mapping.
	 * Handles single to multiple value transitions by creating lists.
	 * Automatically resizes when size limit is exceeded.
	 * @param key The key to map
	 * @param value The value to associate with key
	 * @return true if the value was added, false only if it matched the existing inline single value
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
		if(vCell==value){return false;}
		//CLEVER [verified]: single->multi transition. vCell>=0 means one inline value; promote it to a side
		//LongList and store the negative code in its place. Thereafter vCell<0 (a code) just appends to that list.
		//No size++ here -- multi-value adds add a VALUE, not a KEY.
		if(vCell>=0){
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
//		if(list.contains(value)){return false;}//Slow and not really needed for indexing.
		list.add(value);
//		assert(verify()); //123
		return true;
	}

	/**
	 * Removes the specified key and all its values from the hash map.
	 * For multi-value keys, nulls out the list reference but may leave holes.
	 * Performs rehashing to maintain hash table integrity.
	 * @param key The key to remove
	 * @return The removed value or list index code, NOTPRESENT if not found
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
		if(value<0){
			//multi-value removal: null the side list (leaves a hole unless it was the last, then shrink). Safe:
			//only the LAST listIndex is ever shrunk, and no other key holds a higher index, so no dangling code.
			//A non-last removal leaves a null hole (documented; fine for indexing -- the key is gone, unreachable).
			int idx=(int)(-value-OFFSET);
			multivalues.set(idx, null);
			if(idx+1==multivalues.size()){multivalues.remove(idx);}
		}

		rehashFrom(cell);
//		assert(verify());
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
	 * Returns a string representation showing only the keys.
	 * Delegates to toStringListView for compact output.
	 * @return String representation of keys
	 */
	@Override
	public String toString(){
		return toStringListView();
	}

	/**
	 * Returns a list of occupied cells with cell index, key and stored value/code.
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
	 * Returns mapped keys whose stored inline value/code is at or above thresh.
	 * Multi-value keys store negative side-list codes, so negative thresholds
	 * test the code itself, not the side-list contents.
	 * @param thresh minimum stored value/code to include
	 * @return qualifying keys
	 */
	public long[] toArray(long thresh){
		//FIXED [map/LongLongHashMapHybrid#002] - returns the OCCUPIED keys whose (single) value>=thresh. Previously
		//asserted the counting invariant (value==0 iff empty, key>=0), which false-failed on this Hybrid's value-0
		//singles, multi-value codes, and negative keys. Now asserts only the real invariant (empty cell => value 0)
		//and guards on keys[i]!=invalid so an empty cell is never emitted (the old code could emit the sentinel for
		//thresh<=0). Multi-value keys carry a negative code, so they are naturally excluded for any thresh>=0.
		int len=0;
//		assert(verify());
		for(int i=0; i<values.length; i++){
			assert(keys[i]!=invalid || values[i]==0) : i+", "+values[i]+", "+keys[i]+", "+invalid+"\n"+toStringSetView();
			if(keys[i]!=invalid && values[i]>=thresh){len++;}
		}
		long[] x=KillSwitch.allocLong1D(len);
		for(int i=0, j=0; j<len; i++){
			if(keys[i]!=invalid && values[i]>=thresh){
				x[j]=keys[i];
				j++;
			}
		}
		return x;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Checks Hybrid-specific table invariants.
	 * @return true if all occupied keys are findable and every stored side-list code resolves to a live list
	 */
	public boolean verify(){
		//FIXED [map/LongLongHashMapHybrid#002] - verify() previously carried the COUNTING-map invariant copied from
		//LongLongHashMap (present=>value>=1, value==0 iff empty, key>=0), which is FALSE for this Hybrid: single
		//values may be 0 (put allows >=0), multi-value keys store a NEGATIVE code (<=-OFFSET), and keys may be any
		//long (incl. negative). It now checks the Hybrid's REAL invariants: empty cell (key==invalid) => value==0;
		//occupied cell => value is a single (>=0) OR a valid multi-value code (decodes to a live multivalues[] entry);
		//and each key is findable at its own cell. The same counting-invariant bug was fixed in toArray(thresh) above.
		//Dev/test checker (assert(verify()) pattern; the //123 call sites are commented out) but correct when invoked.
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
				//occupied: single value>=0, OR a multi-value code (<=-OFFSET) decoding to a live multivalues[] entry.
				//-1=NOTPRESENT is never a stored value (its code would decode to idx=-1, caught below).
				if(value<0){
					final int idx=(int)(-value-OFFSET);
					if(idx<0 || idx>=multivalues.size() || multivalues.get(idx)==null){
						assert(false) : i+", "+key+", "+value+", idx="+idx+", lists="+multivalues.size();
						return false;
					}
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
				//FIXED [map/LongLongHashMapHybrid#001]: re-insert via putInner, NOT put. put() asserts value>=0,
				//but multi-value keys store a NEGATIVE code (-(listIndex+OFFSET)) in values[]; re-inserting that
				//code through put() threw AssertionError ("Unsupported negative value") on any resize of a map
				//holding >=1 multi-value key (-ea is always on). putInner copies the code verbatim (multivalues[]
				//is untouched by resize), so it still decodes to the same list. Output-neutral for single-value
				//maps. Empirically proven: crash before / correct decode after. (0 callers today -> latent until used.)
				putInner(key, oldValues[i]);
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

	/**
	 * Returns the borrowed backing key array; empty cells contain invalid().
	 * @return borrowed backing key array
	 */
	public long[] keys(){return keys;}

	/**
	 * Returns the borrowed backing value/code array parallel to keys().
	 * @return borrowed backing value/code array
	 */
	public long[] values(){return values;}

	/**
	 * Returns the current negative empty-cell sentinel.
	 * @return empty-cell sentinel
	 */
	public long invalid(){return invalid;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Side lists for keys that have multiple mapped values. */
	private ArrayList<LongList> multivalues=new ArrayList<LongList>(4);

	/** Keys array; empty cells hold invalid. */
	private long[] keys;
	/** Values or negative side-list codes parallel to keys; empty cells hold zero. */
	private long[] values;
	/** Number of occupied cells. */
	private int size=0;
	/** Sentinel value marking empty cells in the hash table */
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
	/** Return code for missing keys. */
	public static final long NOTPRESENT=-1;
	/** Offset between zero-based side-list indexes and negative stored codes. */
	private static final long OFFSET=2;

	/** Extra cells after the prime modulus to simplify wrapped probing. */
	private static final int extra=10;

}
