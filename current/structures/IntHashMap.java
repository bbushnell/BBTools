package structures;

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;

/**
 * @author Brian Bushnell
 * @date June 7, 2017
 *
 */
public final class IntHashMap extends AbstractIntHashMap implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 5525726007591609843L;

	public static void main(String[] args){
//		IntHashMap set=new IntHashMap(20, 0.7f);
//		test(set);
		int size=args.length>0 ? Integer.parseInt(args[0]) : 100000000;
		bench(size);
	}
	
	/**
	 * Benchmarks IntHashMap performance against HashMap<Integer,Integer>.
	 * Creates maps, populates them with sequential key-value pairs, and measures timing.
	 * @param size Number of key-value pairs to insert during benchmark
	 */
	private static void bench(int size){
		System.gc();
		Timer t=new Timer();
		
		{
			System.err.println("\nIntHashMap:");
			Shared.printMemory();
			t.start();
			IntHashMap map=new IntHashMap();
			for(int i=0; i<size; i++){
				map.put(i, 2*i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			map=null;
			System.gc();
		}
		
		{
			System.err.println("\nHashMap<Integer, Integer>:");
			Shared.printMemory();
			t.start();
			HashMap<Integer, Integer> map=new HashMap<Integer, Integer>();
			for(int i=0; i<size; i++){
				map.put(i, 2*i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			map=null;
			System.gc();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates an IntHashMap with default initial size of 256 */
	public IntHashMap(){
		this(256);
	}
	
	/** Creates an IntHashMap with specified initial size and default load factor.
	 * @param initialSize Initial capacity for the hash table */
	public IntHashMap(int initialSize){
		this(initialSize, 0.7f);
	}
	
	/**
	 * Creates an IntHashMap with specified initial size and load factor.
	 * Load factor is clamped to the range [0.25, 0.90] for performance.
	 * @param initialSize Initial capacity for the hash table
	 * @param loadFactor_ Target load factor (ratio of occupied to total slots)
	 */
	public IntHashMap(int initialSize, float loadFactor_){
		invalid=randy.nextInt()|MINMASK;
		assert(invalid<0);
		assert(initialSize>0);
		assert(loadFactor_>0 && loadFactor_<1);
		loadFactor=Tools.mid(0.25f, loadFactor_, 0.90f);
		resize(initialSize);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
		size=0;
	}
	
	@Override
	public int get(int key){
		int cell=findCell(key);
		return cell<0 ? -1 : values[cell];
	}
	
	@Override
	public int put(int key, int value){return set(key, value);}
	
	/** Copies all key-value mappings from the specified map to this map.
	 * @param map IntHashMap whose mappings are to be copied */
	public void putAll(IntHashMap map) {
		for(int i=0; i<map.keys.length; i++) {
			if(map.keys[i]!=map.invalid) {
				put(map.keys[i], map.values[i]);
			}
		}
	}
	
	@Override
	public int set(int key, int value){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final int oldV=values[cell];
		values[cell]=value;
		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
//		assert(get(key)==value);//TODO: slow
		return oldV;
	}
	
	@Override
	public int increment(int key){
		return increment(key, 1);
	}
	
	@Override
	public int increment(int key, int incr){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final int oldV=values[cell];
		final int value=oldV+incr;
		values[cell]=value;
		values[cell]=Tools.min(Integer.MAX_VALUE, value);
		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
//		assert(get(key)==value);//TODO: slow
		return value;
	}
	
	/**
	 * Increments all keys in this map by the corresponding values from another map.
	 * For each key-value pair in the source map, adds that value to the current value.
	 * @param map IntHashMap containing increment values for each key
	 */
	public void incrementAll(IntHashMap map) {
		for(int i=0; i<map.keys.length; i++) {
			if(map.keys[i]!=map.invalid) {
				increment(map.keys[i], map.values[i]);
			}
		}
	}
	
	/**
	 * Sets each key to the maximum of its current value and the value from another map.
	 * Only updates keys that exist in the source map.
	 * @param map IntHashMap containing values to compare against
	 */
	public void setToMax(IntHashMap map) {
		for(int i=0; i<map.keys.length; i++) {
			final int key=map.keys[i];
			if(key!=map.invalid) {
				put(key, Tools.max(map.values[i], get(i)));
			}
		}
	}
	
	@Override
	public boolean remove(int key){
		if(key==invalid){return false;}
		final int cell=findCell(key);
		if(cell<0){return false;}
		assert(keys[cell]==key);
		keys[cell]=invalid;
		values[cell]=0;
		size--;
		
		rehashFrom(cell);
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Rehashes all entries following a deletion to maintain clustering properties.
	 * Scans from the deletion point forward, then wraps around to maintain probe sequences.
	 * @param initial Cell position where rehashing should begin
	 */
	private void rehashFrom(int initial){
		if(size<1){return;}
		final int limit=keys.length;
		for(int cell=initial+1; cell<limit; cell++){
			final int key=keys[cell];
			if(key==invalid){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final int key=keys[cell];
			if(key==invalid){return;}
			rehashCell(cell);
		}
	}
	
	/**
	 * Attempts to move a key-value pair to its optimal position in the hash table.
	 * Used during deletion cleanup to maintain efficient probe sequences.
	 * @param cell Index of the cell to rehash
	 * @return true if the cell was moved to a different position, false otherwise
	 */
	private boolean rehashCell(final int cell){
		final int key=keys[cell];
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
	
	/**
	 * Generates a new invalid marker value when the current one conflicts with a key.
	 * Scans existing keys to ensure the new invalid value doesn't collide.
	 * Updates all empty cells to use the new invalid marker.
	 */
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;}
		assert(x<0);
		invalid=x;
		Vector.changeAll(keys, old, x);
	}

	@Override
	int findCell(final int key){
//		if(key==invalid){return -1;}
		final int initial=((key&MASK)%modulus);
		return Vector.findKeyScalar(keys, key, initial, invalid);
	}
	
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		final int initial=(int)((key&MASK)%modulus);
		return Vector.findKeyOrInvalid(keys, key, initial, invalid);
	}
	
	/**
	 * Doubles the hash table capacity when the load factor threshold is exceeded
	 */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}
	
	/**
	 * Resizes the hash table to accommodate more entries.
	 * Finds the next suitable prime number and rehashes all existing entries.
	 * @param size2 Target minimum capacity for the new hash table
	 */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		long newPrime=Primes.primeAtLeast(size2);
		if(newPrime+extra>Integer.MAX_VALUE){
			newPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra);
		}
		assert(newPrime>modulus) : "Overflow: "+size+", "+size2+", "+modulus+", "+newPrime;
		modulus=(int)newPrime;
		
		final int size3=(int)(newPrime+extra);
		final float lf=(size3<0x60000000 ? loadFactor : Math.max(loadFactor, 0.80f));
		sizeLimit=(int)(modulus*lf);
		final int[] oldK=keys;
		final int[] oldV=values;
		keys=KillSwitch.allocInt1D(size3);
		values=KillSwitch.allocInt1D(size3);
		Arrays.fill(keys, invalid);
		
//		System.err.println("Resizing "+(old==null ? "null" : ""+old.length)+" to "+size3);
		
		if(size<1){return;}
		
		size=0;
		for(int i=0; i<oldK.length; i++){
			final int k=oldK[i], v=oldV[i];
			if(k!=invalid){
				set(k, v);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int[] toArray(){
		int[] x=KillSwitch.allocInt1D(size);
		int i=0;
		for(int key : keys){
			if(key!=invalid){
				x[i]=key;
				i++;
			}
		}
		return x;
	}
	
	@Override
	public int[] keys(){return keys;}
	
	@Override
	public int[] values(){return values;}
	
	@Override
	public int invalid(){return invalid;}
	
	@Override
	public int size(){return size;}
	
	@Override
	public boolean isEmpty(){return size==0;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Array storing hash table keys; empty slots contain the invalid marker */
	private int[] keys;
	/** Array storing hash table values corresponding to keys array */
	private int[] values;
	/** Number of key-value pairs currently stored in the map */
	private int size=0;
	/** Value for empty cells */
	private int invalid;
	/** Prime number used for hash function modulo operation */
	private int modulus;
	/** Maximum entries before triggering a resize (capacity * loadFactor) */
	private int sizeLimit;
	/** Target ratio of occupied to total slots for resize triggering */
	private final float loadFactor;
	
	/** Random number generator for creating invalid marker values */
	private static final Random randy=new Random(1);
	
}
