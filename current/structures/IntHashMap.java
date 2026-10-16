package structures;

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import shared.Random;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;

/**
 * Hash map implementation with integer keys and integer values.
 * Uses open addressing with linear probing for collision resolution.
 * Provides better memory efficiency than {@code HashMap<Integer,Integer>} for integer mappings.
 * Absent lookups return -1; empty cells are marked by a private negative
 * sentinel that is regenerated if a caller inserts that sentinel value as a key.
 *
 * @author Brian Bushnell
 * @date June 7, 2017
 */
public final class IntHashMap extends AbstractIntHashMap implements Serializable{
	
	private static final long serialVersionUID=5525726007591609843L;

	public static void main(String[] args){
//		IntHashMap set=new IntHashMap(20, 0.7f);
//		test(set);
		int size=args.length>0 ? Integer.parseInt(args[0]) : 100000000;
		bench(size);
	}
	
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
	
	/**
	 * Constructs a map with the default initial capacity.
	 */
	public IntHashMap(){
		this(256);
	}
	
	/**
	 * Constructs a map with the requested initial capacity and the default load
	 * factor.
	 * @param initialSize Initial slot target; must be positive.
	 */
	public IntHashMap(int initialSize){
		this(initialSize, 0.7f);
	}
	
	/**
	 * Constructs a map with the requested initial capacity and load factor.
	 * The load factor is clamped to the supported resize range.
	 *
	 * @param initialSize Initial slot target; must be positive.
	 * @param loadFactor_ Requested load factor; must be strictly between 0 and 1.
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
	
	/**
	 * Removes all keys and values while retaining the current backing arrays.
	 */
	@Override
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
		size=0;
	}
	
	/**
	 * Returns the value mapped to a key.
	 * @param key Key to search for.
	 * @return Stored value, or -1 if the key is absent or equals the current
	 * invalid sentinel.
	 */
	@Override
	public int get(int key){
		int cell=findCell(key);
		return cell<0 ? -1 : values[cell];
	}
	
	@Override
	public int put(int key, int value){return set(key, value);}
	
	/**
	 * Copies every live mapping from another map into this map, overwriting
	 * existing values for matching keys.
	 * @param map Source map.
	 */
	public void putAll(IntHashMap map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				put(map.keys[i], map.values[i]);
			}
		}
	}
	
	/**
	 * Associates a value with a key, overwriting any previous value.
	 *
	 * @param key Key to map.
	 * @param value Value to store.
	 * @return Previous value in the selected cell, or 0 for a newly inserted key.
	 */
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
	
	/**
	 * Increments a key's value by one.
	 * @param key Key to increment.
	 * @return Updated value.
	 */
	@Override
	public int increment(int key){
		return increment(key, 1);
	}
	
	/**
	 * Increments a key's value with upper saturation at {@link Integer#MAX_VALUE}.
	 * Absent keys start from zero and are inserted with the incremented value.
	 *
	 * @param key Key to increment.
	 * @param incr Increment amount.
	 * @return Updated value after saturation.
	 */
	@Override
	public int increment(int key, int incr){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final long oldV=values[cell];
		final int value=(int)Tools.min(Integer.MAX_VALUE, oldV+incr);
		values[cell]=value;
		if(keys[cell]==invalid){
			keys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
//		assert(get(key)==value);//TODO: slow
		return value;
	}
	
	/**
	 * Adds every live count from another map into this one.
	 * @param map Source map.
	 */
	public void incrementAll(IntHashMap map){
		for(int i=0; i<map.keys.length; i++){
			if(map.keys[i]!=map.invalid){
				increment(map.keys[i], map.values[i]);
			}
		}
	}
	
	/**
	 * Replaces each overlapping value with the larger of the two map values.
	 * For absent keys, inserts the larger of the source value and the -1
	 * missing-value sentinel returned by {@link #get(int)}.
	 * @param map Source map.
	 */
	public void setToMax(IntHashMap map){
		for(int i=0; i<map.keys.length; i++){
			final int key=map.keys[i];
			if(key!=map.invalid){
				put(key, Tools.max(map.values[i], get(key)));
			}
		}
	}
	
	/**
	 * Removes one key from the map.
	 * @param key Key to remove.
	 * @return true if the key was present.
	 */
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
	
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;}
		assert(x<0);
		invalid=x;
		Vector.changeAll(keys, old, x);
	}

	/**
	 * Locates the cell containing the specified key using linear probing.
	 * Searches from the hash position forward, then wraps to the beginning if needed.
	 * @param key The key to search for
	 * @return Cell index if found, -1 if not present
	 */
	@Override
	int findCell(final int key){
		if(key==invalid){return -1;}
		final int initial=((key&MASK)%modulus);
		return Vector.findKeyScalar(keys, key, initial, invalid);
	}
	
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		final int initial=(int)((key&MASK)%modulus);
		return Vector.findKeyOrInvalid(keys, key, initial, invalid);
	}
	
	private final void resize(){
		assert(size>=sizeLimit);
		resize(keys.length*2L+1);
	}
	
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
	
	private int[] keys;
	private int[] values;
	private int size=0;
	/** Sentinel value marking empty cells in the keys array */
	private int invalid;
	private int modulus;
	private int sizeLimit;
	private final float loadFactor;
	
	private static final Random randy=shared.Shared.random(1);
	
}
