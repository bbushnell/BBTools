package structures;

import java.util.Arrays;
import shared.Random;

import shared.KillSwitch;
import shared.Tools;

/**
 * Int-to-int hash map with power-of-two primary table sizes.
 * Uses open addressing with linear probing and hashes with {@code key&modulus}.
 * Empty cells are marked by a private negative sentinel that is regenerated if
 * a caller inserts that sentinel value as a key.
 *
 * @author Brian Bushnell
 * @date June 8, 2017
 */
public final class IntHashMapBinary extends AbstractIntHashMap{

	public static void main(String[] args){
		IntHashMapBinary set=new IntHashMapBinary(32, 0.7f);
		test(set);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs a map with the default initial capacity. */
	public IntHashMapBinary(){
		this(256);
	}

	/**
	 * Constructs a map with the requested initial capacity and the default load
	 * factor.
	 * @param initialSize Initial primary slot target; zero asserts, while non-powers
	 * round up to the next power of two, capped at {@code 1<<30} and then
	 * clamped to at least 2.
	 */
	public IntHashMapBinary(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Constructs a map with the requested initial capacity and load factor.
	 * @param initialSize Initial primary slot target; zero asserts, while non-powers
	 * round up to the next power of two, capped at {@code 1<<30} and then
	 * clamped to at least 2.
	 * @param loadFactor_ Requested load factor; must be strictly between 0 and 1,
	 * then clamped to {@code [0.25,0.90]} before use.
	 */
	public IntHashMapBinary(int initialSize, float loadFactor_){
		if(Integer.bitCount(initialSize)>1){
			int zeros=Integer.numberOfLeadingZeros(initialSize);
			if(zeros<2){zeros=2;}
			initialSize=1<<(32-zeros);
		}
		invalid=randy.nextInt()|MINMASK;
		assert(invalid<0);
		assert(initialSize>0);
		assert(loadFactor_>0 && loadFactor_<1);
		loadFactor=Tools.mid(0.25f, loadFactor_, 0.90f);
		//[structures/IntHashMapBinary#001 FIXED 2026-06-21] clamp to >=2 (matches resize()'s own Tools.max(2,modulus+1) growth guard):
		//un-clamped resize(1) -> newModulus=size2-1=0 -> assert(newModulus>0) crash under -ea (always). Reachable via
		//AbstractIntHashMap.toCountHistogram() -> new IntHashMapBinary(Tools.mid(1,size(),64))==1 on a <=1-entry map (testformat.sh, tiny/
		//empty ZMW). initialSize==0 still asserts loud above. Twin of the empirically-proven map/IntHashMapBinary#001 (Eru); structures/
		//was V1-skim so it slipped. Cross-package lead, verified at source + applied by Furina.
		resize(Tools.max(2, initialSize));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Removes all mappings while retaining the current backing arrays. */
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
	 * Removes one key from the map, then rehashes the following probe run.
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

	/** Rehashes cells after a removed slot until the next empty cell. */
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

	/** Moves one live cell to its first valid probe position when needed. */
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
	 * Chooses a fresh empty-cell sentinel and rewrites all old empty cells.
	 * Called before inserting the live sentinel as a user key.
	 */
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<keys.length; i++){
			if(keys[i]==old){
				keys[i]=invalid;
//				assert(volues[i]==0); //TODO: slow
			}
		}
	}

	/**
	 * Locates a key by walking its wrapped linear-probe chain.
	 * @param key Key to search for.
	 * @return Cell index, or -1 when absent or equal to the current invalid sentinel.
	 */
	@Override
	int findCell(final int key){
		if(key==invalid){return -1;}

		final int limit=keys.length, initial=key&modulus;
		for(int cell=initial; cell<limit; cell++){
			final int x=keys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		for(int cell=0; cell<initial; cell++){
			final int x=keys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		return -1;
	}

	/**
	 * Locates a key or the first empty cell in its wrapped linear-probe chain.
	 * @param key Key to search for; must not be the current invalid sentinel.
	 * @return Cell index containing the key or available for insertion.
	 */
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";

		final int limit=keys.length, initial=key&modulus;
		for(int cell=initial; cell<limit; cell++){
			final int x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final int x=keys[cell];
			if(x==key || x==invalid){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}

	/** Doubles the primary table range after the load limit is exceeded. */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(Tools.max(2, modulus+1)*2L);
	}

	/**
	 * Replaces storage with a power-of-two primary range plus overflow cells.
	 * @param size2 Requested primary cell count; must be a power of two.
	 */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		assert(Long.bitCount(size2)==1) : size+", "+size2;
		long newModulus=size2-1;
		assert(newModulus>0) : newModulus;
		assert(newModulus+extra<Integer.MAX_VALUE) : "Overflow";
		assert(newModulus>modulus);

//		System.err.println("size="+size+", modulus="+modulus+" -> size2="+size2+", modulus2="+newModulus);

		modulus=(int)newModulus;

		final int size3=(int)(newModulus+extra);
		sizeLimit=(int)(newModulus*loadFactor);
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

	/**
	 * Copies all live keys into a dense array.
	 * @return Array of live keys in backing-table order.
	 */
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

	/** Returns the mutable backing keys, including invalid sentinels. */
	@Override
	public int[] keys(){return keys;}

	/** Returns the mutable backing values, including empty-cell values. */
	@Override
	public int[] values(){return values;}

	/** Returns the current empty-cell sentinel. */
	@Override
	public int invalid(){return invalid;}

	/** Returns the number of live mappings. */
	@Override
	public int size(){return size;}

	/** Returns true when the map has no live mappings. */
	@Override
	public boolean isEmpty(){return size==0;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Backing keys; empty cells contain {@link #invalid}. */
	private int[] keys;
	/** Backing values parallel to {@link #keys}; empty cells contain zero. */
	private int[] values;
	/** Number of live mappings. */
	private int size=0;
	/** Current empty-cell sentinel for {@link #keys}. */
	private int invalid;
	/** Mask for routing keys into the power-of-two primary table range. */
	private int modulus;
	/** Resize threshold derived from {@link #modulus} and {@link #loadFactor}. */
	private int sizeLimit;
	/** Maximum occupancy fraction for the primary table range. */
	private final float loadFactor;

	/** Deterministic source for regenerating {@link #invalid} sentinels. */
	private static final Random randy=shared.Shared.random(1);

}
