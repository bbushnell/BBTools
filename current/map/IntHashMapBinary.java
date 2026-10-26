package map;

import java.util.Arrays;
import shared.Random;

import shared.KillSwitch;
import shared.Tools;

/**
 * Int-to-int map that sends each key to a power-of-two home bucket using raw key bits.
 * Collision resolution uses a linear probe across the raw table, including physical
 * overflow cells above the logical mask range, and then wraps to cell zero. Values
 * default to -1 on lookup miss and to 0 as the previous value for a fresh insertion.
 * Positive increment overflow is capped at {@link Integer#MAX_VALUE}; negative
 * underflow wraps through signed int truncation. Instances are unsynchronized.
 *
 * @author Brian Bushnell
 * @date June 8, 2017
 */
public final class IntHashMapBinary extends AbstractIntHashMap{

	/** Serialization version. */
	private static final long serialVersionUID=1L;

	/**
	 * Runs the inherited destructive benchmark on a binary int map.
	 * @param args Ignored
	 */
	public static void main(String[] args){
		IntHashMapBinary set=new IntHashMapBinary(32, 0.7f);
		test(set);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs a map with the default initial size and load factor. */
	public IntHashMapBinary(){
		this(256);
	}

	/**
	 * Constructs a map with default load factor 0.7.
	 * @param initialSize Positive requested logical size; rounded up to a power of two, capped at 2^30, and at least two
	 */
	public IntHashMapBinary(int initialSize){
		this(initialSize, 0.7f);
	}

	/**
	 * Constructs a map with caller-selected initial size and load factor.
	 * The load factor is clamped to [0.25,0.90], and the initial size is rounded
	 * up to a power of two, capped at 2^30, and clamped to at least two cells.
	 *
	 * @param initialSize Positive requested logical minimum size
	 * @param loadFactor_ Requested load in (0,1), clamped to [0.25,0.90] after validation
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
		// Clamp one-cell count-histogram instances before resize(): size2=1 would
		// make newModulus=0 and trip assert(newModulus>0). Zero still crashes
		// loud via assert(initialSize>0). [map/IntHashMapBinary#001 FIXED]
		resize(Tools.max(2, initialSize));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Removes all mappings, resets values to zero, and preserves the allocated raw arrays. */
	@Override
	public void clear(){
		if(size<1){return;}
		Arrays.fill(keys, invalid);
		Arrays.fill(values, 0);
		size=0;
	}

	/**
	 * Returns the value mapped from a key.
	 * @param key Key to look up
	 * @return Mapped value, or -1 if absent
	 */
	@Override
	public int get(int key){
		int cell=findCell(key);
		return cell<0 ? -1 : values[cell];
	}

	/**
	 * Stores a key/value pair.
	 *
	 * @param key Key to insert or overwrite
	 * @param value Value to store
	 * @return Previous mapped value, or 0 if absent
	 */
	@Override
	public int put(int key, int value){return set(key, value);}

	/**
	 * Stores a key/value pair, rerolling the sentinel if it collides with {@code key}.
	 *
	 * @param key Key to insert or overwrite
	 * @param value Value to store
	 * @return Previous mapped value, or 0 if absent
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
	 * Adds one to a key's value, inserting 1 if absent.
	 * @param key Key to insert or update
	 * @return New value, capped at {@link Integer#MAX_VALUE} on positive overflow
	 */
	@Override
	public int increment(int key){
		return increment(key, 1);
	}

	/**
	 * Adds a signed delta to a key's value, inserting the delta if absent.
	 * Positive overflow caps at {@link Integer#MAX_VALUE}; negative underflow wraps.
	 *
	 * @param key Key to insert or update
	 * @param incr Signed increment to add
	 * @return New value after overflow handling
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
	 * Removes a key and closes the following probe cluster.
	 * @param key Key to remove
	 * @return true if a live mapping was removed
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

	/**
	 * Rehashes the occupied cells after a deletion until the next empty cell.
	 * @param initial Cell that was just emptied
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
	 * Reinserts one occupied cell if its key can move earlier in its probe sequence.
	 * @param cell Occupied raw cell to test
	 * @return true if the mapping moved to a different cell
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

	/** Rerolls the empty-cell sentinel and rewrites all cells holding the old sentinel. */
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<keys.length; i++){
			if(keys[i]==old){
				keys[i]=invalid;
//				assert(values[i]==0); //TODO: slow
			}
		}
	}

	/**
	 * Locates a key by probing from {@code key & modulus}.
	 * @param key Key to search for
	 * @return Raw cell index if found, or -1 if absent
	 */
	@Override
	int findCell(final int key){
		if(key==invalid){return -1;}

		// `modulus` here is a power-of-two mask, capacity-1. This 2017
		// implementation probes raw key bits and wraps [initial,limit)->[0,initial);
		// IntHashMap2/3 later added Tools.hash32plus mixing.
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
	 * Locates a key or the empty cell that would receive it.
	 * @param key Non-sentinel key to search for
	 * @return Raw cell containing {@code key}, or a raw cell holding {@link #invalid}
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

	/** Doubles the logical power-of-two table size and rehashes live mappings. */
	private final void resize(){
		assert(size>=sizeLimit);
		resize(Tools.max(2, modulus+1)*2L);
	}

	/**
	 * Replaces raw storage with a specified logical power-of-two size and rehashes live mappings.
	 * @param size2 New logical power-of-two table size
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
	 * Copies live keys into a dense new array.
	 * @return Owned key snapshot in raw table iteration order
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

	/**
	 * Borrows the raw key array.
	 * @return Borrowed raw key array with {@link #invalid()} in empty cells
	 */
	@Override
	public int[] keys(){return keys;}

	/**
	 * Borrows the raw value array.
	 * @return Borrowed raw value array, index-aligned with {@link #keys()}
	 */
	@Override
	public int[] values(){return values;}

	/**
	 * Returns the current empty-cell sentinel.
	 * @return Current negative empty-cell sentinel
	 */
	@Override
	public int invalid(){return invalid;}

	/**
	 * Returns the number of live mappings.
	 * @return Number of live mappings
	 */
	@Override
	public int size(){return size;}

	/**
	 * Tests whether the map has no live mappings.
	 * @return true when the map has no live mappings
	 */
	@Override
	public boolean isEmpty(){return size==0;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Raw key table; empty cells contain {@link #invalid}. */
	private int[] keys;
	/** Raw int values, index-aligned with {@link #keys}. */
	private int[] values;
	/** Number of live mappings. */
	private int size=0;
	/** Negative empty-cell sentinel absent from the live key set. */
	private int invalid;
	/** Power-of-two logical mask used as {@code key & modulus}. */
	private int modulus;
	/** Resize threshold derived from the logical mask and load factor. */
	private int sizeLimit;
	/** Clamped requested load factor. */
	private final float loadFactor;

	/** Reproducible source for negative empty-cell sentinel candidates. */
	private static final Random randy=shared.Shared.random(1);

}
