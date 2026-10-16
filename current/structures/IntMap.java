package structures;
import java.util.Arrays;

import shared.Tools;


/**
 * Array-backed integer-to-integer map for a dense inclusive key range.
 * Keys must be within {@link #min} and {@link #max}; unmapped entries store
 * {@link Integer#MIN_VALUE}, so that value cannot be stored as a payload.
 * This class is unsynchronized.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class IntMap{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Empty placeholder entry point.
	 * @param args Unused command-line arguments */
	public static void main(String[] args){

	}

	/** Constructs an IntMap for the inclusive key range {@code [from,to]}.
	 * @param from Minimum supported key, inclusive
	 * @param to Maximum supported key, inclusive */
	public IntMap(int from, int to){
		reset(from, to);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the value associated with a mapped key.
	 * @param key Key in the current supported range
	 * @return Associated value, or {@link Integer#MIN_VALUE} if absent */
	public int get(int key){
		assert(key>=min && key<=max);
		return array[key-min];
	}

	/** Tests whether a key has a non-sentinel mapping.
	 * @param key Key in the current supported range
	 * @return {@code true} if {@code key} is present */
	public boolean containsKey(int key){
		assert(key>=min && key<=max);
		return array[key-min]!=INVALID;
	}

	/** Adds 1 to the existing value for a key.
	 * @param key Present key in the current supported range
	 * @return New value after incrementing */
	public int increment(int key){
		return increment(key, 1);
	}

	/**
	 * Adds an increment to the existing value for a key.
	 * Only positive overflow saturates at {@link Integer#MAX_VALUE}.
	 * @param key Present key in the current supported range
	 * @param value Amount to add
	 * @return New value after incrementing
	 */
	public int increment(int key, int value){
		assert(key>=min && key<=max);
		int index=key-min;
		int old=array[index];
		int v2=(int)Tools.min(Integer.MAX_VALUE, old+(long)value);
		assert(array[index]!=INVALID);
		array[index]=v2;
		return v2;
	}

	/** Associates a value with a key.
	 * @param key Key in the current supported range
	 * @param value Value to store; must not be {@link Integer#MIN_VALUE}
	 * @return Previous value, or {@link Integer#MIN_VALUE} if absent */
	public int put(int key, int value){
		assert(key>=min && key<=max);
		assert(value!=INVALID);
		int index=key-min;
		int old=array[index];
		array[index]=value;
		return old;
	}

	/** Removes the mapping for a key.
	 * @param key Key in the current supported range
	 * @return Previous value, or {@link Integer#MIN_VALUE} if absent */
	public int remove(int key){
		assert(key>=min && key<=max);
		int index=key-min;
		int old=array[index];
		array[index]=INVALID;
		return old;
	}

	/**
	 * Returns the number of non-sentinel values in the retained backing capacity,
	 * which may be wider than the current logical range after a shrinking reset.
	 * @return Count of mapped entries */
	public int size(){
		int sum=0;
		for(int i=0; i<array.length; i++){
			if(array[i]!=INVALID){sum++;}
		}
		return sum;
	}

	/**
	 * Returns the keys for all non-sentinel entries in the retained backing capacity,
	 * which may be wider than the current logical range after a shrinking reset.
	 * @return Keys in ascending backing-index order */
	public int[] keys(){
		int[] r=new int[size()];
		for(int i=0, j=0; j<r.length; i++){
			if(array[i]!=INVALID){
				r[j]=(min+i);
				j++;
			}
		}
		return r;
	}

	/**
	 * Returns the values for all non-sentinel entries in the retained backing capacity,
	 * which may be wider than the current logical range after a shrinking reset.
	 * @return Values corresponding to {@link #keys()} order */
	public int[] values(){
		int[] r=new int[size()];
		for(int i=0, j=0; j<r.length; i++){
			if(array[i]!=INVALID){
				r[j]=array[i];
				j++;
			}
		}
		return r;
	}

	/** Clears every backing slot while retaining the current key range and capacity. */
	public void clear(){
		Arrays.fill(array, INVALID);
	}

	/**
	 * Resets this map to a new inclusive key range.
	 * Grows the backing array only when the retained capacity is too small, then
	 * clears every backing slot, including retained slots outside the new range.
	 * @param from New minimum supported key, inclusive
	 * @param to New maximum supported key, inclusive
	 */
	public void reset(int from, int to){
		min=from;
		max=to;
		assert(max>=min);
		assert(((long)max)-((long)min)<Integer.MAX_VALUE);

		int size=max-min+1;
		if(array==null || array.length<size){
			array=new int[size];
		}
		clear();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Minimum key value supported by the current logical range. */
	public int min;
	/** Maximum key value supported by the current logical range. */
	public int max;
	/** Backing array storing values indexed by {@code key-min}. */
	public int[] array;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Sentinel value indicating an unset or removed entry. */
	private static final int INVALID=Integer.MIN_VALUE;

}
