package map;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import shared.Random;

import shared.Shared;
import shared.Timer;
import structures.IntList;

/**
 * Hash set that tracks newly added integers for optimized bulk clearing.
 * Extends IntHashSet with an internal list to record added elements, enabling
 * faster clearing when the set is sparse relative to its capacity.
 * Not thread-safe; external synchronization is required for concurrent access.
 *
 * Add/clear only: the tracking list is maintained by add() and reset by clear(). remove() is NOT
 * supported -- it would desync the tracking list (which cannot be kept in sync without an O(n) scan),
 * so it is overridden to throw UnsupportedOperationException rather than silently corrupt
 * clear()/toArray()/verify(). SketchIndex uses the add/toArray/clear lifecycle.
 * See [map/IntHashSetList#001].
 *
 * @author Brian Bushnell
 * @date September 13, 2017
 */
public class IntHashSetList extends IntHashSet{

	/** Runs legacy assertion checks and a large fixed-size add/clear stress benchmark.
	 * @param args Unused */
	public static void main(String[] args){
		Random randy2=Shared.threadLocalRandom();
		IntHashSetList set=new IntHashSetList(20, 0.7f);
		HashSet<Integer> set2=new HashSet<Integer>(20, 0.7f);
		ArrayList<Integer> list=new ArrayList<Integer>();
		ArrayList<Integer> list2=new ArrayList<Integer>();
		for(int i=0; i<1000; i++){
			assert(!set.contains(i));
			assert(!set2.contains(i));
			list.add(Integer.valueOf(i));
		}
		for(int i=0; i<1000; i++){
			int r=randy2.nextInt();
			list2.add(r);
		}

		for(int x : list){
			set.add(x);
			set2.add(x);
			assert(set.contains(x));
			assert(set2.contains(x));
		}

		for(int x : list){
			assert(set.contains(x));
			assert(set2.contains(x));
		}
		set.clear();//IntHashSetList is add/clear-only (remove() throws) -> validate via clear, not per-element remove
		set2.clear();
		assert(set.isEmpty());
		assert(set2.isEmpty());

		for(int x : list2){
			set.add(x);
			set2.add(x);
			assert(set.contains(x));
			assert(set2.contains(x));
		}

		for(int x : list2){
			assert(set.contains(x));
			assert(set2.contains(x));
		}
		set.clear();
		set2.clear();
		assert(set.isEmpty());
		assert(set2.isEmpty());

		int count=4000000;
		int runs=32;
		IntList ll=new IntList(count);
		for(int i=0; i<count; i++){ll.add(randy2.nextInt());}

		Shared.printMemory();
		Timer t=new Timer();
		for(int k=0; k<2; k++){
			System.err.println("IntHashSetList:");
			t.start();
			for(int i=0; i<runs; i++){
				final int[] y=ll.array;
				for(int z=0; z<count; z++){
					final int value=y[z];
					set.add(value);
					set.contains(value);
				}
				set.clear();//add/clear-only: clear each run instead of per-element remove (remove() throws)
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

	/** Creates an empty set with default capacity and tracking list. */
	public IntHashSetList(){
		super();
		list=new IntList(4);
	}

	/** Creates an empty set with specified initial capacity.
	 * @param initialSize Positive minimum prime modulus before adding slack cells */
	public IntHashSetList(int initialSize){
		super(initialSize);
		list=new IntList(4);
	}

	/**
	 * Creates an empty set with specified capacity and load factor.
	 * @param initialSize Positive minimum prime modulus before adding slack cells
	 * @param loadFactor_ Requested fraction in (0,1), clamped to [0.25,0.90]
	 */
	public IntHashSetList(int initialSize, float loadFactor_){
		super(initialSize, loadFactor_);
		list=new IntList(4);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Removes all values and resets insertion tracking without shrinking capacity.
	 * Uses whole-array clearing above one quarter of the load threshold; otherwise
	 * removes tracked values through the base implementation, repairing probe chains. */
	@Override
	public void clear(){
		// Check the add-only tracking invariant before either clear path; a desynced
		// list would corrupt both sparse clear and toArray().
		assert(size()==list.size) : list.size+", "+size()+"\n"+list.toString()+"\n"+Arrays.toString(toArray())+"\n"+Arrays.toString(super.toArray())+"\n";
		if(size()<1){return;}
		if(size()>0.25*sizeLimit()){
			super.clear();
		}else{
			for(int i=0; i<list.size; i++){
				boolean b=super.remove(list.get(i));
				assert(b) : list.get(i)+", "+list.size+", "+size()+"\n"+list.toString()+"\n"+Arrays.toString(toArray())+"\n";
			}
		}
		list.clear();
	}

	/** Adds and records an integer only when it was not already present.
	 * @param value Integer to insert
	 * @return True if newly added; false for a duplicate */
	@Override
	public boolean add(int value){
		// super.add returns true only for a new value, so the tracking list mirrors
		// the set exactly and never gains duplicates.
		if(super.add(value)){
			list.add(value);
			return true;
		}
		return false;
	}

	/** Rejects per-element removal, which would invalidate the add-only tracking list.
	 * @param value Unused; every removal request is rejected, including absent values
	 * @return Never returns normally
	 * @throws UnsupportedOperationException Always; use clear() or an IntHashSet */
	@Override
	public boolean remove(int value){
		// [map/IntHashSetList#001] The inherited remove() would silently desync
		// list. Keeping it synchronized would make remove() O(n), and the sole
		// caller sketch/SketchIndex is add/clear-only, so crash loud instead of
		// corrupting clear()/toArray()/verify(). Brian-approved 2026-06-22.
		throw new UnsupportedOperationException("IntHashSetList does not support remove() (it would desync the "
			+"tracking list); use add/clear cycles, or a plain IntHashSet if you need remove().");
	}

	/** Copies values in successful-insertion order since the last clear.
	 * @return Newly allocated dense array of the tracked values */
	@Override
	public int[] toArray(){
		// Return the tracking list in O(size) instead of scanning the whole backing array.
		int[] r=list.toArray();
		return r;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks tracker size/uniqueness followed by the base table's probe invariants.
	 * Intended for explicit diagnostics, not a constant-time health check.
	 * @return True when the tracking and base-table checks pass */
	@Override
	public boolean verify(){
		// list must mirror the set exactly before the base probe-integrity check.
		if(size()!=list.size){return false;}
		if(list.containsDuplicates()){return false;}
		return super.verify();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Values added since the last clear; add()/clear() keep it synchronized with the set. */
	private IntList list;

}
