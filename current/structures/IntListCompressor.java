package structures;

/**
 * Compression utility for integer lists that periodically restores uniqueness
 * and sorted order. Additions append ordinary unsorted values and may create
 * duplicates until a compression pass sorts and shrinks the live prefix.
 * Concurrent mutation requires caller coordination; the direct callers
 * synchronize on this wrapper, and the public list can bypass its bookkeeping.
 *
 * @author Brian Bushnell
 * @date June 2025
 */
public final class IntListCompressor{

	/**
	 * Adds a value and compresses the backing list when it is full and the exact
	 * predicate {@code lastCompression < 0.75f*list.size()} indicates roughly
	 * one-third growth since the last compression pass.
	 * @param value Value to add to the list
	 */
	public void add(int value){
		list.add(value);
		if(list.freeSpace()==0 && lastCompression<0.75f*list.size()){
			sortAndShrink();
		}
	}

	/**
	 * Sorts, removes duplicate values and shrinks the backing array, unless the
	 * current size has not grown beyond the last compressed size. Callers use
	 * this as a finalization step before reading unique counts or hash codes.
	 */
	public void sortAndShrink(){
		if(lastCompression>=list.size()){return;}
		list.sort();
		list.shrinkToUnique();
		lastCompression=list.size();
	}

	/** The underlying integer list that stores the values. */
	public IntList list=new IntList(4);
	/** Size of the list at the time of the last compression operation. */
	private int lastCompression=0;

}
