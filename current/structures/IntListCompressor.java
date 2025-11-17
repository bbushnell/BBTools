package structures;

/**
 * Compression utility for integer lists that maintains uniqueness and sorted order.
 * Automatically triggers compression when the list becomes full and sufficient growth
 * has occurred since the last compression. Provides memory-efficient storage of
 * large integer collections by eliminating duplicates.
 *
 * @author Brian Bushnell
 * @date June 2025
 */
public final class IntListCompressor {
	
	/**
	 * Adds a value to the list and triggers compression if needed.
	 * Compression occurs when the list is full and has grown by at least 25%
	 * since the last compression operation.
	 * @param value The integer value to add to the list
	 */
	public void add(int value){
		list.add(value);
		if(list.freeSpace()==0 && lastCompression<0.75f*list.size()){
			sortAndShrink();
		}
	}
	
	/**
	 * Sorts the list and removes duplicate values to compress storage.
	 * Only performs compression if the list size has grown since the last
	 * compression operation. Updates the compression tracking counter.
	 */
	public void sortAndShrink(){
		if(lastCompression>=list.size()){return;}
		list.sort();
		list.shrinkToUnique();
		lastCompression=list.size();
	}
	
	/** The underlying integer list that stores the values */
	public IntList list=new IntList(4);
	/** Size of the list at the time of the last compression operation */
	private int lastCompression=0;
	
}
