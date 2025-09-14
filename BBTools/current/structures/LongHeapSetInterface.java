package structures;

/**
 * Maintains a heap of unique values.
 * @author Brian Bushnell
 * @date July 6, 2016
 *
 */
public interface LongHeapSetInterface {
	
	public boolean add(long key);
	
	public int increment(long key, int incr);
	
	public void clear();
	
	public int size();
	
	public int capacity();
	
	public boolean hasRoom();
	
	public LongHeap heap();
	
	public long peek();

	public boolean contains(long key);
	
}

