package structures;

/**
 * Maintains a heap of unique long values with a fixed capacity limit.
 * Combines a min-heap with a hash set to provide efficient top-N selection
 * while ensuring uniqueness. When at capacity, smaller values are automatically
 * evicted to make room for larger ones.
 *
 * @author Brian Bushnell
 * @date July 6, 2016
 */
public class LongHeapSet implements LongHeapSetInterface {
	
	public LongHeapSet(int limit_){
		limit=limit_;
		heap=new LongHeap(limit, true);
		set=new LongHashSet(limit*2);
	}
	
	@Override
	public boolean add(long value){
		assert(heap.size()==set.size());
		if(heap.hasRoom()){
			if(set.add(value)){
				heap.add(value);
				assert(heap.size()==set.size());
				return true;
			}
			assert(heap.size()==set.size());
			return false;
		}
		
		final long bottom=heap.peek();
		if(value>bottom){
			if(set.add(value)){
				set.remove(bottom);
				assert(set.size()<=limit);
				heap.add(value);
				assert(heap.size()<=limit);
				assert(heap.size()==set.size());
				return true;
			}
		}
		assert(heap.size()==set.size());
		return false;
	}
	
	/**
	 * Attempts to increment a key's count, but since this is a set structure,
	 * it simply tries to add the key and returns 1 if successful, 0 otherwise.
	 * Note: This implementation doesn't truly support incrementing.
	 *
	 * @param key The key to increment
	 * @param incr The increment amount (ignored)
	 * @return 1 if key was added, 0 if already present
	 */
	@Override
	public int increment(long key, int incr) {
		return add(key) ? 1 : 0; //Not quite correct...
	}
	
	/** Removes all elements from both the heap and set.
	 * After clearing, both structures will be empty and size() will return 0. */
	@Override
	public void clear(){
		assert(heap.size()==set.size()) : heap.size()+", "+set.size();
		if(heap.size()<1){
			assert(set.size()<1) : heap.size()+", "+set.size();
			return;
		}
		heap.clear();
		set.clear();
		assert(heap.size()==set.size());
	}
	
	public void add(LongHeapSet b){
		assert(heap.size()==set.size());
		final long[] array=b.heap.array();
		final int blen=b.heap.size();
		for(int i=1; i<=blen; i++){
			add(array[i]);
		}
		assert(heap.size()==set.size());
	}
	
	/**
	 * Returns the number of elements currently in the set.
	 * The heap and set sizes are maintained to be equal.
	 * @return Current number of elements
	 */
	@Override
	public int size(){
		assert(heap.size()==set.size());
		return heap.size();
	}
	
	@Override
	public boolean contains(long key){return set.contains(key);}

	/** Returns the maximum number of elements this set can hold */
	@Override
	public int capacity(){return heap.capacity();}
	/** Returns true if the set can accept more elements without eviction */
	@Override
	public boolean hasRoom(){return heap.hasRoom();}
	/** Returns the underlying heap structure */
	@Override
	public LongHeap heap(){return heap;}
	/** Returns the smallest value in the set without removing it */
	@Override
	public long peek(){return heap.peek();}
	
	final int limit;
	public LongHeap heap;
	public LongHashSet set;
	
}
