package structures;

/**
 * Maintains a heap of unique values with associated counts using a size-limited structure.
 * Combines a LongHeap and LongHashMap to efficiently track the largest values while keeping counts.
 * When at capacity, smaller values are evicted to make room for larger ones.
 *
 * @author Brian Bushnell
 * @date August 3, 2017
 */
public class LongHeapMap implements LongHeapSetInterface {
	
	public LongHeapMap(int limit_){
		limit=limit_;
		heap=new LongHeap(limit, true);
		map=new LongHashMap(limit*2);
	}
	
	@Override
	public boolean add(long key){
		return increment(key, 1)>0;
	}
	
	/**
	 * Increments the count for a key, adding it to the heap if necessary.
	 * If at capacity, only accepts keys larger than the current minimum.
	 * When evicting, removes the smallest key to make room for larger ones.
	 *
	 * @param key The value to increment
	 * @param incr The increment amount (must be >= 1)
	 * @return The new count for the key, or 0 if rejected due to capacity limits
	 */
	@Override
	public int increment(long key, int incr){
		assert(incr>=1);
		assert(heap.size()==map.size());
		if(heap.hasRoom()){
			int value=map.increment(key, incr);
			if(value==incr){heap.add(key);}
			assert(heap.size()==map.size());
			return value;
		}
		
		final long bottom=heap.peek();
		if(key>bottom){
			int value=map.increment(key, incr);
			if(value==incr){
				map.remove(bottom);
				assert(map.size()<=limit);
				heap.add(key);
				assert(heap.size()<=limit);
				assert(heap.size()==map.size());
				return value;
			}
		}
		assert(heap.size()==map.size());
		return 0;
	}
	
	/** Removes all elements from both the heap and map structures.
	 * Maintains size consistency between internal data structures. */
	@Override
	public void clear(){
		assert(heap.size()==map.size()) : heap.size()+", "+map.size();
		if(heap.size()<1){
			assert(map.size()<1) : heap.size()+", "+map.size();
			return;
		}
		heap.clear();
		map.clear();
		assert(heap.size()==map.size());
	}
	
	public void add(LongHeapMap b){
		assert(heap.size()==map.size());
		final long[] keys=b.map.keys();
		final int[] values=b.map.values();
		final long invalid=b.map.invalid();
		
		for(int i=0; i<keys.length; i++){
			final long key=keys[i];
			final int value=values[i];
			assert((key==invalid)==(value==0));
			if(key!=invalid){
				increment(key, value);
			}
		}
		assert(heap.size()==map.size());
	}
	
	/** Returns the number of unique keys currently stored.
	 * @return The current size of the heap map */
	@Override
	public int size(){
		assert(heap.size()==map.size());
		return heap.size();
	}
	
	@Override
	public boolean contains(long key){return map.contains(key);}

	/** Returns the maximum number of unique keys this heap map can hold.
	 * @return The capacity limit */
	@Override
	public int capacity(){return heap.capacity();}
	/** Checks if the heap map can accept additional keys without eviction.
	 * @return true if below capacity, false if at capacity limit */
	@Override
	public boolean hasRoom(){return heap.hasRoom();}
	/** Returns the internal heap structure for direct access.
	 * @return The underlying LongHeap instance */
	@Override
	public LongHeap heap(){return heap;}
	/** Returns the smallest key in the heap without removing it.
	 * @return The minimum key value currently stored */
	@Override
	public long peek(){return heap.peek();}
	
	final int limit;
	public LongHeap heap;
	public LongHashMap map;
	
}
