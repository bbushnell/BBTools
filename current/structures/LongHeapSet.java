package structures;

import map.LongHashSet;

/**
 * Size-limited set of the largest unique long values.
 * Combines a {@link LongHeap} with a LongHashSet: the heap tracks the eviction
 * floor while the set enforces uniqueness. At capacity, an absent value is
 * accepted only when it is larger than the current minimum.
 *
 * @author Brian Bushnell
 * @date July 6, 2016
 */
public class LongHeapSet implements LongHeapSetInterface{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a set with the requested unique-value limit.
	 * @param limit_ Positive maximum number of retained unique values.
	 * With assertions enabled, nonpositive limits are rejected by the backing
	 * {@link LongHeap}.
	 */
	public LongHeapSet(int limit_){
		limit=limit_;
		heap=new LongHeap(limit, true);
		set=new LongHashSet(limit*2);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds an absent value according to the fixed-capacity largest-value policy.
	 * @param value Value to add.
	 * @return true if the value was newly retained.
	 */
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
	 * Adds the key with set semantics.
	 * This implementation does not track counts and ignores incr.
	 *
	 * @param key Key to add.
	 * @param incr Increment amount, ignored by this set implementation.
	 * @return 1 if key was added, 0 if it was already retained or rejected.
	 */
	@Override
	public int increment(long key, int incr){
		return add(key) ? 1 : 0; //Not quite correct...
	}

	/** Removes all values from both backing structures. */
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

	/**
	 * Merges another heap-set's retained values into this set.
	 * @param b Source heap-set.
	 */
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
	 * @return Current number of retained values.
	 */
	@Override
	public int size(){
		assert(heap.size()==set.size());
		return heap.size();
	}

	/**
	 * Tests whether a value is currently retained.
	 * @param key Value to search for.
	 * @return true if the value is retained.
	 */
	@Override
	public boolean contains(long key){return set.contains(key);}

	/**
	 * Returns the maximum number of unique values this set can retain.
	 * @return Capacity limit.
	 */
	@Override
	public int capacity(){return heap.capacity();}
	/**
	 * Tests whether a new value can be retained without eviction.
	 * @return true if below capacity.
	 */
	@Override
	public boolean hasRoom(){return heap.hasRoom();}
	/**
	 * Returns the internal heap structure for direct access.
	 * @return Underlying {@link LongHeap} instance.
	 */
	@Override
	public LongHeap heap(){return heap;}
	/**
	 * Returns the current eviction floor without removing it.
	 * @return Smallest retained value, or {@link LongHeap#EMPTY} when empty.
	 */
	@Override
	public long peek(){return heap.peek();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum number of unique values retained. */
	final int limit;
	/** Heap that tracks retained values and exposes the current eviction floor. */
	public LongHeap heap;
	/** Backing set enforcing uniqueness of retained values. */
	public LongHashSet set;

}
