package structures;

import map.LongHashMap;

/**
 * Size-limited map of the largest unique long keys and their counts.
 * Combines a {@link LongHeap} with a LongHashMap: the heap tracks the eviction
 * floor while the map tracks retained-key counts. At capacity, an absent key is
 * accepted only when it is larger than the current minimum; retained keys keep
 * accumulating counts.
 *
 * This map is intended for serialized or single-owner access. Every backing map
 * entry for a retained key is expected to have a positive count.
 *
 * @author Brian Bushnell
 * @date August 3, 2017
 */
public class LongHeapMap implements LongHeapSetInterface{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a map with the requested unique-key limit.
	 * @param limit_ Positive maximum number of retained unique keys.
	 * With assertions enabled, nonpositive limits are rejected by the backing
	 * {@link LongHeap}.
	 */
	public LongHeapMap(int limit_){
		limit=limit_;
		heap=new LongHeap(limit, true);
		map=new LongHashMap(limit*2);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Increments a key by one and returns whether its updated count is positive.
	 * Existing retained keys increment even when this map is full; new keys are
	 * retained only if there is room or they exceed the current eviction floor.
	 * Capacity rejection returns false. With assertions disabled, unchecked count
	 * arithmetic can also make this return false for a retained key.
	 * @param key Key to increment.
	 * @return true if the updated count is positive.
	 */
	@Override
	public boolean add(long key){
		return increment(key, 1)>0;
	}

	/**
	 * Increments a key's count, adding it when retained.
	 * Below capacity, absent keys are inserted. At capacity, existing retained keys
	 * still increment, a new key above the floor evicts the current minimum, and a
	 * new key below the floor is rejected.
	 *
	 * @param key Key to increment.
	 * @param incr Positive increment amount.
	 * @return Updated retained count, or 0 on capacity rejection. With assertions
	 * enabled, retained-count overflow is detected after {@link LongHashMap}'s
	 * unchecked addition and the prior count is restored before AssertionError;
	 * with assertions disabled, int wrap can still produce zero or negative
	 * counts.
	 */
	@Override
	public int increment(long key, int incr){
		//TODO: With assertions disabled, unchecked int wrap can corrupt counts and make add() report false for retained keys.
		assert(incr>=1);
		assert(heap.size()==map.size());
		if(heap.hasRoom()){
			int value=map.increment(key, incr);
			assert(value>=incr || rollbackCountOverflow(key, incr)) :
					countOverflowMessage(key, incr, value);
			if(value==incr){heap.add(key);}
			assert(heap.size()==map.size());
			return value;
		}

		final long bottom=heap.peek();
		if(key>bottom){
			int value=map.increment(key, incr);
			assert(value>=incr || rollbackCountOverflow(key, incr)) :
					countOverflowMessage(key, incr, value);
			if(value==incr){//New key above the floor: evict the current minimum to make room
				map.remove(bottom);
				assert(map.size()<=limit);
				heap.add(key);
				assert(heap.size()<=limit);
				assert(heap.size()==map.size());
				return value;
			}
			assert(heap.size()==map.size());
			return value;//Existing member above the floor: count already bumped; return it (was incorrectly 0)
		}else if(key==bottom){//The minimum member must still accrue count when full (was silently dropped)
			int value=map.increment(key, incr);
			assert(value>=incr || rollbackCountOverflow(key, incr)) :
					countOverflowMessage(key, incr, value);
			assert(heap.size()==map.size());
			return value;
		}
		assert(heap.size()==map.size());
		return 0;//key<bottom: a genuinely new, too-small key; reject with no map op (preserves the hot reject path)
	}

	/**
	 * Removes all keys and counts from both backing structures.
	 */
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

	/**
	 * Merges another heap-map's retained keys and counts into this map.
	 * @param b Source heap-map.
	 */
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

	/**
	 * Returns the number of unique keys currently retained.
	 * @return Current unique-key count.
	 */
	@Override
	public int size(){
		assert(heap.size()==map.size());
		return heap.size();
	}

	/**
	 * Tests whether a key is currently retained.
	 * @param key Key to search for.
	 * @return true if the key is retained.
	 */
	@Override
	public boolean contains(long key){return map.contains(key);}

	/**
	 * Returns the maximum number of unique keys this map can retain.
	 * @return Capacity limit.
	 */
	@Override
	public int capacity(){return heap.capacity();}
	/**
	 * Tests whether a new key can be retained without eviction.
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
	 * @return Smallest retained key, or {@link LongHeap#EMPTY} when empty.
	 */
	@Override
	public long peek(){return heap.peek();}

	/**
	 * Rolls back a retained count after an overflowing increment.
	 * @return false so the enclosing assertion still throws.
	 */
	private boolean rollbackCountOverflow(final long key, final int incr){
		map.increment(key, -incr);
		return false;
	}

	/** Explains why a retained count failed the post-increment representability check. */
	private String countOverflowMessage(final long key, final int incr,
			final int value){
		return "key="+key+", wrapped="+value+", incr="+incr+
				"; positive retained counts must fit int after LongHashMap's "+
				"unchecked addition, which also supports LongHeapMap's new-key "+
				"detection.";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum number of unique keys retained. */
	final int limit;
	/** Heap that tracks retained keys and exposes the current eviction floor. */
	public LongHeap heap;
	/** Backing map from retained key to its current count. */
	public LongHashMap map;

}
