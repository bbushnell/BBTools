package structures;

/**
 * Maintains a heap of unique values.
 * @author Brian Bushnell
 * @date July 6, 2016
 *
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
	
	@Override
	public int increment(long key, int incr) {
		return add(key) ? 1 : 0; //Not quite correct...
	}
	
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
	
	@Override
	public int size(){
		assert(heap.size()==set.size());
		return heap.size();
	}
	
	@Override
	public boolean contains(long key){return set.contains(key);}

	@Override
	public int capacity(){return heap.capacity();}
	@Override
	public boolean hasRoom(){return heap.hasRoom();}
	@Override
	public LongHeap heap(){return heap;}
	@Override
	public long peek(){return heap.peek();}
	
	final int limit;
	public LongHeap heap;
	public LongHashSet set;
	
}
