package structures;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * A binary heap that tracks element locations within the heap array.
 * Elements must implement SetLoc to maintain their position information.
 * In rollover mode, a full heap rejects candidates no greater than the current
 * minimum and evicts that minimum before accepting larger candidates.
 *
 * @author Brian Bushnell
 * @param <T> Element type that can track its location in the heap
 */
public final class HeapLoc<T extends SetLoc<? super T>> implements Iterable<T>{

	/**
	 * Constructs a heap with a fixed logical capacity.
	 * A good value for maxSize is {@code 2^N-1}.
	 *
	 * With assertions enabled, rejects nonpositive capacity before allocation
	 * because rollover compares candidates against the current root element.
	 *
	 * @param maxSize Maximum number of elements; must be positive
	 * @param rollover_ If true, a full heap rejects candidates no greater than
	 * the current minimum and evicts the current minimum for larger candidates
	 */
	@SuppressWarnings("unchecked")
	public HeapLoc(int maxSize, boolean rollover_){
		assert(maxSize>0) :
			"HeapLoc capacity must be positive to retain a root element: "+maxSize;
		
		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.
		
		CAPACITY=maxSize;
		array=(T[])new SetLoc[len];
		rollover=rollover_;
//		queue=new PriorityQueue<T>(maxSize);
	}
	
	/**
	 * Returns a replacement heap containing this heap's current elements.
	 * The returned heap owns the transferred elements and this heap is emptied.
	 * Each transferred element is assigned its new location in the returned heap.
	 * With assertions enabled, rejects nonpositive replacement capacity and
	 * capacity smaller than {@link #size()} before mutating recorded locations.
	 * The replacement may shrink spare capacity but must still retain every
	 * current element; rollover mode does not make undersized transfer safe here.
	 *
	 * @param newCapacity Replacement capacity; must retain all current elements
	 * @return New heap containing the transferred elements
	 */
	public HeapLoc<T> resizeNew(int newCapacity){
		assert(newCapacity>=size) :
			"New HeapLoc capacity "+newCapacity+" cannot retain "+size+" elements";
		HeapLoc<T> heap=new HeapLoc<T>(newCapacity, rollover);
		//Technically I could just do a limited array copy,
		//since all the positions will remain the same.
		//TODO Add resizing the current one as an option.
		for(T t : this){
			t.setLoc(-1);
			heap.add(t);
		}
		assert(size()==heap.size());
		this.clear();
		return heap;
	}
	
	/** Adds an element if room exists or rollover accepts it.
	 * A non-rollover heap must not already be full.
	 * @param t Element with {@code loc()<0} before insertion
	 * @return true if the element was accepted */
	public boolean add(T t){
		return addAndReturnLocation(t)>=0;
	}
	
	/** Adds an element, records its array location and returns that location.
	 * A non-rollover heap must not already be full.
	 * @param t Element with {@code loc()<0} before insertion
	 * @return One-based heap-array location, or -1 if rollover rejects {@code t} */
	public int addAndReturnLocation(T t){
		assert(t.loc()<0);
		assert(size==0 || array[size]!=null);
		assert(rollover || size<CAPACITY);
		
		if(size>=CAPACITY){
			
			if(t.compareTo(array[1])<=0){return -1;}
			
			poll(); //Turns into a rolling buffer by removing smallest value.
			
//			{//This is a more efficient alternative to poll() and percDown(), but the result is slightly different.
//				array[1]=t;
//				percUp(1);
//				return true;
//			}
		}
		assert(size<CAPACITY);
		
		//assert(testForDuplicates());
//		assert(queue.size()==size);
//		queue.add(t);
		assert(size==0 || array[size]!=null);
		size++;
		array[size]=t;
		t.setLoc(size);
		return percDown(size);
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
	}
	
	/** Returns the current minimum without removing it.
	 * @return Minimum element, or null if the heap is empty */
	public T peek(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return null;}
//		assert(array[1]==queue.peek()) : size+", "+queue.size()+"\n"+
//			array[1]+"\n"+
//			array[2]+" , "+array[3]+"\n"+
//			array[4]+" , "+array[5]+" , "+array[6]+" , "+array[7]+"\n"+
//			queue.peek()+"\n";
		//assert(testForDuplicates());
		return array[1];
	}
	
	/** Removes the current minimum and marks it outside this heap with loc -1.
	 * @return Minimum element, or null if the heap is empty */
	public T poll(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return null;}
		T t=array[1];
//		assert(t==queue.poll());
		array[1]=array[size];
		array[size]=null;
		size--;
		if(size>0){
			array[1].setLoc(1);
			percUp(1);
		}
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		t.setLoc(-1);
		return t;
	}
	
	/** Restores heap order after an element's sort value changes.
	 * @param t Element still stored at its recorded location
	 * @return New one-based heap-array index for {@code t} */
	public int jiggle(T t){
		final int loc=t.loc();
		assert(array[loc]==t);
		int x=percDown(loc);
		if(x!=loc){return x;}
		return percUp(loc);
	}

	/** Moves an element toward the root while it is smaller than its parent.
	 * @param t Element still stored at its recorded location
	 * @return New one-based heap-array index for {@code t} */
	public int jiggleDown(T t){
		final int loc=t.loc();
		assert(array[loc]==t);
		return percDown(loc);
	}

	/** Moves an element toward the leaves while it is larger than its smaller child.
	 * @param t Element still stored at its recorded location
	 * @return New one-based heap-array index for {@code t} */
	public int jiggleUp(T t){
		final int loc=t.loc();
		assert(array[loc]==t);
		return percUp(loc);
	}

	/** Moves the element at {@code loc} toward the root despite the historical name.
	 * @param loc One-based heap-array index to repair
	 * @return New one-based heap-array index for the moved element */
	private int percDown(int loc){
		//assert(testForDuplicates());
		assert(loc>0);
		if(loc==1){return loc;}
		int next=loc/2;
		T a=array[loc];
		T b=array[next];
		assert(a!=b && a!=null);
		if(a.compareTo(b)<0){
			array[next]=a;
			array[loc]=b;
			a.setLoc(next);
			b.setLoc(loc);
			return percDown(next);
		}
		return loc;
	}

	/** Moves the element at {@code loc} toward the leaves despite the historical name.
	 * @param loc One-based heap-array index to repair
	 * @return New one-based heap-array index for the moved element */
	private int percUp(int loc){
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		int next1=loc*2;
		int next2=next1+1;
		if(next1>size){return loc;}
		T a=array[loc];
		T b=array[next1];
		T c=array[next2];
		assert(a!=b);
		assert(b!=c);
		assert(b!=null && a!=null);
		//assert(testForDuplicates());
		if(c==null || b.compareTo(c)<1){
			if(a.compareTo(b)>0){
				array[next1]=a;
				array[loc]=b;
				a.setLoc(next1);
				b.setLoc(loc);
				//assert(testForDuplicates());
				return percUp(next1);
			}
		}else{
			if(a.compareTo(c)>0){
				array[next2]=a;
				array[loc]=c;
				a.setLoc(next2);
				c.setLoc(loc);
				//assert(testForDuplicates());
				return percUp(next2);
			}
		}
		return loc;
	}
	
	/** Checks whether the heap currently has no elements.
	 * @return true if this heap is empty */
	public boolean isEmpty(){
//		assert((size==0) == queue.isEmpty());
		return size==0;
	}
	
	/** Checks whether a non-rollover add can append without evicting.
	 * @return true if {@link #size()} is less than the logical capacity */
	public boolean hasRoom(){
		return size<CAPACITY;
	}
	
	/**
	 * Removes all elements from this heap.
	 * This does not rewrite the recorded location of former elements.
	 */
	public void clear(){
//		queue.clear();
		for(int i=1; i<=size; i++){array[i]=null;}
		size=0;
	}
	
	/** Returns the number of elements currently stored.
	 * @return Current heap size */
	public int size(){
		return size;
	}
	
	/** Returns the zero-based binary-tree tier for a one-based heap-array index.
	 * @param x One-based heap-array index
	 * @return Zero-based tree tier; -1 for index 0 */
	public static int tier(int x){
		int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}
	
	/** Runs a pairwise backing-array diagnostic over references and recorded locations.
	 * @return true if the checked slots have no duplicate non-null references
	 * or recorded-location mismatch */
	public boolean testForDuplicates(){
		for(int i=0; i<array.length; i++){
			for(int j=i+1; j<array.length; j++){
				if(array[i]!=null && array[i]==array[j]){return false;}
				if(array[i]!=null && array[i].loc()!=i){return false;}
			}
		}
		return true;
	}
	
	/** Drains this heap into a sorted list by repeatedly polling the current minimum.
	 * Removed elements are assigned loc -1.
	 * @return Elements in ascending heap order */
	public ArrayList<T> toList(){
		ArrayList<T> list=new ArrayList<T>(size);
		for(int i=0, lim=size; i<lim; i++){
			list.add(poll());
		}
		assert(isEmpty());
		return list;
	}
	
	/** Iterates over the occupied one-based backing-array slots in array order.
	 * @return Iterator over current backing-array contents */
	@Override
	public Iterator<T> iterator(){
		return new HeapIterator();
	}
	
	private class HeapIterator implements Iterator<T>{

		/** Checks if there are more elements to iterate over.
		 * @return true if more elements remain, false otherwise */
		@Override
		public boolean hasNext(){return loc<=size;}

		/** Returns the next element in iteration order and advances the iterator.
		 * @return Next element in the heap array */
		@Override
		public T next(){return array[loc++];}
		
		int loc=1;
		
	}
	
	/** One-based heap array; slot 0 is unused. */
	private final T[] array;
	/** Maximum number of elements accepted at once. */
	public final int CAPACITY;
	/** True to evict the current minimum when adding a larger element to a full heap. */
	public final boolean rollover;
	/** Number of occupied slots in {@link #array}. */
	private int size=0;
	
//	private PriorityQueue<T> queue;
	
}
