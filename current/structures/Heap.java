package structures;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * A binary min-heap with optional rollover behavior.
 *
 * @param <T> Element type ordered by the heap
 */
public final class Heap<T extends Comparable<? super T>> implements Iterable<T>{

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
	public Heap(int maxSize, boolean rollover_){
		assert(maxSize>0) : "Heap capacity must be positive to retain a root element: "+
			maxSize;
		
		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.
		
		CAPACITY=maxSize;
		array=(T[])new Comparable[len];
		rollover=rollover_;
//		queue=new PriorityQueue<T>(maxSize);
	}
	
	/** Adds an element if room exists or rollover accepts it.
	 * A non-rollover heap must not already be full.
	 * @param t Element to add
	 * @return true if the element was accepted */
	public boolean add(T t){return addAndReturnLocation(t)>=0;}
	
	/** Adds an element and returns its backing-array location.
	 * A non-rollover heap must not already be full.
	 * @param t Element to add
	 * @return One-based heap-array location, or -1 if rollover rejects {@code t} */
	public int addAndReturnLocation(T t){
		
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
	
	/** Removes and returns the current minimum.
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
		if(size>0){percUp(1);}
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return t;
	}
	
	/** Restores heap order after an element's sort value changes.
	 * @param t Element still stored at {@code loc}
	 * @param loc Current one-based heap-array index for {@code t}
	 * @return New one-based heap-array index for {@code t} */
	public int jiggle(T t, int loc){
		assert(array[loc]==t);
		int x=percDown(loc);
		if(x!=loc){return x;}
		return percUp(loc);
	}

	/** Moves an element toward the root while it is smaller than its parent.
	 * @param t Element still stored at {@code loc}
	 * @param loc Current one-based heap-array index for {@code t}
	 * @return New one-based heap-array index for {@code t} */
	public int jiggleDown(T t, int loc){
		assert(array[loc]==t);
		return percDown(loc);
	}

	/** Moves an element toward the leaves while it is larger than its smaller child.
	 * @param t Element still stored at {@code loc}
	 * @param loc Current one-based heap-array index for {@code t}
	 * @return New one-based heap-array index for {@code t} */
	public int jiggleUp(T t, int loc){
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
				//assert(testForDuplicates());
				return percUp(next1);
			}
		}else{
			if(a.compareTo(c)>0){
				array[next2]=a;
				array[loc]=c;
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
	
	/** Checks whether any non-null element reference appears twice in the backing array.
	 * @return true if no duplicate non-null references are present */
	public boolean testForDuplicates(){
		for(int i=0; i<array.length; i++){
			for(int j=i+1; j<array.length; j++){
				if(array[i]!=null && array[i]==array[j]){return false;}
			}
		}
		return true;
	}
	
	/** Drains this heap into a sorted list by repeatedly polling the current minimum.
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

		/** Checks if more elements remain in the iteration.
		 * @return true if more elements exist, false otherwise */
		@Override
		public boolean hasNext(){return loc<=size;}

		/** Returns the next element in the iteration and advances the position.
		 * @return Next element in array order */
		@Override
		public T next(){return array[loc++];}
		
		int loc=1;
		
	}
	
	/** One-based heap array; slot 0 is unused. */
	private final T[] array;
	/** Maximum number of elements accepted at once. */
	private final int CAPACITY;
	/** True to evict the current minimum when adding a larger element to a full heap. */
	private final boolean rollover;
	/** Number of occupied slots in {@link #array}. */
	private int size=0;
	
//	private PriorityQueue<T> queue;
	
}
