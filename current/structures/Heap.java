package structures;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * Binary min-heap implementation with optional rollover behavior.
 * Maintains heap property where parent elements are smaller than their children.
 * Supports fixed-size mode or rollover mode that removes smallest element when capacity is exceeded.
 *
 * @param <T> Element type that must be comparable
 * @author Brian Bushnell
 */
public final class Heap<T extends Comparable<? super T>> implements Iterable<T> {
	
	//A good value for maxSize would be (2^N)-1
	/**
	 * Constructs a heap with specified capacity and rollover behavior.
	 * Array size is adjusted to be even for optimization.
	 * @param maxSize Maximum number of elements the heap can hold
	 * @param rollover_ If true, removes smallest element when capacity exceeded; if false, rejects new elements
	 */
	@SuppressWarnings("unchecked")
	public Heap(int maxSize, boolean rollover_){
		
		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.
		
		CAPACITY=maxSize;
		array=(T[])new Comparable[len];
		rollover=rollover_;
//		queue=new PriorityQueue<T>(maxSize);
	}
	
	/**
	 * Adds an element to the heap.
	 * @param t Element to add
	 * @return true if element was added, false if rejected (when at capacity without rollover)
	 */
	public boolean add(T t){return addAndReturnLocation(t)>=0;}
	
	/**
	 * Adds an element to the heap and returns its final position.
	 * In rollover mode, removes smallest element if at capacity and new element is larger.
	 * @param t Element to add
	 * @return Final array index of the added element, or -1 if element was rejected
	 */
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
	
	/** Returns the smallest element without removing it.
	 * @return Smallest element in the heap, or null if heap is empty */
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
	
	/**
	 * Removes and returns the smallest element.
	 * Maintains heap property by moving last element to root and percolating up.
	 * @return Smallest element from heap, or null if heap is empty
	 */
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
	
	/** Returns the new location */
	public int jiggle(T t, int loc){
		assert(array[loc]==t);
		int x=percDown(loc);
		if(x!=loc) {return x;}
		return percUp(loc);
	}

	/** Returns the new location */
	public int jiggleDown(T t, int loc){
		assert(array[loc]==t);
		return percDown(loc);
	}

	/** Returns the new location */
	public int jiggleUp(T t, int loc){
		assert(array[loc]==t);
		return percDown(loc); //Possible bug: Should be percUp(loc)
	}

	/** Returns the new location */
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

	/** Returns the new location */
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
	
	/** Checks if the heap contains no elements.
	 * @return true if heap is empty, false otherwise */
	public boolean isEmpty(){
//		assert((size==0) == queue.isEmpty());
		return size==0;
	}
	
	/** Checks if the heap can accept more elements without rollover.
	 * @return true if size is less than capacity, false otherwise */
	public boolean hasRoom(){
		return size<CAPACITY;
	}
	
	/** Removes all elements from the heap and resets size to zero.
	 * Nulls out array references for garbage collection. */
	public void clear(){
//		queue.clear();
		for(int i=1; i<=size; i++){array[i]=null;}
		size=0;
	}
	
	/** Returns the number of elements currently in the heap.
	 * @return Current heap size */
	public int size(){
		return size;
	}
	
	/**
	 * Calculates the tier (depth level) of a position in the binary heap.
	 * Uses bit manipulation to find the position of the most significant bit.
	 * @param x Position in the heap array
	 * @return Tier level, with root at tier 0
	 */
	public static int tier(int x){
		int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}
	
	/**
	 * Tests if the heap contains duplicate object references.
	 * Used for debugging and validation purposes.
	 * @return true if no duplicate references found, false if duplicates exist
	 */
	public boolean testForDuplicates(){
		for(int i=0; i<array.length; i++){
			for(int j=i+1; j<array.length; j++){
				if(array[i]!=null && array[i]==array[j]){return false;}
			}
		}
		return true;
	}
	
	/**
	 * Converts heap to sorted ArrayList by repeatedly polling elements.
	 * Empties the heap in the process since poll() removes elements.
	 * @return ArrayList containing all heap elements in sorted order
	 */
	public ArrayList<T> toList(){
		ArrayList<T> list=new ArrayList<T>(size);
		for(int i=0, lim=size; i<lim; i++){
			list.add(poll());
		}
		assert(isEmpty());
		return list;
	}
	
	@Override
	public Iterator<T> iterator() {
		return new HeapIterator();
	}
	
	/** Iterator implementation for traversing heap elements in array order.
	 * Does not provide sorted traversal - elements are returned in internal array order. */
	private class HeapIterator implements Iterator<T> {

		@Override
		public boolean hasNext() {return loc<=size;}

		@Override
		public T next() {return array[loc++];}
		
		/** Current position in the heap array for iteration */
		int loc=1;
		
	}
	
	/** Internal array storing heap elements with 1-based indexing */
	private final T[] array;
	/** Maximum number of elements the heap can hold */
	private final int CAPACITY;
	/** Whether to remove smallest element when capacity exceeded */
	private final boolean rollover;
	/** Current number of elements in the heap */
	private int size=0;
	
//	private PriorityQueue<T> queue;
	
}
