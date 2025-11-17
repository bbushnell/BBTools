package structures;

import java.util.ArrayList;
import java.util.Iterator;

/** A heap that tracks element location */
public final class HeapLoc<T extends SetLoc<? super T>> implements Iterable<T> {
	
	//A good value for maxSize would be (2^N)-1
	/**
	 * Constructs a new heap with specified capacity and rollover behavior.
	 * Array size is adjusted to be even for implementation efficiency.
	 * @param maxSize Maximum number of elements the heap can hold
	 * @param rollover_ If true, adding beyond capacity removes smallest element
	 */
	@SuppressWarnings("unchecked")
	public HeapLoc(int maxSize, boolean rollover_){
		
		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.
		
		CAPACITY=maxSize;
		array=(T[])new SetLoc[len];
		rollover=rollover_;
//		queue=new PriorityQueue<T>(maxSize);
	}
	
	/**
	 * Creates a new heap with different capacity and transfers all elements.
	 * Elements are removed from this heap and added to the new one.
	 * This heap is cleared after the transfer.
	 *
	 * @param newCapacity Capacity of the new heap
	 * @return New heap containing all elements from this heap
	 */
	public HeapLoc<T> resizeNew(int newCapacity){
		HeapLoc<T> heap=new HeapLoc<T>(newCapacity, rollover);
		//Technically I could just do a limited array copy,
		//since all the positions will remain the same.
		//TODO Add resizing the current one as an option.
		for(T t : this) {
			t.setLoc(-1);
			heap.add(t);
		}
		assert(size()==heap.size());
		this.clear();
		return heap;
	}
	
	/**
	 * Adds an element to the heap if there is room or rollover is enabled.
	 * @param t Element to add (must have location set to -1)
	 * @return true if element was added, false if heap is full and rollover disabled
	 */
	public boolean add(T t){
		return addAndReturnLocation(t)>=0;
	}
	
	/**
	 * Adds an element to the heap and returns its final location.
	 * In rollover mode, removes smallest element if heap is full and new element
	 * is larger than the current minimum.
	 *
	 * @param t Element to add (must have location set to -1)
	 * @return Final location of the added element, or -1 if not added
	 */
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
	
	/** Returns the minimum element without removing it from the heap.
	 * @return Minimum element, or null if heap is empty */
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
	 * Removes and returns the minimum element from the heap.
	 * Restructures the heap to maintain heap property after removal.
	 * @return Minimum element, or null if heap is empty
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
	
	/** Returns the new location */
	public int jiggle(T t){
		final int loc=t.loc();
		assert(array[loc]==t);
		int x=percDown(loc);
		if(x!=loc) {return x;}
		return percUp(loc);
	}

	/** Returns the new location */
	public int jiggleDown(T t){
		final int loc=t.loc();
		assert(array[loc]==t);
		return percDown(loc);
	}

	/** Returns the new location */
	public int jiggleUp(T t){
		final int loc=t.loc();
		assert(array[loc]==t);
		return percDown(loc); //Possible bug: should this call percUp() instead?
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
			a.setLoc(next);
			b.setLoc(loc);
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
	 * Sets all array elements to null to prevent memory leaks. */
	public void clear(){
//		queue.clear();
		for(int i=1; i<=size; i++){array[i]=null;}
		size=0;
	}
	
	/** Returns the number of elements currently in the heap.
	 * @return Current number of elements */
	public int size(){
		return size;
	}
	
	/**
	 * Calculates the tier (binary log) of a number.
	 * Returns 31 minus the number of leading zeros in the binary representation.
	 * @param x Input number
	 * @return Tier value (floor of log base 2)
	 */
	public static int tier(int x){
		int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}
	
	/**
	 * Tests the heap integrity by checking for duplicate elements and verifying
	 * that each element's stored location matches its actual array position.
	 * @return true if no duplicates found and all locations are correct
	 */
	public boolean testForDuplicates(){
		for(int i=0; i<array.length; i++){
			for(int j=i+1; j<array.length; j++){
				if(array[i]!=null && array[i]==array[j]){return false;}
				if(array[i]!=null && array[i].loc()!=i){return false;}
			}
		}
		return true;
	}
	
	/**
	 * Converts the heap to a sorted ArrayList by polling all elements.
	 * This empties the heap as a side effect.
	 * @return ArrayList containing all elements in sorted order
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
	
	/** Iterator implementation that traverses heap elements in array order.
	 * Does not provide sorted traversal. */
	private class HeapIterator implements Iterator<T> {

		@Override
		public boolean hasNext() {return loc<=size;}

		@Override
		public T next() {return array[loc++];}
		
		/** Current position in the iteration (1-based indexing) */
		int loc=1;
		
	}
	
	/** Internal array storing heap elements with 1-based indexing */
	private final T[] array;
	/** Maximum number of elements the heap can hold */
	public final int CAPACITY;
	/** Whether to remove smallest element when adding beyond capacity */
	public final boolean rollover;
	/** Current number of elements in the heap */
	private int size=0;
	
//	private PriorityQueue<T> queue;
	
}
