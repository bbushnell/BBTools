package structures;

import shared.KillSwitch;

/**
 * @author Brian Bushnell
 * @date June 30, 2016
 *
 */
public final class LongHeap {
	
	/** Creates a LongHeap with rollover enabled.
	 * @param maxSize Maximum number of elements the heap can contain */
	public LongHeap(int maxSize){this(maxSize, true);}
	
	/**
	 * Creates a LongHeap with specified capacity and rollover behavior.
	 * Array size is automatically adjusted to be even for optimization.
	 * @param maxSize Maximum number of elements the heap can contain
	 * @param rollover_ If true, removes smallest element when capacity exceeded
	 */
	public LongHeap(int maxSize, boolean rollover_){
		
		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.
		
		CAPACITY=maxSize;
		array=KillSwitch.allocLong1D(len);
		rollover=rollover_;
//		queue=new PriorityQueue<T>(maxSize);
	}
	
	/**
	 * Adds an element to the heap, maintaining min-heap property.
	 * If heap is at capacity and rollover is enabled, removes smallest element first.
	 * Elements smaller than or equal to the current minimum are rejected when at capacity.
	 *
	 * @param t The value to add
	 * @return true if element was added, false if rejected due to capacity constraints
	 */
	public boolean add(long t){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
//		queue.add(t);
		assert(size==0 || array[size]!=EMPTY);
		assert(rollover || size<CAPACITY);
		
		if(size>=CAPACITY){
			
			if(t<=array[1]){return false;}
			
			poll(); //Turns into a rolling buffer by removing smallest value.
			
//			{//This is a more efficient alternative to poll() and percDown(), but the result is slightly different.
//				array[1]=t;
//				percUp(1);
//				return true;
//			}
		}
		assert(size<CAPACITY);
		
		size++;
		array[size]=t;
		percDown(size);
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return true;
	}
	
	/** Returns the minimum element without removing it.
	 * @return The smallest element in the heap, or EMPTY if heap is empty */
	public long peek(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return EMPTY;}
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
	 * Restores heap property by moving last element to root and percolating up.
	 * @return The smallest element that was removed, or EMPTY if heap was empty
	 */
	public long poll(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return EMPTY;}
		long t=array[1];
//		assert(t==queue.poll());
		array[1]=array[size];
		array[size]=EMPTY;
		size--;
		if(size>0){percUp(1);}
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return t;
	}
	
//	private void percDownRecursive(int loc){
//		//assert(testForDuplicates());
//		assert(loc>0);
//		if(loc==1){return;}
//		int next=loc/2;
//		long a=array[loc];
//		long b=array[next];
//		assert(a!=b);
//		if(a.compareTo(b)<0){
//			array[next]=a;
//			array[loc]=b;
//			percDown(next);
//		}
//	}
//
//	private void percDown_old(int loc){
//		//assert(testForDuplicates());
//		assert(loc>0);
//
//		final long a=array[loc];
//
//		while(loc>1){
//			int next=loc/2;
//			long b=array[next];
//			assert(a!=b);
//			if(a.compareTo(b)<0){
//				array[next]=a;
//				array[loc]=b;
//				loc=next;
//			}else{return;}
//		}
//	}
	
	/**
	 * Percolates an element down towards the root to maintain min-heap property.
	 * Used when inserting new elements. Continues until proper position is found.
	 * @param loc The array index of the element to percolate down
	 */
	private void percDown(int loc){
		//assert(testForDuplicates());
		assert(loc>0);
		if(loc==1){return;}

		int next=loc/2;
		final long a=array[loc];
		long b=array[next];
		
//		while(loc>1 && (a.site<b.site || (a.site==b.site && a.column<b.column))){
		while(loc>1 && a<b){
			array[loc]=b;
			loc=next;
			next=next/2;
			b=array[next];
		}
			
		array[loc]=a;
	}
	
	/**
	 * Percolates an element up towards leaves to maintain min-heap property.
	 * Used after removing root element. Recursively swaps with smaller child.
	 * @param loc The array index of the element to percolate up
	 */
	private void percUp(int loc){
		//assert(testForDuplicates());
//		assert(loc>0 && loc<=size+1) : loc+", "+size; //Allows use of more-efficient sketch creation, but gives different result...
		assert(loc>0 && loc<=size) : loc+", "+size;
		int next1=loc*2;
		int next2=next1+1;
		if(next1>size){return;}
		long a=array[loc];
		long b=array[next1];
		long c=array[next2];
		assert(a!=b);
		assert(b!=c);
		assert(b!=EMPTY);
		//assert(testForDuplicates());
		if(c==EMPTY || b<=c){
			if(a>b){
//			if((a.site>b.site || (a.site==b.site && a.column>b.column))){
				array[next1]=a;
				array[loc]=b;
				//assert(testForDuplicates());
				percUp(next1);
			}
		}else{
			if(a>c){
//			if((a.site>c.site || (a.site==c.site && a.column>c.column))){
				array[next2]=a;
				array[loc]=c;
				//assert(testForDuplicates());
				percUp(next2);
			}
		}
	}
	
	/**
	 * Iterative version of percUp for better performance.
	 * Percolates element towards leaves using loop instead of recursion.
	 * @param loc The array index of the element to percolate up
	 */
	private void percUpIter(int loc){
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		final long a=array[loc];
		//assert(testForDuplicates());

		int next1=loc*2;
		int next2=next1+1;
		
		while(next1<=size){
			
			long b=array[next1];
			long c=array[next2];
			assert(a!=b);
			assert(b!=c);
			assert(b!=EMPTY);
			
			if(c==EMPTY || b<=c){
//			if(c==EMPTY || (b.site<c.site || (b.site==c.site && b.column<c.column))){
				if(a>b){
//				if((a.site>b.site || (a.site==b.site && a.column>b.column))){
//					array[next1]=a;
					array[loc]=b;
					loc=next1;
				}else{
					break;
				}
			}else{
				if(a>c){
//				if((a.site>c.site || (a.site==c.site && a.column>c.column))){
//					array[next2]=a;
					array[loc]=c;
					loc=next2;
				}else{
					break;
				}
			}
			next1=loc*2;
			next2=next1+1;
		}
		array[loc]=a;
	}
	
	/** Returns true if the heap contains no elements */
	public boolean isEmpty(){
//		assert((size==0) == queue.isEmpty());
		return size==0;
	}
	
	/** Returns true if the heap is at maximum capacity */
	public boolean isFull(){
		return size==CAPACITY;
	}
	
	/** Returns true if the heap can accept more elements without rollover */
	public boolean hasRoom(){
		return size<CAPACITY;
	}
	
	/** Removes all elements from the heap, resetting size to zero */
	public void clear(){
//		queue.clear();
//		for(int i=1; i<=size; i++){array[i]=EMPTY;}
		size=0;
	}
	
	/** Returns the current number of elements in the heap */
	public int size(){
		return size;
	}
	
	/**
	 * Calculates the tier level of a value based on bit position.
	 * Uses leading zero count to determine magnitude tier.
	 * @param x The value to analyze
	 * @return The tier level (0-31) representing value magnitude
	 */
	public static int tier(int x){
		int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}
	
	/**
	 * Debug method that checks for duplicate values in the heap.
	 * Scans entire array to ensure no non-EMPTY values are duplicated.
	 * @return true if no duplicates found, false otherwise
	 */
	public boolean testForDuplicates(){
		for(int i=0; i<array.length; i++){
			for(int j=i+1; j<array.length; j++){
				if(array[i]!=EMPTY && array[i]==array[j]){return false;}
			}
		}
		return true;
	}
	
	/** Returns the underlying array used for heap storage */
	public long[] array(){return array;}
	
	/**
	 * Converts heap contents to a LongList by polling all elements.
	 * This operation empties the heap and returns elements in min-heap order.
	 * @return LongList containing all heap elements in sorted order
	 */
	public LongList toList(){
		LongList list=new LongList(size);
		for(int i=0, lim=size; i<lim; i++){
			list.add(poll());
		}
		assert(isEmpty());
		return list;
	}
	
	/** Returns the maximum capacity of the heap */
	public int capacity(){return CAPACITY;}
	
	/** Internal array storage for heap elements using 1-based indexing */
	private final long[] array;
	/** Maximum number of elements the heap can contain */
	private final int CAPACITY;
	/** If true, removes smallest element when adding to full heap */
	private final boolean rollover;
	/** Current number of elements in the heap */
	private int size=0;
	
	/** Sentinel value representing empty/null elements in the heap */
	public static final long EMPTY=Long.MIN_VALUE;
	
}
