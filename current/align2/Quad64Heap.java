package align2;

/**
 * High-performance binary min-heap implementation for Quad64 objects providing efficient
 * priority queue operations in alignment algorithms. Manages prioritized Quad64 alignment
 * candidates using 1-indexed binary heap with cache-optimized memory layout.
 *
 * Uses traditional heap layout starting at index 1 with efficient parent/child navigation
 * (parent=i/2, left=2i, right=2i+1). Forces even array length for optimal cache line
 * alignment and memory access patterns.
 *
 * @author Brian Bushnell
 * @date December 19, 2013
 */
public final class Quad64Heap {
	
	/**
	 * Constructs a new Quad64Heap with specified maximum capacity.
	 * Forces even array length for optimal cache line alignment.
	 * @param maxSize Maximum number of elements the heap can contain
	 */
	public Quad64Heap(int maxSize){
		
		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.
		
		CAPACITY=maxSize;
		array=new Quad64[len];
//		queue=new PriorityQueue<T>(maxSize);
	}
	
	/**
	 * Inserts element into heap maintaining min-heap property through percolate-up.
	 * O(log n) insertion operation.
	 * @param t The Quad64 element to add to the heap
	 * @return Always true (heap dynamically manages capacity)
	 */
	public boolean add(Quad64 t){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
//		queue.add(t);
		assert(size==0 || array[size]!=null);
		size++;
		array[size]=t;
		percDown(size);
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return true;
	}
	
	/**
	 * Returns minimum element without removing it from heap.
	 * O(1) operation accessing root element at index 1.
	 * @return Minimum Quad64 element, or null if heap is empty
	 */
	public Quad64 peek(){
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
	 * Removes and returns minimum element from heap. O(log n) removal operation
	 * with last-element replacement and percolate-down to maintain heap property.
	 * @return Minimum Quad64 element, or null if heap is empty
	 */
	public Quad64 poll(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return null;}
		Quad64 t=array[1];
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
	
//	private void percDownRecursive(int loc){
//		//assert(testForDuplicates());
//		assert(loc>0);
//		if(loc==1){return;}
//		int next=loc/2;
//		Quad64 a=array[loc];
//		Quad64 b=array[next];
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
//		final Quad64 a=array[loc];
//
//		while(loc>1){
//			int next=loc/2;
//			Quad64 b=array[next];
//			assert(a!=b);
//			if(a.compareTo(b)<0){
//				array[next]=a;
//				array[loc]=b;
//				loc=next;
//			}else{return;}
//		}
//	}
	
	/**
	 * Percolates element upward in heap to maintain min-heap property after insertion.
	 * Optimized upward percolation using while loop, moving smaller elements toward root.
	 * Used during add() operations.
	 * @param loc Index position of element to percolate up
	 */
	private void percDown(int loc){
		//assert(testForDuplicates());
		assert(loc>0);
		if(loc==1){return;}

		int next=loc/2;
		final Quad64 a=array[loc];
		Quad64 b=array[next];
		
//		while(loc>1 && (a.site<b.site || (a.site==b.site && a.column<b.column))){
		while(loc>1 && a.compareTo(b)<0){
			array[loc]=b;
			loc=next;
			next=next/2;
			b=array[next];
		}
			
		array[loc]=a;
	}
	
	/**
	 * Percolates element downward in heap to maintain min-heap property after removal.
	 * Recursively compares with children and swaps with smaller child if needed.
	 * Used during poll() operations.
	 * @param loc Index position of element to percolate down
	 */
	private void percUp(int loc){
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		int next1=loc*2;
		int next2=next1+1;
		if(next1>size){return;}
		Quad64 a=array[loc];
		Quad64 b=array[next1];
		Quad64 c=array[next2];
		assert(a!=b);
		assert(b!=c);
		assert(b!=null);
		//assert(testForDuplicates());
		if(c==null || b.compareTo(c)<1){
			if(a.compareTo(b)>0){
//			if((a.site>b.site || (a.site==b.site && a.column>b.column))){
				array[next1]=a;
				array[loc]=b;
				//assert(testForDuplicates());
				percUp(next1);
			}
		}else{
			if(a.compareTo(c)>0){
//			if((a.site>c.site || (a.site==c.site && a.column>c.column))){
				array[next2]=a;
				array[loc]=c;
				//assert(testForDuplicates());
				percUp(next2);
			}
		}
	}
	
	/**
	 * Iterative alternative to recursive percolate-down operation for performance
	 * optimization. Provides same functionality as percUp() but uses while loop
	 * instead of recursion to avoid stack overhead.
	 * @param loc Index position of element to percolate down iteratively
	 */
	private void percUpIter(int loc){
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		final Quad64 a=array[loc];
		//assert(testForDuplicates());

		int next1=loc*2;
		int next2=next1+1;
		
		while(next1<=size){
			
			Quad64 b=array[next1];
			Quad64 c=array[next2];
			assert(a!=b);
			assert(b!=c);
			assert(b!=null);
			
			if(c==null || b.compareTo(c)<1){
//			if(c==null || (b.site<c.site || (b.site==c.site && b.column<c.column))){
				if(a.compareTo(b)>0){
//				if((a.site>b.site || (a.site==b.site && a.column>b.column))){
//					array[next1]=a;
					array[loc]=b;
					loc=next1;
				}else{
					break;
				}
			}else{
				if(a.compareTo(c)>0){
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
	
	/** Checks if heap contains no elements.
	 * @return True if heap is empty, false otherwise */
	public boolean isEmpty(){
//		assert((size==0) == queue.isEmpty());
		return size==0;
	}
	
	/** Removes all elements from heap without array traversal.
	 * Resets size to 0 for efficient heap reset operation. */
	public void clear(){
//		queue.clear();
//		for(int i=1; i<=size; i++){array[i]=null;}
		size=0;
	}
	
	/** Returns current number of elements in heap.
	 * @return Number of elements currently stored in heap */
	public int size(){
		return size;
	}
	
	/**
	 * Calculates tier level based on bit position of highest set bit.
	 * Uses Integer.numberOfLeadingZeros for efficient bit manipulation.
	 * @param x Input integer value
	 * @return Tier level (31 - leading zeros count)
	 */
	public static int tier(int x){
		int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}
	
	/**
	 * Validation method for heap integrity verification during development.
	 * Checks for duplicate object references in the heap array.
	 * @return True if no duplicate references found, false if duplicates exist
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
	 * 1-indexed binary heap array storing Quad64 elements with even-sized allocation
	 */
	private final Quad64[] array;
	/** Maximum number of elements the heap can contain */
	private final int CAPACITY;
	/** Current number of elements stored in the heap */
	private int size=0;
	
}
