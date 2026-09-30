package align2;

/**
 * Fixed-array binary min-heap for worker-owned Quad64 hit-list cursors.
 * Orders by compareTo (site, then column), never equals. Entries must be nonnull,
 * distinct object references; equal coordinate values in separate objects are
 * allowed. Poll a cursor before changing its key and reinserting it.
 * Indexing is 1-based: parent=i/2, children=2i and 2i+1. Even backing length
 * leaves a sibling slot for the final left child; it does not guarantee cache
 * alignment. clear retains references for pooled reuse. Not thread-safe.
 *
 * @author Brian Bushnell
 * @date December 19, 2013
 */
public final class Quad64Heap{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs a new Quad64Heap with specified maximum capacity.
	 * Rounds maxSize+1 to an even backing length, including unused index0.
	 * @param maxSize Maximum number of elements the heap can contain
	 * @throws IllegalArgumentException if the backing-array length cannot be represented
	 */
	public Quad64Heap(final int maxSize){
		//One unused slot plus even rounding must fit before allocating the array.
		if(maxSize<0 || maxSize>Integer.MAX_VALUE-2){
			throw new IllegalArgumentException("Quad64Heap capacity cannot fit its 1-based, even-length backing array: "+maxSize);
		}

		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.

		CAPACITY=maxSize;
		array=new Quad64[len];
//		queue=new PriorityQueue<T>(maxSize);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Inserts element into heap maintaining min-heap property through percolate-up.
	 * O(log n) insertion operation.
	 * @param t Nonnull cursor not already present by reference
	 * @return true after insertion
	 * @throws IllegalStateException if full; the heap remains unchanged
	 * @throws NullPointerException if the cursor is null; the heap remains unchanged
	 */
	public boolean add(final Quad64 t){
		//Reject invalid insertions before changing size or any active cursor.
		if(t==null){throw new NullPointerException("Quad64Heap requires a nonnull hit-list cursor");}
		if(size>=CAPACITY){throw new IllegalStateException("Quad64Heap is full: size="+size+", capacity="+CAPACITY);}
		//assert(testForDuplicates());
//		assert(queue.size()==size);
//		queue.add(t);
		assert(size==0 || array[size]!=null) : "Every active heap slot must contain an index cursor; size="+size;
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
		final Quad64 t=array[1];
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
		//Despite the name, this sifts TOWARD THE ROOT (index 1); used by add(). Standard sift-up.
		//(QuadHeap's javadocs call this direction "down toward root" — same operation, opposite wording.)
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
		//Despite the name, this sifts TOWARD THE LEAVES; used by poll(). Bounds-safe: only called
		//post-decrement from poll(), so next2<=size+1 stays within the even-length array and
		//array[size+1] is the freshly-vacated null.
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		//Check for a leaf before doubling loc; large valid indices can overflow.
		if(loc>size/2){return;}
		final int next1=loc*2;
		final int next2=next1+1;
		final Quad64 a=array[loc];
		final Quad64 b=array[next1];
		final Quad64 c=array[next2];
		assert(a!=b) : "BBIndex5 supplies distinct pooled cursors; parent/left alias at "+loc;
		assert(b!=c) : "BBIndex5 supplies distinct pooled cursors; siblings alias at "+loc;
		assert(b!=null) : "Active left child must contain a cursor; index="+next1+", size="+size;
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
		//Dead code: never called (poll() uses percUp). Iterative sift-down kept as an unused alternative.
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		final Quad64 a=array[loc];
		//assert(testForDuplicates());

		//Only parents have children; bounding loc first prevents index overflow.
		while(loc<=size/2){
			final int next1=loc*2;
			final int next2=next1+1;

			final Quad64 b=array[next1];
			final Quad64 c=array[next2];
			assert(a!=b) : "BBIndex5 supplies distinct pooled cursors; parent/left alias at "+loc;
			assert(b!=c) : "BBIndex5 supplies distinct pooled cursors; siblings alias at "+loc;
			assert(b!=null) : "Active left child must contain a cursor; index="+next1+", size="+size;

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
	 * Resets size only; inactive backing slots retain their previous references. */
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
	public static int tier(final int x){
		final int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}

	/**
	 * Validation method for heap integrity verification during development.
	 * Checks only active slots for duplicate references; clear() retains inactive slots.
	 * This O(size²) diagnostic is not called by normal heap operations.
	 * @return True if no duplicate references found, false if duplicates exist
	 */
	public boolean testForDuplicates(){
		assert(size>=0 && size<=CAPACITY) : "Active cursor range must fit the declared heap capacity: "+size+" / "+CAPACITY;
		for(int i=1; i<=size; i++){
			for(int j=i+1; j<=size; j++){
				if(array[i]!=null && array[i]==array[j]){return false;}
			}
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * 1-indexed binary heap array storing Quad64 elements with even-sized allocation
	 */
	private final Quad64[] array;
	/** Maximum number of elements the heap can contain */
	private final int CAPACITY;
	/** Current number of elements stored in the heap */
	private int size=0;

}
