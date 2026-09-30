package align2;

/**
 * Fixed-array min-heap of worker-owned hit-list cursors, ordered by site then column.
 * Uses Quad.compareTo and reference identity, never value equality. Entries must
 * be nonnull, distinct object references; equal values in separate cursors are
 * legal. Poll before changing a cursor's ordering key, then reinsert it.
 * The index callers reuse cursor pools; clear resets size but retains references.
 * Unsynchronized and not resizable. Index0 is unused; backing length is rounded
 * even to leave a sibling slot for the final left child, not for cache alignment.
 * @author Brian Bushnell
 */
public final class QuadHeap{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs a QuadHeap with specified maximum capacity.
	 * Allocates an internal array with size maxSize+1 (rounded to even number).
	 * @param maxSize Maximum number of Quad objects this heap can hold
	 */
	public QuadHeap(final int maxSize){
		//One unused slot plus even rounding must fit before allocating the array.
		if(maxSize<0 || maxSize>Integer.MAX_VALUE-2){
			throw new IllegalArgumentException("QuadHeap capacity cannot fit its 1-based, even-length backing array: "+maxSize);
		}

		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.

		CAPACITY=maxSize;
		array=new Quad[len];
//		queue=new PriorityQueue<T>(maxSize);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds a cursor, sifting toward the root despite the percDown method name.
	 * @param t Nonnull cursor not already present by reference
	 * @return true after insertion
	 * @throws IllegalStateException if full; the heap remains unchanged
	 * @throws NullPointerException if the cursor is null; the heap remains unchanged
	 */
	public boolean add(final Quad t){
		//Reject invalid insertions before changing size or any active cursor.
		if(t==null){throw new NullPointerException("QuadHeap requires a nonnull hit-list cursor");}
		if(size>=CAPACITY){throw new IllegalStateException("QuadHeap is full: size="+size+", capacity="+CAPACITY);}
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
	 * Returns the minimum Quad object without removing it from the heap.
	 * The minimum is always at the root position (index 1).
	 * @return The minimum Quad object, or null if heap is empty
	 */
	public Quad peek(){
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
	 * Removes and returns the minimum Quad object from the heap.
	 * Replaces root with the last element and sifts toward the leaves iteratively.
	 * @return The minimum Quad object, or null if heap is empty
	 */
	public Quad poll(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return null;}
		final Quad t=array[1];
//		assert(t==queue.poll());
		array[1]=array[size];
		array[size]=null;
		size--;
		if(size>0){percUpIter(1);}
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
//		Quad a=array[loc];
//		Quad b=array[next];
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
//		final Quad a=array[loc];
//
//		while(loc>1){
//			int next=loc/2;
//			Quad b=array[next];
//			assert(a!=b);
//			if(a.compareTo(b)<0){
//				array[next]=a;
//				array[loc]=b;
//				loc=next;
//			}else{return;}
//		}
//	}

	/**
	 * Percolates element at given location down toward the root.
	 * Used when adding elements to maintain min-heap property by moving
	 * smaller elements closer to the root position.
	 * @param loc Index of element to percolate down
	 */
	private void percDown(int loc){
		//Despite the name, this sifts TOWARD THE ROOT (index 1); used by add(). Standard sift-up;
		//the loc>1 guard short-circuits the harmless array[0] read at the top.
		//assert(testForDuplicates());
		assert(loc>0) : "Heap uses 1-based slots; insertion index="+loc;
		if(loc==1){return;}

		int next=loc/2;
		final Quad a=array[loc];
		Quad b=array[next];

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
	 * Percolates element at given location up toward leaves recursively.
	 * Compares with children and swaps with smaller child to maintain heap order.
	 * Recursively continues down the affected subtree.
	 * @param loc Index of element to percolate up
	 */
	private void percUp(int loc){
		//Unused recursive alternative; poll uses percUpIter. Bounds-safe when called
		//post-decrement from poll(), so next2<=size+1 stays within the even-length array and
		//array[size+1] is the freshly-vacated null. Assumes distinct object refs (BBIndex Quad pool).
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		//Check for a leaf before doubling loc; large valid indices can overflow.
		if(loc>size/2){return;}
		final int next1=loc*2;
		final int next2=next1+1;
		final Quad a=array[loc];
		final Quad b=array[next1];
		final Quad c=array[next2];
		assert(a!=b) : "Index pools provide distinct cursors; parent/left alias at "+loc;
		assert(b!=c) : "Index pools provide distinct cursors; sibling references alias at "+loc;
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
	 * Iterative version of percolate up operation.
	 * Moves element down the tree iteratively rather than recursively,
	 * which may provide better performance for deep trees.
	 * @param loc Index of element to percolate up iteratively
	 */
	private void percUpIter(int loc){
		//Iterative sift-down used by poll(); preserves the recursive version's child and tie choices.
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		final Quad a=array[loc];
		//assert(testForDuplicates());

		//Only parents have children; bounding loc first prevents index overflow.
		while(loc<=size/2){
			final int next1=loc*2;
			final int next2=next1+1;

			final Quad b=array[next1];
			final Quad c=array[next2];
			assert(a!=b) : "Index pools provide distinct cursors; parent/left alias at "+loc;
			assert(b!=c) : "Index pools provide distinct cursors; sibling references alias at "+loc;
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

	public boolean isEmpty(){
//		assert((size==0) == queue.isEmpty());
		return size==0;
	}

	/** Resets the active size without clearing references retained by the backing array. */
	public void clear(){
//		queue.clear();
//		for(int i=1; i<=size; i++){array[i]=null;}
		size=0;
	}

	/** Returns the active cursor count, not the allocated backing length. */
	public int size(){
		return size;
	}

	/**
	 * Calculates the tier level of an integer based on bit position.
	 * Returns 31 minus the number of leading zeros, effectively
	 * computing floor(log2(x)) for the highest set bit position.
	 *
	 * @param x Integer to analyze
	 * @return Highest set bit index; -1 for zero, 31 for negative ints
	 */
	public static int tier(final int x){
		final int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}

	/**
	 * Tests heap integrity by checking for duplicate object references.
	 * Performs O(size²) comparison of active array elements to detect
	 * reference duplicates which would indicate heap corruption.
	 * @return true if no duplicate references found, false if duplicates exist
	 */
	public boolean testForDuplicates(){
		//clear() retains references outside 1..size; they are not heap members.
		assert(size>=0 && size<=CAPACITY) : "Active cursor range must fit the declared heap capacity: "+size+" / "+CAPACITY;
		for(int i=1; i<=size; i++){
			for(int j=i+1; j<=size; j++){
				if(array[i]!=null && array[i]==array[j]){return false;}
			}
		}
		return true;
	}

	/**
	 * Returns string representation of heap contents.
	 * Shows all elements from index 1 to size in array order,
	 * formatted as comma-separated list enclosed in brackets.
	 * @return String representation of heap elements
	 */
	@Override
	public String toString(){
		final StringBuilder sb=new StringBuilder();
		sb.append("[");
		for(int i=1; i<=size; i++){
			sb.append((i==1 ? "" : ", ")+array[i]);
		}
		sb.append("]");
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Quad[] array;
	private final int CAPACITY;
	private int size=0;

}
