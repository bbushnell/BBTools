package align2;

/**
 * Binary min-heap implementation for long values using 1-based array indexing.
 * Values are compared as signed longs; no fields are decoded during comparison.
 * This is not the active BBMapS walk heap. No callers were found in align2 at
 * revision 6e87b382. The legacy implementation reserves -1 as an empty sentinel
 * and its existing assertions reject some duplicate values; it is not a general
 * long priority queue without resolving the annotated limitations below.
 * Mutable and unsynchronized; callers must provide exclusive ownership.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public final class PackedHeap{

	//NOTE: PackedHeap is currently UNUSED (no instantiations anywhere in the tree as of 2026-06-20).
	//Kept as a long-value min-heap primitive. See [align2/PackedHeap#001] re: the value-duplicate asserts.

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs a PackedHeap with specified maximum capacity.
	 * Rounds the backing length to even so the last left-child slot has a right
	 * sibling slot. Array index 0 is unused. This is not an alignment guarantee.
	 * @param maxSize Maximum number of elements the heap can hold
	 */
	public PackedHeap(final int maxSize){
		//TODO: Probable bug - maxSize is unchecked; maxSize+1 and even rounding
		// can overflow. Reject negative or unrepresentable capacities if revived.

		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.

		CAPACITY=maxSize;
		array=new long[len];
//		queue=new PriorityQueue<T>(maxSize);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds a long value to the heap, maintaining min-heap property.
	 * The legacy percDown name means movement toward the root (standard sift-up).
	 * @param t Signed value; -1 is reserved and duplicate values can trip assertions
	 * @return true; this implementation does not return false when full
	 */
	public boolean add(final long t){
		//TODO: Probable bug - CAPACITY is never checked. Even rounding can admit
		// an extra value; an eventual bounds exception occurs after size is raised.
		// Reject a full heap before changing size if this dormant class is revived.
		//assert(testForDuplicates());
//		assert(queue.size()==size);
//		queue.add(t);
		assert(size==0 || array[size]!=-1L);
		size++;
		array[size]=t;
		percDown(size);
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return true;
	}

	/** Returns the minimum without removing it, or -1 if empty. */
	public long peek(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return -1L;}
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
	 * Replaces root with the last element and sifts toward the leaves. The cleared
	 * last slot supplies the -1 sentinel when the final parent has only one child.
	 * @return The minimum long value, or -1L if heap is empty
	 */
	public long poll(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return -1L;}
		final long t=array[1];
//		assert(t==queue.poll());
		array[1]=array[size];
		array[size]=-1L;
		size--;
		if(size>0){percUp(1);}
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return t;
	}

	/**
	 * Sifts toward the root after insertion, despite the historical method name.
	 * The loop may read unused array[0], but loc>1 prevents using it as a parent.
	 * @param loc Starting position for percolation (1-based index)
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
	 * Sifts toward the leaves recursively after poll decrements size.
	 * The missing right child is the freshly cleared array[size+1] sentinel.
	 * @param loc Starting position for percolation (1-based index)
	 */
	private void percUp(int loc){
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		final int next1=loc*2;
		final int next2=next1+1;
		if(next1>size){return;}
		final long a=array[loc];
		final long b=array[next1];
		final long c=array[next2];
		//WARNING [align2/PackedHeap#001]: a,b,c are long VALUES here (not object refs as in QuadHeap),
		//so assert(a!=b)/assert(b!=c) forbid DUPLICATE VALUES — a copy-paste artifact from the object heaps.
		//A value min-heap legitimately holds duplicates; if PackedHeap is ever revived with duplicate
		//scores these crash under -ea. Latent (class unused). //TODO: confirm intended usage before removing.
		assert(a!=b);
		assert(b!=c);
		assert(b!=-1L);
		//assert(testForDuplicates());
		if(c==-1L || b<=c){
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
	 * Unused iterative alternative to percUp with the same sentinel and duplicate
	 * assumptions. No performance advantage is established for this implementation.
	 * @param loc Starting position for percolation (1-based index)
	 */
	private void percUpIter(int loc){
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		final long a=array[loc];
		//assert(testForDuplicates());

		int next1=loc*2;
		int next2=next1+1;

		while(next1<=size){

			final long b=array[next1];
			final long c=array[next2];
			assert(a!=b);
			assert(b!=c);
			assert(b!=-1L);

			if(c==-1L || b<=c){
//			if(c==-1L || (b.site<c.site || (b.site==c.site && b.column<c.column))){
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

	/** Returns whether there are no active values. */
	public boolean isEmpty(){
//		assert((size==0) == queue.isEmpty());
		return size==0;
	}

	/** Resets the active size; backing values remain and are overwritten on reuse. */
	public void clear(){
//		queue.clear();
//		for(int i=1; i<=size; i++){array[i]=-1L;}
		size=0;
	}

	/** Returns the active value count, not the allocated backing length. */
	public int size(){
		return size;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Calculates the tier (floor of log2) of an integer value.
	 * Defined as the highest set-bit index: zero returns -1, negative ints return 31.
	 * @param x The integer value
	 * @return The tier value (31 - number of leading zeros)
	 */
	public static int tier(final int x){
		final int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}

	/**
	 * Debugging method that checks for duplicate values in the heap array.
	 * Performs O(array.length squared) comparison, including inactive storage.
	 * @return true if no duplicates found, false otherwise
	 */
	public boolean testForDuplicates(){
		//TODO: Probable bug - includes unused index 0, zero-initialized slots and
		// stale values after clear(). It can report duplicates in an empty heap.
		// Scan active indices 1..size only if this diagnostic is revived.
		for(int i=0; i<array.length; i++){
			for(int j=i+1; j<array.length; j++){
				if(array[i]!=-1L && array[i]==array[j]){return false;}
			}
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final long[] array;
	private final int CAPACITY;
	private int size=0;

}
