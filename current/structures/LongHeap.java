package structures;

import shared.KillSwitch;

/**
 * Fixed-capacity, 1-indexed min-heap for primitive long values.
 * When rollover is enabled, adding to a full heap replaces the current minimum
 * only if the new value is larger; without rollover, callers must stop before
 * capacity. Removing the current minimum marks the vacated last slot
 * {@link #EMPTY}; percolation treats an {@link #EMPTY} right child as absent.
 * Empty {@link #peek()} and {@link #poll()} calls also return {@link #EMPTY},
 * so callers cannot distinguish an empty heap from a stored
 * {@code Long.MIN_VALUE} using those methods alone. Retained payloads are
 * expected to be unique and not {@link #EMPTY}; {@link LongHeapSet} and
 * {@link LongHeapMap} enforce uniqueness for their heap use, and downward
 * percolation asserts distinct adjacent values plus a non-{@link #EMPTY} left
 * child while treating {@link #EMPTY} as a vacant right child.
 *
 * @author Brian Bushnell
 * @date June 30, 2016
 */
public final class LongHeap{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a rollover heap with the requested fixed capacity.
	 * @param maxSize Positive maximum number of stored values.
	 */
	public LongHeap(int maxSize){this(maxSize, true);}

	/**
	 * Creates a heap with the requested fixed capacity.
	 * @param maxSize Positive maximum number of stored values.
	 * @param rollover_ If true, additions to a full heap may evict the minimum.
	 * If false, callers must avoid adding when the heap is full.
	 */
	public LongHeap(int maxSize, boolean rollover_){
		assert(maxSize>0) :
				"LongHeap.add insertion and full-heap rollover require "+
				"positive maxSize, got "+maxSize;
		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.

		CAPACITY=maxSize;
		array=KillSwitch.allocLong1D(len);
		rollover=rollover_;
//		queue=new PriorityQueue<T>(maxSize);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds a value to the heap.
	 * @param t Value to add; must be distinct from retained values and not
	 * {@link #EMPTY}.
	 * @return true if stored, or false if rollover was enabled and the full heap
	 * already held a minimum greater than or equal to t.
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

	/**
	 * Returns the current minimum without removing it.
	 * @return Minimum value, or {@link #EMPTY} when this heap is empty.
	 */
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
	 * Removes and returns the current minimum.
	 * @return Removed minimum, or {@link #EMPTY} when this heap is empty.
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
	 * Moves the value at a 1-indexed location toward the root until it reaches
	 * the root or sits below a parent no larger than itself.
	 * @param loc Occupied 1-indexed heap location.
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
	 * Moves the value at a 1-indexed location toward the leaves until it is no
	 * larger than both children.
	 * @param loc Occupied 1-indexed heap location.
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
	 * Iterative variant of {@link #percUp(int)} for comparing implementations.
	 * @param loc Occupied 1-indexed heap location.
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

	/** Returns true when this heap has no live values. */
	public boolean isEmpty(){
//		assert((size==0) == queue.isEmpty());
		return size==0;
	}

	/** Returns true when this heap has reached its fixed capacity. */
	public boolean isFull(){
		return size==CAPACITY;
	}

	/** Returns true when a new value can be added without rollover eviction. */
	public boolean hasRoom(){
		return size<CAPACITY;
	}

	/** Clears the logical heap size without clearing the backing array. */
	public void clear(){
//		queue.clear();
//		for(int i=1; i<=size; i++){array[i]=EMPTY;}
		size=0;
	}

	/** Returns the current number of live heap values. */
	public int size(){
		return size;
	}

	/**
	 * Returns the zero-based tier containing a positive integer in a binary heap.
	 * @param x Positive integer.
	 * @return Highest set-bit index.
	 */
	public static int tier(int x){
		int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}

	/**
	 * Slow diagnostic uniqueness check over the full backing array.
	 * @return false if any non-{@link #EMPTY} value appears twice.
	 */
	public boolean testForDuplicates(){
		for(int i=0; i<array.length; i++){
			for(int j=i+1; j<array.length; j++){
				if(array[i]!=EMPTY && array[i]==array[j]){return false;}
			}
		}
		return true;
	}

	/**
	 * Exposes the 1-indexed backing heap array.
	 * Live values occupy indexes {@code 1..size()}.
	 * @return Internal backing array.
	 */
	public long[] array(){return array;}

	/**
	 * Destructively drains this heap to an ascending list.
	 * @return Sorted values formerly held by this heap.
	 */
	public LongList toList(){
		LongList list=new LongList(size);
		for(int i=0, lim=size; i<lim; i++){
			list.add(poll());
		}
		assert(isEmpty());
		return list;
	}

	/**
	 * Returns the fixed value capacity.
	 * @return Maximum number of stored values.
	 */
	public int capacity(){return CAPACITY;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** 1-indexed heap storage; index 0 and non-live slots are not part of the heap. */
	private final long[] array;
	/** Maximum number of live values this heap can hold. */
	private final int CAPACITY;
	/** True to replace the current minimum when a larger value is added at capacity. */
	private final boolean rollover;
	/** Current number of live values in {@link #array}. */
	private int size=0;

	/** Sentinel for empty results and absent right children; not a safe payload value. */
	public static final long EMPTY=Long.MIN_VALUE;

}
