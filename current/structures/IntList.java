package structures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;

import map.IntHashSet;
import shared.Random;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;

public final class IntList{
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Benchmark comparing IntList vs ArrayList vs LinkedList performance.
	 * Tests creation, shuffling, and sorting operations with timing measurements.
	 * @param args Optional array length argument (default: 100,000,000)
	 */
	public static void main(String[] args){
		int length=args.length>0 ? Integer.parseInt(args[0]) : 100000000;
		benchmark(length);
	}
	
	/**
	 * Executes comprehensive performance comparison between list implementations.
	 * Tests IntList, ArrayList, and LinkedList with add, shuffle, sort operations.
	 * Measures execution time and memory usage for each implementation.
	 * @param length Number of elements to test with
	 */
	private static void benchmark(final int length){
		Timer t=new Timer();
		System.gc();
		
		{
			System.err.println("\nIntList:");
			Shared.printMemory();
			t.start();
			IntList list=new IntList();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			list=null;
			System.gc();
		}
		
		{
			System.err.println("\nArrayList:");
			Shared.printMemory();
			t.start();
			ArrayList<Integer> list=new ArrayList<Integer>();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			list=null;
			System.gc();
		}
		
		{
			System.err.println("\nLinkedList:");
			Shared.printMemory();
			t.start();
			LinkedList<Integer> list=new LinkedList<Integer>();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			list=null;
			System.gc();
		}
		
		{
			System.err.println("\nIntList:");
			Shared.printMemory();
			t.start();
			IntList list=new IntList();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time:      \t");
			t.start();
			System.gc();
			t.stop("GC Time:   \t");
			Shared.printMemory();
			t.start();
			list.shuffle();
			t.stop("Shuf Time:  \t");
			t.start();
			list.sort();
			t.stop("Sort Time: \t");
			list=null;
			System.gc();
		}
		
		{
			System.err.println("\nArrayList:");
			Shared.printMemory();
			t.start();
			ArrayList<Integer> list=new ArrayList<Integer>();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time:      \t");
			t.start();
			System.gc();
			t.stop("GC Time:   \t");
			Shared.printMemory();
			t.start();
			Collections.shuffle(list);
			t.stop("Shuf Time:  \t");
			t.start();
			Collections.sort(list);
			t.stop("Sort Time: \t");
			list=null;
			System.gc();
		}
		
		{
			System.err.println("\nLinkedList:");
			Shared.printMemory();
			t.start();
			LinkedList<Integer> list=new LinkedList<Integer>();
			for(int i=0; i<length; i++){
				list.add(i);
			}
			t.stop("Time:      \t");
			t.start();
			System.gc();
			t.stop("GC Time:   \t");
			Shared.printMemory();
			t.start();
			Collections.shuffle(list);
			t.stop("Shuf Time:  \t");
			t.start();
			Collections.sort(list);
			t.stop("Sort Time: \t");
			list=null;
			System.gc();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates IntList with default initial capacity of 256 elements */
	public IntList(){this(256);}
	
	/**
	 * Creates IntList with specified initial capacity.
	 * Minimum capacity is enforced to be at least 1.
	 * @param initial Initial array size (will be at least 1)
	 */
	public IntList(int initial){
		initial=Tools.max(initial, 1);
		array=KillSwitch.allocInt1D(initial);
	}
	
	private IntList(IntList x){
		size=x.size;
		array=Arrays.copyOf(x.array, size);
	}
	
	/**
	 * Creates a deep copy of this IntList with identical contents.
	 * The new IntList has the same size and element values.
	 * @return New IntList containing copies of all elements
	 */
	public IntList copy() {
		return new IntList(this);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Clears the list by setting size to 0.
	 * Does not modify array contents for performance.
	 * @return This IntList for method chaining
	 */
	public IntList clear(){size=0; return this;}
	
	/** Clears list and zeros all array elements for security.
	 * More thorough than clear() but slower performance. */
	public void clearFull(){
		Arrays.fill(array, 0);
		size=0;
	}
	
	public void fill(int value){
		Arrays.fill(array, 0, size, value);
	}
	
	/**
	 * Sets value at specified index, expanding array if necessary.
	 * Updates list size to include the set position.
	 * @param loc Index to set (will expand array if needed)
	 * @param value Value to store at location
	 */
	public final void set(int loc, int value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	/**
	 * Sets the last element to specified value.
	 * Requires list to be non-empty.
	 * @param value New value for the last element
	 */
	public final void setLast(int value){
		assert(size>0);
		array[size-1]=value;
	}
	
	/** Increments value at location by 1, expanding array if necessary.
	 * @param loc Index to increment */
	public final void increment(int loc){increment(loc, 1);}
	
	/**
	 * Increments value at location by specified amount.
	 * Expands array and updates size as needed.
	 * @param loc Index to increment
	 * @param value Amount to add to existing value
	 */
	public final void increment(int loc, int value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	/**
	 * Subtracts all elements from specified value (value - element).
	 * Modifies elements in-place using formula: array[i] = value - array[i].
	 * @param value Value to subtract elements from
	 */
	public void subtractFrom(int value){
		for(int i=0; i<size; i++){
			array[i]=value-array[i];
		}
	}
	
	/**
	 * Appends element to end of list, expanding capacity if needed.
	 * Doubles array size when expansion is required.
	 * @param x Value to add to the list
	 */
	public final void add(int x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	/**
	 * Adds element only if different from last element in list.
	 * Provides pseudo-set behavior to avoid consecutive duplicates.
	 * @param x Value to conditionally add
	 */
	public void addIfNotEqualToLast(int x) {
		if(size<1 || x!=array[size-1]) {add(x);}
	}
	
	/**
	 * Adds element without bounds checking for maximum performance.
	 * Caller must ensure sufficient capacity exists.
	 * @param x Value to add to the list
	 */
	public final void addUnchecked(int x){
		array[size++]=x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Bulk Operations       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Appends all elements from another IntList to this list.
	 * Elements are added in order from the source list.
	 * @param counts Source IntList to copy elements from
	 */
	public void addAll(IntList counts) {
		final int[] array2=counts.array;
		final int size2=counts.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}
	
	/** Sorts elements in ascending order using an efficient algorithm */
	public void sort() {
		if(size>1){Shared.sort(array, 0, size);}
	}
	
	/** Randomly shuffles all elements using Fisher-Yates algorithm */
	public void shuffle() {
		if(size<2){return;}
		Random randy=Shared.threadLocalRandom();
		for(int i=0; i<size; i++){
			int j=randy.nextInt(size);
			int temp=array[i];
			array[i]=array[j];
			array[j]=temp;
		}
	}
	
	/** Reverses the order of all elements in-place */
	public void reverse() {
		if(size>1){Tools.reverseInPlace(array, 0, size);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Resizing           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void ensureCapacity(int extra) {
		int free=array.length-size;
		if(free>=extra) {return;}
		resize(Math.max(size+extra, 2L*array.length+1));
	}
	
	/**
	 * Expands internal array to accommodate at least size2 elements.
	 * Ensures new capacity is within maximum array length limits.
	 * @param size2 New minimum capacity required
	 */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		final int size3=(int)Tools.min(Shared.MAX_ARRAY_LEN, size2);
		assert(size3>size) : "Overflow: "+size+", "+size2+" -> "+size3;
		array=KillSwitch.copyOf(array, size3);
	}
	
	/**
	 * Sets the logical size of the list, expanding array if necessary.
	 * Does not initialize new elements to any particular value.
	 * @param size2 New logical size of the list
	 */
	public final void setSize(final int size2) {
		if(size2>array.length){resize(size2);}
		size=size2;
	}
	
	/**
	 * Shrinks internal array to match current size exactly.
	 * Reduces memory usage by eliminating unused capacity.
	 * @return This IntList for method chaining
	 */
	public final IntList shrink(){
		if(size==array.length){return this;}
		array=KillSwitch.copyOf(array, size);
		return this;
	}
	
	public int maxIdx() {
		if(size<1) {return -1;}
		int max=array[0];
		int maxIdx=0;
		for(int i=1; i<size; i++) {
			if(array[i]>max) {
				max=array[i];
				maxIdx=i;
			}
			max=max(max, array[i]);
		}
		return maxIdx;
	}
	
	public int max() {
		int max=-Integer.MAX_VALUE;
		for(int i=0; i<size; i++) {max=max(max, array[i]);}
		return max;
	}
	
	public int min() {
		int min=Integer.MAX_VALUE;
		for(int i=1; i<size; i++) {min=min(min, array[i]);}
		return min;
	}
	
	public final float stdev(){
		if(size<2){return 0;}
		double sum=sum();
		double avg=sum/size;
		double sumdev2=0;
		for(int i=0; i<size; i++){
			double x=array[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return (float)Math.sqrt(sumdev2/size);
	}
	
	public final double mean(){
		return size<1 ? 0 : sum()/size;
	}
	
	public final double median(){
		if(size<1){return 0;}
		int idx=percentileIndex(0.5f);
		return array[idx];
	}
	
	public final int mode(){
		assert(sorted());
		return mode(array, size);
	}
	
	public final int modeUnsorted(){
		int[] copy=toArray();
		Shared.sort(copy);
		return mode(copy, copy.length);
	}
	
	public static final int mode(int[] array, int size){
		if(size<1){return 0;}
		int streak=1, bestStreak=0;
		int prev=array[0];
		int best=prev;
		for(int i=0; i<size; i++){
			int x=array[i];
			if(x==prev){streak++;}
			else{
				if(streak>bestStreak){
					bestStreak=streak;
					best=prev;
				}
				streak=1;
				prev=x;
			}
		}
		if(streak>bestStreak){
			bestStreak=streak;
			best=prev;
		}
		return best;
	}
	
	/**
	 * Calculates sum of all elements as long to prevent integer overflow.
	 * Uses long arithmetic for large sums that exceed int range.
	 * @return Sum of all elements as long
	 */
	public final long sumLong(){
		long sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	/**
	 * Calculates sum of all elements as double.
	 * Provides floating-point precision for calculations.
	 * @return Sum of all elements as double
	 */
	public final double sum(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	/**
	 * Finds percentile value assuming list is sorted by value frequency.
	 * Uses cumulative sum approach rather than position-based percentile.
	 * @param fraction Percentile as fraction (0.0 to 1.0)
	 * @return Value at specified percentile
	 */
	public double percentile(double fraction){
		if(size<1){return 0;}
		int idx=percentileIndex(fraction);
		return array[idx];
	}
	
	/**
	 * Finds index where cumulative sum reaches target percentile.
	 * Assumes sorted data and uses value-weighted percentile calculation.
	 * @param fraction Percentile as fraction (0.0 to 1.0)
	 * @return Index where percentile threshold is reached
	 */
	public int percentileIndex(double fraction){
		if(size<2){return size-1;}
		assert(sorted());
		double target=(sum()*fraction);
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
			if(sum>=target){
				return i;
			}
		}
		return size-1;
	}
	
	/** Removes duplicate elements and shrinks array to fit exactly */
	public final void shrinkToUnique(){
		condense();
		shrink();
	}
	
	/**
	 * Removes duplicate elements in-place from sorted list.
	 * Maintains sorted order while keeping only unique values.
	 * Optimized algorithm skips initial ascending sequence.
	 */
	public final void condense(){
		if(size<=1){return;}
		
		int i=0, j=1;
		for(; j<size && array[i]<array[j]; i++, j++){}//skip while strictly ascending 
		
		int dupes=0;
		for(; j<size; j++){//This only enters at the first nonascending pair
			int a=array[i], b=array[j];
			assert(a<=b) : "Unsorted: "+i+", "+j+", "+a+", "+b;
			if(b>a){
				i++;
				array[i]=b;
			}else{
				//do nothing
				dupes++;
				assert(a==b);
			}
		}
		assert(dupes==(size-(i+1)));
		assert(size>=(i+1));
		size=i+1;
	}
	
	/**
	 * Keeps only values appearing at least minCopies times in sorted list.
	 * Clever algorithm retains the Nth copy when threshold is reached.
	 * More efficient than first-copy retention for large datasets.
	 * @param minCopies Minimum occurrences required to retain value
	 */
	public final void condenseMinCopies(int minCopies){
		if(minCopies <= 1) { condense(); return; }
		if(size <= 1) { size = 0; return; }

		int writePos = 0;
		int currentCount = 1;

		for(int readPos = 1; readPos < size; readPos++) {
			if(array[readPos] == array[readPos-1]) {
				currentCount++;
				if(currentCount == minCopies) {
					array[writePos++] = array[readPos];
				}
			} else {
				currentCount = 1; // Reset for new value
			}
		}

		size = writePos;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Gets value at specified index with bounds checking.
	 * Returns 0 for out-of-bounds access instead of throwing exception.
	 * @param loc Index to retrieve
	 * @return Value at index, or 0 if index >= size
	 */
	public final int get(int loc){
		return(loc>=size ? 0 : array[loc]);
	}
	
	/**
	 * Removes and returns the last element from the list.
	 * Decreases list size by 1.
	 * @return Value of the removed last element
	 */
	public final int pop() {
		size--;
		return array[size];
	}
	
	/**
	 * Returns last element value without removing it.
	 * Requires list to be non-empty.
	 * @return Value of the last element
	 */
	public int lastElement() {
		assert(size>0);
		return array[size-1];
	}
	
	/**
	 * Returns last element without bounds checking for performance.
	 * Caller must ensure list is non-empty.
	 * @return Value of the last element
	 */
	public final int lastElementUnchecked() {
		return array[size-1];
	}
	
	/**
	 * Checks if list contains any duplicate values using O(n²) algorithm.
	 * Inefficient for large lists - consider sorting first for better performance.
	 * @return true if any duplicate values exist, false otherwise
	 */
	public boolean containsDuplicates(){
		for(int i=0; i<size; i++){
			for(int j=i+1; j<size; j++){
				if(array[i]==array[j]){return true;}
			}
		}
		return false;
	}
	
	/**
	 * Checks if list contains specified value using linear search.
	 * @param x Value to search for
	 * @return true if value is found, false otherwise
	 */
	public boolean contains(int x) {
		for(int i=0; i<size; i++){
			if(array[i]==x){return true;}
		}
		return false;
	}
	
	/**
	 * Creates new array containing copy of all list elements.
	 * Returned array length equals list size, not capacity.
	 * @return New int array with copies of all elements
	 */
	public int[] toArray(){
		return KillSwitch.copyOf(array, size);
	}
	
	/**
	 * Extracts unique values and their occurrence counts from sorted list.
	 * Modifies this list to contain only unique values in-place.
	 * Populates counts list with occurrence frequency for each unique value.
	 * @param counts Output list to receive occurrence counts for each unique value
	 */
	public void getUniqueCounts(IntList counts) {
		counts.size=0;
		if(size<=0){return;}

		int unique=1;
		int count=1;
		
		for(int i=1; i<size; i++){
			assert(array[i]>=array[i-1]);
			if(array[i]==array[i-1]){
				count++;
			}else{
				array[unique]=array[i];
				unique++;
				counts.add(count);
				count=1;
			}
		}
		if(count>0){
			counts.add(count);
		}
		size=unique;
		assert(counts.size==size);
	}
	
	/** Checks if list is sorted in ascending order using O(n) scan.
	 * @return true if elements are in non-decreasing order, false otherwise */
	public boolean sorted(){
		for(int i=1; i<size; i++){
			if(array[i]<array[i-1]){return false;}
		}
		return true;
	}
	
	/**
	 * Checks if all elements are unique using hash set for O(n) performance.
	 * More efficient than containsDuplicates() for large lists.
	 * @return true if no duplicate values exist, false otherwise
	 */
	public boolean unique(){
		if(size<2) {return true;}
		IntHashSet set=new IntHashSet(size*2);
		for(int i=0; i<size; i++) {
			int x=array[i];
			if(set.contains(x)) {return false;}
		}
		return true;
	}
	
	/** Returns current number of elements in the list */
	public int size() {
		return size;
	}
	
	/** Returns true if list contains no elements */
	public boolean isEmpty() {
		return size<1;
	}
	
	/** Returns current internal array capacity (total allocated space) */
	public int capacity() {
		return array.length;
	}
	
	/** Returns unused array capacity (capacity - size) */
	public int freeSpace() {
		return array.length-size;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns string representation showing all elements in list format */
	@Override
	public String toString(){
		return toStringListView();
	}
	
	/**
	 * Returns string showing non-zero elements as (index, value) pairs.
	 * Useful for sparse data visualization.
	 * @return String representation as set view with index-value pairs
	 */
	public String toStringSetView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<size; i++){
			if(array[i]!=0){
				sb.append(comma+"("+i+", "+array[i]+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	/**
	 * Returns string showing all elements in sequential order.
	 * Standard list representation format: [elem1, elem2, elem3].
	 * @return String representation as ordered list
	 */
	public String toStringListView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<size; i++){
				sb.append(comma+array[i]);
				comma=", ";
		}
		sb.append(']');
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns smaller of two integers */
	private static final int min(int x, int y){return x<y ? x : y;}
	/** Returns larger of two integers */
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Backing array for element storage, may have unused capacity */
	public int[] array;
	/** Current number of elements in the list (logical size) */
	public int size=0;
}