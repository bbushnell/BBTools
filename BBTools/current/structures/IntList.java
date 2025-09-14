package structures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Random;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * High-performance integer list implementation optimized for memory efficiency.
 * Uses primitive int arrays instead of boxed Integer objects, providing significant
 * performance and memory benefits for large datasets.
 * 
 * Key advantages over ArrayList<Integer>:
 * - No boxing/unboxing overhead
 * - Lower memory footprint (4 bytes vs ~16 bytes per element)
 * - Reduced garbage collection pressure
 * - Specialized operations for sorted data
 * 
 * @author Brian Bushnell
 * @contributor Isla Winglet
 * @date Sep 20, 2014
 */
public final class IntList{
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Benchmark comparing IntList vs ArrayList vs LinkedList performance.
	 * Accepts an optional argument for list length. 
	 */
	public static void main(String[] args){
		int length=args.length>0 ? Integer.parseInt(args[0]) : 100000000;
		benchmark(length);
	}
	
	/** 
	 * Performance benchmark testing add, shuffle, and sort operations.
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
	
	/** Default constructor with initial capacity of 256 */
	public IntList(){this(256);}
	
	/** 
	 * Constructor with specified initial capacity.
	 * @param initial Initial array size (minimum 1)
	 */
	public IntList(int initial){
		initial=Tools.max(initial, 1);
		array=KillSwitch.allocInt1D(initial);
	}
	
	/** 
	 * Creates a deep copy of this IntList.
	 * @return New IntList with identical contents
	 */
	public IntList copy() {
		IntList copy=new IntList(size);
		copy.addAll(this);
		return copy;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/

	/** Clears the list (sets size to 0) */
	public IntList clear(){size=0; return this;}
	
	/** Clears list and zeros all array elements */
	public void clearFull(){
		Arrays.fill(array, 0);
		size=0;
	}
	
	/** 
	 * Sets value at specified location, expanding array if necessary.
	 * @param loc Index to set
	 * @param value Value to store
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
	 * @param value New value for last element
	 */
	public final void setLast(int value){
		assert(size>0);
		array[size-1]=value;
	}
	
	/** Increments value at location by 1 */
	public final void increment(int loc){increment(loc, 1);}
	
	/** 
	 * Increments value at location by specified amount.
	 * @param loc Index to increment
	 * @param value Amount to add
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
	 * @param value Value to subtract elements from
	 */
	public void subtractFrom(int value){
		for(int i=0; i<size; i++){
			array[i]=value-array[i];
		}
	}
	
	/** 
	 * Adds element to end of list, expanding if necessary.
	 * @param x Value to add
	 */
	public final void add(int x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	/** 
	 * Adds element only if different from last element (pseudo-set behavior).
	 * @param x Value to conditionally add
	 */
	public void addIfNotEqualToLast(int x) {
		if(size<1 || x!=array[size-1]) {add(x);}
	}
	
	/** 
	 * Adds element without bounds checking (for performance).
	 * @param x Value to add
	 */
	public final void addUnchecked(int x){
		array[size]=x;
		size++;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Bulk Operations       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Adds all elements from another IntList.
	 * @param counts Source IntList to copy from
	 */
	public void addAll(IntList counts) {
		final int[] array2=counts.array;
		final int size2=counts.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}
	
	/** Sorts the list in ascending order */
	public void sort() {
		if(size>1){Shared.sort(array, 0, size);}
	}
	
	/** Randomly shuffles the list elements */
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
	
	/** Reverses the order of elements */
	public void reverse() {
		if(size>1){Tools.reverseInPlace(array, 0, size);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Resizing           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * Expands internal array to accommodate more elements.
	 * @param size2 New minimum capacity
	 */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		final int size3=(int)Tools.min(Shared.MAX_ARRAY_LEN, size2);
		assert(size3>size) : "Overflow: "+size+", "+size2+" -> "+size3;
		array=KillSwitch.copyOf(array, size3);
	}
	
	/** 
	 * Sets the logical size of the list.
	 * @param size2 New size
	 */
	public final void setSize(final int size2) {
		if(size2>array.length){resize(size2);}
		size=size2;
	}
	
	/** 
	 * Shrinks internal array to match current size.
	 * @return This IntList for chaining
	 */
	public final IntList shrink(){
		if(size==array.length){return this;}
		array=KillSwitch.copyOf(array, size);
		return this;
	}
	
	/** 
	 * Calculates sum of all elements as long to prevent overflow.
	 * @return Sum of all elements
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
	 * @return Sum of all elements
	 */
	public final double sum(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	/** 
	 * Finds percentile value (assumes sorted data).
	 * @param fraction Percentile as fraction (0.0 to 1.0)
	 * @return Value at specified percentile
	 */
	public double percentile(double fraction){
		if(size<1){return 0;}
		int idx=percentileIndex(fraction);
		return array[idx];
	}
	
	/** 
	 * Finds index of percentile position.
	 * @param fraction Percentile as fraction (0.0 to 1.0)
	 * @return Index of percentile position
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
	
	/** Removes duplicates and shrinks to fit */
	public final void shrinkToUnique(){
		condense();
		shrink();
	}
	
	/** 
	 * Removes duplicate elements in-place (assumes sorted input).
	 * Maintains sorted order while keeping only unique values.
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
	 * Keeps only values that appear at least minCopies times.
	 * Clever algorithm that keeps the Nth copy, not the first copy.
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
	 * Gets value at specified location.
	 * @param loc Index to retrieve
	 * @return Value at index, or 0 if out of bounds
	 */
	public final int get(int loc){
		return(loc>=size ? 0 : array[loc]);
	}
	
	/** 
	 * Removes and returns last element.
	 * @return Last element value
	 */
	public final int pop() {
		size--;
		return array[size];
	}
	
	/** 
	 * Returns last element value without removing it.
	 * @return Last element value
	 */
	public int lastElement() {
		assert(size>0);
		return array[size-1];
	}
	
	/** 
	 * Returns last element without bounds checking.
	 * @return Last element value
	 */
	public final int lastElementUnchecked() {
		return array[size-1];
	}
	
	/** 
	 * Checks for duplicate values (slow O(n²) algorithm).
	 * @return True if duplicates exist
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
	 * Checks if list contains specified value.
	 * @param x Value to search for
	 * @return True if value is found
	 */
	public boolean contains(int x) {
		for(int i=0; i<size; i++){
			if(array[i]==x){return true;}
		}
		return false;
	}
	
	/** 
	 * Creates array copy of current contents.
	 * @return New array with list contents
	 */
	public int[] toArray(){
		return KillSwitch.copyOf(array, size);
	}
	
	/** 
	 * Extracts unique values and their occurrence counts.
	 * Assumes this list is sorted. Modifies this list to contain only unique values.
	 * @param counts Output list to store occurrence counts
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
	
	/** 
	 * Checks if list is sorted in ascending order (slow).
	 * @return True if sorted
	 */
	public boolean sorted(){
		for(int i=1; i<size; i++){
			if(array[i]<array[i-1]){return false;}
		}
		return true;
	}
	
	/** 
	 * Checks if all elements are unique (slow).
	 * @return True if no duplicates exist
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
	
	/** Returns current number of elements */
	public int size() {
		return size;
	}
	
	/** Returns true if list is empty */
	public boolean isEmpty() {
		return size<1;
	}
	
	/** Returns current array capacity */
	public int capacity() {
		return array.length;
	}
	
	/** Returns unused array capacity */
	public int freeSpace() {
		return array.length-size;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	/** 
	 * Returns string showing non-zero elements as (index, value) pairs.
	 * @return String representation as set view
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
	 * Returns string showing all elements in order.
	 * @return String representation as list view
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
	
	/** Returns minimum of two integers */
	private static final int min(int x, int y){return x<y ? x : y;}
	/** Returns maximum of two integers */
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Backing array for element storage */
	public int[] array;
	/** Current number of elements in the list */
	public int size=0;
}