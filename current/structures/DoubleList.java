package structures;

import java.util.ArrayList;
import java.util.LinkedList;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;



/**
 * Dynamic array of double values with resizing capability.
 * Provides memory-efficient storage and common operations like sorting,
 * statistical calculations, and deduplication for double precision numbers.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public final class DoubleList{
	
	/**
	 * Benchmarks DoubleList performance against ArrayList and LinkedList.
	 * Tests memory usage and insertion speed for large datasets.
	 * @param args Command-line arguments; first argument is list length (default 100,000,000)
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		int length=args.length>0 ? Integer.parseInt(args[0]) : 100000000;
		
		System.gc();
		
		{
			System.err.println("\nDoubleList:");
			Shared.printMemory();
			t.start();
			DoubleList list=new DoubleList();
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
			ArrayList<Double> list=new ArrayList<Double>();
			for(int i=0; i<length; i++){
				list.add((double)i);
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
			LinkedList<Double> list=new LinkedList<Double>();
			for(int i=0; i<length; i++){
				list.add((double)i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			list=null;
			System.gc();
		}
	}
	
	/** Creates a DoubleList with initial capacity of 256 elements */
	public DoubleList(){this(256);}
	
	/** Creates a DoubleList with specified initial capacity.
	 * @param initial Initial array capacity; minimum value is 1 */
	public DoubleList(int initial){
//		assert(initial>0) : initial+"\n"+this;
		initial=Tools.max(initial, 1);
		array=KillSwitch.allocDouble1D(initial);
	}
	
	/** Creates a deep copy of this DoubleList.
	 * @return New DoubleList containing same elements in same order */
	public DoubleList copy() {
		DoubleList copy=new DoubleList(size);
		copy.addAll(this);
		return copy;
	}
	
	/** Removes all elements by resetting size to 0 without deallocating array */
	public void clear(){size=0;}
	
	/**
	 * Sets value at specified index, expanding array if necessary.
	 * Updates list size to accommodate the index if it's beyond current size.
	 * @param loc Index position to set
	 * @param value Value to store
	 */
	public final void set(int loc, double value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	/** Sets the last element to specified value.
	 * @param value New value for last element */
	public final void setLast(double value){
		assert(size>0);
		array[size-1]=value;
	}
	
	/**
	 * Adds specified value to element at given index, expanding array if necessary.
	 * Updates list size to accommodate the index if it's beyond current size.
	 * @param loc Index position to increment
	 * @param value Amount to add
	 */
	public final void increment(int loc, double value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	/** Adds 1 to element at given index, expanding array if necessary.
	 * @param loc Index position to increment by 1 */
	public final void increment(int loc){
		increment(loc, 1);
	}
	
	/**
	 * Element-wise addition of another DoubleList to this list.
	 * For each index i, adds b[i] to this[i].
	 * @param b DoubleList to add element-wise
	 */
	public final void incrementBy(DoubleList b){
		for(int i=b.size-1; i>=0; i--){
			increment(i, b.get(i));
		}
	}
	
	/**
	 * Element-wise addition of double array to this list.
	 * For each index i, adds b[i] to this[i].
	 * @param b Array to add element-wise
	 */
	public final void incrementBy(double[] b){
		for(int i=b.length-1; i>=0; i--){
			increment(i, b[i]);
		}
	}
	
	/** Appends all elements from another DoubleList to end of this list.
	 * @param b DoubleList whose elements to append */
	public final void append(DoubleList b){
		for(int i=0; i<b.size; i++){
			add(b.get(i));
		}
	}
	
	/** Appends all elements from double array to end of this list.
	 * @param b Array whose elements to append */
	public final void append(double[] b){
		for(int i=0; i<b.length; i++){
			add(b[i]);
		}
	}
	
	/**
	 * Subtracts each element from specified value, replacing elements with results.
	 * For each element e, sets e = value - e.
	 * @param value Value to subtract elements from
	 */
	public void subtractFrom(double value){
		for(int i=0; i<size; i++){
			array[i]=value-array[i];
		}
	}
	
	/**
	 * Gets value at specified index.
	 * Returns 0 if index is beyond list size.
	 * @param loc Index to retrieve
	 * @return Value at index, or 0 if index >= size
	 */
	public final double get(int loc){
		return(loc>=size ? 0 : array[loc]);//TODO: Shouldn't this crash instead of returning 0?
	}
	
	/** Gets the last element in the list.
	 * @return Last element value */
	public double lastElement() {
		assert(size>0);
		return array[size-1];
	}
	
	/** Appends value to end of list, expanding array if necessary.
	 * @param x Value to add */
	public final void add(double x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	//Slow; for validation
	/**
	 * Checks if list contains any duplicate values using brute force comparison.
	 * Time complexity O(n²) - intended for validation purposes only.
	 * @return true if duplicates found, false otherwise
	 */
	public boolean containsDuplicates(){
		for(int i=0; i<size; i++){
			for(int j=i+1; j<size; j++){
				if(array[i]==array[j]){return true;}
			}
		}
		return false;
	}
	
	/** Appends all elements from another DoubleList to this list.
	 * @param counts DoubleList whose elements to add */
	public void addAll(DoubleList counts) {
		final double[] array2=counts.array;
		final int size2=counts.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}
	
	/**
	 * Checks if list contains specified value using linear search.
	 * @param x Value to search for
	 * @return true if value found, false otherwise
	 */
	public boolean contains(double x) {
		for(int i=0; i<size; i++){
			if(array[i]==x){return true;}
		}
		return false;
	}
	
	public final void setSize(final int size2) {
		if(size2<array.length){resize(size2);}
		size=size2;
	}
	
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		final int size3=(int)Tools.min(Shared.MAX_ARRAY_LEN, size2);
		assert(size3>size) : "Overflow: "+size+", "+size2+" -> "+size3;
		array=KillSwitch.copyOf(array, size3);
	}
	
	/**
	 * Calculates population standard deviation of all elements.
	 * Returns 0 for lists with fewer than 2 elements.
	 * @return Population standard deviation
	 */
	public final double stdev(){
		if(size<2){return 0;}
		double sum=sum();
		double avg=sum/size;
		double sumdev2=0;
		for(int i=0; i<size; i++){
			double x=array[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return Math.sqrt(sumdev2/size);
	}
	
	/** Calculates sum of all elements in the list.
	 * @return Sum of all elements */
	public final double sumLong(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	/** Calculates sum of all elements in the list.
	 * @return Sum of all elements */
	public final double sum(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	/**
	 * Calculates arithmetic mean of all elements.
	 * Returns 0 for empty lists.
	 * @return Mean value of elements, or 0 if empty
	 */
	public final double mean(){
		return size<1 ? 0 : sum()/size;
	}
	
	/** Assumes list is sorted */
	public final double median(){
		if(size<1){return 0;}
		int idx=percentileIndex(0.5);
		return array[idx];
	}
	
	/** Assumes list is sorted */
	public final double mode(){
		if(size<1){return 0;}
		assert(sorted());
		int streak=1, bestStreak=0;
		double prev=array[0];
		double best=prev;
		for(int i=0; i<size; i++){
			double x=array[i];
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
	 * Calculates percentile value for given fraction using sum-weighted approach.
	 * @param fraction Percentile fraction (0.0 to 1.0)
	 * @return Value at specified percentile
	 */
	public double percentile(double fraction){
		if(size<1){return 0;}
		int idx=percentileIndex(fraction);
		return array[idx];
	}
	
	/**
	 * Finds index of element at specified percentile using sum-weighted approach.
	 * Assumes list is sorted and uses cumulative sum to determine percentile position.
	 * @param fraction Percentile fraction (0.0 to 1.0)
	 * @return Index of element at percentile position
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
	
	/**
	 * Reduces array capacity to exactly match current size, freeing unused memory
	 */
	public final void shrink(){
		if(size==array.length){return;}
		array=KillSwitch.copyOf(array, size);
	}
	

	
	/** Removes duplicate values and shrinks array to minimal size */
	public final void shrinkToUnique(){
		condense();
		shrink();
	}
	
	//In-place.
	//Assumes sorted.
	/** Removes duplicate consecutive values in-place, assuming list is sorted.
	 * Preserves one copy of each unique value while reducing list size. */
	public final void condense(){
		if(size<=1){return;}
		
		int i=0, j=1;
		for(; j<size && array[i]<array[j]; i++, j++){}//skip while strictly ascending 
		
		int dupes=0;
		for(; j<size; j++){//This only enters at the first nonascending pair
			double a=array[i], b=array[j];
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
	
	/** Creates copy of list elements as double array.
	 * @return New array containing list elements */
	public double[] toArray(){
		return KillSwitch.copyOf(array, size);
	}
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	/**
	 * Returns string showing only non-zero elements with their indices.
	 * Format: [(index, value), ...]
	 * @return Set view string representation
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
	 * Format: [value1, value2, ...]
	 * @return List view string representation
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
	
	/** Assumes this is sorted.
	 * Reduces the list to a set of unique values;
	 * stores their counts in a second list. */
	public void getUniqueCounts(DoubleList counts) {
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
	
	/** Sorts elements in ascending order using Arrays.sort */
	public void sort() {
		if(size>1){Shared.sort(array, 0, size);}
	}
	
	/** Reverses order of all elements in-place */
	public void reverse() {
		if(size>1){Tools.reverseInPlace(array, 0, size);}
	}
	
	/** Checks if list is sorted in ascending order.
	 * @return true if sorted, false otherwise */
	public boolean sorted(){
		for(int i=1; i<size; i++){
			if(array[i]<array[i-1]){return false;}
		}
		return true;
	}
	
	/** Gets current number of elements in list.
	 * @return Current size */
	public int size() {
		return size;
	}
	
	/** Checks if list contains no elements.
	 * @return true if empty, false otherwise */
	public boolean isEmpty() {
		return size<1;
	}
	
	/** Gets current array capacity.
	 * @return Maximum elements that can be stored without resizing */
	public int capacity() {
		return array.length;
	}
	
	/** Gets number of unused array slots.
	 * @return Available capacity before next resize */
	public int freeSpace() {
		return array.length-size;
	}
	
	/**
	 * Returns smaller of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return Minimum value
	 */
	private static final int min(int x, int y){return x<y ? x : y;}
	/**
	 * Returns larger of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return Maximum value
	 */
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/** Internal array storing the double values */
	public double[] array;
	/** Current number of elements in the list */
	public int size=0;
	
}
