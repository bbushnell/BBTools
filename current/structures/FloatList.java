package structures;

import java.util.ArrayList;
import java.util.LinkedList;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;



/**
 * Dynamic array implementation for float values with optimized memory usage.
 * Provides array-like operations without boxing overhead of ArrayList&lt;Float&gt;.
 * Includes statistical functions, sorting, and specialized bioinformatics operations.
 *
 * @author Brian Bushnell
 * @date 2014
 */
public final class FloatList{
	
	/**
	 * Performance comparison benchmark between FloatList, ArrayList, and LinkedList.
	 * Tests memory usage and execution time for adding elements.
	 * @param args Command-line arguments; args[0] specifies number of elements (default 100M)
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		int length=args.length>0 ? Integer.parseInt(args[0]) : 100000000;
		
		System.gc();
		
		{
			System.err.println("\nFloatList:");
			Shared.printMemory();
			t.start();
			FloatList list=new FloatList();
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
			ArrayList<Float> list=new ArrayList<Float>();
			for(int i=0; i<length; i++){
				list.add((float)i);
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
			LinkedList<Float> list=new LinkedList<Float>();
			for(int i=0; i<length; i++){
				list.add((float)i);
			}
			t.stop("Time: \t");
			System.gc();
			Shared.printMemory();
			list=null;
			System.gc();
		}
	}
	
	/** Creates a FloatList with initial capacity of 256 elements */
	public FloatList(){this(256);}
	
	/** Creates a FloatList with specified initial capacity.
	 * @param initial Initial capacity in elements (minimum 1) */
	public FloatList(int initial){
//		assert(initial>0) : initial+"\n"+this;
		initial=Tools.max(initial, 1);
		array=KillSwitch.allocFloat1D(initial);
	}
	
	/** Creates a deep copy of this FloatList.
	 * @return New FloatList containing all elements from this list */
	public FloatList copy() {
		FloatList copy=new FloatList(size);
		copy.addAll(this);
		return copy;
	}
	
	//TODO: Should zero out also, to prevent data retention with later increment ops 
	/**
	 * Removes all elements by setting size to zero; does not clear array contents
	 */
	public void clear(){size=0;}
	
	/**
	 * Sets value at specified index, expanding array if necessary.
	 * Updates size to include this index if beyond current size.
	 * @param loc Index to set
	 * @param value Float value to store
	 */
	public final void set(int loc, float value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	/** Sets the value of the last element in the list.
	 * @param value New value for the last element */
	public final void setLast(float value){
		assert(size>0);
		array[size-1]=value;
	}
	
	/**
	 * Adds value to the element at specified index, expanding array if necessary.
	 * Updates size to include this index if beyond current size.
	 * @param loc Index to increment
	 * @param value Amount to add
	 */
	public final void increment(int loc, float value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	/** Increments element at specified index by 1.
	 * @param loc Index to increment */
	public final void increment(int loc){
		increment(loc, 1);
	}
	
	/**
	 * Adds corresponding elements from another FloatList to this list.
	 * Element at index i in this list is increased by element at index i in list b.
	 * @param b FloatList whose values will be added to this list
	 */
	public final void incrementBy(FloatList b){
		for(int i=b.size-1; i>=0; i--){
			increment(i, b.get(i));
		}
	}
	
	/**
	 * Adds corresponding elements from float array to this list.
	 * Element at index i in this list is increased by element at index i in array b.
	 * @param b Float array whose values will be added to this list
	 */
	public final void incrementBy(float[] b){
		for(int i=b.length-1; i>=0; i--){
			increment(i, b[i]);
		}
	}
	
	/** Appends all elements from another FloatList to the end of this list.
	 * @param b FloatList whose elements will be appended */
	public final void append(FloatList b){
		assert(b!=this);
		for(int i=0; i<b.size; i++){
			add(b.get(i));
		}
	}
	
	/** Appends all elements from float array to the end of this list.
	 * @param b Float array whose elements will be appended */
	public final void append(float[] b){
		for(int i=0; i<b.length; i++){
			add(b[i]);
		}
	}
	
	/**
	 * Subtracts each element from the specified value (value - element).
	 * Modifies all elements in-place.
	 * @param value Value from which each element will be subtracted
	 */
	public void subtractFrom(float value){
		for(int i=0; i<size; i++){
			array[i]=value-array[i];
		}
	}
	
	/**
	 * Gets value at specified index.
	 * Returns 0 for indices beyond current size instead of throwing exception.
	 * @param loc Index to retrieve
	 * @return Value at index, or 0 if index exceeds size
	 */
	public final float get(int loc){
		return(loc>=size ? 0 : array[loc]);//TODO: Shouldn't this crash instead of returning 0?
	}
	
	/** Gets the last element in the list.
	 * @return Value of the last element */
	public float lastElement() {
		assert(size>0);
		return array[size-1];
	}
	
	/** Adds element to the end of the list, expanding array if necessary.
	 * @param x Float value to add */
	public final void add(float x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	//Slow; for validation
	/**
	 * Checks if list contains duplicate values using O(n²) comparison.
	 * Marked as slow and intended for validation purposes only.
	 * @return true if any value appears more than once
	 */
	public boolean containsDuplicates(){
		for(int i=0; i<size; i++){
			for(int j=i+1; j<size; j++){
				if(array[i]==array[j]){return true;}
			}
		}
		return false;
	}
	
	/** Appends all elements from another FloatList to this list.
	 * @param counts FloatList whose elements will be added */
	public void addAll(FloatList counts) {
		final float[] array2=counts.array;
		final int size2=counts.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}
	
	/**
	 * Checks if list contains the specified value using linear search.
	 * @param x Value to search for
	 * @return true if value is found in the list
	 */
	public boolean contains(float x) {
		for(int i=0; i<size; i++){
			if(array[i]==x){return true;} //Possible bug: Float equality comparison may fail due to precision
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
	
	/** Finds the index of the maximum value in the list.
	 * @return Index of maximum value, or -1 if list is empty */
	public int maxIdx() {
		if(size<1) {return -1;}
		float max=array[0];
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
	
	/** Finds the maximum value in the list.
	 * @return Maximum value, or 0 if list is empty */
	public float max() {
		float max=-Float.MAX_VALUE;
		for(int i=0; i<size; i++) {max=max(max, array[i]);}
		return max;
	}
	
	public float min() {
		float min=Float.MAX_VALUE;
		for(int i=1; i<size; i++) {min=min(min, array[i]);}
		return min;
	}
	
	/** Calculates population standard deviation of all elements.
	 * @return Standard deviation, or 0 if fewer than 2 elements */
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
	
	/** Calculates sum of all elements using double precision.
	 * @return Sum of all elements as double */
	public final double sumLong(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	/** Calculates sum of all elements using double precision.
	 * @return Sum of all elements as double */
	public final double sum(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	/** Calculates arithmetic mean of all elements.
	 * @return Average value, or 0 if list is empty */
	public final double mean(){
		return size<1 ? 0 : sum()/size;
	}
	
	/** Assumes list is sorted */
	public final double median(){
		if(size<1){return 0;}
		int idx=percentileIndex(0.5f);
		return array[idx];
	}
	
	/** Assumes list is sorted */
	public final float mode(){
		if(size<1){return 0;}
		assert(sorted());
		int streak=1, bestStreak=0;
		float prev=array[0];
		float best=prev;
		for(int i=0; i<size; i++){
			float x=array[i];
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
	 * Calculates percentile value using weighted cumulative sum method.
	 * @param fraction Percentile as fraction (0.0 to 1.0)
	 * @return Value at the specified percentile
	 */
	public float percentile(double fraction){
		if(size<1){return 0;}
		int idx=percentileIndex(fraction);
		return array[idx];
	}
	
	/**
	 * Finds index corresponding to percentile using weighted cumulative sum.
	 * Assumes list is sorted and uses sum-based weighting rather than count-based.
	 * @param fraction Percentile as fraction (0.0 to 1.0)
	 * @return Index at the specified percentile
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
	
	/** Reduces array capacity to match current size, freeing unused memory.
	 * @return This FloatList for method chaining */
	public final FloatList shrink(){
		if(size==array.length){return this;}
		array=KillSwitch.copyOf(array, size);
		return this;
	}
	

	
	/** Removes duplicates and shrinks array to minimal size in one operation */
	public final void shrinkToUnique(){
		condense();
		shrink();
	}
	
	//In-place.
	//Assumes sorted.
	/** Removes duplicate elements in-place, assuming list is sorted.
	 * Maintains ascending order and reduces size to unique element count. */
	public final void condense(){
		if(size<=1){return;}
		
		int i=0, j=1;
		for(; j<size && array[i]<array[j]; i++, j++){}//skip while strictly ascending 
		
		int dupes=0;
		for(; j<size; j++){//This only enters at the first nonascending pair
			float a=array[i], b=array[j];
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
	
	/** Creates copy of elements as standard float array.
	 * @return New float array containing all elements */
	public float[] toArray(){
		return KillSwitch.copyOf(array, size);
	}
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	/**
	 * Returns string representation showing only non-zero elements with indices.
	 * Useful for sparse data visualization.
	 * @return String in format [(index1, value1), (index2, value2), ...]
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
	
	/** Returns string representation showing all elements in list format.
	 * @return String in format [element1, element2, ...] */
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
	public void getUniqueCounts(FloatList counts) {
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
	
	public FloatList sort() {
		if(size>1){Shared.sort(array, 0, size);}
		return this;
	}
	
	public FloatList reverse() {
		if(size>1){Tools.reverseInPlace(array, 0, size);}
		return this;
	}
	
	/** Checks if list is sorted in ascending order.
	 * @return true if elements are in non-decreasing order */
	public boolean sorted(){
		for(int i=1; i<size; i++){
			if(array[i]<array[i-1]){return false;}
		}
		return true;
	}
	
	/** Gets the number of elements in the list.
	 * @return Current size */
	public int size() {
		return size;
	}
	
	/** Checks if list contains no elements.
	 * @return true if size is zero */
	public boolean isEmpty() {
		return size<1;
	}
	
	/** Gets the current array capacity.
	 * @return Maximum elements that can be stored without resizing */
	public int capacity() {
		return array.length;
	}
	
	/** Gets unused capacity in the internal array.
	 * @return Number of elements that can be added without resizing */
	public int freeSpace() {
		return array.length-size;
	}
	
	/**
	 * Returns minimum of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return Smaller value
	 */
	private static final int min(int x, int y){return Math.min(x, y);}
	/**
	 * Returns maximum of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return Larger value
	 */
	private static final int max(int x, int y){return Math.max(x, y);}
	private static final float min(float x, float y){return Math.min(x, y);}
	/**
	 * Returns maximum of two floats.
	 * @param x First float
	 * @param y Second float
	 * @return Larger value
	 */
	private static final float max(float x, float y){return Math.max(x, y);}
	
	/** Internal storage array for float elements */
	public float[] array;
	/** Number of elements currently stored in the list */
	public int size=0;
	
}
