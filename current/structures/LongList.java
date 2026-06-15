package structures;

import java.util.Arrays;

import shared.KillSwitch;
import shared.Shared;
import shared.Tools;



public final class LongList{
	
	public static void main(String[] args){
		LongList list=new LongList();
		list.add(3);
		list.add(1);
		list.add(2);
		list.add(5);
		list.add(2);
		list.add(1);
		list.add(7);
		list.add(3);
		list.add(3);
		System.err.println(list);
		list.sort();
		System.err.println(list);
		list.condense();
		System.err.println(list);
		list.condense();
		System.err.println(list);
	}
	
	public LongList(){this(256);}

	/**
	 * Creates LongList with specified initial capacity.
	 * Minimum capacity is enforced to be at least 1.
	 * @param initial Initial array size (will be at least 1)
	 */
	public LongList(int initial){
		array=KillSwitch.allocLong1D(Math.max(1, initial));
	}

	/** Clears the list by setting size to 0; array contents are left intact */
	public void clear(){size=0;}
	/** Clears list and zeros all array elements for security.
	 * More thorough than clear() but slower. */
	public void clearFull(){
		Arrays.fill(array, 0);
		size=0;
	}
	
	/**
	 * Sets value at specified index, expanding array if necessary.
	 * Updates list size to include the set position.
	 * @param loc Index to set (will expand array if needed)
	 * @param value Value to store at location
	 */
	public final void set(int loc, long value){
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
	public final void setLast(long value){
		assert(size>0);
		array[size-1]=value;
	}

	/**
	 * Increments value at location by specified amount.
	 * Expands array and updates size as needed.
	 * @param loc Index to increment
	 * @param value Amount to add to existing value
	 */
	public final void increment(int loc, long value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}

	/** Increments value at location by 1, expanding array if necessary.
	 * @param loc Index to increment */
	public final void increment(int loc){
		increment(loc, 1);
	}

	/** Adds each element of b to this list at the same index, expanding as needed.
	 * @param b Source list whose values are added position-by-position */
	public final void incrementBy(LongList b){
		for(int i=b.size-1; i>=0; i--){
			increment(i, b.get(i));
		}
	}

	/** Adds each element of b to this list at the same index, expanding as needed.
	 * @param b Source array whose values are added position-by-position */
	public final void incrementBy(long[] b){
		for(int i=b.length-1; i>=0; i--){
			increment(i, b[i]);
		}
	}

	/** Appends all elements of b to the end of this list in order.
	 * @param b Source list to copy elements from */
	public final void append(LongList b){
		for(int i=0; i<b.size; i++){
			add(b.get(i));
		}
	}

	/** Appends all elements of b to the end of this list in order.
	 * @param b Source array to copy elements from */
	public final void append(long[] b){
		for(int i=0; i<b.length; i++){
			add(b[i]);
		}
	}

	/**
	 * Gets value at specified index with bounds checking.
	 * Returns 0 for out-of-bounds access instead of throwing exception.
	 * @param loc Index to retrieve
	 * @return Value at index, or 0 if index >= size
	 */
	public final long get(int loc){
		return(loc>=size ? 0 : array[loc]);
	}

	/**
	 * Appends element to end of list, expanding capacity if needed.
	 * Doubles array size when expansion is required.
	 * @param x Value to add to the list
	 */
	public final void add(long x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}

	/**
	 * Checks if list contains specified value using linear search.
	 * @param x Value to search for
	 * @return true if value is found, false otherwise
	 */
	public boolean contains(long x) {
		for(int i=0; i<size; i++) {
			if(array[i]==x) {return true;}
		}
		return false;
	}

	/**
	 * Appends all elements from another LongList to this list.
	 * Elements are added in order from the source list.
	 * @param list2 Source LongList to copy elements from
	 */
	public void addAll(LongList list2) {
		final long[] array2=list2.array;
		final int size2=list2.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}

	/**
	 * Appends a range of elements from a source array via System.arraycopy.
	 * @param array2 Source array to copy from
	 * @param from Inclusive start index in source
	 * @param to Exclusive end index in source
	 */
	public void add(long[] array2, int from, int to) {
		int len=to-from;
		expand(len);
		System.arraycopy(array2, from, array, size, len);
		size+=len;
	}

	/** Resizes if appending extra elements would exceed capacity.
	 * @param extra Number of additional elements to make room for */
	private final void expand(final long extra) {
		if(size+extra>=array.length){
			resize(size*2L+1);
		}
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
	 * Shrinks internal array to match current size exactly.
	 * Reduces memory usage by eliminating unused capacity.
	 * @return This LongList for method chaining
	 */
	public final LongList shrink(){
		if(size==array.length){return this;}
		array=KillSwitch.copyOf(array, size);
		return this;
	}

	/** Computes population standard deviation of the elements.
	 * @return Standard deviation, or 0 if fewer than 2 elements */
	public final double stdev(){
		if(size<2){return 0;}
		double sum=sum();
		double avg=sum/size;
		double sumdev2=0;
		for(int i=0; i<size; i++){
			long x=array[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return Math.sqrt(sumdev2/size);
	}
	
	/** Mean absolute difference between x and each element.
	 * @param x Reference value
	 * @return Average of |x - element| over all elements (0 if empty) */
	public final double avgDif(final double x){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=Tools.absdif(x, array[i]);
		}
		return sum/(Tools.max(1, size));
	}

	/** Root-mean-square difference between x and each element.
	 * @param x Reference value
	 * @return sqrt(mean of (x - element)^2) over all elements (0 if empty) */
	public final double rmsDif(final double x){
		double sum=0;
		for(int i=0; i<size; i++){
			double dif=Tools.absdif(x, array[i]);
			sum+=dif*dif;
		}
		return Math.sqrt(sum/(Tools.max(1, size)));
	}

	/**
	 * Calculates sum of all elements as long.
	 * @return Sum of all elements as long
	 */
	public final long sumLong(){
		long sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}

	/** Histogram-weighted sum: each element is multiplied by its index.
	 * @return Sum of array[i]*i over all elements */
	public final long sumHist(){
		long sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i]*i;
		}
		return sum;
	}

	/**
	 * Calculates sum of all elements as double.
	 * @return Sum of all elements as double
	 */
	public final double sum(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}

	/** Arithmetic mean of all elements.
	 * @return sum/size, or 0 if empty */
	public final double mean(){
		return size<1 ? 0 : sum()/size;
	}

	/** Histogram mean: index-weighted sum divided by element sum.
	 * Treats the list as a histogram of counts indexed by position.
	 * @return sumHist()/sum(), or 0 if empty */
	public final double meanHist(){
		return size<1 ? 0 : sumHist()/sum();
	}

	/** Harmonic mean of elements, ignoring any element less than 1.
	 * @return Harmonic mean of positive elements */
	//Ignores elements below 1
	public final double harmonicMean(){
		double sum=0;
		int count=0;
		for(int i=0; i<size; i++){
			if(array[i]>0){
				sum+=1.0/array[i];
				count++;
			}
		}
		double avg=sum/Tools.max(1, count);
		return 1.0/avg;
	}
	
	/** Geometric mean of elements, ignoring any element less than 1.
	 * @return Geometric mean of positive elements */
	//Ignores elements below 1
	public final double geometricMean(){
		double sum=0;
		int count=0;
		for(int i=0; i<size; i++){
			if(array[i]>0){
				sum+=Math.log(array[i]);
				count++;
			}
		}
		double avg=sum/Tools.max(1, count);
		return Math.exp(avg);
	}
	
	/**
	 * Calculates a weighted average that emphasizes values near the median position.
	 * Assumes the list is sorted. Uses position-based weighting with higher weights for elements closer to the center.
	 * @return Median-weighted average value
	 */
	public final double medianWeightedAverage(){
		if(size<1){return 0;}
		int half=size/2;
		long count=0;
		double sum=0;
		for(int i=0, j=size-1; i<half; i++, j--){
			int mult=i+1;
			double incr=((double)array[i]+(double)array[j])*mult;
			sum+=incr;
			count+=2*mult;
		}
		if((size&1)==1){//odd length
			int mult=half+1;
			double incr=((double)array[half])*mult;
			sum+=incr;
			count+=mult;
		}
		return sum/count;
	}
	
	/**
	 * Returns the median value from the list.
	 * Assumes the list is sorted.
	 * @return Median value, or 0 if the list is empty
	 */
	public final long median(){
		if(size<1){return 0;}
		int idx=percentileIndex(0.5);
		return array[idx];
	}
	
	/**
	 * Finds the minimum value in the list.
	 * Works with unsorted lists.
	 * @return Minimum value, or 0 if the list is empty
	 */
	public final long min(){
		if(size<1){return 0;}
		long x=array[0];
		for(int i=1; i<size; i++){
			x=Tools.min(x, array[i]);
		}
		return x;
	}
	
	/**
	 * Finds the maximum value in the list.
	 * Works with unsorted lists.
	 * @return Maximum value, or 0 if the list is empty
	 */
	public final long max(){
		if(size<1){return 0;}
		long x=array[0];
		for(int i=1; i<size; i++){
			x=Tools.max(x, array[i]);
		}
		return x;
	}
	
	/**
	 * Finds the most frequently occurring value in the list.
	 * Assumes the list is sorted for efficient consecutive value counting.
	 * @return Most frequent value, or 0 if the list is empty
	 */
	public final long mode(){
		if(size<1){return 0;}
		assert(sorted());
		int streak=1, bestStreak=0;
		long prev=array[0];
		long best=prev;
		for(int i=0; i<size; i++){
			long x=array[i];
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
	 * Finds the value where the cumulative sum reaches the given fraction.
	 * Assumes the list is sorted; uses value-weighted (not position) percentile.
	 * @param fraction Percentile as fraction (0.0 to 1.0)
	 * @return Value at specified percentile, or 0 if empty
	 */
	public long percentile(double fraction){
		if(size<1){return 0;}
		int idx=percentileIndex(fraction);
		return array[idx];
	}

	/**
	 * Finds index where cumulative sum reaches the target percentile.
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
	
//	//TODO: This could be done in-place.
//	public final void shrinkToUnique(){
//		//Assumes sorted.
//		if(size<=0){
//			shrink();
//			return;
//		}
//		
//		int unique=1;
//		
//		for(int i=1; i<size; i++){
//			assert(array[i]>=array[i-1]);
//			if(array[i]!=array[i-1]){unique++;}
//		}
//		if(unique==array.length){return;}
//		long[] alt=KillSwitch.allocLong1D(unique);
//		
//		alt[0]=array[0];
//		for(int i=1, j=1; j<unique; i++){
//			if(array[i]!=array[i-1]){
//				alt[j]=array[i];
//				j++;
//			}
//		}
//		
//		array=alt;
//		size=alt.length;
//	}
	
	/** Removes duplicate elements and shrinks array to fit exactly */
	public final void shrinkToUnique(){
		condense();
		shrink();
	}

	/**
	 * Removes duplicate elements in-place from a sorted list.
	 * Maintains sorted order while keeping only unique values.
	 * Skips the initial strictly-ascending run before collapsing duplicates.
	 */
	//In-place.
	//Assumes sorted.
	public final void condense(){
		if(size<=1){return;}
		
		int i=0, j=1;
		for(; j<size && array[i]<array[j]; i++, j++){}//skip while strictly ascending 
		
		int dupes=0;
		for(; j<size; j++){//This only enters at the first nonascending pair
			long a=array[i], b=array[j];
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
	
	/** Removes the element at index i by shifting later elements left.
	 * @param i Index of the element to remove */
	public void removeElementAt(int i) {
		for(int j=i+1; j<size; i++, j++) {
			array[i]=array[j];
		}
		size--;
	}
	
	/** Returns a string representation of the list in list format.
	 * @return String representation showing all elements */
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

	/**
	 * Creates new array containing a copy of all list elements.
	 * Returned array length equals list size, not capacity.
	 * @return New long array with copies of all elements
	 */
	public long[] toArray(){
		long[] x=KillSwitch.allocLong1D(size);
		for(int i=0; i<x.length; i++){
			x[i]=array[i];
		}
		return x;
	}

	/** Sorts elements in ascending order using Shared.sort (parallel-capable) */
	public void sort() {
		if(size>1){Shared.sort(array, 0, size);}
	}

	/** Sorts elements in ascending order using single-threaded Arrays.sort */
	public void sortSerial() {
		if(size>1){Arrays.sort(array, 0, size);}
	}

	/** Reverses the order of all elements in-place */
	public void reverse() {
		if(size>1){Tools.reverseInPlace(array, 0, size);}
	}

	/** Checks if list is sorted in ascending order using O(n) scan.
	 * @return true if elements are in non-decreasing order, false otherwise */
	public boolean sorted(){
		for(int i=1; i<size; i++){
			if(array[i]<array[i-1]){return false;}
		}
		return true;
	}
	
	public int size() {
		return size;
	}
	
	public int capacity() {
		return array.length;
	}
	
	public int freeSpace() {
		return array.length-size;
	}
	
	public boolean isEmpty() {
		return size==0;
	}
	
	/** Finds the first index containing the specified value via linear search.
	 * @param x Value to search for
	 * @return Index of first match, or -1 if not found */
	public int findIndex(long x) {
		for(int i=0; i<size; i++) {
			if(array[i]==x) {return i;}
		}
		return -1;
	}
	
	/**
	 * Finds the first index containing a value greater than the specified value.
	 * Assumes the list is sorted in ascending order.
	 * @param x Value to compare against
	 * @return Index of first value greater than x, or size if none found
	 */
	public int findIndexAfter(long x) {
		for(int i=0; i<size; i++) {
			if(array[i]>x) {return i;}
		}
		return size;
	}
	
	/** Collapses the tail [max, size) into index max by summing those elements,
	 * then truncates size to max+1. Treats the list as a histogram with a cap.
	 * @param max Index of the final retained bin holding the accumulated tail */
	public void capHist(final int max) {
		if(size<=max+1) {return;}
		//size=2, max=0 are the lowest values to enter
		long sum=0;
		for(int i=size-1; i>=max; i--) {
			sum+=array[i];
		}
		array[max]=sum;
		size=max+1;
	}

	/** Returns smaller of two longs */
	private static final long min(long x, long y){return x<y ? x : y;}
	/** Returns larger of two longs */
	private static final long max(long x, long y){return x>y ? x : y;}

	/** Returns smaller of two integers */
	private static final int min(int x, int y){return x<y ? x : y;}
	/** Returns larger of two integers */
	private static final int max(int x, int y){return x>y ? x : y;}

	/** Backing array for element storage, may have unused capacity */
	public long[] array;
	/** Current number of elements in the list (logical size) */
	public int size=0;
	
}
