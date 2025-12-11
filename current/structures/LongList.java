package structures;

import java.util.Arrays;

import shared.KillSwitch;
import shared.Shared;
import shared.Tools;



/**
 * Resizable array of long values with statistical and mathematical operations.
 * Provides efficient storage and manipulation of long sequences with automatic
 * resizing, sorting, and various statistical computations including means,
 * standard deviation, and percentiles.
 *
 * @author Brian Bushnell
 */
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
	
	/** Creates a new LongList with default initial capacity of 256 elements */
	public LongList(){this(256);}
	
	/** Creates a new LongList with specified initial capacity.
	 * @param initial Initial capacity for the internal array */
	public LongList(int initial){
		array=KillSwitch.allocLong1D(Math.max(1, initial));
	}

	/** Resets the size to 0 without clearing array contents */
	public void clear(){size=0;}
	/** Resets the size to 0 and fills the entire array with zeros */
	public void clearFull(){
		Arrays.fill(array, 0);
		size=0;
	}
	
	/**
	 * Sets the value at the specified location, expanding array if necessary.
	 * Automatically adjusts size if location is beyond current size.
	 * @param loc Index position to set
	 * @param value Long value to store
	 */
	public final void set(int loc, long value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	/** Sets the value of the last element in the list.
	 * @param value New value for the last element */
	public final void setLast(long value){
		assert(size>0);
		array[size-1]=value;
	}
	
	/**
	 * Increments the value at the specified location by the given amount.
	 * Expands array if necessary and adjusts size accordingly.
	 * @param loc Index position to increment
	 * @param value Amount to add to the current value
	 */
	public final void increment(int loc, long value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	/** Increments the value at the specified location by 1.
	 * @param loc Index position to increment */
	public final void increment(int loc){
		increment(loc, 1);
	}
	
	/**
	 * Increments corresponding positions by values from another LongList.
	 * Processes from end to beginning for efficiency.
	 * @param b LongList containing values to add to this list
	 */
	public final void incrementBy(LongList b){
		for(int i=b.size-1; i>=0; i--){
			increment(i, b.get(i));
		}
	}
	
	/**
	 * Increments corresponding positions by values from a long array.
	 * Processes from end to beginning for efficiency.
	 * @param b Array containing values to add to this list
	 */
	public final void incrementBy(long[] b){
		for(int i=b.length-1; i>=0; i--){
			increment(i, b[i]);
		}
	}
	
	/** Appends all elements from another LongList to the end of this list.
	 * @param b LongList whose elements will be appended */
	public final void append(LongList b){
		for(int i=0; i<b.size; i++){
			add(b.get(i));
		}
	}
	
	/** Appends all elements from a long array to the end of this list.
	 * @param b Array whose elements will be appended */
	public final void append(long[] b){
		for(int i=0; i<b.length; i++){
			add(b[i]);
		}
	}
	
	/**
	 * Gets the value at the specified location.
	 * Returns 0 if the location is beyond the current size.
	 * @param loc Index position to retrieve
	 * @return Value at the specified location, or 0 if out of bounds
	 */
	public final long get(int loc){
		return(loc>=size ? 0 : array[loc]);
	}
	
	/** Adds a value to the end of the list, expanding array if necessary.
	 * @param x Value to add */
	public final void add(long x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	/**
	 * Checks if the list contains the specified value.
	 * @param x Value to search for
	 * @return true if the value is found, false otherwise
	 */
	public boolean contains(long x) {
		for(int i=0; i<size; i++) {
			if(array[i]==x) {return true;}
		}
		return false;
	}
	
	/** Adds all elements from another LongList to this list.
	 * @param list2 LongList whose elements will be added */
	public void addAll(LongList list2) {
		final long[] array2=list2.array;
		final int size2=list2.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}
	
	/**
	 * Resizes the internal array to accommodate more elements.
	 * Ensures the new size is larger than current size and within array limits.
	 * @param size2 Target size for the array
	 */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		final int size3=(int)Tools.min(Shared.MAX_ARRAY_LEN, size2);
		assert(size3>size) : "Overflow: "+size+", "+size2+" -> "+size3;
		array=KillSwitch.copyOf(array, size3);
	}
	
	/**
	 * Shrinks the internal array to match the current size exactly.
	 * Reduces memory usage by removing unused capacity.
	 * @return This LongList for method chaining
	 */
	public final LongList shrink(){
		if(size==array.length){return this;}
		array=KillSwitch.copyOf(array, size);
		return this;
	}
	
	/**
	 * Calculates the standard deviation of all values in the list.
	 * Returns 0 if the list has fewer than 2 elements.
	 * @return Standard deviation of the values
	 */
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
	
	/**
	 * Calculates the average absolute difference between a target value and all list values.
	 * @param x Target value to compare against
	 * @return Average absolute difference
	 */
	public final double avgDif(final double x){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=Tools.absdif(x, array[i]);
		}
		return sum/(Tools.max(1, size));
	}
	
	/**
	 * Calculates the root mean square difference between a target value and all list values.
	 * @param x Target value to compare against
	 * @return Root mean square difference
	 */
	public final double rmsDif(final double x){
		double sum=0;
		for(int i=0; i<size; i++){
			double dif=Tools.absdif(x, array[i]);
			sum+=dif*dif;
		}
		return Math.sqrt(sum/(Tools.max(1, size)));
	}
	
	/** Calculates the sum of all values as a long.
	 * @return Sum of all values in the list */
	public final long sumLong(){
		long sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	/**
	 * Calculates the weighted sum treating this as a histogram.
	 * Each value is multiplied by its index position before summing.
	 * @return Weighted sum with index weights
	 */
	public final long sumHist(){
		long sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i]*i;
		}
		return sum;
	}
	
	/** Calculates the sum of all values as a double.
	 * @return Sum of all values in the list */
	public final double sum(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	/** Calculates the arithmetic mean of all values.
	 * @return Mean value, or 0 if the list is empty */
	public final double mean(){
		return size<1 ? 0 : sum()/size;
	}
	
	/**
	 * Calculates the weighted mean treating this as a histogram.
	 * Uses index-weighted values divided by total sum.
	 * @return Weighted mean, or 0 if the list is empty
	 */
	public final double meanHist(){
		return size<1 ? 0 : sumHist()/sum();
	}
	
	//Ignores elements below 1
	/**
	 * Calculates the harmonic mean of positive values in the list.
	 * Ignores elements with values less than 1.
	 * @return Harmonic mean of positive values
	 */
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
	
	//Ignores elements below 1
	/**
	 * Calculates the geometric mean of positive values in the list.
	 * Ignores elements with values less than 1.
	 * @return Geometric mean of positive values
	 */
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
	
	/** Assumes list is sorted */
	public final double medianWeightedAverage(){
		if(size<1){return 0;}
		int half=size/2;
		long count=0;
		double sum=0;
		for(int i=0, j=size-1; i<half; i++, j--){
			int mult=i+1;
			double incr=(array[i]+array[j])*mult;
			sum+=incr;
			count+=2*mult;
		}
		if((size&1)==1){//odd length
			int mult=half+1;
			double incr=(array[half])*mult;
			sum+=incr;
			count+=2*mult;
		}
		return sum/count;
	}
	
	/** Assumes list is sorted */
	public final long median(){
		if(size<1){return 0;}
		int idx=percentileIndex(0.5);
		return array[idx];
	}
	
	/** Allows unsorted list */
	public final long min(){
		if(size<1){return 0;}
		long x=array[0];
		for(int i=1; i<size; i++){
			x=Tools.min(x, array[i]);
		}
		return x;
	}
	
	/** Allows unsorted list */
	public final long max(){
		if(size<1){return 0;}
		long x=array[0];
		for(int i=1; i<size; i++){
			x=Tools.max(x, array[i]);
		}
		return x;
	}
	
	/** Assumes list is sorted */
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
	 * Returns the value at the specified percentile fraction.
	 * @param fraction Percentile as a fraction (0.0 to 1.0)
	 * @return Value at the specified percentile
	 */
	public long percentile(double fraction){
		if(size<1){return 0;}
		int idx=percentileIndex(fraction);
		return array[idx];
	}
	
	/**
	 * Returns the index position at the specified percentile fraction.
	 * Uses cumulative sum approach to find the target position.
	 * Assumes the list is sorted.
	 *
	 * @param fraction Percentile as a fraction (0.0 to 1.0)
	 * @return Index position at the specified percentile
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
	
	/** Removes duplicate values and shrinks the list to contain only unique elements.
	 * Calls condense() followed by shrink() for efficiency. */
	public final void shrinkToUnique(){
		condense();
		shrink();
	}
	
	//In-place.
	//Assumes sorted.
	/**
	 * Removes duplicate consecutive values in-place from a sorted list.
	 * Assumes the list is already sorted and compacts it by removing duplicates.
	 * More efficient than creating a new array.
	 */
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
	
	/** Removes the element at the specified index by shifting subsequent elements.
	 * @param i Index of element to remove */
	public void removeElementAt(int i) {
		for(int j=i+1; j<size; i++, j++) {
			array[i]=array[j];
		}
		size--;
	}
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	/**
	 * Returns a string representation showing only non-zero elements with their indices.
	 * Format: [(index, value), (index, value), ...]
	 * @return String representation in set view format
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
	 * Returns a string representation showing all elements in list format.
	 * Format: [value1, value2, value3, ...]
	 * @return String representation in list view format
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
	
	/** Creates a new long array containing all elements from this list.
	 * @return New array with copies of all list elements */
	public long[] toArray(){
		long[] x=KillSwitch.allocLong1D(size);
		for(int i=0; i<x.length; i++){
			x[i]=array[i];
		}
		return x;
	}
	
	/**
	 * Sorts the list elements in ascending order using parallel sort if available
	 */
	public void sort() {
		if(size>1){Shared.sort(array, 0, size);}
	}
	
	/**
	 * Sorts the list elements in ascending order using single-threaded Arrays.sort
	 */
	public void sortSerial() {
		if(size>1){Arrays.sort(array, 0, size);}
	}
	
	/** Reverses the order of elements in the list */
	public void reverse() {
		if(size>1){Tools.reverseInPlace(array, 0, size);}
	}
	
	/** Checks if the list is sorted in ascending order.
	 * @return true if sorted, false otherwise */
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
	
	/**
	 * Finds the first index containing the specified value.
	 * @param x Value to search for
	 * @return Index of first occurrence, or -1 if not found
	 */
	public int findIndex(long x) {
		for(int i=0; i<size; i++) {
			if(array[i]==x) {return i;}
		}
		return -1;
	}
	
	/** Assumes sorted ascending */
	public int findIndexAfter(long x) {
		for(int i=0; i<size; i++) {
			if(array[i]>x) {return i;}
		}
		return size;
	}
	
	/**
	 * Caps the histogram by consolidating all values beyond the specified index.
	 * Sums all values from index max to the end and places the sum at index max.
	 * Useful for limiting histogram size while preserving total counts.
	 * @param max Maximum index to retain; values beyond this are summed into this position
	 */
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
	
	/**
	 * Returns the smaller of two long values.
	 * @param x First value
	 * @param y Second value
	 * @return Minimum of x and y
	 */
	private static final long min(long x, long y){return x<y ? x : y;}
	/**
	 * Returns the larger of two long values.
	 * @param x First value
	 * @param y Second value
	 * @return Maximum of x and y
	 */
	private static final long max(long x, long y){return x>y ? x : y;}
	
	/**
	 * Returns the smaller of two int values.
	 * @param x First value
	 * @param y Second value
	 * @return Minimum of x and y
	 */
	private static final int min(int x, int y){return x<y ? x : y;}
	/**
	 * Returns the larger of two int values.
	 * @param x First value
	 * @param y Second value
	 * @return Maximum of x and y
	 */
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/** Internal array storing the long values */
	public long[] array;
	/** Highest occupied index plus 1, i.e., lowest unoccupied index */
	public int size=0;
	
}
