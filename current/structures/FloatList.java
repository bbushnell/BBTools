package structures;

import java.util.ArrayList;
import java.util.LinkedList;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;


/**
 * Growable list of primitive float values backed by a public array.
 * The logical size limits normal reads and aggregate operations; clear() resets
 * only the logical size and leaves the backing array available for reuse.
 */
public final class FloatList{

	/*--------------------------------------------------------------*/
	/*----------------            Main            ----------------*/
	/*--------------------------------------------------------------*/

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

	/*--------------------------------------------------------------*/
	/*----------------       Initialization       ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs a copy of x's live prefix. */
	public FloatList(FloatList x){
		this(x.size);
		addAll(x);
	}
	
	public FloatList(){this(256);}
	
	public FloatList(int initial){
//		assert(initial>0) : initial+"\n"+this;
		initial=Tools.max(initial, 1);
		array=KillSwitch.allocFloat1D(initial);
	}
	
	/** Returns a new list containing this list's live prefix. */
	public FloatList copy(){
		FloatList copy=new FloatList(size);
		copy.addAll(this);
		return copy;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Mutation          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Clears the logical contents without clearing the backing array.
	 * Subsequent add() calls overwrite the live prefix; increment() on an unwritten
	 * index can combine with a stale backing value from before clear().
	 */
	//TODO: increment after logical clear can reuse stale values from old elements.
	public void clear(){size=0;}
	
	/**
	 * Stores value at loc, expanding the backing array when loc reaches or exceeds current capacity.
	 * If loc is beyond the current logical end, the logical size grows to loc+1.
	 */
	public final void set(int loc, float value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	/** Replaces the last live value; requires at least one element. */
	public final void setLast(float value){
		assert(size>0);
		array[size-1]=value;
	}
	
	/**
	 * Adds value at loc, growing the logical size to include loc when needed.
	 */
	public final void increment(int loc, float value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	/** Increments loc by one, expanding the logical size to include loc when needed. */
	public final void increment(int loc){
		increment(loc, 1);
	}
	
	/** Adds b's live values into matching indices of this list. */
	public final void incrementBy(FloatList b){
		for(int i=b.size-1; i>=0; i--){
			increment(i, b.get(i));
		}
	}
	
	/** Adds b's full array values into matching indices of this list. */
	public final void incrementBy(float[] b){
		for(int i=b.length-1; i>=0; i--){
			increment(i, b[i]);
		}
	}
	
	/**
	 * Appends all logical values from b.
	 * With assertions enabled, self-append is rejected because nonempty self-append
	 * would grow the loop bound while it is being read.
	 */
	public final void append(FloatList b){
		assert(b!=this);
		for(int i=0; i<b.size; i++){
			add(b.get(i));
		}
	}
	
	/** Appends every value from b. */
	public final void append(float[] b){
		for(int i=0; i<b.length; i++){
			add(b[i]);
		}
	}
	
	/** Replaces each live value x with value-x. */
	public void subtractFrom(float value){
		for(int i=0; i<size; i++){
			array[i]=value-array[i];
		}
	}
	
	/**
	 * Returns the logical value at loc, or 0 for nonnegative locs outside the logical size.
	 * Negative locs are rejected by the backing array rather than normalized here.
	 */
	public final float get(int loc){
		return(loc>=size ? 0 : array[loc]);//TODO: Shouldn't this crash instead of returning 0?
	}
	
	/** Returns the last live value; requires at least one element. */
	public float lastElement(){
		assert(size>0);
		return array[size-1];
	}
	
	/** Appends x to the live prefix, expanding capacity when needed. */
	public final void add(float x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	/** Returns true if any two live values compare equal. Slow; for validation. */
	public boolean containsDuplicates(){
		for(int i=0; i<size; i++){
			for(int j=i+1; j<size; j++){
				if(array[i]==array[j]){return true;}
			}
		}
		return false;
	}
	
	/** Appends all live values from counts. */
	public void addAll(FloatList counts){
		final float[] array2=counts.array;
		final int size2=counts.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}
	
	/**
	 * Tests the live prefix for an exact {@code ==} match.
	 * This uses primitive equality: signed zeroes compare equal and NaN does not
	 * match any stored value, including another NaN.
	 */
	public boolean contains(float x){
		for(int i=0; i<size; i++){
			if(array[i]==x){return true;} //Possible bug: Float equality comparison may fail due to precision
		}
		return false;
	}
	
	/**
	 * Sets the logical size to a supported nonnegative size within the array-length cap.
	 * Existing backing slots keep their values; newly allocated slots are zero-filled.
	 */
	public final void setSize(final int size2){
		if(size2>array.length){resize(size2);}
		size=size2;
	}
	
	/** Resizes the backing array to min(size2, Shared.MAX_ARRAY_LEN). */
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		final int size3=(int)Tools.min(Shared.MAX_ARRAY_LEN, size2);
		assert(size3>size) : "Overflow: "+size+", "+size2+" -> "+size3;
		array=KillSwitch.copyOf(array, size3);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Statistics         ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns -1 when empty; otherwise starts at index 0 and replaces the index on strict increases until a NaN is encountered. */
	public int maxIdx(){
		if(size<1){return -1;}
		float max=array[0];
		int maxIdx=0;
		for(int i=1; i<size; i++){
			if(array[i]>max){
				max=array[i];
				maxIdx=i;
			}
			max=max(max, array[i]);
		}
		return maxIdx;
	}
	
	/** Folds Math.max over live values starting from -Float.MAX_VALUE. */
	public float max(){
		float max=-Float.MAX_VALUE;
		for(int i=0; i<size; i++){max=max(max, array[i]);}
		return max;
	}
	
	/** Folds Math.min over live values starting from Float.MAX_VALUE. */
	public float min(){
		float min=Float.MAX_VALUE;
		for(int i=0; i<size; i++){min=min(min, array[i]);}
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
	
	/** Returns the sum of all live values as a double accumulator. */
	public final double sumLong(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	/** Returns the sum of all live values. */
	public final double sum(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	/** Returns the arithmetic mean, or 0 for an empty list. */
	public final double mean(){
		return size<1 ? 0 : sum()/size;
	}
	
	/**
	 * Calculates median value assuming list is sorted.
	 * Uses percentileIndex(0.5) to select the median element.
	 * @return Median value, or 0 if list is empty
	 */
	public final double median(){
		if(size<1){return 0;}
		int idx=percentileIndex(0.5f);
		return array[idx];
	}
	
	/**
	 * Finds the most frequently occurring value assuming list is sorted.
	 * Returns the value with the longest consecutive sequence.
	 * @return Most frequent value, or 0 if list is empty
	 */
	public final float mode(){
		if(size<1){return 0;}
		assert(sorted());
		int streak=1, bestStreak=0;
		float prev=array[0];
		float best=prev;
		for(int i=1; i<size; i++){
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
	 * Returns the sorted value at the sum-weighted percentile index for fraction.
	 * This is not a positional percentile; percentileIndex accumulates values until
	 * the running sum reaches {@code sum()*fraction}. Fraction is not clamped; an
	 * empty list returns 0.
	 */
	public float percentile(double fraction){
		if(size<1){return 0;}
		int idx=percentileIndex(fraction);
		return array[idx];
	}
	
	/**
	 * Returns the first sorted index whose cumulative sum reaches
	 * {@code sum()*fraction}, or the last live index when the target is not reached.
	 * An empty list returns -1 and a singleton returns 0.
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

	/*--------------------------------------------------------------*/
	/*----------------       Deduplication        ----------------*/
	/*--------------------------------------------------------------*/

	/** Shrinks backing capacity to the logical size and returns this list. */
	public final FloatList shrink(){
		if(size==array.length){return this;}
		array=KillSwitch.copyOf(array, size);
		return this;
	}
	

	
	/** Condenses sorted, NaN-free duplicates in place, then shrinks capacity to the new size. */
	public final void shrinkToUnique(){
		condense();
		shrink();
	}
	
	/** Condenses equal runs in a sorted, NaN-free live prefix, preserving one copy of each value. */
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
	
	/** Returns a compact copy of the live prefix. */
	public float[] toArray(){
		return KillSwitch.copyOf(array, size);
	}

	/*--------------------------------------------------------------*/
	/*----------------          ToString          ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns string representation in list format showing all elements.
	 * @return String in format [element1, element2, ...] */
	@Override
	public String toString(){
		return toStringListView();
	}
	
	/** Returns a sparse index/value view containing only nonzero live entries. */
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
	
	/** Returns a dense list view of all live entries. */
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
	/*----------------          Counting          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reduces this sorted list to distinct values and writes their frequencies.
	 * This list receives one copy of each distinct input value, and counts is
	 * cleared and filled with the matching occurrence count for each distinct value.
	 * counts must be a different list from this input list.
	 * @param counts Output list to replace with occurrence counts for the distinct values
	 */
	public void getUniqueCounts(FloatList counts){
		assert(counts!=this) : "Aliased counts would erase input before it is read.";
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

	/*--------------------------------------------------------------*/
	/*----------------         Utilities          ----------------*/
	/*--------------------------------------------------------------*/

	/** Sorts the live prefix in ascending order and returns this list. */
	public FloatList sort(){
		if(size>1){Shared.sort(array, 0, size);}
		return this;
	}
	
	/** Reverses the live prefix in place and returns this list. */
	public FloatList reverse(){
		if(size>1){Tools.reverseInPlace(array, 0, size);}
		return this;
	}
	
	/** Returns true when no adjacent pair is descending; pairs containing NaN do not fail this test. */
	public boolean sorted(){
		for(int i=1; i<size; i++){
			if(array[i]<array[i-1]){return false;}
		}
		return true;
	}
	
	/** Returns the logical element count. */
	public int size(){
		return size;
	}
	
	/** Returns true when this list has no live elements. */
	public boolean isEmpty(){
		return size<1;
	}
	
	/** Returns the current backing-array length. */
	public int capacity(){
		return array.length;
	}
	
	/** Returns capacity minus logical size. */
	public int freeSpace(){
		return array.length-size;
	}
	
	private static final int min(int x, int y){return Math.min(x, y);}
	private static final int max(int x, int y){return Math.max(x, y);}
	private static final float min(float x, float y){return Math.min(x, y);}
	private static final float max(float x, float y){return Math.max(x, y);}

	/*--------------------------------------------------------------*/
	/*----------------           Fields           ----------------*/
	/*--------------------------------------------------------------*/

	/** Backing storage; only indices below size are live. */
	public float[] array;
	/** Number of live values in array. */
	public int size=0;
	
}
