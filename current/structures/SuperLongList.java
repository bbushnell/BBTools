package structures;

import shared.Tools;

/**
 * Hybrid histogram for long-valued distributions.
 *
 * <p>Values in the half-open range {@code [0, limit)} are counted by value in
 * {@link #array}; values greater than or equal to {@code limit} are kept as
 * individual entries in {@link #list}.  The overflow list preserves insertion
 * order until {@link #sort()} sorts it for percentile and mode calculations.</p>
 *
 * @author Brian Bushnell
 */
public class SuperLongList {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public SuperLongList(){
		this(100000);
	}

	public SuperLongList(int limit_){
		limit=limit_;
		array=new long[limit];
		list=new LongList();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns the live count array for values below {@link #limit}; direct edits
	 * do not update {@link #count} or {@link #sum} and require caller
	 * coordination.
	 */
	public long[] array(){return array;}

	/**
	 * Returns the live list of individual values at or above {@link #limit};
	 * direct edits do not update {@link #count} or {@link #sum} and require
	 * caller coordination.
	 */
	public LongList list(){return list;}

	/**
	 * Adds this histogram into {@code ca}, clamping values beyond its last
	 * index into that last cell.
	 */
	public void addTo(long[] ca){
		final int max=ca.length-1;
		{
			for(int i=0; i<array.length; i++){
				ca[Tools.min(i, max)]+=array[i];
			}
		}
		{
			final int listSize=list.size;
			final long[] listArray=list.array;
			for(int i=0; i<listSize; i++){
				long value=listArray[i];
				ca[(int)Tools.min(value, max)]++;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds each live value from {@code ll} as one sample in this histogram.
	 */
	public void add(LongList ll){
		for(int i=0; i<ll.size; i++){
			add(ll.get(i));
		}
	}

	/**
	 * Adds histogram counts where each array index is the represented value and
	 * each cell is the number of samples at that value.
	 */
	public void addCounts(long[] counts){
		for(int i=0; i<counts.length; i++){
			long x=counts[i];
			assert(x>=0);
			if(x>0){increment(i, x);}
		}
	}

	/** Adds one sample with value {@code x}. */
	public void add(long x){
		if(x<limit){array[(int)x]++;}
		else{list.add(x);}
		sum+=x;
		count++;
	}
	/** Adds one sample with value {@code x}. */
	public void increment(long x){add(x);}

	/** Adds {@code amt} samples with value {@code x}. */
	public void increment(long x, long amt){
		assert(amt>=0) : "SLL does not support decrements.";
		if(x<limit){array[(int)x]+=amt;}
		else{
			for(int i=0; i<amt; i++){list.add(x);}
		}
		sum+=x*amt;
		count+=amt;
	}

	/**
	 * Adds another histogram into this one, reclassifying values when its
	 * {@code limit} differs and otherwise appending its overflow list directly.
	 */
	public void add(SuperLongList sllT){
		if(array.length==sllT.array.length){//Fast, expected case
			assert(array.length==sllT.array.length) : "Array lengths must match.";
			for(int i=0; i<sllT.array.length; i++){
				array[i]+=sllT.array[i];
			}
			list.append(sllT.list);
			count+=sllT.count;
			sum+=sllT.sum;
		}else{//Slower generic case of unequal SLLs
			addCounts(sllT.array);
			add(sllT.list);
		}
	}
	public void incrementBy(SuperLongList sllT){add(sllT);}

	/** Sorts only the overflow list; in-array counts are already index ordered. */
	public void sort(){
		list.sort();
//		sorted=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the population standard deviation of all values added. */
	public double stdev(){
		final long div=Tools.max(1, count);
		double avg=sum/(double)div;
		double sumdev2=0;
		for(int i=0; i<array.length; i++){
			double dev=avg-i;
			double dev2=dev*dev;
			sumdev2+=(array[i]*dev2);
		}
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			double dev=avg-x;
			double dev2=dev*dev;
			sumdev2+=dev2;
		}
		return Math.sqrt(sumdev2/div);
	}

	/**
	 * Returns the value at a count percentile, after the overflow list is sorted.
	 */
	public long percentileValueByCount(double percentile){
//		assert(sorted);
		long thresh=(long)(count*percentile);
		long currentSum=0;
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			currentSum+=(x*i);
			currentCount+=x;
			if(currentCount>=thresh){return i;}
		}
		long prev=-1;
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			assert(x>=prev) : "Needs to be sorted ascending.";
			currentSum+=x;
			currentCount++;
			if(currentCount>=thresh){return x;}
			prev=x;
		}
		assert(false) : percentile+", "+count+", "+sum;
		return 0;
	}

	/**
	 * Returns the value that crosses a cumulative-sum percentile, after the
	 * overflow list is sorted.
	 */
	public long percentileValueBySum(double percentile){
//		assert(sorted);
		long thresh=(long)(sum*percentile);
		long currentSum=0;
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			currentSum+=(x*i);
			currentCount+=i;
			if(currentSum>=thresh){return i;}
		}
		long prev=-1;
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			assert(x>=prev) : "Needs to be sorted ascending.";
			currentSum+=x;
			currentCount++;
			if(currentSum>=thresh){return x;}
			prev=x;
		}
		assert(false) : percentile+", "+count+", "+sum;
		return 0;
	}

	/**
	 * Returns the sum of the lowest values through a count percentile, after the
	 * overflow list is sorted.
	 */
	public long percentileSumByCount(double percentile){
//		assert(sorted);
		long thresh=(long)(count*percentile);
		long currentSum=0;
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			currentSum+=(x*i);
			currentCount+=x;
			if(currentCount>=thresh){
				currentSum-=(x*i);
				currentCount-=x;
				while(currentCount<thresh){
					currentSum+=i;
					currentCount++;
				}
				return currentSum;
			}
		}
		long prev=-1;
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			assert(x>=prev) : "Needs to be sorted ascending.";
			currentSum+=x;
			currentCount++;
			if(currentCount>=thresh){return currentSum;}
			prev=x;
		}
		assert(false) : percentile+", "+count+", "+sum;
		return 0;
	}

	/**
	 * Returns how many lowest values are needed to cross a cumulative-sum
	 * percentile, after the overflow list is sorted.
	 */
	public long percentileCountBySum(double percentile){
//		assert(sorted);
		long thresh=(long)(sum*percentile);
		long currentSum=0;
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			currentSum+=(x*i);
			currentCount+=x;
			if(currentSum>=thresh){
				currentSum-=(x*i);
				currentCount-=x;
				while(currentSum<thresh){
					currentSum+=i;
					currentCount++;
				}
				return currentCount;
			}
		}
		long prev=-1;
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			assert(x>=prev) : "Needs to be sorted ascending.";
			currentSum+=x;
			currentCount++;
			if(currentSum>=thresh){return currentCount;}
			prev=x;
		}
		assert(false) : percentile+", "+count+", "+sum;
		return 0;
	}

	/** Returns the highest present value by scanning the overflow list if needed. */
	// Slow, avoid using.
	public long max(){
		if(list.size>0){return list.max();}
		for(int i=array.length-1; i>=0; i--){
			if(array[i]>0){return i;}
		}
		return 0;
	}

	/** Returns the arithmetic mean, or 0 when no values have been added. */
	public double mean(){
		return sum/Tools.max(1.0, count);
	}

	/**
	 * Returns the median value by count, using the lower of the two middle values
	 * for even counts, after the overflow list is sorted.
	 */
	public long median(){
		return countRankValue(count/2+(count&1));
	}

	/**
	 * Returns the value that crosses the 1-based count-rank target. A target of 0
	 * preserves empty-histogram behavior, and overflow values must be sorted so
	 * the rank walk sees them in ascending order.
	 */
	private long countRankValue(long target){
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			currentCount+=array[i];
			if(currentCount>=target){return i;}
		}

		long prev=-1;
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			assert(x>=prev) : "Needs to be sorted ascending.";
			currentCount++;
			if(currentCount>=target){return x;}
			prev=x;
		}
		assert(false) : "Rank not reached for consistent nonempty histogram: "
				+target+", "+count+", "+sum;
		return 0;
	}

	/** Returns the most frequent value, after the overflow list is sorted. */
	public long mode(){
		long maxCount=0;
		long maxValue=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			if(x>maxCount){
				maxCount=x;
				maxValue=i;
			}
		}

		long prev=-1;
		long currentCount=0;
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			if(x==prev){
				currentCount++;
				if(currentCount>maxCount){
					maxCount=currentCount;
					maxValue=x;
				}
			}else{
				assert(x>prev) : "Needs to be sorted ascending.";
				prev=x;
				currentCount=1;
				if(currentCount>maxCount){
					maxCount=currentCount;
					maxValue=x;
				}
			}
		}
		return maxValue;
	}

	/*--------------------------------------------------------------*/
	/*----------------           toString           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns all values in expanded histogram form, with overflow values in
	 * their current list order.
	 */
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		bb.append('[');
		String comma="";
		for(int i=0; i<array.length; i++){
			long value=array[i];
			for(long j=0; j<value; j++){
				bb.append(comma).append(i);
				comma=", ";
			}
		}
		for(int i=0; i<list.size; i++){
			bb.append(comma).append(list.get(i));
			comma=", ";
		}
		bb.append(']');
		return bb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public long count(){return count;}
	public long sum(){return sum;}

	private long count;
	private long sum;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	final long[] array;
	final LongList list;
	final int limit;

}
