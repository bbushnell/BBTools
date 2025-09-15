package structures;

import shared.Tools;

/**
 * Holds counts of numbers for histograms.
 * Small numbers are stored as counts in an array.
 * Large numbers are stored individually in a LongList.
 * @author Brian Bushnell
 *
 */
public class SuperLongList {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public SuperLongList(){
		this(100000);
	}
	
	/**
	 * @param limit_ Length of array used to store small numbers.
	 */
	public SuperLongList(int limit_){
		limit=limit_;
		array=new long[limit];
		list=new LongList();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	public long[] array(){return array;}
	public LongList list(){return list;}
	
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
	
	/** Add or increment these keys. */
	public void add(LongList ll){
		for(int i=0; i<ll.size; i++){
			add(ll.get(i));
		}
	}
	
	/** Add or increment these counts. */
	public void addCounts(long[] counts){
		for(int i=0; i<counts.length; i++){
			long x=counts[i];
			assert(x>=0);
			if(x>0){increment(i, x);}
		}
	}
	
	/** Add or increment this key.
	 * May result in duplicate copies appearing, which is fine. */
	public void add(long x){
		if(x<limit){array[(int)x]++;}
		else{list.add(x);}
		sum+=x;
		count++;
	}
	public void increment(long x){add(x);}
	
	public void increment(long x, long amt){
		assert(amt>=0) : "SLL does not support decrements.";
		if(x<limit){array[(int)x]+=amt;}
		else{
			for(int i=0; i<amt; i++){list.add(x);}
		}
		sum+=x*amt;
		count+=amt;
	}
	
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
	public void incrementBy(SuperLongList sllT) {add(sllT);}
	
	public void sort() {
		list.sort();
//		sorted=true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	/** Returns value such that percentile of values are below that value */
	public long percentileValueByCount(double percentile){
//		assert(sorted);
		long thresh=(long)(count*percentile);
		long currentSum=0;
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			currentSum+=(x*i);
			currentCount+=i;
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
	
	/** Returns value such that percentile of sum of values are below that value */
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

	/** Returns the sum of the lower percentile of values */
	public long percentileSumByCount(double percentile){
//		assert(sorted);
		long thresh=(long)(count*percentile);
		long currentSum=0;
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			currentSum+=(x*i);
			currentCount+=i;
			if(currentCount>=thresh){
				currentSum-=(x*i);
				currentCount-=i;
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

	/** Returns the number of lower values needed to sum to this percentile of the total sum */
	public long percentileCountBySum(double percentile){
//		assert(sorted);
		long thresh=(long)(sum*percentile);
		long currentSum=0;
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			currentSum+=(x*i);
			currentCount+=i;
			if(currentSum>=thresh){
				currentSum-=(x*i);
				currentCount-=i;
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
	
	//Slow, avoid using
	public long max(){
		if(list.size>0){return list.max();}
		for(int i=array.length-1; i>=0; i--){
			if(array[i]>0){return i;}
		}
		return 0;
	}
	
	public double mean(){
		return sum/Tools.max(1.0, count);
	}
	
	public long median(){
		return percentileValueByCount(0.5);
	}
	
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
			}
		}
		return maxValue;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           toString           ----------------*/
	/*--------------------------------------------------------------*/
	
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

	public long count() {return count;}
	public long sum() {return sum;}
	
	private long count;
	private long sum;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final long[] array;
	final LongList list;
	final int limit;
	
}
