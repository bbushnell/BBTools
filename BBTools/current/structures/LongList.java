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
	
	public LongList(int initial){
		assert(initial>0) : initial;
		array=KillSwitch.allocLong1D(initial);
	}

	public void clear(){size=0;}
	public void clearFull(){
		Arrays.fill(array, 0);
		size=0;
	}
	
	public final void set(int loc, long value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	public final void setLast(long value){
		assert(size>0);
		array[size-1]=value;
	}
	
	public final void increment(int loc, long value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	public final void increment(int loc){
		increment(loc, 1);
	}
	
	public final void incrementBy(LongList b){
		for(int i=b.size-1; i>=0; i--){
			increment(i, b.get(i));
		}
	}
	
	public final void incrementBy(long[] b){
		for(int i=b.length-1; i>=0; i--){
			increment(i, b[i]);
		}
	}
	
	public final void append(LongList b){
		for(int i=0; i<b.size; i++){
			add(b.get(i));
		}
	}
	
	public final void append(long[] b){
		for(int i=0; i<b.length; i++){
			add(b[i]);
		}
	}
	
	public final long get(int loc){
		return(loc>=size ? 0 : array[loc]);
	}
	
	public final void add(long x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	public boolean contains(long x) {
		for(int i=0; i<size; i++) {
			if(array[i]==x) {return true;}
		}
		return false;
	}
	
	public void addAll(LongList list2) {
		final long[] array2=list2.array;
		final int size2=list2.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}
	
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		final int size3=(int)Tools.min(Shared.MAX_ARRAY_LEN, size2);
		assert(size3>size) : "Overflow: "+size+", "+size2+" -> "+size3;
		array=KillSwitch.copyOf(array, size3);
	}
	
	public final LongList shrink(){
		if(size==array.length){return this;}
		array=KillSwitch.copyOf(array, size);
		return this;
	}
	
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
	
	public final double avgDif(final double x){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=Tools.absdif(x, array[i]);
		}
		return sum/(Tools.max(1, size));
	}
	
	public final double rmsDif(final double x){
		double sum=0;
		for(int i=0; i<size; i++){
			double dif=Tools.absdif(x, array[i]);
			sum+=dif*dif;
		}
		return Math.sqrt(sum/(Tools.max(1, size)));
	}
	
	public final long sumLong(){
		long sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	public final long sumHist(){
		long sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i]*i;
		}
		return sum;
	}
	
	public final double sum(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	public final double mean(){
		return size<1 ? 0 : sum()/size;
	}
	
	public final double meanHist(){
		return size<1 ? 0 : sumHist()/sum();
	}
	
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
	
	public long percentile(double fraction){
		if(size<1){return 0;}
		int idx=percentileIndex(fraction);
		return array[idx];
	}
	
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
	
	public final void shrinkToUnique(){
		condense();
		shrink();
	}
	
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
	
	public long[] toArray(){
		long[] x=KillSwitch.allocLong1D(size);
		for(int i=0; i<x.length; i++){
			x[i]=array[i];
		}
		return x;
	}
	
	public void sort() {
		if(size>1){Shared.sort(array, 0, size);}
	}
	
	public void sortSerial() {
		if(size>1){Arrays.sort(array, 0, size);}
	}
	
	public void reverse() {
		if(size>1){Tools.reverseInPlace(array, 0, size);}
	}
	
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
	
	private static final long min(long x, long y){return x<y ? x : y;}
	private static final long max(long x, long y){return x>y ? x : y;}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public long[] array;
	/** Highest occupied index plus 1, i.e., lowest unoccupied index */
	public int size=0;
	
}
