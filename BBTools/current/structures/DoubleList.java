package structures;

import java.util.ArrayList;
import java.util.LinkedList;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;



public final class DoubleList{
	
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
	
	public DoubleList(){this(256);}
	
	public DoubleList(int initial){
//		assert(initial>0) : initial+"\n"+this;
		initial=Tools.max(initial, 1);
		array=KillSwitch.allocDouble1D(initial);
	}
	
	public DoubleList copy() {
		DoubleList copy=new DoubleList(size);
		copy.addAll(this);
		return copy;
	}
	
	public void clear(){size=0;}
	
	public final void set(int loc, double value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	public final void setLast(double value){
		assert(size>0);
		array[size-1]=value;
	}
	
	public final void increment(int loc, double value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	public final void increment(int loc){
		increment(loc, 1);
	}
	
	public final void incrementBy(DoubleList b){
		for(int i=b.size-1; i>=0; i--){
			increment(i, b.get(i));
		}
	}
	
	public final void incrementBy(double[] b){
		for(int i=b.length-1; i>=0; i--){
			increment(i, b[i]);
		}
	}
	
	public final void append(DoubleList b){
		for(int i=0; i<b.size; i++){
			add(b.get(i));
		}
	}
	
	public final void append(double[] b){
		for(int i=0; i<b.length; i++){
			add(b[i]);
		}
	}
	
	public void subtractFrom(double value){
		for(int i=0; i<size; i++){
			array[i]=value-array[i];
		}
	}
	
	public final double get(int loc){
		return(loc>=size ? 0 : array[loc]);//TODO: Shouldn't this crash instead of returning 0?
	}
	
	public double lastElement() {
		assert(size>0);
		return array[size-1];
	}
	
	public final void add(double x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	//Slow; for validation
	public boolean containsDuplicates(){
		for(int i=0; i<size; i++){
			for(int j=i+1; j<size; j++){
				if(array[i]==array[j]){return true;}
			}
		}
		return false;
	}
	
	public void addAll(DoubleList counts) {
		final double[] array2=counts.array;
		final int size2=counts.size;
		for(int i=0; i<size2; i++){add(array2[i]);}
	}
	
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
	
	public final double sumLong(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
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
	
	public double percentile(double fraction){
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
	
	public final void shrink(){
		if(size==array.length){return;}
		array=KillSwitch.copyOf(array, size);
	}
	

	
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
	
	public double[] toArray(){
		return KillSwitch.copyOf(array, size);
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
	
	public void sort() {
		if(size>1){Shared.sort(array, 0, size);}
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
	
	public boolean isEmpty() {
		return size<1;
	}
	
	public int capacity() {
		return array.length;
	}
	
	public int freeSpace() {
		return array.length-size;
	}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public double[] array;
	public int size=0;
	
}
