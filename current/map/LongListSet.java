package map;

import java.util.Arrays;

import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;

/**
 * A set of LongLists designed to increase LongList capacity beyond 2B elements.
 * Auto-condenses to avoid representing multiple copies of the same value.
 * Uses modular hashing to distribute elements across multiple LongList instances
 * for scalability and automatic memory management.
 *
 * @author Brian Bushnell
 * @date January 8, 2021
 */
public class LongListSet{
	
	/** Simple demo that builds a set, adds values, sorts, and condenses. */
	public static void main(String[] args){
		LongListSet set=new LongListSet();
		set.add(1);
		set.add(2);
		set.add(3);
		set.add(4);
		set.add(5);
		set.add(2);
		set.add(2);
		set.add(5);
		System.err.println(set);
		set.sort();
		set.condense();
		System.err.println(set);
	}
	
	/** Returns a comma-separated string representation of all elements in iteration order.
	 * @return String form of the set */
	public String toString(){
		LongListSetIterator iter=iterator();
		ByteBuilder bb=new ByteBuilder();
		bb.append('[');
		while(iter.hasMore()){
			long x=iter.next();
			bb.append(x);
			bb.append(',');
		}
		if(bb.endsWith(',')){bb.setLength(bb.length-1);}
		bb.append(']');
		return bb.toString();
	}
	
	/**
	 * Constructs a new LongListSet with default mod partitions and condense thresholds.
	 */
	public LongListSet(){
		array=new LongList[mod];
		for(int i=0; i<mod; i++){
			array[i]=new LongList(32);
		}
		nextCondense=new int[mod];
		Arrays.fill(nextCondense, 64);
	}
	
	/** Adds a value, triggering sort/condense when the target partition hits its threshold.
	 * @param x Value to add */
	public void add(long x){
		int y=(int)((x&Long.MAX_VALUE)%mod);
		LongList list=array[y];
		list.add(x);
		if(list.size>=nextCondense[y]){
			list.sort();
			list.condense();
			nextCondense[y]=(int)Tools.mid(nextCondense[y], list.size*2L, Shared.MAX_ARRAY_LEN);
		}else{
			sorted=false;
		}
	}
	
	public void sort(){
		if(sorted){return;}
		for(LongList list : array){list.sort();}
		sorted=true;
	}
	
	public void condense(){
		assert(sorted) : "Sort first.";
		for(LongList list : array){list.condense();}
	}
	
	public void shrinkToUnique(){
		for(LongList list : array){list.shrinkToUnique();}
	}
	
	/** Returns an iterator over all elements across partitions.
	 * @return Iterator over set contents */
	public LongListSetIterator iterator(){
		return new LongListSetIterator();
	}
	
	private boolean sorted=false;
	
	public final LongList[] array;
	public final int[] nextCondense;
	
	public static final int mod=3;
	
	public class LongListSetIterator{
		
		//Assumes hasMore() has already been called and returned true
		public long next(){
			long x=array[i].get(j);
			j++;
			return x;
		}
		
		public boolean hasMore(){
			return findNextValid();
		}
		
		boolean advance(){
			j++;
			return findNextValid();
		}
		
		boolean findNextValid(){
			if(i<mod && j<array[i].size){return true;}//common case
			while(i<mod){
				if(j<array[i].size){return true;}
				i++;
				j=0;
			}
			return false;
		}
		/** Current element index within the current LongList */
		private int i=0, j=0;
		
	}
	
}
