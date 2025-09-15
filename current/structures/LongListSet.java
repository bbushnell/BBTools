package structures;

import java.util.Arrays;

import shared.Shared;
import shared.Tools;

/**
 * A set of LongLists, designed to increase LongList capacity beyond 2B.
 * Auto-condenses; e.g., not intended to represent multiple copies of a value.
 * 
 * @author Brian Bushnell
 * @date January 8, 2021
 *
 */
public class LongListSet{
	
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
	
	public LongListSet(){
		array=new LongList[mod];
		for(int i=0; i<mod; i++){
			array[i]=new LongList(32);
		}
		nextCondense=new int[mod];
		Arrays.fill(nextCondense, 64);
	}
	
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
		
		/** 
		 * Increment and point to next valid element.
		 * @return true if there is a valid element.
		 */
		boolean advance(){
			j++;
			return findNextValid();
		}
		
		/** 
		 * Point to next valid element.
		 * If already valid, do nothing.
		 * @return true if there is a valid element.
		 */
		boolean findNextValid(){
			if(i<mod && j<array[i].size){return true;}//common case
			while(i<mod){
				if(j<array[i].size){return true;}
				i++;
				j=0;
			}
			return false;
		}
		private int i=0, j=0;
		
	}
	
}
