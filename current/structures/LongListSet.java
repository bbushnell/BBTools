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
	
	/**
	 * Returns string representation of all elements in the set.
	 * Uses iterator to traverse all LongLists and formats as comma-separated
	 * values enclosed in brackets.
	 * @return String representation in format [val1,val2,val3]
	 */
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
	 * Constructs a new LongListSet with default capacity.
	 * Initializes the underlying array of LongLists with initial capacity
	 * and sets up auto-condensing thresholds.
	 */
	public LongListSet(){
		array=new LongList[mod];
		for(int i=0; i<mod; i++){
			array[i]=new LongList(32);
		}
		nextCondense=new int[mod];
		Arrays.fill(nextCondense, 64);
	}
	
	/**
	 * Adds a long value to the appropriate LongList based on hash modulo.
	 * Triggers automatic sorting and condensing when the target list reaches
	 * its condensing threshold to maintain performance and uniqueness.
	 * @param x The long value to add to the set
	 */
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
	
	/** Sorts all LongLists in the set if not already sorted.
	 * Required before condensing to remove duplicates efficiently. */
	public void sort(){
		if(sorted){return;}
		for(LongList list : array){list.sort();}
		sorted=true;
	}
	
	/** Removes duplicate values from all LongLists in the set.
	 * Must be called after sort() to ensure proper functionality. */
	public void condense(){
		assert(sorted) : "Sort first.";
		for(LongList list : array){list.condense();}
	}
	
	/** Shrinks all LongLists to contain only unique values.
	 * Calls shrinkToUnique() on each underlying LongList to minimize memory usage. */
	public void shrinkToUnique(){
		for(LongList list : array){list.shrinkToUnique();}
	}
	
	/**
	 * Returns an iterator for traversing all elements in the set.
	 * The iterator will visit elements from all underlying LongLists in order.
	 * @return A new LongListSetIterator for this set
	 */
	public LongListSetIterator iterator(){
		return new LongListSetIterator();
	}
	
	/** Tracks whether all LongLists in the set are currently sorted */
	private boolean sorted=false;
	
	/** Array of LongLists that store the actual set elements */
	public final LongList[] array;
	/** Size thresholds for triggering automatic condensing on each LongList */
	public final int[] nextCondense;
	
	/** Modulus value used for hash distribution across LongList array */
	public static final int mod=3;
	
	/** Iterator for traversing all elements across all LongLists in the set.
	 * Maintains position state and provides sequential access to all values. */
	public class LongListSetIterator{
		
		//Assumes hasMore() has already been called and returned true
		/**
		 * Returns the next element in the iteration sequence.
		 * Assumes hasMore() has been called and returned true before calling this method.
		 * @return The next long value in the set
		 */
		public long next(){
			long x=array[i].get(j);
			j++;
			return x;
		}
		
		/**
		 * Checks if there are more elements to iterate over.
		 * Advances internal position to the next valid element if necessary.
		 * @return true if more elements are available, false otherwise
		 */
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
		/** Current element index within the current LongList */
		/** Current LongList index in the array */
		private int i=0, j=0;
		
	}
	
}
