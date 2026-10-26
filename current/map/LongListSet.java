package map;

import java.util.Arrays;

import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;

/**
 * Partitioned long accumulator that can exceed a single array's element limit.
 * Values are assigned to three LongLists by their sign-stripped value modulo three.
 * Threshold-triggered sorting/condensing reduces duplicates, but duplicates can
 * remain between compactions: call sort() then condense() for unique contents.
 * Unsynchronized. Public partitions are borrowed mutable storage; external edits
 * must preserve the sorted-state assumptions. Iterators are live, not snapshots,
 * and do not detect concurrent modification.
 *
 * @author Brian Bushnell
 * @date January 8, 2021
 */
public class LongListSet{

	/**
	 * Simple demo that builds a set, adds values, sorts, and condenses.
	 * @param args ignored
	 */
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

	/*--------------------------------------------------------------*/
	/*----------------          Formatting          ----------------*/
	/*--------------------------------------------------------------*/

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

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs an empty accumulator with three partitions and initial condense thresholds of 64.
	 */
	public LongListSet(){
		array=new LongList[mod];
		for(int i=0; i<mod; i++){
			array[i]=new LongList(32);
		}
		nextCondense=new int[mod];
		Arrays.fill(nextCondense, 64);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds a value and compacts its partition when that partition hits its threshold.
	 * @param x Value to add */
	public void add(long x){
		int y=(int)((x&Long.MAX_VALUE)%mod); //partition by sign-stripped value % mod (mod=3 LongLists -> >2B total)
		LongList list=array[y];
		list.add(x);
		if(list.size>=nextCondense[y]){
			//Compacts only the target partition: sort locally, drop duplicates, then move the next threshold
			//toward 2x the post-condense size without lowering it. Duplicate-heavy streams may trigger repeated
			//checks at the same threshold, but only this capped partition is resorted when the threshold is hit.
			//This branch leaves that partition sorted; `sorted` stays as-is because a globally partition-sorted
			//state remains partition-sorted after resorting one partition.
			list.sort();
			list.condense();
			nextCondense[y]=(int)Tools.mid(nextCondense[y], list.size*2L, Shared.MAX_ARRAY_LEN);
		}else{
			sorted=false; //appended an unsorted tail -> the set is no longer globally sorted
		}
	}

	/** Sorts each partition independently; iteration still walks partitions in modulo order. */
	public void sort(){
		if(sorted){return;}
		for(LongList list : array){list.sort();}
		sorted=true;
	}

	/** Condenses duplicate values out of each sorted partition; call sort() first.
	 * @throws AssertionError If the internal sorted flag is false with assertions enabled */
	public void condense(){
		assert(sorted) : "Sort first.";
		for(LongList list : array){list.condense();}
	}

	/** Condenses each already-sorted partition and releases its spare capacity.
	 * Call sort() first: this method does not sort or verify the sorted flag itself. */
	public void shrinkToUnique(){
		for(LongList list : array){list.shrinkToUnique();}
	}

	/** Creates a live traversal of partitions in modulo order, not global value order.
	 * Do not mutate partitions during traversal; call hasMore() before each next().
	 * @return Newly allocated iterator borrowing this accumulator's storage */
	public LongListSetIterator iterator(){
		return new LongListSetIterator();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** True after every partition has been sorted since its last unsorted append. */
	private boolean sorted=false;

	/** Partitioned LongList storage. */
	public final LongList[] array;
	/** Per-partition size threshold that triggers opportunistic condense. */
	public final int[] nextCondense;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of modulo partitions. */
	public static final int mod=3;

	/*--------------------------------------------------------------*/
	/*----------------          Inner Classes       ----------------*/
	/*--------------------------------------------------------------*/

	/** Iterator that walks live partitions in array order and list elements in index order. */
	public class LongListSetIterator{

		/** Constructs an iterator at the first partition and element. */
		public LongListSetIterator(){}

		/**
		 * Returns the current value and advances within the current partition.
		 * Caller must have just observed {@link #hasMore()} true; this method does
		 * not seek to a populated partition.
		 * @return next value
		 */
		public long next(){
			long x=array[i].get(j);
			j++;
			return x;
		}

		/**
		 * Advances to the next populated slot if needed.
		 * @return true if next() can read another value
		 */
		public boolean hasMore(){
			return findNextValid();
		}

		/**
		 * Advances one slot before seeking the next populated slot.
		 * @return true if a later value is available
		 */
		boolean advance(){
			j++;
			return findNextValid();
		}

		//Positions (i,j) at the next populated slot, skipping exhausted/empty partitions; false when done.
		//Idempotent once positioned (the common-case fast path returns true without advancing), so calling
		//hasMore() repeatedly before a single next() is safe. Contract: call hasMore() before each next().
		/**
		 * Moves to the next partition that still has an unread element.
		 * @return true if (i,j) points at a readable element
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
		/** Current partition index and current element index within that partition. */
		private int i=0, j=0;

	}

}
