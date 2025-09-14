package structures;

import java.util.Arrays;

import shared.KillSwitch;


/** 
 * Similar to an ArrayList<int[]> but intended to hold sets.
 * Faster insert could be handled via binary search if sorted,
 * or binary search for the last filled position if 
 * all elements are unique. */
public final class IntList2{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Re-call with default initial size. */
	public IntList2(){this(256);}
	
	/** Construct an IntList3 with this initial size.*/
	public IntList2(int initialSize){
		assert(initialSize>0);
		entries=new int[initialSize][];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Add this entry to the end of the list */
	public final void add(int[] entry){
		if(size>=entries.length){
			resize(max(size*2, 1));
		}
		entries[size]=entry;
		size++;
	}
	
	/** Added for better IntList3 compatibility */
	public final void add(int[] entry, int len){
		if(size>=entries.length){
			resize(max(size*2, 1));
		}
		entries[size]=entry;
		size++;
	}
	
	/** Set this location to specified entry */
	public final void set(int loc, int[] entry){
		if(loc>=entries.length){
			resize((loc+1)*2);
		}
		entries[loc]=entry;
		size=max(size, loc+1);
	}
	
	/** 
	 * Add this value to the specified location. 
	 * If an entry exists, insert the value, enlarging if necessary.
	 * Otherwise, create a new entry. */
	public final int insertIntoList(final int v, final int loc){
		
		if(loc>=size){
			assert(loc==size);
			add(null);
		}
			
		int[] entry=get(loc);
		if(entry==null){
			set(loc, new int[] {v, INVALID});
			return 1;
		}
		
		for(int i=0; i<entry.length; i++){//This is the slow bit; accelerate by hopping to the middle
			if(entry[i]==v){return 0;}
			if(entry[i]==INVALID){entry[i]=v;return 1;}
		}
		
		final int oldSize=entry.length;
		entry=KillSwitch.copyAndFill(entry, oldSize*2L, INVALID);
		set(loc, entry);
		
		//Quick add
		assert(entry[oldSize]==INVALID);
		entry[oldSize]=v;
		return 1;
		
		//Old code
//		final int oldSize=entry.length;
//		final int newSize=(int)Tools.min(Shared.MAX_ARRAY_LEN, oldSize*2L);
//		assert(newSize>entry.length) : "Overflow.";
//		entry=KillSwitch.copyOf(entry, newSize);
//		entry[oldSize]=v;
//		Arrays.fill(entry, oldSize+1, newSize, -1);
//		set(loc, entry);
//		return 1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Resizing           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Resize the array of entries */
	public final void resize(int size2){
		assert(size2>size);
		entries=KillSwitch.copyOf(entries, size2);
	}
	
	/** Compress the list by eliminating the unused trailing space */
	public final void shrink(){
		if(size==entries.length){return;}
		entries=KillSwitch.copyOf(entries, size);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Get the entry at the location */
	public final int[] get(int loc){
		return(loc>=size ? null : entries[loc]);
	}
	
	/** Added for better IntList3 compatibility */
	public int getLen(int i) {
		return i>=size ? 0 : entries[i]==null ? 0 : -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<size; i++){
			if(entries[i]!=null){
				sb.append(comma+"("+i+", "+Arrays.toString(entries[i])+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Holds entries.  Each entry is a sets of numbers in an int[]. 
	 * Leftmost values are valid, rightmost values are invalid. */
	private int[][] entries;
	/** Number of entries in the primary array. */
	public int size=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This could be made mutable per instance, but it would be a lot of work. */
	public static final int INVALID=-1;
	
}
