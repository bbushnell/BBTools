package structures;

import java.util.Arrays;

import shared.KillSwitch;


/**
 * Dynamic array for holding sets of integers, similar to ArrayList<int[]>.
 * Optimized for scenarios where each entry contains a set of unique values.
 * Each entry is an int[] where valid values are at the beginning and
 * invalid positions are marked with INVALID (-1).
 *
 * @author Brian Bushnell
 */
public final class IntList2{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a new IntList2 with default initial capacity of 256 entries */
	public IntList2(){this(256);}
	
	/** Creates a new IntList2 with specified initial capacity.
	 * @param initialSize Initial number of entries the list can hold */
	public IntList2(int initialSize){
		assert(initialSize>0);
		entries=new int[initialSize][];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Appends an entry to the end of the list.
	 * Automatically resizes the underlying array if necessary.
	 * @param entry The int array to add as a new entry
	 */
	public final void add(int[] entry){
		if(size>=entries.length){
			resize(max(size*2, 1));
		}
		entries[size]=entry;
		size++;
	}
	
	/**
	 * Appends an entry to the end of the list with length parameter for IntList3 compatibility.
	 * The length parameter is currently ignored but maintained for interface compatibility.
	 * @param entry The int array to add as a new entry
	 * @param len Length parameter (currently unused)
	 */
	public final void add(int[] entry, int len){
		if(size>=entries.length){
			resize(max(size*2, 1));
		}
		entries[size]=entry;
		size++;
	}
	
	/**
	 * Sets the entry at the specified location.
	 * Automatically resizes the underlying array and updates size if necessary.
	 * @param loc The index where to place the entry
	 * @param entry The int array to place at the specified location
	 */
	public final void set(int loc, int[] entry){
		if(loc>=entries.length){
			resize((loc+1)*2);
		}
		entries[loc]=entry;
		size=max(size, loc+1);
	}
	
	/**
	 * Inserts a value into the set at the specified location.
	 * If the location doesn't exist, creates a new entry.
	 * If the entry exists, attempts to insert the value if not already present.
	 * Automatically grows the entry array if it becomes full.
	 *
	 * @param v The value to insert into the set
	 * @param loc The location/index where the set resides
	 * @return 1 if the value was inserted (new), 0 if already present
	 */
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
	
	/** Resizes the underlying entries array to accommodate more entries.
	 * @param size2 The new size for the entries array, must be larger than current size */
	public final void resize(int size2){
		assert(size2>size);
		entries=KillSwitch.copyOf(entries, size2);
	}
	
	/** Compresses the entries array by removing unused trailing capacity */
	public final void shrink(){
		if(size==entries.length){return;}
		entries=KillSwitch.copyOf(entries, size);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Retrieves the entry at the specified location.
	 * @param loc The index of the entry to retrieve
	 * @return The int array at the specified location, or null if out of bounds
	 */
	public final int[] get(int loc){
		return(loc>=size ? null : entries[loc]);
	}
	
	/**
	 * Gets the length of an entry for IntList3 compatibility.
	 * Returns 0 if index is out of bounds or entry is null, -1 if entry exists.
	 * @param i The index to check
	 * @return 0 if out of bounds or null entry, -1 if entry exists
	 */
	public int getLen(int i) {
		return i>=size ? 0 : entries[i]==null ? 0 : -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Returns a string representation showing all non-null entries with their indices.
	 * Format: [(index, [values]), (index, [values]), ...]
	 * @return String representation of the list contents
	 */
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
	
	private int[][] entries;
	public int size=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Sentinel value (-1) used to mark invalid/empty positions in entry arrays */
	public static final int INVALID=-1;
	
}
