package structures;

import java.util.Arrays;

import shared.KillSwitch;


/**
 * Dynamic array-like data structure for holding sets of long values.
 * Similar to ArrayList&lt;long[]&gt; but optimized for set operations.
 * Each entry is a long array where leftmost values are valid and rightmost values are INVALID (-1).
 * Provides fast insertion and automatic resizing capabilities.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public final class LongList2{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructs a LongList2 with default initial size of 256 entries */
	public LongList2(){this(256);}
	
	/** Constructs a LongList2 with specified initial capacity.
	 * @param initialSize Initial number of entry slots to allocate */
	public LongList2(int initialSize){
		assert(initialSize>0);
		entries=new long[initialSize][];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds an entry to the end of the list.
	 * Automatically resizes the internal array if necessary.
	 * @param entry The long array to add
	 */
	public final void add(long[] entry){
		if(size>=entries.length){
			resize(max(size*2, 1));
		}
		entries[size]=entry;
		size++;
	}
	
	/**
	 * Adds an entry to the end of the list for IntList3 compatibility.
	 * The len parameter is ignored but maintained for interface consistency.
	 * @param entry The long array to add
	 * @param len Length parameter (ignored)
	 */
	public final void add(long[] entry, int len){
		if(size>=entries.length){
			resize(max(size*2, 1));
		}
		entries[size]=entry;
		size++;
	}
	
	/**
	 * Sets the entry at the specified location.
	 * Automatically resizes if location exceeds current capacity.
	 * Updates size to accommodate the new location.
	 *
	 * @param loc Index position to set
	 * @param entry The long array to place at this position
	 */
	public final void set(int loc, long[] entry){
		if(loc>=entries.length){
			resize((loc+1)*2);
		}
		entries[loc]=entry;
		size=max(size, loc+1);
	}
	
	/**
	 * Inserts a value into the set at the specified location.
	 * If no entry exists at location, creates a new entry.
	 * If entry exists, adds value if not already present, expanding entry if needed.
	 * Uses INVALID (-1) as sentinel value for unused slots.
	 *
	 * @param v The long value to insert
	 * @param loc The list location where the set resides
	 * @return 1 if value was added, 0 if value already existed
	 */
	public final int insertIntoList(final long v, final int loc){
		
		if(loc>=size){
			assert(loc==size);
			add(null);
		}
			
		long[] entry=get(loc);
		if(entry==null){
			set(loc, new long[] {v, INVALID});
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
	
	/**
	 * Resizes the internal entries array to the specified size.
	 * New size must be larger than current size.
	 * @param size2 New array capacity
	 */
	public final void resize(int size2){
		assert(size2>size);
		entries=KillSwitch.copyOf(entries, size2);
	}
	
	/** Compresses the list by removing unused trailing capacity */
	public final void shrink(){
		if(size==entries.length){return;}
		entries=KillSwitch.copyOf(entries, size);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Retrieves the entry at the specified location.
	 * @param loc Index position to retrieve
	 * @return The long array at this position, or null if location is out of bounds
	 */
	public final long[] get(int loc){
		return(loc>=size ? null : entries[loc]);
	}
	
	/**
	 * Gets the length of the entry at specified index for IntList3 compatibility.
	 * Returns 0 if index is out of bounds or entry is null, -1 if entry exists.
	 * @param i Index to check
	 * @return -1 if entry exists, 0 if out of bounds or null
	 */
	public int getLen(int i) {
		return i>=size ? 0 : entries[i]==null ? 0 : -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a string representation showing non-null entries with their indices.
	 * Format: [(index, [array_contents]), ...]
	 * @return String representation of the list
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
	
	private long[][] entries;
	public int size=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Sentinel value (-1) used to mark invalid/unused slots in entry arrays */
	public static final int INVALID=-1;
	
}
