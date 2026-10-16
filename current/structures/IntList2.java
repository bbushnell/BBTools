package structures;

import java.util.Arrays;

import shared.KillSwitch;


/**
 * Dynamic table of live int arrays that encode sets using INVALID-terminated
 * spare capacity.
 * Null rows are empty; non-null rows are scanned until a matching value, an
 * INVALID sentinel, or the end of a full row is found.
 * Access to mutation and the live arrays returned by get() must be coordinated
 * by the caller.
 *
 * @author Brian Bushnell
 */
public final class IntList2{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs an IntList2 with the default initial entry-table capacity. */
	public IntList2(){this(256);}

	/**
	 * Constructs an IntList2 with the specified initial entry-table capacity.
	 * @param initialSize Initial number of entry slots to allocate.
	 */
	public IntList2(int initialSize){
		assert(initialSize>0);
		entries=new int[initialSize][];
	}

	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Appends a live entry array at the next location.
	 * @param entry Integer array to store.
	 */
	public final void add(int[] entry){
		if(size>=entries.length){
			resize(max(size*2, 1));
		}
		entries[size]=entry;
		size++;
	}

	/**
	 * Appends a live entry array and ignores len for IntList3 compatibility.
	 * @param entry Integer array to store.
	 * @param len Ignored compatibility parameter.
	 */
	public final void add(int[] entry, int len){
		if(size>=entries.length){
			resize(max(size*2, 1));
		}
		entries[size]=entry;
		size++;
	}

	/**
	 * Stores a live entry array at loc, growing the entry table when loc is
	 * beyond the current capacity.
	 * @param loc Target location index.
	 * @param entry Integer array to store.
	 */
	public final void set(int loc, int[] entry){
		if(loc>=entries.length){
			resize((loc+1)*2);
		}
		entries[loc]=entry;
		size=max(size, loc+1);
	}

	/**
	 * Inserts a value into the entry at loc, creating that entry if needed.
	 * Existing entries are scanned linearly until v, INVALID, or the end of a
	 * full row is found.
	 * @param v Integer value to insert.
	 *        Must differ from INVALID.
	 * @param loc Target location index; loc may append exactly one new slot.
	 * @return 1 if value was added, or 0 if already present before INVALID.
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
		final long newSize=(oldSize>0 ? oldSize*2L : 2);
		entry=KillSwitch.copyAndFill(entry, newSize, INVALID);
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
	 * Resizes the entry table.
	 * @param size2 New capacity; must exceed the current logical size.
	 */
	public final void resize(int size2){
		assert(size2>size);
		entries=KillSwitch.copyOf(entries, size2);
	}

	/** Removes unused trailing slots from the entry table. */
	public final void shrink(){
		if(size==entries.length){return;}
		entries=KillSwitch.copyOf(entries, size);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Gets the live entry array at loc.
	 * @param loc Location index.
	 * @return Entry array, or null if loc is past the logical table extent.
	 */
	public final int[] get(int loc){
		return(loc>=size ? null : entries[loc]);
	}

	/**
	 * Gets the IntList3-compatible length marker for one entry.
	 * @param i Location index.
	 * @return -1 if an entry exists at i, or 0 if i is past the logical table
	 * extent or the entry is null.
	 */
	public int getLen(int i){
		return i>=size ? 0 : entries[i]==null ? 0 : -1;
	}

	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns a string representation showing all non-null entries with their indices.
	 * The complete backing arrays are printed, including INVALID spare slots.
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

	/** Sentinel value used to mark invalid or empty positions in entry arrays. */
	public static final int INVALID=-1;

}
