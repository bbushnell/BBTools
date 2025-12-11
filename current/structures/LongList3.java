package structures;

import java.util.Arrays;

import shared.KillSwitch;


/**
 * A dynamic list of long arrays where each entry maintains its own size independently.
 * Similar to a list of LongList but more memory efficient. Designed for use with
 * HashArrayHybridFast where each entry represents a set of values with specific
 * ordering constraints (ascending or unique additions).
 * @author Brian Bushnell
 * @date Dec 10, 2018
 */
public final class LongList3{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates a LongList3 with default initial size and mode. */
	public LongList3(){this(defaultInitialSize, defaultMode);}
	
	/**
	 * Creates a LongList3 with specified initial capacity and ordering mode.
	 * @param initialSize Initial capacity for the entries array
	 * @param mode_ Ordering constraint mode (ASCENDING or UNIQUE)
	 */
	public LongList3(int initialSize, int mode_){
		assert(initialSize>0);
		entries=new long[initialSize][];
		sizes=new int[initialSize];
		
		mode=mode_;
		assert(mode==ASCENDING || mode==UNIQUE) : "Unsupported mode: "+mode;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds an entry to the end of the list.
	 * @param entry The long array to add
	 * @param len Number of valid elements in the entry
	 */
	public final void add(long[] entry, int len){
		set(size, entry, len);
	}
	
	/**
	 * Sets the entry at the specified location, resizing if necessary.
	 * @param loc The location to set
	 * @param entry The long array to store
	 * @param len Number of valid elements in the entry
	 */
	public final void set(int loc, long[] entry, int len){
		assert((entry==null && len==0) || (entry.length>0 && len<=entry.length)) : len+", "+(entry==null ? entry : ""+entry.length);
		
		if(loc>=entries.length){//Resize by doubling if necessary
			resize(max(size*2, loc+1));
		}
		entries[loc]=entry;
		sizes[loc]=len;
		size=max(size, loc+1);
	}
	
	/**
	 * Inserts a value into the entry at the specified location. Creates a new entry
	 * if none exists, or adds to existing entry if value is not already present.
	 * Optimized for ascending or unique insertion modes.
	 * @param v The value to insert
	 * @param loc The location to insert into
	 * @return 1 if value was added, 0 if already present
	 */
	public final int insertIntoList(final long v, final int loc){
		
		//If the location is empty
		if(loc>=size || entries[loc]==null){
			//Add a new entry and return
			set(loc, new long[] {v, INVALID}, 1);
			return 1;
		}
			
		long[] entry=get(loc);
		final int oldSize=sizes[loc];
		
		//Scan the latest entries to see if v is already present.
		for(int i=oldSize-1, lim=max(0, oldSize-slowAddLimit); i>=lim; i--){
			if(entry[i]==v){return 0;}
			assert(entry[i]<v || entry[i]!=v); //Ascending, the assumption for Seal; unique, the assumption for Sketch
			assert(entry[i]!=INVALID);
		}
		//At this point the element was not found because it was not present or the size is too big
		
		//If the entry is full, resize it
		if(oldSize>=entry.length){
			assert(oldSize==entry.length);
			entry=KillSwitch.copyAndFill(entry, oldSize*2L, INVALID);
			set(loc, entry, sizes[loc]);
		}
		
		//Quick add
		assert(entry[oldSize]==INVALID);
		entry[oldSize]=v;
		sizes[loc]++;
		return 1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Resizing           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Resizes the entries and sizes arrays to the specified capacity.
	 * @param size2 New capacity (must be greater than current size) */
	public final void resize(int size2){
		assert(size2>size);
		entries=KillSwitch.copyOf(entries, size2);//TODO: Safe copy to prevent memory exception
		sizes=KillSwitch.copyOf(sizes, size2);
	}
	
	/** Compresses the list by eliminating unused trailing space in the arrays. */
	public final void shrink(){
		if(size==entries.length){return;}
		entries=KillSwitch.copyOf(entries, size);
		sizes=KillSwitch.copyOf(sizes, size);
	}
	
	/** Makes each entry unique and minimal-length by removing duplicates.
	 * Currently untested and not expected to be used in production. */
	public final void shrinkToUnique(){
		assert(false) : "TODO: This function has not been tested and is not expected to be used.";
		assert(!shrunk);
		assert(mode==DISORDERED);//Not really necessary
		for(int i=0; i<size; i++){
			if(sizes[i]>slowAddLimit){//Under this limit everything is already unique 
				shrinkToUnique(i);
			}
		}
		shrunk=true;
		assert(readyToUse());
	}
	
	/**
	 * Removes duplicates from a specific entry by sorting and condensing.
	 * Creates a new array if duplicates are found.
	 * @param loc The location of the entry to make unique
	 */
	private void shrinkToUnique(int loc){
		final long[] entry=entries[loc];
		final int oldSize=sizes[loc];
		assert(oldSize>1);
		Arrays.sort(entry, 0, oldSize);
		int unique=0;
		long prev=INVALID;
		for(int i=0; i<oldSize; i++){
			assert(entry[i]>=0);
			long current=entry[i];
			if(current!=prev){
				unique++;
			}
			prev=current;
		}
		if(unique==oldSize){return;} //No need to condense aside from saving a little RAM
		
		long[] set2=new long[unique];
		unique=0;
		prev=INVALID;
		for(int i=0; i<oldSize; i++){
			long current=entry[i];
			if(current!=prev){
				set2[unique]=current;
				unique++;
			}
			prev=current;
		}
		set(loc, set2, unique);
	}
	
	private boolean readyToUse(){
		return shrunk || mode==ASCENDING || mode==UNIQUE;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final long[] get(int loc){
		return(loc>=size ? null : entries[loc]);
	}
	
	/**
	 * Returns the length of the entry at the specified location.
	 * Added for better IntList2 compatibility.
	 * @param i The location to check
	 * @return The length of the entry, or 0 if out of bounds
	 */
	public int getLen(int i) {
		return i>=size ? 0 : sizes[i];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Returns a string representation showing all non-null entries with their indices.
	 */
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<size; i++){
			if(entries[i]!=null){//Could be improved to use sizes
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
	
	/** Holds entries where each entry is a set of numbers in a long array.
	 * Leftmost values are valid, rightmost values are invalid. */
	private long[][] entries;
	private int[] sizes;
	public int size=0;
	
	private boolean shrunk=false;
	
	/** Preconditions for adding values (ASCENDING, UNIQUE, or DISORDERED). */
	private final int mode;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Marker value for invalid/unused array positions. */
	public static final int INVALID=-1;
	
	/**
	 * Only scan back this far for duplicates when adding values to optimize performance.
	 */
	public static final int slowAddLimit=4;
	
	/*--------------------------------------------------------------*/
	/*----------------            Modes             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Mode requiring all adds to an entry must be nondescending. */
	public static final int ASCENDING=1;
	/** Mode requiring all adds to an entry must be unique. */
	public static final int UNIQUE=2;
	/** Mode with no requirements for adds. To ensure set functionality,
	 * shrinkToUnique should be called before use. */
	public static final int DISORDERED=3;
	
	/** Default mode for new instances. Should be set prior to creation. */
	public static int defaultMode=ASCENDING;
	/** Default initial capacity (in entry slots) for new LongList3 instances. */
	public static int defaultInitialSize=256;
	
}
