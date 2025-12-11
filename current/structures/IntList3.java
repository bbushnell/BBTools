package structures;

import java.util.Arrays;

import shared.KillSwitch;


/**
 * Dynamic list of integer arrays with independent size tracking for each array.
 * More efficient than a list of IntList objects, designed for use with HashArrayHybridFast.
 * Assumes each entry is a set with additions either ascending (Seal) or unique (Sketch).
 * Does not extend IntList2 to avoid virtual function overhead.
 *
 * @author Brian Bushnell
 * @date Dec 10, 2018
 */
public final class IntList3{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates an IntList3 with default initial size and mode */
	public IntList3(){this(defaultInitialSize, defaultMode);}
	
	/**
	 * Creates an IntList3 with specified initial size and mode.
	 * @param initialSize Initial capacity for the entries array
	 * @param mode_ Addition mode (ASCENDING or UNIQUE)
	 */
	public IntList3(int initialSize, int mode_){
		assert(initialSize>0);
		entries=new int[initialSize][];
		sizes=new int[initialSize];
		
		mode=mode_;
		assert(mode==ASCENDING || mode==UNIQUE) : "Unsupported mode: "+mode;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds an entry to the end of the list.
	 * @param entry Integer array to add
	 * @param len Number of valid elements in the array
	 */
	public final void add(int[] entry, int len){
		set(size, entry, len);
	}
	
	/**
	 * Sets the entry at the specified location.
	 * Resizes the list if necessary to accommodate the location.
	 *
	 * @param loc Target location index
	 * @param entry Integer array to store
	 * @param len Number of valid elements in the array
	 */
	public final void set(int loc, int[] entry, int len){
		assert((entry==null && len==0) || (entry.length>0 && len<=entry.length)) : len+", "+(entry==null ? entry : ""+entry.length);
		
		if(loc>=entries.length){//Resize by doubling if necessary
			resize(max(size*2, loc+1));
		}
		entries[loc]=entry;
		sizes[loc]=len;
		size=max(size, loc+1);
	}
	
	/**
	 * Inserts a value into the list at the specified location.
	 * Creates a new entry if location is empty, otherwise adds to existing entry.
	 * Scans recent entries to avoid duplicates based on mode assumptions.
	 *
	 * @param v Value to insert
	 * @param loc Target location index
	 * @return 1 if value was added, 0 if already present
	 */
	public final int insertIntoList(final int v, final int loc){
		
		//If the location is empty
		if(loc>=size || entries[loc]==null){
			//Add a new entry and return
			set(loc, new int[] {v, INVALID}, 1);
			return 1;
		}
			
		int[] entry=get(loc);
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
	
	/** Resizes the entries and sizes arrays to accommodate more elements.
	 * @param size2 New capacity (must be larger than current size) */
	public final void resize(int size2){
		assert(size2>size);
		entries=KillSwitch.copyOf(entries, size2);//TODO: Safe copy to prevent memory exception
		sizes=KillSwitch.copyOf(sizes, size2);
	}
	
	/** Compresses the list by eliminating unused trailing space */
	public final void shrink(){
		if(size==entries.length){return;}
		entries=KillSwitch.copyOf(entries, size);
		sizes=KillSwitch.copyOf(sizes, size);
	}
	
	/** Makes each entry unique and minimal-length for disordered additions.
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
	
	/** Makes the entry at the specified location unique by sorting and deduplicating.
	 * @param loc Location index to process */
	private void shrinkToUnique(int loc){
		final int[] entry=entries[loc];
		final int oldSize=sizes[loc];
		assert(oldSize>1);
		Arrays.sort(entry, 0, oldSize);
		int unique=0;
		int prev=INVALID;
		for(int i=0; i<oldSize; i++){
			assert(entry[i]>=0);
			int current=entry[i];
			if(current!=prev){
				unique++;
			}
			prev=current;
		}
		if(unique==oldSize){return;} //No need to condense aside from saving a little RAM
		
		int[] set2=new int[unique];
		unique=0;
		prev=INVALID;
		for(int i=0; i<oldSize; i++){
			int current=entry[i];
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
	
	public final int[] get(int loc){
		return(loc>=size ? null : entries[loc]);
	}
	
	/**
	 * Gets the length of the entry at the specified location.
	 * Added for better IntList2 compatibility.
	 * @param i Location index
	 * @return Number of valid elements in the entry, or 0 if out of bounds
	 */
	public int getLen(int i) {
		return i>=size ? 0 : sizes[i];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns a string representation showing all non-null entries with their indices.
	 * @return String representation of the list */
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
	
	/**
	 * Holds entries as integer arrays; leftmost values are valid, rightmost invalid
	 */
	private int[][] entries;
	/** Number of valid values in each entry */
	private int[] sizes;
	/** Number of entries in the primary array */
	public int size=0;
	
	/** True after shrinkToUnique has been called; currently unused */
	private boolean shrunk=false;
	
	/** Preconditions for adding values (ASCENDING, UNIQUE, or DISORDERED) */
	private final int mode;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Sentinel value indicating unused array positions */
	public static final int INVALID=-1;
	
	/** Maximum number of recent entries to scan for duplicates when adding */
	public static final int slowAddLimit=4;
	
	/*--------------------------------------------------------------*/
	/*----------------            Modes             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Mode requiring all additions to an entry be nondescending */
	public static final int ASCENDING=1;
	/** Mode requiring all additions to an entry be unique */
	public static final int UNIQUE=2;
	/**
	 * Mode with no ordering requirements; shrinkToUnique should be called before use
	 */
	public static final int DISORDERED=3;
	
	/** Default mode for new IntList3 instances; should be set by Seal or Sketch */
	public static int defaultMode=ASCENDING;
	/** Default initial capacity (in entry slots) for new IntList3 instances */
	public static int defaultInitialSize=256;
	
}
