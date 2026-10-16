package structures;

import java.util.Arrays;

import shared.KillSwitch;


/**
 * Dynamic list of live long arrays with independent valid-length tracking.
 * Tracks only the leftmost valid values in each entry; spare capacity should
 * contain INVALID when callers grow an entry through insertIntoList().
 * Assumes each entry is a set with additions either ascending or unique,
 * because only the last few values are scanned for duplicates.
 * Access to mutation and the live arrays returned by get() must be coordinated
 * by the caller.
 * @author Brian Bushnell
 * @date Dec 10, 2018
 */
public final class LongList3{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a LongList3 with the current default initial capacity and mode. */
	public LongList3(){this(defaultInitialSize, defaultMode);}

	/**
	 * Creates a LongList3 with the specified initial entry capacity and addition mode.
	 * @param initialSize Initial capacity of the entry table.
	 * @param mode_ Addition mode, ASCENDING or UNIQUE with assertions enabled.
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
	 * Appends a live entry array at the next location.
	 * @param entry Long array to store.
	 * @param len Number of valid elements in entry.
	 */
	public final void add(long[] entry, int len){
		set(size, entry, len);
	}

	/**
	 * Stores a live entry array and valid length at loc, growing the entry table
	 * when loc is beyond the current capacity.
	 * @param loc Target location index.
	 * @param entry Long array to store.
	 * @param len Number of valid elements in entry.
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
	 * Inserts a value into the entry at loc, creating that entry if needed.
	 * Only the newest slowAddLimit values are scanned; older values are assumed
	 * not to duplicate v by the selected addition mode.
	 * @param v Value to insert.
	 * @param loc Target location index.
	 * @return 1 if value was added, or 0 if found in the scanned suffix.
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

	/**
	 * Resizes the entries and sizes arrays to accommodate more elements.
	 * @param size2 New capacity; must exceed the current logical size.
	 */
	public final void resize(int size2){
		assert(size2>size);
		entries=KillSwitch.copyOf(entries, size2);//TODO: Safe copy to prevent memory exception
		sizes=KillSwitch.copyOf(sizes, size2);
	}

	/** Removes unused trailing slots from the entry table. */
	public final void shrink(){
		if(size==entries.length){return;}
		entries=KillSwitch.copyOf(entries, size);
		sizes=KillSwitch.copyOf(sizes, size);
	}

	/**
	 * Untested DISORDERED-mode helper, unconditionally blocked by assertion.
	 * If assertions are disabled, entries above slowAddLimit are sorted and
	 * deduplicated; short entries and already-unique backing arrays are not
	 * reallocated.
	 */
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
	 * Makes one entry unique by sorting its valid prefix and deduplicating it.
	 * @param loc Location index to process.
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

	/**
	 * Gets the live entry array at loc.
	 * @return Entry array, or null if loc is past the logical table extent.
	 */
	public final long[] get(int loc){
		return(loc>=size ? null : entries[loc]);
	}

	/**
	 * Gets the valid length of one entry.
	 * @param i Location index.
	 * @return Number of valid elements in the entry, or 0 if i is past the
	 * logical table extent.
	 */
	public int getLen(int i){
		return i>=size ? 0 : sizes[i];
	}

	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns a string representation showing all non-null entries with indices.
	 * The complete backing arrays are printed, including INVALID spare slots.
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

	/** Live entry arrays; values below sizes[i] are valid. */
	private long[][] entries;

	/** Number of valid values in each entry. */
	private int[] sizes;

	/** Logical table extent, including null or empty slots below size. */
	public int size=0;

	/** Set by shrinkToUnique() and read by readyToUse() assertion checks. */
	private boolean shrunk=false;

	/** Preconditions for adding values. */
	private final int mode;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Sentinel value indicating unused array positions. */
	public static final int INVALID=-1;

	/** Maximum number of recent entries to scan for duplicates when adding. */
	public static final int slowAddLimit=4;

	/*--------------------------------------------------------------*/
	/*----------------            Modes             ----------------*/
	/*--------------------------------------------------------------*/

	/** Mode requiring all additions to an entry be nondescending. */
	public static final int ASCENDING=1;

	/** Mode requiring all additions to an entry be unique. */
	public static final int UNIQUE=2;

	/**
	 * Mode with no ordering requirements; shrinkToUnique should be called before use.
	 */
	public static final int DISORDERED=3;

	/** Default mode for new LongList3 instances. */
	public static int defaultMode=ASCENDING;

	/** Default initial capacity, in entry slots, for new LongList3 instances. */
	public static int defaultInitialSize=256;

}
