package structures;

import java.util.Arrays;

import shared.KillSwitch;


/** 
 * Like IntList2 but has a size for each array independent of its length.
 * This is similar to a list of IntList but more efficient.
 * Designed for use with HashArrayHybridFast.
 * Assumes each entry is a set with all adds either ascending (Seal)
 * or unique (Sketch).
 * 
 * IntList3 does not extend IntList2 because virtual functions might be slower,
 * but it would make the code more concise.
 * 
 * @author Brian Bushnell
 * @date Dec 10, 2018
 *
 */
public final class IntList3{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Re-call with default initial size and mode. */
	public IntList3(){this(defaultInitialSize, defaultMode);}
	
	/** Construct an IntList3 with this initial size.
	 * Setting the proper mode is the caller's responsibility. */
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
	
	/** Add this entry to the end of the list */
	public final void add(int[] entry, int len){
		set(size, entry, len);
	}
	
	/** Set this location to specified entry */
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
	 * Add this value to the specified location. 
	 * If an entry exists, insert the value, enlarging if necessary.
	 * Otherwise, create a new entry. */
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
	
	/** Resize the array of entries */
	public final void resize(int size2){
		assert(size2>size);
		entries=KillSwitch.copyOf(entries, size2);//TODO: Safe copy to prevent memory exception
		sizes=KillSwitch.copyOf(sizes, size2);
	}
	
	/** Compress the list by eliminating the unused trailing space */
	public final void shrink(){
		if(size==entries.length){return;}
		entries=KillSwitch.copyOf(entries, size);
		sizes=KillSwitch.copyOf(sizes, size);
	}
	
	/**
	 * Provided for supporting disordered additions; not currently used.  
	 * Make each entry unique and minimal-length. 
	 * */
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
	 * Provided for supporting disordered additions; not currently used. 
	 * I could redo this so that the sets become packed and the size changes without creating a new array. 
	 * */
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
	
	/** Get the entry at the location */
	public final int[] get(int loc){
		return(loc>=size ? null : entries[loc]);
	}
	
	/** Added for better IntList2 compatibility */
	public int getLen(int i) {
		return i>=size ? 0 : sizes[i];
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
	
	/** Holds entries.  Each entry is a sets of numbers in an int[]. 
	 * Leftmost values are valid, rightmost values are invalid. */
	private int[][] entries;
	/** Number of values in each entry. */
	private int[] sizes;
	/** Number of entries in the primary array. */
	public int size=0;
	
	/** True after shrinkToUnique has been called.
	 * Currently unused. */
	private boolean shrunk=false;
	
	/** Preconditions for adding values */
	private final int mode;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This could be made mutable per instance, but it would be a lot of work. */
	public static final int INVALID=-1;
	
	/** Only scan back this far for duplicates when adding values */
	public static final int slowAddLimit=4;
	
	/*--------------------------------------------------------------*/
	/*----------------            Modes             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** All adds to an entry must be nondescending */
	public static final int ASCENDING=1;
	/** All adds to an entry must be unique */
	public static final int UNIQUE=2;
	/** No requirements for adds.
	 * To ensure set functionality, shrinkToUnique should be called before use. */
	public static final int DISORDERED=3;
	
	/** Should be set prior to creation, e.g. by Seal or Sketch */
	public static int defaultMode=ASCENDING;
	public static int defaultInitialSize=256;
	
}
