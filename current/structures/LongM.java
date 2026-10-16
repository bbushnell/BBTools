package structures;

/**
 * Mutable long wrapper for code paths that need to replace or increment a
 * numeric value through a reference.
 * <p>
 * Instances start mutable unless constructed otherwise. lock() prevents later
 * changes, and iCopy() returns an immutable snapshot when the receiver is still
 * mutable.
 * Mutation is unsynchronized; callers sharing an instance coordinate externally.
 * The mutable value is the equality, hash and sort key, so callers must not
 * mutate instances while they are in hash-based or sorted collections.
 *
 * @author Brian Bushnell
 * @date Feb 8, 2013
 */
public class LongM implements Comparable<LongM>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a mutable wrapper with value 0. */
	public LongM(){this(0L);}

	/** Creates a mutable wrapper with the specified value. */
	public LongM(long v){value=v;}

	/**
	 * Creates a wrapper with the specified value and mutability.
	 * @param v Initial value.
	 * @param mut True if future mutations are allowed.
	 */
	public LongM(long v, boolean mut){
		value=v;
		mutable=mut;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Mutation           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns an immutable snapshot of this value.
	 * Immutable receivers are already snapshots and return themselves.
	 */
	public LongM iCopy(){
		if(!mutable){return this;}
		return new LongM(value, false);
	}

	/** Returns the current wrapped value. */
	public long value(){return value;}
//	public long longValue(){return value;}

	/** Prevents future changes to this wrapper. */
	public void lock(){mutable=false;}

	/**
	 * Replaces and returns the wrapped value.
	 * @throws RuntimeException if this wrapper has been locked.
	 */
	public long set(long v){
		if(!mutable){throw new RuntimeException("Mutating a locked LongM");}
		return (value=v);
	}

	/** Increments the value by 1 and returns the updated value. */
	public long increment(){return set(value+1);}

	/** Adds x to the value and returns the updated value. */
	public long increment(long x){return set(value+x);}

	/*--------------------------------------------------------------*/
	/*----------------           Equality           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int hashCode(){
		return (int)((value^(value>>>32))&0xFFFFFFFFL);
	}

	/** Compares values in ascending order. */
	@Override
	public int compareTo(LongM b){
		return value==b.value ? 0 : value<b.value ? -1 : 1;
	}

	/**
	 * Tests whether b has the same wrapped value.
	 * @param b LongM to compare; may be null.
	 */
	public boolean equals(LongM b){
		return b!=null && value==b.value;
	}

	@Override
	public boolean equals(Object b){
		return b instanceof LongM && equals((LongM)b);
	}

	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public String toString(){return Long.toString(value);}

	/** Returns the value formatted in unsigned hexadecimal. */
	public String toHexString(){return Long.toHexString(value);}

	/** Returns the value formatted in unsigned binary. */
	public String toBinaryString(){return Long.toBinaryString(value);}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** True if set() is allowed to replace the wrapped value. */
	private boolean mutable=true;

	/** Wrapped value. */
	private long value;
}
