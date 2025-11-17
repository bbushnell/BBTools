package structures;

/**
 * A mutable long object
 * @author Brian Bushnell
 * @date Feb 8, 2013
 *
 */
public class LongM implements Comparable<LongM> {
	/** Creates a new mutable LongM with value 0 */
	public LongM(){this(0L);}
	/** Creates a new mutable LongM with the specified value.
	 * @param v Initial value for the long wrapper */
	public LongM(long v){value=v;}

	/**
	 * @param v Value
	 * @param mut Mutable
	 */
	public LongM(long v, boolean mut) {
		value=v;
		mutable=mut;
	}
	
	/**
	 * Creates an immutable copy of this LongM.
	 * If already immutable, returns this instance to avoid unnecessary allocation.
	 * @return Immutable copy with the same value
	 */
	public LongM iCopy(){
		if(!mutable){return this;}
		return new LongM(value, false);
	}
	
	/** Gets the current long value */
	public long value(){return value;}
//	public long longValue(){return value;}
	/** Permanently locks this LongM to prevent further mutations */
	public void lock(){mutable=false;}
	
	/**
	 * Sets the value if this LongM is mutable.
	 * @param v New value to set
	 * @return The new value
	 * @throws RuntimeException if this LongM has been locked
	 */
	public long set(long v){
		if(!mutable){throw new RuntimeException("Mutating a locked LongM");}
		return (value=v);
	}
	/**
	 * Increments the value by 1 if mutable.
	 * @return The new incremented value
	 * @throws RuntimeException if this LongM has been locked
	 */
	public long increment(){return set(value+1);}
	/**
	 * Increments the value by the specified amount if mutable.
	 * @param x Amount to add to current value
	 * @return The new incremented value
	 * @throws RuntimeException if this LongM has been locked
	 */
	public long increment(long x){return set(value+x);}
	
	@Override
	public int hashCode(){
		return (int)((value^(value>>>32))&0xFFFFFFFFL);
	}
	
	@Override
	public int compareTo(LongM b){
		return value==b.value ? 0 : value<b.value ? -1 : 1;
	}
	
	/**
	 * Tests equality with another LongM based on long values.
	 * @param b LongM to compare against
	 * @return true if both have the same long value
	 */
	public boolean equals(LongM b){
		return value==b.value;
	}
	
	@Override
	public boolean equals(Object b){
		return equals((LongM)b); //Possible bug: Unchecked cast may throw ClassCastException if b is not LongM
	}
	@Override
	public String toString(){return Long.toString(value);}
	/** Returns hexadecimal string representation of the long value */
	public String toHexString(){return Long.toHexString(value);}
	/** Returns binary string representation of the long value */
	public String toBinaryString(){return Long.toBinaryString(value);}
	
	/** Controls whether this LongM can be modified after construction */
	private boolean mutable=true;
	/** The wrapped long value */
	private long value;
}
