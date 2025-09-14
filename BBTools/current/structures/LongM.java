package structures;

/**
 * A mutable long object
 * @author Brian Bushnell
 * @date Feb 8, 2013
 *
 */
public class LongM implements Comparable<LongM> {
	public LongM(){this(0L);}
	public LongM(long v){value=v;}

	/**
	 * @param v Value
	 * @param mut Mutable
	 */
	public LongM(long v, boolean mut) {
		value=v;
		mutable=mut;
	}
	
	public LongM iCopy(){
		if(!mutable){return this;}
		return new LongM(value, false);
	}
	
	public long value(){return value;}
//	public long longValue(){return value;}
	public void lock(){mutable=false;}
	
	public long set(long v){
		if(!mutable){throw new RuntimeException("Mutating a locked LongM");}
		return (value=v);
	}
	public long increment(){return set(value+1);}
	public long increment(long x){return set(value+x);}
	
	@Override
	public int hashCode(){
		return (int)((value^(value>>>32))&0xFFFFFFFFL);
	}
	
	@Override
	public int compareTo(LongM b){
		return value==b.value ? 0 : value<b.value ? -1 : 1;
	}
	
	public boolean equals(LongM b){
		return value==b.value;
	}
	
	@Override
	public boolean equals(Object b){
		return equals((LongM)b);
	}
	@Override
	public String toString(){return Long.toString(value);}
	public String toHexString(){return Long.toHexString(value);}
	public String toBinaryString(){return Long.toBinaryString(value);}
	
	private boolean mutable=true;
	private long value;
}