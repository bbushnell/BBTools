package structures;

/**
 * Abstract base class for bit manipulation structures with variable bits per element.
 * Provides factory methods and common operations for bit set implementations.
 * Supports both single-bit (RawBitSet) and multi-bit (MultiBitSet) storage modes.
 * @author Brian Bushnell
 */
public abstract class AbstractBitSet {
	
	/**
	 * Factory method to create appropriate bit set implementation based on bits per element.
	 * Creates RawBitSet for 1-bit elements or MultiBitSet for 2-bit elements.
	 *
	 * @param elements Number of elements to store
	 * @param bitsPerElement Bits per element (must be 1 or 2)
	 * @return RawBitSet for 1-bit elements, MultiBitSet for 2-bit elements
	 * @throws RuntimeException if bitsPerElement is not 1 or 2
	 */
	public static AbstractBitSet make(int elements, int bitsPerElement){
		assert(bitsPerElement==1 || bitsPerElement==2) : bitsPerElement;
		assert(Integer.bitCount(bitsPerElement)==1) : bitsPerElement;
		assert(Integer.bitCount(1+Integer.numberOfTrailingZeros(bitsPerElement))==1) : bitsPerElement;
		//Can also assert
		if(bitsPerElement==1){
			return new RawBitSet(elements);
		}else if(bitsPerElement==2){
			return new MultiBitSet(elements);
		}else{
			throw new RuntimeException(""+bitsPerElement);
		}
	}
	
//	public final void set(int x){increment(x);}
	/** Increments the count at position x by 1.
	 * @param x Position to increment */
	public final void increment(int x){increment(x, 1);}
	/**
	 * Increments the count at position x by the specified amount.
	 * @param x Position to increment
	 * @param incr Amount to increment by
	 */
	public abstract void increment(int x, int incr);
	
//	public final boolean get(int x){return getCount(x)>0;}
	/**
	 * Returns the count value at position x.
	 * @param x Position to query
	 * @return Count value at position x
	 */
	public abstract int getCount(int x);
	
	/** Clears the input BitSet */
	public final void add(AbstractBitSet bs){
		if(bs.getClass()==RawBitSet.class){add((RawBitSet)bs);}
		else if(bs.getClass()==MultiBitSet.class){add((MultiBitSet)bs);}
		else{throw new RuntimeException("Bad class: "+bs.getClass());}
	}
	
	/** Clears the input BitSet */
	public final void add(RawBitSet bs){
		assert(this.getClass()==bs.getClass()) : this.getClass()+", "+bs.getClass();
		RawBitSet bs2=(RawBitSet)this;
		assert(capacity()==bs.capacity()) : capacity()+", "+bs.capacity();
		final int[] rbsArray=bs.array();
		final int[] rbs2Array=bs2.array();
		final int rbsLength=bs.length();
		for(int i=0; i<rbsLength; i++){
			final int value=rbsArray[i];
//			if(value!=0){bs2.addToCell(i, value);}
			rbs2Array[i]|=value;
			rbsArray[i]=0;
		}
	}
	
	/** Clears the input BitSet */
	public final void add(MultiBitSet bs){
		assert(this.getClass()==bs.getClass()) : this.getClass()+", "+bs.getClass();
		MultiBitSet bs2=(MultiBitSet)this;
		assert(bits()==bs.bits());
		assert(capacity()==bs.capacity()) : capacity()+", "+bs.capacity();
		final int[] rbsArray=bs.array();
		final int rbsLength=bs.length();
		for(int i=0; i<rbsLength; i++){
			final int value=rbsArray[i];
			if(value!=0){bs2.addToCell(i, value);}
			rbsArray[i]=0;
		}
	}
	
	/**
	 * Sets each position to the maximum of this BitSet and the input BitSet.
	 * Dispatches to type-specific setToMax methods based on input type.
	 * @param bs BitSet to compare against for maximum values
	 * @throws RuntimeException if bs is not a recognized BitSet type
	 */
	public final void setToMax(AbstractBitSet bs){
		if(bs.getClass()==RawBitSet.class){setToMax((RawBitSet)bs);}
		else if(bs.getClass()==MultiBitSet.class){setToMax((MultiBitSet)bs);}
		else{throw new RuntimeException("Bad class: "+bs.getClass());}
	}
	
	/**
	 * Sets each position to the maximum of this BitSet and the RawBitSet.
	 * For RawBitSet this is equivalent to the add operation.
	 * @param bs RawBitSet to compare against for maximum values
	 */
	public void setToMax(RawBitSet bs) {
		add(bs);
	}
	
	/**
	 * Sets each position to the maximum of this BitSet and the MultiBitSet.
	 * Both BitSets must be the same type, have same bits per element, and same capacity.
	 * @param bs MultiBitSet to compare against for maximum values
	 */
	public void setToMax(MultiBitSet bs) {
		assert(this.getClass()==bs.getClass()) : this.getClass()+", "+bs.getClass();
		assert(bits()==bs.bits());
		assert(capacity()==bs.capacity()) : capacity()+", "+bs.capacity();
		final int[] rbsArray=bs.array();
		final int rbsLength=bs.length();
		for(int i=0; i<rbsLength; i++){
			final int value=rbsArray[i];
			if(value!=0){setToMax(i, value);}
		}
	}

	public abstract void addToCell(final int cell, final int mask);
	public abstract void setToMax(final int cell, final int mask);
	
	/** Clears all values in the BitSet, resetting all positions to zero */
	public abstract void clear();
	/**
	 * Sets the capacity of the BitSet with optional extra space.
	 * @param capacity Target capacity in elements
	 * @param extra Additional elements to allocate beyond capacity
	 */
	public abstract void setCapacity(long capacity, int extra);
	/** Returns the number of non-zero positions in the BitSet.
	 * @return Count of positions with non-zero values */
	public abstract long cardinality();
	/** Returns the maximum number of elements this BitSet can store.
	 * @return Maximum element capacity */
	public abstract long capacity();
	/** Returns the length of the underlying storage array.
	 * @return Storage array length */
	public abstract int length();
	/** Returns the number of bits used per element.
	 * @return Bits per element (typically 1 or 2) */
	public abstract int bits(); //per element
	
	@Override
	public final String toString(){
		
		StringBuilder sb=new StringBuilder();
		
		final long cap=capacity();
		String spacer="";
		sb.append("{");
		for(long i=0; i<cap; i++){
			int x=getCount((int)i);
			if(x>0){
				sb.append(spacer);
				sb.append("("+i+","+x+")");
				spacer=", ";
			}
		}
		sb.append("}");
		
		return sb.toString();
	}
	
//	public final RawBitSet toRaw(){
//		if(this.getClass()==RawBitSet.class){return (RawBitSet)this;}
//		final int cap=(int)capacity();
//		RawBitSet rbs=new RawBitSet(cap, 0);
//		for(int i=0; i<cap; i++){
//			if(get(i)){rbs.set(i);}
//		}
//		return rbs;
//	}
	
}
