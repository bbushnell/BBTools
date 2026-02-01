package map;

/**
 * Abstract base for bit sets with 1- or 2-bit elements.
 * Provides factories and shared operations for RawBitSet and MultiBitSet implementations.
 * @author Brian Bushnell
 */
public abstract class AbstractBitSet {
	
	/**
	 * Factory that creates a RawBitSet (1-bit) or MultiBitSet (2-bit) based on bitsPerElement.
	 *
	 * @param elements Number of elements to store
	 * @param bitsPerElement Bits per element (1 or 2)
	 * @return Concrete AbstractBitSet implementation
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
	public final void increment(int x){increment(x, 1);}
	/**
	 * Increments the count at position x by a specified amount.
	 * @param x Position to increment
	 * @param incr Amount to add
	 */
	public abstract void increment(int x, int incr);
	
//	public final boolean get(int x){return getCount(x)>0;}
	/**
	 * Returns the count value at position x.
	 * @param x Position to query
	 * @return Count at the position
	 */
	public abstract int getCount(int x);
	
	/**
	 * Adds values from another AbstractBitSet of the same size, dispatching by implementation.
	 * Clears the input set after addition.
	 * @param bs BitSet to add from (will be cleared)
	 * @throws RuntimeException if the class is unsupported
	 */
	public final void add(AbstractBitSet bs){
		if(bs.getClass()==RawBitSet.class){add((RawBitSet)bs);}
		else if(bs.getClass()==MultiBitSet.class){add((MultiBitSet)bs);}
		else{throw new RuntimeException("Bad class: "+bs.getClass());}
	}
	
	/** Bitwise-ORs values from a RawBitSet of the same capacity into this set, then clears the input.
	 * @param bs RawBitSet to add and clear */
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
	
	/** Adds cell values from a MultiBitSet of matching shape into this set, then clears the input.
	 * @param bs MultiBitSet to add and clear */
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
	
	/** Sets each position to the maximum of this set and the input, dispatching by implementation type.
	 * @param bs BitSet to compare against */
	public final void setToMax(AbstractBitSet bs){
		if(bs.getClass()==RawBitSet.class){setToMax((RawBitSet)bs);}
		else if(bs.getClass()==MultiBitSet.class){setToMax((MultiBitSet)bs);}
		else{throw new RuntimeException("Bad class: "+bs.getClass());}
	}
	
	/** Sets each position to the maximum of this set and a RawBitSet (equivalent to add()).
	 * @param bs RawBitSet to merge */
	public void setToMax(RawBitSet bs) {
		add(bs);
	}
	
	/** Sets each position to the maximum of this set and a MultiBitSet of matching shape.
	 * @param bs MultiBitSet to merge */
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

	/**
	 * Adds a masked value directly to a storage cell.
	 * @param cell Cell index
	 * @param mask Masked value to add
	 */
	public abstract void addToCell(final int cell, final int mask);
	/**
	 * Sets a storage cell to the maximum of its current value and the masked input.
	 * @param cell Cell index
	 * @param mask Masked value to compare
	 */
	public abstract void setToMax(final int cell, final int mask);
	
	/** Clears all values, resetting the bit set to zero. */
	public abstract void clear();
	/**
	 * Sets the capacity with optional extra space, reallocating if necessary.
	 * @param capacity Target capacity in elements
	 * @param extra Additional elements to allocate
	 */
	public abstract void setCapacity(long capacity, int extra);
	/** Returns the number of positions with non-zero values.
	 * @return Count of non-zero elements */
	public abstract long cardinality();
	/** Returns the maximum number of elements this bit set can store.
	 * @return Element capacity */
	public abstract long capacity();
	public abstract int length();
	/** Returns the number of bits used per element (1 or 2).
	 * @return Bits per element */
	public abstract int bits(); //per element
	
	/** Returns a string of non-zero positions and counts in the form {(pos,count), ...}.
	 * @return String representation of non-zero elements */
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
