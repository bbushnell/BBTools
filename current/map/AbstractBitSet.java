package map;

/**
 * Packed, saturating counters with one or two bits per logical element.
 * Element indices and storage-cell indices are different: a cell packs 32 one-bit
 * elements or 16 two-bit elements. Cell operations must leave padding beyond the
 * logical capacity zero if cardinality is to count only logical elements.
 * Whole-set merges support matching RawBitSet or MultiBitSet instances, require
 * exclusive access, and do not support AtomicBitSet. Source ownership differs
 * between addition and maximum merging; see the individual overloads.
 * Capacity changes are not preserving resizes: replacement storage starts empty,
 * while reused storage retains its contents, including temporarily hidden cells.
 * @author Brian Bushnell
 */
public abstract class AbstractBitSet{

	/** Initializes the stateless base; concrete implementations allocate their own storage. */
	public AbstractBitSet(){}

	/*--------------------------------------------------------------*/
	/*----------------           Factory            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Factory that creates a RawBitSet (1-bit) or MultiBitSet (2-bit) based on bitsPerElement.
	 *
	 * @param elements Nonnegative logical element capacity
	 * @param bitsPerElement Bits per element (1 or 2)
	 * @return Concrete AbstractBitSet implementation
	 * @throws AssertionError If assertions are enabled and the width is unsupported
	 * @throws RuntimeException If the width is unsupported with assertions disabled
	 */
	public static AbstractBitSet make(int elements, int bitsPerElement){
		assert(bitsPerElement==1 || bitsPerElement==2) : bitsPerElement;
		//Retained guards from the intended {1,2,4,8} generalization; redundant for {1,2}.
		assert(Integer.bitCount(bitsPerElement)==1) : bitsPerElement;
		assert(Integer.bitCount(1+Integer.numberOfTrailingZeros(bitsPerElement))==1) : bitsPerElement;
		if(bitsPerElement==1){
			return new RawBitSet(elements);
		}else if(bitsPerElement==2){
			return new MultiBitSet(elements);
		}else{
			throw new RuntimeException(""+bitsPerElement);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Element Access         ----------------*/
	/*--------------------------------------------------------------*/

	/** Increments one logical element, saturating at the implementation's maximum.
	 * @param x Logical index in [0, capacity()), representable as a nonnegative int */
	public final void increment(int x){increment(x, 1);}
	/**
	 * Increments the count at position x by a specified amount.
	 * Changes only the selected element; one-bit implementations set it to one.
	 * @param x Logical index in [0, capacity()), representable as a nonnegative int
	 * @param incr Positive amount to add; one-bit implementations simply set the bit
	 */
	public abstract void increment(int x, int incr);

	/**
	 * Returns the count value at position x.
	 * @param x Logical index in [0, capacity()), representable as a nonnegative int
	 * @return Count in [0, (1&lt;&lt;bits())-1]
	 */
	public abstract int getCount(int x);

	/*--------------------------------------------------------------*/
	/*----------------        Whole-set Merges       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds counts from a set of the same concrete type and logical capacity, with saturation.
	 * Clears the distinct input set after addition. Source and destination must not
	 * be the same object, and callers must provide exclusive access to both sets.
	 * @param bs Non-null, distinct source; its active storage cells are cleared
	 * @throws RuntimeException if the class is unsupported
	 */
	public final void add(AbstractBitSet bs){
		if(bs.getClass()==RawBitSet.class){add((RawBitSet)bs);}
		else if(bs.getClass()==MultiBitSet.class){add((MultiBitSet)bs);}
		else{throw new RuntimeException("Bad class: "+bs.getClass());}
	}

	/** ORs a distinct RawBitSet into this RawBitSet, then clears its active cells.
	 * Both sets must have equal logical capacity and be exclusively accessible.
	 * @param bs Non-null source of matching capacity; consumed by this merge */
	public final void add(RawBitSet bs){
		assert(this.getClass()==bs.getClass()) : this.getClass()+", "+bs.getClass();
		RawBitSet bs2=(RawBitSet)this;
		assert(capacity()==bs.capacity()) : capacity()+", "+bs.capacity();
		final int[] rbsArray=bs.array();
		final int[] rbs2Array=bs2.array();
		final int rbsLength=bs.length();
		//1-bit elements pack 32/word, so set-union == bitwise OR of the words (no per-element saturation
		//needed, unlike MultiBitSet's addToCell). Each source word is OR'd into this, then zeroed (move-merge).
		for(int i=0; i<rbsLength; i++){
			final int value=rbsArray[i];
			rbs2Array[i]|=value;
			rbsArray[i]=0;
		}
	}

	/** Adds a distinct MultiBitSet into this MultiBitSet, saturating each count at three.
	 * Both sets must have equal logical capacity and be exclusively accessible.
	 * @param bs Non-null source of matching capacity; its active cells are cleared */
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

	/** Merges elementwise maxima from a set of matching concrete type and capacity.
	 * Requires exclusive access. A RawBitSet source is cleared and must be distinct;
	 * a MultiBitSet source is preserved and may be this set itself.
	 * @param bs Non-null RawBitSet or MultiBitSet of the same type and capacity
	 * @throws RuntimeException If the source type is unsupported */
	public final void setToMax(AbstractBitSet bs){
		if(bs.getClass()==RawBitSet.class){setToMax((RawBitSet)bs);}
		else if(bs.getClass()==MultiBitSet.class){setToMax((MultiBitSet)bs);}
		else{throw new RuntimeException("Bad class: "+bs.getClass());}
	}

	/** Sets each position to the maximum of this set and a RawBitSet. For 1-bit elements max==OR==add,
	 * so this delegates to add() — which ALSO ZEROES the input bs as a side effect.
	 * Both sets must be exclusively accessible and have equal logical capacity.
	 * @param bs Non-null, distinct RawBitSet source, cleared by this merge */
	public void setToMax(RawBitSet bs){
		//TODO [map/AbstractBitSet#001]: source ownership differs from the preserving MultiBitSet
		//overload. Retained pending a contract decision; Sketch.mergeBitSets discards its source.
		add(bs);
	}

	/** Sets each position to the maximum of this set and a MultiBitSet of matching shape.
	 * Unlike setToMax(RawBitSet), this does NOT clear the input (see #001 asymmetry note above).
	 * This must be a MultiBitSet of equal logical capacity. Both sets require exclusive
	 * access; passing this set itself leaves its values unchanged.
	 * @param bs Non-null MultiBitSet to merge, preserved */
	public void setToMax(MultiBitSet bs){
		assert(this.getClass()==bs.getClass()) : this.getClass()+", "+bs.getClass();
		assert(bits()==bs.bits());
		assert(capacity()==bs.capacity()) : capacity()+", "+bs.capacity();
		final int[] rbsArray=bs.array();
		final int rbsLength=bs.length();
		for(int i=0; i<rbsLength; i++){
			final int value=rbsArray[i];
			if(value!=0){setToMax(i, value);} //per-cell max; input NOT zeroed (cf. the Raw path)
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Storage        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds packed counts to one active storage cell, saturating each element separately.
	 * The caller must leave unused elements in the last cell zero.
	 * @param cell Storage-cell index in [0, length())
	 * @param mask Packed counts in this implementation's element width
	 */
	public abstract void addToCell(final int cell, final int mask);
	/**
	 * Takes the elementwise maximum within one active storage cell, not the integer maximum.
	 * The caller must leave unused elements in the last cell zero.
	 * @param cell Storage-cell index in [0, length())
	 * @param mask Packed counts in this implementation's element width
	 */
	public abstract void setToMax(final int cell, final int mask);

	/** Zeroes active storage cells, including their padding; inactive spare cells are untouched. */
	public abstract void clear();
	/**
	 * Sets logical capacity, reusing storage when the previous reuse bound is sufficient.
	 * Reused cells are not cleared; a replacement allocation discards all prior counts.
	 * Extra cells do not extend the reuse bound and are ignored when storage is reused.
	 * Rounded cell count plus extra must fit a Java array; callers validate these bounds.
	 * @param capacity Nonnegative logical element capacity
	 * @param extra Nonnegative extra storage cells for a replacement array
	 */
	public abstract void setCapacity(long capacity, int extra);
	/** Counts nonzero elements in active cells, including any nonzero last-cell padding.
	 * @return Number of nonzero packed elements, not the sum of their count values */
	public abstract long cardinality();
	/** Returns current logical capacity, which can be smaller than reusable storage.
	 * @return Logical element capacity */
	public abstract long capacity();
	/** Returns the number of active 32-bit storage cells.
	 * @return Rounded-up cell count for the logical capacity, excluding extra cells */
	public abstract int length();
	/** Returns the number of bits used per element (1 or 2).
	 * @return Bits per element */
	public abstract int bits();

	/*--------------------------------------------------------------*/
	/*----------------          Formatting          ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns a string of non-zero positions and counts in the form {(pos,count), ...}.
	 * Scans the logical capacity, so this is diagnostic rather than a sparse traversal.
	 * Capacity must fit a nonnegative int for the element-access API used here.
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

}
