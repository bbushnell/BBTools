package structures;

/**
 * Mutable little-endian cursor over a borrowed byte array.
 * Positions are absolute array indices, including when wrapping a subrange.
 * The logical limit controls remaining-data queries; reads and cursor changes do
 * not enforce it. Callers must provide valid ranges and enough bytes for each read.
 * Java array bounds checks still apply, and a failed scalar read may advance the
 * cursor before throwing. Instances require externally serialized access.
 * 
 * @author Isla
 * @date October 30, 2025
 */
public class BinaryByteWrapperLE{
	
	/** Borrows the whole array and starts at index zero.
	 * @param array_ Non-null backing storage; no copy is made */
	public BinaryByteWrapperLE(byte[] array_){
		wrap(array_);
	}
	
	/** Borrows an array range and starts at its absolute offset.
	 * @param array_ Non-null backing storage; no copy is made
	 * @param offset Nonnegative first readable array index
	 * @param length Nonnegative byte count fitting within the backing array */
	public BinaryByteWrapperLE(byte[] array_, int offset, int length){
		wrap(array_, offset, length);
	}
	
	/** Replaces the borrowed storage and resets the cursor to zero.
	 * @param array_ Non-null array whose length becomes the logical limit */
	public void wrap(byte[] array_){
		array=array_;
		position=0;
		limit=array.length;
	}
	
	/** Replaces the borrowed storage and records a range without validating it.
	 * @param array_ Non-null backing storage
	 * @param offset Nonnegative absolute initial cursor position
	 * @param length Nonnegative byte count; {@code offset+length} must fit the array */
	public void wrap(byte[] array_, int offset, int length){
		array=array_;
		position=offset;
		limit=offset+length;
	}
	
	/** Returns the absolute array index of the next read. */
	public int position(){return position;}
	
	/** Sets the absolute cursor position without changing or checking the limit.
	 * @param newPosition Array index for the next read; caller ensures valid bounds
	 * @return This cursor */
	public BinaryByteWrapperLE position(int newPosition){
		position=newPosition;
		return this;
	}
	
	/** Reports whether the cursor is below the logical limit; does not check read width. */
	public boolean hasRemaining(){return position<limit;}
	
	/** Returns logical limit minus cursor position; may be negative after unchecked movement. */
	public int remaining(){return limit-position;}
	
	/** Moves the cursor by a signed byte count without validating the result.
	 * @param skip Byte displacement; a negative value moves backward
	 * @return true if the resulting cursor is below the logical limit
	 */
	public boolean skip(int skip){
		position+=skip;
		return position<limit;
	}
	
	/** Reads one signed byte and advances the cursor by one. */
	public byte get(){
		return array[position++];
	}
	
	/** Copies enough bytes to fill the destination and advances by its length.
	 * @param dst Destination array
	 * @return This cursor */
	public BinaryByteWrapperLE get(byte[] dst){
		return get(dst, 0, dst.length);
	}
	
	/** Copies a byte range and advances the cursor after a successful copy.
	 * The backing-array bounds apply, but the logical limit is not checked.
	 * @param dst Destination array
	 * @param offset First destination index
	 * @param length Nonnegative byte count fitting both array ranges
	 * @return This cursor */
	public BinaryByteWrapperLE get(byte[] dst, int offset, int length){
		System.arraycopy(array, position, dst, offset, length);
		position+=length;
		return this;
	}
	
	/** Reads a signed little-endian short and advances by two bytes. */
	public short getShort(){
		int b0=array[position++]&0xFF;
		int b1=array[position++]&0xFF;
		return (short)((b1<<8)|b0);
	}
	
	/** Reads a signed little-endian int and advances by four bytes. */
	public int getInt(){
		int b0=array[position++]&0xFF;
		int b1=array[position++]&0xFF;
		int b2=array[position++]&0xFF;
		int b3=array[position++]&0xFF;
		return (b3<<24)|(b2<<16)|(b1<<8)|b0;
	}
	
	/** Reads a signed little-endian long and advances by eight bytes. */
	public long getLong(){
		long b0=array[position++]&0xFFL;
		long b1=array[position++]&0xFFL;
		long b2=array[position++]&0xFFL;
		long b3=array[position++]&0xFFL;
		long b4=array[position++]&0xFFL;
		long b5=array[position++]&0xFFL;
		long b6=array[position++]&0xFFL;
		long b7=array[position++]&0xFFL;
		return (b7<<56)|(b6<<48)|(b5<<40)|(b4<<32)|(b3<<24)|(b2<<16)|(b1<<8)|b0;
	}
	
	/** Reads four little-endian bytes as IEEE float bits and advances by four. */
	public float getFloat(){
		return Float.intBitsToFloat(getInt());
	}
	
	/** Reads eight little-endian bytes as IEEE double bits and advances by eight. */
	public double getDouble(){
		return Double.longBitsToDouble(getLong());
	}
	
	/** Borrowed backing storage; direct replacement does not reset the cursor or limit. */
	public byte[] array;
	/** Absolute index of the next byte to read. */
	private int position;
	/** Exclusive logical end used by queries; reads do not enforce this boundary. */
	private int limit;
}
