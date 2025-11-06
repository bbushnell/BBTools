package structures;

/**
 * Fast little-endian binary reader wrapping a byte array.
 * Drop-in replacement for ByteBuffer with fixed little-endian byte order.
 * Optimized for sequential reading without bounds checking overhead.
 * 
 * @author Isla
 * @date October 30, 2025
 */
public class BinaryByteWrapperLE {
	
	public BinaryByteWrapperLE(byte[] array_){
		wrap(array_);
	}
	
	public BinaryByteWrapperLE(byte[] array_, int offset, int length){
		wrap(array_, offset, length);
	}
	
	public void wrap(byte[] array_) {
		array=array_;
		position=0;
		limit=array.length;
	}
	
	public void wrap(byte[] array_, int offset, int length) {
		array=array_;
		position=offset;
		limit=offset+length;
	}
	
	/** Returns current read position */
	public int position(){return position;}
	
	/** Sets read position */
	public BinaryByteWrapperLE position(int newPosition){
		position=newPosition;
		return this;
	}
	
	/** Returns true if data remains between position and limit */
	public boolean hasRemaining(){return position<limit;}
	
	/** Returns number of bytes remaining */
	public int remaining(){return limit-position;}
	
	/** 
	 * Skip this many bytes 
	 * @return true if there is remaining data 
	 */
	public boolean skip(int skip){
		position+=skip;
		return position<limit;
	}
	
	/** Reads single byte as signed value */
	public byte get(){
		return array[position++];
	}
	
	/** Reads multiple bytes into destination array */
	public BinaryByteWrapperLE get(byte[] dst){
		return get(dst, 0, dst.length);
	}
	
	/** Reads multiple bytes into destination array at offset */
	public BinaryByteWrapperLE get(byte[] dst, int offset, int length){
		System.arraycopy(array, position, dst, offset, length);
		position+=length;
		return this;
	}
	
	/** Reads 2-byte little-endian signed short */
	public short getShort(){
		int b0=array[position++]&0xFF;
		int b1=array[position++]&0xFF;
		return (short)((b1<<8)|b0);
	}
	
	/** Reads 4-byte little-endian signed int */
	public int getInt(){
		int b0=array[position++]&0xFF;
		int b1=array[position++]&0xFF;
		int b2=array[position++]&0xFF;
		int b3=array[position++]&0xFF;
		return (b3<<24)|(b2<<16)|(b1<<8)|b0;
	}
	
	/** Reads 8-byte little-endian signed long */
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
	
	/** Reads 4-byte little-endian float */
	public float getFloat(){
		return Float.intBitsToFloat(getInt());
	}
	
	/** Reads 8-byte little-endian double */
	public double getDouble(){
		return Double.longBitsToDouble(getLong());
	}
	
	/** Underlying byte array */
	public byte[] array;
	/** Current read position in array */
	private int position;
	/** Maximum readable position (exclusive) */
	private int limit;
}
