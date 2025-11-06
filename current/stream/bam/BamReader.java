package stream.bam;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;

import structures.BinaryByteWrapperLE;

/**
 * Reads BAM binary structures with proper little-endian handling.
 * All multi-byte integers in BAM format are little-endian.
 * Optimized version using BinaryByteWrapperLE instead of ByteBuffer.
 *
 * @author Brian Bushnell, Chloe, Isla
 * @date November 5, 2025
 */
public class BamReader {

	public BamReader(InputStream in){
		this.in=in;
		this.temp=new byte[8];
		this.wrapper=new BinaryByteWrapperLE(temp);
	}

	/**
	 * Read a 32-bit signed integer (little-endian).
	 */
	public int readInt32() throws IOException {
		readFully(temp, 0, 4);
		wrapper.position(0);
		return wrapper.getInt();
	}

	/**
	 * Read a 32-bit unsigned integer as long (little-endian).
	 */
	public long readUint32() throws IOException {
		readFully(temp, 0, 4);
		wrapper.position(0);
		return wrapper.getInt()&0xFFFFFFFFL;
	}

	/**
	 * Read a 16-bit signed integer (little-endian).
	 */
	public short readInt16() throws IOException {
		readFully(temp, 0, 2);
		wrapper.position(0);
		return wrapper.getShort();
	}

	/**
	 * Read a 16-bit unsigned integer as int (little-endian).
	 */
	public int readUint16() throws IOException {
		readFully(temp, 0, 2);
		wrapper.position(0);
		return wrapper.getShort()&0xFFFF;
	}

	/**
	 * Read an 8-bit unsigned integer.
	 */
	public int readUint8() throws IOException {
		int b=in.read();
		if(b<0){
			throw new EOFException();
		}
		return b;
	}

	/**
	 * Read exactly n bytes.
	 */
	public byte[] readBytes(int n) throws IOException {
		byte[] result=new byte[n];
		readFully(result, 0, n);
		return result;
	}

	/**
	 * Read n bytes and interpret as ASCII string (no NUL terminator included).
	 */
	public String readString(int n) throws IOException {
		byte[] bytes=readBytes(n);
		return new String(bytes, 0, n, java.nio.charset.StandardCharsets.US_ASCII);
	}

	/**
	 * Read exactly n bytes into array at offset.
	 */
	private void readFully(byte[] array, int offset, int n) throws IOException {
		int total=0;
		while(total<n){
			int bytesRead=in.read(array, offset+total, n-total);
			if(bytesRead<0){
				throw new EOFException("Expected "+n+" bytes, got "+total);
			}
			total+=bytesRead;
		}
	}

	private final InputStream in;
	private final byte[] temp;
	private final BinaryByteWrapperLE wrapper;
}