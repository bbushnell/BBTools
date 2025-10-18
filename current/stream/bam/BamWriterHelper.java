package stream.bam;

import java.io.IOException;
import java.io.OutputStream;

/**
 * Helper class for writing BAM binary structures with little-endian byte order.
 * All multi-byte integers in BAM format are little-endian.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BamWriterHelper {

	public BamWriterHelper(OutputStream out) {
		this.out = out;
	}

	/**
	 * Write a 32-bit signed integer in little-endian format.
	 */
	public void writeInt32(int val) throws IOException {
		out.write(val & 0xFF);
		out.write((val >> 8) & 0xFF);
		out.write((val >> 16) & 0xFF);
		out.write((val >> 24) & 0xFF);
	}

	/**
	 * Write a 32-bit unsigned integer (stored as long) in little-endian format.
	 */
	public void writeUint32(long val) throws IOException {
		writeInt32((int)val);
	}

	/**
	 * Write a 16-bit signed integer in little-endian format.
	 */
	public void writeInt16(int val) throws IOException {
		out.write(val & 0xFF);
		out.write((val >> 8) & 0xFF);
	}

	/**
	 * Write a 16-bit unsigned integer in little-endian format.
	 */
	public void writeUint16(int val) throws IOException {
		writeInt16(val);
	}

	/**
	 * Write an 8-bit unsigned integer.
	 */
	public void writeUint8(int val) throws IOException {
		out.write(val & 0xFF);
	}

	/**
	 * Write a byte array.
	 */
	public void writeBytes(byte[] data) throws IOException {
		out.write(data);
	}

	/**
	 * Write a byte array subset.
	 */
	public void writeBytes(byte[] data, int off, int len) throws IOException {
		out.write(data, off, len);
	}

	/**
	 * Write a string as US-ASCII bytes (no null terminator).
	 */
	public void writeString(String s) throws IOException {
		out.write(s.getBytes("US-ASCII"));
	}

	/**
	 * Write a 32-bit IEEE float in little-endian format.
	 */
	public void writeFloat(float val) throws IOException {
		int bits = Float.floatToIntBits(val);
		writeInt32(bits);
	}

	private final OutputStream out;
}
