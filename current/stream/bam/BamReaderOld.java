package stream.bam;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Reads BAM binary structures with proper little-endian handling.
 * All multi-byte integers in BAM format are little-endian.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BamReaderOld {//TODO: May need to rewrite this faster, not using ByteBuffer

	public BamReaderOld(InputStream in) {
		this.in = in;
		this.buffer = ByteBuffer.allocate(8).order(ByteOrder.LITTLE_ENDIAN);
	}

	/**
	 * Read a 32-bit signed integer (little-endian).
	 */
	public int readInt32() throws IOException {
		readIntoBuffer(4);
		buffer.flip();
		return buffer.getInt();
	}

	/**
	 * Read a 32-bit unsigned integer as long (little-endian).
	 */
	public long readUint32() throws IOException {
		readIntoBuffer(4);
		buffer.flip();
		return buffer.getInt() & 0xFFFFFFFFL;
	}

	/**
	 * Read a 16-bit signed integer (little-endian).
	 */
	public short readInt16() throws IOException {
		readIntoBuffer(2);
		buffer.flip();
		return buffer.getShort();
	}

	/**
	 * Read a 16-bit unsigned integer as int (little-endian).
	 */
	public int readUint16() throws IOException {
		readIntoBuffer(2);
		buffer.flip();
		return buffer.getShort() & 0xFFFF;
	}

	/**
	 * Read an 8-bit unsigned integer.
	 */
	public int readUint8() throws IOException {
		int b = in.read();
		if (b < 0) {
			throw new EOFException();
		}
		return b;
	}

	/**
	 * Read exactly n bytes.
	 */
	public byte[] readBytes(int n) throws IOException {
		byte[] result = new byte[n];
		int offset = 0;
		while (offset < n) {
			int bytesRead = in.read(result, offset, n - offset);
			if (bytesRead < 0) {
				throw new EOFException("Expected " + n + " bytes, got " + offset);
			}
			offset += bytesRead;
		}
		return result;
	}

	/**
	 * Read n bytes and interpret as ASCII string (no NUL terminator included).
	 */
	public String readString(int n) throws IOException {
		byte[] bytes = readBytes(n);
		return new String(bytes, 0, n, java.nio.charset.StandardCharsets.US_ASCII);
	}

	/**
	 * Read bytes into the internal buffer.
	 */
	private void readIntoBuffer(int n) throws IOException {
		buffer.clear();
		buffer.limit(n);
		while (buffer.hasRemaining()) {
			int bytesRead = in.read(buffer.array(), buffer.position(), buffer.remaining());
			if (bytesRead < 0) {
				throw new EOFException();
			}
			buffer.position(buffer.position() + bytesRead);
		}
	}

	private final InputStream in;
	private final ByteBuffer buffer;
}
