package stream.bam;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * Reads BGZF (Blocked GZIP Format) compressed data.
 * BGZF is a variant of gzip with concatenated blocks, each max 64KB uncompressed.
 * Used by BAM files for indexing support.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BgzfInputStream extends InputStream {

	public BgzfInputStream(InputStream in) {
		this.in = in;
		this.inflater = new Inflater(true); // true = nowrap mode for raw deflate
	}

	/**
	 * Return the current BGZF virtual offset.
	 * Upper 48 bits represent the compressed block offset; lower 16 bits represent the in-block position.
	 */
	public long getVirtualOffset() {
		return (blockCompressedStart << 16) | (bufferPos & 0xFFFFL);
	}

	@Override
	public int read() throws IOException {
		if (bufferPos >= bufferLimit) {
			if (!readBlock()) {
				return -1;
			}
		}
		return uncompressedBuffer[bufferPos++] & 0xFF;
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException {
		if (b == null) {
			throw new NullPointerException();
		} else if (off < 0 || len < 0 || len > b.length - off) {
			throw new IndexOutOfBoundsException();
		} else if (len == 0) {
			return 0;
		}

		int totalRead = 0;
		while (totalRead < len) {
			if (bufferPos >= bufferLimit) {
				if (!readBlock()) {
					return totalRead == 0 ? -1 : totalRead;
				}
			}

			int available = bufferLimit - bufferPos;
			int toRead = Math.min(available, len - totalRead);
			System.arraycopy(uncompressedBuffer, bufferPos, b, off + totalRead, toRead);
			bufferPos += toRead;
			totalRead += toRead;
		}

		return totalRead;
	}

	/**
	 * Read and decompress the next BGZF block.
	 * @return true if block was read, false on EOF
	 * @throws IOException on read error or malformed BGZF
	 */
	private boolean readBlock() throws IOException {
		// Read gzip header (minimum 10 bytes)
		byte[] header = new byte[12];
		long blockStart = filePointer;
		int bytesRead = readFully(header, 0, 10);
		if (bytesRead == 0) {
			return false; // EOF
		}
		if (bytesRead < 10) {
			throw new EOFException("Truncated BGZF block header");
		}
		blockCompressedStart = blockStart;

		// Verify gzip signature
		if ((header[0] & 0xFF) != 31 || (header[1] & 0xFF) != 139) {
			throw new IOException("Not a gzip file");
		}

		// Check compression method (should be 8 = DEFLATE)
		if (header[2] != 8) {
			throw new IOException("Unsupported compression method: " + header[2]);
		}

		// Check flags
		int flags = header[3] & 0xFF;
		boolean fextra = (flags & 0x04) != 0;

		if (!fextra) {
			throw new IOException("BGZF block missing FEXTRA flag");
		}

		// Skip MTIME (4 bytes), XFL, OS
		// Already read, just validate positions

		// Read XLEN (2 bytes, little-endian)
		if (readFully(header, 0, 2) < 2) {
			throw new EOFException("Truncated XLEN");
		}
		int xlen = ((header[1] & 0xFF) << 8) | (header[0] & 0xFF);

		// Read extra field and find BC subfield
		byte[] extra = new byte[xlen];
		if (readFully(extra, 0, xlen) < xlen) {
			throw new EOFException("Truncated extra field");
		}

		int bsize = findBsizeInExtra(extra, xlen);
		if (bsize < 0) {
			throw new IOException("BGZF block missing BC subfield");
		}

		// Calculate compressed data length
		// bsize is total block size - 1
		// Already read: 10 (header) + 2 (xlen) + xlen (extra)
		// Still need: compressed_data + 8 (CRC32 + ISIZE)
		int alreadyRead = 10 + 2 + xlen;
		int remaining = (bsize + 1) - alreadyRead;

		if (remaining < 8) {
			throw new IOException("Invalid BSIZE: " + bsize);
		}

		int compressedSize = remaining - 8; // Subtract CRC32 and ISIZE
		byte[] compressed = new byte[compressedSize];
		if (readFully(compressed, 0, compressedSize) < compressedSize) {
			throw new EOFException("Truncated compressed data");
		}

		// Read CRC32 and ISIZE (8 bytes total)
		byte[] trailer = new byte[8];
		if (readFully(trailer, 0, 8) < 8) {
			throw new EOFException("Truncated block trailer");
		}

		ByteBuffer bb = ByteBuffer.wrap(trailer).order(ByteOrder.LITTLE_ENDIAN);
		long crc32 = bb.getInt() & 0xFFFFFFFFL;
		int isize = bb.getInt();

		// Decompress
		inflater.reset();
		inflater.setInput(compressed);

		try {
			bufferLimit = inflater.inflate(uncompressedBuffer);
		} catch (DataFormatException e) {
			throw new IOException("Decompression failed", e);
		}

		if (bufferLimit != isize) {
			throw new IOException("Uncompressed size mismatch: expected " + isize + ", got " + bufferLimit);
		}

		// Verify CRC32
		CRC32 crc = new CRC32();
		crc.update(uncompressedBuffer, 0, bufferLimit);
		if (crc.getValue() != crc32) {
			throw new IOException("CRC32 mismatch");
		}

		bufferPos = 0;
		return true;
	}

	/**
	 * Find BC subfield in gzip extra field and extract BSIZE.
	 * @param extra Extra field bytes
	 * @param xlen Length of extra field
	 * @return BSIZE value, or -1 if not found
	 */
	private int findBsizeInExtra(byte[] extra, int xlen) {
		int pos = 0;
		while (pos + 4 <= xlen) {
			int si1 = extra[pos] & 0xFF;
			int si2 = extra[pos + 1] & 0xFF;
			int slen = ((extra[pos + 3] & 0xFF) << 8) | (extra[pos + 2] & 0xFF);

			if (si1 == 66 && si2 == 67) { // 'B' 'C'
				if (slen == 2 && pos + 6 <= xlen) {
					return ((extra[pos + 5] & 0xFF) << 8) | (extra[pos + 4] & 0xFF);
				}
			}

			pos += 4 + slen;
		}
		return -1;
	}

	/**
	 * Read exactly n bytes from input stream.
	 * @return number of bytes read (0 on immediate EOF, n on success)
	 * @throws EOFException if EOF reached before n bytes read (but after some bytes read)
	 */
	private int readFully(byte[] b, int off, int len) throws IOException {
		int total = 0;
		while (total < len) {
			int n = in.read(b, off + total, len - total);
			if (n < 0) {
				return total;
			}
			filePointer += n;
			total += n;
		}
		return total;
	}

	@Override
	public void close() throws IOException {
		inflater.end();
		in.close();
	}

	private final InputStream in;
	private final Inflater inflater;
	private final byte[] uncompressedBuffer = new byte[65536];
	private int bufferPos = 0;
	private int bufferLimit = 0;
	private long filePointer = 0L;
	private long blockCompressedStart = 0L;
}
