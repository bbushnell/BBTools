package stream.bam;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.GZIPInputStream;
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
	//TODO: Possible bug [stream/bam/BgzfInputStream#001] - the 16-bit in-block uoffset `bufferPos & 0xFFFF`
	//WRAPS to 0 when bufferPos==65536 (a full 65536-byte uncompressed block, exhausted). The voffset is NOT
	//normalized to the next block when a block is fully consumed, so for a block <65536 it returns
	//(thisBlock<<16)|isize - which is seek-EQUIVALENT to (nextBlock<<16)|0 (a reader at uoffset==isize
	//advances to the next block on read), so that case is fine. But at EXACTLY 65536 the &0xFFFF collapses
	//isize→0, yielding (thisBlock<<16)|0 = the START of this block. If a record boundary falls exactly at the
	//end of a 65536-byte block, BamIndexWriter captures this wrapped voffset as the record's .bai chunk offset
	//-> a query seeks to the block START and returns the wrong (earlier) record. NARROW + entirely contingent
	//on a 65536-byte block, which only BBTools' own over-sized writer produces (spec/samtools cap at 0xff00) -
	//see the BGZF block-size lead in leads.md. LOW/latent. Fix: cap uncompressed blocks <=0xff00 (root cause),
	//and/or normalize the voffset to (nextBlock<<16)|0 when bufferPos==bufferLimit.
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
		while (true) {
			if (plainGzipMode) {
				if (fillPlainGzipBuffer()) {
					bufferPos = 0;
					return true;
				}
				exitPlainGzipMode();
				continue;
			}

			// Read gzip header (minimum 10 bytes)
			byte[] header = new byte[10];
			long blockStart = filePointer;
			int bytesRead = readFully(header, 0, header.length);
			if (bytesRead == 0) {
				return false; // EOF
			}
			if (bytesRead < header.length) {
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
				enterPlainGzip(header, null, null);
				continue;
			}

			// Read XLEN (2 bytes, little-endian)
			byte[] xlenBytes = new byte[2];
			if (readFully(xlenBytes, 0, 2) < 2) {
				throw new EOFException("Truncated XLEN");
			}
			int xlen = ((xlenBytes[1] & 0xFF) << 8) | (xlenBytes[0] & 0xFF);

			// Read extra field and find BC subfield
			byte[] extra = new byte[xlen];
			if (readFully(extra, 0, xlen) < xlen) {
				throw new EOFException("Truncated extra field");
			}

			int bsize = findBsizeInExtra(extra, xlen);
			if (bsize < 0) {
				enterPlainGzip(header, xlenBytes, extra);
				continue;
			}

			// Calculate compressed data length
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
				//Single inflate into the 65536 buffer suffices for a valid BGZF block (isize<=64KB fits). A
				//malformed block claiming isize>65536 fills the buffer (returns 65536, finished()==false), but
				//the isize check below catches the mismatch and crashes loud - no silent truncation.
				bufferLimit = inflater.inflate(uncompressedBuffer);
			} catch (DataFormatException e) {
				throw new IOException("Decompression failed", e);
			}

			//Two crash-loud integrity guards (BGZF correctness contract): the decompressed size MUST equal the
			//block's stored ISIZE, and the CRC32 of the decompressed bytes MUST match the stored CRC. Either
			//mismatch throws - a corrupt/truncated BGZF block fails LOUD here, never returns wrong bytes.
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
		try {
			if (plainGzipStream != null) {
				plainGzipStream.close();
			}
		} finally {
			plainGzipStream = null;
			in.close();
		}
	}

	private final InputStream in;
	private final Inflater inflater;
	private final byte[] uncompressedBuffer = new byte[65536];
	private int bufferPos = 0;
	private int bufferLimit = 0;
	private long filePointer = 0L;
	private long blockCompressedStart = 0L;
	private boolean plainGzipMode = false;
	private GZIPInputStream plainGzipStream = null;

	//Clever (verified): a non-BGZF plain gzip (no FEXTRA, or FEXTRA without a BC subfield) is handled by
	//replaying the ALREADY-CONSUMED header/xlen/extra bytes through a PrefixedInputStream in front of the live
	//`in`, then handing the reconstructed stream to a stock GZIPInputStream. So the reader transparently
	//accepts both BGZF and ordinary .gz without losing the bytes it already read to detect the format. (Note:
	//getVirtualOffset is not meaningful in plain-gzip mode - but plain .gz isn't block-indexable anyway, and
	//BAM is always BGZF, so BamIndexWriter never hits this path.)
	private void enterPlainGzip(byte[] header, byte[] xlenBytes, byte[] extra) throws IOException {
		if (plainGzipMode) {
			return;
		}
		byte[] prefix = buildPrefix(header, xlenBytes, extra);
		plainGzipStream = new GZIPInputStream(new PrefixedInputStream(prefix, in), uncompressedBuffer.length);
		plainGzipMode = true;
		bufferPos = 0;
		bufferLimit = 0;
	}

	private boolean fillPlainGzipBuffer() throws IOException {
		if (plainGzipStream == null) {
			return false;
		}
		int n = 0;
		while (n == 0) {
			n = plainGzipStream.read(uncompressedBuffer, 0, uncompressedBuffer.length);
			if (n < 0) {
				return false;
			}
		}
		bufferLimit = n;
		return true;
	}

	private void exitPlainGzipMode() throws IOException {
		if (plainGzipStream != null) {
			plainGzipStream.close();
		}
		plainGzipStream = null;
		plainGzipMode = false;
	}

	private byte[] buildPrefix(byte[] header, byte[] xlenBytes, byte[] extra) {
		int prefixLen = header.length;
		if (xlenBytes != null) {
			prefixLen += xlenBytes.length;
		}
		if (extra != null) {
			prefixLen += extra.length;
		}

		byte[] prefix = new byte[prefixLen];
		int pos = 0;
		System.arraycopy(header, 0, prefix, pos, header.length);
		pos += header.length;
		if (xlenBytes != null) {
			System.arraycopy(xlenBytes, 0, prefix, pos, xlenBytes.length);
			pos += xlenBytes.length;
		}
		if (extra != null && extra.length > 0) {
			System.arraycopy(extra, 0, prefix, pos, extra.length);
		}
		return prefix;
	}

	private static final class PrefixedInputStream extends InputStream {
		private final byte[] prefix;
		private int position = 0;
		private final InputStream tail;

		PrefixedInputStream(byte[] prefix, InputStream tail) {
			this.prefix = prefix;
			this.tail = tail;
		}

		@Override
		public int read() throws IOException {
			if (position < prefix.length) {
				return prefix[position++] & 0xFF;
			}
			return tail.read();
		}

		@Override
		public int read(byte[] b, int off, int len) throws IOException {
			if (position < prefix.length) {
				int toCopy = Math.min(len, prefix.length - position);
				System.arraycopy(prefix, position, b, off, toCopy);
				position += toCopy;
				return toCopy;
			}
			return tail.read(b, off, len);
		}

		@Override
		public void close() {
			// Do not close the tail stream; caller manages lifecycle.
		}
	}
}
