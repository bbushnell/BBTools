package stream.bam;

import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.CRC32;
import java.util.zip.Deflater;

/**
 * Writes BGZF (Blocked GZIP Format) compressed data.
 * BGZF is a variant of gzip with concatenated blocks, each max 64KB uncompressed.
 * Used by BAM files for indexing support.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BgzfOutputStream extends OutputStream {

	private static final boolean DEBUG = false;

	public BgzfOutputStream(OutputStream out) {
		this(out, 6); // Default compression level 6
	}

	public BgzfOutputStream(OutputStream out, int compressionLevel) {
		this.out = out;
		this.deflater = new Deflater(compressionLevel, true); // true = nowrap mode for raw deflate
		this.crc = new CRC32();
	}

	@Override
	public void write(int b) throws IOException {
		buffer[bufferPos++] = (byte)b;
		if (bufferPos >= MAX_BLOCK_SIZE) {
			flushBlock();
		}
	}

	@Override
	public void write(byte[] b, int off, int len) throws IOException {
		while (len > 0) {
			int available = MAX_BLOCK_SIZE - bufferPos;
			int toWrite = Math.min(available, len);
			System.arraycopy(b, off, buffer, bufferPos, toWrite);
			bufferPos += toWrite;
			off += toWrite;
			len -= toWrite;

			if (bufferPos >= MAX_BLOCK_SIZE) {
				flushBlock();
			}
		}
	}

	@Override
	public void flush() throws IOException {
		if (bufferPos > 0) {
			flushBlock();
		}
		out.flush();
	}

	/**
	 * Write the BGZF EOF marker (28-byte empty block).
	 * This is a standard BGZF block with zero uncompressed data.
	 */
	public void writeEOF() throws IOException {
		// Standard 28-byte EOF marker from SAMv1.pdf page 14
		byte[] eof = new byte[]{
			0x1f, (byte)0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
			0x00, (byte)0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
			0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
			0x00, 0x00, 0x00, 0x00
		};
		out.write(eof);
	}

	@Override
	//Contract: close() deliberately does NOT write the BGZF EOF marker - the caller owns it. For BAM,
	//BamOutputStream.close() calls stOut.writeEOF() THEN stOut.close() (so EOF appears exactly once);
	//MT/MT2 close() DO write EOF internally (asymmetric contract, correctly handled per-engine by
	//BamOutputStream). Gap: a direct generic-BGZF .gz via this ST path that calls only close() ends up
	//WITHOUT the EOF sentinel (still a valid gzip-block concatenation, just no truncation-detection
	//marker) - LOW, interop-only; data intact.
	public void close() throws IOException {
		flush();
		deflater.end();
		out.close();
	}

	/**
	 * Compress and write the current buffer as a BGZF block.
	 */
	private void flushBlock() throws IOException {
		if (bufferPos == 0) {
			return;
		}

		// Compress the buffer
		deflater.reset();
		deflater.setInput(buffer, 0, bufferPos);
		deflater.finish();

		//[stream/bam/BgzfOutputStream#001] FIXED 2026-06-20 (greenlit). Two coupled mechanisms (both
		//empirically reproduced on random data): (a) the OLD single deflate() into a no-headroom buffer left
		//finished()==false on incompressible data -> dropped the deflate tail -> truncated block; (b) bsize
		//overflowed the 16-bit BSIZE -> wrapped -> corrupt framing. Fix: cap uncompressed at MAX_BLOCK_SIZE
		//(0xff00, above) so BSIZE always fits, AND give the compressed buffer +1024 headroom + a finish-loop
		//(mirrors the MT engine) so the tail is never dropped even on an incompressible block. assert(bsize)
		//is the loud backstop.
		byte[] compressed = new byte[MAX_BLOCK_SIZE + 1024];
		int compressedSize = 0;
		while(!deflater.finished()){
			int n = deflater.deflate(compressed, compressedSize, compressed.length - compressedSize);
			if(n == 0 && deflater.needsInput()){break;}
			compressedSize += n;
			if(compressedSize == compressed.length && !deflater.finished()){
				compressed = java.util.Arrays.copyOf(compressed, compressed.length * 2); //extreme safety net
			}
		}

		// Calculate CRC32 of uncompressed data
		crc.reset();
		crc.update(buffer, 0, bufferPos);
		long crcValue = crc.getValue();

		// Calculate BSIZE: total block size minus 1
		// Total block size = header(10) + XLEN(2) + BC subfield(6) + compressed data + footer(8)
		// = 10 + 2 + 6 + compressedSize + 8 = 26 + compressedSize
		// BSIZE = 26 + compressedSize - 1 = 25 + compressedSize
		int bsize = 25 + compressedSize;
		assert(bsize>=27 && bsize<=65535) : "BSIZE overflow: "+bsize+" (uncompressed block must be <=0xff00; "+
			"with MAX_BLOCK_SIZE=0xff00 this can't fire on valid input). [#001 loud backstop]";

		// Write gzip header with BC subfield
		writeGzipHeader(bsize);

		// Write compressed data
		out.write(compressed, 0, compressedSize);

		// Write footer: CRC32 (4 bytes) + ISIZE (4 bytes)
		writeInt32((int)crcValue);
		int uncompressedSize = bufferPos;
		writeInt32(uncompressedSize);

		if (DEBUG) {
			System.err.println("BGZF block: uncompressed=" + uncompressedSize + ", compressed=" + compressedSize + ", bsize=" + bsize);
		}

		// Reset buffer
		bufferPos = 0;
	}

	/**
	 * Write BGZF-compliant gzip header with BC subfield.
	 */
	private void writeGzipHeader(int bsize) throws IOException {
		out.write(31);          // ID1
		out.write(139);         // ID2
		out.write(8);           // CM (compression method = DEFLATE)
		out.write(4);           // FLG (FEXTRA flag set)
		writeInt32(0);          // MTIME (4 bytes)
		out.write(0);           // XFL
		out.write(255);         // OS (unknown)
		writeInt16(6);          // XLEN (length of extra field)
		out.write(66);          // SI1 'B'
		out.write(67);          // SI2 'C'
		writeInt16(2);          // SLEN (BC data length)
		writeInt16(bsize);      // BSIZE (already calculated as block size minus 1)
	}

	/**
	 * Write a 16-bit integer in little-endian format.
	 */
	private void writeInt16(int val) throws IOException {
		out.write(val & 0xFF);
		out.write((val >> 8) & 0xFF);
	}

	/**
	 * Write a 32-bit integer in little-endian format.
	 */
	private void writeInt32(int val) throws IOException {
		out.write(val & 0xFF);
		out.write((val >> 8) & 0xFF);
		out.write((val >> 16) & 0xFF);
		out.write((val >> 24) & 0xFF);
	}

	private final OutputStream out;
	private final Deflater deflater;
	private final CRC32 crc;
	private final byte[] buffer = new byte[MAX_BLOCK_SIZE];
	private int bufferPos = 0;

	//[stream/bam/BgzfOutputStream#001 + block-size lead FIXED 2026-06-20 (greenlit)] 0xff00 (65280), NOT
	//65536: caps the UNCOMPRESSED block so compressed+overhead fits the 16-bit BSIZE (see BgzfOutputStreamMT).
	private static final int MAX_BLOCK_SIZE = 65280; // 0xff00 - BGZF/samtools uncompressed-block cap
}
