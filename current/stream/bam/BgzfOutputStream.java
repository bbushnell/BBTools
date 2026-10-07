package stream.bam;

import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.CRC32;
import java.util.zip.Deflater;

/**
 * Writes BGZF (Blocked GZIP Format) compressed data.
 * BGZF is a variant of gzip with concatenated blocks, each max 64KB uncompressed.
 * Used by BAM files for indexing support. This implementation buffers at most
 * 65,280 input bytes per block and owns the supplied output stream on close.
 * Mutable compression state requires a single writer or external synchronization.
 * Flush pending data before calling {@link #writeEOF()}; {@link #close()} flushes
 * and closes the destination but deliberately does not write that marker.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BgzfOutputStream extends OutputStream{

	/** Enables diagnostic block-size printing when compiled true. */
	private static final boolean DEBUG=false;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a raw-deflate BGZF writer using compression level six.
	 * @param out Destination stream, closed by this writer's close method
	 */
	public BgzfOutputStream(OutputStream out){
		this(out, 6); // Default compression level 6
	}

	/**
	 * Creates a writer with a fresh raw deflater and checksum accumulator.
	 * @param out Destination stream, closed by this writer's close method
	 * @param compressionLevel Level forwarded unchanged to Deflater
	 */
	public BgzfOutputStream(OutputStream out, int compressionLevel){
		this.out=out;
		this.deflater=new Deflater(compressionLevel, true); // true = nowrap mode for raw deflate
		this.crc=new CRC32();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Buffers the low eight bits, emitting a block when its input buffer fills.
	 * @param b Byte value to append
	 * @throws IOException If block output fails
	 */
	@Override
	public void write(int b) throws IOException{
		buffer[bufferPos++]=(byte)b;
		if(bufferPos>=MAX_BLOCK_SIZE){flushBlock();}
	}

	/**
	 * Copies a supplied range into successive blocks, retaining a final partial block.
	 * The caller supplies a valid range; this method does not prevalidate it and
	 * performs no work when len is nonpositive.
	 * @param b Source array
	 * @param off First source index
	 * @param len Number of bytes to append
	 * @throws IOException If block output fails
	 */
	@Override
	public void write(byte[] b, int off, int len) throws IOException{
		while(len>0){
			int available=MAX_BLOCK_SIZE-bufferPos;
			int toWrite=Math.min(available, len);
			System.arraycopy(b, off, buffer, bufferPos, toWrite);
			bufferPos+=toWrite;
			off+=toWrite;
			len-=toWrite;

			if(bufferPos>=MAX_BLOCK_SIZE){flushBlock();}
		}
	}

	/**
	 * Emits any pending input as one block, then flushes the destination.
	 * Does not emit an EOF marker.
	 * @throws IOException If block output or destination flushing fails
	 */
	@Override
	public void flush() throws IOException{
		if(bufferPos>0){flushBlock();}
		out.flush();
	}

	/**
	 * Writes the 28-byte empty-block EOF marker directly to the destination.
	 * Does not flush pending input, mark this writer closed, or prevent duplicate
	 * markers. Call flush first when the marker must follow all written data.
	 * @throws IOException If marker output fails
	 */
	public void writeEOF() throws IOException{
		// Fixed caller ordering, STR208 / stream/bam/BgzfOutputStream#002:
		// BamOutputStream now flushes this buffer before calling writeEOF on close.
		// Direct callers still own flush-before-marker ordering and marker uniqueness.
		// Standard 28-byte EOF marker from SAMv1.pdf page 14
		byte[] eof=new byte[]{
			0x1f, (byte)0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
			0x00, (byte)0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
			0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
			0x00, 0x00, 0x00, 0x00
		};
		out.write(eof);
	}

	/**
	 * Flushes pending data, releases the deflater and closes the supplied destination.
	 * Deliberately omits the EOF marker: callers own its emission. For a marker last,
	 * call flush, then writeEOF, then close, without intervening writes.
	 * This method has no closed-state guard or finally-based cleanup.
	 * @throws IOException If flushing or destination close fails
	 */
	@Override
	//Historical rationale distinguishes this ST close from MT/MT2 close, which own
	//marker emission. BamOutputStream explicitly flushes, writes EOF, then closes ST.
	public void close() throws IOException{
		flush();
		deflater.end();
		out.close();
	}

	/**
	 * Compresses nonempty pending input and writes its header, payload and footer.
	 * Computes CRC32 from the input, encodes total block size minus one, and clears
	 * the input count only after all writes complete. Empty input is a no-op.
	 * @throws IOException If block output fails
	 */
	private void flushBlock() throws IOException{
		if(bufferPos==0){return;}

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
		byte[] compressed=new byte[MAX_BLOCK_SIZE+1024];
		int compressedSize=0;
		while(!deflater.finished()){
			int n=deflater.deflate(compressed, compressedSize, compressed.length-compressedSize);
			if(n==0 && deflater.needsInput()){break;}
			compressedSize+=n;
			if(compressedSize==compressed.length && !deflater.finished()){
				compressed=java.util.Arrays.copyOf(compressed, compressed.length*2); //extreme safety net
			}
		}

		// Calculate CRC32 of uncompressed data
		crc.reset();
		crc.update(buffer, 0, bufferPos);
		long crcValue=crc.getValue();

		// Calculate BSIZE: total block size minus 1
		// Total block size = header(10) + XLEN(2) + BC subfield(6) + compressed data + footer(8)
		// = 10 + 2 + 6 + compressedSize + 8 = 26 + compressedSize
		// BSIZE = 26 + compressedSize - 1 = 25 + compressedSize
		int bsize=25+compressedSize;
		assert(bsize>=27 && bsize<=65535) : "BSIZE overflow: "+bsize+" (uncompressed block must be <=0xff00; "+
			"with MAX_BLOCK_SIZE=0xff00 this can't fire on valid input). [#001 loud backstop]";

		// Write gzip header with BC subfield
		writeGzipHeader(bsize);

		// Write compressed data
		out.write(compressed, 0, compressedSize);

		// Write footer: CRC32 (4 bytes) + ISIZE (4 bytes)
		writeInt32((int)crcValue);
		int uncompressedSize=bufferPos;
		writeInt32(uncompressedSize);

		if(DEBUG){System.err.println("BGZF block: uncompressed="+uncompressedSize+", compressed="+compressedSize+", bsize="+bsize);}

		// Reset buffer
		bufferPos=0;
	}

	/**
	 * Writes the fixed 18-byte gzip header and BC subfield.
	 * @param bsize Total block size minus one, checked by the caller
	 * @throws IOException If header output fails
	 */
	private void writeGzipHeader(int bsize) throws IOException{
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
	 * Writes the low 16 bits in little-endian order.
	 * @param val Value to encode
	 * @throws IOException If either byte write fails
	 */
	private void writeInt16(int val) throws IOException{
		out.write(val&0xFF);
		out.write((val>>8)&0xFF);
	}

	/**
	 * Writes all 32 bits in little-endian order.
	 * @param val Value to encode
	 * @throws IOException If a byte write fails
	 */
	private void writeInt32(int val) throws IOException{
		out.write(val&0xFF);
		out.write((val>>8)&0xFF);
		out.write((val>>16)&0xFF);
		out.write((val>>24)&0xFF);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Destination owned by this writer. */
	private final OutputStream out;
	/** Reusable raw-deflate compressor. */
	private final Deflater deflater;
	/** Reusable checksum accumulator for the current input block. */
	private final CRC32 crc;
	/** Uncompressed input for the next block. */
	private final byte[] buffer=new byte[MAX_BLOCK_SIZE];
	/** Number of pending input bytes. */
	private int bufferPos=0;

	//[stream/bam/BgzfOutputStream#001 + block-size lead FIXED 2026-06-20 (greenlit)] 0xff00 (65280), NOT
	//65536: caps the UNCOMPRESSED block so compressed+overhead fits the 16-bit BSIZE (see BgzfOutputStreamMT).
	/** Maximum uncompressed input bytes per emitted block. */
	private static final int MAX_BLOCK_SIZE=65280; // 0xff00 - BGZF/samtools uncompressed-block cap
}
