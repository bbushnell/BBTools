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
 * Reads BGZF blocks through a reusable 65536-byte decompression buffer.
 * Falls back to GZIPInputStream when the parsed header has no BGZF size subfield.
 * Closing this reader closes its underlying input. Instances contain mutable
 * decompression state and are intended for use by one caller at a time.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BgzfInputStream extends InputStream{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a reader without consuming input bytes.
	 * @param in Compressed input owned and closed by this reader
	 */
	public BgzfInputStream(InputStream in){
		this.in=in;
		this.inflater=new Inflater(true); // true = nowrap mode for raw deflate
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Combines the current compressed block start and low 16 bits of the buffer position.
	 * This value is not meaningful in plain-gzip mode and does not normalize an
	 * exhausted block to the next block's start.
	 * @return Encoded block offset and in-block position; initially zero
	 */
	//TODO: Possible bug [stream/bam/BgzfInputStream#001] - bufferPos==65536 wraps the
	//low 16 bits to zero while retaining the same block start. A caller saving this
	//position could mistake the exhausted block for its beginning. Earlier review
	//proposed smaller writer blocks or normalization to the next block; the reader
	//still has this arithmetic boundary. No index-query outcome is established here.
	//The earlier claim that only BBTools writers can produce this size is unverified.
	public long getVirtualOffset(){
		return (blockCompressedStart<<16)|(bufferPos&0xFFFFL);
	}

	/**
	 * Returns one buffered byte, loading a block when the buffer is exhausted.
	 * An empty BGZF block currently exposes a buffer byte; see TODO BgzfInputStream#002.
	 * @return Unsigned byte value, or -1 when no further block is available
	 * @throws IOException If reading or decompression fails
	 */
	//TODO: Probable bug [stream/bam/BgzfInputStream#002] - readBlock can return true
	//with zero decompressed bytes. This one-shot refill then reads a buffer byte
	//despite bufferLimit==0; the bulk overload instead loops past empty blocks.
	@Override
	public int read() throws IOException{
		if(bufferPos>=bufferLimit){
			if(!readBlock()){return -1;}
		}
		return uncompressedBuffer[bufferPos++]&0xFF;
	}

	/**
	 * Copies up to len bytes, crossing blocks until the request is filled or input ends.
	 * @param b Destination array
	 * @param off First destination index
	 * @param len Maximum bytes to copy; zero returns zero without reading input
	 * @return Copied byte count, or -1 if input ends before any byte is copied
	 * @throws IOException If reading or decompression fails
	 * @throws NullPointerException If b is null
	 * @throws IndexOutOfBoundsException If the destination range is invalid
	 */
	@Override
	public int read(byte[] b, int off, int len) throws IOException{
		if(b==null){throw new NullPointerException();}else if(off<0 || len<0 || len>b.length-off){throw new IndexOutOfBoundsException();}else if(len==0){return 0;}

		int totalRead=0;
		while(totalRead<len){
			if(bufferPos>=bufferLimit){
				if(!readBlock()){return totalRead==0 ? -1 : totalRead;}
			}

			int available=bufferLimit-bufferPos;
			int toRead=Math.min(available, len-totalRead);
			System.arraycopy(uncompressedBuffer, bufferPos, b, off+totalRead, toRead);
			bufferPos+=toRead;
			totalRead+=toRead;
		}

		return totalRead;
	}

	/**
	 * Loads a BGZF block or a positive chunk from the plain-gzip fallback.
	 * A BGZF block with zero decompressed bytes still returns true. Stored size
	 * and CRC are checked against the decompressed buffer before it is exposed.
	 * @return true after loading a block or chunk, false at the end of input
	 * @throws IOException If input cannot be read, inflated, or fails a size or CRC check
	 */
	private boolean readBlock() throws IOException{
		while(true){
			if(plainGzipMode){
				if(fillPlainGzipBuffer()){
					bufferPos=0;
					return true;
				}
				exitPlainGzipMode();
				continue;
			}

			// Read gzip header (minimum 10 bytes)
			byte[] header=new byte[10];
			long blockStart=filePointer;
			int bytesRead=readFully(header, 0, header.length);
			if(bytesRead==0){
				return false; // EOF
			}
			if(bytesRead<header.length){throw new EOFException("Truncated BGZF block header");}
			blockCompressedStart=blockStart;

			// Verify gzip signature
			if((header[0]&0xFF)!=31 || (header[1]&0xFF)!=139){throw new IOException("Not a gzip file");}

			// Check compression method (should be 8 = DEFLATE)
			if(header[2]!=8){throw new IOException("Unsupported compression method: "+header[2]);}

			// Check flags
			int flags=header[3]&0xFF;
			boolean fextra=(flags&0x04)!=0;

			if(!fextra){
				enterPlainGzip(header, null, null);
				continue;
			}

			// Read XLEN (2 bytes, little-endian)
			byte[] xlenBytes=new byte[2];
			if(readFully(xlenBytes, 0, 2)<2){throw new EOFException("Truncated XLEN");}
			int xlen=((xlenBytes[1]&0xFF)<<8)|(xlenBytes[0]&0xFF);

			// Read extra field and find BC subfield
			byte[] extra=new byte[xlen];
			if(readFully(extra, 0, xlen)<xlen){throw new EOFException("Truncated extra field");}

			int bsize=findBsizeInExtra(extra, xlen);
			if(bsize<0){
				enterPlainGzip(header, xlenBytes, extra);
				continue;
			}

			// Calculate compressed data length
			int alreadyRead=10+2+xlen;
			int remaining=(bsize+1)-alreadyRead;

			if(remaining<8){throw new IOException("Invalid BSIZE: "+bsize);}

			int compressedSize=remaining-8; // Subtract CRC32 and ISIZE
			byte[] compressed=new byte[compressedSize];
			if(readFully(compressed, 0, compressedSize)<compressedSize){throw new EOFException("Truncated compressed data");}

			// Read CRC32 and ISIZE (8 bytes total)
			byte[] trailer=new byte[8];
			if(readFully(trailer, 0, 8)<8){throw new EOFException("Truncated block trailer");}

			ByteBuffer bb=ByteBuffer.wrap(trailer).order(ByteOrder.LITTLE_ENDIAN);
			long crc32=bb.getInt()&0xFFFFFFFFL;
			int isize=bb.getInt();

			// Decompress
			inflater.reset();
			inflater.setInput(compressed);

			try{
				//Inflate once into the fixed-size buffer. The stored ISIZE is compared with
				//the returned byte count below; this is not an exhaustive format validator.
				bufferLimit=inflater.inflate(uncompressedBuffer);
			}catch(DataFormatException e){
				throw new IOException("Decompression failed", e);
			}

			//Reject a mismatch between the returned byte count and stored ISIZE.
			//The following CRC check compares the bytes actually produced.
			if(bufferLimit!=isize){throw new IOException("Uncompressed size mismatch: expected "+isize+", got "+bufferLimit);}

			// Verify CRC32
			CRC32 crc=new CRC32();
			crc.update(uncompressedBuffer, 0, bufferLimit);
			if(crc.getValue()!=crc32){throw new IOException("CRC32 mismatch");}

			bufferPos=0;
			return true;
		}
	}

	/**
	 * Scans gzip extra subfields for a two-byte BC payload and decodes BSIZE.
	 * @param extra Extra field bytes
	 * @param xlen Number of extra-field bytes to examine
	 * @return Unsigned little-endian BSIZE, or -1 if no matching complete payload is found
	 */
	private int findBsizeInExtra(byte[] extra, int xlen){
		int pos=0;
		while(pos+4<=xlen){
			int si1=extra[pos]&0xFF;
			int si2=extra[pos+1]&0xFF;
			int slen=((extra[pos+3]&0xFF)<<8)|(extra[pos+2]&0xFF);

			if(si1==66 && si2==67){ // 'B' 'C'
				if(slen==2 && pos+6<=xlen){return ((extra[pos+5]&0xFF)<<8)|(extra[pos+4]&0xFF);}
			}

			pos+=4+slen;
		}
		return -1;
	}

	/**
	 * Reads until len bytes have arrived or the underlying input reaches EOF.
	 * Advances filePointer for each byte obtained through this helper.
	 * @param b Destination array
	 * @param off First destination index
	 * @param len Requested byte count
	 * @return Number of bytes read, including a short count on partial EOF
	 * @throws IOException If the underlying input throws while reading
	 */
	private int readFully(byte[] b, int off, int len) throws IOException{
		int total=0;
		while(total<len){
			int n=in.read(b, off+total, len-total);
			if(n<0){return total;}
			filePointer+=n;
			total+=n;
		}
		return total;
	}

	/**
	 * Releases the inflater and closes the optional gzip wrapper and underlying input.
	 * The underlying close is attempted even if closing the wrapper throws.
	 * @throws IOException If an input close fails
	 */
	@Override
	public void close() throws IOException{
		inflater.end();
		try{
			if(plainGzipStream!=null){plainGzipStream.close();}
		}finally{
			plainGzipStream=null;
			in.close();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Gzip Header Replay     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Replays an already-consumed header before delegating to GZIPInputStream.
	 * Headers without FEXTRA, or with no usable BC subfield, take this path.
	 * Virtual offsets are not maintained for bytes consumed by the gzip wrapper.
	 * Does nothing if plain-gzip mode is already active.
	 * @param header Consumed fixed header bytes
	 * @param xlenBytes Consumed XLEN bytes, or null when absent
	 * @param extra Consumed extra-field bytes, or null when absent
	 * @throws IOException If the gzip wrapper cannot read its header
	 */
	private void enterPlainGzip(byte[] header, byte[] xlenBytes, byte[] extra) throws IOException{
		if(plainGzipMode){return;}
		byte[] prefix=buildPrefix(header, xlenBytes, extra);
		plainGzipStream=new GZIPInputStream(new PrefixedInputStream(prefix, in), uncompressedBuffer.length);
		plainGzipMode=true;
		bufferPos=0;
		bufferLimit=0;
	}

	/**
	 * Fills the shared buffer from the gzip wrapper, retrying zero-byte reads.
	 * @return true for a positive chunk, false if no wrapper exists or input ends
	 * @throws IOException If the gzip wrapper cannot supply decoded bytes
	 */
	private boolean fillPlainGzipBuffer() throws IOException{
		if(plainGzipStream==null){return false;}
		int n=0;
		while(n==0){
			n=plainGzipStream.read(uncompressedBuffer, 0, uncompressedBuffer.length);
			if(n<0){return false;}
		}
		bufferLimit=n;
		return true;
	}

	/**
	 * Closes the gzip wrapper and clears its reference and mode after a successful close.
	 * The prefix adapter leaves the underlying input open.
	 * @throws IOException If closing the gzip wrapper fails
	 */
	private void exitPlainGzipMode() throws IOException{
		if(plainGzipStream!=null){plainGzipStream.close();}
		plainGzipStream=null;
		plainGzipMode=false;
	}

	/**
	 * Copies consumed header components into one array in their original order.
	 * @param header Fixed header bytes
	 * @param xlenBytes Optional XLEN bytes, or null
	 * @param extra Optional extra-field bytes, or null
	 * @return Newly allocated concatenation of the supplied components
	 */
	private byte[] buildPrefix(byte[] header, byte[] xlenBytes, byte[] extra){
		int prefixLen=header.length;
		if(xlenBytes!=null){prefixLen+=xlenBytes.length;}
		if(extra!=null){prefixLen+=extra.length;}

		byte[] prefix=new byte[prefixLen];
		int pos=0;
		System.arraycopy(header, 0, prefix, pos, header.length);
		pos+=header.length;
		if(xlenBytes!=null){
			System.arraycopy(xlenBytes, 0, prefix, pos, xlenBytes.length);
			pos+=xlenBytes.length;
		}
		if(extra!=null && extra.length>0){System.arraycopy(extra, 0, prefix, pos, extra.length);}
		return prefix;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Replays a borrowed byte-array prefix, then reads from a borrowed input stream. */
	private static final class PrefixedInputStream extends InputStream{
		/**
		 * Retains the prefix and tail without copying or consuming either.
		 * @param prefix Header bytes to replay
		 * @param tail Remaining compressed input
		 */
		PrefixedInputStream(byte[] prefix, InputStream tail){
			this.prefix=prefix;
			this.tail=tail;
		}

		/** Reports unread prefix and tail bytes, saturating the sum to avoid overflow.
		 * Older GZIPInputStream implementations use available() to find concatenated members;
		 * inherited zero can discard a buffered next-member header and cause false corruption.
		 * @return Unread prefix length plus tail availability, capped at Integer.MAX_VALUE
		 * @throws IOException If querying the tail fails
		 */
		@Override
		public int available() throws IOException{
			assert(position>=0 && position<=prefix.length) : "Gzip header replay position outside prefix: "+position+" of "+prefix.length;
			return (int)Math.min(Integer.MAX_VALUE, (long)(prefix.length-position)+tail.available());
		}

		/**
		 * Returns the next unsigned prefix byte, or delegates when the prefix is exhausted.
		 * @return Unsigned byte value or the tail's EOF result
		 * @throws IOException If reading the tail fails
		 */
		@Override
		public int read() throws IOException{
			if(position<prefix.length){return prefix[position++]&0xFF;}
			return tail.read();
		}

		/**
		 * Copies a prefix chunk or delegates the whole request to the tail.
		 * A single call does not cross from a remaining prefix into the tail.
		 * @param b Destination array
		 * @param off First destination index
		 * @param len Maximum number of bytes to copy
		 * @return Prefix bytes copied, or the tail's read result
		 * @throws IOException If reading the tail fails
		 */
		@Override
		public int read(byte[] b, int off, int len) throws IOException{
			if(position<prefix.length){
				int toCopy=Math.min(len, prefix.length-position);
				System.arraycopy(prefix, position, b, off, toCopy);
				position+=toCopy;
				return toCopy;
			}
			return tail.read(b, off, len);
		}

		/** Does nothing; ownership of the tail stays with the outer reader. */
		@Override
		public void close(){
			// Do not close the tail stream; caller manages lifecycle.
		}

		/** Header bytes replayed before reading the tail. */
		private final byte[] prefix;
		/** Number of prefix bytes already returned. */
		private int position=0;
		/** Borrowed source following the prefix; this adapter does not close it. */
		private final InputStream tail;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Owned source of compressed bytes. */
	private final InputStream in;
	/** Reusable raw-DEFLATE inflater for BGZF blocks. */
	private final Inflater inflater;
	/** Reusable storage for decoded BGZF blocks or plain-gzip chunks. */
	private final byte[] uncompressedBuffer=new byte[65536];
	/** Index of the next byte to return from the decoded buffer. */
	private int bufferPos=0;
	/** Number of decoded bytes in the buffer. */
	private int bufferLimit=0;
	/** Compressed bytes consumed through readFully; excludes plain-gzip wrapper reads. */
	private long filePointer=0L;
	/** Compressed position recorded for the most recently parsed header. */
	private long blockCompressedStart=0L;
	/** Whether reads are delegated to the ordinary gzip wrapper. */
	private boolean plainGzipMode=false;
	/** Optional gzip decoder over a replayed header and the remaining input. */
	private GZIPInputStream plainGzipStream=null;

}
