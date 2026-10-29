package tax;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import shared.Timer;

/**
 * Disk-backed GI-to-TaxID table using memory-mapped files.
 * Replaces GiToTaxid's in-memory int[][] table with a flat binary file:
 * a 16-byte header containing [long maxGi][long reserved], followed by
 * maxGi+1 big-endian int taxIDs indexed directly by GI number. Lookups use
 * segmented read-only mappings because each MappedByteBuffer is limited to an
 * int-sized byte range.
 *
 * @author Brian Bushnell, Chloe
 */
public class DiskGiTable{

	/*--------------------------------------------------------------*/
	/*----------------        Construction         ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens a memory-mapped GI table and maps every body segment read-only.
	 * The instance owns the random-access file and channel until {@link #close()}.
	 * @param path Disk table built by {@link #build(String, int[][], long)}
	 * @throws IOException If the table cannot be opened or mapped */
	public DiskGiTable(String path) throws IOException{
		final File file_=new File(path);
		if(!file_.exists()){throw new IOException("DiskGiTable file not found: "+path);}

		RandomAccessFile raf_=null;
		FileChannel channel_=null;
		try{
			raf_=new RandomAccessFile(file_, "r");
			channel_=raf_.getChannel();

			final long maxGi_=raf_.readLong();
			raf_.readLong(); //reserved

			final long dataOffset_=HEADER_SIZE;
			final long dataLength_=(maxGi_+1)*4L;
			final MappedByteBuffer[] segments_=mapSegments(channel_, dataOffset_, dataLength_);

			file=file_;
			raf=raf_;
			channel=channel_;
			maxGi=maxGi_;
			dataOffset=dataOffset_;
			dataLength=dataLength_;
			segments=segments_;
			raf_=null;
			channel_=null;
		}catch(IOException|RuntimeException|Error e){
			try{
				closeBoth(channel_, raf_);
			}catch(Throwable closeError){
				addSuppressed(e, closeError);
			}
			throw e;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Lookup              ----------------*/
	/*--------------------------------------------------------------*/

	/** Looks up one GI number by direct index into the mapped flat table.
	 * @param gi GI number to resolve
	 * @return TaxID, -1 for negative GI numbers, or -2 past the table maximum */
	public int get(long gi){
		if(gi<0 || gi>maxGi){return gi<0 ? -1 : -2;}
		long bytePos=gi*4L;
		int segIdx=(int)(bytePos/SEGMENT_SIZE);
		int segOff=(int)(bytePos%SEGMENT_SIZE);
		return segments[segIdx].getInt(segOff);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Building             ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds a disk GI table from the in-memory GiToTaxid arrays.
	 * Missing slices or cells are written as zeroes, matching GiToTaxid's dense
	 * default for unmapped GI numbers.
	 * @param outPath Output binary table path
	 * @param array Two-dimensional GI table loaded by GiToTaxid
	 * @param maxGi Largest GI number to write
	 * @throws IOException If the output table cannot be written */
	public static void build(String outPath, int[][] array, long maxGi) throws IOException{
		Timer t=new Timer();
		System.err.println("Building DiskGiTable: maxGi="+maxGi+", file="+outPath);

		try(RandomAccessFile raf=new RandomAccessFile(outPath, "rw")){
			raf.writeLong(maxGi);
			raf.writeLong(0); //reserved

			final int SHIFT=30;
			final long LOWERMASK=(1L<<SHIFT)-1;
			final int BUFSIZE=16384; //16KB write buffer
			byte[] buf=new byte[BUFSIZE];
			int pos=0;

			for(long gi=0; gi<=maxGi; gi++){
				final long upper=gi>>>SHIFT;
				final int lower=(int)(gi&LOWERMASK);
				int tid=0;
				if(upper<array.length){
					int[] slice=array[(int)upper];
					if(slice!=null && lower<slice.length){
						tid=slice[lower];
					}
				}
				buf[pos]=(byte)(tid>>>24);
				buf[pos+1]=(byte)(tid>>>16);
				buf[pos+2]=(byte)(tid>>>8);
				buf[pos+3]=(byte)(tid);
				pos+=4;
				if(pos>=BUFSIZE){
					raf.write(buf, 0, pos);
					pos=0;
				}
			}
			if(pos>0){raf.write(buf, 0, pos);}
		}

		t.stop();
		long fileSize=new File(outPath).length();
		System.err.println("DiskGiTable built: "+fileSize/(1024*1024)+"MB in "+t);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Mapping              ----------------*/
	/*--------------------------------------------------------------*/

	/** Maps a file range in chunks that stay below MappedByteBuffer's int limit.
	 * @param channel Channel to map
	 * @param offset First byte of the logical table body
	 * @param length Number of bytes to map
	 * @return Read-only segment mappings covering the requested range
	 * @throws IOException If a segment cannot be mapped */
	private static MappedByteBuffer[] mapSegments(FileChannel channel, long offset, long length)
			throws IOException{
		int numSegments=(int)((length+SEGMENT_SIZE-1)/SEGMENT_SIZE);
		MappedByteBuffer[] segs=new MappedByteBuffer[numSegments];
		for(int i=0; i<numSegments; i++){
			long segStart=offset+((long)i*SEGMENT_SIZE);
			long segLen=Math.min(SEGMENT_SIZE, offset+length-segStart);
			segs[i]=channel.map(FileChannel.MapMode.READ_ONLY, segStart, segLen);
		}
		return segs;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Closing              ----------------*/
	/*--------------------------------------------------------------*/

	/** Closes the backing channel and file handle.
	 * This releases owned backing handles, but does not explicitly unmap segment buffers.
	 * @throws IOException If either close fails with an IOException */
	public void close() throws IOException{
		closeBoth(channel, raf);
	}

	/** Closes both file resources, preserving both failures if necessary.
	 * RuntimeException and Error close failures propagate unchecked.
	 * @param channel File channel to close; may be null
	 * @param raf Random-access file to close; may be null
	 * @throws IOException If a close failure was an IOException */
	private static void closeBoth(FileChannel channel, RandomAccessFile raf) throws IOException{
		Throwable thrown=null;
		if(channel!=null){
			try{
				channel.close();
			}catch(Throwable e){
				thrown=addSuppressed(thrown, e);
			}
		}
		if(raf!=null){
			try{
				raf.close();
			}catch(Throwable e){
				thrown=addSuppressed(thrown, e);
			}
		}
		if(thrown instanceof IOException){throw (IOException)thrown;}
		if(thrown instanceof RuntimeException){throw (RuntimeException)thrown;}
		if(thrown instanceof Error){throw (Error)thrown;}
	}

	/** Adds a secondary close failure unless it is already the primary failure.
	 * @param thrown Primary failure, or null if no failure has been captured
	 * @param e Newly caught failure
	 * @return Primary failure to keep propagating */
	private static Throwable addSuppressed(Throwable thrown, Throwable e){
		if(thrown==null){return e;}
		if(thrown!=e){thrown.addSuppressed(e);}
		return thrown;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Fields              ----------------*/
	/*--------------------------------------------------------------*/

	/** Backing binary table path. */
	private final File file;
	/** Owned random-access file handle. */
	private final RandomAccessFile raf;
	/** Channel used to create read-only mappings. */
	private final FileChannel channel;
	/** Read-only mappings covering the table body. */
	private final MappedByteBuffer[] segments;
	/** Largest GI number represented by this table. */
	private final long maxGi;
	/** Byte offset of GI 0. */
	private final long dataOffset;
	/** Bytes in the flat GI-to-taxID body. */
	private final long dataLength;

	/** File header length in bytes. */
	static final int HEADER_SIZE=16;
	/** Mapping length below the two-gigabyte buffer ceiling, aligned to four bytes. */
	static final long SEGMENT_SIZE=Integer.MAX_VALUE-(Integer.MAX_VALUE%4); //~2GB, aligned to 4 bytes

}
