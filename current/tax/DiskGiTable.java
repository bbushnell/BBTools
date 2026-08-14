package tax;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;

import shared.Timer;

/**
 * Disk-backed GI-to-TaxID table using memory-mapped files.
 * Replaces the in-memory int[][] array in GiToTaxid with a flat
 * binary file where index=GI and value=taxID (4 bytes per entry).
 *
 * File format:
 *   Header (16 bytes): [long maxGi] [long reserved]
 *   Body: (maxGi+1) x 4-byte ints, indexed by GI number
 *
 * MappedByteBuffer is limited to 2GB per mapping, so the file
 * is mapped in segments.
 *
 * @author Brian Bushnell, Chloe
 */
public class DiskGiTable {

	/*--------------------------------------------------------------*/
	/*----------------        Construction         ----------------*/
	/*--------------------------------------------------------------*/

	public DiskGiTable(String path) throws IOException {
		file=new File(path);
		if(!file.exists()){throw new IOException("DiskGiTable file not found: "+path);}

		raf=new RandomAccessFile(file, "r");
		channel=raf.getChannel();

		maxGi=raf.readLong();
		raf.readLong(); //reserved

		dataOffset=HEADER_SIZE;
		dataLength=(maxGi+1)*4L;

		segments=mapSegments(channel, dataOffset, dataLength);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Lookup              ----------------*/
	/*--------------------------------------------------------------*/

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

	/**
	 * Builds a disk GI table from the in-memory GiToTaxid arrays.
	 * Call after GiToTaxid.initialize() has loaded the data.
	 */
	public static void build(String outPath, int[][] array, long maxGi) throws IOException {
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

	private static MappedByteBuffer[] mapSegments(FileChannel channel, long offset, long length) throws IOException {
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

	public void close() throws IOException {
		channel.close();
		raf.close();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Fields              ----------------*/
	/*--------------------------------------------------------------*/

	private final File file;
	private final RandomAccessFile raf;
	private final FileChannel channel;
	private final MappedByteBuffer[] segments;
	private final long maxGi;
	private final long dataOffset;
	private final long dataLength;

	static final int HEADER_SIZE=16;
	static final long SEGMENT_SIZE=Integer.MAX_VALUE-(Integer.MAX_VALUE%4); //~2GB, aligned to 4 bytes

}
