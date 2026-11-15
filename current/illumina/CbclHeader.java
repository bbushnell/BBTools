package illumina;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Parses a CBCL header containing quality remaps, ordered tile cluster counts,
 * per-tile compressed and uncompressed block sizes, and the non-PF-excluded
 * flag. Only version 1 headers with 2-bit bases and 2-bit quality bins are
 * supported.
 *
 * @author Chloe
 * @date October 15, 2025
 */
public class CbclHeader{

	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parse a CBCL header from a local file.
	 * @param filename path to a CBCL file
	 * @throws IOException if the header cannot be read or is malformed
	 */
	public CbclHeader(String filename) throws IOException{
		readHeader(filename);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Populate this object from the header at the start of {@code filename}.
	 * This method closes its own stream before returning or throwing.
	 */
	private void readHeader(String filename) throws IOException{
		FileInputStream fis=new FileInputStream(filename);
		try{

			//Read first 12 bytes to get version, header size, and basic params
			byte[] initialBytes=new byte[12];
			readFully(fis, initialBytes);

			ByteBuffer bb=ByteBuffer.wrap(initialBytes);
			bb.order(ByteOrder.LITTLE_ENDIAN);

			version=bb.getShort();
			headerSize=bb.getInt();
			bitsPerBasecall=bb.get();
			bitsPerQscore=bb.get();
			numQscoreBins=bb.getInt();
			if(version!=1){
				throw new IOException("Unsupported CBCL version: "+version);
			}
			if(bitsPerBasecall!=2 || bitsPerQscore!=2){
				throw new IOException("Unsupported CBCL encoding: "+bitsPerBasecall+
						" base bits and "+bitsPerQscore+" quality bits");
			}
			if(numQscoreBins<0){
				throw new IOException("Negative CBCL qscore-bin count: "+numQscoreBins);
			}

			//Read rest of header based on headerSize
			if(headerSize<12){
				throw new IOException("Invalid header size: "+headerSize);
			}
			int remainingHeaderBytes=headerSize-12;

			byte[] restOfHeader=new byte[remainingHeaderBytes];
			readFully(fis, restOfHeader);

			//Parse quality bin information
			bb=ByteBuffer.wrap(restOfHeader);
			bb.order(ByteOrder.LITTLE_ENDIAN);
			long qscoreBytes=8L*numQscoreBins;
			if(qscoreBytes+4>bb.remaining()){
				throw new IOException("CBCL header too short for "+numQscoreBins+
						" qscore-bin records");
			}

			qscoreBins=new int[numQscoreBins];
			qscoreRemap=new int[numQscoreBins];

			for(int i=0; i<numQscoreBins; i++){
				qscoreBins[i]=bb.getInt();
				qscoreRemap[i]=bb.getInt();
			}

			//Read number of tiles
			int numTiles=bb.getInt();
			if(numTiles<0){
				throw new IOException("Negative CBCL tile count: "+numTiles);
			}

			final long expectedRemaining=numTiles*16L+1;
			if(bb.remaining()<expectedRemaining){
				throw new IOException("CBCL header too short for "+numTiles+" tile records");
			}
			tileMetadata=new LinkedHashMap<>();
			tileUncompressedSizes=new LinkedHashMap<>();
			tileCompressedSizes=new LinkedHashMap<>();
			for(int i=0; i<numTiles; i++){
				int tileNum=bb.getInt();
				int clusterCount=bb.getInt();
				int uncompressedSize=bb.getInt();
				int compressedSize=bb.getInt();
				if(clusterCount<0){
					throw new IOException("Negative CBCL cluster count for tile "+tileNum);
				}
				if(uncompressedSize<0 || compressedSize<0){
					throw new IOException("Negative CBCL block size for tile "+tileNum);
				}
				tileMetadata.put(tileNum, clusterCount);
				tileUncompressedSizes.put(tileNum, uncompressedSize);
				tileCompressedSizes.put(tileNum, compressedSize);
			}
			nonPfExcluded=bb.get()!=0;

			//Store position where compressed data starts
			compressedDataOffset=headerSize;

		}finally{
			fis.close();
		}
	}

	/**
	 * Fill {@code buffer} exactly or throw if the stream ends first.
	 */
	private static void readFully(FileInputStream fis, byte[] buffer) throws IOException{
		int offset=0;
		while(offset<buffer.length){
			int read=fis.read(buffer, offset, buffer.length-offset);
			if(read<0){
				throw new IOException("Could not read full CBCL header");
			}
			offset+=read;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Return the tile-to-cluster-count map in header order.
	 * Keys are Illumina tile numbers and values are total cluster counts.
	 * The returned map is this header's live metadata map; callers should
	 * treat it as read-only.
	 */
	public Map<Integer, Integer> getTileMetadata(){
		return tileMetadata;
	}

	/**
	 * Return the tile-to-compressed-block-size map in header order.
	 * Values are compressed gzip block byte counts, not cumulative offsets.
	 * The returned map is this header's live metadata map; callers should
	 * treat it as read-only.
	 */
	public Map<Integer, Integer> getTileCompressedSizes(){
		return tileCompressedSizes;
	}

	/**
	 * Return the tile-to-uncompressed-block-size map in header order.
	 * Values are decompressed packed-block byte counts and should equal
	 * {@code ceil(clusterCount/2)} for each tile.
	 * The returned map is this header's live metadata map; callers should
	 * treat it as read-only.
	 */
	public Map<Integer, Integer> getTileUncompressedSizes(){
		return tileUncompressedSizes;
	}

	public int[] getQscoreBins(){
		return qscoreBins;
	}

	public int[] getQscoreRemap(){
		return qscoreRemap;
	}

	public int getCompressedDataOffset(){
		return compressedDataOffset;
	}

	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append("CBCL Header:\n");
		sb.append("  Version: ").append(version).append("\n");
		sb.append("  Header size: ").append(headerSize).append("\n");
		sb.append("  Bits per basecall: ").append(bitsPerBasecall).append("\n");
		sb.append("  Bits per qscore: ").append(bitsPerQscore).append("\n");
		sb.append("  Num qscore bins: ").append(numQscoreBins).append("\n");
		sb.append("  Tiles: ").append(tileMetadata.size()).append("\n");
		for(Map.Entry<Integer, Integer> e : tileMetadata.entrySet()){
			sb.append("    Tile ").append(e.getKey())
			  .append(": ").append(e.getValue()).append(" clusters\n");
		}
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** CBCL header version; only version 1 is supported. */
	public short version;
	/** Total header byte count from the file. */
	public int headerSize;
	/** Bits per basecall; only 2 is supported. */
	public byte bitsPerBasecall;
	/** Bits per quality bin; only 2 is supported. */
	public byte bitsPerQscore;
	/** Number of quality-bin boundary/remap pairs. */
	public int numQscoreBins;
	/** Quality-bin boundaries in file order. */
	public int[] qscoreBins;
	/** ASCII-offset-free Phred remaps in file order. */
	public int[] qscoreRemap;
	/** Tile number to cluster count map, preserving file order. */
	public Map<Integer, Integer> tileMetadata;
	/** Tile number to uncompressed packed-block byte count, preserving file order. */
	public Map<Integer, Integer> tileUncompressedSizes;
	/** Tile number to compressed gzip-block byte count, preserving file order. */
	public Map<Integer, Integer> tileCompressedSizes;
	/** Absolute byte offset where the first compressed tile block starts. */
	public int compressedDataOffset;
	/** True when the CBCL blocks omit non-PF clusters. */
	public boolean nonPfExcluded;

	/*--------------------------------------------------------------*/
	/*----------------        Test Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws Exception{
		if(args.length<1){
			System.err.println("Usage: CbclHeader <cbcl_file>");
			System.exit(1);
		}

		CbclHeader header=new CbclHeader(args[0]);
		System.out.println(header);
	}
}
