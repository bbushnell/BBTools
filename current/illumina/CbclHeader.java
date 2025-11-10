package illumina;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Parses CBCL file headers containing metadata and tile information.
 * Header format:
 * - Version (2 bytes)
 * - Header size (4 bytes)
 * - Bits per basecall (1 byte)
 * - Bits per qscore (1 byte)
 * - Number of quality bins (4 bytes)
 * - Quality bin information
 * - Tile metadata (tile number + cluster count per tile)
 *
 * @author Chloe
 * @date October 15, 2025
 */
public class CbclHeader {

	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

	public CbclHeader(String filename) throws IOException {
		readHeader(filename);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private void readHeader(String filename) throws IOException {
		FileInputStream fis=new FileInputStream(filename);

		//Read first 12 bytes to get version, header size, and basic params
		byte[] initialBytes=new byte[12];
		int read=fis.read(initialBytes);
		if(read<12){
			throw new IOException("File too short to contain valid header");
		}

		ByteBuffer bb=ByteBuffer.wrap(initialBytes);
		bb.order(ByteOrder.LITTLE_ENDIAN);

		version=bb.getShort();
		headerSize=bb.getInt();
		bitsPerBasecall=bb.get();
		bitsPerQscore=bb.get();
		numQscoreBins=bb.getInt();

		//Read rest of header based on headerSize
		int remainingHeaderBytes=headerSize-12;
		if(remainingHeaderBytes<0){
			throw new IOException("Invalid header size: " + headerSize);
		}

		byte[] restOfHeader=new byte[remainingHeaderBytes];
		read=fis.read(restOfHeader);
		if(read<remainingHeaderBytes){
			throw new IOException("Could not read full header");
		}

		//Parse quality bin information
		bb=ByteBuffer.wrap(restOfHeader);
		bb.order(ByteOrder.LITTLE_ENDIAN);

		//Quality bins: numQscoreBins values for bin boundaries
		//Then numQscoreBins values for remapping
		qscoreBins=new int[numQscoreBins];
		qscoreRemap=new int[numQscoreBins];

		for(int i=0; i<numQscoreBins; i++){
			qscoreBins[i]=bb.getInt();
		}
		for(int i=0; i<numQscoreBins; i++){
			qscoreRemap[i]=bb.getInt();
		}

		//Read number of tiles
		int numTiles=bb.getInt();

		//Read tile metadata (tile number + cluster count)
		tileMetadata=new LinkedHashMap<>();
		for(int i=0; i<numTiles; i++){
			int tileNum=bb.getInt();
			int clusterCount=bb.getInt();
			tileMetadata.put(tileNum, clusterCount);
		}

		//Store position where compressed data starts
		compressedDataOffset=headerSize;

		fis.close();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public Map<Integer, Integer> getTileMetadata() {
		return tileMetadata;
	}

	public int[] getQscoreBins() {
		return qscoreBins;
	}

	public int[] getQscoreRemap() {
		return qscoreRemap;
	}

	public int getCompressedDataOffset() {
		return compressedDataOffset;
	}

	@Override
	public String toString() {
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

	public short version;
	public int headerSize;
	public byte bitsPerBasecall;
	public byte bitsPerQscore;
	public int numQscoreBins;
	public int[] qscoreBins;
	public int[] qscoreRemap;
	public Map<Integer, Integer> tileMetadata;
	public int compressedDataOffset;

	/*--------------------------------------------------------------*/
	/*----------------        Test Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws Exception {
		if(args.length<1){
			System.err.println("Usage: CbclHeader <cbcl_file>");
			System.exit(1);
		}

		CbclHeader header=new CbclHeader(args[0]);
		System.out.println(header);
	}
}
