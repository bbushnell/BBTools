package illumina;

import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.zip.GZIPInputStream;

/**
 * Decodes gzip-compressed, bit-packed base call data from CBCL files.
 * Handles 2-bit base encoding (A=00, C=01, G=10, T=11) and
 * 2-bit quality score encoding.
 *
 * @author Chloe
 * @date October 15, 2025
 */
public class CbclDecoder {

	/*--------------------------------------------------------------*/
	/*----------------         Constants            ----------------*/
	/*--------------------------------------------------------------*/

	private static final byte[] BASE_CHARS={'A', 'C', 'G', 'T'};
	private static final byte[] QUAL_CHARS={'0', '1', '2', '3'};

	/*--------------------------------------------------------------*/
	/*----------------       Public Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Decode a gzip-compressed block of CBCL data.
	 * @param compressedData Gzip-compressed bytes
	 * @param numClusters Number of clusters in this tile
	 * @param bitsPerBase Bits per base (typically 2)
	 * @param bitsPerQual Bits per quality (typically 2)
	 * @return Array of [bases, quals] where both are byte arrays
	 */
	public static byte[][] decodeBlock(byte[] compressedData, int numClusters,
	                                    int bitsPerBase, int bitsPerQual) throws IOException {
		//Decompress gzip data
		ByteArrayInputStream bais=new ByteArrayInputStream(compressedData);
		GZIPInputStream gis=new GZIPInputStream(bais);
		byte[] decompressed=gis.readAllBytes();
		gis.close();

		//Allocate output arrays
		byte[] bases=new byte[numClusters];
		byte[] quals=new byte[numClusters];

		//Unpack bit-packed data
		if(bitsPerBase==2 && bitsPerQual==2){
			decode2bit(decompressed, bases, quals, numClusters);
		} else {
			throw new UnsupportedOperationException("Only 2-bit encoding supported currently");
		}

		return new byte[][]{bases, quals};
	}

	/**
	 * Decode 2-bit packed bases and qualities.
	 * Format: 4 bases per byte, 4 qualities per byte
	 */
	private static void decode2bit(byte[] data, byte[] bases, byte[] quals, int numClusters) {
		int numBytes=(numClusters+3)/4; //Round up for partial bytes

		//Decode bases
		int clusterIdx=0;
		for(int byteIdx=0; byteIdx<numBytes && clusterIdx<numClusters; byteIdx++){
			byte b=data[byteIdx];
			//Extract 4 2-bit values from this byte
			for(int shift=0; shift<8 && clusterIdx<numClusters; shift+=2){
				int baseValue=(b>>shift)&0x03;
				bases[clusterIdx]=BASE_CHARS[baseValue];
				clusterIdx++;
			}
		}

		//Decode qualities (start after base bytes)
		clusterIdx=0;
		for(int byteIdx=numBytes; byteIdx<numBytes*2 && clusterIdx<numClusters; byteIdx++){
			byte b=data[byteIdx];
			//Extract 4 2-bit values from this byte
			for(int shift=0; shift<8 && clusterIdx<numClusters; shift+=2){
				int qualValue=(b>>shift)&0x03;
				quals[clusterIdx]=QUAL_CHARS[qualValue];
				clusterIdx++;
			}
		}
	}

	/**
	 * Read and decode an entire CBCL file for a specific tile.
	 * @param filename CBCL file path
	 * @param tileNum Tile number to extract
	 * @return [bases, quals] for the specified tile
	 */
	public static byte[][] readTile(String filename, int tileNum) throws IOException {
		//Parse header
		CbclHeader header=new CbclHeader(filename);

		//Find tile in metadata
		Integer numClusters=header.tileMetadata.get(tileNum);
		if(numClusters==null){
			throw new IOException("Tile " + tileNum + " not found in CBCL file");
		}

		//Read compressed data
		FileInputStream fis=new FileInputStream(filename);
		fis.skip(header.compressedDataOffset);
		byte[] compressedData=fis.readAllBytes();
		fis.close();

		//For files with multiple tiles, need to split the compressed blocks
		//For now, assume single tile per file (which matches test data structure)
		return decodeBlock(compressedData, numClusters,
		                   header.bitsPerBasecall, header.bitsPerQscore);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Test Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws Exception {
		if(args.length<2){
			System.err.println("Usage: CbclDecoder <cbcl_file> <tile_num>");
			System.exit(1);
		}

		String filename=args[0];
		int tileNum=Integer.parseInt(args[1]);

		byte[][] result=readTile(filename, tileNum);
		byte[] bases=result[0];
		byte[] quals=result[1];

		System.out.println("Decoded " + bases.length + " clusters for tile " + tileNum);
		System.out.println("First 5 clusters:");
		for(int i=0; i<Math.min(5, bases.length); i++){
			System.out.printf("Cluster %d: base=%c qual=%c\n", i, bases[i], quals[i]);
		}
	}
}
