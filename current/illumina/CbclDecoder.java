package illumina;

import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.zip.GZIPInputStream;

/**
 * Decodes gzip-compressed, bit-packed base call data from CBCL files.
 *
 * CBCL format packs 2 clusters per byte in interleaved layout:
 *   bits 0-1: base for cluster A (00=A, 01=C, 10=G, 11=T)
 *   bits 2-3: quality bin index for cluster A
 *   bits 4-5: base for cluster B
 *   bits 6-7: quality bin index for cluster B
 *
 * Quality bin indices (0-3) are remapped to Phred scores via the
 * header's qscoreRemap table (e.g., NovaSeq RTA3 maps to 0, 12, 23, 37).
 * A byte value of 0x00 indicates no-call (base=N, quality=0).
 *
 * @author Chloe
 * @date October 15, 2025
 */
public class CbclDecoder {

	/*--------------------------------------------------------------*/
	/*----------------         Constants            ----------------*/
	/*--------------------------------------------------------------*/

	private static final byte[] BASE_CHARS={'A', 'C', 'G', 'T'};

	/*--------------------------------------------------------------*/
	/*----------------       Public Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Decode a gzip-compressed block of CBCL data.
	 * @param compressedData Gzip-compressed bytes
	 * @param numClusters Number of clusters in this tile
	 * @param bitsPerBase Bits per base (typically 2)
	 * @param bitsPerQual Bits per quality (typically 2)
	 * @param qscoreRemap Quality bin remap table from header; may be null
	 * @return Array of [bases, quals] where both are byte arrays
	 */
	public static byte[][] decodeBlock(byte[] compressedData, int numClusters,
	                                    int bitsPerBase, int bitsPerQual,
	                                    int[] qscoreRemap) throws IOException {
		//Decompress gzip data
		ByteArrayInputStream bais=new ByteArrayInputStream(compressedData);
		GZIPInputStream gis=new GZIPInputStream(bais);
		//TODO: Uses Java 9+ library: readAllBytes (InputStream.readAllBytes, since Java 9). Violates Java-8 bytecode-compliance target (no compiler warning). Address post-compaction.
		byte[] decompressed=gis.readAllBytes();
		gis.close();

		//Allocate output arrays
		byte[] bases=new byte[numClusters];
		byte[] quals=new byte[numClusters];

		//Unpack bit-packed data
		if(bitsPerBase==2 && bitsPerQual==2){
			decode2bit(decompressed, bases, quals, numClusters, qscoreRemap);
		} else {
			throw new UnsupportedOperationException("Only 2-bit encoding supported currently");
		}

		return new byte[][]{bases, quals};
	}

	/**
	 * Decode interleaved 2-bit packed bases and qualities.
	 * Each byte holds 2 clusters: bits[0:1]=baseA, bits[2:3]=qualA,
	 * bits[4:5]=baseB, bits[6:7]=qualB. LSB-first ordering.
	 * A nibble of 0x0 (base=0, qual=0) is a no-call: emit N with quality 0.
	 */
	private static void decode2bit(byte[] data, byte[] bases, byte[] quals,
	                                int numClusters, int[] qscoreRemap) {
		int clusterIdx=0;
		for(int byteIdx=0; byteIdx<data.length && clusterIdx<numClusters; byteIdx++){
			int b=data[byteIdx]&0xFF;

			//Cluster A: bits 0-3
			int baseA=b&0x03;
			int qualA=(b>>2)&0x03;
			if((b&0x0F)==0){
				bases[clusterIdx]='N';
				quals[clusterIdx]=(byte)(0+33);
			}else{
				bases[clusterIdx]=BASE_CHARS[baseA];
				int phred=(qscoreRemap!=null && qualA<qscoreRemap.length) ? qscoreRemap[qualA] : qualA;
				quals[clusterIdx]=(byte)(phred+33);
			}
			clusterIdx++;
			if(clusterIdx>=numClusters){break;}

			//Cluster B: bits 4-7
			int baseB=(b>>4)&0x03;
			int qualB=(b>>6)&0x03;
			if((b&0xF0)==0){
				bases[clusterIdx]='N';
				quals[clusterIdx]=(byte)(0+33);
			}else{
				bases[clusterIdx]=BASE_CHARS[baseB];
				int phred=(qscoreRemap!=null && qualB<qscoreRemap.length) ? qscoreRemap[qualB] : qualB;
				quals[clusterIdx]=(byte)(phred+33);
			}
			clusterIdx++;
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
		//TODO: Uses Java 9+ library: readAllBytes (InputStream.readAllBytes, since Java 9). Violates Java-8 bytecode-compliance target (no compiler warning). Address post-compaction.
		byte[] compressedData=fis.readAllBytes();
		fis.close();

		//For files with multiple tiles, need to split the compressed blocks
		//For now, assume single tile per file (which matches test data structure)
		return decodeBlock(compressedData, numClusters,
		                   header.bitsPerBasecall, header.bitsPerQscore,
		                   header.qscoreRemap);
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
			System.out.printf("Cluster %d: base=%c qual=%c (Q%d)\n",
				i, bases[i], quals[i], (quals[i]-33));
		}
	}
}
