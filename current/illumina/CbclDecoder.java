package illumina;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Map.Entry;
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
 * Quality bin indices are remapped to Phred scores via the header's
 * qscoreRemap table. A 4-bit cluster nibble of zero indicates a no-call for
 * that cluster and emits base N with Phred 0.
 *
 * @author Chloe
 * @date October 15, 2025
 */
public class CbclDecoder{

	/*--------------------------------------------------------------*/
	/*----------------         Constants            ----------------*/
	/*--------------------------------------------------------------*/

	private static final byte[] BASE_CHARS={'A', 'C', 'G', 'T'};

	/*--------------------------------------------------------------*/
	/*----------------       Public Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Decode a gzip-compressed block of CBCL data.
	 * The inflated packed block must contain exactly {@code ceil(numClusters/2)}
	 * bytes. Returned bases are ASCII bases; returned qualities are ASCII
	 * Phred+33 bytes.
	 * @param compressedData Gzip-compressed bytes for one tile
	 * @param numClusters number of clusters in this tile
	 * @param bitsPerBase bits per base; only 2 is supported
	 * @param bitsPerQual bits per quality bin; only 2 is supported
	 * @param qscoreRemap quality-bin remap table from the header; may be null
	 * @return array of {@code [bases, quals]} where both entries are byte arrays
	 */
	public static byte[][] decodeBlock(byte[] compressedData, int numClusters,
	                                    int bitsPerBase, int bitsPerQual,
	                                    int[] qscoreRemap) throws IOException{
		//Decompress gzip data
		ByteArrayInputStream bais=new ByteArrayInputStream(compressedData);
		GZIPInputStream gis=new GZIPInputStream(bais);
		byte[] decompressed=readAll(gis);
		gis.close();
		int expectedSize=(int)((numClusters+1L)/2);
		if(decompressed.length!=expectedSize){
			throw new IOException("CBCL block expands to "+decompressed.length+
					" bytes but "+numClusters+" clusters require "+expectedSize);
		}

		//Allocate output arrays
		byte[] bases=new byte[numClusters];
		byte[] quals=new byte[numClusters];

		//Unpack bit-packed data
		if(bitsPerBase==2 && bitsPerQual==2){
			decode2bit(decompressed, bases, quals, numClusters, qscoreRemap);
		}else{
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
	                                int numClusters, int[] qscoreRemap){
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
	public static byte[][] readTile(String filename, int tileNum) throws IOException{
		//Parse header
		CbclHeader header=new CbclHeader(filename);

		//Find tile in metadata
		Integer numClusters=header.tileMetadata.get(tileNum);
		if(numClusters==null){
			throw new IOException("Tile "+tileNum+" not found in CBCL file");
		}

		//Read compressed data
		FileInputStream fis=new FileInputStream(filename);
		byte[] compressedData;
		try{
			long offset=header.compressedDataOffset;
			Integer compressedSize=header.tileCompressedSizes.get(tileNum);
			for(Entry<Integer, Integer> entry : header.tileCompressedSizes.entrySet()){
				if(entry.getKey()==tileNum){break;}
				offset+=entry.getValue();
			}
			skipFully(fis, offset);
			compressedData=readFully(fis, compressedSize);
		}finally{
			fis.close();
		}
		int uncompressedSize=header.tileUncompressedSizes.get(tileNum);
		int expectedSize=(int)((numClusters+1L)/2);
		if(uncompressedSize!=expectedSize){
			throw new IOException("Tile "+tileNum+" declares uncompressed size "+uncompressedSize+
					" but "+numClusters+" clusters require "+expectedSize+" bytes");
		}
		return decodeBlock(compressedData, numClusters,
				header.bitsPerBasecall, header.bitsPerQscore, header.qscoreRemap);
	}

	private static byte[] readFully(InputStream is, int length) throws IOException{
		byte[] data=new byte[length];
		int offset=0;
		while(offset<length){
			int read=is.read(data, offset, length-offset);
			if(read<0){throw new IOException("Unexpected end of compressed CBCL block");}
			offset+=read;
		}
		return data;
	}

	private static byte[] readAll(InputStream is) throws IOException{
		ByteArrayOutputStream baos=new ByteArrayOutputStream();
		byte[] buffer=new byte[32768];
		for(int read=is.read(buffer); read>=0; read=is.read(buffer)){
			if(read>0){baos.write(buffer, 0, read);}
		}
		return baos.toByteArray();
	}

	private static void skipFully(InputStream is, long length) throws IOException{
		long skipped=0;
		while(skipped<length){
			long current=is.skip(length-skipped);
			if(current>0){
				skipped+=current;
			}else if(is.read()<0){
				throw new IOException("Unexpected end of CBCL file");
			}else{
				skipped++;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Test Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws Exception{
		if(args.length<2){
			System.err.println("Usage: CbclDecoder <cbcl_file> <tile_num>");
			System.exit(1);
		}

		String filename=args[0];
		int tileNum=Integer.parseInt(args[1]);

		byte[][] result=readTile(filename, tileNum);
		byte[] bases=result[0];
		byte[] quals=result[1];

		System.out.println("Decoded "+bases.length+" clusters for tile "+tileNum);
		System.out.println("First 5 clusters:");
		for(int i=0; i<Math.min(5, bases.length); i++){
			System.out.printf("Cluster %d: base=%c qual=%c (Q%d)\n",
				i, bases[i], quals[i], (quals[i]-33));
		}
	}
}
