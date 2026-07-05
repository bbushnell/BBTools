package illumina;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Reads Illumina .locs files containing cluster X,Y positions.
 * Format: 12-byte header, then 8 bytes per cluster (2 floats).
 *
 * @author Chloe
 * @date October 15, 2025
 */
public class LocsReader {

	/**
	 * Read cluster positions from .locs file.
	 * @param fname Path to s.locs file
	 * @return Array of [x,y] float pairs indexed by cluster
	 */
	public static float[][] readPositions(String fname) throws IOException {
		FileInputStream fis=new FileInputStream(fname);

		//Read header (12 bytes)
		byte[] headerBytes=new byte[12];
		//TODO: Possible bug [illumina/LocsReader#001] - return value of fis.read(headerBytes) is IGNORED; FileInputStream.read(byte[])
		//may return fewer than 12 bytes even on a valid file (not guaranteed to fill), leaving trailing 0s -> numClusters parsed from
		//garbage. The per-cluster reads below ARE checked (read!=8 throws); the header read is not. Latent-LOW (local files usually fill).
		//Fix: loop until 12 bytes read, or throw on short read (as CbclHeader.readHeader does).
		fis.read(headerBytes);

		//Last 4 bytes = cluster count (little-endian)
		ByteBuffer bb=ByteBuffer.wrap(headerBytes, 8, 4);
		bb.order(ByteOrder.LITTLE_ENDIAN);
		int numClusters=bb.getInt();

		System.err.println("Locs file: " + numClusters + " clusters");

		//Read positions (8 bytes per cluster: 2 floats)
		float[][] positions=new float[numClusters][2];
		byte[] posBytes=new byte[8];

		for(int i=0; i<numClusters; i++){
			int read=fis.read(posBytes);
			if(read!=8){
				throw new IOException("Unexpected end of file at cluster " + i);
			}

			ByteBuffer pb=ByteBuffer.wrap(posBytes);
			pb.order(ByteOrder.LITTLE_ENDIAN);
			positions[i][0]=pb.getFloat(); //X
			positions[i][1]=pb.getFloat(); //Y
		}

		fis.close();
		return positions;
	}

	/**
	 * Test reading a .locs file.
	 */
	public static void main(String[] args) throws Exception {
		if(args.length<1){
			System.err.println("Usage: LocsReader <locs_file>");
			System.exit(1);
		}

		float[][] positions=readPositions(args[0]);
		System.out.println("Read " + positions.length + " cluster positions");
		System.out.println("First 5 clusters:");
		for(int i=0; i<Math.min(5, positions.length); i++){
			System.out.printf("Cluster %d: X=%.2f, Y=%.2f\n", i, positions[i][0], positions[i][1]);
		}
	}
}
