package illumina;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Reads Illumina .filter files containing pass-filter flags.
 * Format: Header bytes + 1 byte per cluster (0=fail, 1=pass).
 *
 * @author Chloe
 * @date October 15, 2025
 */
public class FilterReader {

	/**
	 * Read pass-filter flags from .filter file.
	 * @param fname Path to filter file (e.g., s_1_1101.filter)
	 * @return Array of boolean flags (true=pass filter)
	 */
	public static boolean[] readFilters(String fname) throws IOException {
		FileInputStream fis=new FileInputStream(fname);

		//Read all bytes
		byte[] data=fis.readAllBytes();
		fis.close();

		//First few bytes are header - need to determine header size
		//From Picard: header is first 12 bytes, then cluster count, then flags
		ByteBuffer bb=ByteBuffer.wrap(data, 0, Math.min(data.length, 12));
		bb.order(ByteOrder.LITTLE_ENDIAN);

		//Try reading cluster count from bytes 8-12
		if(data.length>=12){
			bb.position(8);
			int numClusters=bb.getInt();
			System.err.println("Filter file: " + numClusters + " clusters (from header)");

			//Flags start after header (12 bytes)
			boolean[] filters=new boolean[numClusters];
			for(int i=0; i<numClusters && (i+12)<data.length; i++){
				filters[i]=(data[i+12]!=0);
			}
			return filters;
		}

		//Fallback: assume all bytes after first 4 are flags
		int numClusters=data.length-4;
		boolean[] filters=new boolean[numClusters];
		for(int i=0; i<numClusters; i++){
			filters[i]=(data[i+4]!=0);
		}
		return filters;
	}

	/**
	 * Test reading a .filter file.
	 */
	public static void main(String[] args) throws Exception {
		if(args.length<1){
			System.err.println("Usage: FilterReader <filter_file>");
			System.exit(1);
		}

		boolean[] filters=readFilters(args[0]);
		System.out.println("Read " + filters.length + " filter flags");

		int passing=0;
		for(boolean pf : filters){
			if(pf){passing++;}
		}
		System.out.println("Passing filter: " + passing + " (" + (100.0*passing/filters.length) + "%)");
	}
}
