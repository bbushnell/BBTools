package stream.bam;

import java.io.FileInputStream;
import java.io.IOException;

/**
 * Compares total decoded byte counts from BgzfInputStream and BgzfInputStreamMT.
 * Does not compare byte contents; equal lengths can conceal different decoded data.
 */
public class TestSimpleBgzfRead{
	/**
	 * Reads the input twice, closing each reader and file stream after its pass.
	 * Constructs BgzfInputStreamMT with thread argument one. Prints both totals
	 * and a match or mismatch message; a mismatch alone does not throw or exit.
	 * @param args Optional input path, defaulting to mapped.bam; later arguments ignored
	 * @throws IOException If opening, reading or closing a stream fails
	 */
	public static void main(String[] args) throws IOException{
		String file=args.length>0 ? args[0] : "mapped.bam";

		// Test single-threaded
		long bytesS;
		try(FileInputStream fis=new FileInputStream(file);
			BgzfInputStream in=new BgzfInputStream(fis)){
			bytesS=countBytes(in);
		}

		// Test multithreaded
		long bytesMT;
		try(FileInputStream fis=new FileInputStream(file);
			BgzfInputStreamMT in=new BgzfInputStreamMT(fis, 1)){
			bytesMT=countBytes(in);
		}

		System.out.println("Single-threaded: "+bytesS+" bytes");
		System.out.println("Multithreaded:   "+bytesMT+" bytes");

		if(bytesS==bytesMT){System.out.println("✅ Match!");}else{
			System.out.println("❌ MISMATCH! Difference: "+(bytesS-bytesMT));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Consumes the input to EOF and sums returned byte counts using an 8192-byte buffer.
	 * Zero-byte reads add nothing and are retried; this helper does not close the input.
	 * @param in Open input to consume
	 * @return Total bytes reported by nonnegative reads
	 * @throws IOException If reading fails
	 */
	private static long countBytes(java.io.InputStream in) throws IOException{
		byte[] buf=new byte[8192];
		long total=0;
		int n;
		while((n=in.read(buf))>=0){
			total+=n;
		}
		return total;
	}
}
