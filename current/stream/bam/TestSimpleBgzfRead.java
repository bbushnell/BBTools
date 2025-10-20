package stream.bam;

import java.io.FileInputStream;
import java.io.IOException;

/**
 * Simple test to compare single-threaded vs multithreaded BGZF reading.
 */
public class TestSimpleBgzfRead {
	public static void main(String[] args) throws IOException {
		String file = args.length > 0 ? args[0] : "mapped.bam";

		// Test single-threaded
		long bytesS;
		try (FileInputStream fis = new FileInputStream(file);
		     BgzfInputStream in = new BgzfInputStream(fis)) {
			bytesS = countBytes(in);
		}

		// Test multithreaded
		long bytesMT;
		try (FileInputStream fis = new FileInputStream(file);
		     BgzfInputStreamMT in = new BgzfInputStreamMT(fis, 1)) {
			bytesMT = countBytes(in);
		}

		System.out.println("Single-threaded: " + bytesS + " bytes");
		System.out.println("Multithreaded:   " + bytesMT + " bytes");

		if (bytesS == bytesMT) {
			System.out.println("✅ Match!");
		} else {
			System.out.println("❌ MISMATCH! Difference: " + (bytesS - bytesMT));
		}
	}

	private static long countBytes(java.io.InputStream in) throws IOException {
		byte[] buf = new byte[8192];
		long total = 0;
		int n;
		while ((n = in.read(buf)) >= 0) {
			total += n;
		}
		return total;
	}
}
