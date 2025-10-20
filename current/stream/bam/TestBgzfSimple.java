package stream.bam;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Simple BGZF test with known data.
 * Creates simple repeating pattern, compresses/decompresses, verifies correctness.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class TestBgzfSimple {

	public static void main(String[] args) throws IOException {
		System.out.println("BGZF Simple Data Test");
		System.out.println("=====================\n");

		// Create simple test data
		String testFile = "test_simple.txt";
		String bgzfFile = "test_simple.txt.bgz";

		createTestData(testFile);

		// Test 1: Compress with our compressor
		System.out.println("Test 1: Compressing with BgzfOutputStreamMT...");
		compressFile(testFile, bgzfFile);
		System.out.println("  Output: " + bgzfFile);
		System.out.println("  ✅ Compression complete\n");

		// Test 2: Decompress with our decompressor
		System.out.println("Test 2: Decompressing with BgzfInputStreamMT...");
		String decompressed = "test_simple_decompressed.txt";
		decompressFile(bgzfFile, decompressed);

		// Verify
		if (filesMatch(testFile, decompressed)) {
			System.out.println("  ✅ Round-trip successful! Files match perfectly.");
		} else {
			System.out.println("  ❌ MISMATCH! Files differ.");
			System.exit(1);
		}

		System.out.println("\n✅ All tests passed!");
		System.out.println("\nNow test with command-line gzip:");
		System.out.println("  gunzip -c " + bgzfFile + " > test_gzip_decompressed.txt");
		System.out.println("  diff " + testFile + " test_gzip_decompressed.txt");
	}

	private static void createTestData(String filename) throws IOException {
		System.out.println("Creating test data: " + filename);

		try (FileOutputStream fos = new FileOutputStream(filename)) {
			// Write "AAAAAAAAAAAAAAAAAAAA\n" 1000 times
			String line = "AAAAAAAAAAAAAAAAAAAA\n";
			byte[] lineBytes = line.getBytes();

			for (int i = 0; i < 1000; i++) {
				fos.write(lineBytes);
			}
		}

		System.out.println("  Created: 1000 lines of 'AAAAAAAAAAAAAAAAAAAA'");
		System.out.println("  Size: " + (21 * 1000) + " bytes\n");
	}

	private static void compressFile(String input, String output) throws IOException {
		// Enable debug mode
		System.setProperty("bgzf.debug", "true");

		try (FileInputStream fis = new FileInputStream(input);
		     FileOutputStream fos = new FileOutputStream(output);
		     BgzfOutputStreamMT bgzf = new BgzfOutputStreamMT(fos, 1, 6)) {

			byte[] buffer = new byte[8192];
			int n;
			long total = 0;

			while ((n = fis.read(buffer)) >= 0) {
				bgzf.write(buffer, 0, n);
				total += n;
			}

			System.out.println("  Compressed " + total + " bytes");
			// close() will automatically create lastJob and write EOF marker
		}
	}

	private static void decompressFile(String input, String output) throws IOException {
		try (FileInputStream fis = new FileInputStream(input);
		     BgzfInputStreamMT bgzf = new BgzfInputStreamMT(fis, 1);
		     FileOutputStream fos = new FileOutputStream(output)) {

			byte[] buffer = new byte[8192];
			int n;
			long total = 0;

			while ((n = bgzf.read(buffer)) >= 0) {
				fos.write(buffer, 0, n);
				total += n;
			}

			System.out.println("  Decompressed " + total + " bytes");
		}
	}

	private static boolean filesMatch(String file1, String file2) throws IOException {
		try (FileInputStream fis1 = new FileInputStream(file1);
		     FileInputStream fis2 = new FileInputStream(file2)) {

			byte[] buf1 = new byte[8192];
			byte[] buf2 = new byte[8192];
			long pos = 0;

			while (true) {
				int n1 = fis1.read(buf1);
				int n2 = fis2.read(buf2);

				if (n1 != n2) {
					System.err.println("  Read size mismatch at position " + pos + ": " + n1 + " vs " + n2);
					return false;
				}

				if (n1 < 0) break; // EOF

				for (int i = 0; i < n1; i++) {
					if (buf1[i] != buf2[i]) {
						System.err.println("  Byte mismatch at position " + (pos + i));
						return false;
					}
				}

				pos += n1;
			}

			return true;
		}
	}
}
