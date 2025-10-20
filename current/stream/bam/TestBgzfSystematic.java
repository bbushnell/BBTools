package stream.bam;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Systematic test with increasing block counts.
 * Tests 1 block, 2 blocks, 4 blocks, 8 blocks of simple 'A' data.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class TestBgzfSystematic {

	private static final int BLOCK_SIZE = 65536; // 64KB

	public static void main(String[] args) throws IOException {
		System.out.println("BGZF Systematic Block Test");
		System.out.println("==========================\n");

		// Enable debug output
		System.setProperty("bgzf.debug", "true");

		int[] blockCounts = {1, 2, 4, 8};

		for (int blocks : blockCounts) {
			System.out.println("Testing " + blocks + " block(s)...");
			System.out.println("----------------------------");

			boolean success = testBlocks(blocks);

			if (!success) {
				System.out.println("❌ FAILED at " + blocks + " blocks!");
				System.exit(1);
			}

			System.out.println("✅ " + blocks + " block(s) passed\n");
		}

		System.out.println("\n✅ All tests passed!");
	}

	private static boolean testBlocks(int blockCount) throws IOException {
		String input = "test_" + blockCount + "_blocks.txt";
		String compressed = "test_" + blockCount + "_blocks.bgz";
		String decompressed = "test_" + blockCount + "_blocks_out.txt";

		// Create input file
		int totalBytes = blockCount * BLOCK_SIZE;
		createFile(input, totalBytes);
		System.out.println("  Created: " + totalBytes + " bytes (" + blockCount + " blocks)");

		// Compress
		compressFile(input, compressed);

		// Decompress
		decompressFile(compressed, decompressed);

		// Verify
		boolean match = verifyFiles(input, decompressed, totalBytes);

		return match;
	}

	private static void createFile(String filename, int totalBytes) throws IOException {
		try (FileOutputStream fos = new FileOutputStream(filename)) {
			byte[] block = new byte[8192];
			for (int i = 0; i < block.length; i++) {
				block[i] = 'A';
			}

			int remaining = totalBytes;
			while (remaining > 0) {
				int toWrite = Math.min(block.length, remaining);
				fos.write(block, 0, toWrite);
				remaining -= toWrite;
			}
		}
	}

	private static void compressFile(String input, String output) throws IOException {
		long bytesWritten = 0;

		try (FileInputStream fis = new FileInputStream(input);
		     FileOutputStream fos = new FileOutputStream(output);
		     BgzfOutputStreamMT bgzf = new BgzfOutputStreamMT(fos, 1, 6)) {

			byte[] buffer = new byte[8192];
			int n;

			while ((n = fis.read(buffer)) >= 0) {
				bgzf.write(buffer, 0, n);
				bytesWritten += n;
			}

			// Don't call writeEOF() - close() handles it automatically!
		}

		System.out.println("  Compressed: " + bytesWritten + " bytes");
	}

	private static void decompressFile(String input, String output) throws IOException {
		long bytesRead = 0;

		try (FileInputStream fis = new FileInputStream(input);
		     BgzfInputStreamMT bgzf = new BgzfInputStreamMT(fis, 1);
		     FileOutputStream fos = new FileOutputStream(output)) {

			byte[] buffer = new byte[8192];
			int n;

			while ((n = bgzf.read(buffer)) >= 0) {
				fos.write(buffer, 0, n);
				bytesRead += n;
			}
		}

		System.out.println("  Decompressed: " + bytesRead + " bytes");
	}

	private static boolean verifyFiles(String file1, String file2, int expectedSize) throws IOException {
		try (FileInputStream fis1 = new FileInputStream(file1);
		     FileInputStream fis2 = new FileInputStream(file2)) {

			byte[] buf1 = new byte[8192];
			byte[] buf2 = new byte[8192];
			long totalBytes = 0;

			while (true) {
				int n1 = fis1.read(buf1);
				int n2 = fis2.read(buf2);

				if (n1 != n2) {
					System.err.println("  ❌ Read size mismatch at position " + totalBytes);
					System.err.println("     Expected more data, got EOF");
					return false;
				}

				if (n1 < 0) break;

				for (int i = 0; i < n1; i++) {
					if (buf1[i] != buf2[i]) {
						System.err.println("  ❌ Byte mismatch at position " + (totalBytes + i));
						return false;
					}
				}

				totalBytes += n1;
			}

			if (totalBytes != expectedSize) {
				System.err.println("  ❌ Size mismatch: expected " + expectedSize + ", got " + totalBytes);
				return false;
			}

			System.out.println("  Verified: " + totalBytes + " bytes match perfectly");
			return true;
		}
	}
}
