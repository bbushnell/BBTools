package stream.bam;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Test program for multithreaded BGZF streams.
 *
 * Tests both BgzfInputStreamMT and BgzfOutputStreamMT by:
 * 1. Reading a BGZF file (BAM)
 * 2. Writing it to a new file
 * 3. Reading both files and comparing byte-by-byte
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class TestBgzfMultithreaded {

	public static void main(String[] args) throws IOException {
		if (args.length < 1) {
			System.err.println("Usage: java stream.bam.TestBgzfMultithreaded <input.bam> [threads]");
			System.err.println("Example: java stream.bam.TestBgzfMultithreaded mapped.bam 1");
			System.exit(1);
		}

		String inputFile = args[0];
		int threads = args.length > 1 ? Integer.parseInt(args[1]) : 1;
		String outputFile = "test_mt_output.bam";

		System.out.println("BGZF Multithreaded Round-Trip Test");
		System.out.println("===================================");
		System.out.println("Input file:   " + inputFile);
		System.out.println("Output file:  " + outputFile);
		System.out.println("Threads:      " + threads);
		System.out.println();

		// Test 1: Read with single-threaded, write with multithreaded
		System.out.println("Test 1: Single-threaded read -> Multithreaded write");
		testReadWrite(inputFile, outputFile, false, true, threads);

		// Test 2: Read with multithreaded, compare byte-by-byte
		System.out.println("\nTest 2: Multithreaded read -> Byte comparison");
		testRoundTrip(inputFile, outputFile, threads);

		// Test 3: Full multithreaded round-trip
		System.out.println("\nTest 3: Full multithreaded round-trip");
		String output2 = "test_mt_output2.bam";
		testReadWrite(outputFile, output2, true, true, threads);
		testRoundTrip(outputFile, output2, threads);

		System.out.println("\n✅ All tests passed!");
		System.out.println("Multithreaded BGZF streams working correctly with " + threads + " thread(s)");

		// Cleanup
		new File(outputFile).delete();
		new File(output2).delete();
	}

	/**
	 * Test reading and writing with specified stream types.
	 */
	private static void testReadWrite(String input, String output,
	                                   boolean mtRead, boolean mtWrite, int threads) throws IOException {
		long startTime = System.currentTimeMillis();
		long bytesRead = 0;

		try (FileInputStream fis = new FileInputStream(input);
		     FileOutputStream fos = new FileOutputStream(output)) {

			// Create input stream
			Object inStream = mtRead ?
				new BgzfInputStreamMT(fis, threads) :
				new BgzfInputStream(fis);

			// Create output stream
			Object outStream = mtWrite ?
				new BgzfOutputStreamMT(fos, threads, 6) :
				new BgzfOutputStream(fos, 6);

			// Copy data
			byte[] buffer = new byte[8192];
			int n;

			if (mtRead) {
				BgzfInputStreamMT in = (BgzfInputStreamMT)inStream;
				if (mtWrite) {
					BgzfOutputStreamMT out = (BgzfOutputStreamMT)outStream;
					while ((n = in.read(buffer)) >= 0) {
						if (n > 0) {
							out.write(buffer, 0, n);
							bytesRead += n;
						}
					}
					out.close();
				} else {
					BgzfOutputStream out = (BgzfOutputStream)outStream;
					while ((n = in.read(buffer)) >= 0) {
						if (n > 0) {
							out.write(buffer, 0, n);
							bytesRead += n;
						}
					}
					out.writeEOF();
					out.close();
				}
				in.close();
			} else {
				BgzfInputStream in = (BgzfInputStream)inStream;
			if (mtWrite) {
				BgzfOutputStreamMT out = (BgzfOutputStreamMT)outStream;
				while ((n = in.read(buffer)) >= 0) {
					if (n > 0) {
						out.write(buffer, 0, n);
						bytesRead += n;
					}
				}
				out.close();
			} else {
					BgzfOutputStream out = (BgzfOutputStream)outStream;
					while ((n = in.read(buffer)) >= 0) {
						if (n > 0) {
							out.write(buffer, 0, n);
							bytesRead += n;
						}
					}
					out.writeEOF();
					out.close();
				}
				in.close();
			}
		}

		long elapsed = System.currentTimeMillis() - startTime;
		double mbps = (bytesRead / 1024.0 / 1024.0) / (elapsed / 1000.0);

		System.out.println("  Bytes transferred: " + bytesRead);
		System.out.println("  Time:             " + elapsed + " ms");
		System.out.println("  Throughput:       " + String.format("%.2f", mbps) + " MB/s");
	}

	/**
	 * Test round-trip correctness by comparing decompressed data.
	 */
	private static void testRoundTrip(String file1, String file2, int threads) throws IOException {
		long startTime = System.currentTimeMillis();

		try (FileInputStream fis1 = new FileInputStream(file1);
		     FileInputStream fis2 = new FileInputStream(file2);
		     BgzfInputStreamMT in1 = new BgzfInputStreamMT(fis1, threads);
		     BgzfInputStreamMT in2 = new BgzfInputStreamMT(fis2, threads)) {

			byte[] buf1 = new byte[8192];
			byte[] buf2 = new byte[8192];
			long totalBytes = 0;
			int blockNum = 0;

			while (true) {
				int n1 = in1.read(buf1);
				int n2 = in2.read(buf2);

				if (n1 != n2) {
					throw new IOException("Read size mismatch at block " + blockNum +
						": file1=" + n1 + ", file2=" + n2);
				}

				if (n1 < 0) {
					break; // EOF
				}

				// Compare bytes
				for (int i = 0; i < n1; i++) {
					if (buf1[i] != buf2[i]) {
						throw new IOException("Byte mismatch at block " + blockNum +
							", offset " + i + ": file1=" + (buf1[i] & 0xFF) +
							", file2=" + (buf2[i] & 0xFF));
					}
				}

				totalBytes += n1;
				blockNum++;
			}

			long elapsed = System.currentTimeMillis() - startTime;
			System.out.println("  Bytes compared:   " + totalBytes);
			System.out.println("  Blocks compared:  " + blockNum);
			System.out.println("  Time:             " + elapsed + " ms");
			System.out.println("  ✅ Files match perfectly!");
		}
	}
}
