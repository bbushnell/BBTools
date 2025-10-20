package stream.bam;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Stress-test BGZF multithreaded streams using a simple repetitive dataset.
 *
 * Generates a plain-text file containing 40,000 lines of "AAAAAAAAAAAAAAAAAAAA\n"
 * (~840 KB), repeatedly compresses it with {@link BgzfOutputStreamMT}, then
 * decompresses it with {@link BgzfInputStreamMT}, verifying byte-for-byte
 * equality after each round-trip. The goal is to exercise queue ordering and
 * synchronization under repeated runs.
 *
 * Usage: {@code java stream.bam.TestBgzfA80k [iterations] [threads] [blockSize]}
 *  - iterations (optional): number of repetitions per thread count (default 20)
 *  - threads (optional): maximum thread count to test (default 4)
 *  - blockSize (optional): target BGZF block size (default 65536)
 *
 * The test always exercises thread counts from 1 up to the requested maximum.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class TestBgzfA80k {

	private static final String INPUT_FILE = "A80k.txt";
	private static final String COMPRESSED_FILE = "A40k.txt.gz";
	private static final String OUTPUT_FILE = "A80k.roundtrip.txt";
	private static final byte[] LINE_BYTES = "AAAAAAAAAAAAAAAAAAAA\n".getBytes();
	private static final int LINE_REPETITIONS = 40_000;

	public static void main(String[] args) throws IOException {
		int iterations = args.length > 0 ? Integer.parseInt(args[0]) : 20;
		int maxThreads = args.length > 1 ? Integer.parseInt(args[1]) : 4;
		int blockSize = args.length > 2 ? Integer.parseInt(args[2]) : 65536;
		if (iterations < 1) iterations = 1;
		if (maxThreads < 1) maxThreads = 1;
		if (blockSize < 256) blockSize = 256;

		System.out.println("BGZF A80k Stress Test");
		System.out.println("=====================");
		System.out.println("Iterations per thread count: " + iterations);
		System.out.println("Thread counts: 1.." + maxThreads);
		System.out.println("Block size:   " + blockSize + " bytes");
		System.out.println();

		generateInputFile();

		for (int threads = 1; threads <= maxThreads; threads++) {
			System.out.println("Thread count: " + threads);
			for (int iter = 1; iter <= iterations; iter++) {
				System.out.print("  Iteration " + iter + "/" + iterations + " ... ");
				runRoundTrip(threads, blockSize);
				System.out.println("OK");
			}
			System.out.println();
		}

		System.out.println("✅ Completed BGZF round-trip stress test.");

		// Leave artifacts for inspection, but ensure they exist.
		printFileSizes();
	}

	private static void generateInputFile() throws IOException {
		File file = new File(INPUT_FILE);
		if (file.exists()) {
			return;
		}

		try (FileOutputStream fos = new FileOutputStream(file)) {
			for (int i = 0; i < LINE_REPETITIONS; i++) {
				fos.write(LINE_BYTES);
			}
		}
	}

	private static void runRoundTrip(int threads, int blockSize) throws IOException {
		compressFile(INPUT_FILE, COMPRESSED_FILE, threads, blockSize);
		decompressFile(COMPRESSED_FILE, OUTPUT_FILE, threads);

		long mismatch = Files.mismatch(Path.of(INPUT_FILE), Path.of(OUTPUT_FILE));
		if (mismatch != -1) {
			throw new IOException("Round-trip mismatch at byte position " + mismatch);
		}
	}

	private static void compressFile(String input, String output, int threads, int blockSize) throws IOException {
		try (FileInputStream fis = new FileInputStream(input);
		     FileOutputStream fos = new FileOutputStream(output);
		     BgzfOutputStreamMT bgzf = new BgzfOutputStreamMT(fos, threads, 6, blockSize)) {

			byte[] buffer = new byte[32 * 1024];
			int n;
			while ((n = fis.read(buffer)) >= 0) {
				if (n > 0) {
					bgzf.write(buffer, 0, n);
				}
			}
			// close() writes EOF marker and flushes workers
			bgzf.close();
		}
	}

	private static void decompressFile(String input, String output, int threads) throws IOException {
		try (FileInputStream fis = new FileInputStream(input);
		     BgzfInputStreamMT bgzf = new BgzfInputStreamMT(fis, threads);
		     FileOutputStream fos = new FileOutputStream(output)) {

			byte[] buffer = new byte[32 * 1024];
			int n;
			while ((n = bgzf.read(buffer)) >= 0) {
				if (n > 0) {
					fos.write(buffer, 0, n);
				}
			}
		}
	}

	private static void printFileSizes() throws IOException {
		Path input = Path.of(INPUT_FILE);
		Path compressed = Path.of(COMPRESSED_FILE);
		Path output = Path.of(OUTPUT_FILE);

		System.out.println("Artifacts:");
		System.out.println("  " + INPUT_FILE + "       " + Files.size(input) + " bytes");
		System.out.println("  " + COMPRESSED_FILE + "  " + Files.size(compressed) + " bytes");
		System.out.println("  " + OUTPUT_FILE + "      " + Files.size(output) + " bytes");
	}
}
