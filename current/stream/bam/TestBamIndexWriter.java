package stream.bam;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;

/**
 * Regression test for BamIndexWriter: verifies that a generated BAI file
 * matches structural expectations and is internally consistent.
 *
 * Usage:
 *   java stream.bam.TestBamIndexWriter [bamFile] [indexFile]
 *
 * If indexFile is omitted, a temporary path is used.
 */
public final class TestBamIndexWriter {

	public static void main(String[] args) throws Exception {
		String bamPath = args.length > 0 ? args[0] : "../Chloe/mapped.bam";
		String indexPath = args.length > 1 ? args[1] : createTemporaryIndexPath();

		try {
			BamIndexWriter.writeIndex(bamPath, indexPath);
			validateIndex(indexPath);
			System.out.println("BamIndexWriter test passed for " + bamPath);
		} finally {
			if (args.length <= 1) {
				deleteQuietly(indexPath);
			}
		}
	}

	private static String createTemporaryIndexPath() throws IOException {
		File temp = File.createTempFile("bam-index-", ".bai");
		String path = temp.getAbsolutePath();
		if (!temp.delete()) {
			throw new IOException("Failed to prepare temporary index file: " + path);
		}
		return path;
	}

	private static void deleteQuietly(String path) {
		try {
			Files.deleteIfExists(new File(path).toPath());
		} catch (IOException ignored) {
			// Best effort clean up
		}
	}

	private static void validateIndex(String indexPath) throws IOException {
		byte[] data = Files.readAllBytes(new File(indexPath).toPath());
		if (data.length < 16) {
			throw new IOException("BAI file too short: " + indexPath);
		}

		ByteBuffer bb = ByteBuffer.wrap(data).order(ByteOrder.LITTLE_ENDIAN);

		int magic = bb.getInt();
		int expectedMagic = ('B') | ('A' << 8) | ('I' << 16) | (1 << 24);
		if (magic != expectedMagic) {
			throw new IOException("Invalid BAI magic number");
		}

		int nRef = bb.getInt();
		if (nRef < 0) {
			throw new IOException("Negative reference count in index");
		}

		boolean sawIndexedContent = false;
		boolean sawPseudoBin = false;

		for (int ref = 0; ref < nRef; ref++) {
			long nBin = Integer.toUnsignedLong(bb.getInt());
			for (long b = 0; b < nBin; b++) {
				int bin = bb.getInt();
				long nChunk = Integer.toUnsignedLong(bb.getInt());
				if (nChunk > 0) {
					sawIndexedContent = true;
				}
				for (long c = 0; c < nChunk; c++) {
					long beg = bb.getLong();
					long end = bb.getLong();
					if (beg < 0 || end < beg) {
						throw new IOException("Invalid chunk offsets for bin " + bin);
					}
				}
				if (bin == 37450) {
					sawPseudoBin = true;
					long mapped = bb.getLong();
					long unmapped = bb.getLong();
					if (mapped < 0 || unmapped < 0) {
						throw new IOException("Negative read totals in pseudo-bin");
					}
				}
			}

			long nIntv = Integer.toUnsignedLong(bb.getInt());
			for (long i = 0; i < nIntv; i++) {
				long offset = bb.getLong();
				if (offset < 0) {
					throw new IOException("Negative linear index offset encountered");
				}
			}
		}

		long nNoCoord = bb.getLong();
		if (nNoCoord < 0) {
			throw new IOException("Negative n_no_coor value in index");
		}

		if (bb.hasRemaining()) {
			throw new IOException("Unexpected trailing bytes in index");
		}

		if (!sawIndexedContent) {
			throw new IOException("Generated index contains no chunks");
		}
		if (!sawPseudoBin) {
			throw new IOException("Generated index missing pseudo-bin 37450");
		}
	}
}
