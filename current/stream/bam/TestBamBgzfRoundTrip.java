package stream.bam;

import java.io.File;
import java.io.IOException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;

import shared.Shared;
import stream.Read;
import stream.SamLine;
import stream.SamReadInputStream;

/**
 * Utility for verifying BAM read/write compatibility with configurable BGZF modes.
 */
public final class TestBamBgzfRoundTrip {

	public static void main(String[] args) throws Exception {
		String input = args.length > 0 ? args[0] : "../Chloe/mapped.bam";
		File inputFile = new File(input);
		if (!inputFile.exists()) {
			throw new IllegalArgumentException("Input BAM not found: " + input);
		}

		String singleOut = inputFile.getName() + ".single.bam";
		String mtOut = inputFile.getName() + ".mt.bam";

		// Single-threaded run
		configure(false);
		copyBam(input, singleOut);
		String digestOriginal = digestForBam(input);
		String digestSingle = digestForBam(singleOut);
		System.out.println("Single-thread digest match: " + digestOriginal.equals(digestSingle));

		// Multithreaded run
		configure(true);
		copyBam(input, mtOut);
		String digestMulti = digestForBam(mtOut);
		System.out.println("MT digest match: " + digestOriginal.equals(digestMulti));

		// Clean up temporary files
		new File(singleOut).delete();
		new File(mtOut).delete();
	}

	private static void configure(boolean useMultithreaded) {
		BgzfSettings.USE_MULTITHREADED_BGZF = useMultithreaded;
		if (useMultithreaded) {
			int threads = Math.max(1, Shared.threads());
			BgzfSettings.READ_THREADS = threads;
			BgzfSettings.WRITE_THREADS = threads;
		} else {
			BgzfSettings.READ_THREADS = 1;
			BgzfSettings.WRITE_THREADS = 1;
		}
	}

	private static void copyBam(String input, String output) throws IOException {
		BamReadInputStreamST reader = new BamReadInputStreamST(input, true, false, true);
		reader.start();

		List<byte[]> headerCopy = null;
		ArrayList<byte[]> header = SamReadInputStream.getSharedHeader(true);
		if (header != null) {
			headerCopy = new ArrayList<byte[]>(header.size());
			for (byte[] line : header) {
				headerCopy.add(line.clone());
			}
		}

		BamWriter writer = new BamWriter(output);
		if (headerCopy != null) {
			writer.writeHeader(headerCopy);
		}

		ArrayList<Read> reads;
		while ((reads = reader.nextList()) != null && reads.size() > 0) {
			ArrayList<SamLine> lines = new ArrayList<SamLine>(reads.size());
			for (Read r : reads) {
				SamLine sl = r.samline;
				if (sl != null) {
					lines.add(sl);
				}
			}
			writer.writeLines(lines);
		}

		writer.close();
		reader.close();
	}

	private static String digestForBam(String path) throws IOException, NoSuchAlgorithmException {
		MessageDigest md = MessageDigest.getInstance("SHA-256");
		BamReadInputStreamST reader = new BamReadInputStreamST(path, true, false, true);
		reader.start();
		ArrayList<Read> reads;
		while ((reads = reader.nextList()) != null && reads.size() > 0) {
			for (Read r : reads) {
				SamLine sl = r.samline;
				if (sl != null) {
					byte[] bytes = sl.toText().toBytes();
					md.update(bytes, 0, bytes.length);
				}
			}
		}
		reader.close();
		byte[] digest = md.digest();
		StringBuilder sb = new StringBuilder(digest.length * 2);
		for (byte b : digest) {
			sb.append(String.format("%02x", b));
		}
		return sb.toString();
	}

	private TestBamBgzfRoundTrip() {
	}
}
