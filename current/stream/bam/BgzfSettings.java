package stream.bam;

import shared.Shared;

/**
 * Global switches for enabling multithreaded BGZF streams during BAM IO.
 * Users may toggle these flags at runtime prior to constructing readers/writers.
 */
public final class BgzfSettings {

	/** Toggle to enable multithreaded BGZF input/output. */
	public static boolean USE_MULTITHREADED_BGZF = false;

	/**
	 * Number of worker threads to use when decompressing BGZF blocks.
	 * Only consulted when {@link #USE_MULTITHREADED_BGZF} is true.
	 */
	public static int READ_THREADS = Math.max(1, Shared.threads());

	/** Number of worker threads to use when compressing BGZF blocks. */
	public static int WRITE_THREADS = Math.max(1, Shared.threads());

	/** Maximum uncompressed BGZF block size used for writers. */
	public static int WRITE_BLOCK_SIZE = BgzfOutputStreamMT.DEFAULT_BLOCK_SIZE;

	/** Compression level (0-9) used for writers. */
	public static int WRITE_COMPRESSION_LEVEL = 6;

	private BgzfSettings() {
		// Utility class
	}
}
