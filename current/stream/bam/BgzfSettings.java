package stream.bam;

import shared.Shared;
import shared.Tools;

/**
 * Global switches for enabling multithreaded BGZF streams during BAM IO.
 * Users may toggle these flags at runtime prior to constructing readers/writers.
 */
public final class BgzfSettings {

	/** Toggle to enable multithreaded BGZF input/output. */
	public static boolean USE_MULTITHREADED_BGZF = true;

	/**
	 * Number of worker threads to use when decompressing BGZF blocks.
	 * Only consulted when {@link #USE_MULTITHREADED_BGZF} is true.
	 */
	public static int READ_THREADS = Tools.mid(1, Shared.threads(), 8);

	/** Number of worker threads to use when compressing BGZF blocks. */
	public static int WRITE_THREADS = Tools.mid(1, Shared.threads(), 32);

	/** Maximum uncompressed BGZF block size used for writers. */
	public static int WRITE_BLOCK_SIZE = BgzfOutputStreamMT.DEFAULT_BLOCK_SIZE;

	/** Compression level (0-9) used for writers. */
	public static int WRITE_COMPRESSION_LEVEL = 6;

	private BgzfSettings() {
		// Utility class
	}
}
