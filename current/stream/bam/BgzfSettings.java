package stream.bam;

import shared.Shared;
import shared.Tools;

/**
 * Mutable defaults and routing preferences for participating BGZF callers.
 * Used by BAM streams and generic native BGZF paths in ReadWrite.
 * Callers coordinate changes before opening streams; this holder supplies no
 * synchronization or automatic reconfiguration of existing streams.
 * Explicit stream constructors may bypass these defaults or routing switches.
 * Thread defaults are evaluated once at class initialization, not recomputed
 * when Shared's thread count changes.
 */
public final class BgzfSettings{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Prevents construction of this static configuration holder. */
	private BgzfSettings(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Enables multithreaded selection in participating native BGZF factories.
	 * BamOutputStream additionally requires more than one requested thread.
	 * Direct construction of a multithreaded stream need not consult this flag.
	 */
	public static boolean USE_MULTITHREADED_BGZF=true;

	/**
	 * Selects the generic ReadWrite native multithreaded output engine:
	 * true uses BgzfOutputStreamMT2 (OrderedQueueSystem2), false uses
	 * BgzfOutputStreamMT (its own queues). Only consulted on that MT output path.
	 * BamOutputStream's MT branch uses BgzfOutputStreamMT independently of this preference.
	 * ReadWrite's native MT input branch uses BgzfInputStreamMT2 independently too.
	 */
	public static boolean USE_BGZFOS_MT2=true;

	/**
	 * Default decompression worker count, initially Shared.threads() bounded to 1–8.
	 * ReadWrite applies its own bounds and input-dependent selection; the direct
	 * BgzfInputStreamMT2 convenience constructor reads this without checking
	 * USE_MULTITHREADED_BGZF. Explicit thread arguments may bypass this default.
	 */
	public static int READ_THREADS=Tools.mid(1, Shared.threads(), 8);
	// Historical tuning note: peaks at 20; not remeasured, and not the default cap.

	/**
	 * Default compression workers for the BamOutputStream convenience constructor,
	 * initially Shared.threads() bounded to 1–32. Explicit arguments and ReadWrite's
	 * derived thread counts may bypass this default.
	 */
	public static int WRITE_THREADS=Tools.mid(1, Shared.threads(), 32);

	/** Requested uncompressed block size for BamOutputStream's MT writer. */
	public static int WRITE_BLOCK_SIZE=BgzfOutputStreamMT.DEFAULT_BLOCK_SIZE;

	/** Default compression level (0–9) for the BamOutputStream convenience constructor. */
	public static int WRITE_COMPRESSION_LEVEL=6;
}
