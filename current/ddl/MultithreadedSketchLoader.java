package ddl;

import java.util.ArrayList;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/**
 * Utility class for loading a CardinalityTracker (sketch) from a sequence file
 * using multiple worker threads.
 *
 * Threading pattern mirrors LogLogWrapper: one Streamer feeds read batches to
 * N worker threads, each holding its own per-thread tracker.  After all threads
 * finish, the per-thread trackers are merged into a single result via add().
 *
 * When maxThreads == 1 the Streamer + thread overhead is skipped and sequences
 * are hashed directly on the calling thread.
 *
 * @author Brian Bushnell, Noire
 * @date May 2026
 */
public class MultithreadedSketchLoader {

	/*--------------------------------------------------------------*/
	/*----------------           Public API          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads a sequence file and returns a populated CardinalityTracker.
	 *
	 * @param fname       Path to the input sequence file (FASTA or FASTQ).
	 * @param trackerType Tracker type string accepted by CardinalityTracker.makeTracker()
	 *                    (e.g. "DDL", "DynamicDemiLog", "BBLog").  Null uses the
	 *                    current default (Parser.loglogType).
	 * @param buckets     Number of buckets; will be rounded up to the next power of 2.
	 * @param k           K-mer length for hashing.
	 * @param seed        Hash seed.  Use 0 for the library default, 12345L to match
	 *                    DDLCompare's current behaviour.
	 * @param minProb     Minimum per-base probability for quality filtering; 0 disables.
	 * @param makeCounts  If true and tracker is DDL, allocate count/GC arrays for
	 *                    composition analysis.
	 * @param maxThreads  Maximum worker threads.  1 skips thread overhead entirely.
	 * @return A fully populated CardinalityTracker ready for cardinality() or comparison.
	 */
	public static CardinalityTracker loadTrackerFromSequence(
			final String fname,
			final String trackerType,
			final int buckets,
			final int k,
			final long seed,
			final float minProb,
			final boolean makeCounts,
			final int maxThreads) {

		final FileFormat ff = FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);

		if (maxThreads <= 1) {
			return loadSingleThreaded(ff, trackerType, buckets, k, seed, minProb, makeCounts);
		}
		return loadMultiThreaded(ff, trackerType, buckets, k, seed, minProb, makeCounts, maxThreads);
	}

	/**
	 * Convenience overload: seed=0, minProb=0, makeCounts=false.
	 */
	public static CardinalityTracker loadTrackerFromSequence(
			final String fname,
			final String trackerType,
			final int buckets,
			final int maxThreads) {
		return loadTrackerFromSequence(fname, trackerType, buckets, 31, 0L, 0f, false, maxThreads);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Implementation         ----------------*/
	/*--------------------------------------------------------------*/

	/** Single-threaded path: avoids Streamer + Thread overhead for maxThreads=1. */
	private static CardinalityTracker loadSingleThreaded(
			final FileFormat ff,
			final String trackerType,
			final int buckets,
			final int k,
			final long seed,
			final float minProb,
			final boolean makeCounts) {

		final CardinalityTracker tracker = createTracker(trackerType, buckets, k, seed, minProb, makeCounts);

		final Streamer cris = StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		cris.start();

		ListNum<Read> ln = cris.nextList();
		ArrayList<Read> reads = (ln != null ? ln.list : null);
		while (ln != null && reads != null && reads.size() > 0) {
			for (final Read r : reads) {
				tracker.hash(r);
			}
			cris.returnList(ln);
			ln = cris.nextList();
			reads = (ln != null ? ln.list : null);
		}
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);
		return tracker;
	}

	/** Multi-threaded path: N worker threads share one Streamer. */
	private static CardinalityTracker loadMultiThreaded(
			final FileFormat ff,
			final String trackerType,
			final int buckets,
			final int k,
			final long seed,
			final float minProb,
			final boolean makeCounts,
			final int maxThreads) {

		final CardinalityTracker merged = createTracker(trackerType, buckets, k, seed, minProb, makeCounts);

		final Streamer cris = StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		cris.start();

		final SketchThread[] workers = new SketchThread[maxThreads];
		for (int t = 0; t < maxThreads; t++) {
			workers[t] = new SketchThread(
					createTracker(trackerType, buckets, k, seed, minProb, makeCounts),
					cris);
		}
		for (final SketchThread w : workers) { w.start(); }
		for (final SketchThread w : workers) {
			while (w.getState() != Thread.State.TERMINATED) {
				try { w.join(); } catch (final InterruptedException e) { e.printStackTrace(); }
			}
			merged.add(w.tracker);
		}

		ReadWrite.closeStreams(cris);
		return merged;
	}

	/** Creates a tracker, routing through DynamicDemiLog.create when makeCounts is needed. */
	private static CardinalityTracker createTracker(
			final String trackerType, final int buckets, final int k,
			final long seed, final float minProb, final boolean makeCounts) {
		if(makeCounts && ("DDL".equalsIgnoreCase(trackerType) || "DynamicDemiLog".equalsIgnoreCase(trackerType))){
			return DynamicDemiLog.create(buckets, k, seed, minProb, true);
		}
		return CardinalityTracker.makeTracker(trackerType, buckets, k, seed, minProb);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Worker thread: pulls batches from the shared Streamer and hashes k-mers
	 *  into its own per-thread CardinalityTracker. */
	private static class SketchThread extends Thread {

		SketchThread(final CardinalityTracker tracker_, final Streamer cris_) {
			tracker = tracker_;
			cris = cris_;
		}

		@Override
		public void run() {
			ListNum<Read> ln = cris.nextList();
			ArrayList<Read> reads = (ln != null ? ln.list : null);
			while (ln != null && reads != null && reads.size() > 0) {
				for (final Read r : reads) {
					tracker.hash(r);
				}
				cris.returnList(ln);
				ln = cris.nextList();
				reads = (ln != null ? ln.list : null);
			}
			cris.returnList(ln);
		}

		/** Per-thread tracker; merged into the combined result after join(). */
		final CardinalityTracker tracker;
		private final Streamer cris;
	}

}
