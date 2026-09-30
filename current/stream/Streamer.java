package stream;

import structures.ListNum;

/** Common batch interface for sequence readers with implementation-specific threading.
 * Some readers perform parsing in the calling thread; others use background input
 * or conversion workers. Configure before start and consult the concrete reader for
 * restart, shutdown, counter units, sampling and error-reporting guarantees.
 * BBTools reads= limits conventionally count fragments (single reads or pairs),
 * while displayed read totals count individual reads. These concepts must remain
 * distinct when an implementation's counters are used to prepare output statistics.
 * @author Brian Bushnell, Shinobu
 * @date October 30, 2025
 */
public interface Streamer{

	/** Returns the reader's input name.
	 * @return Filename or standard-input name reported by the implementation
	 */
	public String fname();

	/** Initializes input processing, opening input or starting workers as implemented.
	 * Call before consuming batches. Repeated-start and restart behavior is not uniform;
	 * implementations may ignore a repeated call or require a fresh reader instance.
	 */
	public void start();

	/** Requests input closure, including early/emergency shutdown where supported.
	 * Prefer normal consumption to terminal input. This interface does not guarantee
	 * worker joins, final counter aggregation or identical post-close behavior across
	 * readers; backend closure can itself wait. Consult the implementation's contract.
	 */
	public void close();

	/** Reports whether the reader processes linked paired/interleaved input.
	 * This is a reader-level mode, not a substitute for checking an individual Read's mate.
	 * @return Pairing mode reported by the implementation
	 */
	public boolean paired();

	/** Returns the configured input-side pair marker.
	 * Paired/interleaved streams generally report zero while assigning markers to each mate.
	 * @return 0 for R1 or paired input, 1 for separate R2 input
	 */
	public int pairnum();

	/** Returns the implementation's reported processing count.
	 * Units currently differ: FastqStreamer counts converted individual reads, including
	 * both mates, while FastqScanStreamer divides its record count by two when interleaved.
	 * Sampling inclusion, aggregation timing and visibility also depend on the reader.
	 * Do not assume this value can be printed as an individual-read total without conversion.
	 * @return Current implementation-specific read or paired-entry count
	 */
	public long readsProcessed();

	/** Returns the implementation's reported base count.
	 * Sampling inclusion, aggregation timing and visibility depend on the reader;
	 * terminal output does not universally imply that background aggregation has finished.
	 * @return Current reported base count
	 */
	public long basesProcessed();

	/** Configures subsampling before input processing begins.
	 * Algorithms and restrictions vary: some readers use sampleKeep, others use a PRNG,
	 * and some reject fractional sampling of interleaved input. Equal seeds alone do not
	 * establish matching subsets across different implementations. No uniform rate
	 * validation or concurrent-reconfiguration guarantee is supplied by this interface.
	 * @param rate Requested retained fraction, normally in [0,1]
	 * @param seed Sampling seed, with interpretation defined by the concrete reader
	 */
	public void setSampleRate(final float rate, final long seed);

	/**
	 * Returns the next read batch under the implementation's ordering policy.
	 * May block for input. Empty nonterminal batches are possible, for example after
	 * FastqStreamer sampling; use null to detect terminal input and follow the reader's
	 * error-reporting contract to distinguish successful exhaustion from failure.
	 *
	 * The existing source audit dated 2026-09-04 recorded safe concurrent nextList calls
	 * across every implementation then reachable via StreamerFactory: host-driven classes
	 * synchronized read/advance, and worker-driven classes used queue handoffs with
	 * terminal/poison reinjection so each caller could observe EOF. This retains that
	 * historical audit, not a new certification of every implementation or of concurrent
	 * lifecycle/configuration changes. See template/A_SampleStreamerMT.java for the
	 * intended multiple-consumer pattern and the concrete reader for current restrictions.
	 * @return Next batch, possibly empty, or null at terminal input
	 */
	public ListNum<Read> nextList();

	/** Returns the next SamLine batch when the reader supports that representation.
	 * Unsupported readers may return null or throw a runtime exception; do not assume
	 * every unsupported format uses the same exception subtype.
	 * @return Next SamLine batch, or null at terminal input or for an unsupported reader
	 * @throws RuntimeException If the implementation rejects this representation
	 */
	public ListNum<SamLine> nextLines();

	/** Returns a non-consuming availability hint for pre-allocation and similar uses.
	 * False positives are possible, including after forced shutdown in queue-backed
	 * readers. Use nextList's terminal result rather than this hint to control consumption.
	 * @return Whether the implementation currently reports that more data may be available
	 */
	public boolean hasMore();

	/** Reports errors detected by the implementation so far.
	 * This query is not a universal join or final-status barrier. Some readers also fail
	 * directly from their input/consumption paths; absence of this flag is not validation
	 * of every record, especially records skipped before decoding.
	 * @return Currently reported error state
	 */
	public boolean errorState();

	//TODO:  Remove these eventually
	/** Compatibility hook whose default implementation does nothing.
	 * No recycling or ownership transfer is performed by this default method.
	 * @param ln Batch offered to an implementation that overrides the hook
	 */
	public default void returnList(final ListNum<Read> ln){}
	/** Compatibility hook whose default implementation ignores both arguments.
	 * @param id Batch identifier for an implementation that overrides the hook
	 * @param b Legacy flag interpreted only by an overriding implementation
	 */
	public default void returnList(final long id, final boolean b){}

	/**
	 * Deterministic positional subsampler: the keep/drop decision for the record at position
	 * recordNum (0-based in the file) under the given seed and rate. A pure function of its
	 * arguments (SplitMix64-style mix), so it is thread-safe with no shared state, reproducible
	 * across runs and thread counts, and mate-safe: R1 and R2 streamers sampling the same
	 * positions with the same seed and rate keep exactly the same subset. Replaces the shared-PRNG
	 * scheme, whose keep-decisions occurred in worker-scheduling order: under MT that desynced
	 * twin files, broke mate pairing, and killed the consumer on PairStreamer's numericID
	 * assert, hanging the JVM (replicated via stream.sh samplerate=0.5 threadsin=4, 2026-09-05).
	 * This helper performs no rate validation; it compares the top 24 mixed bits with
	 * the integer threshold obtained from the rate scaled by 2^24.
	 * @param recordNum Input position, normally zero-based; use pair positions for paired sampling
	 * @param seed Shared seed for calls that must make the same positional decision
	 * @param rate Retained fraction, normally in [0,1]; zero keeps none and one keeps all
	 * @return true when the mixed position falls below the scaled-rate threshold
	 */
	public static boolean sampleKeep(final long recordNum, final long seed, final float rate){
		long x=recordNum*0x9E3779B97F4A7C15L+seed;
		x=(x^(x>>>30))*0xBF58476D1CE4E5B9L;
		x=(x^(x>>>27))*0x94D049BB133111EBL;
		x^=(x>>>31);
		return (x>>>40)<(long)(rate*0x1p24f);//top 24 bits vs rate scaled to 2^24
	}

	/** Resolves a user-supplied seed; a negative argument requests a fresh random long.
	 * Resolve once and reuse the result directly with sampleKeep. Re-resolving a negative
	 * result requests another random seed, including through setters that call this helper;
	 * do not assume those setters preserve a shared negative argument.
	 * @param seed User seed, returned unchanged when nonnegative
	 * @return Supplied nonnegative seed, or a random long that may itself be negative
	 */
	public static long resolveSampleSeed(final long seed){
		return seed>=0 ? seed : new java.util.Random().nextLong();
	}

}
