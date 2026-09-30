package stream;

import java.util.ArrayList;

import structures.ListNum;

/**
 * Common submission and completion API for sequence-file writers.
 * Implementations differ in ordering, producer concurrency and where formatting occurs.
 * Some use ordered queues and worker threads; others format on the calling thread or
 * ignore batch IDs. Follow the concrete writer's configuration and ownership contract.
 *
 * Historical verification note, 2026-09-04: the then-factory-reachable implementations
 * were reported to accept concurrent producers in any arrival order through
 * OrderedQueueSystem2 or JobQueue, which reordered by dense IDs. That observation is
 * not a universal guarantee for every implementation or configuration. In particular,
 * directly constructed ST/ZT writers are not all safe for multiple producers, and
 * unordered writers need not use IDs for ordering.
 *
 * @author Brian Bushnell
 * @date October 30, 2025
 */
public interface Writer{

	/** Prepares this implementation for submissions.
	 * May start background threads or only initialize local state. Permitted call order,
	 * implicit startup and repeated-start behavior depend on the concrete writer.
	 */
	public void start();

//	/** Emergency shutdown - prefer poisonAndWait() for clean termination */
//	@Deprecated
//	public void close();

	/** Reports the implementation's individual-read output count.
	 * Counters may advance during formatting or after writes; they are not proof of durable
	 * output or successful finalization. Consult the concrete writer for observation timing.
	 * @return Current reported individual-read count
	 */
	public long readsWritten();

	/** Reports the implementation's output base count.
	 * Observation timing follows the concrete writer; this is not a completion barrier.
	 * @return Current reported base count
	 */
	public long basesWritten();

	/** Submits reads using the concrete writer's ordering and concurrency policy.
	 * Ordered queue implementations require IDs forming the expected dense ascending
	 * sequence, even when producer arrival order differs. Follow their starting-ID and
	 * submission rules. Unordered implementations may ignore IDs entirely.
	 * Keep the list and Read payloads stable until the implementation finishes using them:
	 * some copy formatted bytes before returning, while others retain input references.
	 * Submission may block for queue capacity or perform output on the calling thread.
	 * The historical factory verification is recorded in the class documentation.
	 * @param reads Batch to submit; null handling depends on the implementation
	 * @param id Batch ID interpreted according to the concrete writer's ordering contract
	 */
	public void add(ArrayList<Read> reads, long id);

	/** Submits a numbered batch under the same ordering and ownership rules as add.
	 * Producer concurrency, null handling and use of wrapper flags depend on the concrete
	 * writer. A wrapper ID alone does not make an unordered writer preserve that order.
	 * @param reads Batch wrapper to submit
	 * @see #add(ArrayList, long)
	 */
	public void addReads(ListNum<Read> reads);

	/** Submits SAM records for output or conversion supported by the concrete writer.
	 * SAM/BAM writers may serialize these directly; some FASTA/FASTQ writers convert them
	 * to Reads. Other implementations may reject this operation. Ownership, null handling
	 * and preservation of SAM metadata depend on the implementation.
	 * @param lines SAM record batch to submit
	 */
	public void addLines(ListNum<SamLine> lines);

	/** Signals that no more submissions are intended through the implementation's normal path.
	 * Whether this waits for producers, queues a marker or performs finalization varies.
	 * Coordinate with outstanding producers according to the concrete writer's contract.
	 */
	public void poison();

	/** Waits or finalizes according to the concrete implementation's completion policy.
	 * Some implementations also initiate termination; others require poison first.
	 * Interruption/error behavior and whether a return establishes completion are
	 * implementation-specific. Use the documented normal completion sequence.
	 * @return Reported error state; true indicates an error
	 */
	public boolean waitForFinish();

	/** Provides the implementation's normal termination-and-wait convenience operation.
	 * Follow its lifecycle preconditions and inspect the returned error state.
	 * @return Reported error state; true indicates an error
	 */
	public boolean poisonAndWait();

	/** Requests abandonment after an external caller failure instead of draining normal output.
	 * This operation must return promptly and must not wait for the pending backlog to be
	 * written. Queued work may be abandoned. It must set errorState and be idempotent,
	 * including calls after normal finalization. A return is not a success indication.
	 *
	 * Historical rationale, 2026-09-03: this method was added after SamWriter worker/output
	 * failure tests through stream.sh were reported to show internal OQS2 recovery without
	 * hanging. Those tests concerned internal failures; this API exists for the distinct
	 * case of a caller abandoning an otherwise healthy writer after an external failure.
	 * The dated observation does not certify every current implementation's recovery.
	 */
	public void finishError();

	/** Returns the implementation's output name or diagnostic identifier.
	 * @return Output name or identifier
	 */
	public String fname();

//	//Can never be unset
//	public void setErrorState(boolean b);

	/** Reports the implementation's current error indication.
	 * This query alone does not establish that output or cleanup has completed.
	 * @return true when the writer reports an error
	 */
	public boolean errorState();

	/** Reports the implementation's current success indication.
	 * Criteria may involve completion flags and cached errors; use the concrete writer's
	 * finalization contract before treating this as a result for the whole output.
	 * @return Whether the implementation currently reports success
	 */
	public boolean finishedSuccessfully();

}
