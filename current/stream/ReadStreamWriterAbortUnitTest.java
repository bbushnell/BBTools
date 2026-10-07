package stream;

import java.util.ArrayList;

import fileIO.FileFormat;

/** Assertion-based diagnostic for abortNow() with a full, unstarted one-slot writer queue.
 * Requires -ea for the timing, error-state and termination checks; the printed PASS
 * alone does not establish those conditions when assertions are disabled. */
public final class ReadStreamWriterAbortUnitTest{

	/** Queues an empty list for /dev/null, measures abortNow(), then starts and joins the writer.
	 * Asserts measured abort time below 1000 ms, errorState, and termination after a 2000 ms join.
	 * Arguments are ignored. This method must be explicitly run to produce runtime evidence.
	 * @throws Exception If setup or the declared writer/thread operations fail */
	public static void main(String[] args) throws Exception{
		FileFormat ff=FileFormat.testOutput("/dev/null", FileFormat.FASTQ, null, true, true, false, true);
		ReadStreamByteWriter writer=new ReadStreamByteWriter(ff, null, true, 1, null, false);
		writer.addList(new ArrayList<Read>(0)); //Fill the one-slot queue; no consumer is started yet.

		final long start=System.nanoTime();
		writer.abortNow();
		final long elapsed=(System.nanoTime()-start)/1000000L;
		assert(elapsed<1000) : "abortNow took "+elapsed+" ms";
		assert(writer.errorState());

		writer.start();
		writer.join(2000);
		assert(!writer.isAlive()) : "aborted writer did not terminate";
		System.out.println("PASS ReadStreamWriterAbortUnitTest: abortNow returned in "+elapsed+" ms and writer terminated.");
	}
}
