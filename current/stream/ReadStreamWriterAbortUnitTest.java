package stream;

import java.util.ArrayList;

import fileIO.FileFormat;

/** Verifies that abortNow() does not block when the bounded input queue is full. */
public final class ReadStreamWriterAbortUnitTest {

	public static void main(String[] args) throws Exception {
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
