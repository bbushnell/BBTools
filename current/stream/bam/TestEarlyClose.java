package stream.bam;

import java.io.FileInputStream;

/**
 * Timing diagnostic that opens a BGZF MT reader, waits 150 ms, then closes it without reading.
 * The sleep does not prove queues are full. The printed completion line only reports that
 * the two close calls returned; no payload or helper-termination assertion is made.
 */
public class TestEarlyClose{
	//Historical purpose: exercise early closure while producer/workers may have queued work.
	//This remains a main-only timing diagnostic, not a general shutdown correctness oracle.

	/** Opens the positional input with an optional worker count (default one), sleeps, and times closure.
	 * Missing filename exits 1; no input bytes are read by the calling thread.
	 * @throws Exception If opening, sleeping or closing through the declared APIs fails */
	public static void main(String[] args) throws Exception{
		if(args.length<1){
			System.err.println("Usage: java stream.bam.TestEarlyClose <bgzf-file> [threads]");
			System.exit(1);
		}
		final String path=args[0];
		final int threads=(args.length>1 ? Integer.parseInt(args[1]) : 1);

		FileInputStream fis=new FileInputStream(path);
		BgzfInputStreamMT bgzf=new BgzfInputStreamMT(fis, threads);
		// Do not read; let producer/workers run a bit
		Thread.sleep(150);
		long t0=System.currentTimeMillis();
		bgzf.close();
		fis.close();
		long dt=System.currentTimeMillis()-t0;
		System.out.println("Closed cleanly in "+dt+" ms");
	}
}
