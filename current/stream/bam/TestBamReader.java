package stream.bam;

import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/**
 * Simple test program for Streamer.
 * Usage: java -cp build stream.TestBamReader <bamfile> [maxReads]
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class TestBamReader {

	public static void main(String[] args) {
		if (args.length < 1) {
			System.err.println("Usage: java stream.TestBamReader <bamfile> [maxReads]");
			System.exit(1);
		}

		String bamFile = args[0];
		long maxReads = args.length > 1 ? Long.parseLong(args[1]) : 10;

		System.err.println("Testing Streamer on: " + bamFile);
		System.err.println("Reading up to " + maxReads + " records");

		Streamer streamer = StreamerFactory.makeSamOrBamStreamer(bamFile, 2, true, true, maxReads, false);
		streamer.start();

		long count = 0;
		ListNum<SamLine> list;

		while ((list = streamer.nextLines()) != null) {
			for (SamLine line : list) {
				if (line != null) {
					System.out.println(line.toText());
					count++;
					if (count >= maxReads) {
						break;
					}
				}
			}
			if (count >= maxReads) {
				break;
			}
		}

		System.err.println("\nTotal reads processed: " + count);
		System.err.println("Test completed successfully!");
	}
}
