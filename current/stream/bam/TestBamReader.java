package stream.bam;

import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/**
 * Prints records returned by a Streamer selected for the supplied input.
 * This is a diagnostic driver, with no comparison against expected records.
 * Usage: {@code java -ea stream.bam.TestBamReader <bamfile> [maxReads]}.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class TestBamReader{

	/**
	 * Requests ordered output and header retention, then prints nonnull SamLine records.
	 * Forwards the limit unchanged to the selected reader and also breaks the print
	 * loop once the number printed reaches that limit. Defaults to ten; no range
	 * validation is performed. The final count is printed records, not fragments.
	 * Missing input prints usage and exits with status one. Later arguments are ignored.
	 * @param args Input path followed by an optional long-valued limit
	 * @throws NumberFormatException If the supplied limit is not a long integer
	 */
	public static void main(String[] args){
		if(args.length<1){
			System.err.println("Usage: java stream.bam.TestBamReader <bamfile> [maxReads]");
			System.exit(1);
		}

		String bamFile=args[0];
		long maxReads=args.length>1 ? Long.parseLong(args[1]) : 10;

		System.err.println("Testing Streamer on: "+bamFile);
		System.err.println("Reading up to "+maxReads+" records");

		Streamer streamer=StreamerFactory.makeSamOrBamStreamer(bamFile, 2, true, true, maxReads, false);
		streamer.start();

		long count=0;
		ListNum<SamLine> list;

		while((list=streamer.nextLines())!=null){
			for(SamLine line : list){
				if(line!=null){
					System.out.println(line.toText());
					count++;
					if(count>=maxReads){break;}
				}
			}
			if(count>=maxReads){break;}
		}

		System.err.println("\nTotal reads processed: "+count);
		//TODO: Possible bug [stream/bam/TestBamReader#51] - this driver never calls
		//streamer.close(), including after a count-limited break. Cleanup depends on
		//the selected backend; no leak or shutdown outcome has been established here.
		//Keep lifecycle changes separate from this driver's documentation cleanup.
		System.err.println("Test completed successfully!");
	}
}
