package stream.bam;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;

import stream.BamLineStreamer;
import stream.SamLine;
import structures.ListNum;

/**
 * Converts BAM format to SAM format.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class Bam2Sam {

	public static void main(String[] args) {
		if (args.length < 2) {
			System.err.println("Converts BAM format to SAM format.");
			System.err.println("Usage: java stream.Bam2Sam input.bam output.sam");
			System.exit(1);
		}

		String inFile = args[0];
		String outFile = args[1];

		try {
			// Create output writer
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outFile)));

			// Create BAM reader with 4 threads, save header
			BamLineStreamer bls = new BamLineStreamer(inFile, 4, true, false, -1, false);

			// Start reading
			bls.start();

			// Wait for header to be populated
			while (bls.header != null && bls.header.isEmpty()) {
				try {
					Thread.sleep(10);
				} catch (InterruptedException e) {
					break;
				}
			}

			// Make a copy of header to avoid concurrent modification
			java.util.ArrayList<byte[]> headerCopy = null;
			if (bls.header != null) {
				synchronized (bls.header) {
					headerCopy = new java.util.ArrayList<byte[]>(bls.header);
				}
			}

			// Write header lines
			if (headerCopy != null) {
				for (byte[] line : headerCopy) {
					out.println(new String(line));
				}
			}

			// Write alignment records
			long count = 0;
			for (ListNum<SamLine> list = bls.nextLines(); list != null; list = bls.nextLines()) {
				for (SamLine sl : list) {
					out.println(sl.toText());
					count++;
				}
			}

			out.close();

			System.err.println("Converted " + count + " alignments");
			System.err.println("Wrote: " + outFile);

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

}
