package stream.bam;

import java.util.ArrayList;

import stream.BamLineStreamer;
import stream.SamLine;
import stream.SamStreamer;
import stream.Streamer;
import stream.SamReadInputStream;
import structures.ListNum;

public class TestBamRoundtrip {
	public static void main(String[] args) throws Exception {
		String samFile = "test.sam";
		String bamFile = "test2.bam";

		System.err.println("=== Phase 1: Convert SAM to BAM ===");

		// Read SAM
		Streamer sls = SamStreamer.makeStreamer(samFile, 4, true, false, -1, false);
		sls.start();

		// Wait for header to be populated in SHARED_HEADER
		java.util.ArrayList<byte[]> header = SamReadInputStream.getSharedHeader(true);

		// Write BAM
		BamWriter bw = new BamWriter(bamFile);
		bw.writeHeader(header);

		long count = 0;
		for (ListNum<SamLine> list = sls.nextLines(); list != null; list = sls.nextLines()) {
			bw.writeLines(list.list);
			count += list.size();
		}
		bw.close();
		System.err.println("Wrote " + count + " alignments to " + bamFile);

		System.err.println("\n=== Phase 2: Read BAM back ===");

		// Read BAM back
		BamLineStreamer bls = new BamLineStreamer(bamFile, 4, true, true, -1, false);
		bls.start();

		ArrayList<byte[]> header2=SamReadInputStream.getSharedHeader(true);
		System.err.println("Header lines: " + header2.size());
		for (byte[] line : header2) {
			System.err.println(new String(line));
		}

		long count2 = 0;
		for (ListNum<SamLine> list = bls.nextLines(); list != null; list = bls.nextLines()) {
			for (SamLine sl : list.list) {
				System.out.println(sl.toText());
				count2++;
			}
		}
		System.err.println("\nRead " + count2 + " alignments from " + bamFile);
	}
}
