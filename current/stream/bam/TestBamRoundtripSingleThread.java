package stream.bam;

import structures.ListNum;
import fileIO.ByteStreamWriter;
import stream.BamLineStreamer;
import stream.SamLine;
import stream.SamStreamer;
import stream.SamReadInputStream;

public class TestBamRoundtripSingleThread {
	public static void main(String[] args) throws Exception {
		String samFile = "mapped.sam";
		String bamFile = "mapped_st.bam";
		String outFile = "mapped_st_roundtrip.sam";

		System.err.println("=== Phase 1: Convert SAM to BAM (single thread) ===");
		SamStreamer sls = SamStreamer.makeStreamer(samFile, 1, true, true, -1, false);
		sls.start();
		java.util.ArrayList<byte[]> header = SamReadInputStream.getSharedHeader(true);

		BamWriter bw = new BamWriter(bamFile);
		bw.writeHeader(header);

		long count = 0;
		for (ListNum<SamLine> list = sls.nextLines(); list != null; list = sls.nextLines()) {
			bw.writeLines(list.list);
			count += list.size();
		}
		bw.close();
		System.err.println("Wrote " + count + " alignments");

		System.err.println("\n=== Phase 2: Read BAM back (single thread) ===");
		BamLineStreamer bls = new BamLineStreamer(bamFile, 1, true, true, -1, false);
		bls.start();

		ByteStreamWriter bsw = new ByteStreamWriter(outFile, true, false, false);
		bsw.start();

		long count2 = 0;
		for (ListNum<SamLine> list = bls.nextLines(); list != null; list = bls.nextLines()) {
			for (SamLine sl : list.list) {
				bsw.println(sl.toText());
				count2++;
			}
		}
		bsw.poisonAndWait();
		System.err.println("Read " + count2 + " alignments");
	}
}
