package stream.bam;

import stream.SamLine;

/**
 * Quick test of the BAM converter fix
 */
public class QuickTest {
	public static void main(String[] args) throws Exception {
		// Test case from evidence: read with FLAG=83 (reverse strand)
		String originalLine = "read_test\t83\tchr1\t1000\t60\t100M\t=\t1100\t250" +
			"\tTTGAGCCGCTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGACGGGC" +
			"\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

		System.out.println("Original SAM line:");
		System.out.println(originalLine);
		System.out.println();

		// Parse to SamLine
		SamLine sl = new SamLine(originalLine.getBytes());
		System.out.println("Parsed SamLine:");
		System.out.println("  QNAME: " + sl.qname);
		System.out.println("  FLAG: " + sl.flag + " (reverse=" + ((sl.flag & 0x10) != 0) + ")");
		System.out.println("  SEQ (first 20): " + new String(sl.seq, 0, Math.min(20, sl.seq.length)));
		System.out.println();

		// Convert to BAM
		String[] refNames = {"chr1"};
		SamToBamConverter toBam = new SamToBamConverter(refNames);
		byte[] bamRecord = toBam.convertAlignment(sl);
		System.out.println("BAM record size: " + bamRecord.length + " bytes");
		System.out.println();

		// Convert back to SAM
		BamToSamConverter toSam = new BamToSamConverter(refNames);
		// Skip first 4 bytes (block_size)
		byte[] bamAlignment = new byte[bamRecord.length - 4];
		System.arraycopy(bamRecord, 4, bamAlignment, 0, bamAlignment.length);
		byte[] samBytes = toSam.convertAlignment(bamAlignment);
		String roundtrip = new String(samBytes);

		System.out.println("Roundtrip SAM line:");
		System.out.println(roundtrip);
		System.out.println();

		// Extract and compare SEQ fields
		String[] origFields = originalLine.split("\t");
		String[] roundFields = roundtrip.split("\t");

		String origSeq = origFields[9];
		String roundSeq = roundFields[9];

		System.out.println("SEQ comparison:");
		System.out.println("  Original:  " + origSeq);
		System.out.println("  Roundtrip: " + roundSeq);
		System.out.println("  Match: " + origSeq.equals(roundSeq));
		System.out.println();

		if (origSeq.equals(roundSeq)) {
			System.out.println("SUCCESS: Sequence preserved correctly!");
		} else {
			System.out.println("FAILURE: Sequence differs!");
			System.out.println("  Original starts:  " + origSeq.substring(0, 20));
			System.out.println("  Roundtrip starts: " + roundSeq.substring(0, 20));
			System.out.println("  Original ends:    " + origSeq.substring(origSeq.length()-20));
			System.out.println("  Roundtrip ends:   " + roundSeq.substring(roundSeq.length()-20));
		}
	}
}
