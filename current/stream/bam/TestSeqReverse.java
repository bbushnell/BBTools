package stream.bam;

import parse.LineParser1;
import simd.Vector;
import stream.SamLine;

/**
 * Standalone diagnostic printing sequence orientation through SAM/BAM conversion.
 * Constructs one fixed FLAG=83 alignment, prints its stored sequence, and attempts
 * a round trip through the direct BAM-to-SAM text converter. The stored-sequence
 * expectation assumes SamLine.FLIP_ON_LOAD is enabled. Comparisons only print;
 * this driver is not an assertion-based automated test. Its record-framing mismatch
 * is recorded below and must be repaired before relying on round-trip diagnostics.
 */
public class TestSeqReverse{

	//Test harness: round-trips a synthetic FLAG=83 read through SamLine→BAM→SAM and compares
	//SEQ fields. All checks are print-only — no System.exit(1) on mismatch, so CI cannot detect
	//failure automatically. main()-only.
	/**
	 * Prints the fixed alignment, stored sequence, and attempted round-trip comparison.
	 * Arguments are ignored; this method does not configure the global SAM/SIMD options.
	 * A printed mismatch alone does not cause a failing process status.
	 * @param args Ignored command-line arguments
	 * @throws Exception If parsing or conversion fails
	 */
	public static void main(String[] args) throws Exception{
		// Create a SamLine manually with FLAG=83 (reverse strand, paired, mate1)
		// FLAG 83 = 0x53 = 01010011 binary
		// Bit 0 (0x1): paired
		// Bit 1 (0x2): properly paired
		// Bit 4 (0x10): reverse strand ← THIS IS THE KEY BIT
		// Bit 6 (0x40): mate 1

		String originalSeq="ACGTACGTAC";
		String qname="test_read";
		int flag=83; // Reverse strand
		String rname="chr1";
		int pos=1000;
		int mapq=60;
		String cigar="10M";
		String rnext="=";
		int pnext=1100;
		int tlen=200;
		String qual="IIIIIIIIII";

		// Build SAM line
		String samLine=String.join("\t",
			qname,
			String.valueOf(flag),
			rname,
			String.valueOf(pos),
			String.valueOf(mapq),
			cigar,
			rnext,
			String.valueOf(pnext),
			String.valueOf(tlen),
			originalSeq,
			qual
		);

		System.out.println("Original SAM line:");
		System.out.println(samLine);
		System.out.println();

		// Parse into SamLine
		SamLine sl=new SamLine(new LineParser1('\t').set(samLine.getBytes()));

		// Check what's stored in SamLine.seq
		String storedSeq=new String(sl.seq);
		System.out.println("Sequence stored in SamLine.seq:");
		System.out.println(storedSeq);

		// Expected: should be reverse-complemented due to FLIP_ON_LOAD
		String expectedRC=reverseComplement(originalSeq);
		System.out.println("Expected (reverse-complemented):");
		System.out.println(expectedRC);
		System.out.println("Matches stored? "+storedSeq.equals(expectedRC));
		System.out.println();

		// Convert to BAM
		String[] refNames={"chr1"};
		SamToBamConverter toBam=new SamToBamConverter(refNames);
		byte[] bamRecord=toBam.convertAlignment(sl);

		System.out.println("BAM record size: "+bamRecord.length+" bytes");
		System.out.println();

		// Convert back to SAM
		BamToSamConverter toSam=new BamToSamConverter(refNames);
		//TODO: Bug [stream/bam/TestSeqReverse#001] - convertAlignment already returns a body
		//without block_size. The decoder also expects that body; stripping four more bytes
		//here discards refID and shifts the decoder's fields. Pass bamRecord directly in a
		//separate behavior repair. Source-confirmed; not a runtime reproduction (STR-183).
		// Skip first 4 bytes (block_size)
		byte[] bamAlignment=new byte[bamRecord.length-4];
		System.arraycopy(bamRecord, 4, bamAlignment, 0, bamAlignment.length);
		byte[] samBytes=toSam.convertAlignment(bamAlignment);
		String roundtripSam=new String(samBytes);

		System.out.println("Roundtrip SAM line:");
		System.out.println(roundtripSam);
		System.out.println();

		// Extract SEQ field (field 10, 0-indexed field 9)
		String[] fields=roundtripSam.split("\t");
		String roundtripSeq=fields[9];

		System.out.println("Roundtrip SEQ field: "+roundtripSeq);
		System.out.println("Original SEQ field: "+originalSeq);
		System.out.println("Match? "+roundtripSeq.equals(originalSeq));
		System.out.println();

		if(!roundtripSeq.equals(originalSeq)){
			System.out.println("BUG CONFIRMED: Sequences don't match!");
			System.out.println("Roundtrip is reverse of original? "+roundtripSeq.equals(reverse(originalSeq)));
			System.out.println("Roundtrip is RC of original? "+roundtripSeq.equals(expectedRC));
			//TODO: Possible bug [stream/bam/TestSeqReverse#97] - mismatch prints "BUG CONFIRMED" but
			//does not call System.exit(1); automated test runners will see exit code 0 and report PASS
			//even when the seq-reversal bug is present.
		}else{System.out.println("SUCCESS: Sequences match!");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reverse-complements a temporary byte copy using Vector's current dispatch.
	 * Uses the platform charset; this driver's literal input is ASCII DNA.
	 * @param seq Non-null DNA sequence
	 * @return Reverse-complemented text; the input String is unchanged
	 */
	private static String reverseComplement(String seq){
		byte[] bytes=seq.getBytes();
		Vector.reverseComplementInPlace(bytes);
		return new String(bytes);
	}

	/**
	 * Reverses the character order without complementing bases.
	 * @param seq Non-null sequence text
	 * @return Reversed text; the input String is unchanged
	 */
	private static String reverse(String seq){return new StringBuilder(seq).reverse().toString();}
}
