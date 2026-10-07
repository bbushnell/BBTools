package stream.bam;

import java.nio.charset.StandardCharsets;

import parse.LineParser1;
import simd.Vector;
import stream.SamLine;

/**
 * Standalone diagnostic printing sequence orientation through SAM/BAM conversion.
 * Constructs one fixed FLAG=83 alignment and checks a round trip through the direct
 * BAM-to-SAM text converter. Stored-sequence expectations honor SamLine.FLIP_ON_LOAD.
 * A stored-sequence or round-trip sequence mismatch throws even without assertions.
 */
public class TestSeqReverse{

	//Test harness: compares stored orientation and round-trip SEQ for one fixed valid alignment.
	/**
	 * Prints and checks the fixed alignment's stored sequence and round-trip comparison.
	 * Arguments are ignored; this method does not configure the global SAM/SIMD options.
	 * Other SAM parse options must retain the fields needed by this fixture.
	 * @param args Ignored command-line arguments
	 * @throws Exception If parsing, conversion or either sequence comparison fails
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
		final LineParser1 lp=new LineParser1('\t');
		SamLine sl=new SamLine(lp.set(samLine.getBytes(StandardCharsets.US_ASCII)));

		// Check what's stored in SamLine.seq
		String storedSeq=new String(sl.seq, StandardCharsets.US_ASCII);
		System.out.println("Sequence stored in SamLine.seq:");
		System.out.println(storedSeq);

		//The parser flips this mapped reverse-strand fixture only when FLIP_ON_LOAD is enabled.
		String expectedRC=reverseComplement(originalSeq);
		String expectedStored=SamLine.FLIP_ON_LOAD ? expectedRC : originalSeq;
		System.out.println("Expected stored sequence (FLIP_ON_LOAD="+SamLine.FLIP_ON_LOAD+"):");
		System.out.println(expectedStored);
		System.out.println("Matches stored? "+storedSeq.equals(expectedStored));
		if(!storedSeq.equals(expectedStored)){throw new IllegalStateException("Stored sequence does not match FLIP_ON_LOAD orientation");}
		System.out.println();

		// Convert to BAM
		String[] refNames={"chr1"};
		SamToBamConverter toBam=new SamToBamConverter(refNames);
		byte[] bamRecord=toBam.convertAlignment(sl);

		System.out.println("BAM record size: "+bamRecord.length+" bytes");
		System.out.println();

		// Convert back to SAM
		BamToSamConverter toSam=new BamToSamConverter(refNames);
		//Resolved #001: both converters exchange a body without block_size; retain its refID bytes.
		byte[] samBytes=toSam.convertAlignment(bamRecord);
		String roundtripSam=new String(samBytes, StandardCharsets.US_ASCII);

		System.out.println("Roundtrip SAM line:");
		System.out.println(roundtripSam);
		System.out.println();

		// Extract SEQ field (field 10, 0-indexed field 9)
		String roundtripSeq=lp.set(samBytes).parseString(9);

		System.out.println("Roundtrip SEQ field: "+roundtripSeq);
		System.out.println("Original SEQ field: "+originalSeq);
		System.out.println("Match? "+roundtripSeq.equals(originalSeq));
		System.out.println();

		if(!roundtripSeq.equals(originalSeq)){
			System.out.println("BUG CONFIRMED: Sequences don't match!");
			System.out.println("Roundtrip is reverse of original? "+roundtripSeq.equals(reverse(originalSeq)));
			System.out.println("Roundtrip is RC of original? "+roundtripSeq.equals(expectedRC));
			//Resolved #97: diagnostic mismatch must not look like successful process completion.
			throw new IllegalStateException("Round-trip sequence differs from the original SAM sequence");
		}else{System.out.println("SUCCESS: Sequences match!");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reverse-complements a temporary byte copy using Vector's current dispatch.
	 * Uses explicit ASCII for this driver's literal DNA fixture.
	 * @param seq Non-null DNA sequence
	 * @return Reverse-complemented text; the input String is unchanged
	 */
	private static String reverseComplement(String seq){
		byte[] bytes=seq.getBytes(StandardCharsets.US_ASCII);
		Vector.reverseComplementInPlace(bytes);
		return new String(bytes, StandardCharsets.US_ASCII);
	}

	/**
	 * Reverses the character order without complementing bases.
	 * @param seq Non-null sequence text
	 * @return Reversed text; the input String is unchanged
	 */
	private static String reverse(String seq){return new StringBuilder(seq).reverse().toString();}
}
