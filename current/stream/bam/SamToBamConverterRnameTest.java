package stream.bam;

import java.nio.charset.StandardCharsets;

import stream.SamLine;

/** Focused regression tests for native-BAM RNAME/RNEXT dictionary resolution. */
public class SamToBamConverterRnameTest{

	/**
	 * Runs four name-resolution cases with byte-backed reference names enabled.
	 * Restores the previous global RNAME_AS_BYTES value in finally, including when a case fails.
	 * Explicit AssertionError checks do not depend on the JVM assertion switch.
	 * @param args Ignored
	 */
	public static void main(String[] args){
		final boolean oldMode=SamLine.RNAME_AS_BYTES;
		SamLine.RNAME_AS_BYTES=true;
		try{
			testLongRnameAlias();
			testLongRnextAlias();
			testDictionaryCollision();
			testMissingReferenceFails();
			System.err.println("SamToBamConverterRnameTest: PASS");
		}finally{
			SamLine.RNAME_AS_BYTES=oldMode;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Test Cases          ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires a full RNAME alias to encode as the first dictionary entry. */
	private static void testLongRnameAlias(){
		final SamToBamConverter converter=new SamToBamConverter(
				new String[]{"chr1"}, new String[]{"chr1 description"});
		final byte[] record=converter.convertAlignment(line("chr1 description", null));
		assertEquals(0, readIntLE(record, 0), "RNAME alias");
	}

	/** Requires independent full RNAME and RNEXT aliases to encode as dictionary IDs zero and one. */
	private static void testLongRnextAlias(){
		final SamToBamConverter converter=new SamToBamConverter(
				new String[]{"chr1", "chr2"},
				new String[]{"chr1 description", "chr2 another description"});
		final byte[] record=converter.convertAlignment(
				line("chr1 description", "chr2 another description"));
		assertEquals(0, readIntLE(record, 0), "RNAME alias");
		assertEquals(1, readIntLE(record, 20), "RNEXT alias");
	}

	/** Requires duplicate dictionary names to be rejected with the expected ambiguity diagnostic. */
	private static void testDictionaryCollision(){
		try{
			new SamToBamConverter(new String[]{"chr1", "chr1"});
			throw new AssertionError("Duplicate BAM dictionary names were accepted");
		}catch(IllegalArgumentException expected){
			if(!expected.getMessage().contains("Ambiguous BAM dictionary")){throw expected;}
		}
	}

	/** Requires a mapped reference absent from the dictionary to be rejected with the expected diagnostic. */
	private static void testMissingReferenceFails(){
		final SamToBamConverter converter=new SamToBamConverter(new String[]{"chr1"});
		try{
			converter.convertAlignment(line("missing", null));
			throw new AssertionError("Missing mapped RNAME was silently encoded as refID=-1");
		}catch(IllegalArgumentException expected){
			if(!expected.getMessage().contains("absent from the BAM header")){throw expected;}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a forward mapped, one-base fixture with ASCII reference names.
	 * @param rname Nonnull reference name, which may deliberately be absent from a test dictionary
	 * @param rnext Mate reference name, or null for no mate reference
	 * @return New alignment with one A base, a one-base exact-match CIGAR and quality 30
	 */
	private static SamLine line(String rname, String rnext){
		final SamLine sl=new SamLine();
		sl.qname="read1";
		sl.flag=0;
		sl.setRname(rname.getBytes(StandardCharsets.US_ASCII));
		sl.setRnext(rnext==null ? null : rnext.getBytes(StandardCharsets.US_ASCII));
		sl.pos=1;
		sl.mapq=60;
		sl.setCigar("1=");
		sl.pnext=0;
		sl.tlen=0;
		sl.setSeq(new byte[]{'A'});
		sl.setQual(new byte[]{30});
		return sl;
	}

	/**
	 * Decodes four little-endian bytes as a signed integer without advancing a cursor.
	 * @param array Source bytes containing the complete field
	 * @param offset First byte of the four-byte field
	 * @return Decoded 32-bit value
	 */
	private static int readIntLE(byte[] array, int offset){
		return (array[offset]&0xff)|((array[offset+1]&0xff)<<8)
				|((array[offset+2]&0xff)<<16)|((array[offset+3]&0xff)<<24);
	}

	/**
	 * Checks an integer result independently of the JVM assertion switch.
	 * @param expected Required value
	 * @param observed Actual value
	 * @param label Case name included in a mismatch diagnostic
	 * @throws AssertionError If the values differ
	 */
	private static void assertEquals(int expected, int observed, String label){
		if(expected!=observed){throw new AssertionError(label+": expected "+expected+", observed "+observed);}
	}
}
